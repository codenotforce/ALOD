// Offline comparison on the same reduced reference FEM operator and RHS family.
// No LOD reconstruction, reference-grid refinement, or production output writes.
#include "alod/batch.hpp"
#include "alod/checkpoint.hpp"
#include "alod/timing.hpp"
#include "mesh/refine.h"
#include <Eigen/CholmodSupport>
#include <Eigen/UmfPackSupport>
#include <unsupported/Eigen/IterativeSolvers>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <regex>

using namespace alod;
using Clock=std::chrono::steady_clock;
using RowMatrix=Eigen::SparseMatrix<Complex,Eigen::RowMajor>;
using LongMatrix=Eigen::SparseMatrix<Complex,Eigen::ColMajor,SuiteSparse_long>;
double elapsed(Clock::time_point start){return std::chrono::duration<double>(Clock::now()-start).count();}

// A strong SPD preconditioner: exact Cholesky of K+k^2*M, applied to real and
// imaginary parts. This is an energy preconditioner, not shifted-complex AMG.
struct EnergyPreconditioner {
    const Sparse* energy=nullptr;
    Eigen::CholmodSupernodalLLT<Sparse> factor;
    template<class Matrix> EnergyPreconditioner& compute(const Matrix&){
        factor.cholmod().nmethods=1;
        factor.cholmod().method[0].ordering=CHOLMOD_METIS;
        factor.compute(*energy);return *this;
    }
    Eigen::ComputationInfo info()const{return factor.info();}
    Eigen::Index rows()const{return energy->rows();}
    Eigen::Index cols()const{return energy->cols();}
    template<class Rhs> ComplexVector solve(const Rhs& b)const{
        ComplexVector result(b.rows());
        result.real()=factor.solve(b.real());result.imag()=factor.solve(b.imag());
        return result;
    }
};

template<class Matrix> Matrix restrict_free(const Matrix& a,const std::vector<int>& index,int n){
    std::vector<Eigen::Triplet<typename Matrix::Scalar>> entries;
    entries.reserve(a.nonZeros());
    for(int c=0;c<a.outerSize();++c)for(typename Matrix::InnerIterator it(a,c);it;++it)
        if(index[it.row()]>=0&&index[it.col()]>=0)
            entries.emplace_back(index[it.row()],index[it.col()],it.value());
    Matrix result(n,n);result.setFromTriplets(entries.begin(),entries.end());return result;
}

int main(int argc,char** argv){try{
    if(argc!=5&&argc!=6)throw std::invalid_argument(
        "usage: benchmark_reference_solvers checkpoint|smoke threads rhs_count direct|ilut|energy|all [internal_tolerance]");
    const int threads=std::stoi(argv[2]), count=std::stoi(argv[3]);
    const std::string mode=argv[4];const double k=16.;
    const double internal_tolerance=argc==6?std::stod(argv[5]):1e-12;
    if(threads<1||count<1||count>48||(mode!="direct"&&mode!="ilut"&&mode!="energy"&&mode!="all"))
        throw std::invalid_argument("invalid benchmark arguments");
    if(!std::isfinite(internal_tolerance)||internal_tolerance<=0||internal_tolerance>=1)
        throw std::invalid_argument("invalid internal tolerance");
    omp_set_dynamic(0);omp_set_num_threads(threads);PhaseTimer::probe(threads);
    std::cout<<std::setprecision(17)<<std::unitbuf;
    auto start=Clock::now();
    lod2d::TriMesh mesh;std::string members;
    if(std::string(argv[1])=="smoke"){
        mesh=lod2d::refine_mesh_nvb(make_problem("E1",k).initial_mesh,8).mesh;
        members="0 train 0.5 0.5 0 0 0\n";
    }else{
        auto snapshot=load_checkpoint(argv[1]);
        // The explicit scope prevents accidental interpretation of an E2 or
        // another-k checkpoint with E1's parameter factory.
        if(!std::regex_search(snapshot.config_json,std::regex("\"problem\"\\s*:\\s*\"E1\"")))
            throw std::invalid_argument("benchmark currently supports E1 k=16 only");
        std::smatch wavenumber;
        if(!std::regex_search(snapshot.config_json,wavenumber,std::regex("\"wavenumber\"\\s*:\\s*([-+0-9.eE]+)"))
           ||std::stod(wavenumber[1])!=k)
            throw std::invalid_argument("benchmark requires k=16");
        mesh=std::move(snapshot.fine.mesh);members=std::move(snapshot.members_text);
    }
    std::vector<Problem> problems;std::istringstream table(members);std::string line;
    while(std::getline(table,line)&&int(problems.size())<count){
        std::istringstream row(line);int id;std::string role;double x,y,c,a,phase;
        if(!(row>>id>>role>>x>>y>>c>>a>>phase)||id!=int(problems.size()))
            throw std::invalid_argument("expected ordered E1 member table starting at nominal 0");
        problems.push_back(lod2d::helmholtz::benchmarks::make_shifted_r1_paper_case(k,{x,y}));
    }
    if(int(problems.size())!=count)throw std::invalid_argument("missing requested RHS members");
    ComplexSparseMatrix matrix;Sparse energy;ComplexMatrix rhs;
    {
        auto ops=lod2d::helmholtz::assemble_helmholtz_operators(mesh,k);
        auto loads=assemble_load_batch(mesh,problems,paper_quadrature("E1"),threads);
        std::vector<int> index(mesh.nodes.size(),0);
        for(int i:ops.dirichlet_nodes)index[i]=-1;
        int n=0;for(int& i:index)if(i!=-1)i=n++;
        matrix=restrict_free(ops.system,index,n);
        Sparse full_energy=ops.stiffness+k*k*ops.mass;energy=restrict_free(full_energy,index,n);
        rhs.resize(n,count);for(int i=0;i<int(index.size());++i)if(index[i]>=0)rhs.row(index[i])=loads.row(i);
    }
    mesh={};
    std::cout<<"{\"kind\":\"prepared\",\"free_dofs\":"<<matrix.rows()<<",\"nonzeros\":"<<matrix.nonZeros()
        <<",\"rhs_count\":"<<count<<",\"threads\":"<<threads
        <<",\"internal_tolerance\":"<<internal_tolerance<<",\"common_prepare_seconds\":"<<elapsed(start)<<"}\n";
    ComplexMatrix reference;
    auto metrics=[&](const ComplexVector& x,int j){
        double residual=(matrix*x-rhs.col(j)).norm()/rhs.col(j).norm();
        double difference=0;
        if(reference.size()){
            ComplexVector d=x-reference.col(j);
            difference=std::sqrt(std::max(0.,std::real(d.dot(energy.cast<Complex>()*d)))/
                std::max(1e-30,std::real(reference.col(j).dot(energy.cast<Complex>()*reference.col(j)))));
        }
        return std::pair{residual,difference};
    };
    // Retain the direct solution as an independent numerical check for each
    // iterative candidate, even in single-candidate mode. Factor only once.
    {
        std::cout<<"{\"kind\":\"start\",\"method\":\"direct\"}\n";start=Clock::now();
        LongMatrix a=matrix;Eigen::UmfPackLU<LongMatrix> factor;
        factor.umfpackControl()[UMFPACK_ORDERING]=UMFPACK_ORDERING_METIS;
        factor.compute(a);if(factor.info()!=Eigen::Success)throw std::runtime_error("direct factorization failed");
        const double setup=elapsed(start);start=Clock::now();
        reference=factor.solve(rhs);const double solve=elapsed(start);
        if(factor.info()!=Eigen::Success||!reference.allFinite())throw std::runtime_error("direct solve failed");
        double residual=0;for(int j=0;j<count;++j)residual=std::max(residual,metrics(reference.col(j),j).first);
        if(residual>1e-9)throw std::runtime_error("direct true residual failed");
        std::cout<<"{\"kind\":\"result\",\"method\":\"direct\",\"setup_seconds\":"<<setup
            <<",\"solve_seconds\":"<<solve<<",\"total_seconds\":"<<setup+solve
            <<",\"maximum_true_residual\":"<<residual<<",\"converged\":true}\n";
    }
    RowMatrix row_matrix=matrix;
    const auto iterative=[&](auto& solver,const char* name){
        std::cout<<"{\"kind\":\"start\",\"method\":\""<<name<<"\"}\n";
        solver.set_restart(80);solver.setMaxIterations(400);solver.setTolerance(internal_tolerance);
        start=Clock::now();solver.compute(row_matrix);const double setup=elapsed(start);
        if(solver.info()!=Eigen::Success)throw std::runtime_error("preconditioner setup failed");
        double solve=0,max_residual=0,max_difference=0;int total_iterations=0;bool converged=true;
        for(int j=0;j<count;++j){
            start=Clock::now();ComplexVector x=solver.solve(rhs.col(j));const double seconds=elapsed(start);solve+=seconds;
            auto [residual,difference]=metrics(x,j);total_iterations+=solver.iterations();
            const bool passed=x.allFinite()&&std::isfinite(residual)&&residual<=1e-9&&difference<=1e-7;
            max_residual=std::max(max_residual,residual);max_difference=std::max(max_difference,difference);
            std::cout<<"{\"kind\":\"member\",\"method\":\""<<name<<"\",\"sample\":"<<j
                <<",\"seconds\":"<<seconds<<",\"iterations\":"<<solver.iterations()
                <<",\"true_residual\":"<<residual<<",\"relative_energy_difference\":"<<difference
                <<",\"converged\":"<<(passed?"true":"false")<<"}\n";
            if(!passed){converged=false;break;}
        }
        std::cout<<"{\"kind\":\"result\",\"method\":\""<<name<<"\",\"setup_seconds\":"<<setup
            <<",\"solve_seconds\":"<<solve<<",\"total_seconds\":"<<setup+solve
            <<",\"iterations\":"<<total_iterations<<",\"maximum_true_residual\":"<<max_residual
            <<",\"maximum_energy_difference\":"<<max_difference
            <<",\"converged\":"<<(converged?"true":"false")<<"}\n";
    };
    if(mode=="ilut"||mode=="all"){
        Eigen::GMRES<RowMatrix,Eigen::IncompleteLUT<Complex>> solver;
        solver.preconditioner().setDroptol(1e-3);solver.preconditioner().setFillfactor(10);
        iterative(solver,"gmres_ilut");
    }
    if(mode=="energy"||mode=="all"){
        Eigen::GMRES<RowMatrix,EnergyPreconditioner> solver;solver.preconditioner().energy=&energy;
        iterative(solver,"gmres_energy");
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
