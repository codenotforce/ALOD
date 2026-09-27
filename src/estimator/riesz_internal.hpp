#pragma once
#include "alod/estimator.hpp"
#include "alod/execution.hpp"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/Cholesky>
namespace alod {
// Batched Schur/saddle implementation extracted from the archived E2 production source.
struct AdditiveKernelRieszContext::Impl {
    struct Factors {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> energy_factor;
        Eigen::MatrixXd inverse_constraints;
        Eigen::LDLT<Eigen::MatrixXd> schur_factor;
        Eigen::SparseLU<ComplexSparseMatrix> saddle_factor;
        bool schur = true;
    };
    struct Group {
        KernelRieszPatch patch;
        std::vector<int> nodes;
        Eigen::SparseMatrix<double> energy;
        std::shared_ptr<Factors> factors;
    };
    std::vector<std::unique_ptr<Group>> groups;
    std::vector<KernelRieszPatch> patches;
    std::vector<std::size_t> gather_offsets;
    std::vector<std::pair<int,int>> gather_entries;
    std::shared_ptr<const LodHierarchyData> hierarchy;
    const lod2d::TriMesh& mesh;
    const ComplexSparseMatrix& system;
    const Eigen::SparseMatrix<double>& interpolation;
    explicit Impl(std::shared_ptr<const LodHierarchyData> h):hierarchy(std::move(h)),
        mesh(hierarchy->coarse),system(hierarchy->operators.system),interpolation(hierarchy->interpolation){}
    std::vector<int> dirichlet;
    int full_size = 0, coarse_size = 0, elements = 0, threads = 1;
    std::string identity,policy_name;
    std::size_t dense_limit=0;
    double preparation_seconds = 0.0;
    std::atomic<std::size_t> calls{0},columns{0};

    struct LocalSolution {ComplexMatrix values,multipliers;};
    // Both estimator and fused defect actions use this exact constrained solve.
    LocalSolution solve(int index,const ComplexMatrix& rhs) const {
        const auto& g=*groups[index];const auto& c=g.patch.constraints;
        const auto& f=*g.factors;const int n=rhs.rows(),m=c.rows();
        LocalSolution result{ComplexMatrix(n,rhs.cols()),ComplexMatrix::Zero(m,rhs.cols())};
        if(m==n){result.values.setZero();return result;}
        auto& x=result.values;auto& multipliers=result.multipliers;
        if(f.schur){
            x.real()=f.energy_factor.solve(rhs.real());
            x.imag()=f.energy_factor.solve(rhs.imag());
            if(m){
                const ComplexMatrix cr=c.cast<Complex>()*x;
                multipliers.real()=f.schur_factor.solve(cr.real());
                multipliers.imag()=f.schur_factor.solve(cr.imag());
                x-=f.inverse_constraints.cast<Complex>()*multipliers;
            }
        }else{
            ComplexMatrix b=ComplexMatrix::Zero(n+m,rhs.cols());b.topRows(n)=rhs;
            const ComplexMatrix solved=f.saddle_factor.solve(b);
            x=solved.topRows(n);multipliers=solved.bottomRows(m);
        }
        if(!x.allFinite())throw std::runtime_error("AS solve is not finite");
        return result;
    }

    template<class F> void parallel(F fn) {
        std::exception_ptr failure;
        std::mutex mutex;
        const int workers=execution_threads(threads);
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1) num_threads(workers)
#endif
        for (int i = 0; i < static_cast<int>(groups.size()); ++i) {
            try { fn(i); }
            catch (...) {
                std::lock_guard<std::mutex> lock(mutex);
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
    }
};

}
