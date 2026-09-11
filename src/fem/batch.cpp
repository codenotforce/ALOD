#include "alod/batch.hpp"
#include <Eigen/UmfPackSupport>
#include <omp.h>
#include <mutex>
#include <chrono>
namespace alod {
namespace {
template<class F> void parallel(int n,int threads,F f){if(threads<1)throw std::invalid_argument("batch threads must be positive");if(n==1){f(0);return;}threads=std::min(n,threads);std::exception_ptr failure;std::mutex mutex;
#pragma omp parallel for schedule(static) num_threads(threads)
    for(int i=0;i<n;++i){try{f(i);}catch(...){std::lock_guard<std::mutex> lock(mutex);if(!failure)failure=std::current_exception();}}
    if(failure)std::rethrow_exception(failure);
}
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
}
ComplexMatrix assemble_load_batch(const lod2d::TriMesh& mesh,const std::vector<Problem>& problems,const lod2d::helmholtz::QuadraturePolicy& q,int threads){
    ComplexMatrix loads(mesh.nodes.size(),problems.size());parallel(problems.size(),threads,[&](int j){loads.col(j)=lod2d::helmholtz::assemble_helmholtz_load(mesh,problems[j].source,q,problems[j].quadrature_context);});return loads;
}
ErrorBatch integrate_error_batch(const lod2d::TriMesh& mesh,const Sparse& E,const ComplexMatrix& values,const std::vector<Problem>& problems,const lod2d::helmholtz::QuadraturePolicy& q,int threads){
    if(values.rows()!=mesh.nodes.size()||values.cols()!=problems.size())throw std::invalid_argument("error batch dimensions invalid");
    ErrorBatch result;int n=problems.size();result.exact_norm.resize(n);result.exact_error.resize(n);result.energy.resize(n);
    parallel(n,threads,[&](int j){auto& p=problems[j];result.exact_error[j]=lod2d::helmholtz::compute_helmholtz_error(mesh,values.col(j),p.wavenumber,p.exact,p.exact_gradient,q,p.quadrature_context).energy;
        result.exact_norm[j]=lod2d::helmholtz::compute_helmholtz_error(mesh,ComplexVector::Zero(values.rows()),p.wavenumber,p.exact,p.exact_gradient,q,p.quadrature_context).energy;
        result.energy[j]=std::sqrt(std::max(0.,values.col(j).dot(E.cast<Complex>()*values.col(j)).real()));});return result;
}
ErrorBatch integrate_audit_batch(const lod2d::TriMesh& mesh,const Sparse& E,const ComplexMatrix& values,const ComplexMatrix& reference,const std::vector<Problem>& problems,const lod2d::helmholtz::QuadraturePolicy& q,int threads,bool fused){
    if(values.rows()!=mesh.nodes.size()||values.cols()!=problems.size()||reference.rows()!=values.rows()||reference.cols()!=values.cols())throw std::invalid_argument("audit integration dimensions invalid");
    ErrorBatch result;int n=problems.size();result.exact_norm.resize(n);result.exact_error.resize(n);result.reference_error.resize(n);result.energy.resize(n);
    if(!fused){result=integrate_error_batch(mesh,E,values,problems,q,threads);result.reference_error.resize(n);parallel(n,threads,[&](int j){auto& p=problems[j];result.reference_error[j]=lod2d::helmholtz::compute_helmholtz_error(mesh,reference.col(j),p.wavenumber,p.exact,p.exact_gradient,q,p.quadrature_context).energy;});return result;}
    parallel(n,threads,[&](int j){auto& p=problems[j];std::array<double,3> energy_sum{};
        for(int element=0;element<static_cast<int>(mesh.elems.size());++element){auto tri=mesh.elems[element];auto a=mesh.nodes[tri[0]],b=mesh.nodes[tri[1]],c=mesh.nodes[tri[2]];
            double det=(b.x()-a.x())*(c.y()-a.y())-(b.y()-a.y())*(c.x()-a.x());if(std::abs(det)<=1e-15)throw std::invalid_argument("degenerate audit triangle");
            std::array<Eigen::Vector2d,3> gradients{Eigen::Vector2d(b.y()-c.y(),c.x()-b.x())/det,Eigen::Vector2d(c.y()-a.y(),a.x()-c.x())/det,Eigen::Vector2d(a.y()-b.y(),b.x()-a.x())/det};
            std::array<Eigen::Vector2cd,3> discrete_gradient{Eigen::Vector2cd::Zero(),Eigen::Vector2cd::Zero(),Eigen::Vector2cd::Zero()};
            for(int k=0;k<3;++k){discrete_gradient[0]+=values(tri[k],j)*gradients[k].cast<Complex>();discrete_gradient[1]+=reference(tri[k],j)*gradients[k].cast<Complex>();}
            std::array<double,3> l2{},gradient{};
            for(const auto& point:lod2d::helmholtz::triangle_quadrature_points(mesh,element,q,p.quadrature_context)){
                std::array<Complex,3> discrete{};for(int k=0;k<3;++k){discrete[0]+=point.barycentric[k]*values(tri[k],j);discrete[1]+=point.barycentric[k]*reference(tri[k],j);}
                const Complex exact=p.exact(point.point);const Eigen::Vector2cd exact_gradient=p.exact_gradient(point.point);
                for(int k=0;k<3;++k){l2[k]+=point.weight*std::norm(exact-discrete[k]);gradient[k]+=point.weight*(exact_gradient-discrete_gradient[k]).squaredNorm();}
            }
            for(int k=0;k<3;++k)energy_sum[k]+=gradient[k]+p.wavenumber*p.wavenumber*l2[k];
        }
        result.exact_error[j]=std::sqrt(std::max(0.,energy_sum[0]));result.reference_error[j]=std::sqrt(std::max(0.,energy_sum[1]));result.exact_norm[j]=std::sqrt(std::max(0.,energy_sum[2]));
        result.energy[j]=std::sqrt(std::max(0.,values.col(j).dot(E.cast<Complex>()*values.col(j)).real()));
    });return result;
}
struct ReferenceFemContext::Impl {
    ComplexSparseMatrix matrix;std::vector<int> free;int full_size=0;
    Eigen::UmfPackLU<ComplexSparseMatrix> factor;
    double factor_time=0,solve_time=0,residual=0;
};
ReferenceFemContext::ReferenceFemContext(const lod2d::helmholtz::HelmholtzOperators& ops):impl_(std::make_unique<Impl>()){
    auto& p=*impl_;auto start=Clock::now();p.full_size=ops.system.rows();std::vector<int> map(p.full_size,-1);std::vector<bool> d(p.full_size,false);
    for(int i:ops.dirichlet_nodes){if(i<0||i>=p.full_size)throw std::invalid_argument("FEM boundary out of range");d[i]=true;}
    for(int i=0;i<p.full_size;++i)if(!d[i]){map[i]=p.free.size();p.free.push_back(i);}
    std::vector<Eigen::Triplet<Complex>> entries;for(int c=0;c<ops.system.outerSize();++c)for(ComplexSparseMatrix::InnerIterator it(ops.system,c);it;++it)if(map[it.row()]>=0&&map[c]>=0)entries.emplace_back(map[it.row()],map[c],it.value());
    p.matrix.resize(p.free.size(),p.free.size());p.matrix.setFromTriplets(entries.begin(),entries.end());p.factor.compute(p.matrix);if(p.factor.info()!=Eigen::Success)throw std::runtime_error("reference FEM factorization failed");p.factor_time=seconds(start);
}
ReferenceFemContext::~ReferenceFemContext()=default;
ComplexMatrix ReferenceFemContext::solve(const ComplexMatrix& loads){auto& p=*impl_;if(loads.rows()!=p.full_size||loads.cols()<1||!loads.allFinite())throw std::invalid_argument("reference RHS invalid");auto start=Clock::now();ComplexMatrix rhs(p.free.size(),loads.cols());for(int j=0;j<rhs.rows();++j)rhs.row(j)=loads.row(p.free[j]);ComplexMatrix x=p.factor.solve(rhs);
    if(p.factor.info()!=Eigen::Success||!x.allFinite())throw std::runtime_error("reference FEM solve failed");p.residual=(p.matrix*x-rhs).norm()/std::max(1e-30,rhs.norm());if(p.residual>1e-9)throw std::runtime_error("reference FEM residual gate failed");ComplexMatrix full=ComplexMatrix::Zero(p.full_size,loads.cols());for(int j=0;j<x.rows();++j)full.row(p.free[j])=x.row(j);p.solve_time+=seconds(start);return full;}
double ReferenceFemContext::factor_seconds()const{return impl_->factor_time;}
double ReferenceFemContext::solve_seconds()const{return impl_->solve_time;}
double ReferenceFemContext::relative_residual()const{return impl_->residual;}
}
