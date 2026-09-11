#include "alod/regional.hpp"
#include "../lod/fingerprint.hpp"
#include <Eigen/UmfPackSupport>
#include <numeric>
#include <stdexcept>
namespace alod {
struct AdjointTestCache::Impl {
    std::string key;
    std::vector<int> free;
    ComplexSparseMatrix reduced;
    Eigen::UmfPackLU<ComplexSparseMatrix> factor;
    std::size_t factors=0,columns=0;
    double residual=0;
};
AdjointTestCache::AdjointTestCache():impl_(std::make_unique<Impl>()){}
AdjointTestCache::~AdjointTestCache()=default;
ComplexMatrix AdjointTestCache::solve(const lod2d::helmholtz::HelmholtzOperators& ops,const ComplexMatrix& rhs) {
    auto& p=*impl_;const auto& A=ops.system;int n=A.rows();
    if(n!=A.cols()||rhs.rows()!=n||!rhs.allFinite())throw std::invalid_argument("AOT dimensions or RHS invalid");
    std::vector<int> map(n,-1),free;std::vector<bool> boundary(n,false);
    for(int i:ops.dirichlet_nodes){if(i<0||i>=n)throw std::invalid_argument("invalid AOT Dirichlet DOF");boundary[i]=true;}
    for(int i=0;i<n;++i)if(!boundary[i]){map[i]=free.size();free.push_back(i);}
    FingerprintBuilder h;h.add_i64(n);h.add_double(ops.wavenumber);h.add_double(ops.boundary_beta);
    h.add_u64(ops.diffusion.size());for(double x:ops.diffusion)h.add_double(x);
    h.add_u64(ops.refractive_index.size());for(double x:ops.refractive_index)h.add_double(x);
    h.add_u64(free.size());for(int i:free)h.add_i64(i);
    h.add_u64(A.nonZeros());
    for(int c=0;c<A.outerSize();++c)for(ComplexSparseMatrix::InnerIterator it(A,c);it;++it){h.add_i64(it.row());h.add_i64(it.col());h.add_complex(it.value());}
    auto key=h.finish();
    if(p.key!=key){
        p.key.clear();std::vector<Eigen::Triplet<Complex>> entries;
        for(int c=0;c<A.outerSize();++c)for(ComplexSparseMatrix::InnerIterator it(A,c);it;++it)
            if(map[it.row()]>=0&&map[c]>=0)entries.emplace_back(map[c],map[it.row()],std::conj(it.value()));
        p.reduced.resize(free.size(),free.size());p.reduced.setFromTriplets(entries.begin(),entries.end());p.reduced.makeCompressed();
        p.factor.compute(p.reduced);if(p.factor.info()!=Eigen::Success)throw std::runtime_error("AOT factorization failed");
        p.free=free;p.key=key;++p.factors;
    }
    ComplexMatrix b(free.size(),rhs.cols());for(int j=0;j<static_cast<int>(free.size());++j)b.row(j)=rhs.row(free[j]);
    ComplexMatrix x=p.factor.solve(b);p.columns+=rhs.cols();
    if(p.factor.info()!=Eigen::Success||!x.allFinite())throw std::runtime_error("AOT solve failed");
    p.residual=(p.reduced*x-b).norm()/std::max(1e-30,b.norm());
    if(p.residual>1e-10)throw std::runtime_error("AOT residual gate failed");
    ComplexMatrix out=ComplexMatrix::Zero(n,rhs.cols());for(int j=0;j<static_cast<int>(free.size());++j)out.row(free[j])=x.row(j);
    return out;
}
std::string AdjointTestCache::identity()const{return impl_->key;}
std::size_t AdjointTestCache::factorizations()const{return impl_->factors;}
std::size_t AdjointTestCache::solved_columns()const{return impl_->columns;}
double AdjointTestCache::relative_residual()const{return impl_->residual;}
}
