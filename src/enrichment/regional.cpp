#include "alod/regional.hpp"
#include <Eigen/SparseLU>
#include <Eigen/SVD>
#include <stdexcept>
#include <cmath>
namespace alod {
namespace {
double norm(const Sparse& E,const ComplexVector& x){return std::sqrt(std::max(0.,x.dot(E.cast<Complex>()*x).real()));}
void append(ComplexMatrix& a,const ComplexVector& v){auto k=a.cols();a.conservativeResize(Eigen::NoChange,k+1);a.col(k)=v;}
// Reused by every growing/POD evaluation in one state.
struct Coupled {
    const LodSpace& s;AdditiveKernelRieszContext& riesz;AdjointTestCache& aot;
    ComplexSparseMatrix E,B,C,A;
    Eigen::SparseLU<ComplexSparseMatrix> coarse,trial_energy,test_energy;
    bool reuse;
    ComplexMatrix last_loads,last_f0,last_base,root_phi,root_ephi,root_aphi,root_raw;
    ComplexMatrix prepared_phi,prepared_raw,prepared_tests,prepared_response,prepared_as0;
    Eigen::FullPivLU<ComplexMatrix> prepared_lu;
    double prepared_residual=0;
    static bool same(const ComplexMatrix& a,const ComplexMatrix& b){return a.rows()==b.rows()&&a.cols()==b.cols()&&(a.array()==b.array()).all();}
    Coupled(const LodSpace& space,AdditiveKernelRieszContext& r,AdjointTestCache& cache,bool reuse_=true)
        :s(space),riesz(r),aot(cache),E(s.energy().cast<Complex>()),B(s.trial()),C(s.test()),A(s.operators().system),reuse(reuse_){
        ComplexSparseMatrix a=C.adjoint()*A*B,b=B.adjoint()*E*B,c=C.adjoint()*E*C;
        coarse.compute(a);trial_energy.compute(b);test_energy.compute(c);
        if(coarse.info()!=Eigen::Success||trial_energy.info()!=Eigen::Success||test_energy.info()!=Eigen::Success)
            throw std::runtime_error("regional base factorization failed");
    }
    RegionalEvaluation evaluate(const ComplexMatrix& loads,const ComplexMatrix& phi,const std::vector<int>& mask,const ComplexMatrix& transform={}){
        if(loads.rows()!=A.rows()||loads.cols()<1||!loads.allFinite()
            ||(phi.cols()&&(phi.rows()!=A.rows()||!phi.allFinite())))
            throw std::invalid_argument("regional evaluation dimensions or values invalid");
        RegionalEvaluation out;
        ComplexMatrix f0,base;
        if(reuse&&same(loads,last_loads)){f0=last_f0;base=last_base;}
        else {f0=C.adjoint()*loads;base=coarse.solve(f0);if(reuse){last_loads=loads;last_f0=f0;last_base=base;}}
        if(phi.cols()==0){out.values=B*base;out.tests.resize(B.rows(),0);out.raw_tests=out.tests;}
        else {
            if(!(reuse&&same(phi,prepared_phi))){
                ComplexMatrix ephi,aphi,raw;
                if(reuse&&transform.size()&&root_phi.cols()==transform.rows()
                    &&(root_phi*transform-phi).norm()<1e-13*std::max(1.,phi.norm())){
                    ephi=root_ephi*transform;aphi=root_aphi*transform;raw=root_raw*transform;
                }else {
                    int prefix=0;
                    if(reuse&&root_phi.rows()==phi.rows()&&root_phi.cols()<=phi.cols()
                        &&same(root_phi,phi.leftCols(root_phi.cols())))prefix=root_phi.cols();
                    int extra=phi.cols()-prefix;
                    ephi.resize(phi.rows(),phi.cols());aphi=ephi;raw=ephi;
                    if(prefix){ephi.leftCols(prefix)=root_ephi;aphi.leftCols(prefix)=root_aphi;raw.leftCols(prefix)=root_raw;}
                    if(extra){ephi.rightCols(extra)=E*phi.rightCols(extra);aphi.rightCols(extra)=A*phi.rightCols(extra);
                        raw.rightCols(extra)=aot.solve(s.operators(),ephi.rightCols(extra));}
                    if(reuse){root_phi=phi;root_ephi=ephi;root_aphi=aphi;root_raw=raw;}
                }
                // Every reused linear combination is checked against the actual
                // free-DOF adjoint equation, not just a cache label.
                ComplexMatrix defect=A.adjoint()*raw-ephi,rhs=ephi;
                for(int node:s.operators().dirichlet_nodes){defect.row(node).setZero();rhs.row(node).setZero();}
                prepared_residual=defect.norm()/std::max(1e-30,rhs.norm());
                if(prepared_residual>1e-10)throw std::runtime_error("cached AOT residual gate failed");
                prepared_raw=raw;
                prepared_tests=raw-C*test_energy.solve(C.adjoint()*E*raw);
                ComplexMatrix a0s=C.adjoint()*aphi;
                prepared_as0=prepared_tests.adjoint()*A*B;
                prepared_response=coarse.solve(a0s);
                prepared_lu.compute(prepared_tests.adjoint()*aphi-prepared_as0*prepared_response);
                if(!prepared_lu.isInvertible())throw std::runtime_error("regional coupled Schur block is singular");
                prepared_phi=phi;
            }
            out.raw_tests=prepared_raw;out.tests=prepared_tests;out.aot_residual=prepared_residual;
            out.raw_base_block=(out.raw_tests.adjoint()*A*B).norm()/std::max(1.,B.norm());
            out.raw_dictionary_block=(out.raw_tests.adjoint()*A*phi-ComplexMatrix::Identity(phi.cols(),phi.cols())).norm();
            // The projected lower-left block remains in the actual PG solve.
            ComplexMatrix d=prepared_lu.solve(out.tests.adjoint()*loads-prepared_as0*base);
            out.values=B*(base-prepared_response*d)+phi*d;
        }
        ComplexMatrix residual=loads-A*out.values;
        for(int j=0;j<loads.cols();++j){
            double numerator=(C.adjoint()*residual.col(j)).squaredNorm()+(out.tests.adjoint()*residual.col(j)).squaredNorm();
            double denominator=f0.col(j).squaredNorm()+(out.tests.adjoint()*loads.col(j)).squaredNorm();
            out.pg_residual=std::max(out.pg_residual,std::sqrt(numerator/std::max(1e-60,denominator)));
        }
        if(out.pg_residual>1e-9)throw std::runtime_error("regional coupled PG residual gate failed");
        auto local=riesz.apply_selected(residual,mask,false);
        if(local.selected_identity_relative_error>1e-7||local.constraint_relative_residual>1e-8)
            throw std::runtime_error("regional AS identity or kernel gate failed");
        out.eta=local.selected_eta;out.seeds=local.selected_values;return out;
    }
    bool add(ComplexVector raw,ComplexMatrix& psi,ComplexMatrix& phi){
        ComplexVector v=raw-B*trial_energy.solve(B.adjoint()*E*raw);
        for(int pass=0;pass<2;++pass)for(int j=0;j<phi.cols();++j){
            Complex alpha=phi.col(j).dot(E*v);v-=phi.col(j)*alpha;raw-=psi.col(j)*alpha;
        }
        double length=norm(s.energy(),v);if(!(length>1e-12))return false;
        append(phi,(v/length).eval());append(psi,(raw/length).eval());return true;
    }
};
}
struct RegionalEvaluator::Impl { Coupled model;Impl(const LodSpace& s,AdditiveKernelRieszContext& r,AdjointTestCache& a,bool reuse):model(s,r,a,reuse){} };
RegionalEvaluator::RegionalEvaluator(const LodSpace& s,AdditiveKernelRieszContext& r,AdjointTestCache& a,bool reuse):impl_(std::make_unique<Impl>(s,r,a,reuse)){}
RegionalEvaluator::~RegionalEvaluator()=default;
std::size_t RegionalEvaluator::dense_cache_bytes()const{
    const auto& m=impl_->model;std::size_t entries=0;
    for(const auto* matrix:{&m.last_loads,&m.last_f0,&m.last_base,&m.root_phi,&m.root_ephi,&m.root_aphi,&m.root_raw,
        &m.prepared_phi,&m.prepared_raw,&m.prepared_tests,&m.prepared_response,&m.prepared_as0})entries+=matrix->size();
    return sizeof(Complex)*(entries+(m.prepared_phi.cols()?m.prepared_lu.matrixLU().size():0));
}
RegionalEvaluation RegionalEvaluator::evaluate(const ComplexMatrix& loads,const ComplexMatrix& phi,const std::vector<int>& mask){return impl_->model.evaluate(loads,phi,mask);}
RegionalEvaluation evaluate_regional(const LodSpace& s,AdditiveKernelRieszContext& r,AdjointTestCache& aot,
    const ComplexMatrix& loads,const ComplexMatrix& phi,const std::vector<int>& mask){
    return Coupled(s,r,aot).evaluate(loads,phi,mask);
}
RegionalResult train_regional(LodSpace& s,AdditiveKernelRieszContext& r,AdjointTestCache& aot,
    const ComplexMatrix& loads,const RegionalConfig& cfg,const ComplexMatrix& incoming,bool changed){
    if(!std::isfinite(cfg.radius)||cfg.radius<0||cfg.rank_cap<1||cfg.rank_cap>24||loads.cols()<1||loads.rows()!=s.energy().rows()||!loads.allFinite())
        throw std::invalid_argument("invalid regional configuration or loads");
    if(static_cast<std::size_t>(s.energy().rows())*std::max<int>(cfg.rank_cap,loads.cols())>s.limits().maximum_dense_entries)
        throw std::runtime_error("regional dictionary exceeds maximum-dense-entries");
    const auto mask=r.regional_mask(cfg.radius);RegionalResult out;
    for(int selected:mask)out.selected_patches+=selected!=0;
    Coupled model(s,r,aot);ComplexMatrix psi(s.energy().rows(),0),phi=psi;
    auto base=model.evaluate(loads,phi,mask);out.targets=(.1*base.eta.array()).max(1e-12);
    if(cfg.inherit&&incoming.cols()){
        if(incoming.rows()!=s.energy().rows()||!incoming.allFinite())throw std::invalid_argument("incoming kernel must be injected to the current reference");
        ComplexMatrix raw=changed?r.apply_selected(model.E*incoming,mask,false).selected_values:incoming;
        double kernel=(s.interpolation().cast<Complex>()*raw).norm()/std::max(1e-30,s.interpolation().norm()*raw.norm());
        if(kernel>1e-8)throw std::runtime_error("incoming representative is outside the current kernel");
        for(int j=0;j<raw.cols()&&phi.cols()<cfg.rank_cap;++j)model.add(raw.col(j),psi,phi);
    }
    out.inherited_rank=phi.cols();double previous=0;int small_gains=0;RegionalEvaluation full;
    for(;;){
        full=phi.cols()?model.evaluate(loads,phi,mask):base;++out.evaluations;
        Eigen::Index worst=0;double score=(full.eta.array()/out.targets.array()).maxCoeff(&worst);
        out.training_met=score<=1;
        if(out.selected_patches==0){out.stop="empty_region";break;}
        if(out.training_met){out.stop="budget_met";break;}
        if(phi.cols()>=cfg.rank_cap){out.stop="rank_cap";break;}
        if(previous>0)small_gains=(previous-score)/previous<.01?small_gains+1:0;
        previous=score;
        if(small_gains>=2){out.stop="low_marginal_gain";break;}
        ComplexVector seed=full.seeds.col(worst);int added=0;
        for(int k=0;k<2&&phi.cols()<cfg.rank_cap;++k){
            if(!model.add(seed,psi,phi))break;
            ++added;
            if(k==0&&phi.cols()<cfg.rank_cap)
                seed=r.apply_selected(model.A*phi.col(phi.cols()-1),mask,false).selected_values.col(0);
        }
        if(!added){out.stop="dependent_seed";break;}
    }
    out.working_rank=phi.cols();out.full_phi=phi;out.full_values=full.values;out.full_eta=full.eta;
    out.accepted=full;
    if(phi.cols()){
        ComplexMatrix snapshots=full.values-base.values;
        for(int j=0;j<loads.cols();++j)snapshots.col(j)/=std::max(1e-12,norm(s.energy(),full.values.col(j)));
        ComplexMatrix coordinates=phi.adjoint()*model.E*snapshots;
        Eigen::JacobiSVD<ComplexMatrix> svd(coordinates,Eigen::ComputeThinU);
        for(int rank=1;rank<=svd.matrixU().cols();++rank){
            ++out.compression_trials;ComplexMatrix transform=svd.matrixU().leftCols(rank);
            ComplexMatrix candidate=phi*transform;
            auto evaluated=model.evaluate(loads,candidate,mask,transform);bool pass=true;
            for(int j=0;j<loads.cols();++j)
                if(norm(s.energy(),(full.values.col(j)-evaluated.values.col(j)).eval())>.2*out.targets[j]
                    ||evaluated.eta[j]>std::max(1.1*out.targets[j],1.03*full.eta[j]))pass=false;
            if(pass){out.compression_accepted=true;out.accepted=std::move(evaluated);psi=(psi*transform).eval();phi=std::move(candidate);break;}
        }
    }
    out.phi=phi;out.raw_kernel=psi;
    std::vector<bool> support(s.energy().rows(),false);
    for(int z=0;z<static_cast<int>(mask.size());++z)if(mask[z])for(int dof:r.patches()[z].discrete_dofs)support[dof]=true;
    double outside=0;for(int i=0;i<psi.rows();++i)if(!support[i])outside+=psi.row(i).squaredNorm();
    if(std::sqrt(outside)>1e-8*std::max(1e-30,psi.norm()))throw std::runtime_error("kernel representative escaped selected patch support");
    out.kernel_residual=(s.interpolation().cast<Complex>()*psi).norm()/std::max(1e-30,s.interpolation().norm()*psi.norm());
    out.base_orthogonality=(model.B.adjoint()*model.E*phi).norm()/std::max(1.,model.B.norm());
    out.gram_residual=(phi.adjoint()*model.E*phi-ComplexMatrix::Identity(phi.cols(),phi.cols())).norm();
    if(out.kernel_residual>1e-8||out.base_orthogonality>1e-8||out.gram_residual>1e-8)throw std::runtime_error("regional dictionary invariant failed");
    return out;
}
}
