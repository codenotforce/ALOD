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
    Coupled(const LodSpace& space,AdditiveKernelRieszContext& r,AdjointTestCache& cache)
        :s(space),riesz(r),aot(cache),E(s.energy().cast<Complex>()),B(s.trial()),C(s.test()),A(s.operators().system){
        ComplexSparseMatrix a=C.adjoint()*A*B,b=B.adjoint()*E*B,c=C.adjoint()*E*C;
        coarse.compute(a);trial_energy.compute(b);test_energy.compute(c);
        if(coarse.info()!=Eigen::Success||trial_energy.info()!=Eigen::Success||test_energy.info()!=Eigen::Success)
            throw std::runtime_error("regional base factorization failed");
    }
    RegionalEvaluation evaluate(const ComplexMatrix& loads,const ComplexMatrix& phi,const std::vector<int>& mask){
        if(loads.rows()!=A.rows()||loads.cols()<1||!loads.allFinite()
            ||(phi.cols()&&(phi.rows()!=A.rows()||!phi.allFinite())))
            throw std::invalid_argument("regional evaluation dimensions or values invalid");
        RegionalEvaluation out;
        ComplexMatrix f0=C.adjoint()*loads,base=coarse.solve(f0);
        if(phi.cols()==0){out.values=B*base;out.tests.resize(B.rows(),0);out.raw_tests=out.tests;}
        else {
            out.raw_tests=aot.solve(s.operators(),E*phi);out.aot_residual=aot.relative_residual();
            out.raw_base_block=(out.raw_tests.adjoint()*A*B).norm()/std::max(1.,B.norm());
            out.raw_dictionary_block=(out.raw_tests.adjoint()*A*phi-ComplexMatrix::Identity(phi.cols(),phi.cols())).norm();
            out.tests=out.raw_tests-C*test_energy.solve(C.adjoint()*E*out.raw_tests);
            // Keep the projected lower-left block: projection generally makes
            // it nonzero, even though the RAW optimal-test block vanishes.
            ComplexMatrix a0s=C.adjoint()*A*phi,as0=out.tests.adjoint()*A*B;
            ComplexMatrix response=coarse.solve(a0s),schur=out.tests.adjoint()*A*phi-as0*response;
            Eigen::FullPivLU<ComplexMatrix> lu(schur);
            if(!lu.isInvertible())throw std::runtime_error("regional coupled Schur block is singular");
            ComplexMatrix d=lu.solve(out.tests.adjoint()*loads-as0*base);
            out.values=B*(base-response*d)+phi*d;
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
            auto evaluated=model.evaluate(loads,candidate,mask);bool pass=true;
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
