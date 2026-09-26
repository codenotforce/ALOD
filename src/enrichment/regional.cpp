#include "alod/regional.hpp"
#include "alod/timing.hpp"
#include <Eigen/SparseLU>
#include <Eigen/SVD>
#include <stdexcept>
#include <cmath>
#include <cstdlib>
namespace alod {
namespace {
double norm(const Sparse& E,const ComplexVector& x){return std::sqrt(std::max(0.,x.dot(E.cast<Complex>()*x).real()));}
void append(ComplexMatrix& a,const ComplexVector& v){auto k=a.cols();a.conservativeResize(Eigen::NoChange,k+1);a.col(k)=v;}
// Reused by every growing/POD evaluation in one state.
struct Coupled {
    const LodSpace& s;AdditiveKernelRieszContext* riesz;AdjointTestCache& aot;
    ComplexSparseMatrix E;
    const ComplexSparseMatrix &B,&C,&A;
    Eigen::SparseLU<ComplexSparseMatrix> coarse,trial_energy,test_energy;
    bool reuse; EnrichmentTests policy;
    bool independent_factor=std::getenv("ALOD_REFERENCE_EXECUTION")!=nullptr;
    ComplexMatrix correction_coordinates;
    bool last_estimate_direct=true;
    ComplexMatrix riesz_base,riesz_phi,riesz_images;
    std::vector<int> riesz_mask;
    std::size_t training_budget=256ULL*1024*1024;
    void clear_riesz(){riesz_base.resize(0,0);riesz_phi.resize(0,0);riesz_images.resize(0,0);riesz_mask.clear();}
    ComplexMatrix last_loads,last_f0,last_base,root_phi,root_ephi,root_aphi,root_raw;
    ComplexMatrix root_as0,root_response,root_test_block;
    ComplexMatrix prepared_phi,prepared_raw,prepared_tests,prepared_response,prepared_as0;
    Eigen::FullPivLU<ComplexMatrix> prepared_lu;
    double prepared_residual=0,prepared_raw_base=0,prepared_raw_dictionary=0;
    bool trial_energy_ready=false;
    ComplexMatrix orth_phi,orth_images;
    static bool same(const ComplexMatrix& a,const ComplexMatrix& b){return a.rows()==b.rows()&&a.cols()==b.cols()&&(a.array()==b.array()).all();}
    Coupled(const LodSpace& space,AdditiveKernelRieszContext& r,AdjointTestCache& cache,bool reuse_=true,EnrichmentTests policy_=EnrichmentTests::KernelLift)
        :Coupled(space,&r,cache,reuse_,policy_){}
    Coupled(const LodSpace& space,AdditiveKernelRieszContext* r,AdjointTestCache& cache,bool reuse_,EnrichmentTests policy_)
        :s(space),riesz(r),aot(cache),E(s.energy().cast<Complex>()),B(s.trial()),C(s.test()),A(s.operators().system),reuse(reuse_),policy(policy_){
        if(const char* value=std::getenv("ALOD_TRAIN_RIESZ_BYTES")){
            std::string text(value);if(text.empty()||text.find_first_not_of("0123456789")!=std::string::npos)throw std::invalid_argument("invalid ALOD_TRAIN_RIESZ_BYTES");
            training_budget=std::stoull(text);
        }
        if(independent_factor){
            ComplexSparseMatrix a=C.adjoint()*A*B;coarse.compute(a);
        }
        if(independent_factor && coarse.info()!=Eigen::Success)
            throw std::runtime_error("regional base factorization failed");
        if(policy==EnrichmentTests::ArchivedAdjoint){
            ComplexSparseMatrix c=C.adjoint()*E*C;test_energy.compute(c);
            if(test_energy.info()!=Eigen::Success)throw std::runtime_error("regional test energy factorization failed");
        }
    }
    ComplexMatrix solve_base(const ComplexMatrix& rhs){
        if(independent_factor)return coarse.solve(rhs);
        return s.solve_reduced(rhs);
    }
    RegionalEvaluation evaluate(const ComplexMatrix& loads,const ComplexMatrix& phi,const std::vector<int>& mask,const ComplexMatrix& transform={},bool estimate=true,bool direct_riesz=false){
        if(loads.rows()!=A.rows()||loads.cols()<1||!loads.allFinite()
            ||(phi.cols()&&(phi.rows()!=A.rows()||!phi.allFinite())))
            throw std::invalid_argument("regional evaluation dimensions or values invalid");
        RegionalEvaluation out;
        if(!same(loads,last_loads)||mask!=riesz_mask)clear_riesz();
        ComplexMatrix f0,base;
        if(reuse&&same(loads,last_loads)){f0=last_f0;base=last_base;}
        else {f0=C.adjoint()*loads;base=solve_base(f0);if(reuse){last_loads=loads;last_f0=f0;last_base=base;}}
        if(phi.cols()==0){out.values=B*base;out.tests.resize(B.rows(),0);out.raw_tests=out.tests;}
        else {
            if(!(reuse&&same(phi,prepared_phi))){
                ComplexMatrix ephi,aphi,raw;int appended_from=0;bool transformed=false;
                if(reuse&&transform.size()&&root_phi.cols()==transform.rows()
                    &&(root_phi*transform-phi).norm()<1e-13*std::max(1.,phi.norm())){
                    ephi=root_ephi*transform;aphi=root_aphi*transform;raw=root_raw*transform;transformed=true;
                }else {
                    int prefix=0;
                    if(reuse&&root_phi.rows()==phi.rows()&&root_phi.cols()<=phi.cols()
                        &&same(root_phi,phi.leftCols(root_phi.cols())))prefix=root_phi.cols();
                    appended_from=prefix;int extra=phi.cols()-prefix;
                    ephi.resize(phi.rows(),phi.cols());aphi=ephi;raw=ephi;
                    if(prefix){ephi.leftCols(prefix)=root_ephi;aphi.leftCols(prefix)=root_aphi;raw.leftCols(prefix)=root_raw;}
                    if(extra){ephi.rightCols(extra)=E*phi.rightCols(extra);aphi.rightCols(extra)=A*phi.rightCols(extra);
                        if(policy==EnrichmentTests::ArchivedAdjoint)raw.rightCols(extra)=aot.solve(s.operators(),ephi.rightCols(extra));
                        else {
                            ComplexMatrix interpolated=s.interpolation().cast<Complex>()*phi.rightCols(extra);
                            ComplexMatrix coefficients(B.cols(),extra);
                            for(int i=0;i<B.cols();++i)coefficients.row(i)=interpolated.row(s.coarse_nodes()[i]);
                            raw.rightCols(extra)=phi.rightCols(extra)-B*coefficients;
                        }}
                    if(reuse){root_phi=phi;root_ephi=ephi;root_aphi=aphi;root_raw=raw;}
                }
                // Check the actual defining equation for the selected test policy.
                ComplexMatrix defect,rhs;
                if(policy==EnrichmentTests::ArchivedAdjoint){defect=A.adjoint()*raw-ephi;rhs=ephi;}
                else {defect=s.interpolation().cast<Complex>()*raw;rhs=s.interpolation().cast<Complex>()*phi;}
                if(policy==EnrichmentTests::ArchivedAdjoint)
                for(int node:s.operators().dirichlet_nodes){defect.row(node).setZero();rhs.row(node).setZero();}
                prepared_residual=defect.norm()/std::max(1e-30,policy==EnrichmentTests::ArchivedAdjoint?rhs.norm():s.interpolation().norm()*raw.norm());
                if(prepared_residual>1e-10)throw std::runtime_error("enrichment test residual gate failed");
                prepared_raw=raw;
                if(policy==EnrichmentTests::ArchivedAdjoint)prepared_tests=raw-C*test_energy.solve(C.adjoint()*E*raw);
                else prepared_tests=raw; // Y(E) = Y0 + (I - T_ell I_H) E.
                ComplexMatrix test_block;
                if(policy==EnrichmentTests::KernelLift&&reuse&&!independent_factor){
                    if(transformed&&root_as0.rows()==transform.rows()){
                        prepared_as0=transform.adjoint()*root_as0;
                        prepared_response=root_response*transform;
                        test_block=transform.adjoint()*root_test_block*transform;
                    }else{
                        if(root_as0.rows()!=appended_from||root_response.cols()!=appended_from)appended_from=0;
                        const int extra=phi.cols()-appended_from;
                        prepared_as0.resize(phi.cols(),B.cols());prepared_response.resize(B.cols(),phi.cols());
                        if(appended_from){
                            prepared_as0.topRows(appended_from)=root_as0;
                            prepared_response.leftCols(appended_from)=root_response;
                        }
                        if(extra){
                            // Only new test rows and coarse responses are needed
                            // when the dictionary grows. No rank-by-fine-grid
                            // test-action matrix is retained across evaluations.
                            ComplexMatrix new_test_A=prepared_tests.rightCols(extra).adjoint()*A;
                            prepared_as0.bottomRows(extra)=new_test_A*B;
                            prepared_response.rightCols(extra)=solve_base(C.adjoint()*aphi.rightCols(extra));
                        }
                        test_block=prepared_tests.adjoint()*aphi;
                        if(!transformed){root_as0=prepared_as0;root_response=prepared_response;root_test_block=test_block;}
                    }
                    prepared_raw_base=prepared_as0.norm()/std::max(1.,B.norm());
                    prepared_raw_dictionary=(test_block-ComplexMatrix::Identity(phi.cols(),phi.cols())).norm();
                }else{
                    ComplexMatrix test_A=prepared_tests.adjoint()*A;
                    prepared_as0=test_A*B;
                    if(policy==EnrichmentTests::KernelLift){
                        prepared_raw_base=prepared_as0.norm()/std::max(1.,B.norm());
                        prepared_raw_dictionary=(test_A*phi-ComplexMatrix::Identity(phi.cols(),phi.cols())).norm();
                    }else{
                        prepared_raw_base=(raw.adjoint()*A*B).norm()/std::max(1.,B.norm());
                        prepared_raw_dictionary=(raw.adjoint()*A*phi-ComplexMatrix::Identity(phi.cols(),phi.cols())).norm();
                    }
                    prepared_response=solve_base(C.adjoint()*aphi);
                    test_block=prepared_tests.adjoint()*aphi;
                }
                prepared_lu.compute(test_block-prepared_as0*prepared_response);
                if(!prepared_lu.isInvertible())throw std::runtime_error("regional coupled Schur block is singular");
                prepared_phi=phi;
            }
            out.raw_tests=prepared_raw;out.tests=prepared_tests;if(policy==EnrichmentTests::ArchivedAdjoint)out.aot_residual=prepared_residual;
            else out.kernel_lift_residual=prepared_residual;
            out.raw_base_block=prepared_raw_base;
            out.raw_dictionary_block=prepared_raw_dictionary;
            // The projected lower-left block remains in the actual PG solve.
            ComplexMatrix d=prepared_lu.solve(out.tests.adjoint()*loads-prepared_as0*base);
            correction_coordinates=d;
            out.values=B*(base-prepared_response*d)+phi*d;
        }
        ComplexMatrix residual=loads-A*out.values;
        for(int j=0;j<loads.cols();++j){
            double numerator=(C.adjoint()*residual.col(j)).squaredNorm()+(out.tests.adjoint()*residual.col(j)).squaredNorm();
            double denominator=f0.col(j).squaredNorm()+(out.tests.adjoint()*loads.col(j)).squaredNorm();
            out.pg_residual=std::max(out.pg_residual,std::sqrt(numerator/std::max(1e-60,denominator)));
        }
        if(out.pg_residual>1e-9)throw std::runtime_error("regional coupled PG residual gate failed");
        if(!estimate)return out;
        if(!riesz)throw std::logic_error("regional estimation requires a Riesz context");
        last_estimate_direct=true;
        const bool cache_riesz=reuse&&!independent_factor&&loads.cols()>2&&!transform.size()
            &&static_cast<long double>(A.rows())*(loads.cols()+phi.cols())*3*sizeof(Complex)<=training_budget;
        if(!cache_riesz)clear_riesz();
        if(cache_riesz&&!direct_riesz&&phi.cols()&&riesz_base.cols()==loads.cols()){
            const bool prefix=riesz_phi.cols()<=phi.cols()
                &&(!riesz_phi.cols()||same(riesz_phi,phi.leftCols(riesz_phi.cols())));
            if(prefix){
                PhaseTimer timer("training_riesz_coordinates",-1);
                const int first=riesz_phi.cols();riesz_images.conservativeResize(A.rows(),phi.cols());
                for(int col=first;col<phi.cols();col+=loads.cols()){
                    int count=std::min<int>(loads.cols(),phi.cols()-col);
                    ComplexMatrix rhs=A*(phi.middleCols(col,count)-B*prepared_response.middleCols(col,count));
                    auto image=riesz->apply_selected(rhs,mask,false);
                    if(image.selected_identity_relative_error>1e-7||image.constraint_relative_residual>1e-8)
                        throw std::runtime_error("training residual image gate failed");
                    riesz_images.middleCols(col,count)=image.selected_values;
                }
                riesz_phi=phi;
                out.seeds=riesz_base-riesz_images*correction_coordinates;
                out.eta.resize(loads.cols());bool cancellation=false;
                for(int j=0;j<loads.cols();++j){
                    const Complex action=residual.col(j).dot(out.seeds.col(j));
                    const double magnitude=residual.col(j).norm()*out.seeds.col(j).norm();
                    cancellation|=!std::isfinite(action.real())||!std::isfinite(action.imag())
                        ||action.real()<=1e-6*magnitude||std::abs(action.imag())>1e-8*std::max(1e-30,action.real());
                    out.eta[j]=std::sqrt(std::max(0.,action.real()));
                }
                if(!cancellation){last_estimate_direct=false;return out;}
            }else clear_riesz();
        }
        auto local=riesz->apply_selected(residual,mask,false);
        if(local.selected_identity_relative_error>1e-7||local.constraint_relative_residual>1e-8)
            throw std::runtime_error("regional AS identity or kernel gate failed");
        out.eta=local.selected_eta;out.seeds=local.selected_values;
        if(cache_riesz&&!phi.cols()){riesz_base=out.seeds;riesz_mask=mask;}
        return out;
    }
    bool add(ComplexVector raw,ComplexMatrix& psi,ComplexMatrix& phi){
        // Frozen evaluations and empty training regions never need this factor.
        if(!trial_energy_ready){
            ComplexSparseMatrix b=B.adjoint()*E*B;trial_energy.compute(b);
            if(trial_energy.info()!=Eigen::Success)throw std::runtime_error("regional trial energy factorization failed");
            trial_energy_ready=true;
        }
        ComplexVector v=raw-B*trial_energy.solve(B.adjoint()*E*raw);
        const auto original_v=v,original_raw=raw;
        const bool incremental=std::getenv("ALOD_REFERENCE_EXECUTION")==nullptr;
        if(incremental){
            int prefix=orth_phi.rows()==phi.rows()&&orth_phi.cols()<=phi.cols()
                &&same(orth_phi,phi.leftCols(orth_phi.cols()))?orth_phi.cols():0;
            orth_images.conservativeResize(phi.rows(),phi.cols());
            if(phi.cols()>prefix)orth_images.rightCols(phi.cols()-prefix)=E*phi.rightCols(phi.cols()-prefix);
            orth_phi=phi;
            ComplexVector ev=E*v;
            for(int pass=0;pass<2;++pass){
                for(int j=0;j<phi.cols();++j){Complex alpha=phi.col(j).dot(ev);
                    v-=phi.col(j)*alpha;raw-=psi.col(j)*alpha;ev-=orth_images.col(j)*alpha;}
                ev=E*v; // Bound accumulated drift before the next pass and gate.
            }
            const double scale=std::max(1e-30,std::sqrt(std::max(0.,v.dot(ev).real())));
            if(!v.allFinite()||(phi.adjoint()*ev).norm()>1e-10*scale){
                v=original_v;raw=original_raw;
                for(int pass=0;pass<2;++pass)for(int j=0;j<phi.cols();++j){
                    Complex alpha=phi.col(j).dot(E*v);v-=phi.col(j)*alpha;raw-=psi.col(j)*alpha;}
            }
        }else for(int pass=0;pass<2;++pass)for(int j=0;j<phi.cols();++j){
            Complex alpha=phi.col(j).dot(E*v);v-=phi.col(j)*alpha;raw-=psi.col(j)*alpha;
        }
        double length=norm(s.energy(),v);if(!(length>1e-12))return false;
        append(phi,(v/length).eval());append(psi,(raw/length).eval());return true;
    }
};
}
struct RegionalEvaluator::Impl { Coupled model;Impl(const LodSpace& s,AdditiveKernelRieszContext& r,AdjointTestCache& a,bool reuse,EnrichmentTests policy):model(s,r,a,reuse,policy){}
    Impl(const LodSpace& s,AdjointTestCache& a,bool reuse,EnrichmentTests policy):model(s,nullptr,a,reuse,policy){} };
RegionalEvaluator::RegionalEvaluator(const LodSpace& s,AdditiveKernelRieszContext& r,AdjointTestCache& a,bool reuse,EnrichmentTests policy):impl_(std::make_unique<Impl>(s,r,a,reuse,policy)){}
RegionalEvaluator::RegionalEvaluator(const LodSpace& s,AdjointTestCache& a,bool reuse,EnrichmentTests policy):impl_(std::make_unique<Impl>(s,a,reuse,policy)){}
RegionalEvaluation RegionalEvaluator::evaluate_frozen(const ComplexMatrix& loads,const ComplexMatrix& phi){return impl_->model.evaluate(loads,phi,{},{},false);}
RegionalEvaluator::~RegionalEvaluator()=default;
std::size_t RegionalEvaluator::dense_cache_bytes()const{
    const auto& m=impl_->model;std::size_t entries=0;
    for(const auto* matrix:{&m.last_loads,&m.last_f0,&m.last_base,&m.root_phi,&m.root_ephi,&m.root_aphi,&m.root_raw,
        &m.prepared_phi,&m.prepared_raw,&m.prepared_tests,&m.prepared_response,&m.prepared_as0,&m.orth_phi,&m.orth_images})entries+=matrix->size();
    return sizeof(Complex)*(entries+(m.prepared_phi.cols()?m.prepared_lu.matrixLU().size():0));
}
RegionalEvaluation RegionalEvaluator::evaluate(const ComplexMatrix& loads,const ComplexMatrix& phi,const std::vector<int>& mask){return impl_->model.evaluate(loads,phi,mask);}
RegionalEvaluation evaluate_regional(const LodSpace& s,AdditiveKernelRieszContext& r,AdjointTestCache& aot,
    const ComplexMatrix& loads,const ComplexMatrix& phi,const std::vector<int>& mask,EnrichmentTests policy){
    return Coupled(s,r,aot,true,policy).evaluate(loads,phi,mask);
}
RegionalResult train_regional(LodSpace& s,AdditiveKernelRieszContext& r,AdjointTestCache& aot,
    const ComplexMatrix& loads,const RegionalConfig& cfg,const ComplexMatrix& incoming,bool changed){
    if(!std::isfinite(cfg.radius)||cfg.radius<0||cfg.rank_cap<1||cfg.rank_cap>24||loads.cols()<1||loads.rows()!=s.energy().rows()||!loads.allFinite())
        throw std::invalid_argument("invalid regional configuration or loads");
    if(static_cast<std::size_t>(s.energy().rows())*std::max<int>(cfg.rank_cap,loads.cols())>s.limits().maximum_dense_entries)
        throw std::runtime_error("regional dictionary exceeds maximum-dense-entries");
    const auto mask=r.regional_mask(cfg.radius);RegionalResult out;
    for(int selected:mask)out.selected_patches+=selected!=0;
    Coupled model(s,r,aot,true,cfg.tests);ComplexMatrix psi(s.energy().rows(),0),phi=psi;
    auto base=model.evaluate(loads,phi,mask);out.targets=(.1*base.eta.array()).max(1e-12);
    if(cfg.inherit&&incoming.cols()){
        if(incoming.rows()!=s.energy().rows()||!incoming.allFinite())throw std::invalid_argument("incoming kernel must be injected to the current reference");
        ComplexMatrix raw=changed?r.apply_selected(model.E*incoming,mask,false).selected_values:incoming;
        double kernel=(s.interpolation().cast<Complex>()*raw).norm()/std::max(1e-30,s.interpolation().norm()*raw.norm());
        if(kernel>1e-8)throw std::runtime_error("incoming representative is outside the current kernel");
        for(int j=0;j<raw.cols()&&phi.cols()<cfg.rank_cap;++j)model.add(raw.col(j),psi,phi);
    }
    out.inherited_rank=phi.cols();double previous=0;int small_gains=0;RegionalEvaluation full;bool full_verified=false;
    for(;;){
        full=phi.cols()?model.evaluate(loads,phi,mask):base;++out.evaluations;
        full_verified=!phi.cols()||model.last_estimate_direct;
        Eigen::Index worst=0;double score=(full.eta.array()/out.targets.array()).maxCoeff(&worst);
        bool close_decision=std::abs(score-1.)<=1e-6;
        if(previous>0)close_decision|=std::abs((previous-score)/previous-.01)<=1e-6;
        for(int j=0;j<loads.cols();++j)if(j!=worst)
            close_decision|=std::abs(full.eta[j]/out.targets[j]-score)<=1e-8*std::max(1.,score);
        if(!full_verified&&close_decision){
            full=model.evaluate(loads,phi,mask,{},true,true);full_verified=true;
            score=(full.eta.array()/out.targets.array()).maxCoeff(&worst);
        }
        out.training_met=score<=1;
        if(out.selected_patches==0){out.stop="empty_region";break;}
        if(out.training_met||phi.cols()>=cfg.rank_cap){
            if(!full_verified)full=model.evaluate(loads,phi,mask,{},true,true);full_verified=true;
            score=(full.eta.array()/out.targets.array()).maxCoeff(&worst);out.training_met=score<=1;
            if(out.training_met){out.stop="budget_met";break;}
            if(phi.cols()>=cfg.rank_cap){out.stop="rank_cap";break;}
        }
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
    // Publish accepted estimators and seeds through the original local-energy
    // path, including low-gain/dependent-seed termination.
    if(phi.cols()&&!full_verified)full=model.evaluate(loads,phi,mask,{},true,true);
    out.training_met=(full.eta.array()/out.targets.array()).maxCoeff()<=1;
    out.working_rank=phi.cols();out.full_phi=phi;out.full_values=full.values;out.full_eta=full.eta;
    out.accepted=full;
    if(phi.cols()){
        ComplexMatrix snapshots=full.values-base.values;
        for(int j=0;j<loads.cols();++j)snapshots.col(j)/=std::max(1e-12,norm(s.energy(),full.values.col(j)));
        ComplexMatrix coordinates=phi.adjoint()*model.E*snapshots;
        Eigen::JacobiSVD<ComplexMatrix> svd(coordinates,Eigen::ComputeThinU);
        // POD candidates live in the fixed correction space W = phi-B*response.
        // Build its residual Gram once, instead of solving every local Riesz
        // problem for every candidate and every training member.
        ComplexMatrix residual_gram;
        const int members=loads.cols(), width=members+phi.cols();
        std::size_t pod_budget=256ULL*1024*1024;
        if(const char* value=std::getenv("ALOD_POD_GRAM_BYTES")){
            std::string text(value);if(text.empty()||text.find_first_not_of("0123456789")!=std::string::npos)throw std::invalid_argument("invalid ALOD_POD_GRAM_BYTES");
            pod_budget=std::stoull(text);
        }
        const bool reduced_pod=!std::getenv("ALOD_REFERENCE_EXECUTION")
            && svd.matrixU().cols()>2
            && static_cast<long double>(s.energy().rows())*width*3*sizeof(Complex)<=pod_budget
            && static_cast<std::size_t>(s.energy().rows())*width*3<=s.limits().maximum_dense_entries;
        if(reduced_pod){
            PhaseTimer timer("pod_residual_gram",-1);
            ComplexMatrix directions(loads.rows(),width);
            directions.leftCols(members)=loads-model.A*base.values;
            directions.rightCols(phi.cols())=model.A*(phi-model.B*model.prepared_response);
            for(int node:s.operators().dirichlet_nodes)directions.row(node).setZero();
            residual_gram.resize(width,width);
            if(model.riesz_base.cols()==members&&Coupled::same(model.riesz_phi,phi)){
                residual_gram.leftCols(members)=directions.adjoint()*model.riesz_base;
                residual_gram.rightCols(phi.cols())=directions.adjoint()*model.riesz_images;
            }else {
                // Bound local patch scratch using the same batch width as training.
                for(int first=0;first<width;first+=members){
                    int count=std::min(members,width-first);
                    auto image=r.apply_selected(directions.middleCols(first,count),mask,false);
                    if(image.selected_identity_relative_error>1e-7||image.constraint_relative_residual>1e-8)
                        throw std::runtime_error("POD residual basis Riesz gate failed");
                    residual_gram.middleCols(first,count)=directions.adjoint()*image.selected_values;
                }
            }
            residual_gram=(.5*(residual_gram+residual_gram.adjoint())).eval();
        }
        model.clear_riesz();
        for(int rank=1;rank<=svd.matrixU().cols();++rank){
            ++out.compression_trials;ComplexMatrix transform=svd.matrixU().leftCols(rank);
            ComplexMatrix candidate=phi*transform;
            auto evaluated=model.evaluate(loads,candidate,mask,transform,!reduced_pod);bool pass=true;
            bool direct=!reduced_pod;
            if(reduced_pod){
                ComplexMatrix coordinates=ComplexMatrix::Zero(width,members);
                coordinates.topRows(members).setIdentity();
                coordinates.bottomRows(phi.cols())=-transform*model.correction_coordinates;
                evaluated.eta.resize(members);
                bool cancellation=false;
                for(int j=0;j<members;++j){
                    const auto c=coordinates.col(j);
                    const double sq=std::real(c.dot(residual_gram*c));
                    const double magnitude=(c.cwiseAbs().transpose()*residual_gram.cwiseAbs()*c.cwiseAbs())(0,0);
                    const double threshold=std::max(1.1*out.targets[j],1.03*full.eta[j]);
                    cancellation|=!std::isfinite(sq)||sq<=1e-6*magnitude
                        ||std::abs(sq-threshold*threshold)<=1e-8*std::max(magnitude,threshold*threshold);
                    evaluated.eta[j]=std::sqrt(std::max(0.,sq));
                }
                if(cancellation){evaluated=model.evaluate(loads,candidate,mask,transform);direct=true;}
            }
            for(int j=0;j<loads.cols();++j)
                if(norm(s.energy(),(full.values.col(j)-evaluated.values.col(j)).eval())>.2*out.targets[j]
                    ||evaluated.eta[j]>std::max(1.1*out.targets[j],1.03*full.eta[j]))pass=false;
            // Verify any accepted candidate through the original Riesz path.
            if(pass&&!direct){
                evaluated=model.evaluate(loads,candidate,mask,transform);
                for(int j=0;j<members;++j)
                    if(evaluated.eta[j]>std::max(1.1*out.targets[j],1.03*full.eta[j]))pass=false;
            }
            if(pass){out.compression_accepted=true;out.accepted=std::move(evaluated);psi=(psi*transform).eval();phi=std::move(candidate);break;}
        }
    }
    out.phi=phi;out.raw_kernel=psi;
    if(cfg.tests==EnrichmentTests::KernelLift &&
       (out.accepted.tests-psi).norm()>1e-8*std::max(1e-30,psi.norm()))
        throw std::runtime_error("kernel lifting does not recover regional representatives");
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
