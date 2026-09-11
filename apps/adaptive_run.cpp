#include "adaptive_run.hpp"
#include "fixed_support.hpp"
#include "alod/adaptive.hpp"
#include "alod/regional.hpp"
#include "alod/mesh_state.hpp"
#include <chrono>
#include <numeric>
namespace {
using namespace alod;
using namespace lod2d;
using namespace lod2d::helmholtz;
struct Input {
    fixed::Input fixed;
    int cycles=1,m=2,state_limit=0,gap=2,rank=24;
    std::size_t patch_entries=8000000,dense_entries=8000000;
    bool inherit=true,audit=false,emit=false;
    double radius=.6;
    std::vector<double> ref_theta{.3,.2};
    std::string method="ALOD";
    EllPolicy policy;
    std::set<int> force;
};
std::vector<double> numbers(const std::string& text){
    std::istringstream stream(text);std::string word;std::vector<double> out;
    while(std::getline(stream,word,',')){std::size_t n;double v=std::stod(word,&n);if(n!=word.size()||!std::isfinite(v))throw std::invalid_argument("invalid numeric list");out.push_back(v);}
    if(!text.empty()&&text.back()==',')throw std::invalid_argument("invalid trailing comma");return out;
}
Input parse(int argc,char** argv){
    Input in;std::map<std::string,std::string> opts{{"cycles","1"},{"m-ref","2"},{"state-limit","0"},{"minimum-gap","2"},
        {"rank-cap","24"},{"radius","0.6"},{"inherit","1"},{"audit","0"},{"emit-solution","0"},{"method","ALOD"},
        {"ell-mode","lazy"},{"maximum-ell","4"},{"extra-checks",""},{"force-promotions",""},{"reference-theta","0.3,0.2"},
        {"maximum-patch-entries","8000000"},{"maximum-dense-entries","8000000"}};
    std::vector<char*> forward{argv[0]};if(argc>1)forward.push_back(argv[1]);std::set<std::string> seen;
    for(int i=2;i<argc;++i){std::string arg=argv[i];auto eq=arg.find('=');auto key=arg.substr(2,eq-2);
        if(arg.starts_with("--")&&eq!=std::string::npos&&opts.contains(key)){
            if(!seen.insert(key).second)throw std::invalid_argument("duplicate adaptive option");opts[key]=arg.substr(eq+1);
        }else forward.push_back(argv[i]);
    }
    in.fixed=fixed::parse(forward.size(),forward.data());in.policy.problem=in.fixed.problem;
    auto integer=[&](std::string key,int lo,int hi){auto v=numbers(opts.at(key));if(v.size()!=1||v[0]!=std::floor(v[0])||v[0]<lo||v[0]>hi)throw std::invalid_argument("invalid "+key);return static_cast<int>(v[0]);};
    in.cycles=integer("cycles",0,1000);in.m=integer("m-ref",1,16);in.state_limit=integer("state-limit",0,10001);
    in.gap=integer("minimum-gap",0,8);in.rank=integer("rank-cap",1,24);in.inherit=integer("inherit",0,1);
    in.audit=integer("audit",0,1);in.emit=integer("emit-solution",0,1);in.policy.maximum=integer("maximum-ell",in.fixed.ell,4);
    in.patch_entries=integer("maximum-patch-entries",1,1000000000);in.dense_entries=integer("maximum-dense-entries",1,1000000000);
    auto radius=numbers(opts.at("radius"));if(radius.size()!=1||radius[0]<0)throw std::invalid_argument("invalid radius");in.radius=radius[0];
    in.ref_theta=numbers(opts.at("reference-theta"));if(in.ref_theta.size()!=static_cast<std::size_t>(in.m))throw std::invalid_argument("reference-theta length must equal m_ref");
    for(double v:in.ref_theta)if(v<=0||v>1)throw std::invalid_argument("invalid reference theta");
    in.method=opts.at("method");if(in.method!="ALOD"&&in.method!="AFEM")throw std::invalid_argument("unknown adaptive method");
    if(in.method=="AFEM"&&in.m!=1)throw std::invalid_argument("AFEM requires m_ref=1");
    auto mode=opts.at("ell-mode");if(mode!="lazy"&&mode!="every"&&mode!="fixed")throw std::invalid_argument("invalid ell mode");
    in.policy.mode=mode=="lazy"?EllMode::Lazy:mode=="every"?EllMode::EveryState:EllMode::Fixed;
    auto ids=[&](std::string key,std::set<int>& out){for(double v:numbers(opts.at(key))){if(v<0||v!=std::floor(v)||v>10000||!out.insert(v).second)throw std::invalid_argument("invalid event state IDs");}};
    ids("extra-checks",in.policy.extra_checks);ids("force-promotions",in.force);
    if(!in.force.empty()&&(in.fixed.problem!="E2"||in.method!="ALOD"||in.policy.mode==EllMode::Fixed))throw std::invalid_argument("forced promotions require the E2 adaptive diagnostic");
    bool nominal=false;for(const auto& m:in.fixed.members)nominal|=m.id==0;
    if(!nominal)throw std::invalid_argument("adaptive localization requires nominal member 0");return in;
}
}
int adaptive_main(int argc,char** argv){
 try {
    auto config=parse(argc,argv);auto& input=config.fixed;
    auto H=fixed::coarse_mesh(input);auto reference=refine_mesh_nvb(H,config.method=="AFEM"?0:input.gap);
    const auto quad=paper_quadrature(input.problem);int ell=input.ell,event_id=0,revision=0;
    int horizon=1+config.cycles*config.m,states=config.state_limit?std::min(horizon,config.state_limit):horizon;
    std::vector<FamilyMember> members;std::vector<int> training;int nominal=0;
    for(int j=0;j<static_cast<int>(input.members.size());++j){auto& m=input.members[j];if(m.id==0)nominal=j;
        members.push_back({m.id,m.role=="train"?MemberRole::Train:m.role=="test"?MemberRole::Test:MemberRole::Shift});}
    for(int id:input.ids)for(int j=0;j<static_cast<int>(members.size());++j)if(members[j].id==id)training.push_back(j);
    AdjointTestCache cache;ComplexMatrix incoming,warm_full;AdaptiveCursor cursor;bool mesh_changed=false;
    LodLimits limits;limits.maximum_reference_nodes=input.cap;limits.threads=input.threads;
    limits.maximum_patch_entries=config.patch_entries;limits.maximum_dense_entries=config.dense_entries;
    std::cout<<std::setprecision(17);
    for(int state=0;state<states;++state){
        const auto start=std::chrono::steady_clock::now();bool terminal=state+1==states;
        if(reference.mesh.nodes.size()*members.size()>config.dense_entries)
            throw std::runtime_error("adaptive load block exceeds maximum-dense-entries");
        ComplexMatrix loads(reference.mesh.nodes.size(),members.size());
        for(int j=0;j<loads.cols();++j)loads.col(j)=assemble_helmholtz_load(reference.mesh,input.members[j].problem.source,quad,input.members[j].problem.quadrature_context);
        ComplexMatrix values;Sparse E;Eigen::VectorXd eta;Eigen::MatrixXd coarse_mass,strong;
        std::unique_ptr<LodSpace> space;std::unique_ptr<AdditiveKernelRieszContext> riesz;
        RegionalResult regional;double pg=0;std::optional<double> accepted_theta;int theta_ell=-1;
        bool check=config.policy.due(state,ell,terminal)||config.force.contains(state);
        if(config.method=="AFEM"){
            auto ops=assemble_helmholtz_operators(H,16);E=ops.stiffness+256*ops.mass;
            values=solve_helmholtz_fem_batch(ops,loads,HelmholtzFemSolverKind::Umfpack);
            strong.resize(H.elems.size(),loads.cols());eta.resize(loads.cols());
            for(int j=0;j<loads.cols();++j){auto e=adaptive::diagnostics::estimate_conforming_p1_residual(H,ops,values.col(j),loads.col(j),input.members[j].problem.source,quad,input.members[j].problem.quadrature_context);
                for(int i=0;i<strong.rows();++i)strong(i,j)=e.element_squared[i];eta[j]=e.eta;
                ComplexVector residual=ops.system*values.col(j)-loads.col(j),rhs=loads.col(j);for(int n:ops.dirichlet_nodes){residual[n]=0;rhs[n]=0;}pg=std::max(pg,residual.norm()/std::max(1e-30,rhs.norm()));}
            coarse_mass=strong;
        }else for(;;){
            riesz.reset();space=std::make_unique<LodSpace>(H,reference,16,ell,input.policy=="area"?InterpolationPolicy::ManuscriptAreaWeighted:InterpolationPolicy::ArchivedArithmetic,limits);
            riesz=std::make_unique<AdditiveKernelRieszContext>(*space,input.riesz_policy=="n2"?RieszPatchPolicy::ManuscriptN2:RieszPatchPolicy::ArchivedSupportExpanded);
            E=space->energy();
            if(input.problem=="E2"){
                ComplexMatrix train_loads(loads.rows(),training.size());for(int j=0;j<train_loads.cols();++j)train_loads.col(j)=loads.col(training[j]);
                regional=train_regional(*space,*riesz,cache,train_loads,{config.radius,config.rank,config.inherit},incoming,mesh_changed);
                auto evaluated=evaluate_regional(*space,*riesz,cache,loads,regional.phi,riesz->regional_mask(config.radius));values=evaluated.values;pg=evaluated.pg_residual;
                incoming=regional.raw_kernel;mesh_changed=false;
            }else {auto solved=space->solve(loads);values=solved.values;pg=solved.pg_relative_residual;}
            auto estimate=riesz->estimate(loads,values,input.theta);eta=estimate.eta;coarse_mass=estimate.element_eta_squared;
            if(!check)break;
            LocalizationEigenConfig eig;eig.relative_tolerance=input.tolerance;eig.maximum_iterations=input.iterations;eig.dense_cross_check_max_dimension=input.dense;
            std::string warm_transport="cold";
            if(warm_full.rows()==static_cast<int>(H.nodes.size())&&warm_full.cols()){
                ComplexMatrix block(space->coarse_nodes().size(),warm_full.cols());
                for(int j=0;j<block.rows();++j)block.row(j)=warm_full.row(space->coarse_nodes()[j]);
                if(block.norm()>0){eig.warm_start=LocalizationWarmStart{space->identity()+":"+riesz->patch_policy_name()+":adjoint",block};warm_transport="coarse_nodal_injection";}
            }
            auto localization=localization_theta(*space,*riesz,eig);accepted_theta=localization.theta;theta_ell=ell;
            warm_full=ComplexMatrix::Zero(H.nodes.size(),localization.warm_start.block.cols());
            for(int j=0;j<static_cast<int>(space->coarse_nodes().size());++j)warm_full.row(space->coarse_nodes()[j])=localization.warm_start.block.row(j);
            auto ratio=localization_ratio(localization.theta,eta[nominal]);std::string action=config.policy.decision(localization.theta,eta[nominal],ell);
            if(config.force.contains(state)&&ell<config.policy.maximum)action="forced_promotion";
            bool promote=action=="promote"||action=="forced_promotion";
            std::cout<<"{\"kind\":\"ell_check\",\"ell_event_id\":"<<event_id++<<",\"state_id\":"<<state<<",\"ell\":"<<ell
                <<",\"next_ell\":"<<ell+promote<<",\"action\":\""<<action<<"\",\"theta\":"<<localization.theta<<",\"eta\":"<<eta[nominal]
                <<",\"ratio_status\":\""<<ratio.status<<"\",\"ratio\":";if(ratio.value)std::cout<<*ratio.value;else std::cout<<"null";
            std::cout<<",\"mesh_revision\":"<<revision<<",\"ritz_iterations\":"<<localization.spectrum.iterations<<",\"ritz_residual\":"<<localization.spectrum.relative_residual
                <<",\"warm_transport\":\""<<warm_transport<<"\",\"identity\":\""<<localization.identity<<"\"}\n";
            config.policy.last_check=state;
            if(!promote)break;
            ++ell;
            if(input.problem=="E1")check=false;
        }
        if(config.method!="AFEM"){
            strong.resize(reference.mesh.elems.size(),loads.cols());
            for(int j=0;j<loads.cols();++j){auto e=reference_strong_residual(*space,values.col(j),loads.col(j),input.members[j].problem.source,quad,input.members[j].problem.quadrature_context);
                for(int i=0;i<strong.rows();++i)strong(i,j)=e.fine.element_squared[i];}
        }
        Eigen::VectorXd energy(loads.cols());for(int j=0;j<loads.cols();++j)energy[j]=std::sqrt(std::max(0.,values.col(j).dot(E.cast<Complex>()*values.col(j)).real()));
        auto cm=mark_family(coarse_mass,energy,members,input.ids,input.theta);
        auto next=cursor;next.advance(config.m);
        auto fm=mark_family(strong,energy,members,input.ids,config.ref_theta[next.reference_sweep-1]);
        std::cout<<"{\"kind\":\"accepted\",\"state_id\":"<<state<<",\"coarse_cycle\":"<<cursor.coarse_cycle<<",\"reference_sweep\":"<<cursor.reference_sweep
            <<",\"cycle_complete\":"<<(cursor.cycle_complete(config.m)?"true":"false")<<",\"mesh_revision\":"<<revision<<",\"ell\":";
        if(config.method=="AFEM")std::cout<<"null";else std::cout<<ell;
        std::cout<<",\"ell_last_check\":"<<config.policy.last_check<<",\"theta\":";if(accepted_theta&&theta_ell==ell)std::cout<<*accepted_theta;else std::cout<<"null";
        std::cout<<",\"coarse_free\":"<<H.nodes.size()-H.dirichlet.size()<<",\"reference_free\":"<<reference.mesh.nodes.size()-reference.mesh.dirichlet.size()
            <<",\"coarse_nodes\":"<<H.nodes.size()<<",\"reference_nodes\":"<<reference.mesh.nodes.size()<<",\"coarse_fingerprint\":\""<<mesh_fingerprint(H)
            <<"\",\"reference_fingerprint\":\""<<mesh_fingerprint(reference.mesh)<<"\",\"pg_residual\":"<<pg<<",\"eta\":";fixed::vector(eta);
        std::cout<<",\"training_ids\":";fixed::vector(input.ids);std::cout<<",\"frozen_scales\":";fixed::vector(cm.frozen_scales);
        std::cout<<",\"coarse_marks\":";fixed::vector(cm.marked_elements);std::cout<<",\"reference_marks\":";fixed::vector(fm.marked_elements);
        std::cout<<",\"worst_member_id\":"<<cm.worst_member_id<<",\"rank\":"<<regional.phi.cols()<<",\"working_rank\":"<<regional.working_rank
            <<",\"inherited_rank\":"<<regional.inherited_rank<<",\"training_stop\":\""<<regional.stop<<"\",\"training_met\":"<<(regional.training_met?"true":"false")
            <<",\"compression_accepted\":"<<(regional.compression_accepted?"true":"false")<<",\"selected_patches\":"<<regional.selected_patches
            <<",\"training_evaluations\":"<<regional.evaluations<<",\"compression_trials\":"<<regional.compression_trials
            <<",\"kernel_residual\":"<<regional.kernel_residual<<",\"base_orthogonality\":"<<regional.base_orthogonality<<",\"gram_residual\":"<<regional.gram_residual
            <<",\"aot_factorizations\":"<<cache.factorizations()<<",\"aot_solved_columns\":"<<cache.solved_columns()<<",\"aot_identity\":\""<<cache.identity()<<"\",\"targets\":";fixed::vector(regional.targets);
        std::cout<<",\"next_phase\":\""<<(next.reference_sweep==1?"coarse_and_reference":"reference_only")<<"\",\"complete_horizon\":"<<(state+1==horizon?"true":"false");
        if(config.emit){std::cout<<",\"solution\":";fixed::complex_matrix(values);}
        if(config.audit){
            std::cout<<",\"audit\":[";
            for(int j=0;j<loads.cols();++j){auto& p=input.members[j].problem;
                auto error=compute_helmholtz_error(reference.mesh,values.col(j),16,p.exact,p.exact_gradient,quad,p.quadrature_context);
                auto exact=compute_helmholtz_error(reference.mesh,ComplexVector::Zero(values.rows()),16,p.exact,p.exact_gradient,quad,p.quadrature_context);
                if(j)std::cout<<',';std::cout<<"{\"sample\":"<<members[j].id<<",\"role\":\""<<input.members[j].role<<"\",\"energy\":"<<energy[j]<<",\"exact_norm\":"<<exact.energy<<",\"exact_error\":"<<error.energy<<'}';
            }std::cout<<']';
        }
        std::cout<<",\"wall_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n"<<std::flush;
        if(terminal)break;
        if(config.method=="AFEM"){
            auto refined=bisect_newest_vertex(H,cm.marked_elements);H=std::move(refined.mesh);
            if(H.nodes.size()>static_cast<std::size_t>(input.cap))throw std::runtime_error("AFEM node limit exceeded");
            reference=refine_mesh_nvb(H,0);
        }else {
            bool commit_coarse=next.reference_sweep==1;
            auto transition=refine_pair(H,reference,cm.marked_elements,fm.marked_elements,config.gap,ell+1,input.cap,commit_coarse);
            if(incoming.cols())incoming=(transition.injection.cast<Complex>()*incoming).eval();
            if(warm_full.cols()&&commit_coarse)warm_full=(build_nested_mesh_embedding(H,transition.coarse).P_node.cast<Complex>()*warm_full).eval();
            H=std::move(transition.coarse);reference=std::move(transition.reference);mesh_changed=true;
        }
        ++revision;cursor=next;
    }
    return 0;
 }catch(const std::exception& e){std::cerr<<"alod_run adaptive: "<<e.what()<<'\n';return 1;}
}
