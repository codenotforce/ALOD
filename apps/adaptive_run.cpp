#include "alod/patch_cache.hpp"
#include "alod/timing.hpp"
#include "adaptive_run.hpp"
#include "fixed_support.hpp"
#include "alod/adaptive.hpp"
#include "alod/regional.hpp"
#include "alod/mesh_state.hpp"
#include "alod/checkpoint.hpp"
#include "alod/batch.hpp"
#include "../src/lod/fingerprint.hpp"
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
    std::string checkpoint_dir,resume,config_file,pause_phase="accepted";
    int pause_state=-1;
    EllPolicy policy;
    EnrichmentTests tests=EnrichmentTests::KernelLift;
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
        {"enrichment-tests","kernel_lift"},{"ell-mode","lazy"},{"ell-ratio-mode","raw"},{"ell-threshold","0"},{"ell-absolute-threshold","-1"},{"maximum-ell","4"},{"extra-checks",""},{"force-promotions",""},{"reference-theta","0.3,0.2"},
        {"maximum-patch-entries","8000000"},{"maximum-dense-entries","8000000"},
        {"checkpoint-dir",""},{"resume",""},{"config-file",""},{"pause-state","-1"},{"pause-phase","accepted"}};
    std::vector<char*> forward{argv[0]};if(argc>1)forward.push_back(argv[1]);std::set<std::string> seen;
    for(int i=2;i<argc;++i){std::string arg=argv[i];auto eq=arg.find('=');auto key=arg.substr(2,eq-2);
        if(arg.starts_with("--")&&eq!=std::string::npos&&opts.contains(key)){
            if(!seen.insert(key).second)throw std::invalid_argument("duplicate adaptive option");opts[key]=arg.substr(eq+1);
        }else forward.push_back(argv[i]);
    }
    in.fixed=fixed::parse(forward.size(),forward.data());in.policy.problem=in.fixed.problem;
    auto integer=[&](std::string key,int lo,int hi){auto v=numbers(opts.at(key));if(v.size()!=1||v[0]!=std::floor(v[0])||v[0]<lo||v[0]>hi)throw std::invalid_argument("invalid "+key);return static_cast<int>(v[0]);};
    in.cycles=integer("cycles",0,1000);in.m=integer("m-ref",1,16);in.state_limit=integer("state-limit",0,10001);
    in.checkpoint_dir=opts["checkpoint-dir"];in.resume=opts["resume"];in.config_file=opts["config-file"];in.pause_state=integer("pause-state",-1,10000);in.pause_phase=opts["pause-phase"];
    if(in.pause_phase!="accepted"&&in.pause_phase!="training"&&in.pause_phase!="ell")throw std::invalid_argument("invalid checkpoint pause phase");
    if(in.pause_state>=0&&in.checkpoint_dir.empty())throw std::invalid_argument("pause requires checkpoint-dir");
    in.gap=integer("minimum-gap",0,8);in.rank=integer("rank-cap",1,24);in.inherit=integer("inherit",0,1);
    in.audit=integer("audit",0,1);in.emit=integer("emit-solution",0,1);in.policy.maximum=integer("maximum-ell",in.fixed.ell,4);
    auto entry=numbers(opts.at("maximum-patch-entries"));if(entry.size()!=1||entry[0]!=std::floor(entry[0])||entry[0]<1||entry[0]>64000000000.)throw std::invalid_argument("invalid maximum-patch-entries");in.patch_entries=static_cast<std::size_t>(entry[0]);in.dense_entries=integer("maximum-dense-entries",1,1000000000);
    auto radius=numbers(opts.at("radius"));if(radius.size()!=1||radius[0]<0)throw std::invalid_argument("invalid radius");in.radius=radius[0];
    in.ref_theta=numbers(opts.at("reference-theta"));if(in.ref_theta.size()!=static_cast<std::size_t>(in.m))throw std::invalid_argument("reference-theta length must equal m_ref");
    for(double v:in.ref_theta)if(v<=0||v>1)throw std::invalid_argument("invalid reference theta");
    in.method=opts.at("method");if(in.method!="ALOD"&&in.method!="AFEM")throw std::invalid_argument("unknown adaptive method");
    if(in.method=="AFEM"&&in.m!=1)throw std::invalid_argument("AFEM requires m_ref=1");
    auto mode=opts.at("ell-mode");if(mode!="lazy"&&mode!="every"&&mode!="fixed")throw std::invalid_argument("invalid ell mode");
    in.policy.mode=mode=="lazy"?EllMode::Lazy:mode=="every"?EllMode::EveryState:EllMode::Fixed;
    auto ratio_mode=opts.at("ell-ratio-mode");
    if(ratio_mode!="raw"&&ratio_mode!="solution_scaled")throw std::invalid_argument("invalid ell ratio mode");
    in.policy.solution_scaled=ratio_mode=="solution_scaled";
    auto tests=opts.at("enrichment-tests");
    if(tests!="kernel_lift"&&tests!="adjoint")throw std::invalid_argument("invalid enrichment tests");
    in.tests=tests=="kernel_lift"?EnrichmentTests::KernelLift:EnrichmentTests::ArchivedAdjoint;
    auto threshold=numbers(opts.at("ell-threshold"));
    if(threshold.size()!=1||threshold[0]<0)throw std::invalid_argument("invalid ell threshold");
    in.policy.ratio_threshold=threshold[0];
    auto absolute=numbers(opts.at("ell-absolute-threshold"));
    if(absolute.size()!=1||(absolute[0]<0&&absolute[0]!=-1))throw std::invalid_argument("invalid manual tau(k)");
    in.policy.absolute_threshold=absolute[0];
    auto ids=[&](std::string key,std::set<int>& out){for(double v:numbers(opts.at(key))){if(v<0||v!=std::floor(v)||v>10000||!out.insert(v).second)throw std::invalid_argument("invalid event state IDs");}};
    ids("extra-checks",in.policy.extra_checks);ids("force-promotions",in.force);
    if(!in.force.empty()&&(in.fixed.problem!="E2"||in.method!="ALOD"||in.policy.mode==EllMode::Fixed))throw std::invalid_argument("forced promotions require the E2 adaptive diagnostic");
    bool nominal=false;for(const auto& m:in.fixed.members)nominal|=m.id==0;
    if(!nominal)throw std::invalid_argument("adaptive localization requires nominal member 0");return in;
}
std::string read_text(const std::string& path){std::ifstream in(path);if(!in)throw std::runtime_error("cannot read configuration or member table");return std::string(std::istreambuf_iterator<char>(in),{});}
std::string mathematics_key(const Input& c,const std::string& members){
    FingerprintBuilder f;f.add_string("ALOD-adaptive-schema-1-p56");auto& p=c.fixed;for(const auto& s:{p.problem,p.policy,p.riesz_policy,c.method,members})f.add_string(s);
    for(int n:{p.level,p.gap,p.ell,int(p.graded),p.iterations,p.dense,c.m,c.gap,c.rank,int(c.inherit),c.policy.maximum,int(c.policy.mode)})f.add_i64(n);
    if(c.policy.absolute_threshold>=0){f.add_string("lazy-absolute-theta-v1");f.add_double(c.policy.absolute_threshold);}
    // Preserve archived keys only for the archived raw policy.
    if(c.policy.solution_scaled||c.policy.ratio_threshold!=0){f.add_string("solution-scaled-ell-v2-relative-Ritz-repeat");f.add_i64(c.policy.solution_scaled);f.add_double(c.policy.ratio_threshold);}
    if(c.tests==EnrichmentTests::KernelLift)f.add_string("kernel-lift-v4");
    // Preserve the default-kappa checkpoint key; nondefault operators get a distinct key.
    if(p.wavenumber!=16){f.add_string("wavenumber");f.add_double(p.wavenumber);}
    for(double x:{p.theta,p.tolerance,c.radius})f.add_double(x);for(double x:c.ref_theta)f.add_double(x);
    f.add_u64(p.ids.size());for(int i:p.ids)f.add_i64(i);f.add_u64(c.policy.extra_checks.size());for(int i:c.policy.extra_checks)f.add_i64(i);f.add_u64(c.force.size());for(int i:c.force)f.add_i64(i);return f.finish();
}
struct BufferedEvents {
    std::ostringstream buffer;std::streambuf* original;
    explicit BufferedEvents(const std::string& pending):original(std::cout.rdbuf(buffer.rdbuf())){std::cout<<pending;}
    ~BufferedEvents(){std::cout.rdbuf(original);}
    void publish(){std::cout.rdbuf(original);std::cout<<buffer.str()<<std::flush;}
};
}
int adaptive_main(int argc,char** argv){
 try {
    auto config=parse(argc,argv);auto& input=config.fixed;PhaseTimer::probe(input.threads);
    PhaseTimer total_timer("adaptive_process",-1);
    auto H=fixed::coarse_mesh(input);auto reference=refine_mesh_nvb(H,config.method=="AFEM"?0:input.gap);
    MeshState coarse_history(make_problem(input.problem).initial_mesh);
    coarse_history.adopt(refine_mesh_nvb(coarse_history.mesh,input.level));
    if(input.graded)coarse_history.refine({0});
    MeshState fine_history=coarse_history;fine_history.adopt(reference);
    const auto quad=paper_quadrature(input.problem);int ell=input.ell,event_id=0,revision=0;
    int horizon=1+config.cycles*config.m,states=config.state_limit?std::min(horizon,config.state_limit):horizon;
    // Production evaluates only training members plus the nominal localization
    // denominator. The complete frozen table remains in every checkpoint.
    input.members.erase(std::remove_if(input.members.begin(),input.members.end(),[&](const auto& m){
        return m.id!=0&&std::find(input.ids.begin(),input.ids.end(),m.id)==input.ids.end();
    }),input.members.end());
    std::vector<FamilyMember> members;std::vector<int> training;int nominal=0;
    for(int j=0;j<static_cast<int>(input.members.size());++j){auto& m=input.members[j];if(m.id==0)nominal=j;
        members.push_back({m.id,m.role=="train"?MemberRole::Train:m.role=="test"?MemberRole::Test:MemberRole::Shift});}
    for(int id:input.ids)for(int j=0;j<static_cast<int>(members.size());++j)if(members[j].id==id)training.push_back(j);
    AdjointTestCache cache;ComplexMatrix incoming,warm_full;AdaptiveCursor cursor;bool mesh_changed=false;
    std::string members_text=read_text(input.members_path),config_json=config.config_file.empty()?"{}":read_text(config.config_file);
    std::string key=mathematics_key(config,members_text),pending_events;std::uint64_t committed_lines=0,prefix_hash=14695981039346656037ULL;
    bool resumed_current=false,resumed_check=false;int begin_state=0;
    LodLimits limits;limits.maximum_reference_nodes=input.cap;limits.threads=input.threads;
    limits.maximum_patch_entries=config.patch_entries;limits.maximum_dense_entries=config.dense_entries;
    if(!std::getenv("ALOD_REFERENCE_EXECUTION")){
        std::size_t cap=0; // Opt in only when patch reuse offsets identity/copy costs.
        if(const char* value=std::getenv("ALOD_PATCH_CACHE_BYTES")){std::string text(value);if(text.empty()||text.find_first_not_of("0123456789")!=std::string::npos)throw std::invalid_argument("invalid ALOD_PATCH_CACHE_BYTES");cap=std::stoull(text);}
        if(cap)limits.patch_cache=std::make_shared<LodPatchCache>(cap);
    }
    auto advance=[&](const std::vector<int>& coarse_marks,const std::vector<int>& reference_marks){
        PhaseTimer timing("mesh_transition",cursor.state_id);
        auto next=cursor;next.advance(config.m);
        if(config.method=="AFEM"){
            auto refined=bisect_newest_vertex(H,coarse_marks);coarse_history.adopt(refined);H=refined.mesh;
            if(H.nodes.size()>static_cast<std::size_t>(input.cap))throw std::runtime_error("AFEM node limit exceeded");
            reference=refine_mesh_nvb(H,0);fine_history=coarse_history;
        }else {
            bool commit_coarse=next.reference_sweep==1;
            auto transition=refine_pair(H,reference,coarse_marks,reference_marks,config.gap,ell+1,input.cap,commit_coarse);
            if(incoming.cols())incoming=(transition.injection.cast<Complex>()*incoming).eval();
            if(warm_full.cols()&&commit_coarse)warm_full=(transition.coarse_injection.cast<Complex>()*warm_full).eval();
            coarse_history.adopt({transition.coarse,{},transition.coarse_parent,{}});
            fine_history.adopt({transition.reference.mesh,{},transition.reference_parent,{}});
            H=std::move(transition.coarse);reference=std::move(transition.reference);mesh_changed=true;
        }
        ++revision;cursor=next;
    };
    if(!config.resume.empty()){
        auto saved=load_checkpoint(config.resume);
        if(saved.mathematics_key!=key)throw std::invalid_argument("checkpoint mathematical configuration or members mismatch");
        H=saved.coarse.mesh;reference=saved.reference();coarse_history=std::move(saved.coarse);fine_history=std::move(saved.fine);
        incoming=saved.raw_kernel;warm_full=saved.warm_full;cursor=saved.cursor;ell=saved.ell;revision=saved.revision;event_id=saved.next_event;
        config.policy.last_check=saved.last_check;mesh_changed=saved.mesh_changed;committed_lines=saved.committed_lines;prefix_hash=saved.journal_hash;
        if(saved.phase==CheckpointPhase::Accepted){
            committed_lines+=std::count(saved.journal.begin(),saved.journal.end(),'\n');prefix_hash=journal_hash(saved.journal,prefix_hash);
            if(cursor.state_id+1>=states)return 0;
            advance(saved.coarse_marks,saved.reference_marks);
        }else {pending_events=saved.journal;resumed_current=true;resumed_check=saved.check_pending;}
        begin_state=cursor.state_id;
    }
    auto snapshot=[&](CheckpointPhase phase,const std::string& journal,bool check,const ComplexMatrix& values,const ComplexMatrix& phi,
                      const std::vector<int>& cm,const std::vector<int>& fm,const std::string& identity,const LodSpace* accepted_space=nullptr){
        if(config.checkpoint_dir.empty())return;
        PhaseTimer timing(phase==CheckpointPhase::Accepted?"checkpoint_accepted":"checkpoint_restart",cursor.state_id);
        Checkpoint s;s.coarse=coarse_history;s.fine=fine_history;s.P_node=reference.P_node;s.P_elem=reference.P_elem;s.P_dg=reference.P_dg;
        s.raw_kernel=incoming;s.phi=phi;s.values=values;s.warm_full=warm_full;s.coarse_marks=cm;s.reference_marks=fm;
        for(const auto& m:input.members)s.computed_ids.push_back(m.id);
        s.cursor=cursor;s.ell=ell;s.last_check=config.policy.last_check;s.next_event=event_id;s.revision=revision;s.phase=phase;s.mesh_changed=mesh_changed;s.check_pending=check;
        s.committed_lines=committed_lines;s.journal_hash=prefix_hash;s.journal=journal;s.mathematics_key=key;s.config_json=config_json;s.members_text=members_text;s.space_identity=identity;
        save_checkpoint(config.checkpoint_dir,s,accepted_space?&accepted_space->trial():nullptr,
            accepted_space&&!std::getenv("ALOD_REFERENCE_EXECUTION")?&accepted_space->reduced():nullptr,
            !std::getenv("ALOD_REFERENCE_EXECUTION"));
    };
    std::cout<<std::setprecision(17);
    for(int state=begin_state;state<states;++state){
        PhaseTimer state_timer("adaptive_state",state);
        BufferedEvents events(pending_events);pending_events.clear();
        const auto start=std::chrono::steady_clock::now();bool terminal=state+1==states;
        if(reference.mesh.nodes.size()*members.size()>config.dense_entries)
            throw std::runtime_error("adaptive load block exceeds maximum-dense-entries");
        std::vector<Problem> problems;for(const auto& m:input.members)problems.push_back(m.problem);
        double prepare_seconds=0,training_seconds=0,estimate_seconds=0,theta_seconds=0,strong_seconds=0;
        auto phase_start=std::chrono::steady_clock::now();
        auto phase_time=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-phase_start).count();};
        ComplexMatrix values;Sparse E;Eigen::VectorXd eta;Eigen::MatrixXd coarse_mass,strong;
        std::unique_ptr<LodSpace> space;std::unique_ptr<AdditiveKernelRieszContext> riesz;
        RegionalResult regional;double pg=0;std::optional<double> accepted_theta;int theta_ell=-1;
        bool was_resumed=resumed_current;
        bool check=resumed_current?resumed_check:(config.policy.due(state,ell,terminal)||config.force.contains(state));resumed_current=false;
        if(!was_resumed)snapshot(CheckpointPhase::BeforeTraining,events.buffer.str(),check,{},{},{},{},"");
        if(config.pause_state==state&&config.pause_phase=="training")return 0;
        auto load_start=std::chrono::steady_clock::now();
        std::vector<SourceMomentData> source_moments;
        ComplexMatrix loads;{PhaseTimer timing("load",state);loads=assemble_load_batch(reference.mesh,problems,quad,input.threads,
            std::getenv("ALOD_REFERENCE_EXECUTION")||config.method=="AFEM"?nullptr:&source_moments);}
        double load_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-load_start).count();

        if(config.method=="AFEM"){
            auto ops=assemble_helmholtz_operators(H,input.wavenumber);E=ops.stiffness+input.wavenumber*input.wavenumber*ops.mass;
            values=solve_helmholtz_fem_batch(ops,loads,HelmholtzFemSolverKind::Umfpack);
            strong.resize(H.elems.size(),loads.cols());eta.resize(loads.cols());
            for(int j=0;j<loads.cols();++j){auto e=adaptive::diagnostics::estimate_conforming_p1_residual(H,ops,values.col(j),loads.col(j),input.members[j].problem.source,quad,input.members[j].problem.quadrature_context);
                for(int i=0;i<strong.rows();++i)strong(i,j)=e.element_squared[i];eta[j]=e.eta;
                ComplexVector residual=ops.system*values.col(j)-loads.col(j),rhs=loads.col(j);for(int n:ops.dirichlet_nodes){residual[n]=0;rhs[n]=0;}pg=std::max(pg,residual.norm()/std::max(1e-30,rhs.norm()));}
            coarse_mass=strong;
        }else for(;;){
            phase_start=std::chrono::steady_clock::now();
            {PhaseTimer timing("lod_construct",state);
            if(space)space=std::make_unique<LodSpace>(*space,ell);
            else space=std::make_unique<LodSpace>(H,reference,input.wavenumber,ell,input.policy=="area"?InterpolationPolicy::ManuscriptAreaWeighted:InterpolationPolicy::ArchivedArithmetic,limits);
            }
            {PhaseTimer timing("riesz_prepare",state);if(!riesz)riesz=std::make_unique<AdditiveKernelRieszContext>(*space,input.riesz_policy=="n2"?RieszPatchPolicy::ManuscriptN2:RieszPatchPolicy::ArchivedSupportExpanded);
            }
            if(riesz->reference_identity()!=space->reference_identity())throw std::runtime_error("ell-only Riesz cache identity mismatch");
            E=space->energy();prepare_seconds+=phase_time();phase_start=std::chrono::steady_clock::now();
            {PhaseTimer timing("training_or_lod_solve",state);
            if(input.problem=="E2"){
                ComplexMatrix train_loads(loads.rows(),training.size());for(int j=0;j<train_loads.cols();++j)train_loads.col(j)=loads.col(training[j]);
                regional=train_regional(*space,*riesz,cache,train_loads,{config.radius,config.rank,config.inherit,config.tests},incoming,mesh_changed);
                values.resize(loads.rows(),loads.cols());pg=regional.accepted.pg_residual;
                for(int j=0;j<static_cast<int>(training.size());++j)values.col(training[j])=regional.accepted.values.col(j);
                if(std::find(training.begin(),training.end(),nominal)==training.end()){
                    auto evaluated=evaluate_regional(*space,*riesz,cache,loads.col(nominal),regional.phi,riesz->regional_mask(config.radius),config.tests);
                    values.col(nominal)=evaluated.values.col(0);pg=std::max(pg,evaluated.pg_residual);
                }
                incoming=std::move(regional.raw_kernel);mesh_changed=false;
                regional.full_phi.resize(0,0);regional.full_values.resize(0,0);regional.full_eta.resize(0);
                regional.accepted.values.resize(0,0);regional.accepted.seeds.resize(0,0);
                regional.accepted.tests.resize(0,0);regional.accepted.raw_tests.resize(0,0);
            }else {auto solved=space->solve(loads);values=solved.values;pg=solved.pg_relative_residual;}
            }
            training_seconds+=phase_time();phase_start=std::chrono::steady_clock::now();
            {PhaseTimer timing("riesz_estimate",state);auto estimate=riesz->estimate(loads,values,input.theta);eta=estimate.eta;coarse_mass=estimate.element_eta_squared;}estimate_seconds+=phase_time();phase_start=std::chrono::steady_clock::now();
            if(!check)break;
            LocalizationEigenConfig eig;eig.relative_tolerance=input.tolerance;eig.maximum_iterations=input.iterations;eig.dense_cross_check_max_dimension=input.dense;
            eig.eigenvalue_relative_residual=config.policy.solution_scaled;
            std::string warm_transport="cold";
            if(warm_full.rows()==static_cast<int>(H.nodes.size())&&warm_full.cols()){
                ComplexMatrix block(space->coarse_nodes().size(),warm_full.cols());
                for(int j=0;j<block.rows();++j)block.row(j)=warm_full.row(space->coarse_nodes()[j]);
                if(block.norm()>0){eig.warm_start=LocalizationWarmStart{space->identity()+":"+riesz->patch_policy_name()+":adjoint",block};warm_transport="coarse_nodal_injection";}
            }
            auto localization=[&]{PhaseTimer timing("theta",state);return localization_theta(*space,*riesz,eig);}();accepted_theta=localization.theta;theta_ell=ell;
            warm_full=ComplexMatrix::Zero(H.nodes.size(),localization.warm_start.block.cols());
            for(int j=0;j<static_cast<int>(space->coarse_nodes().size());++j)warm_full.row(space->coarse_nodes()[j])=localization.warm_start.block.row(j);
            theta_seconds+=phase_time();
            double solution_scale=std::max(1e-12,std::sqrt(std::max(0.,values.col(nominal).dot(E.cast<Complex>()*values.col(nominal)).real())));
            auto raw_ratio=localization_ratio(localization.theta,eta[nominal]);
            auto scaled_ratio=localization_ratio(solution_scale*localization.theta,eta[nominal]);
            auto ratio=config.policy.solution_scaled?scaled_ratio:raw_ratio;
            std::string action=config.policy.decision(localization.theta,eta[nominal],ell,solution_scale);
            if(config.force.contains(state)&&ell<config.policy.maximum)action="forced_promotion";
            bool promote=action=="promote"||action=="forced_promotion";
            bool absolute_exceeded=config.policy.absolute_threshold>=0&&localization.theta>config.policy.absolute_threshold;
            bool ratio_exceeded=ratio.status=="zero_denominator"||(ratio.value&&*ratio.value>config.policy.threshold());
            std::cout<<"{\"kind\":\"ell_check\",\"ell_event_id\":"<<event_id++<<",\"state_id\":"<<state<<",\"ell\":"<<ell
                <<",\"next_ell\":"<<ell+promote<<",\"action\":\""<<action<<"\",\"theta\":"<<localization.theta<<",\"eta\":"<<eta[nominal]
                <<",\"ratio_status\":\""<<ratio.status<<"\",\"ratio\":";if(ratio.value)std::cout<<*ratio.value;else std::cout<<"null";
            std::cout<<",\"absolute_threshold\":";if(config.policy.absolute_threshold>=0)std::cout<<config.policy.absolute_threshold;else std::cout<<"null";
            std::cout<<",\"absolute_exceeded\":"<<(absolute_exceeded?"true":"false")<<",\"ratio_exceeded\":"<<(ratio_exceeded?"true":"false");
            std::cout<<",\"ratio_mode\":\""<<(config.policy.solution_scaled?"solution_scaled":"raw")<<"\",\"solution_scale\":"<<solution_scale<<",\"scaled_theta\":"<<solution_scale*localization.theta<<",\"threshold\":"<<config.policy.threshold();
            std::cout<<",\"raw_ratio\":";if(raw_ratio.value)std::cout<<*raw_ratio.value;else std::cout<<"null";
            std::cout<<",\"scaled_ratio\":";if(scaled_ratio.value)std::cout<<*scaled_ratio.value;else std::cout<<"null";
            std::cout<<",\"criterion_satisfied\":"<<(!absolute_exceeded&&!ratio_exceeded?"true":"false");
            std::cout<<",\"ritz_eigenvalue_relative\":"<<(eig.eigenvalue_relative_residual?"true":"false");
            std::cout<<",\"mesh_revision\":"<<revision<<",\"ritz_iterations\":"<<localization.spectrum.iterations<<",\"ritz_residual\":"<<localization.spectrum.relative_residual
                <<",\"warm_transport\":\""<<warm_transport<<"\",\"identity\":\""<<localization.identity<<"\"}\n";
            config.policy.last_check=state;
            if(!promote)break;
            ++ell;
            if(input.problem=="E1"&&!config.policy.solution_scaled&&config.policy.absolute_threshold<0)check=false; // archived behavior only
            else check=true; // promotion never consumes the new ell check
            snapshot(CheckpointPhase::AfterEll,events.buffer.str(),check,{},{},{},{},"");
            if(config.pause_state==state&&config.pause_phase=="ell")return 0;
        }
        phase_start=std::chrono::steady_clock::now();
        {PhaseTimer timing("strong_residual",state);
        if(config.method!="AFEM"){
            strong.resize(reference.mesh.elems.size(),loads.cols());
            lod2d::helmholtz::adaptive::diagnostics::ResidualMeshContext residual_mesh(space->fine());
            for(int j=0;j<loads.cols();++j){auto e=residual_mesh.estimate(space->operators(),values.col(j),loads.col(j),input.members[j].problem.source,quad,input.members[j].problem.quadrature_context,source_moments.empty()?nullptr:&source_moments[j]);
                for(int i=0;i<strong.rows();++i)strong(i,j)=e.element_squared[i];}
        }
        }
        if(limits.patch_cache){
            PhaseTimer::counter("patch_cache_hits",limits.patch_cache->hits(),state);
            PhaseTimer::counter("patch_cache_misses",limits.patch_cache->misses(),state);
            PhaseTimer::counter("patch_cache_bytes",limits.patch_cache->bytes(),state);
        }
        PhaseTimer::counter("source_moment_members",source_moments.size(),state);
        source_moments.clear();source_moments.shrink_to_fit();
        strong_seconds=phase_time();
        Eigen::VectorXd energy(loads.cols());for(int j=0;j<loads.cols();++j)energy[j]=std::sqrt(std::max(0.,values.col(j).dot(E.cast<Complex>()*values.col(j)).real()));
        auto cm=mark_family(coarse_mass,energy,members,input.ids,input.theta);
        auto next=cursor;next.advance(config.m);
        auto fm=mark_family(strong,energy,members,input.ids,config.ref_theta[next.reference_sweep-1]);
        auto captured=[](const Eigen::VectorXd& mass,const std::vector<int>& marks){double value=0;for(int e:marks)value+=mass[e];return value;};
        std::cout<<"{\"kind\":\"accepted\",\"state_id\":"<<state<<",\"coarse_cycle\":"<<cursor.coarse_cycle<<",\"reference_sweep\":"<<cursor.reference_sweep
            <<",\"cycle_complete\":"<<(cursor.cycle_complete(config.m)?"true":"false")<<",\"mesh_revision\":"<<revision<<",\"ell\":";
        if(config.method=="AFEM")std::cout<<"null";else std::cout<<ell;
        std::cout<<",\"ell_last_check\":"<<config.policy.last_check<<",\"theta\":";if(accepted_theta&&theta_ell==ell)std::cout<<*accepted_theta;else std::cout<<"null";
        std::cout<<",\"coarse_free\":"<<H.nodes.size()-H.dirichlet.size()<<",\"reference_free\":"<<reference.mesh.nodes.size()-reference.mesh.dirichlet.size()
            <<",\"coarse_nodes\":"<<H.nodes.size()<<",\"reference_nodes\":"<<reference.mesh.nodes.size()<<",\"coarse_fingerprint\":\""<<mesh_fingerprint(H)
            <<"\",\"reference_fingerprint\":\""<<mesh_fingerprint(reference.mesh)<<"\",\"pg_residual\":"<<pg<<",\"eta\":";fixed::vector(eta);
        std::cout<<",\"patch_residual\":";if(space)std::cout<<space->patch_residual();else std::cout<<"null";
        std::cout<<",\"constraint_residual\":";if(space)std::cout<<space->constraint_residual();else std::cout<<"null";
        std::cout<<",\"dictionary_hash\":"<<json_string(matrix_hash(regional.phi))<<",\"kernel_hash\":"<<json_string(matrix_hash(incoming));
        std::cout<<",\"threads\":"<<input.threads<<",\"wavenumber\":"<<input.wavenumber<<",\"H_max\":"<<mesh_diameter(H)<<",\"h_max\":"<<mesh_diameter(reference.mesh)<<",\"coarse_elements\":"<<H.elems.size()<<",\"reference_elements\":"<<reference.mesh.elems.size();
        std::cout<<",\"maximum_budget_ratio\":";if(regional.targets.size())std::cout<<(regional.accepted.eta.array()/regional.targets.array()).maxCoeff();else std::cout<<"null";
        std::cout<<",\"coarse_mean_total\":"<<cm.mean_squared.sum()<<",\"coarse_worst_total\":"<<cm.worst_squared.sum()<<",\"coarse_mean_captured\":"<<captured(cm.mean_squared,cm.marked_elements)<<",\"coarse_worst_captured\":"<<captured(cm.worst_squared,cm.marked_elements);
        std::cout<<",\"reference_mean_total\":"<<fm.mean_squared.sum()<<",\"reference_worst_total\":"<<fm.worst_squared.sum()<<",\"reference_mean_captured\":"<<captured(fm.mean_squared,fm.marked_elements)<<",\"reference_worst_captured\":"<<captured(fm.worst_squared,fm.marked_elements);
        std::cout<<",\"computed_ids\":[";for(int j=0;j<static_cast<int>(input.members.size());++j){if(j)std::cout<<',';std::cout<<input.members[j].id;}std::cout<<']';
        std::cout<<",\"prepare_seconds\":"<<prepare_seconds<<",\"training_seconds\":"<<training_seconds<<",\"estimate_seconds\":"<<estimate_seconds<<",\"theta_seconds\":"<<theta_seconds<<",\"strong_seconds\":"<<strong_seconds;
        std::cout<<",\"coarse_mean_bulk\":"<<cm.mean_bulk<<",\"coarse_worst_bulk\":"<<cm.worst_bulk<<",\"reference_mean_bulk\":"<<fm.mean_bulk<<",\"reference_worst_bulk\":"<<fm.worst_bulk;
        std::cout<<",\"load_seconds\":"<<load_seconds<<",\"training_ids\":";fixed::vector(input.ids);std::cout<<",\"frozen_scales\":";fixed::vector(cm.frozen_scales);
        std::cout<<",\"coarse_marks\":";fixed::vector(cm.marked_elements);std::cout<<",\"reference_marks\":";fixed::vector(fm.marked_elements);
        std::cout<<",\"worst_member_id\":"<<cm.worst_member_id<<",\"reference_worst_member_id\":"<<fm.worst_member_id<<",\"rank\":"<<regional.phi.cols()<<",\"working_rank\":"<<regional.working_rank
            <<",\"inherited_rank\":"<<regional.inherited_rank<<",\"training_stop\":\""<<regional.stop<<"\",\"training_met\":"<<(regional.training_met?"true":"false")
            <<",\"compression_accepted\":"<<(regional.compression_accepted?"true":"false")<<",\"selected_patches\":"<<regional.selected_patches
            <<",\"training_evaluations\":"<<regional.evaluations<<",\"compression_trials\":"<<regional.compression_trials
            <<",\"kernel_residual\":"<<regional.kernel_residual<<",\"base_orthogonality\":"<<regional.base_orthogonality<<",\"gram_residual\":"<<regional.gram_residual
            <<",\"enrichment_tests\":\""<<(config.tests==EnrichmentTests::KernelLift?"kernel_lift":"adjoint")<<"\",\"kernel_lift_residual\":"<<regional.accepted.kernel_lift_residual<<",\"aot_factorizations\":"<<cache.factorizations()<<",\"aot_solved_columns\":"<<cache.solved_columns()<<",\"aot_identity\":\""<<cache.identity()<<"\",\"targets\":";fixed::vector(regional.targets);
        std::cout<<",\"next_phase\":\""<<(next.reference_sweep==1?"coarse_and_reference":"reference_only")<<"\",\"complete_horizon\":"<<(state+1==horizon?"true":"false");
        if(config.emit){std::cout<<",\"solution\":";fixed::complex_matrix(values);}
        std::cout<<",\"wall_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n"<<std::flush;
        snapshot(CheckpointPhase::Accepted,events.buffer.str(),false,values,regional.phi,cm.marked_elements,fm.marked_elements,space?space->identity():"AFEM",space.get());
        auto journal=events.buffer.str();committed_lines+=std::count(journal.begin(),journal.end(),'\n');prefix_hash=journal_hash(journal,prefix_hash);
        events.publish();
        if(config.pause_state==state&&config.pause_phase=="accepted")return 0;
        if(terminal)break;
        riesz.reset();space.reset();regional=RegionalResult{};
        loads.resize(0,0);values.resize(0,0);E=Sparse{};coarse_mass.resize(0,0);strong.resize(0,0);
        advance(cm.marked_elements,fm.marked_elements);
    }
    return 0;
 }catch(const std::exception& e){std::cerr<<"alod_run adaptive: "<<e.what()<<'\n';return 1;}
}
