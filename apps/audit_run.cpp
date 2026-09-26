#include "alod/timing.hpp"
#include "fixed_support.hpp"
#include "alod/checkpoint.hpp"
#include "alod/batch.hpp"
#include "alod/regional.hpp"
#include <chrono>
using namespace alod;
int audit_main(int argc,char** argv){try{
    PhaseTimer process_timer("audit_process",-1);
    std::map<std::string,std::string> options{{"reuse-basis","1"},{"enrichment-tests","adjoint"},{"refinement-steps","0"},{"ell-override","0"},{"quadrature-boost","0"},{"checkpoint",""},{"batch-size","8"},{"fresh","0"},{"rank-zero","1"},{"audit-ids",""},{"method","ALOD"},{"radius",".6"},{"rank-cap","24"},{"maximum-patch-entries","8000000"},{"maximum-dense-entries","8000000"}};
    std::vector<char*> forward{argv[0]};if(argc>1)forward.push_back(argv[1]);std::set<std::string> seen;
    for(int j=2;j<argc;++j){std::string arg=argv[j];auto eq=arg.find('=');auto key=arg.substr(2,eq-2);
        if(arg.starts_with("--")&&eq!=std::string::npos&&options.contains(key)){if(!seen.insert(key).second)throw std::invalid_argument("duplicate audit option");options[key]=arg.substr(eq+1);}else forward.push_back(argv[j]);}
    auto policy=options["enrichment-tests"];
    if(policy!="kernel_lift"&&policy!="adjoint")throw std::invalid_argument("invalid enrichment tests");
    auto tests=policy=="kernel_lift"?EnrichmentTests::KernelLift:EnrichmentTests::ArchivedAdjoint;
    auto input=fixed::parse(forward.size(),forward.data());PhaseTimer::probe(input.threads);
    auto snapshot=[&]{PhaseTimer timing("checkpoint_load",-1);return load_checkpoint(options["checkpoint"]);}();
    PhaseTimer state_timer("audit_state",snapshot.cursor.state_id);
    if(snapshot.phase!=CheckpointPhase::Accepted)throw std::invalid_argument("audit requires an accepted checkpoint");
    // The accepted journal binds the test policy, including direct CLI audits.
    const bool lifted=snapshot.journal.find("\"enrichment_tests\":\"kernel_lift\"")!=std::string::npos;
    if(snapshot.phi.cols() && lifted!=(tests==EnrichmentTests::KernelLift))
        throw std::invalid_argument("audit enrichment policy differs from accepted checkpoint");
    auto integer=[&](std::string key,int lo,int hi){std::size_t n;int v=std::stoi(options[key],&n);if(n!=options[key].size()||v<lo||v>hi)throw std::invalid_argument("invalid audit "+key);return v;};
    int batch=integer("batch-size",1,50),fresh=integer("fresh",0,1),rankzero=integer("rank-zero",0,1),rankcap=integer("rank-cap",1,24);
    int refinement_steps=integer("refinement-steps",0,3);
    if(refinement_steps && input.problem!="E1")throw std::invalid_argument("refinement diagnostic is E1 only");
    int ell_override=integer("ell-override",0,4),quadrature_boost=integer("quadrature-boost",0,8);
    if(ell_override && (input.problem!="E1" || snapshot.phi.cols()!=0 || options["method"]!="ALOD"))
        throw std::invalid_argument("ell override requires an E1 rank-zero ALOD checkpoint");
    std::size_t consumed;double radius=std::stod(options["radius"],&consumed);if(consumed!=options["radius"].size()||!std::isfinite(radius)||radius<0)throw std::invalid_argument("invalid audit radius");
    bool afem=options["method"]=="AFEM";if(!afem&&options["method"]!="ALOD")throw std::invalid_argument("invalid audit method");
    if(afem!=(snapshot.space_identity=="AFEM"))throw std::invalid_argument("audit method does not match checkpoint");
    std::set<int> ids;std::istringstream id_stream(options["audit-ids"]);std::string word;
    while(std::getline(id_stream,word,',')){std::size_t n;int i=std::stoi(word,&n);if(n!=word.size()||!ids.insert(i).second)throw std::invalid_argument("invalid audit sample IDs");}
    std::vector<fixed::Member> selected;for(const auto& m:input.members)if(ids.empty()||ids.contains(m.id))selected.push_back(m);
    if(selected.empty()||(!ids.empty()&&selected.size()!=ids.size()))throw std::invalid_argument("unknown audit sample");
    std::map<int,std::string> parameter_ids;{std::ifstream table(input.members_path);std::string row;while(std::getline(table,row)){std::istringstream parsed(row);int id;if(parsed>>id)parameter_ids[id]=std::to_string(journal_hash(row));}}
    auto start=std::chrono::steady_clock::now();auto seconds=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    LodLimits limits;limits.maximum_reference_nodes=input.cap;limits.threads=input.threads;{std::size_t n;auto v=std::stoull(options["maximum-patch-entries"],&n);if(n!=options["maximum-patch-entries"].size()||v<1||v>64000000000ULL)throw std::invalid_argument("invalid patch entry bound");limits.maximum_patch_entries=v;}limits.maximum_dense_entries=integer("maximum-dense-entries",1,1000000000);
    if(snapshot.fine.mesh.nodes.size()*static_cast<std::size_t>(batch)>limits.maximum_dense_entries)throw std::runtime_error("audit batch exceeds dense allocation bound");
    std::unique_ptr<LodSpace> space;std::unique_ptr<AdditiveKernelRieszContext> riesz;AdjointTestCache aot;std::unique_ptr<RegionalEvaluator> coupled;
    lod2d::helmholtz::HelmholtzOperators afem_operators;
    const bool restored_basis=integer("reuse-basis",0,1)&&!afem&&snapshot.lod_trial.cols()&&(!ell_override||ell_override==snapshot.ell);
    if(!afem){
        PhaseTimer timing(restored_basis?"lod_restore":"lod_rebuild",snapshot.cursor.state_id);
        if(restored_basis)space=std::make_unique<LodSpace>(snapshot.coarse.mesh,snapshot.reference(),input.wavenumber,snapshot.ell,input.policy=="area"?InterpolationPolicy::ManuscriptAreaWeighted:InterpolationPolicy::ArchivedArithmetic,limits,snapshot.lod_trial,std::getenv("ALOD_REFERENCE_EXECUTION")?nullptr:&snapshot.lod_reduced);
        else space=std::make_unique<LodSpace>(snapshot.coarse.mesh,snapshot.reference(),input.wavenumber,snapshot.ell,input.policy=="area"?InterpolationPolicy::ManuscriptAreaWeighted:InterpolationPolicy::ArchivedArithmetic,limits);
        if(space->identity()!=snapshot.space_identity)throw std::invalid_argument("checkpoint/audit operator or space identity mismatch");
        if(ell_override && ell_override!=snapshot.ell)space=std::make_unique<LodSpace>(*space,ell_override);
        if(input.problem=="E2"){if(fresh)riesz=std::make_unique<AdditiveKernelRieszContext>(*space,input.riesz_policy=="n2"?RieszPatchPolicy::ManuscriptN2:RieszPatchPolicy::ArchivedSupportExpanded);coupled=std::make_unique<RegionalEvaluator>(*space,aot,true,tests);}}
    // Only accepted solutions/dictionary and fine geometry are used below.
    snapshot.lod_trial=ComplexSparseMatrix{};snapshot.lod_reduced=ComplexSparseMatrix{};snapshot.P_node=Sparse{};snapshot.P_elem=Sparse{};snapshot.P_dg=Sparse{};
    snapshot.raw_kernel.resize(0,0);snapshot.warm_full.resize(0,0);
    snapshot.coarse=MeshState(lod2d::TriMesh{});snapshot.fine.elements.clear();snapshot.fine.elements.shrink_to_fit();snapshot.fine.ancestry.clear();snapshot.fine.ancestry.shrink_to_fit();
    if(afem)afem_operators=lod2d::helmholtz::assemble_helmholtz_operators(snapshot.fine.mesh,input.wavenumber);
    const auto& audit_mesh=afem?snapshot.fine.mesh:space->fine();
    if(!afem)snapshot.fine=MeshState(lod2d::TriMesh{});
    const auto& operators=afem?afem_operators:space->operators();
    Sparse afem_energy;if(afem)afem_energy=operators.stiffness+input.wavenumber*input.wavenumber*operators.mass;
    const Sparse& E=afem?afem_energy:space->energy();
    const auto reference_mesh_hash=mesh_fingerprint(audit_mesh),dictionary_hash=matrix_hash(snapshot.phi);
    auto reference_ptr=[&]{PhaseTimer timing("reference_factor",snapshot.cursor.state_id);return std::make_unique<ReferenceFemContext>(operators);}();auto& reference=*reference_ptr;double prepare_seconds=seconds(),load_seconds=0,error_seconds=0,coupled_seconds=0,fresh_seconds=0;
    auto q=paper_quadrature(input.problem);q.base_triangle_order+=quadrature_boost;q.gaussian_triangle_order+=quadrature_boost;q.singular_triangle_order+=quadrature_boost;std::cout<<std::setprecision(17);int completed=0;
    std::unique_ptr<AuditIntegrationGeometry> integration_geometry;
    if(!std::getenv("ALOD_REFERENCE_EXECUTION")){PhaseTimer timing("audit_geometry",snapshot.cursor.state_id);
        integration_geometry=std::make_unique<AuditIntegrationGeometry>(audit_mesh,input.threads);}
    for(int begin=0;begin<static_cast<int>(selected.size());begin+=batch){int count=std::min(batch,static_cast<int>(selected.size())-begin);std::vector<Problem> problems;for(int j=0;j<count;++j)problems.push_back(selected[begin+j].problem);
        double t=seconds();auto loads=[&]{PhaseTimer timing("audit_load",snapshot.cursor.state_id);return assemble_load_batch(audit_mesh,problems,q,input.threads);}();load_seconds+=seconds()-t;
        auto ref=[&]{PhaseTimer timing("reference_solve",snapshot.cursor.state_id);return reference.solve(loads);}();ComplexMatrix values,base;double pg=0;
        t=seconds();{PhaseTimer timing("frozen_solve",snapshot.cursor.state_id);if(afem)values=ref;else if(coupled){auto result=coupled->evaluate_frozen(loads,snapshot.phi);values=std::move(result.values);pg=result.pg_residual;if(rankzero)base=space->solve(loads).values;}
        else {auto result=space->solve(loads);values=std::move(result.values);pg=result.pg_relative_residual;}}coupled_seconds+=seconds()-t;
        for(int step=0;step<refinement_steps;++step){
            ComplexMatrix residual=loads-operators.system*ref;
            for(int node:operators.dirichlet_nodes)residual.row(node).setZero();
            ref+=reference.solve(residual);
            if(afem)values=ref;
            else {residual=loads-operators.system*values;
                for(int node:operators.dirichlet_nodes)residual.row(node).setZero();
                values+=space->solve(residual).values;}
        }
        ComplexMatrix ref_defect=operators.system*ref-loads,ref_rhs=loads;
        for(int node:operators.dirichlet_nodes){ref_defect.row(node).setZero();ref_rhs.row(node).setZero();}
        double audited_reference_residual=ref_defect.norm()/std::max(1e-30,ref_rhs.norm());
        ref_defect.resize(0,0);ref_rhs.resize(0,0);
        if(!afem && refinement_steps)pg=(space->test().adjoint()*(operators.system*values-loads)).norm()/std::max(1e-30,(space->test().adjoint()*loads).norm());
        if(refinement_steps && std::max(pg,audited_reference_residual)>1e-10)throw std::runtime_error("strict refined residual gate failed");
        t=seconds();auto errors=[&]{PhaseTimer timing("error_integral",snapshot.cursor.state_id);return integrate_audit_batch(audit_mesh,E,values,ref,problems,q,input.threads,true,integration_geometry.get());}();
        const auto& floors=errors.reference_error;error_seconds+=seconds()-t;
        PhaseTimer metrics_timer("sample_metrics",snapshot.cursor.state_id);
        for(int j=0;j<count;++j){const auto& member=selected[begin+j];double n=errors.exact_norm[j];auto norm=[&](const ComplexVector& x){return std::sqrt(std::max(0.,x.dot(E.cast<Complex>()*x).real()));};
            double gap=afem?0:norm((ref.col(j)-values.col(j)).eval());
            std::cout<<"{\"kind\":\"sample\",\"state_id\":"<<snapshot.cursor.state_id<<",\"sample\":"<<member.id<<",\"role\":"<<json_string(member.role)<<",\"nominal\":"<<(member.id==0?"true":"false")
                <<",\"marking_member\":"<<(std::find(input.ids.begin(),input.ids.end(),member.id)!=input.ids.end()?"true":"false")<<",\"ell\":";if(afem)std::cout<<"null";else std::cout<<(ell_override?ell_override:snapshot.ell);std::cout<<",\"rank\":"<<snapshot.phi.cols()<<",\"energy\":"<<errors.energy[j]
                <<",\"exact_norm\":"<<n<<",\"exact_error\":"<<errors.exact_error[j]<<",\"e\":"<<errors.exact_error[j]<<",\"f\":";if(afem)std::cout<<"null";else std::cout<<floors[j];
            std::cout<<",\"g\":";if(afem)std::cout<<"null";else std::cout<<gap;
            auto ratio=[&](double x){if(n>1e-12)std::cout<<x/n;else std::cout<<"null";};
            std::cout<<",\"E\":";ratio(errors.exact_error[j]);std::cout<<",\"E_ref\":";if(afem)std::cout<<"null";else ratio(floors[j]);std::cout<<",\"F\":";if(afem)std::cout<<"null";else ratio(floors[j]);std::cout<<",\"G\":";if(afem)std::cout<<"null";else ratio(gap);
            auto saved=std::find(snapshot.computed_ids.begin(),snapshot.computed_ids.end(),member.id);
            std::cout<<",\"accepted_solution_distance\":";
            if(saved==snapshot.computed_ids.end())std::cout<<"null";
            else std::cout<<norm((values.col(j)-snapshot.values.col(saved-snapshot.computed_ids.begin())).eval());
            std::cout<<",\"wavenumber\":"<<input.wavenumber<<",\"source_ell\":"<<snapshot.ell<<",\"quadrature_boost\":"<<quadrature_boost;
            std::cout<<",\"reference_residual\":"<<audited_reference_residual<<",\"refinement_steps\":"<<refinement_steps;
            std::cout<<",\"parameter_id\":"<<json_string(parameter_ids.at(member.id))<<",\"reference_mesh_hash\":"<<json_string(reference_mesh_hash)<<",\"dictionary_hash\":"<<json_string(dictionary_hash)<<",\"solution_hash\":"<<json_string(matrix_hash(values.col(j)));
            std::cout<<",\"ratio_status\":\""<<(n>1e-12?"finite":"near_zero_exact_norm")<<"\",\"reference_status\":\""<<(afem?"not_applicable":"available")<<"\",\"PG_residual\":"<<(afem?reference.relative_residual():pg);
            if(rankzero&&!afem)std::cout<<",\"base_gap\":"<<(coupled?norm((ref.col(j)-base.col(j)).eval()):gap)<<",\"as_correction_energy\":"<<(coupled?norm((values.col(j)-base.col(j)).eval()):0.);
            if(fresh&&!afem){t=seconds();ComplexVector u;int rank=0;if(input.problem=="E2"){auto trained=train_regional(*space,*riesz,aot,loads.col(j),{radius,rankcap,false,tests});u=trained.accepted.values.col(0);rank=trained.phi.cols();}else u=space->solve(loads.col(j)).values.col(0);
                double error=lod2d::helmholtz::compute_helmholtz_error(audit_mesh,u,input.wavenumber,member.problem.exact,member.problem.exact_gradient,q,member.problem.quadrature_context).energy;fresh_seconds+=seconds()-t;
                std::cout<<",\"fresh_rank\":"<<rank<<",\"fresh_e\":"<<error<<",\"fresh_g\":"<<norm((ref.col(j)-u).eval())<<",\"shared_fresh_signed_percent\":";if(error>1e-12)std::cout<<100*(errors.exact_error[j]-error)/error;else std::cout<<"null";}
            std::cout<<"}\n";++completed;
        }
    }
    std::cout<<"{\"kind\":\"audit_complete\",\"state_id\":"<<snapshot.cursor.state_id<<",\"samples\":"<<completed<<",\"batch_size\":"<<batch<<",\"reference_factorizations\":1,\"aot_factorizations\":"<<aot.factorizations()
        <<",\"lod_basis_reused\":"<<(restored_basis?"true":"false")<<",\"lod_patch_rebuilds\":"<<(!afem&&!restored_basis?1:0)<<",\"threads\":"<<input.threads<<",\"prepare_seconds\":"<<prepare_seconds<<",\"load_seconds\":"<<load_seconds<<",\"error_seconds\":"<<error_seconds<<",\"reference_factor_seconds\":"<<reference.factor_seconds()<<",\"reference_solve_seconds\":"<<reference.solve_seconds()
        <<",\"coupled_seconds\":"<<coupled_seconds<<",\"fresh_seconds\":"<<fresh_seconds<<",\"wall_seconds\":"<<seconds()<<"}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<"alod_audit: "<<e.what()<<'\n';return 1;}}
