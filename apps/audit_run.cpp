#include "audit_run.hpp"
#include "audit_diagnostics.hpp"
#include "mesh_export.hpp"
#include "alod/timing.hpp"
#include "alod/execution.hpp"
#include <functional>
#include "fixed_support.hpp"
#include "alod/checkpoint.hpp"
#include "alod/batch.hpp"
#include "alod/regional.hpp"
#include "alod/localization.hpp"
#include <chrono>
#include <bit>
#include "mesh/refine.h"
using namespace alod;
int audit_main(int argc,char** argv,AuditWorkerState* worker,std::ostream* output,std::ostream* error){
    auto& out=output?*output:std::cout;auto& err=error?*error:std::cerr;try{
    PhaseTimer process_timer("audit_process",-1);
    std::map<std::string,std::string> options{{"localization","0"},{"drop-enrichment","0"},{"retrain-shared","0"},{"reference-levels","0"},{"region-radius","-1"},{"two-level","0"},{"mesh-output",""},{"audit-mode","full"},{"reuse-basis","1"},{"enrichment-tests","adjoint"},{"refinement-steps","0"},{"ell-override","0"},{"quadrature-boost","0"},{"checkpoint",""},{"batch-size","8"},{"fresh","0"},{"rank-zero","1"},{"audit-ids",""},{"method","ALOD"},{"radius",".6"},{"rank-cap","24"},{"maximum-patch-entries","8000000"},{"maximum-dense-entries","8000000"}};
    std::vector<char*> forward{argv[0]};if(argc>1)forward.push_back(argv[1]);std::set<std::string> seen;
    for(int j=2;j<argc;++j){std::string arg=argv[j];auto eq=arg.find('=');auto key=arg.substr(2,eq-2);
        if(arg.starts_with("--")&&eq!=std::string::npos&&options.contains(key)){if(!seen.insert(key).second)throw std::invalid_argument("duplicate audit option");options[key]=arg.substr(eq+1);}else forward.push_back(argv[j]);}
    const bool exact_only=options["audit-mode"]=="exact";if(!exact_only&&options["audit-mode"]!="full")throw std::invalid_argument("invalid audit mode");
    auto policy=options["enrichment-tests"];
    if(policy!="kernel_lift"&&policy!="adjoint")throw std::invalid_argument("invalid enrichment tests");
    auto tests=policy=="kernel_lift"?EnrichmentTests::KernelLift:EnrichmentTests::ArchivedAdjoint;
    auto input=fixed::parse(forward.size(),forward.data());PhaseTimer::probe(input.threads);
    const bool shared_snapshot=worker&&worker->prepared&&worker->prepared_path==options["checkpoint"];
    auto snapshot=[&]{PhaseTimer timing(shared_snapshot?"checkpoint_transfer":"checkpoint_load",-1);
        if(shared_snapshot){auto owned=std::move(worker->prepared);return std::move(*owned);}
        return load_checkpoint(options["checkpoint"]);}();
    PhaseTimer state_timer("audit_state",snapshot.cursor.state_id);
    if(snapshot.phase!=CheckpointPhase::Accepted)throw std::invalid_argument("audit requires an accepted checkpoint");
    // The accepted journal binds the test policy, including direct CLI audits.
    const bool lifted=snapshot.journal.find("\"enrichment_tests\":\"kernel_lift\"")!=std::string::npos;
    if(snapshot.phi.cols() && lifted!=(tests==EnrichmentTests::KernelLift))
        throw std::invalid_argument("audit enrichment policy differs from accepted checkpoint");
    auto integer=[&](std::string key,int lo,int hi){std::size_t n;int v=std::stoi(options[key],&n);if(n!=options[key].size()||v<lo||v>hi)throw std::invalid_argument("invalid audit "+key);return v;};
    int batch=integer("batch-size",1,50),fresh=integer("fresh",0,1),rankzero=integer("rank-zero",0,1),rankcap=integer("rank-cap",1,24);
    int diagnostic_theta=integer("localization",0,1);
    int drop_enrichment=integer("drop-enrichment",0,1);
    int retrain=integer("retrain-shared",0,1),deepening=integer("reference-levels",0,3),two_level=integer("two-level",0,1);
    int refinement_steps=integer("refinement-steps",0,3);
    if(refinement_steps && input.problem!="E1")throw std::invalid_argument("refinement diagnostic is E1 only");
    int ell_override=integer("ell-override",0,4),quadrature_boost=integer("quadrature-boost",0,8);
    if(ell_override && ((snapshot.phi.cols()!=0&&!retrain&&!drop_enrichment) || options["method"]!="ALOD"))
        throw std::invalid_argument("ell override requires rank-zero or shared retraining on an ALOD checkpoint");
    if(exact_only&&(fresh||refinement_steps))throw std::invalid_argument("exact audit does not compute reference/fresh diagnostics");
    std::size_t consumed;double radius=std::stod(options["radius"],&consumed);if(consumed!=options["radius"].size()||!std::isfinite(radius)||radius<0)throw std::invalid_argument("invalid audit radius");
    std::size_t region_used;double region_radius=std::stod(options["region-radius"],&region_used);if(region_used!=options["region-radius"].size()||!std::isfinite(region_radius)||(region_radius<=0&&region_radius!=-1))throw std::invalid_argument("invalid physical region radius");
    bool afem=options["method"]=="AFEM";if(!afem&&options["method"]!="ALOD")throw std::invalid_argument("invalid audit method");
    if(afem!=(snapshot.space_identity=="AFEM"))throw std::invalid_argument("audit method does not match checkpoint");
    if(afem&&(retrain||deepening||two_level||diagnostic_theta))throw std::invalid_argument("LOD diagnostics require ALOD");
    if(retrain&&input.problem!="E2")throw std::invalid_argument("regional retraining requires E2");
    std::set<int> ids;std::istringstream id_stream(options["audit-ids"]);std::string word;
    while(std::getline(id_stream,word,',')){std::size_t n;int i=std::stoi(word,&n);if(n!=word.size()||!ids.insert(i).second)throw std::invalid_argument("invalid audit sample IDs");}
    std::vector<fixed::Member> selected;for(const auto& m:input.members)if(ids.empty()||ids.contains(m.id))selected.push_back(m);
    if(selected.empty()||(!ids.empty()&&selected.size()!=ids.size()))throw std::invalid_argument("unknown audit sample");
    std::map<int,std::string> parameter_ids;{std::ifstream table(input.members_path);std::string row;while(std::getline(table,row)){std::istringstream parsed(row);int id;if(parsed>>id)parameter_ids[id]=std::to_string(journal_hash(row));}}
    auto start=std::chrono::steady_clock::now();auto seconds=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    LodLimits limits;limits.maximum_reference_nodes=input.cap;limits.threads=input.threads;{std::size_t n;auto v=std::stoull(options["maximum-patch-entries"],&n);if(n!=options["maximum-patch-entries"].size()||v<1||v>64000000000ULL)throw std::invalid_argument("invalid patch entry bound");limits.maximum_patch_entries=v;}limits.maximum_dense_entries=integer("maximum-dense-entries",1,1000000000);
    if(static_cast<std::size_t>(snapshot.values.rows())*static_cast<std::size_t>(batch)>limits.maximum_dense_entries)throw std::runtime_error("audit batch exceeds dense allocation bound");
    std::unique_ptr<LodSpace> space;std::unique_ptr<AdditiveKernelRieszContext> riesz;AdjointTestCache aot;std::unique_ptr<RegionalEvaluator> coupled;
    lod2d::helmholtz::HelmholtzOperators afem_operators;
    const bool restored_basis=integer("reuse-basis",0,1)&&!afem&&snapshot.lod_trial.cols()&&(!ell_override||ell_override==snapshot.ell);
    if(shared_snapshot&&(!restored_basis||deepening)&&snapshot.fine.mesh.nodes.empty())
        snapshot=load_checkpoint(options["checkpoint"]); // Explicit rebuild diagnostics need full geometry.
    if(!afem){
        PhaseTimer timing(restored_basis?"lod_restore":"lod_rebuild",snapshot.cursor.state_id);
        const std::string hierarchy_key=snapshot.geometry_file+":"+std::to_string(std::bit_cast<std::uint64_t>(input.wavenumber))+":"+input.policy;
        const bool shared_hierarchy=restored_basis&&worker&&worker->hierarchy&&!snapshot.geometry_file.empty()&&worker->hierarchy_key==hierarchy_key;
        if(shared_hierarchy){
            space=std::make_unique<LodSpace>(worker->hierarchy,snapshot.ell,input.policy=="area"?InterpolationPolicy::ManuscriptAreaWeighted:InterpolationPolicy::ArchivedArithmetic,limits,snapshot.lod_trial,&snapshot.lod_reduced);
            PhaseTimer::counter("audit_hierarchy_hits",1,snapshot.cursor.state_id);
        }
        else if(restored_basis)space=std::make_unique<LodSpace>(snapshot.coarse.mesh,snapshot.reference(),input.wavenumber,snapshot.ell,input.policy=="area"?InterpolationPolicy::ManuscriptAreaWeighted:InterpolationPolicy::ArchivedArithmetic,limits,snapshot.lod_trial,std::getenv("ALOD_REFERENCE_EXECUTION")?nullptr:&snapshot.lod_reduced);
        else space=std::make_unique<LodSpace>(snapshot.coarse.mesh,snapshot.reference(),input.wavenumber,snapshot.ell,input.policy=="area"?InterpolationPolicy::ManuscriptAreaWeighted:InterpolationPolicy::ArchivedArithmetic,limits);
        if(worker){
            worker->hierarchy.reset();worker->hierarchy_key.clear();
            const auto h=space->hierarchy();std::size_t bytes=0;
            auto sparse_bytes=[](const auto& a){return std::size_t(a.nonZeros())*(sizeof(typename std::decay_t<decltype(a)>::Scalar)+sizeof(int))+std::size_t(a.cols()+1)*sizeof(int);};
            for(const Sparse* a:{&h->interpolation,&h->coarse_basis,&h->energy,&h->reference.P_node,&h->reference.P_elem,&h->reference.P_dg,&h->operators.stiffness,&h->operators.mass,&h->operators.boundary_mass})bytes+=sparse_bytes(*a);
            bytes+=sparse_bytes(h->operators.system)+h->operators.element_blocks.size()*sizeof(Eigen::Matrix3cd);
            bytes+=(h->coarse.nodes.size()+h->reference.mesh.nodes.size())*64+(h->coarse.elems.size()+h->reference.mesh.elems.size())*64;
            std::size_t cap=64ULL*1024*1024;
            if(const char* value=std::getenv("ALOD_AUDIT_HIERARCHY_BYTES")){std::string text(value);if(text.empty()||text.find_first_not_of("0123456789")!=std::string::npos)throw std::invalid_argument("invalid audit hierarchy budget");cap=std::stoull(text);}
            if(bytes<=cap&&!snapshot.geometry_file.empty()){worker->hierarchy=h;worker->hierarchy_key=hierarchy_key;}
            PhaseTimer::counter("audit_retained_hierarchy_bytes",worker->hierarchy?bytes:0,snapshot.cursor.state_id);
        }
        if(space->identity()!=snapshot.space_identity)throw std::invalid_argument("checkpoint/audit operator or space identity mismatch");
        if(ell_override && ell_override!=snapshot.ell)space=std::make_unique<LodSpace>(*space,ell_override);
        if(deepening) audit_diagnostics::deepen_reference(
            space, snapshot, deepening, ell_override ? ell_override : snapshot.ell,
            input.wavenumber, input.policy == "area" ? InterpolationPolicy::ManuscriptAreaWeighted
                                                    : InterpolationPolicy::ArchivedArithmetic, limits);
        if(two_level){
            const int current=ell_override?ell_override:snapshot.ell;
            const double difference = audit_diagnostics::two_level_corrector(*space,current);
            out<<"{\"kind\":\"two_level_corrector\",\"ell\":"<<current<<",\"next_ell\":"<<current+1<<",\"basis_energy_frobenius\":"<<difference<<"}\n";
        }
        if(drop_enrichment)snapshot.phi.resize(space->fine().nodes.size(),0);
        if(retrain){
            AdditiveKernelRieszContext training_riesz(*space,input.riesz_policy=="n2"?RieszPatchPolicy::ManuscriptN2:RieszPatchPolicy::ArchivedSupportExpanded);
            std::vector<Problem> training;for(int id:input.ids)for(const auto& m:input.members)if(m.id==id)training.push_back(m.problem);
            auto training_loads=assemble_load_batch(space->fine(),training,paper_quadrature(input.problem),input.threads);
            auto trained=train_regional(*space,training_riesz,aot,training_loads,{radius,rankcap,false,tests});snapshot.phi=std::move(trained.phi);
        }
        if(input.problem=="E2"){if(fresh)riesz=std::make_unique<AdditiveKernelRieszContext>(*space,input.riesz_policy=="n2"?RieszPatchPolicy::ManuscriptN2:RieszPatchPolicy::ArchivedSupportExpanded);coupled=std::make_unique<RegionalEvaluator>(*space,aot,true,tests);}}
    if(diagnostic_theta){
        AdditiveKernelRieszContext theta_riesz(*space,input.riesz_policy=="n2"?RieszPatchPolicy::ManuscriptN2:RieszPatchPolicy::ArchivedSupportExpanded);
        LocalizationEigenConfig eig;eig.relative_tolerance=input.tolerance;eig.maximum_iterations=input.iterations;eig.dense_cross_check_max_dimension=input.dense;eig.eigenvalue_relative_residual=true;
        auto theta=localization_theta(*space,theta_riesz,eig);
        out<<std::setprecision(17)<<"{\"kind\":\"diagnostic_theta\",\"theta\":"<<theta.theta<<",\"ritz_residual\":"<<theta.spectrum.relative_residual<<",\"ritz_iterations\":"<<theta.spectrum.iterations<<",\"relative_tolerance\":"<<input.tolerance<<"}\n";
    }
    // Only accepted solutions/dictionary and fine geometry are used below.
    snapshot.lod_trial=ComplexSparseMatrix{};snapshot.lod_reduced=ComplexSparseMatrix{};snapshot.P_node=Sparse{};snapshot.P_elem=Sparse{};snapshot.P_dg=Sparse{};
    snapshot.raw_kernel.resize(0,0);snapshot.warm_full.resize(0,0);
    snapshot.coarse=MeshState(lod2d::TriMesh{});snapshot.fine.elements.clear();snapshot.fine.elements.shrink_to_fit();snapshot.fine.ancestry.clear();snapshot.fine.ancestry.shrink_to_fit();
    if(afem)afem_operators=lod2d::helmholtz::assemble_helmholtz_operators(snapshot.fine.mesh,input.wavenumber);
    const auto& audit_mesh=afem?snapshot.fine.mesh:space->fine();
    if(!afem)snapshot.fine=MeshState(lod2d::TriMesh{});
    // Select once for the entire member family, rather than copying and filtering
    // the same fine mesh again for every RHS. Disabled diagnostics allocate nothing.
    lod2d::TriMesh region;
    if(region_radius>0) region=audit_diagnostics::centroid_region(audit_mesh,region_radius);
    export_mesh_pair(options["mesh-output"],afem?audit_mesh:space->coarse(),audit_mesh);
    const auto& operators=afem?afem_operators:space->operators();
    Sparse afem_energy;if(afem)afem_energy=operators.stiffness+input.wavenumber*input.wavenumber*operators.mass;
    const Sparse& E=afem?afem_energy:space->energy();
    const auto reference_mesh_hash=mesh_fingerprint(audit_mesh),dictionary_hash=matrix_hash(snapshot.phi);
    const bool accepted_only=exact_only&&!retrain&&!drop_enrichment&&!deepening&&!ell_override&&!quadrature_boost
        &&std::all_of(selected.begin(),selected.end(),[&](const auto& m){return std::find(snapshot.computed_ids.begin(),snapshot.computed_ids.end(),m.id)!=snapshot.computed_ids.end();});
    double accepted_pg=0.;
    if(accepted_only){auto position=snapshot.journal.find("\"pg_residual\":");if(position==std::string::npos)throw std::runtime_error("accepted PG evidence missing");accepted_pg=std::stod(snapshot.journal.substr(position+14));if(!std::isfinite(accepted_pg)||accepted_pg>1e-8||accepted_pg<0)throw std::runtime_error("invalid accepted PG evidence");}
    std::unique_ptr<ReferenceFemContext> reference;
    if(!exact_only||(afem&&!accepted_only)){PhaseTimer timing("reference_factor",snapshot.cursor.state_id);reference=std::make_unique<ReferenceFemContext>(operators);}double prepare_seconds=seconds(),load_seconds=0,error_seconds=0,coupled_seconds=0,fresh_seconds=0;
    auto q=paper_quadrature(input.problem);q.base_triangle_order+=quadrature_boost;q.gaussian_triangle_order+=quadrature_boost;q.singular_triangle_order+=quadrature_boost;out<<std::setprecision(17);int completed=0;
    std::unique_ptr<AuditIntegrationGeometry> integration_geometry;
    if(!std::getenv("ALOD_REFERENCE_EXECUTION")){PhaseTimer timing("audit_geometry",snapshot.cursor.state_id);
        integration_geometry=std::make_unique<AuditIntegrationGeometry>(audit_mesh,input.threads);}
    // A bounded three-stage pipeline: one prefetched load, one solve, and one
    // pending integration. Factors stay on the foreground thread; only immutable
    // batch values reach the integrator. Fresh-training diagnostics remain serial.
    const bool pipeline=asynchronous_execution(input.threads)&&input.threads>=3&&!fresh&&!accepted_only
        &&selected.size()>static_cast<std::size_t>(batch)
        &&asynchronous_buffers_fit(std::size_t(audit_mesh.nodes.size())*batch*sizeof(Complex)*5);
    // Quadrature stages share one persistent team rather than competing for
    // memory bandwidth. The foreground exclusively owns the sparse factors.
    const int integral_workers=pipeline?input.threads-1:input.threads;
    ExecutionScope solve_scope(thread_budget(pipeline?1:input.threads));
    auto integrate_budget=thread_budget(integral_workers);
    auto load_batch=[&](int begin){
        auto start=seconds();PhaseTimer timer("audit_load",snapshot.cursor.state_id);
        std::vector<Problem> problems;
        for(int j=begin;j<std::min<int>(begin+batch,selected.size());++j)problems.push_back(selected[j].problem);
        auto loads=accepted_only?ComplexMatrix{}:assemble_load_batch(audit_mesh,problems,q,input.threads);
        return std::make_pair(std::move(loads),seconds()-start);
    };
    std::unique_ptr<TaskLane> quadrature_lane;
    if(pipeline)quadrature_lane=std::make_unique<TaskLane>(integrate_budget);
    std::future<std::pair<ComplexMatrix,double>> next_load;
    if(pipeline)next_load=quadrature_lane->submit([&]{return load_batch(0);});
    std::function<void(ErrorBatch)> emit_metrics;
    // The lane drains on exceptional unwinding before borrowed mesh/solver state
    // is destroyed. Batch payloads remain owned by their queued tasks.
    std::future<std::pair<ErrorBatch,double>> pending_errors;
    auto drain=[&]{if(pending_errors.valid()){
        PhaseTimer wait("audit_integral_wait",snapshot.cursor.state_id);
        auto result=pending_errors.get();error_seconds+=result.second;
        emit_metrics(std::move(result.first));emit_metrics={};
    }};
    for(int begin=0;begin<static_cast<int>(selected.size());begin+=batch){int count=std::min(batch,static_cast<int>(selected.size())-begin);std::vector<Problem> problems;for(int j=0;j<count;++j)problems.push_back(selected[begin+j].problem);
        auto loaded=next_load.valid()?next_load.get():load_batch(begin);
        auto loads=std::move(loaded.first);load_seconds+=loaded.second;double t=seconds();
        if(pipeline&&begin+batch<static_cast<int>(selected.size()))
            next_load=quadrature_lane->submit([&,next=begin+batch]{return load_batch(next);});
        auto ref=[&]{PhaseTimer timing("reference_solve",snapshot.cursor.state_id);return reference?reference->solve(loads):ComplexMatrix::Zero(loads.rows(),loads.cols()).eval();}();ComplexMatrix values,base;double pg=0;
        t=seconds();{PhaseTimer timing("frozen_solve",snapshot.cursor.state_id);if(accepted_only){values.resize(audit_mesh.nodes.size(),count);for(int j=0;j<count;++j){auto found=std::find(snapshot.computed_ids.begin(),snapshot.computed_ids.end(),selected[begin+j].id);values.col(j)=snapshot.values.col(found-snapshot.computed_ids.begin());}pg=accepted_pg;}else if(afem)values=ref;else if(coupled){auto result=coupled->evaluate_frozen(loads,snapshot.phi);values=std::move(result.values);pg=result.pg_residual;if(rankzero&&!exact_only)base=space->solve(loads).values;}
        else {auto result=space->solve(loads);values=std::move(result.values);pg=result.pg_relative_residual;}}coupled_seconds+=seconds()-t;
        for(int step=0;step<refinement_steps;++step){
            ComplexMatrix residual=loads-operators.system*ref;
            for(int node:operators.dirichlet_nodes)residual.row(node).setZero();
            ref+=reference->solve(residual);
            if(afem)values=ref;
            else {residual=loads-operators.system*values;
                for(int node:operators.dirichlet_nodes)residual.row(node).setZero();
                values+=space->solve(residual).values;}
        }
        double audited_reference_residual=0.;
        if(reference){ComplexMatrix ref_defect=operators.system*ref-loads,ref_rhs=loads;
        for(int node:operators.dirichlet_nodes){ref_defect.row(node).setZero();ref_rhs.row(node).setZero();}
        audited_reference_residual=ref_defect.norm()/std::max(1e-30,ref_rhs.norm());
        ref_defect.resize(0,0);ref_rhs.resize(0,0);}
        if(!afem && refinement_steps)pg=(space->test().adjoint()*(operators.system*values-loads)).norm()/std::max(1e-30,(space->test().adjoint()*loads).norm());
        if(refinement_steps && std::max(pg,audited_reference_residual)>1e-10)throw std::runtime_error("strict refined residual gate failed");
        std::string exact_family=snapshot.members_text;for(int j=0;j<count;++j)exact_family+=":"+std::to_string(selected[begin+j].id);
        // Finish the previous batch before reusing its integration cache or
        // emitting this batch. Ordered output is independent of task completion.
        drain();
        struct BatchValues {ComplexMatrix values,ref,base,loads;std::vector<Problem> problems;};
        auto owned=std::make_shared<BatchValues>(BatchValues{std::move(values),std::move(ref),std::move(base),std::move(loads),std::move(problems)});
        auto integrate=[&,owned,exact_family]{
            auto start=seconds();PhaseTimer timing("error_integral",snapshot.cursor.state_id);
            auto errors=exact_only?integrate_error_batch(audit_mesh,E,owned->values,owned->problems,q,input.threads,worker?&worker->exact_reuse:nullptr,exact_family)
                :integrate_audit_batch(audit_mesh,E,owned->values,owned->ref,owned->problems,q,input.threads,true,integration_geometry.get(),worker?&worker->exact_reuse:nullptr,exact_family);
            return std::make_pair(std::move(errors),seconds()-start);
        };
        const double reference_pg=reference?reference->relative_residual():0.;
        emit_metrics=[&,owned,begin,count,pg,audited_reference_residual,reference_pg](ErrorBatch errors){
        const auto& values=owned->values;const auto& ref=owned->ref;const auto& base=owned->base;const auto& loads=owned->loads;
        const auto& floors=errors.reference_error;double t=seconds();
        PhaseTimer metrics_timer("sample_metrics",snapshot.cursor.state_id);
        for(int j=0;j<count;++j){const auto& member=selected[begin+j];double n=errors.exact_norm[j];auto norm=[&](const ComplexVector& x){return std::sqrt(std::max(0.,x.dot(E.cast<Complex>()*x).real()));};
            double gap=(afem||exact_only)?0:norm((ref.col(j)-values.col(j)).eval());
            out<<"{\"kind\":\"sample\",\"state_id\":"<<snapshot.cursor.state_id<<",\"sample\":"<<member.id<<",\"role\":"<<json_string(member.role)<<",\"nominal\":"<<(member.id==0?"true":"false")
                <<",\"marking_member\":"<<(std::find(input.ids.begin(),input.ids.end(),member.id)!=input.ids.end()?"true":"false")<<",\"ell\":";if(afem)out<<"null";else out<<(ell_override?ell_override:snapshot.ell);out<<",\"rank\":"<<snapshot.phi.cols()<<",\"energy\":"<<errors.energy[j]
                <<",\"exact_norm\":"<<n<<",\"exact_error\":"<<errors.exact_error[j]<<",\"e\":"<<errors.exact_error[j]<<",\"f\":";if(afem||exact_only)out<<"null";else out<<floors[j];
            out<<",\"g\":";if(afem||exact_only)out<<"null";else out<<gap;
            auto ratio=[&](double x){if(n>1e-12)out<<x/n;else out<<"null";};
            out<<",\"E\":";ratio(errors.exact_error[j]);out<<",\"E_ref\":";if(afem||exact_only)out<<"null";else ratio(floors[j]);out<<",\"F\":";if(afem||exact_only)out<<"null";else ratio(floors[j]);out<<",\"G\":";if(afem||exact_only)out<<"null";else ratio(gap);
            auto saved=std::find(snapshot.computed_ids.begin(),snapshot.computed_ids.end(),member.id);
            out<<",\"accepted_solution_distance\":";
            if(saved==snapshot.computed_ids.end())out<<"null";
            else out<<norm((values.col(j)-snapshot.values.col(saved-snapshot.computed_ids.begin())).eval());
            out<<",\"wavenumber\":"<<input.wavenumber<<",\"source_ell\":"<<snapshot.ell<<",\"quadrature_boost\":"<<quadrature_boost;
            out<<",\"reference_residual\":";if(!reference)out<<"null";else out<<audited_reference_residual;out<<",\"refinement_steps\":"<<refinement_steps;
            out<<",\"parameter_id\":"<<json_string(parameter_ids.at(member.id))<<",\"reference_mesh_hash\":"<<json_string(reference_mesh_hash)<<",\"dictionary_hash\":"<<json_string(dictionary_hash)<<",\"solution_hash\":"<<json_string(matrix_hash(values.col(j)));
            out<<",\"ratio_status\":\""<<(n>1e-12?"finite":"near_zero_exact_norm")<<"\",\"reference_status\":\""<<(afem?"not_applicable":exact_only?"not_computed":"available")<<"\",\"PG_residual\":"<<(afem&&reference?reference_pg:pg);
            if(rankzero&&!afem&&!exact_only)out<<",\"base_gap\":"<<(coupled?norm((ref.col(j)-base.col(j)).eval()):gap)<<",\"as_correction_energy\":"<<(coupled?norm((values.col(j)-base.col(j)).eval()):0.);
            if(fresh&&!afem){t=seconds();ComplexVector u;int rank=0;if(input.problem=="E2"){auto trained=train_regional(*space,*riesz,aot,loads.col(j),{radius,rankcap,false,tests});u=trained.accepted.values.col(0);rank=trained.phi.cols();}else u=space->solve(loads.col(j)).values.col(0);
                double error=lod2d::helmholtz::compute_helmholtz_error(audit_mesh,u,input.wavenumber,member.problem.exact,member.problem.exact_gradient,q,member.problem.quadrature_context).energy;fresh_seconds+=seconds()-t;
                out<<",\"fresh_rank\":"<<rank<<",\"fresh_e\":"<<error<<",\"fresh_g\":"<<norm((ref.col(j)-u).eval())<<",\"shared_fresh_signed_percent\":";if(error>1e-12)out<<100*(errors.exact_error[j]-error)/error;else out<<"null";}
            if(region_radius>0){
                double region_error=0.,region_norm=0.;
                if(!region.elems.empty()){
                    region_error=lod2d::helmholtz::compute_helmholtz_error(region,values.col(j),input.wavenumber,member.problem.exact,member.problem.exact_gradient,q,member.problem.quadrature_context).energy;
                    region_norm=lod2d::helmholtz::compute_helmholtz_error(region,ComplexVector::Zero(values.rows()),input.wavenumber,member.problem.exact,member.problem.exact_gradient,q,member.problem.quadrature_context).energy;
                }
                out<<",\"region_radius\":"<<region_radius<<",\"region_elements\":"<<region.elems.size()<<",\"region_e\":"<<region_error<<",\"region_exact_norm\":"<<region_norm<<",\"region_E\":";if(region_norm>1e-12)out<<region_error/region_norm;else out<<"null";
            }
            out<<"}\n";++completed;
        }
        };
        if(pipeline)pending_errors=quadrature_lane->submit(std::move(integrate));
        else {auto result=integrate();error_seconds+=result.second;emit_metrics(std::move(result.first));emit_metrics={};}
    }
    integrate_budget->store(input.threads);drain();
    out<<"{\"kind\":\"audit_complete\",\"state_id\":"<<snapshot.cursor.state_id<<",\"samples\":"<<completed<<",\"batch_size\":"<<batch<<",\"reference_factorizations\":"<<(reference?1:0)<<",\"audit_mode\":"<<json_string(options["audit-mode"])<<",\"aot_factorizations\":"<<aot.factorizations()
        <<",\"batch_pipeline\":"<<(pipeline?"true":"false")<<",\"accepted_values_reused\":"<<(accepted_only?"true":"false")<<",\"snapshot_transferred\":"<<(shared_snapshot?"true":"false")<<",\"lod_basis_reused\":"<<(restored_basis?"true":"false")<<",\"lod_patch_rebuilds\":"<<(!afem&&!restored_basis?1:0)<<",\"threads\":"<<input.threads<<",\"prepare_seconds\":"<<prepare_seconds<<",\"load_seconds\":"<<load_seconds<<",\"error_seconds\":"<<error_seconds<<",\"reference_factor_seconds\":"<<(reference?reference->factor_seconds():0.)<<",\"reference_solve_seconds\":"<<(reference?reference->solve_seconds():0.)
        <<",\"coupled_seconds\":"<<coupled_seconds<<",\"fresh_seconds\":"<<fresh_seconds<<",\"wall_seconds\":"<<seconds()<<"}\n";
    return 0;
}catch(const std::exception& e){err<<"alod_audit: "<<e.what()<<'\n';return 1;}}
