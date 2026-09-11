#include "fixed_support.hpp"
#include "alod/checkpoint.hpp"
#include "alod/batch.hpp"
#include "alod/regional.hpp"
#include <chrono>
using namespace alod;
int main(int argc,char** argv){try{
    std::map<std::string,std::string> options{{"checkpoint",""},{"batch-size","8"},{"fresh","0"},{"rank-zero","1"},{"audit-ids",""},{"method","ALOD"},{"radius",".6"},{"rank-cap","24"},{"maximum-patch-entries","8000000"},{"maximum-dense-entries","8000000"}};
    std::vector<char*> forward{argv[0]};if(argc>1)forward.push_back(argv[1]);std::set<std::string> seen;
    for(int j=2;j<argc;++j){std::string arg=argv[j];auto eq=arg.find('=');auto key=arg.substr(2,eq-2);
        if(arg.starts_with("--")&&eq!=std::string::npos&&options.contains(key)){if(!seen.insert(key).second)throw std::invalid_argument("duplicate audit option");options[key]=arg.substr(eq+1);}else forward.push_back(argv[j]);}
    auto input=fixed::parse(forward.size(),forward.data());auto snapshot=load_checkpoint(options["checkpoint"]);
    if(snapshot.phase!=CheckpointPhase::Accepted)throw std::invalid_argument("audit requires an accepted checkpoint");
    auto integer=[&](std::string key,int lo,int hi){std::size_t n;int v=std::stoi(options[key],&n);if(n!=options[key].size()||v<lo||v>hi)throw std::invalid_argument("invalid audit "+key);return v;};
    int batch=integer("batch-size",1,50),fresh=integer("fresh",0,1),rankzero=integer("rank-zero",0,1),rankcap=integer("rank-cap",1,24);
    std::size_t consumed;double radius=std::stod(options["radius"],&consumed);if(consumed!=options["radius"].size()||!std::isfinite(radius)||radius<0)throw std::invalid_argument("invalid audit radius");
    bool afem=options["method"]=="AFEM";if(!afem&&options["method"]!="ALOD")throw std::invalid_argument("invalid audit method");
    if(afem!=(snapshot.space_identity=="AFEM"))throw std::invalid_argument("audit method does not match checkpoint");
    std::set<int> ids;std::istringstream id_stream(options["audit-ids"]);std::string word;
    while(std::getline(id_stream,word,',')){std::size_t n;int i=std::stoi(word,&n);if(n!=word.size()||!ids.insert(i).second)throw std::invalid_argument("invalid audit sample IDs");}
    std::vector<fixed::Member> selected;for(const auto& m:input.members)if(ids.empty()||ids.contains(m.id))selected.push_back(m);
    if(selected.empty()||(!ids.empty()&&selected.size()!=ids.size()))throw std::invalid_argument("unknown audit sample");
    std::map<int,std::string> parameter_ids;{std::ifstream table(input.members_path);std::string row;while(std::getline(table,row)){std::istringstream parsed(row);int id;if(parsed>>id)parameter_ids[id]=std::to_string(journal_hash(row));}}
    auto start=std::chrono::steady_clock::now();auto seconds=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    LodLimits limits;limits.maximum_reference_nodes=input.cap;limits.threads=input.threads;limits.maximum_patch_entries=integer("maximum-patch-entries",1,1000000000);limits.maximum_dense_entries=integer("maximum-dense-entries",1,1000000000);
    if(snapshot.fine.mesh.nodes.size()*static_cast<std::size_t>(batch)>limits.maximum_dense_entries)throw std::runtime_error("audit batch exceeds dense allocation bound");
    std::unique_ptr<LodSpace> space;std::unique_ptr<AdditiveKernelRieszContext> riesz;AdjointTestCache aot;std::unique_ptr<RegionalEvaluator> coupled;
    lod2d::helmholtz::HelmholtzOperators afem_operators;
    if(!afem){space=std::make_unique<LodSpace>(snapshot.coarse.mesh,snapshot.reference(),16,snapshot.ell,input.policy=="area"?InterpolationPolicy::ManuscriptAreaWeighted:InterpolationPolicy::ArchivedArithmetic,limits);
        if(space->identity()!=snapshot.space_identity)throw std::invalid_argument("checkpoint/audit operator or space identity mismatch");
        if(input.problem=="E2"){riesz=std::make_unique<AdditiveKernelRieszContext>(*space,input.riesz_policy=="n2"?RieszPatchPolicy::ManuscriptN2:RieszPatchPolicy::ArchivedSupportExpanded);coupled=std::make_unique<RegionalEvaluator>(*space,*riesz,aot);}}
    if(afem)afem_operators=lod2d::helmholtz::assemble_helmholtz_operators(snapshot.fine.mesh,16);
    const auto& operators=afem?afem_operators:space->operators();Sparse E=operators.stiffness+256*operators.mass;
    ReferenceFemContext reference(operators);double prepare_seconds=seconds(),load_seconds=0,error_seconds=0,coupled_seconds=0,fresh_seconds=0;
    auto q=paper_quadrature(input.problem);std::cout<<std::setprecision(17);int completed=0;
    for(int begin=0;begin<static_cast<int>(selected.size());begin+=batch){int count=std::min(batch,static_cast<int>(selected.size())-begin);std::vector<Problem> problems;for(int j=0;j<count;++j)problems.push_back(selected[begin+j].problem);
        double t=seconds();auto loads=assemble_load_batch(snapshot.fine.mesh,problems,q,input.threads);load_seconds+=seconds()-t;
        auto ref=reference.solve(loads);ComplexMatrix values,base;double pg=0;
        t=seconds();if(afem)values=ref;else if(coupled){auto result=coupled->evaluate(loads,snapshot.phi,riesz->regional_mask(radius));values=result.values;pg=result.pg_residual;if(rankzero)base=space->solve(loads).values;}
        else {auto result=space->solve(loads);values=result.values;pg=result.pg_relative_residual;if(rankzero)base=values;}coupled_seconds+=seconds()-t;
        t=seconds();auto errors=integrate_audit_batch(snapshot.fine.mesh,E,values,ref,problems,q,input.threads);
        const auto& floors=errors.reference_error;error_seconds+=seconds()-t;
        for(int j=0;j<count;++j){const auto& member=selected[begin+j];double n=errors.exact_norm[j];auto norm=[&](const ComplexVector& x){return std::sqrt(std::max(0.,x.dot(E.cast<Complex>()*x).real()));};
            double gap=afem?0:norm((ref.col(j)-values.col(j)).eval());
            std::cout<<"{\"kind\":\"sample\",\"state_id\":"<<snapshot.cursor.state_id<<",\"sample\":"<<member.id<<",\"role\":"<<json_string(member.role)<<",\"nominal\":"<<(member.id==0?"true":"false")
                <<",\"marking_member\":"<<(std::find(input.ids.begin(),input.ids.end(),member.id)!=input.ids.end()?"true":"false")<<",\"ell\":";if(afem)std::cout<<"null";else std::cout<<snapshot.ell;std::cout<<",\"rank\":"<<snapshot.phi.cols()<<",\"energy\":"<<errors.energy[j]
                <<",\"exact_norm\":"<<n<<",\"exact_error\":"<<errors.exact_error[j]<<",\"e\":"<<errors.exact_error[j]<<",\"f\":";if(afem)std::cout<<"null";else std::cout<<floors[j];
            std::cout<<",\"g\":";if(afem)std::cout<<"null";else std::cout<<gap;
            auto ratio=[&](double x){if(n>1e-12)std::cout<<x/n;else std::cout<<"null";};
            std::cout<<",\"E\":";ratio(errors.exact_error[j]);std::cout<<",\"E_ref\":";if(afem)std::cout<<"null";else ratio(floors[j]);std::cout<<",\"F\":";if(afem)std::cout<<"null";else ratio(floors[j]);std::cout<<",\"G\":";if(afem)std::cout<<"null";else ratio(gap);
            std::cout<<",\"parameter_id\":"<<json_string(parameter_ids.at(member.id))<<",\"reference_mesh_hash\":"<<json_string(mesh_fingerprint(snapshot.fine.mesh))<<",\"dictionary_hash\":"<<json_string(matrix_hash(snapshot.phi))<<",\"solution_hash\":"<<json_string(matrix_hash(values.col(j)));
            std::cout<<",\"ratio_status\":\""<<(n>1e-12?"finite":"near_zero_exact_norm")<<"\",\"reference_status\":\""<<(afem?"not_applicable":"available")<<"\",\"PG_residual\":"<<(afem?reference.relative_residual():pg);
            if(rankzero&&!afem)std::cout<<",\"base_gap\":"<<norm((ref.col(j)-base.col(j)).eval())<<",\"as_correction_energy\":"<<norm((values.col(j)-base.col(j)).eval());
            if(fresh&&!afem){t=seconds();ComplexVector u;int rank=0;if(input.problem=="E2"){auto trained=train_regional(*space,*riesz,aot,loads.col(j),{radius,rankcap,false});u=trained.accepted.values.col(0);rank=trained.phi.cols();}else u=space->solve(loads.col(j)).values.col(0);
                double error=lod2d::helmholtz::compute_helmholtz_error(snapshot.fine.mesh,u,16,member.problem.exact,member.problem.exact_gradient,q,member.problem.quadrature_context).energy;fresh_seconds+=seconds()-t;
                std::cout<<",\"fresh_rank\":"<<rank<<",\"fresh_e\":"<<error<<",\"fresh_g\":"<<norm((ref.col(j)-u).eval())<<",\"shared_fresh_signed_percent\":";if(error>1e-12)std::cout<<100*(errors.exact_error[j]-error)/error;else std::cout<<"null";}
            std::cout<<"}\n";++completed;
        }
    }
    std::cout<<"{\"kind\":\"audit_complete\",\"state_id\":"<<snapshot.cursor.state_id<<",\"samples\":"<<completed<<",\"batch_size\":"<<batch<<",\"reference_factorizations\":1,\"aot_factorizations\":"<<aot.factorizations()
        <<",\"prepare_seconds\":"<<prepare_seconds<<",\"load_seconds\":"<<load_seconds<<",\"error_seconds\":"<<error_seconds<<",\"reference_factor_seconds\":"<<reference.factor_seconds()<<",\"reference_solve_seconds\":"<<reference.solve_seconds()
        <<",\"coupled_seconds\":"<<coupled_seconds<<",\"fresh_seconds\":"<<fresh_seconds<<",\"wall_seconds\":"<<seconds()<<"}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<"alod_audit: "<<e.what()<<'\n';return 1;}}
