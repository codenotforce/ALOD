#include "alod/patch_cache.hpp"
#include "alod/batch.hpp"
#include "alod/afem.hpp"
#include "alod/checkpoint.hpp"
#include "mesh/refine.h"
#include <fstream>
#include <iostream>
#include <chrono>
#include <atomic>
using namespace alod;
using namespace lod2d;
using namespace lod2d::helmholtz;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main(){try{
    omp_set_num_threads(2);
    for(const auto* id:{"E1","E2"})for(double k:{8.,128.}){
        auto problem=make_problem(id,k);auto mesh=refine_mesh_nvb(problem.initial_mesh,4).mesh;
        auto q=paper_quadrature(id);std::vector<Problem> problems{problem,problem};
        std::vector<SourceMomentData> moments;
        auto load=assemble_load_batch(mesh,problems,q,2,&moments);
        require(moments.size()==2&&moments[0].mass==moments[1].mass,"source moments not shared");
        require(load==assemble_load_batch(mesh,problems,q,2),"moment load changed");
        auto ops=assemble_helmholtz_operators(mesh,k);
        adaptive::diagnostics::ResidualMeshContext context(mesh);
        ComplexVector value(mesh.nodes.size());for(int i=0;i<value.size();++i)value[i]=problem.exact(mesh.nodes[i]);
        auto expected=context.estimate(ops,value,load.col(0),problem.source,q,problem.quadrature_context);
        auto actual=context.estimate(ops,value,load.col(0),problem.source,q,problem.quadrature_context,&moments[0]);
        require(std::abs(expected.eta-actual.eta)<1e-10*std::max(1.,expected.eta),"moment estimator changed");
        for(int e=0;e<mesh.elems.size();++e)require(std::abs(expected.element_squared[e]-actual.element_squared[e])<1e-10*std::max(1.,expected.element_squared[e]),"moment element changed");
        moments[0].mesh_identity="wrong";bool rejected=false;
        try{context.estimate(ops,value,load.col(0),problem.source,q,problem.quadrature_context,&moments[0]);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"wrong moment geometry accepted");
        std::atomic<int> calls{0};auto constant=problem;constant.source=[&](const Point2&){++calls;return Complex(-k*k,0);};
        std::vector<Problem> constants{constant};auto rhs=assemble_load_batch(mesh,constants,q,2,&moments);calls=0;
        auto tiny=context.estimate(ops,ComplexVector::Ones(value.size()),rhs.col(0),constant.source,q,constant.quadrature_context,&moments[0]);
        require(calls>0,"cancellation did not use direct quadrature");
    }
    auto problem=make_problem("E1");auto H=refine_mesh_nvb(problem.initial_mesh,3).mesh;auto h=refine_mesh_nvb(H,2);
    LodLimits limits;limits.threads=2;limits.patch_cache=std::make_shared<LodPatchCache>(8*1024*1024);
    HelmholtzPatchSystem signature;signature.local_vertices={0,1};signature.wavenumber=16;
    signature.helmholtz.resize(2,2);signature.helmholtz.setIdentity();signature.constraints=Eigen::MatrixXd::Ones(1,2);signature.rhs=ComplexMatrix::Ones(2,3);
    const auto original_key=LodPatchCache::key(signature,h.mesh);
    signature.constraints(0,0)=2;require(LodPatchCache::key(signature,h.mesh)!=original_key,"constraint change not invalidated");signature.constraints(0,0)=1;
    signature.rhs(0,0)=2;require(LodPatchCache::key(signature,h.mesh)!=original_key,"RHS change not invalidated");signature.rhs(0,0)=1;
    signature.helmholtz.coeffRef(0,0)=2;require(LodPatchCache::key(signature,h.mesh)!=original_key,"operator change not invalidated");
    LodPatchCache disabled(0);require(!disabled.fits(signature),"zero cache budget ignored");
    LodSpace first(H,h,16,1,InterpolationPolicy::ManuscriptAreaWeighted,limits);
    auto misses=limits.patch_cache->misses();
    LodSpace second(H,h,16,1,InterpolationPolicy::ManuscriptAreaWeighted,limits);
    require(limits.patch_cache->hits()>0&&limits.patch_cache->misses()==misses,"identical patch cache miss");
    require((first.trial()-second.trial()).norm()==0.,"cached correctors changed");
    LodSpace changed(H,h,17,1,InterpolationPolicy::ManuscriptAreaWeighted,limits);
    require(limits.patch_cache->misses()>misses&&limits.patch_cache->bytes()<=8*1024*1024,"cache invalidation/budget failed");
    LodSpace restored(H,h,16,1,InterpolationPolicy::ManuscriptAreaWeighted,limits,first.trial(),&first.reduced());
    require((restored.reduced()-first.reduced()).norm()==0.,"restored reduced matrix changed");
    auto bad=first.reduced();bad.coeffRef(0,0)+=1.;bool rejected=false;
    try{LodSpace wrong(H,h,16,1,InterpolationPolicy::ManuscriptAreaWeighted,limits,first.trial(),&bad);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"wrong reduced operator accepted");
    auto root=std::filesystem::temp_directory_path()/("alod-reuse-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    Checkpoint state;state.coarse=MeshState(H);state.fine=MeshState(h.mesh);state.P_node=h.P_node;state.P_elem=h.P_elem;state.P_dg=h.P_dg;
    state.values=ComplexMatrix::Zero(h.mesh.nodes.size(),1);state.computed_ids={0};state.config_json="{}";
    auto file=save_checkpoint(root,state,&first.trial(),&first.reduced(),true);
    auto saved=load_checkpoint(file);require(saved.format_version==3&&!saved.geometry_file.empty(),"shared checkpoint missing geometry");
    state.cursor.state_id=1;auto next=save_checkpoint(root,state,&first.trial(),&first.reduced(),true);
    int blobs=0;for(const auto& entry:std::filesystem::directory_iterator(root/"meshes"))++blobs;
    require(blobs==1,"duplicate geometry written");
    auto blob=root/"meshes"/saved.geometry_file;
    {std::fstream stream(blob,std::ios::in|std::ios::out|std::ios::binary);stream.seekg(20);char c;stream.read(&c,1);c^=1;stream.seekp(20);stream.write(&c,1);}
    rejected=false;try{inspect_checkpoint(next);}catch(const std::exception&){rejected=true;}
    require(rejected,"corrupt shared geometry accepted");
    state.cursor.state_id=2;auto repaired=save_checkpoint(root,state,&first.trial(),&first.reduced(),true);
    require(load_checkpoint(repaired).values.rows()==h.mesh.nodes.size(),"geometry repair failed");
    require(std::filesystem::exists(root/"meshes/recovery"),"corrupt geometry evidence lost");
    saved=load_checkpoint(file);auto portable=save_checkpoint(root/"portable",saved);
    require(load_checkpoint(portable).geometry_file.empty(),"portable checkpoint has dependencies");
    std::filesystem::remove_all(root);
    std::cout<<"Patch cache, source cancellation, reduced operator and shared geometry validated\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
