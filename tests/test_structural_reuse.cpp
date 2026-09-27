#include "alod/regional.hpp"
#include "alod/patch_cache.hpp"
#include "alod/batch.hpp"
#include "alod/checkpoint.hpp"
#include <cstdlib>
#include <chrono>
#include <iostream>
using namespace alod;
void require(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
void reference(bool yes){
#ifdef _WIN32
    _putenv_s("ALOD_REFERENCE_EXECUTION",yes?"1":"");
#else
    if(yes)setenv("ALOD_REFERENCE_EXECUTION","1",1);else unsetenv("ALOD_REFERENCE_EXECUTION");
#endif
}
int main(int argc,char** argv){try{
    const double wavenumber=argc>1?std::stod(argv[1]):32.;
    auto problem=make_problem("E2",wavenumber);auto H=lod2d::refine_mesh_nvb(problem.initial_mesh,5).mesh;
    auto h=lod2d::refine_mesh_nvb(H,2);LodLimits limits;limits.threads=2;
    limits.patch_cache=std::make_shared<LodPatchCache>(32*1024*1024);
    auto owner=std::make_unique<LodSpace>(H,h,wavenumber,1,InterpolationPolicy::ManuscriptAreaWeighted,limits);
    auto hierarchy=owner->hierarchy();AdditiveKernelRieszContext riesz(*owner);
    LodSpace promoted(*owner,2);require(promoted.hierarchy()==hierarchy,"ell rebuild copied immutable hierarchy");
    owner.reset(); // The estimator must remain valid after an ell replacement.
    ComplexMatrix loads(h.mesh.nodes.size(),6);
    for(int j=0;j<loads.cols();++j)for(int i=0;i<loads.rows();++i)loads(i,j)=Complex(std::sin((i+1.)*(j+1)),std::cos((i+2.)*(j+3)));
    auto solved=promoted.solve(loads);
    require((solved.coefficients-promoted.solve_reduced(promoted.test().adjoint()*loads)).norm()==0,"shared coarse solve changed");
    AdjointTestCache cache;
    const auto initial_columns=riesz.applied_columns();
    reference(true);auto before=train_regional(promoted,riesz,cache,loads,{10.,12,true});
    const auto before_columns=riesz.applied_columns()-initial_columns;
    reference(false);auto after=train_regional(promoted,riesz,cache,loads,{10.,12,true});
    const auto after_columns=riesz.applied_columns()-initial_columns-before_columns;
    if(wavenumber==32)require(after_columns<before_columns,"training coordinates did not reduce Riesz RHS work");
    require(before.working_rank==after.working_rank&&before.phi.cols()==after.phi.cols()&&before.stop==after.stop&&before.compression_trials==after.compression_trials,"POD/training decision changed");
    require((before.accepted.values-after.accepted.values).norm()<1e-8*std::max(1.,before.accepted.values.norm()),"POD solution changed");
    require((before.accepted.eta-after.accepted.eta).norm()<1e-8*std::max(1.,before.accepted.eta.norm()),"POD estimator changed");
    auto hits=limits.patch_cache->hits();LodSpace identical(H,h,wavenumber,2,InterpolationPolicy::ManuscriptAreaWeighted,limits);
    require(limits.patch_cache->hits()>hits,"early dependency cache did not hit");
    require((identical.trial()-promoted.trial()).norm()==0,"early hit changed basis");
    auto tight=limits;tight.maximum_dense_entries=1;bool limited=false;
    try{LodSpace rejected(H,h,wavenumber,2,InterpolationPolicy::ManuscriptAreaWeighted,tight);}catch(const std::runtime_error&){limited=true;}
    require(limited,"cache hit bypassed current resource ceiling");
    auto changed=refine_pair(H,h,{0,2},{0,4},2,3,20000);
    LodSpace cached(changed.coarse,changed.reference,wavenumber,2,InterpolationPolicy::ManuscriptAreaWeighted,limits);
    LodLimits cold_limits=limits;cold_limits.patch_cache.reset();
    LodSpace cold(changed.coarse,changed.reference,wavenumber,2,InterpolationPolicy::ManuscriptAreaWeighted,cold_limits);
    // Matching includes arithmetic order, so reuse preserves the original bits.
    const double reuse_error=(cached.trial()-cold.trial()).norm()/std::max(1.,cold.trial().norm());
    require(reuse_error==0,"refinement dependency invalidation failed");
    std::cout<<"Corrector relative difference="<<reuse_error<<'\n';
    require(limits.patch_cache->bytes()<=32*1024*1024,"dependency storage exceeded bound");
    std::vector<Problem> problems;
    for(int j=0;j<6;++j)problems.push_back(lod2d::helmholtz::benchmarks::make_parameterized_boundary_gaussian_s_paper_case(8.*(j+1),.1*j,.1*j,.3*j,80.,lod2d::Point2(-.4-.03*j,.4+.02*j)));
    auto shared_source=lod2d::helmholtz::benchmarks::shared_source_evaluator(problems);
    auto shared_jet=lod2d::helmholtz::benchmarks::shared_jet_evaluator(problems);
    std::vector<Complex> sources(6);std::vector<lod2d::helmholtz::benchmarks::ExactJet> jets(6);
    lod2d::Point2 point(-.37,.43);shared_source(point,sources.data());shared_jet(point,jets.data());
    for(int j=0;j<6;++j){require(sources[j]==problems[j].source(point),"shared source changed");auto jet=problems[j].exact_jet(point);require(jet.first==jets[j].first&&jet.second==jets[j].second,"shared jet changed");}
    auto custom=problems;custom[0].source=[](const lod2d::Point2&){return Complex(7.,0);};
    require(!lod2d::helmholtz::benchmarks::shared_source_evaluator(custom),"custom callback incorrectly shared");
    auto q=paper_quadrature("E2");
    AuditIntegrationGeometry geometry(promoted.fine(),2);
    auto a=integrate_audit_batch(promoted.fine(),promoted.energy(),after.accepted.values,solved.values,problems,q,2);
    auto b=integrate_audit_batch(promoted.fine(),promoted.energy(),after.accepted.values,solved.values,problems,q,2,true,&geometry);
    require(a.exact_error==b.exact_error&&a.reference_error==b.reference_error,"geometry reuse changed integration");
    MeshState coarse(H),fine(h.mesh);CheckpointGeometryView view{coarse,fine,h.P_node,h.P_elem,h.P_dg};
    Checkpoint cp;cp.config_json="{}";cp.values=ComplexMatrix::Zero(h.mesh.nodes.size(),1);cp.computed_ids={0};auto root=std::filesystem::temp_directory_path()/("alod-structural-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    auto file=save_checkpoint(root,cp,nullptr,nullptr,true,&view);auto read=load_checkpoint(file);
    require(read.fine.mesh.nodes.size()==h.mesh.nodes.size()&&!view.object_name.empty(),"borrowed geometry was not published");
    cp.cursor.state_id=1;auto file2=save_checkpoint(root,cp,nullptr,nullptr,true,&view);
    require(load_checkpoint(file2).geometry_file==read.geometry_file,"geometry identity not reused");
    std::filesystem::remove_all(root);
    std::cout<<"Structural reuse passed; working rank="<<after.working_rank<<", POD trials="<<after.compression_trials<<", Riesz columns="<<before_columns<<" -> "<<after_columns<<'\n';
}catch(const std::exception& e){reference(false);std::cerr<<e.what()<<'\n';return 1;}}
