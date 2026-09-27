#include "alod/batch.hpp"
#include "mesh/refine.h"
#include <atomic>
#include <iostream>
using namespace alod;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main(){try{
    auto problem=make_problem("E1");auto mesh=lod2d::refine_mesh_nvb(problem.initial_mesh,5).mesh;
    auto refined=lod2d::bisect_newest_vertex(mesh,{0});auto q=paper_quadrature("E1");
    std::atomic<int> evaluations=0;auto source=problem.source;
    problem.source=[&](const lod2d::Point2& x){++evaluations;return source(x);};
    std::vector<Problem> family{problem};std::vector<lod2d::helmholtz::SourceMomentData> moments,cold;
    SourceMomentReuse reuse;
    assemble_load_batch(mesh,family,q,2,&moments,&reuse,"frozen-family-1");evaluations=0;
    auto reused=assemble_load_batch(refined.mesh,family,q,2,&moments,&reuse,"frozen-family-1");int reused_calls=evaluations;evaluations=0;
    auto direct=assemble_load_batch(refined.mesh,family,q,2,&cold);int direct_calls=evaluations;
    require(reused_calls<direct_calls,"local refinement did not reuse integrals");
    require((reused-direct).norm()==0,"incremental load changed arithmetic");
    for(int e=0;e<static_cast<int>(refined.mesh.elems.size());++e){
        require(moments[0].squared[e]==cold[0].squared[e],"source squared integral changed");
        require((*moments[0].mass)[e]==(*cold[0].mass)[e],"quadrature mass changed");
    }
    ExactIntegrationReuse exact;
    auto op=lod2d::helmholtz::assemble_helmholtz_operators(mesh,16.);
    Sparse energy=op.stiffness+256.*op.mass;ComplexMatrix values=ComplexMatrix::Zero(mesh.nodes.size(),1);
    integrate_error_batch(mesh,energy,values,family,q,2,&exact,"frozen-family-1");
    auto fine_op=lod2d::helmholtz::assemble_helmholtz_operators(refined.mesh,16.);energy=fine_op.stiffness+256.*fine_op.mass;
    values=ComplexMatrix::Constant(refined.mesh.nodes.size(),1,Complex(.17,.09));
    auto cached=integrate_error_batch(refined.mesh,energy,values,family,q,2,&exact,"frozen-family-1");
    auto uncached=integrate_error_batch(refined.mesh,energy,values,family,q,2);
    require(exact.entries.hits()>0,"unchanged triangles did not reuse exact jets");
    require((cached.exact_error-uncached.exact_error).norm()==0,"cached jets reused stale discrete values");
    require((cached.exact_norm-uncached.exact_norm).norm()==0,"cached exact norm changed");
    auto hits=exact.entries.hits();auto boosted=q;boosted.base_triangle_order+=2;
    integrate_error_batch(refined.mesh,energy,values,family,boosted,2,&exact,"frozen-family-1");
    require(exact.entries.hits()==hits,"changed quadrature reused stale jets");
    std::cout<<"Source evaluations "<<direct_calls<<" -> "<<reused_calls<<"; exact-jet hits "<<hits<<"\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
