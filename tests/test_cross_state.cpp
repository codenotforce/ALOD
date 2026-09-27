#include "alod/patch_cache.hpp"
#include "alod/local_factor_cache.hpp"
#include "alod/adaptive.hpp"
#include "alod/problems.hpp"
#include "alod/estimator.hpp"
#include <iostream>
#include <future>
using namespace alod;
void require(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
int main(){try{
    auto problem=make_problem("E1");auto H=lod2d::refine_mesh_nvb(problem.initial_mesh,6).mesh;
    auto h=lod2d::refine_mesh_nvb(H,2);
    LodLimits limits;limits.threads=2;
    limits.patch_cache=std::make_shared<LodPatchCache>(64ULL*1024*1024);
    limits.riesz_cache=std::make_shared<LocalFactorCache>(64ULL*1024*1024);
    LodSpace first(H,h,16,1,InterpolationPolicy::ManuscriptAreaWeighted,limits);
    AdditiveKernelRieszContext first_riesz(first);
    const auto start_hits=limits.patch_cache->hits(),start_riesz=limits.riesz_cache->hits();
    // A genuine local mesh transition, not reconstruction of the same state.
    auto next=refine_pair(H,h,{0},{0},2,3,20000);
    LodSpace cached(next.coarse,next.reference,16,1,InterpolationPolicy::ManuscriptAreaWeighted,limits);
    AdditiveKernelRieszContext cached_riesz(cached);
    const auto patch_hits=limits.patch_cache->hits()-start_hits,riesz_hits=limits.riesz_cache->hits()-start_riesz;
    require(patch_hits>0,"no corrector reuse across local refinement");
    require(riesz_hits>0,"no Riesz factor reuse across local refinement");
    auto cold_limits=limits;cold_limits.patch_cache.reset();cold_limits.riesz_cache.reset();
    LodSpace cold(next.coarse,next.reference,16,1,InterpolationPolicy::ManuscriptAreaWeighted,cold_limits);
    AdditiveKernelRieszContext cold_riesz(cold);
    require((cached.trial()-cold.trial()).norm()<1e-11*std::max(1.,cold.trial().norm()),"cached refined basis differs");
    ComplexMatrix rhs(cold.fine().nodes.size(),3);
    for(int j=0;j<3;++j)for(int i=0;i<rhs.rows();++i)rhs(i,j)=Complex(std::sin((i+1.)*(j+1)),std::cos((i+2.)*(j+3)));
    auto a=cached_riesz.apply(rhs),b=cold_riesz.apply(rhs);
    require((a.values-b.values).norm()<1e-11*std::max(1.,b.values.norm()),"reused Riesz factors changed action");
    AdditiveKernelRieszContext concurrent_riesz(cached);
    auto left=std::async(std::launch::async,[&]{return cached_riesz.apply(rhs);});
    auto right=std::async(std::launch::async,[&]{return concurrent_riesz.apply(rhs);});
    require((left.get().values-b.values).norm()<1e-11*std::max(1.,b.values.norm())
        &&(right.get().values-b.values).norm()<1e-11*std::max(1.,b.values.norm()),"concurrent immutable factor solves changed");
    // Reverse global fine-node numbering, preserving physical elements.
    auto renumbered=next.reference;const int n=renumbered.mesh.nodes.size();
    std::reverse(renumbered.mesh.nodes.begin(),renumbered.mesh.nodes.end());
    for(auto& e:renumbered.mesh.elems)for(int& v:e)v=n-1-v;
    for(int& v:renumbered.mesh.dirichlet)v=n-1-v;
    for(auto& e:renumbered.mesh.boundary_edges){for(int& v:e.nodes)v=n-1-v;std::sort(e.nodes.begin(),e.nodes.end());}
    std::vector<Eigen::Triplet<double>> entries;
    for(int j=0;j<next.reference.P_node.outerSize();++j)for(Sparse::InnerIterator it(next.reference.P_node,j);it;++it)entries.emplace_back(n-1-it.row(),it.col(),it.value());
    renumbered.P_node.setZero();renumbered.P_node.setFromTriplets(entries.begin(),entries.end());
    const auto before=limits.patch_cache->hits();
    LodSpace shuffled(next.coarse,renumbered,16,1,InterpolationPolicy::ManuscriptAreaWeighted,limits);
    require(limits.patch_cache->hits()>before,"renumbering invalidated every corrector");
    ComplexMatrix ordered=ComplexMatrix(shuffled.trial()).colwise().reverse().eval();
    require((ordered-ComplexMatrix(cached.trial())).norm()<1e-10*std::max(1.,cached.trial().norm()),"cached row mapping is wrong");
    const auto renumbering_hits=limits.patch_cache->hits()-before;
    auto translated_H=next.coarse;auto translated_h=next.reference;
    for(auto& p:translated_H.nodes)p+=lod2d::Point2(4.,4.);
    for(auto& p:translated_h.mesh.nodes)p+=lod2d::Point2(4.,4.);
    const auto translation_start=limits.patch_cache->hits();
    LodSpace translated(translated_H,translated_h,16,1,InterpolationPolicy::ManuscriptAreaWeighted,limits);
    LodSpace translated_cold(translated_H,translated_h,16,1,InterpolationPolicy::ManuscriptAreaWeighted,cold_limits);
    require(limits.patch_cache->hits()>translation_start,"equivalent translated patch not reused");
    require((translated.trial()-translated_cold.trial()).norm()<1e-10*std::max(1.,translated_cold.trial().norm()),"translation reused a different local problem");
    require(limits.patch_cache->bytes()<=64ULL*1024*1024&&limits.riesz_cache->bytes()<=64ULL*1024*1024,"cache budget exceeded");
    std::cout<<"{\"refinement_corrector_hits\":"<<patch_hits<<",\"refinement_riesz_hits\":"<<riesz_hits<<",\"renumbering_hits\":"<<renumbering_hits<<",\"translation_hits\":"<<limits.patch_cache->hits()-translation_start<<"}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
