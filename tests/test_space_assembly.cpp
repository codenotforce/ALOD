#include "alod/lod.hpp"
#include "alod/problems.hpp"
#include "helmholtz/patch_system.h"
#include "helmholtz/boundary.h"
#include "lod/patches.h"
#include <set>
#include <iostream>

int main(){try{
    using namespace alod;
    using namespace lod2d;
    using namespace lod2d::helmholtz;
    for(const auto* id:{"E1","E2"})for(bool implicit:{false,true}){
        auto p=make_problem(id);
        auto H=refine_mesh_nvb(p.initial_mesh,3).mesh;
        H=bisect_newest_vertex(H,{0}).mesh;
        auto h=refine_mesh_nvb(H,2);
        LodLimits limits;limits.threads=2;
        LodSpace s(H,h,16,2,InterpolationPolicy::ManuscriptAreaWeighted,limits);
        LodSpace restored(H,h,16,2,InterpolationPolicy::ManuscriptAreaWeighted,limits,s.trial());
        if(restored.identity()!=s.identity()||(restored.trial()-s.trial()).norm()!=0.)
            throw std::runtime_error("restored LOD basis differs");
        if(implicit)h.mesh.boundary_edges.clear();
        const auto [edges,boundary]=compute_edges(h.mesh);
        std::set<Edge> natural;
        for(std::size_t e=0;e<edges.size();++e)if(boundary[e]){
            auto tag=boundary_tag(h.mesh,edges[e]);
            if(tag==BoundaryTag::Robin||tag==BoundaryTag::Neumann)natural.insert(edges[e]);
        }
        auto patches=build_patches(H,2);
        const std::vector<TriMesh> no_meshes;
        const std::vector<Sparse> no_prolongations;
        HelmholtzPatchAssembler assembler(H,h.mesh,h.P_elem,h.P_dg,s.interpolation(),
            patches,no_meshes,no_prolongations,no_prolongations,s.operators());
        for(int target=0;target<static_cast<int>(H.elems.size());++target){
            auto system=assembler.assemble(target);
            auto compact=assembler.assemble(target,false);
            if((compact.helmholtz-system.helmholtz).norm()!=0.||compact.rhs!=system.rhs||compact.constraints!=system.constraints
                ||compact.stiffness.size()!=0||compact.mass.size()!=0||compact.robin.size()!=0)
                throw std::runtime_error("compact patch system differs");
            bool touches=false;
            for(int e:system.patch_elements)for(int j=0;j<3;++j)
                touches=touches||natural.contains(canonical_edge(h.mesh.elems[e][j],h.mesh.elems[e][(j+1)%3]));
            if(touches!=system.touches_physical_boundary)throw std::runtime_error("cached patch boundary changed");
            std::vector<int> local(h.mesh.nodes.size(),-1);
            for(int i=0;i<static_cast<int>(system.local_vertices.size());++i)local[system.local_vertices[i]]=i;
            ComplexMatrix rhs=ComplexMatrix::Zero(system.rhs.rows(),3);
            for(Sparse::InnerIterator it(h.P_elem,target);it;++it)if(it.value()!=0.){
                int e=it.row();
                for(int i=0;i<3;++i){
                    int row=local[h.mesh.elems[e][i]];if(row<0)continue;
                    for(int c=0;c<3;++c){
                        Complex value=0.;
                        for(int j=0;j<3;++j)value+=s.operators().element_blocks[e](i,j)*h.P_dg.coeff(3*e+j,3*target+c);
                        rhs(row,c)+=value;
                    }
                }
            }
            if(rhs!=system.rhs)throw std::runtime_error("local embedding reuse changed patch RHS");
        }
    }
    std::cout<<"Patch boundary and RHS match independent assembly on graded E1/E2 meshes\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
