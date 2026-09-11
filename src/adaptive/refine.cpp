#include "alod/adaptive.hpp"
#include <cmath>
#include <stdexcept>
namespace alod {
using namespace lod2d;
MeshTransition refine_pair(const TriMesh& H,const RefineOutput& h,
    const std::vector<int>& hm,const std::vector<int>& fm,int gap,int layers,int cap,bool commit_coarse) {
    if(gap<0||layers<0||cap<1)throw std::invalid_argument("invalid refinement limits");
    auto coarse=bisect_newest_vertex(H,hm);
    auto fine=bisect_newest_vertex(h.mesh,fm);
    auto check=[&]{if(fine.mesh.nodes.size()>static_cast<std::size_t>(cap))throw std::runtime_error("reference node limit exceeded");};
    check();
    Sparse injection=fine.P_node;
    Sparse reference_parent=fine.P_elem;
    Sparse old_parent=fine.P_elem*h.P_elem;
    auto parent=fine_element_parents(coarse.P_elem,coarse.mesh.elems.size(),H.elems.size());
    std::vector<int> count(H.elems.size());for(int p:parent)++count[p];
    std::vector<char> active(parent.size(),false);
    for(int i=0;i<static_cast<int>(parent.size());++i)active[i]=count[parent[i]]>1;
    // Legacy element rings share a vertex, including conformity closure seeds.
    for(int ring=0;ring<layers;++ring){
        std::vector<char> nodes(coarse.mesh.nodes.size(),false);
        for(int i=0;i<static_cast<int>(active.size());++i)if(active[i])for(int v:coarse.mesh.elems[i])nodes[v]=true;
        for(int i=0;i<static_cast<int>(active.size());++i)for(int v:coarse.mesh.elems[i])if(nodes[v])active[i]=true;
    }
    auto refine=[&](const std::vector<int>& marks){
        // NVB closure can allocate more than the requested bisections. Also
        // check the actual result before any operator or dense block is built.
        if(fine.mesh.nodes.size()+marks.size()/2>static_cast<std::size_t>(cap))throw std::runtime_error("reference refinement budget exceeded");
        auto next=bisect_newest_vertex(fine.mesh,marks);
        injection=next.P_node*injection;old_parent=next.P_elem*old_parent;
        reference_parent=next.P_elem*reference_parent;
        fine=std::move(next);check();
    };
    RefineOutput embedding;int rounds=0;
    for(;;){
        auto old=fine_element_parents(old_parent,fine.mesh.elems.size(),H.elems.size());
        try {embedding=update_nested_mesh_embedding_after_parent_refinement(H,coarse,fine.mesh,old);break;}
        catch(const std::invalid_argument& e){
            if(std::string(e.what())!="refined parent mesh is not contained in the fixed child mesh")throw;
        }
        if(++rounds>32)throw std::runtime_error("nested closure failed after 32 rounds");
        std::vector<int> marks;for(int i=0;i<static_cast<int>(old.size());++i)if(count[old[i]]>1)marks.push_back(i);
        if(marks.empty())throw std::runtime_error("nested closure has no active children");
        refine(marks);
    }
    const auto Harea=compute_area(coarse.mesh);
    for(int round=0;;++round){
        const auto farea=compute_area(fine.mesh);
        const auto p=fine_element_parents(embedding.P_elem,fine.mesh.elems.size(),coarse.mesh.elems.size());
        std::vector<int> marks;
        for(int i=0;i<static_cast<int>(p.size());++i)
            if(active[p[i]]&&std::llround(std::log2(Harea[p[i]]/farea[i]))<gap)marks.push_back(i);
        if(marks.empty())break;
        if(round>=32)throw std::runtime_error("reference gap closure failed");
        refine(marks);++rounds;
        embedding=update_nested_mesh_embedding_after_parent_refinement(H,coarse,fine.mesh,
            fine_element_parents(old_parent,fine.mesh.elems.size(),H.elems.size()));
    }
    if(!commit_coarse){
        // E1's archived reference-only sweep closes over the proposed coarse
        // mesh, then discards that proposal. Retain its reference closure.
        embedding.P_node=embedding.P_node*coarse.P_node;
        embedding.P_elem=embedding.P_elem*coarse.P_elem;
        embedding.P_dg=embedding.P_dg*coarse.P_dg;
        Sparse identity_nodes(H.nodes.size(),H.nodes.size()),identity_elements(H.elems.size(),H.elems.size());
        identity_nodes.setIdentity();identity_elements.setIdentity();
        return {H,std::move(embedding),std::move(injection),rounds,std::move(identity_nodes),std::move(identity_elements),std::move(reference_parent)};
    }
    return {std::move(coarse.mesh),std::move(embedding),std::move(injection),rounds,std::move(coarse.P_node),std::move(coarse.P_elem),std::move(reference_parent)};
}
}
