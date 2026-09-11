#pragma once
#include "mesh/refine.h"
#include <bit>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace alod {
// Identities refer to accepted refinement transactions. Closure can bisect a
// parent more than once; generation records its total NVB depth increment.
struct ElementIdentity { std::uint64_t id, parent; int generation; };
struct MeshState {
    lod2d::TriMesh mesh;
    std::vector<ElementIdentity> elements, ancestry;
    std::uint64_t next_id=0, version=0;
    explicit MeshState(lod2d::TriMesh initial):mesh(std::move(initial)) {
        for(std::size_t i=0;i<mesh.elems.size();++i) elements.push_back({next_id++,UINT64_MAX,0});
        ancestry=elements;
    }
    void refine(const std::vector<int>& marks) {
        auto out=lod2d::bisect_newest_vertex(mesh,marks);
        adopt(std::move(out));
    }
    void adopt(lod2d::RefineOutput out) {
        auto before=lod2d::compute_area(mesh), after=lod2d::compute_area(out.mesh);
        std::vector<int> parents(out.mesh.elems.size(),-1);
        for(int c=0;c<out.P_elem.outerSize();++c)
            for(Eigen::SparseMatrix<double>::InnerIterator it(out.P_elem,c);it;++it) {
                if(parents[it.row()]!=-1 || it.value()!=1) throw std::runtime_error("invalid NVB parent map");
                parents[it.row()]=c;
            }
        std::vector<ElementIdentity> identities;
        for(std::size_t i=0;i<parents.size();++i) {
            int p=parents[i]; if(p<0) throw std::runtime_error("missing NVB parent");
            if(out.mesh.elems[i]==mesh.elems[p]) identities.push_back(elements[p]);
            else {
                int depth=static_cast<int>(std::llround(std::log2(before[p]/after[i])));
                if(depth<=0 || std::abs(std::ldexp(after[i],depth)-before[p])>1e-12*before[p])
                    throw std::runtime_error("invalid child area");
                ElementIdentity id{next_id++,elements[p].id,elements[p].generation+depth};
                identities.push_back(id);ancestry.push_back(id);
            }
        }
        if(out.mesh.elems!=mesh.elems) ++version;
        mesh=std::move(out.mesh);elements=std::move(identities);
    }
};

// Portable diagnostic fingerprint (FNV-1a; not a cryptographic checksum).
inline std::string mesh_fingerprint(const lod2d::TriMesh& mesh) {
    std::uint64_t h=14695981039346656037ULL;
    auto word=[&](std::uint64_t x){for(int b=0;b<8;++b){h^=(x>>(8*b))&255;h*=1099511628211ULL;}};
    word(mesh.nodes.size());word(mesh.elems.size());word(mesh.boundary_edges.size());
    for(const auto& x:mesh.nodes){word(std::bit_cast<std::uint64_t>(x.x()));word(std::bit_cast<std::uint64_t>(x.y()));}
    for(const auto& t:mesh.elems)for(int v:t)word(v);
    for(const auto& e:mesh.boundary_edges){word(e.nodes[0]);word(e.nodes[1]);word(static_cast<int>(e.tag));}
    std::ostringstream out;out<<std::hex<<std::setw(16)<<std::setfill('0')<<h;return out.str();
}
}
