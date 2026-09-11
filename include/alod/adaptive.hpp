#pragma once
#include "alod/localization.hpp"
#include <set>

namespace alod {
struct AdaptiveCursor {
    int state_id=0, coarse_cycle=0, reference_sweep=0;
    void advance(int m_ref);
    bool cycle_complete(int m_ref) const;
};
enum class EllMode { Lazy, EveryState, Fixed };
struct EllPolicy {
    std::string problem="E1";
    EllMode mode=EllMode::Lazy;
    int maximum=4, last_check=-1;
    std::set<int> extra_checks;
    bool due(int state_id, int ell, bool terminal) const;
    // E1 promotes once per scheduled check. E2 repeats at the new ell.
    std::string decision(double theta, double eta, int ell) const;
};
lod2d::RefineOutput build_nested_mesh_embedding(const lod2d::TriMesh&,const lod2d::TriMesh&);
std::vector<int> fine_element_parents(const Sparse&,int,int);
lod2d::RefineOutput update_nested_mesh_embedding_after_parent_refinement(
    const lod2d::TriMesh&,const lod2d::RefineOutput&,const lod2d::TriMesh&,const std::vector<int>&);
struct MeshTransition {
    lod2d::TriMesh coarse;
    lod2d::RefineOutput reference;
    Sparse injection; // Old reference values to the new reference.
    int closure_rounds=0;
};
// Both sets of marks must come from the same accepted solution and frozen scales.
MeshTransition refine_pair(const lod2d::TriMesh&,const lod2d::RefineOutput&,
    const std::vector<int>& coarse_marks,const std::vector<int>& reference_marks,
    int minimum_gap,int active_layers,int maximum_nodes,bool commit_coarse=true);
}
