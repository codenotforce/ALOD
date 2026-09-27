#include "audit_diagnostics.hpp"
#include "mesh/refine.h"

namespace audit_diagnostics {
using namespace alod;

void deepen_reference(std::unique_ptr<LodSpace>& space, Checkpoint& snapshot,
                      int levels, int ell, double wavenumber,
                      InterpolationPolicy policy, const LodLimits& limits) {
    long double bound = space->fine().nodes.size();
    long double triangles = space->fine().elems.size();
    for (int level = 0; level < levels; ++level) {
        bound += triangles;
        triangles *= 2;
    }
    if (bound > limits.maximum_reference_nodes)
        throw std::runtime_error("reference diagnostic exceeds conservative node limit");
    auto fine = lod2d::refine_mesh_nvb(space->fine(), levels);
    if (fine.mesh.nodes.size() > limits.maximum_reference_nodes)
        throw std::runtime_error("reference diagnostic exceeds node limit");

    // Borrow the old embeddings until construction completes; copying them can
    // be much larger than the diagnostic's new coarse matrix.
    const auto& previous = space->hierarchy()->reference;
    lod2d::RefineOutput lifted{fine.mesh, fine.P_node * previous.P_node,
                              fine.P_elem * previous.P_elem, fine.P_dg * previous.P_dg};
    snapshot.phi = (fine.P_node.cast<Complex>() * snapshot.phi).eval();
    snapshot.values = (fine.P_node.cast<Complex>() * snapshot.values).eval();
    ComplexSparseMatrix basis = fine.P_node.cast<Complex>() * space->trial();
    // Injection preserves the frozen space. Shared-family retraining is a
    // separate explicit diagnostic and must not be hidden in mesh deepening.
    space = std::make_unique<LodSpace>(space->coarse(), std::move(lifted),
                                      wavenumber, ell, policy, limits, basis);
}

double two_level_corrector(const LodSpace& space, int ell) {
    if (ell >= 4) throw std::invalid_argument("two-level diagnostic needs ell < 4");
    LodSpace other(space, ell + 1);
    ComplexMatrix difference = ComplexMatrix(other.trial() - space.trial());
    double square = std::real(
        (difference.adjoint() * space.energy().cast<Complex>() * difference).trace());
    return std::sqrt(std::max(0., square));
}

lod2d::TriMesh centroid_region(const lod2d::TriMesh& mesh, double radius) {
    lod2d::TriMesh region;
    region.nodes = mesh.nodes;
    for (auto tri : mesh.elems) {
        auto center = (mesh.nodes[tri[0]] + mesh.nodes[tri[1]] + mesh.nodes[tri[2]]) / 3.;
        if (center.norm() <= radius) region.elems.push_back(tri);
    }
    return region;
}
}
