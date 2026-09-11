#include "alod/meshes.hpp"
#include "helmholtz/boundary.h"
#include "mesh/refine.h"
#include <cmath>
namespace lod2d::helmholtz {
TriMesh make_helmholtz_unit_square_mesh() {
    TriMesh mesh;
    mesh.nodes = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    mesh.elems = {{0, 1, 3}, {2, 3, 1}};
    // The generic unit-square mesh and the R2 cases use a pure impedance
    // boundary. R1 replaces these tags with its D/N/R partition in the paper
    // case factory. Keep the classification on the mesh itself so every
    // downstream assembly path consumes the same explicit edge contract.
    tag_all_boundary_edges(mesh, BoundaryTag::Robin);
    return mesh;
}

TriMesh make_helmholtz_l_shape_mesh() {
    TriMesh mesh;
    mesh.nodes = {
        {-1, -1}, {0, -1},
        {-1, 0}, {0, 0}, {1, 0},
        {-1, 1}, {0, 1}, {1, 1}};
    // Three unit squares, each using the same compatible two-triangle NVB
    // pattern as make_helmholtz_unit_square_mesh(). The lower-right square is
    // removed.
    mesh.elems = {
        {0, 1, 2}, {3, 2, 1},
        {2, 3, 5}, {6, 5, 3},
        {3, 4, 6}, {7, 6, 4}};
    const auto [edges, boundary] = compute_edges(mesh);
    for (std::size_t index = 0; index < edges.size(); ++index) {
        if (!boundary[index]) continue;
        const Edge edge = edges[index];
        const Point2 midpoint = 0.5 * (mesh.nodes[edge[0]] + mesh.nodes[edge[1]]);
        const bool reentrant_edge =
            (std::abs(midpoint.x()) < 1e-14 && midpoint.y() < 0.0)
            || (std::abs(midpoint.y()) < 1e-14 && midpoint.x() > 0.0);
        mesh.boundary_edges.push_back({
            edge,
            reentrant_edge ? BoundaryTag::Dirichlet : BoundaryTag::Robin});
    }
    synchronize_dirichlet_nodes(mesh);
    validate_boundary_tags(mesh);
    return mesh;
}
}
