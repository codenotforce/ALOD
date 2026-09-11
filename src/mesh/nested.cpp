#include "alod/adaptive.hpp"

#include "helmholtz/boundary.h"
#include "lod/quasi_interp.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace alod {
using namespace lod2d;
using namespace lod2d::helmholtz;
namespace {

Eigen::SparseMatrix<double> identity_sparse(int size) {
    Eigen::SparseMatrix<double> result(size, size);
    result.setIdentity();
    return result;
}

Eigen::SparseMatrix<double> cg_to_dg(const TriMesh &mesh) {
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * mesh.elems.size());
    for (int element = 0; element < static_cast<int>(mesh.elems.size()); ++element) {
        for (int local = 0; local < 3; ++local)
            triplets.emplace_back(3 * element + local, mesh.elems[element][local], 1.0);
    }
    Eigen::SparseMatrix<double> result(
        3 * static_cast<int>(mesh.elems.size()),
        static_cast<int>(mesh.nodes.size()));
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

int refinement_increment(double parent_area, double child_area) {
    if (!(parent_area > 0.0) || !(child_area > 0.0))
        throw std::runtime_error("NVB hierarchy encountered a nonpositive element area");
    const double raw = std::log2(parent_area / child_area);
    const int increment = static_cast<int>(std::llround(raw));
    if (increment < 0 || std::abs(raw - increment) > 1e-9)
        throw std::runtime_error("NVB child area is not a dyadic fraction of its parent");
    return increment;
}

double cross(const Point2 &first, const Point2 &second) {
    return first.x() * second.y() - first.y() * second.x();
}

std::array<double, 3> barycentric_coordinates(
    const Point2 &point,
    const Point2 &first,
    const Point2 &second,
    const Point2 &third) {
    const Point2 edge_one = second - first;
    const Point2 edge_two = third - first;
    const double denominator = cross(edge_one, edge_two);
    if (std::abs(denominator) <= 1e-15)
        throw std::runtime_error("nested mesh contains a degenerate parent triangle");
    const Point2 relative = point - first;
    const double second_weight = cross(relative, edge_two) / denominator;
    const double third_weight = cross(edge_one, relative) / denominator;
    return {1.0 - second_weight - third_weight, second_weight, third_weight};
}

bool inside_triangle(const std::array<double, 3> &weights, double tolerance) {
    return weights[0] >= -tolerance
        && weights[1] >= -tolerance
        && weights[2] >= -tolerance
        && weights[0] <= 1.0 + tolerance
        && weights[1] <= 1.0 + tolerance
        && weights[2] <= 1.0 + tolerance;
}

std::array<double, 3> clean_weights(std::array<double, 3> weights) {
    for (double &weight : weights) {
        if (std::abs(weight) <= 1e-13) weight = 0.0;
        if (std::abs(weight - 1.0) <= 1e-13) weight = 1.0;
    }
    const double sum = weights[0] + weights[1] + weights[2];
    if (std::abs(sum) <= 1e-15)
        throw std::runtime_error("nested mesh interpolation produced zero barycentric weight");
    for (double &weight : weights) weight /= sum;
    return weights;
}

double total_area(const TriMesh &mesh) {
    double result = 0.0;
    for (double area : compute_area(mesh)) result += std::abs(area);
    return result;
}

double triangle_diameter(const TriMesh &mesh, const Triangle &triangle) {
    return std::max({
        (mesh.nodes[triangle[0]] - mesh.nodes[triangle[1]]).norm(),
        (mesh.nodes[triangle[1]] - mesh.nodes[triangle[2]]).norm(),
        (mesh.nodes[triangle[2]] - mesh.nodes[triangle[0]]).norm()});
}

void validate_indices(
    const std::vector<int> &indices,
    int upper_bound,
    const char *description) {
    std::set<int> unique;
    for (int index : indices) {
        if (index < 0 || index >= upper_bound)
            throw std::out_of_range(std::string(description) + " index is out of range");
        if (!unique.insert(index).second)
            throw std::invalid_argument(std::string(description) + " contains a duplicate index");
    }
}

} // namespace

RefineOutput build_nested_mesh_embedding(
    const TriMesh &parent_mesh,
    const TriMesh &child_mesh) {
    if (parent_mesh.nodes.empty() || parent_mesh.elems.empty()
        || child_mesh.nodes.empty() || child_mesh.elems.empty()) {
        throw std::invalid_argument("nested mesh embedding requires two nonempty meshes");
    }
    validate_boundary_tags(parent_mesh);
    validate_boundary_tags(child_mesh);

    const double parent_area = total_area(parent_mesh);
    const double child_area = total_area(child_mesh);
    const double area_scale = std::max({1.0, parent_area, child_area});
    if (std::abs(parent_area - child_area) > 1e-10 * area_scale)
        throw std::invalid_argument("nested meshes do not cover the same domain");
    if (!parent_mesh.boundary_edges.empty()) {
        if (child_mesh.boundary_edges.empty())
            throw std::invalid_argument("nested child mesh lost explicit boundary tags");
        for (BoundaryTag tag : {
                 BoundaryTag::Dirichlet,
                 BoundaryTag::Neumann,
                 BoundaryTag::Robin}) {
            const double parent_measure = boundary_measure(parent_mesh, tag);
            const double child_measure = boundary_measure(child_mesh, tag);
            const double scale = std::max({1.0, parent_measure, child_measure});
            if (std::abs(parent_measure - child_measure) > 1e-10 * scale)
                throw std::invalid_argument("nested meshes have inconsistent boundary tags");
        }
    }

    const int parent_nodes = static_cast<int>(parent_mesh.nodes.size());
    const int parent_elements = static_cast<int>(parent_mesh.elems.size());
    const int child_nodes = static_cast<int>(child_mesh.nodes.size());
    const int child_elements = static_cast<int>(child_mesh.elems.size());
    constexpr double tolerance = 2e-11;

    std::vector<int> element_parents(child_elements, -1);
    std::vector<std::array<std::array<double, 3>, 3>> element_weights(child_elements);
    for (int child = 0; child < child_elements; ++child) {
        const Triangle &child_triangle = child_mesh.elems[child];
        const Point2 centroid = (
            child_mesh.nodes[child_triangle[0]]
            + child_mesh.nodes[child_triangle[1]]
            + child_mesh.nodes[child_triangle[2]]) / 3.0;
        for (int parent = 0; parent < parent_elements; ++parent) {
            const Triangle &parent_triangle = parent_mesh.elems[parent];
            const auto centroid_weights = barycentric_coordinates(
                centroid,
                parent_mesh.nodes[parent_triangle[0]],
                parent_mesh.nodes[parent_triangle[1]],
                parent_mesh.nodes[parent_triangle[2]]);
            if (!inside_triangle(centroid_weights, tolerance)) continue;

            bool all_inside = true;
            std::array<std::array<double, 3>, 3> weights{};
            for (int local = 0; local < 3; ++local) {
                weights[local] = barycentric_coordinates(
                    child_mesh.nodes[child_triangle[local]],
                    parent_mesh.nodes[parent_triangle[0]],
                    parent_mesh.nodes[parent_triangle[1]],
                    parent_mesh.nodes[parent_triangle[2]]);
                if (!inside_triangle(weights[local], tolerance)) {
                    all_inside = false;
                    break;
                }
                weights[local] = clean_weights(weights[local]);
            }
            if (!all_inside) continue;
            if (element_parents[child] >= 0)
                throw std::invalid_argument("nested child triangle has more than one parent");
            element_parents[child] = parent;
            element_weights[child] = weights;
        }
        if (element_parents[child] < 0)
            throw std::invalid_argument("child mesh is not a conforming refinement of parent mesh");
    }

    std::vector<Eigen::Triplet<double>> element_triplets;
    std::vector<Eigen::Triplet<double>> dg_triplets;
    element_triplets.reserve(child_elements);
    dg_triplets.reserve(9 * child_elements);
    std::vector<std::map<int, double>> nodal_rows(child_nodes);
    for (int child = 0; child < child_elements; ++child) {
        const int parent = element_parents[child];
        const Triangle &parent_triangle = parent_mesh.elems[parent];
        const Triangle &child_triangle = child_mesh.elems[child];
        element_triplets.emplace_back(child, parent, 1.0);
        for (int child_local = 0; child_local < 3; ++child_local) {
            std::map<int, double> row;
            for (int parent_local = 0; parent_local < 3; ++parent_local) {
                const double weight = element_weights[child][child_local][parent_local];
                if (std::abs(weight) <= 1e-14) continue;
                row[parent_triangle[parent_local]] += weight;
                dg_triplets.emplace_back(
                    3 * child + child_local,
                    3 * parent + parent_local,
                    weight);
            }
            const int child_node = child_triangle[child_local];
            if (nodal_rows[child_node].empty()) {
                nodal_rows[child_node] = row;
            } else {
                std::set<int> columns;
                for (const auto &[column, value] : nodal_rows[child_node]) {
                    (void)value;
                    columns.insert(column);
                }
                for (const auto &[column, value] : row) {
                    (void)value;
                    columns.insert(column);
                }
                for (int column : columns) {
                    const double old_value = nodal_rows[child_node].contains(column)
                        ? nodal_rows[child_node].at(column) : 0.0;
                    const double new_value = row.contains(column) ? row.at(column) : 0.0;
                    if (std::abs(old_value - new_value) > 2e-10)
                        throw std::invalid_argument(
                            "nested nodal interpolation is inconsistent across a parent edge");
                }
            }
        }
    }

    std::vector<Eigen::Triplet<double>> node_triplets;
    node_triplets.reserve(3 * child_nodes);
    for (int child_node = 0; child_node < child_nodes; ++child_node) {
        if (nodal_rows[child_node].empty())
            throw std::invalid_argument("nested child mesh contains an unused node");
        Point2 reconstructed = Point2::Zero();
        for (const auto &[parent_node, weight] : nodal_rows[child_node]) {
            node_triplets.emplace_back(child_node, parent_node, weight);
            reconstructed += weight * parent_mesh.nodes[parent_node];
        }
        if ((reconstructed - child_mesh.nodes[child_node]).norm() > 2e-10)
            throw std::runtime_error("nested nodal prolongation does not reproduce coordinates");
    }

    Eigen::SparseMatrix<double> node_prolongation(child_nodes, parent_nodes);
    node_prolongation.setFromTriplets(node_triplets.begin(), node_triplets.end());
    Eigen::SparseMatrix<double> element_prolongation(child_elements, parent_elements);
    element_prolongation.setFromTriplets(element_triplets.begin(), element_triplets.end());
    Eigen::SparseMatrix<double> dg_prolongation(3 * child_elements, 3 * parent_elements);
    dg_prolongation.setFromTriplets(dg_triplets.begin(), dg_triplets.end());
    return {
        child_mesh,
        std::move(node_prolongation),
        std::move(element_prolongation),
        std::move(dg_prolongation)};
}

RefineOutput update_nested_mesh_embedding_after_parent_refinement(
    const TriMesh &old_parent_mesh,
    const RefineOutput &parent_refinement,
    const TriMesh &child_mesh,
    const std::vector<int> &old_child_parent_elements) {
    const TriMesh &new_parent_mesh = parent_refinement.mesh;
    if (old_parent_mesh.elems.empty() || new_parent_mesh.elems.empty()
        || child_mesh.elems.empty()) {
        throw std::invalid_argument(
            "incremental nested embedding requires nonempty meshes");
    }
    if (old_child_parent_elements.size() != child_mesh.elems.size()) {
        throw std::invalid_argument(
            "incremental nested embedding parent map has the wrong size");
    }
    const std::vector<int> new_parent_old_parents = fine_element_parents(
        parent_refinement.P_elem,
        static_cast<int>(new_parent_mesh.elems.size()),
        static_cast<int>(old_parent_mesh.elems.size()));
    std::vector<std::vector<int>> descendants(old_parent_mesh.elems.size());
    for (int element = 0;
         element < static_cast<int>(new_parent_old_parents.size()); ++element) {
        descendants[new_parent_old_parents[element]].push_back(element);
    }

    constexpr double tolerance = 2e-11;
    const int parent_nodes = static_cast<int>(new_parent_mesh.nodes.size());
    const int parent_elements = static_cast<int>(new_parent_mesh.elems.size());
    const int child_nodes = static_cast<int>(child_mesh.nodes.size());
    const int child_elements = static_cast<int>(child_mesh.elems.size());
    std::vector<int> element_parents(child_elements, -1);
    std::vector<std::array<std::array<double, 3>, 3>> element_weights(
        child_elements);
    for (int child = 0; child < child_elements; ++child) {
        const int old_parent = old_child_parent_elements[child];
        if (old_parent < 0
            || old_parent >= static_cast<int>(descendants.size())) {
            throw std::out_of_range(
                "incremental nested embedding old parent is out of range");
        }
        const Triangle &child_triangle = child_mesh.elems[child];
        const Point2 centroid = (
            child_mesh.nodes[child_triangle[0]]
            + child_mesh.nodes[child_triangle[1]]
            + child_mesh.nodes[child_triangle[2]]) / 3.0;
        for (int parent : descendants[old_parent]) {
            const Triangle &parent_triangle = new_parent_mesh.elems[parent];
            if (!inside_triangle(
                    barycentric_coordinates(
                        centroid,
                        new_parent_mesh.nodes[parent_triangle[0]],
                        new_parent_mesh.nodes[parent_triangle[1]],
                        new_parent_mesh.nodes[parent_triangle[2]]),
                    tolerance)) {
                continue;
            }
            std::array<std::array<double, 3>, 3> weights{};
            bool all_inside = true;
            for (int local = 0; local < 3; ++local) {
                weights[local] = barycentric_coordinates(
                    child_mesh.nodes[child_triangle[local]],
                    new_parent_mesh.nodes[parent_triangle[0]],
                    new_parent_mesh.nodes[parent_triangle[1]],
                    new_parent_mesh.nodes[parent_triangle[2]]);
                if (!inside_triangle(weights[local], tolerance)) {
                    all_inside = false;
                    break;
                }
                weights[local] = clean_weights(weights[local]);
            }
            if (!all_inside) continue;
            if (element_parents[child] >= 0) {
                throw std::invalid_argument(
                    "incremental nested child triangle has more than one parent");
            }
            element_parents[child] = parent;
            element_weights[child] = weights;
        }
        if (element_parents[child] < 0) {
            throw std::invalid_argument(
                "refined parent mesh is not contained in the fixed child mesh");
        }
    }

    std::vector<Eigen::Triplet<double>> element_triplets;
    std::vector<Eigen::Triplet<double>> dg_triplets;
    element_triplets.reserve(child_elements);
    dg_triplets.reserve(9 * child_elements);
    std::vector<std::map<int, double>> nodal_rows(child_nodes);
    for (int child = 0; child < child_elements; ++child) {
        const int parent = element_parents[child];
        const Triangle &parent_triangle = new_parent_mesh.elems[parent];
        const Triangle &child_triangle = child_mesh.elems[child];
        element_triplets.emplace_back(child, parent, 1.0);
        for (int child_local = 0; child_local < 3; ++child_local) {
            std::map<int, double> row;
            for (int parent_local = 0; parent_local < 3; ++parent_local) {
                const double weight =
                    element_weights[child][child_local][parent_local];
                if (std::abs(weight) <= 1e-14) continue;
                row[parent_triangle[parent_local]] += weight;
                dg_triplets.emplace_back(
                    3 * child + child_local,
                    3 * parent + parent_local,
                    weight);
            }
            const int child_node = child_triangle[child_local];
            if (nodal_rows[child_node].empty()) {
                nodal_rows[child_node] = row;
            } else {
                std::set<int> columns;
                for (const auto &[column, value] : nodal_rows[child_node]) {
                    (void)value;
                    columns.insert(column);
                }
                for (const auto &[column, value] : row) {
                    (void)value;
                    columns.insert(column);
                }
                for (int column : columns) {
                    const double old_value = nodal_rows[child_node].contains(column)
                        ? nodal_rows[child_node].at(column) : 0.0;
                    const double new_value = row.contains(column)
                        ? row.at(column) : 0.0;
                    if (std::abs(old_value - new_value) > 2e-10) {
                        throw std::invalid_argument(
                            "incremental nodal interpolation is inconsistent across an edge");
                    }
                }
            }
        }
    }

    std::vector<Eigen::Triplet<double>> node_triplets;
    node_triplets.reserve(3 * child_nodes);
    for (int child_node = 0; child_node < child_nodes; ++child_node) {
        if (nodal_rows[child_node].empty()) {
            throw std::invalid_argument(
                "incremental nested child mesh contains an unused node");
        }
        Point2 reconstructed = Point2::Zero();
        for (const auto &[parent_node, weight] : nodal_rows[child_node]) {
            node_triplets.emplace_back(child_node, parent_node, weight);
            reconstructed += weight * new_parent_mesh.nodes[parent_node];
        }
        if ((reconstructed - child_mesh.nodes[child_node]).norm() > 2e-10) {
            throw std::runtime_error(
                "incremental nodal prolongation does not reproduce coordinates");
        }
    }

    Eigen::SparseMatrix<double> node_prolongation(child_nodes, parent_nodes);
    node_prolongation.setFromTriplets(node_triplets.begin(), node_triplets.end());
    Eigen::SparseMatrix<double> element_prolongation(
        child_elements, parent_elements);
    element_prolongation.setFromTriplets(
        element_triplets.begin(), element_triplets.end());
    Eigen::SparseMatrix<double> dg_prolongation(
        3 * child_elements, 3 * parent_elements);
    dg_prolongation.setFromTriplets(dg_triplets.begin(), dg_triplets.end());
    return {
        child_mesh,
        std::move(node_prolongation),
        std::move(element_prolongation),
        std::move(dg_prolongation)};
}

std::vector<int> fine_element_parents(
    const Eigen::SparseMatrix<double> &prolongation,
    int fine_element_count,
    int coarse_element_count) {
    if (prolongation.rows() != fine_element_count
        || prolongation.cols() != coarse_element_count) {
        throw std::invalid_argument("element prolongation dimensions do not match the meshes");
    }
    std::vector<int> parents(fine_element_count, -1);
    for (int coarse = 0; coarse < prolongation.outerSize(); ++coarse) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(prolongation, coarse); it; ++it) {
            if (std::abs(it.value()) <= 1e-14) continue;
            if (std::abs(it.value() - 1.0) > 1e-12 || parents[it.row()] >= 0)
                throw std::runtime_error("element prolongation is not a unique parent map");
            parents[it.row()] = coarse;
        }
    }
    if (std::find(parents.begin(), parents.end(), -1) != parents.end())
        throw std::runtime_error("element prolongation leaves a fine element without a parent");
    return parents;
}


} // namespace alod
