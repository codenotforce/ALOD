#include "alod/problems.hpp"
#include "helmholtz/boundary.h"
#include "helmholtz/operators.h"
#include "mesh/refine.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--help") {
            std::cout << "Usage: alod_fem_smoke E1|E2 [NVB-level 0..8]\n"
                      << "Outputs a nominal fixed-grid FEM regression as JSON.\n";
            return 0;
        }
        if (argc < 2 || argc > 3) throw std::invalid_argument("expected E1|E2 and optional level");
        const std::string id = argv[1];
        int level = 6;
        if (argc == 3) {
            std::size_t used = 0;
            level = std::stoi(argv[2], &used);
            if (used != std::string(argv[2]).size()) throw std::invalid_argument("invalid level");
        }
        if (level < 0 || level > 8) throw std::invalid_argument("smoke level must be in 0..8");
        const auto problem = alod::make_problem(id);
        const auto quadrature = alod::paper_quadrature(id);
        using namespace lod2d;
        using namespace lod2d::helmholtz;
        const auto refined = refine_mesh_nvb(problem.initial_mesh, level);
        const auto& mesh = refined.mesh;
        validate_boundary_tags(mesh);
        Eigen::VectorXd coarse_linear(problem.initial_mesh.nodes.size());
        for (int i = 0; i < coarse_linear.size(); ++i) {
            const auto& x = problem.initial_mesh.nodes[i];
            coarse_linear[i] = 1 + 2*x.x() - 3*x.y();
        }
        const Eigen::VectorXd injected = refined.P_node * coarse_linear;
        double injection_error = 0;
        for (int i = 0; i < injected.size(); ++i) {
            const auto& x = mesh.nodes[i];
            injection_error = std::max(injection_error, std::abs(injected[i] - (1+2*x.x()-3*x.y())));
        }
        const auto operators = assemble_helmholtz_operators(mesh, problem.wavenumber);
        const auto load = assemble_helmholtz_load(mesh, problem.source, quadrature, problem.quadrature_context);
        const auto solution = solve_helmholtz_fem(operators, load, HelmholtzFemSolverKind::Umfpack);
        ComplexVector residual = operators.system * solution - load;
        ComplexVector free_load = load;
        for (const int i : operators.dirichlet_nodes) { residual[i] = 0; free_load[i] = 0; }
        const double relative_residual = residual.norm() / std::max(1e-12, free_load.norm());
        const auto error = compute_helmholtz_error(mesh, solution, problem.wavenumber,
            problem.exact, problem.exact_gradient, quadrature, problem.quadrature_context);
        const auto norm = compute_helmholtz_error(mesh, ComplexVector::Zero(mesh.nodes.size()),
            problem.wavenumber, problem.exact, problem.exact_gradient, quadrature, problem.quadrature_context);
        const auto areas = compute_area(mesh);
        double area = 0;
        for (double a : areas) { if (!(a > 0)) throw std::runtime_error("nonpositive element area"); area += a; }
        if (!solution.allFinite() || !std::isfinite(error.energy) || !std::isfinite(relative_residual)
            || relative_residual > 1e-9 || injection_error > 1e-13)
            throw std::runtime_error("FEM or nested injection invariant failed");
        std::cout << std::setprecision(17)
            << "{\"problem\":\"" << id << "\",\"level\":" << level
            << ",\"backend\":\"UMFPACK\",\"nodes\":" << mesh.nodes.size()
            << ",\"elements\":" << mesh.elems.size()
            << ",\"free_dof\":" << mesh.nodes.size()-operators.dirichlet_nodes.size()
            << ",\"area\":" << area
            << ",\"boundary_dirichlet\":" << boundary_measure(mesh, BoundaryTag::Dirichlet)
            << ",\"boundary_neumann\":" << boundary_measure(mesh, BoundaryTag::Neumann)
            << ",\"boundary_robin\":" << boundary_measure(mesh, BoundaryTag::Robin)
            << ",\"stiffness_norm\":" << operators.stiffness.norm()
            << ",\"mass_norm\":" << operators.mass.norm()
            << ",\"load_norm\":" << load.norm()
            << ",\"energy_error\":" << error.energy
            << ",\"l2_error\":" << error.l2
            << ",\"exact_energy_norm\":" << norm.energy
            << ",\"relative_energy_error\":" << error.energy/norm.energy
            << ",\"relative_residual\":" << relative_residual
            << ",\"injection_error\":" << injection_error << ",\"solution\":[";
        for (int i = 0; i < solution.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << '[' << solution[i].real() << ',' << solution[i].imag() << ']';
        }
        std::cout << "]}\n";
    } catch (const std::exception& e) {
        std::cerr << "alod_fem_smoke: " << e.what() << '\n';
        return 1;
    }
}
