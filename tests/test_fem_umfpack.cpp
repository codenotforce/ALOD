#include "alod/problems.hpp"
#include "mesh/refine.h"
#include "helmholtz/operators.h"
#include <iostream>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
int main() {
    try {
        for (const auto* id : {"E1", "E2"}) for (double k : {8., 16., 128.}) {
            auto p = alod::make_problem(id, k);
            auto mesh = refine_mesh_nvb(p.initial_mesh, 6).mesh;
            auto ops = assemble_helmholtz_operators(mesh, k);
            auto b = assemble_helmholtz_load(mesh, p.source, alod::paper_quadrature(id), p.quadrature_context);
            auto reference = solve_helmholtz_fem(ops, b, HelmholtzFemSolverKind::SparseLu);
            auto actual = solve_helmholtz_fem(ops, b, HelmholtzFemSolverKind::Umfpack);
            require((actual-reference).norm() < 1e-9*std::max(1., reference.norm()), "UMFPACK single RHS differs from SparseLU");
            ComplexMatrix loads(b.size(), 2);
            loads.col(0) = b;
            loads.col(1) = Complex(.3, -.7)*b;
            auto batch = solve_helmholtz_fem_batch(ops, loads, HelmholtzFemSolverKind::Umfpack);
            require((batch.col(0)-actual).norm() < 1e-10*std::max(1., actual.norm()), "UMFPACK block first RHS differs");
            require((batch.col(1)-Complex(.3, -.7)*actual).norm() < 1e-10*std::max(1., actual.norm()), "UMFPACK complex block RHS differs");
            for (int i : ops.dirichlet_nodes) require(batch.row(i).norm()==0., "Dirichlet constraint lost");
            ComplexMatrix residual = ops.system*batch-loads;
            for (int i : ops.dirichlet_nodes) residual.row(i).setZero();
            require(residual.norm()/std::max(1e-30, loads.norm()) < 1e-10, "UMFPACK free-row residual failed");
        }
        std::cout << "64-bit UMFPACK single and block FEM solves passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
