#include "helmholtz/patch_solver.h"




#include <Eigen/LU>
#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <memory>
#include <string>
#include <cstdint>
#include <vector>

namespace lod2d::helmholtz {
namespace {

using ComplexSparseLu = Eigen::SparseLU<ComplexSparseMatrix>;

struct SparseLuPatternCache {
    ComplexSparseLu solver;
    std::vector<int> outer_indices;
    std::vector<int> inner_indices;
    int rows = -1;
    int cols = -1;
    std::vector<Complex> values;
    std::uint64_t last_use = 0;

    bool matches(const ComplexSparseMatrix &matrix) const {
        if (rows != matrix.rows() || cols != matrix.cols()) return false;
        if (outer_indices.size() != static_cast<std::size_t>(matrix.outerSize() + 1)
            || inner_indices.size() != static_cast<std::size_t>(matrix.nonZeros()))
            return false;
        return std::equal(
                   outer_indices.begin(), outer_indices.end(), matrix.outerIndexPtr())
            && std::equal(
                   inner_indices.begin(), inner_indices.end(), matrix.innerIndexPtr());
    }

    bool values_match(const ComplexSparseMatrix &matrix) const {
        return matches(matrix)
            && values.size() == static_cast<std::size_t>(matrix.nonZeros())
            && std::equal(values.begin(), values.end(), matrix.valuePtr());
    }

    ComplexSparseLu &factorize(
        const ComplexSparseMatrix &matrix,
        bool reuse_identical,
        bool &symbolic_reused,
        bool &factorization_reused,
        const char *description) {
        symbolic_reused = matches(matrix);
        factorization_reused = reuse_identical && values_match(matrix);
        if (!symbolic_reused) {
            solver.analyzePattern(matrix);
            rows = matrix.rows();
            cols = matrix.cols();
            outer_indices.assign(
                matrix.outerIndexPtr(), matrix.outerIndexPtr() + matrix.outerSize() + 1);
            inner_indices.assign(
                matrix.innerIndexPtr(), matrix.innerIndexPtr() + matrix.nonZeros());
        }
        if (factorization_reused) return solver;
        solver.factorize(matrix);
        if (solver.info() != Eigen::Success)
            throw std::runtime_error(std::string(description) + " factorization failed");
        values.assign(matrix.valuePtr(), matrix.valuePtr() + matrix.nonZeros());
        return solver;
    }
};

struct SparseLuCachePool {
    std::vector<std::unique_ptr<SparseLuPatternCache>> entries;
    std::uint64_t clock = 0;

    ComplexSparseLu &factorize(
        const ComplexSparseMatrix &matrix,
        int requested_slots,
        bool reuse_identical,
        bool &symbolic_reused,
        bool &factorization_reused,
        const char *description) {
        const int slots = std::max(1, requested_slots);
        if (static_cast<int>(entries.size()) > slots)
            entries.resize(slots);

        SparseLuPatternCache *selected = nullptr;
        if (reuse_identical) {
            for (const auto &entry : entries) {
                if (entry->values_match(matrix)) {
                    selected = entry.get();
                    break;
                }
            }
        }
        if (!selected && !reuse_identical) {
            for (const auto &entry : entries) {
                if (entry->matches(matrix)) {
                    selected = entry.get();
                    break;
                }
            }
        }
        if (!selected && static_cast<int>(entries.size()) < slots) {
            entries.push_back(std::make_unique<SparseLuPatternCache>());
            selected = entries.back().get();
        }
        if (!selected) {
            selected = entries.front().get();
            for (const auto &entry : entries) {
                if (entry->last_use < selected->last_use)
                    selected = entry.get();
            }
        }
        selected->last_use = ++clock;
        return selected->factorize(
            matrix, reuse_identical, symbolic_reused,
            factorization_reused, description);
    }
};

thread_local SparseLuCachePool helmholtz_cache;

struct SchurInputs {
    ComplexMatrix z;
    ComplexMatrix y;
    HelmholtzPatchSolveDiagnostics diagnostics;
};

SchurInputs solve_direct_helmholtz_blocks(
    const HelmholtzPatchSystem &system,
    const HelmholtzPatchSolverConfig &config) {
    const int n = system.helmholtz.rows();
    const int m = static_cast<int>(system.constraints.rows());
    ComplexMatrix combined(n, m + system.rhs.cols());
    if (m > 0)
        combined.leftCols(m) = system.constraints.transpose().cast<Complex>();
    combined.rightCols(system.rhs.cols()) = system.rhs;

    bool symbolic_reused = false;
    bool factorization_reused = false;
    ComplexSparseLu &solver = helmholtz_cache.factorize(
        system.helmholtz, config.symbolic_cache_slots,
        config.reuse_identical_factorization,
        symbolic_reused, factorization_reused,
        "Helmholtz patch A block");
    const ComplexMatrix solved = solver.solve(combined);
    if (solver.info() != Eigen::Success || !solved.allFinite())
        throw std::runtime_error("Helmholtz patch A-block solve failed");

    SchurInputs result;
    result.z = solved.leftCols(m);
    result.y = solved.rightCols(system.rhs.cols());
    result.diagnostics.symbolic_reused = symbolic_reused;
    result.diagnostics.factorization_reused = factorization_reused;
    return result;
}

HelmholtzPatchSolveResult recover_schur_solution(
    const HelmholtzPatchSystem &system,
    SchurInputs inputs) {
    const int m = static_cast<int>(system.constraints.rows());
    HelmholtzPatchSolveResult result;
    result.diagnostics = inputs.diagnostics;
    if (m == 0) {
        result.corrector = std::move(inputs.y);
        result.multipliers.resize(0, system.rhs.cols());
        return result;
    }

    const ComplexMatrix constraints = system.constraints.cast<Complex>();
    const ComplexMatrix schur = constraints * inputs.z;
    const ComplexMatrix schur_rhs = constraints * inputs.y;
    Eigen::FullPivLU<ComplexMatrix> factorization(schur);
    result.diagnostics.schur_rcond = factorization.rcond();
    if (!factorization.isInvertible())
        throw std::runtime_error("Helmholtz patch Schur complement is singular");
    result.multipliers = factorization.solve(schur_rhs);
    if (!result.multipliers.allFinite())
        throw std::runtime_error("Helmholtz patch Schur solve returned non-finite values");
    result.corrector = inputs.y - inputs.z * result.multipliers;
    result.diagnostics.schur_residual =
        (schur * result.multipliers - schur_rhs).norm()
        / std::max(1.0, schur_rhs.norm());
    return result;
}

void compute_final_residuals(
    const HelmholtzPatchSystem &system,
    HelmholtzPatchSolveResult &result) {
    ComplexMatrix equation_residual =
        system.helmholtz * result.corrector - system.rhs;
    if (system.constraints.rows() > 0) {
        equation_residual.noalias() +=
            system.constraints.transpose().cast<Complex>() * result.multipliers;
    }
    const double rhs_scale = std::max(1.0, system.rhs.norm());
    result.diagnostics.primal_residual = equation_residual.norm() / rhs_scale;
    result.diagnostics.constraint_residual = system.constraints.rows() > 0
        ? (system.constraints.cast<Complex>() * result.corrector).norm() / rhs_scale
        : 0.0;

    ComplexMatrix adjoint_residual = system.helmholtz.adjoint()
        * result.corrector.conjugate() - system.rhs.conjugate();
    if (system.constraints.rows() > 0) {
        adjoint_residual.noalias() +=
            system.constraints.transpose().cast<Complex>()
            * result.multipliers.conjugate();
    }
    result.diagnostics.adjoint_residual = adjoint_residual.norm() / rhs_scale;
}

} // namespace

HelmholtzPatchSolveResult solve_helmholtz_patch(
    const HelmholtzPatchSystem &system,
    const HelmholtzPatchSolverConfig &config) {
    if (system.helmholtz.rows() <= 0
        || system.helmholtz.rows() != system.helmholtz.cols())
        throw std::invalid_argument("Helmholtz patch operator must be square and nonempty");
    if (system.rhs.rows() != system.helmholtz.rows())
        throw std::invalid_argument("Helmholtz patch right-hand side has the wrong size");
    if (system.constraints.cols() != system.helmholtz.cols())
        throw std::invalid_argument("Helmholtz patch constraints have the wrong size");
    if (config.symbolic_cache_slots <= 0)
        throw std::invalid_argument(
            "Helmholtz symbolic cache slots must be positive");

    auto result = recover_schur_solution(system, solve_direct_helmholtz_blocks(system, config));
    compute_final_residuals(system, result);
    return result;
}
}
