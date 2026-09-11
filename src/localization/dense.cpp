#include "alod/localization.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
namespace alod {
ComplexMatrix whitened_matrix(
    const ComplexMatrix &numerator,
    const ComplexMatrix &denominator,
    ComplexMatrix &inverse_lower) {
    if (numerator.rows() != numerator.cols()
        || denominator.rows() != denominator.cols()
        || numerator.rows() != denominator.rows()) {
        throw std::invalid_argument("generalized spectrum dimensions disagree");
    }
    const ComplexMatrix hermitian_denominator =
        0.5 * (denominator + denominator.adjoint());
    Eigen::LLT<ComplexMatrix> factorization(hermitian_denominator);
    if (factorization.info() != Eigen::Success)
        throw std::runtime_error("coarse energy matrix is not positive definite");
    const ComplexMatrix lower = factorization.matrixL();
    inverse_lower = lower.triangularView<Eigen::Lower>().solve(
        ComplexMatrix::Identity(lower.rows(), lower.cols()));
    const ComplexMatrix hermitian_numerator =
        0.5 * (numerator + numerator.adjoint());
    ComplexMatrix result = inverse_lower * hermitian_numerator
        * inverse_lower.adjoint();
    return 0.5 * (result + result.adjoint());
}

LocalizationSpectrum largest_generalized_eigenvalue_dense(
    const ComplexMatrix &numerator,
    const ComplexMatrix &denominator,
    const LocalizationEigenConfig &config) {
    if (config.maximum_iterations <= 0
        || !(config.relative_tolerance > 0.0)
        || config.dense_cross_check_max_dimension < 0
        || config.dense_fallback_max_dimension < 0) {
        throw std::invalid_argument("localization eigen configuration is invalid");
    }
    ComplexMatrix inverse_lower;
    const ComplexMatrix whitened = whitened_matrix(
        numerator, denominator, inverse_lower);
    const int dimension = whitened.rows();
    LocalizationSpectrum result;
    if (dimension == 0) {
        result.converged = true;
        return result;
    }

    const ComplexVector warm_vector=config.warm_start?ComplexVector(config.warm_start->block.col(config.warm_start->block.cols()-1)):ComplexVector();
    ComplexVector iterate;
    if (warm_vector.size() == dimension
        && warm_vector.allFinite()
        && warm_vector.norm() > 0.0) {
        result.used_warm_start = true;
        const ComplexMatrix lower_inverse_adjoint =
            inverse_lower.adjoint();
        // y=L^*x and x=L^{-H}y.
        iterate = lower_inverse_adjoint.fullPivLu().solve(
            warm_vector);
    } else {
        iterate = ComplexVector::Ones(dimension);
        for (int index = 0; index < dimension; ++index)
            iterate(index) += Complex(0.0, 0.125 * (index + 1));
    }
    iterate.normalize();
    const double matrix_scale = std::max(1.0, whitened.norm());
    for (int iteration = 1; iteration <= config.maximum_iterations; ++iteration) {
        const ComplexVector applied = whitened * iterate;
        const double norm = applied.norm();
        if (norm <= std::numeric_limits<double>::epsilon() * matrix_scale) {
            result.lambda_max = 0.0;
            result.relative_residual = 0.0;
            result.iterations = iteration;
            result.converged = true;
            break;
        }
        iterate = applied / norm;
        const ComplexVector next = whitened * iterate;
        result.lambda_max = std::max(
            0.0, std::real(iterate.dot(next)));
        result.relative_residual =
            (next - result.lambda_max * iterate).norm() / matrix_scale;
        result.iterations = iteration;
        if (result.relative_residual <= config.relative_tolerance) {
            result.converged = true;
            break;
        }
    }
    if (!result.converged) {
        if (dimension > config.dense_fallback_max_dimension) {
            throw std::runtime_error(
                "localization largest-eigenvalue iteration did not converge: "
                "dimension=" + std::to_string(dimension)
                + ", iterations=" + std::to_string(result.iterations)
                + ", relative_residual="
                + std::to_string(result.relative_residual));
        }
        Eigen::SelfAdjointEigenSolver<ComplexMatrix> dense_solver(whitened);
        if (dense_solver.info() != Eigen::Success)
            throw std::runtime_error("localization dense eigen fallback failed");
        iterate = dense_solver.eigenvectors().col(dimension - 1);
        result.lambda_max = std::max(
            0.0, dense_solver.eigenvalues()(dimension - 1));
        const ComplexVector applied = whitened * iterate;
        result.relative_residual =
            (applied - result.lambda_max * iterate).norm() / matrix_scale;
        if (!std::isfinite(result.relative_residual)
            || result.relative_residual > config.relative_tolerance) {
            throw std::runtime_error(
                "localization dense eigen fallback residual is too large");
        }
        result.converged = true;
        result.used_dense_fallback = true;
        result.dense_cross_checked = true;
        result.dense_lambda_max = result.lambda_max;
        result.dense_relative_difference = 0.0;
    }

    result.dominant_vector = inverse_lower.adjoint() * iterate;
    const double energy_norm = std::sqrt(std::max(
        0.0, std::real(result.dominant_vector.dot(
            denominator * result.dominant_vector))));
    if (energy_norm > 0.0) result.dominant_vector /= energy_norm;
    result.dominant_subspace.resize(dimension, 1);
    result.dominant_subspace.col(0) = result.dominant_vector;

    if (!result.used_dense_fallback
        && dimension <= config.dense_cross_check_max_dimension) {
        Eigen::SelfAdjointEigenSolver<ComplexMatrix> dense_solver(whitened);
        if (dense_solver.info() != Eigen::Success)
            throw std::runtime_error("dense localization eigen cross-check failed");
        result.dense_cross_checked = true;
        result.dense_lambda_max = std::max(
            0.0, dense_solver.eigenvalues().maxCoeff());
        result.dense_relative_difference = std::abs(
            result.lambda_max - result.dense_lambda_max)
            / std::max(1.0, result.dense_lambda_max);
    }
    return result;
}

}
