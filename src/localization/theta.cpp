#include "alod/localization.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>
#include <iomanip>
#include <sstream>
namespace alod {
LocalizationSpectrum largest_generalized_eigenvalue_dense(const ComplexMatrix&,const ComplexMatrix&,const LocalizationEigenConfig&);
namespace {
ComplexVector multiply_real_sparse(
    const Eigen::SparseMatrix<double> &matrix,
    const ComplexVector &vector) {
    ComplexVector result(vector.size());
    result.real() = matrix * vector.real();
    result.imag() = matrix * vector.imag();
    return result;
}

ComplexMatrix multiply_real_sparse(
    const Eigen::SparseMatrix<double> &matrix,
    const ComplexMatrix &vectors) {
    ComplexMatrix result(vectors.rows(), vectors.cols());
    result.real() = matrix * vectors.real();
    result.imag() = matrix * vectors.imag();
    return result;
}

ComplexVector solve_real_sparse(
    const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &factorization,
    const ComplexVector &rhs) {
    ComplexVector result(rhs.size());
    result.real() = factorization.solve(rhs.real());
    result.imag() = factorization.solve(rhs.imag());
    return result;
}

ComplexMatrix solve_real_sparse(
    const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &factorization,
    const ComplexMatrix &rhs) {
    ComplexMatrix result(rhs.rows(), rhs.cols());
    result.real() = factorization.solve(rhs.real());
    result.imag() = factorization.solve(rhs.imag());
    return result;
}


LocalizationSpectrum reference_defect_spectrum_matrix_free(
    AdditiveKernelRieszContext& riesz,const ComplexSparseMatrix& defect_rhs,
    const Sparse& denominator,const LocalizationEigenConfig& config,std::size_t dense_entry_limit) {
    const int dimension=denominator.rows();
    const ComplexMatrix warm_block=config.warm_start?config.warm_start->block:ComplexMatrix();
    const ComplexVector warm_vector;
    const auto apply_block=[&](const ComplexMatrix& vectors)->ComplexMatrix {
        return defect_rhs.adjoint()*riesz.apply(defect_rhs*vectors).values;
    };
    const auto dense_spectrum=[&](bool fallback) {
        if(static_cast<std::size_t>(dimension)*dimension>dense_entry_limit)
            throw std::runtime_error("localization dense Gram resource limit exceeded");
        ComplexMatrix gram(dimension,dimension);
        for(int first=0;first<dimension;first+=16) {
            int count=std::min(16,dimension-first);
            ComplexMatrix units=ComplexMatrix::Zero(dimension,count);
            for(int j=0;j<count;++j)units(first+j,j)=1;
            gram.middleCols(first,count)=apply_block(units);
        }
        gram=(0.5*(gram+gram.adjoint())).eval();
        LocalizationSpectrum result=largest_generalized_eigenvalue_dense(gram,ComplexMatrix(denominator),config);
        result.used_dense_fallback=result.used_dense_fallback||fallback;
        return result;
    };
    if(dimension<=config.dense_cross_check_max_dimension)return dense_spectrum(false);
    Eigen::SparseMatrix<double> energy = 0.5 * (
        denominator
        + Eigen::SparseMatrix<double>(denominator.transpose()));
    energy.makeCompressed();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factorization;
    factorization.compute(energy);
    if (factorization.info() != Eigen::Success) {
        throw std::runtime_error(
            "reference corrector coarse energy factorization failed");
    }

    LocalizationSpectrum result;
    result.used_sparse_generalized_solver = true;
    const bool valid_warm_block = warm_block.rows() == dimension
        && warm_block.cols() > 0
        && warm_block.allFinite()
        && warm_block.norm() > 0.0;
    const bool valid_warm_vector = warm_vector.size() == dimension
        && warm_vector.allFinite()
        && warm_vector.norm() > 0.0;
    result.used_warm_start = valid_warm_block || valid_warm_vector;
    const auto append_orthonormal = [&](std::vector<ComplexVector> &basis,
                                        ComplexVector candidate) {
        for (int pass = 0; pass < 2; ++pass) {
            for (const ComplexVector &vector : basis) {
                candidate -= vector * vector.dot(
                    multiply_real_sparse(energy, candidate));
            }
        }
        const double squared = std::real(candidate.dot(
            multiply_real_sparse(energy, candidate)));
        if (!(squared > 256.0 * std::numeric_limits<double>::epsilon()))
            return false;
        candidate /= std::sqrt(squared);
        basis.push_back(std::move(candidate));
        return true;
    };
    const auto deterministic_vector = [&](int seed) {
        ComplexVector vector(dimension);
        for (int index = 0; index < dimension; ++index) {
            const double phase = static_cast<double>(
                (index + 1) * (seed + 1));
            vector(index) = Complex(
                std::sin(phase), 0.5 * std::cos(phase * 1.6180339887498948));
        }
        return vector;
    };
    const auto as_matrix = [&](const std::vector<ComplexVector> &columns) {
        ComplexMatrix matrix(dimension, columns.size());
        for (int column = 0; column < static_cast<int>(columns.size()); ++column)
            matrix.col(column) = columns[column];
        return matrix;
    };

    // Four vectors cover the clustered dominant modes seen in E1 without the
    // memory and RHS traffic of the earlier eight-vector prototype.
    const int block_size = std::min(4, dimension);
    std::vector<ComplexVector> initial;
    initial.reserve(block_size);
    if (valid_warm_block) {
        for (int column = 0;
             column < warm_block.cols()
             && static_cast<int>(initial.size()) < block_size; ++column) {
            (void)append_orthonormal(
                initial, warm_block.col(column));
        }
    }
    if (valid_warm_vector && static_cast<int>(initial.size()) < block_size)
        (void)append_orthonormal(initial, warm_vector);
    for (int seed = 0;
         static_cast<int>(initial.size()) < block_size
         && seed < 4 * dimension; ++seed) {
        (void)append_orthonormal(initial, deterministic_vector(seed));
    }
    if (static_cast<int>(initial.size()) != block_size) {
        throw std::runtime_error(
            "reference corrector block eigensolve could not initialize an energy-orthonormal basis");
    }
    ComplexMatrix iterate = as_matrix(initial);
    ComplexMatrix applied = apply_block(iterate);
    ComplexMatrix search_direction(dimension, 0);
    for (int iteration = 1; iteration <= config.maximum_iterations;
         ++iteration) {
        ComplexMatrix energy_times = multiply_real_sparse(energy, iterate);
        ComplexMatrix projected = iterate.adjoint() * applied;
        projected = (0.5 * (projected + projected.adjoint())).eval();
        Eigen::SelfAdjointEigenSolver<ComplexMatrix> eigensolver(projected);
        if (eigensolver.info() != Eigen::Success) {
            throw std::runtime_error(
                "reference corrector matrix-free block Ritz solve failed");
        }
        iterate = iterate * eigensolver.eigenvectors();
        applied = applied * eigensolver.eigenvectors();
        energy_times = energy_times * eigensolver.eigenvectors();
        const Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
        result.lambda_max = std::max(0.0, eigenvalues(block_size - 1));
        ComplexMatrix residual = applied;
        for (int column = 0; column < block_size; ++column)
            residual.col(column) -=
                eigenvalues(column) * energy_times.col(column);
        const ComplexMatrix inverse_residual = solve_real_sparse(
            factorization, residual);
        const ComplexVector top_residual = residual.col(block_size - 1);
        const ComplexVector top_inverse_residual =
            inverse_residual.col(block_size - 1);
        const double dual_squared = std::max(
            0.0, std::real(top_residual.dot(top_inverse_residual)));
        result.relative_residual = std::sqrt(dual_squared)
            / std::max(1.0, std::abs(result.lambda_max));
        result.iterations = iteration;
        if (result.relative_residual <= config.relative_tolerance) {
            result.converged = true;
            result.dominant_vector = iterate.col(block_size - 1);
            result.dominant_subspace = iterate;
            break;
        }

        std::vector<ComplexVector> basis;
        basis.reserve(2 * block_size);
        for (int column = 0; column < block_size; ++column)
            (void)append_orthonormal(basis, iterate.col(column));
        for (int column = 0; column < block_size; ++column)
            (void)append_orthonormal(basis, inverse_residual.col(column));
        for (int column = 0; column < search_direction.cols(); ++column)
            (void)append_orthonormal(basis, search_direction.col(column));
        for (int seed = iteration * block_size;
             static_cast<int>(basis.size()) == block_size
             && seed < (iteration + 4) * block_size; ++seed) {
            (void)append_orthonormal(basis, deterministic_vector(seed));
        }
        if (static_cast<int>(basis.size()) == block_size) {
            if (dimension <= config.dense_fallback_max_dimension) {
                LocalizationSpectrum fallback = dense_spectrum(true);
                return fallback;
            }
            throw std::runtime_error(
                "reference corrector block Ritz iteration stagnated without an independent residual direction");
        }
        const ComplexMatrix subspace = as_matrix(basis);
        const ComplexMatrix subspace_action = apply_block(subspace);
        ComplexMatrix small = subspace.adjoint() * subspace_action;
        small = (0.5 * (small + small.adjoint())).eval();
        Eigen::SelfAdjointEigenSolver<ComplexMatrix> subspace_solver(small);
        if (subspace_solver.info() != Eigen::Success) {
            throw std::runtime_error(
                "reference corrector matrix-free block subspace solve failed");
        }
        const ComplexMatrix selected =
            subspace_solver.eigenvectors().rightCols(block_size);
        ComplexMatrix direction_coefficients = selected;
        direction_coefficients.topRows(block_size).setZero();
        search_direction = subspace * direction_coefficients;
        iterate = subspace * selected;
        applied = subspace_action * selected;
    }
    if (!result.converged) {
        if (dimension <= config.dense_fallback_max_dimension) {
            LocalizationSpectrum fallback = dense_spectrum(true);
            return fallback;
        }
        std::ostringstream detail;
        detail << std::setprecision(17)
               << "reference corrector matrix-free eigensolve did not converge: "
               << "dimension=" << dimension
               << ", iterations=" << result.iterations
               << ", relative_residual=" << result.relative_residual
               << ", tolerance=" << config.relative_tolerance;
        throw std::runtime_error(detail.str());
    }
    result.dominant_vector = iterate.col(block_size - 1);
    result.dominant_subspace = iterate;
    return result;
}

} // namespace
LocalizationResult localization_theta(const LodSpace& space,AdditiveKernelRieszContext& riesz,
    const LocalizationEigenConfig& config,bool adjoint) {
    if(riesz.reference_identity()!=space.reference_identity())
        throw std::invalid_argument("localization Riesz context belongs to another reference space");
    const std::string identity=space.identity()+":"+riesz.patch_policy_name()+(adjoint?":adjoint":":primal");
    if(config.maximum_iterations<1 || !std::isfinite(config.relative_tolerance) || config.relative_tolerance<=0
        || config.dense_cross_check_max_dimension<0 || config.dense_fallback_max_dimension<0)
        throw std::invalid_argument("invalid localization eigen configuration");
    if(config.warm_start && (config.warm_start->identity!=identity
        || config.warm_start->block.rows()!=space.trial().cols() || config.warm_start->block.cols()<1
        || !config.warm_start->block.allFinite() || config.warm_start->block.norm()==0))
        throw std::invalid_argument("localization warm start identity/dimensions/values mismatch");
    const auto& basis=space.coarse_basis();
    Sparse denominator=basis.transpose()*space.energy()*basis;
    denominator=0.5*(denominator+Sparse(denominator.transpose()));denominator.makeCompressed();
    ComplexSparseMatrix defect=adjoint?ComplexSparseMatrix(space.operators().system.adjoint()*space.test())
        :ComplexSparseMatrix(space.operators().system*space.trial());
    LocalizationResult result;result.identity=identity;
    result.spectrum=reference_defect_spectrum_matrix_free(riesz,defect,denominator,config,space.limits().maximum_dense_entries);
    result.theta=std::sqrt(result.spectrum.lambda_max);
    result.warm_start={identity,result.spectrum.dominant_subspace};
    return result;
}
LocalizationRatio localization_ratio(double theta,double eta) {
    if(!std::isfinite(theta) || !std::isfinite(eta) || theta<0 || eta<0)
        throw std::invalid_argument("invalid localization ratio inputs");
    if(eta==0)return {theta==0?"zero_defect":"zero_denominator",std::nullopt};
    return {"finite",theta/eta};
}
}
