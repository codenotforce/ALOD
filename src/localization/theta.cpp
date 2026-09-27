#include "alod/localization.hpp"
#include "alod/timing.hpp"
#include "alod/kernel_defect.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>
#include <iomanip>
#include <sstream>
#include <chrono>
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
    std::size_t applications=0,columns=0;
    std::unique_ptr<KernelDefectOperator> defect_operator;
    if((config.fused_defect||config.parallel_defect)&&dimension>config.dense_cross_check_max_dimension){
        const auto bytes=std::size_t(defect_rhs.nonZeros())*(sizeof(Complex)+sizeof(int))
            +std::size_t(defect_rhs.rows()+1)*sizeof(int);
        if(bytes<=dense_entry_limit*sizeof(Complex))
            defect_operator=std::make_unique<KernelDefectOperator>(riesz,defect_rhs,config.fused_defect);
        else PhaseTimer::counter("theta_parallel_memory_fallback",1,-1);
    }
    double operator_seconds=0;
    const auto measured=[&](LocalizationSpectrum result){
        result.operator_applications=applications;result.operator_columns=columns;
        result.operator_seconds=operator_seconds;return result;
    };
    const auto apply_block=[&](const ComplexMatrix& vectors,bool original=false)->ComplexMatrix {
        const auto start=std::chrono::steady_clock::now();
        ComplexMatrix result;
        if(defect_operator&&config.fused_defect&&!original)result=defect_operator->apply(vectors);
        else if(defect_operator&&config.parallel_defect)result=defect_operator->apply_global(vectors);
        else {
        ComplexMatrix rhs=defect_rhs*vectors;
        ComplexMatrix values=config.lightweight_riesz?riesz.apply_action(rhs):riesz.apply(rhs).values;
        result=defect_rhs.adjoint()*values;
        }
        ++applications;columns+=vectors.cols();
        operator_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        return result;
    };
    const auto dense_spectrum=[&](bool fallback) {
        if(static_cast<std::size_t>(dimension)*dimension>dense_entry_limit)
            throw std::runtime_error("localization dense Gram resource limit exceeded");
        ComplexMatrix gram(dimension,dimension);
        for(int first=0;first<dimension;first+=16) {
            int count=std::min(16,dimension-first);
            ComplexMatrix units=ComplexMatrix::Zero(dimension,count);
            for(int j=0;j<count;++j)units(first+j,j)=1;
            gram.middleCols(first,count)=apply_block(units,true);
        }
        gram=(0.5*(gram+gram.adjoint())).eval();
        LocalizationSpectrum result=largest_generalized_eigenvalue_dense(gram,ComplexMatrix(denominator),config);
        result.used_dense_fallback=result.used_dense_fallback||fallback;
        return measured(result);
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
                                        ComplexVector candidate,
                                        std::vector<ComplexVector>* images=nullptr,
                                        ComplexVector candidate_image=ComplexVector()) {
        {
            const double norm_squared=std::real(candidate.dot(multiply_real_sparse(energy,candidate)));
            if(!(norm_squared>0)||!std::isfinite(norm_squared))return false;
            candidate/=std::sqrt(norm_squared);
            if(images)candidate_image/=std::sqrt(norm_squared);
        }
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t j=0;j<basis.size();++j) {
                const Complex projection=basis[j].dot(multiply_real_sparse(energy,candidate));
                candidate-=basis[j]*projection;
                if(images)candidate_image-=(*images)[j]*projection;
            }
        }
        const double squared = std::real(candidate.dot(
            multiply_real_sparse(energy, candidate)));
        if (!(squared > 256.0 * std::numeric_limits<double>::epsilon()))
            return false;
        candidate /= std::sqrt(squared);
        if(images){candidate_image/=std::sqrt(squared);images->push_back(std::move(candidate_image));}
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
             && static_cast<int>(initial.size()) < (config.eigenvalue_relative_residual?std::max(1,block_size-1):block_size); ++column) {
            (void)append_orthonormal(
                initial, warm_block.col(config.eigenvalue_relative_residual?warm_block.cols()-1-column:column));
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
    ComplexMatrix search_applied(dimension, 0);
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
            / std::max(config.eigenvalue_relative_residual?1e-30:1.0, std::abs(result.lambda_max));
        result.iterations = iteration;
        PhaseTimer::counter("theta_ritz_iteration",iteration,-1);
        if((config.reuse_operator_actions||defect_operator) && result.relative_residual<=config.relative_tolerance){
            // Cached linear combinations may accumulate roundoff. Accept only
            // after a fresh action certifies the original generalized residual.
            applied=apply_block(iterate,true);
            const ComplexVector vector=iterate.col(block_size-1);
            result.lambda_max=std::max(0.,std::real(vector.dot(applied.col(block_size-1))));
            ComplexVector fresh_residual=applied.col(block_size-1)-result.lambda_max*energy_times.col(block_size-1);
            result.relative_residual=std::sqrt(std::max(0.,std::real(fresh_residual.dot(
                solve_real_sparse(factorization,fresh_residual)))))
                /std::max(config.eigenvalue_relative_residual?1e-30:1.,std::abs(result.lambda_max));
            if(result.relative_residual>config.relative_tolerance){
                search_direction.resize(dimension,0);search_applied.resize(dimension,0);
                continue;
            }
        }
        if (result.relative_residual <= config.relative_tolerance) {
            result.converged = true;
            result.dominant_vector = iterate.col(block_size - 1);
            result.dominant_subspace = iterate;
            break;
        }

        std::vector<ComplexVector> basis;
        std::vector<ComplexVector> images;
        basis.reserve(2 * block_size);
        // Only new preconditioned residual directions require patch solves.
        // Carry AX/AP through exactly the same orthogonalization as X/P.
        ComplexMatrix inverse_applied;
        if(config.reuse_operator_actions)inverse_applied=apply_block(inverse_residual);
        const auto append=[&](const ComplexMatrix& vectors,const ComplexMatrix& actions){
            for(int column=0;column<vectors.cols();++column){
                if(config.reuse_operator_actions)append_orthonormal(basis,vectors.col(column),&images,actions.col(column));
                else append_orthonormal(basis,vectors.col(column));
            }
        };
        append(iterate,applied);append(inverse_residual,inverse_applied);append(search_direction,search_applied);
        for (int seed = iteration * block_size;
             static_cast<int>(basis.size()) == block_size
             && seed < (iteration + 4) * block_size; ++seed) {
            ComplexVector candidate=deterministic_vector(seed);
            if(config.reuse_operator_actions){
                ComplexMatrix candidate_block=candidate;
                append_orthonormal(basis,candidate,&images,apply_block(candidate_block).col(0));
            }else append_orthonormal(basis,candidate);
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
        const ComplexMatrix subspace_action = config.reuse_operator_actions?as_matrix(images):apply_block(subspace);
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
        if(config.reuse_operator_actions)search_applied=subspace_action*direction_coefficients;
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
               << ", tolerance=" << config.relative_tolerance
               << ", lambda=" << result.lambda_max;
        throw std::runtime_error(detail.str());
    }
    result.dominant_vector = iterate.col(block_size - 1);
    result.dominant_subspace = iterate;
    return measured(result);
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
    Sparse denominator;
    {PhaseTimer timer("theta_denominator",-1);
    denominator=basis.transpose()*space.energy()*basis;
    denominator=0.5*(denominator+Sparse(denominator.transpose()));denominator.makeCompressed();
    }
    ComplexSparseMatrix defect;
    {PhaseTimer timer("theta_defect",-1);
    defect=adjoint?ComplexSparseMatrix(space.operators().system.adjoint()*space.test())
        :ComplexSparseMatrix(space.operators().system*space.trial());
    }
    if(config.relative_tolerance<1e-6){
        PhaseTimer timer("theta_cancellation_projection",-1);
        // R annihilates range(I_H^*). Remove that component before applying R
        // and before the dual dot product to avoid cancellation at small Theta.
        // I_H P_H = I on free coarse nodes, so this leaves D^* R D unchanged.
        ComplexSparseMatrix coarse_action=space.prolongation().transpose().cast<Complex>()*defect;
        defect-=space.interpolation().transpose().cast<Complex>()*coarse_action;
        defect.makeCompressed();
    }
    LocalizationResult result;result.identity=identity;
    {PhaseTimer timer("theta_ritz",-1);
    result.spectrum=reference_defect_spectrum_matrix_free(riesz,defect,denominator,config,space.limits().maximum_dense_entries);}
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
