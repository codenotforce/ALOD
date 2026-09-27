// Fine-grid inf-sup diagnostic. No LOD spaces, loads, or audit solves are built.
#include "alod/problems.hpp"
#include "helmholtz/operators.h"
#include "mesh/refine.h"
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <Eigen/SparseCholesky>
#include <Eigen/UmfPackSupport>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>

using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;
// Fine-grid calibration can exceed the workspace-index range of UMFPACK's
// 32-bit interface even when the input matrix itself has fewer than 2^31 entries.
using Sparse = Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor, SuiteSparse_long>;
using Clock = std::chrono::steady_clock;

int main(int argc, char** argv) {
    try {
        if (argc != 8) throw std::invalid_argument(
            "usage: alod_stability E1|E2 k level block seed iterations tolerance");
        const std::string problem = argv[1];
        const double k = std::stod(argv[2]), tolerance = std::stod(argv[7]);
        const int level = std::stoi(argv[3]), block = std::stoi(argv[4]);
        const int seed = std::stoi(argv[5]), maximum_iterations = std::stoi(argv[6]);
        if (!(k > 0) || !std::isfinite(k) || level < 0 || level > 21
            || block < 1 || block > 16 || maximum_iterations < 1
            || !(tolerance > 0) || !std::isfinite(tolerance))
            throw std::invalid_argument("invalid stability parameters");
        const auto start = Clock::now();
        const auto seconds = [&] { return std::chrono::duration<double>(Clock::now()-start).count(); };
        auto mesh = alod::make_problem(problem, k).initial_mesh;
        // Discard each prolongation immediately; this diagnostic only needs V_h.
        for (int i=0; i<level; ++i) mesh = lod2d::refine_nvb(mesh).mesh;
        Sparse A, E;
        int nodes = mesh.nodes.size();
        const double h = alod::mesh_diameter(mesh);
        {
            auto op = lod2d::helmholtz::assemble_helmholtz_operators(mesh, k);
            std::vector<int> free(nodes, 0);
            for (int i : op.dirichlet_nodes) free[i] = -1;
            int n=0;
            for (int& i : free) if (i != -1) i=n++;
            auto restrict_matrix = [&](const Sparse& full) {
                std::vector<Eigen::Triplet<std::complex<double>>> entries;
                entries.reserve(full.nonZeros());
                for (int j=0; j<full.outerSize(); ++j)
                    for (Sparse::InnerIterator it(full,j); it; ++it)
                        if (free[it.row()]>=0 && free[it.col()]>=0)
                            entries.emplace_back(free[it.row()],free[it.col()],it.value());
                Sparse result(n,n); result.setFromTriplets(entries.begin(),entries.end());
                return result;
            };
            A = restrict_matrix(op.system);
            E = restrict_matrix((op.stiffness+k*k*op.mass).cast<std::complex<double>>());
        }
        mesh = {};
        const int n=A.rows();
        if (block>n) throw std::invalid_argument("block exceeds free dimension");
        const double assembly_seconds=seconds();
        // A is complex symmetric, not Hermitian. Its adjoint inverse is obtained
        // by conjugation, so forward and adjoint actions share one factorization.
        if ((A-Sparse(A.transpose())).norm()>1e-13*A.norm())
            throw std::runtime_error("adjoint reuse requires complex symmetry");
        std::cerr<<"factorizing n="<<n<<" nonzeros="<<A.nonZeros()<<'\n';
        Eigen::UmfPackLU<Sparse> lu;
        if(n>256) lu.umfpackControl()[UMFPACK_ORDERING]=UMFPACK_ORDERING_METIS;
        lu.compute(A);
        if (lu.info()!=Eigen::Success) {
            lu.umfpackControl()[UMFPACK_PRL]=3;
            lu.printUmfpackStatus(); lu.printUmfpackInfo();
            throw std::runtime_error("Helmholtz factorization failed");
        }
        const double factor_seconds=seconds()-assembly_seconds;
        auto apply = [&](const Matrix& x)->Matrix {
            Matrix rhs=E*x;
            Matrix conjugate_rhs=rhs.conjugate();
            Matrix adjoint=lu.solve(conjugate_rhs).conjugate();
            Matrix result=lu.solve(Matrix(E*adjoint));
            if (lu.info()!=Eigen::Success || !result.allFinite())
                throw std::runtime_error("inverse action failed");
            return result;
        };
        // T=A^{-1} E A^{-*} E is self-adjoint in the E inner product;
        // lambda_max(T)=gamma_h^{-2}. Never form an inverse or normal matrix.
        auto orthonormalize = [&](const Matrix& candidates) {
            Matrix q(n,candidates.cols()), eq(n,candidates.cols()); int count=0;
            for (int j=0; j<candidates.cols(); ++j) {
                Vector v=candidates.col(j), ev=E*v;
                double original=std::sqrt(std::max(0.,std::real(v.dot(ev))));
                if (!(original>0)) continue;
                v/=original; ev/=original;
                for (int pass=0; pass<2; ++pass) {
                    ev=E*v;
                    for (int i=0; i<count; ++i) {
                        auto coefficient=q.col(i).dot(ev);
                        v-=coefficient*q.col(i); ev-=coefficient*eq.col(i);
                    }
                }
                ev=E*v;
                double norm=std::sqrt(std::max(0.,std::real(v.dot(ev))));
                if (norm<1e-10) continue;
                q.col(count)=v/norm; eq.col(count)=ev/norm; ++count;
            }
            return Matrix(q.leftCols(count));
        };
        std::mt19937 rng(seed); std::normal_distribution<double> normal;
        Matrix random(n,block);
        for (int j=0;j<block;++j) for(int i=0;i<n;++i) random(i,j)={normal(rng),normal(rng)};
        Matrix x=orthonormalize(random), previous(n,0), tx=apply(x);
        double lambda=0, residual=1; int iteration=0;
        std::cerr<<std::setprecision(10)<<"assembled n="<<n<<" factor_seconds="<<factor_seconds<<'\n';
        for (iteration=1; iteration<=maximum_iterations; ++iteration) {
            // Refresh actions: propagating them through nearly dependent trial
            // bases can otherwise accumulate a false residual floor.
            tx=apply(x);
            Matrix small=x.adjoint()*(E*tx);
            small=(.5*(small+small.adjoint())).eval();
            Eigen::SelfAdjointEigenSolver<Matrix> ritz(small);
            if (ritz.info()!=Eigen::Success) throw std::runtime_error("Ritz solve failed");
            x=(x*ritz.eigenvectors()).eval(); tx=(tx*ritz.eigenvectors()).eval();
            lambda=ritz.eigenvalues().tail(1)(0);
            Vector r=tx.col(block-1)-lambda*x.col(block-1);
            residual=std::sqrt(std::max(0.,std::real(r.dot(E*r))))/lambda;
            if(iteration==1 || iteration%10==0)
                std::cerr<<"iteration="<<iteration<<" gamma="<<1/std::sqrt(lambda)<<" residual="<<residual<<'\n';
            if(residual<=tolerance) break;
            Matrix candidates(n,2*block+previous.cols()); candidates<<x,tx,previous;
            Matrix basis=orthonormalize(candidates), tb=apply(basis);
            small=basis.adjoint()*(E*tb); small=(.5*(small+small.adjoint())).eval();
            Eigen::SelfAdjointEigenSolver<Matrix> expanded(small);
            if(expanded.info()!=Eigen::Success) throw std::runtime_error("expanded Ritz solve failed");
            Matrix selected=expanded.eigenvectors().rightCols(block);
            previous=x; x=basis*selected; tx=tb*selected;
        }
        if (residual>tolerance || !(lambda>0)) throw std::runtime_error("inverse eigensolve did not converge");
        // One inverse polish damps high-frequency contamination before checking
        // the original problem, whose small eigenvalue amplifies that component.
        Matrix polished=apply(Matrix(x.rightCols(1)));
        polished/=std::sqrt(std::real(polished.col(0).dot(E*polished.col(0))));
        const double eigen_seconds=seconds()-assembly_seconds-factor_seconds;
        // Validate against the original generalized eigenproblem, independently
        // of the inverse iteration's stopping test.
        Eigen::SimplicialLDLT<Sparse> energy; energy.compute(E);
        if(energy.info()!=Eigen::Success) throw std::runtime_error("energy factorization failed");
        Vector v=polished.col(0), av=A*v, z=energy.solve(av);
        const double gamma_squared=1/lambda;
        Vector r=A.adjoint()*z-gamma_squared*(E*v), er=energy.solve(r);
        double original_residual=std::sqrt(std::max(0.,std::real(r.dot(er))))/gamma_squared;
        double dense_gamma=0, dense_difference=0;
        if(n<=256) {
            Eigen::LLT<Matrix> chol{Matrix(E)};
            Matrix lower=chol.matrixL();
            Matrix inverse_lower=lower.triangularView<Eigen::Lower>().solve(Matrix::Identity(n,n));
            Matrix whitened=inverse_lower*Matrix(A)*inverse_lower.adjoint();
            Eigen::JacobiSVD<Matrix> svd(whitened);
            dense_gamma=svd.singularValues().tail(1)(0);
            dense_difference=std::abs(dense_gamma-std::sqrt(gamma_squared))/dense_gamma;
            if(dense_difference>1e-7) throw std::runtime_error("dense SVD cross-check failed");
        }
        if(!std::isfinite(original_residual) || original_residual>1e-5)
            throw std::runtime_error("original generalized residual failed: "+std::to_string(original_residual));
        std::cout<<std::setprecision(17)<<"{\"problem\":\""<<problem<<"\",\"k\":"<<k
          <<",\"level\":"<<level<<",\"nodes\":"<<nodes<<",\"free_dofs\":"<<n<<",\"h\":"<<h
          <<",\"block\":"<<block<<",\"seed\":"<<seed<<",\"gamma\":"<<std::sqrt(gamma_squared)
          <<",\"index_bits\":"<<8*sizeof(SuiteSparse_long)<<",\"ordering\":\""<<(n>256?"metis":"amd")<<"\""
          <<",\"inverse_residual\":"<<residual<<",\"original_residual\":"<<original_residual
          <<",\"iterations\":"<<iteration<<",\"dense_gamma\":"<<dense_gamma<<",\"dense_relative_difference\":"<<dense_difference
          <<",\"assembly_seconds\":"<<assembly_seconds<<",\"factor_seconds\":"<<factor_seconds
          <<",\"eigen_seconds\":"<<eigen_seconds<<",\"wall_seconds\":"<<seconds()<<"}\n";
    } catch(const std::exception& e) { std::cerr<<"alod_stability: "<<e.what()<<'\n'; return 1; }
}
