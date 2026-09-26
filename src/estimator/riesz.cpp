#include "alod/estimator.hpp"
#include "../lod/fingerprint.hpp"
#include "helmholtz/boundary.h"
#include <Eigen/QR>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <chrono>
#include <mutex>
#include <unordered_map>
#include <omp.h>
namespace alod {
using namespace lod2d;
using namespace lod2d::helmholtz;
using lod2d::helmholtz::adaptive::mark_doerfler;
namespace {
constexpr double kSupportTolerance=1e-13, kRankTolerance=2e-12;
struct ConstraintReduction {
    std::vector<int> active_rows,independent_rows;
    Eigen::MatrixXd matrix;
};
std::vector<int> sorted_unique(std::vector<int> values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

std::vector<std::vector<int>> coarse_element_graph(const TriMesh &mesh) {
    std::vector<std::vector<int>> incident(mesh.nodes.size());
    for (int element = 0; element < static_cast<int>(mesh.elems.size()); ++element) {
        for (int node : mesh.elems[element]) incident[node].push_back(element);
    }
    std::vector<std::vector<int>> graph(mesh.elems.size());
    for (const std::vector<int> &star : incident) {
        for (int first : star) {
            graph[first].insert(graph[first].end(), star.begin(), star.end());
        }
    }
    for (std::vector<int> &neighbors : graph) neighbors = sorted_unique(std::move(neighbors));
    return graph;
}

std::vector<std::vector<int>> coarse_hat_support(const TriMesh &mesh) {
    std::vector<std::vector<int>> result(mesh.nodes.size());
    for (int element = 0; element < static_cast<int>(mesh.elems.size()); ++element) {
        for (int node : mesh.elems[element]) result[node].push_back(element);
    }
    for (std::vector<int> &support : result) support = sorted_unique(std::move(support));
    return result;
}

std::vector<int> expand_elements(
    const std::vector<std::vector<int>> &graph,
    const std::vector<int> &initial,
    int layers) {
    if (layers < 0) throw std::invalid_argument("patch layer count must be nonnegative");
    std::vector<char> selected(graph.size(), false);
    std::vector<int> frontier;
    for (int element : initial) {
        if (element < 0 || element >= static_cast<int>(graph.size()))
            throw std::out_of_range("patch contains an out-of-range coarse element");
        if (!selected[element]) {
            selected[element] = true;
            frontier.push_back(element);
        }
    }
    for (int layer = 0; layer < layers && !frontier.empty(); ++layer) {
        std::vector<int> next;
        for (int element : frontier) {
            for (int neighbor : graph[element]) {
                if (selected[neighbor]) continue;
                selected[neighbor] = true;
                next.push_back(neighbor);
            }
        }
        frontier = std::move(next);
    }
    std::vector<int> result;
    for (int element = 0; element < static_cast<int>(selected.size()); ++element) {
        if (selected[element]) result.push_back(element);
    }
    return result;
}

ConstraintReduction reduce_constraints(
    const Eigen::SparseMatrix<double> &interpolation,
    const std::vector<int> &discrete_dofs) {
    // Only a small number of coarse interpolation rows touch a local patch.
    // The former implementation nevertheless allocated and cleared an
    // interpolation.rows() by patch_dofs dense matrix for every patch.  On a
    // deep adaptive hierarchy this turned Gram-structure preparation into an
    // O(N_H * sum patch_dofs) memory-bandwidth operation.  Gather the touched
    // rows first and keep all dense work patch-local.
    std::vector<int> candidate_rows;
    candidate_rows.reserve(discrete_dofs.size());
    for (int local = 0; local < static_cast<int>(discrete_dofs.size()); ++local) {
        const int discrete = discrete_dofs[local];
        if (discrete < 0 || discrete >= interpolation.cols())
            throw std::out_of_range("local Riesz degree of freedom is out of range");
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 interpolation, discrete); it; ++it) {
            candidate_rows.push_back(it.row());
        }
    }
    std::sort(candidate_rows.begin(), candidate_rows.end());
    candidate_rows.erase(
        std::unique(candidate_rows.begin(), candidate_rows.end()),
        candidate_rows.end());

    ConstraintReduction result;
    if (candidate_rows.empty()) {
        result.matrix.resize(0, discrete_dofs.size());
        return result;
    }
    Eigen::MatrixXd candidates = Eigen::MatrixXd::Zero(
        candidate_rows.size(), discrete_dofs.size());
    for (int local = 0; local < static_cast<int>(discrete_dofs.size()); ++local) {
        const int discrete = discrete_dofs[local];
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 interpolation, discrete); it; ++it) {
            const auto found = std::lower_bound(
                candidate_rows.begin(), candidate_rows.end(), it.row());
            candidates(
                static_cast<int>(found - candidate_rows.begin()), local) =
                it.value();
        }
    }
    std::vector<int> active_candidate_rows;
    active_candidate_rows.reserve(candidate_rows.size());
    for (int row = 0; row < candidates.rows(); ++row) {
        if (candidates.row(row).norm() <= kSupportTolerance) continue;
        active_candidate_rows.push_back(row);
        result.active_rows.push_back(candidate_rows[row]);
    }
    if (result.active_rows.empty()) {
        result.matrix.resize(0, discrete_dofs.size());
        return result;
    }

    Eigen::MatrixXd active(active_candidate_rows.size(), discrete_dofs.size());
    for (int row = 0; row < static_cast<int>(result.active_rows.size()); ++row)
        active.row(row) = candidates.row(active_candidate_rows[row]);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(active.transpose());
    qr.setThreshold(kRankTolerance);
    const int rank = qr.rank();
    const Eigen::VectorXi permutation = qr.colsPermutation().indices();
    result.matrix.resize(rank, discrete_dofs.size());
    for (int row = 0; row < rank; ++row) {
        const int active_row = permutation(row);
        const int original_row = result.active_rows[active_row];
        result.independent_rows.push_back(original_row);
        result.matrix.row(row) = active.row(active_row);
    }
    return result;
}

Eigen::SparseMatrix<double> energy_matrix(const HelmholtzOperators &operators) {
    Eigen::SparseMatrix<double> result = operators.stiffness;
    result += operators.wavenumber * operators.wavenumber * operators.mass;
    result.makeCompressed();
    return result;
}

Eigen::SparseMatrix<double> restrict_sparse_matrix(
    const Eigen::SparseMatrix<double> &matrix,
    const std::vector<int> &dofs) {
    // Patch restrictions must scale with the patch, not the global reference
    // space.  Avoid allocating and clearing matrix.rows() integers for every
    // patch.  Kernel-patch DOFs are sorted; retain a map fallback for the few
    // generic callers that may provide another ordering.
    const bool sorted_unique = std::is_sorted(dofs.begin(), dofs.end())
        && std::adjacent_find(dofs.begin(), dofs.end()) == dofs.end();
    std::unordered_map<int, int> local_index;
    if (!sorted_unique) local_index.reserve(2 * dofs.size() + 1);
    for (int local = 0; local < static_cast<int>(dofs.size()); ++local) {
        if (dofs[local] < 0 || dofs[local] >= matrix.rows())
            throw std::out_of_range("matrix restriction degree of freedom is out of range");
        if (!sorted_unique) local_index[dofs[local]] = local;
    }
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(dofs.size() * 7);
    for (int local_column = 0; local_column < static_cast<int>(dofs.size()); ++local_column) {
        const int global_column = dofs[local_column];
        for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, global_column); it; ++it) {
            int local_row = -1;
            if (sorted_unique) {
                const auto found = std::lower_bound(
                    dofs.begin(), dofs.end(), it.row());
                if (found != dofs.end() && *found == it.row())
                    local_row = static_cast<int>(found - dofs.begin());
            } else {
                const auto found = local_index.find(it.row());
                if (found != local_index.end()) local_row = found->second;
            }
            if (local_row >= 0)
                triplets.emplace_back(local_row, local_column, it.value());
        }
    }
    Eigen::SparseMatrix<double> result(dofs.size(), dofs.size());
    result.setFromTriplets(triplets.begin(), triplets.end());
    result.makeCompressed();
    return result;
}


std::vector<KernelRieszPatch> build_kernel_riesz_patches(const LodSpace& space,RieszPatchPolicy policy) {
    const auto& coarse=space.coarse();const auto& fine=space.fine();
    const auto graph=coarse_element_graph(coarse), stars=coarse_hat_support(coarse);
    std::vector<std::vector<int>> children(coarse.elems.size());
    const Eigen::SparseMatrix<double,Eigen::RowMajor> parents(space.element_prolongation());
    for(int e=0;e<parents.rows();++e) {
        int count=0;
        for(decltype(parents)::InnerIterator it(parents,e);it;++it) {
            if(it.value()!=1)throw std::invalid_argument("non-integral reference parent map");
            children[it.col()].push_back(e);++count;
        }
        if(count!=1)throw std::invalid_argument("reference element needs exactly one parent");
    }
    std::vector<int> incidence(fine.nodes.size(),0);
    for(const auto& element:fine.elems)for(int node:element)++incidence[node];
    std::vector<bool> dirichlet(fine.nodes.size(),false);
    for(int node:dirichlet_nodes(fine))dirichlet[node]=true;
    int layers=1;
    if(policy==RieszPatchPolicy::ArchivedSupportExpanded) {
        std::vector<std::vector<int>> node_parents(fine.nodes.size()),support(coarse.nodes.size());
        for(int parent=0;parent<static_cast<int>(children.size());++parent)
            for(int e:children[parent])for(int node:fine.elems[e])node_parents[node].push_back(parent);
        for(auto& v:node_parents)v=sorted_unique(std::move(v));
        for(int node=0;node<space.interpolation().outerSize();++node)
            for(Sparse::InnerIterator it(space.interpolation(),node);it;++it)if(std::abs(it.value())>kSupportTolerance)
                support[it.row()].insert(support[it.row()].end(),node_parents[node].begin(),node_parents[node].end());
        int propagation=0;
        for(int z=0;z<static_cast<int>(support.size());++z) {
            support[z]=sorted_unique(std::move(support[z]));
            int distance=0;auto reached=stars[z];
            while(!std::includes(reached.begin(),reached.end(),support[z].begin(),support[z].end())) {
                if(++distance>static_cast<int>(coarse.elems.size()))throw std::runtime_error("unreachable interpolation support");
                reached=expand_elements(graph,reached,1);
            }
            propagation=std::max(propagation,distance);
        }
        layers=propagation+1;
    } else if(policy!=RieszPatchPolicy::ManuscriptN2)throw std::invalid_argument("invalid Riesz patch policy");
    std::vector<KernelRieszPatch> patches(coarse.nodes.size());
    std::size_t total_entries=0;
    // A star is N^1(z); one vertex-adjacency expansion is N^2(z).
    // Include geometric Dirichlet vertices as patch centers as well.
    for(int z=0;z<static_cast<int>(coarse.nodes.size());++z) {
        auto& patch=patches[z];patch.coarse_node=z;patch.coarse_hat_support=stars[z];
        patch.coarse_elements=expand_elements(graph,stars[z],layers);
        std::vector<int> nodes;
        for(int parent:patch.coarse_elements)for(int e:children[parent])
            for(int node:fine.elems[e])nodes.push_back(node);
        std::sort(nodes.begin(),nodes.end());
        for(std::size_t first=0;first<nodes.size();) {
            std::size_t last=first+1;
            while(last<nodes.size() && nodes[last]==nodes[first])++last;
            int node=nodes[first];
            if(!dirichlet[node] && last-first==static_cast<std::size_t>(incidence[node]))patch.discrete_dofs.push_back(node);
            first=last;
        }
        // Conservative bound before creating the local dense constraint matrix.
        std::vector<int> touched_rows;
        for(int node:patch.discrete_dofs)for(Sparse::InnerIterator it(space.interpolation(),node);it;++it)touched_rows.push_back(it.row());
        touched_rows=sorted_unique(std::move(touched_rows));
        total_entries+=touched_rows.size()*patch.discrete_dofs.size();
        if(total_entries>space.limits().maximum_patch_entries)
            throw std::runtime_error("Riesz patch storage resource limit exceeded");
        auto reduced=reduce_constraints(space.interpolation(),patch.discrete_dofs);
        patch.constraints=std::move(reduced.matrix);
        patch.active_constraint_rows=std::move(reduced.active_rows);
        patch.independent_constraint_rows=std::move(reduced.independent_rows);
    }
    return patches;
}
} // namespace
// Batched Schur/saddle implementation extracted from the archived E2 production source.
struct AdditiveKernelRieszContext::Impl {
    struct Group {
        KernelRieszPatch patch;
        std::vector<int> nodes;
        Eigen::SparseMatrix<double> energy;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> energy_factor;
        Eigen::MatrixXd inverse_constraints;
        Eigen::LDLT<Eigen::MatrixXd> schur_factor;
        Eigen::SparseLU<ComplexSparseMatrix> saddle_factor;
        bool schur = true;
    };
    std::vector<std::unique_ptr<Group>> groups;
    std::vector<KernelRieszPatch> patches;
    std::shared_ptr<const LodHierarchyData> hierarchy;
    const TriMesh& mesh;
    const ComplexSparseMatrix& system;
    const Eigen::SparseMatrix<double>& interpolation;
    explicit Impl(std::shared_ptr<const LodHierarchyData> h):hierarchy(std::move(h)),
        mesh(hierarchy->coarse),system(hierarchy->operators.system),interpolation(hierarchy->interpolation){}
    std::vector<int> dirichlet;
    int full_size = 0, coarse_size = 0, elements = 0, threads = 1;
    std::string identity,policy_name;
    std::size_t dense_limit=0;
    double preparation_seconds = 0.0;
    std::size_t calls = 0,columns = 0;

    template<class F> void parallel(F fn) {
        std::exception_ptr failure;
        std::mutex mutex;
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
        for (int i = 0; i < static_cast<int>(groups.size()); ++i) {
            try { fn(i); }
            catch (...) {
                std::lock_guard<std::mutex> lock(mutex);
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
    }
};

AdditiveKernelRieszContext::AdditiveKernelRieszContext(
    const LodSpace &space,RieszPatchPolicy policy)
    : impl_(std::make_unique<Impl>(space.hierarchy())) {
    const auto start = std::chrono::steady_clock::now();
    auto &p = *impl_;
    const auto& operators=space.operators();
    const int maximum_threads=space.limits().threads;
    p.identity=space.reference_identity();
    p.policy_name=policy==RieszPatchPolicy::ManuscriptN2?"manuscript-n2":"archived-support-expanded";
    p.dense_limit=space.limits().maximum_dense_entries;
    p.full_size = operators.system.rows();
    p.coarse_size = space.coarse().nodes.size();
    p.elements = space.coarse().elems.size();
    p.dirichlet = operators.dirichlet_nodes;
    p.patches = build_kernel_riesz_patches(space,policy);
    const auto& energy = space.energy();
    std::unordered_map<std::string, std::vector<int>> lookup;
    for (const auto &patch : p.patches) {
        if (patch.discrete_dofs.empty()) continue;
        FingerprintBuilder fp;
        for (int d : patch.discrete_dofs) fp.add_i64(d);
        fp.add_i64(patch.constraints.rows());
        fp.add_i64(patch.constraints.cols());
        for (int j = 0; j < patch.constraints.size(); ++j)
            fp.add_double(patch.constraints.data()[j]);
        auto &candidates = lookup[fp.finish()];
        int found = -1;
        for (int index : candidates) {
            const auto &other = p.groups[index]->patch;
            if (other.discrete_dofs == patch.discrete_dofs
                && other.constraints.rows() == patch.constraints.rows()
                && other.constraints.cols() == patch.constraints.cols()
                && (other.constraints.array() == patch.constraints.array()).all()) {
                found = index; break;
            }
        }
        if (found < 0) {
            found = p.groups.size();
            auto group = std::make_unique<Impl::Group>();
            group->patch = patch;
            p.groups.push_back(std::move(group));
            candidates.push_back(found);
        }
        p.groups[found]->nodes.push_back(patch.coarse_node);
    }
#ifdef _OPENMP
    p.threads = omp_get_max_threads();
#endif
    if (maximum_threads > 0) p.threads = std::min(p.threads, maximum_threads);
    p.threads = std::max(1, std::min(p.threads, static_cast<int>(p.groups.size())));
    p.parallel([&](int index) {
        auto &g = *p.groups[index];
        g.energy = restrict_sparse_matrix(energy, g.patch.discrete_dofs);
        const auto &c = g.patch.constraints;
        const int n = g.energy.rows(), m = c.rows();
        if(m==n)return; // The constrained space is exactly zero; no factor is used.
        g.schur = m <= 256;
        if (g.schur) {
            g.energy_factor.compute(g.energy);
            if (g.energy_factor.info() != Eigen::Success)
                throw std::runtime_error("AS energy factorization failed");
            if (m) {
                g.inverse_constraints = g.energy_factor.solve(c.transpose());
                g.schur_factor.compute(c * g.inverse_constraints);
                if (g.schur_factor.info() != Eigen::Success)
                    throw std::runtime_error("AS Schur factorization failed");
            }
        } else {
            std::vector<ComplexTriplet> entries;
            for (int j = 0; j < n; ++j)
                for (Eigen::SparseMatrix<double>::InnerIterator it(g.energy,j); it; ++it)
                    entries.emplace_back(it.row(),it.col(),it.value());
            for (int i = 0; i < m; ++i) for (int j = 0; j < n; ++j) {
                if (c(i,j) == 0.0) continue;
                entries.emplace_back(n+i,j,c(i,j));
                entries.emplace_back(j,n+i,c(i,j));
            }
            ComplexSparseMatrix saddle(n+m,n+m);
            saddle.setFromTriplets(entries.begin(),entries.end());
            g.saddle_factor.compute(saddle);
            if (g.saddle_factor.info() != Eigen::Success)
                throw std::runtime_error("AS saddle factorization failed");
        }
    });
    p.preparation_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
}

AdditiveKernelRieszContext::~AdditiveKernelRieszContext() = default;

std::vector<int> AdditiveKernelRieszContext::regional_mask(double radius) const {
    if(!std::isfinite(radius) || radius<0)throw std::invalid_argument("invalid region radius");
    std::vector<int> mask(impl_->coarse_size,0);
    for(const auto &patch:impl_->patches) {
        bool inside=!patch.discrete_dofs.empty();
        for(int el:patch.coarse_elements)
            for(int j=0;j<3;++j)
                if(impl_->mesh.nodes[impl_->mesh.elems[el][j]].norm()>radius+1e-13) inside=false;
        if(patch.coarse_elements.empty()) inside=false;
        mask[patch.coarse_node]=inside?1:0;
    }
    return mask;
}
AdditiveKernelRieszContext::Result AdditiveKernelRieszContext::apply(const ComplexMatrix &input) {
    return apply_selected(input,std::vector<int>(impl_->coarse_size,1),true);
}
AdditiveKernelRieszContext::Result AdditiveKernelRieszContext::apply_selected(
    const ComplexMatrix &input,const std::vector<int> &mask,bool full_estimator) {
    const auto start = std::chrono::steady_clock::now();
    auto &p = *impl_;
    if (input.rows() != p.full_size || input.cols()<1 || !input.allFinite()
        || static_cast<std::size_t>(input.size())>p.dense_limit)
        throw std::invalid_argument("AS residual dimension/value mismatch");
    if(mask.size()!=static_cast<std::size_t>(p.coarse_size)) throw std::invalid_argument("AS mask size");
    for(int v:mask) if(v!=0 && v!=1) throw std::invalid_argument("AS mask value");
    std::size_t local_entries=0;
    for(const auto& g:p.groups)local_entries+=g->patch.discrete_dofs.size()*input.cols();
    if(local_entries>p.dense_limit)throw std::runtime_error("Riesz batch workspace resource limit exceeded");
    ++p.calls;p.columns+=input.cols();
    ComplexMatrix rhs = input;
    for (int node : p.dirichlet) rhs.row(node).setZero();
    Result result;
    if(full_estimator)result.values = ComplexMatrix::Zero(p.full_size,input.cols());
    result.node_eta_squared = Eigen::MatrixXd::Zero(p.coarse_size,input.cols());
    result.selected_values = ComplexMatrix::Zero(p.full_size,input.cols());
    result.selected_eta = Eigen::VectorXd::Zero(input.cols());
    std::vector<int> selected_counts(p.groups.size(),0);
    for(int k=0;k<static_cast<int>(p.groups.size());++k)
        for(int node:p.groups[k]->nodes) selected_counts[k]+=mask[node];
    std::vector<ComplexMatrix> local_values(p.groups.size());
    std::vector<double> local_errors(p.groups.size(),0.0);
    std::vector<double> kernel_errors(p.groups.size(),0.0);
    p.parallel([&](int index) {
        auto &g = *p.groups[index];
        if(!full_estimator && selected_counts[index]==0) return;
        const auto &dofs = g.patch.discrete_dofs;
        const auto &c = g.patch.constraints;
        const int n = dofs.size(), m = c.rows();
        // Constraints are already reduced to independent rows. When m == n,
        // the admissible space is exactly {0}; its Riesz solution is zero.
        // Do not subtract Schur terms and test relative roundoff against zero.
        if (m == n) {
            local_values[index] = ComplexMatrix::Zero(n,input.cols());
            return;
        }
        ComplexMatrix r(n,input.cols()), x(n,input.cols());
        for (int i=0;i<n;++i) r.row(i)=rhs.row(dofs[i]);
        ComplexMatrix multipliers = ComplexMatrix::Zero(m,input.cols());
        if (g.schur) {
            x.real() = g.energy_factor.solve(r.real());
            x.imag() = g.energy_factor.solve(r.imag());
            if (m) {
                const ComplexMatrix cr = c.cast<Complex>() * x;
                multipliers.real() = g.schur_factor.solve(cr.real());
                multipliers.imag() = g.schur_factor.solve(cr.imag());
                x -= g.inverse_constraints.cast<Complex>() * multipliers;
            }
        } else {
            ComplexMatrix b = ComplexMatrix::Zero(n+m,input.cols());
            b.topRows(n)=r;
            const ComplexMatrix sol = g.saddle_factor.solve(b);
            x=sol.topRows(n); multipliers=sol.bottomRows(m);
        }
        if (!x.allFinite()) throw std::runtime_error("AS solve is not finite");
        const ComplexMatrix ex = g.energy.cast<Complex>()*x;
        kernel_errors[index]=(c.cast<Complex>()*x).norm()/std::max(1e-30,c.norm()*x.norm());
        local_errors[index]=(ex+c.transpose().cast<Complex>()*multipliers-r).norm()
            /std::max(1e-30,r.norm());
        for (int j=0;j<input.cols();++j) {
            const double sq=std::max(0.0,std::real(x.col(j).dot(ex.col(j))));
            for (int node:g.nodes) if(full_estimator || mask[node]) result.node_eta_squared(node,j)=sq;
        }
        local_values[index]=std::move(x);
    });
    // Deterministic reduction: no overlapping parallel writes or atomics.
    for (int k=0;k<static_cast<int>(p.groups.size());++k) {
        const auto &g=*p.groups[k];
        if(!full_estimator && selected_counts[k]==0) continue;
        for (int i=0;i<static_cast<int>(g.patch.discrete_dofs.size());++i) {
            if(full_estimator)result.values.row(g.patch.discrete_dofs[i]) +=
                static_cast<double>(g.nodes.size()) * local_values[k].row(i);
            result.selected_values.row(g.patch.discrete_dofs[i]) +=
                static_cast<double>(selected_counts[k]) * local_values[k].row(i);
        }
        result.local_relative_residual=std::max(result.local_relative_residual,local_errors[k]);
        result.constraint_relative_residual=std::max(result.constraint_relative_residual,kernel_errors[k]);
    }
    // For selected-only training both fields represent exactly the same sum.
    // Scatter once in the original deterministic order, then copy contiguously.
    if(!full_estimator)result.values=result.selected_values;
    result.eta.resize(input.cols());
    for (int j=0;j<input.cols();++j) {
        const double sq=result.node_eta_squared.col(j).sum();
        double selected_sq=0;
        for(int node=0;node<p.coarse_size;++node) if(mask[node]) selected_sq+=result.node_eta_squared(node,j);
        result.selected_eta(j)=std::sqrt(std::max(0.0,selected_sq));
        result.selected_identity_relative_error=std::max(result.selected_identity_relative_error,
            std::abs(rhs.col(j).dot(result.selected_values.col(j))-Complex(selected_sq,0))/std::max(1e-30,selected_sq));
        result.eta(j)=std::sqrt(std::max(0.0,sq));
        const Complex action=rhs.col(j).dot(result.values.col(j));
        result.identity_relative_error=std::max(result.identity_relative_error,
            std::abs(action-Complex(sq,0.0))/std::max(1e-30,sq));
    }
    result.constraint_relative_residual=std::max(result.constraint_relative_residual,
        (p.interpolation.cast<Complex>()*result.values).norm()/std::max(1e-30,p.interpolation.norm()*result.values.norm()));
    result.solve_seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
    return result;
}

ResidualRieszBatch AdditiveKernelRieszContext::estimate(
    const ComplexMatrix &loads,const ComplexMatrix &values,double theta) {
    auto &p=*impl_;
    ResidualRieszBatch result;
    if(loads.rows()!=p.full_size || values.rows()!=p.full_size || loads.cols()!=values.cols()
        || !loads.allFinite() || !values.allFinite())throw std::invalid_argument("Riesz estimate dimensions/values");
    result.global_residuals=loads-p.system*values;
    auto applied=apply(result.global_residuals);
    result.node_eta_squared=std::move(applied.node_eta_squared);
    result.eta=std::move(applied.eta);
    result.element_eta_squared=Eigen::MatrixXd::Zero(p.elements,loads.cols());
    for (const auto &patch:p.patches)
        for (int element:patch.coarse_hat_support)
            result.element_eta_squared.row(element)+=result.node_eta_squared.row(patch.coarse_node)
                /static_cast<double>(patch.coarse_hat_support.size());
    result.marked_elements.resize(loads.cols());
    for (int j=0;j<loads.cols();++j) {
        std::vector<double> mass(p.elements);
        for(int i=0;i<p.elements;++i) mass[i]=result.element_eta_squared(i,j);
        result.marked_elements[j]=mark_doerfler(mass,theta);
    }
    result.parallel_threads=p.threads;
    result.patch_solve_seconds=applied.solve_seconds;
    return result;
}
double AdditiveKernelRieszContext::factorization_seconds() const {return impl_->preparation_seconds;}
std::size_t AdditiveKernelRieszContext::factorizations() const {return impl_->groups.size();}
std::size_t AdditiveKernelRieszContext::applications() const {return impl_->calls;}
std::size_t AdditiveKernelRieszContext::applied_columns() const {return impl_->columns;}

const std::vector<KernelRieszPatch>& AdditiveKernelRieszContext::patches()const{return impl_->patches;}
const std::string& AdditiveKernelRieszContext::reference_identity()const{return impl_->identity;}
const std::string& AdditiveKernelRieszContext::patch_policy_name()const{return impl_->policy_name;}
}
