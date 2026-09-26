#include "alod/batch.hpp"
#include "alod/afem.hpp"
#include "mesh/refine.h"
#include <omp.h>
#include <iostream>
#include <stdexcept>
#include <atomic>

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void close(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    require((a-b).norm() <= 2e-13 * std::max(1.0, b.norm()), "integration changed");
}
int main() {
    using namespace alod;
    using namespace lod2d::helmholtz;
    for (const auto* id : {"E1", "E2"}) for (double k : {8., 128.}) {
        auto p = make_problem(id, k);
        auto mesh = lod2d::refine_mesh_nvb(p.initial_mesh, 5).mesh;
        auto q = paper_quadrature(id);
        std::vector<PhysicalTriangleQuadraturePoint> scratch;
        for (int element = 0; element < static_cast<int>(mesh.elems.size()); ++element) {
            auto expected = triangle_quadrature_points(mesh, element, q, p.quadrature_context);
            triangle_quadrature_points_into(mesh, element, q, p.quadrature_context, scratch);
            require(expected.size() == scratch.size(), "scratch point count changed");
            for (std::size_t i = 0; i < scratch.size(); ++i) {
                require(expected[i].point == scratch[i].point && expected[i].weight == scratch[i].weight
                    && expected[i].barycentric == scratch[i].barycentric, "scratch point changed");
                const auto jet = p.exact_jet(scratch[i].point);
                require(jet.first == p.exact(scratch[i].point)
                    && jet.second == p.exact_gradient(scratch[i].point), "joint analytic evaluation changed");
            }
        }
        auto ops = assemble_helmholtz_operators(mesh, k);
        Sparse energy = ops.stiffness + k*k*ops.mass;
        ComplexMatrix values(mesh.nodes.size(), 2), reference(mesh.nodes.size(), 2);
        for (int i = 0; i < values.rows(); ++i) {
            values(i,0) = p.exact(mesh.nodes[i]);
            values(i,1) = Complex(0.2,0.1) * values(i,0);
            reference(i,0) = Complex(0.9,0.05) * values(i,0);
            reference(i,1) = 0.;
        }
        std::vector<Problem> problems{p,p};
        // Separate scalar passes remain independent of the joint callback.
        auto expected = integrate_audit_batch(mesh,energy,values,reference,problems,q,1,false);
        for (int threads : {1, 2, 32}) {
            omp_set_num_threads(threads);
            auto actual = integrate_audit_batch(mesh,energy,values,reference,problems,q,threads);
            close(actual.exact_norm,expected.exact_norm);
            close(actual.exact_error,expected.exact_error);
            close(actual.reference_error,expected.reference_error);
            close(actual.energy,expected.energy);
            auto ordinary = integrate_error_batch(mesh,energy,values,problems,q,threads);
            close(ordinary.exact_error,expected.exact_error);
            close(ordinary.exact_norm,expected.exact_norm);
            require(ordinary.reference_error.size()==0,"ordinary error schema changed");
        }
        for (auto& problem : problems) problem.exact_jet = {};
        auto fallback = integrate_audit_batch(mesh,energy,values,reference,problems,q,2);
        close(fallback.exact_error,expected.exact_error);
        close(fallback.reference_error,expected.reference_error);
        auto load = assemble_load_batch(mesh,problems,q,2);
        omp_set_num_threads(1);
        require(load == assemble_load_batch(mesh,problems,q,1),"load thread reproducibility");
        // Compare the shared geometry context with independently rebuilt contexts.
        lod2d::helmholtz::adaptive::diagnostics::ResidualMeshContext geometry(mesh);
        for(int j=0;j<2;++j){
            auto cached=geometry.estimate(ops,values.col(j),load.col(j),p.source,q,p.quadrature_context);
            auto separate=lod2d::helmholtz::adaptive::diagnostics::estimate_conforming_p1_residual(mesh,ops,values.col(j),load.col(j),p.source,q,p.quadrature_context);
            require(cached.element_squared==separate.element_squared&&cached.eta==separate.eta,
                    "shared residual geometry changed the estimator");
        }
        // Distinct quadrature features must not be merged by the batch path.
        if(std::string(id)=="E2"){
            auto mixed=problems;mixed[1].quadrature_context.feature_point+=lod2d::Point2(0.1,0.05);
            auto scalar=integrate_audit_batch(mesh,energy,values,reference,mixed,q,1,false);
            auto batch=integrate_audit_batch(mesh,energy,values,reference,mixed,q,32);
            close(scalar.exact_error,batch.exact_error);
            auto mixed_load=assemble_load_batch(mesh,mixed,q,32);
            for(int j=0;j<2;++j)require(mixed_load.col(j)==assemble_helmholtz_load(mesh,mixed[j].source,q,mixed[j].quadrature_context),"mixed-context load changed");
        }
        std::atomic<unsigned long long> load_workers{0},error_workers{0};
        for(auto& problem:problems){
            auto source=problem.source;problem.source=[&,source](const lod2d::Point2& x){load_workers.fetch_or(1ULL<<omp_get_thread_num());return source(x);};
            auto jet=p.exact_jet;problem.exact_jet=[&,jet](const lod2d::Point2& x){error_workers.fetch_or(1ULL<<omp_get_thread_num());return jet(x);};
        }
        omp_set_dynamic(0);
        require(load==assemble_load_batch(mesh,problems,q,32),"32-thread load reproducibility");
        auto tiled=integrate_audit_batch(mesh,energy,values,reference,problems,q,32);
        close(tiled.exact_error,expected.exact_error);
        require(load_workers==0xffffffffULL&&error_workers==0xffffffffULL,"small RHS batches must use all 32 workers");
    }
    std::cout << "joint and separate integrals agree; reusable quadrature and thread reproducibility passed\n";
}
