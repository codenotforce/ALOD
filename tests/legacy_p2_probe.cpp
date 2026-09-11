// Developer-only probe, linked exclusively to the archived numerical library.
#include "../apps/fixed_support.hpp"
#include "helmholtz/model.h"
#include "helmholtz/adaptive/hierarchy.h"
#include "helmholtz/adaptive/kernel_residual.h"
#include "helmholtz/adaptive/certificates.h"
#include "helmholtz/adaptive/estimator.h"

int main(int argc,char** argv) {
 try {
    auto input=fixed::parse(argc,argv);
    if(input.policy!="arithmetic"||input.riesz_policy!="archive")throw std::invalid_argument("archived oracle requires arithmetic interpolation and archive patches");
    using namespace lod2d;using namespace lod2d::helmholtz;using namespace adaptive;
    auto initial=alod::make_problem(input.problem).initial_mesh;
    auto coarse=fixed::coarse_mesh(input);auto fine=refine_mesh_nvb(coarse,input.gap).mesh;
    double root_area=compute_area(initial)[0];
    auto levels=[&](const TriMesh& mesh){std::vector<int> out;for(double a:compute_area(mesh))out.push_back(std::lround(std::log2(root_area/a)));return out;};
    ReferenceEpochHierarchy hierarchy(initial,coarse,levels(coarse),fine,levels(fine));
    HelmholtzProblemConfig config;config.H=input.level;config.h=input.level+input.gap;
    config.initial_mesh=initial;config.ell=input.ell;config.wavenumber=16;
    config.patch_solver.kind=HelmholtzPatchSolverKind::DirectSchur;config.quadrature=alod::paper_quadrature(input.problem);
    auto model=HelmholtzLodModel::build_adaptive(config,hierarchy);
    ComplexMatrix loads(fine.nodes.size(),input.members.size());
    for(int j=0;j<loads.cols();++j)loads.col(j)=assemble_helmholtz_load(fine,input.members[j].problem.source,config.quadrature,input.members[j].problem.quadrature_context);
    auto values=model.solve_loads(loads);
    AdditiveKernelRieszContext riesz(hierarchy,model.operators(),input.threads);
    auto estimate=riesz.estimate(loads,values,input.theta);
    Eigen::MatrixXd strong(fine.elems.size(),loads.cols());
    for(int j=0;j<loads.cols();++j){auto e=diagnostics::estimate_conforming_p1_residual(fine,model.operators(),values.col(j),loads.col(j),input.members[j].problem.source,config.quadrature,input.members[j].problem.quadrature_context);
        for(int i=0;i<strong.rows();++i)strong(i,j)=e.element_squared[i];}
    std::vector<int> nodes;std::set<int> boundary(coarse.dirichlet.begin(),coarse.dirichlet.end());
    for(int i=0;i<static_cast<int>(coarse.nodes.size());++i)if(!boundary.contains(i))nodes.push_back(i);
    LocalizationEigenConfig eigen;eigen.relative_tolerance=input.tolerance;eigen.maximum_iterations=input.iterations;
    eigen.dense_cross_check_max_dimension=input.dense;eigen.dense_fallback_max_dimension=0;
    auto theta=build_reference_corrector_certificate(hierarchy,model.operators(),model.test_basis(),nodes,KernelRieszSolver::SaddlePoint,eigen);
    Eigen::SparseMatrix<double> energy=model.operators().stiffness+256*model.operators().mass;
    Eigen::VectorXd scales(loads.cols());for(int j=0;j<loads.cols();++j)scales(j)=std::max(1e-12,std::sqrt(std::real(values.col(j).dot(energy.cast<Complex>()*values.col(j)))));
    std::cout<<std::setprecision(17)<<"{\"coarse_nodes\":"<<coarse.nodes.size()<<",\"reference_nodes\":"<<fine.nodes.size()<<",\"theta\":"<<theta.theta_loc;
    std::cout<<",\"eta\":";fixed::vector(estimate.eta);
    std::cout<<",\"node_squared\":";fixed::matrix(estimate.node_eta_squared);
    std::cout<<",\"element_squared\":";fixed::matrix(estimate.element_eta_squared);
    std::cout<<",\"reference_strong_squared\":";fixed::matrix(strong);
    auto patches=build_kernel_riesz_patches(hierarchy,KernelRieszSpace::ReferenceResidual,kernel_riesz_patch_policy(hierarchy,KernelRieszSpace::ReferenceResidual));
    std::cout<<",\"patch_elements\":[";for(int z=0;z<static_cast<int>(patches.size());++z){if(z)std::cout<<',';fixed::vector(patches[z].coarse_elements);}std::cout<<']';
    std::cout<<",\"solution\":";fixed::complex_matrix(values);
    std::cout<<",\"all_scales\":";fixed::vector(scales);std::cout<<"}\n";
 } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
