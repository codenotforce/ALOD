#include "fixed_support.hpp"
#include "alod/localization.hpp"
#include <chrono>

int main(int argc,char** argv) {
 try {
    if(argc==2&&std::string(argv[1])=="--help") {
        std::cout<<"Usage: alod_fixed E1|E2 --members=TABLE --training-ids=0,1,... [--level=2 --gap=3 --ell=1 --graded=0 --interpolation=area|arithmetic --riesz-patches=n2|archive --threads=1 --maximum-nodes=20000 --theta=0.15 --ritz-tolerance=0.0001 --ritz-iterations=750 --dense-threshold=64]\n";return 0;
    }
    const auto start=std::chrono::steady_clock::now();
    auto input=fixed::parse(argc,argv);auto coarse=fixed::coarse_mesh(input);
    auto fine=lod2d::refine_mesh_nvb(coarse,input.gap);
    alod::LodLimits limits;limits.threads=input.threads;limits.maximum_reference_nodes=input.cap;
    alod::LodSpace space(coarse,std::move(fine),input.wavenumber,input.ell,input.policy=="area"
        ?alod::InterpolationPolicy::ManuscriptAreaWeighted:alod::InterpolationPolicy::ArchivedArithmetic,limits);
    using namespace alod;
    const auto quad=paper_quadrature(input.problem);
    ComplexMatrix loads(space.fine().nodes.size(),input.members.size());
    std::vector<FamilyMember> members;
    for(int j=0;j<loads.cols();++j) {
        const auto& m=input.members[j];
        loads.col(j)=lod2d::helmholtz::assemble_helmholtz_load(space.fine(),m.problem.source,quad,m.problem.quadrature_context);
        members.push_back({m.id,m.role=="train"?MemberRole::Train:m.role=="test"?MemberRole::Test:m.role=="shift"?MemberRole::Shift:MemberRole::Pure});
    }
    auto solution=space.solve(loads);AdditiveKernelRieszContext riesz(space,input.riesz_policy=="n2"?RieszPatchPolicy::ManuscriptN2:RieszPatchPolicy::ArchivedSupportExpanded);
    auto estimate=riesz.estimate(loads,solution.values,input.theta);
    Eigen::VectorXd energies(loads.cols());
    for(int j=0;j<loads.cols();++j)energies(j)=std::sqrt(std::max(0.0,std::real(solution.values.col(j).dot(space.energy().cast<Complex>()*solution.values.col(j)))));
    auto marks=mark_family(estimate.element_eta_squared,energies,members,input.ids,input.theta);
    Eigen::MatrixXd strong(space.fine().elems.size(),loads.cols());
    for(int j=0;j<loads.cols();++j) {
        auto e=reference_strong_residual(space,solution.values.col(j),loads.col(j),input.members[j].problem.source,quad,input.members[j].problem.quadrature_context);
        for(int i=0;i<strong.rows();++i)strong(i,j)=e.fine.element_squared[i];
    }
    LocalizationEigenConfig cfg;cfg.relative_tolerance=input.tolerance;cfg.maximum_iterations=input.iterations;cfg.dense_cross_check_max_dimension=input.dense;
    auto theta=localization_theta(space,riesz,cfg);
    int nominal=-1;for(int j=0;j<loads.cols();++j)if(input.members[j].id==0)nominal=j;
    auto ratio=nominal>=0?localization_ratio(theta.theta,estimate.eta(nominal)):LocalizationRatio{"nominal_missing",std::nullopt};
    std::cout<<std::setprecision(17)<<"{\"problem\":\""<<input.problem<<"\",\"interpolation\":\""<<input.policy
        <<"\",\"riesz_patches\":\""<<input.riesz_policy<<"\",\"identity\":\""<<space.identity()<<"\",\"coarse_nodes\":"<<coarse.nodes.size()<<",\"reference_nodes\":"<<space.fine().nodes.size()
        <<",\"pg_residual\":"<<solution.pg_relative_residual<<",\"patch_residual\":"<<space.patch_residual()<<",\"constraint_residual\":"<<space.constraint_residual()
        <<",\"theta\":"<<theta.theta<<",\"ritz_residual\":"<<theta.spectrum.relative_residual<<",\"ritz_iterations\":"<<theta.spectrum.iterations
        <<",\"ritz_converged\":true,\"ritz_used_warm_start\":"<<(theta.spectrum.used_warm_start?"true":"false")
        <<",\"ritz_dense\":"<<(theta.spectrum.dense_cross_checked?"true":"false")<<",\"ritz_identity\":\""<<theta.identity<<"\""
        <<",\"riesz_factorizations\":"<<riesz.factorizations()<<",\"riesz_geometric_vertices\":"<<riesz.patches().size()
        <<",\"worst_member_id\":"<<marks.worst_member_id<<",\"mean_bulk\":"<<marks.mean_bulk<<",\"worst_bulk\":"<<marks.worst_bulk;
    std::cout<<",\"training_ids\":";fixed::vector(marks.member_ids);
    std::cout<<",\"ratio_status\":\""<<ratio.status<<"\",\"ratio\":";
    if(ratio.value)std::cout<<*ratio.value;else std::cout<<"null";
    std::cout<<",\"nominal_eta\":";if(nominal>=0)std::cout<<estimate.eta(nominal);else std::cout<<"null";
    std::cout<<",\"ritz_kind\":\""<<(theta.spectrum.dense_cross_checked?"dense_power":"block_ritz")<<"\",\"ritz_side\":\"adjoint\"";
    std::cout<<",\"ritz_dense_lambda\":";if(theta.spectrum.dense_cross_checked)std::cout<<theta.spectrum.dense_lambda_max;else std::cout<<"null";
    std::cout<<",\"ritz_dense_difference\":";if(theta.spectrum.dense_cross_checked)std::cout<<theta.spectrum.dense_relative_difference;else std::cout<<"null";
    std::cout<<",\"marked_elements\":";fixed::vector(marks.marked_elements);
    std::cout<<",\"frozen_scales\":";fixed::vector(marks.frozen_scales);
    std::cout<<",\"eta\":";fixed::vector(estimate.eta);
    std::cout<<",\"node_squared\":";fixed::matrix(estimate.node_eta_squared);
    std::cout<<",\"element_squared\":";fixed::matrix(estimate.element_eta_squared);
    std::cout<<",\"reference_strong_squared\":";fixed::matrix(strong);
    std::cout<<",\"patch_elements\":[";for(int z=0;z<static_cast<int>(riesz.patches().size());++z){if(z)std::cout<<',';fixed::vector(riesz.patches()[z].coarse_elements);}std::cout<<']';
    std::cout<<",\"solution\":";fixed::complex_matrix(solution.values);
    std::cout<<",\"wall_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 } catch(const std::exception& e){std::cerr<<"alod_fixed: "<<e.what()<<'\n';return 1;}
}
