#include "alod/localization.hpp"
#include "alod/problems.hpp"
#include <iostream>
#include <iomanip>
using namespace alod;
int main(){try{
    auto p=make_problem("E1",16);
    auto H=lod2d::refine_mesh_nvb(p.initial_mesh,6).mesh;
    auto h=lod2d::refine_mesh_nvb(H,3);
    LodLimits limits;limits.threads=2;limits.maximum_reference_nodes=200000;
    limits.maximum_patch_entries=1000000000;limits.maximum_dense_entries=100000000;
    LodSpace s(H,h,16,3,InterpolationPolicy::ManuscriptAreaWeighted,limits);
    AdditiveKernelRieszContext r(s,RieszPatchPolicy::ManuscriptN2);
    LocalizationEigenConfig cfg;cfg.eigenvalue_relative_residual=true;cfg.dense_cross_check_max_dimension=0;
    auto insufficient=cfg;insufficient.maximum_iterations=1;
    bool rejected=false;try{(void)localization_theta(s,r,insufficient);}catch(const std::runtime_error&){rejected=true;}
    if(!rejected)throw std::runtime_error("unconverged nontrivial Ritz problem was accepted");
    auto cold=localization_theta(s,r,cfg);
    auto fused_config=cfg;fused_config.fused_defect=true;
    auto fused_result=localization_theta(s,r,fused_config);
    if(fused_result.spectrum.relative_residual>1e-9||std::abs(fused_result.theta/cold.theta-1)>1e-8)
        throw std::runtime_error("fused eigenproblem differs from parallel global path");
    auto reference_config=cfg;
    reference_config.reuse_operator_actions=false;reference_config.lightweight_riesz=false;
    reference_config.fused_defect=false;reference_config.parallel_defect=false;
    auto reference=localization_theta(s,r,reference_config);
    if(std::abs(cold.theta/reference.theta-1)>1e-8)
        throw std::runtime_error("cached/lightweight Theta differs from reference");
    if(cold.spectrum.operator_columns>=reference.spectrum.operator_columns)
        throw std::runtime_error("operator action reuse did not reduce solved columns");
    for(bool light:{false,true}){
        auto variant=cfg;variant.lightweight_riesz=light;variant.reuse_operator_actions=!light;
        variant.fused_defect=false;variant.parallel_defect=false;
        auto result=localization_theta(s,r,variant);
        if(result.spectrum.relative_residual>1e-9||std::abs(result.theta/reference.theta-1)>1e-8)
            throw std::runtime_error("individual Theta optimization changed result");
    }
    if(cold.spectrum.relative_residual>1e-9)throw std::runtime_error("cold residual exceeds 1e-9");
    cfg.warm_start=cold.warm_start;
    auto warm=localization_theta(s,r,cfg);
    if(warm.spectrum.relative_residual>1e-9)throw std::runtime_error("warm residual exceeds 1e-9");
    cfg.dense_cross_check_max_dimension=1024;
    auto dense=localization_theta(s,r,cfg);
    cfg.relative_tolerance=1e-4;
    auto unprojected=localization_theta(s,r,cfg);
    if(std::abs(dense.theta/unprojected.theta-1)>1e-7)throw std::runtime_error("dual kernel projection changed Theta");
    double err=std::abs(cold.theta/dense.theta-1),err2=std::abs(warm.theta/dense.theta-1);
    std::cout<<std::setprecision(17)<<"{\"cold\":"<<cold.theta<<",\"warm\":"<<warm.theta<<",\"dense\":"<<dense.theta
             <<",\"cold_iterations\":"<<cold.spectrum.iterations<<",\"warm_iterations\":"<<warm.spectrum.iterations
             <<",\"relative_difference\":"<<std::max(err,err2)<<"}\n";
    if(err>1e-6||err2>1e-6)throw std::runtime_error("relative matrix-free Theta differs from dense largest eigenvalue");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
