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
    auto cold=localization_theta(s,r,cfg);
    cfg.warm_start=cold.warm_start;
    auto warm=localization_theta(s,r,cfg);
    cfg.dense_cross_check_max_dimension=1024;
    auto dense=localization_theta(s,r,cfg);
    double err=std::abs(cold.theta/dense.theta-1),err2=std::abs(warm.theta/dense.theta-1);
    std::cout<<std::setprecision(17)<<"{\"cold\":"<<cold.theta<<",\"warm\":"<<warm.theta<<",\"dense\":"<<dense.theta
             <<",\"cold_iterations\":"<<cold.spectrum.iterations<<",\"warm_iterations\":"<<warm.spectrum.iterations
             <<",\"relative_difference\":"<<std::max(err,err2)<<"}\n";
    if(err>1e-6||err2>1e-6)throw std::runtime_error("relative matrix-free Theta differs from dense largest eigenvalue");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
