#include "../apps/fixed_support.hpp"
#include "alod/batch.hpp"
#include "alod/regional.hpp"
#include "alod/checkpoint.hpp"
#include <chrono>

// Opt-in phase probe. Use the same source against before/after core libraries.
int main(int argc,char** argv){try{
    using namespace alod;
    using Clock=std::chrono::steady_clock;
    auto input=fixed::parse(argc,argv);
    auto H=fixed::coarse_mesh(input);
    auto h=lod2d::refine_mesh_nvb(H,input.gap);
    LodLimits limits;limits.threads=input.threads;limits.maximum_reference_nodes=input.cap;
    limits.maximum_patch_entries=64000000000ULL;limits.maximum_dense_entries=1000000000ULL;
    auto start=Clock::now();
    auto lap=[&](){auto now=Clock::now();double t=std::chrono::duration<double>(now-start).count();start=now;return t;};
    LodSpace space(H,h,input.wavenumber,input.ell,
        input.policy=="area"?InterpolationPolicy::ManuscriptAreaWeighted:InterpolationPolicy::ArchivedArithmetic,limits);
    const double lod=lap();
    AdditiveKernelRieszContext riesz(space);const double riesz_time=lap();
    std::vector<Problem> problems;for(const auto& member:input.members)problems.push_back(member.problem);
    auto loads=assemble_load_batch(h.mesh,problems,paper_quadrature(input.problem),input.threads);
    const double load_time=lap();
    AdjointTestCache aot;
    auto trained=train_regional(space,riesz,aot,loads,{.9,24,false});
    const double train_time=lap();
    RegionalEvaluator frozen(space,aot);const double frozen_setup=lap();
    auto cold=frozen.evaluate_frozen(loads,trained.phi);const double cold_time=lap();
    auto hot=frozen.evaluate_frozen(loads,trained.phi);const double hot_time=lap();
    const double difference=(hot.values-trained.accepted.values).norm()/std::max(1.,trained.accepted.values.norm());
    if(difference>1e-10||(cold.values-hot.values).norm()>1e-12*std::max(1.,hot.values.norm()))
        throw std::runtime_error("frozen evaluation differs from training");
    std::cout<<std::setprecision(17)<<"{\"lod_seconds\":"<<lod<<",\"riesz_seconds\":"<<riesz_time
        <<",\"load_seconds\":"<<load_time<<",\"training_seconds\":"<<train_time
        <<",\"frozen_setup_seconds\":"<<frozen_setup<<",\"frozen_cold_seconds\":"<<cold_time
        <<",\"frozen_hot_seconds\":"<<hot_time<<",\"reference_nodes\":"<<loads.rows()
        <<",\"coarse_columns\":"<<space.trial().cols()<<",\"rank\":"<<trained.phi.cols()
        <<",\"working_rank\":"<<trained.working_rank<<",\"evaluations\":"<<trained.evaluations
        <<",\"compression_trials\":"<<trained.compression_trials<<",\"training_stop\":\""<<trained.stop
        <<"\",\"dictionary_hash\":\""<<matrix_hash(trained.phi)<<"\",\"solution_hash\":\""<<matrix_hash(trained.accepted.values)
        <<"\",\"kernel_hash\":\""<<matrix_hash(trained.raw_kernel)<<"\",\"eta\":";
    fixed::vector(trained.accepted.eta);
    std::cout<<",\"pg_residual\":"<<trained.accepted.pg_residual<<",\"dense_cache_bytes\":"<<frozen.dense_cache_bytes()<<"}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
