// Offline same-checkpoint benchmark. Restoration is timed separately from Theta.
#include "alod/checkpoint.hpp"
#include "alod/localization.hpp"
#include "alod/timing.hpp"
#include "alod/kernel_defect.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sys/resource.h>

using namespace alod;
using Clock=std::chrono::steady_clock;
int main(int argc,char** argv){try{
    if(argc!=6)throw std::invalid_argument(
        "usage: benchmark_theta checkpoint k threads baseline|light|reuse|optimized|parallel|fused|probe cold|warm");
    const double k=std::stod(argv[2]);const int threads=std::stoi(argv[3]);
    const std::string mode=argv[4],temperature=argv[5];
    if(!(k>0)||threads<1||(temperature!="cold"&&temperature!="warm")
       ||(mode!="baseline"&&mode!="light"&&mode!="reuse"&&mode!="optimized"&&mode!="fused"&&mode!="parallel"&&mode!="probe"))
        throw std::invalid_argument("invalid benchmark arguments");
    omp_set_dynamic(0);omp_set_num_threads(threads);PhaseTimer::probe(threads);
    const auto start=Clock::now();
    auto snapshot=load_checkpoint(argv[1]);
    if(snapshot.phase!=CheckpointPhase::Accepted||!snapshot.lod_trial.cols())
        throw std::invalid_argument("benchmark requires an accepted checkpoint with a LOD basis");
    LodLimits limits;limits.threads=threads;limits.maximum_reference_nodes=4000000;
    limits.parallel_riesz_gather=mode=="parallel"||mode=="fused"||mode=="probe";
    limits.maximum_patch_entries=64000000000ULL;limits.maximum_dense_entries=1000000000;
    // This benchmark targets the current manuscript's area interpolation/N2 policy.
    LodSpace space(snapshot.coarse.mesh,snapshot.reference(),k,snapshot.ell,
        InterpolationPolicy::ManuscriptAreaWeighted,limits,snapshot.lod_trial,&snapshot.lod_reduced);
    const int state=snapshot.cursor.state_id,ell=snapshot.ell;
    LocalizationEigenConfig config;config.eigenvalue_relative_residual=true;
    config.reuse_operator_actions=mode=="reuse"||mode=="optimized"||mode=="fused"||mode=="parallel";
    config.lightweight_riesz=mode=="light"||mode=="optimized"||mode=="fused"||mode=="parallel";
    config.fused_defect=mode=="fused";
    config.parallel_defect=mode=="fused"||mode=="parallel";
    if(temperature=="warm"){
        if(snapshot.warm_full.cols()==0)throw std::runtime_error("checkpoint has no warm block");
        ComplexMatrix block(space.coarse_nodes().size(),snapshot.warm_full.cols());
        for(int j=0;j<block.rows();++j)block.row(j)=snapshot.warm_full.row(space.coarse_nodes()[j]);
        config.warm_start=LocalizationWarmStart{space.identity()+":manuscript-n2:adjoint",std::move(block)};
    }
    snapshot=Checkpoint{};
    AdditiveKernelRieszContext riesz(space);
    const double restore=std::chrono::duration<double>(Clock::now()-start).count();
    if(mode=="probe"){
        ComplexSparseMatrix defect=space.operators().system.adjoint()*space.test();
        const ComplexSparseMatrix coarse=space.prolongation().transpose().cast<Complex>()*defect;
        defect-=space.interpolation().transpose().cast<Complex>()*coarse;defect.makeCompressed();
        KernelDefectOperator action(riesz,defect);
        ComplexMatrix block(defect.cols(),4);
        for(int i=0;i<block.rows();++i)for(int j=0;j<block.cols();++j)
            block(i,j)=Complex(std::sin((i+1.)*(j+1.)),std::cos((i+1.)*(j+2.)));
        block/=block.norm();ComplexMatrix expected;
        for(const auto* name:{"global","fused","parallel","parallel","fused","global"}){
            const auto begin=Clock::now();ComplexMatrix result;
            if(std::string(name)=="fused")result=action.apply(block);
            else if(std::string(name)=="parallel")result=action.apply_global(block);
            else result=defect.adjoint()*riesz.apply_action(defect*block);
            const double elapsed=std::chrono::duration<double>(Clock::now()-begin).count();
            if(expected.size()==0)expected=result;
            const double difference=(result-expected).norm()/std::max(1e-30,expected.norm());
            if(difference>1e-9)throw std::runtime_error("defect probe failed numerical equivalence");
            std::cout<<std::setprecision(17)<<"{\"mode\":\""<<name<<"\",\"seconds\":"<<elapsed
                <<",\"relative_difference\":"<<difference<<",\"defect_nonzeros\":"<<defect.nonZeros()
                <<",\"restore_seconds\":"<<restore<<"}\n"<<std::flush;
        }
        return 0;
    }
    const auto begin=Clock::now();auto result=localization_theta(space,riesz,config);
    const double wall=std::chrono::duration<double>(Clock::now()-begin).count();
    struct rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<std::setprecision(17)<<"{\"state\":"<<state<<",\"ell\":"<<ell
        <<",\"mode\":\""<<mode<<"\",\"temperature\":\""<<temperature<<"\",\"threads\":"<<threads
        <<",\"theta\":"<<result.theta<<",\"residual\":"<<result.spectrum.relative_residual
        <<",\"iterations\":"<<result.spectrum.iterations
        <<",\"operator_applications\":"<<result.spectrum.operator_applications
        <<",\"operator_columns\":"<<result.spectrum.operator_columns
        <<",\"operator_seconds\":"<<result.spectrum.operator_seconds
        <<",\"restore_seconds\":"<<restore<<",\"theta_seconds\":"<<wall
        <<",\"process_peak_rss_kib\":"<<usage.ru_maxrss<<"}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
