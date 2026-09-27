#include "alod/kernel_defect.hpp"
#include "alod/execution.hpp"
#include "alod/problems.hpp"
#include <iostream>
using namespace alod;
int main(){try{
    for(const auto* name:{"E1","E2"})for(int gap:{0,2}){
        auto coarse=lod2d::refine_mesh_nvb(make_problem(name).initial_mesh,3).mesh;
        LodLimits limits;limits.threads=4;
        LodSpace space(coarse,lod2d::refine_mesh_nvb(coarse,gap),16,1,
                       InterpolationPolicy::ManuscriptAreaWeighted,limits);
        AdditiveKernelRieszContext context(space);
        limits.parallel_riesz_gather=false;
        LodSpace serial_space(coarse,lod2d::refine_mesh_nvb(coarse,gap),16,1,
                              InterpolationPolicy::ManuscriptAreaWeighted,limits);
        AdditiveKernelRieszContext serial(serial_space);
        std::vector<Eigen::Triplet<Complex>> entries;
        for(int i=0;i<static_cast<int>(space.fine().nodes.size());++i)
            for(int j=0;j<7;++j)if((i+3*j)%4==0)
                entries.emplace_back(i,j,Complex(std::sin(i+j),std::cos(2*i+j)));
        ComplexSparseMatrix defect(space.fine().nodes.size(),7);
        defect.setFromTriplets(entries.begin(),entries.end());
        KernelDefectOperator fused(context,defect);
        ComplexMatrix block(7,3);
        for(int i=0;i<7;++i)for(int j=0;j<3;++j)block(i,j)=Complex(std::cos(i+j),std::sin(i-j));
        const ComplexMatrix expected=defect.adjoint()*serial.apply_action(defect*block);
        for(int threads:{1,4}){
            ExecutionScope scope(thread_budget(threads));
            const ComplexMatrix rhs=defect*block;
            std::vector<int> mask(coarse.nodes.size());
            for(std::size_t i=0;i<mask.size();++i)mask[i]=i%2;
            for(bool full:{false,true}){
                const auto a=context.apply_selected(rhs,mask,full);
                const auto b=serial.apply_selected(rhs,mask,full);
                if((a.values-b.values).norm()!=0 ||
                   (a.selected_values-b.selected_values).norm()!=0 ||
                   (a.node_eta_squared-b.node_eta_squared).norm()!=0)
                    throw std::runtime_error("parallel Riesz gather changed ordered reduction");
            }
            auto actual=fused.apply(block);
            const auto global=fused.apply_global(block);
            if((global-expected).norm()>1e-12*std::max(1.,expected.norm()))
                throw std::runtime_error("parallel sparse products changed action");
            if((actual-expected).norm()>1e-12*std::max(1.,expected.norm()))
                throw std::runtime_error("fused patch action differs from global assembly");
            const ComplexMatrix gram=block.adjoint()*actual;
            if((gram-gram.adjoint()).norm()>1e-11*std::max(1.,gram.norm()))
                throw std::runtime_error("fused operator lost Hermitian symmetry");
            if(gram.diagonal().real().minCoeff()<-1e-12)
                throw std::runtime_error("fused operator lost positivity");
        }
    }
    std::cout<<"Fused/global equivalence, zero kernels, Hermitian action and thread parity passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
