#include "alod/localization.hpp"
#include "alod/problems.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace alod;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
int main(){try{
    for(const auto id:{"E1","E2"}) for(int threads:{1,4}) {
        const auto problem=make_problem(id);
        const auto coarse=lod2d::refine_mesh_nvb(problem.initial_mesh,4).mesh;
        // V_h == V_H: every local kernel is zero, independently of the RHS.
        LodLimits limits;limits.threads=threads;
        LodSpace space(coarse,lod2d::refine_mesh_nvb(coarse,0),16,2,
                       InterpolationPolicy::ManuscriptAreaWeighted,limits);
        AdditiveKernelRieszContext context(space);
        int nonempty=0;
        for(const auto& patch:context.patches()) {
            const int n=patch.discrete_dofs.size();
            require(patch.constraints.rows()==n,"fixture must have zero-dimensional kernels");
            if(n>0)++nonempty;
        }
        require(nonempty>0,"fixture must exercise positive-size square constraint matrices");
        ComplexMatrix rhs(space.fine().nodes.size(),3);
        for(int i=0;i<rhs.rows();++i)for(int j=0;j<rhs.cols();++j)
            rhs(i,j)=Complex(std::sin((i+1)*(j+1)),std::cos((i+1)*(j+2)));
        std::vector<int> mask(coarse.nodes.size());
        for(int i=0;i<static_cast<int>(mask.size());++i)mask[i]=i%2;
        for(double scale:{1.0,1e-12}) for(bool full:{false,true}) {
            auto r=context.apply_selected(scale*rhs,mask,full);
            require(r.values.norm()==0 && r.selected_values.norm()==0,
                    "zero-dimensional Riesz map must return exact zero");
            require(r.eta.norm()==0 && r.selected_eta.norm()==0 && r.node_eta_squared.norm()==0,
                    "zero-dimensional estimator contributions must be exact zero");
            require(r.constraint_relative_residual==0 && r.local_relative_residual==0,
                    "vacuous local problem must pass validation");
        }
        std::cout<<id<<" threads="<<threads<<" nonempty_zero_kernels="<<nonempty<<" passed\n";
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
