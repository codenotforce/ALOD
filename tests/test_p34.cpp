#include "alod/adaptive.hpp"
#include "alod/regional.hpp"
#include "alod/problems.hpp"
#include <iostream>
#include <stdexcept>
using namespace alod;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
int main(){try{
    for(int m:{1,2,3}){AdaptiveCursor c;for(int s=1;s<=3*m;++s){c.advance(m);require(c.state_id==s&&c.coarse_cycle==(s-1)/m+1&&c.reference_sweep==(s-1)%m+1,"cursor progression");require(c.cycle_complete(m)==(s%m==0),"partial cycle");}}
    EllPolicy p;p.maximum=4;p.last_check=0;
    require(!p.due(3,2,false)&&p.due(4,2,false)&&p.due(1,2,true),"lazy cadence");
    p.extra_checks={2};require(p.due(2,2,false),"extra check");
    require(p.decision(.13,.1,2)=="promote","E1 ratio must take priority over absolute theta");
    require(p.decision(.12,.1,2)=="absolute_accept","E1 strict ratio threshold");
    p.problem="E2";require(p.decision(.011,.1,2)=="promote"&&p.decision(.009,.1,2)=="keep","E2 ratio threshold");
    require(p.decision(1,0,4)=="ell_cap_reached"&&p.decision(0,0,2)=="keep","zero denominator");
    p.mode=EllMode::Fixed;require(p.decision(1,0,2)=="fixed","fixed ell observation");
    auto problem=make_problem("E2");auto H=lod2d::refine_mesh_nvb(problem.initial_mesh,5).mesh;
    auto h=lod2d::refine_mesh_nvb(H,2);LodSpace s(H,h,16,1);AdditiveKernelRieszContext r(s);AdjointTestCache cache;
    ComplexMatrix loads(h.mesh.nodes.size(),2);
    loads.col(0)=lod2d::helmholtz::assemble_helmholtz_load(h.mesh,problem.source,paper_quadrature("E2"),problem.quadrature_context);
    auto second=lod2d::helmholtz::benchmarks::make_parameterized_boundary_gaussian_s_paper_case(16,.7,.6,.4,80,{-.5,.5});
    loads.col(1)=lod2d::helmholtz::assemble_helmholtz_load(h.mesh,second.source,paper_quadrature("E2"),second.quadrature_context);
    auto empty=train_regional(s,r,cache,loads,{0,4,true});require(empty.stop=="empty_region"&&empty.phi.cols()==0&&empty.selected_patches==0,"empty region");
    auto mask=r.regional_mask(.9),all=r.regional_mask(10);int count=0;for(int v:mask)count+=v;
    require(count>0&&count<static_cast<int>(mask.size()),"partial patch mask");
    auto seed=r.apply_selected(loads,mask,true);
    require(seed.selected_identity_relative_error<1e-8&&seed.identity_relative_error<1e-8,"additive AS identities");
    auto full_seed=r.apply_selected(loads,all,true);require((full_seed.values-full_seed.selected_values).norm()<1e-10,"full mask equals full sum");
    auto result=train_regional(s,r,cache,loads,{.9,6,true});
    require(result.phi.cols()>0&&result.working_rank<=6,"nonempty training and rank bound");
    require(result.kernel_residual<1e-8&&result.gram_residual<1e-8&&result.base_orthogonality<1e-8,"dictionary invariants");
    require(result.accepted.pg_residual<1e-9&&result.accepted.raw_base_block<1e-8&&result.accepted.raw_dictionary_block<1e-8,"AOT and coupled block identities");
    ComplexMatrix trial(s.trial().rows(),s.trial().cols()+result.phi.cols());trial<<ComplexMatrix(s.trial()),result.phi;
    ComplexMatrix test(trial.rows(),trial.cols());test<<ComplexMatrix(s.test()),result.accepted.tests;
    ComplexMatrix a=test.adjoint()*s.operators().system*trial;
    ComplexMatrix direct=trial*a.fullPivLu().solve(test.adjoint()*loads);
    require((direct-result.accepted.values).norm()/std::max(1.,direct.norm())<1e-9,"actual coupled PG solution");
    if(result.compression_accepted)for(int j=0;j<loads.cols();++j){ComplexVector delta=result.full_values.col(j)-result.accepted.values.col(j);
        require(std::sqrt(std::max(0.,delta.dot(s.energy().cast<Complex>()*delta).real()))<=.2*result.targets[j]*(1+1e-8),"POD actual solution distance");
        require(result.accepted.eta[j]<=std::max(1.1*result.targets[j],1.03*result.full_eta[j])*(1+1e-8),"POD regional residual");}
    auto factors=cache.factorizations();auto inherited=train_regional(s,r,cache,loads,{.9,6,true},result.raw_kernel,false);
    require(inherited.inherited_rank>0&&cache.factorizations()==factors,"same operator reuse");
    auto reset=train_regional(s,r,cache,loads,{.9,6,false},result.raw_kernel,false);require(reset.inherited_rank==0,"reset arm");
    require((reset.accepted.values-result.accepted.values).norm()<1e-9,"reset ignores inherited representatives");
    auto capped=train_regional(s,r,cache,loads,{.9,1,false});
    require(capped.working_rank<=1&&(!capped.training_met?capped.stop=="rank_cap":capped.stop=="budget_met"),"rank stop remains separate from budget success");
    ComplexMatrix transferred=result.raw_kernel;factors=cache.factorizations();
    for(int ell:{2,3,4}){
        LodSpace promoted(H,h,16,ell);AdditiveKernelRieszContext promoted_r(promoted);
        auto trained=train_regional(promoted,promoted_r,cache,loads,{.9,6,true},transferred,false);
        require(trained.inherited_rank>0&&trained.kernel_residual<1e-8,"ell-only representative inheritance");
        require(cache.factorizations()==factors,"ell-only promotion reuses the reference AOT factor");
        transferred=trained.raw_kernel;
    }
    auto next=refine_pair(H,h,{0,2},{0,4},2,2,20000);
    require((next.reference.P_node*Eigen::VectorXd::Ones(next.coarse.nodes.size())-Eigen::VectorXd::Ones(next.reference.mesh.nodes.size())).norm()<1e-10,"nested constant preservation");
    require((next.injection*Eigen::VectorXd::Ones(h.mesh.nodes.size())-Eigen::VectorXd::Ones(next.reference.mesh.nodes.size())).norm()<1e-10,"reference injection");
    LodSpace changed(next.coarse,next.reference,16,1);AdditiveKernelRieszContext changed_r(changed);
    ComplexMatrix changed_loads=next.injection.cast<Complex>()*loads;
    auto repaired=train_regional(changed,changed_r,cache,changed_loads,{.9,6,true},next.injection.cast<Complex>()*result.raw_kernel,true);
    require(repaired.kernel_residual<1e-8&&cache.factorizations()>factors,"mesh-change kernel repair and AOT invalidation");
    AdjointTestCache probe;auto ops=s.operators();probe.solve(ops,loads);factors=probe.factorizations();probe.solve(ops,2*loads);
    require(probe.factorizations()==factors,"RHS-independent AOT factor");
    ops.diffusion[0]*=1.01;probe.solve(ops,loads);require(probe.factorizations()==++factors,"coefficient invalidation");
    ops.system.coeffRef(0,ops.system.cols()-1)+=1e-9;probe.solve(ops,loads);require(probe.factorizations()==++factors,"pattern invalidation");
    ops.system.coeffRef(1,1)+=1e-9;probe.solve(ops,loads);require(probe.factorizations()==++factors,"numeric invalidation");
    int new_boundary=0;while(std::find(ops.dirichlet_nodes.begin(),ops.dirichlet_nodes.end(),new_boundary)!=ops.dirichlet_nodes.end())++new_boundary;
    ops.dirichlet_nodes.push_back(new_boundary);probe.solve(ops,loads);require(probe.factorizations()==++factors,"free DOF invalidation");
    std::cout<<"P3/P4 scheduling, nested transfer, regional AS, POD, inheritance, coupled PG and AOT cache passed\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
