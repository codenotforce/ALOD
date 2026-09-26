#include "alod/regional.hpp"
#include "alod/adaptive.hpp"
#include "alod/problems.hpp"
#include <iostream>
using namespace alod;
void require(bool v,const char* s){if(!v)throw std::runtime_error(s);}
int main(){try{
    auto p=make_problem("E2");auto H=lod2d::refine_mesh_nvb(p.initial_mesh,5).mesh;
    auto h=lod2d::refine_mesh_nvb(H,2);LodSpace s(H,h,16,1);
    AdditiveKernelRieszContext r(s);AdjointTestCache cache;
    ComplexMatrix loads(h.mesh.nodes.size(),2);
    loads.col(0)=lod2d::helmholtz::assemble_helmholtz_load(h.mesh,p.source,paper_quadrature("E2"),p.quadrature_context);
    for(int i=0;i<loads.rows();++i)loads(i,1)=Complex(std::sin(i+1.),std::cos(i+2.));
    auto a=train_regional(s,r,cache,loads,{.9,6,true});
    require(a.phi.cols()>0,"fixture must enrich");
    require(cache.factorizations()==0&&cache.solved_columns()==0,"kernel lifting must not solve global adjoint systems");
    require((a.accepted.tests-a.raw_kernel).norm()<1e-8*std::max(1.,a.raw_kernel.norm()),"J phi must recover supported kernel representatives after POD");
    require(a.accepted.kernel_lift_residual<1e-10,"lift must belong to ker I_H");
    ComplexMatrix X(s.trial().rows(),s.trial().cols()+a.phi.cols());X<<ComplexMatrix(s.trial()),a.phi;
    ComplexMatrix Y(X.rows(),X.cols());Y<<ComplexMatrix(s.test()),a.raw_kernel;
    ComplexMatrix block=Y.adjoint()*s.operators().system*X;
    ComplexMatrix direct=X*block.fullPivLu().solve(Y.adjoint()*loads);
    require((direct-a.accepted.values).norm()<1e-8*std::max(1.,direct.norm()),"coupled PG must equal independent kernel-test block solve");
    RegionalEvaluator reused(s,r,cache),cold(s,r,cache,false);
    RegionalEvaluator frozen(s,cache);
    auto frozen_result=frozen.evaluate_frozen(loads,a.phi);
    require((frozen_result.values-a.accepted.values).norm()<1e-10*std::max(1.,a.accepted.values.norm()),"frozen audit needs no estimator factors");
    auto mask=r.regional_mask(.9);
    auto b=reused.evaluate(loads,a.phi,mask);
    auto check_blocks=[&](const RegionalEvaluation& v,const ComplexMatrix& phi){
        double base=(v.raw_tests.adjoint()*s.operators().system*s.trial()).norm()/std::max(1.,s.trial().norm());
        double dictionary=(v.raw_tests.adjoint()*s.operators().system*phi-ComplexMatrix::Identity(phi.cols(),phi.cols())).norm();
        require(std::abs(v.raw_base_block-base)<1e-13*std::max(1.,base),"cached lower-left diagnostic");
        require(std::abs(v.raw_dictionary_block-dictionary)<1e-13*std::max(1.,dictionary),"cached dictionary diagnostic");
    };
    check_blocks(b,a.phi);
    auto scaled=reused.evaluate((2.*loads).eval(),a.phi,mask);
    require((scaled.values-2.*b.values).norm()<1e-10*std::max(1.,b.values.norm()),"changed RHS with prepared dictionary");
    check_blocks(scaled,a.phi);
    ComplexMatrix mixed=a.phi; mixed.col(0)*=Complex(.7,.4);
    auto c=reused.evaluate(loads,mixed,mask),d=cold.evaluate(loads,mixed,mask);
    check_blocks(c,mixed);
    require((b.values-c.values).norm()<1e-8*std::max(1.,b.values.norm()),"test space invariant under basis scaling");
    require((c.values-d.values).norm()<1e-10*std::max(1.,d.values.norm()),"cache invalidation");
    auto inherited=train_regional(s,r,cache,loads,{.9,6,true},a.raw_kernel,false);
    require(inherited.inherited_rank>0&&cache.factorizations()==0,"kernel lift inheritance without AOT");
    auto next=refine_pair(H,h,{0,2},{0,4},2,2,20000);
    LodSpace ns(next.coarse,next.reference,16,2);AdditiveKernelRieszContext nr(ns);
    auto repaired=train_regional(ns,nr,cache,next.injection.cast<Complex>()*loads,{.9,6,true},next.injection.cast<Complex>()*a.raw_kernel,true);
    require(repaired.kernel_residual<1e-8&&cache.factorizations()==0,"mesh repair and promotion");
    require((repaired.accepted.tests-repaired.raw_kernel).norm()<1e-8*std::max(1.,repaired.raw_kernel.norm()),"repaired lift consistency");
    EllPolicy policy;policy.solution_scaled=true;
    require(policy.threshold()==.3,"v4 normalized threshold");
    require(policy.decision(.2,1.,2,2.)=="promote"&&policy.decision(.2,10.,2,20.)=="promote","scale invariant localization");
    require(policy.decision(.1,1.,2,2.)=="absolute_accept","absolute gate follows normalized balance");
    require(policy.decision(.1,0.,4,2.)=="ell_cap_reached","zero estimator cap reporting");
    policy.absolute_threshold=.2;
    require(policy.decision(.25,100.,2,1.)=="promote","absolute defect alone promotes");
    require(policy.decision(.2,100.,2,1.)!="promote","absolute gate is strict");
    require(policy.decision(.25,100.,4,1.)=="ell_cap_reached","absolute defect cap reporting");
    std::cout<<"Kernel lifting, POD, coupled PG, cache, inheritance, repair and normalized policy passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
