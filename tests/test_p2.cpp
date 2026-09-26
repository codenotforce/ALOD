#include "alod/localization.hpp"
#include "alod/problems.hpp"
#include "helmholtz/boundary.h"
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <cmath>
#include <set>

using namespace alod;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F> void rejects(F fn){bool rejected=false;try{fn();}catch(const std::exception&){rejected=true;}require(rejected,"invalid input was accepted");}
void close(double a,double b,double tolerance=1e-9){if(std::abs(a-b)>tolerance*std::max({1.0,std::abs(a),std::abs(b)}))throw std::runtime_error("numerical comparison failed: "+std::to_string(a)+" vs "+std::to_string(b));}
void run(const std::string& id,bool graded) {
    auto problem=make_problem(id);auto coarse=lod2d::refine_mesh_nvb(problem.initial_mesh,2).mesh;
    if(graded)coarse=lod2d::bisect_newest_vertex(coarse,{0}).mesh;
    LodSpace space(coarse,lod2d::refine_mesh_nvb(coarse,3),16,1);
    AdditiveKernelRieszContext riesz(space);
    const int n=space.fine().nodes.size();
    ComplexMatrix residual(n,3);
    for(int i=0;i<n;++i)for(int j=0;j<3;++j)residual(i,j)=Complex(std::sin((i+1)*(j+1)),std::cos((i+1)*(j+2)));
    for(int node:space.operators().dirichlet_nodes)residual.row(node).setZero();
    auto actual=riesz.apply(residual);
    ComplexMatrix map=ComplexMatrix::Zero(n,n);
    Eigen::MatrixXd local_mass=Eigen::MatrixXd::Zero(coarse.nodes.size(),3);
    require(riesz.patches().size()==coarse.nodes.size(),"geometric boundary vertices omitted");
    for(const auto& patch:riesz.patches()) {
        require(patch.coarse_node>=0,"missing patch center");
        int size=patch.discrete_dofs.size();if(size==0)continue;
        Eigen::MatrixXd kernel;
        if(patch.constraints.rows())kernel=patch.constraints.fullPivLu().kernel();
        else kernel=Eigen::MatrixXd::Identity(size,size);
        Eigen::MatrixXd energy(size,size);
        for(int i=0;i<size;++i)for(int j=0;j<size;++j)energy(i,j)=space.energy().coeff(patch.discrete_dofs[i],patch.discrete_dofs[j]);
        Eigen::MatrixXd local=Eigen::MatrixXd::Zero(size,size);
        if(kernel.norm()>1e-12) {
            Eigen::MatrixXd gram=kernel.transpose()*energy*kernel;
            local=kernel*gram.ldlt().solve(kernel.transpose());
        }
        ComplexMatrix embedded=ComplexMatrix::Zero(n,n);
        for(int i=0;i<size;++i)for(int j=0;j<size;++j)embedded(patch.discrete_dofs[i],patch.discrete_dofs[j])=local(i,j);
        require((space.interpolation().cast<Complex>()*embedded).norm()<1e-9,"raw representative left the kernel");
        // Every active fine hat has its complete support in N^2(z).
        std::set<int> allowed(patch.coarse_elements.begin(),patch.coarse_elements.end());
        for(int e=0;e<static_cast<int>(space.fine().elems.size());++e) {
            bool touches=false;for(int node:space.fine().elems[e])
                if(std::binary_search(patch.discrete_dofs.begin(),patch.discrete_dofs.end(),node))touches=true;
            if(touches)for(int parent=0;parent<space.element_prolongation().cols();++parent)
                if(space.element_prolongation().coeff(e,parent)!=0)require(allowed.contains(parent),"patch support escaped N^2");
        }
        map+=embedded;
        ComplexMatrix values=embedded*residual;
        for(int j=0;j<3;++j)local_mass(patch.coarse_node,j)=std::real(residual.col(j).dot(values.col(j)));
    }
    require((actual.values-map*residual).norm()<1e-8,"batched Schur differs from independent dense kernel inverse");
    require((actual.node_eta_squared-local_mass).norm()<1e-8,"local patch energies disagree");
    require(actual.identity_relative_error<1e-9&&actual.constraint_relative_residual<1e-9&&actual.local_relative_residual<1e-9,"Riesz identity/residual failed");
    std::vector<int> mask(coarse.nodes.size(),0);
    require(riesz.apply_selected(residual,mask).values.norm()==0,"empty selected mask");
    for(int i=0;i<static_cast<int>(mask.size());++i)mask[i]=i%2;
    auto partial=riesz.apply_selected(residual,mask,true),selected=riesz.apply_selected(residual,mask);
    require((selected.values-partial.selected_values).norm()<1e-10,"selected multiplicity inconsistent");
    require((partial.values-actual.values).norm()<1e-10,"selected mask changed full estimator");
    require(selected.identity_relative_error<1e-9,"selected identity failed");
    for(int j=0;j<3;++j){double sum=0;for(int i=0;i<static_cast<int>(mask.size());++i)if(mask[i])sum+=actual.node_eta_squared(i,j);close(sum,selected.eta(j)*selected.eta(j));}
    std::fill(mask.begin(),mask.end(),1);require((riesz.apply_selected(residual,mask).values-actual.values).norm()<1e-10,"full mask multiplicity");
    require(riesz.factorizations()<riesz.patches().size(),"fixture did not exercise duplicate patch factors");
    for(int j=0;j<3;++j)require((riesz.apply(residual.col(j)).values-actual.values.col(j)).norm()<1e-10,"batch/scalar disagreement");
    auto zero=riesz.apply(ComplexMatrix::Zero(n,2));require(zero.eta.norm()==0&&zero.values.norm()==0,"zero residual");
    auto region=riesz.regional_mask(.8);
    for(const auto& patch:riesz.patches()) {
        bool inside=!patch.discrete_dofs.empty()&&!patch.coarse_elements.empty();
        for(int e:patch.coarse_elements)for(int node:coarse.elems[e])if(coarse.nodes[node].norm()>.8+1e-13)inside=false;
        require(region[patch.coarse_node]==int(inside),"region must contain complete patch");
    }
    rejects([&]{riesz.regional_mask(-1);});rejects([&]{riesz.apply_selected(residual,{});});mask[0]=2;rejects([&]{riesz.apply_selected(residual,mask);});
    auto quad=paper_quadrature(id);
    ComplexVector load=lod2d::helmholtz::assemble_helmholtz_load(space.fine(),problem.source,quad,problem.quadrature_context);
    ComplexMatrix loads(n,2);loads.col(0)=load;loads.col(1)=Complex(2,-1)*load;
    auto solution=space.solve(loads);require((solution.values.col(1)-Complex(2,-1)*solution.values.col(0)).norm()<1e-10,"PG block linearity");
    auto estimate=riesz.estimate(loads,solution.values,.15);
    for(int j=0;j<2;++j){close(estimate.node_eta_squared.col(j).sum(),estimate.element_eta_squared.col(j).sum());close(estimate.eta(j)*estimate.eta(j),estimate.element_eta_squared.col(j).sum());}
    for(int e=0;e<static_cast<int>(coarse.elems.size());++e)for(int j=0;j<2;++j) {
        double allocated=0;
        for(int z:coarse.elems[e])allocated+=estimate.node_eta_squared(z,j)/riesz.patches()[z].coarse_hat_support.size();
        close(allocated,estimate.element_eta_squared(e,j));
    }
    auto strong=reference_strong_residual(space,solution.values.col(0),load,problem.source,quad,problem.quadrature_context);
    require(strong.fine.element_squared.size()==space.fine().elems.size(),"reference indicator is on the wrong mesh");
    require(strong.fine.algebraic_relative_difference<1e-8,"strong residual reconstruction");
    require(std::abs(strong.fine.eta-estimate.eta(0))>1e-5,"strong and kernel estimator were mixed");
    LocalizationEigenConfig exact;exact.relative_tolerance=1e-10;
    auto theta=localization_theta(space,riesz,exact);
    ComplexMatrix defect=space.operators().system.adjoint()*space.test();
    ComplexMatrix gram=defect.adjoint()*map*defect;
    gram=(.5*(gram+gram.adjoint())).eval();
    ComplexMatrix denominator=space.coarse_basis().transpose()*space.energy()*space.coarse_basis();
    Eigen::GeneralizedSelfAdjointEigenSolver<ComplexMatrix> spectrum(gram,denominator);
    close(theta.theta*theta.theta,spectrum.eigenvalues().tail(1)(0),1e-8);
    close(theta.spectrum.lambda_max,localization_theta(space,riesz,exact,false).spectrum.lambda_max);
    LocalizationEigenConfig iterative;iterative.dense_cross_check_max_dimension=0;iterative.relative_tolerance=1e-8;
    auto cold=localization_theta(space,riesz,iterative);close(cold.theta,theta.theta,1e-7);
    iterative.warm_start=cold.warm_start;auto warm=localization_theta(space,riesz,iterative);
    require(warm.spectrum.used_warm_start,"block warm start not used");close(warm.theta,theta.theta,1e-7);
    auto dense_warm=exact;dense_warm.warm_start=cold.warm_start;
    close(localization_theta(space,riesz,dense_warm).spectrum.lambda_max,theta.spectrum.lambda_max,1e-7);
    iterative.warm_start->identity+="stale";rejects([&]{localization_theta(space,riesz,iterative);});
    iterative.warm_start.reset();iterative.maximum_iterations=0;iterative.relative_tolerance=1e-30;
    rejects([&]{localization_theta(space,riesz,iterative);});
    if(theta.theta>0)close(*localization_ratio(theta.theta,estimate.eta(1)).value*std::sqrt(5.0),*localization_ratio(theta.theta,estimate.eta(0)).value);
    std::cout<<id<<" graded="<<graded<<" theta="<<theta.theta<<" eta="<<estimate.eta(0)<<" checked\n";
}
void marking() {
    Eigen::MatrixXd supplemental(4,2);supplemental<<10,0, 0,6, 0,5, 0,5;
    auto supplement=mark_family(supplemental,Eigen::VectorXd::Ones(2),{{0,MemberRole::Train},{1,MemberRole::Train}},{0,1},.3);
    require(supplement.marked_elements==std::vector<int>({0,1}),"worst-member supplement missing");
    Eigen::MatrixXd mass(4,3);mass<<9,0,1000, 1,0,1000, 0,0,1000, 0,16,1000;
    Eigen::VectorXd energies(3);energies<<1,1,1;
    std::vector<FamilyMember> members{{10,MemberRole::Train},{20,MemberRole::Train},{0,MemberRole::Test}};
    auto family=mark_family(mass,energies,members,{20,10},.7);
    require(family.worst_member_id==20&&family.mean_bulk>=.7&&family.worst_bulk>=.7,"mean+worst bulk");
    auto nominal=mark_family(mass,energies,members,{10},.7);
    require(nominal.marked_elements==std::vector<int>{0},"explicit nominal member");
    mass.col(2).setConstant(1e100);require(mark_family(mass,energies,members,{10,20},.7).marked_elements==family.marked_elements,"audit leakage");
    mass.col(0).swap(mass.col(1));std::swap(members[0],members[1]);
    require(mark_family(mass,energies,members,{20,10},.7).marked_elements==family.marked_elements,"column order affected IDs");
    rejects([&]{mark_family(mass,energies,members,{0},.7);});rejects([&]{mark_family(mass,energies,members,{99},.7);});rejects([&]{mark_family(mass,energies,members,{10,10},.7);});
    for(auto role:{MemberRole::Shift,MemberRole::Pure}){members[2].role=role;rejects([&]{mark_family(mass,energies,members,{0},.7);});}
    require(!localization_ratio(1,0).value&&localization_ratio(1,0).status=="zero_denominator","zero ratio denominator");
    require(!localization_ratio(0,0).value&&localization_ratio(0,0).status=="zero_defect","zero defect status");
    require(localization_ratio(0,1).value==0,"zero numerator has a defined zero ratio");
}
int main(){try{marking();for(auto id:{"E1","E2"})for(bool graded:{false,true})run(id,graded);}
catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
