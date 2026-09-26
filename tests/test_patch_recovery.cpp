#include "helmholtz/patch_solver.h"
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/LU>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <omp.h>
using namespace lod2d::helmholtz;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
double residual(const HelmholtzPatchSolveResult& x){return std::max({x.diagnostics.primal_residual,x.diagnostics.adjoint_residual,x.diagnostics.constraint_residual});}
void check(double eps,bool complex){
    HelmholtzPatchSystem s;ComplexMatrix A(2,2);A<<.5+eps,.5,.5,.5+eps;
    if(complex)A(1,1)+=Complex(0,.2);
    s.helmholtz=A.sparseView();s.constraints.resize(1,2);s.constraints<<1,0;
    s.rhs.resize(2,3);s.rhs<<1.,Complex(.2,.3),-.7,.7,1.,Complex(.1,.8);
    auto x=solve_helmholtz_patch(s);
    ComplexMatrix expected=ComplexMatrix::Zero(2,3);expected.row(1)=s.rhs.row(1)/A(1,1);
    require((x.corrector-expected).norm()<1e-10,"constrained solution differs from analytic kernel solve");
    require(residual(x)<1e-10,"residual unchanged gate failed");
    if(eps==0&&!complex)require(x.diagnostics.saddle_fallback,"singular unconstrained block did not use saddle solve");
}
int main(int argc,char**argv){try{
    std::cout<<std::setprecision(17);
    for(double e:{1.,1e-8,1e-12,0.})for(bool c:{false,true})check(e,c);
    // Full-rank constraints imply an exactly zero corrector; indefinite A is allowed.
    HelmholtzPatchSystem z;z.helmholtz.resize(2,2);z.helmholtz.insert(0,0)=1;z.helmholtz.insert(1,1)=-1;
    z.constraints=Eigen::MatrixXd::Identity(2,2);z.rhs=ComplexMatrix::Ones(2,3);
    require(solve_helmholtz_patch(z).corrector.norm()<1e-12,"zero kernel incorrect");
    z.constraints.resize(0,2);z.helmholtz.setZero();bool rejected=false;
    try{solve_helmholtz_patch(z);}catch(const std::runtime_error&){rejected=true;}
    require(rejected,"truly singular constrained problem accepted");
    int checked=0;
    if(argc==2)for(const auto& f:std::filesystem::directory_iterator(argv[1]))if(f.path().filename().string().ends_with("-A.mtx")) {
        auto base=f.path().string();base.resize(base.size()-6);HelmholtzPatchSystem s;ComplexSparseMatrix C,b;
        require(Eigen::loadMarket(s.helmholtz,base+"-A.mtx"),"read A");
        require(Eigen::loadMarket(C,base+"-C.mtx"),"read C");require(Eigen::loadMarket(b,base+"-rhs.mtx"),"read rhs");
        s.constraints=ComplexMatrix(C).real();s.rhs=ComplexMatrix(b);
        auto x=solve_helmholtz_patch(s);require(x.diagnostics.saddle_fallback,"actual failed patch did not use repair");require(residual(x)<1e-10,"actual patch recovery residual");
        ComplexMatrix values[2];
        #pragma omp parallel for num_threads(32)
        for(int j=0;j<2;++j)values[j]=solve_helmholtz_patch(s).corrector;
        require((values[0]-x.corrector).norm()<1e-12&&(values[1]-x.corrector).norm()<1e-12,"parallel patch recovery differs");
        std::cout<<base<<" original="<<x.diagnostics.original_residual<<" repaired="<<residual(x)<<" n="<<s.helmholtz.rows()<<" m="<<s.constraints.rows()<<std::endl;++checked;
    }
    if(argc==2)require(checked>0,"no actual failed patches checked");
    std::cout<<"PASS analytic cases, singular gate, parallel recovery; actual patches="<<checked<<std::endl;
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
