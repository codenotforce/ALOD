#include "alod/bordered_solve.hpp"
#include <iostream>
using namespace alod;
int main(){try{
    ComplexMatrix a=ComplexMatrix::Random(24,24);a+=8.*ComplexMatrix::Identity(24,24);
    BorderedSolve solver;
    for(int n=2;n<=24;n+=2){
        ComplexMatrix m=a.topLeftCorner(n,n), rhs=ComplexMatrix::Random(n,16);
        solver.compute(m);auto x=solver.solve(rhs);
        if((x-m.fullPivLu().solve(rhs)).norm()>1e-11)throw std::runtime_error("bordered family mismatch");
    }
    if(solver.updates()!=11||solver.rebuilds()!=1)throw std::runtime_error("borders were not reused");
    a(0,0)+=1.;solver.compute(a);
    if(solver.rebuilds()!=2)throw std::runtime_error("changed basis did not invalidate");
    ComplexMatrix rhs=ComplexMatrix::Random(24,8);
    if((a*solver.solve(rhs)-rhs).norm()>1e-11)throw std::runtime_error("changed basis solve failed");
    BorderedSolve unsafe;ComplexMatrix b=ComplexMatrix::Identity(1,1);unsafe.compute(b);
    ComplexMatrix c=ComplexMatrix::Identity(3,3);c(1,1)=1e-13;
    unsafe.compute(c);
    if(unsafe.rebuilds()!=2||unsafe.updates()!=0)throw std::runtime_error("unsafe border not rebuilt");
    ComplexMatrix crhs=ComplexMatrix::Random(3,4);
    if((c*unsafe.solve(crhs)-crhs).norm()>1e-11)throw std::runtime_error("fallback solve failed");
    std::cout<<"bordered multi-RHS solve verified\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
