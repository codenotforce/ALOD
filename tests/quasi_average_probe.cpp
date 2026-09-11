#include "lod/quasi_interp.h"
#include "mesh/refine.h"
#include "helmholtz/boundary.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>

// Compare the archived arithmetic operator with the manuscript's area average.
// Every column is a continuous fine P1 nodal function, not arbitrary DG data.
double probe(double depth) {
    using namespace lod2d;
    using Sparse = Eigen::SparseMatrix<double>;
    using Triplet = Eigen::Triplet<double>;
    TriMesh coarse;
    coarse.nodes={{0,0},{1,0},{0,1},{0,-depth}};
    coarse.elems={{0,1,2},{0,3,1}};
    helmholtz::tag_all_boundary_edges(coarse,BoundaryTag::Robin);
    auto f=refine_mesh_nvb(coarse,2);
    int ne=f.mesh.elems.size(), nc=coarse.elems.size();
    std::vector<Triplet> cg, mass, inv, average;
    auto fa=compute_area(f.mesh), ca=compute_area(coarse);
    Eigen::VectorXd sums=Eigen::VectorXd::Zero(coarse.nodes.size());
    for(int e=0;e<nc;++e) for(int i=0;i<3;++i) sums[coarse.elems[e][i]]+=ca[e];
    for(int e=0;e<ne;++e) for(int i=0;i<3;++i) {
        cg.emplace_back(3*e+i,f.mesh.elems[e][i],1);
        for(int j=0;j<3;++j) mass.emplace_back(3*e+i,3*e+j,fa[e]/12*(i==j?2:1));
    }
    for(int e=0;e<nc;++e) for(int i=0;i<3;++i) {
        average.emplace_back(coarse.elems[e][i],3*e+i,ca[e]/sums[coarse.elems[e][i]]);
        for(int j=0;j<3;++j) inv.emplace_back(3*e+i,3*e+j,(i==j?9:-3)/ca[e]);
    }
    Sparse C(3*ne,f.mesh.nodes.size()),M(3*ne,3*ne),B(3*nc,3*nc),E(coarse.nodes.size(),3*nc);
    C.setFromTriplets(cg.begin(),cg.end()); M.setFromTriplets(mass.begin(),mass.end());
    B.setFromTriplets(inv.begin(),inv.end()); E.setFromTriplets(average.begin(),average.end());
    Sparse expected=E*B*f.P_dg.transpose()*M*C;
    Sparse manuscript=build_quasi_interp(coarse,f.mesh,f.P_dg,C,f.mesh.nodes.size(),coarse.nodes.size(),QuasiInterpolationPolicy::ManuscriptAreaWeighted);
    if((manuscript-expected).norm()>1e-12)throw std::runtime_error("manuscript interpolation assembly mismatch");
    Sparse actual=build_quasi_interp(coarse,f.mesh,f.P_dg,C,f.mesh.nodes.size(),coarse.nodes.size());
    Sparse identity(coarse.nodes.size(),coarse.nodes.size());identity.setIdentity();
    Sparse check=actual*f.P_node-identity, check2=expected*f.P_node-identity;
    if(check.norm()>1e-12 || check2.norm()>1e-12) throw std::runtime_error("projection identity failed");
    return Sparse(actual-expected).norm();
}
int main() {
    double uniform=probe(1), graded=probe(2);
    std::cout<<std::setprecision(17)<<"{\"uniform_difference\":"<<uniform
             <<",\"graded_difference\":"<<graded<<",\"both_reproduce_coarse_P1\":true}\n";
    return uniform<1e-12 && graded>1e-3 ? 0 : 1;
}
