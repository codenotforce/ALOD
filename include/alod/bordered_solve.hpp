#pragma once
#include "alod/regional.hpp"
#include <Eigen/LU>
#include <vector>
#include <stdexcept>
namespace alod {
// Block elimination for a growing, nonsymmetric complex Schur matrix.
// No explicit inverse; changed bases and unsafe updates rebuild the full LU.
class BorderedSolve {
    struct Border { ComplexMatrix response, lower; Eigen::FullPivLU<ComplexMatrix> lu; };
    ComplexMatrix matrix_;
    Eigen::FullPivLU<ComplexMatrix> base_;
    std::vector<Border> borders_;
    int rebuilds_=0, updates_=0;
    ComplexMatrix apply(const ComplexMatrix& rhs, std::size_t count) const {
        if(!count)return base_.solve(rhs);
        const auto& b=borders_[count-1];
        const int n=b.lower.cols(), s=b.lower.rows();
        ComplexMatrix x=apply(rhs.topRows(n),count-1);
        ComplexMatrix y=b.lu.solve(rhs.bottomRows(s)-b.lower*x);
        ComplexMatrix out(rhs.rows(),rhs.cols());
        out.topRows(n)=x-b.response*y;out.bottomRows(s)=y;return out;
    }
    void rebuild(){
        borders_.clear();base_.compute(matrix_);++rebuilds_;
        if(!base_.isInvertible())throw std::runtime_error("regional coupled Schur block is singular");
    }
public:
    void compute(const ComplexMatrix& next,bool allow_append=true){
        const int n=matrix_.rows(), extra=next.rows()-n;
        // Exact matching prevents an update across a changed POD basis.
        if(allow_append && n>0 && extra>0
           && (next.topLeftCorner(n,n).array()==matrix_.array()).all()){
            Border b;b.lower=next.bottomLeftCorner(extra,n);
            b.response=solve(next.topRightCorner(n,extra));
            b.lu.compute(next.bottomRightCorner(extra,extra)-b.lower*b.response);
            if(b.lu.isInvertible() && b.lu.rcond()>1e-12 && b.response.allFinite()){
                borders_.push_back(std::move(b));matrix_=next;++updates_;return;
            }
        }
        matrix_=next;rebuild();
    }
    ComplexMatrix solve(const ComplexMatrix& rhs){
        ComplexMatrix x=apply(rhs,borders_.size());
        auto good=[&]{return x.allFinite() && (matrix_*x-rhs).norm()<=1e-11*std::max(1e-30,rhs.norm());};
        if(!good()&&!borders_.empty()){rebuild();x=apply(rhs,0);}
        if(!good())throw std::runtime_error("regional Schur solve residual gate failed");
        return x;
    }
    std::size_t entries()const{
        std::size_t n=matrix_.size()+base_.matrixLU().size();
        for(const auto& b:borders_)n+=b.response.size()+b.lower.size()+b.lu.matrixLU().size();
        return n;
    }
    int updates()const{return updates_;}
    int rebuilds()const{return rebuilds_;}
};
}
