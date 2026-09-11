#include "alod/estimator.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

namespace alod {
FamilyMarking mark_family(const Eigen::MatrixXd& mass,const Eigen::VectorXd& energy,
    const std::vector<FamilyMember>& members,const std::vector<int>& ids,double theta) {
    if(mass.cols()!=static_cast<int>(members.size()) || energy.size()!=mass.cols()
        || mass.rows()<1 || ids.empty() || !mass.allFinite() || !energy.allFinite()
        || (mass.array()<0).any() || (energy.array()<0).any() || !std::isfinite(theta) || theta<0 || theta>1)
        throw std::invalid_argument("invalid family marking inputs");
    std::map<int,int> columns;
    for(int i=0;i<static_cast<int>(members.size());++i)
        if(members[i].id<0 || !columns.emplace(members[i].id,i).second)
            throw std::invalid_argument("family member IDs must be unique and nonnegative");
    FamilyMarking result;result.member_ids=ids;
    std::sort(result.member_ids.begin(),result.member_ids.end());
    if(std::adjacent_find(result.member_ids.begin(),result.member_ids.end())!=result.member_ids.end())
        throw std::invalid_argument("duplicate training ID");
    result.frozen_scales.resize(ids.size());result.mean_squared=Eigen::VectorXd::Zero(mass.rows());
    double worst=-1;
    for(int j=0;j<static_cast<int>(ids.size());++j) {
        auto found=columns.find(result.member_ids[j]);
        if(found==columns.end() || members[found->second].role!=MemberRole::Train)
            throw std::invalid_argument("marking requires explicit existing training members");
        int col=found->second;double scale=std::max(energy(col),1e-12);
        result.frozen_scales(j)=scale;
        for(int e=0;e<mass.rows();++e)result.mean_squared(e)+=mass(e,col)/(ids.size()*scale*scale);
        double proxy=std::sqrt(mass.col(col).sum())/scale;
        if(proxy>worst) {
            worst=proxy;result.worst_member_id=result.member_ids[j];
            result.worst_squared=mass.col(col)/(scale*scale);
        }
    }
    if(!result.mean_squared.allFinite() || !result.worst_squared.allFinite())
        throw std::invalid_argument("normalized family mass overflow");
    std::vector<double> aggregate(result.mean_squared.data(),result.mean_squared.data()+mass.rows());
    result.marked_elements=lod2d::helmholtz::adaptive::mark_doerfler(aggregate,theta);
    std::vector<bool> selected(mass.rows(),false);
    double covered=0,total=result.worst_squared.sum();
    for(int e:result.marked_elements){selected[e]=true;covered+=result.worst_squared(e);}
    std::vector<int> order;
    for(int e=0;e<mass.rows();++e)if(!selected[e])order.push_back(e);
    std::stable_sort(order.begin(),order.end(),[&](int a,int b){
        if(result.worst_squared(a)!=result.worst_squared(b))return result.worst_squared(a)>result.worst_squared(b);
        return a<b;
    });
    // Preserve the archived minimal-supplement tolerance and stable index ties.
    for(int e:order) {
        if(covered+1e-14*std::max(1.0,total)>=theta*total)break;
        result.marked_elements.push_back(e);covered+=result.worst_squared(e);
    }
    std::sort(result.marked_elements.begin(),result.marked_elements.end());
    double mean_covered=0;
    for(int e:result.marked_elements)mean_covered+=result.mean_squared(e);
    result.mean_bulk=result.mean_squared.sum()>0?mean_covered/result.mean_squared.sum():1;
    result.worst_bulk=total>0?covered/total:1;
    return result;
}
ReferenceStrongResidual reference_strong_residual(const LodSpace& space,
    const ComplexVector& value,const ComplexVector& load,
    const lod2d::helmholtz::ComplexFunction& source,
    const lod2d::helmholtz::QuadraturePolicy& quadrature,
    const lod2d::helmholtz::QuadratureContext& context) {
    return {lod2d::helmholtz::adaptive::diagnostics::estimate_conforming_p1_residual(
        space.fine(),space.operators(),value,load,source,quadrature,context)};
}
}
