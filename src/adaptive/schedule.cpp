#include "alod/adaptive.hpp"
#include <cmath>
#include <stdexcept>
namespace alod {
void AdaptiveCursor::advance(int m) {
    if(m<1)throw std::invalid_argument("m_ref must be positive");
    ++state_id;
    if(coarse_cycle==0||reference_sweep==m){++coarse_cycle;reference_sweep=1;}
    else ++reference_sweep;
}
bool AdaptiveCursor::cycle_complete(int m) const {return state_id==0||reference_sweep==m;}
bool EllPolicy::due(int state,int ell,bool terminal) const {
    if(ell<1||ell>maximum||state<0)throw std::invalid_argument("invalid ell schedule state");
    return mode==EllMode::EveryState||state==0||terminal||extra_checks.contains(state)
        ||last_check<0||state-last_check>=(1<<ell);
}
std::string EllPolicy::decision(double theta,double eta,int ell) const {
    auto ratio=localization_ratio(theta,eta);
    if(mode==EllMode::Fixed)return "fixed";
    bool promote=ratio.status=="zero_denominator"||(ratio.value&&*ratio.value>(problem=="E1"?1.2:.1));
    if(promote)return ell<maximum?"promote":"ell_cap_reached";
    if(problem=="E1"&&theta<.2)return "absolute_accept";
    return "keep";
}
}
