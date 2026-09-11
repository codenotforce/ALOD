#pragma once
#include "alod/problems.hpp"
#include "mesh/refine.h"
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <omp.h>
namespace fixed {
using namespace lod2d;
using namespace lod2d::helmholtz;
struct Member { int id; std::string role; alod::Problem problem; };
struct Input {
    std::string problem,policy,riesz_policy,members_path;
    int level=2,gap=3,ell=1,threads=1,cap=20000;
    bool graded=false;
    double theta=.15,tolerance=1e-4;
    int iterations=750,dense=64;
    std::vector<int> ids;
    std::vector<Member> members;
};
inline Input parse(int argc,char** argv) {
    if(argc<2)throw std::invalid_argument("expected E1 or E2; see --help");
    Input in;in.problem=argv[1];(void)alod::make_problem(in.problem);
    std::map<std::string,std::string> options{{"level","2"},{"gap","3"},{"ell","1"},{"graded","0"},
        {"threads","1"},{"maximum-nodes","20000"},{"interpolation","area"},{"members",""},
        {"training-ids","0"},{"theta","0.15"},{"ritz-tolerance","0.0001"},{"ritz-iterations","750"},{"dense-threshold","64"},{"riesz-patches","n2"}};
    std::set<std::string> seen;
    for(int i=2;i<argc;++i){std::string arg=argv[i];auto equal=arg.find('=');
        if(!arg.starts_with("--")||equal==std::string::npos)throw std::invalid_argument("expected --name=value");
        auto key=arg.substr(2,equal-2);
        if(!options.contains(key)||!seen.insert(key).second)throw std::invalid_argument("unknown or duplicate fixed-state option");
        options[key]=arg.substr(equal+1);
    }
    auto integer=[&](std::string key){std::size_t n;int v=std::stoi(options.at(key),&n);if(n!=options.at(key).size())throw std::invalid_argument("invalid integer");return v;};
    auto real=[&](std::string key){std::size_t n;double v=std::stod(options.at(key),&n);if(n!=options.at(key).size()||!std::isfinite(v))throw std::invalid_argument("invalid number");return v;};
    in.level=integer("level");in.gap=integer("gap");in.ell=integer("ell");in.threads=integer("threads");in.cap=integer("maximum-nodes");
    int graded=integer("graded");in.graded=graded;in.policy=options.at("interpolation");in.members_path=options.at("members");
    in.riesz_policy=options.at("riesz-patches");
    in.theta=real("theta");in.tolerance=real("ritz-tolerance");in.iterations=integer("ritz-iterations");in.dense=integer("dense-threshold");
    if(in.level<0||in.level>12||in.gap<1||in.gap>8||in.ell<1||in.ell>4||in.cap<1||in.cap>2000000
        ||in.threads<1||in.threads>64||(graded!=0&&graded!=1)||(in.policy!="area"&&in.policy!="arithmetic")
        ||(in.riesz_policy!="n2"&&in.riesz_policy!="archive")||in.theta<=0||in.theta>1||in.tolerance<=0||in.iterations<1||in.dense<0)
        throw std::invalid_argument("fixed-state option out of range");
    std::istringstream ids(options.at("training-ids"));std::string word;
    while(std::getline(ids,word,',')){std::size_t n;int id=std::stoi(word,&n);if(n!=word.size())throw std::invalid_argument("invalid member ID");in.ids.push_back(id);}
    if(in.ids.empty()||options.at("training-ids").back()==',')throw std::invalid_argument("empty training ID");
    std::ifstream file(in.members_path);
    if(!file)throw std::invalid_argument("cannot read explicit member table");
    std::string line;
    while(std::getline(file,line)) {
        std::istringstream row(line);int id;std::string role,extra;double x,y,c,a,phase;
        if(!(row>>id>>role>>x>>y>>c>>a>>phase)||row>>extra)throw std::invalid_argument("invalid member row");
        if(role!="train"&&role!="test"&&role!="shift"&&role!="pure")throw std::invalid_argument("invalid member role");
        if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(c)||!std::isfinite(a)||!std::isfinite(phase))throw std::invalid_argument("nonfinite member");
        auto p=in.problem=="E1"?benchmarks::make_shifted_r1_paper_case(16,{x,y})
            :benchmarks::make_parameterized_boundary_gaussian_s_paper_case(16,c,a,phase,80,{x,y});
        in.members.push_back({id,role,std::move(p)});
    }
    if(in.members.empty()||in.members.size()>50)throw std::invalid_argument("member count outside 1..50");
    std::set<int> member_ids,training_ids;
    for(const auto& m:in.members)if(m.id<0||!member_ids.insert(m.id).second)throw std::invalid_argument("duplicate or negative member ID");
    for(int id:in.ids){
        if(!training_ids.insert(id).second)throw std::invalid_argument("duplicate training ID");
        bool found=false;for(const auto& m:in.members)if(m.id==id&&m.role=="train")found=true;
        if(!found)throw std::invalid_argument("selected member is missing or audit-only");
    }
    std::sort(in.ids.begin(),in.ids.end());omp_set_num_threads(in.threads);
    return in;
}
inline TriMesh coarse_mesh(const Input& in) {
    auto initial=alod::make_problem(in.problem).initial_mesh;
    if(std::ldexp(static_cast<double>(initial.elems.size()),in.level+in.gap+(in.graded?2:0))>in.cap)
        throw std::runtime_error("preallocation reference mesh limit exceeded");
    auto mesh=refine_mesh_nvb(initial,in.level).mesh;
    if(in.graded)mesh=bisect_newest_vertex(mesh,{0}).mesh;
    return mesh;
}
template<class Values> void vector(const Values& v) {
    std::cout<<'[';for(int i=0;i<static_cast<int>(v.size());++i){if(i)std::cout<<',';std::cout<<v[i];}std::cout<<']';
}
inline void matrix(const Eigen::MatrixXd& m) {
    std::cout<<'[';for(int i=0;i<m.rows();++i){if(i)std::cout<<',';vector(m.row(i).transpose().eval());}std::cout<<']';
}
inline void complex_matrix(const ComplexMatrix& m) {
    std::cout<<'[';for(int i=0;i<m.rows();++i){if(i)std::cout<<',';std::cout<<'[';
        for(int j=0;j<m.cols();++j){if(j)std::cout<<',';std::cout<<'['<<m(i,j).real()<<','<<m(i,j).imag()<<']';}std::cout<<']';}std::cout<<']';
}
}
