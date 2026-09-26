#include "alod/problems.hpp"
#include "mesh/refine.h"
#include <iostream>
#include <iomanip>
#include <cmath>
int main(int argc,char** argv){try{
    bool paper=argc==2&&std::string(argv[1])=="--paper-v4";
    if(argc>1&&!paper)throw std::invalid_argument("expected --paper-v4");
    std::cout<<std::setprecision(17);
    for(double k:std::vector<double>{8.,16.,32.,64.,128.}){
        if(k==128&&!paper)continue;
        auto mesh=alod::make_problem("E1",k).initial_mesh;int level=0;
        while(alod::mesh_diameter(mesh)>std::acos(-1.)/k){mesh=lod2d::refine_mesh_nvb(mesh,1).mesh;++level;}
        int gap=paper?static_cast<int>(std::log2(k)):4;
        auto fine=lod2d::refine_mesh_nvb(mesh,gap).mesh;
        double H=alod::mesh_diameter(mesh),h=alod::mesh_diameter(fine);
        std::cout<<"{\"wavenumber\":"<<k<<",\"level\":"<<level<<",\"reference_level\":"<<level+gap
          <<",\"gap\":"<<gap<<",\"H0\":"<<H<<",\"h0\":"<<h<<",\"kH0\":"<<k*H<<",\"kh0\":"<<k*h<<",\"H0_over_h0\":"<<H/h
          <<",\"coarse_nodes\":"<<mesh.nodes.size()<<",\"reference_nodes\":"<<fine.nodes.size()
          <<",\"coarse_free\":"<<mesh.nodes.size()-mesh.dirichlet.size()<<",\"reference_free\":"<<fine.nodes.size()-fine.dirichlet.size()<<"}\n";
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
