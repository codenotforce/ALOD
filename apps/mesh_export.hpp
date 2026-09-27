#pragma once
#include "mesh/types.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>
inline void export_mesh_pair(const std::string& path,const lod2d::TriMesh& coarse,const lod2d::TriMesh& fine){
    if(path.empty())return;
    std::ofstream out(path);if(!out)throw std::runtime_error("cannot open mesh export");out<<std::setprecision(17);
    auto emit=[&](const lod2d::TriMesh& mesh){
        out<<"{\"nodes\":[";
        for(int i=0;i<static_cast<int>(mesh.nodes.size());++i){if(i)out<<',';out<<'['<<mesh.nodes[i].x()<<','<<mesh.nodes[i].y()<<']';}
        out<<"],\"elements\":[";
        for(int i=0;i<static_cast<int>(mesh.elems.size());++i){if(i)out<<',';auto t=mesh.elems[i];out<<'['<<t[0]<<','<<t[1]<<','<<t[2]<<']';}
        out<<"]}";
    };
    out<<"{\"coarse\":";emit(coarse);out<<",\"reference\":";emit(fine);out<<"}\n";
    out.flush();if(!out)throw std::runtime_error("mesh export failed");
}
