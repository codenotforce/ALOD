// Opt-in server test: large serialization only, no Helmholtz solve.
#include "alod/checkpoint.hpp"
#include "alod/problems.hpp"
#include <chrono>
#include <iostream>
#include <numeric>
using namespace alod;
int main(int argc,char** argv){try{
    if(argc!=2)throw std::invalid_argument("provide a new output directory");
    std::filesystem::path directory=argv[1];
    if(!std::filesystem::create_directory(directory))throw std::invalid_argument("output directory already exists");
    auto start=std::chrono::steady_clock::now();
    auto elapsed=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    auto H=lod2d::refine_mesh_nvb(make_problem("E1").initial_mesh,2).mesh;
    auto ref=lod2d::refine_mesh_nvb(H,18);
    Checkpoint s;s.coarse=MeshState(std::move(H));s.fine=MeshState(std::move(ref.mesh));
    s.P_node=std::move(ref.P_node);s.P_elem=std::move(ref.P_elem);s.P_dg=std::move(ref.P_dg);
    const auto nodes=s.fine.mesh.nodes.size();
    s.values=ComplexMatrix::Constant(nodes,48,Complex(.125,-.25));
    s.computed_ids.resize(48);std::iota(s.computed_ids.begin(),s.computed_ids.end(),0);
    s.config_json="{}";s.mathematics_key="large-checkpoint-serialization-test";
    const double prepared=elapsed();auto file=save_checkpoint(directory,s);double written=elapsed();
    auto bytes=std::filesystem::file_size(file);
    if(bytes<=1073741824ULL)throw std::runtime_error("test did not exceed the old one GiB limit");
    s=Checkpoint{};auto metadata=inspect_checkpoint(file);double inspected=elapsed();
    auto restored=load_checkpoint(file);double loaded=elapsed();
    if(restored.values.rows()!=nodes||restored.values.cols()!=48||restored.values(0,0)!=Complex(.125,-.25)
       ||restored.values(nodes-1,47)!=Complex(.125,-.25)||metadata!=checkpoint_metadata(restored))
        throw std::runtime_error("large checkpoint round trip failed");
    std::cout<<"{\"bytes\":"<<bytes<<",\"nodes\":"<<nodes<<",\"columns\":48,\"prepare_seconds\":"<<prepared
             <<",\"write_verify_seconds\":"<<written-prepared<<",\"inspect_seconds\":"<<inspected-written
             <<",\"load_seconds\":"<<loaded-inspected<<",\"passed\":true}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
