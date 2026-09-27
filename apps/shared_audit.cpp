#include "shared_audit.hpp"
#include "audit_run.hpp"
#include "alod/timing.hpp"
#include <atomic>
#include <iostream>
#include <algorithm>
#include <bit>
#include <fstream>
#include <map>
#include <mutex>
#include <thread>
#ifdef __linux__
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>
#include <cstring>
namespace {
using namespace alod;
std::size_t budget(const char* name,std::size_t fallback){
    const char* value=std::getenv(name);if(!value)return fallback;
    std::string text(value);if(text.empty()||text.find_first_not_of("0123456789")!=std::string::npos)throw std::invalid_argument("invalid shared audit budget");
    return std::stoull(text);
}
std::size_t hierarchy_bytes(const LodHierarchyData& h){
    std::size_t bytes=(h.coarse.nodes.size()+h.reference.mesh.nodes.size())*64+(h.coarse.elems.size()+h.reference.mesh.elems.size())*64;
    auto sparse=[&](const auto& a){bytes+=a.nonZeros()*(sizeof(typename std::decay_t<decltype(a)>::Scalar)+sizeof(int))+(a.cols()+1)*sizeof(int);};
    for(const Sparse* a:{&h.interpolation,&h.coarse_basis,&h.energy,&h.reference.P_node,&h.reference.P_elem,&h.reference.P_dg,&h.operators.stiffness,&h.operators.mass,&h.operators.boundary_mass})sparse(*a);
    sparse(h.operators.system);return bytes+h.operators.element_blocks.size()*sizeof(Eigen::Matrix3cd);
}
struct Snapshot {
    std::unique_ptr<Checkpoint> state;
    std::shared_ptr<const LodHierarchyData> hierarchy;
    std::string metadata,key;
    std::shared_ptr<std::atomic<std::size_t>> retained;
    std::size_t bytes=0;
    ~Snapshot(){if(retained)retained->fetch_sub(bytes);}
};
void receive(int fd,char* data,std::size_t size){
    while(size){auto n=recv(fd,data,size,0);if(n<=0)throw std::runtime_error("audit connection closed");data+=n;size-=n;}
}
std::string field(int fd){
    unsigned char bytes[4];receive(fd,reinterpret_cast<char*>(bytes),4);
    std::uint32_t n=0;for(int i=0;i<4;++i)n|=std::uint32_t(bytes[i])<<(8*i);
    if(n>1024*1024)throw std::runtime_error("audit request field too large");
    std::string text(n,'\0');receive(fd,text.data(),n);return text;
}
void reply(int fd,std::string text){
    text+='\n';const char* data=text.data();std::size_t size=text.size();
    while(size){auto n=send(fd,data,size,MSG_NOSIGNAL);if(n<=0)throw std::runtime_error("audit reply connection closed");data+=n;size-=n;}
}
class Server {
    int listener=-1;
    std::string path;
    std::thread acceptor;
    std::vector<std::thread> handlers;
    std::mutex mutex;
    std::map<std::string,std::shared_ptr<Snapshot>> snapshots;
    std::vector<int> connections;
    std::atomic<bool> stopped=false;
    std::atomic<int> clients=0;
    std::shared_ptr<std::atomic<std::size_t>> retained=std::make_shared<std::atomic<std::size_t>>(0);
    std::size_t cap=budget("ALOD_SHARED_SNAPSHOT_BYTES",256ULL*1024*1024);
    int slots=budget("ALOD_SHARED_AUDIT_WORKERS",2);
    void serve(int fd){
        AuditWorkerState worker;
        try{
            reply(fd,"{\"ready\":1,\"pid\":"+std::to_string(getpid())+"}");
            while(!stopped){
                auto command=field(fd);
                if(command=="quit")break;
                if(command=="shutdown"){stopped=true;reply(fd,"{\"ok\":true}");break;}
                if(command=="prepare"){
                    worker.prepared.reset();worker.snapshot_owner.reset();worker.prepared_path=field(fd);
                    std::shared_ptr<Snapshot> saved;
                    {std::lock_guard lock(mutex);auto it=snapshots.find(worker.prepared_path);
                        if(it!=snapshots.end()){saved=it->second;snapshots.erase(it);}}
                    std::string metadata;
                    if(saved){
                        // The in-memory state never bypasses immutable disk integrity.
                        metadata=inspect_checkpoint(worker.prepared_path);
                        if(metadata!=saved->metadata)throw std::runtime_error("shared audit snapshot differs from durable checkpoint");
                        worker.prepared=std::move(saved->state);worker.hierarchy=saved->hierarchy;
                        worker.hierarchy_key=saved->key;worker.snapshot_owner=saved;
                    }else{
                        worker.prepared=std::make_unique<Checkpoint>(load_checkpoint(worker.prepared_path));
                        metadata=checkpoint_metadata(*worker.prepared);
                        metadata.insert(metadata.size()-1,",\"geometry_file\":"+json_string(worker.prepared->geometry_file));
                    }
                    reply(fd,"{\"ok\":true,\"memory_snapshot\":"+std::string(saved?"true":"false")+",\"metadata\":"+json_string(metadata)+"}");
                }else if(command=="audit"){
                    auto output=field(fd),error=field(fd),timing=field(fd),count_text=field(fd);
                    std::size_t used;int count=std::stoi(count_text,&used);
                    if(used!=count_text.size()||count<2||count>128)throw std::runtime_error("invalid audit arguments");
                    std::vector<std::string> args;for(int i=0;i<count;++i)args.push_back(field(fd));
                    std::vector<char*> argv;for(auto& arg:args)argv.push_back(arg.data());
                    std::ofstream out(output),err(error);if(!out||!err)throw std::runtime_error("cannot open audit output");
                    PhaseTimer::thread_file=timing;
                    PhaseTimer::counter("audit_memory_snapshot",worker.snapshot_owner?1:0,-1);
                    int status=audit_main(argv.size(),argv.data(),&worker,&out,&err);
                    out.flush();err.flush();worker.prepared.reset();worker.prepared_path.clear();
                    // No shared state may accumulate behind the bounded queue.
                    worker.hierarchy.reset();worker.hierarchy_key.clear();worker.snapshot_owner.reset();
                    PhaseTimer::thread_file.clear();
                    reply(fd,"{\"ok\":true,\"returncode\":"+std::to_string(status)+"}");
                }else throw std::runtime_error("unknown audit request");
            }
        }catch(const std::exception& e){try{reply(fd,"{\"ok\":false,\"error\":"+json_string(e.what())+"}");}catch(...){} }
        {std::lock_guard lock(mutex);connections.erase(std::remove(connections.begin(),connections.end(),fd),connections.end());}
        close(fd);--clients;
    }
public:
    explicit Server(std::string name):path(std::move(name)){
        if(slots<1||slots>128)throw std::invalid_argument("invalid shared audit slots");
        sockaddr_un addr{};addr.sun_family=AF_UNIX;
        if(path.size()>=sizeof(addr.sun_path))throw std::invalid_argument("audit socket path too long");
        std::memcpy(addr.sun_path,path.c_str(),path.size()+1);
        listener=socket(AF_UNIX,SOCK_STREAM,0);
        if(listener<0||bind(listener,reinterpret_cast<sockaddr*>(&addr),sizeof(addr))||listen(listener,slots+1))throw std::runtime_error("cannot create audit socket");
        acceptor=std::thread([this]{
            while(!stopped){
                int fd=accept(listener,nullptr,nullptr);if(fd<0)break;
                // One spare connection permits orderly shutdown after queue drain.
                if(clients.fetch_add(1)>=slots+1){--clients;close(fd);continue;}
                std::lock_guard lock(mutex);connections.push_back(fd);handlers.emplace_back([this,fd]{serve(fd);});
            }
        });
    }
    void publish(const std::filesystem::path& file,Checkpoint&& state,const LodSpace& space,const CheckpointGeometryView& geometry){
        if(!cap)return;
        const auto hierarchy=space.hierarchy();
        std::size_t bytes=hierarchy_bytes(*hierarchy)+sizeof(Complex)*(state.phi.size()+state.values.size())
            +(space.trial().nonZeros()+space.reduced().nonZeros())*(sizeof(Complex)+sizeof(int))
            +(space.trial().cols()+space.reduced().cols()+2)*sizeof(int)+state.journal.size()+state.config_json.size()+state.members_text.size()+4096;
        std::lock_guard lock(mutex);
        while(!snapshots.empty()&&(snapshots.size()>=2||retained->load()+bytes>cap))snapshots.erase(snapshots.begin());
        if(retained->load()+bytes>cap){PhaseTimer::counter("audit_snapshot_spilled",1,state.cursor.state_id);return;}
        auto saved=std::make_shared<Snapshot>();saved->bytes=bytes;saved->retained=retained;retained->fetch_add(bytes);
        saved->metadata=checkpoint_metadata(state,true,3,&geometry);
        state.format_version=3;state.geometry_file=geometry.object_name;
        saved->metadata.insert(saved->metadata.size()-1,",\"geometry_file\":"+json_string(state.geometry_file));
        state.raw_kernel.resize(0,0);state.warm_full.resize(0,0);
        state.lod_trial=space.trial();state.lod_reduced=space.reduced();
        saved->state=std::make_unique<Checkpoint>(std::move(state));saved->hierarchy=hierarchy;
        saved->key=saved->state->geometry_file+":"+std::to_string(std::bit_cast<std::uint64_t>(space.operators().wavenumber))+":"+
            (hierarchy->policy==InterpolationPolicy::ManuscriptAreaWeighted?"area":"arithmetic");
        snapshots[file.string()]=std::move(saved);
        PhaseTimer::counter("audit_retained_snapshot_bytes",retained->load(),PhaseTimer::current_state);
    }
    void drain(){
        {std::ofstream out(path+".done.tmp");out<<"0\n";}
        std::filesystem::rename(path+".done.tmp",path+".done");
        while(!stopped)std::this_thread::sleep_for(std::chrono::milliseconds(20));
    }
    ~Server(){
        stopped=true;shutdown(listener,SHUT_RDWR);close(listener);if(acceptor.joinable())acceptor.join();
        {std::lock_guard lock(mutex);for(int fd:connections)shutdown(fd,SHUT_RDWR);}
        for(auto& thread:handlers)if(thread.joinable())thread.join();unlink(path.c_str());
    }
};
Server* active=nullptr;
}
int with_shared_audits(const std::function<int()>& adaptive){
    const char* socket=std::getenv("ALOD_AUDIT_SOCKET");if(!socket||!*socket)return adaptive();
    try{
        Server server(socket);active=&server;int status=adaptive();active=nullptr;
        if(status){std::cout.flush();std::cerr.flush();std::_Exit(status);}
        server.drain();return status;
    }catch(const std::exception& e){std::cerr<<"shared audit server: "<<e.what()<<'\n';return 1;}
}
void publish_audit_snapshot(const std::filesystem::path& file,alod::Checkpoint&& state,
    const alod::LodSpace& space,const alod::CheckpointGeometryView& geometry){
    if(active)active->publish(file,std::move(state),space,geometry);
}
#else
int with_shared_audits(const std::function<int()>& adaptive){return adaptive();}
void publish_audit_snapshot(const std::filesystem::path&,alod::Checkpoint&&,const alod::LodSpace&,const alod::CheckpointGeometryView&){}
#endif
