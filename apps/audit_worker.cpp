#include "audit_run.hpp"
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif
#ifdef __linux__
#include <malloc.h>
#endif
namespace {
std::string field(){
    unsigned char bytes[4];if(!std::cin.read(reinterpret_cast<char*>(bytes),4))throw std::runtime_error("truncated worker request");
    std::uint32_t n=0;for(int i=0;i<4;++i)n|=std::uint32_t(bytes[i])<<(8*i);
    if(n>1024*1024)throw std::runtime_error("worker field exceeds limit");
    std::string s(n,'\0');if(!std::cin.read(s.data(),n))throw std::runtime_error("truncated worker field");return s;
}
void timing_file(const std::string& path){
#ifdef _WIN32
    _putenv_s("ALOD_TIMING_FILE",path.c_str());
#else
    setenv("ALOD_TIMING_FILE",path.c_str(),1);
#endif
}
struct Redirect {
    std::streambuf *out,*err;
    Redirect(std::ostream& a,std::ostream& b):out(std::cout.rdbuf(a.rdbuf())),err(std::cerr.rdbuf(b.rdbuf())){}
    ~Redirect(){std::cout.flush();std::cerr.flush();std::cout.rdbuf(out);std::cerr.rdbuf(err);}
};
}
int audit_worker_main(){
#ifdef _WIN32
    _setmode(_fileno(stdin),_O_BINARY);
#endif
    AuditWorkerState worker;
    std::cout<<"{\"ready\":1}\n"<<std::flush;
    while(std::cin.peek()!=EOF){
        try{
            const auto command=field();
            if(command=="quit")return 0;
            if(command=="prepare"){
                worker.prepared.reset();worker.prepared_path=field();
                worker.prepared=std::make_unique<alod::Checkpoint>(alod::load_checkpoint(worker.prepared_path));
                auto metadata=alod::checkpoint_metadata(*worker.prepared);
                metadata.insert(metadata.size()-1,",\"geometry_file\":"+alod::json_string(worker.prepared->geometry_file));
                std::cout<<"{\"ok\":true,\"metadata\":"<<alod::json_string(metadata)<<"}\n"<<std::flush;
            }else if(command=="audit"){
                const auto output=field(),error=field(),timing=field();
                const auto count_text=field();std::size_t used;int count=std::stoi(count_text,&used);
                if(used!=count_text.size()||count<2||count>128)throw std::runtime_error("invalid worker argument count");
                std::vector<std::string> args;for(int i=0;i<count;++i)args.push_back(field());
                std::vector<char*> pointers;for(auto& arg:args)pointers.push_back(arg.data());
                std::ofstream out(output),err(error);if(!out||!err)throw std::runtime_error("cannot open audit output");
                timing_file(timing);int status;
                {Redirect redirect(out,err);status=audit_main(pointers.size(),pointers.data(),&worker);}
                worker.prepared.reset();worker.prepared_path.clear();
                if(status){worker.hierarchy.reset();worker.hierarchy_key.clear();}
#ifdef __linux__
                malloc_trim(0); // Long-lived workers must release obsolete state arenas.
#endif
                std::cout<<"{\"ok\":true,\"returncode\":"<<status<<"}\n"<<std::flush;
            }else throw std::runtime_error("unknown audit worker command");
        }catch(const std::exception& e){
            std::cout<<"{\"ok\":false,\"error\":"<<alod::json_string(e.what())<<"}\n"<<std::flush;
            // Do not attempt to resynchronize a possibly partial binary request.
            return 1;
        }
    }
    return 0;
}
