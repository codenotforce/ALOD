#pragma once
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <mutex>
#include <omp.h>

// A separate, append-only runtime log; mathematical journals remain unchanged.
struct PhaseTimer {
    using Clock=std::chrono::steady_clock;
    inline static thread_local int current_state=-1;
    inline static thread_local std::string thread_file;
    inline static std::mutex output_mutex;
    const char* phase; int state,previous_state; Clock::time_point start=Clock::now();
    static void record(const char* event,const char* phase,int state,double seconds,int workers=0){
        const char* file=thread_file.empty()?std::getenv("ALOD_TIMING_FILE"):thread_file.c_str();if(!file||!*file)return;
        std::lock_guard<std::mutex> lock(output_mutex);
        std::ofstream out(file,std::ios::app);
        out<<std::setprecision(17)<<"{\"event\":\""<<event<<"\",\"phase\":\""<<phase
           <<"\",\"state_id\":"<<state<<",\"seconds\":"<<seconds
           <<",\"unix_seconds\":"<<std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count()
           <<",\"workers\":"<<workers<<"}\n";
    }
    static void counter(const char* name,std::size_t value,int state){
        const char* file=thread_file.empty()?std::getenv("ALOD_TIMING_FILE"):thread_file.c_str();if(!file||!*file)return;
        std::lock_guard<std::mutex> lock(output_mutex);
        std::ofstream out(file,std::ios::app);
        out<<"{\"event\":\"counter\",\"phase\":\""<<name<<"\",\"state_id\":"<<state
           <<",\"value\":"<<value<<"}\n";
    }
    PhaseTimer(const char* name,int id):phase(name),state(id<0?current_state:id),previous_state(current_state){current_state=state;record("start",phase,state,0);}
    ~PhaseTimer(){record("end",phase,state,std::chrono::duration<double>(Clock::now()-start).count());current_state=previous_state;}
    static void probe(int requested){
        int observed=0;
        #pragma omp parallel num_threads(requested)
        {
            #pragma omp single
            observed=omp_get_num_threads();
        }
        record("runtime","openmp_team",-1,0,observed);
        if(observed!=requested)throw std::runtime_error("OpenMP runtime prevented the requested thread team");
    }
};
