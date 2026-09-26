#pragma once
#include <cstdlib>
#ifdef __linux__
#include <sys/prctl.h>
#include <unistd.h>
#include <signal.h>
#endif
// Prevent orphan workers from writing into a run after its supervisor dies.
inline void supervised_parent(){
#ifdef __linux__
    if(const char* parent=std::getenv("ALOD_PARENT_PID")){
        prctl(PR_SET_PDEATHSIG,SIGTERM);
        if(getppid()!=std::strtol(parent,nullptr,10))std::exit(125);
    }
#endif
}
