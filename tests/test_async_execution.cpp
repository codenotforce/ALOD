#include "alod/execution.hpp"
#include "alod/estimator.hpp"
#include "alod/problems.hpp"
#include <barrier>
#include <iostream>
using namespace alod;
void require(bool value,const char* message) {if(!value)throw std::runtime_error(message);}
int observed_team() {
    const int requested=execution_threads(4);
    int observed=0;
    #pragma omp parallel num_threads(requested)
    {
        #pragma omp master
        observed=omp_get_num_threads();
    }
    return observed;
}
void check_completion_handoff() {
    // Background wins: a live foreground observes a larger actual team.
    for(bool fail:{false,true}) {
        ConcurrentBudgetPair pair(4,2,"test_background_wins");
        std::promise<void> release;auto gate=release.get_future();
        auto task=background_task(pair.budget(0),[&,pair]{
            auto done=pair.completion(0);gate.wait();
            if(fail)throw std::runtime_error("injected paired task failure");
        });
        auto foreground_done=pair.completion(1);
        int before,after;bool caught=false;
        {ExecutionScope scope(pair.budget(1));before=observed_team();release.set_value();
            try{task.get();}catch(const std::runtime_error&){caught=true;}
            after=observed_team();}
        foreground_done.finish();foreground_done.finish();
        require(before==2&&after==4,"background completion did not expand foreground team");
        require(caught==fail,"paired task exception lost");
        require(pair.budget(0)->load()==0&&pair.budget(1)->load()==0,"completed shares were not released");
    }
    // Foreground wins: the waiting child observes the returned share.
    ConcurrentBudgetPair pair(4,2,"test_foreground_wins");
    std::promise<void> ready,release;auto gate=release.get_future();
    auto task=background_task(pair.budget(0),[&,pair]{
        auto done=pair.completion(0);const int before=observed_team();
        ready.set_value();gate.wait();return std::pair{before,observed_team()};
    });
    bool foreground_failed=false;
    try {
        auto done=pair.completion(1);ExecutionScope scope(pair.budget(1));ready.get_future().wait();
        throw std::runtime_error("injected foreground failure");
    }catch(const std::runtime_error&){foreground_failed=true;}
    release.set_value();const auto teams=task.get();
    require(foreground_failed&&teams.first==2&&teams.second==4,"foreground exception did not expand background team");
    ConcurrentBudgetPair serial(1,0,"test_no_background");
    {auto done=serial.completion(0);}
    require(serial.budget(1)->load()==1,"absent task changed serial budget");
}
int main(){try{
    check_completion_handoff();
    const auto coarse=lod2d::refine_mesh_nvb(make_problem("E2").initial_mesh,3).mesh;
    LodLimits limits;limits.threads=4;
    LodSpace space(coarse,lod2d::refine_mesh_nvb(coarse,2),16,1,
                   InterpolationPolicy::ManuscriptAreaWeighted,limits);
    AdditiveKernelRieszContext riesz(space);
    ComplexMatrix rhs(space.fine().nodes.size(),3);
    for(int i=0;i<rhs.rows();++i)for(int j=0;j<rhs.cols();++j)
        rhs(i,j)=Complex(std::sin(i+2*j),std::cos(2*i+j));
    const auto expected=riesz.apply(rhs);
    const auto calls=riesz.applications(),columns=riesz.applied_columns();
    std::barrier rendezvous(2);
    auto budget=thread_budget(2);
    auto task=background_task(budget,[&]{rendezvous.arrive_and_wait();
        ComplexMatrix result;
        for(int i=0;i<8;++i)result=riesz.apply_action(rhs);
        return result;
    });
    {ExecutionScope foreground(thread_budget(2));rendezvous.arrive_and_wait();
        for(int i=0;i<8;++i){const auto actual=riesz.apply(rhs);
            require((actual.values-expected.values).norm()<1e-12,"concurrent estimator changed");
            require((actual.eta-expected.eta).norm()<1e-12,"concurrent indicators changed");
        }
    }
    budget->store(4);
    require((task.get()-expected.values).norm()<1e-12,"concurrent action changed");
    require(riesz.applications()==calls+16&&riesz.applied_columns()==columns+48,"concurrent counters lost updates");
    const int previous=omp_get_max_threads();
    {ExecutionScope scope(budget);budget->store(2);require(execution_threads(32)==2,"budget split ignored");
        budget->store(4);require(execution_threads(32)==4,"budget return ignored");}
    require(omp_get_max_threads()==previous,"OpenMP scope leaked");
    auto failing=background_task(thread_budget(1),[]{throw std::runtime_error("injected");});
    bool caught=false;try{failing.get();}catch(const std::runtime_error&){caught=true;}
    require(caught,"background exception lost");
    std::future<int> tail;
    std::atomic<int> completed{0};
    {TaskLane lane(thread_budget(2));
        const auto worker=lane.submit([]{return std::this_thread::get_id();}).get();
        for(int i=0;i<8;++i)tail=lane.submit([&,worker]{
            require(std::this_thread::get_id()==worker,"pipeline recreated its worker");
            return ++completed;
        });
    }
    require(tail.get()==8&&completed==8,"lane destruction did not drain queued work");
    std::cout<<"Concurrent Riesz, budget handoff and exception propagation passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
