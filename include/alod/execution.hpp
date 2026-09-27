#pragma once
#include "alod/timing.hpp"
#include <algorithm>
#include <atomic>
#include <future>
#include <memory>
#include <utility>
#include <condition_variable>
#include <functional>
#include <thread>

namespace alod {
// A task owns its OpenMP budget. A waiting parent may return its share to the
// child; the next parallel region observes it without changing another thread's ICV.
using ThreadBudget = std::shared_ptr<std::atomic<int>>;
inline thread_local ThreadBudget execution_budget;
inline int execution_threads(int requested) {
    if (!execution_budget) return requested;
    const int available = std::max(1, execution_budget->load());
    omp_set_num_threads(available);
    return std::max(1, std::min(requested, available));
}
class ExecutionScope {
    ThreadBudget previous_ = execution_budget;
    int previous_threads_ = omp_get_max_threads();
    int previous_dynamic_ = omp_get_dynamic();
public:
    explicit ExecutionScope(ThreadBudget budget) {
        execution_budget = std::move(budget);
        omp_set_dynamic(0);
        execution_threads(previous_threads_);
    }
    ~ExecutionScope() {
        execution_budget = std::move(previous_);
        omp_set_num_threads(previous_threads_);
        omp_set_dynamic(previous_dynamic_);
    }
    ExecutionScope(const ExecutionScope&) = delete;
    ExecutionScope& operator=(const ExecutionScope&) = delete;
};
inline ThreadBudget thread_budget(int threads) {
    return std::make_shared<std::atomic<int>>(std::max(1, threads));
}
inline bool asynchronous_execution(int threads) {
    const char* serial = std::getenv("ALOD_ASYNC_DISABLE");
    return threads > 1 && !(serial && std::string(serial) == "1");
}
inline bool asynchronous_buffers_fit(std::size_t extra_bytes) {
    // Additional dense buffers only; existing factor/patch ceilings still apply.
    std::size_t limit=1024ULL*1024*1024;
    if(const char* value=std::getenv("ALOD_ASYNC_BUFFER_BYTES")) {
        const std::string text(value);
        if(text.empty()||text.find_first_not_of("0123456789")!=std::string::npos)
            throw std::invalid_argument("invalid asynchronous buffer budget");
        limit=std::stoull(text);
    }
    return extra_bytes<=limit;
}
template<class F> auto background_task(ThreadBudget budget, F task) {
    // All captures are owned except the explicitly borrowed immutable state.
    // The future must be destroyed/joined before that state, also on exceptions.
    return std::async(std::launch::async,
        [budget=std::move(budget), task=std::move(task),
         file=PhaseTimer::thread_file, state=PhaseTimer::current_state]() mutable {
            ExecutionScope scope(std::move(budget));
            PhaseTimer::thread_file=std::move(file);
            PhaseTimer::current_state=state;
            return task();
        });
}
// One persistent OpenMP team for bandwidth-heavy pipeline stages. Keeping load
// and integration on this lane avoids competing quadrature teams and repeated
// thread-pool creation. Capacity is one running task plus one queued task.
class TaskLane {
    std::mutex mutex_;
    std::condition_variable changed_;
    std::function<void()> next_;
    bool stopping_=false;
    std::thread worker_;
public:
    explicit TaskLane(ThreadBudget budget):worker_(
        [this,budget=std::move(budget),file=PhaseTimer::thread_file,state=PhaseTimer::current_state]{
            ExecutionScope scope(budget);
            PhaseTimer::thread_file=file;PhaseTimer::current_state=state;
            for(;;){
                std::function<void()> task;
                {std::unique_lock lock(mutex_);
                    changed_.wait(lock,[&]{return stopping_||bool(next_);});
                    if(!next_)return;
                    task=std::move(next_);next_={};
                }
                changed_.notify_all();task();
            }
        }) {}
    ~TaskLane(){
        {std::lock_guard lock(mutex_);stopping_=true;}
        changed_.notify_all();worker_.join();
    }
    TaskLane(const TaskLane&)=delete;
    TaskLane& operator=(const TaskLane&)=delete;
    template<class F> auto submit(F task){
        auto packaged=std::make_shared<std::packaged_task<std::invoke_result_t<F>()>>(std::move(task));
        auto future=packaged->get_future();
        {std::unique_lock lock(mutex_);
            changed_.wait(lock,[&]{return stopping_||!next_;});
            if(stopping_)throw std::runtime_error("submission to a stopped task lane");
            next_=[packaged]{(*packaged)();};
        }
        changed_.notify_all();return future;
    }
};
}
