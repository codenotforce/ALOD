#pragma once
#include <memory>
#include <mutex>
#include <list>
#include <unordered_map>
#include <string>
namespace alod {
// Run-owned bounded retention. Active contexts retain their own shared owners;
// eviction never invalidates a live factor. Keys contain exact matrix bytes.
class LocalFactorCache {
    struct Entry {std::shared_ptr<void> value;std::size_t bytes;std::list<std::string>::iterator position;};
    std::mutex mutex_;
    std::unordered_map<std::string,Entry> entries_;
    std::list<std::string> order_;
    std::size_t limit_,bytes_=0,hits_=0,misses_=0;
public:
    explicit LocalFactorCache(std::size_t limit):limit_(limit){}
    std::shared_ptr<void> find(const std::string& key){
        std::lock_guard lock(mutex_);auto it=entries_.find(key);
        if(it==entries_.end()){++misses_;return {};}
        ++hits_;order_.splice(order_.end(),order_,it->second.position);return it->second.value;
    }
    void insert(std::string key,std::shared_ptr<void> value,std::size_t bytes){
        bytes+=2*key.capacity()+256;
        std::lock_guard lock(mutex_);if(bytes>limit_||entries_.contains(key))return;
        while(bytes_+bytes>limit_&&!order_.empty()){
            auto it=entries_.find(order_.front());bytes_-=it->second.bytes;entries_.erase(it);order_.pop_front();
        }
        order_.push_back(key);entries_.emplace(std::move(key),Entry{std::move(value),bytes,std::prev(order_.end())});bytes_+=bytes;
    }
    std::size_t hits(){std::lock_guard lock(mutex_);return hits_;}
    std::size_t misses(){std::lock_guard lock(mutex_);return misses_;}
    std::size_t bytes(){std::lock_guard lock(mutex_);return bytes_;}
};
}
