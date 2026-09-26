#pragma once
#include "helmholtz/patch_solver.h"
#include <bit>
#include <list>
#include <memory>
#include <mutex>
#include <unordered_map>

namespace alod {
// Full byte keys, not hash-only identities. Cached local solutions are sufficient
// to avoid both factorization and solve, and cheaper to retain than LU factors.
class LodPatchCache {
    using Result=lod2d::helmholtz::HelmholtzPatchSolveResult;
    struct Entry {std::shared_ptr<const Result> result;std::size_t bytes;std::list<std::string>::iterator position;};
    mutable std::mutex mutex_;
    std::unordered_map<std::string,Entry> entries_;
    std::list<std::string> order_,previous_order_;
    std::unordered_map<std::string,Entry> previous_entries_;
    std::size_t previous_bytes_=0;
    std::size_t limit_,bytes_=0,hits_=0,misses_=0;
    std::vector<std::string> dependencies_;
    std::vector<std::uint64_t> versions_;
    std::uint64_t clock_=0;
    std::size_t dependency_bytes_=0;
    std::size_t generation_budget() const {return (limit_-dependency_bytes_)/2;}
public:
    // Called between parallel constructions. Tokens change on exact byte
    // inequality; hashes never decide correctness. Reindexing is conservative.
    std::vector<std::uint64_t> update_dependencies(std::vector<std::string> records){
        std::size_t size=records.capacity()*sizeof(std::string)+records.size()*sizeof(std::uint64_t);
        for(const auto& row:records)size+=row.capacity();
        if(size>limit_/4){std::vector<std::string>().swap(dependencies_);std::vector<std::uint64_t>().swap(versions_);dependency_bytes_=0;return {};}
        std::vector<std::uint64_t> versions(records.size());
        for(std::size_t i=0;i<records.size();++i)
            versions[i]=i<dependencies_.size()&&records[i]==dependencies_[i]?versions_[i]:++clock_;
        dependencies_=std::move(records);versions_=versions;dependency_bytes_=size;
        // A larger dependency table must not violate the total storage bound.
        while(previous_bytes_>generation_budget()&&!previous_order_.empty()){
            auto old=previous_entries_.find(previous_order_.front());previous_bytes_-=old->second.bytes;
            previous_entries_.erase(old);previous_order_.pop_front();
        }
        return versions;
    }
    explicit LodPatchCache(std::size_t bytes=256ULL*1024*1024):limit_(bytes){}
    bool fits(const lod2d::helmholtz::HelmholtzPatchSystem& s)const{
        const long double key=64.L+16.L*s.local_vertices.size()+32.L*s.helmholtz.nonZeros()+24.L*(s.constraints.array()!=0.).count()+16.L*s.rhs.size();
        return 2*key+16.L*(s.rhs.size()+s.constraints.rows()*s.rhs.cols())+512<=generation_budget();
    }
    void begin_generation(){
        std::lock_guard lock(mutex_);
        previous_entries_.clear();previous_order_.clear();
        previous_entries_=std::move(entries_);previous_order_=std::move(order_);
        previous_bytes_=bytes_;entries_.clear();order_.clear();bytes_=0;
    }
    static std::string key(const lod2d::helmholtz::HelmholtzPatchSystem& s,const lod2d::TriMesh& mesh){
        std::string out;
        auto word=[&](std::uint64_t x){for(int i=0;i<8;++i)out.push_back(char(x>>(8*i)));};
        auto real=[&](double x){word(std::bit_cast<std::uint64_t>(x));};
        word(s.local_vertices.size());real(s.wavenumber);
        for(int v:s.local_vertices){real(mesh.nodes[v].x());real(mesh.nodes[v].y());}
        word(s.helmholtz.rows());word(s.helmholtz.cols());word(s.helmholtz.nonZeros());
        for(int c=0;c<s.helmholtz.outerSize();++c)
            for(lod2d::helmholtz::ComplexSparseMatrix::InnerIterator it(s.helmholtz,c);it;++it){word(it.row());word(c);real(it.value().real());real(it.value().imag());}
        word(s.constraints.rows());word(s.constraints.cols());
        for(int c=0;c<s.constraints.cols();++c)for(int r=0;r<s.constraints.rows();++r)if(s.constraints(r,c)!=0.){word(c);word(r);real(s.constraints(r,c));}
        word(UINT64_MAX);
        word(s.rhs.rows());word(s.rhs.cols());
        for(int c=0;c<s.rhs.cols();++c)for(int r=0;r<s.rhs.rows();++r){real(s.rhs(r,c).real());real(s.rhs(r,c).imag());}
        return out;
    }
    std::shared_ptr<const Result> find_shared(const std::string& key){
        std::lock_guard lock(mutex_);auto found=entries_.find(key);
        if(found!=entries_.end()){
            ++hits_;order_.splice(order_.end(),order_,found->second.position);return found->second.result;
        }
        auto old=previous_entries_.find(key);
        if(old==previous_entries_.end()){++misses_;return {};}
        ++hits_;return old->second.result;
    }
    void insert_shared(std::string key,std::shared_ptr<const Result> result){
        // Account for both retained key copies, numerical buffers and metadata.
        const std::size_t size=2*key.capacity()+sizeof(Entry)+256+
            sizeof(lod2d::helmholtz::Complex)*(result->corrector.size()+result->multipliers.size());
        std::lock_guard lock(mutex_);
        if(size>generation_budget()||entries_.contains(key))return;
        while(bytes_+size>generation_budget()&&!order_.empty()){
            auto old=entries_.find(order_.front());bytes_-=old->second.bytes;entries_.erase(old);order_.pop_front();
        }
        order_.push_back(key);entries_.emplace(std::move(key),Entry{std::move(result),size,std::prev(order_.end())});bytes_+=size;
    }
    bool find(const std::string& key,Result& result){auto found=find_shared(key);if(!found)return false;result=*found;return true;}
    void insert(std::string key,const Result& result){insert_shared(std::move(key),std::make_shared<Result>(result));}
    std::size_t hits()const{std::lock_guard lock(mutex_);return hits_;}
    std::size_t misses()const{std::lock_guard lock(mutex_);return misses_;}
    std::size_t bytes()const{std::lock_guard lock(mutex_);return bytes_+previous_bytes_+dependency_bytes_;}
};
}
