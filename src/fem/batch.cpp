#include "alod/batch.hpp"
#include "alod/mesh_state.hpp"
#include <cstdlib>
#include <bit>
#include <unordered_map>
#include "alod/timing.hpp"
#include <Eigen/UmfPackSupport>
#include <omp.h>
#include <mutex>
#include <chrono>
namespace alod {
namespace {
template<class F> void parallel(int n,int threads,F f){if(threads<1)throw std::invalid_argument("batch threads must be positive");if(n==0)return;if(n==1){f(0);return;}threads=std::min(n,threads);std::exception_ptr failure;std::mutex mutex;
#pragma omp parallel for schedule(static) num_threads(threads)
    for(int i=0;i<n;++i){try{f(i);}catch(...){std::lock_guard<std::mutex> lock(mutex);if(!failure)failure=std::current_exception();}}
    if(failure)std::rethrow_exception(failure);
}
// Only share points when every quadrature context is exactly identical.
// Distinct localized features must retain their own adaptive quadrature.
bool shared_quadrature(const std::vector<Problem>& problems){
    if(problems.empty())return false;
    const auto& a=problems.front().quadrature_context;
    for(const auto& p:problems){const auto& b=p.quadrature_context;
        if(a.integrand_class!=b.integrand_class||a.feature_point!=b.feature_point||a.feature_scale!=b.feature_scale)return false;}
    return true;
}
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
}
AuditIntegrationGeometry::AuditIntegrationGeometry(const lod2d::TriMesh& mesh,int threads):mesh_(&mesh),gradients_(mesh.elems.size()){
    parallel(static_cast<int>(mesh.elems.size()),threads,[&](int e){
        auto tri=mesh.elems[e];auto a=mesh.nodes[tri[0]],b=mesh.nodes[tri[1]],c=mesh.nodes[tri[2]];
        const double det=(b.x()-a.x())*(c.y()-a.y())-(b.y()-a.y())*(c.x()-a.x());
        if(std::abs(det)<=1e-15)throw std::invalid_argument("degenerate audit triangle");
        gradients_[e]={Eigen::Vector2d(b.y()-c.y(),c.x()-b.x())/det,Eigen::Vector2d(c.y()-a.y(),a.x()-c.x())/det,Eigen::Vector2d(a.y()-b.y(),b.x()-a.x())/det};
    });
}
ComplexMatrix assemble_load_batch(const lod2d::TriMesh& mesh,const std::vector<Problem>& problems,const lod2d::helmholtz::QuadraturePolicy& q,int threads,std::vector<lod2d::helmholtz::SourceMomentData>* moments,SourceMomentReuse* reuse,const std::string& family){
    if(threads<1)throw std::invalid_argument("batch threads must be positive");
    const int n=problems.size(),elements=mesh.elems.size();
    ComplexMatrix loads=ComplexMatrix::Zero(mesh.nodes.size(),n);
    auto batch_source=std::getenv("ALOD_REFERENCE_EXECUTION")?std::function<void(const lod2d::Point2&,Complex*)>{}:lod2d::helmholtz::benchmarks::shared_source_evaluator(problems);
    if(moments){
        std::vector<lod2d::helmholtz::SourceMomentData> previous;
        std::vector<int> parents(elements,-1);std::vector<std::string> keys;
        if(reuse&&!family.empty()&&!std::getenv("ALOD_REFERENCE_EXECUTION")&&!std::getenv("ALOD_INTEGRATION_COLD")){
            keys.resize(elements);
            auto word=[](std::string& key,std::uint64_t value){for(int i=0;i<8;++i)key.push_back(char(value>>(8*i)));};
            for(int e=0;e<elements;++e)for(int v:mesh.elems[e])for(double x:{mesh.nodes[v].x(),mesh.nodes[v].y()})word(keys[e],std::bit_cast<std::uint64_t>(x));
            bool match=reuse->family==family&&moments->size()==problems.size();
            for(int j=0;match&&j<n;++j){auto& m=(*moments)[j];auto& c=problems[j].quadrature_context;
                match=m.policy==q&&m.context.integrand_class==c.integrand_class&&m.context.feature_point==c.feature_point&&m.context.feature_scale==c.feature_scale;}
            if(match){
                std::unordered_map<std::string,int> old;for(int e=0;e<static_cast<int>(reuse->elements.size());++e)old.emplace(reuse->elements[e],e);
                for(int e=0;e<elements;++e){auto it=old.find(keys[e]);if(it!=old.end())parents[e]=it->second;}
                previous=std::move(*moments);
            }
            reuse->family=family;reuse->elements=std::move(keys);
        }else if(reuse){reuse->family.clear();reuse->elements.clear();}
        moments->clear();const bool shared=shared_quadrature(problems);
        std::size_t budget=768ULL*1024*1024;
        if(const char* value=std::getenv("ALOD_SOURCE_MOMENT_BYTES")){std::string text(value);if(text.empty()||text.find_first_not_of("0123456789")!=std::string::npos)throw std::invalid_argument("invalid ALOD_SOURCE_MOMENT_BYTES");budget=std::stoull(text);}
        const long double needed=static_cast<long double>(elements)*(56.L*n+72.L*(shared?1:n));
        if(!previous.empty()&&2*needed>budget){previous.clear();std::fill(parents.begin(),parents.end(),-1);}
        if(n&&needed<=budget){
            moments->resize(n);const auto identity=mesh_fingerprint(mesh);
            for(int j=0;j<n;++j){auto& m=(*moments)[j];if(!problems[j].source)throw std::invalid_argument("empty batch source");
                m.mesh_identity=identity;m.policy=q;m.context=problems[j].quadrature_context;
                m.squared.resize(elements);m.load.resize(elements);
                m.mass=shared&&j?(*moments)[0].mass:std::make_shared<std::vector<std::array<double,9>>>(elements);
            }
            parallel(threads,threads,[&](int part){
                std::vector<lod2d::helmholtz::PhysicalTriangleQuadraturePoint> points;
                std::vector<Complex> source_values(n);
                for(int e=static_cast<long long>(elements)*part/threads;e<static_cast<long long>(elements)*(part+1)/threads;++e){
                    if(parents[e]>=0){
                        for(int j=0;j<n;++j){auto& m=(*moments)[j];auto& old=previous[j];int pe=parents[e];
                            m.squared[e]=old.squared.at(pe);m.load[e]=old.load.at(pe);
                            if(!shared||j==0)(*m.mass)[e]=old.mass->at(pe);}
                        continue;
                    }
                    if(shared&&batch_source){
                        lod2d::helmholtz::triangle_quadrature_points_into(mesh,e,q,problems[0].quadrature_context,points);
                        auto& mass=(*(*moments)[0].mass)[e];
                        for(const auto& point:points){
                            for(int a=0;a<3;++a)for(int b=0;b<3;++b)mass[3*a+b]+=point.weight*point.barycentric[a]*point.barycentric[b];
                            batch_source(point.point,source_values.data());
                            for(int j=0;j<n;++j){auto& m=(*moments)[j];const auto value=source_values[j];
                                m.squared[e]+=point.weight*std::norm(value);
                                for(int i=0;i<3;++i)m.load[e][i]+=point.weight*value*point.barycentric[i];}
                        }
                        continue;
                    }
                    for(int j=0;j<n;++j){auto& m=(*moments)[j];
                        if(!j||!shared){
                            lod2d::helmholtz::triangle_quadrature_points_into(mesh,e,q,m.context,points);
                            auto& mass=(*m.mass)[e];
                            for(const auto& point:points)for(int a=0;a<3;++a)for(int b=0;b<3;++b)
                                mass[3*a+b]+=point.weight*point.barycentric[a]*point.barycentric[b];
                        }
                        for(const auto& point:points){const auto value=problems[j].source(point.point);
                            m.squared[e]+=point.weight*std::norm(value);
                            for(int i=0;i<3;++i)m.load[e][i]+=point.weight*value*point.barycentric[i];}
                    }
                }
            });
            PhaseTimer::counter("source_reused_elements",std::count_if(parents.begin(),parents.end(),[](int p){return p>=0;}),PhaseTimer::current_state);
            parallel(n,threads,[&](int j){for(int e=0;e<elements;++e)for(int i=0;i<3;++i)
                loads(mesh.elems[e][i],j)+=(*moments)[j].load[e][i];});
            return loads;
        }
    }
    if(n>1&&mesh.elems.size()>=static_cast<std::size_t>(threads)&&shared_quadrature(problems)){
        for(const auto& problem:problems)if(!problem.source)throw std::invalid_argument("empty batch source");
        std::vector<std::array<Complex,3>> local(static_cast<std::size_t>(n)*elements);
        parallel(threads,threads,[&](int part){
            std::vector<lod2d::helmholtz::PhysicalTriangleQuadraturePoint> points;
            std::vector<Complex> source_values(n);
            for(int e=static_cast<long long>(elements)*part/threads;e<static_cast<long long>(elements)*(part+1)/threads;++e){
                lod2d::helmholtz::triangle_quadrature_points_into(mesh,e,q,problems.front().quadrature_context,points);
                for(const auto& point:points){
                    if(batch_source)batch_source(point.point,source_values.data());
                    for(int j=0;j<n;++j){auto& values=local[static_cast<std::size_t>(j)*elements+e];
                        const Complex value=batch_source?source_values[j]:problems[j].source(point.point);
                        for(int i=0;i<3;++i)values[i]+=point.weight*value*point.barycentric[i];}
                }
            }
        });
        parallel(n,threads,[&](int j){for(int e=0;e<elements;++e)for(int i=0;i<3;++i)
            loads(mesh.elems[e][i],j)+=local[static_cast<std::size_t>(j)*elements+e][i];});
        return loads;
    }
    if(n>=threads||n==0){parallel(n,threads,[&](int j){loads.col(j)=lod2d::helmholtz::assemble_helmholtz_load(mesh,problems[j].source,q,problems[j].quadrature_context);});return loads;}
    const int blocks=(threads+n-1)/n;
    std::vector<std::array<Complex,3>> local(static_cast<std::size_t>(n)*elements);
    parallel(n*blocks,threads,[&](int task){
        int j=task/blocks,part=task%blocks;
        if(!problems[j].source)throw std::invalid_argument("empty batch source");
        std::vector<lod2d::helmholtz::PhysicalTriangleQuadraturePoint> points;
        for(int e=static_cast<long long>(elements)*part/blocks;e<static_cast<long long>(elements)*(part+1)/blocks;++e){
            lod2d::helmholtz::triangle_quadrature_points_into(mesh,e,q,problems[j].quadrature_context,points);
            auto& values=local[static_cast<std::size_t>(j)*elements+e];
            for(const auto& point:points){const Complex value=problems[j].source(point.point);
                for(int i=0;i<3;++i)values[i]+=point.weight*value*point.barycentric[i];}
        }
    });
    parallel(n,threads,[&](int j){for(int e=0;e<elements;++e)for(int i=0;i<3;++i)
        loads(mesh.elems[e][i],j)+=local[static_cast<std::size_t>(j)*elements+e][i];});
    return loads;
}
static ErrorBatch integrate_error_separate(const lod2d::TriMesh& mesh,const Sparse& E,const ComplexMatrix& values,const std::vector<Problem>& problems,const lod2d::helmholtz::QuadraturePolicy& q,int threads){
    if(values.rows()!=mesh.nodes.size()||values.cols()!=problems.size())throw std::invalid_argument("error batch dimensions invalid");
    ErrorBatch result;int n=problems.size();result.exact_norm.resize(n);result.exact_error.resize(n);result.energy.resize(n);
    parallel(n,threads,[&](int j){auto& p=problems[j];result.exact_error[j]=lod2d::helmholtz::compute_helmholtz_error(mesh,values.col(j),p.wavenumber,p.exact,p.exact_gradient,q,p.quadrature_context).energy;
        result.exact_norm[j]=lod2d::helmholtz::compute_helmholtz_error(mesh,ComplexVector::Zero(values.rows()),p.wavenumber,p.exact,p.exact_gradient,q,p.quadrature_context).energy;
        result.energy[j]=std::sqrt(std::max(0.,values.col(j).dot(E.cast<Complex>()*values.col(j)).real()));});return result;
}
ErrorBatch integrate_error_batch(const lod2d::TriMesh& mesh,const Sparse& E,const ComplexMatrix& values,const std::vector<Problem>& problems,const lod2d::helmholtz::QuadraturePolicy& q,int threads,ExactIntegrationReuse* reuse,const std::string& family){
    auto result=integrate_audit_batch(mesh,E,values,values,problems,q,threads,true,nullptr,reuse,family);
    result.reference_error.resize(0);
    return result;
}
ErrorBatch integrate_audit_batch(const lod2d::TriMesh& mesh,const Sparse& E,const ComplexMatrix& values,const ComplexMatrix& reference,const std::vector<Problem>& problems,const lod2d::helmholtz::QuadraturePolicy& q,int threads,bool fused,const AuditIntegrationGeometry* geometry,ExactIntegrationReuse* reuse,const std::string& family){
    if(values.rows()!=mesh.nodes.size()||values.cols()!=problems.size()||reference.rows()!=values.rows()||reference.cols()!=values.cols())throw std::invalid_argument("audit integration dimensions invalid");
    if(geometry&&geometry->mesh()!=&mesh)throw std::invalid_argument("audit geometry belongs to a different mesh");
    ErrorBatch result;int n=problems.size();result.exact_norm.resize(n);result.exact_error.resize(n);result.reference_error.resize(n);result.energy.resize(n);
    if(!fused){result=integrate_error_separate(mesh,E,values,problems,q,threads);result.reference_error.resize(n);parallel(n,threads,[&](int j){auto& p=problems[j];result.reference_error[j]=lod2d::helmholtz::compute_helmholtz_error(mesh,reference.col(j),p.wavenumber,p.exact,p.exact_gradient,q,p.quadrature_context).energy;});return result;}
    if(threads<1)throw std::invalid_argument("batch threads must be positive");
    if(n==0)return result;
    auto batch_jet=std::getenv("ALOD_REFERENCE_EXECUTION")?std::function<void(const lod2d::Point2&,lod2d::helmholtz::benchmarks::ExactJet*)>{}:lod2d::helmholtz::benchmarks::shared_jet_evaluator(problems);
    struct Samples {std::vector<lod2d::helmholtz::PhysicalTriangleQuadraturePoint> points;std::vector<lod2d::helmholtz::benchmarks::ExactJet> jets;};
    auto word=[](std::string& key,std::uint64_t value){for(int i=0;i<8;++i)key.push_back(char(value>>(8*i)));};
    std::string sample_key;word(sample_key,family.size());sample_key+=family;
    for(int v:{q.base_triangle_order,q.gaussian_triangle_order,q.singular_triangle_order,q.max_recursive_subdivisions})word(sample_key,v);
    for(const auto& p:problems){word(sample_key,std::bit_cast<std::uint64_t>(p.wavenumber));word(sample_key,int(p.quadrature_context.integrand_class));for(double x:{p.quadrature_context.feature_point.x(),p.quadrature_context.feature_point.y(),p.quadrature_context.feature_scale})word(sample_key,std::bit_cast<std::uint64_t>(x));}
    bool cache_samples=reuse&&!family.empty()&&!std::getenv("ALOD_REFERENCE_EXECUTION")&&!std::getenv("ALOD_INTEGRATION_COLD");
    // Admit only working sets likely to fit. Streaming through an undersized LRU
    // evicts unchanged cells before the next refinement state can reuse them.
    if(cache_samples&&!mesh.elems.empty()){
        std::vector<lod2d::helmholtz::PhysicalTriangleQuadraturePoint> probe;
        std::size_t total=0,count=std::min<std::size_t>(16,mesh.elems.size());
        for(std::size_t i=0;i<count;++i){
            lod2d::helmholtz::triangle_quadrature_points_into(mesh,i*mesh.elems.size()/count,q,problems.front().quadrature_context,probe);
            total+=probe.size();
        }
        const double per_element=double(total)/count*(sizeof(lod2d::helmholtz::PhysicalTriangleQuadraturePoint)+n*sizeof(lod2d::helmholtz::benchmarks::ExactJet))+2*sample_key.size()+1024;
        cache_samples=2.*mesh.elems.size()*per_element<=reuse->budget;
        if(!cache_samples)PhaseTimer::counter("exact_jet_cache_bypassed",1,PhaseTimer::current_state);
    }
    const bool shared=n>1&&mesh.elems.size()>=static_cast<std::size_t>(threads)&&shared_quadrature(problems);
    const int blocks=shared?threads:std::max(1,(threads+n-1)/n),elements=mesh.elems.size();
    const bool staged=shared||blocks>1;
    std::vector<std::array<double,3>> element_energy(staged?static_cast<std::size_t>(n)*elements:0);
    auto finish=[&](int j,const std::array<double,3>& energy_sum){
        result.exact_error[j]=std::sqrt(std::max(0.,energy_sum[0]));result.reference_error[j]=std::sqrt(std::max(0.,energy_sum[1]));result.exact_norm[j]=std::sqrt(std::max(0.,energy_sum[2]));
        result.energy[j]=std::sqrt(std::max(0.,values.col(j).dot(E.cast<Complex>()*values.col(j)).real()));
    };
    parallel(shared?blocks:n*blocks,threads,[&](int task){int first=shared?0:task/blocks,last=shared?n:first+1,part=task%blocks;std::array<double,3> energy_sum{};
        std::vector<lod2d::helmholtz::PhysicalTriangleQuadraturePoint> points;
        auto scratch=std::make_shared<Samples>();
        std::vector<std::array<Eigen::Vector2cd,2>> discrete_gradient(last-first);
        std::vector<std::array<double,3>> l2(last-first),gradient(last-first);
        std::vector<lod2d::helmholtz::benchmarks::ExactJet> jets(shared&&batch_jet?n:0);
        for(int element=static_cast<long long>(elements)*part/blocks;element<static_cast<long long>(elements)*(part+1)/blocks;++element){auto tri=mesh.elems[element];
            std::array<Eigen::Vector2d,3> gradients;
            if(geometry)gradients=geometry->gradients(element);
            else {
                auto a=mesh.nodes[tri[0]],b=mesh.nodes[tri[1]],c=mesh.nodes[tri[2]];
                double det=(b.x()-a.x())*(c.y()-a.y())-(b.y()-a.y())*(c.x()-a.x());if(std::abs(det)<=1e-15)throw std::invalid_argument("degenerate audit triangle");
                gradients={Eigen::Vector2d(b.y()-c.y(),c.x()-b.x())/det,Eigen::Vector2d(c.y()-a.y(),a.x()-c.x())/det,Eigen::Vector2d(a.y()-b.y(),b.x()-a.x())/det};
            }
            std::shared_ptr<Samples> saved=cache_samples?nullptr:scratch;std::string key;
            if(cache_samples){
                key=sample_key;word(key,first);word(key,last);
                for(int v:tri)for(double x:{mesh.nodes[v].x(),mesh.nodes[v].y()})word(key,std::bit_cast<std::uint64_t>(x));
                saved=std::static_pointer_cast<Samples>(reuse->entries.find(key));
            }
            if(!saved||!cache_samples){
                if(cache_samples)saved=std::make_shared<Samples>();
                lod2d::helmholtz::triangle_quadrature_points_into(mesh,element,q,problems[first].quadrature_context,saved->points);
                if(cache_samples){
                    saved->jets.resize(saved->points.size()*(last-first));
                    for(int point=0;point<static_cast<int>(saved->points.size());++point){
                        auto* target=saved->jets.data()+point*(last-first);const auto& x=saved->points[point].point;
                        if(shared&&batch_jet)batch_jet(x,target);
                        else for(int j=first;j<last;++j){auto& p=problems[j];target[j-first]=p.exact_jet?p.exact_jet(x):std::make_pair(p.exact(x),p.exact_gradient(x));}
                    }
                    reuse->entries.insert(key,saved,saved->points.capacity()*sizeof(lod2d::helmholtz::PhysicalTriangleQuadraturePoint)+saved->jets.capacity()*sizeof(lod2d::helmholtz::benchmarks::ExactJet)+128);
                }
            }

            std::fill(l2.begin(),l2.end(),std::array<double,3>{});
            std::fill(gradient.begin(),gradient.end(),std::array<double,3>{});
            for(int j=first;j<last;++j){auto& dg=discrete_gradient[j-first];dg[0].setZero();dg[1].setZero();
                for(int k=0;k<3;++k){dg[0]+=values(tri[k],j)*gradients[k].cast<Complex>();dg[1]+=reference(tri[k],j)*gradients[k].cast<Complex>();}}
            std::size_t point_index=0;
            for(const auto& point:saved->points){
                if(!cache_samples&&shared&&batch_jet)batch_jet(point.point,jets.data());
                for(int j=first;j<last;++j){auto& p=problems[j];int index=j-first;
                    std::array<Complex,3> discrete{};for(int k=0;k<3;++k){discrete[0]+=point.barycentric[k]*values(tri[k],j);discrete[1]+=point.barycentric[k]*reference(tri[k],j);}
                    const auto jet=cache_samples?saved->jets[point_index*(last-first)+j-first]:shared&&batch_jet?jets[j]:(p.exact_jet?p.exact_jet(point.point):std::make_pair(p.exact(point.point),p.exact_gradient(point.point)));
                    for(int k=0;k<3;++k){
                        l2[index][k]+=point.weight*std::norm(jet.first-discrete[k]);
                        gradient[index][k]+=point.weight*(k<2?(jet.second-discrete_gradient[index][k]).squaredNorm():jet.second.squaredNorm());
                    }
                }
                ++point_index;
            }
            for(int j=first;j<last;++j)for(int k=0;k<3;++k){double v=gradient[j-first][k]+problems[j].wavenumber*problems[j].wavenumber*l2[j-first][k];
                if(!staged)energy_sum[k]+=v;else element_energy[static_cast<std::size_t>(j)*elements+element][k]=v;}
        }
        if(!staged)finish(first,energy_sum);
    });
    if(staged)parallel(n,threads,[&](int j){std::array<double,3> sum{};
        for(int e=0;e<elements;++e)for(int k=0;k<3;++k)sum[k]+=element_energy[static_cast<std::size_t>(j)*elements+e][k];finish(j,sum);});
    if(cache_samples){PhaseTimer::counter("exact_jet_cache_hits",reuse->entries.hits(),PhaseTimer::current_state);PhaseTimer::counter("exact_jet_cache_bytes",reuse->entries.bytes(),PhaseTimer::current_state);}
    return result;
}
struct ReferenceFemContext::Impl {
    ComplexSparseMatrix matrix;std::vector<int> free;int full_size=0;
    Eigen::UmfPackLU<ComplexSparseMatrix> factor;
    double factor_time=0,solve_time=0,residual=0;
};
ReferenceFemContext::ReferenceFemContext(const lod2d::helmholtz::HelmholtzOperators& ops):impl_(std::make_unique<Impl>()){
    auto& p=*impl_;auto start=Clock::now();p.full_size=ops.system.rows();std::vector<int> map(p.full_size,-1);std::vector<bool> d(p.full_size,false);
    for(int i:ops.dirichlet_nodes){if(i<0||i>=p.full_size)throw std::invalid_argument("FEM boundary out of range");d[i]=true;}
    for(int i=0;i<p.full_size;++i)if(!d[i]){map[i]=p.free.size();p.free.push_back(i);}
    std::vector<Eigen::Triplet<Complex>> entries;for(int c=0;c<ops.system.outerSize();++c)for(ComplexSparseMatrix::InnerIterator it(ops.system,c);it;++it)if(map[it.row()]>=0&&map[c]>=0)entries.emplace_back(map[it.row()],map[c],it.value());
    p.matrix.resize(p.free.size(),p.free.size());p.matrix.setFromTriplets(entries.begin(),entries.end());p.factor.compute(p.matrix);if(p.factor.info()!=Eigen::Success)throw std::runtime_error("reference FEM factorization failed");p.factor_time=seconds(start);
}
ReferenceFemContext::~ReferenceFemContext()=default;
ComplexMatrix ReferenceFemContext::solve(const ComplexMatrix& loads){auto& p=*impl_;if(loads.rows()!=p.full_size||loads.cols()<1||!loads.allFinite())throw std::invalid_argument("reference RHS invalid");auto start=Clock::now();ComplexMatrix rhs(p.free.size(),loads.cols());for(int j=0;j<rhs.rows();++j)rhs.row(j)=loads.row(p.free[j]);ComplexMatrix x=p.factor.solve(rhs);
    if(p.factor.info()!=Eigen::Success||!x.allFinite())throw std::runtime_error("reference FEM solve failed");p.residual=(p.matrix*x-rhs).norm()/std::max(1e-30,rhs.norm());if(p.residual>1e-9)throw std::runtime_error("reference FEM residual gate failed");ComplexMatrix full=ComplexMatrix::Zero(p.full_size,loads.cols());for(int j=0;j<x.rows();++j)full.row(p.free[j])=x.row(j);p.solve_time+=seconds(start);return full;}
double ReferenceFemContext::factor_seconds()const{return impl_->factor_time;}
double ReferenceFemContext::solve_seconds()const{return impl_->solve_time;}
double ReferenceFemContext::relative_residual()const{return impl_->residual;}
}
