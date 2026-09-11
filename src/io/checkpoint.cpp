#include "alod/checkpoint.hpp"
#include "helmholtz/boundary.h"
#include "../lod/fingerprint.hpp"
#include <bit>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <set>
#include <array>
#include <cstdio>
#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif
namespace alod {
namespace {
constexpr std::uint64_t magic=0x3154504b43444f4cULL;
struct Writer {
    std::ofstream out;std::uint64_t hash=14695981039346656037ULL;
    explicit Writer(const std::filesystem::path& p):out(p,std::ios::binary){if(!out)throw std::runtime_error("cannot create checkpoint temporary");}
    void u(std::uint64_t v){for(int i=0;i<8;++i){unsigned char c=(v>>(8*i))&255;out.put(c);hash=(hash^c)*1099511628211ULL;}if(!out)throw std::runtime_error("checkpoint write failed");}
    void real(double v){if(!std::isfinite(v))throw std::invalid_argument("checkpoint contains nonfinite value");u(std::bit_cast<std::uint64_t>(v));}
    void str(const std::string& s){u(s.size());for(unsigned char c:s){out.put(c);hash=(hash^c)*1099511628211ULL;}if(!out)throw std::runtime_error("checkpoint write failed");}
    void ints(const std::vector<int>& a){u(a.size());for(int v:a)u(static_cast<std::uint64_t>(static_cast<std::int64_t>(v)));}
    void dense(const ComplexMatrix& a){u(a.rows());u(a.cols());for(int j=0;j<a.cols();++j)for(int i=0;i<a.rows();++i){real(a(i,j).real());real(a(i,j).imag());}}
    void sparse(const Sparse& a){u(a.rows());u(a.cols());u(a.nonZeros());for(int j=0;j<a.outerSize();++j)for(Sparse::InnerIterator it(a,j);it;++it){u(it.row());u(it.col());real(it.value());}}
    void mesh(const MeshState& state){const auto& m=state.mesh;u(m.nodes.size());for(auto p:m.nodes){real(p.x());real(p.y());}u(m.elems.size());for(auto t:m.elems)for(int v:t)u(v);ints(m.dirichlet);u(m.boundary_edges.size());for(auto e:m.boundary_edges){u(e.nodes[0]);u(e.nodes[1]);u(static_cast<int>(e.tag));}
        u(state.next_id);u(state.version);for(const auto* ids:{&state.elements,&state.ancestry}){u(ids->size());for(auto id:*ids){u(id.id);u(id.parent);u(id.generation);}}}
};
struct Reader {
    std::ifstream in;std::uint64_t remaining,hash=14695981039346656037ULL;
    Reader(const std::filesystem::path& p,std::uint64_t cap):in(p,std::ios::binary){auto size=std::filesystem::file_size(p);if(!in||size<24||size>cap)throw std::runtime_error("checkpoint size outside allowed bounds");remaining=size;}
    unsigned char byte(){if(!remaining--)throw std::runtime_error("truncated checkpoint");int c=in.get();if(c==EOF)throw std::runtime_error("truncated checkpoint");hash=(hash^static_cast<unsigned char>(c))*1099511628211ULL;return c;}
    std::uint64_t u(){if(remaining<8)throw std::runtime_error("truncated checkpoint word");std::uint64_t v=0;for(int i=0;i<8;++i)v|=std::uint64_t(byte())<<(i*8);return v;}
    int integer(){auto v=static_cast<std::int64_t>(u());if(v<INT32_MIN||v>INT32_MAX)throw std::runtime_error("checkpoint integer out of range");return v;}
    double real(){double v=std::bit_cast<double>(u());if(!std::isfinite(v))throw std::runtime_error("nonfinite checkpoint value");return v;}
    std::size_t count(std::uint64_t width){auto n=u();if(n>remaining/width||n>100000000)throw std::runtime_error("checkpoint allocation bound exceeded");return n;}
    std::string str(){auto n=count(1);std::string s(n,' ');for(auto& c:s)c=byte();return s;}
    std::vector<int> ints(){std::vector<int> a(count(8));for(auto& v:a)v=integer();return a;}
    ComplexMatrix dense(){auto n=u(),m=u();if(n>2000000||m>2000000||(m&&n>remaining/(16*m)))throw std::runtime_error("checkpoint dense dimensions invalid");ComplexMatrix a(n,m);for(int j=0;j<a.cols();++j)for(int i=0;i<a.rows();++i){double re=real(),im=real();a(i,j)={re,im};}return a;}
    Sparse sparse(){auto n=u(),m=u();auto count_=count(24);if(n>30000000||m>30000000)throw std::runtime_error("checkpoint sparse dimensions invalid");std::vector<Eigen::Triplet<double>> entries;entries.reserve(count_);std::pair<int,int> previous{-1,-1};
        for(std::size_t j=0;j<count_;++j){int row=integer(),col=integer();double value=real();if(row<0||col<0||row>=n||col>=m||std::pair{col,row}<=previous)throw std::runtime_error("invalid checkpoint sparse coordinate order");previous={col,row};entries.emplace_back(row,col,value);}Sparse a(n,m);a.setFromTriplets(entries.begin(),entries.end());return a;}
    MeshState mesh(){lod2d::TriMesh m;m.nodes.resize(count(16));for(auto& p:m.nodes){double x=real(),y=real();p={x,y};}m.elems.resize(count(24));for(auto& t:m.elems)for(auto& v:t){v=integer();if(v<0||v>=m.nodes.size())throw std::runtime_error("checkpoint triangle node out of range");}m.dirichlet=ints();m.boundary_edges.resize(count(24));for(auto& e:m.boundary_edges){e.nodes={integer(),integer()};e.tag=static_cast<lod2d::BoundaryTag>(integer());}
        lod2d::helmholtz::validate_boundary_tags(m);MeshState state(std::move(m));state.next_id=u();state.version=u();
        for(auto* ids:{&state.elements,&state.ancestry}){ids->resize(count(24));for(auto& id:*ids){id.id=u();id.parent=u();id.generation=integer();}}
        if(state.elements.size()!=state.mesh.elems.size())throw std::runtime_error("checkpoint lineage size mismatch");std::set<std::uint64_t> seen;for(auto id:state.ancestry){if(id.id>=state.next_id||!seen.insert(id.id).second||id.generation<0||(id.parent!=UINT64_MAX&&!seen.contains(id.parent)))throw std::runtime_error("invalid checkpoint ancestry");}
        std::set<std::uint64_t> leaves;for(auto id:state.elements)if(!seen.contains(id.id)||!leaves.insert(id.id).second)throw std::runtime_error("invalid checkpoint leaf IDs");return state;}
};
void sync_file(const std::filesystem::path& path){
#ifndef _WIN32
    int fd=::open(path.c_str(),O_RDONLY);if(fd<0)throw std::runtime_error("cannot open checkpoint for sync");int rc=::fsync(fd);::close(fd);if(rc)throw std::runtime_error("checkpoint sync failed");
#endif
}
void atomic_text(const std::filesystem::path& target,const std::string& text){auto temp=target;temp+=".tmp";{std::ofstream out(temp,std::ios::binary);out<<text;out.flush();if(!out)throw std::runtime_error("checkpoint pointer write failed");}sync_file(temp);std::filesystem::rename(temp,target);sync_file(target.parent_path());}
}
std::uint64_t journal_hash(const std::string& text,std::uint64_t hash){for(unsigned char c:text)hash=(hash^c)*1099511628211ULL;return hash;}
std::string json_string(const std::string& text){std::ostringstream out;out<<'"';for(unsigned char c:text){if(c=='"'||c=='\\')out<<'\\'<<c;else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;else out<<c;}out<<'"';return out.str();}
std::string matrix_hash(const ComplexMatrix& matrix){FingerprintBuilder hash;hash.add_i64(matrix.rows());hash.add_i64(matrix.cols());for(int j=0;j<matrix.cols();++j)for(int i=0;i<matrix.rows();++i)hash.add_complex(matrix(i,j));return hash.finish();}
std::string checkpoint_metadata(const Checkpoint& s){std::ostringstream out;out<<"{\"format\":1,\"state_id\":"<<s.cursor.state_id<<",\"coarse_cycle\":"<<s.cursor.coarse_cycle<<",\"reference_sweep\":"<<s.cursor.reference_sweep<<",\"phase\":"<<int(s.phase)<<",\"ell\":"<<s.ell<<",\"committed_lines\":"<<s.committed_lines<<",\"journal_hash\":"<<json_string(std::to_string(s.journal_hash))<<",\"journal\":"<<json_string(s.journal)<<",\"mathematics_key\":"<<json_string(s.mathematics_key)<<",\"config\":"<<s.config_json<<",\"members_text\":"<<json_string(s.members_text)<<",\"space_identity\":"<<json_string(s.space_identity)
        <<",\"coarse_mesh_hash\":"<<json_string(mesh_fingerprint(s.coarse.mesh))<<",\"reference_mesh_hash\":"<<json_string(mesh_fingerprint(s.fine.mesh))
        <<",\"kernel_hash\":"<<json_string(matrix_hash(s.raw_kernel))<<",\"dictionary_hash\":"<<json_string(matrix_hash(s.phi))
        <<",\"solution_hash\":"<<json_string(matrix_hash(s.values))<<",\"warm_hash\":"<<json_string(matrix_hash(s.warm_full))
        <<",\"config_hash\":"<<json_string(std::to_string(journal_hash(s.config_json)))<<",\"members_hash\":"<<json_string(std::to_string(journal_hash(s.members_text)))
        <<",\"format_identity\":\"ALOD-checkpoint-le-ieee754-v1\"}";return out.str();}
std::filesystem::path save_checkpoint(const std::filesystem::path& directory,const Checkpoint& s){
    std::filesystem::create_directories(directory);std::ostringstream name;name<<"state-"<<std::setw(6)<<std::setfill('0')<<s.cursor.state_id<<"-phase-"<<int(s.phase)<<"-ell-"<<s.ell<<".bin";
    auto file=directory/name.str(),temp=file;temp+=".tmp";
    if(std::filesystem::exists(file))throw std::runtime_error("refusing to overwrite an immutable checkpoint");
    {Writer w(temp);w.u(magic);w.str(checkpoint_metadata(s));w.str(s.mathematics_key);w.str(s.config_json);w.str(s.members_text);w.str(s.space_identity);w.str(s.journal);
        for(int v:{s.cursor.state_id,s.cursor.coarse_cycle,s.cursor.reference_sweep,s.ell,s.last_check,s.next_event,s.revision,int(s.phase),int(s.mesh_changed),int(s.check_pending)})w.u(static_cast<std::uint64_t>(static_cast<std::int64_t>(v)));
        w.u(s.committed_lines);w.u(s.journal_hash);w.mesh(s.coarse);w.mesh(s.fine);w.sparse(s.P_node);w.sparse(s.P_elem);w.sparse(s.P_dg);
        for(const auto* m:{&s.raw_kernel,&s.phi,&s.values,&s.warm_full})w.dense(*m);
        w.ints(s.coarse_marks);w.ints(s.reference_marks);w.ints(s.computed_ids);auto checksum=w.hash;w.u(checksum);w.out.flush();if(!w.out)throw std::runtime_error("checkpoint flush failed");}
    // Verify in constant memory instead of decoding a second dense snapshot.
    {Reader r(temp,1073741824);while(r.remaining>8)r.byte();auto checksum=r.hash;if(r.u()!=checksum||r.remaining)throw std::runtime_error("checkpoint verification failed");}
    sync_file(temp);std::filesystem::rename(temp,file);atomic_text(directory/"latest",name.str()+"\n");return file;
}
Checkpoint load_checkpoint(const std::filesystem::path& file,std::uint64_t cap){
    Reader r(file,cap);if(r.u()!=magic)throw std::runtime_error("unsupported checkpoint format");auto metadata=r.str();Checkpoint s;
    s.mathematics_key=r.str();s.config_json=r.str();s.members_text=r.str();s.space_identity=r.str();s.journal=r.str();
    s.cursor={r.integer(),r.integer(),r.integer()};s.ell=r.integer();s.last_check=r.integer();s.next_event=r.integer();s.revision=r.integer();int phase=r.integer();s.mesh_changed=r.integer();s.check_pending=r.integer();s.phase=static_cast<CheckpointPhase>(phase);s.committed_lines=r.u();s.journal_hash=r.u();
    s.coarse=r.mesh();s.fine=r.mesh();s.P_node=r.sparse();s.P_elem=r.sparse();s.P_dg=r.sparse();s.raw_kernel=r.dense();s.phi=r.dense();s.values=r.dense();s.warm_full=r.dense();s.coarse_marks=r.ints();s.reference_marks=r.ints();s.computed_ids=r.ints();
    auto checksum=r.hash;if(r.u()!=checksum||r.remaining)throw std::runtime_error("checkpoint checksum or trailing-data mismatch");
    if(phase<0||phase>2||s.cursor.state_id<0||s.ell<1||s.ell>4||s.P_node.rows()!=s.fine.mesh.nodes.size()||s.P_node.cols()!=s.coarse.mesh.nodes.size()||s.P_elem.rows()!=s.fine.mesh.elems.size()||s.P_elem.cols()!=s.coarse.mesh.elems.size()||s.P_dg.rows()!=3*s.fine.mesh.elems.size()||s.P_dg.cols()!=3*s.coarse.mesh.elems.size())throw std::runtime_error("checkpoint state dimensions invalid");
    const auto nh=s.fine.mesh.nodes.size(),nH=s.coarse.mesh.nodes.size();
    if(s.cursor.coarse_cycle<0||s.cursor.reference_sweep<0||s.cursor.reference_sweep>16||s.last_check>s.cursor.state_id||s.next_event<0||s.revision<0
        ||(s.raw_kernel.cols()&&(s.raw_kernel.rows()!=nh||s.raw_kernel.cols()>24))
        ||(s.phi.cols()&&(s.phi.rows()!=nh||s.phi.cols()>24))
        ||(s.warm_full.cols()&&(s.warm_full.rows()!=nH||s.warm_full.cols()>16))
        ||(phase==0&&(s.values.rows()!=nh||s.values.cols()!=s.computed_ids.size())))throw std::runtime_error("checkpoint matrix/cursor contract invalid");
    for(int mark:s.coarse_marks)if(mark<0||mark>=s.coarse.mesh.elems.size())throw std::runtime_error("checkpoint coarse mark invalid");
    for(int mark:s.reference_marks)if(mark<0||mark>=s.fine.mesh.elems.size())throw std::runtime_error("checkpoint reference mark invalid");
    if(metadata!=checkpoint_metadata(s))throw std::runtime_error("checkpoint metadata mismatch");return s;
}
}
