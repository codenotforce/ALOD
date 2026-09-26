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
constexpr std::uint64_t magic_basis=0x3254504b43444f4cULL;
constexpr std::uint64_t magic_shared=0x3354504b43444f4cULL;
constexpr std::uint64_t magic_mesh=0x314853454d444f4cULL;
struct Writer {
    bool writing=false;
    Writer()=default;
    std::ofstream out;std::uint64_t hash=14695981039346656037ULL;
    explicit Writer(const std::filesystem::path& p):writing(true),out(p,std::ios::binary){if(!out)throw std::runtime_error("cannot create checkpoint temporary");}
    void u(std::uint64_t v){char bytes[8];for(int i=0;i<8;++i){unsigned char c=(v>>(8*i))&255;bytes[i]=c;hash=(hash^c)*1099511628211ULL;}if(writing)out.write(bytes,8);if(writing&&!out)throw std::runtime_error("checkpoint write failed");}
    void real(double v){if(!std::isfinite(v))throw std::invalid_argument("checkpoint contains nonfinite value");u(std::bit_cast<std::uint64_t>(v));}
    void str(const std::string& s){u(s.size());if(writing)out.write(s.data(),s.size());for(unsigned char c:s){hash=(hash^c)*1099511628211ULL;}if(writing&&!out)throw std::runtime_error("checkpoint write failed");}
    void ints(const std::vector<int>& a){u(a.size());for(int v:a)u(static_cast<std::uint64_t>(static_cast<std::int64_t>(v)));}
    void dense(const ComplexMatrix& a){u(a.rows());u(a.cols());for(int j=0;j<a.cols();++j)for(int i=0;i<a.rows();++i){real(a(i,j).real());real(a(i,j).imag());}}
    void sparse(const Sparse& a){u(a.rows());u(a.cols());u(a.nonZeros());for(int j=0;j<a.outerSize();++j)for(Sparse::InnerIterator it(a,j);it;++it){u(it.row());u(it.col());real(it.value());}}
    void sparse_complex(const ComplexSparseMatrix& a){u(a.rows());u(a.cols());u(a.nonZeros());for(int j=0;j<a.outerSize();++j)for(ComplexSparseMatrix::InnerIterator it(a,j);it;++it){u(it.row());u(it.col());real(it.value().real());real(it.value().imag());}}
    void mesh(const MeshState& state){const auto& m=state.mesh;u(m.nodes.size());for(auto p:m.nodes){real(p.x());real(p.y());}u(m.elems.size());for(auto t:m.elems)for(int v:t)u(v);ints(m.dirichlet);u(m.boundary_edges.size());for(auto e:m.boundary_edges){u(e.nodes[0]);u(e.nodes[1]);u(static_cast<int>(e.tag));}
        u(state.next_id);u(state.version);for(const auto* ids:{&state.elements,&state.ancestry}){u(ids->size());for(auto id:*ids){u(id.id);u(id.parent);u(id.generation);}}}
};
struct Reader {
    std::ifstream in;std::uint64_t remaining,hash=14695981039346656037ULL;
    std::array<unsigned char,65536> buffer;std::size_t pos=0,available=0;
    Reader(const std::filesystem::path& p,std::uint64_t cap):in(p,std::ios::binary){auto size=std::filesystem::file_size(p);if(!in||size<24||size>cap)throw std::runtime_error("checkpoint size outside allowed bounds");remaining=size;}
    unsigned char byte(){
        if(!remaining)throw std::runtime_error("truncated checkpoint");
        if(pos==available){in.read(reinterpret_cast<char*>(buffer.data()),std::min<std::uint64_t>(buffer.size(),remaining));available=in.gcount();pos=0;if(!available)throw std::runtime_error("truncated checkpoint");}
        --remaining;unsigned char c=buffer[pos++];hash=(hash^c)*1099511628211ULL;return c;
    }
    std::uint64_t u(){if(remaining<8)throw std::runtime_error("truncated checkpoint word");std::uint64_t v=0;for(int i=0;i<8;++i)v|=std::uint64_t(byte())<<(i*8);return v;}
    int integer(){auto v=static_cast<std::int64_t>(u());if(v<INT32_MIN||v>INT32_MAX)throw std::runtime_error("checkpoint integer out of range");return v;}
    double real(){double v=std::bit_cast<double>(u());if(!std::isfinite(v))throw std::runtime_error("nonfinite checkpoint value");return v;}
    std::size_t count(std::uint64_t width){auto n=u();if(n>remaining/width||n>INT32_MAX)throw std::runtime_error("checkpoint allocation bound exceeded");return n;}
    std::string str(){auto n=count(1);std::string s(n,' ');for(auto& c:s)c=byte();return s;}
    std::vector<int> ints(){std::vector<int> a(count(8));for(auto& v:a)v=integer();return a;}
    ComplexMatrix dense(){auto n=u(),m=u();if(n>INT32_MAX||m>INT32_MAX||(m&&n>remaining/(16*m)))throw std::runtime_error("checkpoint dense dimensions invalid");ComplexMatrix a(n,m);for(int j=0;j<a.cols();++j)for(int i=0;i<a.rows();++i){double re=real(),im=real();a(i,j)={re,im};}return a;}
    Sparse sparse(){auto n=u(),m=u();auto count_=count(24);if(n>INT32_MAX||m>INT32_MAX)throw std::runtime_error("checkpoint sparse dimensions invalid");std::vector<Eigen::Triplet<double>> entries;entries.reserve(count_);std::pair<int,int> previous{-1,-1};
        for(std::size_t j=0;j<count_;++j){int row=integer(),col=integer();double value=real();if(row<0||col<0||row>=n||col>=m||std::pair{col,row}<=previous)throw std::runtime_error("invalid checkpoint sparse coordinate order");previous={col,row};entries.emplace_back(row,col,value);}Sparse a(n,m);a.setFromTriplets(entries.begin(),entries.end());return a;}
    ComplexSparseMatrix sparse_complex(){auto n=u(),m=u();auto count_=count(32);if(n>INT32_MAX||m>INT32_MAX)throw std::runtime_error("checkpoint basis dimensions invalid");
        ComplexSparseMatrix a(n,m);a.reserve(count_);std::pair<int,int> previous{-1,-1};int column=0;if(m)a.startVec(0);
        for(std::size_t j=0;j<count_;++j){int row=integer(),col=integer();double re=real(),im=real();if(row<0||col<0||row>=n||col>=m||std::pair{col,row}<=previous)throw std::runtime_error("invalid checkpoint basis coordinate order");while(column<col)a.startVec(++column);a.insertBack(row,col)=Complex(re,im);previous={col,row};}
        while(column+1<m)a.startVec(++column);a.finalize();return a;}
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
std::string geometry_name(std::uint64_t hash){std::ostringstream out;out<<"mesh-"<<std::hex<<std::setw(16)<<std::setfill('0')<<hash<<".bin";return out.str();}
std::filesystem::path geometry_path(const std::filesystem::path& checkpoint,const std::string& name){
    if(name.size()!=25||!name.starts_with("mesh-")||!name.ends_with(".bin")
       ||name.substr(5,16).find_first_not_of("0123456789abcdef")!=std::string::npos)
        throw std::runtime_error("invalid shared checkpoint geometry name");
    return checkpoint.parent_path()/"meshes"/name;
}
void write_geometry(Writer& w,const Checkpoint& s){w.u(magic_mesh);w.mesh(s.coarse);w.mesh(s.fine);w.sparse(s.P_node);w.sparse(s.P_elem);w.sparse(s.P_dg);}
std::uint64_t finish_read(Reader& r){auto checksum=r.hash;if(r.u()!=checksum||r.remaining)throw std::runtime_error("checkpoint checksum or trailing-data mismatch");return checksum;}
void verify_geometry(const std::filesystem::path& file,const std::string& name,std::uint64_t cap){
    Reader r(file,cap);if(r.u()!=magic_mesh)throw std::runtime_error("invalid shared geometry format");
    while(r.remaining>8)r.byte();if(geometry_name(finish_read(r))!=name)throw std::runtime_error("shared geometry identity mismatch");
}
std::string save_geometry(const std::filesystem::path& directory,const Checkpoint& s){
    Writer fingerprint;write_geometry(fingerprint,s);const auto name=geometry_name(fingerprint.hash);
    const auto folder=directory/"meshes";std::filesystem::create_directories(folder);
    const auto file=folder/name;
    if(std::filesystem::exists(file)){
        try{verify_geometry(file,name,UINT64_MAX);return name;}
        catch(const std::exception&){
            // Preserve corruption evidence; the supplied live mesh regenerates
            // exactly the content-addressed object needed by older snapshots.
            const auto recovery=folder/"recovery";std::filesystem::create_directories(recovery);
            auto old=recovery/name;int attempt=0;
            while(std::filesystem::exists(old))old=recovery/(name+"-"+std::to_string(++attempt));
            std::filesystem::rename(file,old);sync_file(folder);
        }
    }
    auto temp=file;temp+=".tmp";
    {Writer w(temp);write_geometry(w,s);auto checksum=w.hash;w.u(checksum);w.out.flush();if(!w.out)throw std::runtime_error("geometry flush failed");}
    verify_geometry(temp,name,UINT64_MAX);sync_file(temp);std::filesystem::rename(temp,file);sync_file(folder);return name;
}
void atomic_text(const std::filesystem::path& target,const std::string& text){auto temp=target;temp+=".tmp";{std::ofstream out(temp,std::ios::binary);out<<text;out.flush();if(!out)throw std::runtime_error("checkpoint pointer write failed");}sync_file(temp);std::filesystem::rename(temp,target);sync_file(target.parent_path());}
}
std::uint64_t journal_hash(const std::string& text,std::uint64_t hash){for(unsigned char c:text)hash=(hash^c)*1099511628211ULL;return hash;}
std::string json_string(const std::string& text){std::ostringstream out;out<<'"';for(unsigned char c:text){if(c=='"'||c=='\\')out<<'\\'<<c;else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;else out<<c;}out<<'"';return out.str();}
std::string matrix_hash(const ComplexMatrix& matrix){FingerprintBuilder hash;hash.add_i64(matrix.rows());hash.add_i64(matrix.cols());for(int j=0;j<matrix.cols();++j)for(int i=0;i<matrix.rows();++i)hash.add_complex(matrix(i,j));return hash.finish();}
std::string checkpoint_metadata(const Checkpoint& s,bool has_basis,int format){has_basis=has_basis||s.lod_trial.cols();if(!format)format=s.format_version?s.format_version:(has_basis?2:1);std::ostringstream out;out<<"{\"format\":"<<format<<",\"state_id\":"<<s.cursor.state_id<<",\"coarse_cycle\":"<<s.cursor.coarse_cycle<<",\"reference_sweep\":"<<s.cursor.reference_sweep<<",\"phase\":"<<int(s.phase)<<",\"ell\":"<<s.ell<<",\"committed_lines\":"<<s.committed_lines<<",\"journal_hash\":"<<json_string(std::to_string(s.journal_hash))<<",\"journal\":"<<json_string(s.journal)<<",\"mathematics_key\":"<<json_string(s.mathematics_key)<<",\"config\":"<<s.config_json<<",\"members_text\":"<<json_string(s.members_text)<<",\"space_identity\":"<<json_string(s.space_identity)
        <<",\"coarse_mesh_hash\":"<<json_string(mesh_fingerprint(s.coarse.mesh))<<",\"reference_mesh_hash\":"<<json_string(mesh_fingerprint(s.fine.mesh))
        <<",\"kernel_hash\":"<<json_string(matrix_hash(s.raw_kernel))<<",\"dictionary_hash\":"<<json_string(matrix_hash(s.phi))
        <<",\"solution_hash\":"<<json_string(matrix_hash(s.values))<<",\"warm_hash\":"<<json_string(matrix_hash(s.warm_full))
        <<",\"config_hash\":"<<json_string(std::to_string(journal_hash(s.config_json)))<<",\"members_hash\":"<<json_string(std::to_string(journal_hash(s.members_text)))
        <<",\"format_identity\":\"ALOD-checkpoint-le-ieee754-v"<<format<<"\"}";return out.str();}
std::filesystem::path save_checkpoint(const std::filesystem::path& directory,const Checkpoint& s,const ComplexSparseMatrix* accepted_trial,const ComplexSparseMatrix* accepted_reduced,bool share_geometry){
    const auto& basis=accepted_trial?*accepted_trial:s.lod_trial;
    const auto& reduced=accepted_reduced?*accepted_reduced:s.lod_reduced;
    const int format=share_geometry||reduced.cols()?3:(basis.cols()?2:1);
    std::filesystem::create_directories(directory);std::ostringstream name;name<<"state-"<<std::setw(6)<<std::setfill('0')<<s.cursor.state_id<<"-phase-"<<int(s.phase)<<"-ell-"<<s.ell<<".bin";
    auto file=directory/name.str(),temp=file;temp+=".tmp";
    if(std::filesystem::exists(file))throw std::runtime_error("refusing to overwrite an immutable checkpoint");
    const std::string geometry=share_geometry?save_geometry(directory,s):"";
    {Writer w(temp);w.u(format==3?magic_shared:basis.cols()?magic_basis:magic);w.str(checkpoint_metadata(s,basis.cols()>0,format));
        if(format==3){w.str(geometry);w.u(basis.cols()>0);w.u(reduced.cols()>0);}
        w.str(s.mathematics_key);w.str(s.config_json);w.str(s.members_text);w.str(s.space_identity);w.str(s.journal);
        for(int v:{s.cursor.state_id,s.cursor.coarse_cycle,s.cursor.reference_sweep,s.ell,s.last_check,s.next_event,s.revision,int(s.phase),int(s.mesh_changed),int(s.check_pending)})w.u(static_cast<std::uint64_t>(static_cast<std::int64_t>(v)));
        w.u(s.committed_lines);w.u(s.journal_hash);if(geometry.empty()){w.mesh(s.coarse);w.mesh(s.fine);w.sparse(s.P_node);w.sparse(s.P_elem);w.sparse(s.P_dg);}
        for(const auto* m:{&s.raw_kernel,&s.phi,&s.values,&s.warm_full})w.dense(*m);
        w.ints(s.coarse_marks);w.ints(s.reference_marks);w.ints(s.computed_ids);if(basis.cols())w.sparse_complex(basis);if(format==3&&reduced.cols())w.sparse_complex(reduced);auto checksum=w.hash;w.u(checksum);w.out.flush();if(!w.out)throw std::runtime_error("checkpoint flush failed");}
    // Verify in constant memory instead of decoding a second dense snapshot.
    {Reader r(temp,UINT64_MAX);while(r.remaining>8)r.byte();auto checksum=r.hash;if(r.u()!=checksum||r.remaining)throw std::runtime_error("checkpoint verification failed");}
    sync_file(temp);std::filesystem::rename(temp,file);atomic_text(directory/"latest",name.str()+"\n");return file;
}
std::string inspect_checkpoint(const std::filesystem::path& file,std::uint64_t cap){
    Reader r(file,cap);auto signature=r.u();if(signature!=magic&&signature!=magic_basis&&signature!=magic_shared)throw std::runtime_error("unsupported checkpoint format");
    auto metadata=r.str();const auto geometry=signature==magic_shared?r.str():"";
    while(r.remaining>8)r.byte();finish_read(r);
    if(!geometry.empty())verify_geometry(geometry_path(file,geometry),geometry,cap-std::filesystem::file_size(file));
    if(signature==magic_shared)metadata.insert(metadata.size()-1,",\"geometry_file\":"+json_string(geometry));
    return metadata;
}
Checkpoint load_checkpoint(const std::filesystem::path& file,std::uint64_t cap){
    Reader r(file,cap);auto signature=r.u();if(signature!=magic&&signature!=magic_basis&&signature!=magic_shared)throw std::runtime_error("unsupported checkpoint format");auto metadata=r.str();Checkpoint s;
    s.format_version=signature==magic_shared?3:signature==magic_basis?2:1;
    bool has_basis=signature==magic_basis,has_reduced=false;
    if(signature==magic_shared){s.geometry_file=r.str();auto b=r.u(),a=r.u();if(b>1||a>1||a>b)throw std::runtime_error("invalid checkpoint cache flags");has_basis=b;has_reduced=a;}
    s.mathematics_key=r.str();s.config_json=r.str();s.members_text=r.str();s.space_identity=r.str();s.journal=r.str();
    s.cursor={r.integer(),r.integer(),r.integer()};s.ell=r.integer();s.last_check=r.integer();s.next_event=r.integer();s.revision=r.integer();int phase=r.integer();s.mesh_changed=r.integer();s.check_pending=r.integer();s.phase=static_cast<CheckpointPhase>(phase);s.committed_lines=r.u();s.journal_hash=r.u();
    if(s.geometry_file.empty()){s.coarse=r.mesh();s.fine=r.mesh();s.P_node=r.sparse();s.P_elem=r.sparse();s.P_dg=r.sparse();}
    else {
        Reader g(geometry_path(file,s.geometry_file),cap-std::filesystem::file_size(file));
        if(g.u()!=magic_mesh)throw std::runtime_error("invalid geometry format");
        s.coarse=g.mesh();s.fine=g.mesh();s.P_node=g.sparse();s.P_elem=g.sparse();s.P_dg=g.sparse();
        if(geometry_name(finish_read(g))!=s.geometry_file)throw std::runtime_error("shared geometry identity mismatch");
    }
    s.raw_kernel=r.dense();s.phi=r.dense();s.values=r.dense();s.warm_full=r.dense();s.coarse_marks=r.ints();s.reference_marks=r.ints();s.computed_ids=r.ints();if(has_basis)s.lod_trial=r.sparse_complex();if(has_reduced)s.lod_reduced=r.sparse_complex();
    auto checksum=r.hash;if(r.u()!=checksum||r.remaining)throw std::runtime_error("checkpoint checksum or trailing-data mismatch");
    if(phase<0||phase>2||s.cursor.state_id<0||s.ell<1||s.ell>4||s.P_node.rows()!=s.fine.mesh.nodes.size()||s.P_node.cols()!=s.coarse.mesh.nodes.size()||s.P_elem.rows()!=s.fine.mesh.elems.size()||s.P_elem.cols()!=s.coarse.mesh.elems.size()||s.P_dg.rows()!=3*s.fine.mesh.elems.size()||s.P_dg.cols()!=3*s.coarse.mesh.elems.size())throw std::runtime_error("checkpoint state dimensions invalid");
    const auto nh=s.fine.mesh.nodes.size(),nH=s.coarse.mesh.nodes.size();
    if(has_basis&&(phase!=0||s.lod_trial.rows()!=nh||s.lod_trial.cols()!=nH-lod2d::helmholtz::dirichlet_nodes(s.coarse.mesh).size()))throw std::runtime_error("checkpoint LOD basis contract invalid");
    if(has_reduced&&(s.lod_reduced.rows()!=s.lod_trial.cols()||s.lod_reduced.cols()!=s.lod_trial.cols()))throw std::runtime_error("checkpoint reduced operator dimensions invalid");
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
