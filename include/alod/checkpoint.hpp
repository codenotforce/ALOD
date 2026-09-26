#pragma once
#include "alod/adaptive.hpp"
#include "alod/mesh_state.hpp"
#include <filesystem>
namespace alod {
// Checkpoints at a training boundary or after an ell decision are resumable,
// but only Accepted snapshots may be used as published audit states.
enum class CheckpointPhase { Accepted=0, BeforeTraining=1, AfterEll=2 };
struct Checkpoint {
    MeshState coarse{lod2d::TriMesh{}}, fine{lod2d::TriMesh{}};
    Sparse P_node,P_elem,P_dg;
    ComplexSparseMatrix lod_trial,lod_reduced;
    int format_version=0;
    std::string geometry_file;
    ComplexMatrix raw_kernel,phi,values,warm_full;
    std::vector<int> coarse_marks,reference_marks,computed_ids;
    AdaptiveCursor cursor;
    int ell=1,last_check=-1,next_event=0,revision=0;
    CheckpointPhase phase=CheckpointPhase::Accepted;
    bool mesh_changed=false,check_pending=false;
    std::uint64_t committed_lines=0;
    std::uint64_t journal_hash=14695981039346656037ULL;
    std::string journal,mathematics_key,config_json,members_text,space_identity;
    lod2d::RefineOutput reference() const {return {fine.mesh,P_node,P_elem,P_dg};}
};
// Borrowed immutable geometry for synchronous publication. The owner must not
// mutate these objects during the view's lifetime; create a new view on refinement.
struct CheckpointGeometryView {
    const MeshState &coarse,&fine;
    const Sparse &P_node,&P_elem,&P_dg;
    mutable std::string object_name,coarse_hash,fine_hash;
};
std::string json_string(const std::string&);
std::uint64_t journal_hash(const std::string&,std::uint64_t prefix=14695981039346656037ULL);
std::string matrix_hash(const ComplexMatrix&);
std::string checkpoint_metadata(const Checkpoint&, bool has_basis=false,int format=0,const CheckpointGeometryView* geometry=nullptr);
// Immutable, checksummed little-endian files. A temporary is verified before
// rename; latest is replaced only after the new checkpoint is committed.
std::filesystem::path save_checkpoint(const std::filesystem::path& directory,const Checkpoint&, const ComplexSparseMatrix* accepted_trial=nullptr,
    const ComplexSparseMatrix* accepted_reduced=nullptr,bool share_geometry=false,const CheckpointGeometryView* geometry=nullptr);
// Check the complete checksum while reading only the metadata into memory.
std::string inspect_checkpoint(const std::filesystem::path&,std::uint64_t maximum_bytes=UINT64_MAX);
Checkpoint load_checkpoint(const std::filesystem::path& file,std::uint64_t maximum_bytes=UINT64_MAX);
}
