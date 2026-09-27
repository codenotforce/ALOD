#include "alod/kernel_defect.hpp"
#include "riesz_internal.hpp"
#include <algorithm>

namespace alod {
struct KernelDefectOperator::Impl {
    using Rows=Eigen::SparseMatrix<Complex,Eigen::RowMajor>;
    AdditiveKernelRieszContext::Impl& riesz;
    AdditiveKernelRieszContext& owner;
    const ComplexSparseMatrix& original;
    Rows defect;
    std::vector<std::vector<int>> columns;
    // Each coarse row gathers patch contributions in the original group order.
    std::vector<std::vector<std::pair<int,int>>> contributions;
    std::size_t output_rows=0;
    bool local_ready=false;
    Impl(AdditiveKernelRieszContext& context,const ComplexSparseMatrix& matrix)
        :riesz(*context.impl_),owner(context),original(matrix){}
};

KernelDefectOperator::KernelDefectOperator(AdditiveKernelRieszContext& context,
    const ComplexSparseMatrix& defect,bool prepare_local):impl_(std::make_unique<Impl>(context,defect)) {
    auto& p=*impl_;auto& r=p.riesz;
    if(defect.rows()!=r.full_size||defect.cols()<1)
        throw std::invalid_argument("kernel defect dimensions invalid");
    const auto row_storage=std::size_t(defect.nonZeros())*(sizeof(Complex)+sizeof(int))
        +std::size_t(defect.rows()+1)*sizeof(int);
    if(row_storage>r.dense_limit*sizeof(Complex))
        throw std::runtime_error("kernel defect row storage resource limit exceeded");
    PhaseTimer timer("theta_fused_prepare",-1);
    p.defect=defect;p.defect.makeCompressed();
    PhaseTimer::counter("theta_defect_nonzeros",defect.nonZeros(),-1);
    PhaseTimer::counter("theta_defect_row_storage_bytes",row_storage,-1);
    if(!prepare_local)return;
    p.columns.resize(r.groups.size());p.contributions.resize(defect.cols());
    r.parallel([&](int index){
        const auto& group=*r.groups[index];
        if(group.patch.constraints.rows()==static_cast<int>(group.patch.discrete_dofs.size()))return;
        auto& columns=p.columns[index];
        std::vector<unsigned char> seen(defect.cols(),0);
        for(int row:group.patch.discrete_dofs)
            for(Impl::Rows::InnerIterator it(p.defect,row);it;++it)
                if(!seen[it.col()]){seen[it.col()]=1;columns.push_back(it.col());}
        std::sort(columns.begin(),columns.end());
    });
    for(int group=0;group<static_cast<int>(p.columns.size());++group){
        const auto& columns=p.columns[group];p.output_rows+=columns.size();
        for(int j=0;j<static_cast<int>(columns.size());++j)
            p.contributions[columns[j]].emplace_back(group,j);
    }
    PhaseTimer::counter("theta_fused_support_entries",p.output_rows,-1);
    p.local_ready=true;
}
KernelDefectOperator::~KernelDefectOperator()=default;

ComplexMatrix KernelDefectOperator::apply_global(const ComplexMatrix& block) const {
    const auto& p=*impl_;auto& r=p.riesz;
    if(block.rows()!=p.defect.cols()||block.cols()<1||!block.allFinite())
        throw std::invalid_argument("parallel defect block dimensions/values invalid");
    if(std::size_t(p.defect.rows())*block.cols()>r.dense_limit)
        throw std::runtime_error("parallel defect RHS resource limit exceeded");
    const int workers=execution_threads(r.threads);
    ComplexMatrix rhs=ComplexMatrix::Zero(p.defect.rows(),block.cols());
    {PhaseTimer timer("theta_global_rhs",-1);
    #pragma omp parallel for schedule(dynamic,256) num_threads(workers)
    for(int row=0;row<p.defect.rows();++row)
        for(Impl::Rows::InnerIterator it(p.defect,row);it;++it)
            for(int j=0;j<block.cols();++j)rhs(row,j)+=it.value()*block(it.col(),j);
    }
    auto values=p.owner.apply_action(rhs);
    rhs.resize(0,0);
    ComplexMatrix result=ComplexMatrix::Zero(block.rows(),block.cols());
    {PhaseTimer timer("theta_global_dual",-1);
    #pragma omp parallel for schedule(dynamic,4) num_threads(workers)
    for(int row=0;row<p.original.cols();++row)
        for(ComplexSparseMatrix::InnerIterator it(p.original,row);it;++it)
            for(int j=0;j<block.cols();++j)result(row,j)+=std::conj(it.value())*values(it.row(),j);
    }
    return result;
}

ComplexMatrix KernelDefectOperator::apply(const ComplexMatrix& block) const {
    const auto& p=*impl_;auto& r=p.riesz;
    if(!p.local_ready)throw std::logic_error("local defect action was not prepared");
    if(block.rows()!=p.defect.cols()||block.cols()<1||!block.allFinite())
        throw std::invalid_argument("kernel defect block dimensions/values invalid");
    if(p.output_rows*std::size_t(block.cols())>r.dense_limit)
        throw std::runtime_error("kernel defect projection resource limit exceeded");
    ++r.calls;r.columns+=block.cols();
    std::vector<ComplexMatrix> projected(r.groups.size());
    {PhaseTimer timer("theta_fused_patches",-1);
    r.parallel([&](int index){
        const auto& columns=p.columns[index];if(columns.empty())return;
        const auto& group=*r.groups[index];const auto& dofs=group.patch.discrete_dofs;
        if(std::size_t(dofs.size())*block.cols()>r.dense_limit)
            throw std::runtime_error("kernel defect local RHS resource limit exceeded");
        ComplexMatrix rhs=ComplexMatrix::Zero(dofs.size(),block.cols());
        // Stream shared CSR rows; never duplicate D once per overlapping patch.
        for(int i=0;i<static_cast<int>(dofs.size());++i)
            for(Impl::Rows::InnerIterator it(p.defect,dofs[i]);it;++it)
                for(int j=0;j<block.cols();++j)rhs(i,j)+=it.value()*block(it.col(),j);
        auto solved=r.solve(index,rhs);
        auto& output=projected[index];output=ComplexMatrix::Zero(columns.size(),block.cols());
        std::vector<int> positions(block.rows(),-1);
        for(int j=0;j<static_cast<int>(columns.size());++j)positions[columns[j]]=j;
        for(int i=0;i<static_cast<int>(dofs.size());++i)
            for(Impl::Rows::InnerIterator it(p.defect,dofs[i]);it;++it){
                const int target=positions[it.col()];
                for(int j=0;j<block.cols();++j)output(target,j)+=std::conj(it.value())*solved.values(i,j);
            }
        // Identical patches were grouped when factors were prepared.
        output*=static_cast<double>(group.nodes.size());
    });}
    ComplexMatrix result=ComplexMatrix::Zero(block.rows(),block.cols());
    {PhaseTimer timer("theta_fused_gather",-1);
    const int workers=execution_threads(r.threads);
    #pragma omp parallel for schedule(static) num_threads(workers)
    for(int row=0;row<result.rows();++row)
        for(const auto& [group,local]:p.contributions[row])
            for(int j=0;j<block.cols();++j)result(row,j)+=projected[group](local,j);
    }
    return result;
}
}
