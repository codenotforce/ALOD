#pragma once
#include "alod/estimator.hpp"
namespace alod {
// D* R_AS D without global fine-grid dense intermediates. Borrows immutable
// factors and CSC defect; both must outlive the operator and all applications.
class KernelDefectOperator {
public:
    KernelDefectOperator(AdditiveKernelRieszContext&,const ComplexSparseMatrix&,bool prepare_local=true);
    KernelDefectOperator(AdditiveKernelRieszContext&,ComplexSparseMatrix&&,bool=true)=delete;
    ~KernelDefectOperator();
    KernelDefectOperator(const KernelDefectOperator&)=delete;
    KernelDefectOperator& operator=(const KernelDefectOperator&)=delete;
    ComplexMatrix apply(const ComplexMatrix&) const;
    // Independent global Riesz/scatter path, with row-owned parallel products.
    ComplexMatrix apply_global(const ComplexMatrix&) const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
