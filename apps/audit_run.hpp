#pragma once
#include "alod/batch.hpp"
#include "alod/checkpoint.hpp"
#include <memory>
#include <ostream>
struct AuditWorkerState {
    alod::ExactIntegrationReuse exact_reuse;
    std::unique_ptr<alod::Checkpoint> prepared;
    std::string prepared_path, hierarchy_key;
    std::shared_ptr<void> snapshot_owner;
    std::shared_ptr<const alod::LodHierarchyData> hierarchy;
};
int audit_main(int argc,char** argv,AuditWorkerState* worker=nullptr,std::ostream* output=nullptr,std::ostream* error=nullptr);
int audit_worker_main();
