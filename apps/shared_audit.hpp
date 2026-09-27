#pragma once
#include "alod/checkpoint.hpp"
#include <functional>
// Linux in-process transport; other platforms retain standalone native workers.
int with_shared_audits(const std::function<int()>& adaptive);
void publish_audit_snapshot(const std::filesystem::path&,alod::Checkpoint&&,
    const alod::LodSpace&,const alod::CheckpointGeometryView&);
