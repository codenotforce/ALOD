# P0 acceptance evidence

The original E1 production deployment was retrieved read-only from its archived host. Its runner SHA-256 is `5b4867a4f1a0e7dd6b88bda92e597c56a9645ed6581baa37c62684b6bb26e807`, workflow include SHA-256 is `186f89fbac07e39106d304198ac246ad14bb39cd9dc92aea776991beca10eda8`, executable SHA-256 is `f9a92ccd2807a3d6c4e05914b001052f0deb0ca7c627fd9bee0ec0d2c4e12d4b`, and core archive SHA-256 is `6488b977a316992f48a573c63b80e194e42bc2cbf2822efb2f2e94084e741941`. These match the archived production checksum list. Full original/replay source hashes are recorded separately for E1 and E2.

Both source snapshots were rebuilt from source in isolated local directories with no reused legacy libraries. The old production host guards require more available memory than the local WSL machine has. Explicit resource-only patches set the process RSS limit to 4 GiB and host reserve to 2 GiB; formulas and acceptance decisions are unchanged. Configurations also bound horizon, workload, thread count and observation indices. E1 uses the historically validated H2/h4 27-state smoke configuration; E2 preserves H6/h10 and the production build-hash switches over two accepted states. OpenMP=2 and BLAS=1.

E1: 170 archived family audit rows compared at every configured audit point, 120 coarse marking entries compared exactly, and all 27 accepted states replayed. Maximum compared numerical difference: 3.06e-10, within 1e-12 + 1e-8*abs(expected). E2: 100 accepted rows over two states, maximum difference 4.22e-15. The E2 shortened endpoint produces extra regional diagnostic partitions; the comparison report identifies these explicitly. Training-only eta fields remain missing for audit-only members, as in the archive. Both replays stopped at their configured horizon. Observed peaks were approximately 153 MiB (E1) and 642 MiB (E2), not performance benchmarks.

`*-build-commands.json` records each translation unit and final link command with symbolic roots. `*-link-closure.json` lists archive members actually pulled in by GNU ld, unlike the conservative include inventories in the parent directory. The source hash inventory preserves headers and included workflow sources as well. `<CXX>`, `<e1_SOURCE>`, `<e2_SOURCE>`, corresponding build roots and `<DEPENDENCY>` are operator-supplied locations, not literal filesystem paths.

`*-branch-inventory.json` records source locations of build-hash token tests and environment reads. Enabled flags preserve archived production workflows; controls are separately documented in `docs/experiments.md`. E1 and E2 require different snapshots, including the production-only selected-value kernel API. Do not combine the main worktree header with a production library. The AOT optimization snapshot is recorded as a later integration input, not falsely substituted for these accepted-source replays.

Dependency versions are in `dependencies.json`. Archived compiler/CMake evidence is available; package versions were additionally read from the same host on the observation date. No historical package lock exists, so bitwise toolchain reconstruction is not claimed. Numerical replay is the acceptance evidence.

The manuscript inventory, frozen samples, four control audits, provenance hashes, effective configurations, old numerical replays and actual link closure complete P0's source-discovery gate. The area-average discrepancy is an explicit scientific finding with a reproducer, not a claim of adaptive manuscript equivalence. Resolve it before P2 acceptance.

To reproduce the archived smoke runs on Linux/WSL, supply the original snapshot roots explicitly:

```sh
python3 tools/replay_legacy.py --problem E1 --source "$E1_SOURCE_ROOT" --work _local/replay-e1
python3 tools/replay_legacy.py --problem E2 --source "$E2_SOURCE_ROOT" --work _local/replay-e2
```

The work directories must be new. The tool copies the sources, applies only the recorded resource changes, builds the old executable, runs the bounded configuration, and compares the results with the frozen evidence. For baseline probes, use `tools/build_baseline_oracles.py --help`; its `--native` flag must match the supplied legacy library's architecture flags. Normal ALOD builds do not invoke any of these developer tools.
