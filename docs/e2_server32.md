# E2 server deployment (2026-09-22)

**Historical deployment record.** All configuration, source-match, cache-default, process and progress statements below refer to the September 22 deployment only. Current E2 uses 28 cycles and later cache/runtime defaults; see [the production protocol](production_campaign_20260927.md). This record must not be used as live server status.

The isolated deployment directory is `e2-32-20260922` under the server project root. The tmux session is `alod-e2-32-20260922`. Existing experiments are preserved.

`configs/adaptive/e2_server32.json` retains every mathematical setting from `e2_main.json` and sets `threads` to 32. The run uses 25 cycles, up to 51 accepted states, 16 training members, all 48 audit members, radius 0.6, kernel lifting, and lazy localization checks. One asynchronous audit worker uses 32 threads while adaptation independently uses 32 threads.

The deployed executable matched the numerical sources at deployment time and passed all 43 tests. Its SHA-256 is `54000ff7c8a699a98b9ffeb699a903d5e38b7177dc7c8052f4f8a44e1293ee48`. An unused legacy timing header in the build source directory was confirmed absent from the executable dependency files. The source archive and executable checksums are recorded in `deployment.sha256` on the server.

Effective defaults include source moments with a 768 MiB budget, incremental enrichment products, audit restoration of the LOD basis and reduced operator, and shared checkpoint geometry. The optional local corrector cache is disabled (`ALOD_PATCH_CACHE_BYTES=0`) because earlier benchmarks did not show a benefit. `ALOD_REFERENCE_EXECUTION` is unset. OpenMP thread limits and affinity overrides are cleared; OpenMP dynamic teams are disabled.

The supervisor runs `tools/run_adaptive.py` with `--audit-workers 1 --audit-threads 32 --timeout 604800`. Output is in `results/e2-main`. Adaptive phase timings are in `results/e2-main/timings.jsonl`; audit phase timings are in its `audits` subdirectories. Resource samples are in `results/e2-resources.jsonl`, the supervisor log is `results/e2-supervisor.log`, and final process status is written to `results/e2-exit-code`.

Initial verification observed 25 accepted states, 17 completed audits, one running audit and seven pending audits. Both numerical processes had 32 threads, unrestricted CPU eligibility across cores 0-127, and zero swap usage. Completed audits reused the saved LOD basis. This is a deployment observation, not a claim that the experiment has completed.
