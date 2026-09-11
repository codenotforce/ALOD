"""Run one fixed-state profile on explicitly selected physical Linux cores.

This bounded P2 probe is not a production scheduler. It requires a fresh output
directory and enforces affinity before starting the numerical process.
"""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import resource
import subprocess
import time
from run_fixed import ROOT, validate, write_members, arguments

PROC = Path(os.sep) / "proc"
SYS_CPU = Path(os.sep) / "sys/devices/system/cpu"


def fields(path):
    return dict(line.split(":", 1) for line in path.read_text().splitlines() if ":" in line)


def kib(value):
    return int(value.split()[0])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--cpus", required=True, help="comma-separated logical CPUs, one per physical core")
    p.add_argument("--reserve-gib", type=float, default=80)
    p.add_argument("--maximum-rss-gib", type=float, default=4)
    p.add_argument("--maximum-seconds", type=int, default=1800)
    a = p.parse_args()
    config = json.loads(a.config.read_text())
    rows = validate(config)
    cpus = [int(c) for c in a.cpus.split(",")]
    if len(set(cpus)) != len(cpus) or len(cpus) != config["threads"] or not set(cpus) <= os.sched_getaffinity(0):
        raise ValueError("CPU set must be allowed, unique and match configured threads")
    physical = []
    for cpu in cpus:
        topology = SYS_CPU / f"cpu{cpu}" / "topology"
        physical.append((int((topology / "physical_package_id").read_text()), int((topology / "core_id").read_text())))
    if len(set(physical)) != len(physical):
        raise ValueError("selected CPUs share physical cores")
    mem = fields(PROC / "meminfo")
    reserve = max(a.reserve_gib * 2**20, .2 * kib(mem["MemTotal"]))
    if not math.isfinite(a.maximum_rss_gib) or not math.isfinite(a.reserve_gib) or a.maximum_rss_gib <= 0 or a.maximum_seconds <= 0 or a.reserve_gib < 0:
        raise ValueError("invalid profile resource limits")
    if kib(mem["MemAvailable"]) - 1.5 * a.maximum_rss_gib * 2**20 < reserve:
        raise RuntimeError("insufficient memory headroom for the conservative 1.5x job allowance")
    a.output.mkdir(parents=True, exist_ok=False)
    write_members(a.output / "members.txt", rows)
    report = dict(config=config, cpus=cpus, physical_cores=physical, reserve_kib=reserve,
                  initial_available_kib=kib(mem["MemAvailable"]), initial_swap_free_kib=kib(mem["SwapFree"]),
                  peak_group_rss_kib=0, peak_group_swap_kib=0, minimum_available_kib=kib(mem["MemAvailable"]),
                  maximum_rss_gib=a.maximum_rss_gib, maximum_seconds=a.maximum_seconds,
                  worker_affinity={}, worker_cpu_ticks={}, status="running")
    os.sched_setaffinity(0, cpus)
    env = {**os.environ, "OMP_NUM_THREADS": str(config["threads"]), "OPENBLAS_NUM_THREADS": "1",
           "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1", "OMP_PROC_BIND": "false"}
    env.pop("OMP_PLACES", None)
    start = time.monotonic()
    child = None
    try:
        with (a.output / "state.json").open("w") as out, (a.output / "stderr.log").open("w") as err:
            child = subprocess.Popen([str(ROOT / "build/alod_fixed"), *arguments(config, (a.output / "members.txt").resolve())],
                                     env=env, stdout=out, stderr=err, start_new_session=True)
            while True:
                rss = swap = 0
                for process in PROC.iterdir():
                    if not process.name.isdigit():
                        continue
                    try:
                        stat = (process / "stat").read_text().rsplit(")", 1)[1].split()
                        if int(stat[2]) != child.pid:
                            continue
                        status = fields(process / "status")
                        rss += kib(status.get("VmRSS", "0"))
                        swap += kib(status.get("VmSwap", "0"))
                        for worker in (process / "task").iterdir():
                            affinity = fields(worker / "status")["Cpus_allowed_list"].strip()
                            report["worker_affinity"][worker.name] = affinity
                            stat = (worker / "stat").read_text().rsplit(")", 1)[1].split()
                            report["worker_cpu_ticks"][worker.name] = int(stat[11]) + int(stat[12])
                            if not os.sched_getaffinity(int(worker.name)) <= set(cpus):
                                raise RuntimeError("worker affinity escaped the selected physical cores")
                    except (FileNotFoundError, ProcessLookupError):
                        continue
                mem = fields(PROC / "meminfo")
                report["peak_group_rss_kib"] = max(report["peak_group_rss_kib"], rss)
                report["peak_group_swap_kib"] = max(report["peak_group_swap_kib"], swap)
                report["minimum_available_kib"] = min(report["minimum_available_kib"], kib(mem["MemAvailable"]))
                if rss > a.maximum_rss_gib * 2**20 or swap > 256*1024 or kib(mem["MemAvailable"]) < reserve:
                    raise RuntimeError("profile exceeded RSS, swap or available-memory limit")
                if time.monotonic() - start > a.maximum_seconds:
                    raise RuntimeError("profile exceeded wall-clock limit")
                if child.poll() is not None:
                    break
                time.sleep(.1)
            report["returncode"] = child.returncode
            if child.returncode:
                raise RuntimeError("fixed-state program failed; see stderr.log")
            report["participating_workers"] = sum(ticks > 0 for ticks in report["worker_cpu_ticks"].values())
            if report["participating_workers"] < config["threads"]:
                raise RuntimeError("profile did not observe all configured workers performing CPU work")
            report["status"] = "complete"
    except BaseException as error:
        report["status"] = "failed"
        report["reason"] = str(error)
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
        raise
    finally:
        report["wall_seconds"] = time.monotonic() - start
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        report["child_kernel_peak_rss_kib"] = usage.ru_maxrss
        report["child_user_seconds"] = usage.ru_utime
        report["child_system_seconds"] = usage.ru_stime
        (a.output / "profile.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
