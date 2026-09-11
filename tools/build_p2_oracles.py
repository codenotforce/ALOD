"""Developer-only creation of new P2 fixtures from the archived numerical library.

Existing fixture files are never replaced. Runtime paths are not stored in fixtures.
"""
import argparse
import hashlib
import json
import os
import shlex
from pathlib import Path
import subprocess
from run_fixed import ROOT, DEFAULT, validate, write_members, arguments


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--legacy-build", type=Path, required=True)
    p.add_argument("--work", type=Path, default=ROOT / "_local/p2-oracle")
    p.add_argument("--fixtures", type=Path, default=ROOT / "tests/fixtures/p2")
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    a.fixtures.mkdir(parents=True, exist_ok=True)
    names = [f"{problem}_{case}" for problem in ("E1", "E2") for case in ("uniform", "graded", "ritz")]
    if any((a.fixtures / (n + ".json")).exists() for n in names):
        raise ValueError("refusing to replace frozen P2 fixtures")
    executable = a.work / "legacy_probe"
    eigen_flags = shlex.split(subprocess.check_output(["pkg-config", "--cflags", "eigen3"], text=True))
    subprocess.run(["c++", "-O3", "-DNDEBUG", "-std=c++20", "-I" + str(a.source / "include"),
                    "-I" + str(ROOT / "include"), *eigen_flags,
                    str(ROOT / "tests/legacy_p2_probe.cpp"), str(ROOT / "src/problems/problems.cpp"),
                    str(a.legacy_build / "liblod2d_core.a"), "-fopenmp", "-lumfpack", "-o", str(executable)], check=True)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    provenance = dict(library_sha256=sha(a.legacy_build / "liblod2d_core.a"),
                      probe_sha256=sha(ROOT / "tests/legacy_p2_probe.cpp"),
                      model_sha256=sha(a.source / "src/helmholtz/model.cpp"),
                      riesz_sha256=sha(a.source / "src/helmholtz/adaptive/kernel_residual.cpp"),
                      additive_sha256=sha(a.source / "src/helmholtz/adaptive/additive_kernel_riesz.inc"),
                      certificates_sha256=sha(a.source / "src/helmholtz/adaptive/certificates.cpp"),
                      architecture="portable GCC flags matching ALOD; no native architecture option",
                      scope="Complete archived model, Riesz and spectrum kernels; shared input parsing and serialization only")
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    for name in names:
        problem, case = name.split("_")
        config = dict(DEFAULT, problem=problem, interpolation="arithmetic", riesz_patches="archive",
                      graded=case == "graded", level=4 if case == "ritz" else 2,
                      dense_threshold=0 if case == "ritz" else 64)
        table = a.work / f"{problem}.txt"
        write_members(table, validate(config))
        result = subprocess.run([str(executable.resolve()), *arguments(config, table.resolve())],
                                env=env, capture_output=True, text=True, check=True, timeout=180)
        output = json.loads(result.stdout)
        # The archived controller aggregates in member order and supplements the
        # worst member using the same stable comparison and 1e-14 mass tolerance.
        ids = sorted(config["training_ids"])
        columns = [config["member_ids"].index(i) for i in ids]
        scales = output["all_scales"]
        masses = output["element_squared"]
        mean = [sum(row[j] / (len(ids)*scales[j]*scales[j]) for j in columns) for row in masses]
        worst = max(columns, key=lambda j: output["eta"][j] / scales[j])
        mass = [row[worst] / (scales[worst]*scales[worst]) for row in masses]
        marked, covered = [], 0.0
        for i in sorted(range(len(mean)), key=lambda i: (-mean[i], i)):
            marked.append(i)
            covered += mean[i]
            if covered >= config["theta"] * sum(mean):
                break
        covered, total = sum(mass[i] for i in marked), sum(mass)
        for i in sorted((i for i in range(len(mean)) if i not in marked), key=lambda i: (-mass[i], i)):
            if covered + 1e-14 * max(1, total) >= config["theta"] * total:
                break
            marked.append(i)
            covered += mass[i]
        output.update(marked_elements=sorted(marked), worst_member_id=config["member_ids"][worst],
                      frozen_scales=[scales[j] for j in columns])
        del output["all_scales"]
        (a.fixtures / f"{name}.json").write_text(json.dumps(dict(provenance=provenance, config=config, output=output), separators=(",", ":")) + "\n")
        print(f"{name}: archived solution, indicators and Theta captured", flush=True)


if __name__ == "__main__":
    main()
