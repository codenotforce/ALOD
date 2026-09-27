"""Explicit native/historical adapters and portable scientific-data validation."""

import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET


def read(path):
    path = Path(path)
    return json.loads(
        gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    )


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def adaptive(data, *, complete=True):
    states = data["states"]
    audits = data["audits"]
    seen = {}
    for row in audits:
        key = (row["state_id"], row["sample"])
        if key in seen:
            raise ValueError("duplicate state/member audit")
        # Validate every member, including held-out samples. A good nominal
        # solve must not mask a failed or non-finite family audit.
        residual = float(row.get("PG_residual", 0))
        if not math.isfinite(residual) or not 0 <= residual <= 1e-8:
            raise ValueError("audit residual gate failed")
        if not math.isfinite(float(row["E"])) or float(row["E"]) <= 0:
            raise ValueError("invalid member error")
        seen[key] = row
    config = data.get("config", {})
    expected = config.get("member_ids", list(range(48)))
    curves = []
    for state in states:
        sid = state["state_id"]
        if {m for s, m in seen if s == sid} != set(expected):
            raise ValueError("incomplete frozen audit coverage")
        audit = seen[sid, 0]
        row = dict(
            state=sid,
            N_on=state["coarse_free"] + state.get("rank", 0),
            E=float(audit["E"]),
            ell=state.get("ell"),
            cycle_complete=state.get("cycle_complete", True),
            rank=state.get("rank", 0),
            eta_relative=state["eta"][0] / float(audit["exact_norm"]),
        )
        tests = [
            float(a["E"]) for (s, m), a in seen.items() if s == sid and a.get("role") == "test"
        ]
        if tests:
            import numpy as np

            row.update(
                test_count=len(tests),
                test_p50=float(np.quantile(tests, 0.5)),
                test_p90=float(np.quantile(tests, 0.9)),
                test_max=max(tests),
            )
        curves.append(row)
    return dict(
        rows=curves,
        config=config,
        complete=complete,
        promotions=[
            {k: c[k] for k in ("state_id", "ell", "next_ell")}
            for c in data.get("ell_checks", [])
            if c["action"] == "promote"
        ],
        audits=[
            {k: a.get(k) for k in ("state_id", "sample", "role", "E", "F", "G", "exact_norm")}
            for a in audits
        ],
    )


def native(folder, executable=None):
    folder = Path(folder)
    manifest = read(folder / "run.json")
    if (folder / "events.jsonl").exists():
        events = [json.loads(line) for line in (folder / "events.jsonl").read_text().splitlines()]
        states = [r for r in events if r["kind"] == "accepted"]
        out = adaptive(
            dict(
                states=states,
                audits=[a for r in states for a in r.get("audit", [])],
                ell_checks=[e for e in events if e["kind"] == "ell_check"],
                config=manifest["config"],
            ),
            complete=bool(manifest.get("audit_complete") and manifest.get("solver_completed")),
        )
    else:
        rows = [json.loads(line) for line in (folder / "states.jsonl").read_text().splitlines()]
        out = dict(
            rows=[
                dict(
                    state=r["state"],
                    N_on=r["free_dof"],
                    E=r["relative_energy_error"],
                    ell=r.get("ell"),
                    cycle_complete=True,
                )
                for r in rows
            ],
            config=manifest["config"],
            complete=manifest["status"] == "completed",
        )
    out.update(
        origin="native_run",
        source_sha256=sha(folder / "run.json"),
        source_data_sha256=sha(
            folder / ("events.jsonl" if (folder / "events.jsonl").exists() else "states.jsonl")
        ),
        executable_sha256=manifest.get("executable_sha256"),
    )
    # Audit folders can contain intermediate or deepened meshes. Their directory
    # order does not identify the terminal production mesh.
    if (folder / "mesh.json").exists():
        out["mesh"] = read(folder / "mesh.json")
    elif executable is not None and (folder / "checkpoints/latest").exists():
        import subprocess, tempfile
        from checkpoint_io import resolve_checkpoint

        with tempfile.TemporaryDirectory() as temp:
            file = Path(temp) / "mesh.json"
            subprocess.run(
                [
                    str(Path(executable).resolve()),
                    "checkpoint-mesh",
                    str(resolve_checkpoint(folder / "checkpoints/latest").resolve()),
                    str(file),
                ],
                check=True,
            )
            out["mesh"] = read(file)
    return out


def mesh(path):
    path = Path(path)
    if ".vtu" not in path.name:
        x = read(path)
        return dict(nodes=x["nodes"], elements=x.get("elements", x.get("triangles")))
    root = ET.fromstring(
        gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    )
    import numpy as np

    points = root.find(".//Points/DataArray")
    if points.get("format", "ascii") != "ascii":
        raise ValueError("VTU importer requires explicit ASCII arrays")
    nodes = np.fromstring(points.text, sep=" ").reshape(-1, 3)[:, :2].tolist()
    cells = {
        a.get("Name"): np.fromstring(a.text, sep=" ", dtype=int)
        for a in root.findall(".//Cells/DataArray")
    }
    if any(cells["types"] != 5):
        raise ValueError("VTU importer requires triangular cells")
    return dict(nodes=nodes, elements=cells["connectivity"].reshape(-1, 3).tolist())


def import_archives(e1, e2, e3):
    """Named adapter for the supplied q03/kernel-lift/k8--128 exported schemas."""
    e1, e2, e3 = map(Path, (e1, e2, e3))
    runs = {}
    for exp, source in [("e1", e1 / "data"), ("e2", e2 / "data")]:
        file = source / f"{exp}-production.json"
        data = read(file)
        runs[f"{exp}-ALOD"] = adaptive(data)
        runs[f"{exp}-ALOD"].update(
            origin="historical_v4_export",
            source_sha256=sha(file),
            mesh={
                key: mesh(source / f"{exp}-final-{suffix}.json.gz")
                for key, suffix in [("coarse", "coarse"), ("reference", "reference")]
            },
        )
        baselines = e1 / "data" / f"{exp}-baselines.json"
        for method, record in read(baselines)["runs"].items():
            name = "LOD" if method == "SLOD" else method
            runs[f"{exp}-{name}"] = dict(
                display_state_limit=(
                    {"AFEM": 77, "UFEM": 13, "LOD": 10}[name]
                    if exp == "e2"
                    else len(record["trajectory"])
                ),
                rows=record["trajectory"],
                config={},
                complete=True,
                origin="historical_v4_export",
                source_sha256=sha(baselines),
            )
        runs[f"{exp}-AFEM"]["mesh"] = {"coarse": mesh(e1 / "data" / f"{exp}-afem-final.vtu.gz")}
    controls = e1 / "data/e1-controls.json"
    for name, data in read(controls).items():
        runs["control-" + name] = adaptive(data)
        runs["control-" + name].update(origin="historical_v4_export", source_sha256=sha(controls))
    control_csv = e1 / "controls/data/afem/state_statistics.csv"
    with control_csv.open(newline="", encoding="utf-8-sig") as f:
        statistics = list(csv.DictReader(f))
    for arm in ("nominal-afem", "family-afem"):
        rows = [
            dict(
                state=int(r["state"]),
                N_on=int(r["N_online"]),
                E=float(r["nominal_E"]),
                test_count=24,
                test_p50=float(r["test_E_p50"]),
                test_p90=float(r["test_E_p90"]),
                test_max=float(r["test_E_max"]),
                cycle_complete=True,
            )
            for r in statistics
            if r["arm"] == arm
        ]
        runs["control-" + arm] = dict(
            rows=rows,
            config={},
            complete=True,
            origin="historical_summary_export",
            source_sha256=sha(control_csv),
        )
    summary = read(e1 / "controls/analysis.json")["e1"]["budgets"]
    for arm, budgets in summary.items():
        runs["control-" + arm]["budget_statistics"] = budgets["2000"]
    for r in read(e3 / "data.json")["runs"]:
        name = "LOD" if r["method"] == "SLOD" else r["method"]
        runs[f'k{r["k"]}-{name}'] = dict(
            rows=[
                dict(
                    state=x.get("state_id", x.get("state")),
                    N_on=x.get("coarse_free", x.get("free_dof")),
                    E=x["relative_energy_error"],
                    ell=x.get("ell"),
                    cycle_complete=x.get("cycle_complete", True),
                )
                for x in r["rows"]
            ],
            config=r["config"],
            complete=r["status"] in ("complete", "target_reached"),
            resource_stopped=r.get("guard_reason") is not None,
            origin="historical_v4_export",
            source_sha256=sha(e3 / "data.json"),
        )
    return runs


def clean(runs):
    """Whitelist portable scientific data; original archives may contain host paths."""
    allowed = (
        "problem",
        "method",
        "wavenumber",
        "level",
        "gap",
        "initial_level",
        "ell",
        "maximum_ell",
        "ell_threshold",
        "ell_absolute_threshold",
        "ritz_tolerance",
        "enrichment_tests",
        "theta",
        "reference_theta",
        "m_ref",
        "training_ids",
        "member_ids",
    )
    result = {}
    for name, run in runs.items():
        if not name or any(
            c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_"
            for c in name
        ):
            raise ValueError("unsafe run name")
        rows = []
        for index, row in enumerate(run["rows"]):
            r = {
                k: row[k]
                for k in (
                    "state",
                    "N_on",
                    "E",
                    "eta_relative",
                    "ell",
                    "rank",
                    "cycle_complete",
                    "test_count",
                    "test_p50",
                    "test_p90",
                    "test_max",
                )
                if k in row
            }
            r.setdefault("state", index)
            if (
                type(r["N_on"]) not in (int, float)
                or not math.isfinite(r["N_on"])
                or r["N_on"] <= 0
            ):
                raise ValueError("invalid online dimension")
            if type(r["E"]) not in (int, float) or not math.isfinite(r["E"]) or r["E"] <= 0:
                raise ValueError("invalid error")
            rows.append(r)
        if not rows or [r["state"] for r in rows] != list(range(len(rows))):
            raise ValueError("trajectory has gaps or duplicates")
        result[name] = {
            k: run[k]
            for k in (
                "origin",
                "source_sha256",
                "source_data_sha256",
                "executable_sha256",
                "complete",
                "resource_stopped",
                "mesh",
                "audits",
                "display_state_limit",
                "budget_statistics",
            )
            if k in run
        }
        if "promotions" in run:
            result[name]["promotions"] = [
                {k: int(c[k]) for k in ("state_id", "ell", "next_ell")} for c in run["promotions"]
            ]
        result[name].update(
            rows=rows, config={k: v for k, v in run.get("config", {}).items() if k in allowed}
        )
    return result
