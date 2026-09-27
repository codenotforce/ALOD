"""Orchestrate validation, rendering and a self-contained delivery package."""

import gzip
import json
import shutil
from pathlib import Path
from .data import clean, sha
from .selection import METHODS, KAPPAS
from .figures import render_figures, STYLE_PROVENANCE
from .tables import write_tables

FIGURES = ("e1_final_meshes", "e1_error_dof", "e2_final_meshes", "e2_error_dof", "kappa_error_dof")
TABLES = (
    "adaptive-config",
    "e1-terminal-results",
    "e1-family-control",
    "e2-terminal-results",
    "kappa-target-summary",
)


def deliver(runs, output, *, allow_partial=False):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    for name in ("figures", "tables", "data"):
        (output / name).mkdir()
    runs = clean(runs)
    required = (
        [f"{e}-{m}" for e in ("e1", "e2") for m in METHODS]
        + [f"k{k}-{m}" for k in KAPPAS for m in ("ALOD", "AFEM", "LOD")]
        + ["control-" + x for x in ("nominal-alod", "family-alod", "nominal-afem", "family-afem")]
    )
    missing = [name for name in required if name not in runs]
    issues = []
    for name in required:
        if name not in runs:
            continue
        r = runs[name]
        if name in ("e1-ALOD", "e2-ALOD"):
            if len(r["rows"]) < 51:
                issues.append(name + ": fewer than 51 accepted states")
            if {a["sample"] for a in r.get("audits", [])} != set(range(48)):
                issues.append(name + ": incomplete 48-member audit")
        if not r.get("complete") and not (name.startswith("k") and r.get("resource_stopped")):
            issues.append(name + ": incomplete run")
    render_figures(runs, output / "figures", issues)
    metrics, inventory = write_tables(runs, output, issues)
    # Canonical input can regenerate the package without any original directories.
    (output / "data/runs.json.gz").write_bytes(
        gzip.compress(
            json.dumps(dict(schema=1, runs=runs), separators=(",", ":")).encode(),
            compresslevel=4,
            mtime=0,
        )
    )
    report = dict(
        schema=1,
        manuscript="ALOD version4",
        manuscript_sha256=inventory["manuscript_sha256"],
        figures=list(FIGURES),
        tables=list(TABLES),
        complete=not missing and not issues,
        missing_runs=missing,
        issues=issues,
        scientific_notes=[
            "The manuscript family-control text says first N >= 2000, but its supplied published table values are log-budget interpolation at N = 2000. Both estimands are exported separately; this package does not silently reconcile the manuscript."
        ],
        fit_policy="E1: ALOD five completed cycles, AFEM/UFEM eight states, LOD six states. E2: ALOD ten completed cycles, AFEM/UFEM ten displayed states, LOD six displayed states. Historical E2 controls use the archived display limits; native controls end at the first crossing of the final ALOD error.",
        origins={k: r.get("origin") for k, r in runs.items()},
        metrics=metrics,
        plot_style=STYLE_PROVENANCE,
    )
    (output / "tools").mkdir()
    (output / "configs").mkdir()
    project = Path(__file__).resolve().parents[2]
    shutil.copyfile(project / "tools/paper_delivery.py", output / "tools/paper_delivery.py")
    (output / "tools/paper").mkdir()
    for module in (project / "tools/paper").glob("*.py"):
        shutil.copyfile(module, output / "tools/paper" / module.name)
    shutil.copyfile(
        project / "configs/paper_v4_inventory.json", output / "configs/paper_v4_inventory.json"
    )
    shutil.copyfile(project / "requirements-paper.txt", output / "requirements-paper.txt")
    (output / "README.md").write_text(
        "# ALOD version 4 result delivery\n\n"
        + (
            "Complete inventory.\n"
            if report["complete"]
            else "Partial inventory: see delivery.json for missing evidence.\n"
        )
        + "\nHistorical imports are not measurements of the current executable. Tables report actual crossings; resource-stopped endpoints are explicitly labelled. All five figures and five table generators are included. Rebuild from data/runs.json.gz with tools/paper_delivery.py --canonical.\n",
        encoding="utf-8",
    )
    report["files"] = {
        str(p.relative_to(output)).replace("\\", "/"): sha(p)
        for p in sorted(output.rglob("*"))
        if p.is_file()
    }
    (output / "delivery.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if not report["complete"] and not allow_partial:
        raise ValueError(
            "paper evidence incomplete; see delivery.json (use --allow-partial for an explicitly partial package)"
        )
    return report
