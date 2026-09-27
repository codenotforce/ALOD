"""Scientific table calculations, independent of Matplotlib styling."""

import csv
import json
import math
from pathlib import Path
import numpy as np
from .data import read
from .selection import METHODS, KAPPAS, displayed, completed_cycles


def write_table(output, name, head, rows):
    with (output / "tables" / f"{name}.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(head)
        w.writerows(rows)

    def tex(x):
        return str(x).replace("_", r"\_").replace("%", r"\%")

    lines = [
        r"\begin{tabular}{" + "l" * len(head) + "}",
        " & ".join(map(tex, head)) + r" \\ \hline",
    ]
    lines += [" & ".join(map(tex, row)) + r" \\" for row in rows]
    lines += [r"\end{tabular}"]
    (output / "tables" / f"{name}.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def terminal_rows(exp, runs):
    rows = []
    for method in METHODS:
        run = runs.get(f"{exp}-{method}")
        if not run:
            continue
        run = displayed(exp, method, run, runs)
        rr = run["rows"]
        last = rr[-1]
        window = (
            6
            if method == "LOD"
            else 5 if exp == "e1" and method == "ALOD" else 8 if exp == "e1" else 10
        )
        tail = completed_cycles(run)[-window:]
        p = r2 = None
        if len(tail) >= 3 and len({r["N_on"] for r in tail}) >= 3:
            x = np.log([r["N_on"] for r in tail])
            y = np.log([r["E"] for r in tail])
            a, b = np.polyfit(x, y, 1)
            p = -float(a)
            r2 = (
                float(1 - np.sum((y - a * x - b) ** 2) / np.sum((y - y.mean()) ** 2))
                if np.var(y) > 0
                else None
            )
        crossings = [
            next((r["N_on"] for r in rr if r["E"] <= t), "not reached") for t in (0.02, 0.01)
        ]
        rows.append([method, *crossings, last["N_on"], last["E"], p, r2])
    return rows


def family_crossing_rows(runs, issues):
    rows = []
    for arm in ("nominal-alod", "family-alod", "nominal-afem", "family-afem"):
        run = runs.get("control-" + arm)
        r = next((r for r in run["rows"] if r["N_on"] >= 2000), None) if run else None
        if not r or r.get("test_count") != 24:
            issues.append(
                "control-" + arm + ": missing 2000-DOF crossing or 24-member held-out audit"
            )
            continue
        rows.append(
            [
                arm,
                r["state"],
                r["N_on"],
                *[100 * r[k] for k in ("E", "test_p50", "test_p90", "test_max")],
            ]
        )
    return rows


def family_budget_rows(runs, issues):
    budget_rows = []
    for arm in ("nominal-alod", "family-alod", "nominal-afem", "family-afem"):
        run = runs.get("control-" + arm)
        if not run:
            continue
        stats = run.get("budget_statistics")
        if not stats:
            rr = completed_cycles(run)
            hi = next((i for i, r in enumerate(rr) if r["N_on"] >= 2000), None)
            if hi is not None and hi > 0:
                lo, up = rr[hi - 1], rr[hi]
                w = math.log(2000 / lo["N_on"]) / math.log(up["N_on"] / lo["N_on"])
                lower = {
                    a["sample"]: a for a in run.get("audits", []) if a["state_id"] == lo["state"]
                }
                upper = {
                    a["sample"]: a for a in run.get("audits", []) if a["state_id"] == up["state"]
                }
                if lower.keys() == upper.keys() and 0 in lower:
                    errors = {
                        m: math.exp(
                            (1 - w) * math.log(float(lower[m]["E"]))
                            + w * math.log(float(upper[m]["E"]))
                        )
                        for m in lower
                    }
                    held = [v for m, v in errors.items() if lower[m].get("role") == "test"]
                    if len(held) == 24:
                        stats = dict(
                            nominal=errors[0],
                            p50=float(np.quantile(held, 0.5)),
                            p90=float(np.quantile(held, 0.9)),
                            maximum=max(held),
                        )
        if stats:
            budget_rows.append(
                [arm, 2000, *[100 * stats[k] for k in ("nominal", "p50", "p90", "maximum")]]
            )
        else:
            issues.append("control-" + arm + ": missing member-level interpolation evidence")
    return budget_rows


def kappa_rows(runs, issues):
    rows = []
    for k in KAPPAS:
        for method in ("ALOD", "AFEM", "LOD"):
            run = runs.get(f"k{k}-{method}")
            if not run:
                continue
            crossing = next((r for r in run["rows"] if r["E"] <= 0.01), None)
            r = crossing or run["rows"][-1]
            status = (
                "target reached"
                if crossing
                else "resource stopped" if run.get("resource_stopped") else "target not reached"
            )
            if not crossing and not run.get("resource_stopped"):
                issues.append(f"k{k}-{method}: target not reached without resource-stop evidence")
            rows.append([k, method, r["state"], r["N_on"], 100 * r["E"], status])
    return rows


def write_tables(runs, output, issues):
    metrics = {}
    for exp in ("e1", "e2"):
        rows = terminal_rows(exp, runs)
        write_table(
            output,
            exp + "-terminal-results",
            [
                "method",
                "first N at 2%",
                "first N at 1%",
                "endpoint N",
                "endpoint E",
                "tail p",
                "tail R2",
            ],
            rows,
        )
        metrics[exp] = rows
    write_table(
        output,
        "e1-family-first-crossing",
        ["method", "state", "N", "nominal %", "test median %", "test p90 %", "test max %"],
        family_crossing_rows(runs, issues),
    )
    write_table(
        output,
        "e1-family-control",
        ["method", "interpolated N", "nominal %", "test median %", "test p90 %", "test max %"],
        family_budget_rows(runs, issues),
    )
    write_table(
        output,
        "kappa-target-summary",
        ["k", "method", "state", "N", "error %", "status"],
        kappa_rows(runs, issues),
    )
    settings = [
        [name, key, json.dumps(value)]
        for name in ("e1-ALOD", "e2-ALOD", *[f"k{k}-ALOD" for k in KAPPAS])
        for key, value in runs.get(name, {}).get("config", {}).items()
    ]
    inventory = read(Path(__file__).resolve().parents[2] / "configs/paper_v4_inventory.json")
    write_table(output, "adaptive-config", ["setting", "E1", "E2", "E3"], inventory["settings"])
    write_table(output, "observed-run-config", ["run", "setting", "value"], settings)
    return metrics, inventory
