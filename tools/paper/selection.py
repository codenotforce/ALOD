"""Keep each manuscript's sampling convention explicit and shared with tables."""

METHODS = ("ALOD", "AFEM", "UFEM", "LOD")
KAPPAS = (8, 16, 32, 64, 128)


def completed_cycles(run):
    """E1/E2 omit intermediate reference sweeps with repeated online dimension."""
    return [r for r in run["rows"] if r.get("cycle_complete", True) or r["state"] == 0]


def target_prefix(run, target=0.01):
    """E3 keeps *all* accepted states through the first crossing, including sweeps."""
    rows = run["rows"]
    end = next((i + 1 for i, row in enumerate(rows) if row["E"] <= target), len(rows))
    return rows[:end]


def displayed(exp, method, run, runs):
    rows = run["rows"]
    if run.get("display_state_limit"):
        rows = rows[: run["display_state_limit"]]
    elif exp == "e2" and method != "ALOD" and "e2-ALOD" in runs:
        rows = target_prefix(run, runs["e2-ALOD"]["rows"][-1]["E"])
    return dict(run, rows=rows)


def promotions(run):
    if "promotions" in run:
        return run["promotions"]
    # Older canonical packages did not preserve ell-check events. Inspect every
    # accepted state before mapping a promotion to its next plotted H-step.
    return [
        dict(state_id=after["state"], ell=before["ell"], next_ell=after["ell"])
        for before, after in zip(run["rows"], run["rows"][1:])
        if before.get("ell") and after.get("ell", before["ell"]) > before["ell"]
    ]
