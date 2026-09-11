"""Exercise complete/partial sweeps, fixed controls and same-state E2 promotion."""
import json
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from run_adaptive import DEFAULT, run

exe = Path(sys.argv[1])
with tempfile.TemporaryDirectory() as temp:
    root = Path(temp)
    def check(name, **changes):
        config = dict(DEFAULT, level=1, gap=2, member_ids=[0, 1], training_ids=[0, 1],
                      emit_solution=True, **changes)
        manifest = run(config, root/name, exe, 180)
        events = [json.loads(line) for line in (root/name/"events.jsonl").read_text().splitlines()]
        accepted = [e for e in events if e["kind"] == "accepted"]
        assert all(e["pg_residual"] < 1e-8 for e in accepted)
        return config, manifest, events, accepted
    for m in (1, 2, 3):
        _, _, _, states = check(f"complete-{m}", m_ref=m, reference_theta=[.2]*m, ell_mode="fixed")
        assert len(states) == 1+m and states[-1]["cycle_complete"]
        assert [s["reference_sweep"] for s in states] == [0, *range(1, m+1)]
        for s in states[2:]:
            assert s["coarse_fingerprint"] == states[1]["coarse_fingerprint"]
        if m > 1:
            _, manifest, _, partial = check(f"partial-{m}", m_ref=m, reference_theta=[.2]*m,
                                             ell_mode="fixed", state_limit=m)
            assert manifest["interrupted_cycle"] and not partial[-1]["complete_horizon"]
            assert partial[-1]["next_phase"] == "reference_only"
            assert [s["solution"] for s in partial] == [s["solution"] for s in states[:m]]
    _, _, events, states = check("forced", problem="E2", ell=2, maximum_ell=4,
                                  force_promotions=[0], cycles=0, m_ref=1, reference_theta=[.2])
    checks = [e for e in events if e["kind"] == "ell_check"]
    assert [e["ell"] for e in checks] == [2, 3, 4]
    assert len(states) == 1 and states[0]["ell"] == 4
    assert len({e["mesh_revision"] for e in events}) == 1
    assert checks[1]["warm_transport"] == "coarse_nodal_injection"
    for mode in ("lazy", "every"):
        _, _, events, states = check(f"cadence-{mode}", problem="E2", ell=2, maximum_ell=2,
                                      ell_mode=mode, cycles=2, m_ref=1, reference_theta=[.2])
        assert [e["state_id"] for e in events if e["kind"] == "ell_check"] == ([0, 2] if mode == "lazy" else [0, 1, 2])
        if mode == "lazy":
            assert states[1]["theta"] is None
    for method in ("ALOD", "AFEM"):
        for family in (False, True):
            config = dict(DEFAULT, level=1, gap=2, method=method, cycles=1,
                          ell_mode="fixed", ell=3, maximum_ell=3, audit=True,
                          m_ref=2 if method == "ALOD" else 1,
                          reference_theta=[.3, .2] if method == "ALOD" else [.15],
                          member_ids=list(range(48)), training_ids=list(range(16)) if family else [0])
            name=f"{method}-{family}"
            run(config, root/name, exe, 180)
            states = [json.loads(l) for l in (root/name/"events.jsonl").read_text().splitlines() if json.loads(l)["kind"] == "accepted"]
            assert len(states) == 1+config["m_ref"]
            assert all(s["training_ids"] == config["training_ids"] and len(s["audit"]) == 48 for s in states)
            assert all([a["sample"] for a in s["audit"]] == list(range(48)) for s in states)
print("P3/P4 CLI trajectories and four E1 control contracts passed")
