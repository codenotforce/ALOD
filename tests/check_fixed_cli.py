"""Exercise P2 artifact status, explicit IDs and failure reporting end to end."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools.run_fixed import DEFAULT, write_members, validate, arguments

exe = Path(sys.argv[1]).resolve()
env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    for problem in ("E1", "E2"):
        config = dict(DEFAULT, problem=problem, graded=True, member_ids=[0, 1, 16, 40], training_ids=[0])
        cfg = root / f"{problem}.json"
        cfg.write_text(json.dumps(config))
        output = root / problem
        subprocess.run([sys.executable, str(ROOT / "tools/run_fixed.py"), "--config", str(cfg),
                        "--output", str(output), "--executable", str(exe)], check=True, env=env)
        manifest = json.loads((output / "run.json").read_text())
        state = json.loads((output / "state.json").read_text())
        assert manifest["status"] == "complete" and state["training_ids"] == [0]
        assert state["interpolation"] == "area" and state["riesz_patches"] == "n2"
        assert state["worst_member_id"] == 0 and state["ritz_converged"]
        assert state["ratio_status"] == "finite" and abs(state["ratio"]*state["nominal_eta"]-state["theta"]) < 1e-10
    table = root / "members.txt"
    write_members(table, validate(DEFAULT))
    base = arguments(DEFAULT, table)
    for replacement in ("--theta=nan", "--training-ids=16", "--training-ids=0,0", "--maximum-nodes=1",
                        "--riesz-patches=unknown", "--threads=0"):
        key = replacement.split("=", 1)[0] + "="
        args = [arg for arg in base if not arg.startswith(key)] + [replacement]
        run = subprocess.run([str(exe), *args], capture_output=True, text=True, env=env)
        assert run.returncode and "alod_fixed:" in run.stderr and not run.stdout
    failing = dict(DEFAULT, problem="E2", ritz_iterations=1, ritz_tolerance=1e-30, dense_threshold=0)
    cfg = root / "failure.json"
    cfg.write_text(json.dumps(failing))
    output = root / "failure"
    run = subprocess.run([sys.executable, str(ROOT / "tools/run_fixed.py"), "--config", str(cfg),
                          "--output", str(output), "--executable", str(exe)], capture_output=True, env=env)
    assert run.returncode and json.loads((output / "run.json").read_text())["status"] == "failed"
    assert "did not converge" in (output / "stderr.log").read_text()
print("P2 CLI artifacts, nominal selection and explicit failures passed")
