"""Compare a P2 fixed state with frozen archived-kernel evidence."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
from run_fixed import validate, write_members, arguments


def compare(expected, actual, path=""):
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"shape mismatch at {path}")
        for i, (x, y) in enumerate(zip(expected, actual)):
            compare(x, y, f"{path}[{i}]")
    elif type(expected) is int:
        if expected != actual:
            raise AssertionError(f"integer mismatch at {path}: {expected} != {actual}")
    else:
        # Theta is a square root; an exactly vanishing defect has a roundoff
        # floor of order sqrt(epsilon) in the archived dense Gram solver.
        atol = 1e-7 if path == "theta" else 1e-10
        if not math.isfinite(actual) or abs(expected-actual) > atol + 1e-8 * max(abs(expected), abs(actual)):
            raise AssertionError(f"numerical mismatch at {path}: {expected} != {actual}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--executable", type=Path, required=True)
    p.add_argument("--fixture", type=Path, required=True)
    p.add_argument("--threads", type=int)
    a = p.parse_args()
    fixture = json.loads(a.fixture.read_text())
    config = dict(fixture["config"])
    if a.threads:
        config["threads"] = a.threads
    env = {**os.environ, "OMP_NUM_THREADS": str(config["threads"]), "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    with tempfile.TemporaryDirectory() as tmp:
        table = Path(tmp) / "members.txt"
        write_members(table, validate(config))
        result = subprocess.run([str(a.executable.resolve()), *arguments(config, table)],
                                env=env, capture_output=True, text=True, timeout=180, check=True)
        actual = json.loads(result.stdout)
    for key, expected in fixture["output"].items():
        compare(expected, actual[key], key)
    if not actual["ritz_converged"] or actual["ritz_residual"] > config["ritz_tolerance"]:
        raise AssertionError("localization did not converge")
    print(f"{a.fixture.stem}: archived P2 comparison passed")


if __name__ == "__main__":
    main()
