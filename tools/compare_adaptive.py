"""Compare a bounded new trajectory with immutable archived P3/P4 evidence."""
import argparse
import json
import math
from pathlib import Path
import tempfile
from run_adaptive import run


def near(actual, expected, relative=1e-8, absolute=1e-12):
    if isinstance(expected, list):
        assert len(actual) == len(expected)
        for a, b in zip(actual, expected):
            near(a, b, relative, absolute)
    else:
        assert math.isfinite(actual) and abs(actual-expected) <= absolute+relative*abs(expected), (actual, expected)


def compare(fixture, events):
    states = [e for e in events if e["kind"] == "accepted"]
    if "states" in fixture:
        for expected in fixture["states"]:
            s = states[expected["state_id"]]
            for key in ("state_id", "ell", "coarse_free", "reference_free", "coarse_marks"):
                assert s[key] == expected[key], (key, s[key], expected[key])
            near(s["eta"], expected["eta"])
            near(s["audit"][0]["exact_error"], expected["exact_error"])
            near(s["audit"][0]["exact_norm"], expected["exact_norm"])
        checks = [e for e in events if e["kind"] == "ell_check"]
        for old in fixture["checks"]:
            new = next(c for c in checks if c["state_id"] == old["state_id"])
            assert new["ell"] == old["ell"]
            near(new["theta"], old["theta"], absolute=1e-7)
            near(new["eta"], old["eta"])
        return dict(scope="E1 archive equality", compared_states=len(fixture["states"]))
    s = states[0]
    assert s["rank"] == fixture["rank"]
    near(s["targets"], fixture["targets"])
    near(s["theta"], fixture["theta"], absolute=1e-7)
    near([a["exact_norm"] for a in s["audit"]], fixture["exact_norm"])
    # The explicit corrected bootstrap is numerically different. This 0.01%
    # bound is a proximity guard, never an archive-equivalence claim.
    near(s["eta"], fixture["eta"], relative=1e-4)
    errors = [a["exact_error"] for a in s["audit"]]
    near(errors, fixture["exact_error"], relative=1e-4)
    return dict(scope="E2 equivalent targets/Theta and enriched proximity only",
                eta_max_absolute_difference=max(abs(a-b) for a,b in zip(s["eta"], fixture["eta"])),
                exact_error_max_absolute_difference=max(abs(a-b) for a,b in zip(errors, fixture["exact_error"])))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--executable", type=Path, required=True)
    p.add_argument("--fixture", type=Path, required=True)
    args = p.parse_args()
    fixture = json.loads(args.fixture.read_text())
    with tempfile.TemporaryDirectory() as temp:
        out = Path(temp)/"run"
        run(fixture["config"], out, args.executable, 180)
        print(json.dumps(compare(fixture, [json.loads(line) for line in (out/"events.jsonl").read_text().splitlines()])))


if __name__ == "__main__":
    main()
