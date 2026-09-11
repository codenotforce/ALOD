"""Compare a fresh nominal FEM solve with the frozen legacy-kernel output."""
import argparse
import json
import math
import subprocess
from pathlib import Path


def compare(expected, actual, path="root"):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise ValueError(f"{path}: key mismatch")
        return max((compare(value, actual[key], f"{path}.{key}") for key, value in expected.items()), default=0)
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise ValueError(f"{path}: length mismatch")
        return max((compare(x, y, f"{path}[{i}]") for i, (x, y) in enumerate(zip(expected, actual))), default=0)
    if isinstance(expected, (int, float)) and not isinstance(expected, bool):
        if not isinstance(actual, (int, float)) or not math.isfinite(actual):
            raise ValueError(f"{path}: nonfinite or nonnumeric output")
        delta = abs(expected-actual)
        tolerance = 0 if isinstance(expected, int) else 1e-12 + 1e-8*abs(expected)
        if delta > tolerance:
            raise ValueError(f"{path}: delta {delta:.3g} exceeds {tolerance:.3g}")
        return delta
    if expected != actual:
        raise ValueError(f"{path}: value mismatch")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', required=True, type=Path)
    parser.add_argument('--fixture', required=True, type=Path)
    args = parser.parse_args()
    fixture = json.loads(args.fixture.read_text(encoding='utf-8'))
    expected = fixture['output']
    result = subprocess.run([str(args.executable.resolve()), expected['problem'], str(expected['level'])],
                            check=True, text=True, capture_output=True, timeout=90)
    delta = compare(expected, json.loads(result.stdout))
    print(f"{expected['problem']}: legacy FEM comparison passed; max absolute delta={delta:.3g}")


if __name__ == '__main__':
    main()
