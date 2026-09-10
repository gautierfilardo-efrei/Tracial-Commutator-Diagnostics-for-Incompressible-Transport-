"""Compare a fresh run with the packaged numerical reference, without fitting."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
RTOL = 1e-9
ATOL = 1e-10


def verify(generated, reference=ROOT):
    generated, reference = Path(generated).resolve(), Path(reference).resolve()
    if generated == reference:
        raise ValueError("Generate into a separate directory before verification")
    reference_results = reference / "results"
    generated_results = generated / "results"
    source_metadata = json.loads((generated_results / "environment.json").read_text())
    for name, expected in source_metadata["sha256"].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise AssertionError(f"Generated metadata does not match source: {name}")
    checked = []
    for extension in ("*.csv", "*.npy"):
        expected_names = {p.name for p in reference_results.glob(extension)}
        actual_names = {p.name for p in generated_results.glob(extension)}
        if not expected_names or expected_names != actual_names:
            raise AssertionError(f"Missing or unexpected {extension} files")
        for name in sorted(expected_names):
            left, right = reference_results / name, generated_results / name
            if left.suffix == ".npy":
                np.testing.assert_allclose(np.load(right, allow_pickle=False),
                                           np.load(left, allow_pickle=False),
                                           rtol=RTOL, atol=ATOL, equal_nan=False,
                                           err_msg=name)
            else:
                with left.open(newline="") as stream:
                    reader = csv.DictReader(stream)
                    expected_fields, expected_rows = reader.fieldnames, list(reader)
                with right.open(newline="") as stream:
                    reader = csv.DictReader(stream)
                    fields, rows = reader.fieldnames, list(reader)
                if fields != expected_fields or len(rows) != len(expected_rows):
                    raise AssertionError(f"Different table structure: {name}")
                for index, (actual, expected) in enumerate(zip(rows, expected_rows)):
                    for field in fields:
                        a, b = actual[field], expected[field]
                        try:
                            a_number, b_number = float(a), float(b)
                        except (TypeError, ValueError):
                            if a != b:
                                raise AssertionError(f"Different label: {name}, row {index}, {field}")
                        else:
                            if not (np.isfinite(a_number) and np.isfinite(b_number)
                                    and np.isclose(a_number, b_number, rtol=RTOL, atol=ATOL)):
                                raise AssertionError(f"Different numerical value: {name}, row {index}, {field}")
            checked.append(name)
    for name in ("cutoff_limits.pdf", "small_time.pdf", "shear_blindness.pdf"):
        figure = generated / "figures" / name
        if not figure.is_file() or figure.stat().st_size == 0:
            raise AssertionError(f"Missing figure: {name}")
    if (generated / "benchmark_table.tex").read_text() != (reference / "benchmark_table.tex").read_text():
        raise AssertionError("The rounded manuscript table differs from the reference")
    return {"status": "passed", "numerical_files_checked": checked,
            "rtol": RTOL, "atol": ATOL,
            "figure_check": "existence only; PDF metadata are not compared",
            "source_metadata": "matched current implementation"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generated", required=True, type=Path)
    parser.add_argument("--reference", type=Path, default=ROOT)
    args = parser.parse_args()
    print(json.dumps(verify(args.generated, args.reference), indent=2))


if __name__ == "__main__":
    main()
