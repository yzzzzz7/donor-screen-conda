import csv
import json
from pathlib import Path
run = Path("result")
report = json.loads((run / "results.json").read_text())
assert report["cohort"]["offspring_count"] == 1
assert report["cohort"]["het_sites"] > 0
assert [r["k"] for r in report["results"]] == [0, 1]
assert len(report["results"][1]["donors"]) == 1
with (run / "results.csv").open(newline="") as generated, (Path(__file__).with_name("results.csv")).open(newline="") as expected:
    generated_csv = csv.DictReader(generated)
    expected_csv = csv.DictReader(expected)
    assert generated_csv.fieldnames == expected_csv.fieldnames
    assert len(list(generated_csv)) == len(list(expected_csv)) == 2
print("Full pipeline smoke test passed")
