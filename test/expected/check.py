import json
from pathlib import Path
run = Path("result")
report = json.loads((run / "results.json").read_text())
assert report["cohort"]["offspring_count"] == 1
assert report["cohort"]["het_sites"] > 0
assert [r["k"] for r in report["results"]] == [0, 1]
assert len(report["results"][1]["donors"]) == 1
assert (run / "results.csv").stat().st_size > 0
print("Full pipeline smoke test passed")
