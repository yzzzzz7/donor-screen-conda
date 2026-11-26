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

from zipfile import ZipFile
from xml.etree import ElementTree
for folder, count in (("01_selection", 2), ("02_pairings", 2), ("03_combinations", 3)):
    with ZipFile(run / folder / "results.xlsx") as book:
        assert book.testzip() is None
        sheet = ElementTree.fromstring(book.read("xl/worksheets/sheet1.xml"))
        assert len(sheet.find("{*}sheetData")) == count + 1
assert len(list(run.glob("0*/*.svg"))) == 4
assert (run / "report.html").read_text(encoding="utf-8").count("<svg ") == 8
