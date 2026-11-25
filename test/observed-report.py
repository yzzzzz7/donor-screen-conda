"""Create an HTML and SVG demonstration from the observed FROH output."""
from html import escape
from pathlib import Path
import csv
import sys

root = Path(sys.argv[1])
with (root / "froh.tsv").open(encoding="utf-8", newline="") as stream:
    rows = list(csv.DictReader(stream, delimiter="\t"))
columns = ("FROH_100kb", "FROH_gt1Mb")
numbers = [[100 * float(row[key]) for key in columns] for row in rows]
scale = max((max(row) for row in numbers), default=0) or 1
height = 90 + 30 * len(rows)
svg = [
    f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 740 {height}" role="img" aria-label="Observed FROH per sample">',
    f'<rect width="740" height="{height}" fill="white"/>',
    '<g font-family="system-ui,sans-serif" font-size="13" fill="#233245">',
    '<text x="16" y="25" font-size="17" font-weight="600">Observed FROH by sample (%)</text>',
    '<line x1="16" x2="35" y1="45" y2="45" stroke="#245d91" stroke-width="6"/><text x="42" y="50">Total ≥100 kb</text>',
    '<line x1="176" x2="195" y1="45" y2="45" stroke="#be6c35" stroke-width="6"/><text x="202" y="50">Long &gt;1 Mb</text>',
]
for i, (row, values) in enumerate(zip(rows, numbers)):
    y = 72 + i * 30
    svg.append(f'<text x="16" y="{y + 6}">{escape(row["sample"])}</text>')
    for j, (value, color) in enumerate(zip(values, ("#245d91", "#be6c35"))):
        offset = 12 * j
        width = 480 * value / scale
        svg.append(f'<rect x="110" y="{y + offset - 7}" width="{width:.2f}" height="6" fill="{color}"/>')
        svg.append(f'<text x="{116 + width:.2f}" y="{y + offset}" font-size="11">{value:.5g}%</text>')
figure = "".join(svg) + "</g></svg>"
(root / "froh.svg").write_text(figure, encoding="utf-8")
table = "".join(
    f'<tr><th scope="row">{escape(row["sample"])}</th><td>{values[0]:.5g}%</td><td>{values[1]:.5g}%</td></tr>'
    for row, values in zip(rows, numbers)
)
page = f"""<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Observed FROH example</title>
<style>
body{{font:15px/1.5 system-ui,sans-serif;color:#233245;max-width:850px;margin:36px auto;padding:0 20px}}
table{{border-collapse:collapse}} th,td{{padding:6px 16px;text-align:left;border-bottom:1px solid #cbd5df}}
svg{{width:100%;max-width:740px;display:block}} a{{color:#245d91}}
</style>
<h1>Observed FROH example</h1>
<p>Demo using 12 public 1000 Genomes samples, approximately 1 Mb of markers on chr22 and the whole-chromosome denominator of 51,304,566 bp. Whole-chromosome and genome-wide estimates require full reference coverage. These samples are used to exercise the software.</p>
<p><a href="froh.xlsx">Excel</a> · <a href="froh.tsv">TSV</a> · <a href="roh.hom">PLINK segments</a> · <a href="froh.svg">SVG chart</a></p>
{figure}
<table><thead><tr><th>Sample</th><th>Total FROH (≥100 kb)</th><th>Long FROH (&gt;1 Mb)</th></tr></thead>
<tbody>{table}</tbody></table>
</html>"""
(root / "report.html").write_text(page, encoding="utf-8")
