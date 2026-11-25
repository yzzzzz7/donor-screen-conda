"""Write Excel tables, SVG plots and a self-contained HTML run report."""
from html import escape
from itertools import chain, islice
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

XMLNS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
RELS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG = "http://schemas.openxmlformats.org/package/2006"
MAX_DATA_ROWS = 1_048_575


def label(value):
    if value is None:
        return ""
    if isinstance(value, list):
        return ", ".join(map(str, value))
    return str(value)


def cell(value):
    if value is None:
        return "<c/>"
    if isinstance(value, bool):
        return f'<c t="b"><v>{int(value)}</v></c>'
    if isinstance(value, (float, int)):
        return f"<c><v>{value}</v></c>"
    return f'<c t="inlineStr"><is><t xml:space="preserve">{escape(label(value))}</t></is></c>'


def write_xlsx(path, columns, records):
    headers = tuple(columns.values())
    widths = "".join(f'<col min="{i}" max="{i}" width="{min(44, max(16, len(title) + 2))}" customWidth="1"/>'
                     for i, title in enumerate(headers, 1))
    rows = iter(records)
    names = []
    with ZipFile(path, "w", ZIP_DEFLATED) as book:
        while (first := next(rows, None)) is not None or not names:
            index = len(names) + 1
            names.append(f'<sheet name="Results {index}" sheetId="{index}" r:id="rId{index}"/>')
            with book.open(f"xl/worksheets/sheet{index}.xml", "w") as sheet:
                sheet.write(f'<worksheet xmlns="{XMLNS}"><cols>{widths}</cols><sheetData>'.encode())
                sheet.write(("<row>" + "".join(map(cell, headers)) + "</row>").encode())
                chunk = () if first is None else chain((first,), islice(rows, MAX_DATA_ROWS - 1))
                for record in chunk:
                    cells = "".join(cell(record.get(key)) for key in columns)
                    sheet.write(f"<row>{cells}</row>".encode())
                sheet.write(b"</sheetData></worksheet>")
        sheet_type = "application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"
        book_type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"
        indices = range(1, len(names) + 1)
        overrides = "".join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="{sheet_type}"/>' for i in indices)
        links = "".join(f'<Relationship Id="rId{i}" Type="{RELS}/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in indices)
        book.writestr("[Content_Types].xml", f'<Types xmlns="{PKG}/content-types">'
                      f'<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                      f'<Default Extension="xml" ContentType="application/xml"/>'
                      f'<Override PartName="/xl/workbook.xml" ContentType="{book_type}"/>'
                      f'<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>{overrides}</Types>')
        book.writestr("_rels/.rels", f'<Relationships xmlns="{PKG}/relationships"><Relationship Id="rId1" '
                      f'Type="{RELS}/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        book.writestr("xl/workbook.xml", f'<workbook xmlns="{XMLNS}" xmlns:r="{RELS}"><sheets>{"".join(names)}</sheets></workbook>')
        book.writestr("xl/_rels/workbook.xml.rels", f'<Relationships xmlns="{PKG}/relationships">{links}'
                      f'<Relationship Id="styles" Type="{RELS}/styles" Target="styles.xml"/></Relationships>')
        book.writestr("xl/styles.xml", f'<styleSheet xmlns="{XMLNS}">'
                      '<fonts count="1"><font><sz val="11"/><color rgb="FF000000"/><name val="Times New Roman"/></font></fonts>'
                      '<fills count="2"><fill><patternFill patternType="none"/></fill>'
                      '<fill><patternFill patternType="gray125"/></fill></fills>'
                      '<borders count="1"><border/></borders>'
                      '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
                      '<cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>'
                      '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
                      '</styleSheet>')


SUMMARY = {
    "k": "Donor count", "donors": "Selected donors", "rare": "Expected new rare alleles",
    "het": "Predicted offspring heterozygosity", "hbd": "Predicted recent HBD",
    "cap": "Recent HBD limit", "eligible_count": "Eligible combinations",
    "possible_count": "Combinations evaluated", "reference_donors": "High-heterozygosity reference donors",
    "reference_rare": "Reference expected new rare alleles", "reference_het": "Reference offspring heterozygosity",
    "gain_over_reference": "Expected rare alleles gained over reference", "runner_up": "Runner-up",
    "runner_gap": "Expected rare allele gap over runner-up", "reason": "Selection rationale",
}
PAIRS = {
    "k": "Donor count", "parent1": "Parent 1", "parent2": "Parent 2", "n_offspring": "Offspring",
    "het": "Predicted offspring heterozygosity", "hbd": "Predicted recent HBD",
}
COMBINATIONS = {
    "k": "Donor count", "donors": "Candidate donors", "pass": "Eligible", "rare": "Expected new rare alleles",
    "het": "Predicted offspring heterozygosity", "hbd": "Predicted recent HBD", "cap": "Recent HBD limit",
}


def plot(title, rows, first, second, legend, percentage=False):
    """Draw paired bars against a shared zero-based scale, including missing values."""
    scale = 100 if percentage else 1
    values = [(f'{r["k"]} donor(s)' if r["k"] else "Control", r.get(first), r.get(second)) for r in rows]
    maximum = max((x for _, a, b in values for x in (a, b) if x is not None), default=0)
    width = 760
    height = 96 + 49 * len(values)
    color = ("#245d91", "#be6c35")
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="{escape(title)}">',
           f'<rect width="{width}" height="{height}" fill="#fff"/>',
           '<g font-family="system-ui,sans-serif" font-size="13" fill="#233245">',
           f'<text x="16" y="25" font-size="17" font-weight="600">{escape(title)}</text>']
    for i, name in enumerate(legend):
        svg.append(f'<rect x="{16 + i * 230}" y="43" width="12" height="12" fill="{color[i]}"/>')
        svg.append(f'<text x="{34 + i * 230}" y="54">{escape(name)}</text>')
    for i, (name, *numbers) in enumerate(values):
        y = 75 + i * 49
        svg.append(f'<text x="16" y="{y + 20}">{escape(name)}</text>')
        for j, number in enumerate(numbers):
            ybar = y + j * 19
            bar = 500 * number / maximum if maximum and number is not None else 0
            shown = "—" if number is None else f"{number * scale:.6g}" + ("%" if percentage else "")
            svg.append(f'<rect x="119" y="{ybar}" width="{bar:.3f}" height="13" fill="{color[j]}"/>')
            svg.append(f'<text x="{124 + bar:.3f}" y="{ybar + 11}">{shown}</text>')
    return "".join(svg) + "</g></svg>"


def html_table(columns, rows):
    header = "".join(f"<th scope='col'>{escape(name)}</th>" for name in columns.values())
    body = []
    for row in rows:
        cells = []
        for key in columns:
            value = row.get(key)
            if isinstance(value, float):
                value = f"{value:.6g}"
            cells.append(f"<td>{escape(label(value)) or '—'}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f"<div class='table'><table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"


def write_report(out):
    data = json.loads((out / "results.json").read_text(encoding="utf-8"))
    results = data["results"]
    cohort = data["cohort"]
    pairs = [dict(pair, k=row["k"]) for row in results for pair in row["pairs"]]
    groups = [
        ("01_selection", "Selection summary", SUMMARY, results),
        ("02_pairings", "Selected pairings", PAIRS, pairs),
        ("03_combinations", "Candidate combinations", COMBINATIONS, data["alternatives"]),
    ]
    for folder, _, columns, rows in groups:
        (out / folder).mkdir()
        write_xlsx(out / folder / "results.xlsx", columns, rows)

    plots = [
        ("01_selection", "rare-alleles.svg", "Expected new rare alleles", "rare", "reference_rare",
         ("Selected plan", "High-heterozygosity reference"), False),
        ("01_selection", "offspring-heterozygosity.svg", "Predicted offspring heterozygosity", "het", "reference_het",
         ("Selected plan", "High-heterozygosity reference"), True),
        ("02_pairings", "recent-hbd.svg", "Mean predicted recent HBD across selected pairings", "hbd", "cap",
         ("Selected plan", "HBD limit"), True),
        ("03_combinations", "eligible-combinations.svg", "Candidate combinations", "eligible_count", "possible_count",
         ("Eligible", "Evaluated"), False),
    ]
    graphics = {}
    for folder, filename, title, first, second, legend, percentage in plots:
        image = plot(title, results, first, second, legend, percentage)
        (out / folder / filename).write_text(image, encoding="utf-8")
        graphics.setdefault(folder, []).append(image)

    sections = []
    for folder, title, columns, rows in groups:
        if folder == "01_selection":
            keys = ("k", "donors", "rare", "het", "hbd", "reference_donors", "gain_over_reference")
            table = html_table({key: columns[key] for key in keys}, rows)
        elif folder == "02_pairings":
            table = html_table(columns, rows)
        else:
            table = f"<p>{len(rows):,} candidate combinations evaluated. The workbook lists each combination.</p>"
        sections.append(f'<section><h2>{title}</h2><p><a href="{folder}/results.xlsx">Download Excel</a></p>'
                        f'{table}{"".join(graphics[folder])}</section>')
    page = f"""<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Donor-screen report</title>
<style>
body{{font:15px/1.5 system-ui,sans-serif;color:#233245;max-width:1024px;margin:32px auto;padding:0 20px}}
section{{border-top:1px solid #cbd5df;margin-top:26px}} .table{{overflow:auto}} table{{border-collapse:collapse;width:100%}}
th,td{{text-align:left;padding:8px 10px;border-bottom:1px solid #e0e6ec}} th{{background:#f0f4f8}}
svg{{display:block;width:100%;max-width:760px;margin:12px 0}} a{{color:#245d91}}
@media print{{section{{break-inside:avoid}}}}
</style>
<h1>Donor-screen report</h1><p>Run {escape(cohort["name"])} · {cohort["offspring_count"]} simulated F1 offspring per plan ·
{len(cohort["recipients"])} recipients · {len(cohort["candidates"])} candidate donors</p>
<p>Heterozygosity and HBD are proportions in the tables and percentages in the plots.
An em dash indicates an unavailable or inapplicable value. The control introduces zero donors.</p>
<p><a href="results.csv">Original results CSV</a> · <a href="results.json">Original results JSON</a></p>
{"".join(sections)}
</html>"""
    (out / "report.html").write_text(page, encoding="utf-8")
