"""Localize static report text while leaving numeric results and sample IDs intact."""
import re

ZH = {
    'Selected plan profile': '入选供体及方案预测',
    'Expected F1 SNP heterozygosity': '预计 F1 SNP 杂合率',
    'Planned matings': '计划配对',
    'Donor': '供体',
    'Recipient': '受体',
    'F1 slots': '计划 F1 数',
    'Feasible donor combinations': '可行供体组合',
    'Eligible / evaluated': '达标组合数 / 评估组合数',
    'Total FROH': '总纯合片段覆盖率',
    'Long FROH (>1 Mb)': '长纯合片段覆盖率（>1 Mb）',
    'The selected plan shows expected distinct new rare alleles across its planned offspring. Recent HBD is a proportion; heterozygosity, R and FROH are percentages in the profile. An em dash indicates an unavailable or inapplicable value. The control introduces zero donors.': '预计新增稀有等位基因数按整个繁育方案去重。近期 HBD 用比例表示；方案图的杂合率、R 和纯合片段覆盖率用百分比表示。横线表示无相应数据。对照为不引入供体。',

    "Donor-screen report":"供体筛选报告",
    "Run":"运行",
    "Selection summary":"筛选汇总",
    "Selected pairings":"入选配对",
    "Candidate combinations":"候选组合",
    "Recipient–donor pairwise results":"受体与供体逐对结果",
    "Download Excel":"下载 Excel",
    "Original results CSV":"原始 CSV 结果",
    "Original results JSON":"原始 JSON 结果",
    "Donor count":"供体数",
    "Selected donors":"入选供体",
    "Expected new rare alleles":"预计新增稀有等位基因数",
    "Predicted offspring heterozygosity":"预测后代杂合度",
    "Predicted recent HBD":"预测近期 HBD",
    "High-heterozygosity reference donors":"高杂合度对照供体",
    "Expected rare alleles gained over reference":"相对对照预计多获得的稀有等位基因数",
    "Parent 1":"亲本 1","Parent 2":"亲本 2","Offspring":"后代数",
    "Predicted homozygous missense proportion R":"预测纯合错义比例 R",
    "Predicted R":"预测 R","Reference R":"对照 R","R limit":"R 上限",
    "Selected plan":"入选方案","High-heterozygosity reference":"高杂合度对照",
    "Control":"无引入对照","HBD limit":"HBD 上限","Eligible":"达标","Evaluated":"已评估",
    "Mean predicted recent HBD across selected pairings":"入选配对的平均预测近期 HBD",
    "Predicted offspring FROH":"后代 FROH 预测",
    "Plan":"方案","Donors":"供体","Selected":"入选方案","No introduction":"无引入对照",
    "Mean total FROH (≥100 kb)":"平均总 FROH（≥100 kb）",
    "Mean long FROH (>1 Mb)":"平均长片段 FROH（>1 Mb）",
    "Total FROH 2.5th percentile":"总 FROH 第 2.5 百分位",
    "Total FROH 97.5th percentile":"总 FROH 第 97.5 百分位",
    "Long FROH 2.5th percentile":"长片段 FROH 第 2.5 百分位",
    "Long FROH 97.5th percentile":"长片段 FROH 第 97.5 百分位",
    "Excel":"Excel","Simulated F1 results":"模拟 F1 逐次结果",
    "Heterozygosity and HBD are proportions in the tables and percentages in the plots. An em dash indicates an unavailable or inapplicable value. The control introduces zero donors.":"表中的杂合度和 HBD 为比例，图中以百分比表示。长破折号表示没有结果或不适用。无引入对照的供体数为零。",
    "Each row is one recipient–donor pair. Heterozygosity and rare-allele gain describe one F1. With FROH prediction enabled, every recipient–donor pair is simulated, including pairs that fail the pairing filter. R requires CDS annotations.":"每行对应一个受体与供体配对。杂合度和新增稀有等位基因数均按单个 F1 计算。启用 FROH 预测时模拟每一个受体与供体配对，包括未通过筛选的配对。R 需要 CDS 注释。",
    "R is the mean across predicted offspring of the expected homozygous ALT missense count divided by the expected homozygous ALT CDS count. R is reported alongside the three-metric choice.":"逐个后代计算预期纯合 ALT 错义位点数与预期纯合 ALT CDS 位点数之比后，取方案均值。R 与原三指标筛选结果并列展示。",
    "R is the mean across predicted offspring of the expected homozygous ALT missense count divided by the expected homozygous ALT CDS count. R is capped at the high-heterozygosity reference plan.":"逐个后代计算预期纯合 ALT 错义位点数与预期纯合 ALT CDS 位点数之比后，取方案均值。R 上限取高杂合度对照方案值。",
    "R is the mean across predicted offspring of the expected homozygous ALT missense count divided by the expected homozygous ALT CDS count. R is capped at the no-introduction control.":"逐个后代计算预期纯合 ALT 错义位点数与预期纯合 ALT CDS 位点数之比后，取方案均值。R 上限取无引入对照值。",
}


def translate(match):
    raw = match.group(1)
    text = " ".join(raw.split())
    if text.replace("&gt;", ">") in ZH:
        return ">" + ZH[text.replace("&gt;", ">")].replace(">", "&gt;") + "<"
    if (m := re.fullmatch(r"(\d+) donor\(s\)", text)):
        return f">{m[1]} 个供体<"
    if (m := re.fullmatch(r"Run (.*?) · (\d+) simulated F1 offspring per plan · (\d+) recipients · (\d+) candidate donors", text)):
        return f">运行 {m[1]} · 每个方案 {m[2]} 个 F1 · {m[3]} 个受体 · {m[4]} 个候选供体<"
    if (m := re.fullmatch(r"([\d,]+) candidate combinations evaluated\. The workbook lists each combination\.", text)):
        return f">评估了 {m[1]} 个候选组合。Excel 列出了每个组合。<"
    if (m := re.match(r"(\d+) simulated offspring per pair, seed (\d+)\.", text)):
        return (f">每对亲本模拟 {m[1]} 个后代，随机种子为 {m[2]}。按恒定 1 cM/Mb 重组，"
                "沿用实测样本的 PLINK ROH 参数。第 2.5 和第 97.5 百分位描述模拟方案之间的波动。"
                "FROH 用于比较方案。<")
    return match.group(0)


def bilingual(page):
    head, content = page.split("<h1>", 1)
    body = "<h1>" + content.removesuffix("</html>")
    zh = re.sub(r">([^<>]+)<", translate, body)
    return head + """<nav aria-label="Report language"><button type="button" id="language" aria-controls="english chinese">中文</button></nav>
<div id="english" lang="en">""" + body + """</div>
<div id="chinese" lang="zh-CN" hidden>""" + zh + """</div>
<script>
const button = document.getElementById('language');
button.addEventListener('click', () => {
    const english = document.getElementById('english');
    const chinese = document.getElementById('chinese');
    english.hidden = !english.hidden;
    chinese.hidden = !english.hidden;
    const isChinese = !chinese.hidden;
    document.documentElement.lang = isChinese ? 'zh-CN' : 'en';
    document.title = isChinese ? '供体筛选报告' : 'Donor-screen report';
    button.textContent = isChinese ? 'English' : '中文';
});
</script></html>"""
