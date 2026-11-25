# Linux 供体筛选

## Conda 安装（Linux）

解压后进入软件目录，运行：

```bash
bash install-conda.sh
conda activate donor-screen
donor-screen joint.vcf.gz roles.tsv --out run1
```

安装脚本创建独立环境，从作者网站下载原分析的 Beagle jar，并从 CRAN 编译 RZooRoH 0.4.3。需要网络和已安装的 Conda。若同名环境已存在，可用 `bash install-conda.sh donor-screen-new` 指定新名称。安装完成会检查 RZooRoH 版本和 CLI 帮助。

## 手动安装

安装 Python ≥3.10、bcftools、PLINK 2、Java、Beagle 5.4（原分析 jar：beagle.25Mar22.4f6.jar）、R 和 RZooRoH 0.4.3，前四个命令及 Rscript 置于 PATH。将本包解压到同一目录，运行 chmod +x donor-screen，然后 export BEAGLE_JAR=/path/to/beagle.25Mar22.4f6.jar。

运行： ./donor-screen joint.vcf.gz roles.tsv --out run1

roles.tsv 用制表符，示例：

    #autosomes=1,2
    #reference_bp=123456789
    #king_panel=joint
    sample	role
    R01	recipient
    R02	recipient
    R03	recipient
    R04	recipient
    R05	recipient
    R06	recipient
    D01	donor
    D02	donor
    D03	donor
    F01	reference

请换成实际常染色体范围和现有个体杂合度使用的参考分母。所有列出的样本组成频率参照和亲缘计算队列，reference 仅作参照。使用质控完成、按基因组位置排序的联合二倍体双等位 SNP VCF。默认每方案 3 个 F1、最多 3 个供体，至少需 6 只受体及相应候选。可以用 --offspring 和 --max-donors 指定配额。运行目录须尚不存在。

#king_panel=joint 表示明确采用联合 VCF 的指定常染色体 SNP 做 KING。复现原论文的独立 LD 面板时，用 --king-vcf adopted_king.vcf.gz，并提供原频率参照和亲缘分簇队列。程序不会另设 LD 或缺失阈值；联合 SNP 面板与原独立面板可产生不同的亲缘结果。先对基因型、参考范围和队列进行核对。

保持 KING 0.177、Beagle seed 50483、RZooRoH 六速率和 1 cM/Mb 假设、近期 HBD Rc≤25、原有最大化首次传递方案。HBD 不收敛时停止。输出 results.json 和 results.csv，保留日志及原生 HBD 批次。只提供两个输入文件时计算前三项。可按下节传入 CDS 注释计算第四项 R。


安装 R 包时可用 Rscript -e 'remotes::install_version("RZooRoH", version="0.4.3")'（先安装 remotes，必要时使用系统的编译工具）。

## Output report

Open `report.html` in a browser. It contains the summary, selected pairings, and four charts (five when R is available). Its Excel links work when the report stays alongside the numbered folders. The browser can print the report to PDF.

| Folder | Contents |
| --- | --- |
| `01_selection/` | Selected plans and high-heterozygosity reference plans, with rare-allele and heterozygosity charts |
| `02_pairings/` | Selected parents and predicted offspring metrics, with a recent-HBD chart |
| `03_combinations/` | All candidate combinations and their eligibility, with a count chart |

Each folder contains a `results.xlsx` workbook. The original `results.csv` and `results.json` remain in the run root. Tables store the original proportions; charts label heterozygosity and HBD as percentages.

## 第四项指标 R

```bash
donor-screen joint.vcf.gz roles.tsv --cds-sites coding.tsv --out run-r
```

`--cds-sites` 接受 TSV 或 TSV.GZ，列名为 `chrom`、`pos`、`ref`、`alt`、`missense`。请提供分析范围内完整的、经人工或既定注释流程确认的 CDS SNP 面板，错义位点标 1，其他纳入分母的 CDS 位点标 0。坐标和两个等位基因必须与 VCF 一致。程序使用原始观测基因型的 ALT 剂量。若原分析有独立 CDS 面板，用 `--cds-panel` 传入同样五列加受体和供体各样本的 ALT 剂量列（0、1、2）。纳入位点缺失或配对分母为零时停止。

每个后代的 R 为预期 ALT 纯合错义位点数除以预期 ALT 纯合 CDS 位点数，方案按后代取平均。默认 `--r-policy report` 计算并展示 R，保留原三指标入选方案。`--r-policy reference` 或 `baseline` 分别采用原筛选器内置的高杂合度对照或无引入对照 R 上限。R 写入 CSV、JSON、Excel 和 HTML 报告，附 R 对比图。
