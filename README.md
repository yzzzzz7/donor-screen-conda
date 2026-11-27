# donor-screen

Screen donor–recipient pairings for genetic rescue from a joint SNP VCF. The software applies kinship and predicted recent HBD limits, then ranks eligible plans by their expected gain of rare alleles absent from the recipient population.

## Install on Linux

The public Conda channel provides **donor-screen 0.2.10** for Linux x86-64. Run this in the Conda environment where you want to use the tool.

```bash
conda install -c https://raw.githubusercontent.com/yzzzzz7/donor-screen-conda/main \
  -c conda-forge -c bioconda donor-screen=0.2.10
donor-screen --help
```

Installation fetches PLINK, bcftools, Java and R dependencies, along with the required Beagle and RZooRoH versions. KING is provided by PLINK 2. The [official Bioconda submission](https://github.com/bioconda/bioconda-recipes/pull/72231) is under review. See the [Linux installation test](https://github.com/yzzzzz7/donor-screen-conda/actions/runs/45801557968).

## Copy-and-run synthetic example

After installing donor-screen 0.2.10, run this command from the repository root. It uses one synthetic VCF and plans **2 donors and 3 F1 slots**, with **50 simulated F1 per pair** for FROH.

```bash
donor-screen test/input/synthetic-contrast.vcf.gz test/input/roles.tsv --froh-bed test/input/autosomes.bed --cds-sites test/input/synthetic-cds.tsv --froh-chr-set 1 -k 2 -n 3 -r 50 -t 2 --java-memory 2g --out demo-two-donors
```

**[Test VCF](test/input/synthetic-contrast.vcf.gz)** · **[Complete 18-pair Excel workbook](https://yzzzzz7.github.io/donor-screen-conda/05_pairwise/results.xlsx)** · [Bilingual report](https://yzzzzz7.github.io/donor-screen-conda/)

Across the 18 pairs, expected F1 SNP heterozygosity ranges from **14.7% to 23.4%**. Total FROH is above zero in **18 pairs** and long FROH (>1 Mb) is above zero in **6 pairs**. The selected two-donor plan adds **650.0** expected rare alleles versus **147.5** for the high-heterozygosity reference.

The fixture uses artificial genotypes and CDS labels. The complete pairwise workbook preserves all 18 pairs, with KING in the first column and both FROH columns filled.

## Run on your data

Use a sorted, quality-filtered joint VCF of diploid biallelic SNPs and a tab-separated role sheet.

| Input | What to provide |
| --- | --- |
| `joint.vcf.gz` | Genotypes for recipients, donor candidates and any additional reference samples |
| `roles.tsv` | VCF sample IDs and their roles, plus the autosome list, reference length and KING panel setting. See the [example sample sheet](test/input/roles.tsv). |
| `autosomes.bed` | Whole-autosome coordinates for optional offspring FROH prediction, passed with `--froh-bed` |
| `coding.tsv` | Curated allele-matched CDS annotations for optional missense ratio R, passed with `--cds-sites` |

For the basic screen, which calculates expected rare-allele gain, offspring SNP heterozygosity and recent HBD:

```bash
donor-screen joint.vcf.gz roles.tsv --out run1
```

To include FROH and R on your own data, replace the test paths in the copy-and-run command above with your files.

| Option | What it controls |
| --- | --- |
| `-k 2` / `--donors 2` | Evaluate **exactly 2 donors** and the no-introduction control. Use `--max-donors 2` to compare the 1- and 2-donor scenarios alongside the control. |
| `-n 3` / `--offspring 3` | Allocate **3 F1 slots per plan**. Default: **3**. |
| `-r 50` / `--froh-replicates 50` | Simulate **50 F1 per parental pair** for FROH. Default: **100**. Requires `--froh-bed`. |
| `-t 8` / `--threads 8` | Use **8 threads** for PLINK 2 KING and Beagle. Default: **4**. FROH batches run sequentially. |

Donor count and F1 slots both support values from 1 to 3, with donor count at most the number of F1 slots. The example plans for three F1 and uses 50 simulations per pair to estimate FROH. Choose a new output directory. See the [full input formats and options](usage.md).

## Open the results

| Output | Contents |
| --- | --- |
| `report.html` | Summary tables, charts and Excel links. English by default, with a **中文 / English** switch. |
| `01_selection/`, `02_pairings/`, `03_combinations/` | Selection summaries, recommended pairings and candidate combinations |
| `05_pairwise/results.xlsx` | English workbook with a row for each recipient–donor pair. The first column is the KING coefficient. When FROH is enabled, both FROH columns are filled for every pair. |
| `04_predicted_froh/` | Optional total and long FROH predictions for comparing the selected plan, the high-heterozygosity reference and the no-introduction control |

For FROH measured from sequenced individuals, use `donor-screen froh --help`.

### Reading the screening figures

For each feasible nonzero donor budget, `01_selection/selected-plan-k*.svg` identifies the selected donors and reports expected distinct new rare alleles across the planned offspring, mean F1 SNP heterozygosity, recent HBD, R and optional FROH estimates. `02_pairings/planned-matings-k*.svg` links each donor to the assigned recipient and gives the number of planned F1. `03_combinations/eligible-combinations.svg` shows how many evaluated combinations passed the constraints. The HTML report embeds the figures. The Excel tables retain the high-heterozygosity reference and all pairwise estimates.
