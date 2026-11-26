# donor-screen

Screen donor–recipient pairings for genetic rescue from a joint SNP VCF. The software applies kinship and predicted recent HBD limits, then ranks eligible plans by their expected gain of rare alleles absent from the recipient population.

## Install on Linux

The public Conda channel provides **donor-screen 0.2.9** for Linux x86-64. Run this in the Conda environment where you want to use the tool.

```bash
conda install -c https://raw.githubusercontent.com/yzzzzz7/donor-screen-conda/main \
  -c conda-forge -c bioconda donor-screen=0.2.9
donor-screen --help
```

Installation fetches PLINK, bcftools, Java and R dependencies, along with the required Beagle and RZooRoH versions. KING is provided by PLINK 2. The [official Bioconda submission](https://github.com/bioconda/bioconda-recipes/pull/72231) is under review. See the [Linux installation test](https://github.com/yzzzzz7/donor-screen-conda/actions/runs/45801557968).

## Copy-and-run example: two donors, three F1 slots

After installing donor-screen 0.2.9, open a terminal in this repository's root and copy the entire command. It uses the files in [test/input/](test/input/) and creates a new output directory named `demo-two-donors`.

```bash
donor-screen test/input/fragment.vcf.gz test/input/roles.tsv --froh-bed test/input/autosomes.bed --cds-sites test/input/synthetic-cds.tsv --froh-chr-set 22 -k 2 -n 3 -r 50 -t 2 --java-memory 2g --out demo-two-donors
```

The example evaluates exactly **2 donors** for **3 planned F1** and simulates **50 F1 per parental pair** to estimate FROH. All **18 recipient–donor pairs** have both FROH columns filled in the Excel workbook.

| Plan | Donors | Expected new rare alleles | F1 SNP heterozygosity | Recent HBD (proportion) | Total FROH | Long FROH (>1 Mb) | R |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Selected | NA12812, NA12815 | 54.500 | 7.315% | 1.63e-09 | 0.486% | 0.000% | 52.582% |
| High-heterozygosity reference | NA12812, NA12815 | 54.500 | 7.315% | 1.63e-09 | 0.486% | 0.000% | 52.582% |
| No introduction | — | 0.000 | 6.347% | 2.97e-09 | 0.509% | 0.000% | 50.946% |

[Open the English / Chinese report](https://yzzzzz7.github.io/donor-screen-conda/04_two_donors/report.html) · [Download the 18-pair Excel workbook](https://yzzzzz7.github.io/donor-screen-conda/04_two_donors/05_pairwise/results.xlsx) · [Download the plan comparison Excel workbook](https://yzzzzz7.github.io/donor-screen-conda/04_two_donors/04_predicted_froh/results.xlsx) · [View the Linux test run](https://github.com/yzzzzz7/donor-screen-conda/actions/runs/45806224958)

This example uses a short public chr22 VCF fragment and synthetic CDS annotations. The results demonstrate the workflow.

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

## Other public examples

[Open the bilingual example report](https://yzzzzz7.github.io/donor-screen-conda/02_screen_with_R/report.html) · [Download the 18-pair Excel workbook](https://yzzzzz7.github.io/donor-screen-conda/02_screen_with_R/05_pairwise/results.xlsx)

From a checkout of this repository, run `bash test/run-all.sh demo-run` after installation. The [test directory](test/) contains the inputs and commands. The example uses a short public chromosome 22 VCF fragment and synthetic CDS labels to demonstrate the workflow.
