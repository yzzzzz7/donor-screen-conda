# Donor-screen usage

Install the Linux package by following the [repository README](README.md).

```bash
donor-screen joint.vcf.gz roles.tsv --out run1
```

Use a quality-filtered, sorted, diploid, biallelic SNP VCF. The tab-delimited sample sheet starts with `#autosomes=1,2` and `#reference_bp=123456789` settings followed by a `sample` and `role` header. Roles are `recipient`, `donor`, or `reference`. Set `#king_panel=joint` to use the joint VCF for KING or supply the original independent panel with `--king-vcf`. Use `--offspring` and `--max-donors` to set quotas of 1–3. Each F1 requires two distinct recipients in the no-introduction control. Choose a new output directory.

## Common options

| Option | Meaning |
| --- | --- |
| -r, --froh-replicates | Simulated F1 offspring per parental pair for FROH, default 100. Supply --froh-bed to enable prediction. |
| -t, --threads | Threads for PLINK 2 KING and Beagle phasing, default 4. FROH batches run sequentially. |
| -k, --donors | Evaluate exactly this many donors (1–3) and the no-introduction control. |
| -n, --offspring | Planned F1 offspring per scenario (1–3), default 3. |
| --max-donors | Compare all donor counts from 1 to the specified maximum, alongside the no-introduction control. Use either this option or -k. |

Choose at most one donor per F1 slot. With neither donor-count option, the program evaluates all counts up to the smaller of the F1 quota and the number of candidate donors. To select two donors for a plan with three F1 slots and predict FROH from 50 simulated F1 per parental pair:

    donor-screen joint.vcf.gz roles.tsv --froh-bed autosomes.bed -k 2 -n 3 -r 50 -t 8 --out run1

The simulation count applies to FROH per pair. It does not increase the three F1 slots in the plan.

Supply a curated TSV or TSV.GZ with `chrom`, `pos`, `ref`, `alt`, and `missense` columns using `--cds-sites` to calculate R. A zero flag indicates another CDS SNP included in the denominator. A one flag indicates a missense SNP. Coordinates and alleles must match the VCF. Alternatively, use `--cds-panel` with the same columns plus observed ALT dosages by sample. The default `--r-policy report` reports R. The `reference` and `baseline` policies impose the high-heterozygosity and no-introduction R limits, respectively.

Use `--froh-bed` with whole-autosome BED rows (chromosome, zero, length) to predict F1 FROH for every recipient–donor pair and for the pairings used by the comparison plans. Defaults are 100 simulated F1 per pair, seed 50483, and 1 cM/Mb. PLINK 1.9 identifies total ROH at least 100 kb and long ROH over 1 Mb. Use `--froh-replicates`, `--froh-seed` and `--froh-chr-set` to specify these run settings. For observed FROH, run `donor-screen froh --help`.

Open `report.html` locally and use the **中文 / English** button to switch languages. The result tables and Excel workbooks are in English. The pairwise workbook at `05_pairwise/results.xlsx` includes one row per recipient and candidate donor. Its heterozygosity and rare-allele columns describe **one F1**, with rare-allele gain measured against the recipient population. At plan level, an allele counts once across all proposed F1. Pairwise R requires CDS input. When FROH prediction is enabled, both FROH columns are filled for every recipient–donor pair, including disallowed pairs. Screening criteria are unchanged.

Each pairwise row includes the original signed KING kinship coefficient for the two parents. The run root also contains the complete native PLINK 2 `king.kin0` table. KING is included in PLINK 2; no separate KING installation is needed.

KING is the first column of the pairwise workbook. The program retains the original selected and reference plan simulations and simulates the remaining recipient–donor pairs in deterministic batches. Run time increases with the number of pairs.

### Reading the screening figures

For each feasible nonzero donor budget, `01_selection/selected-plan-k*.svg` identifies the selected donors and reports expected distinct new rare alleles across the planned offspring, mean F1 SNP heterozygosity, recent HBD, R and optional FROH estimates. `02_pairings/planned-matings-k*.svg` links each donor to the assigned recipient and gives the number of planned F1. `03_combinations/eligible-combinations.svg` shows how many evaluated combinations passed the constraints. The HTML report embeds the figures. The Excel tables retain the high-heterozygosity reference and all pairwise estimates.
