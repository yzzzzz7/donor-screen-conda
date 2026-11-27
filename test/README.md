# Synthetic public example

[Synthetic VCF](input/synthetic-contrast.vcf.gz) · [Complete pairwise Excel](expected/outputs/05_pairwise/results.xlsx) · [Bilingual report](https://yzzzzz7.github.io/donor-screen-conda/)

The repository ships one example VCF and one expected Excel workbook. The independent artificial dataset has six recipients, three candidate donors, three additional reference samples, and 8,680 markers on one simulated chromosome. All genotypes and CDS labels are synthetic.

Run this command from the repository root after installing donor-screen 0.2.10.

```bash
donor-screen test/input/synthetic-contrast.vcf.gz test/input/roles.tsv --froh-bed test/input/autosomes.bed --cds-sites test/input/synthetic-cds.tsv --froh-chr-set 1 -k 2 -n 3 -r 50 -t 2 --java-memory 2g --out demo-two-donors
```

The role sheet, reference BED and CDS labels are in [input/](input/). The generator [make-contrast-vcf.py](make-contrast-vcf.py) defines a fixed seed, low-heterozygosity background markers, and constructed homozygous runs.

To regenerate the fixture with the same seed, run the following from the repository root.

```bash
python test/make-contrast-vcf.py | bgzip -c > test/input/synthetic-contrast.vcf.gz
python test/make-contrast-vcf.py --cds > test/input/synthetic-cds.tsv
```
