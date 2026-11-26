# Public VCF examples

[View the reports online](https://yzzzzz7.github.io/donor-screen-conda/)

From the repository root after installing donor-screen 0.2.9:

    bash test/run-all.sh demo-run

Use `bash test/run.sh demo-run` for the original screening example alone.

- [input/](input/) contains the example VCF, sample roles and additional inputs for R and FROH.
- [expected/outputs/](expected/outputs/) contains the HTML reports, result tables and figures. Screening reports support English and Chinese and link to the English pairwise workbook. The five-metric report uses 16 simulated F1 per mating to exercise predicted FROH. The original CSV is retained at `expected/results.csv` for the smoke check.
- `fetch-input.sh` regenerates the VCF from the public callset. Run `bash test/fetch-input.sh regenerated.vcf.gz`.

The example uses 12 public 1000 Genomes Phase 1 samples and 8,576 biallelic SNPs from hg19 chr22:20,000,000–21,000,000. Sample roles and coding/missense annotations are assigned for software testing. The observed FROH example divides detected lengths in the fragment by the chr22 reference length. Treat these outputs as software examples. [IGSR data policy](https://www.internationalgenome.org/IGSR_disclaimer/).

The pairwise workbook starts with the signed KING kinship coefficient for each recipient–donor pair.

The five-metric example provides predicted FROH for all 18 recipient–donor pairs.

The joint VCF is named `donor-recipient.vcf.gz`. Its samples use role-based English IDs (Donor_01, Recipient_01, and Reference_01). [sample-names.tsv](input/sample-names.tsv) maps each ID to the public source sample. Variant records and genotype calls are unchanged.
