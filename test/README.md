# Public VCF examples

[View all three reports online](https://yzzzzz7.github.io/donor-screen-conda/)

From the repository root after installing donor-screen 0.2.3:

    bash test/run-all.sh demo-run

Use `bash test/run.sh demo-run` for the original screening example alone.

- [input/](input/) contains the example VCF, sample roles and additional inputs for R and FROH.
- [expected/outputs/](expected/outputs/) contains the three HTML reports, result tables and figures. The original CSV is retained at `expected/results.csv` for the smoke check.
- `fetch-input.sh` regenerates the VCF from the public callset. Run `bash test/fetch-input.sh regenerated.vcf.gz`.

The example uses 12 public 1000 Genomes Phase 1 samples and 8,576 biallelic SNPs from hg19 chr22:20,000,000–21,000,000. Sample roles and coding/missense annotations are assigned for software testing. The observed FROH example divides detected lengths in the fragment by the chr22 reference length. Treat these outputs as software examples. [IGSR data policy](https://www.internationalgenome.org/IGSR_disclaimer/).
