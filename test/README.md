# Public VCF test

From the repository root, with donor-screen installed:

    bash test/run.sh test-run

- input/ holds the sample role sheet and example VCF.
- expected/ holds the automated checks and reference output.
- fetch-input.sh regenerates the VCF from the public upstream callset. Run it from a working directory containing roles.tsv.

The example uses 12 public 1000 Genomes Phase 1 samples and 8,576 biallelic SNPs from hg19 chromosome 22 (20,000,000-21,000,000). The roles are assigned for software testing. The test runs one F1 comparison against the no-introduction control, then checks the generated CSV and JSON. See [IGSR data-use policy](https://www.internationalgenome.org/data-use-policy/).
