# donor-screen

Public Conda channel for donor-screen v0.2.3 on Linux x86_64.

```bash
conda install -c https://raw.githubusercontent.com/yzzzzz7/donor-screen-conda/main -c conda-forge -c bioconda donor-screen
donor-screen --help
donor-screen joint.vcf.gz roles.tsv --out run1
```

Installation requires access to Conda channels, the Beagle author's site and CRAN.
Use `conda update donor-screen` with the same channel to upgrade an existing installation. With optional curated CDS input, the report includes the fourth metric, R. Each new run generates a local HTML report with Excel workbooks and charts. See [usage.md](usage.md) for the VCF and sample-sheet format.
The native Linux workflow passed an [end-to-end smoke test](https://github.com/yzzzzz7/donor-screen-conda/actions/runs/45628034716) on the public example below (12 samples and 8,576 SNPs). This test runs the workflow from source; it is separate from Conda package validation.

For sequenced parents or offspring, run `donor-screen froh --help` to measure total and long FROH.

## Public VCF example

Run `bash test/run.sh test-run` from the repository root after installation. See [test files and instructions](test/README.md).
