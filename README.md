# donor-screen

Public Conda channel for donor-screen v0.2.0 on Linux x86_64.

```bash
conda install -c https://raw.githubusercontent.com/yzzzzz7/donor-screen-conda/main -c conda-forge -c bioconda donor-screen
donor-screen --help
donor-screen joint.vcf.gz roles.tsv --out run1
```

Installation requires access to Conda channels, the Beagle author's site and CRAN.
See [usage.md](usage.md) for the VCF and sample-sheet format.
This is a test release. Native end-to-end Linux validation is pending.
