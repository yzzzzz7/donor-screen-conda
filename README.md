# donor-screen

Public Conda channel for donor-screen v0.2.0 on Linux x86_64.

```bash
conda install -c https://raw.githubusercontent.com/yzzzzz7/donor-screen-conda/main -c conda-forge -c bioconda donor-screen
donor-screen --help
donor-screen joint.vcf.gz roles.tsv --out run1
```

Installation requires access to Conda channels, the Beagle author's site and CRAN.
See [usage.md](usage.md) for the VCF and sample-sheet format.
The native Linux workflow passed an [end-to-end smoke test](https://github.com/yzzzzz7/donor-screen-conda/actions/runs/45628034716) on the public example below (12 samples and 8,576 SNPs). This test runs the workflow from source; it is separate from Conda package validation.

## Public VCF example

Run `bash example-low-pass.sh example-run` after installation. This downloads a 1 Mb region on chromosome 22 at runtime from the public 1000 Genomes Phase 1 integrated callset (hg19). Phase 1 includes low-coverage whole-genome sequencing. The script uses 12 published sample IDs and arbitrary test roles and expects a successful one-offspring comparison. It demonstrates software execution; the short fragment is unsuitable for interpreting kinship or donor suitability. Genotype files remain in the run directory.
