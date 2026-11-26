#!/usr/bin/env bash
set -euo pipefail
url=https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz
region=22:20000000-21000000
here="$(cd -- "$(dirname -- "$0")" && pwd)"
out="${1:?Usage: bash test/fetch-input.sh OUTPUT_VCF.gz}"
names="$here/input/sample-names.tsv"
bcftools view --no-version -r "$region" -S <(tail -n +2 "$names" | cut -f1) -m2 -M2 -v snps -Oz "$url" |
  bcftools reheader --samples <(tail -n +2 "$names") --output "$out"
