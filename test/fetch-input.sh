#!/usr/bin/env bash
set -euo pipefail
url=https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz
region=22:20000000-21000000
cut -f1 roles.tsv | tail -n +5 > samples.txt
bcftools view --no-version -r "$region" -S samples.txt -m2 -M2 -v snps -Oz -o fragment.vcf.gz "$url"
rm -- samples.txt
