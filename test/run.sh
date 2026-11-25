#!/usr/bin/env bash
set -euo pipefail
# Public Phase 1 integrated genotypes. Roles below are arbitrary software fixtures.
url=https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz
region=22:20000000-21000000
out="${1:?Usage: bash test/run.sh NEW_OUTPUT_DIRECTORY}"
tool="${DONOR_SCREEN:-donor-screen}"
here="$(cd -- "$(dirname -- "$0")" && pwd)"
mkdir -- "$out"
cd -- "$out"
cp "$here/input/roles.tsv" roles.tsv
if [[ -f "$here/input/fragment.vcf.gz" ]]; then
  cp "$here/input/fragment.vcf.gz" .
else
  bash "$here/fetch-input.sh"
fi
printf 'Samples: '; bcftools query -l fragment.vcf.gz | wc -l
printf 'SNPs: '; bcftools view -H fragment.vcf.gz | wc -l
"$tool" fragment.vcf.gz roles.tsv --out result --offspring 1 --threads 2 --java-memory 2g
python "$here/expected/check.py"
