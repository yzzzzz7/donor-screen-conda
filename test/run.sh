#!/usr/bin/env bash
set -euo pipefail
# Public Phase 1 integrated genotypes. Roles below are arbitrary software fixtures.
out="${1:?Usage: bash test/run.sh NEW_OUTPUT_DIRECTORY}"
tool="${DONOR_SCREEN:-donor-screen}"
here="$(cd -- "$(dirname -- "$0")" && pwd)"
mkdir -- "$out"
cd -- "$out"
cp "$here/input/roles.tsv" roles.tsv
cp "$here/input/fragment.vcf.gz" .
printf 'Samples: '; bcftools query -l fragment.vcf.gz | wc -l
printf 'SNPs: '; bcftools view -H fragment.vcf.gz | wc -l
"$tool" fragment.vcf.gz roles.tsv --out result --offspring 1 --threads 2 --java-memory 2g
python "$here/expected/check.py"
