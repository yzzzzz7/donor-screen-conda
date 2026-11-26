#!/usr/bin/env bash
set -euo pipefail
# Public Phase 1 integrated genotypes. Roles below are arbitrary software fixtures.
out="${1:?Usage: bash test/run.sh NEW_OUTPUT_DIRECTORY}"
tool="${DONOR_SCREEN:-donor-screen}"
here="$(cd -- "$(dirname -- "$0")" && pwd)"
mkdir -- "$out"
cd -- "$out"
cp "$here/input/roles.tsv" roles.tsv
cp "$here/input/donor-recipient.vcf.gz" .
printf 'Samples: '; bcftools query -l donor-recipient.vcf.gz | wc -l
printf 'SNPs: '; bcftools view -H donor-recipient.vcf.gz | wc -l
"$tool" donor-recipient.vcf.gz roles.tsv --out result --offspring 1 --threads 2 --java-memory 2g
python "$here/expected/check.py"
