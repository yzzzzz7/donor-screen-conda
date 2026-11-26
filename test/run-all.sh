#!/usr/bin/env bash
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
out="${1:?Usage: bash test/run-all.sh NEW_OUTPUT_DIRECTORY}"
tool="${DONOR_SCREEN:-donor-screen}"
bash "$here/run.sh" "$out"
"$tool" "$out/fragment.vcf.gz" "$out/roles.tsv" --out "$out/r-result" \
  --offspring 1 --threads 2 --java-memory 2g --cds-sites "$here/input/synthetic-cds.tsv" \
  --froh-bed "$here/input/autosomes.bed" --froh-chr-set 22 --froh-replicates 16
bcftools index "$out/fragment.vcf.gz"
"$tool" froh "$out/fragment.vcf.gz" "$here/input/observed-samples.txt" "$here/input/autosomes.bed" \
  --chr-set 22 --out "$out/observed-roh"
python "$here/observed-report.py" "$out/observed-roh"
