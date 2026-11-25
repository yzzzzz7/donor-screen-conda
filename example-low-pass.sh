#!/usr/bin/env bash
set -euo pipefail
# Public Phase 1 integrated genotypes. Roles below are arbitrary software fixtures.
url=https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz
region=22:20000000-21000000
out="${1:?Usage: bash example-low-pass.sh NEW_OUTPUT_DIRECTORY}"
tool="${DONOR_SCREEN:-donor-screen}"
mkdir -- "$out"
cd -- "$out"
cat > roles.tsv <<'EOF'
#autosomes=22
#reference_bp=1000001
#king_panel=joint
sample	role
NA06984	recipient
NA06986	recipient
NA07000	recipient
NA07048	recipient
NA07051	recipient
NA07347	recipient
NA12812	donor
NA12814	donor
NA12815	donor
NA10847	reference
NA07357	reference
NA10851	reference
EOF
cut -f1 roles.tsv | tail -n +5 > samples.txt
bcftools view -r "$region" -S samples.txt -m2 -M2 -v snps -Oz -o fragment.vcf.gz "$url"
printf 'Samples: '; bcftools query -l fragment.vcf.gz | wc -l
printf 'SNPs: '; bcftools view -H fragment.vcf.gz | wc -l
"$tool" fragment.vcf.gz roles.tsv --out result --offspring 1 --threads 2 --java-memory 2g
python - <<'PY'
import json
from pathlib import Path
run = Path("result")
report = json.loads((run / "results.json").read_text())
assert report["cohort"]["offspring_count"] == 1
assert report["cohort"]["het_sites"] > 0
assert [r["k"] for r in report["results"]] == [0, 1]
assert len(report["results"][1]["donors"]) == 1
assert (run / "results.csv").stat().st_size > 0
print("Full pipeline smoke test passed")
PY
