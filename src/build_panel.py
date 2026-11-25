"""Build donor-screening inputs from a jointly aligned, QC-complete VCF (Python 3.10+).

Use the screening JSON with two additional explicit lists: reference_samples
(frequency denominator) and chromosomes (retained diploid autosomes). A site must
have complete calls for all reference individuals, recipients and candidates.
This program does not infer chromosome scope, call variants or estimate HBD.
"""
import argparse
import csv
import gzip
import json
import os
from contextlib import ExitStack
from itertools import combinations, product
from pathlib import Path
from tempfile import NamedTemporaryFile


DIPLOID = {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2,
           "0|0": 0, "0|1": 1, "1|0": 1, "1|1": 2}
MISSING = {".", "./.", ".|.", "./0", "./1", "0/.", "1/.", ".|0", ".|1", "0|.", "1|."}
PRESENT = bytes.maketrans(bytes((0, 1, 2)), bytes((0, 1, 1)))
DOUBLE = bytes.maketrans(bytes((0, 1, 2)), bytes((0, 0, 1)))


def names(values, label):
    if (not isinstance(values, list) or not values or
            not all(isinstance(x, str) and x and not any(ch.isspace() for ch in x) for x in values)
            or len(set(values)) != len(values)):
        raise ValueError(f"{label} needs a nonempty list of unique IDs")
    return values


class PairTotals:
    """Exact sums over SNPs using two bit arrays per individual for each chunk."""
    def __init__(self, n_rec, n_don):
        self.columns = [bytearray() for _ in range(n_rec + n_don)]
        self.pairs = list(combinations(range(n_rec), 2))
        self.pairs += list(product(range(n_rec), range(n_rec, n_rec + n_don)))
        self.twice_het = [0] * len(self.pairs)
        self.sites = 0

    def consume(self, genotypes):
        for column, gt in zip(self.columns, genotypes):
            column.append(gt)
        self.sites += 1
        if len(self.columns[0]) == 8192:
            self.flush()

    def flush(self):
        if not self.columns[0]:
            return
        bits = [(int.from_bytes(c.translate(PRESENT), "little"),
                 int.from_bytes(c.translate(DOUBLE), "little")) for c in self.columns]
        copies = [p.bit_count() + d.bit_count() for p, d in bits]
        for k, (i, j) in enumerate(self.pairs):
            # g_i*g_j = (I[g_i>=1]+I[g_i==2])*(I[g_j>=1]+I[g_j==2]).
            both = sum((a & b).bit_count() for a in bits[i] for b in bits[j])
            self.twice_het[k] += copies[i] + copies[j] - both
        for column in self.columns:
            column.clear()


def build_panel(config_path, vcf_path):
    cfg = json.loads(config_path.read_text(encoding="utf-8-sig"))
    reference = names(cfg["reference_samples"], "reference_samples")
    chroms = set(names(cfg["chromosomes"], "chromosomes"))
    individuals = cfg["individuals"]
    names([r["sample"] for r in individuals], "individuals")
    if any(r["role"] not in ("recipient", "candidate") for r in individuals):
        raise ValueError("Each individual needs role recipient or candidate")
    rec, don = (sorted(r["sample"] for r in individuals if r["role"] == role)
                for role in ("recipient", "candidate"))
    if len(rec) < 2 or not don or set(don) & {"chrom","pos","ref","alt","minor"}:
        raise ValueError("Need >=2 recipients, >=1 candidate, and unambiguous output column names")
    selected = rec + don
    needed = list(dict.fromkeys(reference + selected))
    idx = [needed.index(x) for x in selected]
    rare_path = config_path.parent / cfg["rare_file"]
    het_path = config_path.parent / cfg["het_file"]
    outputs = [rare_path, het_path]
    if rare_path.suffix not in (".gz", ".tsv") or het_path.suffix != ".tsv":
        raise ValueError("rare_file must be TSV (.tsv or .gz) and het_file must be plain .tsv")
    if len({p.resolve() for p in [config_path, vcf_path, *outputs]}) != 4:
        raise ValueError("Input and output paths must be distinct")
    for out in outputs:
        if out.exists() or out.is_symlink() or not out.parent.is_dir():
            raise FileExistsError(f"Use an existing directory and a new output filename: {out}")
    pair_counts = PairTotals(len(rec), len(don))
    seen_chroms = set()
    current_chrom = None
    previous_pos = 0
    at_position = set()
    missing_count = rare_count = total = 0
    with ExitStack() as cleanup:
        staged = []
        for output in outputs:
            with NamedTemporaryFile(dir=output.parent, prefix=output.name+".", suffix=".tmp", delete=False) as f:
                path = Path(f.name)
            cleanup.callback(path.unlink, missing_ok=True)
            staged.append(path)
        read = gzip.open if vcf_path.suffix in (".gz", ".bgz", ".bgzf") else open
        write = gzip.open if rare_path.suffix == ".gz" else open
        with read(vcf_path, "rt", encoding="utf-8-sig") as vcf, write(staged[0], "wt", encoding="utf-8", newline="") as f:
            sites = csv.writer(f, delimiter="\t", lineterminator="\n")
            sites.writerow(["chrom", "pos", "ref", "alt", "minor", *don])
            for line in vcf:
                if line.startswith("##"):
                    continue
                header = line.rstrip("\r\n").split("\t")
                if header[:9] != ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]:
                    raise ValueError("Missing VCF header or invalid fixed columns")
                break
            else:
                raise ValueError("Missing VCF header")
            samples = names(header[9:], "VCF sample IDs")
            absent = set(needed) - set(samples)
            if absent:
                raise ValueError(f"Samples absent from VCF: {sorted(absent)}")
            sample_cols = [samples.index(x)+9 for x in needed]
            for line in vcf:
                fields = line.rstrip("\r\n").split("\t")
                chrom = fields[0]
                if chrom not in chroms:
                    continue
                if len(fields) != len(header):
                    raise ValueError(f"Incorrect number of VCF fields: {fields[:2]}")
                pos, ref, alt = int(fields[1]), fields[3], fields[4]
                if ref not in {"A","C","G","T"} or alt not in {"A","C","G","T"} or ref == alt:
                    raise ValueError(f"Expected a biallelic SNP at {chrom}:{pos}")
                if current_chrom != chrom:
                    if chrom in seen_chroms:
                        raise ValueError(f"VCF is not grouped by chromosome: {chrom}")
                    seen_chroms.add(chrom)
                    current_chrom, previous_pos = chrom, 0
                if pos <= 0 or pos < previous_pos:
                    raise ValueError(f"Invalid VCF position/order: {chrom}:{pos}")
                if previous_pos != pos:
                    at_position.clear()
                if (ref,alt) in at_position:
                    raise ValueError(f"Duplicate VCF allele at {chrom}:{pos}")
                at_position.add((ref,alt))
                previous_pos = pos
                if fields[8].split(":",1)[0] != "GT":
                    raise ValueError(f"Expected standard GT-first VCF FORMAT at {chrom}:{pos}")
                total += 1
                g = []
                complete = True
                for col in sample_cols:
                    gt = fields[col].split(":",1)[0]
                    if gt in MISSING:
                        complete = False
                        g.append(None)
                    elif gt in DIPLOID:
                        g.append(DIPLOID[gt])
                    else:
                        raise ValueError(f"Expected diploid biallelic GT, got {gt} for {header[col]} at {chrom}:{pos}")
                if not complete:
                    missing_count += 1
                    continue
                calls = [g[i] for i in idx]
                pair_counts.consume(calls)
                an = 2 * len(reference)
                ac = sum(g[:len(reference)])
                mac = min(ac, an-ac)
                if mac == 0 or 20*mac > an:
                    continue
                minor_is_ref = 2*ac > an
                if minor_is_ref:
                    calls = [2-d for d in calls]
                if max(calls[:len(rec)]) == 0 and max(calls[len(rec):]) > 0:
                    sites.writerow([chrom,pos,ref,alt,ref if minor_is_ref else alt,*calls[len(rec):]])
                    rare_count += 1
        if seen_chroms != chroms or pair_counts.sites == 0:
            raise ValueError(f"No complete SNPs or missing requested chromosomes: {sorted(chroms-seen_chroms)}")
        pair_counts.flush()
        with staged[1].open("w", encoding="utf-8", newline="") as f:
            pairs = csv.writer(f, delimiter="\t", lineterminator="\n")
            pairs.writerow(["sample1","sample2","n_sites","expected_het"])
            for (i,j), count in zip(pair_counts.pairs,pair_counts.twice_het):
                pairs.writerow([selected[i],selected[j],pair_counts.sites,format(count/(2*pair_counts.sites),".17g")])
        # Same-directory hard links publish complete files without overwriting on a concurrent run.
        os.link(staged[0], outputs[0])
        try:
            os.link(staged[1], outputs[1])
        except BaseException:
            if outputs[0].exists() and outputs[0].samefile(staged[0]):
                outputs[0].unlink()
            raise
    return {"in_scope_snp_records": total, "incomplete_skipped": missing_count,
            "complete_snp_records": pair_counts.sites, "recipient_new_rare_sites": rare_count,
            "pairs": len(pair_counts.pairs), "recipients": rec, "candidates": don}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--vcf", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(build_panel(args.config, args.vcf), ensure_ascii=False, indent=2))
