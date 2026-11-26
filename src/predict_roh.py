"""Simulate offspring genotypes and estimate their FROH with the observed ROH caller."""
from argparse import Namespace
import csv
import math
from pathlib import Path
from random import Random
from statistics import mean, quantiles
import subprocess
from tempfile import TemporaryDirectory

from observed_roh import run as observed_roh
from report import html_table, write_xlsx

COLUMNS = {
    "plan": "Plan", "k": "Donor count", "donors": "Donors",
    "rare": "Expected new rare alleles", "het": "Predicted offspring heterozygosity",
    "hbd": "Predicted recent HBD", "r": "Predicted R",
    "total": "Mean total FROH (≥100 kb)", "long": "Mean long FROH (>1 Mb)",
    "total_low": "Total FROH 2.5th percentile", "total_high": "Total FROH 97.5th percentile",
    "long_low": "Long FROH 2.5th percentile", "long_high": "Long FROH 97.5th percentile",
}


def pair_key(pair):
    return tuple(sorted((pair["parent1"], pair["parent2"])))


def plans_for(rows):
    for row in rows:
        if row["pairs"]:
            yield dict(row, plan="Selected" if row["k"] else "No introduction")
        if "reference_pairs" in row:
            yield dict(plan="High-heterozygosity reference", k=row["k"],
                       donors=row["reference_donors"], pairs=row["reference_pairs"],
                       rare=row["reference_rare"], het=row["reference_het"],
                       hbd=row["cap"], r=row.get("reference_r"))


def simulate(genotypes, names, pairs, count, bed, out, seed):
    rng = Random(seed)
    parents = [(names.index(a), names.index(b)) for a, b in pairs for _ in range(count)]
    samples = [f"F1_{i + 1:08d}" for i in range(len(parents))]
    ref_lengths = {chrom: int(end) for chrom, start, end, *_ in (
        line.split() for line in bed.read_text(encoding="utf-8-sig").splitlines() if line.strip()
    )}
    prev_chrom = None
    prev_pos = 0
    observed_chromosomes = set()
    with (out / "simulate.log").open("w", encoding="utf-8") as log:
        with subprocess.Popen(["bcftools", "view", "-Oz", "-o", str(out / "simulated.vcf.gz")],
                              stdin=subprocess.PIPE, stderr=log, text=True, encoding="utf-8") as process:
            stream = process.stdin
            print("##fileformat=VCFv4.2", file=stream)
            for chrom, length in ref_lengths.items():
                print(f"##contig=<ID={chrom},length={length}>", file=stream)
            print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=stream)
            print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                  *samples, sep="\t", file=stream)
            for (chrom, pos, ref, alt), calls in genotypes:
                bp = int(pos)
                if chrom not in ref_lengths or not 1 <= bp <= ref_lengths[chrom]:
                    raise ValueError(f"ROH BED does not cover site {chrom}:{pos}")
                if chrom != prev_chrom:
                    if chrom in observed_chromosomes:
                        raise ValueError("VCF chromosome records must be contiguous")
                    observed_chromosomes.add(chrom)
                    maternal = [rng.randrange(2) for _ in parents]
                    paternal = [rng.randrange(2) for _ in parents]
                else:
                    if bp < prev_pos:
                        raise ValueError("VCF positions must be sorted")
                    # Haldane mapping gives the probability of an odd number of crossovers.
                    switch = -math.expm1(-2 * (bp - prev_pos) / 100_000_000) / 2
                    maternal = [state ^ (rng.random() < switch) for state in maternal]
                    paternal = [state ^ (rng.random() < switch) for state in paternal]
                diplotypes = (calls[a][2 * first] + "/" + calls[b][2 * second]
                             for (a, b), first, second in zip(parents, maternal, paternal))
                print(chrom, pos, f"{chrom}:{pos}:{ref}:{alt}", ref, alt, ".", "PASS", ".", "GT",
                      *diplotypes, sep="\t", file=stream)
                prev_chrom, prev_pos = chrom, bp
            stream.close()
            if process.wait():
                raise ValueError("bcftools could not write simulated genotypes")
    if observed_chromosomes != set(ref_lengths):
        raise ValueError("The ROH BED must cover exactly the analyzed chromosomes")
    (out / "samples.txt").write_text("\n".join(samples) + "\n", encoding="utf-8")
    subprocess.run(["bcftools", "index", str(out / "simulated.vcf.gz")], check=True)
    return samples


def predict(data, genotypes, names, bed, chr_set, out, count, seed):
    plans = list(plans_for(data["results"]))
    pairs = sorted({pair_key(pair) for plan in plans for pair in plan["pairs"]})
    out.mkdir()
    with TemporaryDirectory(dir=out, prefix="simulation-") as tmp:
        work = Path(tmp)
        samples = simulate(genotypes, names, pairs, count, bed, work, seed)
        called = work / "called"
        observed_roh(Namespace(vcf=work / "simulated.vcf.gz", samples=work / "samples.txt",
                               bed=bed, chr_set=chr_set, out=called), announce=False)
        with (called / "froh.tsv").open(encoding="utf-8", newline="") as source:
            values = {row["sample"]: (float(row["FROH_100kb"]), float(row["FROH_gt1Mb"]))
                      for row in csv.DictReader(source, delimiter="\t")}
    lookup = {pair: [values[s] for s in samples[i * count:(i + 1) * count]]
              for i, pair in enumerate(pairs)}
    with (out / "simulated-f1-froh.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("parent1", "parent2", "replicate", "FROH_100kb", "FROH_gt1Mb"))
        for pair, rows in lookup.items():
            for i, row in enumerate(rows, 1):
                writer.writerow((*pair, i, *row))
    summary = []
    for plan in plans:
        matings = [lookup[pair_key(pair)] for pair in plan["pairs"]]
        totals, longs = zip(*(tuple(map(mean, zip(*replicate))) for replicate in zip(*matings)))
        low, high = (quantiles(values, n=40, method="inclusive") for values in (totals, longs))
        summary.append(dict(plan, total=mean(totals), long=mean(longs),
                            total_low=low[0], total_high=low[-1],
                            long_low=high[0], long_high=high[-1]))
    write_xlsx(out / "results.xlsx", COLUMNS, summary)
    return summary


def report_section(rows, replicates, seed):
    return (
        '<section><h2>Predicted offspring FROH</h2>'
        f'<p>{replicates} simulated offspring per pair, seed {seed}. Parental haplotypes '
        'were recombined at a constant 1 cM/Mb and scored using the same PLINK ROH '
        'settings as observed samples. The 2.5th and 97.5th percentiles describe variation '
        'among simulated plans. FROH did not enter the ranking.</p>'
        '<p><a href="04_predicted_froh/results.xlsx">Excel</a> · '
        '<a href="04_predicted_froh/simulated-f1-froh.tsv">Simulated F1 results</a></p>'
        + html_table(COLUMNS, rows) + '</section>'
    )
