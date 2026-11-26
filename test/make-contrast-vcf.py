"""Generate a synthetic VCF to demonstrate donor-ranking trade-offs."""
import random
import sys


SEED = 170621
NAMES = [f"Recipient_{n:02d}" for n in range(1, 7)] + [
    f"Donor_{n:02d}" for n in range(1, 4)
] + [f"Reference_{n:02d}" for n in range(1, 4)]
DONOR_START = 6
CHROMOSOME_BP = 4_000_000
SITE_COUNTS = {"background": 2300, "donor_1": 160, "donor_2": 90, "donor_3": 760}


def dosage(rng, p):
    return sum(rng.random() < p for _ in range(2))


def emit_vcf(output):
    rng = random.Random(SEED)
    classes = [kind for kind, count in SITE_COUNTS.items() for _ in range(count)]
    rng.shuffle(classes)
    print("##fileformat=VCFv4.2", file=output)
    print(f"##contig=<ID=1,length={CHROMOSOME_BP}>", file=output)
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=output)
    print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
          *NAMES, sep="\t", file=output)
    for index, kind in enumerate(classes, 1):
        if kind == "background":
            p = rng.uniform(0.22, 0.46)
            while True:
                calls = [dosage(rng, p) for _ in NAMES]
                for j in (DONOR_START, DONOR_START + 1):
                    calls[j] = dosage(rng, 0.50)
                # Give this donor lower heterozygosity and more contrasting genotypes.
                calls[DONOR_START + 2] = dosage(rng, 0.94)
                if 0 < sum(calls) < 2 * len(calls):
                    break
        else:
            calls = [0] * len(NAMES)
            calls[DONOR_START + int(kind[-1]) - 1] = 1
        position = index * (CHROMOSOME_BP // (len(classes) + 1))
        genotypes = [("0/0", "0/1", "1/1")[i] for i in calls]
        print("1", position, f"synthetic_{index}", "A", "G", ".", "PASS", ".", "GT",
              *genotypes, sep="\t", file=output)


if __name__ == "__main__":
    emit_vcf(sys.stdout)
