"""Generate the single artificial donor-screen fixture and its CDS labels."""
import argparse
import random
import sys


SEED = 170621
RECIPIENTS = [f"Recipient_{n:02d}" for n in range(1, 7)]
DONORS = [f"Donor_{n:02d}" for n in range(1, 4)]
REFERENCES = [f"Reference_{n:02d}" for n in range(1, 4)]
NAMES = RECIPIENTS + DONORS + REFERENCES
DONOR_START = len(RECIPIENTS)
REFERENCE_START = DONOR_START + len(DONORS)
CHROMOSOME_BP = 32_000_000
PRIVATE_SITES = (160, 90, 760)
LONG_BLOCK_SAMPLES = {"Recipient_01", "Recipient_02", "Recipient_03", "Donor_01", "Donor_03"}
SHORT_BLOCK_SAMPLES = {"Recipient_04", "Recipient_05", "Donor_02"}


def dosage(rng, p):
    return sum(rng.random() < p for _ in range(2))


def synthetic_sites():
    rng = random.Random(SEED)
    block_rng = random.Random(SEED + 1)
    classes = ["background"] * 2300
    for donor, count in enumerate(PRIVATE_SITES, 1):
        classes.extend([f"donor_{donor}"] * count)
    rng.shuffle(classes)
    # Keep the original donor-ranking contrast, then add lower-heterozygosity markers.
    backgrounds = []
    for kind in classes:
        if kind == "background":
            p = rng.uniform(0.22, 0.46)
            while True:
                calls = [dosage(rng, p) for _ in NAMES]
                for j in (DONOR_START, DONOR_START + 1):
                    calls[j] = dosage(rng, 0.50)
                calls[DONOR_START + 2] = dosage(rng, 0.94)
                if 0 < sum(calls) < 2 * len(calls):
                    break
        else:
            calls = [0] * len(NAMES)
            calls[DONOR_START + int(kind[-1]) - 1] = 1
        backgrounds.append(calls)
    for _ in range(4900):
        calls = [0] * len(NAMES)
        for j in block_rng.sample(range(REFERENCE_START, len(NAMES)), 2):
            calls[j] = 1
        backgrounds.append(calls)
    block_rng.shuffle(backgrounds)

    def homozygous_block(count, target):
        for _ in range(count):
            calls = [0 if name in target else dosage(block_rng, 0.34) for name in NAMES]
            # Keep every marker polymorphic without interrupting the intended tract.
            calls[REFERENCE_START] = 1
            yield calls

    block_long = list(homozygous_block(360, LONG_BLOCK_SAMPLES))
    block_short = list(homozygous_block(110, SHORT_BLOCK_SAMPLES))
    calls = backgrounds[:2150] + block_long + backgrounds[2150:5700] + block_short + backgrounds[5700:]
    step = CHROMOSOME_BP // (len(calls) + 1)
    assert (len(block_long) - 1) * step > 1_000_000
    for index, sample_dosages in enumerate(calls, 1):
        yield index, index * step, sample_dosages


def emit_vcf(output):
    print("##fileformat=VCFv4.2", file=output)
    print(f"##contig=<ID=1,length={CHROMOSOME_BP}>", file=output)
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=output)
    print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
          *NAMES, sep="\t", file=output)
    for index, position, calls in synthetic_sites():
        print("1", position, f"synthetic_{index}", "A", "G", ".", "PASS", ".", "GT",
              *(("0/0", "0/1", "1/1")[i] for i in calls), sep="\t", file=output)


def emit_cds(output):
    print("chrom", "pos", "ref", "alt", "missense", sep="\t", file=output)
    for index, position, _ in synthetic_sites():
        print("1", position, "A", "G", index % 2, sep="\t", file=output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cds", action="store_true", help="Emit the synthetic CDS labels instead of a VCF")
    arguments = parser.parse_args()
    (emit_cds if arguments.cds else emit_vcf)(sys.stdout)
