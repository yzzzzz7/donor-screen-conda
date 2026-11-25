"""Measure FROH in sequenced individuals with the study's PLINK 1.9 settings."""
import argparse
import csv
from pathlib import Path
import subprocess

from report import write_xlsx


COLUMNS = {
    "sample": "Sample",
    "FROH_100kb": "Total FROH (ROH at least 100 kb)",
    "FROH_gt1Mb": "Long FROH (ROH over 1 Mb)",
    "reference_bp": "Reference length (bp)",
}


def run(args):
    source = args.vcf.resolve(strict=True)
    samples = args.samples.resolve(strict=True)
    bed = args.bed.resolve(strict=True)
    names = samples.read_text(encoding="utf-8-sig").split()
    if not names or len(set(names)) != len(names):
        raise ValueError("Provide one unique sample ID per line")
    chromosomes = {}
    for line in bed.read_text(encoding="utf-8-sig").splitlines():
        if not line.strip():
            continue
        chrom, start, end, *_ = line.split()
        if chrom in chromosomes or int(start) != 0 or int(end) <= 0:
            raise ValueError("The BED must contain unique whole chromosomes starting at zero")
        chromosomes[chrom] = int(end)
    if not chromosomes:
        raise ValueError("The BED contains no chromosomes")
    if not 1 <= args.chr_set <= 95:
        raise ValueError("--chr-set must be within 1 to 95")
    version = subprocess.check_output(["plink", "--version"], text=True)
    if "PLINK v1.9" not in version:
        raise ValueError("Install PLINK 1.9 and expose it as plink")
    output = args.out.resolve()
    output.mkdir(exist_ok=False)
    (output / "plink-version.txt").write_text(version, encoding="utf-8")
    denominator = sum(chromosomes.values())
    plink_flags = ("--allow-extra-chr", "--allow-no-sex", "--chr-set", str(args.chr_set))

    with (output / "commands.log").open("w", encoding="utf-8") as log:
        def command(*tokens):
            command_line = [str(token) for token in tokens]
            print(repr(command_line), file=log, flush=True)
            subprocess.run(command_line, cwd=output, stdout=log, stderr=subprocess.STDOUT, check=True)

        command("bcftools", "view", "-S", samples, "-R", bed, "-m2", "-M2", "-v", "snps",
                source, "-Oz", "-o", "observed.vcf.gz")
        command("plink", "--vcf", "observed.vcf.gz", "--double-id", *plink_flags,
                "--make-bed", "--out", "observed")
        observed = [line.split()[1] for line in (output / "observed.fam").read_text().splitlines()]
        if len(observed) != len(names) or set(observed) != set(names):
            raise ValueError("VCF samples do not match the requested samples")
        command("plink", "--bfile", "observed", *plink_flags, "--homozyg",
                "--homozyg-kb", "100", "--homozyg-snp", "20", "--homozyg-window-snp", "20",
                "--homozyg-density", "50", "--homozyg-gap", "1000", "--homozyg-window-het", "1",
                "--homozyg-window-missing", "5", "--homozyg-window-threshold", "0.05", "--out", "roh")
    lengths = {sample: [] for sample in names}
    with (output / "roh.hom").open() as segments:
        header = segments.readline().split()
        if not {"KB", "IID"}.issubset(header):
            raise ValueError("Unexpected PLINK .hom columns")
        for line in segments:
            row = dict(zip(header, line.split(), strict=True))
            kb = float(row["KB"])
            if kb < 100:
                raise ValueError("PLINK returned a segment below the 100 kb threshold")
            lengths[row["IID"]].append(kb)
    records = [
        dict(sample=sample, FROH_100kb=sum(kb) * 1000 / denominator,
             FROH_gt1Mb=sum(x for x in kb if x > 1000) * 1000 / denominator,
             reference_bp=denominator)
        for sample, kb in lengths.items()
    ]
    with (output / "froh.tsv").open("w", encoding="utf-8", newline="") as handle:
        table = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        table.writeheader()
        table.writerows(records)
    write_xlsx(output / "froh.xlsx", COLUMNS, records)
    print(output / "froh.xlsx")


def main(argv):
    parser = argparse.ArgumentParser(prog="donor-screen froh", description=__doc__)
    parser.add_argument("vcf", type=Path, help="Indexed, filtered VCF with observed genotypes")
    parser.add_argument("samples", type=Path, help="One sequenced sample ID per line")
    parser.add_argument("bed", type=Path, help="Whole-autosome BED for the reference denominator")
    parser.add_argument("--chr-set", type=int, required=True, help="PLINK autosome count for the reference")
    parser.add_argument("--out", type=Path, required=True)
    try:
        run(parser.parse_args(argv))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"donor-screen froh: {error}\n")
