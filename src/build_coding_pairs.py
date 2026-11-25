"""Summarize a complete, allele-matched CDS panel into per-mating expected counts."""
import argparse
import csv
import gzip
import json
from itertools import combinations,product
from pathlib import Path
from screen import output_file


def build(sites, config, output):
    if output.exists():
        raise FileExistsError(output)
    cfg=json.loads(config.read_text(encoding="utf-8-sig"))
    rec,don=([r["sample"] for r in cfg["individuals"] if r["role"]==role]
             for role in ("recipient","candidate"))
    people=rec+don
    if len(people)!=len(set(people)) or not rec or not don:
        raise ValueError("Invalid recipients/candidates")
    coords=list(combinations(range(len(rec)),2))+list(product(range(len(rec)),range(len(rec),len(people))))
    totals=[[0,0] for _ in coords]
    seen=set()
    opener=gzip.open if sites.suffix==".gz" else open
    with opener(sites,"rt",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        expected={"chrom","pos","ref","alt","missense",*people}
        if not expected.issubset(reader.fieldnames or []):
            raise ValueError("CDS panel requires allele identity, curated missense flag and ALT dosages")
        for r in reader:
            key=(r["chrom"],int(r["pos"]),r["ref"],r["alt"])
            if key in seen or key[0] not in cfg["chromosomes"] or key[1]<1 or key[2]==key[3]:
                raise ValueError(f"Invalid/duplicate CDS allele: {key}")
            seen.add(key)
            gt=[int(r[x]) for x in people]
            if any(x not in (0,1,2) for x in gt) or r["missense"] not in ("0","1"):
                raise ValueError(f"Missing or invalid genotype/annotation: {key}")
            flag=int(r["missense"])
            for idx,(i,j) in enumerate(coords):
                copies=gt[i]*gt[j]
                totals[idx][0]+=flag*copies
                totals[idx][1]+=copies
    if any(c<=0 for _,c in totals):
        raise ValueError("At least one pair has an undefined R (zero expected CDS count)")
    with output_file(output,"utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(["sample1","sample2","expected_missense","expected_cds"])
        for (i,j),(m,c) in zip(coords,totals):
            w.writerow([people[i],people[j],m/4,c/4])




def coding_panel(sites, genotypes, people, output):
    """Attach observed ALT dosages to a curated set of CDS alleles."""
    opener = gzip.open if sites.suffix == ".gz" else open
    flags = {}
    with opener(sites, "rt", encoding="utf-8-sig", newline="") as source:
        rows = csv.DictReader(source, delimiter="\t")
        fields = ("chrom", "pos", "ref", "alt", "missense")
        if not set(fields).issubset(rows.fieldnames or []):
            raise ValueError("CDS annotation requires chrom, pos, ref, alt and missense")
        for row in rows:
            key = (row["chrom"], str(int(row["pos"])), row["ref"], row["alt"])
            if key in flags or row["missense"] not in ("0", "1"):
                raise ValueError(f"Invalid or duplicate CDS annotation: {key}")
            flags[key] = row["missense"]
    if not flags:
        raise ValueError("The CDS annotation is empty")
    with output_file(output, "utf-8", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow((*fields, *people))
        matched = set()
        for site, gt in genotypes:
            if site not in flags:
                continue
            if site in matched:
                raise ValueError(f"Duplicate CDS allele in VCF: {site}")
            matched.add(site)
            writer.writerow((*site, flags[site], *gt))
        if len(matched) != len(flags):
            example = next(key for key in flags if key not in matched)
            raise ValueError(f"CDS annotation has {len(flags) - len(matched)} unmatched VCF alleles, including {example}")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("sites","config","output"):p.add_argument(arg,type=Path)
    args=p.parse_args()
    build(args.sites,args.config,args.output)
