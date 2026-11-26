"""Recipient-specific donor screening, v0.1. Python 3.9+; standard library only."""
import argparse
import csv
import gzip
import json
import math
import os
from collections import Counter
from contextlib import contextmanager
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from tempfile import NamedTemporaryFile


@contextmanager
def output_file(path,encoding,newline=None):
    with NamedTemporaryFile("w",dir=path.parent,encoding=encoding,newline=newline,delete=False) as file:
        tmp = Path(file.name)
        try:
            yield file
            file.close()
            os.replace(tmp,path)
        finally:
            file.close()
            tmp.unlink(missing_ok=True)


def pair_key(a, b):
    return tuple(sorted((a, b)))


def prob(value):
    value = float(value)
    if not 0 <= value <= 1:
        raise ValueError("Probability outside [0,1]")
    return value


def load_panel(path, donors):
    patterns = Counter()
    seen = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as f:
        rows = csv.DictReader(f, delimiter="\t")
        required = {"chrom", "pos", "ref", "alt", "minor", *donors}
        if not required.issubset(rows.fieldnames or []):
            raise ValueError("Missing allele or sample columns")
        for r in rows:
            key = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
            if key in seen or key[1] < 1 or r["minor"] not in key[2:]:
                raise ValueError(f"Invalid or duplicate allele: {key}")
            seen.add(key)
            gt = tuple(int(r[d]) for d in donors)
            if not all(x in (0,1,2) for x in gt):
                raise ValueError("Minor-allele dosages must be complete diploid calls")
            patterns[gt] += 1
    return patterns


def read_pairs(cfg, root):
    by_id = {}
    for r in cfg["pairs"]:
        key = pair_key(r["a"],r["b"])
        if key in by_id or key[0] == key[1] or not isinstance(r["allowed"],bool):
            raise ValueError(f"Invalid pair {key}")
        by_id[key] = {"allowed":r["allowed"], "hbd":prob(r["hbd"])}
    if cfg.get("export_pairwise"):
        for r in cfg["pairs"]:
            if "king" in r:
                by_id[pair_key(r["a"], r["b"])]["king"] = float(r["king"])
    sites = set()
    with (root / cfg["het_file"]).open(encoding="utf-8-sig") as f:
        for row in csv.DictReader(f,delimiter="\t"):
            key = pair_key(row["sample1"],row["sample2"])
            if key in by_id:
                if "het" in by_id[key]:
                    raise ValueError(f"Duplicate heterozygosity for {key}")
                by_id[key]["het"] = prob(row["expected_het"])
                sites.add(int(row["n_sites"]))
    if len(sites) != 1 or min(sites) <= 0:
        raise ValueError("A common, nonempty callable SNP panel is required")
    return by_id,sites.pop()


def best_baseline(recipients, n, pairs):
    edges = [p for p in combinations(recipients,2) if pairs[p]["allowed"]]
    return min((math.fsum(pairs[e]["hbd"] for e in lineup)/n,lineup)
               for lineup in combinations(edges,n) if len({s for e in lineup for s in e}) == 2*n)


def match_donors(slots, recipients, pairs):
    # DP over donor slots. Each recipient may be used once. Cost ties use sample IDs.
    states = {0:(0.0,tuple([None]*len(slots)))}
    for a in recipients:
        next_states = states.copy()
        for mask,(_,assignment) in states.items():
            for j,d in enumerate(slots):
                bit = 1 << j
                if mask & bit or not pairs[pair_key(a,d)]["allowed"]:
                    continue
                full = assignment[:j] + ((d,a),) + assignment[j+1:]
                cost = math.fsum(pairs[pair_key(*e)]["hbd"] for e in full if e is not None)
                option = (cost,full)
                if mask|bit not in next_states or option < next_states[mask|bit]:
                    next_states[mask|bit] = option
        states = next_states
    return states.get((1 << len(slots))-1)


def run(config_path, *, missense_policy=None):
    cfg = json.loads(config_path.read_text(encoding="utf-8-sig"))
    people = cfg["individuals"]
    info = {r["sample"]:r for r in people}
    if len(info) != len(people) or any(r["role"] not in ("candidate","recipient") for r in people):
        raise ValueError("Each individual needs a unique ID and a single role")
    donors = sorted(s for s in info if info[s]["role"] == "candidate")
    recipients = sorted(s for s in info if info[s]["role"] == "recipient")
    n, kmax = cfg["offspring_count"],cfg["max_donors"]
    exact = cfg.get("donor_count")
    if exact is not None and (type(exact) is not int or not 1 <= exact <= kmax):
        raise ValueError("Invalid exact donor count")
    if not 1 <= kmax <= n <= 3 or len(recipients) < 2*n or len(donors)<kmax:
        raise ValueError("v0.1 requires 1–3 F1, enough candidates and at least 2 recipients per F1")
    for s in donors:
        if not info[s]["cluster"]:
            raise ValueError(f"Missing donor cluster for {s}")
        prob(info[s]["h"])
    pairs,nsites = read_pairs(cfg,config_path.parent)
    r_file = cfg.get("missense_file")
    policy = missense_policy if missense_policy is not None else cfg.get("missense_policy","report")
    if policy not in ("reference","baseline","report") or (("missense_policy" in cfg or missense_policy is not None) and not r_file):
        raise ValueError("Specify missense_file and a valid missense_policy")
    if r_file:
        from missense_screen import attach_pair_counts, annotate, compare, enrich_group_options, minimum_r_for_original_score, decision_intervals
        attach_pair_counts(config_path.parent/r_file,recipients,donors,pairs)
    for a,b in (*combinations(recipients,2),*((a,d) for a in recipients for d in donors)):
        if "het" not in pairs[pair_key(a,b)]:
            raise ValueError(f"Missing heterozygosity for {a}, {b}")
    patterns = load_panel(config_path.parent / cfg["rare_file"],donors)
    for j,d in enumerate(donors):
        info[d]["novel_carried"] = sum(c for gt,c in patterns.items() if gt[j])
    for a in recipients:
        info[a]["novel_carried"] = 0
    cohort = dict(name=cfg["cohort"],offspring_count=n,het_sites=nsites,rare_sites=sum(patterns.values()),
                  candidates=donors,recipients=recipients,notes=cfg.get("notes",[]))

    def metrics(lineup):
        details = []
        for a,b in lineup:
            p = pairs[pair_key(a,b)]
            details.append(dict(parent1=a,parent2=b,n_offspring=1,hbd=p["hbd"],het=p["het"]))
        return dict(pairs=details,hbd=math.fsum(p["hbd"] for p in details)/n,
                    het=math.fsum(p["het"] for p in details)/n)

    try:
        _, base = best_baseline(recipients,n,pairs)
    except ValueError as e:
        raise ValueError("No recipient-only baseline satisfies the supplied pairing limits") from e
    baseline = dict(k=0, donors=[],rare=0,cap=None,reference_donors=[],
                    reason="The no-introduction control uses the feasible recipient-only pairings with the lowest recent HBD.",**metrics(base))
    if r_file:
        annotate(baseline,pairs)
        compare(baseline,baseline.copy(),None,policy)
        cohort["missense_policy"] = policy
    rows = [baseline]
    alternatives = []
    donor_counts = (exact,) if exact is not None else range(1, kmax + 1)
    for k in donor_counts:
        groups = list(combinations(donors,k))
        by_group = {}
        for group in groups:
            if len({info[d]["cluster"] for d in group}) < k:
                continue
            options = []
            for slots in combinations_with_replacement(group,n):
                if len(set(slots)) != k:
                    continue
                match = match_donors(slots,recipients,pairs)
                if match is None:
                    continue
                idx = tuple(donors.index(d) for d in slots)
                missed = sum(c * math.prod(2-gt[j] for j in idx) for gt,c in patterns.items())
                rare = sum(patterns.values()) - missed/(2**n)
                options.append(dict(k=k,donors=list(group),rare=rare,**metrics(match[1])))
            if options:
                by_group[group] = options
        def order(r):
            return (-r["rare"],r["hbd"],r["donors"],[(p["parent1"],p["parent2"]) for p in r["pairs"]])
        if not by_group:
            rows.append(dict(k=k,donors=[],pairs=[],hbd=None,het=None,rare=None,cap=None,
                             eligible_count=0,possible_count=len(groups),
                             reason="No combination meets the donor-cluster and pairing criteria."))
            alternatives.extend(dict(k=k,donors=list(g),cap=None,rare=None,hbd=None,het=None,
                                     **{"pass":False}) for g in groups)
            if r_file:
                rows[-1]["r"] = None
                compare(rows[-1],rows[-1].copy(),None,policy)
            continue
        min_hbd = {g:min(r["hbd"] for r in opts) for g,opts in by_group.items()}
        reference = min(by_group,key=lambda g:(-math.fsum(info[d]["h"] for d in g),min_hbd[g],g))
        cap = min_hbd[reference]
        admitted = []
        for g in groups:
            opts = by_group.get(g,[])
            allowed = [r for r in opts if r["hbd"] <= cap]
            rec = dict(k=k,donors=list(g),cap=cap,**{"pass":bool(allowed)})
            if allowed:
                best = min(allowed,key=order)
                admitted.append(best)
                rec.update(rare=best["rare"],hbd=best["hbd"],het=best["het"])
            else:
                rec.update(rare=None,hbd=min_hbd.get(g),het=None)
            alternatives.append(rec)
        admitted.sort(key=order)
        winner = admitted[0].copy()
        runner = admitted[1] if len(admitted)>1 else None
        reference_plan = min(by_group[reference],key=lambda r:(r["hbd"],order(r)))
        winner.update(cap=cap,eligible_count=len(admitted),possible_count=len(groups),reference_donors=list(reference),
                      reference_pairs=reference_plan["pairs"],
                      reference_rare=reference_plan["rare"],reference_het=reference_plan["het"],
                      gain_over_reference=winner["rare"]-reference_plan["rare"],
                      runner_up=runner["donors"] if runner else [],runner_gap=winner["rare"]-runner["rare"] if runner else None)
        winner["reason"] = f"Of {len(groups)} combinations, {len(admitted)} qualified. The selected plan maximizes expected distinct new rare alleles."
        if runner and winner["runner_gap"]==0:
            winner["reason"] += " Ties are resolved by lower recent HBD, then sample identifiers."
        if r_file:
            annotate(winner,pairs)
            annotate(reference_plan,pairs)
            old = winner.copy()
            r_cap = {"reference":reference_plan["r"],"baseline":baseline["r"],"report":None}[policy]
            admitted4 = []
            for rec in alternatives[-len(groups):]:
                key = tuple(rec["donors"])
                opts = by_group.get(key,[])
                if policy=="report":
                    valid = [annotate(o,pairs) for o in opts if o["hbd"] <= cap]
                else:
                    valid = enrich_group_options(opts,recipients,pairs,cap,r_cap,metrics)
                rec["r_cap"] = r_cap
                rec["pass3"] = rec["pass"]
                rec["pass"] = bool(valid)
                if valid:
                    chosen = min(valid,key=order)
                    admitted4.append(chosen)
                    rec.update({key:chosen[key] for key in ("rare","hbd","het","r")})
                else:
                    rec.update(rare=None,het=None,r=None)
            admitted4.sort(key=order)
            if admitted4:
                winner=admitted4[0].copy()
                runner=admitted4[1] if len(admitted4)>1 else None
                winner.update(runner_up=runner["donors"] if runner else [],
                              runner_gap=winner["rare"]-runner["rare"] if runner else None,
                              gain_over_reference=winner["rare"]-reference_plan["rare"])
            else:
                winner=dict(k=k,donors=[],pairs=[],rare=None,hbd=None,het=None,r=None)
            winner.update(cap=cap,eligible_count=len(admitted4),possible_count=len(groups),
                          reference_donors=list(reference),reference_pairs=reference_plan["pairs"],reference_rare=reference_plan["rare"],
                          reference_het=reference_plan["het"],reference_r=reference_plan["r"])
            compare(winner,old,r_cap,policy)
            winner["min_r_for_original_score"]=minimum_r_for_original_score(
                by_group[tuple(old["donors"])],old,recipients,pairs,cap)
            if policy=="report":
                winner["r_intervals"],winner["conditional"]=decision_intervals(by_group,old,recipients,pairs,cap,
                    {"reference":reference_plan["r"],"baseline":baseline["r"]},metrics)
            winner["reason"]=winner["change"]
        rows.append(winner)
    result = dict(version="0.2" if r_file else "0.1",cohort=cohort,individuals=people,results=rows,alternatives=alternatives)
    if cfg.get("export_pairwise"):
        gain = {d: sum(c * gt[j] for gt, c in patterns.items()) / 2
                for j, d in enumerate(donors)}
        result["pairwise"] = [
            dict(recipient=a, donor=d, het=pairs[pair_key(a, d)]["het"],
                 rare=gain[d], hbd=pairs[pair_key(a, d)]["hbd"] if pairs[pair_key(a, d)]["allowed"] else None,
                 king=pairs[pair_key(a, d)].get("king"), r=pairs[pair_key(a, d)].get("r"), allowed=pairs[pair_key(a, d)]["allowed"])
            for a in recipients for d in donors
        ]
    return result


def main():
    folder = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Recipient-specific donor screening for genetic rescue")
    default_config=folder/"20261112_R_inputs.json"
    if not default_config.exists():
        default_config=folder/"20261112_Snow_inputs.json"
    parser.add_argument("config",nargs="?",type=Path,default=default_config)
    parser.add_argument("output",nargs="?",type=Path,default=folder/"results.json")
    parser.add_argument("--csv",type=Path,help="Write a CSV summary")
    parser.add_argument("--r-policy",choices=["report","reference","baseline"],help="Override the configured R screening policy")
    args = parser.parse_args()
    paths = [p.resolve() for p in (args.config,args.output,args.csv) if p is not None]
    if len(set(paths)) != len(paths):
        parser.error("Input and output paths must differ.")
    result = run(args.config,missense_policy=args.r_policy)
    with output_file(args.output,"utf-8") as f:
        json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
    if args.csv:
        with output_file(args.csv,"utf-8-sig",newline="") as f:
            writer=csv.writer(f)
            writer.writerow(["Donor count","Selected donors","Expected new rare alleles","Predicted offspring SNP heterozygosity","Predicted recent HBD","Predicted homozygous missense proportion R","Original three-metric selection","Selection with high-heterozygosity R limit","Selection with no-introduction R limit","Rationale"])
            for r in result["results"]:
                cases = [("+".join(r["conditional"][key]["donors"]) or "No eligible combination") if "conditional" in r else ""
                         for key in ("reference","baseline")]
                writer.writerow([r["k"],"+".join(r["donors"]),r["rare"],r["het"],r["hbd"],r.get("r"),"+".join(r.get("old_donors",[])),*cases,r["reason"]])
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
