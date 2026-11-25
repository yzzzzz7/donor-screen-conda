"""Add the CDS homozygous-missense proportion to the original paired screen."""
import csv
import math
from itertools import groupby, permutations


def attach_pair_counts(path, recipients, donors, pairs):
    required = {tuple(sorted((a,b))) for a in recipients for b in recipients+donors if a!=b}
    with path.open(encoding="utf-8-sig",newline="") as f:
        for row in csv.DictReader(f,delimiter="\t"):
            key = tuple(sorted((row["sample1"],row["sample2"])))
            if key not in required:
                continue
            pair = pairs[key]
            if "r" in pair:
                raise ValueError(f"Repeated coding counts: {key}")
            m,c = float(row["expected_missense"]),float(row["expected_cds"])
            if not (math.isfinite(m) and math.isfinite(c) and 0 <= m <= c and c > 0):
                raise ValueError(f"Invalid predicted coding counts: {key}")
            pair.update(r=m/c,expected_missense=m,expected_cds=c)
    missing = [k for k in required if "r" not in pairs[k]]
    if missing:
        raise ValueError(f"Missing coding counts for {len(missing)} pairs, including {missing[0]}")


def annotate(plan, pairs):
    for p in plan["pairs"]:
        score = pairs[tuple(sorted((p["parent1"],p["parent2"])))]
        for key in ("r","expected_missense","expected_cds"):
            p[key] = score[key]
    plan["r"] = math.fsum(p["r"] for p in plan["pairs"])/len(plan["pairs"])
    return plan


def feasible_matings(option, recipients, pairs, cap):
    slots=tuple(p["parent1"] for p in option["pairs"])
    n=len(slots)
    if option["hbd"]>cap:
        return
    for rs in permutations(recipients,n):
        if any(slots[j]==slots[j-1] and rs[j]<rs[j-1] for j in range(1,n)):
            continue
        lineup=tuple(zip(slots,rs))
        cells=[pairs[tuple(sorted(e))] for e in lineup]
        if all(p["allowed"] for p in cells):
            hbd=math.fsum(p["hbd"] for p in cells)/n
            if hbd<=cap:
                yield hbd,math.fsum(p["r"] for p in cells)/n,lineup


def enrich_group_options(options, recipients, pairs, cap, r_cap, metrics):
    """For each offspring allocation seek the lowest HBD among *all* jointly eligible matings."""
    choices = []
    for old in options:
        best=min(((h,line) for h,r,line in feasible_matings(old,recipients,pairs,cap) if r<=r_cap),default=None)
        if best is not None:
            choice = dict(k=old["k"],donors=old["donors"],rare=old["rare"],**metrics(best[1]))
            choices.append(annotate(choice,pairs))
    return choices


def compare(row, old, r_cap, policy):
    row.update(old_donors=old["donors"],old_pairs=old["pairs"],old_r=old.get("r"),
               old_rare=old.get("rare"),r_cap=r_cap,r_policy=policy)
    if row["k"] == 0:
        row["change"]="无引入对照，沿用相同配对。"
    elif policy=="report":
        row["change"]="已计算R；未指定上限，四指标名次待定。"
    elif not old["donors"]:
        row["change"]="原亲缘和近交限制下已无可行供体。"
    elif not row["donors"]:
        row["change"]="原方案的R超过上限，本名额没有同时达标的组合。"
    elif old["donors"] != row["donors"]:
        row["change"]="原组合在双上限下无法保住原排序，新组合的预期新增变异最多。"
    elif [(p["parent1"],p["parent2"]) for p in old["pairs"]] != [(p["parent1"],p["parent2"]) for p in row["pairs"]]:
        row["change"]="供体相同，调整受体配对后同时满足R和近交上限。"
    else:
        row["change"]="原首选的R已达标，供体和配对均相同。"
    return row


def minimum_r_for_original_score(options, old, recipients, pairs, cap):
    """Smallest R attainable by the original donors without sacrificing their original score."""
    values=[]
    for o in options:
        if o["rare"]!=old["rare"] or o["hbd"]>cap:
            continue
        values.extend(r for _,r,_ in feasible_matings(o,recipients,pairs,cap))
    return min(values)


def decision_intervals(by_group, old, recipients, pairs, cap, refs, metrics):
    """Exact donor changes as a function of the user's cap, from one set of feasible matings."""
    points=[]
    for group,options in by_group.items():
        for o in options:
            for h,r,line in feasible_matings(o,recipients,pairs,cap):
                points.append((r,(-o["rare"],h,group,line)))
    points.sort()
    best=None
    intervals=[]
    optimal_at={}
    for boundary, same in groupby(points,key=lambda p:p[0]):
        rank=min(p[1] for p in same)
        if best is not None and best<=rank:
            continue
        if best is None or best[2]!=rank[2] or best[0]!=rank[0]:
            intervals.append(dict(min_cap=boundary,donors=list(rank[2]),rare=-rank[0]))
        best=rank
        for name,cutoff in refs.items():
            if boundary<=cutoff:
                optimal_at[name]=rank
    cases={}
    for name,cutoff in refs.items():
        rank=optimal_at.get(name)
        if rank is None:
            row=dict(k=old["k"],donors=[],pairs=[],r=None,hbd=None,het=None,rare=None)
        else:
            row=dict(k=old["k"],donors=list(rank[2]),rare=-rank[0],**metrics(rank[3]))
            annotate(row,pairs)
        cases[name]=compare(row,old,cutoff,name)
    return intervals,cases
