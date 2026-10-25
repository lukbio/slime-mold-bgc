"""Tabulate antiSMASH regions per genome and draw a heatmap of cluster types."""

import csv
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

with open(snakemake.params.genomes) as fh:  # noqa: F821
    meta = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}

regions = []
for path in snakemake.input:  # noqa: F821
    acc = Path(path).stem
    data = json.loads(Path(path).read_text())
    for rec in data["records"]:
        for i, area in enumerate(rec.get("areas", []), 1):
            regions.append({
                "accession": acc,
                "species": meta[acc]["name"],
                "group": meta[acc]["group"],
                "record": rec["id"],
                "region": i,
                "start": area["start"],
                "end": area["end"],
                "length": area["end"] - area["start"],
                "products": "+".join(area["products"]),
            })

with open(snakemake.output.regions, "w", newline="") as fh:  # noqa: F821
    w = csv.DictWriter(fh, fieldnames=list(regions[0]), delimiter="\t")
    w.writeheader()
    w.writerows(regions)

# Count each product type once per region (hybrids count toward each type).
counts: dict[str, Counter] = defaultdict(Counter)
for r in regions:
    for p in set(r["products"].split("+")):
        counts[r["accession"]][p] += 1
types = sorted({t for c in counts.values() for t in c}, key=lambda t: -sum(c[t] for c in counts.values()))
order = list(meta)  # config order: dictyostelids, other amoebae, fungi, bacteria

with open(snakemake.output.counts, "w", newline="") as fh:  # noqa: F821
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "species", "group", "total_regions", *types])
    for acc in order:
        total = sum(1 for r in regions if r["accession"] == acc)
        w.writerow([acc, meta[acc]["name"], meta[acc]["group"], total, *(counts[acc][t] for t in types)])

# Per-group ranges and medians, so README numbers come straight from here.
with open(snakemake.output.groups, "w", newline="") as fh:  # noqa: F821
    w = csv.writer(fh, delimiter="	")
    w.writerow(["group", "genomes", "measure", "min", "median", "max"])
    groups: dict[str, list[str]] = defaultdict(list)
    for acc in order:
        groups[meta[acc]["group"].split(" (")[0]].append(acc)
    for g, accs in groups.items():
        totals = [sum(1 for r in regions if r["accession"] == a) for a in accs]
        for name, vals in [("all regions", totals), *((t, [counts[a][t] for a in accs]) for t in types[:5])]:
            w.writerow([g, len(accs), name, min(vals), statistics.median(vals), max(vals)])

fig, ax = plt.subplots(figsize=(1.0 + 0.55 * len(types), 1.2 + 0.42 * len(order)))
grid = [[counts[a][t] for t in types] for a in order]
im = ax.imshow(grid, cmap="Blues", aspect="auto")
for i, row in enumerate(grid):
    for j, v in enumerate(row):
        if v:
            ax.text(j, i, v, ha="center", va="center", fontsize=8,
                    color="white" if v > max(map(max, grid)) / 2 else "black")
ax.set_xticks(range(len(types)), types, rotation=60, ha="right", fontsize=8)
ax.set_yticks(range(len(order)), [f"{meta[a]['name']}  ({meta[a]['group']})" for a in order], fontsize=8,
              style="italic")
ax.set_title("Biosynthetic gene clusters by type (antiSMASH)", fontsize=10)
fig.colorbar(im, ax=ax, shrink=0.6, label="regions")
fig.tight_layout()
fig.savefig(snakemake.output.heatmap, dpi=200)  # noqa: F821
