"""Drop gene features whose exons overlap (malformed annotation antiSMASH rejects).

A few annotations (e.g. Acanthamoeba castellanii) contain one or two genes
whose CDS join has overlapping exons. antiSMASH stops on the first one, so
remove those genes (and their mRNA/CDS) and log what was dropped.
"""

import sys

from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")  # noqa: F821


def overlapping(feature) -> bool:
    parts = sorted((int(p.start), int(p.end)) for p in feature.location.parts)
    return any(b_start < a_end for (_, a_end), (b_start, _) in zip(parts, parts[1:]))


def tag(feature) -> str:
    q = feature.qualifiers
    return (q.get("locus_tag") or q.get("gene") or [""])[0]


records, dropped = [], 0
for rec in SeqIO.parse(snakemake.input[0], "genbank"):  # noqa: F821
    bad = {tag(f) for f in rec.features if f.type in ("CDS", "mRNA") and overlapping(f)}
    bad.discard("")
    keep = [f for f in rec.features if not (f.type in ("gene", "mRNA", "CDS") and tag(f) in bad)]
    dropped += len(bad)
    for t in sorted(bad):
        print(f"{rec.id}: dropped {t} (overlapping exons)", file=sys.stderr)
    rec.features = keep
    records.append(rec)
SeqIO.write(records, snakemake.output[0], "genbank")  # noqa: F821
print(f"dropped {dropped} genes", file=sys.stderr)
