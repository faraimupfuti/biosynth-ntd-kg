#!/usr/bin/env python3
"""
BioSynth — Cross-Domain Target Overlap Detector
==================================================

The single highest-signal repurposing flag available to a company working
across NTDs, rare diseases, AND cancer: a target validated in one domain
that also shows up in another. A drug already approved for an oncology
target that happens to be druggable in a trypanosome, or vice versa, is a
much stronger lead than either target's association score alone suggests.

This does NOT re-run any API calls itself — it reads the JSON/CSV output
of two separate ntd_kg_pipeline.py runs (which is already domain-agnostic;
nothing in that script is NTD-specific except its default disease list) and
finds the intersection.

USAGE
-----
    # Run 1: the NTD graph (as always)
    python ntd_kg_pipeline.py --out data/ntd_graph.json --csv data/ntd_candidates.csv

    # Run 2: same script, a different disease list — e.g. cancers/rare diseases
    python ntd_kg_pipeline.py --diseases-csv "Breast cancer,Pancreatic cancer,Cystic fibrosis" \
        --out data/other_graph.json --csv data/other_candidates.csv

    # Then find the overlap
    python cross_domain_overlap.py \
        --graph-a data/ntd_graph.json --candidates-a data/ntd_candidates.csv --label-a NTD \
        --graph-b data/other_graph.json --candidates-b data/other_candidates.csv --label-b Cancer/Rare \
        --out data/cross_domain_overlap.md

NOTE ON THE DEFAULT DISEASE LIST BELOW
---------------------------------------
DEFAULT_OTHER_DISEASES is an illustrative starter list only — a mix of
common cancers and well-known rare diseases to make this runnable out of
the box. It is NOT BioSynth's actual chosen focus areas; that's a business
decision for your team, not something to infer from a generic example list.
Replace it with --diseases-csv on the second pipeline run once you know
which specific cancers/rare diseases you're prioritising.
"""

import argparse
import json
import sys
from collections import defaultdict

DEFAULT_OTHER_DISEASES = (
    "Breast cancer,Lung cancer,Pancreatic cancer,Colorectal cancer,"
    "Acute myeloid leukemia,Melanoma,Cystic fibrosis,"
    "Duchenne muscular dystrophy,Huntington's disease,Gaucher disease"
)


def load_graph(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def load_candidates(path: str) -> list[dict]:
    import csv
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def target_ids(graph: dict) -> set[str]:
    return {n["id"] for n in graph["nodes"] if n.get("type") == "target"}


def drugs_for_target(candidates: list[dict], target: str) -> list[dict]:
    return [c for c in candidates if c.get("target") == target]


def diseases_for_target(candidates: list[dict], target: str) -> set[str]:
    return {c["disease"] for c in candidates if c.get("target") == target}


def find_overlap(graph_a, cand_a, label_a, graph_b, cand_b, label_b) -> list[dict]:
    shared = target_ids(graph_a) & target_ids(graph_b)
    overlaps = []
    for target in sorted(shared):
        drugs_a = drugs_for_target(cand_a, target)
        drugs_b = drugs_for_target(cand_b, target)
        if not drugs_a or not drugs_b:
            continue  # target overlaps but no known drug on one side — nothing actionable yet
        overlaps.append({
            "target": target,
            f"diseases_{label_a}": sorted(diseases_for_target(cand_a, target)),
            f"diseases_{label_b}": sorted(diseases_for_target(cand_b, target)),
            f"drugs_{label_a}": sorted({d["drug"] for d in drugs_a}),
            f"drugs_{label_b}": sorted({d["drug"] for d in drugs_b}),
            # A drug already known on one side, appearing on the other side's
            # target list too, is the actual repurposing lead — surface it
            # explicitly rather than making the reader cross-reference lists.
            "cross_domain_drug_candidates": sorted(
                {d["drug"] for d in drugs_a} | {d["drug"] for d in drugs_b}
            ),
        })
    return overlaps


def generate_report(overlaps: list[dict], label_a: str, label_b: str) -> str:
    lines = [
        f"# Cross-Domain Target Overlap — {label_a} × {label_b}",
        "",
        "> **Research lead, not a recommendation.** A shared target across "
        "domains means the underlying biology is druggable in both contexts "
        "— it does NOT mean any specific drug is safe or effective for the "
        "new indication. Every entry below needs the same pharmacology "
        "review as any other candidate in the main report.",
        "",
    ]
    if not overlaps:
        lines.append(
            f"No targets with known drugs on both sides were found between "
            f"the {label_a} and {label_b} runs. This can mean either that "
            f"there's genuinely no overlap in this data, or that one/both "
            f"runs only covered a handful of diseases (increase "
            f"--top-targets or the disease list and re-run)."
        )
        return "\n".join(lines)

    lines.append(f"Found **{len(overlaps)} shared target(s)** with known drug candidates on both sides.")
    lines.append("")
    for o in overlaps:
        lines.append(f"## Target: {o['target']}")
        lines.append("")
        lines.append(f"- **{label_a} disease(s):** {', '.join(o[f'diseases_{label_a}'])}")
        lines.append(f"- **{label_a} drug(s) known against this target:** {', '.join(o[f'drugs_{label_a}']) or '—'}")
        lines.append(f"- **{label_b} disease(s):** {', '.join(o[f'diseases_{label_b}'])}")
        lines.append(f"- **{label_b} drug(s) known against this target:** {', '.join(o[f'drugs_{label_b}']) or '—'}")
        lines.append(f"- **Cross-domain candidates worth a literature check:** "
                     f"{', '.join(o['cross_domain_drug_candidates'])}")
        lines.append("")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graph-a", required=True)
    ap.add_argument("--candidates-a", required=True)
    ap.add_argument("--label-a", default="DomainA")
    ap.add_argument("--graph-b", required=True)
    ap.add_argument("--candidates-b", required=True)
    ap.add_argument("--label-b", default="DomainB")
    ap.add_argument("--out", default="cross_domain_overlap.md")
    args = ap.parse_args()

    try:
        graph_a, cand_a = load_graph(args.graph_a), load_candidates(args.candidates_a)
        graph_b, cand_b = load_graph(args.graph_b), load_candidates(args.candidates_b)
    except FileNotFoundError as e:
        print(f"error: {e} — run both pipeline passes first.", file=sys.stderr)
        sys.exit(1)

    overlaps = find_overlap(graph_a, cand_a, args.label_a, graph_b, cand_b, args.label_b)
    report = generate_report(overlaps, args.label_a, args.label_b)
    with open(args.out, "w") as f:
        f.write(report)
    print(f"wrote {args.out}: {len(overlaps)} shared target(s) with cross-domain drug candidates")


if __name__ == "__main__":
    main()
