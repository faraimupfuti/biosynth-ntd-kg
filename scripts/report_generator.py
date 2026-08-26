#!/usr/bin/env python3
"""
BioSynth — NTD Repurposing Candidate Research Brief Generator
================================================================

Turns ntd_candidates.csv into a readable Markdown report: review flags up
front (boxed warnings, zero-trial candidates), a ranked summary table, then
full detail per top candidate.

Deliberately called a "research brief," not a "clinical report" — the
latter has a specific regulatory meaning (an ICH Clinical Study Report tied
to an actual trial) that this output does not meet and should never be
represented as meeting. This is a triage document for a human reviewer.

USAGE
-----
    python report_generator.py --csv data/ntd_candidates.csv --out data/ntd_report.md
"""

import argparse
import csv
import sys

DISCLAIMER = (
    "> **Not a clinical or investment recommendation.** Every score below "
    "reflects computational association evidence (target-disease "
    "association strength plus the regulatory maturity of the drug's "
    "*original* indication) plus, where available, real-world safety and "
    "trial signals. None of this substitutes for expert pharmacology "
    "review, preclinical studies, or clinical trials. This report is "
    "generated automatically — verify every flagged item independently "
    "before acting on it."
)

TRUTHY = {"true", "1", "yes"}


def is_true(v) -> bool:
    return str(v).strip().lower() in TRUTHY


def fnum(v, default=0.0) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def load_candidates(path: str) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def generate_report(candidates: list[dict], top_n: int = 15) -> str:
    lines = [
        "# BioSynth NTD Repurposing — Candidate Research Brief",
        "",
        DISCLAIMER,
        "",
        f"Generated from {len(candidates)} scored drug-target-disease triples.",
        "",
    ]

    # Review flags first — this is the most important section, so it goes
    # at the top rather than being buried under the ranked table.
    boxed = [c for c in candidates if is_true(c.get("boxed_warning"))]
    no_trials = [c for c in candidates if c.get("disease_trial_count") == "0"]
    has_flags = boxed or no_trials
    if has_flags:
        lines.append("## Review flags")
        lines.append("")
        if boxed:
            names = ", ".join(f"{c['drug']} ({c['disease']})" for c in boxed)
            lines.append(f"- **{len(boxed)} candidate(s) carry an FDA boxed warning** "
                         f"— requires pharmacology sign-off before further "
                         f"consideration: {names}")
        if no_trials:
            names = ", ".join(f"{c['drug']} ({c['disease']})" for c in no_trials)
            lines.append(f"- **{len(no_trials)} candidate(s) have zero registered "
                         f"trials** for their NTD indication — early-stage, "
                         f"unvalidated hypotheses: {names}")
        lines.append("")

    ranked = sorted(candidates, key=lambda c: fnum(c.get("priority_score")), reverse=True)[:top_n]

    lines.append(f"## Top {len(ranked)} candidates by priority score")
    lines.append("")
    lines.append("| Rank | Drug | NTD | Target | Score | Max phase | Trials for this NTD | Boxed warning |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for i, c in enumerate(ranked, 1):
        lines.append(
            f"| {i} | {c.get('drug','')} | {c.get('disease','')} | {c.get('target','')} | "
            f"{fnum(c.get('priority_score')):.1f} | {c.get('max_phase','—')} | "
            f"{c.get('disease_trial_count','—')} | {c.get('boxed_warning','—')} |"
        )
    lines.append("")

    lines.append("## Candidate detail")
    lines.append("")
    for i, c in enumerate(ranked, 1):
        lines.append(f"### {i}. {c.get('drug','')} → {c.get('disease','')}")
        lines.append("")
        lines.append(f"- **Target:** {c.get('target','')}")
        lines.append(f"- **Mechanism / action type:** {c.get('action_type','unknown')}")
        lines.append(f"- **Open Targets association score:** {c.get('association_score','—')}")
        lines.append(f"- **ChEMBL max clinical phase (original indication):** {c.get('max_phase','—')}")
        lines.append(f"- **Priority score:** {fnum(c.get('priority_score')):.1f}/100")
        if c.get("disease_trial_count", "") not in ("", None):
            lines.append(
                f"- **Registered trials for this NTD:** {c.get('disease_trial_count')} "
                f"(max phase reached: {c.get('disease_max_trial_phase','—')}, "
                f"results posted: {c.get('disease_trial_has_results','—')})"
            )
        if c.get("boxed_warning", "") not in ("", None):
            lines.append(f"- **FDA boxed warning:** {c.get('boxed_warning')}")
        if c.get("contraindications"):
            lines.append(f"- **Contraindications (label excerpt):** {c.get('contraindications')}")
        if c.get("serious_ae_reports", "") not in ("", None):
            lines.append(f"- **Serious adverse event reports (FAERS):** {c.get('serious_ae_reports')}")
        lines.append("")

    lines.append("---")
    lines.append(
        "*Recommended next step for any candidate above: literature review by "
        "a domain pharmacologist, followed by in vitro confirmation if not "
        "already published, before any further investment of time or capital.*"
    )
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default="data/ntd_candidates.csv")
    ap.add_argument("--out", default="data/ntd_report.md")
    ap.add_argument("--top", type=int, default=15, help="how many top candidates to detail")
    args = ap.parse_args()

    try:
        candidates = load_candidates(args.csv)
    except FileNotFoundError:
        print(f"error: {args.csv} not found — run the pipeline first.", file=sys.stderr)
        sys.exit(1)

    if not candidates:
        print("No candidates found in CSV — nothing to report.", file=sys.stderr)
        sys.exit(1)

    report = generate_report(candidates, top_n=args.top)
    with open(args.out, "w") as f:
        f.write(report)
    print(f"wrote report: {args.out} ({len(candidates)} candidates, top {min(args.top, len(candidates))} detailed)")


if __name__ == "__main__":
    main()
