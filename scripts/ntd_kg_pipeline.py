#!/usr/bin/env python3
"""
BioSynth — NTD Drug Repurposing Knowledge Graph Pipeline
==========================================================

Builds a real, queryable knowledge graph of neglected tropical diseases (NTDs),
their pathogens/targets, and candidate drugs — pulling live from public
bio-data APIs instead of a hand-typed dataset:

  - Open Targets Platform (GraphQL) — disease -> target association scores,
    the same evidence Open Targets uses for its own target prioritisation.
    https://api.platform.opentargets.org/api/v4/graphql

  - ChEMBL (REST) — target -> known/candidate drugs, mechanism of action,
    clinical phase, and approval status.
    https://www.ebi.ac.uk/chembl/api/data

  - ChEMBL-NTD (via the main ChEMBL REST API) — phenotypic whole-organism
    screening hits specifically donated for NTD drug discovery (GSK TCAMS,
    Novartis-GNF, St Jude, DNDi lead-optimisation sets, the MMV Pathogen Box,
    etc). ChEMBL-NTD is "largely a subset of ChEMBL" per its own docs, so we
    reach it by filtering the standard /assay and /activity endpoints by
    pathogen organism rather than a separate API. Sets not yet folded into a
    ChEMBL release only exist as per-set downloads — see
    CHEMBL_NTD_DEPOSITED_SETS below and chembl.gitbook.io/chembl-ntd.
    https://www.ebi.ac.uk/chembl/api/data

  - TDR Targets (tdrtargets.org) — a chemogenomics resource curated
    specifically for NTD pathogens, with real complementary (not
    duplicate) coverage vs ChEMBL: its own paper reports only ~20% overlap.
    It does NOT publish a documented public REST/GraphQL API — target lists
    are built and exported through its web query builder. This script
    treats it as a manual-export source: run a target-prioritisation query
    on tdrtargets.org, export to CSV/TSV, and pass it via --tdr-csv so it
    gets merged into the graph. See load_tdr_targets_csv() below.

Output: a JSON graph (nodes/links) in the same schema the BioSynth demo
front-end (biosynth-ntd-knowledge-graph.html) consumes, plus a CSV of
scored drug-target-disease triples for review in a spreadsheet.

USAGE
-----
    pip install requests networkx
    python ntd_kg_pipeline.py --out ntd_graph.json --csv ntd_candidates.csv \
        --include-ntd-screens --tdr-csv path/to/tdrtargets_export.csv

This machine's sandbox has outbound network access disabled, so this script
has not been run against the live APIs here — run it locally / on a server
with internet access. Endpoints and query shapes were verified against
current API docs as of Aug 2026.

NOTE ON SCOPE
-------------
This pulls target-disease association evidence and known drug-target links.
It does NOT do de novo structure-based target discovery (docking, AlphaFold-based
pocket detection, etc.) — that is a separate, heavier pipeline stage. See the
"NEXT STAGES" notes at the bottom of this file.
"""

import argparse
import csv
import json
import sys
import time
from dataclasses import dataclass, field, asdict
from typing import Optional

import requests

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_URL = "https://www.ebi.ac.uk/chembl/api/data"

# The 21 WHO-listed NTDs (2024 roadmap list). Free-text names — resolved to
# Open Targets disease IDs (EFO/MONDO/Orphanet) at runtime via search, since
# hardcoding IDs silently goes stale when ontologies get remapped.
NTD_NAMES = [
    "Buruli ulcer", "Chagas disease", "Dengue", "Chikungunya",
    "Dracunculiasis", "Echinococcosis", "Foodborne trematodiases",
    "Human African trypanosomiasis", "Leishmaniasis", "Leprosy",
    "Lymphatic filariasis", "Mycetoma", "Onchocerciasis", "Rabies",
    "Scabies", "Schistosomiasis", "Snakebite envenoming",
    "Soil-transmitted helminthiasis", "Taeniasis", "Trachoma", "Yaws",
]

HEADERS = {"User-Agent": "BioSynth-NTD-KG/0.2 (research pipeline)"}

# NTD -> causative organism, needed because ChEMBL-style phenotypic screening
# data (and TDR Targets) are indexed by pathogen, not by WHO disease name.
# Diseases with no single causative organism (snakebite envenoming, scabies'
# ectoparasite is not typically the screening target) are omitted here.
ORGANISM_MAP = {
    "Buruli ulcer": "Mycobacterium ulcerans",
    "Chagas disease": "Trypanosoma cruzi",
    "Dengue": "Dengue virus",
    "Chikungunya": "Chikungunya virus",
    "Echinococcosis": "Echinococcus granulosus",
    "Human African trypanosomiasis": "Trypanosoma brucei",
    "Leishmaniasis": "Leishmania donovani",
    "Leprosy": "Mycobacterium leprae",
    "Lymphatic filariasis": "Wuchereria bancrofti",
    "Mycetoma": "Madurella mycetomatis",
    "Onchocerciasis": "Onchocerca volvulus",
    "Schistosomiasis": "Schistosoma mansoni",
    "Soil-transmitted helminthiasis": "Ascaris lumbricoides",
    "Taeniasis": "Taenia solium",
    "Trachoma": "Chlamydia trachomatis",
    "Yaws": "Treponema pallidum",
}

# Pointer list for the deposited sets on ChEMBL-NTD that are NOT (yet) folded
# into a main ChEMBL release, so they can't be reached via the REST API at
# all — only as direct per-set downloads. Kept short and disease-tagged;
# full list at https://chembl.gitbook.io/chembl-ntd/downloads
CHEMBL_NTD_DEPOSITED_SETS = [
    {"set": 27, "disease": "Chagas disease", "title": "T. cruzi LAPTc inhibitor screen (2024)"},
    {"set": 26, "disease": "Schistosomiasis", "title": "S. mansoni life-cycle-stage screen (2023)"},
    {"set": 25, "disease": "Schistosomiasis", "title": "Single-dose juvenile/adult schistosome actives (2021)"},
    {"set": 21, "disease": "multiple", "title": "MMV Pathogen Box — 400 actives across 12 NTDs"},
    {"set": 15, "disease": "multiple", "title": "DNDi: antiprotozoal activity profiling of approved drugs (repositioning)"},
    {"set": 14, "disease": "Leishmaniasis / Chagas / HAT", "title": "GSK TCAKS kinetoplastid screen"},
    {"set": 10, "disease": "Chagas disease", "title": "DNDi/Epichem T. cruzi HTS hit optimisation"},
    {"set": 9, "disease": "Tuberculosis (comparator)", "title": "GSK TCAMS M. tuberculosis screen"},
]


# --------------------------------------------------------------------------
# Data model — mirrors the node/edge shape the front-end graph expects
# --------------------------------------------------------------------------

@dataclass
class TargetHit:
    disease_name: str
    disease_id: str
    target_symbol: str
    target_id: str
    association_score: float  # Open Targets overall association score, 0-1


@dataclass
class DrugHit:
    drug_name: str
    molecule_chembl_id: str
    target_chembl_id: str
    target_name: str
    max_phase: Optional[float]      # 4 = approved somewhere
    first_approval: Optional[int]
    action_type: Optional[str]
    mechanism_desc: Optional[str]


@dataclass
class ScoredCandidate:
    drug: str
    disease: str
    target: str
    association_score: float   # 0-1, from Open Targets (target validation proxy)
    max_phase: float           # 0-4 (maturity/safety proxy — already approved = fast path)
    action_type: str
    priority_score: float = 0.0


# --------------------------------------------------------------------------
# Open Targets: resolve disease name -> ID, then pull associated targets
# --------------------------------------------------------------------------

def ot_query(query: str, variables: dict) -> dict:
    resp = requests.post(OT_URL, json={"query": query, "variables": variables},
                          headers=HEADERS, timeout=30)
    resp.raise_for_status()
    payload = resp.json()
    if "errors" in payload:
        raise RuntimeError(f"Open Targets GraphQL error: {payload['errors']}")
    return payload["data"]


def ot_search_disease(name: str) -> Optional[str]:
    q = """
    query search($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 1}) {
        hits { id name entity }
      }
    }
    """
    data = ot_query(q, {"q": name})
    hits = data["search"]["hits"]
    return hits[0]["id"] if hits else None


def ot_targets_for_disease(disease_id: str, disease_name: str, top_n: int = 8) -> list[TargetHit]:
    q = """
    query assoc($efoId: String!) {
      disease(efoId: $efoId) {
        associatedTargets(page: {index: 0, size: 25}) {
          rows {
            score
            target { id approvedSymbol }
          }
        }
      }
    }
    """
    data = ot_query(q, {"efoId": disease_id})
    rows = data["disease"]["associatedTargets"]["rows"]
    rows.sort(key=lambda r: r["score"], reverse=True)
    out = []
    for r in rows[:top_n]:
        out.append(TargetHit(
            disease_name=disease_name, disease_id=disease_id,
            target_symbol=r["target"]["approvedSymbol"],
            target_id=r["target"]["id"],
            association_score=r["score"],
        ))
    return out


# --------------------------------------------------------------------------
# ChEMBL: for each target, pull known/candidate drugs and their mechanism
# --------------------------------------------------------------------------

def chembl_target_id_for_symbol(symbol: str) -> Optional[str]:
    r = requests.get(f"{CHEMBL_URL}/target.json",
                      params={"target_synonym__icontains": symbol, "limit": 1},
                      headers=HEADERS, timeout=30)
    r.raise_for_status()
    targets = r.json().get("targets", [])
    return targets[0]["target_chembl_id"] if targets else None


def chembl_drugs_for_target(target_chembl_id: str, target_symbol: str) -> list[DrugHit]:
    r = requests.get(f"{CHEMBL_URL}/mechanism.json",
                      params={"target_chembl_id": target_chembl_id, "limit": 25},
                      headers=HEADERS, timeout=30)
    r.raise_for_status()
    mechanisms = r.json().get("mechanisms", [])
    out = []
    for m in mechanisms:
        mol_id = m.get("molecule_chembl_id")
        if not mol_id:
            continue
        mol = requests.get(f"{CHEMBL_URL}/molecule/{mol_id}.json",
                            headers=HEADERS, timeout=30).json()
        name = (mol.get("pref_name") or mol_id)
        out.append(DrugHit(
            drug_name=name,
            molecule_chembl_id=mol_id,
            target_chembl_id=target_chembl_id,
            target_name=target_symbol,
            max_phase=mol.get("max_phase"),
            first_approval=mol.get("first_approval"),
            action_type=m.get("action_type"),
            mechanism_desc=m.get("mechanism_of_action"),
        ))
        time.sleep(0.05)  # be polite to the free EBI endpoint
    return out


# --------------------------------------------------------------------------
# ChEMBL-NTD: phenotypic whole-organism screening hits, reached by filtering
# the main ChEMBL API by pathogen organism rather than a separate service.
# These are hits against the whole parasite/pathogen (no confirmed human
# target yet) — a different, earlier-stage evidence type than the
# mechanism-based drug hits above, valuable precisely because it surfaces
# compounds before anyone has mapped them to a target.
# --------------------------------------------------------------------------

def chembl_ntd_screening_hits(organism: str, max_ic50_nm: float = 10000,
                               limit: int = 20) -> list[dict]:
    """Potent (<= max_ic50_nm) phenotypic actives against a given pathogen
    organism, drawn from whole-cell/whole-organism assays in ChEMBL —
    the same donor sets (GSK, Novartis, St Jude, DNDi, MMV) that make up
    ChEMBL-NTD, insofar as they've been folded into a main ChEMBL release."""
    r = requests.get(f"{CHEMBL_URL}/activity.json", params={
        "target_organism": organism,
        "standard_type": "IC50",
        "standard_value__lte": max_ic50_nm,
        "standard_units": "nM",
        "limit": limit,
    }, headers=HEADERS, timeout=30)
    r.raise_for_status()
    out = []
    for a in r.json().get("activities", []):
        out.append({
            "molecule_chembl_id": a.get("molecule_chembl_id"),
            "compound_name": a.get("molecule_pref_name") or a.get("molecule_chembl_id"),
            "assay_description": a.get("assay_description"),
            "ic50_nm": a.get("standard_value"),
            "organism": organism,
            "document_chembl_id": a.get("document_chembl_id"),
        })
    return out


# --------------------------------------------------------------------------
# TDR Targets: no documented public API, so this ingests a manual CSV/TSV
# export from the tdrtargets.org query builder (Search for target genes ->
# prioritise -> export). Expected columns (rename to match your export):
#   gene_id, gene_symbol, organism, druggability_score, essentiality
# --------------------------------------------------------------------------

def load_tdr_targets_csv(path: str) -> list[dict]:
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "target_symbol": row.get("gene_symbol") or row.get("gene_id"),
                "organism": row.get("organism"),
                "tdr_druggability": _safe_float(row.get("druggability_score")),
                "tdr_essential": row.get("essentiality"),
            })
    return rows


def _safe_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------
# Scoring — same weighted-factor model as the demo front-end, but fed by
# live association scores + regulatory phase instead of hand-typed numbers.
# Swap in DisGeNET / patent-status / ATC-based affordability lookups here
# once you wire those sources in.
# --------------------------------------------------------------------------

def priority_score(assoc: float, max_phase: float,
                    w_target=0.4, w_maturity=0.6) -> float:
    """0-100 composite. assoc is 0-1 (Open Targets), max_phase is 0-4 (ChEMBL)."""
    return round((assoc * w_target * 100) + ((max_phase / 4.0) * w_maturity * 100), 1)


# --------------------------------------------------------------------------
# Pipeline
# --------------------------------------------------------------------------

def build_graph(disease_names: list[str], top_targets_per_disease: int = 5,
                 include_ntd_screens: bool = False, tdr_csv: Optional[str] = None):
    nodes: dict[str, dict] = {}
    links: list[dict] = []
    candidates: list[ScoredCandidate] = []

    def add_node(node_id, **kw):
        if node_id not in nodes:
            nodes[node_id] = {"id": node_id, **kw}
        else:
            nodes[node_id].update({k: v for k, v in kw.items() if v is not None})

    # Optional: merge a manually-exported TDR Targets prioritisation list.
    # Keyed by target_symbol so it enriches the same target nodes Open
    # Targets produces, rather than creating a parallel disconnected graph.
    tdr_by_target = {}
    if tdr_csv:
        print(f"[tdr] loading manual export: {tdr_csv}", file=sys.stderr)
        for row in load_tdr_targets_csv(tdr_csv):
            if row["target_symbol"]:
                tdr_by_target[row["target_symbol"]] = row

    for name in disease_names:
        print(f"[disease] resolving: {name}", file=sys.stderr)
        efo_id = ot_search_disease(name)
        if not efo_id:
            print(f"  no Open Targets match for {name}, skipping", file=sys.stderr)
            continue
        add_node(name, type="disease")

        targets = ot_targets_for_disease(efo_id, name, top_n=top_targets_per_disease)
        for t in targets:
            tdr_hit = tdr_by_target.get(t.target_symbol)
            add_node(t.target_symbol, type="target",
                     tdr_druggability=tdr_hit["tdr_druggability"] if tdr_hit else None,
                     tdr_essential=tdr_hit["tdr_essential"] if tdr_hit else None)
            links.append({"source": name, "target": t.target_symbol,
                          "kind": "associated", "score": t.association_score})

            chembl_id = chembl_target_id_for_symbol(t.target_symbol)
            if not chembl_id:
                continue
            drugs = chembl_drugs_for_target(chembl_id, t.target_symbol)
            for d in drugs:
                add_node(d.drug_name, type="drug",
                         max_phase=d.max_phase, action=d.action_type)
                links.append({"source": d.drug_name, "target": t.target_symbol,
                              "kind": "modulates", "action_type": d.action_type})
                links.append({"source": d.drug_name, "target": name, "kind": "candidate"})

                mp = d.max_phase or 0
                candidates.append(ScoredCandidate(
                    drug=d.drug_name, disease=name, target=t.target_symbol,
                    association_score=t.association_score, max_phase=mp,
                    action_type=d.action_type or "unknown",
                    priority_score=priority_score(t.association_score, mp),
                ))

        # ChEMBL-NTD-style phenotypic hits: compounds active against the
        # whole pathogen with no confirmed target yet. Added as a distinct
        # edge kind ("phenotypic_hit") so the front-end / analysts can tell
        # these apart from mechanism-confirmed candidates above — they're
        # earlier-stage, higher-risk, but sometimes the only lead that exists.
        if include_ntd_screens and name in ORGANISM_MAP:
            organism = ORGANISM_MAP[name]
            print(f"  [ntd-screen] {organism}", file=sys.stderr)
            try:
                hits = chembl_ntd_screening_hits(organism)
            except requests.RequestException as e:
                print(f"    screen query failed: {e}", file=sys.stderr)
                hits = []
            for h in hits:
                add_node(h["compound_name"], type="drug", phenotypic_only=True)
                links.append({"source": h["compound_name"], "target": name,
                              "kind": "phenotypic_hit", "ic50_nm": h["ic50_nm"],
                              "document": h["document_chembl_id"]})

    candidates.sort(key=lambda c: c.priority_score, reverse=True)
    return {"nodes": list(nodes.values()), "links": links}, candidates


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="ntd_graph.json", help="graph JSON output path")
    ap.add_argument("--csv", default="ntd_candidates.csv", help="scored candidates CSV path")
    ap.add_argument("--diseases", nargs="*", default=NTD_NAMES,
                     help="subset of NTD names to process (default: all 21)")
    ap.add_argument("--top-targets", type=int, default=5,
                     help="top-N Open Targets targets to pull per disease")
    ap.add_argument("--include-ntd-screens", action="store_true",
                     help="also pull ChEMBL-NTD-style phenotypic screening hits per pathogen organism")
    ap.add_argument("--tdr-csv", default=None,
                     help="path to a manually-exported TDR Targets CSV (see load_tdr_targets_csv)")
    ap.add_argument("--list-ntd-deposits", action="store_true",
                     help="print the ChEMBL-NTD deposited-set pointers (data not in the main API) and exit")
    args = ap.parse_args()

    if args.list_ntd_deposits:
        for s in CHEMBL_NTD_DEPOSITED_SETS:
            print(f"Set {s['set']:>2}  [{s['disease']}]  {s['title']}")
        print("\nFull list + downloads: https://chembl.gitbook.io/chembl-ntd/downloads")
        return

    graph, candidates = build_graph(
        args.diseases, top_targets_per_disease=args.top_targets,
        include_ntd_screens=args.include_ntd_screens, tdr_csv=args.tdr_csv,
    )

    with open(args.out, "w") as f:
        json.dump(graph, f, indent=2)
    print(f"wrote graph: {args.out} ({len(graph['nodes'])} nodes, {len(graph['links'])} links)")

    with open(args.csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(candidates[0]).keys()) if candidates else
                            ["drug", "disease", "target", "association_score",
                             "max_phase", "action_type", "priority_score"])
        w.writeheader()
        for c in candidates:
            w.writerow(asdict(c))
    print(f"wrote candidates: {args.csv} ({len(candidates)} rows)")


if __name__ == "__main__":
    main()

# --------------------------------------------------------------------------
# NEXT STAGES (not implemented here — roadmap notes)
# --------------------------------------------------------------------------
# 1. Affordability/access factor: pull ATC classification + generic
#    manufacturer count (e.g. via WHO ATC/DDD index, or a patent-status
#    lookup) to replace the placeholder in priority_score().
# 2. Cross-disease target overlap: after running this for both the NTD list
#    and a cancer/rare-disease target list, join on target_id to surface
#    drugs whose target is shared across BioSynth's three focus areas —
#    that overlap is the single highest-signal repurposing flag.
# 3. Literature evidence: layer in a PubMed count/recency query per
#    drug-disease pair (E-utilities API) as the "clinical evidence level"
#    factor instead of relying on ChEMBL max_phase alone.
# 4. Load into Neo4j: MERGE nodes/relationships from the JSON output using
#    the neo4j Python driver, or generate a LOAD CSV / apoc.load.json script,
#    to get real Cypher traversal ("targets shared by >=2 NTDs") instead of
#    scanning the JSON in Python.
# 5. Swap static disease list for the WHO NTD data portal feed so the
#    pipeline stays current as WHO adds/removes conditions from the list.
# 6. The deposited ChEMBL-NTD sets NOT yet in a ChEMBL release (see
#    CHEMBL_NTD_DEPOSITED_SETS / --list-ntd-deposits) still require manual
#    SDF/CSV download and a small loader — same pattern as
#    load_tdr_targets_csv() — since neither has a live query API.
# 7. TDR Targets: if their team confirms/publishes a stable bulk-download
#    format (check https://tdrtargets.org/releases), replace the manual
#    --tdr-csv step with a direct fetch.
