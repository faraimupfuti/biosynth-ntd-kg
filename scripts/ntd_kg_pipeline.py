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

  - openFDA Drug Labels + Adverse Events (--include-clinical-evidence) —
    boxed warnings, contraindications, and dosing pulled from the actual
    FDA-approved label (SPL), plus a serious-adverse-event report count
    from FAERS. https://api.fda.gov/drug/label.json and .../drug/event.json
    openFDA's own docs are explicit that this data is NOT validated for
    clinical use and should not inform medical care decisions — treat it
    the same way here: a flag for a human reviewer to look closer, never a
    number folded silently into the priority score.

  - ClinicalTrials.gov API v2 (--include-clinical-evidence) — whether this
    EXACT drug-disease pair has ever been registered as a trial, at what
    phase, and whether results were posted. This is deliberately kept
    separate from ChEMBL's max_phase, which only reflects the drug's
    ORIGINAL indication — a drug at Phase 4 for cancer can easily have zero
    trials for an NTD, and that gap is exactly the signal a reviewer needs.
    https://clinicaltrials.gov/api/v2/studies

  NOT INTEGRATED: DrugBank. Free access only covers browsing/academic
  downloads; API access and any commercial use require a paid license
  negotiated directly with DrugBank — there is no self-serve commercial
  tier. Worth asking their team about the "internal research use" carve-out
  before assuming a full license is needed.

  - PubMed E-utilities (--include-literature) — publication count/recency
    per drug-disease pair, as a real evidence-level signal that catches
    candidates with genuine preclinical/literature support before they've
    reached a registered trial. https://eutils.ncbi.nlm.nih.gov/entrez/eutils/

  - WHO Global Health Observatory OData API (--include-burden) — disease
    burden (reported case counts / people requiring treatment) per NTD,
    resolved dynamically by searching indicator names rather than
    hardcoding indicator codes that could silently go stale.
    https://ghoapi.azureedge.net/api/ — free, no auth.

  - openFDA storage/handling (bundled into --include-clinical-evidence,
    no extra API call) — the label lookup already fetches this field; we
    just weren't reading it. Cold-chain requirements matter specifically
    for NTDs, where deployment is disproportionately rural/tropical.

  NOT INTEGRATED AS LIVE APIs — patent status, WHO Essential Medicines List
  membership, pathogen drug-resistance data, and manufacturing cost/
  feasibility. None of these have a clean, free, queryable API:
  WIPO PatentScope's programmatic access is a paid SOAP service (2,000 CHF/
  year); resistance surveillance lives in WHO PDF reports and scattered
  literature; manufacturing cost data lives in Global Fund/MSF Access
  Campaign reports read by people, not machines. Building a fake
  integration against any of these would produce false confidence exactly
  where the most caution is needed. Instead, see --expert-annotations-csv
  and load_expert_annotations() below: a human curates these fields, the
  pipeline merges them in, same honest pattern as --tdr-csv.

  ALSO NOT INTEGRATED, DELIBERATELY: ADMET data beyond what's on the
  approved label. There's no reasonable way to curate this speculatively
  per repurposing candidate — it's a genuine wet-lab/computational-chemistry
  question, not a data-integration one. Flagging that plainly here rather
  than building an empty placeholder that would look like coverage.

Output: a JSON graph (nodes/links) in the same schema the BioSynth demo
front-end (biosynth-ntd-knowledge-graph.html) consumes, plus a CSV of
scored drug-target-disease triples for review in a spreadsheet.

USAGE
-----
    pip install requests networkx
    python ntd_kg_pipeline.py --out ntd_graph.json --csv ntd_candidates.csv \
        --include-ntd-screens --include-clinical-evidence \
        --tdr-csv path/to/tdrtargets_export.csv

This machine's sandbox has outbound network access disabled, so this script
has not been run against the live APIs here — run it locally / on a server
with internet access. Endpoints and query shapes were verified against
current API docs as of Aug 2026.

NOTE ON SCOPE
-------------
This pulls target-disease association evidence and known drug-target links,
plus (optionally) real-world safety/trial signals. It does NOT do de novo
structure-based target discovery (docking, AlphaFold-based pocket
detection, etc.), and it is NOT a clinical decision tool — see the
"NEXT STAGES" notes at the bottom of this file for what real clinical-grade
validation would still require (wet-lab work, regulatory review, expert
pharmacology review of every flagged candidate).
"""

import argparse
import csv
import json
import random
import sys
import time
from dataclasses import dataclass, field, asdict
from typing import Optional

import requests

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_URL = "https://www.ebi.ac.uk/chembl/api/data"
OPENFDA_URL = "https://api.fda.gov"
CTGOV_URL = "https://clinicaltrials.gov/api/v2/studies"
PUBMED_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GHO_URL = "https://ghoapi.azureedge.net/api"

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

HEADERS = {"User-Agent": "BioSynth-NTD-KG/0.3 (research pipeline)"}


def request_with_retry(method: str, url: str, *, max_retries: int = 4,
                        base_delay: float = 2.0, **kwargs) -> requests.Response:
    """Shared retry/backoff wrapper for every outbound call in this script.

    The public EBI (ChEMBL) and Open Targets endpoints are free services —
    they occasionally return 500/502/503/504 under load, or 429 if you're
    hitting them fast. Without this, a single transient error kills the
    entire run (that's exactly what happened on CHEMBL227 in the first live
    test — the target has a very large mechanism list and the free server
    timed out on it once). Retries with exponential backoff + jitter turn a
    one-off blip into a few extra seconds instead of a failed run.
    """
    last_exc = None
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.request(method, url, headers=HEADERS, timeout=30, **kwargs)
            if resp.status_code in (429, 500, 502, 503, 504):
                raise requests.exceptions.HTTPError(
                    f"{resp.status_code} on attempt {attempt}", response=resp)
            resp.raise_for_status()
            return resp
        except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError,
                requests.exceptions.Timeout) as e:
            last_exc = e
            if attempt == max_retries:
                break
            delay = base_delay * (2 ** (attempt - 1)) + random.uniform(0, 1)
            print(f"    [retry] {url.split('?')[0]} failed ({e}); "
                  f"retrying in {delay:.1f}s ({attempt}/{max_retries})", file=sys.stderr)
            time.sleep(delay)
    raise last_exc

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
    # --- clinical-evidence fields (only populated with --include-clinical-evidence) ---
    # Deliberately kept OUT of priority_score — a boxed warning or a missing
    # trial is a judgment call for a human reviewer, not something to average
    # into one number where it could get diluted or hidden.
    boxed_warning: Optional[bool] = None
    contraindications: Optional[str] = None
    serious_ae_reports: Optional[int] = None
    disease_trial_count: Optional[int] = None
    disease_max_trial_phase: Optional[int] = None
    disease_trial_has_results: Optional[bool] = None
    # --- literature / burden (--include-literature, --include-burden) ---
    pubmed_result_count: Optional[int] = None
    disease_burden_estimate: Optional[str] = None  # kept as text: WHO reports mix units/years
    # --- storage (bundled into openFDA label lookup above) ---
    storage_requires_cold_chain: Optional[bool] = None
    # --- expert-curated overlay (--expert-annotations-csv) ---
    patent_status: Optional[str] = None
    eml_listed: Optional[bool] = None
    resistance_notes: Optional[str] = None
    manufacturing_notes: Optional[str] = None


# --------------------------------------------------------------------------
# Open Targets: resolve disease name -> ID, then pull associated targets
# --------------------------------------------------------------------------

def ot_query(query: str, variables: dict) -> dict:
    resp = request_with_retry("POST", OT_URL, json={"query": query, "variables": variables})
    payload = resp.json()
    if "errors" in payload:
        raise RuntimeError(f"Open Targets GraphQL error: {payload['errors']}")
    return payload["data"]


def ot_search_disease(name: str) -> Optional[tuple[str, str]]:
    """Returns (efo_id, matched_name) so callers can see and log exactly what
    Open Targets resolved the query to — silently trusting a fuzzy free-text
    match is how a mistyped or malformed disease name turns into a
    completely unrelated graph (see the Chagas/quoting incident)."""
    q = """
    query search($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 1}) {
        hits { id name entity }
      }
    }
    """
    data = ot_query(q, {"q": name})
    hits = data["search"]["hits"]
    if not hits:
        return None
    return hits[0]["id"], hits[0]["name"]


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
    r = request_with_retry("GET", f"{CHEMBL_URL}/target.json",
                            params={"target_synonym__icontains": symbol, "limit": 1})
    targets = r.json().get("targets", [])
    return targets[0]["target_chembl_id"] if targets else None


def chembl_drugs_for_target(target_chembl_id: str, target_symbol: str) -> list[DrugHit]:
    r = request_with_retry("GET", f"{CHEMBL_URL}/mechanism.json",
                            params={"target_chembl_id": target_chembl_id, "limit": 25})
    mechanisms = r.json().get("mechanisms", [])
    out = []
    for m in mechanisms:
        mol_id = m.get("molecule_chembl_id")
        if not mol_id:
            continue
        mol = request_with_retry("GET", f"{CHEMBL_URL}/molecule/{mol_id}.json").json()
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
    r = request_with_retry("GET", f"{CHEMBL_URL}/activity.json", params={
        "target_organism": organism,
        "standard_type": "IC50",
        "standard_value__lte": max_ic50_nm,
        "standard_units": "nM",
        "limit": limit,
    })
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
# openFDA: drug label flags (dosing/contraindications/boxed warnings) and a
# serious-adverse-event report count from FAERS. openFDA's own docs state
# this data is NOT validated for clinical use — we surface it as flags for
# a human reviewer, never as an input to priority_score.
# --------------------------------------------------------------------------

def openfda_label_flags(drug_name: str) -> dict:
    """Best-effort lookup by generic name. Matching is exact-ish (openFDA's
    Lucene search over openfda.generic_name), so brand names, salts, or
    naming mismatches with ChEMBL's pref_name will legitimately miss —
    that's a real coverage gap to be aware of, not a bug to silently paper
    over with fuzzy matching that could attach the wrong drug's label."""
    try:
        r = request_with_retry("GET", f"{OPENFDA_URL}/drug/label.json", params={
            "search": f'openfda.generic_name:"{drug_name}"', "limit": 1,
        })
    except requests.RequestException as e:
        return {"label_found": False, "lookup_error": str(e)}
    results = r.json().get("results", [])
    if not results:
        return {"label_found": False}
    label = results[0]

    def first(field):
        vals = label.get(field)
        return vals[0][:400] if vals else None

    return {
        "label_found": True,
        "boxed_warning": bool(label.get("boxed_warning")),
        "contraindications": first("contraindications"),
        "drug_interactions_noted": bool(label.get("drug_interactions")),
        "dosage_snippet": first("dosage_and_administration"),
        # Cold-chain matters disproportionately for NTDs given where they're
        # deployed — crude heuristic (label mentions refrigeration/freezing),
        # not a substitute for actually reading the storage section.
        "requires_cold_chain": bool(
            label.get("storage_and_handling") and any(
                kw in (label["storage_and_handling"][0] or "").lower()
                for kw in ("refrigerat", "2°c", "2-8", "freez", "do not store above 25")
            )
        ),
    }


def openfda_serious_ae_count(drug_name: str) -> Optional[int]:
    try:
        r = request_with_retry("GET", f"{OPENFDA_URL}/drug/event.json", params={
            "search": f'patient.drug.medicinalproduct:"{drug_name}" AND serious:1',
            "limit": 1,
        })
    except requests.RequestException:
        return None
    return r.json().get("meta", {}).get("results", {}).get("total")


# --------------------------------------------------------------------------
# ClinicalTrials.gov v2: has THIS SPECIFIC drug-disease pair ever been
# registered as a trial? Deliberately separate from ChEMBL's max_phase,
# which only reflects the drug's ORIGINAL indication — a drug at Phase 4
# for cancer can have zero trials for an NTD, and that gap matters.
# --------------------------------------------------------------------------

_CT_PHASE_RANK = {"NA": 0, "EARLY_PHASE1": 1, "PHASE1": 1, "PHASE2": 2,
                  "PHASE3": 3, "PHASE4": 4}


def ct_trials_for_pair(drug_name: str, disease_name: str) -> dict:
    try:
        r = request_with_retry("GET", CTGOV_URL, params={
            "query.term": f"{drug_name} AND {disease_name}",
            "pageSize": 10,
            "fields": "protocolSection.designModule.phases,hasResults",
        })
    except requests.RequestException as e:
        return {"trial_count": None, "lookup_error": str(e)}
    studies = r.json().get("studies", [])
    phases, has_results = [], False
    for s in studies:
        phases.extend(s.get("protocolSection", {}).get("designModule", {}).get("phases", []))
        if s.get("hasResults"):
            has_results = True
    max_phase = max((_CT_PHASE_RANK.get(p, 0) for p in phases), default=0)
    return {
        "trial_count": len(studies),
        "max_trial_phase_for_disease": max_phase,
        "any_results_posted": has_results,
    }


# --------------------------------------------------------------------------
# PubMed E-utilities: publication count for a drug-disease pair, as a real
# evidence-level signal. Catches candidates with genuine literature support
# that hasn't reached a registered trial yet — a gap ClinicalTrials.gov
# alone can't see. No API key required at this call volume; NCBI asks for
# one (free, via an NCBI account) if you scale this up past ~3 req/sec.
# --------------------------------------------------------------------------

def pubmed_result_count(drug_name: str, disease_name: str) -> Optional[int]:
    try:
        r = request_with_retry("GET", f"{PUBMED_URL}/esearch.fcgi", params={
            "db": "pubmed",
            "term": f'("{drug_name}"[Title/Abstract]) AND ("{disease_name}"[Title/Abstract])',
            "retmode": "json",
            "retmax": 0,  # we only want the count, not the actual PMIDs
        })
    except requests.RequestException:
        return None
    try:
        return int(r.json()["esearchresult"]["count"])
    except (KeyError, ValueError):
        return None


# --------------------------------------------------------------------------
# WHO GHO OData API: disease burden, resolved by searching indicator NAMES
# rather than hardcoding indicator codes — codes can be renamed/retired,
# and a hardcoded table would silently go stale with no error to notice.
# Free, no auth: https://ghoapi.azureedge.net/api/
# --------------------------------------------------------------------------

def gho_burden_for_disease(disease_name: str) -> Optional[dict]:
    """Best-effort: finds a GHO indicator whose name mentions this disease
    and reports a case/burden count, then returns its most recent value.
    NTDs are unevenly covered in GHO (case-count indicators exist for most
    of the preventive-chemotherapy diseases, not all 21) — a None/empty
    result here is a real coverage gap, not necessarily a bug."""
    try:
        r = request_with_retry("GET", f"{GHO_URL}/Indicator", params={
            "$filter": f"contains(IndicatorName,'{disease_name}')",
        })
    except requests.RequestException as e:
        return {"lookup_error": str(e)}
    indicators = r.json().get("value", [])
    # Prefer an indicator that looks like a case-count/burden metric over
    # e.g. "status of endemicity" (categorical, not a number worth trending).
    burden_kw = ("case", "number", "requiring", "reported")
    candidates = [i for i in indicators
                  if any(k in i.get("IndicatorName", "").lower() for k in burden_kw)]
    target = (candidates or indicators)
    if not target:
        return None
    code = target[0]["IndicatorCode"]
    try:
        r2 = request_with_retry("GET", f"{GHO_URL}/{code}", params={"$top": 1000})
    except requests.RequestException as e:
        return {"lookup_error": str(e)}
    rows = r2.json().get("value", [])
    if not rows:
        return None
    latest = max(rows, key=lambda row: row.get("TimeDim", 0))
    return {
        "indicator_name": target[0]["IndicatorName"],
        "year": latest.get("TimeDim"),
        "value": latest.get("NumericValue") or latest.get("Value"),
    }


# --------------------------------------------------------------------------
# Expert-curated overlay: patent status, WHO EML membership, resistance
# notes, manufacturing notes. None of these have a clean free API — see the
# module docstring for why. This is a human filling in a spreadsheet, not
# automation pretending otherwise. Merged onto candidates by drug name.
# --------------------------------------------------------------------------

def load_expert_annotations(path: str) -> dict:
    """Expected CSV columns: drug, patent_status, eml_listed, resistance_notes,
    manufacturing_notes. Any column can be blank if not yet researched —
    that's an honest 'not yet reviewed' state, not an error."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            drug = (row.get("drug") or "").strip()
            if not drug:
                continue
            out[drug] = {
                "patent_status": row.get("patent_status") or None,
                "eml_listed": is_true(row.get("eml_listed")) if row.get("eml_listed") else None,
                "resistance_notes": row.get("resistance_notes") or None,
                "manufacturing_notes": row.get("manufacturing_notes") or None,
            }
    return out


def is_true(v) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes")


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
                 include_ntd_screens: bool = False, tdr_csv: Optional[str] = None,
                 include_clinical_evidence: bool = False, include_literature: bool = False,
                 include_burden: bool = False, expert_csv: Optional[str] = None):
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

    # Optional: expert-curated overlay (patent/EML/resistance/manufacturing).
    # Keyed by drug name — see load_expert_annotations() for why this is a
    # human-filled CSV rather than a fake API integration.
    expert_by_drug = {}
    if expert_csv:
        print(f"[expert] loading manual annotations: {expert_csv}", file=sys.stderr)
        expert_by_drug = load_expert_annotations(expert_csv)

    burden_by_disease = {}

    for name in disease_names:
        print(f"[disease] resolving: {name}", file=sys.stderr)
        if include_burden:
            burden = gho_burden_for_disease(name)
            if burden and "lookup_error" not in burden:
                burden_by_disease[name] = burden
                print(f"  [burden] {burden.get('indicator_name')}: "
                      f"{burden.get('value')} ({burden.get('year')})", file=sys.stderr)
        resolved = ot_search_disease(name)
        if not resolved:
            print(f"  no Open Targets match for {name}, skipping", file=sys.stderr)
            continue
        efo_id, matched_name = resolved
        print(f"  matched -> {matched_name} ({efo_id})", file=sys.stderr)
        if matched_name.strip().lower() != name.strip().lower():
            print(f"  [warn] matched name differs from input — verify this is "
                  f"the disease you meant before trusting this branch of the graph",
                  file=sys.stderr)
        add_node(name, type="disease", matched_name=matched_name)

        try:
            targets = ot_targets_for_disease(efo_id, name, top_n=top_targets_per_disease)
        except requests.RequestException as e:
            print(f"  [warn] Open Targets lookup failed for {name}: {e}", file=sys.stderr)
            continue

        for t in targets:
            tdr_hit = tdr_by_target.get(t.target_symbol)
            add_node(t.target_symbol, type="target",
                     tdr_druggability=tdr_hit["tdr_druggability"] if tdr_hit else None,
                     tdr_essential=tdr_hit["tdr_essential"] if tdr_hit else None)
            links.append({"source": name, "target": t.target_symbol,
                          "kind": "associated", "score": t.association_score})

            try:
                chembl_id = chembl_target_id_for_symbol(t.target_symbol)
            except requests.RequestException as e:
                print(f"  [warn] ChEMBL target lookup failed for "
                      f"{t.target_symbol}: {e}", file=sys.stderr)
                continue
            if not chembl_id:
                continue
            try:
                drugs = chembl_drugs_for_target(chembl_id, t.target_symbol)
            except requests.RequestException as e:
                # A single problematic target (e.g. one with an unusually
                # large mechanism list) should not take down the whole run —
                # log it and keep going with everything else.
                print(f"  [warn] ChEMBL lookup failed for target "
                      f"{t.target_symbol} ({chembl_id}): {e}", file=sys.stderr)
                continue
            for d in drugs:
                add_node(d.drug_name, type="drug",
                         max_phase=d.max_phase, action=d.action_type)
                links.append({"source": d.drug_name, "target": t.target_symbol,
                              "kind": "modulates", "action_type": d.action_type})
                links.append({"source": d.drug_name, "target": name, "kind": "candidate"})

                mp = _safe_float(d.max_phase) or 0.0
                candidate = ScoredCandidate(
                    drug=d.drug_name, disease=name, target=t.target_symbol,
                    association_score=t.association_score, max_phase=mp,
                    action_type=d.action_type or "unknown",
                    priority_score=priority_score(t.association_score, mp),
                )

                # Optional real-world safety/trial-evidence lookups. Off by
                # default because this roughly triples the API calls per
                # drug (label + FAERS + trials on top of the ChEMBL calls
                # already made) — enable with --include-clinical-evidence.
                if include_clinical_evidence:
                    label = openfda_label_flags(d.drug_name)
                    candidate.boxed_warning = label.get("boxed_warning")
                    candidate.contraindications = label.get("contraindications")
                    candidate.serious_ae_reports = openfda_serious_ae_count(d.drug_name)

                    ct = ct_trials_for_pair(d.drug_name, name)
                    candidate.disease_trial_count = ct.get("trial_count")
                    candidate.disease_max_trial_phase = ct.get("max_trial_phase_for_disease")
                    candidate.disease_trial_has_results = ct.get("any_results_posted")

                    if candidate.boxed_warning:
                        print(f"    [flag] {d.drug_name} carries an FDA boxed "
                              f"warning — needs pharmacology review before any "
                              f"further consideration", file=sys.stderr)
                    if candidate.disease_trial_count == 0:
                        print(f"    [note] {d.drug_name} has no registered "
                              f"trials for {name} — this is an early-stage, "
                              f"unvalidated repurposing hypothesis", file=sys.stderr)

                    add_node(d.drug_name, type="drug", boxed_warning=candidate.boxed_warning,
                             disease_trial_count=candidate.disease_trial_count)

                if include_literature:
                    candidate.pubmed_result_count = pubmed_result_count(d.drug_name, name)

                if include_burden and name in burden_by_disease:
                    b = burden_by_disease[name]
                    candidate.disease_burden_estimate = f"{b.get('value')} ({b.get('year')}, {b.get('indicator_name')})"

                if include_clinical_evidence:
                    candidate.storage_requires_cold_chain = label.get("requires_cold_chain")

                if d.drug_name in expert_by_drug:
                    exp = expert_by_drug[d.drug_name]
                    candidate.patent_status = exp["patent_status"]
                    candidate.eml_listed = exp["eml_listed"]
                    candidate.resistance_notes = exp["resistance_notes"]
                    candidate.manufacturing_notes = exp["manufacturing_notes"]

                candidates.append(candidate)

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
    ap.add_argument("--diseases", nargs="*", default=None,
                     help="subset of NTD names to process, space-separated "
                          "SINGLE-WORD-SAFE names only (e.g. Dengue Yaws). "
                          "For multi-word names use --diseases-csv instead — "
                          "shell/CI word-splitting will silently break multi-word "
                          "names passed here (see --diseases-csv docstring below).")
    ap.add_argument("--diseases-csv", default=None,
                     help="comma-separated disease names, safe for multi-word names "
                          "and for passing through CI inputs without quoting pitfalls, "
                          "e.g. --diseases-csv \"Chagas disease,Dengue,Yaws\". "
                          "Takes precedence over --diseases if both are given.")
    ap.add_argument("--top-targets", type=int, default=5,
                     help="top-N Open Targets targets to pull per disease")
    ap.add_argument("--include-ntd-screens", action="store_true",
                     help="also pull ChEMBL-NTD-style phenotypic screening hits per pathogen organism")
    ap.add_argument("--include-clinical-evidence", action="store_true",
                     help="also pull openFDA label flags (boxed warnings, contraindications, "
                          "cold-chain storage), FAERS serious-adverse-event counts, and "
                          "ClinicalTrials.gov trial history for this specific drug-disease "
                          "pair. Roughly triples API calls per drug — off by default.")
    ap.add_argument("--include-literature", action="store_true",
                     help="also pull a PubMed publication count per drug-disease pair "
                          "(one extra API call per drug)")
    ap.add_argument("--include-burden", action="store_true",
                     help="also pull WHO GHO disease burden (case counts) per disease "
                          "(one extra API call per disease, not per drug)")
    ap.add_argument("--expert-annotations-csv", default=None,
                     help="path to a human-curated CSV with columns: drug, patent_status, "
                          "eml_listed, resistance_notes, manufacturing_notes. These fields "
                          "have no clean free API — see the module docstring for why — so "
                          "this is where a person's research gets merged in instead.")
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

    if args.diseases_csv:
        disease_list = [d.strip() for d in args.diseases_csv.split(",") if d.strip()]
    elif args.diseases:
        disease_list = args.diseases
    else:
        disease_list = NTD_NAMES

    graph, candidates = build_graph(
        disease_list, top_targets_per_disease=args.top_targets,
        include_ntd_screens=args.include_ntd_screens, tdr_csv=args.tdr_csv,
        include_clinical_evidence=args.include_clinical_evidence,
        include_literature=args.include_literature,
        include_burden=args.include_burden,
        expert_csv=args.expert_annotations_csv,
    )

    with open(args.out, "w") as f:
        json.dump(graph, f, indent=2)
    print(f"wrote graph: {args.out} ({len(graph['nodes'])} nodes, {len(graph['links'])} links)")

    fallback_fields = ["drug", "disease", "target", "association_score", "max_phase",
                        "action_type", "priority_score", "boxed_warning", "contraindications",
                        "serious_ae_reports", "disease_trial_count",
                        "disease_max_trial_phase", "disease_trial_has_results",
                        "pubmed_result_count", "disease_burden_estimate",
                        "storage_requires_cold_chain", "patent_status", "eml_listed",
                        "resistance_notes", "manufacturing_notes"]
    with open(args.csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(candidates[0]).keys()) if candidates
                            else fallback_fields)
        w.writeheader()
        for c in candidates:
            w.writerow(asdict(c))
    print(f"wrote candidates: {args.csv} ({len(candidates)} rows)")
    if args.include_clinical_evidence and candidates:
        flagged = sum(1 for c in candidates if c.boxed_warning)
        unvalidated = sum(1 for c in candidates if c.disease_trial_count == 0)
        print(f"  {flagged} candidate(s) carry an FDA boxed warning — review before proceeding")
        print(f"  {unvalidated} candidate(s) have zero registered trials for their NTD indication")


if __name__ == "__main__":
    main()

# --------------------------------------------------------------------------
# NEXT STAGES (not implemented here — roadmap notes)
# --------------------------------------------------------------------------
# 1. Affordability/access factor: --expert-annotations-csv now carries
#    patent_status and eml_listed, but neither feeds priority_score yet —
#    same deliberate choice as boxed_warning: these are judgment calls for
#    a human reviewer, not something to silently average into one number.
# 2. Cross-disease target overlap: after running this for both the NTD list
#    and a cancer/rare-disease target list, join on target_id to surface
#    drugs whose target is shared across BioSynth's three focus areas —
#    that overlap is the single highest-signal repurposing flag.
# 3. DONE (--include-clinical-evidence): openFDA labels + FAERS + trial
#    history for the specific drug-disease pair, including cold-chain
#    storage flags. DONE (--include-literature): PubMed publication counts.
#    DONE (--include-burden): WHO GHO disease burden. DONE
#    (--expert-annotations-csv): patent status, EML listing, resistance and
#    manufacturing notes as a human-curated overlay — see the module
#    docstring for why these couldn't be live API integrations.
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
# 7. Genuinely still out of scope, not just deferred: ADMET beyond the
#    label, and any form of automated resistance/manufacturing-cost
#    inference. These need a domain expert's actual judgment, not a script
#    guessing at it — see --expert-annotations-csv for the honest version
#    of "integrating" this: a person fills it in.
# 7. TDR Targets: if their team confirms/publishes a stable bulk-download
#    format (check https://tdrtargets.org/releases), replace the manual
#    --tdr-csv step with a direct fetch.
