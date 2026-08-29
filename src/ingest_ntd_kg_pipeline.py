"""
ingest_ntd_kg_pipeline.py
===========================
Bridges scripts/ntd_kg_pipeline.py's output (live Open Targets/ChEMBL/
openFDA/ClinicalTrials.gov/PubMed/WHO-GHO data) into this project's
graph schema (src/schema.py), merging it with the existing hand-curated
seed data rather than replacing it.

Run scripts/ntd_kg_pipeline.py FIRST, from a machine with network
access:

    python3 scripts/ntd_kg_pipeline.py --out ntd_graph.json \\
        --csv ntd_candidates.csv --include-ntd-screens \\
        --include-clinical-evidence --include-literature --include-burden

Then run this to merge its output into data/nodes.csv / data/edges.csv:

    python3 src/ingest_ntd_kg_pipeline.py --graph ntd_graph.json --csv ntd_candidates.csv

IMPORTANT — ground truth vs. hypothesis, read before running:
------------------------------------------------------------
ntd_kg_pipeline.py's "candidate" links (drug known to modulate a target
associated with the disease) and "phenotypic_hit" links (whole-organism
screening actives) are HYPOTHESES, not confirmed treatments. They are
mapped to NEW edge types here -- CfD (candidateFor) and CpH
(phenotypicHitAgainst) -- deliberately kept separate from CtD
(Compound-treats-Disease), which is reserved for actual approved
indications. This matters: src/metapath_features.py and src/gnn_numpy.py
train on CtD as ground-truth positive labels. If candidate/hypothesis
links were mislabeled as CtD, the models would be trained to reproduce
this pipeline's own guesses rather than learn from real approved
indications -- silently corrupting evaluation.

The clinical-evidence fields (boxed_warning, disease_trial_count,
serious_ae_reports, pubmed_result_count, priority_score, etc.) from
ntd_candidates.csv are stored as edge properties on the CfD/CpH edges
(packed into the `evidence` column as a semicolon-separated string,
matching this project's existing edges.csv shape) so nothing is lost,
even though the current CSV schema doesn't have dedicated columns for
them. If you want them queryable as first-class columns, extend
edges.csv's schema and update build_graph.py / deploy_neo4j.py to
match -- that's a reasonable next step once this data is flowing.

Node type mapping:
    JSON "disease" -> Disease
    JSON "target"  -> Gene         (Open Targets targets are human genes;
                                     these are NOT PathogenGene nodes --
                                     that pool is reserved for pathogen-
                                     encoded targets, e.g. from TDR Targets)
    JSON "drug"    -> Compound

Edge type mapping:
    "associated"     (disease -> target)  -> DaG  (Disease-associates-Gene)
    "modulates"       (drug -> target)     -> CbG  (Compound-binds-Gene)
    "candidate"        (drug -> disease)    -> CfD  (NEW: hypothesis, not CtD)
    "phenotypic_hit"  (drug -> disease)    -> CpH  (NEW: hypothesis, not CtD)
"""

import argparse
import csv
import json
import os
import sys

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def load_existing(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def merge_nodes(existing_nodes, api_graph):
    """Merge API-pulled nodes into the existing node list, matching by
    name (case-insensitive) so the same real-world entity doesn't get
    duplicated if it's already in the curated seed data under a
    different id scheme."""
    xref_cols = [c for c in existing_nodes[0].keys() if c not in ("id", "name", "kind")] \
        if existing_nodes else []
    by_name = {n["name"].strip().lower(): n for n in existing_nodes}
    next_ids = {"Disease": 0, "Gene": 0, "Compound": 0}

    kind_map = {"disease": "Disease", "target": "Gene", "drug": "Compound"}
    id_prefix = {"Disease": "D:api", "Gene": "G:api", "Compound": "C:api"}
    name_to_id = {}

    added = 0
    for node in api_graph["nodes"]:
        kind = kind_map.get(node.get("type"))
        if not kind:
            continue
        key = node["id"].strip().lower()
        if key in by_name:
            name_to_id[node["id"]] = by_name[key]["id"]
            continue
        next_ids[kind] += 1
        new_id = f"{id_prefix[kind]}:{next_ids[kind]}"
        row = {"id": new_id, "name": node["id"], "kind": kind}
        for col in xref_cols:
            row[col] = ""
        existing_nodes.append(row)
        by_name[key] = row
        name_to_id[node["id"]] = new_id
        added += 1

    print(f"  merged nodes: {added} new, "
          f"{len(api_graph['nodes']) - added} matched existing entities")
    return existing_nodes, name_to_id


def merge_edges(existing_edges, api_graph, candidates_by_pair, name_to_id):
    kind_map = {
        "associated": "DaG",
        "modulates": "CbG",
        "candidate": "CfD",
        "phenotypic_hit": "CpH",
    }
    added = 0
    for link in api_graph["links"]:
        abbr = kind_map.get(link.get("kind"))
        if not abbr:
            continue
        src_id = name_to_id.get(link["source"])
        tgt_id = name_to_id.get(link["target"])
        if not src_id or not tgt_id:
            continue

        evidence_parts = [f"source=ntd_kg_pipeline_live_api"]
        if link.get("score") is not None:
            evidence_parts.append(f"assoc_score={link['score']:.3f}")

        if abbr == "CfD":
            cand = candidates_by_pair.get((link["source"], link["target"]))
            if cand:
                for field in ("priority_score", "max_phase", "boxed_warning",
                              "disease_trial_count", "serious_ae_reports",
                              "pubmed_result_count"):
                    val = cand.get(field)
                    if val not in (None, "", "None"):
                        evidence_parts.append(f"{field}={val}")

        existing_edges.append({
            "source": src_id, "relation": link["kind"], "target": tgt_id,
            "abbr": abbr, "evidence": ";".join(evidence_parts),
        })
        added += 1
    print(f"  merged edges: {added} new (CfD=hypothesis candidates, "
          f"CpH=phenotypic hits, DaG/CbG=supporting evidence)")
    return existing_edges


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--graph", required=True, help="ntd_graph.json from ntd_kg_pipeline.py")
    ap.add_argument("--csv", required=True, help="ntd_candidates.csv from ntd_kg_pipeline.py")
    args = ap.parse_args()

    with open(args.graph) as f:
        api_graph = json.load(f)

    candidates_by_pair = {}
    with open(args.csv) as f:
        for row in csv.DictReader(f):
            candidates_by_pair[(row["drug"], row["disease"])] = row

    nodes_path = os.path.join(DATA_DIR, "nodes.csv")
    edges_path = os.path.join(DATA_DIR, "edges.csv")
    existing_nodes = load_existing(nodes_path)
    existing_edges = load_existing(edges_path)

    print(f"Starting from {len(existing_nodes)} existing nodes, {len(existing_edges)} existing edges")
    print(f"Merging {len(api_graph['nodes'])} API nodes, {len(api_graph['links'])} API links...")

    existing_nodes, name_to_id = merge_nodes(existing_nodes, api_graph)
    existing_edges = merge_edges(existing_edges, api_graph, candidates_by_pair, name_to_id)

    node_cols = list(existing_nodes[0].keys())
    with open(nodes_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=node_cols)
        w.writeheader()
        w.writerows(existing_nodes)

    with open(edges_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["source", "relation", "target", "abbr", "evidence"])
        w.writeheader()
        w.writerows(existing_edges)

    print(f"\nDone. data/nodes.csv now has {len(existing_nodes)} nodes, "
          f"data/edges.csv now has {len(existing_edges)} edges.")
    print("Re-run src/pipeline.py to retrain on the enlarged graph -- note "
          "CfD/CpH edges are hypothesis-stage and are NOT used as CtD "
          "ground-truth positives (see module docstring).")


if __name__ == "__main__":
    main()
