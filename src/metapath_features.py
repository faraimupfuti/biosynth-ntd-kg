"""
metapath_features.py
=====================
Implements the metapath-based feature extraction approach from
Himmelstein et al. 2017 ("Project Rephetio"): for every Compound-Disease
pair, count the number of paths following each biologically meaningful
metapath template (e.g. Compound-binds-PathogenGene-hasGene_rev-Pathogen
-causes-Disease), then use these path counts as features for a
classifier trained to predict Compound-treats-Disease edges.

This is deliberately the SIMPLE, INTERPRETABLE baseline: every feature
has a direct biological reading ("this drug binds a target in a pathway
shared with a gene associated with this disease"), which matters for
explaining predictions to grant reviewers / collaborators, and for
sanity-checking before layering on less-interpretable embedding/GNN
models (see gnn_numpy.py).
"""

import itertools
import networkx as nx
import numpy as np

# Hand-picked metapath templates (sequences of edge abbreviations to
# traverse from a Compound node to a Disease node). These mirror the
# style of Rephetio's DWPC (degree-weighted path count) metapaths,
# extended with the new Pathogen/PathogenGene/Vector node types.
METAPATH_TEMPLATES = [
    ["CtD"],                                   # direct known treatment (only for train positives; excluded at inference)
    ["CbPG", "PhPG_rev", "PcD"],                # Compound -> target -> pathogen -> disease
    ["CbPG", "PGpPW", "PGpPW_rev", "PhPG_rev", "PcD"],  # shared pathway route
    ["CbPG", "PGoG", "GpPW", "PGpPW_rev", "PhPG_rev", "PcD"],  # via human ortholog + pathway
    ["CbG", "GpPW", "PGpPW_rev", "PhPG_rev", "PcD"],   # human-target compound -> pathway -> pathogen -> disease
    ["PCiC_rev", "CbPG", "PhPG_rev", "PcD"],    # same pharmacologic class, shared pathogen target
]

METAPATH_NAMES = [
    "direct_CtD",
    "target_to_pathogen_to_disease",
    "shared_pathway_route",
    "human_ortholog_pathway_route",
    "human_target_pathway_route",
    "same_drug_class_shared_target",
]


def _dwpc(G, source, target, path_abbrs, damping=0.4):
    """
    Degree-Weighted Path Count between source and target following an
    exact sequence of edge-type abbreviations. Down-weights paths that
    pass through high-degree ("hub") nodes, per Rephetio's DWPC metric,
    so promiscuous nodes (e.g. a pathway with hundreds of genes) don't
    dominate the signal.
    """
    frontier = {source: 1.0}
    for abbr in path_abbrs:
        next_frontier = {}
        for node, weight in frontier.items():
            for _, nbr, data in G.out_edges(node, data=True):
                if data["abbr"] != abbr:
                    continue
                deg = max(G.out_degree(nbr), 1)
                w = weight / (deg ** damping)
                next_frontier[nbr] = next_frontier.get(nbr, 0.0) + w
        frontier = next_frontier
        if not frontier:
            break
    return frontier.get(target, 0.0)


def compound_disease_pairs(G):
    compounds = [n for n, d in G.nodes(data=True) if d["kind"] == "Compound"]
    diseases = [n for n, d in G.nodes(data=True) if d["kind"] == "Disease"]
    return list(itertools.product(compounds, diseases))


def build_feature_matrix(G, pairs, exclude_direct=True):
    """
    For each (compound, disease) pair, compute DWPC along every metapath
    template (excluding the trivial 'direct_CtD' feature when
    exclude_direct=True, since that would leak the label itself).
    Returns X (features), and the ordered feature names used.
    """
    templates = METAPATH_TEMPLATES
    names = METAPATH_NAMES
    if exclude_direct:
        templates = templates[1:]
        names = names[1:]

    X = np.zeros((len(pairs), len(templates)))
    for i, (c, d) in enumerate(pairs):
        for j, path_abbrs in enumerate(templates):
            X[i, j] = _dwpc(G, c, d, path_abbrs)
    return X, names


def label_pairs(G, pairs):
    """Label a (compound, disease) pair 1 if a known CtD edge exists."""
    y = np.zeros(len(pairs), dtype=int)
    for i, (c, d) in enumerate(pairs):
        if G.has_edge(c, d):
            edata = G.get_edge_data(c, d)
            if any(v["abbr"] == "CtD" for v in edata.values()):
                y[i] = 1
    return y
