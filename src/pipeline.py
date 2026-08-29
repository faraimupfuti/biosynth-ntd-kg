"""
pipeline.py
===========
End-to-end run: build graph -> extract features -> train baseline
(metapath + logistic regression) and GNN link-prediction models ->
evaluate on held-out known Compound-treats-Disease edges -> generate
a ranked list of NEW candidate drug-repurposing predictions for NTDs.

Run with:  python3 src/pipeline.py
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_graph import load_graph, summary
from metapath_features import (
    build_feature_matrix,
    compound_disease_pairs,
    label_pairs,
    METAPATH_NAMES,
)
from gnn_numpy import RelationalGNN, build_adjacency, random_search_train

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score

RNG_SEED = 42


def train_test_split_edges(pairs, y, test_frac=0.25, seed=RNG_SEED):
    rng = np.random.default_rng(seed)
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    rng.shuffle(pos_idx)
    n_test_pos = max(1, int(len(pos_idx) * test_frac))
    test_pos = pos_idx[:n_test_pos]
    train_pos = pos_idx[n_test_pos:]

    # balanced negative sampling
    rng.shuffle(neg_idx)
    train_neg = neg_idx[: len(train_pos) * 3]
    test_neg = neg_idx[len(train_pos) * 3: len(train_pos) * 3 + n_test_pos * 3]

    train_idx = np.concatenate([train_pos, train_neg])
    test_idx = np.concatenate([test_pos, test_neg])
    rng.shuffle(train_idx)
    rng.shuffle(test_idx)
    return train_idx, test_idx


def precision_at_k(y_true, scores, k):
    order = np.argsort(-scores)[:k]
    return y_true[order].sum() / k


def main():
    print("=" * 70)
    print("BIOSYNTH NTD KNOWLEDGE GRAPH — DRUG REPURPOSING PIPELINE")
    print("=" * 70)

    G = load_graph()
    kinds, rels = summary(G)
    print(f"\nGraph loaded: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges() // 2} unique edges (+ reverses)")
    print(f"Node types: {dict(kinds)}")

    pairs = compound_disease_pairs(G)
    y = label_pairs(G, pairs)
    print(f"\nTotal Compound-Disease pairs: {len(pairs)} "
          f"({y.sum()} known 'treats' edges = positive class)")

    # ---------------------------------------------------------- BASELINE
    print("\n" + "-" * 70)
    print("MODEL 1: Metapath DWPC features + Logistic Regression (Rephetio-style)")
    print("-" * 70)

    X, feat_names = build_feature_matrix(G, pairs, exclude_direct=True)
    print(f"Feature matrix: {X.shape[0]} pairs x {X.shape[1]} metapath features")
    print(f"Features: {feat_names}")

    train_idx, test_idx = train_test_split_edges(pairs, y)
    clf = LogisticRegression(max_iter=1000, class_weight="balanced")
    clf.fit(X[train_idx], y[train_idx])

    test_scores_baseline = clf.predict_proba(X[test_idx])[:, 1]
    y_test = y[test_idx]

    auc_baseline = roc_auc_score(y_test, test_scores_baseline)
    ap_baseline = average_precision_score(y_test, test_scores_baseline)
    k = min(5, y_test.sum())
    p_at_k_baseline = precision_at_k(y_test, test_scores_baseline, k) if k > 0 else float("nan")

    print(f"\nHeld-out test set: {len(test_idx)} pairs ({y_test.sum()} positive)")
    print(f"  AUROC:          {auc_baseline:.3f}")
    print(f"  Average Prec.:  {ap_baseline:.3f}")
    print(f"  Precision@{k}:    {p_at_k_baseline:.3f}")

    print("\nFeature importance (logistic regression coefficients):")
    for name, coef in sorted(zip(feat_names, clf.coef_[0]), key=lambda x: -abs(x[1])):
        print(f"  {name:35s} {coef:+.3f}")

    # ---------------------------------------------------------------- GNN
    print("\n" + "-" * 70)
    print("MODEL 2: Relational GNN (R-GCN style, numpy implementation)")
    print("-" * 70)

    relations = sorted(set(d["abbr"] for _, _, d in G.edges(data=True)))
    node_ids = list(G.nodes())
    id_to_idx = {n: i for i, n in enumerate(node_ids)}
    adj = build_adjacency(G, relations, id_to_idx)

    base_gnn = RelationalGNN(node_ids, relations, embed_dim=16, seed=RNG_SEED)

    train_pairs = [pairs[i] for i in train_idx]
    train_y = y[train_idx]
    pos_train_pairs = [p for p, lbl in zip(train_pairs, train_y) if lbl == 1]
    neg_train_pairs = [p for p, lbl in zip(train_pairs, train_y) if lbl == 0]

    print(f"Model-selecting over random restarts using {len(pos_train_pairs)} "
          f"positive / {len(neg_train_pairs)} negative training pairs...")
    best_gnn, train_auc = random_search_train(
        base_gnn, adj, pos_train_pairs, neg_train_pairs, n_trials=40
    )
    print(f"Best training AUROC across restarts: {train_auc:.3f}")

    final_embeddings = best_gnn.forward(adj)
    test_pairs = [pairs[i] for i in test_idx]
    test_scores_gnn = best_gnn.score_pairs(final_embeddings, test_pairs)

    auc_gnn = roc_auc_score(y_test, test_scores_gnn)
    ap_gnn = average_precision_score(y_test, test_scores_gnn)
    p_at_k_gnn = precision_at_k(y_test, test_scores_gnn, k) if k > 0 else float("nan")

    print(f"\nHeld-out test set (same split as baseline):")
    print(f"  AUROC:          {auc_gnn:.3f}")
    print(f"  Average Prec.:  {ap_gnn:.3f}")
    print(f"  Precision@{k}:    {p_at_k_gnn:.3f}")

    # ------------------------------------------------- CANDIDATE RANKING
    print("\n" + "-" * 70)
    print("NEW CANDIDATE PREDICTIONS (unknown pairs, ranked by ensemble score)")
    print("-" * 70)

    unknown_mask = y == 0
    unknown_pairs = [p for p, m in zip(pairs, unknown_mask) if m]
    X_unknown, _ = build_feature_matrix(G, unknown_pairs, exclude_direct=True)
    baseline_unknown_scores = clf.predict_proba(X_unknown)[:, 1]
    gnn_unknown_scores = best_gnn.score_pairs(final_embeddings, unknown_pairs)

    # simple ensemble: average of both models' calibrated rank
    def rank_normalize(s):
        order = np.argsort(s)
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(len(s)) / max(len(s) - 1, 1)
        return ranks

    ensemble = 0.5 * rank_normalize(baseline_unknown_scores) + 0.5 * rank_normalize(gnn_unknown_scores)

    top_n = 12
    top_idx = np.argsort(-ensemble)[:top_n]
    node_names = {n: d["name"] for n, d in G.nodes(data=True)}

    print(f"\nTop {top_n} candidate drug-NTD repurposing pairs:\n")
    print(f"{'Rank':<5}{'Compound':<20}{'Disease':<32}{'Baseline':<10}{'GNN':<8}{'Ensemble':<10}")
    for rank, i in enumerate(top_idx, 1):
        c, d = unknown_pairs[i]
        print(f"{rank:<5}{node_names[c]:<20}{node_names[d]:<32}"
              f"{baseline_unknown_scores[i]:<10.3f}{gnn_unknown_scores[i]:<8.3f}{ensemble[i]:<10.3f}")

    print("\n" + "=" * 70)
    print("NOTE: This ranking is generated from a small, hand-curated seed")
    print("graph (see data/build_seed_data.py) meant to validate the")
    print("pipeline end-to-end. Treat these specific rankings as a")
    print("methodology demo, NOT as validated candidates -- rerun once")
    print("bulk sources (DrugBank, TDR Targets, DisGeNET) are ingested,")
    print("and route any promising hits through literature/wet-lab")
    print("validation before including them in a grant narrative.")
    print("=" * 70)

    # ------------------------------------------------------------- SAVE
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output")
    os.makedirs(out_dir, exist_ok=True)
    import csv as _csv
    with open(os.path.join(out_dir, "candidate_predictions.csv"), "w", newline="") as f:
        w = _csv.writer(f)
        w.writerow(["rank", "compound", "disease", "baseline_score", "gnn_score", "ensemble_score"])
        for rank, i in enumerate(sorted(range(len(unknown_pairs)), key=lambda i: -ensemble[i])[:50], 1):
            c, d = unknown_pairs[i]
            w.writerow([rank, node_names[c], node_names[d],
                        round(float(baseline_unknown_scores[i]), 4),
                        round(float(gnn_unknown_scores[i]), 4),
                        round(float(ensemble[i]), 4)])
    with open(os.path.join(out_dir, "metrics.csv"), "w", newline="") as f:
        w = _csv.writer(f)
        w.writerow(["model", "auroc", "average_precision", f"precision_at_{k}"])
        w.writerow(["metapath_logreg_baseline", round(auc_baseline, 4), round(ap_baseline, 4), round(p_at_k_baseline, 4)])
        w.writerow(["relational_gnn_numpy", round(auc_gnn, 4), round(ap_gnn, 4), round(p_at_k_gnn, 4)])
    print(f"\nSaved: output/candidate_predictions.csv (top 50), output/metrics.csv")

    return {
        "baseline_auc": auc_baseline,
        "gnn_auc": auc_gnn,
        "top_candidates": [(node_names[unknown_pairs[i][0]], node_names[unknown_pairs[i][1]], float(ensemble[i]))
                            for i in top_idx],
    }


if __name__ == "__main__":
    main()
