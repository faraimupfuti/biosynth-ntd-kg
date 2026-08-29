"""
build_graph.py
===============
Loads nodes.csv / edges.csv into a networkx MultiDiGraph representing
the Biosynth NTD-extended Hetionet graph, and adds the reverse edge
for every relation (Hetionet's graph is used as undirected/bidirectional
for metapath traversal, since e.g. Gene-interacts-Gene and most
biological relations are queried in both directions).
"""

import os
import csv
import networkx as nx

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def load_graph(nodes_path=None, edges_path=None):
    nodes_path = nodes_path or os.path.join(DATA_DIR, "nodes.csv")
    edges_path = edges_path or os.path.join(DATA_DIR, "edges.csv")

    G = nx.MultiDiGraph()

    with open(nodes_path) as f:
        for row in csv.DictReader(f):
            attrs = {"name": row["name"], "kind": row["kind"]}
            # carry through any populated cross-reference ID columns
            for k, v in row.items():
                if k not in ("id", "name", "kind") and v:
                    attrs[k] = v
            G.add_node(row["id"], **attrs)

    with open(edges_path) as f:
        for row in csv.DictReader(f):
            G.add_edge(
                row["source"],
                row["target"],
                relation=row["relation"],
                abbr=row["abbr"],
                evidence=row.get("evidence", ""),
            )
            # add reverse direction with an inverse abbreviation tag so
            # metapaths can traverse either way (standard Hetionet practice)
            G.add_edge(
                row["target"],
                row["source"],
                relation=row["relation"] + "_rev",
                abbr=row["abbr"] + "_rev",
                evidence=row.get("evidence", ""),
            )
    return G


def summary(G):
    kinds = {}
    for _, data in G.nodes(data=True):
        kinds[data["kind"]] = kinds.get(data["kind"], 0) + 1
    rels = {}
    for _, _, data in G.edges(data=True):
        rels[data["abbr"]] = rels.get(data["abbr"], 0) + 1
    return kinds, rels


if __name__ == "__main__":
    G = load_graph()
    kinds, rels = summary(G)
    print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} directed edges (incl. reverses)\n")
    print("Node types:")
    for k, v in sorted(kinds.items(), key=lambda x: -x[1]):
        print(f"  {k:20s} {v}")
    print("\nEdge (relation) types (forward only shown, reverses mirror these):")
    for k, v in sorted(rels.items(), key=lambda x: -x[1]):
        if not k.endswith("_rev"):
            print(f"  {k:12s} {v}")
