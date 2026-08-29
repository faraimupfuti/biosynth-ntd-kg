"""
deploy_neo4j.py
================
Loads data/nodes.csv + data/edges.csv into a running Neo4j instance,
so the graph can be explored via Neo4j Browser / Bloom instead of
only through the Python/networkx pipeline.

Setup:
    1. Spin up a Neo4j instance -- easiest is Neo4j Aura Free
       (neo4j.com/cloud/aura), a managed cloud instance, no server
       setup required. Alternative: Neo4j Desktop for a self-hosted
       instance on a shared machine.
    2. pip install neo4j
    3. Set the three env vars below (Aura gives you these on
       instance creation; for local Desktop, default user is "neo4j"
       and you set the password on first launch).
    4. python3 src/deploy_neo4j.py

Design notes:
    - Every node gets a label matching its `kind` column (Compound,
      Disease, Pathogen, etc.) plus a shared `:Entity` label, so you
      can query either a specific type (MATCH (c:Compound)) or
      everything (MATCH (n:Entity)).
    - Every edge becomes a relationship whose TYPE is the `abbr`
      column (CtD, CbPG, PcD, ...) so Cypher queries read naturally,
      e.g.:
          MATCH (c:Compound)-[:CtD]->(d:Disease) RETURN c, d
    - Reverse edges (the `_rev` suffixed ones networkx uses internally
      for undirected-style traversal) are NOT loaded into Neo4j --
      Neo4j relationships are natively traversable in both directions
      with MATCH (a)-[:REL]-(b), so they're redundant here and would
      just double the relationship count.
    - Cross-reference IDs (doid, mesh_id, ncbi_gene_id, etc. -- see
      XREF_FIELDS in src/schema.py) are loaded as node properties
      whenever populated, so they're queryable/filterable directly:
          MATCH (d:Disease) WHERE d.mesh_id = 'D014355' RETURN d
"""

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schema import XREF_FIELDS

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

NEO4J_URI = os.environ.get("NEO4J_URI", "neo4j+s://<your-instance>.databases.neo4j.io")
# Aura's exported credentials file uses NEO4J_USERNAME, not NEO4J_USER --
# accept either so you can `source` that file directly without renaming
# anything in it.
NEO4J_USER = os.environ.get("NEO4J_USER") or os.environ.get("NEO4J_USERNAME", "neo4j")
NEO4J_PASSWORD = os.environ.get("NEO4J_PASSWORD", "")
# Aura Free instances get a named database (matching the instance ID),
# not always "neo4j" -- pass it through explicitly rather than relying
# on the driver's default, which can silently point at the wrong database.
NEO4J_DATABASE = os.environ.get("NEO4J_DATABASE")  # None = driver default

ALL_XREF_COLUMNS = sorted({col for cols in XREF_FIELDS.values() for col in cols})


def load_csv_rows(filename):
    with open(os.path.join(DATA_DIR, filename)) as f:
        return list(csv.DictReader(f))


def push_to_neo4j():
    try:
        from neo4j import GraphDatabase
    except ImportError:
        raise ImportError("pip install neo4j  # official Neo4j Python driver")

    nodes = load_csv_rows("nodes.csv")
    edges = load_csv_rows("edges.csv")

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
    session_kwargs = {"database": NEO4J_DATABASE} if NEO4J_DATABASE else {}

    with driver.session(**session_kwargs) as session:
        # constraint for fast MERGE lookups + de-dup safety on re-runs
        session.run("CREATE CONSTRAINT entity_id IF NOT EXISTS FOR (n:Entity) REQUIRE n.id IS UNIQUE")

        print(f"Loading {len(nodes)} nodes...")
        for row in nodes:
            props = {"id": row["id"], "name": row["name"]}
            for col in ALL_XREF_COLUMNS:
                if row.get(col):
                    props[col] = row[col]
            # label = the node's kind (e.g. Compound, Disease), plus shared :Entity label
            session.run(
                f"MERGE (n:Entity:{row['kind']} {{id: $id}}) SET n += $props",
                id=row["id"], props=props,
            )

        print(f"Loading {len(edges)} edges...")
        for row in edges:
            if row["abbr"].endswith("_rev"):
                continue  # skip reverse edges -- redundant in Neo4j, see module docstring
            rel_type = row["abbr"]
            session.run(
                f"""
                MATCH (a:Entity {{id: $source}}), (b:Entity {{id: $target}})
                MERGE (a)-[r:{rel_type}]->(b)
                SET r.evidence = $evidence
                """,
                source=row["source"], target=row["target"], evidence=row.get("evidence", ""),
            )

    driver.close()
    print("\nDone. Open Neo4j Browser and try:")
    print("  MATCH (c:Compound)-[:CtD]->(d:Disease) RETURN c, d LIMIT 25")
    print("  MATCH (c:Compound)-[:CbPG]->(t:PathogenGene)<-[:PhPG]-(p:Pathogen)-[:PcD]->(d:Disease)")
    print("  RETURN c, t, p, d LIMIT 25")


if __name__ == "__main__":
    if not NEO4J_PASSWORD:
        print("Set NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD environment variables first")
        print("(Aura gives you these when you create an instance; for Neo4j")
        print("Desktop, default user is 'neo4j' with the password you set locally).")
        sys.exit(1)
    push_to_neo4j()
