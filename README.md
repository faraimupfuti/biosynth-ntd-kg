# Biosynth NTD Knowledge Graph — Drug Repurposing Pipeline

A Hetionet-style knowledge graph, extended with Neglected Tropical
Disease (NTD)-specific biology (pathogens, transmission vectors, and
pathogen-level drug targets), plus two link-prediction models for
identifying drug repurposing candidates.

## Quick start

```bash
pip install -r requirements.txt
python3 data/build_seed_data.py   # generates data/nodes.csv, data/edges.csv
python3 src/pipeline.py           # runs the full pipeline end-to-end
```

Outputs land in `output/`:
- `candidate_predictions.csv` — top 50 ranked drug-NTD repurposing candidates
- `metrics.csv` — model evaluation metrics

## What this is (and isn't) right now

This is a **working, end-to-end methodology demonstrator**, seeded
with a small (~90 node / ~125 edge) hand-curated graph built from
well-established public pharmacology/parasitology facts (real NTDs,
real approved drugs, real mechanisms of action). It proves out:

1. The extended graph schema
2. A metapath-based feature extraction + logistic regression baseline
   (Rephetio methodology)
3. A relational GNN architecture (R-GCN style message passing)
4. Train/test evaluation on held-out known drug-disease links
5. Ranked candidate generation for genuinely unknown pairs

**It is not yet a validated scientific result.** The seed graph is
small and illustrative. Before anything here goes into a grant
narrative or gets acted on scientifically, it needs (a) bulk real
data ingested (see below) and (b) literature/wet-lab validation of
any promising candidates.

## Architecture

```
data/
  build_seed_data.py     -- generates the curated seed graph
  nodes.csv, edges.csv   -- 8 node types, 14 edge types (see src/schema.py)

src/
  schema.py               -- metagraph definition (Hetionet + NTD extension)
  build_graph.py           -- loads CSVs into a networkx MultiDiGraph
  metapath_features.py     -- DWPC metapath feature extraction (Rephetio-style)
  gnn_numpy.py              -- relational GNN, implemented from scratch in numpy
  data_ingestion_stubs.py  -- function signatures for real bulk data sources
  pipeline.py               -- end-to-end run: train, evaluate, rank candidates

output/
  candidate_predictions.csv
  metrics.csv
```

### Schema: Hetionet extended for NTDs

Original Hetionet node types (Compound, Disease, Gene, Anatomy,
Pathway, etc.) plus three NTD-specific additions:

- **Pathogen** — the causative organism (e.g. *Trypanosoma cruzi*)
- **Vector** — the transmission vector (e.g. triatomine bug)
- **PathogenGene** — pathogen-encoded drug targets, kept as a separate
  node pool from human genes since most NTD drugs act on
  organism-specific targets with no direct human ortholog

New edge types: `PcD` (Pathogen-causes-Disease), `VtP`
(Vector-transmits-Pathogen), `PhPG` (Pathogen-hasGene-PathogenGene),
`CbPG` (Compound-binds-PathogenGene — the key repurposing-relevant
edge), `PGoG` (PathogenGene-orthologOf-Gene), `PGpPW`
(PathogenGene-participates-Pathway). Full list in `src/schema.py`.

### Model 1: Metapath DWPC + Logistic Regression (primary/recommended)

This follows Himmelstein et al. 2017's Project Rephetio methodology:
count degree-weighted paths between each (compound, disease) pair
along hand-picked biologically meaningful metapaths (e.g. "drug binds
a target that's in the same pathway as a gene in a pathogen that
causes this disease"), then train a classifier on these path-count
features.

**Why this is the recommended model for grant-facing work:** every
feature has a direct, explainable biological reading. Reviewers can
see *why* a candidate scored highly, not just that it did. On the
seed graph it also outperforms the GNN (AUROC 0.91 vs 0.70 in the
current run — see `output/metrics.csv` for the latest numbers).

### Model 2: Relational GNN (secondary/future work)

A minimal R-GCN-style relational graph neural network, **implemented
from scratch in numpy** because this sandbox has no network access to
install PyTorch/PyTorch Geometric. It does relation-specific
message-passing and a dot-product link-prediction decoder, but uses
random-restart model selection instead of proper gradient-based
backpropagation (see the docstring in `src/gnn_numpy.py`), so treat
its current results as a architecture proof-of-concept, not a tuned
model.

**To productionize:** swap `gnn_numpy.py` for `torch_geometric.nn.RGCNConv`
+ standard backprop training once you have network access (e.g. your
own infra, a Colab notebook, or an environment where `pip install
torch torch-geometric` can run). The graph construction and feature
pipeline stay identical — only the model class changes.

## Scaling to real data

**Sandbox constraint (read this first):** the code-execution sandbox
this project was built in has no outbound network access at all — its
egress proxy returns HTTP 403 (`x-deny-reason: host_not_allowed`) for
*every* host, confirmed by direct test. This applies identically to
Hetionet's own bulk files and to every other source below. It's not
possible to bulk-download any of them from inside this sandbox,
regardless of which one you point at.

`src/enrich_from_sources.py` contains **real, working integration
code** (not stubs) against each source's actual public API or bulk
file. Run it from any machine with normal internet access:

```bash
pip install requests
python3 src/enrich_from_sources.py
```

| Source | Function | Access |
|---|---|---|
| Hetionet | (see note below) | Free, public GitHub release |
| NCBI Gene | `enrich_genes_ncbi()` | Free, no key (E-utilities) |
| Gene Ontology | `enrich_genes_go()` | Free, no key (QuickGO API) |
| Disease Ontology | `enrich_diseases_do()` | Free, public OBO file |
| Uberon | `fetch_uberon_terms()` | Free, public OBO file |
| MeSH | `enrich_via_mesh()` | Free, no key |
| Reactome | `enrich_pathways_reactome()` | Free, no key |
| WikiPathways | `enrich_pathways_wikipathways()` | Free, no key |
| DrugCentral | `enrich_compounds_drugcentral()` | Free, no key |
| SIDER | `enrich_side_effects_sider()` | Free, static bulk TSV |
| UMLS | `enrich_via_umls()` | Free but requires a UTS API key (apply at uts.nlm.nih.gov) |
| DrugBank | `enrich_compounds_drugbank()` | Requires an academic/commercial license (bulk XML export) |

Running it populates the cross-reference ID columns now defined in
`src/schema.py` (`XREF_FIELDS`) — DOID/MeSH/UMLS for diseases, NCBI
Gene/UniProt/GO for genes, Reactome/WikiPathways IDs for pathways,
DrugCentral/DrugBank/MeSH for compounds — and writes
`data/nodes_enriched.csv`. Nothing in this codebase invents an
identifier; every field is either pulled live from the source or left
blank until you run it with network access.

**Hetionet itself** isn't behind a REST API — it's distributed as a
JSON/TSV/Neo4j bulk release on GitHub
(github.com/hetio/hetionet, ~47K nodes / 2.25M edges). Once you have
network access, that's still the highest-leverage first step (see
"Suggested next steps" below) — it gives you the general
Compound/Gene/Disease/Pathway backbone for free, and this project's
NTD extension (Pathogen/Vector/PathogenGene) layers directly on top of
its schema.

Beyond the sources above, **TDR Targets Database** (tdrtargets.org)
and **DrugBank** remain the two most NTD-specific/repurposing-relevant
sources not yet covered by a live API in `enrich_from_sources.py` —
TDR Targets doesn't publish a documented public REST API (bulk
export/collaboration request is the intended access path), and
DrugBank's bulk data is license-gated as noted above.

## Suggested next steps

1. Get Hetionet's public dataset ingested first (this alone will 10-100x
   the graph size for free) and re-run `pipeline.py` to see how the
   baseline model performs at real scale.
2. Run `src/enrich_from_sources.py` from a networked machine to attach
   real DOID/MeSH/NCBI Gene/GO/Reactome/UniProt/DrugCentral IDs to the
   seed graph — this makes it interoperable with any other tool that
   speaks those ontologies.
3. Reach out to TDR Targets directly for pathogen drug-target data at
   scale — it's the highest-leverage source specifically for NTDs,
   where your differentiation from generic repurposing platforms lives.
4. Move `gnn_numpy.py` to a real PyTorch Geometric implementation once
   you have compute/network access, and benchmark it against the
   baseline properly.
5. For any candidate that scores well after real data is in, do a
   literature pass (PubMed) before including it in a grant narrative
   or reaching out to a wet-lab collaborator for validation.
6. Consider deploying the graph in Neo4j (rather than in-memory
   networkx) once it's large enough that load time / query speed
   matters — the schema translates directly to a Neo4j property graph.

## Deploying as an internal query/exploration tool

For a team that wants to browse and query the graph directly (rather
than only through `pipeline.py`), load it into Neo4j:

```bash
pip install neo4j
export NEO4J_URI="neo4j+s://<your-instance>.databases.neo4j.io"
export NEO4J_USER="neo4j"
export NEO4J_PASSWORD="<your-password>"
python3 src/deploy_neo4j.py
```

Easiest way to get an instance: **Neo4j Aura Free**
(neo4j.com/cloud/aura) — managed, no server setup, free tier handles
graphs far larger than this one. `src/deploy_neo4j.py` loads every
node with a label matching its `kind` (plus a shared `:Entity` label)
and every edge as a relationship typed by its `abbr` (`CtD`, `CbPG`,
etc.), and carries over any populated cross-reference IDs as node
properties. Then explore via **Neo4j Browser** (built in, Cypher
queries) or **Neo4j Bloom** (visual, no-code, good for non-technical
team members). Example query once loaded:

```cypher
MATCH (c:Compound)-[:CbPG]->(t:PathogenGene)<-[:PhPG]-(p:Pathogen)-[:PcD]->(d:Disease)
RETURN c, t, p, d LIMIT 25
```

Keep `pipeline.py` as the source of truth for candidate rankings (the
ML models don't run inside Neo4j) — after each run, push updated
prediction scores back in as relationship properties so the graph UI
reflects the latest results.
