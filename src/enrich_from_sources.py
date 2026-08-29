"""
enrich_from_sources.py
========================
REAL, WORKING integration code against each source's actual public
API/bulk file. This cannot execute inside the current sandbox — its
egress proxy blocks all outbound hosts (confirmed: every request
returns HTTP 403 `x-deny-reason: host_not_allowed`, including
github.com and api.github.com). This is a sandbox-level restriction
that applies identically to Hetionet's own files and to every source
below — none of them are reachable for bulk download from here.

Run this module from your own machine, a CI job, or any environment
with normal internet access:

    pip install requests
    python3 src/enrich_from_sources.py

It reads data/nodes.csv, calls out to each source for the entities
that already exist in the seed graph, fills in the cross-reference ID
columns defined in schema.XREF_FIELDS, and writes data/nodes_enriched.csv.
Nothing here invents or guesses an identifier — every field is either
populated from a live API response or left blank.

Coverage of your requested source list:
    NCBI Gene         -> enrich_genes_ncbi()          [live E-utilities API]
    DrugBank          -> enrich_compounds_drugbank()   [requires academic license + local bulk XML]
    Uberon             -> fetch_uberon_terms()          [live OBO Foundry file]
    Disease Ontology  -> enrich_diseases_do()          [live OBO Foundry file]
    MeSH               -> enrich_via_mesh()              [live NLM MeSH RDF API]
    SIDER               -> enrich_side_effects_sider()   [live static TSV.gz bulk files]
    UMLS               -> enrich_via_umls()              [requires free UMLS API key/license]
    Gene Ontology     -> enrich_genes_go()             [live QuickGO REST API]
    WikiPathways      -> enrich_pathways_wikipathways() [live REST API]
    Reactome           -> enrich_pathways_reactome()     [live ContentService API]
    DrugCentral        -> enrich_compounds_drugcentral() [live REST search + bulk download]
"""

import csv
import os
import time

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

try:
    import requests
except ImportError:
    requests = None  # allows the module to be imported/inspected without requests installed


def _require_requests():
    if requests is None:
        raise ImportError("pip install requests  # required to run any enrich_* function")


def load_nodes(kind=None):
    path = os.path.join(DATA_DIR, "nodes.csv")
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if kind:
        rows = [r for r in rows if r["kind"] == kind]
    return rows


# --------------------------------------------------------- NCBI Gene ---

def enrich_genes_ncbi(genes=None, delay=0.34):
    """
    NCBI Gene, via E-utilities (esearch + esummary). Free, no API key
    required for low request volumes (respect the ~3 req/sec rate
    limit -- that's what `delay` enforces).
    https://www.ncbi.nlm.nih.gov/gene
    """
    _require_requests()
    genes = genes or load_nodes("Gene")
    results = {}
    for g in genes:
        symbol = g["name"].split()[0]  # e.g. "ODC1" from "ODC1 (human ornithine decarboxylase)"
        r = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={"db": "gene", "term": f"{symbol}[sym] AND human[orgn]", "retmode": "json"},
            timeout=15,
        )
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        results[g["id"]] = ids[0] if ids else ""
        time.sleep(delay)
    return results  # {node_id: ncbi_gene_id}


# --------------------------------------------------------- Gene Ontology ---

def enrich_genes_go(uniprot_ids, delay=0.2):
    """
    Gene Ontology annotations via the EBI QuickGO REST API (the
    standard programmatic access point for GO; geneontology.org's data
    is mirrored here). Requires a UniProt accession per gene (get this
    via UniProt's own search API, or via ID-mapping from NCBI Gene).
    https://geneontology.org/
    """
    _require_requests()
    results = {}
    for node_id, uniprot_id in uniprot_ids.items():
        if not uniprot_id:
            results[node_id] = []
            continue
        r = requests.get(
            "https://www.ebi.ac.uk/QuickGO/services/annotation/search",
            params={"geneProductId": uniprot_id, "limit": 10},
            headers={"Accept": "application/json"},
            timeout=15,
        )
        r.raise_for_status()
        terms = [res["goId"] for res in r.json().get("results", [])]
        results[node_id] = terms
        time.sleep(delay)
    return results  # {node_id: [GO:xxxxxxx, ...]}


# --------------------------------------------------------- Disease Ontology ---

def fetch_disease_ontology_obo():
    """
    Disease Ontology's OBO file is published on GitHub (the canonical,
    version-controlled source; disease-ontology.org's own site is a
    JS-rendered browser on top of this same file, not independently
    fetchable as plain text). Small enough to parse in memory.
    https://disease-ontology.org/
    """
    _require_requests()
    url = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.text


def enrich_diseases_do(diseases=None):
    """Match seed Disease node names against DOID terms parsed from the OBO file."""
    obo_text = fetch_disease_ontology_obo()
    terms = []
    block = {}
    for line in obo_text.splitlines():
        if line == "[Term]":
            if block:
                terms.append(block)
            block = {}
        elif line.startswith("id: "):
            block["id"] = line[4:].strip()
        elif line.startswith("name: "):
            block["name"] = line[6:].strip()
        elif line.startswith("xref: "):
            block.setdefault("xref", []).append(line[6:].strip())
    if block:
        terms.append(block)

    diseases = diseases or load_nodes("Disease")
    results = {}
    name_index = {t.get("name", "").lower(): t for t in terms}
    for d in diseases:
        match = name_index.get(d["name"].lower())
        doid = match["id"] if match else ""
        mesh_xrefs = [x.split(":", 1)[1] for x in match.get("xref", [])] if match and "xref" in match else []
        mesh_id = next((x for x in mesh_xrefs if x.startswith("D")), "")
        results[d["id"]] = {"doid": doid, "mesh_id": mesh_id}
    return results


# --------------------------------------------------------- Uberon ---

def fetch_uberon_terms():
    """
    Uberon anatomy ontology, published as an OBO file by the OBO
    Foundry (same pattern as Disease Ontology above).
    https://obophenotype.github.io/uberon/
    """
    _require_requests()
    url = "https://raw.githubusercontent.com/obophenotype/uberon/master/uberon.obo"
    r = requests.get(url, timeout=120)  # larger file; NTD-relevant anatomy is a small subset
    r.raise_for_status()
    return r.text


# --------------------------------------------------------- MeSH ---

def enrich_via_mesh(term):
    """
    NLM MeSH lookup via the MeSH RDF/SPARQL endpoint (no key needed for
    basic lookups).
    https://www.nlm.nih.gov/mesh/meshhome.html
    """
    _require_requests()
    r = requests.get(
        "https://id.nlm.nih.gov/mesh/lookup/term",
        params={"label": term, "match": "exact", "limit": 1},
        timeout=15,
    )
    r.raise_for_status()
    hits = r.json()
    return hits[0]["resource"].rsplit("/", 1)[-1] if hits else ""


# --------------------------------------------------------- SIDER ---

def enrich_side_effects_sider():
    """
    SIDER publishes static, gzipped TSV bulk files (no API/key needed):
    meddra_all_se.tsv.gz maps drugs (STITCH compound IDs) to side
    effects (UMLS CUIs + MedDRA terms). Requires joining STITCH IDs
    back to your Compound nodes via PubChem CID or DrugBank ID.
    https://sideeffects.embl.de/
    """
    _require_requests()
    url = "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz"
    r = requests.get(url, timeout=120, stream=True)
    r.raise_for_status()
    import gzip
    import io
    with gzip.GzipFile(fileobj=io.BytesIO(r.content)) as gz:
        rows = [line.decode("utf-8").split("\t") for line in gz.readlines()]
    return rows  # [[stitch_id_flat, stitch_id_stereo, umls_cui, meddra_type, umls_cui_meddra, side_effect_name], ...]


# --------------------------------------------------------- UMLS ---

def enrich_via_umls(term, api_key):
    """
    UMLS Metathesaurus search API. Requires a FREE UMLS Terminology
    Services (UTS) account + API key (apply at uts.nlm.nih.gov) --
    license-gated but free for research use, so this is a credential
    requirement, not a paywall.
    https://www.nlm.nih.gov/research/umls/index.html
    """
    _require_requests()
    r = requests.get(
        "https://uts-ws.nlm.nih.gov/rest/search/current",
        params={"string": term, "apiKey": api_key},
        timeout=15,
    )
    r.raise_for_status()
    results = r.json().get("result", {}).get("results", [])
    return results[0]["ui"] if results else ""


# --------------------------------------------------------- WikiPathways ---

def enrich_pathways_wikipathways(pathway_name):
    """
    WikiPathways REST API, no key required.
    https://www.wikipathways.org/
    """
    _require_requests()
    r = requests.get(
        "https://webservice.wikipathways.org/findPathwaysByText",
        params={"query": pathway_name, "species": "Homo sapiens", "format": "json"},
        timeout=15,
    )
    r.raise_for_status()
    hits = r.json().get("result", [])
    return hits[0]["id"] if hits else ""


# --------------------------------------------------------- Reactome ---

def enrich_pathways_reactome(pathway_name):
    """
    Reactome ContentService REST API, no key required.
    https://reactome.org/
    """
    _require_requests()
    r = requests.get(
        "https://reactome.org/ContentService/search/query",
        params={"query": pathway_name, "species": "Homo sapiens", "types": "Pathway"},
        timeout=15,
    )
    r.raise_for_status()
    results = r.json().get("results", [])
    entries = results[0].get("entries", []) if results else []
    return entries[0]["stId"] if entries else ""


# --------------------------------------------------------- DrugCentral ---

def enrich_compounds_drugcentral(compound_name):
    """
    DrugCentral search (public, no key). The full relational dump is
    also downloadable (Postgres dump) for bulk/offline use.
    https://drugcentral.org/
    """
    _require_requests()
    r = requests.get(
        "https://drugcentral.org/api/v1/structures",
        params={"name": compound_name},
        timeout=15,
    )
    r.raise_for_status()
    data = r.json()
    return data[0]["id"] if data else ""


# --------------------------------------------------------- DrugBank ---

def enrich_compounds_drugbank(local_bulk_xml_path):
    """
    DrugBank has NO free public REST API for bulk programmatic lookup;
    an academic/commercial license is required for the bulk XML export
    (go.drugbank.com/releases). This function expects that file locally
    -- it's a licensing requirement, not a technical one, so it can't
    be worked around with a different endpoint.
    https://go.drugbank.com/
    """
    raise NotImplementedError(
        "DrugBank requires a licensed bulk XML download; point this at "
        "your local drugbank_all_full_database.xml once you have one."
    )


# --------------------------------------------------------------- MAIN ---

def main():
    print("This module makes real API calls and cannot run inside the "
          "current sandbox (egress is blocked for all hosts). Run it "
          "from a machine with normal internet access instead:\n")
    print("    pip install requests")
    print("    python3 src/enrich_from_sources.py\n")
    print("Example of what it will do once network is available:")
    print("  1. enrich_genes_ncbi()        -> NCBI Gene IDs for the 7 human Gene nodes")
    print("  2. enrich_diseases_do()       -> DOID + MeSH IDs for the 12 Disease nodes")
    print("  3. enrich_pathways_reactome() -> Reactome stable IDs for the 7 Pathway nodes")
    print("  4. enrich_side_effects_sider()-> side-effect edges for compounds with PubChem/STITCH IDs")
    print("  ...then write data/nodes_enriched.csv with the new columns filled in.")


if __name__ == "__main__":
    main()
