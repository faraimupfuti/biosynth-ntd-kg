"""
schema.py
=========
Defines the metagraph (node types + edge types) for the Biosynth NTD
knowledge graph. This is Hetionet's original schema, extended with
node/edge types needed to represent Neglected Tropical Disease (NTD)
biology: causative pathogens, transmission vectors, and pathogen-level
drug targets (which is where a lot of NTD repurposing signal lives,
since many NTD drugs act on parasite/bacterial targets that have no
direct human ortholog in the original Hetionet graph).

NODE TYPES
----------
Original Hetionet:
    Compound (C), Disease (D), Gene (G), Anatomy (A), Pathway (PW),
    Biological Process (BP), Molecular Function (MF),
    Cellular Component (CC), Pharmacologic Class (PC),
    Side Effect (SE), Symptom (S)

NTD extension (new):
    Pathogen (P)   - the causative organism (parasite, bacterium, virus)
    Vector (V)     - the transmission vector (e.g. triatomine bug, blackfly)
    PathogenGene (PG) - a pathogen-encoded gene/protein (distinct pool from
                         human Gene nodes, since most NTD drug targets are
                         organism-specific, not human)

EDGE TYPES (metaedges)
-----------------------
Original Hetionet (subset relevant to repurposing):
    CtD  Compound-treats-Disease
    CpD  Compound-palliates-Disease
    CbG  Compound-binds-Gene
    CuG  Compound-upregulates-Gene
    CdG  Compound-downregulates-Gene
    DaG  Disease-associates-Gene
    DlA  Disease-localizes-Anatomy
    GiG  Gene-interacts-Gene
    GpPW Gene-participates-Pathway
    GpBP Gene-participates-BiologicalProcess
    GpMF Gene-participates-MolecularFunction
    GpCC Gene-participates-CellularComponent
    PCiC PharmacologicClass-includes-Compound

NTD extension (new):
    PcD  Pathogen-causes-Disease
    VtP  Vector-transmits-Pathogen
    PhPG Pathogen-hasGene-PathogenGene   (pathogen genome membership)
    CbPG Compound-binds-PathogenGene     (the key repurposing edge type:
                                          which drugs bind which pathogen
                                          targets)
    PGoG PathogenGene-orthologOf-Gene    (ortholog link back to human/
                                          model-organism gene space, so
                                          human pathway data can still
                                          inform predictions)
    PGpPW PathogenGene-participates-Pathway

This schema is intentionally close to Hetionet's original design so
that Hetionet's published metapath methodology (Himmelstein et al.,
2017, "Systematic integration of biomedical knowledge prioritizes
drugs for repurposing", eLife) transfers directly, with the pathogen
layer added as new metapaths for NTD-specific reasoning.
"""

NODE_TYPES = [
    "Compound",
    "Disease",
    "Gene",
    "Anatomy",
    "Pathway",
    "BiologicalProcess",
    "PharmacologicClass",
    "Pathogen",
    "Vector",
    "PathogenGene",
    "SideEffect",   # NEW: enables SIDER integration (sideeffects.embl.de)
]

EDGE_TYPES = {
    # (source_type, relation, target_type): abbreviation
    ("Compound", "treats", "Disease"): "CtD",
    ("Compound", "binds", "Gene"): "CbG",
    ("Compound", "binds", "PathogenGene"): "CbPG",
    ("Compound", "causes", "SideEffect"): "CcSE",   # NEW: SIDER edge type
    ("Disease", "associates", "Gene"): "DaG",
    ("Disease", "localizes", "Anatomy"): "DlA",
    ("Gene", "interacts", "Gene"): "GiG",
    ("Gene", "participates", "Pathway"): "GpPW",
    ("Gene", "participates", "BiologicalProcess"): "GpBP",
    ("PharmacologicClass", "includes", "Compound"): "PCiC",
    ("Pathogen", "causes", "Disease"): "PcD",
    ("Vector", "transmits", "Pathogen"): "VtP",
    ("Pathogen", "hasGene", "PathogenGene"): "PhPG",
    ("PathogenGene", "orthologOf", "Gene"): "PGoG",
    ("PathogenGene", "participates", "Pathway"): "PGpPW",
}

# Which edge type is the prediction TARGET for repurposing
TARGET_EDGE = ("Compound", "treats", "Disease")
TARGET_ABBR = "CtD"

# Cross-reference ID fields, by node kind, populated by
# src/enrich_from_sources.py once network access to the source is
# available. Empty string until enriched. Keeping these as separate
# optional columns (rather than folding into `id`) means the graph
# structure never depends on any one external ontology being present.
XREF_FIELDS = {
    "Disease": ["doid", "mesh_id", "umls_cui"],
    "Gene": ["ncbi_gene_id", "uniprot_id", "go_ids"],
    "PathogenGene": ["uniprot_id", "go_ids"],
    "Compound": ["drugbank_id", "drugcentral_id", "mesh_id"],
    "Pathway": ["reactome_id", "wikipathways_id"],
    "Anatomy": ["uberon_id"],
    "SideEffect": ["umls_cui", "meddra_id"],
}

