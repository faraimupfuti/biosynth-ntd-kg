"""
build_seed_data.py
===================
Generates a curated SEED knowledge graph (nodes.csv, edges.csv) for the
Biosynth NTD-extended Hetionet graph.

IMPORTANT — read this before using in production:
---------------------------------------------------
This seed graph is hand-curated from well-established public domain
pharmacology/parasitology knowledge (WHO NTD list, standard drug
mechanism-of-action facts). It is meant to (a) prove out the schema,
pipeline and ML methodology end-to-end, and (b) give you a small,
inspectable, human-verifiable graph to sanity check before scaling up.

It is NOT a substitute for bulk data ingestion. For a production
system, replace/augment this module with loaders that pull from:
    - Hetionet's own public data (github.com/hetio/hetionet) for the
      general Compound/Gene/Disease/Pathway backbone
    - DrugBank (compound-target binding, CbG/CbPG)
    - TDR Targets Database (tdrtargets.org) — purpose-built for NTD
      pathogen drug target data, ideal for PathogenGene nodes
    - DisGeNET / OMIM (disease-gene associations, DaG)
    - Reactome / KEGG (pathway membership, GpPW / PGpPW)
    - UniProt + OrthoDB / OrthoMCL (pathogen-gene to human-gene
      ortholog mapping, PGoG)
    - WHO NTD portal + CDC (Pathogen-causes-Disease, Vector-transmits-
      Pathogen ground truth)
See src/data_ingestion_stubs.py for loader function signatures to fill
in once you have API/bulk-file access.
"""

import csv
import os

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- NODES ---

nodes = []

# Union of all XREF_FIELDS values (see src/schema.py) so nodes.csv has
# a stable column set regardless of node kind; blank where not
# applicable or not yet verified.
XREF_COLUMNS = [
    "doid", "mesh_id", "umls_cui", "ncbi_gene_id", "uniprot_id", "go_ids",
    "drugbank_id", "drugcentral_id", "reactome_id", "wikipathways_id",
    "uberon_id", "meddra_id",
]


def add_node(node_id, name, kind, **xrefs):
    row = {"id": node_id, "name": name, "kind": kind}
    for col in XREF_COLUMNS:
        row[col] = xrefs.get(col, "")
    nodes.append(row)


# Diseases (WHO NTD list, subset)
diseases = [
    ("D:Chagas", "Chagas disease"),
    ("D:HAT", "Human African trypanosomiasis"),
    ("D:Leish", "Leishmaniasis"),
    ("D:Schisto", "Schistosomiasis"),
    ("D:LF", "Lymphatic filariasis"),
    ("D:Oncho", "Onchocerciasis"),
    ("D:STH", "Soil-transmitted helminthiasis"),
    ("D:Trachoma", "Trachoma"),
    ("D:Leprosy", "Leprosy"),
    ("D:Taeniasis", "Taeniasis/Cysticercosis"),
    ("D:Echino", "Echinococcosis"),
    ("D:Scabies", "Scabies"),
]
# Real MeSH descriptor IDs verified live via NLM MeSH (id.nlm.nih.gov/
# meshb.nlm.nih.gov) on 2026-08-29. Blank = not yet verified -- run
# src/enrich_from_sources.py's enrich_via_mesh() to fill the rest
# properly and fast via the live API rather than one-by-one search.
disease_mesh_ids = {
    "D:Chagas": "D014355",
    "D:HAT": "D014353",
}
for nid, name in diseases:
    add_node(nid, name, "Disease", mesh_id=disease_mesh_ids.get(nid, ""))

# Pathogens
pathogens = [
    ("P:Tcruzi", "Trypanosoma cruzi", "D:Chagas"),
    ("P:Tbrucei", "Trypanosoma brucei", "D:HAT"),
    ("P:Leishmania", "Leishmania donovani", "D:Leish"),
    ("P:Smansoni", "Schistosoma mansoni", "D:Schisto"),
    ("P:Wbancrofti", "Wuchereria bancrofti", "D:LF"),
    ("P:Ovolvulus", "Onchocerca volvulus", "D:Oncho"),
    ("P:Alumbricoides", "Ascaris lumbricoides", "D:STH"),
    ("P:Ctrachomatis", "Chlamydia trachomatis", "D:Trachoma"),
    ("P:Mleprae", "Mycobacterium leprae", "D:Leprosy"),
    ("P:Tsolium", "Taenia solium", "D:Taeniasis"),
    ("P:Egranulosus", "Echinococcus granulosus", "D:Echino"),
    ("P:Sscabiei", "Sarcoptes scabiei", "D:Scabies"),
]
for nid, name, _ in pathogens:
    add_node(nid, name, "Pathogen")

# Vectors
vectors = [
    ("V:Triatomine", "Triatomine bug", "P:Tcruzi"),
    ("V:Tsetse", "Tsetse fly", "P:Tbrucei"),
    ("V:Sandfly", "Sandfly", "P:Leishmania"),
    ("V:Snail", "Freshwater snail (Biomphalaria)", "P:Smansoni"),
    ("V:CulexMosquito", "Culex mosquito", "P:Wbancrofti"),
    ("V:Blackfly", "Blackfly (Simulium)", "P:Ovolvulus"),
]
for nid, name, _ in vectors:
    add_node(nid, name, "Vector")

# Pathogen-encoded drug targets (genes/proteins)
pathogen_genes = [
    ("PG:TryR", "Trypanothione reductase", ["P:Tcruzi", "P:Tbrucei", "P:Leishmania"]),
    ("PG:NTR", "Type I nitroreductase", ["P:Tcruzi", "P:Tbrucei"]),
    ("PG:ODC1_path", "Ornithine decarboxylase (pathogen)", ["P:Tbrucei"]),
    ("PG:PyrK_path", "Pyruvate kinase (pathogen)", ["P:Tbrucei"]),
    ("PG:TRPMPZQ", "TRPM-PZQ calcium channel", ["P:Smansoni"]),
    ("PG:GluCl", "Glutamate-gated chloride channel", ["P:Ovolvulus", "P:Wbancrofti", "P:Alumbricoides", "P:Sscabiei"]),
    ("PG:BetaTub", "Beta-tubulin (pathogen)", ["P:Wbancrofti", "P:Alumbricoides", "P:Tsolium"]),
    ("PG:DHFR_path", "Dihydrofolate reductase (pathogen)", ["P:Tbrucei"]),
    ("PG:kDNA_topo", "Kinetoplast DNA topoisomerase", ["P:Leishmania", "P:Tbrucei"]),
    ("PG:PhosphoMeta", "Phospholipid metabolism enzyme", ["P:Leishmania"]),
    ("PG:Ergosterol_path", "Ergosterol/sterol biosynthesis", ["P:Leishmania"]),
    ("PG:Ribosome30S", "Bacterial 30S ribosomal subunit", ["P:Ctrachomatis"]),
    ("PG:Ribosome50S", "Bacterial 50S ribosomal subunit", ["P:Ctrachomatis"]),
    ("PG:RNAPol_myco", "Mycobacterial RNA polymerase", ["P:Mleprae"]),
    ("PG:DHPS_myco", "Dihydropteroate synthase (mycobacterial)", ["P:Mleprae"]),
    ("PG:OxPhos_helminth", "Oxidative phosphorylation (helminth)", ["P:Tsolium"]),
    ("PG:ThioredoxinR", "Thioredoxin reductase (protozoal)", ["P:Leishmania"]),
]
for nid, name, _ in pathogen_genes:
    add_node(nid, name, "PathogenGene")

# Human ortholog / pathway-linked genes (small set; enables PGoG + GpPW reasoning)
human_genes = [
    ("G:ODC1", "ODC1 (human ornithine decarboxylase)"),
    ("G:TXNRD1", "TXNRD1 (human thioredoxin reductase 1)"),
    ("G:TUBB", "TUBB (human beta-tubulin)"),
    ("G:GABRA1", "GABRA1 (human GABA-A receptor, Cys-loop family)"),
    ("G:P2RX4", "P2RX4 (human purinergic receptor)"),
    ("G:AKT1", "AKT1 (human, phospholipid/Akt signaling)"),
    ("G:DHFR", "DHFR (human dihydrofolate reductase)"),
]
for nid, name in human_genes:
    add_node(nid, name, "Gene")

pathways = [
    ("PW:Polyamine", "Polyamine biosynthesis pathway"),
    ("PW:TrypanothioneMetab", "Trypanothione/redox metabolism"),
    ("PW:Microtubule", "Microtubule dynamics"),
    ("PW:CysLoopSignaling", "Cys-loop ion channel signaling"),
    ("PW:Purinergic", "Purinergic signaling"),
    ("PW:PI3K_Akt", "PI3K/Akt signaling"),
    ("PW:FolateMetab", "Folate metabolism"),
]
for nid, name in pathways:
    add_node(nid, name, "Pathway")

pharm_classes = [
    ("PC:Nitroheterocyclic", "Nitroheterocyclic antiparasitics"),
    ("PC:Benzimidazole", "Benzimidazole anthelmintics"),
    ("PC:Macrocyclic_lactone", "Macrocyclic lactone antiparasitics"),
    ("PC:Antileishmanial", "Antileishmanial agents"),
]
for nid, name in pharm_classes:
    add_node(nid, name, "PharmacologicClass")

# Side effects (schema now supports SIDER integration; a couple of
# well-known, real side effects seeded here so the CcSE edge type has
# at least illustrative coverage until enrich_side_effects_sider()
# pulls the full SIDER bulk file)
side_effects = [
    ("SE:Nausea", "Nausea"),
    ("SE:PeripheralNeuropathy", "Peripheral neuropathy"),
    ("SE:Hepatotoxicity", "Hepatotoxicity"),
    ("SE:QTprolongation", "QT interval prolongation"),
]
for nid, name in side_effects:
    add_node(nid, name, "SideEffect")

# Compounds (approved / repurposing candidates)
compounds = [
    ("C:Benznidazole", "Benznidazole"),
    ("C:Nifurtimox", "Nifurtimox"),
    ("C:Fexinidazole", "Fexinidazole"),
    ("C:Melarsoprol", "Melarsoprol"),
    ("C:Eflornithine", "Eflornithine (DFMO)"),
    ("C:Suramin", "Suramin"),
    ("C:Pentamidine", "Pentamidine"),
    ("C:Miltefosine", "Miltefosine"),
    ("C:AmphotericinB", "Amphotericin B"),
    ("C:Praziquantel", "Praziquantel"),
    ("C:Ivermectin", "Ivermectin"),
    ("C:DEC", "Diethylcarbamazine"),
    ("C:Albendazole", "Albendazole"),
    ("C:Mebendazole", "Mebendazole"),
    ("C:Doxycycline", "Doxycycline"),
    ("C:Azithromycin", "Azithromycin"),
    ("C:Rifampicin", "Rifampicin"),
    ("C:Dapsone", "Dapsone"),
    ("C:Clofazimine", "Clofazimine"),
    ("C:Niclosamide", "Niclosamide"),
    ("C:Auranofin", "Auranofin"),  # repurposing candidate, not yet NTD-indicated
    ("C:Fluoxetine", "Fluoxetine"),  # repurposing candidate
]
for nid, name in compounds:
    add_node(nid, name, "Compound")

with open(os.path.join(OUT_DIR, "nodes.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id", "name", "kind"] + XREF_COLUMNS)
    w.writeheader()
    w.writerows(nodes)

# ---------------------------------------------------------------- EDGES ---
edges = []


def add_edge(src, rel, dst, abbr, evidence="curated_seed"):
    edges.append({"source": src, "relation": rel, "target": dst, "abbr": abbr, "evidence": evidence})


# Pathogen -> Disease
for nid, name, disease in pathogens:
    add_edge(nid, "causes", disease, "PcD")

# Vector -> Pathogen
for nid, name, pathogen in vectors:
    add_edge(nid, "transmits", pathogen, "VtP")

# Pathogen -> PathogenGene
for nid, name, hosts in pathogen_genes:
    for h in hosts:
        add_edge(h, "hasGene", nid, "PhPG")

# PathogenGene -> ortholog Gene
ortholog_map = {
    "PG:ODC1_path": "G:ODC1",
    "PG:ThioredoxinR": "G:TXNRD1",
    "PG:BetaTub": "G:TUBB",
    "PG:GluCl": "G:GABRA1",
    "PG:DHFR_path": "G:DHFR",
}
for pg, g in ortholog_map.items():
    add_edge(pg, "orthologOf", g, "PGoG")

# PathogenGene -> Pathway / Gene -> Pathway
gene_pathway = [
    ("PG:ODC1_path", "PW:Polyamine"),
    ("G:ODC1", "PW:Polyamine"),
    ("PG:TryR", "PW:TrypanothioneMetab"),
    ("PG:kDNA_topo", "PW:TrypanothioneMetab"),
    ("PG:BetaTub", "PW:Microtubule"),
    ("G:TUBB", "PW:Microtubule"),
    ("PG:GluCl", "PW:CysLoopSignaling"),
    ("G:GABRA1", "PW:CysLoopSignaling"),
    ("G:P2RX4", "PW:Purinergic"),
    ("G:AKT1", "PW:PI3K_Akt"),
    ("PG:PhosphoMeta", "PW:PI3K_Akt"),
    ("PG:DHFR_path", "PW:FolateMetab"),
    ("G:DHFR", "PW:FolateMetab"),
    ("PG:ThioredoxinR", "PW:TrypanothioneMetab"),
]
for g, pw in gene_pathway:
    kind = "PGpPW" if g.startswith("PG:") else "GpPW"
    add_edge(g, "participates", pw, kind)

# PharmacologicClass -> Compound
pc_map = {
    "PC:Nitroheterocyclic": ["C:Benznidazole", "C:Nifurtimox", "C:Fexinidazole"],
    "PC:Benzimidazole": ["C:Albendazole", "C:Mebendazole"],
    "PC:Macrocyclic_lactone": ["C:Ivermectin"],
    "PC:Antileishmanial": ["C:Miltefosine", "C:AmphotericinB"],
}
for pc, comps in pc_map.items():
    for c in comps:
        add_edge(pc, "includes", c, "PCiC")

# Compound -> PathogenGene (CbPG) -- the key repurposing-relevant binding edges
compound_target = [
    ("C:Benznidazole", "PG:TryR"),
    ("C:Benznidazole", "PG:NTR"),
    ("C:Nifurtimox", "PG:NTR"),
    ("C:Nifurtimox", "PG:TryR"),
    ("C:Fexinidazole", "PG:NTR"),
    ("C:Melarsoprol", "PG:PyrK_path"),
    ("C:Melarsoprol", "PG:TryR"),
    ("C:Eflornithine", "PG:ODC1_path"),
    ("C:Suramin", "PG:DHFR_path"),
    ("C:Pentamidine", "PG:kDNA_topo"),
    ("C:Miltefosine", "PG:PhosphoMeta"),
    ("C:AmphotericinB", "PG:Ergosterol_path"),
    ("C:Praziquantel", "PG:TRPMPZQ"),
    ("C:Ivermectin", "PG:GluCl"),
    ("C:Albendazole", "PG:BetaTub"),
    ("C:Mebendazole", "PG:BetaTub"),
    ("C:Doxycycline", "PG:Ribosome30S"),
    ("C:Azithromycin", "PG:Ribosome50S"),
    ("C:Rifampicin", "PG:RNAPol_myco"),
    ("C:Dapsone", "PG:DHPS_myco"),
    ("C:Niclosamide", "PG:OxPhos_helminth"),
    # repurposing-candidate binding evidence (in vitro / other-indication
    # evidence, NOT yet an approved NTD indication -- these are exactly
    # the kind of edges a repurposing pipeline should surface as leads)
    ("C:Auranofin", "PG:ThioredoxinR"),
    ("C:Fluoxetine", "PG:kDNA_topo"),
]
for c, pg in compound_target:
    add_edge(c, "binds", pg, "CbPG", evidence="curated_seed_moa")

# Compound -> Gene (human targets, for repurposing-candidate compounds)
compound_human_target = [
    ("C:Auranofin", "G:TXNRD1"),
    ("C:Fluoxetine", "G:P2RX4"),  # simplified illustrative link
    ("C:Miltefosine", "G:AKT1"),
]
for c, g in compound_human_target:
    add_edge(c, "binds", g, "CbG")

# Compound -> Disease (TREATS, ground-truth approved indications)
# This is the TARGET edge type (CtD) the ML pipeline learns to predict.
compound_disease = [
    ("C:Benznidazole", "D:Chagas"),
    ("C:Nifurtimox", "D:Chagas"),
    ("C:Nifurtimox", "D:HAT"),          # NECT combination therapy
    ("C:Fexinidazole", "D:HAT"),
    ("C:Melarsoprol", "D:HAT"),
    ("C:Eflornithine", "D:HAT"),
    ("C:Suramin", "D:HAT"),
    ("C:Pentamidine", "D:HAT"),
    ("C:Pentamidine", "D:Leish"),
    ("C:Miltefosine", "D:Leish"),
    ("C:AmphotericinB", "D:Leish"),
    ("C:Praziquantel", "D:Schisto"),
    ("C:Praziquantel", "D:Taeniasis"),
    ("C:Ivermectin", "D:Oncho"),
    ("C:Ivermectin", "D:LF"),
    ("C:Ivermectin", "D:STH"),
    ("C:Ivermectin", "D:Scabies"),
    ("C:DEC", "D:LF"),
    ("C:Albendazole", "D:LF"),
    ("C:Albendazole", "D:STH"),
    ("C:Albendazole", "D:Echino"),
    ("C:Mebendazole", "D:STH"),
    ("C:Doxycycline", "D:Oncho"),       # anti-Wolbachia macrofilaricidal effect
    ("C:Doxycycline", "D:LF"),
    ("C:Azithromycin", "D:Trachoma"),
    ("C:Rifampicin", "D:Leprosy"),
    ("C:Dapsone", "D:Leprosy"),
    ("C:Clofazimine", "D:Leprosy"),
    ("C:Niclosamide", "D:Taeniasis"),
]
for c, d in compound_disease:
    add_edge(c, "treats", d, "CtD", evidence="approved_indication")

# Compound -> SideEffect (CcSE) -- a few well-documented examples;
# enrich_side_effects_sider() replaces/extends this with the full
# SIDER bulk mapping once network access is available.
compound_side_effect = [
    ("C:Benznidazole", "SE:PeripheralNeuropathy"),
    ("C:Nifurtimox", "SE:Nausea"),
    ("C:Melarsoprol", "SE:PeripheralNeuropathy"),
    ("C:AmphotericinB", "SE:Hepatotoxicity"),
    ("C:Pentamidine", "SE:QTprolongation"),
]
for c, se in compound_side_effect:
    add_edge(c, "causes", se, "CcSE", evidence="curated_seed_known_ade")

with open(os.path.join(OUT_DIR, "edges.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["source", "relation", "target", "abbr", "evidence"])
    w.writeheader()
    w.writerows(edges)

print(f"Wrote {len(nodes)} nodes -> nodes.csv")
print(f"Wrote {len(edges)} edges -> edges.csv")
