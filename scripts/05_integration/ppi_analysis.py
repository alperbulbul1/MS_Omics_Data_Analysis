#!/usr/bin/env python3
"""ppi_analysis.py — protein-protein interaction analysis for MS_GEO candidates.

Queries STRING (REST API, physical-only) and GeneMANIA (physical interaction
attributes) for the union of our high-confidence candidates:
  - CO7 panel (LXN, SH3BP4, CHL1, CTSZ, RPAP2, PCNP, THRB)
  - Top 30 inverse-concordant genes (nb09)
  - Multi-stratum validated TX (LXN, RPAP2, STAT1, STAT3, TYK2)
  - scRNA brain/CSF/blood validated (nb13)

Outputs:
  PPI/
    ppi_string_physical_full.tsv     full STRING physical edge list
    ppi_string_partners_top10.tsv    each query gene's top-10 partners
    ppi_string_enrichment.tsv        GO/KEGG/Reactome enrichment
    ppi_genemania_physical.tsv       GeneMANIA physical interactions
    ppi_hub_ranking.tsv              degree-ranked hub genes
    ppi_subnetwork_*.png             3-4 community sub-networks
    ppi_network_full.png             full network figure
"""
from pathlib import Path
import time, json
import warnings; warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import requests
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

PROJ = Path("__MS_GEO_ROOT__")
ME   = PROJ / "Methylation" / "results"
OUT  = PROJ / "PPI"; OUT.mkdir(exist_ok=True)
SPECIES = "9606"  # H. sapiens

# Build candidate gene panel
CO7 = ["LXN","SH3BP4","CHL1","CTSZ","RPAP2","PCNP","THRB"]
inv = pd.read_csv(ME / "INVERSE_CONCORDANT_by_gene.tsv", sep="\t")
INV_TOP = inv.sort_values(["n_pairings","best_rna_fdr"], ascending=[False, True]).gene.tolist()[:30]
MULTI_TX = pd.read_csv(PROJ / "Transcriptome" / "results" / "MultiStratum_Validated_R.tsv", sep="\t").gene.tolist()
master = pd.read_csv(ME / "MASTER_4layer_validation.tsv", sep="\t")
# scRNA validated (>=1 layer scRNA sig)
SCRNA_VALID = master[master["scRNA_n_sig"] > 0].gene.tolist()
# Additional curated MS-relevant
ADDITIONAL = ["STAT1","STAT3","TYK2","JCHAIN","SDC1","CHIT1","IGHM","CHI3L1","MBP","GFAP","NEFL",
              "CXCL13","IL17A","IFNG","FOXP3","IL2RA","CD4","CD8A","CD19","ITGAL","ITGAM","ICAM1"]

PANEL = list(dict.fromkeys(CO7 + INV_TOP + MULTI_TX + SCRNA_VALID + ADDITIONAL))
print(f"Query panel: {len(PANEL)} genes")
print("  ", ", ".join(PANEL[:30]), "..." if len(PANEL) > 30 else "")

# =====================================================================
# 1. STRING physical network
# =====================================================================
STRING_BASE = "https://string-db.org/api"
CALLER = "MSGEO_pipeline"

def string_network(genes, network_type="physical"):
    """Query STRING network endpoint. Returns DataFrame of edges."""
    ids = "%0d".join(genes)
    url = (f"{STRING_BASE}/tsv/network"
           f"?identifiers={ids}&species={SPECIES}"
           f"&network_type={network_type}&caller_identity={CALLER}"
           f"&required_score=0")  # we'll filter ourselves
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    from io import StringIO
    return pd.read_csv(StringIO(r.text), sep="\t")

def string_partners(gene, network_type="physical", limit=10):
    """Get top-N partners for one gene."""
    url = (f"{STRING_BASE}/tsv/interaction_partners"
           f"?identifiers={gene}&species={SPECIES}"
           f"&network_type={network_type}&limit={limit}"
           f"&caller_identity={CALLER}")
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    from io import StringIO
    return pd.read_csv(StringIO(r.text), sep="\t")

def string_enrichment(genes):
    """GO/KEGG/Reactome enrichment via STRING."""
    ids = "%0d".join(genes)
    url = (f"{STRING_BASE}/tsv/enrichment"
           f"?identifiers={ids}&species={SPECIES}"
           f"&caller_identity={CALLER}")
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    from io import StringIO
    return pd.read_csv(StringIO(r.text), sep="\t")

print("\n[1/4] STRING physical network (network within panel)")
try:
    net = string_network(PANEL, "physical")
    print(f"  Got {len(net)} edges")
    net.to_csv(OUT / "ppi_string_physical_full.tsv", sep="\t", index=False)
except Exception as e:
    print(f"  ERR: {e}")
    net = pd.DataFrame()

print("\n[2/4] STRING per-gene top-10 physical partners")
partner_rows = []
for i, g in enumerate(PANEL):
    try:
        d = string_partners(g, "physical", limit=10)
        if not d.empty:
            d["query_gene"] = g
            partner_rows.append(d)
        if i % 10 == 0: print(f"  ...{i}/{len(PANEL)} ({g})")
        time.sleep(0.15)  # gentle rate limit
    except Exception as e:
        print(f"  {g}: {e}")
partners = pd.concat(partner_rows, ignore_index=True) if partner_rows else pd.DataFrame()
partners.to_csv(OUT / "ppi_string_partners_top10.tsv", sep="\t", index=False)
print(f"  Total: {len(partners)} partner edges across {partners.query_gene.nunique() if not partners.empty else 0} query genes")

print("\n[3/4] STRING enrichment (functional + pathway)")
try:
    enr = string_enrichment(PANEL)
    enr_sig = enr[enr["fdr"] < 0.05].copy()
    print(f"  Got {len(enr)} terms total ({len(enr_sig)} FDR<0.05)")
    enr.to_csv(OUT / "ppi_string_enrichment_all.tsv", sep="\t", index=False)
    enr_sig.to_csv(OUT / "ppi_string_enrichment_sig.tsv", sep="\t", index=False)
except Exception as e:
    print(f"  ERR: {e}")
    enr = pd.DataFrame()

# =====================================================================
# 2. GeneMANIA physical interactions
# =====================================================================
print("\n[4/4] GeneMANIA physical interactions")
# GeneMANIA REST API: takes gene list, returns network with attributes
GM_BASE = "https://genemania.org/json/network_data"
def genemania_physical(genes):
    """Query GeneMANIA for physical interactions among input genes."""
    params = {
        "organism": "9606",
        "genes": "|".join(genes),
        "weighting": "average",
        "network_categories": "Physical_Interactions",
        "geneThreshold": 0,
        "attrThreshold": 0,
    }
    try:
        r = requests.post("https://genemania.org/json/network_data",
                          data=params, timeout=90)
        if r.status_code == 200:
            return r.json()
    except Exception as e:
        print(f"  GeneMANIA error: {e}")
    return None

# GeneMANIA can be slow; try with a subset
gm_data = genemania_physical(PANEL[:80])  # GeneMANIA recommends <=100 genes
gm_edges = []
if gm_data and "interactions" in gm_data:
    for it in gm_data["interactions"]:
        for e in it.get("elements", []):
            gm_edges.append(dict(
                gene_a=e["start"], gene_b=e["end"],
                weight=e.get("weight", 0.0),
                network=it.get("name", ""), category="Physical_Interactions"))
gm_df = pd.DataFrame(gm_edges)
gm_df.to_csv(OUT / "ppi_genemania_physical.tsv", sep="\t", index=False)
print(f"  GeneMANIA: {len(gm_df)} physical interaction edges")

# =====================================================================
# 3. Build networkx graph + hub ranking
# =====================================================================
print("\n[5/x] Building unified physical PPI graph")
G = nx.Graph()
G.add_nodes_from(PANEL)
# Add STRING physical edges (combined score >= 0.4)
if not net.empty:
    for _, row in net.iterrows():
        a, b = row["preferredName_A"], row["preferredName_B"]
        sc = float(row["score"]) if "score" in row else 0
        ec = float(row.get("escore", 0))  # experimental physical evidence
        dc = float(row.get("dscore", 0))  # database
        if sc < 0.4: continue
        if G.has_edge(a, b):
            G[a][b]["string_score"] = max(G[a][b]["string_score"], sc)
        else:
            G.add_edge(a, b, string_score=sc, escore=ec, dscore=dc, source="STRING")
# Add GeneMANIA edges
for _, row in gm_df.iterrows():
    a, b = row["gene_a"], row["gene_b"]
    if not G.has_node(a): G.add_node(a)
    if not G.has_node(b): G.add_node(b)
    if G.has_edge(a, b):
        G[a][b]["genemania_weight"] = row["weight"]
        G[a][b]["source"] = G[a][b].get("source", "") + ";GeneMANIA"
    else:
        G.add_edge(a, b, genemania_weight=row["weight"], source="GeneMANIA")

# Limit to nodes IN our panel + their direct partners
panel_set = set(PANEL)
in_panel_edges = [(u, v) for u, v in G.edges() if u in panel_set and v in panel_set]
G_panel = G.edge_subgraph(in_panel_edges).copy()
print(f"  Full graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
print(f"  Within-panel: {G_panel.number_of_nodes()} nodes, {G_panel.number_of_edges()} edges")

# Hub ranking (degree)
hub_rows = []
for n in G.nodes():
    deg = G.degree(n)
    deg_in_panel = G_panel.degree(n) if n in G_panel else 0
    hub_rows.append(dict(
        gene=n, degree_total=deg, degree_within_panel=deg_in_panel,
        in_CO7=(n in CO7),
        in_INV_TOP=(n in INV_TOP),
        in_multi_TX=(n in MULTI_TX),
        in_scRNA_valid=(n in SCRNA_VALID)))
hub_df = pd.DataFrame(hub_rows).sort_values(["degree_within_panel","degree_total"],
                                              ascending=[False, False])
hub_df.to_csv(OUT / "ppi_hub_ranking.tsv", sep="\t", index=False)
print(f"  Hub ranking saved. Top 15:")
print(hub_df.head(15).to_string(index=False))

# =====================================================================
# 4. Visualization
# =====================================================================
print("\n[6/x] Network visualization")
def color_for(node):
    if node in CO7:           return "#D62828"   # CO7 = red
    if node in INV_TOP[:15]:  return "#1F4E79"   # top inverse = dark blue
    if node in INV_TOP:       return "#3E92CC"   # inverse rest = light blue
    if node in MULTI_TX:      return "#7B3FA0"   # multi-stratum TX
    if node in SCRNA_VALID:   return "#0F8B5C"   # scRNA validated
    return "#888888"

def size_for(node, scale=80):
    return scale + 30 * G.degree(node)

# Full panel sub-network
fig, ax = plt.subplots(figsize=(16, 14), dpi=150)
pos = nx.spring_layout(G_panel, k=0.6, iterations=80, seed=42)
edges = G_panel.edges(data=True)
edge_colors = ["#2A9D8F" if "STRING" in e.get("source","") and "GeneMANIA" in e.get("source","")
                else "#1F4E79" if "STRING" in e.get("source","")
                else "#F4A261" for _, _, e in edges]
edge_widths = [1.0 + 4*e.get("string_score", 0.5) for _, _, e in edges]
nx.draw_networkx_edges(G_panel, pos, edge_color=edge_colors, width=edge_widths,
                        alpha=0.45, ax=ax)
node_colors = [color_for(n) for n in G_panel.nodes()]
node_sizes  = [size_for(n) for n in G_panel.nodes()]
nx.draw_networkx_nodes(G_panel, pos, node_color=node_colors, node_size=node_sizes,
                        alpha=0.85, linewidths=0.7, edgecolors="black", ax=ax)
nx.draw_networkx_labels(G_panel, pos, font_size=8, font_weight="bold", ax=ax)

legend = [
    mpatches.Patch(color="#D62828", label="CO7 panel (n=7)"),
    mpatches.Patch(color="#1F4E79", label="Inverse-concordant top 15"),
    mpatches.Patch(color="#3E92CC", label="Inverse-concordant 16-30"),
    mpatches.Patch(color="#7B3FA0", label="Multi-stratum TX validated"),
    mpatches.Patch(color="#0F8B5C", label="scRNA-validated (any tissue)"),
    mpatches.Patch(color="#888888", label="Other / curated MS"),
]
ax.legend(handles=legend, loc="upper left", fontsize=9, framealpha=0.92)
ax.set_title(f"Physical PPI network — MS_GEO candidates within-panel "
              f"({G_panel.number_of_nodes()} nodes, {G_panel.number_of_edges()} edges)\n"
              f"STRING (physical) + GeneMANIA; edge width = STRING combined score",
              fontsize=12, fontweight="bold")
ax.set_axis_off()
plt.tight_layout()
plt.savefig(OUT / "ppi_network_full.png", dpi=150, bbox_inches="tight", facecolor="white")
print(f"  Wrote {OUT / 'ppi_network_full.png'}")

# Community detection (Louvain) - sub-networks
try:
    import community as community_louvain  # python-louvain
    parts = community_louvain.best_partition(G_panel, random_state=42)
except ImportError:
    # Fallback: greedy modularity (built into networkx)
    comm_sets = list(nx.community.greedy_modularity_communities(G_panel))
    parts = {n: i for i, com in enumerate(comm_sets) for n in com}

print(f"  Detected {len(set(parts.values()))} communities")
comm_df = pd.DataFrame([dict(gene=n, community=c) for n, c in parts.items()])
comm_df = comm_df.merge(hub_df[["gene","degree_within_panel","in_CO7","in_INV_TOP","in_scRNA_valid"]],
                        on="gene")
comm_df = comm_df.sort_values(["community","degree_within_panel"], ascending=[True, False])
comm_df.to_csv(OUT / "ppi_communities.tsv", sep="\t", index=False)
print(f"  Wrote {OUT / 'ppi_communities.tsv'}")

# Render per-community panel
unique_comms = sorted(set(parts.values()))[:6]  # top 6 communities
ncols = 3; nrows = int(np.ceil(len(unique_comms) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows), dpi=150)
axes = np.atleast_1d(axes).flatten()
for i, c in enumerate(unique_comms):
    ax = axes[i]
    nodes_c = [n for n in G_panel.nodes() if parts.get(n) == c]
    Gc = G_panel.subgraph(nodes_c)
    if Gc.number_of_nodes() == 0:
        ax.set_axis_off(); continue
    pos_c = nx.spring_layout(Gc, seed=42, k=0.5)
    nx.draw_networkx_edges(Gc, pos_c, alpha=0.4, edge_color="#888", ax=ax)
    nx.draw_networkx_nodes(Gc, pos_c,
                            node_color=[color_for(n) for n in Gc.nodes()],
                            node_size=[size_for(n, 60) for n in Gc.nodes()],
                            alpha=0.85, linewidths=0.5, edgecolors="black", ax=ax)
    nx.draw_networkx_labels(Gc, pos_c, font_size=7, font_weight="bold", ax=ax)
    ax.set_title(f"Community {c+1}  ({Gc.number_of_nodes()} nodes, "
                  f"{Gc.number_of_edges()} edges)",
                  fontsize=10, fontweight="bold")
    ax.set_axis_off()
for j in range(i+1, len(axes)): axes[j].set_axis_off()
plt.suptitle("Physical PPI sub-networks (Louvain communities)",
              fontsize=13, fontweight="bold", y=1.01)
plt.tight_layout()
plt.savefig(OUT / "ppi_subnetworks.png", dpi=150, bbox_inches="tight", facecolor="white")
print(f"  Wrote {OUT / 'ppi_subnetworks.png'}")

print(f"\n✓ Done. Outputs in {OUT}")
for f in sorted(OUT.iterdir()):
    print(f"  {f.name}  ({f.stat().st_size/1024:.1f} KB)")
