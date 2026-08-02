#!/usr/bin/env python3
"""Figure 7 (v7): STRING physical-interaction network of the 38-gene panel
(17 named candidates + connected context genes; non-tier isolated genes dropped)
+ g:Profiler pathway over-representation on the 17 named candidates.
v7b: smaller nodes, more spacing so edges are visible."""
import pandas as pd, numpy as np, networkx as nx, textwrap
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings("ignore")
ROOT="__MS_GEO_ROOT__"
OUT=ROOT+"/Poster_v2/figures/string_network_v6.png"   # overwrite the embedded Fig 7
T1=["ITGB2","CD79B","IKZF1"]
T2A=["LXN","SH3BP4","CASP6","CASP8","DGKQ","MX1","IFIT1","NUP210","RUNX3"]
T2N=["CTSZ","CHL1","ICAM1","FOXP3","ITGAL"]
NAMED=set(T1+T2A+T2N)
RED="#B71C1C"; ORANGE="#E65100"; PURPLE="#6A1B9A"; GREY="#9E9E9E"; HUBGREY="#777"
def tcol(g): return RED if g in set(T1) else (ORANGE if g in set(T2A) else (PURPLE if g in set(T2N) else GREY))
plt.rcParams.update({"font.family":"DejaVu Sans","axes.titleweight":"bold"})

# ---- 58-gene panel + 53 physical edges ----
hub=pd.read_csv(ROOT+"/PPI/ppi_hub_ranking.tsv",sep="\t")
gcol=[c for c in hub.columns if c.lower() in ("gene","name","node","symbol")][0]
allg=hub[gcol].astype(str).tolist()
ed=pd.read_csv(ROOT+"/PPI/ppi_string_physical_full.tsv",sep="\t")
G=nx.Graph(); G.add_nodes_from(allg)
for _,r in ed.iterrows(): G.add_edge(r["preferredName_A"],r["preferredName_B"],w=float(r["score"]))
conn=[g for g in allg if G.degree(g)>0]
iso_cand=[g for g in allg if G.degree(g)==0 and g in NAMED]
panel=conn+iso_cand
print(f"panel={len(panel)}  connected={len(conn)}  isolated-candidates={len(iso_cand)}: {iso_cand}")
Gc=G.subgraph(conn)

# ---- g:Profiler on the 17 named candidates ----
ENR_CACHE="/tmp/string17_gprofiler.csv"
try:
    from gprofiler import GProfiler
    gp=GProfiler(return_dataframe=True)
    enr=gp.profile(organism="hsapiens",query=list(NAMED),sources=["GO:BP","KEGG","REAC","WP"],
                   significance_threshold_method="g_SCS",user_threshold=0.05,no_evidences=True)
    enr=enr.sort_values("p_value").head(10)
    if len(enr): enr.to_csv(ENR_CACHE,index=False)
except Exception as e:
    print("gProfiler live failed:",e)
    import os
    enr=pd.read_csv(ENR_CACHE).sort_values("p_value").head(10) if os.path.exists(ENR_CACHE) else pd.DataFrame(columns=["name","p_value","source"])

fig=plt.figure(figsize=(17.5,10.0),dpi=200)
gs=gridspec.GridSpec(1,2,figure=fig,width_ratios=[1.65,1.0],wspace=0.30)

# ===== Panel A — network (smaller nodes, more spacing) =====
axA=fig.add_subplot(gs[0,0]); axA.axis("off")
deg=dict(G.degree())
SZ={n:380+195*deg.get(n,0) for n in panel}            # smaller nodes
pos=nx.spring_layout(Gc,k=4.3,iterations=1200,seed=7,weight="w")
P={n:np.asarray(pos[n],float) for n in conn}
_a=np.array([P[n] for n in conn]); rng=max((_a.max(0)-_a.min(0)).max(),1e-9)
P={n:(P[n]-_a.mean(0))/rng*5.6 for n in conn}         # much wider box
rad={n:0.0098*np.sqrt(SZ[n]) for n in conn}           # effective rendered radius (was under-estimated)
for _ in range(1500):                                  # de-overlap with a clear gap
    moved=False
    for i in range(len(conn)):
        for j in range(i+1,len(conn)):
            a,b=conn[i],conn[j]; dvec=P[a]-P[b]; dist=float(np.hypot(*dvec)); need=rad[a]+rad[b]+0.42
            if dist<need:
                if dist<1e-6: dvec=np.array([np.cos(i*1.3),np.sin(i*1.3)]); dist=1.0
                u=dvec/dist; sh=(need-dist)/2.0; P[a]=P[a]+u*sh; P[b]=P[b]-u*sh; moved=True
    if not moved: break
ymin=min(P[n][1] for n in conn)
xs=np.linspace(-5.0,5.0,len(iso_cand)); yiso=ymin-1.4
for k,g in enumerate(iso_cand): P[g]=np.array([xs[k],yiso])
for (u,v,dd) in G.edges(data=True):
    axA.plot([P[u][0],P[v][0]],[P[u][1],P[v][1]],color="#9AA0A6",lw=0.9+3.8*dd["w"],alpha=0.78,zorder=1)
for n in panel:
    axA.scatter(*P[n],s=SZ[n],color=tcol(n),alpha=0.95,edgecolors="white",lw=1.3,zorder=3)
    axA.text(P[n][0],P[n][1],n,fontsize=6.4,fontstyle="italic",fontweight="bold",ha="center",va="center",color="white",zorder=4)
axA.plot([-5.4,5.4],[yiso+0.6,yiso+0.6],color="#DDD",lw=1,zorder=0)
axA.text(0,yiso+0.63,"named candidates with no high-confidence physical edge within the panel",
         fontsize=9,ha="center",va="bottom",color="#777",style="italic")
axA.set_aspect("equal")
xs_all=[P[n][0] for n in panel]; ys_all=[P[n][1] for n in panel]
axA.set_xlim(min(xs_all)-0.6,max(xs_all)+0.6); axA.set_ylim(min(ys_all)-0.45,max(ys_all)+0.6)
axA.set_title("A",fontsize=27,fontweight="bold",loc="left",pad=6)
axA.legend(handles=[Patch(fc=RED,label="Inverse-concordant Tier-1"),Patch(fc=ORANGE,label="Tier-2 auxiliary"),
                    Patch(fc=PURPLE,label="Tier-2 non-concordant anchor"),Patch(fc=GREY,label="canonical MS marker (context)")],
           fontsize=11.5,loc="upper center",bbox_to_anchor=(0.5,-0.005),ncol=2,frameon=False,handletextpad=0.4,columnspacing=1.4)
hubs=sorted(((g,deg[g]) for g in conn),key=lambda z:-z[1])[:4]
axA.text(0.99,0.99,"hubs: "+", ".join(f"{h} ({d})" for h,d in hubs),transform=axA.transAxes,
         ha="right",va="top",fontsize=12,color=HUBGREY)

# ===== Panel B — g:Profiler =====
axB=fig.add_subplot(gs[0,1])
if len(enr):
    e=enr.iloc[::-1]; y=np.arange(len(e))
    SRCCOL={"GO:BP":"#1565C0","KEGG":"#2E7D32","REAC":"#00838F","WP":"#8E24AA"}
    axB.barh(y,-np.log10(e["p_value"]),color=[SRCCOL.get(s,"#777") for s in e["source"]],alpha=0.88,height=0.66)
    axB.set_yticks(y); axB.set_yticklabels(["\n".join(textwrap.wrap(n,30)) for n in e["name"]],fontsize=10.6)
    axB.set_xlabel("−log₁₀ (g:SCS-corrected P)",fontsize=14.5)
    axB.legend(handles=[Patch(fc=c,label=s) for s,c in SRCCOL.items()],fontsize=12,loc="lower right",frameon=False,title="source",title_fontsize=18)
else:
    axB.text(0.5,0.5,"no enrichment",ha="center")
axB.set_title("B",fontsize=27,fontweight="bold",loc="left",pad=6)
axB.spines[["top","right"]].set_visible(False)

plt.savefig(OUT,dpi=200,bbox_inches="tight",facecolor="white")
from PIL import Image
print("✓ saved",OUT,Image.open(OUT).size,"| panel genes:",len(panel),"| edges:",G.number_of_edges())
