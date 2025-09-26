#!/usr/bin/env python3
import scanpy as sc, pandas as pd, numpy as np, os, matplotlib.pyplot as plt

IN=[("DC2_ileum.h5ad","DC2"),("DP_ileum.h5ad","DP")]
GMT="resources/metabolic_9sets.gmt"
OUT="METAB_ENRICH"; os.makedirs(OUT, exist_ok=True)

def read_gmt(path):
    sets={}
    with open(path) as f:
        for line in f:
            parts=line.strip().split("\t")
            if len(parts)>=3:
                sets[parts[0]]= [g for g in parts[2:] if g]
    return sets

def score_sets(ad, sets):
    for name, genes in sets.items():
        g = [g for g in genes if g in ad.var_names]
        if len(g) < 5: 
            print(f"[WARN] {name}: only {len(g)} genes present; skipping")
            continue
        sc.tl.score_genes(ad, gene_list=g, score_name=f"GS_{name}", use_raw=False)
    return ad

def summarize(ad, tag, sets):
    cols=[f"GS_{k}" for k in sets.keys() if f"GS_{k}" in ad.obs.columns]
    if not cols: return
    key_state = "state_label" if "state_label" in ad.obs else "leiden_sub"
    key_group = "disease_group" if "disease_group" in ad.obs else None

    # mean by state
    m = ad.obs.groupby(key_state)[cols].mean().T
    m.to_csv(f"{OUT}/{tag}_ileum_state_means.tsv", sep="\t")

    # heatmap
    plt.figure(figsize=(0.6*m.shape[1]+3, 0.4*m.shape[0]+2), dpi=200)
    plt.imshow(m, aspect="auto")
    plt.yticks(range(m.shape[0]), m.index, fontsize=7)
    plt.xticks(range(m.shape[1]), m.columns, rotation=45, ha="right", fontsize=7)
    plt.colorbar(label="Module score"); plt.title(f"{tag} metabolic pathways (ileum)")
    plt.tight_layout(); plt.savefig(f"{OUT}/FIG_{tag}_metabolic_heatmap.png"); plt.close()

    # UC vs CD per pathway (pooled across states) — quick box
    if key_group:
        for name in [c.replace("GS_","") for c in cols]:
            col=f"GS_{name}"
            d = ad.obs[[col, key_group]].dropna()
            plt.figure(figsize=(3.2,3), dpi=180)
            for i,g in enumerate(["HC","UC","CD"]):
                y = d.loc[d[key_group]==g, col]
                x = np.random.normal(i, 0.06, size=len(y))
                plt.plot(x, y, "o", ms=2, alpha=0.35, color="black")
                plt.scatter([i], [y.mean() if len(y)>0 else np.nan], s=55, color="#ff7f0e")
            plt.xticks([0,1,2], ["HC","UC","CD"]); plt.ylabel(name)
            plt.title(f"{tag} ileum"); plt.tight_layout()
            plt.savefig(f"{OUT}/FIG_{tag}_{name}_box.png"); plt.close()

def main():
    sets = read_gmt(GMT)
    for h5, tag in IN:
        ad = sc.read(h5)
        ad = score_sets(ad, sets)
        ad.write(f"{tag}_ileum_with_gscores.h5ad")
        summarize(ad, tag, sets)
        print("[OK]", tag, "metabolic scoring done.")

if __name__ == "__main__":
    main()
