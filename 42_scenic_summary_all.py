#!/usr/bin/env python3
import scanpy as sc, pandas as pd, numpy as np, os, matplotlib.pyplot as plt

SETS=[("DC2_subset_reclustered_fullX_states.h5ad","DC2"),
      ("DP_subset_reclustered_fullX_states.h5ad","DP")]
OUT="SCENIC_SUMMARY_ALL"; os.makedirs(OUT, exist_ok=True)

def attach_auc(ad, auc_csv):
    auc = pd.read_csv(auc_csv, index_col=0)   # regulons x cells
    auc = auc.loc[:, auc.columns.intersection(ad.obs_names)]
    auc = auc.loc[:, ad.obs_names]            # reorder to cells
    ad.obsm["AUCell"] = auc.T
    return ad

def heatmap_states(ad, tag):
    if "state_label" not in ad.obs: return
    M = ad.obsm["AUCell"].groupby(ad.obs["state_label"]).mean()
    # top 40 variable regulons
    var = M.var(axis=1).sort_values(ascending=False).head(40).index
    plt.figure(figsize=(12,8), dpi=200)
    plt.imshow(M.loc[var], aspect="auto")
    plt.yticks(range(len(var)), var, fontsize=6)
    plt.xticks(range(M.shape[1]), M.columns, rotation=45, ha="right", fontsize=8)
    plt.colorbar(label="AUCell")
    plt.title(f"{tag}: top regulons by state (all tissues)")
    plt.tight_layout()
    plt.savefig(f"{OUT}/FIG_{tag}_states_TF_heatmap.png"); plt.close()
    M.to_csv(f"{OUT}/{tag}_mean_AUCell_by_state.csv")

def box_by_group(ad, tag, tfs):
    if "disease_group" not in ad.obs: return
    for tf in tfs:
        if tf not in ad.obsm["AUCell"].columns: continue
        s = ad.obsm["AUCell"][tf]
        d = pd.DataFrame({"AUCell": s.values, "group": ad.obs["disease_group"].astype(str).values})
        plt.figure(figsize=(3.6,3), dpi=180)
        for i, g in enumerate(["HC","UC","CD"]):
            y = d.loc[d["group"]==g, "AUCell"].dropna()
            x = np.random.normal(i, 0.06, size=len(y))
            plt.plot(x, y, "o", ms=2, alpha=0.35, color="black")
            plt.scatter([i], [y.mean() if len(y)>0 else np.nan], s=55, color="#ff7f0e")
        plt.xticks([0,1,2], ["HC","UC","CD"])
        plt.ylabel(f"{tf} AUCell")
        plt.title(f"{tag} (all tissues)")
        plt.tight_layout()
        plt.savefig(f"{OUT}/FIG_{tag}_{tf}_AUCell_groups.png"); plt.close()

def main():
    EX_TFS = {"DC2": ["IRF4","RELB","KLF4","NFKB1","STAT1"],
              "DP":  ["AHR","TCF7","LEF1","STAT1","IRF7"]}
    for h5, tag in SETS:
        ad = sc.read(h5)
        ad = attach_auc(ad, f"{tag}_all.auc.csv")
        heatmap_states(ad, tag)
        box_by_group(ad, tag, EX_TFS[tag])
        print("[OK]", tag, "summaries in", OUT)

if __name__ == "__main__":
    main()
