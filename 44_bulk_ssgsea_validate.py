#!/usr/bin/env python3
"""
Runs ssGSEA on ileum bulk for the regulon GMTs, then makes UC vs CD boxplots.
Accepts a full bulk matrix; filters ileum using metadata file.
"""
import os, argparse, pandas as pd, numpy as np, matplotlib.pyplot as plt
from gseapy import ssgsea
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu

def read_table(path):
    df = pd.read_csv(path, sep=None, engine="python", index_col=0)
    return df

def orient_genes_x_samples(df):
    # If rows < cols and columns look like genes, transpose
    if df.shape[0] < df.shape[1]:
        # try: are columns gene-like (uppercase letters)?
        upp = sum([str(c).isupper() for c in df.columns[:50]])
        if upp > 20:
            df = df.T
            print("[INFO] Transposed bulk to genes x samples")
    return df

def filter_ileum(df, meta, id_col="sample_id"):
    dat = df.copy()
    # meta must have rows indexed by sample_id
    meta = meta.copy()
    if id_col not in meta.columns: raise ValueError("meta must have 'sample_id'")
    meta = meta.set_index(id_col)
    shared = dat.columns.intersection(meta.index)
    dat = dat.loc[:, shared]
    meta = meta.loc[shared]
    ile = meta["tissue"].str.contains("ileum", case=False, na=False)
    return dat.loc[:, ile], meta.loc[ile]

def run_ssgsea(expr_gxs, gmt_path, outdir, tag):
    os.makedirs(outdir, exist_ok=True)
    # gseapy needs genes x samples with gene names in index
    pre_res = ssgsea(data=expr_gxs, gene_sets=gmt_path, sample_norm_method="rank", outdir=None, scale=False, processes=4)
    # pre_res results is a DataFrame pathways x samples
    scores = pre_res.res2d
    scores.to_csv(os.path.join(outdir, f"{tag}_ssgsea_scores.tsv"), sep="\t")
    return scores

def boxplots(scores, meta, outdir, tag):
    os.makedirs(outdir, exist_ok=True)
    results=[]
    for pathway in scores.index:
        sc = scores.loc[pathway]
        df = pd.DataFrame({"score": sc, "group": meta["disease_group"]})
        # UC vs CD
        uc = df.loc[df["group"]=="UC","score"]
        cd = df.loc[df["group"]=="CD","score"]
        if len(uc) >= 3 and len(cd) >= 3:
            stat, p = mannwhitneyu(uc, cd, alternative="two-sided")
            results.append((pathway, stat, p, len(uc), len(cd)))
        # plot
        plt.figure(figsize=(3.6,3), dpi=180)
        for i,g in enumerate(["HC","UC","CD"]):
            y = df.loc[df["group"]==g,"score"].dropna()
            x = np.random.normal(i, 0.06, size=len(y))
            plt.plot(x, y, "o", ms=2, alpha=0.35, color="black")
            plt.scatter([i], [y.mean() if len(y)>0 else np.nan], s=55, color="#ff7f0e")
        plt.xticks([0,1,2], ["HC","UC","CD"]); plt.ylabel("ssGSEA NES")
        plt.title(f"{pathway} (ileum)")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{tag}_{pathway}_box.png")); plt.close()
    if results:
        T = pd.DataFrame(results, columns=["pathway","U","p","n_UC","n_CD"]).set_index("pathway")
        T["q"] = multipletests(T["p"], method="fdr_bh")[1]
        T.to_csv(os.path.join(outdir, f"{tag}_UC_vs_CD.tsv"), sep="\t")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bulk", required=True, help="bulk expression TSV (genes x samples or samples x genes)")
    ap.add_argument("--meta", required=True, help="bulk metadata CSV with sample_id,tissue,disease_group")
    ap.add_argument("--gmt_dc2", required=True, help="GMT file for DC2 regulons")
    ap.add_argument("--gmt_dp", required=True, help="GMT file for DP regulons")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    bulk = read_table(args.bulk)
    bulk = orient_genes_x_samples(bulk)
    meta = pd.read_csv(args.meta)

    ile_expr, ile_meta = filter_ileum(bulk, meta)
    print("[INFO] ileum samples:", ile_expr.shape[1])

    # DC2 regulons
    sc_dc2 = run_ssgsea(ile_expr, args.gmt_dc2, os.path.join(args.outdir,"DC2"), "DC2")
    boxplots(sc_dc2, ile_meta, os.path.join(args.outdir,"DC2"), "DC2")

    # DP regulons
    sc_dp  = run_ssgsea(ile_expr, args.gmt_dp,  os.path.join(args.outdir,"DP"),  "DP")
    boxplots(sc_dp,  ile_meta, os.path.join(args.outdir,"DP"),  "DP")

if __name__ == "__main__":
    main()
