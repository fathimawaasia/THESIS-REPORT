#!/usr/bin/env python3
# Make tidy UMAPs for CellTypist outputs (no re-annotation)

import os
import argparse
import scanpy as sc
import matplotlib.pyplot as plt

def parse_args():
    p = argparse.ArgumentParser(description="Plot CellTypist UMAPs from an annotated H5AD.")
    p.add_argument("--adata", default="adata_celltypist_qc.h5ad",
                   help="Annotated H5AD (already contains celltypist_* columns).")
    p.add_argument("--label_col", default="ct_seed",
                   help="Seed label column to plot (default: ct_seed).")
    p.add_argument("--figdir", default="celltypist_qc/figs",
                   help="Directory to write figures.")
    p.add_argument("--pt", type=float, default=0.6,
                   help="UMAP point size (default: 0.6).")
    p.add_argument("--topk", type=int, default=20,
                   help="Number of top groups to label on-map for label_col.")
    return p.parse_args()

def ensure_neighbors_umap(adata):
    """Ensure neighbors/UMAP exist; prefer X_scVI if available."""
    use_rep = "X_scVI" if "X_scVI" in adata.obsm else None
    need = (
        "neighbors" not in adata.uns or
        (use_rep and adata.uns["neighbors"].get("params", {}).get("use_rep") != use_rep) or
        "X_umap" not in adata.obsm
    )
    if need:
        if use_rep:
            sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep)
        else:
            sc.pp.pca(adata, n_comps=50, svd_solver="arpack")
            sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.umap(adata)

def umap_with_side_legend(ad, key, fname, size=1.0):
    if key in ad.obs:
        sc.pl.umap(
            ad,
            color=[key],
            legend_loc="right margin",   # side legend, no text on points
            frameon=False,
            size=size,
            show=False,
            save=f"_{fname}.png",
        )

def umap_ondata_topk(ad, key, k, fname, size=1.0):
    if key in ad.obs:
        vc = ad.obs[key].value_counts()
        topk = list(vc.head(k).index.astype(str))
        sc.pl.umap(
            ad,
            color=[key],
            groups=topk,                 # only top-K get colors/labels
            legend_loc="on data",
            legend_fontsize=6,
            frameon=False,
            size=size,
            show=False,
            save=f"_{fname}.png",
        )

def main():
    args = parse_args()
    os.makedirs(args.figdir, exist_ok=True)
    sc.set_figure_params(dpi=180, figsize=(6, 5))
    sc.settings.figdir = args.figdir

    print(f"[LOAD] {args.adata}")
    ad = sc.read(args.adata)

    ensure_neighbors_umap(ad)

    # 1) Side-legend (all categories), no on-data text
    umap_with_side_legend(ad, "celltypist_majority", "umap_celltypist_majority_legend", size=args.pt)
    umap_with_side_legend(ad, args.label_col,        "umap_ct_seed_legend",              size=args.pt)

    # 2) On-data labels, only top-K of label_col
    umap_ondata_topk(ad, args.label_col, args.topk, "umap_ct_seed_topK_ondata", size=args.pt)

    # 3) Leiden clusters with on-data labels (if present)
    if "leiden" in ad.obs:
        sc.pl.umap(
            ad,
            color=["leiden"],
            legend_loc="on data",
            legend_fontsize=7,
            frameon=False,
            size=args.pt,
            show=False,
            save="_umap_leiden_ondata.png",
        )

    print(f"[DONE] Figures in: {args.figdir}")

if __name__ == "__main__":
    main()
