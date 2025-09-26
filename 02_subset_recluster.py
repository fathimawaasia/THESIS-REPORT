#!/usr/bin/env python3
"""
02_subset_recluster_raw_simple.py

Subsets scANVI object (with .raw full genes) to DC2 and DP,
reclusters using scANVI latent, ranks markers from .raw,
and writes both HVG and full-gene variants.

Inputs:
  - H5AD_PATH (must have .raw with full gene set)

Outputs (per subset):
  - <TAG>_subset_reclustered.h5ad           # .X = HVG; .raw = full (unchanged)
  - <TAG>_subset_reclustered_fullX.h5ad     # .X = full genes copied from .raw
"""

import scanpy as sc
import pandas as pd
import numpy as np

# ---------- config ----------
H5AD_PATH  = "adata_scanvi_withRAW_plusTD.h5ad"
LABEL_COL  = "cell_type_final"
LEIDEN_RES = 0.4
USE_REP    = "X_scANVI"
WRITE_FULLX_VARIANT = True   # also write *_fullX.h5ad with .X replaced by .raw
# ----------------------------

LABELS_DC2 = ["DC2", "cDC2", "Dendritic cell type 2"]
LABELS_DP  = ["Double-positive thymocytes", "Double positive thymocytes", "DoublePositive", "CD4+CD8+"]

def log(*a): print(*a, flush=True)

def pick_embedding(adata: sc.AnnData) -> str:
    for rep in [USE_REP, "X_scVI", "X_pca"]:
        if rep in adata.obsm: return rep
    sc.pp.pca(adata, n_comps=min(50, adata.n_vars))
    return "X_pca"

def rank_markers_from_raw(sub: sc.AnnData, group_key: str = "leiden_sub"):
    """Compute rank_genes_groups from full genes in .raw, copy back to sub.uns."""
    if sub.raw is None:
        log("[WARN] .raw missing; ranking from HVGs only")
        tmp = sub.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        sc.tl.rank_genes_groups(tmp, groupby=group_key, method="wilcoxon", use_raw=False, pts=True)
        sub.uns["rank_genes_groups"] = tmp.uns["rank_genes_groups"]
        return

    tmp = sub.raw.to_adata()   # full gene matrix
    tmp.obs = sub.obs.copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    sc.tl.rank_genes_groups(tmp, groupby=group_key, method="wilcoxon", use_raw=False, pts=True)
    sub.uns["rank_genes_groups"] = tmp.uns["rank_genes_groups"]

def write_with_fullX(sub: sc.AnnData, out_path: str):
    """
    Create an AnnData whose .X and .var are the full-gene matrix from .raw,
    while preserving obs/obsm/uns/obsp/layers entries from the subset.
    """
    if sub.raw is None:
        log(f"[WARN] .raw not present; cannot write fullX variant for {out_path}")
        sub.write_h5ad(out_path)  # fallback
        return
    full = sub.raw.to_adata()          # full .X and .var
    # copy annotations/backed results
    full.obs  = sub.obs.copy()
    full.obsm = sub.obsm.copy()
    full.obsp = sub.obsp.copy() if hasattr(sub, "obsp") else {}
    full.uns  = sub.uns.copy()
    # keep counts if you had them in the parent HVG object
    for lname in list(sub.layers.keys()):
        try:
            full.layers[lname] = sub.layers[lname]
        except Exception:
            pass
    # drop raw inside full to avoid duplication (optional)
    full.raw = None
    full.write_h5ad(out_path)

def subset_and_cluster(adata: sc.AnnData, labels, tag: str):
    mask = adata.obs[LABEL_COL].astype(str).isin(labels)
    sub = adata[mask].copy()
    log(f"[INFO] {tag}: {sub.n_obs:,} cells in subset")
    if sub.n_obs == 0:
        log(f"[WARN] {tag}: no cells; skipping.")
        return None

    rep = pick_embedding(sub)
    sc.pp.neighbors(sub, use_rep=rep, n_neighbors=15)
    sc.tl.umap(sub)
    sc.tl.leiden(sub, resolution=LEIDEN_RES, key_added="leiden_sub")

    # rank markers from full gene space
    rank_markers_from_raw(sub, "leiden_sub")

    # write standard HVG-based .X (raw retained)
    out = f"{tag}_subset_reclustered.h5ad"
    sub.write_h5ad(out)
    log(f"[OK] wrote {out} (.X=HVGs, .raw=full)")

    # optional: write variant with full-gene .X
    if WRITE_FULLX_VARIANT:
        out_full = f"{tag}_subset_reclustered_fullX.h5ad"
        write_with_fullX(sub, out_full)
        log(f"[OK] wrote {out_full} (.X=FULL genes)")

    return sub

def main():
    log(f"[STEP] Reading {H5AD_PATH}")
    ad = sc.read(H5AD_PATH)
    log(f"[OK] adata: {ad.n_obs:,} cells, {ad.n_vars:,} HVGs | raw present: {ad.raw is not None}")

    if ad.raw is None:
        log("[WARN] Your input has no .raw; full-gene outputs will not be produced.")

    if LABEL_COL not in ad.obs.columns:
        raise KeyError(f"'{LABEL_COL}' not in .obs (have: {list(ad.obs.columns)[:25]})")

    subset_and_cluster(ad, LABELS_DC2, "DC2")
    subset_and_cluster(ad, LABELS_DP,  "DP")
    log("[DONE]")

if __name__ == "__main__":
    main()



