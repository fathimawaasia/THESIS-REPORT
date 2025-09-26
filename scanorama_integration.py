#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scanorama

# ---------- Paths ----------
IN_H5AD = "qc_outputs/adata_all_HVGfiltered_meta.h5ad"
OUT_DIR = "qc_outputs"
FIG_DIR = "figures"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR

def is_log1p_like(adata, sample=5000):
    """Heuristic: if values look already log1p-normalized, skip re-log."""
    X = adata.X
    if not isinstance(X, np.ndarray):
        X = X[:sample].toarray() if X.shape[0] > sample else X.toarray()
    else:
        X = X[:sample] if X.shape[0] > sample else X
    mx = float(np.nanmax(X))
    # raw counts often have max >> 50; log1p tends to be ~<15
    return mx < 20

print("▶ Loading:", IN_H5AD)
adata = sc.read(IN_H5AD)

# ---------- 0) Normalize + set .raw (only if not already log1p) ----------
if is_log1p_like(adata):
    print("ℹ️  adata.X appears already log1p-like; keeping as-is")
else:
    print("▶ Normalizing + log1p …")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

# Freeze counts for later CellTypist
adata.raw = adata.copy()
adata.write(os.path.join(OUT_DIR, "adata_all_HVGfiltered_meta_log1p_rawset.h5ad"))

# ---------- 1) PCA/UMAP/Leiden BEFORE Scanorama (QC) ----------
print("▶ PCA/UMAP/Leiden (pre-Scanorama) …")
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_pre")

sc.pl.pca(adata, color="sample_id", save="_pre_scanorama_pca.png", show=False)
sc.pl.umap(adata, color=["sample_id","leiden_pre"], legend_loc="right margin",
           save="_pre_scanorama_umap.png", show=False)

# ---------- 2) Split by sample + run Scanorama ----------
print("▶ Splitting by sample_id and running Scanorama …")
splits, raw_splits = [], []
for sid in adata.obs["sample_id"].unique():
    ad = adata[adata.obs["sample_id"] == sid].copy()
    splits.append(ad)
    raw_splits.append(ad.copy())  # keep per-sample .raw source

print("▶ Scanorama correction …")
corrected = scanorama.correct_scanpy(splits, return_dimred=True)

# Reattach per-sample .raw to preserve normalized counts
for i in range(len(corrected)):
    corrected[i].raw = raw_splits[i]

# ---------- 3) Concatenate corrected sets ----------
print("▶ Concatenating corrected AnnData objects …")
keys = [a.obs["sample_id"][0] for a in corrected]
adata_int = sc.concat(corrected, label="sample_id", keys=keys)
adata_int.obs_names_make_unique()

# ---------- 4) Graph on Scanorama embeddings + UMAP/Leiden ----------
adata_int.obsm["X_pca"] = adata_int.obsm["X_scanorama"]
sc.pp.neighbors(adata_int, use_rep="X_pca", n_neighbors=15, n_pcs=40)
sc.tl.umap(adata_int)
sc.tl.leiden(adata_int, resolution=0.5, key_added="leiden")

# Plots
sc.pl.pca(adata_int, color=["sample_id"], save="_scanorama_pca.png", show=False)
sc.pl.umap(adata_int, color=["sample_id","tissue","donor_id","leiden"],
           ncols=2, legend_loc="right margin", save="_scanorama_umap.png", show=False)

# ---------- 5) Save outputs ----------
out_h5ad = os.path.join(OUT_DIR, "adata_scanorama_integrated.h5ad")
adata_int.write(out_h5ad)

# Optional: export integrated expression matrix (dense)
df_expr = pd.DataFrame(
    adata_int.X.toarray() if not isinstance(adata_int.X, np.ndarray) else adata_int.X,
    index=adata_int.obs_names,
    columns=adata_int.var_names
)
df_expr.to_csv(os.path.join(OUT_DIR, "scanorama_integrated_expression_matrix.csv"))

print("✅ Wrote:", out_h5ad)
print("✅ Figures in:", FIG_DIR)
