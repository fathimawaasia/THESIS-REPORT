#!/usr/bin/env python3
# scanorama_celltypist_pipeline_fixed.py
# - Attaches Group from metadata (GSM_ID -> Group) BEFORE integration
# - Sanitizes .X (drop NaN/Inf; drop empty cells/genes) before PCA
# - Uses use_raw=False when plotting obs colors
# - For CellTypist, re-normalizes a fresh copy to avoid raw-handling pitfalls

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scanorama
import celltypist
import scipy.sparse as sp
import matplotlib.pyplot as plt

# ========== Paths ==========
INPUT_RAW_HVG = "qc_outputs/adata_all_raw_set.h5ad"              # your pre-HVG filtered input (as you used)
NORMED_FILE   = "qc_outputs/adata_all_raw_set_normalized_3.h5ad" # normalized+log1p version (we create/update)
META_CSV      = "sample_metadata_7.csv"                          # uploaded metadata with columns GSM_ID, Group
OUTPUT_DIR    = "qc_outputs"
FIG_DIR       = "figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# ========== Helpers ==========

def sanitize_X(adata):
    """Make .X finite and drop empty cells/genes to keep sklearn happy."""
    X = adata.X
    if sp.issparse(X):
        X = X.tocsr(copy=True)
        d = X.data
        bad = ~np.isfinite(d)
        if bad.any():
            print(f"⚠️ Found {bad.sum()} non-finite entries in sparse .X; setting to 0.")
            d[bad] = 0.0
        adata.X = X
    else:
        n_nan = np.isnan(adata.X).sum()
        n_inf = np.isinf(adata.X).sum()
        if n_nan or n_inf:
            print(f"⚠️ Found {n_nan} NaNs and {n_inf} infs in dense .X; setting to 0.")
        np.nan_to_num(adata.X, copy=False, nan=0.0, posinf=0.0, neginf=0.0)

    # Drop all-zero cells/genes
    cell_sums = np.asarray(adata.X.sum(axis=1)).ravel()
    gene_sums = np.asarray(adata.X.sum(axis=0)).ravel()
    keep_cells = cell_sums > 0
    keep_genes = gene_sums > 0
    n_drop_cells = (~keep_cells).sum()
    n_drop_genes = (~keep_genes).sum()
    if n_drop_cells or n_drop_genes:
        print(f"🧹 Dropping {n_drop_cells} empty cells and {n_drop_genes} empty genes.")
        adata._inplace_subset_obs(keep_cells)
        adata._inplace_subset_var(keep_genes)

def attach_group_from_meta(adata, meta_csv):
    """Merge Group onto .obs using GSM_ID or sample_id if present."""
    if not os.path.exists(meta_csv):
        raise FileNotFoundError(f"Metadata CSV not found: {meta_csv}")
    meta = pd.read_csv(meta_csv)

    # Ensure 'Group' present in metadata
    if "Group" not in meta.columns:
        raise KeyError("Metadata must contain a 'Group' column.")
    # Pick a join key
    join_key = None
    if "GSM_ID" in adata.obs.columns and "GSM_ID" in meta.columns:
        join_key = "GSM_ID"
    elif "sample_id" in adata.obs.columns and "sample_id" in meta.columns:
        join_key = "sample_id"
    else:
        raise KeyError("Neither 'GSM_ID' nor 'sample_id' are present in both AnnData.obs and metadata for merging.")

    meta = meta[[join_key, "Group"]].drop_duplicates()
    before_n = adata.n_obs
    adata.obs = adata.obs.merge(meta.set_index(join_key), left_on=join_key, right_index=True, how="left")
    if adata.n_obs != before_n:
        raise RuntimeError("Unexpected change in number of observations after metadata merge.")
    adata.obs["Group"] = adata.obs["Group"].astype("category")
    print("Group value counts after merge:\n", adata.obs["Group"].value_counts(dropna=False))
    return adata

# ========== Step 0: Normalize & save (if needed) ==========
print("🔹 Loading raw HVG-filtered file...")
adata = sc.read(INPUT_RAW_HVG)

# Normalize/log1p and save normalized file
print("🔹 Normalizing (total=1e4) and log1p ...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()  # keep normalized values as raw for provenance
adata.write(NORMED_FILE)
print(f"✅ Saved log1p-normalized data with .raw attached → {NORMED_FILE}")

# From here on, work with the normalized file
adata_all = sc.read(NORMED_FILE)

# ========== Step 1: Attach Group BEFORE integration ==========
print("🔹 Attaching Group from metadata before integration ...")
adata_all = attach_group_from_meta(adata_all, META_CSV)

# ========== Step 2: PRE-Scanorama PCA/UMAP/Leiden ==========
print("🔹 Running PCA, UMAP, Leiden BEFORE Scanorama...")
sanitize_X(adata_all)

sc.tl.pca(adata_all, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata_all, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata_all)
sc.tl.leiden(adata_all, resolution=0.5, key_added="leiden_pre")

sc.pl.pca(adata_all, use_raw=False, save="_pre_scanorama_pca_3.png", show=False)
sc.pl.umap(adata_all, color=["sample_id"], use_raw=False,
           save="_pre_scanorama_umap_by_sample_3.png", show=False)
sc.pl.umap(adata_all, color=["leiden_pre"], legend_loc="on data", use_raw=False,
           save="_pre_scanorama_umap_by_leiden_3.png", show=False)

# ========== Step 3: Split by sample_id & run Scanorama ==========
print("🔹 Splitting by sample_id and running Scanorama integration...")
if "sample_id" not in adata_all.obs.columns:
    raise KeyError("Expected 'sample_id' in adata_all.obs for per-sample integration.")

adatas_split = []
for sid in adata_all.obs["sample_id"].unique():
    ad = adata_all[adata_all.obs["sample_id"] == sid].copy()
    adatas_split.append(ad)

# Correct/batch integrate; return_dimred gives X_scanorama in .obsm
corrected_adatas = scanorama.correct_scanpy(adatas_split, return_dimred=True)

# ========== Step 4: Concatenate integrated AnnDatas ==========
print("🔹 Concatenating corrected AnnData objects...")
sample_ids = []
for a in corrected_adatas:
    # ensure sample_id is a scalar per split
    sids = a.obs["sample_id"].unique().tolist()
    sample_ids.append(sids[0] if len(sids) == 1 else sids[0])

adata_scanorama = sc.concat(corrected_adatas, label="sample_id", keys=sample_ids)
adata_scanorama.obs_names_make_unique(join="_")

# ========== Step 5: Save integrated expression matrix & AnnData ==========
print("🔹 Saving integrated matrix and AnnData ... (this CSV can be large)")
X_dense = adata_scanorama.X.toarray() if not isinstance(adata_scanorama.X, np.ndarray) else adata_scanorama.X
df_expr = pd.DataFrame(X_dense, index=adata_scanorama.obs_names, columns=adata_scanorama.var_names)
df_expr.to_csv("scanorama_integrated_expression_matrix_3.csv")
adata_scanorama.write(os.path.join(OUTPUT_DIR, "adata_scanorama_integrated_3.h5ad"))
print("✅ Saved: scanorama_integrated_expression_matrix_3.csv")
print("✅ Saved: qc_outputs/adata_scanorama_integrated_3.h5ad")

# ========== Step 6: POST-Scanorama PCA/UMAP/Leiden & plots ==========
print("🔹 Running PCA/UMAP/Leiden AFTER Scanorama...")
adata_scanorama = sc.read(os.path.join(OUTPUT_DIR, "adata_scanorama_integrated_3.h5ad"))

# Assign Scanorama embeddings as PCA basis
if "X_scanorama" not in adata_scanorama.obsm:
    raise KeyError("X_scanorama embedding not found in .obsm after scanorama.")
adata_scanorama.obsm["X_pca"] = adata_scanorama.obsm["X_scanorama"]

# PCA colored by Group (ensure use_raw=False)
sc.pl.pca(
    adata_scanorama,
    color=["Group"],
    legend_loc="right margin",
    legend_fontsize=10,
    use_raw=False,
    save="_scanorama_Group_PCA_3.png",
)
print("PCA plot saved to figures/pca_scanorama_Group_PCA_3.png")

# Neighbors/UMAP/Leiden on X_pca
sc.pp.neighbors(adata_scanorama, use_rep="X_pca", n_neighbors=15, n_pcs=40)
sc.tl.umap(adata_scanorama)
sc.tl.leiden(adata_scanorama, resolution=0.5, flavor="igraph", directed=False)

# UMAPs
sc.set_figure_params(dpi=150)
sc.pl.umap(
    adata_scanorama,
    color=["Group", "sample_id", "leiden", "GSM_ID"],
    legend_loc="right margin",
    legend_fontsize=10,
    ncols=2,
    wspace=0.4,
    use_raw=False,
    save="_scanorama_integrated_UMAP_3.png",
)
sc.pl.umap(
    adata_scanorama,
    color=["leiden", "Group"],
    legend_loc="right margin",
    legend_fontsize=10,
    ncols=2,
    wspace=0.4,
    use_raw=False,
    save="_scanorama_leiden_Group_3.png",
)

# Save clustered object
adata_scanorama.write(os.path.join(OUTPUT_DIR, "adata_scanorama_clustered_umap_3.h5ad"))
print("✅ UMAP & clustering complete. Saved → qc_outputs/adata_scanorama_clustered_umap_3.h5ad")

# ========== Step 7: CellTypist on a re-normalized copy ==========
# ========== Step 7: CellTypist on a cleaned, re-normalized copy ==========
# ========== Step 7: CellTypist on a correctly scaled matrix ==========
# ========== Step 7: CellTypist on a correctly scaled matrix ==========
print("🔹 Preparing data for CellTypist ...")

import numpy as np, scipy.sparse as sp, pandas as pd, celltypist

def _count_nonfinite(A):
    if sp.issparse(A):
        d = A.data
        return np.isnan(d).sum(), np.isinf(d).sum()
    else:
        return np.isnan(A).sum(), np.isinf(A).sum()

def _ensure_csr_f32(adata):
    if sp.issparse(adata.X):
        adata.X = adata.X.tocsr().astype(np.float32)
    else:
        adata.X = sp.csr_matrix(adata.X.astype(np.float32))

# Prefer using .raw if present (should be log1p-normalized to 1e4 from Step 0)
if adata_scanorama.raw is not None:
    print("   Using .raw as CellTypist input (expected log1p to 1e4).")
    adata_ct = adata_scanorama.raw.to_adata()
    adata_ct.obs = adata_scanorama.obs.copy()  # keep obs alignment
else:
    print("   No .raw on integrated object; using .X as-is (assumed log1p to 1e4).")
    adata_ct = adata_scanorama.copy()

# Sanitize BEFORE any scaling checks (replace NaN/Inf, drop empty rows/cols)
n_nan, n_inf = _count_nonfinite(adata_ct.X)
print(f"   Non-finite before clean → NaN: {n_nan:,}  Inf: {n_inf:,}")
sanitize_X(adata_ct)   # <- you already defined this earlier
n_nan, n_inf = _count_nonfinite(adata_ct.X)
print(f"   After initial clean      → NaN: {n_nan:,}  Inf: {n_inf:,}")

# Clamp negatives to 0 (defensive)
if sp.issparse(adata_ct.X):
    d = adata_ct.X.data
    n_neg = (d < 0).sum()
    if n_neg:
        print(f"   Clamping {n_neg:,} negative entries to 0.")
        d[d < 0] = 0.0
else:
    n_neg = (adata_ct.X < 0).sum()
    if n_neg:
        print(f"   Clamping {int(n_neg):,} negative entries to 0.")
        adata_ct.X[adata_ct.X < 0] = 0.0

# Check scale: median total counts after expm1 should be ~1e4
if sp.issparse(adata_ct.X):
    X_tmp = adata_ct.X.tocsr(copy=True)
    X_tmp.data = np.expm1(X_tmp.data)
    totals = np.asarray(X_tmp.sum(axis=1)).ravel()
else:
    totals = np.expm1(adata_ct.X).sum(axis=1)
med_total = float(np.median(totals))
print(f"   Median total counts (expm1) ≈ {med_total:,.1f} (expected ~10,000)")

if not (3_000 <= med_total <= 30_000):
    print("   ⚠️ Totals off-scale; re-normalizing from un-logged values to 1e4 then log1p.")
    if sp.issparse(adata_ct.X):
        X_unlogged = adata_ct.X.tocsr(copy=True)
        X_unlogged.data = np.expm1(X_unlogged.data)
        adata_ct.X = X_unlogged
    else:
        adata_ct.X = np.expm1(adata_ct.X)
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct)
    sanitize_X(adata_ct)
    _ensure_csr_f32(adata_ct)
    # Re-check
    if sp.issparse(adata_ct.X):
        X_tmp = adata_ct.X.tocsr(copy=True)
        X_tmp.data = np.expm1(X_tmp.data)
        totals = np.asarray(X_tmp.sum(axis=1)).ravel()
    else:
        totals = np.expm1(adata_ct.X).sum(axis=1)
    med_total = float(np.median(totals))
    print(f"   Median total counts after fix ≈ {med_total:,.1f}")
else:
    _ensure_csr_f32(adata_ct)

print("🔹 Running CellTypist annotation...")
celltypist.models.download_models()  # no-op if already present
model = "Immune_All_Low.pkl"
pred = celltypist.annotate(adata_ct, model=model, majority_voting=True)

# === Assign predictions back to adata_scanorama.obs ===
pl = pred.predicted_labels
if isinstance(pl, pd.DataFrame):
    pl = pl.iloc[:, 0]
pl = pl.reindex(adata_scanorama.obs_names)
adata_scanorama.obs["celltypist_labels"] = pl.astype(str).values

if hasattr(pred, "majority_voting") and pred.majority_voting is not None:
    mv = pred.majority_voting
    if isinstance(mv, pd.DataFrame):
        mv = mv.iloc[:, 0]
    mv = mv.reindex(adata_scanorama.obs_names)
    adata_scanorama.obs["celltypist_labels_mv"] = mv.astype(str).values

if hasattr(pred, "probability") and pred.probability is not None:
    probs = pred.probability.reindex(index=adata_scanorama.obs_names)
    adata_scanorama.obs["celltypist_confidence"] = probs.max(axis=1).values

# ========== Step 8: Save annotated data (guarded) ==========
print("🔹 Saving CellTypist outputs ...")
annot_h5ad = os.path.join(OUTPUT_DIR, "adata_celltypist_annotated_3.h5ad")
adata_scanorama.write(annot_h5ad)
print(f"✅ Saved {annot_h5ad}")

if "celltypist_labels" in adata_scanorama.obs.columns:
    annot_csv = os.path.join(OUTPUT_DIR, "celltypist_labels_all_3.csv")
    adata_scanorama.obs[["celltypist_labels"]].to_csv(annot_csv)
    print(f"✅ Saved {annot_csv}")
else:
    print("⚠️ 'celltypist_labels' not found in .obs; skipping CSV export.")

# ========== Step 9: UMAPs with CellTypist/Leiden (guarded) ==========
print("🔹 Plotting UMAPs with CellTypist and Leiden labels ...")
if "celltypist_labels" in adata_scanorama.obs.columns:
    sc.pl.umap(adata_scanorama, color=["celltypist_labels"], use_raw=False,
               save="_celltypist_labels_3.png", show=False)
else:
    print("⚠️ Skipping UMAP colored by 'celltypist_labels' (column missing).")

sc.pl.umap(adata_scanorama, color=["leiden"], use_raw=False,
           save="_leiden_clusters_3.png", show=False)

# ========== Step 10: Cluster-majority labels + immune split (guarded) ==========
if "celltypist_labels" in adata_scanorama.obs.columns:
    print("🔹 Computing cluster-level majority CellTypist labels ...")
    cluster_labels = (
        adata_scanorama.obs.groupby("leiden")["celltypist_labels"]
        .agg(lambda x: x.value_counts().index[0])
    )
    adata_scanorama.obs["cluster_identity"] = adata_scanorama.obs["leiden"].map(cluster_labels)

    sc.pl.umap(
        adata_scanorama,
        color="cluster_identity",
        legend_loc="on data",
        legend_fontsize=8,
        size=20,
        use_raw=False,
        save="_cluster_identity_umap_3.png",
        show=False
    )

    print("🔹 Splitting immune vs non-immune cells ...")
    immune_keywords = ["T cell", "B cell", "Monocyte", "Macrophage", "NK", "DC", "Neutrophil"]
    immune_mask = adata_scanorama.obs["celltypist_labels"].str.contains("|".join(immune_keywords), case=False, na=False)

    adata_immune    = adata_scanorama[immune_mask].copy()
    adata_nonimmune = adata_scanorama[~immune_mask].copy()

    adata_immune.write(os.path.join(OUTPUT_DIR, "adata_celltypist_immune_3.h5ad"))
    adata_nonimmune.write(os.path.join(OUTPUT_DIR, "adata_celltypist_nonimmune_3.h5ad"))
    adata_immune.obs[["celltypist_labels"]].to_csv(os.path.join(OUTPUT_DIR, "celltypist_labels_immune_3.csv"))
    adata_nonimmune.obs[["celltypist_labels"]].to_csv(os.path.join(OUTPUT_DIR, "celltypist_labels_nonimmune_3.csv"))
    print("✅ Immune/non-immune splits saved.")
else:
    print("⚠️ 'celltypist_labels' missing; skipping cluster-identity UMAP and immune split.")
