#!/usr/bin/env python3
# CellTypist annotation + cross-dataset QC + BayesPrism reference (CPM means, filtered)

import os, math, itertools
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import celltypist
from celltypist import models

# -------------------- CONFIG --------------------
INP   = "adata_all_scvi_integrated.h5ad"   # same folder as this script
AOUT  = "adata_celltypist_qc.h5ad"
OUTD  = "celltypist_qc"
LOGD  = "logs"

SIG_TSV    = f"{OUTD}/reference_celltypist_meanCPM.filtered.tsv.gz"
REPORT_CSV = f"{OUTD}/qc_report_celltypist.csv"

LABEL_COL    = "ct_seed"          # seed labels to feed scANVI later
CONF_COL     = "celltypist_confidence"
BATCH_KEY    = "dataset_id"       # must exist in .obs
COUNTS_LAYER = "counts"           # must be counts if present

# thresholds
CONF_MIN           = 0.70          # keep cells with conf >= this for seeding & reference
MIN_DATASETS       = 2             # cell type must appear in at least this many datasets
MIN_PER_CT_TOTAL   = 100           # minimum cells across all datasets for a cell type
MIN_PER_CT_DATASET = 30            # minimum cells per dataset to compute a mean
CORR_MIN           = 0.85          # minimum pairwise corr among dataset means for the cell type
CPM_GENE_MIN       = 1.0           # keep genes with >= 1 CPM in at least one cell type
# ------------------------------------------------

os.makedirs(OUTD, exist_ok=True)
os.makedirs(LOGD, exist_ok=True)

# ------------- load -------------
print(f"[LOAD] Reading: {INP}")
ad = sc.read(INP)
print("[LOAD]", ad)

if BATCH_KEY not in ad.obs:
    raise KeyError(f"obs['{BATCH_KEY}'] not found. Add dataset identifiers before running.")

# ------------- CellTypist (expects log1p-normalized X to 1e4) -------------
# ------------- CellTypist (expects log1p-normalized X to 1e4) -------------
from pathlib import Path

# Build a working copy and put counts into .X
work = ad.copy()
if COUNTS_LAYER in work.layers:
    Xc = work.layers[COUNTS_LAYER]
else:
    print("[WARN] No counts layer; assuming .X are counts.")
    Xc = work.X
work.X = Xc  # set counts as X for normalization

# Normalize to 1e4 per cell, then log1p — required by CellTypist
sc.pp.normalize_total(work, target_sum=1e4, inplace=True)
sc.pp.log1p(work)

# Ensure we have cluster labels for majority voting
if "leiden" in ad.obs:
    mv_labels = ad.obs["leiden"].astype(str).values
    print("[MV] Using existing ad.obs['leiden'] for majority voting.")
else:
    print("[MV] ad.obs['leiden'] not found — computing lightweight Leiden on normalized copy.")
    # Prefer X_scVI if available (copied from ad to work)
    use_rep = "X_scVI" if "X_scVI" in work.obsm else None
    if use_rep is None:
        sc.pp.pca(work, n_comps=50, svd_solver="arpack")
        sc.pp.neighbors(work, n_neighbors=15)
    else:
        sc.pp.neighbors(work, n_neighbors=15, use_rep=use_rep)
    sc.tl.leiden(work, resolution=1.0)
    mv_labels = work.obs["leiden"].astype(str).values

# Download models (no-op if already present) and select Immune_All_Low
models.download_models()
model_path = Path(models.models_path) / "Immune_All_Low.pkl"
assert model_path.exists(), f"Model not found at: {model_path}"

pred = celltypist.annotate(
    work,
    model=str(model_path),
    majority_voting=True,
    over_clustering=mv_labels  # MUST be array-like with len == n_cells
)
pred_ad = pred.to_adata()

# Attach predictions back to original AnnData
ad.obs["celltypist_labels"]   = pred_ad.obs["predicted_labels"].reindex(ad.obs_names)
ad.obs[CONF_COL]              = pred_ad.obs["conf_score"].reindex(ad.obs_names).astype(float)
ad.obs["celltypist_majority"] = pred_ad.obs["majority_voting"].reindex(ad.obs_names)

# ---- SAFE string handling: no nullable StringDtype anywhere ----
# 1) Coerce CellTypist string columns to plain Python strings (object) and fill NAs
for col in ["celltypist_labels", "celltypist_majority"]:
    if col in ad.obs:
        ad.obs[col] = ad.obs[col].astype(object).fillna("Unknown")

# 2) Build seed labels with confidence gate, keep plain object -> then categorical
seed = ad.obs["celltypist_majority"].astype(object)
seed = seed.where(ad.obs[CONF_COL] >= CONF_MIN, other="Unknown").fillna("Unknown")
ad.obs[LABEL_COL] = pd.Categorical(seed)   # e.g., 'ct_seed'

# (Optional belt-and-suspenders: demote *any* nullable string columns to object)
import pandas as pd
from pandas.api.types import is_string_dtype
for c in ad.obs.columns:
    if is_string_dtype(ad.obs[c].dtype):   # catches pandas StringDtype
        ad.obs[c] = ad.obs[c].astype(object)

# ------------- write annotated h5ad -------------
ad.write(AOUT, compression="gzip")
print("[WRITE] Annotated:", AOUT)

# ------------- Prepare CPM (no log) -------------
work2 = ad.copy()
if COUNTS_LAYER in work2.layers:
    Xc = work2.layers[COUNTS_LAYER]
else:
    print("[WARN] No counts layer; using .X as counts for CPM.")
    Xc = work2.X

if not sparse.issparse(Xc):
    Xc = sparse.csr_matrix(Xc)

# CPM normalization without densifying the whole matrix
libsize = np.asarray(Xc.sum(axis=1)).ravel()
libsize[libsize == 0] = 1.0
sf = libsize / 1e6
inv_sf = 1.0 / sf
X_cpm = Xc.multiply(inv_sf[:, None]).tocsr()

labels  = ad.obs[LABEL_COL].astype(str).values
batches = ad.obs[BATCH_KEY].astype(str).values
genes   = ad.var_names.to_list()

# ------------- Cross-dataset QC per cell type -------------
def safe_corr(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if np.allclose(a, a[0]) or np.allclose(b, b[0]):
        return np.nan
    return np.corrcoef(a, b)[0, 1]

def mean_cpm_for_indices(idx_bool):
    if idx_bool.sum() == 0:
        return None
    v = X_cpm[idx_bool, :].mean(axis=0)
    v = v.A1 if sparse.issparse(v) else np.asarray(v).ravel()
    return v

qc_rows = []
ctypes = sorted(set(labels) - {"Unknown"})
kept_celltypes = []
dataset_presence = {}

for ct in ctypes:
    idx_ct = (labels == ct)
    n_total = int(idx_ct.sum())
    if n_total < MIN_PER_CT_TOTAL:
        qc_rows.append([ct, n_total, 0, np.nan, False, "total_cells_lt_min"])
        continue

    batches_ct, counts_ct = np.unique(batches[idx_ct], return_counts=True)
    valid_mask = counts_ct >= MIN_PER_CT_DATASET
    valid_batches = list(batches_ct[valid_mask])
    n_valid = len(valid_batches)

    if n_valid < MIN_DATASETS:
        qc_rows.append([ct, n_total, n_valid, np.nan, False, "datasets_lt_min"])
        continue

    means = {}
    for b in valid_batches:
        idx = (labels == ct) & (batches == b)
        means[b] = mean_cpm_for_indices(idx)

    corrs = []
    for (b1, b2) in itertools.combinations(valid_batches, 2):
        r = safe_corr(means[b1], means[b2])
        if not np.isnan(r):
            corrs.append(r)

    min_corr = float(np.nanmin(corrs)) if len(corrs) else np.nan
    pass_corr = (len(corrs) > 0) and (min_corr >= CORR_MIN)
    reason = "ok" if pass_corr else "low_pairwise_corr"

    qc_rows.append([ct, n_total, n_valid, min_corr, pass_corr, reason])
    if pass_corr:
        kept_celltypes.append(ct)
        dataset_presence[ct] = valid_batches

qc_df = pd.DataFrame(qc_rows, columns=[
    "cell_type", "n_cells_total", "n_datasets_valid", "min_pairwise_corr", "pass", "reason"
])
qc_df.to_csv(REPORT_CSV, index=False)
print(f"[QC] Wrote cross-dataset report: {REPORT_CSV}")
print(qc_df.head(20))

# ------------- Build filtered BayesPrism reference -------------
print(f"[REF] Building reference from {len(kept_celltypes)} cell types (passed QC)…")
ref_cols, ref_mat = [], []

for ct in kept_celltypes:
    idx = (labels == ct) & (ad.obs[CONF_COL].values >= CONF_MIN)
    if idx.sum() < MIN_PER_CT_TOTAL:
        continue
    vec = mean_cpm_for_indices(idx)
    if vec is None:
        continue
    ref_cols.append(ct)
    ref_mat.append(vec)

if len(ref_cols) == 0:
    raise RuntimeError("No cell types passed QC to build reference. Relax thresholds or inspect qc_report.")

ref = pd.DataFrame(np.vstack(ref_mat).T, index=genes, columns=ref_cols)
keep_genes = (ref.max(axis=1) >= CPM_GENE_MIN)
ref = ref.loc[keep_genes]

ref.to_csv(SIG_TSV, sep="\t", compression="gzip")
print(f"[REF] Wrote filtered reference: {SIG_TSV}")

print("✅ Done CellTypist+QC.")

# ---------------- UMAP visualization (cell types & clusters) ----------------
import matplotlib.pyplot as plt

FIGD = f"{OUTD}/figs"
os.makedirs(FIGD, exist_ok=True)

def ensure_neighbors_umap(adata):
    """
    Ensure neighbors/UMAP exist; prefer X_scVI if available.
    Recompute if neighbors exist but were built on a different representation.
    """
    use_rep = None
    if "X_scVI" in adata.obsm:
        use_rep = "X_scVI"

    needs_neighbors = (
        "neighbors" not in adata.uns or
        (use_rep is not None and adata.uns["neighbors"].get("params", {}).get("use_rep") != use_rep)
    )
    if needs_neighbors:
        if use_rep is None:
            # Fall back to PCA on .X (assumes integrated/log space adequate for viz)
            sc.pp.pca(adata, n_comps=50, svd_solver="arpack")
            sc.pp.neighbors(adata, n_neighbors=15)
        else:
            sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep)
        sc.tl.umap(adata)

# Make sure UMAP is ready
ensure_neighbors_umap(ad)

# Scanpy figure defaults
sc.set_figure_params(dpi=120, figsize=(5.0, 5.0))
sc.settings.figdir = FIGD  # so sc.pl.umap(..., save="...") writes into FIGD

# Helper: safe color list if some keys are missing
def existing_keys(keys):
    return [k for k in keys if k in ad.obs_keys()]

# 1) UMAP by CellTypist outputs
ct_keys = existing_keys(["celltypist_labels", "celltypist_majority", LABEL_COL])
if ct_keys:
    sc.pl.umap(
        ad,
        color=ct_keys,
        ncols=2,
        wspace=0.35,
        frameon=False,
        legend_loc="on data",
        legend_fontsize=8,
        show=False,
        save="_celltypist_labels.png",
    )

# 2) UMAP by cluster identity (Leiden if present)
if "leiden" in ad.obs:
    sc.pl.umap(
        ad,
        color=["leiden"],
        frameon=False,
        legend_loc="on data",
        legend_fontsize=8,
        show=False,
        save="_leiden.png",
    )

# 3) UMAP by dataset/batch and optional group/condition
to_try = existing_keys([BATCH_KEY, "Group", "condition", "study_id"])
if to_try:
    sc.pl.umap(
        ad,
        color=to_try,
        ncols=2,
        wspace=0.35,
        frameon=False,
        legend_loc="right margin",
        legend_fontsize=8,
        show=False,
        save="_meta.png",
    )

print(f"[PLOTS] Saved UMAPs to: {FIGD}")


