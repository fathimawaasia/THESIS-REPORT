#!/usr/bin/env python3
# scANVI refinement + cross-dataset QC + BayesPrism reference (CPM means, filtered)

import os, math, itertools
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scipy import sparse
from pandas.api.types import is_string_dtype

# -------------------- CONFIG --------------------
INP       = "../cell_typist/adata_celltypist_qc.h5ad"     # parent dir
SCVI_DIR  = "../qc_outputs/scvi_model"                    # parent dir
OUT_AD    = "../qc_outputs/adata_all_scanvi_refined_qc.h5ad"
OUTD      = "../qc_outputs/scanvi_qc"
SIG_TSV   = f"{OUTD}/reference_scanvi_meanCPM.filtered.tsv.gz"
REPORT_CSV= f"{OUTD}/qc_report_scanvi.csv"

LABEL_SEED     = "ct_seed"
FINAL_LABEL    = "cell_type_final"
BATCH_KEY      = "dataset_id"
COUNTS_LAYER   = "counts"

# thresholds
CONF_MIN_SEED  = 0.50
CONF_MIN_FINAL = 0.60
MIN_DATASETS   = 2
MIN_PER_CT_TOTAL   = 100
MIN_PER_CT_DATASET = 30
CORR_MIN       = 0.85
CPM_GENE_MIN   = 1.0
# ------------------------------------------------

os.makedirs(OUTD, exist_ok=True)

print(f"[LOAD] {INP}")
ad = sc.read(INP)
print("[LOAD]", ad)

if LABEL_SEED not in ad.obs:
    raise ValueError(f"obs['{LABEL_SEED}'] not found; run the CellTypist step first.")
if BATCH_KEY not in ad.obs:
    raise KeyError(f"obs['{BATCH_KEY}'] not found in {INP}")

# Gate low-confidence seeds
if "celltypist_confidence" in ad.obs:
    conf = ad.obs["celltypist_confidence"].astype(float).values
    seed = ad.obs[LABEL_SEED].astype(object)
    seed = pd.Series(seed, index=ad.obs_names).where(conf >= CONF_MIN_SEED, other="Unknown").fillna("Unknown")
    ad.obs[LABEL_SEED] = pd.Categorical(seed)

# ---- Load SCVI and build SCANVI
print("[INFO] Loading SCVI model from:", SCVI_DIR)
if not os.path.isdir(SCVI_DIR):
    # fallback if someone trained under scANVI/scvi_model
    alt = "scvi_model"
    if os.path.isdir(alt):
        SCVI_DIR = alt
        print("[WARN] Falling back to:", SCVI_DIR)
    else:
        raise FileNotFoundError(f"SCVI model dir not found at {SCVI_DIR} or {alt}")

scvi_model = scvi.model.SCVI.load(SCVI_DIR, adata=ad)

scanvi = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key=LABEL_SEED,
    unlabeled_category="Unknown",
)

print("[TRAIN] SCANVI refining labels…")
scanvi.train(
    max_epochs=100,
    early_stopping=True,
    early_stopping_patience=15,
    batch_size=4096,
)

# Predictions + confidence
ad.obs[FINAL_LABEL] = scanvi.predict()
soft = scanvi.predict(soft=True)
ad.obs["scanvi_confidence"] = soft.max(axis=1).values

# Latent + UMAP
ad.obsm["X_scANVI"] = scanvi.get_latent_representation()
sc.pp.neighbors(ad, use_rep="X_scANVI", n_neighbors=15)
sc.tl.umap(ad, min_dist=0.3, spread=1.0)

# Sanitize strings before write
for c in ad.obs.columns:
    if is_string_dtype(ad.obs[c].dtype):
        ad.obs[c] = ad.obs[c].astype(object)
ad.obs[FINAL_LABEL] = pd.Categorical(ad.obs[FINAL_LABEL].astype(object))

ad.write(OUT_AD, compression="gzip")
print("[WRITE] SCANVI-refined:", OUT_AD)

# ---- Build CPM-based reference with QC
work = ad.copy()
if COUNTS_LAYER in work.layers:
    Xc = work.layers[COUNTS_LAYER]
else:
    print("[WARN] No counts layer; using .X as counts.")
    Xc = work.X
if not sparse.issparse(Xc):
    Xc = sparse.csr_matrix(Xc)

libsize = np.asarray(Xc.sum(axis=1)).ravel()
libsize[libsize == 0] = 1.0
X_cpm = Xc.multiply((1.0 / (libsize / 1e6))[:, None]).tocsr()

labels = ad.obs[FINAL_LABEL].astype(str).values
batches = ad.obs["dataset_id"].astype(str).values
conf_final = ad.obs["scanvi_confidence"].astype(float).values
genes = ad.var_names.to_list()

def safe_corr(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if np.allclose(a, a[0]) or np.allclose(b, b[0]): return np.nan
    return np.corrcoef(a, b)[0, 1]

def mean_cpm(idx):
    if idx.sum() == 0: return None
    v = X_cpm[idx, :].mean(axis=0)
    return v.A1 if sparse.issparse(v) else np.asarray(v).ravel()

qc_rows, kept = [], []
for ct in sorted(set(labels) - {"Unknown"}):
    idx_ct = (labels == ct)
    n_total = int(idx_ct.sum())
    if n_total < MIN_PER_CT_TOTAL:
        qc_rows.append([ct, n_total, 0, np.nan, False, "total_cells_lt_min"]); continue
    batches_ct, counts_ct = np.unique(batches[idx_ct], return_counts=True)
    valid_batches = list(batches_ct[counts_ct >= MIN_PER_CT_DATASET])
    if len(valid_batches) < MIN_DATASETS:
        qc_rows.append([ct, n_total, len(valid_batches), np.nan, False, "datasets_lt_min"]); continue
    means, corrs = {}, []
    for b in valid_batches:
        idx = (labels == ct) & (batches == b) & (conf_final >= CONF_MIN_FINAL)
        means[b] = mean_cpm(idx)
    for i in range(len(valid_batches)):
        for j in range(i+1, len(valid_batches)):
            r = safe_corr(means[valid_batches[i]], means[valid_batches[j]])
            if not np.isnan(r): corrs.append(r)
    min_corr = float(np.nanmin(corrs)) if corrs else np.nan
    passed = (len(corrs) > 0) and (min_corr >= CORR_MIN)
    qc_rows.append([ct, n_total, len(valid_batches), min_corr, passed, "ok" if passed else "low_pairwise_corr"])
    if passed: kept.append(ct)

pd.DataFrame(qc_rows, columns=["cell_type","n_cells_total","n_datasets_valid","min_pairwise_corr","pass","reason"])\
  .to_csv(REPORT_CSV, index=False)
print(f"[QC] Wrote cross-dataset report: {REPORT_CSV}")

print(f"[REF] Building reference from {len(kept)} cell types…")
ref_cols, ref_mat = [], []
for ct in kept:
    idx = (labels == ct) & (conf_final >= CONF_MIN_FINAL)
    if idx.sum() < MIN_PER_CT_TOTAL: continue
    v = mean_cpm(idx)
    if v is None: continue
    ref_cols.append(ct); ref_mat.append(v)
if not ref_cols:
    raise RuntimeError("No cell types passed QC to build reference.")

ref = pd.DataFrame(np.vstack(ref_mat).T, index=genes, columns=ref_cols)
ref = ref.loc[ref.max(axis=1) >= CPM_GENE_MIN]
ref.to_csv(SIG_TSV, sep="\t", compression="gzip")
print(f"[REF] Wrote filtered reference: {SIG_TSV}")

print("✅ Done scANVI+QC.")


