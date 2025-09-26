# scanorama_celltypist_pipeline.py

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scanorama
import celltypist

# ---------- Paths ----------
in_h5ad   = "qc_outputs/adata_all_HVGfiltered_meta.h5ad"
out_dir   = "qc_outputs"
fig_dir   = "figures"
os.makedirs(out_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)
sc.settings.figdir = fig_dir

# ---------- 0) Load & (re)normalize, set .raw ----------
print("▶ Loading HVG-filtered + meta AnnData …")
adata = sc.read(in_h5ad)

# (Re)normalize to be explicit/consistent for CellTypist
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Preserve a frozen copy for CellTypist later
adata.raw = adata.copy()
adata.write(os.path.join(out_dir, "adata_all_HVGfiltered_meta_log1p_rawset.h5ad"))

# ---------- 1) Quick PCA/UMAP/Leiden BEFORE Scanorama (optional QC) ----------
print("▶ PCA/UMAP/Leiden (pre-Scanorama) …")
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_pre")

sc.pl.pca(adata, color="sample_id", save="_pre_scanorama_pca.png", show=False)
sc.pl.umap(adata, color=["sample_id","leiden_pre"], legend_loc="right margin",
           save="_pre_scanorama_umap.png", show=False)

# ---------- 2) Split by sample and run Scanorama ----------
print("▶ Splitting by sample_id and running Scanorama …")
splits = []
raw_splits = []
for sid in adata.obs["sample_id"].unique():
    ad = adata[adata.obs["sample_id"] == sid].copy()
    splits.append(ad)
    raw_splits.append(ad.copy())  # keep a per-sample raw copy

# Integrate; return_dimred=True gives X_scanorama embeddings in .obsm
corrected = scanorama.correct_scanpy(splits, return_dimred=True)

# Reattach raw to each corrected piece
for i in range(len(corrected)):
    corrected[i].raw = raw_splits[i]

# ---------- 3) Concatenate corrected sets ----------
print("▶ Concatenating corrected AnnData objects …")
sample_keys = [a.obs["sample_id"][0] for a in corrected]
adata_int = sc.concat(corrected, label="sample_id", keys=sample_keys)
adata_int.obs_names_make_unique()

# ---------- 4) Use Scanorama embeddings for graph/UMAP/Leiden ----------
adata_int.obsm["X_pca"] = adata_int.obsm["X_scanorama"]
sc.pp.neighbors(adata_int, use_rep="X_pca", n_neighbors=15, n_pcs=40)
sc.tl.umap(adata_int)
sc.tl.leiden(adata_int, resolution=0.5, key_added="leiden")

# Plots (use columns that exist: sample_id / donor_id / tissue / leiden)
sc.pl.pca(adata_int, color=["sample_id"], save="_scanorama_pca.png", show=False)
sc.pl.umap(adata_int, color=["sample_id","tissue","donor_id","leiden"],
           ncols=2, legend_loc="right margin", save="_scanorama_umap.png", show=False)

# Save integrated object
adata_int.write(os.path.join(out_dir, "adata_scanorama_integrated.h5ad"))

# Optional: export integrated expression matrix (dense)
df_expr = pd.DataFrame(
    adata_int.X.toarray() if not isinstance(adata_int.X, np.ndarray) else adata_int.X,
    index=adata_int.obs_names,
    columns=adata_int.var_names
)
df_expr.to_csv(os.path.join(out_dir, "scanorama_integrated_expression_matrix.csv"))

# ---------- 5) CellTypist on the log1p-normalized data ----------
# Ensure we pass log1p-normalized counts; easiest is to restore from .raw (set at the top).
print("▶ Restoring log1p-normalized counts for CellTypist …")
adata_int.X   = adata_int.raw.X.copy()
adata_int.var = adata_int.raw.var.copy()

print("▶ Running CellTypist … (downloads models if missing)")
# celltypist.models.download_models()  # uncomment once to download all models
model = "Immune_All_Low.pkl"  # pick model; change if you need a different atlas
pred = celltypist.annotate(adata_int, model=model, majority_voting=True)

# Add labels to obs
labels = pred.predicted_labels
if isinstance(labels, pd.DataFrame):
    labels = labels.iloc[:, 0]
adata_int.obs["celltypist_labels"] = labels.astype(str).values

if getattr(pred, "majority_voting", None) is not None:
    mv = pred.majority_voting
    if isinstance(mv, pd.DataFrame):
        mv = mv.iloc[:, 0]
    adata_int.obs["celltypist_labels_mv"] = mv.astype(str).values

if getattr(pred, "probability", None) is not None:
    adata_int.obs["celltypist_confidence"] = pred.probability.max(axis=1).values

# Save annotated object + labels CSV
annot_h5ad = os.path.join(out_dir, "adata_celltypist_annotated.h5ad")
annot_csv  = os.path.join(out_dir, "celltypist_labels_all.csv")
adata_int.write(annot_h5ad)
adata_int.obs[["celltypist_labels"]].to_csv(annot_csv)

# UMAPs colored by CellTypist & Leiden
sc.pl.umap(adata_int, color=["celltypist_labels"], save="_celltypist_labels.png", show=False)
sc.pl.umap(adata_int, color=["leiden"],           save="_leiden.png",            show=False)

# ---------- 6) Optional: majority label per cluster & immune split ----------
# Majority label per Leiden cluster
clu_major = (
    adata_int.obs.groupby("leiden")["celltypist_labels"]
    .agg(lambda x: x.value_counts().index[0])
)
adata_int.obs["cluster_identity"] = adata_int.obs["leiden"].map(clu_major)
sc.pl.umap(adata_int, color="cluster_identity", legend_loc="on data",
           save="_cluster_identity.png", show=False)

# Immune vs non-immune split
immune_terms = ["T cell","B cell","Monocyte","Macrophage","NK","DC","Neutrophil"]
immune_mask  = adata_int.obs["celltypist_labels"].str.contains("|".join(immune_terms),
                                                               case=False, na=False)
adata_int[immune_mask].write(os.path.join(out_dir, "adata_celltypist_immune.h5ad"))
adata_int[~immune_mask].write(os.path.join(out_dir, "adata_celltypist_nonimmune.h5ad"))

print("✅ Scanorama + CellTypist pipeline complete.")
