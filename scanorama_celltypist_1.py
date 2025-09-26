# scanorama_celltypist_pipeline.py

import scanpy as sc
import celltypist
import pandas as pd
import numpy as np
import os
import scanorama

# ========== Paths ==========

# fix_normalization.py
import scanpy as sc

# Load the raw HVG-filtered file
adata = sc.read("qc_outputs/adata_all_raw_set.h5ad")

# Normalize and log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Now save .raw correctly
adata.raw = adata.copy()

# Save fixed version
adata.write("qc_outputs/adata_all_raw_set_normalized.h5ad")

print("✅ Saved log1p-normalized data with .raw attached.")

input_file = "qc_outputs/adata_all_raw_set_normalized.h5ad"
output_dir = "qc_outputs"
figures_dir = "figures"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(figures_dir, exist_ok=True)

adata_all = sc.read("qc_outputs/adata_all_raw_set_normalized.h5ad")

# ========== Step 1: Load and attach .raw ==========
print("🔹 Loading raw HVG-filtered log-normalized data...")
adata_all = sc.read(input_file)
adata_all.raw = adata_all.copy()

# Save a version with .raw for reproducibility
adata_all.write(os.path.join(output_dir, "adata_all_raw_set_with_raw.h5ad"))
# ======== Step 1b: PCA/UMAP/Leiden before Scanorama =========
print("🔹 Running PCA, UMAP, Leiden BEFORE Scanorama...")

# Run PCA
sc.tl.pca(adata_all, n_comps=50, svd_solver='arpack')

# Neighbors, UMAP, Leiden clustering
sc.pp.neighbors(adata_all, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata_all)
sc.tl.leiden(adata_all, resolution=0.5, key_added="leiden_pre")

# Save pre-Scanorama plots
sc.pl.pca(adata_all, save="_pre_scanorama_pca_1.png", show=False)
sc.pl.umap(adata_all, color=["sample_id"], save="_pre_scanorama_umap_by_sample_1.png", show=False)
sc.pl.umap(adata_all, color=["leiden_pre"], legend_loc="on data",
           save="_pre_scanorama_umap_by_leiden_1.png", show=False)


# ========== Step 2: Split by sample_id and preserve .raw ==========
adatas_split = []
raws_split = []

for sid in adata_all.obs["sample_id"].unique():
    ad = adata_all[adata_all.obs["sample_id"] == sid].copy()
    adatas_split.append(ad)
    raws_split.append(ad.copy())  # separately keep raw expression for later

# ========== Step 3: Run Scanorama integration ==========
print("🔹 Running Scanorama integration...")
corrected_adatas = scanorama.correct_scanpy(adatas_split, return_dimred=True)

# Reattach .raw to each corrected AnnData
for i in range(len(corrected_adatas)):
    corrected_adatas[i].raw = raws_split[i]

# ========== Step 4: Concatenate integrated AnnDatas ==========
print("🔹 Concatenating corrected AnnData objects...")
sample_ids = [a.obs["sample_id"][0] for a in corrected_adatas]
adata_scanorama = sc.concat(corrected_adatas, label="sample_id", keys=sample_ids)
adata_scanorama.obs_names_make_unique(join="_")


# ========== Step 5: Save integrated expression matrix ==========
# Convert expression matrix to dense format if needed
df_expr = pd.DataFrame(
    adata_scanorama.X.toarray() if not isinstance(adata_scanorama.X, np.ndarray) else adata_scanorama.X,
    index=adata_scanorama.obs_names,
    columns=adata_scanorama.var_names
)
df_expr.to_csv("scanorama_integrated_expression_matrix_1.csv")

# Save AnnData object itself
adata_scanorama.write("qc_outputs/adata_scanorama_integrated_1.h5ad")

print("\n Scanorama integration complete!")
print("Saved integrated expression matrix: scanorama_integrated_expression_matrix_1.csv")
print("Saved integrated AnnData object: qc_outputs/adata_scanorama_integrated_1.h5ad")

# saving expression matrix with header:
import scanpy as sc
import pandas as pd
import numpy as np

# Load the integrated Scanorama AnnData object
adata = sc.read("qc_outputs/adata_scanorama_integrated_1.h5ad")

# Convert the expression matrix to a dense DataFrame with gene and cell names
df_expr = pd.DataFrame(
    adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X,
    index=adata.obs_names,
    columns=adata.var_names
)

# Save to CSV with row index (cells) and column headers (genes)
df_expr.to_csv("scanorama_integrated_expression_matrix_1.csv", index=True, header=True)

print("Saved: scanorama_integrated_expression_matrix_1.csv with row and column labels")

# Rows are cells and columns are genes.
# Scanorama integrated expression values are -transformed gene expression per cell.
print(df_expr.head())

adata_scanorama = sc.read("qc_outputs/adata_scanorama_integrated_1.h5ad")
print(adata_scanorama.obs.columns.tolist())

# ======== Step 1b: PCA/UMAP/Leiden After Scanorama =========
print("🔹 Running PCA, UMAP, Leiden AFTER Scanorama...")

import scanpy as sc
import matplotlib.pyplot as plt

# Load integrated data
adata_scanorama = sc.read("qc_outputs/adata_scanorama_integrated_1.h5ad")

# ========== Step 1: Assign Scanorama embeddings as PCA ==========
adata_scanorama.obsm["X_pca"] = adata_scanorama.obsm["X_scanorama"]

# PCA plot
sc.pl.pca(
    adata_scanorama,
    color=["Group"],
    legend_loc="right margin",
    legend_fontsize=10,
    save="_scanorama_Group_PCA_1.png"
)

print("PCA plot saved to figures/_scanorama_Group_PCA_1.png")


# ========== Step 2: Neighbors, UMAP ==========
sc.pp.neighbors(adata_scanorama, use_rep="X_pca", n_neighbors=15, n_pcs=40)
sc.tl.umap(adata_scanorama)

# ========== Step 3: Leiden Clustering ==========
sc.tl.leiden(adata_scanorama, resolution=0.5, flavor="igraph", directed=False)

# ========== Step 4: UMAP Plotting ==========
sc.set_figure_params(dpi=150)

# Grouping overview
sc.pl.umap(
    adata_scanorama,
    color=["Group", "sample_id", "leiden", "GSM_ID"],
    legend_loc="right margin",
    legend_fontsize=10,
    ncols=2,
    wspace=0.4,
    save="_scanorama_integrated_UMAP_1.png"
)

# Cluster overview
sc.pl.umap(
    adata_scanorama,
    color=["leiden", "Group"],
    legend_loc="right margin",
    legend_fontsize=10,
    ncols=2,
    wspace=0.4,
    save="_scanorama_leiden_Group_1.png"
)

# ========== Step 5: Save clustered object ==========
adata_scanorama.write("qc_outputs/adata_scanorama_clustered_umap_1.h5ad")

print("UMAP & clustering complete. Plots saved to figures/. AnnData saved to qc_outputs/adata_scanorama_clustered_umap_1.h5ad") 

# ========== Step 6: Restore expression from .raw for CellTypist ==========
print("🔹 Restoring .raw for CellTypist...")
adata_scanorama.X = adata_scanorama.raw.X.copy()
adata_scanorama.var = adata_scanorama.raw.var.copy()

# ========== Step 7: Run CellTypist ==========
print("🔹 Running CellTypist annotation...")
celltypist.models.download_models()  # Downloads all standard models
model = 'Immune_All_Low.pkl'  # Use the model name as a string

predictions = celltypist.annotate(adata_scanorama, model=model, majority_voting=True)
# ========= Step X: Run CellTypist =========
print("🔹 Running CellTypist annotation...")
predictions = celltypist.annotate(adata_scanorama, model=model, majority_voting=True)

# Extract single label column from predictions
pl = predictions.predicted_labels
if isinstance(pl, pd.DataFrame):
    pl = pl.iloc[:, 0]  # first column contains the label strings
adata_scanorama.obs["celltypist_labels"] = pl.astype(str).values

# Majority voting labels
if hasattr(predictions, "majority_voting") and predictions.majority_voting is not None:
    mv = predictions.majority_voting
    if isinstance(mv, pd.DataFrame):
        mv = mv.iloc[:, 0]
    adata_scanorama.obs["celltypist_labels_mv"] = mv.astype(str).values

# Optional: add prediction confidence
if hasattr(predictions, "probability") and predictions.probability is not None:
    adata_scanorama.obs["celltypist_confidence"] = predictions.probability.max(axis=1).values

# ========== Step 8: Save annotated data and labels ==========
annot_h5ad = os.path.join(output_dir, "adata_celltypist_annotated.h5ad")
annot_csv = os.path.join(output_dir, "celltypist_labels_all.csv")
adata_scanorama.write(annot_h5ad)
adata_scanorama.obs[['celltypist_labels']].to_csv(annot_csv)

# ========== Step 9: Plot UMAPs ==========
print("🔹 Plotting UMAPs...")
sc.pl.umap(adata_scanorama, color=['celltypist_labels'], save="_celltypist_labels.png", show=False)
sc.pl.umap(adata_scanorama, color=['leiden'], save="_leiden_clusters.png", show=False)

# ======== Step 9c: Clean UMAP colored by majority CellTypist label per cluster ========
import matplotlib.pyplot as plt

# 1) Get majority CellTypist label per Leiden cluster
cluster_labels = (
    adata_scanorama.obs.groupby("leiden")["celltypist_labels"]
    .agg(lambda x: x.value_counts().index[0])
)

# 2) Map each cell to its cluster's majority label
adata_scanorama.obs["cluster_identity"] = adata_scanorama.obs["leiden"].map(cluster_labels)

# 3) Plot UMAP colored by cluster identity
sc.pl.umap(
    adata_scanorama,
    color="cluster_identity",
    legend_loc="on data",       # puts labels on the map
    legend_fontsize=8,
    size=20,
    save="_cluster_identity_umap.png",
    show=False
)

# 4) Also overlay text manually at cluster centroids (optional, for clarity)
coords   = adata_scanorama.obsm["X_umap"]
clusters = adata_scanorama.obs["leiden"].astype(str).values
centers  = (
    pd.DataFrame(coords, columns=["x", "y"])
      .assign(cluster=clusters)
      .groupby("cluster")[["x", "y"]]
      .median()
)

plt.figure()
sc.pl.umap(adata_scanorama, color="cluster_identity", show=False)
ax = plt.gca()
for cl, row in centers.iterrows():
    label = cluster_labels.loc[cl]
    ax.text(
        row["x"], row["y"], label,
        ha="center", va="center",
        fontsize=8, weight="bold",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", alpha=0.7)
    )

plt.savefig("figures/umap_cluster_identity_with_labels.png", dpi=300, bbox_inches="tight")
plt.close()

print("✅ Saved: figures/umap_cluster_identity_umap.png and umap_cluster_identity_with_labels.png")


# ========== Step 10: Split immune vs non-immune ==========
print("🔹 Splitting immune vs non-immune cells...")
immune_keywords = ["T cell", "B cell", "Monocyte", "Macrophage", "NK", "DC", "Neutrophil"]
immune_mask = adata_scanorama.obs['celltypist_labels'].str.contains('|'.join(immune_keywords), case=False, na=False)

adata_immune = adata_scanorama[immune_mask].copy()
adata_nonimmune = adata_scanorama[~immune_mask].copy()

adata_immune.write(os.path.join(output_dir, "adata_celltypist_immune.h5ad"))
adata_nonimmune.write(os.path.join(output_dir, "adata_celltypist_nonimmune.h5ad"))

adata_immune.obs[['celltypist_labels']].to_csv(os.path.join(output_dir, "celltypist_labels_immune.csv"))
adata_nonimmune.obs[['celltypist_labels']].to_csv(os.path.join(output_dir, "celltypist_labels_nonimmune.csv"))

print("✅ Full Scanorama + CellTypist pipeline complete.")