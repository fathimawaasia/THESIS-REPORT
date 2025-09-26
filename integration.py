#!/usr/bin/env python
# coding: utf-8

"""
GSE125527 Scanorama Integration and Clustering
- Input: HVG-filtered, scaled, metadata-attached h5ad files
- Output: Scanorama-integrated, clustered AnnData with PCA, UMAP, and Leiden plots
"""
import os
import scanpy as sc
import pandas as pd
import numpy as np
import scanorama

# ========== Paths ==========
qc_dir = "qc_outputs"
input_file = f"{qc_dir}/adata_all_HVGfiltered_meta.h5ad"
os.makedirs(qc_dir, exist_ok=True)

# ========== Step 1: Load data ==========
print("▶ Loading HVG filtered and scaled h5ad file...")
adata = sc.read(input_file)

# ========== Step 2: Split by sample for Scanorama ==========
print("▶ Splitting AnnData by sample_id for Scanorama...")
adata_list = [adata[adata.obs['sample_id'] == sid].copy() for sid in adata.obs['sample_id'].unique()]

# ========== Step 3: Run Scanorama integration ==========
print("▶ Running Scanorama integration...")
scanorama.integrate_scanpy(adata_list, dimred=50)

# Concatenate integrated data (using anndata.concat instead of deprecated concatenate)
print("▶ Concatenating Scanorama-integrated datasets...")
adata_scanorama = adata_list[0].concatenate(
    *adata_list[1:],
    batch_key="sample_id",
    index_unique="-",
    join="inner"
)

# Store integrated embedding
adata_scanorama.obsm["X_scanorama"] = np.concatenate([ad.obsm["X_scanorama"] for ad in adata_list])

# ========== Step 4: Save integrated file and expression matrix ==========
data_path = f"{qc_dir}/adata_all_scanorama_integrated.h5ad"
print(f"▶ Saving integrated h5ad to {data_path}...")
adata_scanorama.write(data_path)

# Save expression matrix as CSV
print("▶ Saving expression matrix to CSV...")
X_array = adata_scanorama.X if isinstance(adata_scanorama.X, np.ndarray) else adata_scanorama.X.toarray()
expr_matrix = pd.DataFrame(X_array, index=adata_scanorama.obs_names, columns=adata_scanorama.var_names)
expr_matrix.to_csv(f"{qc_dir}/adata_all_scanorama_expression.csv")

# ========== Step 5: PCA/Neighbors/UMAP/Leiden ==========
print("▶ Performing PCA, UMAP, and Leiden clustering on Scanorama embedding...")
sc.tl.pca(adata_scanorama, svd_solver="arpack")
sc.pl.pca(adata_scanorama, color="sample_id", save="_GSE125527_PCA_scanorama_sample.png")
if "Group" in adata_scanorama.obs.columns:
    sc.pl.pca(adata_scanorama, color="Group", save="_GSE125527_PCA_scanorama_group.png")
sc.pp.neighbors(adata_scanorama, use_rep="X_scanorama")
sc.tl.umap(adata_scanorama)
sc.tl.leiden(adata_scanorama, resolution=0.5)

# Save final clustered UMAP h5ad
adata_scanorama.write(f"{qc_dir}/adata_all_scanorama_clustered_umap.h5ad")

# ========== Step 6: Plots ==========
print("▶ Generating UMAP and clustering plots...")
sc.pl.umap(adata_scanorama, color=["sample_id"], save="_GSE125527_UMAP_scanorama_sample.png")
if "Group" in adata_scanorama.obs.columns:
    sc.pl.umap(adata_scanorama, color=["Group"], save="_GSE125527_UMAP_scanorama_group.png")
sc.pl.umap(adata_scanorama, color="leiden", save="_GSE125527_UMAP_scanorama_leiden.png")

print("✅ Scanorama integration, clustering, and expression export complete. All outputs saved in qc_outputs/")
