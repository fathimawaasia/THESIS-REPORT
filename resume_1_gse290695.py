#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#==============================
# GSE290695_resume_1
#==============================

#==============================
# STEP 01: import libraries////
#==============================

import os
import numpy as np
import pandas as pd
import random as rd
import scanpy as sc
import anndata as ad
import seaborn as sns
import glob
import matplotlib.pyplot as pl
from scipy.sparse import csr_matrix
import scanorama as scrma
import pickle
from scipy.stats import median_abs_deviation
import scrublet as scr
from scipy import sparse
import matplotlib.pyplot as plt

# In[ ]:

# Integration by scanorama run as sbatch with run_scanorama_integration.py
# ========== Step 1: Load preprocessed AnnData ==========
adata_all = sc.read("qc_outputs/adata_all_raw_set.h5ad")
# ========== Step 2: Split AnnData by sample_id ==========
adatas_split = [
    adata_all[adata_all.obs["sample_id"] == sid].copy()
    for sid in adata_all.obs["sample_id"].unique()
]
# ========== Step 3: Run Scanorama integration ==========
import scanpy as sc
import scanorama
corrected_adatas = scanorama.correct_scanpy(adatas_split, return_dimred=True)
# ========== Step 4: Concatenate back into a single AnnData ==========
# Keep track of sample names
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
df_expr.to_csv("scanorama_integrated_expression_matrix.csv")

# Save AnnData object itself
adata_scanorama.write("qc_outputs/adata_scanorama_integrated.h5ad")

print("\n✅ Scanorama integration complete!")
print("Saved integrated expression matrix: scanorama_integrated_expression_matrix.csv")
print("Saved integrated AnnData object: qc_outputs/adata_scanorama_integrated.h5ad")




# In[2]:


# saving expression matrix with header:
import scanpy as sc
import pandas as pd
import numpy as np

# Load the integrated Scanorama AnnData object
adata = sc.read("qc_outputs/adata_scanorama_integrated.h5ad")

# Convert the expression matrix to a dense DataFrame with gene and cell names
df_expr = pd.DataFrame(
    adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X,
    index=adata.obs_names,
    columns=adata.var_names
)

# Save to CSV with row index (cells) and column headers (genes)
df_expr.to_csv("scanorama_integrated_expression_matrix.csv", index=True, header=True)

print("Saved: scanorama_integrated_expression_matrix.csv with row and column labels")



# In[3]:


# Rows are cells and columns are genes.
# Scanorama integrated expression values are -transformed gene expression per cell.
print(df_expr.head())


# In[2]:


import scanpy as sc
adata_scanorama = sc.read("qc_outputs/adata_scanorama_integrated.h5ad")
print(adata_scanorama.obs.columns.tolist())


# In[ ]:


#PCA, UMAP and Leiden clustering for scanorama integrated file was carried on as a sbatch with run_scanorama_umap_clustering.py.

# Load integrated data
adata_scanorama = sc.read("qc_outputs/adata_scanorama_integrated.h5ad")

# ========== Step 1: Assign Scanorama embeddings as PCA ==========
adata_scanorama.obsm["X_pca"] = adata_scanorama.obsm["X_scanorama"]

# PCA plot
sc.pl.pca(
    adata_scanorama,
    color=["Group"],
    legend_loc="right margin",
    legend_fontsize=10,
    save="_scanorama_Group_PCA.png"
)

print("PCA plot saved to figures/_scanorama_Group_PCA.png")

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
    save="_scanorama_integrated_UMAP.png"
)

# Cluster overview
sc.pl.umap(
    adata_scanorama,
    color=["leiden", "Group"],
    legend_loc="right margin",
    legend_fontsize=10,
    ncols=2,
    wspace=0.4,
    save="_scanorama_leiden_Group.png"
)


# In[ ]:


# ========== Step 5: Save clustered object ==========
adata_scanorama.write("qc_outputs/adata_scanorama_clustered_umap.h5ad")

print("UMAP & clustering complete. Plots saved to figures/. AnnData saved to qc_outputs/adata_scanorama_clustered_umap.h5ad") 


# In[ ]:


# Detection and Ranking of Marker genes:
# A: Rank genes by Leiden: Ran as sbatch using run_marker_leiden_detection.py: This is for genes distinguishing clusters.

# B:Rank genes by Group(Condition-UC,UCA,HC): Ran as sbatch using marker_genes_by_group.py. :This gives genes that define condition-specific responses across all cells or per cell type/cluster. Also DEGs between groups.


# A: Rank genes by Leiden: Ran as sbatch using run_marker_leiden_detection.py: This is for genes distinguishing clusters.
# Load clustered data
adata = sc.read("qc_outputs/adata_scanorama_clustered_umap.h5ad")

# Identify marker genes for each Leiden cluster
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
    key_added="rank_genes_leiden"
)

# Create results directory
os.makedirs("results", exist_ok=True)

# Save top marker table
markers_df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_leiden")
markers_df.to_csv("results/scanorama_leiden_markers.csv", index=False)

# Plot top 20 genes per cluster
sc.pl.rank_genes_groups(
    adata,
    n_genes=20,
    sharey=False,
    key="rank_genes_leiden",
    save="_scanorama_leiden_top20.png"
)

print("✅ Marker gene detection complete: results/scanorama_leiden_markers.csv + top20 plot saved")







# In[ ]:


# B:Rank genes by Group(Condition-UC,UCA,HC): Ran as sbatch using marker_genes_by_group.py. :This gives genes that define condition-specific responses across all cells or per cell type/cluster. Also DEGs between groups.
# Load clustered data
adata = sc.read("qc_outputs/adata_scanorama_clustered_umap.h5ad")
# Identify DEGs by condition group (e.g., UC vs CD vs HC)
sc.tl.rank_genes_groups(
    adata,
    groupby="Group",  # must match your adata.obs['Group']
    method="wilcoxon",
    key_added="rank_genes_Group"
)
# Create results folder
os.makedirs("result", exist_ok=True)

# Save all group comparisons
markers_group_df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_Group")
markers_group_df.to_csv("results/condition_DEGs_by_Group.csv", index=False)

# Plot top 20 per group
sc.pl.rank_genes_groups(
    adata,
    n_genes=20,
    sharey=False,
    key="rank_genes_Group",
    save="_top20_Group.png"
)

print(" Marker gene detection by Group complete.")


# In[ ]:


# feature plots
adata_path = "qc_outputs/adata_scanorama_clustered_umap.h5ad"
output_dir = "umap_feature_plots"
os.makedirs(output_dir, exist_ok=True)
# Load AnnData
adata = sc.read(adata_path)
# Set output directory and figure parameters
sc.settings.figdir = output_dir
sc.set_figure_params(dpi=150, format="png")

# Feature genes
genes_to_plot = ["CXCL2", "IL5RA"]

# Plot each gene
for gene in genes_to_plot:
    sc.pl.umap(adata, color=gene, show=False, save=f"_feature_{gene}.png")

print("Feature plots saved in:", output_dir)


# In[ ]:


# Cell type annotation by Cell typist models:
#A: Immune cells annotation by Immune_All_Low.pkl
#B: Nonimmune cells annotation by Cells_Intestinal_Tract.pkl
# This was run as a sbatch with run_celltypist_dual_annotation.py   ####


# # end


