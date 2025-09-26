#!/usr/bin/env python3
# save as export_scenic_tsv.py
import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad

INFILES = [
    ("DC2_subset_reclustered_fullX_states.h5ad", "DC2"),
    ("DP_subset_reclustered_fullX_states.h5ad",  "DP"),
]
OUTDIR = "SCENIC_INPUT"; USE_RAW = True; MIN_CELLS_FRAC = 0.005
os.makedirs(OUTDIR, exist_ok=True)

def boolean_sum_axis0(X):
    return np.asarray((X>0).sum(axis=0)).ravel() if hasattr(X,"tocsr") else (X>0).sum(axis=0)

for path, tag in INFILES:
    A = sc.read_h5ad(path)
    Ar = ad.AnnData(X=A.raw.X, var=A.raw.var.copy(), obs=A.obs.copy()) if (USE_RAW and A.raw is not None) else A.copy()
    Ar.var_names_make_unique()
    if MIN_CELLS_FRAC>0:
        min_cells = int(np.ceil(MIN_CELLS_FRAC*Ar.n_obs)); keep = boolean_sum_axis0(Ar.X)>=min_cells
        Ar = Ar[:, keep].copy()

    X = Ar.X.toarray() if hasattr(Ar.X,"toarray") else Ar.X  # make dense for TSV
    expr = os.path.join(OUTDIR, f"{tag}_expr.tsv.gz")
    genes= os.path.join(OUTDIR, f"{tag}_genes.txt")
    cells= os.path.join(OUTDIR, f"{tag}_cells.txt")

    print(f"[WRITE] {expr}  (cells={Ar.n_obs:,}, genes={Ar.n_vars:,})")
    pd.DataFrame(X, index=Ar.obs_names, columns=Ar.var_names).to_csv(expr, sep="\t", compression="gzip")
    pd.Series(Ar.var_names).to_csv(genes, index=False, header=False)
    pd.Series(Ar.obs_names).to_csv(cells, index=False, header=False)
