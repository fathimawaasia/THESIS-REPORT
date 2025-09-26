#!/usr/bin/env python3
import scanpy as sc, numpy as np
from scipy.sparse import issparse

path = "qc_outputs/adata_all_postqc_preHVG_GSE125527_counts.h5ad"
ad = sc.read(path)

def is_int_like(M):
    if M is None: return False
    if issparse(M):
        d = M.data
        return d.size == 0 or (d.min() >= 0 and np.allclose(d, np.round(d), atol=1e-8))
    arr = np.asarray(M)
    return arr.size == 0 or (arr.min() >= 0 and np.allclose(arr, np.round(arr), atol=1e-8))

src = None
if "counts" in ad.layers and is_int_like(ad.layers["counts"]):
    src = "layers['counts']"
elif ad.raw is not None and is_int_like(ad.raw.X):
    src = ".raw.X"
elif is_int_like(ad.X):
    src = ".X"

print("File:", path)
print("Cells:", ad.n_obs, "| Genes:", ad.n_vars)
print("Integer-like source:", src)
