#!/usr/bin/env python3
"""
choose_scvi_source_and_prepare.py
- Probe multiple candidate .h5ad files to find one with integer-like counts.
- Prefer layers['counts'] > .raw.X > .X.
- Once a source is found, subset genes to the HVGs from a provided HVG meta file,
  copy useful .obs, and write an scVI-ready file with layers['counts'] set.

Usage:
  python3 choose_scvi_source_and_prepare.py \
    --candidates qc_outputs/adata_all_raw_set_with_raw.h5ad \
                 qc_outputs/adata_all_raw_set.h5ad \
                 qc_outputs/adata_all_HVGfiltered.h5ad \
                 qc_outputs/adata_all_postqc_preHVG.h5ad \
                 qc_outputs/adata_preprocessed.pkl\
    --hvg qc_outputs/adata_all_HVGfiltered_meta.h5ad \
    --out qc_outputs/adata_all_HVGfiltered_for_scvi.h5ad
"""
import argparse, sys
from pathlib import Path
import numpy as np
import scanpy as sc
from scipy.sparse import issparse

def is_integer_like_matrix(M) -> bool:
    if issparse(M):
        d = M.data
        if d.size == 0:
            return True
        if np.min(d) < 0:
            return False
        return np.allclose(d, np.round(d), atol=1e-6)
    else:
        arr = np.asarray(M)
        if arr.size == 0 or np.min(arr) < 0:
            return False
        return np.allclose(arr, np.round(arr), atol=1e-6)

def ensure_counts_layer(ad):
    # keep existing counts if present
    if "counts" in ad.layers:
        return True, "existing layers['counts']"
    # prefer .raw.X
    if ad.raw is not None and ad.raw.X is not None and is_integer_like_matrix(ad.raw.X):
        ad.layers["counts"] = ad.raw.X
        return True, ".raw.X → layers['counts']"
    # fallback to .X
    if is_integer_like_matrix(ad.X):
        ad.layers["counts"] = ad.X
        return True, ".X → layers['counts']"
    return False, "no integer-like counts found"

def load_first_scvi_candidate(paths):
    tried = []
    for p in paths:
        P = Path(p)
        tried.append(str(P))
        if not P.exists():
            print(f"[WARN] Missing candidate: {P}", file=sys.stderr)
            continue
        print(f"[INFO] Trying candidate: {P}")
        ad = sc.read(str(P))
        ok, src = ensure_counts_layer(ad)
        print(f"[INFO] {src}")
        if ok:
            print(f"[OK] Using {P} (obs={ad.n_obs}, var={ad.n_vars})")
            return ad, str(P), src
        else:
            print(f"[SKIP] {P} not suitable for scVI (no counts).")
    raise SystemExit(f"[ERROR] None of the candidates had integer-like counts.\nTried:\n- " + "\n- ".join(tried))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", nargs="+", required=True, help="Candidate .h5ad files to probe (order matters).")
    ap.add_argument("--hvg", required=True, help="HVG meta .h5ad to obtain the HVG gene list and optional obs to copy.")
    ap.add_argument("--out", required=True, help="Output scVI-ready .h5ad")
    args = ap.parse_args()

    # 1) load first candidate with counts
    ad_raw, chosen_path, source = load_first_scvi_candidate(args.candidates)

    # 2) load HVG meta to get the HVG set (e.g., 2000 genes) and obs metadata
    hvg_p = Path(args.hvg)
    if not hvg_p.exists():
        raise SystemExit(f"[ERROR] HVG file not found: {hvg_p}")
    print(f"[INFO] Loading HVG meta: {hvg_p}")
    ad_hvg = sc.read(str(hvg_p))
    hvg_genes = ad_hvg.var_names
    print(f"[INFO] HVG genes: {len(hvg_genes)}")

    # 3) subset genes to the intersection of HVGs and candidate
    inter = ad_raw.var_names.intersection(hvg_genes)
    if len(inter) == 0:
        raise SystemExit("[ERROR] No overlap between candidate genes and HVG set.")
    ad = ad_raw[:, inter].copy()
    print(f"[INFO] Subset to HVGs present: {ad.n_vars} genes")

    # counts already ensured in ensure_counts_layer()
    assert "counts" in ad.layers, "counts layer not set unexpectedly."

    # 4) transfer .obs columns from HVG (only for common cells)
    common = ad.obs_names.intersection(ad_hvg.obs_names)
    if len(common) == 0:
        print("[WARN] No overlapping cells between candidate and HVG meta; skipping obs transfer.")
    else:
        if len(common) < ad.n_obs:
            ad = ad[common].copy()
            print(f"[INFO] Restricted to common cells: {ad.n_obs} cells")
        bring_cols = [c for c in ad_hvg.obs.columns if c not in ad.obs.columns]
        if bring_cols:
            ad.obs = ad.obs.join(ad_hvg.obs.loc[ad.obs_names, bring_cols])
            print(f"[INFO] Transferred obs columns: {len(bring_cols)}")

    # 5) write output
    out_p = Path(args.out)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    ad.write(str(out_p), compression="gzip")
    print(f"[DONE] Wrote scVI-ready file: {out_p}")
    print(f"[NOTE] Chosen source: {chosen_path} ({source})")

if __name__ == "__main__":
    main()
