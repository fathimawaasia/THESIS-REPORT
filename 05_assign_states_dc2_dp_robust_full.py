#!/usr/bin/env python3
"""
Robust state-labelling for DC2 & DP using full-gene .X subsets.
- Validates signature coverage (keeps if >= COV_MIN_FRACTION genes present)
- Scores modules on full X, z-scores, penalizes stress signatures slightly
- Assigns labels by per-Leiden median winner
Outputs -> QC_SUMMARY_FULLX/CLUSTER_STATES/*
"""

import os, numpy as np, pandas as pd, scanpy as sc
import matplotlib.pyplot as plt

IN = {
    "DC2": "DC2_subset_reclustered_fullX.h5ad",
    "DP":  "DP_subset_reclustered_fullX.h5ad",
}
OUTDIR = "QC_SUMMARY_FULLX/CLUSTER_STATES"
os.makedirs(OUTDIR, exist_ok=True)

# ---- signatures (gut/IBD-aware) ----
DC2_SIGS = {
  "classical DC2 (IRF4/KLF4/CD1C)": ["CD1C","FCER1A","CLEC10A","IRF4","KLF4","ITGAX","CSF2RA","CSF1R","HLA-DRA","HLA-DPA1","HLA-DPB1"],
  "activated/APC DC2 (CD83/CD40+)": ["CD83","CD40","CCR7","LAMP3","FSCN1","RELB","NFKBIA","CCL19","CCL22","BIRC3","ICAM1"],
  "migratory DC2 (CCR7+/LAMP3+)":   ["CCR7","LAMP3","FSCN1","MARCKSL1","CCL19","CCL21","CCL22","ITGAE"],
  "inflammatory DC2 (TNF/IL1/NFKB)": ["TNF","IL1B","IL6","NFKBIA","REL","RELB","JUN","FOS","PTGS2","CXCL8","CCL20"],
  "interferon-stimulated DC2":      ["ISG15","IFIT1","IFIT2","IFIT3","IFITM1","MX1","OAS1","OASL","RSAD2","IRF7","STAT1"],
  "tolerogenic/regulatory DC2":     ["IL10","IDO1","LGALS9","LILRB1","LILRB2","PDCD1LG2","AHR","ALDH1A2"],
  "cycling DC2":                    ["MKI67","TOP2A","PCNA","UBE2C","CENPF","TYMS","BIRC5","STMN1"],
  "stress/metal":                   ["MT1X","MT1E","MT1F","MT2A","HSPA1A","HSPH1","DDIT3","ATF4","FTL","FTH1"],
}

DP_SIGS = {
  "canonical DP thymocytes":        ["LEF1","TCF7","RAG1","RAG2","PTCRA","CCR9","S1PR1","CD3D","CD3E","CD3G","BCL11B"],
  "cycling DP":                     ["MKI67","TOP2A","PCNA","UBE2C","CENPA","CENPE","TK1","TYMS"],
  "IFN/inflammatory DP":            ["ISG15","IFITM1","IFIT1","MX1","OAS1","CXCL10","IRF7","STAT1"],
  "migratory/egress DP":            ["KLF2","CCR7","S1PR1","SELL","ITGA4","TCF7"],
  "metabolic/oxidative (AHR/mito)": ["AHR","AHRR","CYP1B1","IDO1","PDK1","PDK3","ACADVL","NDUFAF6","COX7A2L"],
  "stress/metal DP":                ["MT1X","MT1E","MT1F","MT2A","HSPA1A","HSPH1","DDIT3","ATF4"],
}

S_GENES   = ["MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","TOP2A","MKI67"]
G2M_GENES = ["HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","AURKB"]

# ---- knobs ----
COV_MIN_FRACTION = 0.40   # keep a signature only if >=40% genes exist in .var_names
STRESS_PENALTY   = 0.15   # subtract 0.15 * stress score from all other program medians
# ----------------

def _validate_signatures(ad, sigs: dict, label_prefix: str):
    """Return (kept_dict, coverage_df)."""
    rows = []
    kept = {}
    genes_set = set(ad.var_names)
    for label, genes in sigs.items():
        present = [g for g in genes if g in genes_set]
        frac = len(present) / max(1, len(genes))
        rows.append({"program": f"{label_prefix}__{label}", "n_total": len(genes), "n_present": len(present), "fraction_present": frac})
        if frac >= COV_MIN_FRACTION and len(present) >= 2:
            kept[label] = present
    cov = pd.DataFrame(rows).sort_values("program")
    return kept, cov

def _score_modules(ad, sigs_kept: dict, prefix: str):
    """Score each signature on full .X, then z-score."""
    for label, genes in sigs_kept.items():
        score_col = f"{prefix}__{label}"
        sc.tl.score_genes(ad, genes, score_name=score_col, use_raw=False)
        s = ad.obs[score_col]
        ad.obs[score_col] = (s - s.mean()) / (s.std(ddof=0) + 1e-8)

def _label_by_cluster_medians(ad, prefix: str, cluster_key: str, penalty_key: str | None = None):
    cols = [c for c in ad.obs.columns if c.startswith(prefix + "__")]
    df = ad.obs[[cluster_key] + cols].copy()
    med = df.groupby(cluster_key, observed=False).median(numeric_only=True)

    if penalty_key:
        pk = f"{prefix}__{penalty_key}"
        if pk in med.columns:
            # apply small penalty to *other* programs
            penal = med[pk] * STRESS_PENALTY
            for col in med.columns:
                if col != pk:
                    med[col] = med[col] - penal

    winners = med.idxmax(axis=1).str.replace(prefix+"__", "", regex=False)
    return med, winners

def _plot_states(ad, winners, title, out_png):
    tmp = ad.obs["leiden_sub"].map(winners).astype("category")
    ad.obs["_state_label_tmp"] = tmp
    fig, ax = plt.subplots(figsize=(9,7))
    sc.pl.umap(ad, color="_state_label_tmp", legend_loc="right", frameon=False,
               title=title, ax=ax, show=False)
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    del ad.obs["_state_label_tmp"]

def _write_by_tissue(ad, winners, out_csv):
    lab = ad.obs["leiden_sub"].map(winners).rename("state")
    tcol = "tissue" if "tissue" in ad.obs.columns else None
    if tcol is None:
        return
    tab = ad.obs.assign(state=lab).groupby([tcol, "state"], observed=False).size().unstack(fill_value=0)
    tab.to_csv(out_csv)

def process_one(tag: str, in_path: str, sigs: dict, prefix: str, stress_label_for_penalty: str):
    if not os.path.exists(in_path):
        print(f"[SKIP] {tag} file missing: {in_path}"); return
    ad = sc.read(in_path)

    # 1) validate signatures against full var_names
    kept, cov = _validate_signatures(ad, sigs, prefix)
    cov_path = os.path.join(OUTDIR, f"{tag}_signature_coverage.csv")
    cov.to_csv(cov_path, index=False)
    kept_names = list(kept.keys())
    dropped = sorted(set(sigs.keys()) - set(kept_names))
    print(f"[{tag}] kept {len(kept)} programs: {kept_names}")
    if dropped:
        print(f"[{tag}] dropped (low coverage): {dropped}")

    # 2) score modules + cell cycle on full X
    _score_modules(ad, kept, prefix)
    s_ok = [g for g in S_GENES if g in ad.var_names]
    g_ok = [g for g in G2M_GENES if g in ad.var_names]
    if s_ok and g_ok:
        sc.tl.score_genes_cell_cycle(ad, s_genes=s_ok, g2m_genes=g_ok, use_raw=False)
        cc = (ad.obs["S_score"] + ad.obs["G2M_score"])
        ad.obs[f"{prefix}__cell_cycle"] = (cc - cc.mean()) / (cc.std(ddof=0) + 1e-8)

    # 3) assign labels per Leiden median (with stress penalty)
    med, winners = _label_by_cluster_medians(ad, prefix, "leiden_sub", penalty_key=stress_label_for_penalty)

    med.to_csv(os.path.join(OUTDIR, f"{tag}_state_medians_by_leiden.csv"))
    winners.rename("state").to_csv(os.path.join(OUTDIR, f"{tag}_state_labels_by_leiden.csv"))

    # 4) plots + tissue tables
    _plot_states(ad, winners, f"{tag} subclusters (state-labeled, fullX)", 
                 os.path.join(OUTDIR, f"umap_{tag}_states_fullX.png"))
    _write_by_tissue(ad, winners, os.path.join(OUTDIR, f"{tag}_states_by_tissue.csv"))
    print(f"[OK] {tag} states assigned (fullX)")

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    process_one("DC2", IN["DC2"], DC2_SIGS, "dc2", "stress/metal")
    process_one("DP",  IN["DP"],  DP_SIGS,  "dp",  "stress/metal DP")
    print("[DONE]")

if __name__ == "__main__":
    main()
