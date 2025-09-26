# label_common_scatter.py
import os, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---- inputs (edit if your paths differ)
merged_csv = "results/compare_UC_vs_HC/UC_vs_HC_enrichment_merged.csv"
out_png    = "results/compare_UC_vs_HC/scatter_common_terms_UC_vs_HC_labeled.png"
pthresh    = 0.05       # significance threshold used before
top_k      = 10         # how many labels to show

df = pd.read_csv(merged_csv)

# keep terms significant in BOTH and marked as Common (created by compare script)
mask = (df["padj_UC"] < pthresh) & (df["padj_HC"] < pthresh)
if "category" in df.columns:
    mask &= (df["category"] == "Common")

common = df.loc[mask].copy()
if common.empty:
    raise SystemExit("No common significant terms found — check paths/thresholds.")

# compute -log10 p and diagonal distance
common["x"] = -np.log10(common["padj_HC"].clip(lower=1e-300))
common["y"] = -np.log10(common["padj_UC"].clip(lower=1e-300))
common["delta"] = common["y"] - common["x"]           # +ve => stronger in UC
common["d_abs"] = common["delta"].abs()                # distance from diagonal

# choose labels:
# - strongest in UC (above diagonal), then (optionally) strongest in HC
top_uc = common.sort_values("delta", ascending=False).head(top_k)
# If you want some that are stronger in HC too, uncomment next line:
# top_hc = common.sort_values("delta", ascending=True).head(max(0, top_k//3))
labels = top_uc # .append(top_hc).drop_duplicates(subset=["Term"])

# ---- plot
plt.figure(figsize=(6,6))
plt.scatter(common["x"], common["y"], s=10, alpha=0.6)

lims = [0, max(common[["x","y"]].to_numpy().max()+0.5, 5)]
plt.plot(lims, lims, ls="--", lw=1, color="gray")
plt.xlim(lims); plt.ylim(lims)

# annotations
for _, r in labels.iterrows():
    plt.scatter([r["x"]], [r["y"]], s=20)  # re-highlight
    # small offset to reduce overlap
    ha = "left" if r["delta"]>=0 else "right"
    dx = 0.05 if ha=="left" else -0.05
    plt.annotate(
        r["Term"],
        (r["x"], r["y"]),
        xytext=(r["x"]+dx, r["y"]+0.15),
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="black", alpha=0.7),
        arrowprops=dict(arrowstyle="-", lw=0.5, color="black"),
        ha=ha, va="bottom"
    )

plt.title("Common enriched terms (labels = most UC-shifted)")
plt.xlabel("-log10(adjP) HC")
plt.ylabel("-log10(adjP) UC")
plt.tight_layout()
os.makedirs(os.path.dirname(out_png), exist_ok=True)
plt.savefig(out_png, dpi=300)
print(f"[OK] wrote {out_png}")
