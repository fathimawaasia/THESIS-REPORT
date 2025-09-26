import scanpy as sc, pandas as pd, json
H5AD = "adata_scanvi_withRAW.h5ad"
GMT = "metabolic_sets.gmt"

# read
adata = sc.read_h5ad(H5AD)
with open(GMT) as f:
    lines = [l.strip().split('\t') for l in f]

rows = []
for l in lines:
    set_name = l[0]
    genes = set(l[2:])
    present = [g for g in genes if g in adata.var_names]
    rows.append({
        "set": set_name,
        "n_genes_total": len(genes),
        "n_in_adata": len(present),
        "frac_in_adata": round(len(present)/max(1,len(genes)), 3),
        "example_present": ",".join(present[:10])
    })

df = pd.DataFrame(rows).sort_values("set")
df.to_csv("metabolic_set_coverage.csv", index=False)
print(df)
