#!/usr/bin/env python3
import argparse, os, sys, re, hashlib, glob
import pandas as pd

def first_header(path):
    with open(path, 'r') as fh:
        return fh.readline().rstrip('\n')

def read_header_list(path):
    hdr = first_header(path).split('\t')
    if len(hdr) < 2:
        raise ValueError(f"Header in {path} doesn’t look like TSV.")
    return hdr

def get_cell_types(signature_tsv):
    hdr = read_header_list(signature_tsv)
    # drop first col (GeneSymbol)
    return hdr[1:]

def get_sample_ids(mixture_tsv):
    hdr = read_header_list(mixture_tsv)
    # drop first col (GeneSymbol)
    return hdr[1:]

def detect_sep_and_read(path):
    """Try TSV first; if it collapses to 1 column, fall back to whitespace."""
    try:
        df = pd.read_csv(path, sep='\t', dtype=str)
        if df.shape[1] == 1:
            raise ValueError("Looks not TSV; retry whitespace")
        # cast numeric where possible
        return df.apply(pd.to_numeric, errors='ignore')
    except Exception:
        df = pd.read_csv(path, sep=r"\s+", engine="python", header=None, dtype=str)
        return df.apply(pd.to_numeric, errors='ignore')

def normalize_colnames(cols):
    def norm(x):
        s = str(x).strip()
        s = s.replace('.', ' ')
        return s
    return [norm(c) for c in cols]

def choose_result_file(result_path, base_dir):
    if result_path and os.path.exists(result_path):
        return result_path
    # otherwise try to auto-pick something sensible from base_dir/cibersortx
    cand = sorted(glob.glob(os.path.join(base_dir, "cibersortx", "*result*.txt"))) \
        + sorted(glob.glob(os.path.join(base_dir, "cibersortx", "CIBERSORTx*.txt"))) \
        + sorted(glob.glob(os.path.join(base_dir, "cibersortx", "*.tsv")))
    if not cand:
        raise FileNotFoundError("Couldn’t find a CIBERSORTx result file. Specify --results.")
    return cand[-1]

def main():
    ap = argparse.ArgumentParser(description="Map CIBERSORTx results to GSM IDs and write wide + tidy outputs.")
    ap.add_argument("--base", default=".", help="Working directory (default: .)")
    ap.add_argument("--signature", default="cibersortx/cibersortx_signature.txt")
    ap.add_argument("--mixture",   default="cibersortx/cibersortx_mixture.txt")
    ap.add_argument("--results",   default="", help="CIBERSORTx result file (.txt/.tsv). If omitted I’ll try to detect.")
    ap.add_argument("--outdir",    default="cibersortx_out")
    args = ap.parse_args()

    base = os.path.abspath(args.base)
    sig  = os.path.join(base, args.signature)
    mix  = os.path.join(base, args.mixture)
    os.makedirs(os.path.join(base, args.outdir), exist_ok=True)

    cell_types  = get_cell_types(sig)
    sample_ids  = get_sample_ids(mix)

    res_path = choose_result_file(args.results, base)
    res = detect_sep_and_read(res_path)

    # Case A: official CIBERSORTx table with header
    # Usually has "Mixture" + each cell type + stats.
    stats_known = ["P-value", "P value", "Pvalue", "P value (perm)",
                   "Correlation", "Pearson Correlation", "Pearson_correlation",
                   "RMSE", "Root Mean Squared Error", "Absolute score", "Absolute_score"]
    if list(res.columns)[0] == "Mixture" or "Mixture" in res.columns:
        res.columns = normalize_colnames(res.columns)
        # harmonize stat names
        col_map = {}
        for c in res.columns:
            cn = c.strip()
            if cn.lower().startswith("mixture"):
                col_map[c] = "Mixture"
            elif cn.lower().startswith("p-value") or cn.lower().startswith("p value") or cn.lower().startswith("pvalue"):
                col_map[c] = "P-value"
            elif "correlation" in cn.lower():
                col_map[c] = "Correlation"
            elif "rmse" in cn.lower():
                col_map[c] = "RMSE"
            elif "absolute" in cn.lower():
                col_map[c] = "Absolute score"
        res = res.rename(columns=col_map)
        # Ensure all cell types exist as columns (some builds change spacing)
        # If not, try exact names from signature
        missing = [ct for ct in cell_types if ct not in res.columns]
        if missing:
            raise ValueError(f"Missing cell-type columns in results: {missing}")
    else:
        # Case B: raw numeric matrix without header.
        # Expect len(cell_types) + 3 or +4 numeric cols.
        n = res.shape[1]
        if n not in (len(cell_types)+3, len(cell_types)+4):
            raise ValueError(f"Unexpected columns ({n}) in raw results; expected "
                             f"{len(cell_types)+3} or {len(cell_types)+4}.")
        cols = cell_types + ["P-value", "Correlation", "RMSE"]
        if n == len(cell_types)+4:
            cols = cell_types + ["P-value", "Correlation", "RMSE", "Absolute score"]
        res.columns = cols
        if len(sample_ids) != res.shape[0]:
            raise ValueError(f"Rows in results ({res.shape[0]}) != samples in mixture header ({len(sample_ids)}).")
        res.insert(0, "Mixture", sample_ids)

    # Order columns neatly: Mixture, cell types…, stats (if present)
    stat_cols = [c for c in ["P-value", "Correlation", "RMSE", "Absolute score"] if c in res.columns]
    res = res[["Mixture"] + cell_types + stat_cols]

    # Save wide and tidy
    out_wide = os.path.join(base, args.outdir, "cibersortx_fractions_wide.tsv")
    res.to_csv(out_wide, sep="\t", index=False)

    id_vars = ["Mixture"] + stat_cols
    tidy = res.melt(id_vars=id_vars, value_vars=cell_types,
                    var_name="cell_type", value_name="fraction")
    tidy = tidy.rename(columns={"Mixture": "sample"})
    out_tidy = os.path.join(base, args.outdir, "cibersortx_fractions_tidy.csv")
    tidy.to_csv(out_tidy, index=False)

    # small provenance summary
    with open(os.path.join(base, args.outdir, "_summary.txt"), "w") as fh:
        fh.write(f"Result file: {res_path}\n")
        fh.write(f"samples={len(tidy['sample'].unique())}  cell_types={len(cell_types)}  rows_tidy={len(tidy)}\n")
        fh.write(f"wide: {out_wide}\n")
        fh.write(f"tidy: {out_tidy}\n")

    print(f"[OK] Wrote:\n  {out_wide}\n  {out_tidy}")

if __name__ == "__main__":
    main()
