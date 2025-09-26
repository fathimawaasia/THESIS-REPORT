#!/usr/bin/env python3
import argparse, pandas as pd, anndata as ad

def main():
    ap = argparse.ArgumentParser(
        description="Export meta (cell_id, state_label, disease_group, tissue) from h5ad; optionally map state labels from a leiden→state CSV."
    )
    ap.add_argument("h5ad", help="Input .h5ad")
    ap.add_argument("out_tsv", help="Output TSV")
    ap.add_argument("--state-col", default="state_label", help="Obs column with state labels if present (default: state_label)")
    ap.add_argument("--disease-col", default="disease_group", help="Obs column with diagnosis (default: disease_group)")
    ap.add_argument("--tissue-col", default="tissue", help="Obs column with tissue (default: tissue)")

    # mapping options (used only if --state-col missing/empty)
    ap.add_argument("--map", default=None, help="CSV mapping file for leiden→state (e.g. *_state_labels_by_leiden.csv)")
    ap.add_argument("--leiden-col", default="leiden_sub", help="Obs column with leiden ids (default: leiden_sub)")
    ap.add_argument("--map-leiden-col", default="leiden_sub", help="Column name in mapping CSV for leiden ids (default: leiden_sub)")
    ap.add_argument("--map-state-col", default="state", help="Column name in mapping CSV for state labels (default: state)")
    args = ap.parse_args()

    A = ad.read_h5ad(args.h5ad)
    obs = A.obs.copy()

    # ensure string dtypes (avoid pandas ‘category’ surprises)
    for c in [args.state_col, args.disease_col, args.tissue_col, args.leiden_col]:
        if c in obs.columns:
            obs[c] = obs[c].astype(str)

    # choose state labels
    use_mapping = False
    if args.state_col in obs.columns and not obs[args.state_col].isna().all():
        state_series = obs[args.state_col].astype(str)
        source = f"obs['{args.state_col}']"
    else:
        if args.map is None:
            raise SystemExit(
                f"[ERR] '{args.state_col}' not in obs (or empty) and no --map provided.\n"
                f"       Either pass an existing state column with --state-col, or provide a mapping CSV with --map."
            )
        if args.leiden_col not in obs.columns:
            raise SystemExit(f"[ERR] obs is missing '{args.leiden_col}' (needed to merge the mapping). Columns: {list(obs.columns)}")

        map_df = pd.read_csv(args.map)
        # be permissive about header names (strip spaces)
        map_df.columns = [c.strip() for c in map_df.columns]
        if args.map_leiden_col not in map_df.columns or args.map_state_col not in map_df.columns:
            raise SystemExit(
                f"[ERR] Mapping CSV must contain columns '{args.map_leiden_col}' and '{args.map_state_col}'. "
                f"Found: {list(map_df.columns)}"
            )
        map_df = map_df[[args.map_leiden_col, args.map_state_col]].copy()
        map_df[args.map_leiden_col] = map_df[args.map_leiden_col].astype(str)

        # left-join obs by leiden
        tmp = obs[[args.leiden_col]].reset_index().rename(columns={"index":"cell_id"})
        tmp[args.leiden_col] = tmp[args.leiden_col].astype(str)
        merged = tmp.merge(map_df, how="left", left_on=args.leiden_col, right_on=args.map_leiden_col)
        state_series = merged[args.map_state_col].astype(str).set_axis(tmp["cell_id"])
        use_mapping = True
        source = f"mapping CSV ({args.map}) via obs['{args.leiden_col}']"

        # warn on unmapped
        n_unmapped = state_series.isna().sum() + (state_series == "nan").sum()
        if n_unmapped > 0:
            unmapped_levels = sorted(set(tmp.loc[(state_series.isna()) | (state_series == "nan"), args.leiden_col]))
            print(f"[WARN] {n_unmapped} cells could not be mapped to a state. Unmapped leiden ids: {unmapped_levels[:10]}{' ...' if len(unmapped_levels)>10 else ''}")

    # build output dataframe
    out = pd.DataFrame({
        "cell_id": obs.index.astype(str),
        "state_label": state_series.values,
        "disease_group": obs[args.disease_col].astype(str) if args.disease_col in obs.columns else "NA",
        "tissue": obs[args.tissue_col].astype(str) if args.tissue_col in obs.columns else "NA",
    })

    out.to_csv(args.out_tsv, sep="\t", index=False)
    n_states = out["state_label"].nunique(dropna=True)
    print(f"[OK] Wrote: {args.out_tsv}")
    print(f"[INFO] state source: {source}; unique states: {n_states}")
    print("[INFO] head:")
    print(out.head())

if __name__ == "__main__":
    main()


