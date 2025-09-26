#!/usr/bin/env python3
"""
Builds GMT files of TF regulon targets from pySCENIC ctx output (*.reg.csv).
Handles common column name variants.
"""
import pandas as pd, argparse, os

DEF_DC2 = ["IRF4","RELB","KLF4","NFKB1"]
DEF_DP  = ["AHR","STAT1","IRF7","TCF7","LEF1"]

def parse_reg(reg_path):
    df = pd.read_csv(reg_path)
    cols = {c.lower(): c for c in df.columns}
    # Try to get TF & target columns
    if "tf" in cols and "target" in cols:
        tfc, tgtc = cols["tf"], cols["target"]
        pairs = df[[tfc, tgtc]].dropna().astype(str)
    elif "tf" in cols and "gene" in cols:
        tfc, tgtc = cols["tf"], cols["gene"]
        pairs = df[[tfc, tgtc]].dropna().astype(str)
    elif "regulon" in cols and "genes" in cols:
        # one row per regulon with comma-separated targets
        regs = []
        for _, r in df.iterrows():
            regname = str(r[cols["regulon"]])
            tf = regname.split('(')[0].strip()
            genes = str(r[cols["genes"]]).split(',')
            for g in genes:
                g = g.strip()
                if g:
                    regs.append((tf,g))
        pairs = pd.DataFrame(regs, columns=["TF","target"])
    else:
        raise ValueError("Could not find TF/target columns in .reg.csv")
    pairs.columns = ["TF","target"]
    return pairs

def build_gmt(pairs, tf_list, tag, outg):
    with open(outg, "w") as f:
        for tf in tf_list:
            tgts = sorted(set(pairs.loc[pairs["TF"]==tf, "target"].astype(str)))
            if len(tgts) < 5: 
                print(f"[WARN] {tag}: {tf} has only {len(tgts)} targets; skipping")
                continue
            name = f"{tag}_{tf}_regulon"
            f.write(name + "\tNA\t" + "\t".join(tgts) + "\n")
    print("[OK] wrote", outg)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dc2_reg", default="DC2_all.reg.csv")
    ap.add_argument("--dp_reg",  default="DP_all.reg.csv")
    ap.add_argument("--dc2_tfs", nargs="*", default=DEF_DC2)
    ap.add_argument("--dp_tfs",  nargs="*", default=DEF_DP)
    ap.add_argument("--out", default="REGULON_GMT")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    dc2_pairs = parse_reg(args.dc2_reg)
    dp_pairs  = parse_reg(args.dp_reg)
    build_gmt(dc2_pairs, args.dc2_tfs, "DC2", os.path.join(args.out, "DC2_regulons.gmt"))
    build_gmt(dp_pairs,  args.dp_tfs,  "DP",  os.path.join(args.out, "DP_regulons.gmt"))

if __name__ == "__main__":
    main()
