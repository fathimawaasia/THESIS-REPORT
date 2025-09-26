#!/usr/bin/env python3
"""
Fetch GSM metadata for one or more GEO Series (GSE*) and filter to biopsy samples
from colon / ileum / rectum (or any tissues you pass via --tissues).

Adds optional splits:
- by site
- by disease
- by inflammation
- cross (site ∩ disease)

Requires: GEOparse  (pip install GEOparse)
"""

import argparse, re, sys
from pathlib import Path
try:
    import GEOparse
except ImportError:
    sys.exit("Please install GEOparse first:  pip install GEOparse")

def norm(x): return (x or "").strip()
def stringify(v): return "; ".join(map(str,v)) if isinstance(v,(list,tuple)) else str(v or "")

def parse_characteristics(meta):
    out = {}
    for k,v in meta.items():
        if k.lower().startswith("characteristics_ch1"):
            vals = v if isinstance(v,list) else [v]
            for entry in vals:
                s = stringify(entry)
                if ":" in s: kk,vv=s.split(":",1)
                elif "=" in s: kk,vv=s.split("=",1)
                else: kk,vv="characteristics",s
                kk=kk.strip().lower().replace(" ","_"); vv=vv.strip()
                out[kk]=out.get(kk,""); out[kk]+=(";"+vv if out[kk] else vv)
    return out

def build_row(gsm):
    meta=gsm.metadata or {}; ch=parse_characteristics(meta)
    row={"GSM":gsm.name,"title":norm(meta.get("title",[""])[0] if "title" in meta else ""),
         "tissue":"","diagnosis":"","inflammation":"","other_characteristics":""}
    for key in ["tissue","biopsy_site","anatomical_site","site","body_site","location"]:
        if key in ch and not row["tissue"]: row["tissue"]=ch[key].lower()
    for key in ["diagnosis","disease","phenotype","condition","group"]:
        if key in ch and not row["diagnosis"]: row["diagnosis"]=ch[key].lower()
    for key in ["inflammation_status","activity","status"]:
        if key in ch and not row["inflammation"]: row["inflammation"]=ch[key].lower()
    if ch: row["other_characteristics"]="; ".join([f"{k}:{v}" for k,v in ch.items()])
    return row

def write_csv(rows,path):
    import csv; Path(path).parent.mkdir(parents=True,exist_ok=True)
    fields=["GSM","title","tissue","diagnosis","inflammation","other_characteristics"]
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for r in rows: w.writerow({k:r.get(k,"") for k in fields})
    print(f"[WRITE] {path}  (n={len(rows)})")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gse",nargs="+",required=True)
    ap.add_argument("--tissues",nargs="+",default=["colon","ileum","rectum"])
    ap.add_argument("--outdir",default="geo_filtered")
    ap.add_argument("--split-by-site",action="store_true")
    ap.add_argument("--split-by-disease",action="store_true")
    ap.add_argument("--split-by-inflammation",action="store_true")
    ap.add_argument("--split-cross",action="store_true",
                   help="Write per-site ∩ per-disease CSVs")
    ap.add_argument("--destdir",default="geo_cache")
    args=ap.parse_args()

    Path(args.outdir).mkdir(exist_ok=True); Path(args.destdir).mkdir(exist_ok=True)
    t_regex=re.compile("|".join([re.escape(t) for t in args.tissues]),re.I)

    for series in args.gse:
        gse=GEOparse.get_GEO(series,destdir=args.destdir,how="full")
        rows=[]
        for gsm in gse.gsms.values():
            r=build_row(gsm); hay=" ".join(r.values()).lower()
            if not t_regex.search(hay): continue
            rows.append(r)
        pooled=Path(args.outdir)/f"{series}_pooled.csv"; write_csv(rows,pooled)

        if args.split_by_site:
            by_site={}
            for r in rows:
                site="other"
                for t in args.tissues:
                    if t.lower() in r.get("tissue",""): site=t.lower()
                by_site.setdefault(site,[]).append(r)
            for site,rs in by_site.items():
                if site=="other": continue
                write_csv(rs,Path(args.outdir)/series/"by_site"/f"{series}_{site}.csv")

        if args.split_by_disease:
            by_dis={}
            for r in rows:
                d=r.get("diagnosis","unknown").lower()
                by_dis.setdefault(d,[]).append(r)
            for d,rs in by_dis.items():
                dtag=re.sub(r"[^a-z0-9]+","-",d).strip("-")
                write_csv(rs,Path(args.outdir)/series/"by_disease"/f"{series}_{dtag}.csv")

        if args.split_by_inflammation:
            by_inf={}
            for r in rows:
                inf=r.get("inflammation","unknown").lower()
                by_inf.setdefault(inf,[]).append(r)
            for inf,rs in by_inf.items():
                itag=re.sub(r"[^a-z0-9]+","-",inf).strip("-")
                write_csv(rs,Path(args.outdir)/series/"by_inflammation"/f"{series}_{itag}.csv")

        if args.split_cross:
            cross={}
            for r in rows:
                site="other"
                for t in args.tissues:
                    if t in r.get("tissue",""): site=t
                dis=r.get("diagnosis","unknown").lower()
                cross.setdefault((site,dis),[]).append(r)
            for (site,dis),rs in cross.items():
                if site=="other": continue
                stag=site.lower(); dtag=re.sub(r"[^a-z0-9]+","-",dis).strip("-")
                write_csv(rs,Path(args.outdir)/series/"by_cross"/f"{series}_{stag}_{dtag}.csv")

if __name__=="__main__":
    main()
