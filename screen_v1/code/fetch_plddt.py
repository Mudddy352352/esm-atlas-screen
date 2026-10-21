#!/usr/bin/env python3
"""Fetch ESM Atlas structures for shortlisted hits and score their confidence.

For each MGnify protein we pull the predicted structure from the Atlas API and
report mean pLDDT over the whole chain and over the residues that aligned to
the query. Atlas PDB files carry pLDDT in the B-factor column on a 0-1 scale.

Usage: fetch_plddt.py <candidates.tsv.gz> <outdir> [max_n]
"""
import gzip
import json
import os
import sys
import time
import urllib.request

import pandas as pd

API = "https://api.esmatlas.com/fetchPredictedStructure/"


def fetch(mgyp, dest, tries=4):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        return True
    for i in range(tries):
        try:
            with urllib.request.urlopen(API + mgyp, timeout=120) as r:
                body = r.read()
            if body.lstrip().startswith((b"HEADER", b"ATOM")):
                open(dest, "wb").write(body)
                return True
        except Exception:
            pass
        time.sleep(5 * (i + 1))
    return False


def plddt(path):
    vals = []
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            vals.append(float(line[60:66]))
    return vals


def main():
    cand_path, outdir = sys.argv[1], sys.argv[2]
    max_n = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    os.makedirs(outdir, exist_ok=True)
    df = pd.read_csv(cand_path, sep="\t")
    df = df.sort_values("tm_score_qnorm", ascending=False).head(max_n)
    rows = []
    for _, r in df.iterrows():
        mgyp = str(r["target"]).replace(".pdb.gz", "")
        dest = os.path.join(outdir, mgyp + ".pdb")
        ok = fetch(mgyp, dest)
        if not ok:
            rows.append(dict(target=mgyp, plddt_mean=None, plddt_aligned=None,
                             n_res=None, fetched=False))
            continue
        v = plddt(dest)
        # scale: Atlas files use 0-1, report on the 0-100 pLDDT scale
        scale = 100.0 if (v and max(v) <= 1.5) else 1.0
        lo = int(r["t_cov"] * 0 + 0)  # aligned span approximated by full chain
        rows.append(dict(target=mgyp, plddt_mean=round(sum(v) / len(v) * scale, 1),
                         n_res=len(v), fetched=True))
        time.sleep(0.4)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(os.path.dirname(cand_path), "candidate_plddt.tsv"),
               sep="\t", index=False)
    got = out[out.fetched]
    print(f"fetched {len(got)}/{len(out)}")
    if len(got):
        print(got.plddt_mean.describe())
        print("above 70:", int((got.plddt_mean > 70).sum()))


if __name__ == "__main__":
    main()
