#!/usr/bin/env python3
"""Score the LABELLED pass-1 hits (AlphaFold DB and PDB) with the same
active-site check the screen applied to the Atlas hits.

The screen's own check (screen_v2/code/active_site_check.py) skips every hit
that is not from mgnify_esm30. The Atlas carries no functional annotation, so
nothing in the paper measures whether the filter separates enzymes that do the
chemistry from same-fold enzymes that do not. The AlphaFold DB and PDB hits of
the same searches do carry annotation, so they give ground truth for free.

Nothing here reimplements the check: the geometry, the alignment walk and the
superposition are imported from the screen's own module, and only the database
filter changes.

Usage: score_labelled_hits.py <screen_v1_dir> <sites_tsv> <out_tsv> [dbs]
"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

V2_CODE = _ROOT + "/screen_v2/code"
sys.path.insert(0, V2_CODE)
from active_site_check import ca_from_pdb, best_fit, walk, group_of  # noqa: E402

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))


def score(v1, sites_tsv, out_tsv, dbs):
    sites = pd.read_csv(sites_tsv, sep="\t")
    rows, skipped = [], []

    for jf in sorted(glob.glob(os.path.join(v1, "foldseek_raw", "*.json"))):
        qname = os.path.basename(jf).replace(".pdb.json", "")
        acc = qname.split("_")[0]
        m = re.search(r"_w(\d+)-(\d+)$", qname)
        offset = int(m.group(1)) if m else 1
        qxyz, qseq = ca_from_pdb(os.path.join(v1, "structures", qname + ".pdb"))

        s = sites[sites.acc == acc].copy()
        s["wpos"] = s.pos - offset + 1
        s = s[(s.wpos >= 1) & (s.wpos <= len(qseq))]
        bad = [(int(r.pos), r.res, qseq[int(r.wpos) - 1]) for r in s.itertuples()
               if qseq[int(r.wpos) - 1] != r.res]
        if bad:
            skipped.append((qname, "residue mismatch in model"))
            continue
        if s.empty:
            skipped.append((qname, "no curated site residues in this window"))
            continue
        s = s.sort_values("wpos")
        wpos = s.wpos.to_numpy()
        q_site_xyz = qxyz[wpos - 1]
        qd = np.linalg.norm(q_site_xyz[:, None] - q_site_xyz[None, :], axis=-1)

        d = json.load(open(jf))
        n_before = len(rows)
        for res in d["results"]:
            if res["db"] not in dbs:
                continue
            al = res["alignments"]
            al = al[0] if al and isinstance(al[0], list) else al
            for a in al:
                try:
                    txyz = np.array([float(v) for v in a["tCa"].split(",")]).reshape(-1, 3)
                except Exception:
                    continue
                pairs = walk(a)
                qp = np.array([qxyz[i - 1] for i in pairs if i <= len(qxyz)])
                tp = np.array([txyz[pairs[i][0] - 1] for i in pairs
                               if i <= len(qxyz) and pairs[i][0] <= len(txyz)])
                if len(qp) < 3 or len(qp) != len(tp):
                    continue
                R, t, tm = best_fit(qp, tp, a["qLen"])
                if R is None:
                    continue

                tseq = a.get("tSeq", "") or ""
                n_cat = n_cat_al = n_cat_id = 0
                n_al = n_id = 0
                devs, detail, t_site_xyz, keep = [], [], [], []
                for k, r in enumerate(s.itertuples()):
                    cat = r.kind in ("Active site", "Modified residue")
                    n_cat += cat
                    hit = pairs.get(int(r.wpos))
                    tag = f"{r.res}{int(r.pos)}"
                    if hit is None or hit[0] > len(txyz):
                        detail.append(tag + ":-")
                        continue
                    tpos, tres = hit
                    n_al += 1
                    same = (tres == r.res)
                    n_id += same
                    if cat:
                        n_cat_al += 1
                        n_cat_id += same
                    moved = R @ txyz[tpos - 1] + t
                    devs.append(float(np.linalg.norm(moved - qxyz[int(r.wpos) - 1])))
                    t_site_xyz.append(txyz[tpos - 1])
                    keep.append(k)
                    sim = (not same) and group_of(tres) == group_of(r.res)
                    detail.append(f'{tag}:{tres}{"=" if same else ("~" if sim else "x")}{devs[-1]:.1f}')

                geom = np.nan
                if len(keep) >= 3:
                    ts = np.array(t_site_xyz)
                    td = np.linalg.norm(ts[:, None] - ts[None, :], axis=-1)
                    sub = np.ix_(keep, keep)
                    a_ = qd[sub][np.triu_indices(len(keep), 1)]
                    b_ = td[np.triu_indices(len(keep), 1)]
                    geom = float(np.mean(np.abs(a_ - b_)))

                target_full = str(a["target"])
                rows.append(dict(
                    query=qname, acc=acc, db=res["db"],
                    target=target_full.split()[0],
                    target_desc=" ".join(target_full.split()[1:]),
                    tax_name=a.get("taxName", ""),
                    seq_id_pct=round(float(a.get("seqId", 0)), 1),
                    evalue=a.get("eval"), tm_score_qnorm=round(float(tm), 4),
                    n_sites=len(s), n_aligned=n_al, n_identical=n_id,
                    n_catalytic=n_cat, n_catalytic_aligned=n_cat_al,
                    n_catalytic_identical=n_cat_id,
                    frac_identical=round(n_id / len(s), 3),
                    site_rmsd=round(float(np.sqrt(np.mean(np.square(devs)))), 2) if devs else np.nan,
                    site_dist_mad=round(geom, 2) if geom == geom else np.nan,
                    site_detail=";".join(detail)))
        print(f"[bench] {qname}: {len(s)} site residues, "
              f"{len(rows) - n_before} labelled hits scored", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(out_tsv, sep="\t", index=False)
    print("wrote", out_tsv, df.shape, flush=True)
    for q, r in skipped:
        print("skipped", q, r, flush=True)


if __name__ == "__main__":
    v1, sites_tsv, out_tsv = sys.argv[1:4]
    dbs = set((sys.argv[4] if len(sys.argv) > 4 else "afdb50,pdb100").split(","))
    score(v1, sites_tsv, out_tsv, dbs)
