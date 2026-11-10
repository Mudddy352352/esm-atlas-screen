#!/usr/bin/env python3
"""How often does the active-site filter pass a candidate by chance?

The screen reports that 5 of 1,643 structure-only relatives keep a complete
annotated site and 31 keep part of it, but nothing in the paper says how often
those thresholds would be met by residues that mean nothing. This builds the
null: the same check, on the same hits, against decoy site definitions drawn
from the query model itself.

Three decoys, increasingly hard:
  random     - positions anywhere in the query, same number of residues
  composed   - same number AND the same residue types as the real site,
               which controls for how common glutamate and aspartate are
  clustered  - composition-matched and spatially clustered like a real site,
               which also controls for the geometry terms

Everything else, the alignment walk, the superposition and the thresholds, is
the screen's own code, imported rather than reimplemented.

Usage: decoy_site_null.py <screen_v1> <sites_tsv> <candidates_tsv> <out_dir> [n_rep]
"""
import json
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, _ROOT + "/screen_v2/code")
from active_site_check import ca_from_pdb, best_fit, walk  # noqa: E402

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

RMSD_OK, MAD_OK = 2.0, 1.5          # screen_v2/code/rank_candidates.py
RNG = np.random.default_rng(20261107)


def decoy_positions(kind, qseq, qxyz, real_wpos, rng):
    n = len(real_wpos)
    L = len(qseq)
    if kind == "random":
        return rng.choice(np.arange(1, L + 1), size=n, replace=False)
    types = [qseq[p - 1] for p in real_wpos]
    pools = {a: np.array([i + 1 for i, c in enumerate(qseq) if c == a]) for a in set(types)}
    if any(len(pools[a]) < types.count(a) for a in set(types)):
        return None
    if kind == "composed":
        out, used = [], set()
        for a in types:
            choices = [p for p in pools[a] if p not in used]
            if not choices:
                return None
            p = int(rng.choice(choices)); used.add(p); out.append(p)
        return np.array(out)
    if kind == "clustered":
        # pick a seed, then take the nearest unused residue of each required
        # type, so the decoy site is as compact as a real one
        for _ in range(60):
            seed = int(rng.integers(1, L + 1))
            centre = qxyz[seed - 1]
            out, used = [], set()
            ok = True
            for a in types:
                cand = [p for p in pools[a] if p not in used]
                if not cand:
                    ok = False
                    break
                d = np.linalg.norm(qxyz[np.array(cand) - 1] - centre, axis=1)
                p = int(cand[int(np.argmin(d))]); used.add(p); out.append(p)
            if ok:
                return np.array(out)
        return None
    raise ValueError(kind)


def site_geometry(qxyz, wpos):
    xyz = qxyz[np.asarray(wpos) - 1]
    d = np.linalg.norm(xyz[:, None] - xyz[None, :], axis=-1)
    return xyz, d


def score_one(a, qxyz, qseq, wpos, residues, qd):
    try:
        txyz = np.array([float(v) for v in a["tCa"].split(",")]).reshape(-1, 3)
    except Exception:
        return None
    pairs = walk(a)
    qp = np.array([qxyz[i - 1] for i in pairs if i <= len(qxyz)])
    tp = np.array([txyz[pairs[i][0] - 1] for i in pairs
                   if i <= len(qxyz) and pairs[i][0] <= len(txyz)])
    if len(qp) < 3 or len(qp) != len(tp):
        return None
    R, t, _ = best_fit(qp, tp, a["qLen"])
    if R is None:
        return None
    n_id = n_al = 0
    devs, t_site, keep = [], [], []
    for k, (p, res) in enumerate(zip(wpos, residues)):
        hit = pairs.get(int(p))
        if hit is None or hit[0] > len(txyz):
            continue
        tpos, tres = hit
        n_al += 1
        n_id += (tres == res)
        moved = R @ txyz[tpos - 1] + t
        devs.append(float(np.linalg.norm(moved - qxyz[int(p) - 1])))
        t_site.append(txyz[tpos - 1])
        keep.append(k)
    if not n_al:
        return None
    rmsd = float(np.sqrt(np.mean(np.square(devs))))
    mad = np.nan
    if len(keep) >= 3:
        ts = np.array(t_site)
        td = np.linalg.norm(ts[:, None] - ts[None, :], axis=-1)
        sub = np.ix_(keep, keep)
        mad = float(np.mean(np.abs(qd[sub][np.triu_indices(len(keep), 1)]
                                   - td[np.triu_indices(len(keep), 1)])))
    complete = (n_id == len(wpos) and rmsd <= RMSD_OK
                and (mad != mad or mad <= MAD_OK))
    partial = (not complete) and (n_id / len(wpos) >= 0.7)
    return dict(target=str(a["target"]).split()[0].replace(".pdb.gz", ""), n_identical=n_id,
                site_rmsd=rmsd, site_dist_mad=mad,
                complete=bool(complete), partial=bool(partial))


def main(v1, sites_tsv, cand_tsv, out_dir, n_rep=10, mode="all", tm_min=0.0):
    n_rep = int(n_rep)
    tm_min = float(tm_min)
    os.makedirs(out_dir, exist_ok=True)
    sites = pd.read_csv(sites_tsv, sep="\t")
    cdf = pd.read_csv(cand_tsv, sep="\t")
    if tm_min:
        cdf = cdf[cdf.tm_score_qnorm >= tm_min]
    cands = set(cdf.target)
    if mode == "catalytic":
        sites = sites[sites.kind.isin(["Active site", "Modified residue"])]
    rows = []

    for jf in sorted(glob.glob(os.path.join(v1, "foldseek_raw", "*.json"))):
        qname = os.path.basename(jf).replace(".pdb.json", "")
        acc = qname.split("_")[0]
        m = re.search(r"_w(\d+)-(\d+)$", qname)
        offset = int(m.group(1)) if m else 1
        qxyz, qseq = ca_from_pdb(os.path.join(v1, "structures", qname + ".pdb"))
        s = sites[sites.acc == acc].copy()
        s["wpos"] = s.pos - offset + 1
        s = s[(s.wpos >= 1) & (s.wpos <= len(qseq))].sort_values("wpos")
        if s.empty or any(qseq[int(r.wpos) - 1] != r.res for r in s.itertuples()):
            continue
        real_wpos = s.wpos.to_numpy()
        real_res = [qseq[p - 1] for p in real_wpos]

        d = json.load(open(jf))
        al = []
        for res in d["results"]:
            if res["db"] != "mgnify_esm30":
                continue
            block = res["alignments"]
            al = block[0] if block and isinstance(block[0], list) else block
        al = [a for a in al
          if str(a["target"]).split()[0].replace(".pdb.gz", "") in cands]
        if not al:
            continue

        jobs = [("real", 0, real_wpos, real_res)]
        for kind in ("random", "composed", "clustered"):
            for rep in range(n_rep):
                pos = decoy_positions(kind, qseq, qxyz, real_wpos, RNG)
                if pos is None:
                    continue
                jobs.append((kind, rep, pos, [qseq[p - 1] for p in pos]))

        for kind, rep, wpos, residues in jobs:
            _, qd = site_geometry(qxyz, wpos)
            n = n_complete = n_partial = 0
            for a in al:
                r = score_one(a, qxyz, qseq, np.asarray(wpos), residues, qd)
                if r is None:
                    continue
                n += 1
                n_complete += r["complete"]
                n_partial += r["partial"]
            rows.append(dict(query=qname, acc=acc, site_kind=kind, replicate=rep,
                             n_site_residues=len(wpos), n_hits=n,
                             n_complete=n_complete, n_partial=n_partial))
        print("[null] %s: %d candidate hits, %d site sets" % (qname, len(al), len(jobs)),
              flush=True)

    df = pd.DataFrame(rows)
    tag = "" if mode == "all" else "_" + mode
    df.to_csv(os.path.join(out_dir, "decoy_site_null%s.tsv" % tag), sep="\t", index=False)
    agg = (df.groupby("site_kind")[["n_hits", "n_complete", "n_partial"]].sum()
             .assign(pct_complete=lambda x: 100 * x.n_complete / x.n_hits,
                     pct_partial=lambda x: 100 * x.n_partial / x.n_hits))
    agg.to_csv(os.path.join(out_dir, "decoy_site_null%s_summary.tsv" % tag), sep="\t")
    print(agg.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:])
