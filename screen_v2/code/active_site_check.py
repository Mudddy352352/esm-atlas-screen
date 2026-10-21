#!/usr/bin/env python3
"""Step 3 of the study plan: active-site check on the pass-1 Foldseek hits.

For every Atlas hit of every query model we take the curated site residues of
the query (results/query_sites.tsv, built by build_sites.py), push them
through the Foldseek structural alignment, and ask three separate questions:

  1. Is the site residue aligned at all, and is it the same amino acid?
  2. After superposing the candidate on the query using the whole alignment,
     how far is each candidate site CA from the query site CA?
  3. Independently of any superposition, do the site residues sit the same
     distances apart in the candidate as in the query? (rotation-invariant,
     so it cannot be flattered by a good global fit)

Question 3 is the "conserved in three-dimensional space, not just in the
sequence alignment" part of the plan: a candidate can keep the residue types
and still place them wrongly, and it will fail here.

Target CA coordinates come from the `tCa` field of the Foldseek server JSON,
so no structure download is needed and all hits can be scored, not only the
300 whose structures were fetched in pass 1.

Usage: active_site_check.py <screen_v1_dir> <sites_tsv> <out_dir>
"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V'}

# conservative substitution groups, used only to label a mismatch as
# "similar" or "different"; nothing is counted as conserved on this basis.
GROUPS = ['DE', 'KR', 'NQ', 'ST', 'ILVM', 'FYW', 'AG', 'C', 'H', 'P']


def ca_from_pdb(path):
    xyz, seq = [], []
    for line in open(path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            seq.append(THREE_TO_ONE.get(line[17:20].strip(), 'X'))
    return np.array(xyz), ''.join(seq)


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def best_fit(qxyz, txyz, L_norm):
    """TM-score style superposition of paired CAs; returns (R, t, tm)."""
    n = len(qxyz)
    d0 = max(1.24 * (L_norm - 15) ** (1 / 3) - 1.8, 0.5) if L_norm > 15 else 0.5
    best = (None, None, -1.0)
    seeds = [np.arange(n)]
    if n >= 16:
        seeds += [np.arange(0, n // 2), np.arange(n // 2, n)]
    for idx in seeds:
        for _ in range(8):
            if len(idx) < 3:
                break
            R, t = kabsch(txyz[idx], qxyz[idx])
            d = np.linalg.norm((R @ txyz.T).T + t - qxyz, axis=1)
            tm = float(np.sum(1 / (1 + (d / d0) ** 2)) / L_norm)
            if tm > best[2]:
                best = (R, t, tm)
            cut = d0 + 1.0
            new = np.where(d < cut)[0]
            while len(new) < 3 and cut <= 20:
                cut += 1.0
                new = np.where(d < cut)[0]
            if len(new) < 3 or np.array_equal(new, idx):
                break
            idx = new
    return best


def group_of(a):
    for g in GROUPS:
        if a in g:
            return g
    return a


def walk(a):
    """query-model position -> (target position, target residue)."""
    qi, ti = a['qStartPos'] - 1, a['dbStartPos'] - 1
    pairs = {}
    for qc, tc in zip(a['qAln'], a['dbAln']):
        if qc != '-' and tc != '-':
            pairs[qi + 1] = (ti + 1, tc)
        if qc != '-':
            qi += 1
        if tc != '-':
            ti += 1
    return pairs


def main():
    v1, sites_tsv, out_dir = sys.argv[1:4]
    os.makedirs(out_dir, exist_ok=True)
    sites = pd.read_csv(sites_tsv, sep='\t')
    rows = []
    skipped = []

    for jf in sorted(glob.glob(os.path.join(v1, 'foldseek_raw', '*.json'))):
        qname = os.path.basename(jf).replace('.pdb.json', '')
        acc = qname.split('_')[0]
        m = re.search(r'_w(\d+)-(\d+)$', qname)
        offset = int(m.group(1)) if m else 1
        qxyz, qseq = ca_from_pdb(os.path.join(v1, 'structures', qname + '.pdb'))

        s = sites[sites.acc == acc].copy()
        s['wpos'] = s.pos - offset + 1
        s = s[(s.wpos >= 1) & (s.wpos <= len(qseq))]
        # sanity: the model really carries the annotated residue there
        bad = [(int(r.pos), r.res, qseq[int(r.wpos) - 1]) for r in s.itertuples()
               if qseq[int(r.wpos) - 1] != r.res]
        if bad:
            skipped.append((qname, 'residue mismatch in model', str(bad)))
            continue
        if s.empty:
            skipped.append((qname, 'no curated site residues inside this window', ''))
            continue
        s = s.sort_values('wpos')
        wpos = s.wpos.to_numpy()
        q_site_xyz = qxyz[wpos - 1]
        qd = np.linalg.norm(q_site_xyz[:, None] - q_site_xyz[None, :], axis=-1)
        iu = np.triu_indices(len(wpos), 1)

        d = json.load(open(jf))
        for res in d['results']:
            if res['db'] != 'mgnify_esm30':
                continue
            al = res['alignments']
            al = al[0] if al and isinstance(al[0], list) else al
            for a in al:
                try:
                    txyz = np.array([float(v) for v in a['tCa'].split(',')]).reshape(-1, 3)
                except Exception:
                    continue
                pairs = walk(a)
                qp = np.array([qxyz[i - 1] for i in pairs if i <= len(qxyz)])
                tp = np.array([txyz[pairs[i][0] - 1] for i in pairs
                               if i <= len(qxyz) and pairs[i][0] <= len(txyz)])
                if len(qp) < 3 or len(qp) != len(tp):
                    continue
                R, t, _ = best_fit(qp, tp, a['qLen'])
                if R is None:
                    continue

                tseq = a.get('tSeq', '') or ''
                moved_all = (R @ txyz.T).T + t

                n_cat = n_cat_al = n_cat_id = 0
                n_bind = n_bind_id = 0
                n_al = n_id = n_sim = 0
                devs, detail = [], []
                t_site_xyz, keep = [], []
                # 3D proximity check, independent of alignment register:
                # after superposition, which candidate residue is physically
                # closest to each query site CA, and is it the same residue?
                n_near_id = 0
                near_d, near_detail = [], []
                for r in s.itertuples():
                    qv = qxyz[int(r.wpos) - 1]
                    dd = np.linalg.norm(moved_all - qv, axis=1)
                    k = int(np.argmin(dd))
                    nres = tseq[k] if k < len(tseq) else 'X'
                    near_d.append(float(dd[k]))
                    ok = (nres == r.res) and dd[k] <= 4.0
                    n_near_id += ok
                    near_detail.append(f'{r.res}{int(r.pos)}:{nres}@{dd[k]:.1f}')

                for k, r in enumerate(s.itertuples()):
                    cat = r.kind in ('Active site', 'Modified residue')
                    n_cat += cat
                    n_bind += (not cat)
                    hit = pairs.get(int(r.wpos))
                    tag = f'{r.res}{int(r.pos)}'
                    if hit is None or hit[0] > len(txyz):
                        detail.append(tag + ':-')
                        continue
                    tpos, tres = hit
                    n_al += 1
                    same = (tres == r.res)
                    n_id += same
                    if cat:
                        n_cat_al += 1
                        n_cat_id += same
                    else:
                        n_bind_id += same
                    sim = (not same) and group_of(tres) == group_of(r.res)
                    n_sim += sim
                    moved = R @ txyz[tpos - 1] + t
                    dev = float(np.linalg.norm(moved - qxyz[int(r.wpos) - 1]))
                    devs.append(dev)
                    t_site_xyz.append(txyz[tpos - 1])
                    keep.append(k)
                    detail.append(f'{tag}:{tres}{"=" if same else ("~" if sim else "x")}{dev:.1f}')

                # rotation-invariant geometry of the site constellation
                geom = np.nan
                if len(keep) >= 3:
                    ts = np.array(t_site_xyz)
                    td = np.linalg.norm(ts[:, None] - ts[None, :], axis=-1)
                    sub = np.ix_(keep, keep)
                    a_ = qd[sub][np.triu_indices(len(keep), 1)]
                    b_ = td[np.triu_indices(len(keep), 1)]
                    geom = float(np.mean(np.abs(a_ - b_)))

                rows.append(dict(
                    query=qname, acc=acc, target=str(a['target']).split()[0],
                    n_sites=len(s), n_aligned=n_al, n_identical=n_id,
                    n_catalytic=n_cat, n_catalytic_aligned=n_cat_al,
                    n_catalytic_identical=n_cat_id,
                    n_binding=n_bind, n_binding_identical=n_bind_id,
                    n_similar=n_sim,
                    frac_identical=round(n_id / len(s), 3),
                    site_rmsd=round(float(np.sqrt(np.mean(np.square(devs)))), 2) if devs else np.nan,
                    site_max_dev=round(float(np.max(devs)), 2) if devs else np.nan,
                    site_dist_mad=round(geom, 2) if geom == geom else np.nan,
                    n_nearest_identical=n_near_id,
                    frac_nearest_identical=round(n_near_id / len(s), 3),
                    nearest_max_dist=round(float(np.max(near_d)), 2) if near_d else np.nan,
                    site_detail=';'.join(detail),
                    nearest_detail=';'.join(near_detail)))
        print(f'[sites] {qname}: {len(s)} site residues in window, '
              f'{sum(1 for r in rows if r["query"] == qname)} hits scored', flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out_dir, 'active_site_all_hits.tsv'), sep='\t', index=False)
    with open(os.path.join(out_dir, 'active_site_skipped.tsv'), 'w') as fh:
        fh.write('query\treason\tdetail\n')
        for r in skipped:
            fh.write('\t'.join(r) + '\n')
    print(df.shape)
    print('\n'.join(f'{a}: {b} | {c}' for a, b, c in skipped))


if __name__ == '__main__':
    main()
