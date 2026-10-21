#!/usr/bin/env python3
"""Turn raw Foldseek server JSON into per-hit tables.

For every alignment we record the Foldseek numbers (sequence identity,
E-value, alignment length, coverage) and a TM-score estimated from the
Foldseek alignment itself: the aligned CA pairs are superimposed with Kabsch
and the superposition is refined iteratively, keeping the best TM-score.
This approximates TM-align from a fixed alignment; it is not TM-align, and
shortlisted hits are re-scored with the server's tmalign mode.

Usage: parse_results.py <raw_dir> <structures_dir> <out_dir>
"""
import glob
import json
import os
import sys

import numpy as np


def ca_from_pdb(path):
    xyz, seq = [], []
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
        'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V'}
    for line in open(path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            seq.append(three_to_one.get(line[17:20].strip(), 'X'))
    return np.array(xyz), ''.join(seq)


def kabsch(P, Q):
    """Rotation+translation taking P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def tm_score(qxyz, txyz, L_norm):
    """TM-score for pre-paired coordinates, normalised by L_norm, with
    iterative refinement of the residue subset used for superposition."""
    n = len(qxyz)
    if n < 3:
        return 0.0
    d0 = max(1.24 * (L_norm - 15) ** (1.0 / 3.0) - 1.8, 0.5) if L_norm > 15 else 0.5
    best = 0.0
    seeds = [np.arange(n)]
    if n >= 16:
        half = n // 2
        seeds += [np.arange(0, half), np.arange(half, n)]
    for idx in seeds:
        for _ in range(8):
            if len(idx) < 3:
                break
            R, t = kabsch(txyz[idx], qxyz[idx])
            moved = (R @ txyz.T).T + t
            d = np.linalg.norm(moved - qxyz, axis=1)
            score = float(np.sum(1.0 / (1.0 + (d / d0) ** 2)) / L_norm)
            if score > best:
                best = score
            cut = d0 + 1.0
            new = np.where(d < cut)[0]
            while len(new) < 3 and cut <= 20:
                cut += 1.0
                new = np.where(d < cut)[0]
            if len(new) < 3 or np.array_equal(new, idx):
                break
            idx = new
    return best


def aligned_pairs(a, qxyz, txyz):
    """Walk the Foldseek alignment strings and return paired CA coordinates."""
    qi = a['qStartPos'] - 1
    ti = a['dbStartPos'] - 1
    qp, tp = [], []
    for qc, tc in zip(a['qAln'], a['dbAln']):
        if qc != '-' and tc != '-':
            if qi < len(qxyz) and ti < len(txyz):
                qp.append(qxyz[qi])
                tp.append(txyz[ti])
        if qc != '-':
            qi += 1
        if tc != '-':
            ti += 1
    return np.array(qp), np.array(tp)


def parse_one(args):
    jf, struct_dir = args
    qname = os.path.basename(jf).replace('.pdb.json', '').replace('.json', '')
    d = json.load(open(jf))
    qxyz, _ = ca_from_pdb(os.path.join(struct_dir, qname + '.pdb'))
    rows = []
    for res in d['results']:
        al = res['alignments']
        al = al[0] if al and isinstance(al[0], list) else al
        for a in al:
            try:
                txyz = np.array([float(v) for v in a['tCa'].split(',')]).reshape(-1, 3)
            except Exception:
                continue
            qp, tp = aligned_pairs(a, qxyz, txyz)
            tm_q = tm_score(qp, tp, a['qLen']) if len(qp) >= 3 else 0.0
            tm_t = tm_score(qp, tp, a['dbLen']) if len(qp) >= 3 else 0.0
            rows.append(dict(
                query=qname, db=res['db'],
                target=str(a['target']).split()[0],
                target_desc=' '.join(str(a['target']).split()[1:])[:120],
                seq_id_pct=a['seqId'], evalue=a['eval'], bits=a['score'],
                aln_len=a['alnLength'], q_len=a['qLen'], t_len=a['dbLen'],
                q_cov=round(len(qp) / a['qLen'], 4),
                t_cov=round(len(qp) / a['dbLen'], 4),
                n_aligned=len(qp),
                tm_score_qnorm=round(tm_q, 4), tm_score_tnorm=round(tm_t, 4),
                tax_name=a.get('taxName', ''), tax_id=a.get('taxId', ''),
                target_seq=a.get('tSeq', '')))
    print(f'[parsed] {qname}: {len(rows)} alignments', flush=True)
    return rows


def main():
    import multiprocessing as mp
    raw_dir, struct_dir, out_dir = sys.argv[1:4]
    os.makedirs(out_dir, exist_ok=True)
    jobs = [(jf, struct_dir) for jf in sorted(glob.glob(os.path.join(raw_dir, '*.json')))]
    with mp.Pool(min(len(jobs), max(1, (os.cpu_count() or 2)))) as pool:
        rows = [r for chunk in pool.map(parse_one, jobs) for r in chunk]
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out_dir, 'all_hits.tsv.gz'), sep='\t', index=False,
              compression='gzip')
    df.drop(columns=['target_seq']).to_csv(
        os.path.join(out_dir, 'all_hits.tsv'), sep='\t', index=False)
    print(df.groupby(['query', 'db']).size())



if __name__ == '__main__':
    main()
