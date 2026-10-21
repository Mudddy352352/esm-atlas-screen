#!/usr/bin/env python3
"""Locate each tier-1 candidate's active site in its own structure.

The pass-2 active-site check knows which query site residues survive, but it
reports them in *query* numbering. To place a docking box we need the
coordinates in the candidate's own structure, so this script repeats the
alignment walk from screen_v2/code/active_site_check.py and writes out, for
each candidate, the residue indices and the centroid of its site CAs.

Usage: map_sites.py <screen_v1_dir> <query_sites.tsv> <tier1.tsv> <struct_dirs...> <out.tsv>
"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V'}


def ca_from_pdb(path):
    """CA coordinates, one-letter sequence, and author residue numbers."""
    xyz, seq, num = [], [], []
    for line in open(path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            seq.append(THREE_TO_ONE.get(line[17:20].strip(), 'X'))
            num.append(int(line[22:26]))
    return np.array(xyz), ''.join(seq), num


def walk(a):
    """query-model position -> (target position, target residue)"""
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
    v1, sites_tsv, tier1_tsv = sys.argv[1:4]
    out_tsv = sys.argv[-1]
    struct_dirs = sys.argv[4:-1]

    sites = pd.read_csv(sites_tsv, sep='\t')
    tier1 = pd.read_csv(tier1_tsv, sep='\t')
    want = {(r.target, r.query) for r in tier1.itertuples()}
    rows = []

    for qname in sorted({q for _, q in want}):
        jf = os.path.join(v1, 'foldseek_raw', qname + '.pdb.json')
        if not os.path.exists(jf):
            sys.exit('missing foldseek json: ' + jf)
        acc = qname.split('_')[0]
        m = re.search(r'_w(\d+)-(\d+)$', qname)
        offset = int(m.group(1)) if m else 1
        qxyz, qseq, _ = ca_from_pdb(os.path.join(v1, 'structures', qname + '.pdb'))

        s = sites[sites.acc == acc].copy()
        s['wpos'] = s.pos - offset + 1
        s = s[(s.wpos >= 1) & (s.wpos <= len(qseq))].sort_values('wpos')
        targets = {t for t, q in want if q == qname}

        d = json.load(open(jf))
        for res in d['results']:
            if res['db'] != 'mgnify_esm30':
                continue
            al = res['alignments']
            al = al[0] if al and isinstance(al[0], list) else al
            for a in al:
                tgt = re.sub(r'\.pdb(\.gz)?$', '', a['target'])
                if tgt not in targets:
                    continue
                path = None
                for sd in struct_dirs:
                    p = os.path.join(sd, tgt + '.pdb')
                    if os.path.exists(p):
                        path = p
                        break
                if path is None:
                    sys.exit('no structure for ' + tgt)
                txyz, tseq_pdb, tnum = ca_from_pdb(path)
                tseq_aln = a.get('tSeq', '') or ''
                # the docking box must be placed on the *downloaded* structure,
                # so refuse to proceed unless its sequence is the one Foldseek
                # aligned against
                agree = (tseq_aln[:len(tseq_pdb)] == tseq_pdb) if tseq_aln else None

                pairs = walk(a)
                idx, detail = [], []
                for r in s.itertuples():
                    hit = pairs.get(int(r.wpos))
                    tag = f'{r.res}{int(r.pos)}'
                    if hit is None or hit[0] > len(txyz):
                        detail.append(tag + ':-')
                        continue
                    tpos, tres = hit
                    detail.append(f'{tag}:{tres}{tnum[tpos - 1]}')
                    idx.append((tpos, tres, r.kind))

                if not idx:
                    continue
                cat = [i for i in idx if i[2] in ('Active site', 'Modified residue')]
                # centre the box on the catalytic residues when they are there,
                # otherwise on all mapped site residues
                use = cat if cat else idx
                cen_cat = txyz[[i[0] - 1 for i in use]].mean(0)
                cen_all = txyz[[i[0] - 1 for i in idx]].mean(0)
                spread = float(np.max(np.linalg.norm(
                    txyz[[i[0] - 1 for i in idx]] - cen_all, axis=1)))
                rows.append(dict(
                    target=tgt, query=qname, pdb=path, t_len=len(txyz),
                    seq_matches_aln=agree,
                    n_site_mapped=len(idx), n_cat_mapped=len(cat),
                    cat_resnums=';'.join(str(tnum[i[0] - 1]) for i in cat),
                    site_resnums=';'.join(str(tnum[i[0] - 1]) for i in idx),
                    cx=round(float(cen_cat[0]), 3), cy=round(float(cen_cat[1]), 3),
                    cz=round(float(cen_cat[2]), 3),
                    site_spread_ang=round(spread, 2),
                    mapping=';'.join(detail)))

    out = pd.DataFrame(rows).sort_values(['query', 'target'])
    out.to_csv(out_tsv, sep='\t', index=False)
    print(out[['target', 'query', 't_len', 'seq_matches_aln', 'n_site_mapped',
               'n_cat_mapped', 'cat_resnums', 'site_spread_ang']].to_string(index=False))
    missing = want - {(r.target, r.query) for r in out.itertuples()}
    if missing:
        print('NOT MAPPED:', missing)


if __name__ == '__main__':
    main()
