#!/usr/bin/env python3
"""Join the active-site check onto the pass-1 hit tables and rank candidates.

The structure-vs-sequence split (phmmer) is recomputed here for every
same-fold Atlas hit, with the same parameters as pass 1
(screen_v1/code/filter_and_compare.py), because pass 1 only stored the
phmmer result for the hits that became candidates. The Foldseek screen
itself is not re-run; it is read from screen_v1/results.

Outputs (in <v2_results>)
  active_site_all_hits_joined.tsv.gz  every Atlas hit with its site columns
  active_site_by_query.tsv            per query, site conservation in
                                      structure-only hits vs hits a
                                      sequence search already finds
  active_site_candidates.tsv          the 1,643 pass-1 candidates, ranked
  tier1_candidates.tsv                the shortlist

Usage: rank_candidates.py <screen_v1_dir> <v2_results>
"""
import os
import sys

import pandas as pd
import pyhmmer

TM_MIN, QCOV_MIN, ID_MAX, PHMMER_E = 0.5, 0.5, 30.0, 1e-3
TM_TIER, ID_TIER = 0.70, 20.0
RMSD_OK = 2.0      # angstrom, site CA deviation after global superposition
MAD_OK = 1.5       # angstrom, mean |difference| of site-site distances


def read_query_seqs(qdir):
    out = {}
    for fn in os.listdir(qdir):
        if not fn.endswith('.fasta'):
            continue
        name, seq = None, []
        for line in open(os.path.join(qdir, fn)):
            if line.startswith('>'):
                name = line[1:].split()[0]
            else:
                seq.append(line.strip())
        out[fn[:-6]] = (name, ''.join(seq))
    return out


def phmmer_hits(query_seq, query_name, targets):
    alpha = pyhmmer.easel.Alphabet.amino()
    q = pyhmmer.easel.TextSequence(name=query_name.encode(),
                                   sequence=query_seq).digitize(alpha)
    seqs = [pyhmmer.easel.TextSequence(name=k.encode(), sequence=v).digitize(alpha)
            for k, v in targets.items() if v]
    if not seqs:
        return {}
    db = pyhmmer.easel.DigitalSequenceBlock(alpha, seqs)
    found = {}
    for top in pyhmmer.hmmer.phmmer([q], db, E=10.0, cpus=0):
        for hit in top:
            name = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
            found[name] = hit.evalue
    return found


def classify(r):
    if not r.n_sites or r.n_sites != r.n_sites or not r.n_sites_aligned:
        return 'not_checkable'
    if (r.n_identical == r.n_sites and r.site_rmsd <= RMSD_OK
            and (r.site_dist_mad != r.site_dist_mad or r.site_dist_mad <= MAD_OK)):
        return 'site_complete'
    if r.frac_identical >= 0.7:
        return 'site_partial'
    return 'site_lost'


def catclass(r):
    """Rank on the catalytic residues alone (UniProt Active site and the
    PLP-lysine), then on the cofactor/substrate-binding residues. For the
    queries whose UniProt record annotates no catalytic residue, only the
    binding set exists and this returns 'no_catalytic_annotation'."""
    if not r.n_sites or r.n_sites != r.n_sites:
        return 'not_checkable'
    if not r.n_catalytic:
        return 'no_catalytic_annotation'
    if r.n_catalytic_identical == r.n_catalytic:
        return 'catalytic_kept'
    if r.n_catalytic_identical:
        return 'catalytic_partial'
    return 'catalytic_lost'


def classify3d(r):
    """Same question asked by 3D proximity instead of alignment register:
    for each query site CA, the physically closest candidate residue after
    superposition must be the same amino acid and within 4 A."""
    if not r.n_sites or r.n_sites != r.n_sites:
        return 'not_checkable'
    if r.frac_nearest_identical == 1.0:
        return 'site_complete_3d'
    if r.frac_nearest_identical >= 0.7:
        return 'site_partial_3d'
    return 'site_lost_3d'


def strip(s):
    return s.str.replace('.pdb.gz', '', regex=False)


def main():
    v1, v2r = sys.argv[1:3]
    v1r = os.path.join(v1, 'results')
    sites = pd.read_csv(os.path.join(v2r, 'active_site_all_hits.tsv'), sep='\t')
    sites['target'] = strip(sites.target)
    sites = sites.rename(columns={'n_aligned': 'n_sites_aligned'}).drop(columns=['acc'])

    allh = pd.read_csv(os.path.join(v1r, 'all_hits.tsv.gz'), sep='\t')
    atlas = allh[allh.db == 'mgnify_esm30'].copy()
    atlas['target'] = strip(atlas.target)
    queries = read_query_seqs(os.path.join(v1, 'queries'))

    fold_all = []
    for qname, sub in atlas.groupby('query'):
        qkey = qname.split('_w')[0]
        match = [k for k in queries if k.startswith(qkey)]
        if not match:
            continue
        qacc, qseq = queries[match[0]]
        fold = sub[(sub.tm_score_qnorm >= TM_MIN) & (sub.q_cov >= QCOV_MIN)].copy()
        ph = phmmer_hits(qseq, qacc, dict(zip(fold.target, fold.target_seq.fillna(''))))
        fold['phmmer_evalue'] = [ph.get(t, float('nan')) for t in fold.target]
        fold['seq_detectable'] = [(t in ph and ph[t] <= PHMMER_E) for t in fold.target]
        fold_all.append(fold)
    fold = pd.concat(fold_all).drop(columns=['target_seq'])

    j = fold.merge(sites, on=['query', 'target'], how='left')
    j['site_class'] = j.apply(classify, axis=1)
    j['site_class_3d'] = j.apply(classify3d, axis=1)
    j['catalytic_class'] = j.apply(catclass, axis=1)
    j.to_csv(os.path.join(v2r, 'active_site_all_hits_joined.tsv.gz'), sep='\t',
             index=False, compression='gzip')

    rows = []
    for q, g in j.groupby('query'):
        so = g[(~g.seq_detectable) & (g.seq_id_pct < ID_MAX)]
        sd = g[g.seq_detectable]
        rows.append(dict(
            query=q,
            site_residues_in_window=int(g.n_sites.max()) if g.n_sites.notna().any() else 0,
            same_fold=len(g),
            structure_only=len(so),
            so_site_complete=int((so.site_class == 'site_complete').sum()),
            so_site_partial=int((so.site_class == 'site_partial').sum()),
            so_site_lost=int((so.site_class == 'site_lost').sum()),
            so_median_frac_id=round(float(so.frac_identical.median()), 3) if len(so) else None,
            so_site_complete_3d=int((so.site_class_3d == 'site_complete_3d').sum()),
            so_catalytic_kept=int((so.catalytic_class == 'catalytic_kept').sum()),
            so_median_frac_id_3d=round(float(so.frac_nearest_identical.median()), 3) if len(so) else None,
            seq_detectable=len(sd),
            sd_site_complete=int((sd.site_class == 'site_complete').sum()),
            sd_median_frac_id=round(float(sd.frac_identical.median()), 3) if len(sd) else None,
            sd_site_complete_3d=int((sd.site_class_3d == 'site_complete_3d').sum()),
            sd_catalytic_kept=int((sd.catalytic_class == 'catalytic_kept').sum()),
            sd_median_frac_id_3d=round(float(sd.frac_nearest_identical.median()), 3) if len(sd) else None))
    by_q = pd.DataFrame(rows).sort_values('query')
    by_q.to_csv(os.path.join(v2r, 'active_site_by_query.tsv'), sep='\t', index=False)

    cand = pd.read_csv(os.path.join(v1r, 'candidates.tsv'), sep='\t')
    cand['target'] = strip(cand.target)
    c = cand.merge(sites, on=['query', 'target'], how='left')
    c['site_class'] = c.apply(classify, axis=1)
    c['site_class_3d'] = c.apply(classify3d, axis=1)
    c['catalytic_class'] = c.apply(catclass, axis=1)
    order = {'site_complete': 0, 'site_partial': 1, 'site_lost': 2, 'not_checkable': 3}
    c = c.assign(_o=c.site_class.map(order)).sort_values(
        ['_o', 'frac_identical', 'tm_score_qnorm'], ascending=[True, False, False])
    cols = ['target', 'query', 'tm_score_qnorm', 'seq_id_pct', 'q_cov', 'n_aligned',
            'phmmer_evalue', 'n_sites', 'n_sites_aligned', 'n_identical',
            'n_similar', 'frac_identical', 'site_rmsd', 'site_max_dev',
            'site_dist_mad', 'site_class', 'n_catalytic',
            'n_catalytic_identical', 'n_binding', 'n_binding_identical',
            'catalytic_class', 'n_nearest_identical',
            'frac_nearest_identical', 'nearest_max_dist', 'site_class_3d',
            'site_detail', 'nearest_detail', 't_len', 'evalue']
    c[cols].to_csv(os.path.join(v2r, 'active_site_candidates.tsv'), sep='\t', index=False)

    t1 = c[(c.catalytic_class == 'catalytic_kept') & (c.site_rmsd <= RMSD_OK)
           & (c.site_dist_mad <= MAD_OK) & (c.tm_score_qnorm >= 0.6)
           & (c.phmmer_evalue.isna())]
    t2 = c[(c.catalytic_class == 'no_catalytic_annotation')
           & (c.site_class == 'site_complete')
           & (c.site_class_3d == 'site_complete_3d')
           & (c.phmmer_evalue.isna())]
    t2.sort_values('tm_score_qnorm', ascending=False)[cols].to_csv(
        os.path.join(v2r, 'tier2_binding_only_candidates.tsv'), sep='\t', index=False)
    t1 = t1.sort_values('tm_score_qnorm', ascending=False)
    t1[cols].to_csv(os.path.join(v2r, 'tier1_candidates.tsv'), sep='\t', index=False)

    print(by_q.to_string(index=False))
    print('\ncandidate site classes by query:')
    print(c.groupby(['query', 'site_class']).size().unstack(fill_value=0).to_string())
    print('\nall candidates, by alignment-register check:')
    print(c.site_class.value_counts().to_string())
    print('\nall candidates, by 3D proximity check:')
    print(c.site_class_3d.value_counts().to_string())
    print('\ncatalytic-residue class:')
    print(c.catalytic_class.value_counts().to_string())
    print('\ncross-tab:')
    print(pd.crosstab(c.site_class, c.site_class_3d).to_string())
    print(f'\ntier 1 (catalytic residues kept, geometry ok, TM>=0.6, phmmer blind): {len(t1)}')
    print(f'tier 2 (no catalytic annotation, full binding set kept both ways): {len(t2)}')
    print(t1[['target', 'query', 'tm_score_qnorm', 'seq_id_pct', 'n_sites',
              'n_catalytic', 'n_catalytic_identical', 'n_binding',
              'n_binding_identical', 'site_rmsd', 'site_dist_mad',
              'site_detail']].head(30).to_string(index=False))


if __name__ == '__main__':
    main()
