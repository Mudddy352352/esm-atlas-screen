#!/usr/bin/env python3
"""Filter Atlas hits and compare structure search with sequence search.

Steps
  1. Keep Atlas (mgnify_esm30) hits that share the query fold:
     TM-score (query-normalised) >= 0.5 and query coverage >= 0.5.
  2. Split them by sequence identity to the query (< 30% = sequence-divergent).
  3. Run phmmer (HMMER3, via pyhmmer) from the query against the set of hit
     sequences, so we can say which of the structural hits a sequence search
     would also have found. Threshold E <= 1e-3.
  4. Write per-query summary counts and a ranked candidate table.

Usage: filter_and_compare.py <results_dir>
"""
import os
import sys

import pandas as pd
import pyhmmer

TM_MIN = 0.5
QCOV_MIN = 0.5
ID_MAX = 30.0
PHMMER_E = 1e-3


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
        out[fn.replace('.fasta', '')] = (name, ''.join(seq))
    return out


def phmmer_hits(query_seq, query_name, targets):
    """targets: dict name -> sequence. Returns {name: evalue} for E <= PHMMER_E."""
    alpha = pyhmmer.easel.Alphabet.amino()
    q = pyhmmer.easel.TextSequence(name=query_name.encode(), sequence=query_seq).digitize(alpha)
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


def main():
    rdir = sys.argv[1]
    base = os.path.dirname(os.path.abspath(rdir))
    df = pd.read_csv(os.path.join(rdir, 'all_hits.tsv.gz'), sep='\t')
    queries = read_query_seqs(os.path.join(base, 'queries'))

    atlas = df[df.db == 'mgnify_esm30'].copy()
    atlas['same_fold'] = (atlas.tm_score_qnorm >= TM_MIN) & (atlas.q_cov >= QCOV_MIN)
    atlas['divergent'] = atlas.seq_id_pct < ID_MAX

    summaries, candidates = [], []
    for qname, sub in atlas.groupby('query'):
        qkey = qname.replace('.pdb', '')
        if qkey not in queries:
            # window-derived queries fall back to the parent accession
            qkey = qkey.split('_w')[0]
            matches = [k for k in queries if k.startswith(qkey)]
            qkey = matches[0] if matches else None
        if qkey is None:
            continue
        qacc, qseq = queries[qkey]
        fold = sub[sub.same_fold]
        targets = dict(zip(fold.target, fold.target_seq.fillna('')))
        ph = phmmer_hits(qseq, qacc, targets)
        fold = fold.assign(
            phmmer_evalue=[ph.get(t, float('nan')) for t in fold.target],
            seq_detectable=[(t in ph and ph[t] <= PHMMER_E) for t in fold.target])
        summaries.append(dict(
            query=qname, atlas_hits_returned=len(sub),
            same_fold=int(fold.shape[0]),
            same_fold_div=int((fold.seq_id_pct < ID_MAX).sum()),
            seq_detectable=int(fold.seq_detectable.sum()),
            structure_only=int((~fold.seq_detectable).sum()),
            structure_only_div=int(((~fold.seq_detectable) & (fold.seq_id_pct < ID_MAX)).sum()),
            min_seq_id=round(float(sub.seq_id_pct.min()), 1),
            max_tm=round(float(sub.tm_score_qnorm.max()), 3)))
        cand = fold[(~fold.seq_detectable) & (fold.seq_id_pct < ID_MAX)]
        candidates.append(cand.sort_values('tm_score_qnorm', ascending=False))

    pd.DataFrame(summaries).to_csv(os.path.join(rdir, 'summary_per_query.tsv'),
                                   sep='\t', index=False)
    cdf = pd.concat(candidates) if candidates else pd.DataFrame()
    if len(cdf):
        cdf.to_csv(os.path.join(rdir, 'candidates.tsv.gz'), sep='\t', index=False,
                   compression='gzip')
        cdf.drop(columns=['target_seq']).to_csv(
            os.path.join(rdir, 'candidates.tsv'), sep='\t', index=False)
    print(pd.DataFrame(summaries).to_string(index=False))
    print(f'\ncandidates (same fold, <{ID_MAX}% id, not found by phmmer): {len(cdf)}')


if __name__ == '__main__':
    main()
