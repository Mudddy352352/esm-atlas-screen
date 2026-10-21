#!/usr/bin/env python3
"""Build the curated active-site table for the query enzymes.

Two sources, both recorded per residue so the provenance stays visible:

  direct      the residue is annotated on the query accession itself in
              UniProtKB (Active site / Binding site / Modified residue
              feature).
  transferred the query accession has no site features, so residues are
              taken from a reviewed, experimentally annotated member of the
              same family and carried across a pairwise alignment. Only
              positions where the aligned residue is identical are kept.

Nothing here is typed in from memory or from the literature; every residue
comes out of a UniProtKB JSON record fetched in this run.

Usage: build_sites.py <uniprot_json_dir> <out_dir>
"""
import json
import os
import sys

from Bio import Align

# query accession -> annotated donor accession, for the two queries that
# carry no site features of their own.
TRANSFER = {
    'P07650': 'P77836',   # E. coli deoA  <- Geobacillus pyrimidine-nucleoside phosphorylase
    'Q7A782': 'O07529',   # S. aureus azo1 <- B. subtilis FMN-dependent NADPH-azoreductase
}

WANTED = ('Active site', 'Binding site', 'Modified residue')


def load(acc, d):
    rec = json.load(open(os.path.join(d, acc + '.json')))
    name = rec['proteinDescription']['recommendedName']['fullName']['value']
    return rec['sequence']['value'], name, rec.get('features', [])


def sites_of(seq, feats):
    """Expand site features to one row per residue."""
    out = []
    for ft in feats:
        if ft['type'] not in WANTED:
            continue
        if ft['type'] == 'Modified residue' and 'pyridoxal' not in ft.get('description', ''):
            continue
        s = ft['location']['start']['value']
        e = ft['location']['end']['value']
        if not s or not e:
            continue
        lig = ft.get('ligand', {}).get('name', '') or ft.get('description', '')
        ev = ';'.join(sorted({x['evidenceCode'] for x in ft.get('evidences', [])}))
        for p in range(s, e + 1):
            out.append(dict(pos=p, res=seq[p - 1], kind=ft['type'], ligand=lig,
                            evidence=ev))
    # one row per position; keep the most specific kind
    best = {}
    rank = {'Active site': 3, 'Modified residue': 2, 'Binding site': 1}
    for r in out:
        cur = best.get(r['pos'])
        if cur is None or rank[r['kind']] > rank[cur['kind']]:
            best[r['pos']] = r
        elif rank[r['kind']] == rank[cur['kind']] and r['ligand'] not in cur['ligand']:
            cur['ligand'] = (cur['ligand'] + '+' + r['ligand']).strip('+')
    return [best[p] for p in sorted(best)]


def align_pair(a, b):
    al = Align.PairwiseAligner(open_gap_score=-11, extend_gap_score=-1,
                               substitution_matrix=Align.substitution_matrices.load('BLOSUM62'))
    al.mode = 'global'
    aln = al.align(a, b)[0]
    # map donor position (1-based) -> query position (1-based)
    qa, qb = aln[0], aln[1]
    ia = ib = 0
    m = {}
    ident = n = 0
    for ca, cb in zip(qa, qb):
        if ca != '-' and cb != '-':
            m[ib + 1] = (ia + 1, ca, cb)
            n += 1
            ident += (ca == cb)
        if ca != '-':
            ia += 1
        if cb != '-':
            ib += 1
    return m, (100.0 * ident / n if n else 0.0)


def main():
    ud, out_dir = sys.argv[1:3]
    os.makedirs(out_dir, exist_ok=True)
    queries = ['P05804', 'P07650', 'P0ABF6', 'P0DTQ4', 'P15888', 'P17117',
               'P19079', 'Q7A782', 'Q8KU07']
    rows = []
    notes = []
    for acc in queries:
        seq, name, feats = load(acc, ud)
        direct = sites_of(seq, feats)
        if direct:
            for r in direct:
                rows.append(dict(acc=acc, protein=name, source='direct',
                                 donor='', **r))
            notes.append(f'{acc}\t{name}\tdirect\t{len(direct)} residues')
            continue
        donor = TRANSFER.get(acc)
        if not donor:
            notes.append(f'{acc}\t{name}\tnone\t0 residues (no site features, no donor)')
            continue
        dseq, dname, dfeats = load(donor, ud)
        dsites = sites_of(dseq, dfeats)
        m, pid = align_pair(seq, dseq)
        kept = dropped = 0
        for r in dsites:
            hit = m.get(r['pos'])
            if hit is None:
                dropped += 1
                continue
            qpos, qres, dres = hit
            if qres != dres:          # only transfer identical residues
                dropped += 1
                continue
            rows.append(dict(acc=acc, protein=name, source='transferred',
                             donor=donor, pos=qpos, res=qres, kind=r['kind'],
                             ligand=r['ligand'], evidence=r['evidence']))
            kept += 1
        notes.append(f'{acc}\t{name}\ttransferred from {donor} ({dname}, '
                     f'{pid:.1f}% id)\t{kept} residues kept, {dropped} dropped')

    import pandas as pd
    df = pd.DataFrame(rows).sort_values(['acc', 'pos'])
    df.to_csv(os.path.join(out_dir, 'query_sites.tsv'), sep='\t', index=False)
    with open(os.path.join(out_dir, 'query_sites_provenance.tsv'), 'w') as fh:
        fh.write('accession\tprotein\tsource\tdetail\n')
        fh.write('\n'.join(notes) + '\n')
    print(df.groupby(['acc', 'source', 'kind']).size())
    print('\n'.join(notes))


if __name__ == '__main__':
    main()
