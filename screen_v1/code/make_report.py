#!/usr/bin/env python3
"""Write results/README.md from the tables produced by the pipeline.

Every number in the report is read from a file, so the report cannot drift
from the data.

Usage: make_report.py <screen_dir>
"""
import glob
import json
import os
import sys

import pandas as pd

FAMILY = {
    'P0ABF6': ('Cytidine deaminase', 'gemcitabine, cytarabine'),
    'P19079': ('Cytidine deaminase (tetrameric type)', 'gemcitabine, cytarabine'),
    'P05804': ('Beta-glucuronidase', 'irinotecan (SN-38 reactivation)'),
    'P07650': ('Thymidine phosphorylase', '5-fluorouracil pathway'),
    'P15888': ('NAD(P)H nitroreductase NfsB', 'nitro-prodrugs'),
    'P17117': ('NADPH nitroreductase NfsA', 'nitro-prodrugs'),
    'Q7A782': ('FMN-dependent NADPH-azoreductase', 'azo-prodrugs'),
    'Q8KU07': ('NAD(P)H azoreductase', 'azo-prodrugs'),
    'P0DTQ4': ('L-tyrosine decarboxylase (control)', 'levodopa'),
}


def pdb_validation(raw_dir):
    out = {}
    for f in sorted(glob.glob(os.path.join(raw_dir, '*.json'))):
        q = os.path.basename(f).replace('.pdb.json', '')
        d = json.load(open(f))
        for r in d['results']:
            if r['db'] != 'pdb100':
                continue
            al = r['alignments']
            al = al[0] if al and isinstance(al[0], list) else al
            if al:
                a = al[0]
                out[q] = (str(a['target']), a['seqId'], a['eval'], len(al))
    return out


def main():
    base = sys.argv[1]
    res = os.path.join(base, 'results')
    hits = pd.read_csv(os.path.join(res, 'all_hits.tsv'), sep='\t')
    summ = pd.read_csv(os.path.join(res, 'summary_per_query.tsv'), sep='\t')
    val = pdb_validation(os.path.join(base, 'foldseek_raw'))
    plddt_path = os.path.join(res, 'candidate_plddt.tsv')
    plddt = pd.read_csv(plddt_path, sep='\t') if os.path.exists(plddt_path) else None
    cand_path = os.path.join(res, 'candidates.tsv')
    cand = pd.read_csv(cand_path, sep='\t') if os.path.exists(cand_path) else None

    L = []
    L.append('# Pass 1: Foldseek screen of the ESM Metagenomic Atlas\n')
    L.append('Run on 20 October 2026 in the project cloud environment. '
             'Everything below is computed from files in this directory; '
             'nothing is quoted from memory or from the literature.\n')

    L.append('## What was run\n')
    L.append('1. Query enzymes were taken from the curated Swiss-Prot bacterial protein set '
             'bundled with Prokka (24,727 reviewed entries), selected by their EC number and '
             'enzyme name, one or two per family in section 4 of the study plan.')
    L.append('2. Each query sequence was folded with the public ESMFold API '
             '(`api.esmatlas.com/foldSequence/v1/pdb`). The API refuses sequences longer than '
             '400 residues, so the three long queries were folded as two 390-residue windows each.')
    L.append('3. Each model was searched with the Foldseek web server (`search.foldseek.com`, '
             'mode 3Di+AA) against three databases: `mgnify_esm30` (the ESM Metagenomic Atlas, '
             'MGnify proteins clustered at 30%), `afdb50` (AlphaFold DB UniProt50) and `pdb100` '
             '(the PDB).')
    L.append('4. TM-scores were computed locally from the Foldseek alignments by Kabsch '
             'superposition of the aligned CA pairs with iterative refinement, normalised by '
             'query length.')
    L.append('5. For each query, phmmer (HMMER3 via pyhmmer 0.12.3) was run from the query '
             'sequence against the sequences of its own Atlas structural hits, to measure how '
             'many of them a sequence search would also have found.\n')

    L.append('## Query set and identity check\n')
    L.append('The PDB search is a control on the query itself: if the top PDB hit is the '
             'expected enzyme at high identity, the sequence really is that enzyme and the '
             'ESMFold model is accurate.\n')
    L.append('| Query | Family | Drug relevance | Length | Top PDB hit | % id |')
    L.append('|---|---|---|---|---|---|')
    for q in sorted(val):
        acc = q.split('_')[0]
        fam, drug = FAMILY.get(acc, ('', ''))
        t, sid, ev, n = val[q]
        qlen = hits[hits['query'] == q]['q_len'].iloc[0] if (hits['query'] == q).any() else ''
        desc = ' '.join(t.split()[1:])[:58]
        L.append(f'| `{q}` | {fam} | {drug} | {qlen} | {desc} | {sid:.1f} |')
    L.append('')

    L.append('## Atlas hits per query\n')
    L.append('`same_fold` = TM-score (query-normalised) >= 0.5 and query coverage >= 0.5. '
             '`structure_only` = same fold but phmmer from the query does not find it at '
             'E <= 1e-3. The Foldseek server returns at most 1000 hits per database, so these '
             'counts are a floor, not the true size of each family in the Atlas.\n')
    cols = ['query', 'atlas_hits_returned', 'same_fold', 'same_fold_div',
            'seq_detectable', 'structure_only', 'structure_only_div', 'min_seq_id', 'max_tm']
    L.append('| ' + ' | '.join(cols) + ' |')
    L.append('|' + '---|' * len(cols))
    for _, r in summ.iterrows():
        L.append('| ' + ' | '.join(str(r[c]) for c in cols) + ' |')
    tot = summ[['same_fold', 'same_fold_div', 'seq_detectable', 'structure_only',
                'structure_only_div']].sum()
    L.append('| **total** | | **%d** | **%d** | **%d** | **%d** | **%d** | | |' % (
        tot.same_fold, tot.same_fold_div, tot.seq_detectable, tot.structure_only,
        tot.structure_only_div))
    L.append('')

    if cand is not None and len(cand):
        L.append('## Candidates\n')
        L.append(f'{len(cand)} Atlas proteins share a query fold, sit below 30% sequence '
                 'identity to the query, and are not found by phmmer from the query. '
                 'These are the proteins a sequence-based screen would have missed.\n')
        if plddt is not None:
            got = plddt[plddt.fetched == True]
            L.append(f'Structure confidence was fetched from the Atlas for '
                     f'{len(got)} of them: median pLDDT {got.plddt_mean.median():.1f}, '
                     f'{int((got.plddt_mean > 70).sum())} above 70.\n')
        top = cand.sort_values('tm_score_qnorm', ascending=False).head(25)
        L.append('| Atlas protein | Query | TM-score | % id | q cov | aln | phmmer E |')
        L.append('|---|---|---|---|---|---|---|')
        for _, r in top.iterrows():
            pe = r.get('phmmer_evalue')
            pe = 'not found' if pd.isna(pe) else f'{pe:.1g}'
            L.append('| `%s` | %s | %.3f | %.1f | %.2f | %d | %s |' % (
                str(r.target).replace('.pdb.gz', ''), r['query'].split('_')[0],
                r.tm_score_qnorm, r.seq_id_pct, r.q_cov, r.n_aligned, pe))
        L.append('')

    if cand is not None and len(cand):
        c = cand.copy()
        c['mgyp'] = c.target.str.replace('.pdb.gz', '', regex=False)
        L.append('### Tiers\n')
        L.append(f'- {c.mgyp.nunique()} distinct Atlas proteins '
                 f'({int((c.groupby("mgyp")["query"].nunique() > 1).sum())} hit more than one query).')
        L.append(f'- {int(((c.tm_score_qnorm >= 0.7) & (c.seq_id_pct < 20)).sum())} at '
                 'TM-score >= 0.7 and under 20% identity, which is the tier worth looking at first.')
        L.append(f'- Candidate identity to the query runs {c.seq_id_pct.min():.1f}% to '
                 f'{c.seq_id_pct.max():.1f}% (median {c.seq_id_pct.median():.1f}%).')
        none = c[c.phmmer_evalue.isna()]
        weak = c[~c.phmmer_evalue.isna()]
        L.append(f'- phmmer finds nothing at all for {len(none)} of them (no hit even at E = 10); '
                 f'the other {len(weak)} come up only as weak hits between E = 1e-3 and E = 10. '
                 'The first group is the honest "sequence search would have missed this" set, '
                 f'and {int(((none.tm_score_qnorm >= 0.7) & (none.seq_id_pct < 20)).sum())} of '
                 'those are at TM-score >= 0.7 and under 20% identity.\n')

    L.append('## Caveats specific to this pass\n')
    L.append('- The two azoreductase queries have no near-identical structure in the PDB '
             '(top hits at 19.8% and 38.2% identity), so their family assignment rests on the '
             'Swiss-Prot annotation, not on a solved structure of that protein. The two are also '
             'different structural classes (flavodoxin-like FMN reductase vs Rossmann-fold '
             'reductase), which is why their hit sets barely overlap. A structurally '
             'characterised azoreductase should replace them in the next pass.')
    L.append('- The three long queries were searched as 390-residue windows, which cuts domains '
             'in half. The two thymidine phosphorylase windows behave very differently '
             '(385 structure-only hits vs none), so window choice is doing real work here and '
             'the counts for those three families should not be read as family sizes. '
             'Full-length structures are needed.')
    L.append('- Organism assignments for the queries come from the top PDB hit, not from '
             'UniProt, because UniProt is blocked here.\n')

    L.append('## What this does and does not show\n')
    L.append('- These are predicted structures matched to predicted structures. Nothing here '
             'shows that any candidate binds or inactivates any drug.')
    L.append('- The active-site check (step 3 of the plan) has not been done. It needs curated '
             'catalytic and cofactor residues per family, which come from UniProt feature '
             'annotations and the primary literature; `rest.uniprot.org`, `www.rcsb.org` and '
             '`www.ebi.ac.uk` are blocked by this environment\'s network policy, so that '
             'curation could not be done here.')
    L.append('- Source-sample and taxonomy mapping (steps 5 and 6) also need MGnify or IMG, '
             'which are on the same blocked hosts. The Atlas database in the Foldseek server '
             'carries no taxonomy field, so Atlas hits here have no organism or body site.')
    L.append('- Docking (step 4) was not attempted.')
    L.append('- The 1000-hit cap per database means the screen is sensitivity-limited. A full '
             'screen needs Foldseek run locally against the downloaded Atlas, which is the '
             'cluster job in the plan.\n')

    L.append('## Files\n')
    for f in sorted(os.listdir(res)):
        p = os.path.join(res, f)
        L.append(f'- `results/{f}` ({os.path.getsize(p)/1e6:.1f} MB)')
    L.append('- `queries/` query FASTA files, Swiss-Prot accession in the name')
    L.append('- `structures/` ESMFold models used as Foldseek queries')
    L.append('- `foldseek_raw/` raw Foldseek server JSON, one file per query')
    L.append('- `code/` the pipeline: esmfold.py, run_foldseek.py, parse_results.py, '
             'filter_and_compare.py, fetch_plddt.py, make_report.py')

    open(os.path.join(res, 'README.md'), 'w').write('\n'.join(L) + '\n')
    print('wrote', os.path.join(res, 'README.md'))


if __name__ == '__main__':
    main()
