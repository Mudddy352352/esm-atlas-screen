# Does a shared fold carry a working active site? Asked of 62 enzyme families

Run 7 November 2026. The screen reported, as a side result, that 83.9% of same-fold hits a
sequence search also finds keep their catalytic residues, against 42.3% of the hits only
structure search reaches. That rested on nine query enzymes and on UniProt's catalytic-site
annotation, which for most of those queries names one or two residues. This repeats the same
measurement on 62 families from the Mechanism and Catalytic Site Atlas, whose catalytic
residues are curated from the literature, and on 49,735 same-fold hits.

Nothing is reimplemented. `code/mcsa_retention_scan.py` imports the alignment walk, the
superposition and the geometry from `screen_v2/code/active_site_check.py` and applies the
screen's own thresholds: site CA deviation 2.0 A, site-distance MAD 1.5 A, TM-score 0.5 or
better for a hit to count as same-fold.

## Method

For each M-CSA entry: take the curated catalytic residues of its reference PDB chain, search
AlphaFold DB with Foldseek in 3Di+AA mode, and for every hit above TM 0.5 measure whether all
of those residues are aligned, identical and still correctly placed. A hit counts as
sequence-detectable if a single-sequence phmmer search of the query against the hit sequences
finds it at E <= 1e-3, and structure-only otherwise. Entries with fewer than two curated
residues, or whose reference chain could not be matched, are skipped and logged in
`results/done.txt`.

## Result

| Group | Hits | Site kept | 95% CI |
| --- | ---: | ---: | ---: |
| Sequence-detectable | 44,120 | 37.4% | 37.0 - 37.9% |
| Structure-only | 5,615 | 10.2% | 9.4 - 11.0% |
| All same-fold hits | 49,735 | 34.4% | 33.9 - 34.8% |

Pooling over hits lets a large family dominate, so the comparison is also made family by
family. Of the 52 families with hits in both groups, 45 retain the site more often among the
sequence-detectable hits; the median within-family difference is 16.4 percentage points
(Wilcoxon signed-rank on the paired family rates, p = 1.3e-07).

**The confound, and whether it explains the result.** Structure-only hits are remote by
construction, so the gap could be a sequence-identity effect and nothing more. It is not:
inside the band below 30% identity, where 22,922 sequence-detectable and 5,605 structure-only
hits sit together, retention is 24.8% against 10.2%. The gap narrows from about 3.7-fold to
about 2.4-fold and survives.

| Identity to query | Sequence-detectable | Structure-only |
| --- | ---: | ---: |
| Below 30% | 24.8% (22,922) | 10.2% (5,605) |
| 30 to 50% | 50.2% (16,899) | 0.0% (10) |
| 50% and above | 54.9% (4,299) | - |

Retention also tracks identity on its own: 21.9% below 30%, 50.2% from 30 to 50%, 54.9% above
50%.

## How this relates to the paper's 83.9% against 42.3%

The direction is the same and now rests on 62 families rather than nine. The absolute levels
are much lower, and the reason is the annotation, not the method: M-CSA curates a median of
five or six catalytic residues per enzyme, while UniProt names one or two for most of the
screen's queries, and the test requires every annotated residue to survive. Retention falls as
the site grows, from 56.4% at two residues to about 30% at eight. So the two numbers answer the
same question at different stringencies, and the general statement to make is the ratio and its
direction, not the percentage.

## Caveats

- Hits come from a search capped at 1,000 per query, so each family's set is a sample of its
  fold neighbourhood, not a census.
- AlphaFold DB models, not experimental structures, on both sides of the comparison. A predicted
  side chain in the right place is weaker evidence than a crystal structure.
- Sequence-detectability is defined by one phmmer search from one query sequence at E <= 1e-3.
  A profile method would find more of the structure-only hits, which would shrink the gap.
- 62 families is a sample of M-CSA's roughly 1,000, taken in entry order and not chosen at
  random, so it is not a stratified sample of enzyme chemistry.
- Two families account for the extremes of the site-size table at n = 9 and n = 10 residues, so
  those rows are family effects rather than a size effect.

## Files

| File | What it holds |
| --- | --- |
| `code/mcsa_retention_scan.py` | The resumable scan: M-CSA, RCSB, Foldseek, phmmer, the screen's geometry |
| `code/analyse_retention.py` | Pooling, the per-family pairing, the bands and the summary |
| `results/mcsa_hits.tsv` | 49,735 scored hits, one row each |
| `results/retention_by_family.tsv` | One row per family, both groups |
| `results/retention_by_identity_band.tsv`, `retention_band_by_detectability.tsv` | The band tables above |
| `results/RETENTION_SUMMARY.txt` | The printed summary as run |
| `results/done.txt` | Every entry attempted and why it was skipped |
