# Pfam-profile comparator

Run on 21 October 2026. This answers the question left open by passes 1 and 2:
the screen compared structure search against phmmer from a single query sequence,
but the study plan called for a comparison against profile methods (BLAST and
Pfam HMMs). phmmer from one sequence is a weak comparator, so a profile search
could recover proteins the screen calls "structure-only". This step runs that
comparison.

## What was run

`hmmscan` (pyhmmer 0.12.3) against the full Pfam-A HMM library
(Pfam-A.hmm, 30,134 models, current release from the EBI FTP server), at Pfam
**gathering thresholds** (the per-family bit-score cutoffs Pfam curators set, the
standard way to call Pfam domains), over:

- the 9 full-length query enzymes, and
- the 1,640 distinct structure-only candidate proteins from pass 1.

A candidate "recovers the link" when it shares at least one Pfam family with its
own query enzyme. Sequences for the candidates were taken from the `tSeq` field of
the pass-1 Foldseek JSON, the same sequences that were screened. Query sequences
are the Swiss-Prot query FASTA from pass 1.

Inputs `queries.faa`, `candidates.faa`. Script `pfam_scan.py`. Outputs
`pfam_query_domains.tsv`, `pfam_candidate_domains.tsv`, `pfam_recovery.tsv` (one
row per candidate x query pair with the shared families and the verdict),
`pfam_summary.txt`.

## Result: a Pfam profile search recovers most of the candidates

| | n | % |
|---|---|---|
| candidate (query, protein) pairs | 1,643 | |
| pairs where the candidate carries any Pfam domain | 1,613 | 98.2 |
| pairs where the candidate shares a Pfam family with its query | 1,519 | 92.5 |
| distinct proteins sharing a Pfam family with any query | 1,516 | 92.4 |

So the 1,643 are not proteins a profile search would miss. A Pfam HMM search
assigns the great majority of them to the same Pfam family as the
drug-metabolising query. The screen's phmmer-from-one-sequence comparator
overstated how invisible these proteins are to sequence methods, and the
manuscript has to say so.

Recovery varies by family, and the variation is informative:

| Query | Candidates | Recovered by Pfam | % | Query's Pfam families |
|---|---|---|---|---|
| P07650 thymidine phosphorylase (w1-390) | 385 | 384 | 99.7 | Glycos_trans_3N, Glycos_transf_3, PYNP_C |
| P17117 nitroreductase NfsA | 318 | 318 | 100 | Nitroreductase, Nitroreductase_C |
| P15888 nitroreductase NfsB | 296 | 295 | 99.7 | Nitroreductase |
| Q8KU07 azoreductase | 255 | 243 | 95.3 | Epimerase, NAD_binding_10, NmrA, TrkA_N |
| P05804 beta-glucuronidase (w1-390) | 55 | 54 | 98.2 | Glyco_hydro_2 (+_C, _N, _N2) |
| P0DTQ4 tyrosine decarboxylase (w231-620) | 36 | 21 | 58.3 | Pyridoxal_deC, tyr_de_CO2_C |
| P05804 beta-glucuronidase (w214-603) | 28 | 26 | 92.9 | Glyco_hydro_2, Glyco_hydro_2_C |
| Q7A782 azoreductase | 261 | 174 | 66.7 | FMN_red, Flavodoxin_2 |
| P0DTQ4 tyrosine decarboxylase (w1-390) | 2 | 2 | 100 | |
| P19079 cytidine deaminase (tetrameric) | 5 | 1 | 20.0 | dCMP_cyt_deam_1, dCMP_cyt_deam_2 |
| P0ABF6 cytidine deaminase | 2 | 1 | 50.0 | dCMP_cyt_deam_1, dCMP_cyt_deam_2 |

The cytidine deaminase family, the gemcitabine-relevant one, is the exception:
Pfam profiles recover only 2 of its 7 structure-only hits. For that family,
structure search does reach proteins a Pfam search does not. For the
irinotecan-relevant GH2 family and the reductase families, Pfam recovers almost
everything.

## The 13-protein shortlist

11 of the 13 are recovered by a Pfam search. Both exceptions are cytidine
deaminase relatives from the tetrameric query, MGYP000574843407 and
MGYP001319153113, for which hmmscan finds no deaminase family (nor any Pfam
domain) at gathering thresholds, and both were scored `site_lost` overall in pass
2. The ten GH2 relatives and the one P0ABF6 cytidine deaminase relative all carry
the query's Pfam families (Glyco_hydro_2 / Glyco_hydro_2_C, or dCMP_cyt_deam).

## What this does and does not change for the manuscript

- "1,643 proteins sequence search misses" is wrong, and so is any softer version
  that keeps the word "invisible". A profile search finds 92.5% of them.
- The Pfam family it assigns is broad. Glyco_hydro_2 covers glucuronidases,
  galactosidases and mannosidases alike, so a Pfam hit places a protein in GH2
  and does not identify it as a glucuronidase. This matches the taxonomy result:
  the nearest labelled relatives of the GH2 shortlist are generic GH2 enzymes.
- So the screen's contribution is not reach into unannotated space. It is
  prioritisation within families that annotation already identifies: ranking Atlas
  proteins by 3D similarity to a specific drug-metabolising structure, and then
  filtering on active-site geometry. A Pfam label gives neither the ranking nor
  the active-site check.
- The one place structure search adds genuine reach is the cytidine deaminase
  family, where Pfam profiles recover 2 of 7. That is worth stating, with the
  caveat that the two shortlisted deaminase relatives Pfam misses have both lost
  the site.

## Caveats

- Shared Pfam family is a generous definition of "recovered". It counts a shared
  broad family (Glyco_hydro_2), not identification as the specific drug-relevant
  enzyme. Read this way, recovery is a lower bound on how much a profile search
  explains, which makes the conclusion stronger not weaker.
- Gathering thresholds are the curated Pfam cutoffs. A more permissive E-value
  cutoff would recover even more, so 92.5% is conservative.
- This compares against Pfam only, not against a full BLAST of a large sequence
  database, which was in the plan and is not done here. BLAST would likely recover
  a similar or larger fraction, since the taxonomy step already found labelled
  AFDB50 relatives at a median 47.7% identity for the shortlist.
