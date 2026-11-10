# How well does the screen actually perform?

Run 7 November 2026. This is the evaluation the manuscript has never had. The paper
reports what the pipeline found; nothing in it says how often the pipeline is right, how
often its active-site filter passes residues that mean nothing, or whether structure
similarity beats plain sequence identity at the job it claims to do.

Two questions are answered here, both using the screen's own code. Nothing is
reimplemented: `score_labelled_hits.py` and `decoy_site_null.py` import the alignment
walk, the superposition and the geometry from `screen_v2/code/active_site_check.py`, and
reuse the thresholds in `screen_v2/code/rank_candidates.py` (site CA deviation 2.0 Å,
site-distance MAD 1.5 Å, TM 0.6 for the shortlist).

## 1. Is the active-site filter better than chance?

The pass-1 searches also returned AlphaFold DB and PDB hits, but the headline question
needs no labels at all: run the same filter against **decoy site definitions** taken from
the query model itself, and see how often a candidate passes by accident. Three decoys,
increasingly hard: residues at random positions; residues matched to the real site in
number and amino-acid type; and residues matched in type and also spatially clustered, so
the geometry terms are controlled too. Ten replicates each, scored against the same 1,643
structure-only candidates.

**Full annotated site, excluding the one query window whose site is a single residue**

| Site definition | Candidate scores | Complete pass | Partial pass |
| --- | ---: | ---: | ---: |
| Real curated site | 1,694 | 0.295% (5) | 1.181% (20) |
| Clustered, composition-matched decoy | 16,864 | 0.000% (0) | 0.030% |
| Composition-matched decoy | 16,920 | 0.000% (0) | 0.035% |
| Random decoy | 16,893 | 0.006% (1) | 0.083% |

Not one of 16,864 compact, composition-matched decoy sites passed the complete-site test.
The filter is doing real work.

**Catalytic residues only, which is the rule the 13-protein shortlist actually uses**
(TM ≥ 0.6 candidates)

| Query | Site size | Real | Clustered decoy | Composed decoy | Random decoy |
| --- | ---: | ---: | ---: | ---: | ---: |
| β-glucuronidase P05804 | 2 | 62.5% | 0.83% | 0.00% | 0.43% |
| Tyrosine decarboxylase P0DTQ4 | 2 | 0.0% | 0.00% | 0.00% | 2.08% |
| Cytidine deaminase P0ABF6 | 1 | 100% | 33.3% | 55.6% | 0.0% |
| Cytidine deaminase P19079 | 1 | 80% | 46.8% | 22.0% | 11.1% |

For the glucuronidase pair of glutamates the rule is about 75 times better than the
hardest null, so among the ten GUS-fold shortlisted proteins fewer than one is expected
to pass by chance. **For the two queries whose UniProt record annotates a single catalytic
residue the rule is close to meaningless**: a decoy residue passes a fifth to a half of
the time. Three of the thirteen shortlisted proteins are cytidine deaminase relatives and
inherit that weakness, and the manuscript should say so.

## 2. Can the pipeline tell the right enzyme from its same-fold neighbours?

The Atlas carries no functional annotation, which is why the screen could never measure
this. The same searches hit AlphaFold DB, and 1,684 of those hits carry a UniProt EC
number, which turns the search into a labelled retrieval problem: among proteins with the
query's fold, find the ones with the query's activity. Positives are hits whose EC matches
the query's at the query's own level of specificity. Negatives are hits with another EC,
most of them same-family enzymes acting on a different sugar, which is the hard case.

For the glucuronidase query the EC-labelled set is 74 true beta-glucuronidases against 195
other glycoside hydrolases, mostly beta-galactosidases and beta-mannosidases.

| Measure | Whole EC-labelled set | Below 30% identity |
| --- | ---: | ---: |
| AUROC, TM-score | 0.810 | 1.000 |
| AUROC, sequence identity | 0.998 | 0.943 |
| AUROC, site retention | 0.712 | 0.851 |
| Screen rule precision | 0.770 | 0.176 |
| Screen rule recall | 0.635 | 1.000 |
| Pfam same-family precision | 0.288 | 0.016 |

Only 3 of the 74 EC-labelled positives sit below 30% identity, so that column rests on
three proteins and is reported for completeness, not as an estimate.

A wider set labels hits by their UniProt protein name instead of a curated EC. Two kinds
of entry carry no label either way and are **excluded rather than counted as negatives**:
records UniProt has deleted (1,116 of the low-identity targets) and names that state no
activity at all ("uncharacterized protein", "glycoside hydrolase family 2",
DUF and domain-containing names). That leaves 107 positives against 251 annotated
negatives, 12 against 247 of them below 30% identity.

### Against the obvious simpler methods

The question a reviewer asks is not what AUROC each feature reaches but whether the
pipeline beats thresholding one number. Each baseline below is given the threshold that
matches the screen rule's own recall, so the precisions compare like for like. The Pfam
baseline has nothing to tune and is shown at its own operating point
(`results/baseline_comparison.tsv`, from `code/baseline_comparison.py`).

| Method | EC set, all | Name set, all | Name set, below 30% |
| --- | ---: | ---: | ---: |
| Screen rule | 0.770 | 0.686 | 0.156 |
| Sequence identity, matched recall | 1.000 | 0.967 | 0.106 |
| TM-score, matched recall | 1.000 | 0.983 | 0.294 |
| Site retention, matched recall | 0.456 | 0.492 | 0.088 |
| Pfam same family | 0.288 | 0.309 | 0.034 |
| Recall all methods are held to | 0.635 | 0.551 | 0.417 |

### What this says, including the parts that do not flatter the paper

- **Above 30% identity the pipeline is accurate and unnecessary.** Precision 1.00 on the
  name-labelled band at and above 30%, and sequence identity alone separates
  glucuronidases from their relatives almost perfectly (AUROC 0.998 on the EC set).
- **Where the labels are good, a single sequence-identity threshold matches or beats the
  screen rule.** Precision 1.000 against 0.770 on the EC-labelled set at the same recall.
  The screen's extra machinery buys nothing in the regime where annotation exists.
- **Below 30% identity, where all 1,643 candidates sit, every method is poor and TM-score
  is the least poor**: 0.294 against 0.106 for sequence identity and 0.034 for Pfam, at
  recall 0.417. Structure does beat sequence identity in the remote regime, which is the
  paper's claim, but on 12 positives, and the absolute precision means roughly three in
  ten flagged proteins would be real. Applied to the ten GUS-fold shortlisted proteins
  that is an expectation of about three real ones.
- **The Pfam baseline cannot discriminate**, as expected: precision 0.288 at recall 1.0 on
  the EC set, because Glyco_hydro_2 covers glucuronidases and galactosidases alike. That
  part of the manuscript's argument survives intact.
- **Site retention on its own is the weakest ranking of the three.** It earns its place in
  the pipeline through the decoy null in part 1, as a specificity filter, not as a way of
  ranking candidates.

### A correction to an earlier version of this file

An earlier run of the name-labelled analysis reported 85 positives against 1,373
negatives below 30% identity, with screen precision 0.256 and TM AUROC 0.742 against
sequence identity 0.797, and concluded that structure does not beat sequence identity in
the remote regime. That construction could not be reproduced from the stored data, and the
only constructions that come close reach those counts by scoring deleted and
unannotated UniProt entries as negatives. Those figures are withdrawn. The table above
replaces them and is produced by a script in `code/`.

## Caveats

- EC numbers in AlphaFold DB entries are mostly automatic annotation: of 9,875 targets
  with a UniProt record, 1,684 carry an EC and only 10 are reviewed. The name-labelled
  analysis is noisier still, since "beta-glucuronidase" in a TrEMBL name is often
  inherited from a homolog.
- The low-identity band is the one that matters and the one with fewest positives. Three
  EC-labelled positives is not a benchmark, and the 12 name-labelled ones are a weak
  estimate: the confidence interval on a precision of 0.156 from 12 positives is wide, and
  nothing here pins it down further.
- Excluding deleted and activity-free UniProt names is the right call but it is not free.
  Of 1,690 low-identity targets with any UniProt record, 1,332 are excluded, so the set
  that remains is the annotated minority and may not represent the fold's neighbourhood.
- Hits come from a search capped at 1,000 per database, so the labelled set is a sample of
  the fold's neighbourhood and not a census.
- Decoy sites can only be drawn from the query model, so the null answers "would these
  thresholds pass meaningless residues", not "would another method do better".

## Files

| File | What it holds |
| --- | --- |
| `code/score_labelled_hits.py` | Scores the AlphaFold DB and PDB hits with the screen's own check |
| `code/fetch_uniprot_labels.py` | EC numbers, Pfam and review status for the AlphaFold targets, cached under `labels/` |
| `code/evaluate_retrieval.py` | AUROC, precision and recall per query, and the Pfam baseline |
| `code/decoy_site_null.py` | The decoy-site null model, all three decoy types |
| `results/labelled_hits_scored.tsv` | 17,320 annotated hits with site scores |
| `results/uniprot_labels.tsv` | 9,875 UniProt records for the AlphaFold targets |
| `code/baseline_comparison.py` | The matched-recall comparison against identity, TM, site retention and Pfam |
| `results/retrieval_by_query.tsv`, `retrieval_by_identity_band.tsv` | Per-query retrieval, and the EC-labelled identity bands |
| `results/baseline_comparison.tsv`, `name_labelled_gus.tsv` | The baseline table above and the set it is computed on |
| `results/decoy_site_null*.tsv` | The null model, per query and pooled |
