# Evaluating the screen (manuscript section 3.10)

This directory holds the parts of the section 3.10 evaluation that run offline from
the files deposited in this repository. It was added after release v1.0.0, when the
section 3.10 analyses turned out not to have been committed with the rest.

## What is here

### `threshold_sweep.py` — Supplementary Note S4

How much of the shortlist is the data and how much is the cut-off. Sweeps the three
geometric thresholds of section 2.6 across TM-score 0.5–0.7, site RMSD 1.5–3.0 Å and
site geometry MAD 1.0–2.5 Å, five values each, and reports the shortlist at every
combination.

Runs from `screen_v1/results/candidates.tsv` and
`screen_v2/results/active_site_all_hits.tsv`. No network, no large intermediates.

```
python3 screen_v2/evaluation/threshold_sweep.py
```

It verifies itself against the manuscript and against the deposited shortlist, and
exits non-zero if any reported number fails to reproduce. As of this commit all of
them do:

| | Reproduced | Section 3.10 |
|---|---|---|
| Candidates keeping every annotated catalytic residue | 30 | 30 |
| Inside the geometric thresholds | 17 | 17 |
| Also TM-score ≥ 0.6 | 15 | 15 |
| Also phmmer-undetected — the shortlist | 13 | 13 |
| Threshold combinations swept | 125 | 125 |
| Shortlist size across the sweep | 3 to 15 | 3 to 15 |
| Of the 13, shortlisted at every setting | 3 | 3 |

The rebuilt shortlist is identical, protein for protein, to `tier1_candidates.tsv`.

The MAD criterion never binds at any value in the swept range, and site RMSD stops
binding above 2.0 Å, as section 3.10 states. The three proteins shortlisted at every
setting are MGYP000574654335, MGYP000574843407 and MGYP001319153113.

Outputs: `threshold_sweep_grid.tsv` (one row per threshold combination) and
`threshold_sweep_per_protein.tsv` (how many settings each shortlisted protein
survives). The per-protein fractions are what section 4.1 quotes: MGYP000751766417
and MGYP000016927400 at four fifths of settings, MGYP000051719573 at 48%.

## What is not here, and why

Two of the section 3.10 analyses cannot be rebuilt from this repository. Neither is a
matter of effort; both need an input that is not deposited.

- **The decoy-site null** (Supplementary Tables S8 and S9) scores decoy site
  definitions against the candidates' own coordinates, which live in the per-residue
  `tCa` fields of the raw pass-1 Foldseek JSON — the ~250 MB intermediate that is not
  shipped (see `EXCLUDED.md` §1). Regenerating that JSON means re-running the pass-1
  search against the Foldseek web server.
- **The 62-family M-CSA comparison** (Figure 7, Supplementary Table S10) needs fresh
  Foldseek searches of 62 Mechanism and Catalytic Site Atlas enzymes against
  AlphaFold DB, and the M-CSA catalytic-residue curation, neither of which is a local
  computation.

The **per-pose contact analysis** behind Supplementary Figures S2 and S3 is partly
recoverable: the 72 pose files, the 13 candidate structures and the catalytic residue
numbers are all deposited, so candidate contacts can be recomputed, but the crystal
control receptors (PDB 3K4D, 6D4O, 6CXS, 193L, 1TIM) are not vendored, so the control
comparison that validates the contact reading needs them fetched from the RCSB PDB
first.
