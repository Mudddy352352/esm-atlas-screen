# Step 4: docking SN-38 glucuronide into the tier-1 shortlist

Run on 21 October 2026 in the project cloud environment. Every number here is
computed from the files in this directory. Nothing from pass 1 or pass 2 was
re-run; this reads `tier1_candidates.tsv`, `query_sites.tsv` and
`bodysite_candidates.tsv` and adds step 4 of the study plan.

## Headline, in two parts

**1. The docking score on its own is uninformative, and no ranking should be
read from it.** The characterised gut bacterial beta-glucuronidases score
between -8.7 and -10.1 kcal/mol. The unrelated proteins score between -7.1 and
-9.4. They overlap: an NAD(P)H azoreductase that has nothing to do with
glucuronides (-9.4) beats two of the four positive controls, and an L-tyrosine
decarboxylase (-8.72) ties the *E. coli* beta-glucuronidase crystal (-8.73) to
within four thousandths of a kcal/mol. The separation between the positive and
negative control means is 1.25 kcal/mol against a negative-control range of
2.3 kcal/mol. Run-to-run noise is small (mean SD across three seeds, 0.03
kcal/mol), so the overlap is not noise; the Vina score simply does not
distinguish a glucuronidase from any other protein with a pocket of the right
size. **The two candidates with the best scores of all, MGYP000051719573
(-10.7) and MGYP001050749229 (-10.6), beat every positive control, which is by
itself a reason to distrust the score rather than to promote them.**

**2. Where the ligand sits is informative, and it does separate the controls.**
The question a score cannot answer is whether the pocket holds the glucuronide
the way a glucuronidase has to hold it, with the sugar's anomeric end against
the catalytic carboxylates rather than the drug end wedged in somewhere
comfortable. Measuring that over all 27 poses per receptor (9 poses × 3 seeds):

| | positive controls (4) | unrelated proteins (3) | off-site box |
|---|---|---|---|
| poses with the sugar within 4 Å of the catalytic residues | 67%, 80%, 89%, 93% | 4%, 15%, 22% | 0% |
| score the ligand gives up to adopt such a pose | 0.00 in all four | 0.61, 0.70, 0.99 | n/a |

In all four positive controls the catalytically placed pose *is* the
best-scoring pose, and most poses are placed that way. In the unrelated proteins a sugar
occasionally lands near a carboxylate — they all have carboxylates — but it is
a minority of poses and it costs 0.6 to 1.0 kcal/mol relative to that
receptor's own preferred pose. So "majority of poses catalytic, at no cost in
score" is a criterion the controls pass and the negatives fail.

### What that criterion says about the 13 candidates

Four of the ten beta-glucuronidase-fold candidates behave like the positive
controls — the catalytic pose is also the best-scoring pose, and it is reached
in 42-56% of poses:

| MGYP | catalytic poses | score penalty | pass-2 active site | habitat |
|---|---|---|---|---|
| MGYP000751766417 | 52% | 0.00 | `site_partial` | marine sediment |
| MGYP000574654335 | 56% | 0.00 | `site_lost` | marine, brackish |
| MGYP000016927400 | 52% | 0.00 | `site_partial` | no metadata |
| MGYP000497510748 | 42% | 0.00 | `site_lost` | no metadata |

Taking pass 2 at its word — a candidate whose active site was scored lost is
not rescued by a favourable docking result — **two candidates are consistent
across both passes: MGYP000751766417 and MGYP000016927400.** The first is from
marine sediment in the Gulf of California, the second has no biome metadata in
any reachable MGnify release. Neither is a gut protein.

**Neither human-gut candidate produces a single catalytically placed pose.**
MGYP001533818318 and MGYP000245278657 were the two proteins the body-site
mapping placed in the human large intestine, and they are the two the
manuscript would want for an irinotecan claim. In 26-27 poses each, the closest
the glucuronide ring ever comes to their catalytic glutamates is 6.8 Å and
6.4 Å. For comparison, the positive controls reach 2.6-2.7 Å, and even the
unrelated proteins get to 3.0-3.8 Å. On this measure the two gut proteins are
worse than hen egg lysozyme.

This is not explained by model confidence. Local pLDDT at the catalytic site is
0.79 for MGYP001533818318 and 0.70 for MGYP000245278657, against 0.66 for
MGYP000051719573, which does produce catalytic poses; and the sites are no more
crowded (47-48 receptor atoms within 6 Å, against 42-53 across the set). See
`pose_scan_with_confidence.tsv`.

Three of the thirteen (MGYP000574843407, MGYP001319153113, MGYP001178215973)
are cytidine deaminase relatives, not glucuronidases. SN-38 glucuronide is not
their substrate, so for them this is a negative control rather than a test, and
the fact that MGYP001319153113 passes the geometric criterion should be read as
a false positive of the criterion, not as a result.

## What this can and cannot show

- These are **predicted structures**, including every candidate. ESMFold models
  a backbone, not a catalytically competent active site, and side-chain
  placement in a pocket is exactly what prediction is least reliable at.
- Docking is a **plausibility check, not evidence of activity.** Nothing here
  measures turnover, and a pose is not a reaction. The only way to know whether
  any of these proteins deglucuronidates SN-38G is to express it and assay it.
- The receptor is **rigid** and there is no explicit solvent, no pH, no
  cofactor and no protonation search. Several GH2 beta-glucuronidases close a
  loop over the substrate; that cannot happen here.
- The scoring function is **not calibrated for glucuronides**, which is most of
  what part 1 of the headline is about.
- The geometric criterion is **our own construction**, validated on four
  positives and four negatives. Eight control receptors is not a benchmark.
  Note also that no candidate reaches the pose share of any positive control:
  the best candidates sit at 42-56% against 67-93%, so they pass the criterion
  without matching the controls. It
  should be read as "this candidate does not fail an obvious test", not as
  "this candidate is a drug-inactivating enzyme".
- **A negative result here is the more trustworthy direction.** That the two
  gut candidates never place the sugar catalytically is a reason not to build a
  gut irinotecan claim on them. It is not proof they are inactive.

## Method

| | |
|---|---|
| Engine | AutoDock Vina 1.2.5 (static linux x86_64 build), Vina scoring function |
| Search | exhaustiveness 24, 9 poses per run, 26 Å cubic box, 4 CPU |
| Seeds | SN-38 glucuronide run at 3 seeds (20261021, 7, 999331) on every receptor; SN-38 and glucuronic acid at 1 seed |
| Runs | 118, none failed |
| Ligand prep | RDKit 2026.09.1 ETKDGv3, 20 conformers, MMFF94-optimised, lowest-energy conformer; Open Babel 3.1.1 to PDBQT at pH 7.4, so the glucuronide carboxylate is ionised |
| Receptor prep | Open Babel 3.1.1, rigid receptor, polar hydrogens at pH 7.4 |

### Ligands

| Name | ChEBI | Formula | Heavy atoms | Rotatable bonds |
|---|---|---|---|---|
| SN-38 glucuronide | CHEBI:8990 | C28H28N2O11 | 41 | 9 |
| SN-38 (the aglycone) | CHEBI:8988 | C22H20N2O5 | 29 | 4 |
| beta-D-glucuronic acid | CHEBI:28860 | C6H10O7 | 13 | 5 |

Structures were taken from the ChEBI API, not written from memory. The
glucuronide differs from SN-38 only by the beta-glucuronide on the 10-hydroxyl,
as it should.

### Receptors

22 in total.

- **13 candidates**, the tier-1 shortlist. Seven structures came from
  `screen_v1/atlas_structures/`; the other six were not saved in pass 1 and were
  fetched from the ESM Atlas API in this run, into `atlas_structures_extra/`.
- **Positive controls.** `3K4D` chain A, *E. coli* beta-glucuronidase with the
  glucaro-1,5-lactam inhibitor bound; `6D4O` chain A, *Eubacterium eligens*
  beta-glucuronidase; `6CXS` chain A, *Clostridium perfringens*
  beta-glucuronidase. All three are gut organisms and all three have a ligand in
  the catalytic site. Plus `P05804_uidA_gus_w214-603`, the ESMFold model of the
  *E. coli* enzyme that the screen used as a query, which is a positive control
  with the **same provenance and the same structure quality as the candidates**.
- **Negative controls.** `193L` hen egg lysozyme and `1TIM` chicken
  triose-phosphate isomerase; the ESMFold query models of *E. coli* cytidine
  deaminase, an NAD(P)H azoreductase and L-tyrosine decarboxylase; and
  `3K4D_A_offsite`, the *E. coli* crystal with the box moved 55 Å to a surface
  patch, which tests whether a score is really just "large protein with
  crevices".

### How the boxes were placed

For a candidate, the pass-2 alignment was walked again
(`code/map_sites.py`) to find where the query's catalytic residues land in the
candidate's **own** numbering, and the box was centred on them. The script
refuses to run unless the sequence of the downloaded structure matches the one
Foldseek aligned, which it did for all 13.

The box definition was checked against ground truth: in `3K4D`, the centroid of
the catalytic Glu413/Glu504 CA pair is 5.9 Å from the centroid of the actual
bound ligand, well inside a 26 Å box. Both definitions were used for that
receptor and give the same answer.

### Setup validation

Redocking the crystal ligand of `3K4D` into its own site recovers the
crystallographic pose at 1.65-1.68 Å (lower-bound RMSD, nearest-atom matching),
identical across three seeds. The pipeline works; part 1 of the headline is a
statement about the scoring function, not about a broken setup.

One snag worth recording: ESM Atlas and ESMFold PDB files left-justify the atom
name from column 13, so Open Babel reads a serine OG as the element Oganesson
and writes an atom type Vina rejects. All 45 runs of the first attempt failed
this way. `code/prep.py` now repairs the columns.

## Files

| File | What it holds |
|---|---|
| `final_table.tsv` | the combined table: pose geometry, scores, pass-2 site class, habitat |
| `pose_scan_sn38g.tsv` | the geometric criterion, per receptor, over all poses and seeds |
| `pose_scan_all_poses.tsv` | one row per pose: affinity and sugar-to-catalytic distance |
| `pose_scan_with_confidence.tsv` | the above plus local pLDDT and site crowding |
| `dock_summary_sn38g.tsv` | SN-38G scores per receptor, averaged over seeds, with the other two ligands alongside |
| `dock_by_receptor_ligand.tsv` | every receptor × ligand, mean/SD/min/max across seeds |
| `dock_raw.tsv` | all 118 runs, one row each, with the box and settings used |
| `control_check.tsv` | the positive-versus-negative comparison that part 1 rests on |
| `redock_check.tsv` | the pose-recovery check |
| `boxes_final.tsv` | every receptor's box centre and how it was chosen |
| `site_boxes.tsv` | the alignment walk: query site residue → candidate residue number |
| `ligands.tsv`, `ligands/` | ligand provenance, 3D SDF and PDBQT |
| `poses/` | every SN-38G pose, as Vina PDBQT, named `<receptor>__SN38G__s<seed>.pdbqt` |
| `atlas_structures_extra/` | the six candidate structures fetched in this run |
| `code/` | `map_sites.py`, `prep.py`, `run_dock.py`, `analyze.py`, `pose_scan.py` |

## Reproducing

```
python3 code/map_sites.py  <screen_v1> <query_sites.tsv> <tier1.tsv> <struct dirs...> work/site_boxes.tsv
python3 code/prep.py       <workdir>
python3 code/run_dock.py   <workdir>      # resumable; skips runs already in dock_raw.tsv
python3 code/analyze.py    <workdir> <out>
python3 code/pose_scan.py  <workdir> <out>
```

`run_dock.py` is resumable because the sandbox restarted mid-run once; it keeps
whatever is already in `dock_raw.tsv` and redoes only what is missing.

## What is still open

- Docking says nothing about turnover. An assay is the only way forward on
  MGYP000751766417 or MGYP000016927400, and neither is a gut protein.
- The two catalytic glutamates this run centres its boxes on are a GH2-wide
  feature. Keeping them, and holding a glucuronide against them in a model,
  places a protein in GH2 with an intact catalytic centre; it does not identify
  it as a glucuronidase rather than some other GH2 glycosidase.
- **Closed on 21 October:** 6D4O and 6CXS now carry mapped catalytic numbering
  (Glu425/Glu516 and Glu412/Glu505, located by global pairwise alignment to
  P05804 at 47.3% and 48.8% identity, both confirmed as glutamate), so the
  geometric criterion rests on four positives rather than two. They are the two
  strongest positives, at 93% and 89% of poses placed catalytically.
- The score and the placement are close to independent: across all 560 poses the
  correlation between affinity and sugar-to-catalytic distance is r = +0.26, and
  +0.26 within the candidates alone. A better score does not mean a better
  placed pose, which is the clearest single statement of why part 1 and part 2
  disagree.
- The score shows no glucuronide specificity either. Positive controls gain
  0.90-0.92 kcal/mol going from SN-38 to SN-38 glucuronide, but the off-site box
  gains 1.13 and the unrelated proteins 0.65-0.74, so the gain tracks ligand
  size rather than recognition.
