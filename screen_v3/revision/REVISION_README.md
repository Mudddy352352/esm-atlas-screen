# Answering the internal peer review of v1.6

Run 9 November 2026. `manuscript/reviews/peer_review_v1_6.md` holds three
in-character referee reports on draft v1.6. Everything in them that could be
settled by measurement rather than by rewording was settled here, in
`code/`, with the results in `results/`. Nothing was transcribed from the
draft: every number below comes from the screen's own output files.

Run order does not matter; the scripts are independent.

```
python3 code/stats_for_revision.py          # intervals, Fisher test, attrition, pLDDT
python3 ../benchmark/code/baseline_comparison.py ../benchmark --sensitivity
python3 code/threshold_sensitivity.py       # the three cut-offs, one at a time and together
python3 code/side_chain_check.py            # side-chain geometry at the curated sites
python3 code/specificity_motif.py           # the residues that read the C6 carboxylate
python3 code/pocket_enclosure.py            # is the site open or closed
python3 code/site_loops.py                  # is there a loop there to close
python3 code/ensemble_contacts.py           # contacts over 27 poses, not one
```

`code/pocket_enclosure.py` needs the four control crystals, which are in
`structures/` (downloaded from RCSB: 3K4D, 6CXS, 193L, 1TIM); 6D4O is reused
from `manuscript/figures/pose_inputs/`.

## Section 3.4 reproduces, and now carries a test

`results/section_3_4_by_query.tsv`. Recomputed from
`screen_v2/results/active_site_all_hits_joined.tsv.gz` with the screen's own
definitions: 1,882 of 2,243 sequence-detectable hits keep every annotated
catalytic residue (83.9%, 95% CI 82.3 to 85.4) against 30 of 71 structure-only
hits (42.3%, 95% CI 31.5 to 53.8). Fisher exact two-sided p = 7.0e-15, odds
ratio 7.1. Four query windows contribute, and all four go the same way.

## Every precision in Table 16 now has an interval, and the choice of label set is tested

`baseline_comparison.py` gained the number of calls, the number of true calls
and a Wilson interval on every row, and a `--sensitivity` mode that repeats the
name-labelled comparison under four rules for which UniProt names carry no
label (`results/baseline_sensitivity.tsv` in the benchmark directory).

**This changes what the paper may claim.** Below 30% identity the screen rule's
precision is 0.156 (0.069 to 0.318) over 32 calls and the TM-score baseline's is
0.294 (0.133 to 0.531) over 17. The intervals overlap heavily, so "structure
beats sequence identity in the remote regime" is a direction the data support,
not a separation the data establish. The ordering of the four methods is stable
under all four exclusion rules; the absolute precisions move by up to four
points (screen rule 0.114 to 0.156).

The expectation for the shortlist follows the rule, not the baseline: at 0.156,
about 1.6 of the ten GUS-fold shortlisted proteins would be real, not three.
v1.6 used the TM-score figure in two places and read in our favour.

## The thresholds matter, and most of the shortlist does not survive moving them

`results/threshold_sensitivity.tsv`, `results/shortlist_stability.tsv`. Sweeping
TM-score 0.5 to 0.7, site RMSD 1.5 to 3.0 A and site MAD 1.0 to 2.5 A:

* the shortlist runs from 3 proteins (all three strict) to 15 (all three loose);
* the MAD cut-off never binds on the candidate set, at any value in that range;
* site RMSD stops binding above 2.0 A, so the published value sits at the edge
  of where it has any effect;
* across all 125 threshold combinations only 3 of the 13 are shortlisted every
  time, and 2 proteins outside the published 13 appear under looser settings;
* benchmark precision below 30% identity rises from 0.156 to 0.25 as the
  cut-offs tighten, at unchanged recall.

## The side chains, not just the backbone

`results/side_chain_check.tsv`. The screen's site test is a register and
backbone test: same amino acid in the alignment, CA close after superposition,
site residues the same distances apart. For the thirteen, whose full-atom models
are on disk, this measures the distance between the functional atoms of the two
side chains (carboxylate midpoint for Asp and Glu, OH for Tyr, NZ for Lys, ring
centroid for Trp) and the chi-1 torsions.

The catalytic residues hold up: 23 catalytic positions, functional atoms a
median 1.54 A apart, worst case 2.97 A. Over all 86 site positions the median is
1.72 A and the worst 9.98 A, and in 15 of 86 the CA is within 2 A while the
functional atoms are more than 3 A apart. Those 15 are at substituted
binding residues, not at the catalytic pair, so the backbone criterion was not
hiding a flipped carboxylate where it matters. The text still has to say that
the published test was backbone-only.

## The glucuronide-specificity residues are gone in all ten

`results/specificity_motif.tsv`, `results/specificity_motif_wide.tsv`. What
separates a GH2 beta-glucuronidase from a GH2 beta-galactosidase is not the
catalytic glutamate pair, which both have, but the residues that stabilise the
carboxyl group at C6 of the sugar: the Asn-x-Lys motif in loop 8, an arginine
upstream in the same loop, and a tyrosine in loop 6 (Lebreton et al. 2025,
Biotechnol Biofuels Bioprod 18:68, doi:10.1186/s13068-025-02669-8, and the
structures cited there). In the E. coli query P05804 those are Asn566, Lys568,
Arg562 and Tyr472, of which UniProt annotates only Tyr472 and Lys568, so the
screen never asked this question separately.

Not one of the ten GUS-fold shortlisted proteins keeps all four. Lys568 and
Tyr472 are kept by six each, Asn566 by one, and **Arg562 by none** — in every
one of the ten it falls in a 5 to 11 residue stretch of loop 8 that the
candidate does not have at all. Five of the ten keep two of the four and five
keep one or none.

This is the strongest structural support the paper has for its central negative
claim, and it is new: the candidates are GH2 proteins with an intact catalytic
centre and without the apparatus for reading a glucuronide.

## The loop-closure explanation does not hold

`results/pocket_enclosure.tsv`, `results/site_loops.tsv`. The referee's
alternative reading of section 3.9 is that the two gut-sampled proteins fail on
pose because a loop that real GH2 glucuronidases close over the substrate is
absent or open in a rigid predicted model. Two measurements, both against the
same control set the docking used:

* **Enclosure.** From the most open point near each docking box centre, 1,000
  rays are cast and the fraction that run into protein within 15 A is counted.
  The positive controls score 0.616 to 0.965 (median 0.879); the candidates
  score 0.808 to 1.000 (median 0.947). The two gut proteins are 1.000 and
  0.940, both above every positive control crystal except 3K4D. Their sites are
  not more open. Enclosure does not predict the pose fraction
  (Spearman rho = -0.14, p = 0.59, n = 18).
* **Is there a loop at all.** Every insertion of five residues or more within
  12 A of the site, in either direction. Ten of thirteen candidates are
  *shorter* than the query there, by 5 to 14 residues; neither gut protein
  carries any extra segment near its site. Three candidates have a
  candidate-only segment, and the largest of those, 157 residues in
  MGYP000574843407, is where the alignment starts late rather than a loop.

So there is no loop in these models to close, and the open-pocket reading is not
supported. What the measurements cannot rule out is that a flexible loop in the
real protein would reposition the ligand; that limit belongs in the text, and it
applies to the controls equally.

## Contacts over the ensemble rather than the best pose

`results/ensemble_contacts.tsv`, `results/ensemble_contacts_per_pose.tsv`. The
same contact code as Figures 9 and 10, run over all 27 poses of each of the four
figure receptors, with each frequency counting poses rather than atom-level
contacts.

The ensemble reading is stronger than the single pose, in the same direction:

* 6D4O, the positive control: both catalytic glutamates hydrogen-bond the sugar
  in 27 of 27 poses; Lys580 makes its salt bridge in 67% and Arg574 in 52%.
* MGYP000751766417: Glu237 in 89% of poses, Glu338 in 7%. Half a site, in every
  pose, not just the best one.
* MGYP000016927400: one catalytic glutamate in 41% of poses, the pocket held
  mostly by Phe141 and Trp366 hydrophobics (85% each).
* MGYP001533818318, from the human large intestine: **neither catalytic residue
  contacts the ligand in any of the 27 poses**, while Arg122 hydrogen-bonds it
  in 27 of 27. The protein binds the drug firmly and by the wrong residues, and
  that was not an artefact of picking one pose.

## Other numbers the referees asked for

From `stats_for_revision.py`:

* **Attrition from 30 to 13.** 30 candidates keep every annotated catalytic
  residue; 17 of those are inside the geometric thresholds; 15 of those are at
  TM-score >= 0.6; 13 of those are phmmer-blind at any threshold.
* **Model confidence for the thirteen** (`results/shortlist_plddt.tsv`):
  whole-model pLDDT median 86.5, range 75.4 to 91.5 (7 proteins have it);
  local confidence at the site, median 0.81, range 0.66 to 0.92, all thirteen.
* **The pose criterion's own false positive rate.** MGYP001319153113, a cytidine
  deaminase relative with no glucuronidase site at all, reaches 0.59 of poses
  catalytic, inside the positive-control range of 0.67 to 0.93 only just below
  it, and above two negative controls. The criterion is not specific on its own.
* **Table 15, tyrosine decarboxylase.** The real site passes 0.0% and the random
  decoy 2.08%. Both are small counts on a two-residue site: 0 of 24 real against
  7 of 332 random. Nothing follows from the real rate being lower, and the draft
  should say so rather than leave it unexplained.

## References checked against Europe PMC

* Reference 24, Takasuna et al. 1996 Cancer Res 56:3752-3757, PMID 8706020:
  authors are Takasuna K, Hagiwara T, Hirohashi M, Kato M, Nomura M, Nagai E,
  Yokoi T, Kamataki T. The draft's spelling is right.
* Reference 26, Pollet et al. 2017: full author list is Pollet RM, D'Agostino EH,
  Walton WG, Xu Y, Little MS, Biernat KA, Pellock SJ, Patterson LM, Creekmore BC,
  Isenberg HN, Bahethi RR, Bhatt AP, Liu J, Gharaibeh RZ, Redinbo MR.
* Reference 27, Bhatt et al. 2020: Bhatt AP, Pellock SJ, Biernat KA, Walton WG,
  Wallace BD, Creekmore BC, Letertre MM, Swann JR, Wilson ID, Roques JR, Darr DB,
  Bailey ST, Montgomery SA, Roach JM, Azcarate-Peril MA, Sartor RB, Gharaibeh RZ,
  Bultman SJ, Redinbo MR.
* The GH2 shared-mechanism placeholder in section 3.5 is answered by Lebreton
  et al. 2025 (above), which states the catalytic pair is common to the family
  and that the NxK motif, the loop-8 arginine and the loop-6 tyrosine are what
  distinguish the beta-hexuronidases. Pellock et al. 2018 (JBC 293:18559-18573,
  doi:10.1074/jbc.RA118.005414) is the matching structural reference: nine
  crystal structures of three GUS enzymes from one gut commensal, reporting that
  active-site features of this family are hard to predict by structural
  modelling, which is also the honest citation for our own limitation.
