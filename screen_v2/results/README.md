# Pass 2: active-site check on the pass-1 candidates

Run on 21 October 2026 in the project cloud environment. Every number below is
computed from the files in this directory and from `screen_v1/`. The Foldseek
screen was **not** re-run; pass 2 reads `screen_v1/results` and adds step 3 of
the study plan.

## What was and was not possible

| Plan step | Status | Why |
|---|---|---|
| 3. Active-site check | **done** | `rest.uniprot.org` is reachable from this environment |
| 5. Body-site mapping | **done on 21 Oct, after the environment was widened** | see `BODYSITE_README.md`. 2 of the 13 are human faecal, 5 environmental/industrial, 6 have no metadata in any reachable release |
| 6. Taxonomy of Atlas hits | **done on 21 Oct, as nearest-relative inference only** | see `TAXONOMY_README.md`. MGnify assigns no taxonomy to MGYP accessions at all, so the source organism is not recoverable; what is reported is the taxonomy of each candidate's closest characterised relatives |
| 4. Docking | not attempted | out of scope for this pass |

Host check made at the start of this run, by direct request:

| Host | Result |
|---|---|
| `rest.uniprot.org` | reachable at 13:20 UTC, **refused at 13:44 UTC** after the environment was switched (see the rerun note below) |
| `rcsb.org`, `ebi.ac.uk` | reachable, but both only redirect to their `www` host |
| `www.rcsb.org`, `files.rcsb.org`, `data.rcsb.org`, `search.rcsb.org` | refused |
| `www.ebi.ac.uk`, `alphafold.ebi.ac.uk`, `ftp.ebi.ac.uk`, `api.mgnify.org` | refused |
| `esmatlas.com`, `api.esmatlas.com`, `search.foldseek.com`, `afesm.foldseek.com` | refused (these worked on 20 October, so pass 1 cannot be repeated or extended in this environment as it now stands) |

The pattern is that allowed entries cover the exact hostname only, not its
subdomains, so the `www.` and `api.` hosts each have to be listed.

### Rerun on 21 October, 13:44 UTC

All three steps were rerun after the project's cloud environment was switched.
Two things came out of it:

- Every host above is now refused, `rest.uniprot.org` included, so the
  active-site step can no longer be run against UniProt from this environment.
  The network policy is applied live rather than being fixed at session start,
  so a switch can take away access that a session already had.
- The step is still reproducible offline, because the nine UniProt JSON records
  it needs were saved to `../sites/` on the first run. Rerunning
  `build_sites.py`, `active_site_check.py` and `rank_candidates.py` from those
  cached records reproduced all nine output files byte for byte. The numbers in
  this README are therefore confirmed by a second, independent execution.

## How the query active sites were obtained

Catalytic and cofactor residues come from UniProtKB feature annotations
(`Active site`, `Binding site`, and the pyridoxal-phosphate `Modified residue`),
fetched per accession in this run. Seven of the nine queries carry such features
themselves. The two that do not were filled from a reviewed, experimentally
annotated member of the same family, carried across a global pairwise alignment
(BLOSUM62), keeping only positions where the aligned residue is identical:

| Query | Source | Residues |
|---|---|---|
| P05804 beta-glucuronidase | own UniProt record | 8 (2 catalytic, 6 binding) |
| P0ABF6 cytidine deaminase | own UniProt record | 7 (1 catalytic, 6 binding) |
| P0DTQ4 L-tyrosine decarboxylase | own UniProt record | 9 (2 catalytic incl. PLP-Lys, 7 binding) |
| P15888 nitroreductase NfsB | own UniProt record | 11 (0 catalytic, 11 binding) |
| P17117 nitroreductase NfsA | own UniProt record | 13 (0 catalytic, 13 binding) |
| P19079 cytidine deaminase (tetrameric) | own UniProt record | 7 (1 catalytic, 6 binding) |
| Q8KU07 NAD(P)H azoreductase | own UniProt record | 13 (0 catalytic, 13 binding) |
| P07650 thymidine phosphorylase | transferred from P77836 (*Geobacillus* pyrimidine-nucleoside phosphorylase, 42.4% id) | 10 kept, 5 dropped |
| Q7A782 FMN-dependent NADPH-azoreductase | transferred from O07529 (*B. subtilis* azoreductase, 32.0% id) | 3 kept, 7 dropped |

Residue-level provenance is in `query_sites.tsv`; the per-query summary is in
`query_sites_provenance.tsv`.

**This is the weak point of pass 2.** UniProt annotates no catalytic residue at
all for the nitroreductases, the azoreductases or thymidine phosphorylase, only
cofactor and substrate contacts. Those four families account for 1,572 of the
1,643 candidates, so for 96% of the candidate list the check tests
cofactor binding, not catalysis. Fixing that needs catalytic residues from solved
structures and the primary literature (M-CSA, PDB), both on blocked hosts.

## How the check was run

For each of the 12,000 Atlas alignments (12 query models × 1,000 hits) the query
site residues were pushed through the Foldseek structural alignment and asked
three separate questions:

1. **Alignment register.** Is the site residue aligned, and is it the same amino
   acid?
2. **Superposed deviation.** After superposing the candidate on the query using
   the whole alignment (TM-score style iterative refinement), how far is each
   candidate site CA from the query site CA? (`site_rmsd`, `site_max_dev`)
3. **Rotation-invariant geometry.** Do the site residues sit the same distances
   apart in the candidate as in the query? (`site_dist_mad`, the mean absolute
   difference of all pairwise site-site CA distances) This cannot be flattered
   by a good global fit, and it is the "conserved in three-dimensional space,
   not just in the sequence alignment" part of the plan.

Question 1 can fail on a one-residue misregistration even when the site is really
conserved, so a fourth, register-free check was added: after superposition, find
the candidate residue whose CA is physically closest to each query site CA and
ask whether it is the same amino acid within 4 Å (`frac_nearest_identical`). It
agrees with question 1 almost everywhere (1,589 of 1,594 `site_lost` candidates
are also lost on proximity), so the result below is not an alignment artefact.

Target CA coordinates come from the `tCa` field of the Foldseek server JSON, so
all 12,000 hits could be scored, not only the 300 whose structures were
downloaded in pass 1.

## The method behaves correctly on close homologs

Restricted to the queries that have a curated catalytic residue, among same-fold
Atlas hits:

| Group | Hits | Catalytic residues all kept |
|---|---|---|
| Found by phmmer (E ≤ 1e-3) | 2,243 | 1,882 (83.9%) |
| Structure-only, under 30% identity | 71 | 30 (42.3%) |

Close homologs come out as site-conserved at a high rate, which is the positive
control this check needed. The gap between the two rows is the actual finding.

## Result: most structure-only candidates lose the active site

Of the 1,643 pass-1 candidates:

| Class (alignment register) | n |
|---|---|
| `site_lost` (under 70% of site residues identical) | 1,594 |
| `site_partial` (70% or more identical, not all) | 31 |
| `site_complete` (all identical, RMSD ≤ 2 Å, geometry MAD ≤ 1.5 Å) | 5 |
| `not_checkable` (no site residue aligned) | 13 |

| Class (catalytic residues only) | n |
|---|---|
| `no_catalytic_annotation` (nitroreductase, azoreductase, TP families) | 1,572 |
| `catalytic_partial` | 36 |
| `catalytic_kept` | 30 |
| `catalytic_lost` | 5 |

Read plainly: a fold match at 15% identity, on its own, mostly does **not** come
with the catalytic machinery. The plan anticipated this and asked for these to be
reported separately as probable non-functional relatives, and that is what the
bulk of the 1,643 are. The honest headline from pass 1 and pass 2 together is not
"1,643 candidate drug-inactivating enzymes" but "1,643 structure-only relatives,
of which a few dozen keep the active site".

## The shortlist

`tier1_candidates.tsv` — 13 candidates whose curated catalytic residues are all
present, with site CA RMSD ≤ 2 Å and site geometry MAD ≤ 1.5 Å, at TM ≥ 0.6,
and which phmmer does not find at all (not even at E = 10). Ten of the thirteen
are beta-glucuronidases, the irinotecan family; the other three are cytidine
deaminases, the gemcitabine family.

| Atlas protein | Query | TM | % id | catalytic kept | binding kept | site RMSD | geometry MAD |
|---|---|---|---|---|---|---|---|
| `MGYP000574843407` | P19079 cytidine deaminase | 0.820 | 18.8 | 1/1 | 2/6 | 0.99 | 0.44 |
| `MGYP001319153113` | P19079 cytidine deaminase | 0.744 | 12.7 | 1/1 | 1/6 | 0.98 | 0.37 |
| `MGYP000751766417` | P05804 beta-glucuronidase | 0.722 | 17.6 | 2/2 | 3/5 | 1.57 | 0.65 |
| `MGYP000574654335` | P05804 beta-glucuronidase | 0.706 | 15.4 | 2/2 | 1/5 | 1.37 | 0.96 |
| `MGYP000968177834` | P05804 beta-glucuronidase | 0.697 | 17.4 | 2/2 | 3/5 | 1.25 | 0.42 |
| `MGYP000245278657` | P05804 beta-glucuronidase | 0.692 | 14.5 | 2/2 | 2/5 | 1.98 | 0.94 |
| `MGYP000729808394` | P05804 beta-glucuronidase | 0.684 | 16.0 | 2/2 | 3/5 | 1.28 | 0.43 |
| `MGYP000016927400` | P05804 beta-glucuronidase | 0.674 | 14.8 | 2/2 | 3/5 | 1.43 | 0.57 |
| `MGYP001533818318` | P05804 beta-glucuronidase | 0.663 | 17.0 | 2/2 | 3/5 | 1.44 | 0.76 |
| `MGYP000497510748` | P05804 beta-glucuronidase | 0.657 | 14.5 | 2/2 | 0/5 | 1.15 | 0.87 |
| `MGYP001050749229` | P05804 beta-glucuronidase | 0.644 | 15.6 | 2/2 | 2/5 | 1.76 | 0.50 |
| `MGYP000051719573` | P05804 beta-glucuronidase | 0.635 | 14.4 | 2/2 | 1/5 | 1.74 | 0.46 |
| `MGYP001178215973` | P0ABF6 cytidine deaminase | 0.614 | 18.4 | 1/1 | 4/6 | 0.76 | 0.25 |

The best-supported single candidate is `MGYP000751766417`: both catalytic
glutamates of beta-glucuronidase (E413 the proton donor, E504 the nucleophile)
are present and within 1.7 Å of their query positions, along with N412, Y472 and
K568 of the glucuronate site, at 17.6% sequence identity and with phmmer blind to
it. Its two substitutions are N466→His and W549→His. `MGYP000968177834` and
`MGYP000729808394` have the same pattern. These are the proteins worth docking
and, later, testing.

`tier2_binding_only_candidates.tsv` — 5 azoreductase-fold candidates
(`MGYP002632877957`, `MGYP003573197260`, `MGYP003574735644`, `MGYP001046310337`,
`MGYP001551561688`) that keep the whole annotated FMN site on both checks. They
are listed separately because that site is only three residues, transferred from
a 32%-identity homolog, so the evidence behind them is much thinner than the
table above. Do not present them as equal to the tier-1 list.

## Caveats specific to this pass

- For 1,572 of 1,643 candidates no catalytic residue is annotated anywhere in
  UniProt, so what was tested is cofactor binding. The nitroreductase,
  azoreductase and thymidine phosphorylase families need catalytic residues
  curated from solved structures before their candidates can be judged.
- The beta-glucuronidase window `w1-390` contains exactly one annotated residue
  (D163), so its 55 candidates are effectively unchecked; 13 of them have no
  site residue aligned at all. Its catalytic residues are both in the other
  window.
- Q7A782's site is 3 transferred residues and P07650's is 10, none catalytic.
  Both rest on a homolog, not on the query's own record.
- A conserved site does not mean activity on any drug. Nothing in pass 1 or
  pass 2 touches a drug molecule.
- Candidate structures are ESMFold predictions, and so are the query models. The
  site geometry measured here is geometry between predicted CA positions. Only
  CA atoms are compared, because the Foldseek server returns CA only; side-chain
  orientation is not tested.
- The pass-1 1,000-hit-per-database cap still applies, so these counts are
  floors.

## Files

| File | What it is |
|---|---|
| `query_sites.tsv` | the curated site residues, one row per residue, with source and evidence code |
| `query_sites_provenance.tsv` | per query: direct or transferred, and from where |
| `active_site_all_hits.tsv` | the site check for all 12,000 Atlas alignments |
| `active_site_all_hits_joined.tsv.gz` | same, joined to the pass-1 Foldseek columns and the recomputed phmmer result, with all three classifications |
| `active_site_by_query.tsv` | per query, site conservation in structure-only hits vs hits a sequence search finds |
| `active_site_candidates.tsv` | the 1,643 pass-1 candidates with their site columns, ranked |
| `tier1_candidates.tsv` | the 13-candidate shortlist |
| `tier2_binding_only_candidates.tsv` | the 5 binding-site-only azoreductase candidates |
| `active_site_skipped.tsv` | query models dropped from the check, with reasons (empty in this run) |

Code in `../code/`: `build_sites.py` (UniProt to site table), `active_site_check.py`
(the three checks), `rank_candidates.py` (join, phmmer recomputation, ranking).
The UniProt JSON records used are in `../sites/`.

## What to do next

1. Add the blocked hosts to the environment so steps 5 and 6 can run:
   `www.ebi.ac.uk` for MGnify, the `rcsb.org` subdomains for solved structures,
   and `esmatlas.com` / `search.foldseek.com` to be able to extend the screen at
   all.
2. Curate catalytic residues for the nitroreductase, azoreductase and thymidine
   phosphorylase families from solved structures, then rerun
   `active_site_check.py`. That is the single biggest improvement available.
3. Replace the two 390-residue beta-glucuronidase windows with a full-length
   structure, so the catalytic pair and the glucuronate site are in one model.
4. Dock SN-38 glucuronide into the ten beta-glucuronidase tier-1 candidates.
