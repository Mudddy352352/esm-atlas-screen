# β-Glucuronidase-fold proteins in the ESM Metagenomic Atlas: fold conservation without active-site conservation

Code and data accompanying the manuscript. This repository holds everything
needed to check the numbers in the paper and to rerun each analysis step, plus
cached copies of the external records the steps depend on so that most of them
reproduce without network access.

> **Please read this first.** The honest claim of this work is narrow, and the
> files here are arranged to support that narrow claim and no wider one. The
> study does **not** report the discovery of new drug-inactivating enzymes. Its
> contribution is **prioritisation within protein families that annotation
> already identifies**: ranking metagenomic proteins by 3D similarity to a
> specific drug-metabolising structure, then filtering on active-site geometry.
> A docking step found **no** supported human-gut candidate for irinotecan
> reactivation. See "What the study shows" below.

## What the study did

Twelve ESMFold models of nine characterised bacterial enzymes (acting on
gemcitabine, irinotecan, 5-fluorouracil, nitro- and azo-prodrugs, and levodopa)
were searched against the ESM Metagenomic Atlas with Foldseek. The pipeline then:

1. **Pass 1 (`screen_v1/`)** — Foldseek structure search, giving 8,637 same-fold
   alignments and 1,643 hits below 30% identity to the query that a
   single-sequence search (phmmer) does not find.
2. **Pass 2, active-site check (`screen_v2/`)** — for each hit, whether the
   query's curated catalytic and binding residues are present and in the right
   place in 3D. Most hits have lost the site; 13 keep every curated catalytic
   residue (the "tier-1" shortlist).
3. **Pfam comparator (`screen_v2/pfam/`)** — a Pfam-A profile search recovers
   92.5% of the 1,643 and 11 of the 13, so these are not proteins sequence
   methods miss.
4. **Body-site mapping (`screen_v2/results/` + `screen_v2/bodysite/`)** — the
   habitat of each shortlisted protein from MGnify. 2 are human-gut, 5
   environmental/industrial, 6 unplaceable.
5. **Taxonomy (`screen_v2/taxonomy/`)** — nearest characterised relatives. MGnify
   assigns no taxonomy to these accessions, so the source organism is not
   recoverable; only nearest-relative inference exists.
6. **Docking (`screen_v2/results/docking/`)** — SN-38 glucuronide docked into the
   shortlist against positive and negative controls.
7. **Pairwise and nr searches (`screen_v3/nr/`)** — each shortlisted protein
   searched against NCBI nr. Twelve of the thirteen have a relative already in
   nr at a median 74% identity, so these proteins were in the databases
   undescribed rather than missing from them.
8. **Evaluation (`screen_v3/benchmark/`, `screen_v3/mcsa/`)** — what the screen
   is worth: a decoy-site null model, a labelled retrieval benchmark built from
   the AlphaFold DB and PDB hits the same searches returned, baselines at
   matched recall, and a repeat of the retention comparison over 62 curated
   M-CSA families and 49,735 same-fold hits.
9. **Revision analyses (`screen_v3/revision/`)** — the pocket-enclosure,
   site-loop, side-chain, specificity-motif, pose-ensemble and
   threshold-sensitivity measurements the manuscript's later sections report.

## What the study shows (and does not)

- A shared fold at ~15% identity mostly does **not** bring the active site with
  it: among queries with a curated catalytic residue, 83.9% of
  sequence-detectable same-fold hits keep all catalytic residues, against 42.3%
  of structure-only hits. This difference is the central result.
- The candidates are **not invisible to sequence search**. A Pfam profile search
  places 92.5% of them, and the nearest labelled AlphaFold-DB relative sits at a
  median 47.7% identity. Phrases like "dark matter" or "invisible to sequence
  search" are not supported by the data and are not used.
- The Pfam family a hit is assigned to is **broad** (e.g. Glyco_hydro_2 covers
  glucuronidases, galactosidases and mannosidases alike), so a hit places a
  protein in a family, not on a specific drug substrate.
- **Docking does not rank the candidates.** Positive and negative controls
  overlap on Vina score; what separates them is pose geometry (whether the sugar
  sits against the catalytic carboxylates). On that measure, **neither of the two
  human-gut candidates places the sugar catalytically**, so the screen supports
  no gut irinotecan candidate.
- **Nothing here has been tested against a drug.** Every structure is predicted,
  and docking is a plausibility check, not evidence of activity.

The per-step READMEs state the caveats for each step in full. They are the
authoritative record; this file is a map.

## Layout

```
screen_v1/                 Pass 1: Foldseek structure search
  code/                    esmfold.py, run_foldseek.py, parse_results.py,
                           filter_and_compare.py, fetch_plddt.py, make_report.py
  queries/                 the 9 query enzyme FASTA files (Swiss-Prot)
  structures/              the 12 ESMFold query models used as Foldseek queries
  atlas_structures/        ESMFold models of the 13 shortlisted Atlas proteins
  results/                 README.md + hit tables (all_hits.tsv.gz,
                           candidates.tsv[.gz], summary_per_query.tsv, pLDDT)
screen_v2/                 Pass 2 and the follow-on steps
  code/                    build_sites.py, active_site_check.py, rank_candidates.py
  sites/                   cached UniProt JSON records for the 9 queries
  bodysite/                build_bodysite.py + cached MGnify records
  taxonomy/                build_taxonomy.py, run_taxonomy_search.py,
                           candidate structures, cached Foldseek + NCBI records
  pfam/                    pfam_scan.py + query/candidate FASTAs and outputs
  results/                 README.md, the site-check tables, both shortlists,
                           body-site, taxonomy and docking outputs
    docking/               step-4 README, code/, all 72 SN-38G pose files,
                           ligands, and the docking result tables
screen_v3/                 nr search and the evaluation of the screen
  nr/                      NR_README.md, code/, the classified nr hit tables
                           (raw BLAST XML not shipped, see EXCLUDED.md)
  benchmark/               BENCHMARK_README.md, code/, labels/, the decoy null,
                           the labelled-hit scores and the baseline comparison
  mcsa/                    MCSA_README.md, code/, the 62 M-CSA query families
                           and the per-family retention tables
  revision/                REVISION_README.md, code/, structures/ and the
                           pocket, loop, side-chain, contact and sweep tables
manuscript/
  figures/                 make_figures.py, make_tables.py, and the rendered
                           figures (300-dpi PNG + vector PDF)
LICENSE-CODE               MIT (code)
LICENSE-DATA               CC-BY-4.0 (data) — confirm before release, see below
CITATION.cff               citation metadata (fill in the full author list)
EXCLUDED.md                large files not shipped, and how to regenerate them
```

## Reproducing each step

All scripts resolve paths relative to the repository root, so a clone works
wherever it sits. To run against a copy elsewhere, set `ESMATLAS_ROOT` to its
path. Python 3.11 was used.

```bash
pip install numpy pandas biopython pyhmmer rdkit matplotlib
# AutoDock Vina 1.2.5 and Open Babel 3.1.1 are needed only for docking.
```

| Step | Command (from the repository root) | Needs the network? |
|---|---|---|
| Query active sites | `python3 screen_v2/code/build_sites.py screen_v2/sites <out>` | no — reads cached UniProt JSON |
| Rank candidates | `python3 screen_v2/code/rank_candidates.py screen_v1 screen_v2/results` | no |
| Habitat mapping | `python3 screen_v2/bodysite/build_bodysite.py` | no — reads cached MGnify records |
| Pfam comparator | `python3 screen_v2/pfam/pfam_scan.py --pfam <Pfam-A.hmm> --queries screen_v2/pfam/queries.faa --candidates screen_v2/pfam/candidates.faa --pairs screen_v1/results/candidates.tsv --out <dir>` | Pfam-A.hmm (see EXCLUDED.md) |
| M-CSA retention | `python3 screen_v3/mcsa/code/analyse_retention.py screen_v3/mcsa` | no |
| Baselines at matched recall | `python3 screen_v3/benchmark/code/baseline_comparison.py screen_v3/benchmark` | no |
| Revision statistics | `python3 screen_v3/revision/code/stats_for_revision.py` | no |
| nr hit tables | `python3 screen_v3/nr/code/analyse_nr.py screen_v3/nr` | raw BLAST XML (see EXCLUDED.md) |
| Figures | `python3 manuscript/figures/make_figures.py --outdir <dir>` | no |
| Tables | `python3 manuscript/figures/make_tables.py > tables.md` | no |

Steps that **cannot** be fully rerun inside this package, by design, because
they need a blocked/heavy external resource or a large intermediate that is not
shipped (see EXCLUDED.md for each):

- **Pass 1 Foldseek search** (`run_foldseek.py`) needs the Foldseek web server
  and the ESMFold API. Its parsed output (`screen_v1/results/all_hits.tsv.gz`)
  is shipped, so everything downstream reproduces from it.
- **The active-site geometry step** (`active_site_check.py`) needs the raw
  Foldseek JSON (`screen_v1/foldseek_raw/`, ~250 MB, not shipped) for the
  per-residue coordinates. Its output (`active_site_all_hits.tsv`) is shipped as
  the checkpoint, and `rank_candidates.py` runs from there.
- **The nearest-relative search** (`build_taxonomy.py`) needs its own raw
  Foldseek output (`screen_v2/taxonomy/raw_foldseek/`, ~172 MB, not shipped).
  Its outputs (`taxonomy_candidates.tsv`, `taxonomy_top_hits.tsv`) are shipped.
  Run without that directory the script warns per protein and then stops with an
  `IndexError`, writing nothing, so the shipped tables survive — see
  `EXCLUDED.md` §4.
- **The nr searches** (`submit_all.py`, `collect_nr.py`) need NCBI's BLAST URL
  API, and the parsing steps need the raw XML it returns (~116 MB, not shipped).
  The classified hit tables in `screen_v3/nr/results/` are shipped, and every
  number in the paper's nr section comes from them. A fresh search will not
  match them, because nr grows between runs.
- **The decoy null, the labelled benchmark and the M-CSA scan**
  (`decoy_site_null.py`, `score_labelled_hits.py`, `mcsa_retention_scan.py`)
  read the same raw Foldseek JSON as the active-site geometry step, and the
  M-CSA scan also fetches curated families from the M-CSA and AlphaFold DB. All
  of their outputs are shipped, and the three analyses that produce the paper's
  numbers from those outputs (`analyse_retention.py`,
  `baseline_comparison.py`, `stats_for_revision.py`) run offline from this
  package.
- **Docking prep and run** (`prep.py`, `run_dock.py`) need AutoDock Vina, Open
  Babel, RDKit and the crystal/ligand structures. The 72 pose files and all
  result tables are shipped, so the analysis (`analyze.py`, `pose_scan.py`) can
  be inspected and the geometry audited; a from-scratch rerun requires the tools
  above and regenerating the `work/receptors/` intermediates with `prep.py`.

Verification done when this package was assembled, all of it from the files
here and with no network access:

- `build_sites.py` reproduced `query_sites.tsv` and
  `query_sites_provenance.tsv` **byte for byte**.
- `rank_candidates.py` reproduced `tier1_candidates.tsv` (the 13-protein
  shortlist), `active_site_candidates.tsv`, `active_site_by_query.tsv` and the
  joined hit table **byte for byte**.
- `build_bodysite.py` reproduced `bodysite_candidates.tsv` and
  `bodysite_biome_rows.tsv` **byte for byte**.
- `make_tables.py` regenerated Tables 2–10 and `make_figures.py` all six
  figures. (Figures are not byte-compared: the rendering shifts slightly with
  the matplotlib version.)
- When `screen_v3/` was added on 10 November 2026,
  `screen_v3/revision/code/stats_for_revision.py` was rerun from this package
  with `ESMATLAS_ROOT` pointing at it and reproduced every table in
  `screen_v3/revision/results/` **byte for byte**, and
  `analyse_retention.py` and `baseline_comparison.py` reprinted the retention
  and baseline numbers the manuscript quotes (37.4% against 10.2%,
  p = 1.3 × 10⁻⁷, and 0.156 against the TM-score baseline's 0.294).

## Licensing

- **Code** (everything under `*/code/`, the `*.py` scripts) is offered under the
  MIT License (`LICENSE-CODE`).
- **Data** (result tables, structures, cached records) is offered under
  CC-BY-4.0 (`LICENSE-DATA`) as a sensible default. **Confirm this before
  release.** Cached third-party records carry their own terms: UniProt
  (CC-BY-4.0), MGnify / EMBL-EBI and RCSB PDB data, and ChEBI structures are
  redistributed here for reproducibility under those providers' terms, which
  take precedence over this repository's data licence for those files.

## Citation

See `CITATION.cff`. The author name and the Zenodo DOI are placeholders until
you fill them in.
