# A structural screen of the ESM Metagenomic Atlas for microbial drug-metabolising enzymes

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
manuscript/
  figures/                 make_figures.py, make_tables.py, and the rendered
                           figures (300-dpi PNG + vector PDF)
LICENSE-CODE               MIT (code)
LICENSE-DATA               CC-BY-4.0 (data) — confirm before release, see below
CITATION.cff               citation metadata (fill in author, DOI)
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
| Pfam comparator | `python3 screen_v2/pfam/pfam_scan.py --pfam <Pfam-A.hmm> --queries screen_v2/pfam/queries.faa --candidates screen_v2/pfam/candidates.faa --pairs screen_v1/results/candidates.tsv --out <dir>` | Pfam-A.hmm (see EXCLUDED.md) |
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
- **Docking prep and run** (`prep.py`, `run_dock.py`) need AutoDock Vina, Open
  Babel, RDKit and the crystal/ligand structures. The 72 pose files and all
  result tables are shipped, so the analysis (`analyze.py`, `pose_scan.py`) can
  be inspected and the geometry audited; a from-scratch rerun requires the tools
  above and regenerating the `work/receptors/` intermediates with `prep.py`.

Verification done when this package was assembled: `build_sites.py`,
`rank_candidates.py`, `make_tables.py` and `make_figures.py` were rerun from the
files here. `build_sites.py` and `rank_candidates.py` reproduced their shipped
outputs byte for byte (including the 13-protein shortlist); all six figures and
all tables regenerated cleanly.

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
