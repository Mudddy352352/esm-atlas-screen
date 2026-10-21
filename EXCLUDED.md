# Files deliberately not included, and how to regenerate them

This package is kept small (~30 MB) so it can live in git. A few large or
third-party inputs are left out. In every case either the derived output they
feed is shipped, or the file is freely downloadable from a public source. None
of the paper's numbers depends on a file that is both excluded and
unregenerable.

## 1. Pass-1 raw Foldseek server JSON — `screen_v1/foldseek_raw/` (~250 MB)

The raw per-query JSON returned by the Foldseek web server, one file per query
model. It holds the per-residue Cα coordinates (`tCa`) that
`screen_v2/code/active_site_check.py` reads to measure site geometry over all
12,000 alignments.

- **Why excluded:** ~250 MB, a bulky intermediate.
- **Downstream output that is shipped:** the parsed hit table
  `screen_v1/results/all_hits.tsv.gz` and the active-site checkpoint
  `screen_v2/results/active_site_all_hits.tsv`. `rank_candidates.py` and all
  figures/tables reproduce from these.
- **To regenerate:** rerun `screen_v1/code/run_foldseek.py` against the Foldseek
  web server (`search.foldseek.com`, mode 3Di+AA, databases `mgnify_esm30`,
  `afdb50`, `pdb100`) with the query models in `screen_v1/structures/`. The
  server returns at most 1,000 hits per database, so results are a floor and may
  shift slightly as the server's databases are updated.

## 2. Pfam-A HMM library — `Pfam-A.hmm` (~1.5 GB uncompressed)

The full Pfam-A profile library that `screen_v2/pfam/pfam_scan.py` scans the
query and candidate sequences against at gathering thresholds.

- **Why excluded:** large, and freely available from EBI; versioned there, not
  by this study.
- **Downstream output that is shipped:** the scan results
  (`pfam_query_domains.tsv`, `pfam_candidate_domains.tsv`, `pfam_recovery.tsv`,
  `pfam_summary.txt`) are all in `screen_v2/pfam/`.
- **To regenerate:** download `Pfam-A.hmm.gz` from the EBI FTP server
  (`ftp.ebi.ac.uk`, Pfam current release), gunzip, `hmmpress` it, then run
  `pfam_scan.py --pfam Pfam-A.hmm ...` as shown in the top-level README. Record
  the Pfam release used; recovery percentages can move with the release.

## 3. Full MGnify biome table — `mgy_biomes.tsv.gz` (5.5 GB gzipped)

The MGnify 2022_05 release's protein-to-biome table, streamed and filtered in
flight by `screen_v2/bodysite/build_bodysite.py` (nothing was written to disk;
941 M rows were read and the stream stopped past the shortlist's highest
accession).

- **Why excluded:** 5.5 GB; it is a public release file, not study data.
- **What is shipped instead:** only the rows matching the 13 shortlisted
  proteins, in `screen_v2/bodysite/mgy_biomes_2022_05.tsv`, plus the cached
  MGnify API records in `screen_v2/bodysite/raw/` and `raw_samples/`, and the
  final `bodysite_candidates.tsv` / `bodysite_biome_rows.tsv` in
  `screen_v2/results/`.
- **To regenerate:** stream `mgy_biomes.tsv.gz` for the MGnify **2022_05**
  release from `ftp.ebi.ac.uk` (the release the ESM Atlas was built on; later
  releases renumber accessions and are not a drop-in substitute — see
  `BODYSITE_README.md`).
- **You do not need it to rerun the habitat step.** `build_bodysite.py` runs
  offline from the shipped subset and the cached API records, and reproduces
  `bodysite_candidates.tsv` and `bodysite_biome_rows.tsv` byte for byte. This
  was verified from a clean checkout of this repository.

## 4. Taxonomy raw Foldseek JSON — `screen_v2/taxonomy/raw_foldseek/` (~172 MB)

Raw Foldseek results of the 13 candidate structures against AFDB50,
AlphaFold/Swiss-Prot and PDB100, used by `build_taxonomy.py`.

- **Why excluded:** ~172 MB of raw server output.
- **What is shipped instead:** the candidate structures
  (`screen_v2/taxonomy/structures/`), the NCBI lineage cache
  (`lineage_cache.json`), the run log, and the final
  `taxonomy_candidates.tsv` / `taxonomy_top_hits.tsv` in `screen_v2/results/`.
- **To regenerate:** rerun `run_taxonomy_search.py` against the Foldseek web
  server with the shipped candidate structures.

> **What happens if you run `build_taxonomy.py` without this directory:** it
> prints `warning: no Foldseek result for <MGYP>` for each of the 13 proteins,
> then stops with `IndexError: list index out of range` (it tries to read the
> column names off an empty result list). It writes nothing, so **the shipped
> tables are left intact** — the traceback is ugly but nothing is lost.
> Regenerate `raw_foldseek/` first if you mean to rerun the step.

## 5. Docking working intermediates — `work/`, `receptors/`

Prepped PDBQT receptors and scratch files produced by
`screen_v2/results/docking/code/prep.py`.

- **Why excluded:** regenerable, and they need Open Babel to rebuild anyway.
- **What is shipped instead:** all 72 SN-38G pose files (`docking/poses/`), the
  prepared ligands (`docking/ligands/`), and every docking result table.
- **To regenerate:** run `prep.py` (needs Open Babel 3.1.1 and RDKit) over a
  workdir seeded with the candidate structures and the crystal/control
  structures, then `run_dock.py` (needs AutoDock Vina 1.2.5). The crystal
  structures (3K4D, 6D4O, 6CXS, 193L, 1TIM) come from the RCSB PDB and the
  ligands from ChEBI; see `docking/README.md` for the exact accessions.

## Note on duplicated structures

The 13 shortlisted ESMFold models appear in more than one place
(`screen_v1/atlas_structures/`, `screen_v2/results/docking/atlas_structures_extra/`
and `screen_v2/taxonomy/structures/`). This is intentional: each step stays
self-contained, at a cost of a few MB. The copies are identical.
