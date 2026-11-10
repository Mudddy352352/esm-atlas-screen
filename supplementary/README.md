# Supplementary material

Everything here is generated from files tracked in this repository. Nothing is
hand-entered, so a clean checkout reproduces every table and figure:

    python3 supplementary/make_supplementary.py
    python3 supplementary/make_supplementary_figures.py

`make_supplementary.py` writes `tables/`, `make_supplementary_figures.py` writes
`figures/`. Each table carries two comment lines naming the script and the
source file it was built from.

## What is here

| Item | File | Built from |
|---|---|---|
| Table S1, hits per query model | `tables/Table_S1_hits_per_query.tsv` | `screen_v1/results/summary_per_query.tsv`, `all_hits.tsv.gz` |
| Table S2, habitat of the shortlist | `tables/Table_S2_habitat.tsv` | `screen_v2/results/bodysite_candidates.tsv` |
| Table S3, Pfam recovery | `tables/Table_S3_pfam_recovery.tsv` | `screen_v2/pfam/pfam_recovery.tsv` |
| Table S6, SN-38 glucuronide docking | `tables/Table_S6_docking_sn38g.tsv` | `screen_v2/results/docking/dock_summary_sn38g.tsv` |
| Table S7, docking by ligand | `tables/Table_S7_docking_by_ligand.tsv` | `screen_v2/results/docking/dock_by_receptor_ligand.tsv` |
| Table S12, site classification | `tables/Table_S12_site_classification.tsv` | `screen_v1/results/candidates.tsv`, `screen_v2/results/active_site_all_hits.tsv` |
| Data S2, query site provenance | `tables/Data_S2_query_site_provenance.tsv` | `screen_v2/results/query_sites_provenance.tsv` |
| Data S11, docking boxes | `tables/Data_S11_docking_boxes.tsv` | `screen_v2/results/docking/boxes_final.tsv` |
| Data S12, ligand geometry | `tables/Data_S12_ligand_geometry.tsv` | `screen_v2/results/docking/ligands.tsv` |
| Note S4, per-protein site confidence | `tables/Note_S4_candidate_plddt.tsv` | `screen_v1/results/candidate_plddt.tsv` |
| Note S4, threshold sweep | `screen_v2/evaluation/` | `candidates.tsv`, `active_site_all_hits.tsv` |
| Figure S1, habitat of the shortlist | `figures/Figure_S1_habitat.{png,pdf}` | `screen_v2/results/bodysite_candidates.tsv` |
| Figure S4, pose placement | `figures/Figure_S4_pose_placement.{png,pdf}` | `screen_v2/results/docking/pose_scan_all_poses.tsv` |

Two consistency checks run as part of the build and are printed: Table S12's
site and catalytic class counts, which must match the counts reported in
Results, and the pLDDT summary behind Note S4.

## What is not here, and why

These items are cited in the manuscript but cannot be rebuilt from this
deposit. In each case the missing input is named rather than approximated.

| Item | Missing input |
|---|---|
| Tables S4 and S5, pairwise sequence search | the per-pair SWORD alignment outputs were never written to a checkpoint |
| Table S8, decoy-site null | per-residue Cα coordinates, which live only in the raw Foldseek JSON (~250 MB, not deposited; see `EXCLUDED.md`) |
| Table S9, per-query decoy pass rates | same as S8 |
| Table S10 and Figure 7, M-CSA comparison | fresh Foldseek searches of 62 M-CSA entries against AlphaFold DB |
| Table S11, specificity-residue retention | the per-pair alignment register, available only alongside the raw Foldseek output |
| Note S1, Figures S2 and S3, contact maps | candidate contacts are recoverable, but the crystal-receptor controls need PDB structures that are not vendored here |
| Note S2, site enclosure and insertions | the enclosure measurements were not written to a checkpoint |
| Note S3, functional-atom agreement | needs the full-atom models of the thirteen shortlisted proteins |

The analyses behind these ran after the v1.0.0 snapshot, in an environment that
no longer exists. Section 2.14 of the manuscript states this directly rather
than implying the whole supplement regenerates.
