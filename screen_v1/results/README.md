# Pass 1: Foldseek screen of the ESM Metagenomic Atlas

Run on 20 October 2026 in the project cloud environment. Everything below is computed from files in this directory; nothing is quoted from memory or from the literature.

## What was run

1. Query enzymes were taken from the curated Swiss-Prot bacterial protein set bundled with Prokka (24,727 reviewed entries), selected by their EC number and enzyme name, one or two per family in section 4 of the study plan.
2. Each query sequence was folded with the public ESMFold API (`api.esmatlas.com/foldSequence/v1/pdb`). The API refuses sequences longer than 400 residues, so the three long queries were folded as two 390-residue windows each.
3. Each model was searched with the Foldseek web server (`search.foldseek.com`, mode 3Di+AA) against three databases: `mgnify_esm30` (the ESM Metagenomic Atlas, MGnify proteins clustered at 30%), `afdb50` (AlphaFold DB UniProt50) and `pdb100` (the PDB).
4. TM-scores were computed locally from the Foldseek alignments by Kabsch superposition of the aligned CA pairs with iterative refinement, normalised by query length.
5. For each query, phmmer (HMMER3 via pyhmmer 0.12.3) was run from the query sequence against the sequences of its own Atlas structural hits, to measure how many of them a sequence search would also have found.

## Query set and identity check

The PDB search is a control on the query itself: if the top PDB hit is the expected enzyme at high identity, the sequence really is that enzyme and the ESMFold model is accurate.

| Query | Family | Drug relevance | Length | Top PDB hit | % id |
|---|---|---|---|---|---|
| `P05804_uidA_gus_w1-390` | Beta-glucuronidase | irinotecan (SN-38 reactivation) | 390 | Crystal structure of E. coli beta-glucuronidase in complex | 98.4 |
| `P05804_uidA_gus_w214-603` | Beta-glucuronidase | irinotecan (SN-38 reactivation) | 390 | Crystal structure of E. coli beta-glucuronidase in complex | 99.4 |
| `P07650_deoA_tp_w1-390` | Thymidine phosphorylase | 5-fluorouracil pathway | 390 | STRUCTURAL AND THEORETICAL STUDIES SUGGEST DOMAIN MOVEMENT | 99.7 |
| `P07650_deoA_tp_w51-440` | Thymidine phosphorylase | 5-fluorouracil pathway | 390 | X-ray structure of the thymidine phosphorylase from Salmon | 94.6 |
| `P0ABF6_cdd_ecoli` | Cytidine deaminase | gemcitabine, cytarabine | 294 | CRYSTAL STRUCTURE OF CYTIDINE DEAMINASE COMPLEXED WITH URI | 100.0 |
| `P0DTQ4_tyrDC_control_w1-390` | L-tyrosine decarboxylase (control) | levodopa | 390 | Crystal structure of a tyrosine decarboxylase from Enteroc | 99.7 |
| `P0DTQ4_tyrDC_control_w231-620` | L-tyrosine decarboxylase (control) | levodopa | 390 | Crystal structure of a tyrosine decarboxylase from Enteroc | 97.9 |
| `P15888_nfsB_nitroreductase` | NAD(P)H nitroreductase NfsB | nitro-prodrugs | 217 | Structure of the Salmonella typhimurium nfnB dihydropterid | 100.0 |
| `P17117_nfsA_nitroreductase` | NADPH nitroreductase NfsA | nitro-prodrugs | 240 | E. coli NfsA bound to NADP+ | 98.7 |
| `P19079_cdd_tetrameric` | Cytidine deaminase (tetrameric type) | gemcitabine, cytarabine | 136 | Bacillus subtilis cytidine deaminase with an Arg56 - Gln s | 99.2 |
| `Q7A782_azo1_azoreductase` | FMN-dependent NADPH-azoreductase | azo-prodrugs | 188 | Crystal Structure of ArsH from Sinorhizobium meliloti | 19.8 |
| `Q8KU07_azoB_azoreductase` | NAD(P)H azoreductase | azo-prodrugs | 286 | polyketide ketoreductase SimC7 - ternary complex with NADP | 38.2 |

## Atlas hits per query

`same_fold` = TM-score (query-normalised) >= 0.5 and query coverage >= 0.5. `structure_only` = same fold but phmmer from the query does not find it at E <= 1e-3. The Foldseek server returns at most 1000 hits per database, so these counts are a floor, not the true size of each family in the Atlas.

| query | atlas_hits_returned | same_fold | same_fold_div | seq_detectable | structure_only | structure_only_div | min_seq_id | max_tm |
|---|---|---|---|---|---|---|---|---|
| P05804_uidA_gus_w1-390 | 1000 | 628 | 607 | 573 | 55 | 55 | 3.5 | 0.928 |
| P05804_uidA_gus_w214-603 | 1000 | 803 | 720 | 775 | 28 | 28 | 8.1 | 0.954 |
| P07650_deoA_tp_w1-390 | 1000 | 927 | 533 | 542 | 385 | 385 | 9.2 | 0.989 |
| P07650_deoA_tp_w51-440 | 1000 | 848 | 171 | 848 | 0 | 0 | 16.9 | 0.991 |
| P0ABF6_cdd_ecoli | 1000 | 109 | 23 | 107 | 2 | 2 | 12.6 | 0.993 |
| P0DTQ4_tyrDC_control_w1-390 | 1000 | 93 | 62 | 91 | 2 | 2 | 6.1 | 0.917 |
| P0DTQ4_tyrDC_control_w231-620 | 1000 | 408 | 383 | 372 | 36 | 36 | 4.4 | 0.964 |
| P15888_nfsB_nitroreductase | 1000 | 995 | 881 | 699 | 296 | 296 | 7.6 | 0.98 |
| P17117_nfsA_nitroreductase | 1000 | 972 | 883 | 654 | 318 | 318 | 10.8 | 0.975 |
| P19079_cdd_tetrameric | 1000 | 994 | 59 | 989 | 5 | 5 | 12.7 | 0.94 |
| Q7A782_azo1_azoreductase | 1000 | 983 | 947 | 722 | 261 | 261 | 8.0 | 0.95 |
| Q8KU07_azoB_azoreductase | 1000 | 877 | 761 | 622 | 255 | 255 | 10.0 | 0.97 |
| **total** | | **8637** | **6030** | **6994** | **1643** | **1643** | | |

## Candidates

1643 Atlas proteins share a query fold, sit below 30% sequence identity to the query, and are not found by phmmer from the query. These are the proteins a sequence-based screen would have missed.

Structure confidence was fetched from the Atlas for 300 of them: median pLDDT 95.2, 300 above 70.

| Atlas protein | Query | TM-score | % id | q cov | aln | phmmer E |
|---|---|---|---|---|---|---|
| `MGYP000010903125` | Q8KU07 | 0.849 | 16.0 | 0.96 | 275 | 0.006 |
| `MGYP001241316458` | Q7A782 | 0.829 | 19.7 | 0.98 | 184 | not found |
| `MGYP001627315251` | Q7A782 | 0.827 | 16.3 | 0.98 | 184 | 0.005 |
| `MGYP000574843407` | P19079 | 0.820 | 18.8 | 0.90 | 123 | not found |
| `MGYP001269804397` | Q7A782 | 0.819 | 20.5 | 0.94 | 176 | not found |
| `MGYP000913334691` | Q7A782 | 0.817 | 16.9 | 0.93 | 174 | 0.002 |
| `MGYP000187952076` | Q8KU07 | 0.813 | 18.4 | 0.97 | 278 | not found |
| `MGYP001070460782` | P15888 | 0.809 | 17.2 | 0.93 | 202 | 0.004 |
| `MGYP001087563302` | Q8KU07 | 0.808 | 18.3 | 0.96 | 276 | not found |
| `MGYP001422630899` | Q7A782 | 0.807 | 15.6 | 0.95 | 178 | not found |
| `MGYP000196079655` | Q8KU07 | 0.805 | 17.8 | 0.96 | 276 | not found |
| `MGYP002640187951` | Q8KU07 | 0.804 | 14.9 | 0.98 | 281 | not found |
| `MGYP001461789676` | Q7A782 | 0.800 | 18.1 | 0.90 | 169 | 0.002 |
| `MGYP001165233942` | P19079 | 0.797 | 18.0 | 0.88 | 120 | not found |
| `MGYP000265226133` | Q8KU07 | 0.794 | 19.4 | 0.95 | 272 | 0.002 |
| `MGYP000063503949` | P19079 | 0.793 | 19.4 | 0.88 | 120 | 0.002 |
| `MGYP000695878256` | Q7A782 | 0.791 | 19.2 | 0.95 | 179 | 0.004 |
| `MGYP000670750219` | Q7A782 | 0.790 | 17.3 | 0.90 | 170 | not found |
| `MGYP001464171282` | Q7A782 | 0.790 | 16.5 | 0.91 | 171 | not found |
| `MGYP000520589075` | Q8KU07 | 0.789 | 20.5 | 0.95 | 271 | not found |
| `MGYP000042736691` | Q8KU07 | 0.787 | 17.1 | 0.97 | 278 | not found |
| `MGYP000234380966` | Q8KU07 | 0.787 | 17.9 | 0.93 | 267 | not found |
| `MGYP001342696119` | Q7A782 | 0.787 | 21.0 | 0.93 | 174 | not found |
| `MGYP000947209926` | Q7A782 | 0.784 | 17.8 | 0.92 | 173 | not found |
| `MGYP002724366964` | P15888 | 0.784 | 22.2 | 0.92 | 200 | not found |

### Tiers

- 1640 distinct Atlas proteins (3 hit more than one query).
- 200 at TM-score >= 0.7 and under 20% identity, which is the tier worth looking at first.
- Candidate identity to the query runs 7.6% to 24.4% (median 15.6%).
- phmmer finds nothing at all for 1423 of them (no hit even at E = 10); the other 220 come up only as weak hits between E = 1e-3 and E = 10. The first group is the honest "sequence search would have missed this" set, and 145 of those are at TM-score >= 0.7 and under 20% identity.

## Caveats specific to this pass

- The two azoreductase queries have no near-identical structure in the PDB (top hits at 19.8% and 38.2% identity), so their family assignment rests on the Swiss-Prot annotation, not on a solved structure of that protein. The two are also different structural classes (flavodoxin-like FMN reductase vs Rossmann-fold reductase), which is why their hit sets barely overlap. A structurally characterised azoreductase should replace them in the next pass.
- The three long queries were searched as 390-residue windows, which cuts domains in half. The two thymidine phosphorylase windows behave very differently (385 structure-only hits vs none), so window choice is doing real work here and the counts for those three families should not be read as family sizes. Full-length structures are needed.
- Organism assignments for the queries come from the top PDB hit, not from UniProt, because UniProt is blocked here.

## What this does and does not show

- These are predicted structures matched to predicted structures. Nothing here shows that any candidate binds or inactivates any drug.
- The active-site check (step 3 of the plan) has not been done. It needs curated catalytic and cofactor residues per family, which come from UniProt feature annotations and the primary literature; `rest.uniprot.org`, `www.rcsb.org` and `www.ebi.ac.uk` are blocked by this environment's network policy, so that curation could not be done here.
- Source-sample and taxonomy mapping (steps 5 and 6) also need MGnify or IMG, which are on the same blocked hosts. The Atlas database in the Foldseek server carries no taxonomy field, so Atlas hits here have no organism or body site.
- Docking (step 4) was not attempted.
- The 1000-hit cap per database means the screen is sensitivity-limited. A full screen needs Foldseek run locally against the downloaded Atlas, which is the cluster job in the plan.

## Files

- `results/README.md` (0.0 MB)
- `results/all_hits.tsv` (4.9 MB)
- `results/all_hits.tsv.gz` (6.7 MB)
- `results/candidate_plddt.tsv` (0.0 MB)
- `results/candidates.tsv` (0.2 MB)
- `results/candidates.tsv.gz` (0.3 MB)
- `results/summary_per_query.tsv` (0.0 MB)
- `queries/` query FASTA files, Swiss-Prot accession in the name
- `structures/` ESMFold models used as Foldseek queries
- `foldseek_raw/` raw Foldseek server JSON, one file per query
- `code/` the pipeline: esmfold.py, run_foldseek.py, parse_results.py, filter_and_compare.py, fetch_plddt.py, make_report.py
