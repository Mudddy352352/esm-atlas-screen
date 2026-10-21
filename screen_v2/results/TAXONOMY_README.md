# Step 6: taxonomy of the tier-1 shortlist

Run on 21 October 2026, after `www.ebi.ac.uk`, `api.esmatlas.com` and
`search.foldseek.com` became reachable again. This is step 6 of the study plan.
Nothing earlier was re-run: the Foldseek screen, the active-site check and the
body-site mapping are unchanged.

Inputs: the 13 proteins in `tier1_candidates.tsv`.
Outputs: `taxonomy_candidates.tsv` (one row per protein) and
`taxonomy_top_hits.tsv` (top 10 hits per protein per database).
Structures, raw Foldseek results, the NCBI lineage cache and both scripts are in
`../taxonomy/`.

## The honest answer to "what organisms are these?" is: not knowable

MGnify assigns **no taxonomy to MGYP accessions at all**. This was checked at
source rather than assumed:

- The 2022_05 protein-database release ships no taxonomy file. Its files cover
  sequences, assemblies, biomes, contigs, clusters, counts and Pfam, and the
  README documents no taxonomic assignment.
- Foldseek's database listing reports `taxonomy: false` for `mgnify_esm30`, the
  ESM Atlas database, against `taxonomy: true` for AFDB50, AlphaFold/Swiss-Prot,
  AlphaFold/Proteome, PDB100, CATH50, BFMD and BFVD.
- A metagenomic assembly is a mixed community, and the contig a protein was
  called from carries no organism label in this resource.

So step 6 cannot report a source organism. What it reports instead is the
taxonomy of each protein's **closest structurally characterised relatives**, from
Foldseek searches of the real ESM Atlas structures (not ProstT5-predicted 3Di)
against AFDB50, AlphaFold/Swiss-Prot and PDB100. That says what each protein
resembles. It is not where it came from, and it must not be written up as though
it were.

## The finding that matters most for the manuscript

**These proteins are not invisible to sequence search. Their link to the
drug-metabolising query enzymes is.**

| | Identity |
|---|---|
| To its query enzyme (the screen's own number) | min 12.7%, **median 15.6%**, max 18.8% |
| To its closest labelled AFDB50 protein | min 19.0%, **median 47.7%**, max 61.7% |

**12 of the 13 have a labelled AlphaFold DB relative at 30% identity or better**,
up to 61.7%. Those relatives are annotated as ordinary carbohydrate-active
enzymes: "Glycoside hydrolase family 2", "Beta-glycosidase",
"Beta-galactosidase", "Beta-mannosidase B", "CMP/dCMP-type deaminase".

A plain sequence search does find these proteins. What it returns is a generic
GH2 glycosidase or a generic deaminase, not beta-glucuronidase or cytidine
deaminase specifically. The screen's contribution is the *link* to a
drug-metabolising enzyme, which structure search makes and sequence search does
not. "Invisible to sequence search" overstates it, and the number 1,643 cannot be
described as 1,643 proteins unknown to sequence databases.

## Body site and taxonomy point at different candidates

No protein on the shortlist is both sampled from a human gut **and**
beta-glucuronidase-like by its labelled relatives.

**The most beta-glucuronidase-like is MGYP000051719573**, and it is not a gut
sample. Its closest AFDB50 relative is a *Phocaeicola vulgatus* beta-glycosidase
at 47.7% identity, its best Swiss-Prot hit is *Bacteroides fragilis*
beta-glucuronidase, and its top PDB hits are characterised gut bacterial
beta-glucuronidases: *Bacteroides uniformis* GUS (5uj6, 6d50) and *Tannerella
forsythia* GUS (8dhl). Taxonomically its relatives are gut Bacteroidota. But it
was sampled from terephthalate-degrading wastewater in North Carolina, and pass 2
scored its active site `site_lost`.

**The two human-gut proteins look like generic GH2 enzymes.** For both
MGYP001533818318 and MGYP000245278657 the closest AFDB50 relatives are
*Anaerolineaceae* and *Verrucomicrobiales* "Glycoside hydrolase family 2" at
41–51% identity, and the best Swiss-Prot hit is fungal beta-mannosidase B. Their
closest solved structures are gut Bacteroidota GH2 enzymes (*Phocaeicola
plebeius* 8u01, *Bacteroides thetaiotaomicron* 2vl4/2vot/2vqu) but only at
14–18% identity.

Read together with the body-site step, the shortlist offers a gut-sampled pair
whose relatives are unremarkable glycosidases, and a GUS-like protein from
wastewater that has lost its site. That is a weaker position than either step
suggests alone, and the write-up should say so rather than pick whichever
candidate suits the sentence.

## Phylum-level attribution is not supportable

Taking the 20 closest AFDB50 relatives per protein and counting phyla, most
proteins spread across 9 to 16 distinct phyla. No phylum call is defensible for
them. Two exceptions, both narrow:

| Protein | Distinct phyla in top 20 | Dominant |
|---|---|---|
| MGYP000574654335 | 5 | Bacteroidota 14/20 |
| MGYP001178215973 | 5 | Streptophyta 12/20 |

At domain level the best hit is Bacteria for 10 of 13, Archaea for 2
(MGYP000574843407, MGYP001319153113, both cytidine-deaminase queries) and
Eukaryota for 1.

**MGYP001178215973 is the odd one.** Its nearest relatives are land plants
(Streptophyta 12 of 20, best hit a cytidine deaminase from *Escallonia rubra* at
19.0% identity), and 19.0% is the weakest best-hit identity on the shortlist. For
a protein from an unplaced metagenomic contig, a plant-dominated relative set is
more likely to reflect what is densely represented in AlphaFold DB than anything
about the protein. Treat it as unresolved.

## Caveats

- Everything here is nearest-relative inference. Identity to the best relative
  runs 19–62%, so the confidence varies a lot across the 13 and a single table
  row flattens that. Use `taxonomy_top_hits.tsv`, not just the best hit.
- Database composition drives the consensus. AlphaFold DB and the PDB are not a
  random sample of life, so a phylum count over 20 hits measures what has been
  sequenced and solved as much as what the protein is related to.
- Annotation names are transferred labels, not assays. "Beta-glucuronidase" on a
  Swiss-Prot entry at 14.4% identity to our protein is a statement about that
  entry, not about our protein.
- The pass-2 caveat is unchanged and still governs the headline number: 1,643
  structure-only relatives, a few dozen of which retain the active site, never
  1,643 candidate drug-inactivating enzymes. For 1,572 of the 1,643 UniProt
  annotates no catalytic residue at all, so for most of the list the active-site
  check tested cofactor binding rather than catalysis, and an absent site is
  partly absent annotation.
- Nothing in this step is evidence that any candidate acts on any drug. That
  still needs assays.

## How it was run

Structures for all 13 came from `screen_v1/atlas_structures/` where present (7)
and from `api.esmatlas.com/fetchPredictedStructure/` otherwise (6). Sequences
extracted from the PDB files match `t_len` in `tier1_candidates.tsv` for all 13,
and match the MGnify API sequence exactly for the 6 proteins where both exist, so
the structures correspond to the proteins that were screened.

Each structure was submitted to `search.foldseek.com` in `3diaa` mode against
AFDB50, AlphaFold/Swiss-Prot and PDB100 (`../taxonomy/run_taxonomy_search.py`),
giving 1000, ~300 and ~170 hits respectively per protein. Hits were re-sorted by
e-value, since the server returns them by score. NCBI lineages for the hit taxIds
came from `rest.uniprot.org/taxonomy/{id}.json` and are cached in
`../taxonomy/lineage_cache.json` (`../taxonomy/build_taxonomy.py`). UniProt labels
the top rank `domain`, not `superkingdom`.
