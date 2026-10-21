# Step 5: body-site / habitat mapping for the tier-1 shortlist

Run on 21 October 2026, after the project's cloud environment was widened so that
`www.ebi.ac.uk` and `ftp.ebi.ac.uk` became reachable. This is step 5 of the study
plan, the one pass 2 could not do. Nothing else was re-run: the Foldseek screen
and the active-site check are unchanged.

Inputs: the 13 proteins in `tier1_candidates.tsv`.
Outputs: `bodysite_candidates.tsv` (one row per protein) and
`bodysite_biome_rows.tsv` (one row per protein × biome, with observation counts).
Cached source records are in `../bodysite/raw/`, `../bodysite/raw_samples/` and
`../bodysite/mgy_biomes_2022_05.tsv`; the build script is
`../bodysite/build_bodysite.py`.

## Headline

**Of the 13 shortlisted proteins, 2 come from a human body site. Both are human
faecal / large intestine. 5 come from environmental or industrial habitats, and 6
cannot be placed at all.**

| Habitat | Proteins |
|---|---|
| Human, large intestine (faecal) | 2 |
| Environmental (marine sediment, marine brackish) | 2 |
| Engineered (activated sludge, terephthalate wastewater, fermented beverages) | 3 |
| No metadata in any reachable release | 6 |

### The two human-gut proteins

Both are beta-glucuronidase (P05804) structural relatives, which is the
irinotecan-relevant family, and both keep their two catalytic glutamates:

| MGYP | TM (q-norm) | Seq id | Site class | Catalytic | Body site |
|---|---|---|---|---|---|
| MGYP001533818318 | 0.6633 | 17.0% | `site_partial` (5/7 identical) | 2/2 kept | Human:Digestive system:Large intestine:Fecal |
| MGYP000245278657 | 0.6918 | 14.5% | `site_lost` / `site_partial_3d` (5/7 by nearest-residue) | 2/2 kept | Human:Digestive system:Large intestine:Fecal |

MGYP000245278657 is traceable to a single sample: ERS396295 (`MH0222`), a Danish
human stool metagenome, assembly ERZ509137, study ERP108254. MGYP001533818318 has
a biome assignment but no assembly row in the 2022_05 release, so it stops at the
biome lineage.

### The best-scoring candidate is not a gut protein

MGYP000751766417, which pass 2 ranked first (TM 0.7223, 17.6% identity, both
catalytic glutamates present), comes from **marine sediment** in the Gulf of
California (sample SRS3033588, assembly ERZ1742239, study ERP109156). It is a
genuine beta-glucuronidase-fold relative with an intact pair of catalytic
glutamates, but it is not a gut organism's protein, so it cannot carry a claim
about drug inactivation in the human gut. If the manuscript wants a gut-resident
example, the two proteins above are the ones to use, and both score slightly
lower than MGYP000751766417.

The remaining placed proteins are from marine brackish water (MGYP000574654335),
activated sludge (MGYP001319153113), terephthalate-degrading wastewater
(MGYP000051719573) and a Chinese fermented-beverage pit (MGYP000968177834).

## How the mapping was done

MGnify assigns a biome to each **assembly**, and the protein database propagates
that to every MGYP found in it, so "body site" here is the biome of the assembly
the sequence was called from, not an independent measurement on the protein.

Two sources were used, and `biome_source` in the output records which one each row
came from:

1. **MGnify 2022_05 `mgy_biomes.tsv.gz`**, the release ESM Atlas was built on, and
   therefore the matching release for these accessions. The file is 5.5 GB
   gzipped and holds one row per MGYP × biome, sorted by accession. It was
   streamed from `ftp.ebi.ac.uk` and filtered in flight, stopping once the stream
   passed the highest accession on the shortlist: 941,039,296 rows read, nothing
   written to disk. 5 of the 13 have a row.
2. **The MGnify Proteins API** (`/metagenomics/proteins/api/v1/protein/{mgyp}`),
   which serves the current release. 6 of the 13 resolve there; the other 7 return
   404. For the 2 proteins it covers that 2022_05 does not
   (MGYP000751766417, MGYP001319153113) its biome was used, and `biome_source`
   says `current_release_api`.

Where both sources have a protein, they agree exactly, biome lineage and
observation count alike, in all 4 overlapping cases.

Sample-level metadata (geography, host, environmental package) was then pulled
from the main MGnify API for the 6 assemblies the Proteins API exposes, and is
cached in `../bodysite/raw_samples/`. Only ERS396295 carries human host
annotation.

## Caveats, and what this does not say

**Six of the thirteen are unplaced, which is not the same as not being human.**
The 2022_05 README states that a subset of sequences in that release have no
associated metadata and therefore no biome and no ERZ assembly. The six are in
that subset: MGYP000574843407, MGYP000729808394, MGYP000016927400,
MGYP000497510748, MGYP001050749229 and MGYP001178215973. This was checked
directly rather than assumed, and MGYP000016927400 is absent from the 2022_05
assembly table as well as the biome table. They are also absent from the current
release. So the honest reading is 2 human, 5 non-human, 6 unknown, and the
denominator for any proportion is 7, not 13.

**Later releases were not usable.** The 2023_02, 2024_04 and 2026_07 biome tables
carry far more rows per protein (a prefix that reaches accession ~23,000,000 in
2022_05 reaches only ~400,000 in 2023_02), so streaming them to the accessions on
this shortlist is not practical. 2026_07 also publishes a 24 GB Parquet version,
which would allow selective reads if this ever needs to be revisited.

**The counts are tiny.** Every placed protein has an observation count of 1, i.e.
each was seen in a single assembly. These are singletons, not abundant or
widespread gut enzymes, and cluster sizes are 1 or 2. Nothing here supports a
claim about prevalence across people or samples.

**The pass-2 caveat still stands and is not weakened by this step.** The correct
description of the pass-1 output is 1,643 structure-only relatives, a few dozen of
which retain the active site, never 1,643 candidate drug-inactivating enzymes. For
1,572 of the 1,643, UniProt annotates no catalytic residue at all, so for most of
the list the active-site check tested cofactor binding rather than catalysis, and
an absent site is partly absent annotation. Body-site mapping narrows the
shortlist further; it adds no evidence that any of these proteins is active on any
drug. That still needs assays.
