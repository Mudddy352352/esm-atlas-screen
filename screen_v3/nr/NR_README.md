# The real nr search

Run 8 November 2026. Section 3.8 of the manuscript reports a BLAST-class aligner run on
local sequence sets with its E-values scaled to the size of a large database, because no
sequence search service was reachable when the screen ran. NCBI's BLAST URL API answers
from this environment now, so the search itself has been run: all thirteen shortlisted
Atlas proteins against nr, with nr's own statistics.

## What was run

`blastp` against nr through NCBI's URL API, one job per protein, hit lists taken to 5,000
so that a weak but real alignment is not buried under thousands of same-family hits at
E = 0. Request ids are written to disk as each job is submitted and results are collected
separately, because a long-running process does not survive this environment pausing; the
searches live on NCBI's side and are fetched whenever the session is awake.

## What nr says about the thirteen

**These proteins are already in the databases, unannotated.** Twelve of the thirteen have a
close relative in nr, median 74.1% identity to the best hit and two above 99%. In twelve of
the thirteen the single best hit carries a name that states no activity: "glycoside
hydrolase family 2", "hypothetical protein", "DUF8487 domain-containing protein". So the
screen is not reaching proteins the sequence databases lack. It is reaching proteins the
sequence databases hold and have not annotated, which is what the manuscript already
claims and can now state as a measurement.

**No shortlisted protein is significantly linked to the characterised query enzyme.** The
only hit to the query enzyme's own organism anywhere in the thirteen deep searches is
*Escherichia coli* cytidine deaminase in MGYP001178215973's results, at E = 1.5, which is
not significant. Section 3.8's central claim survives the real search.

**But twelve of the thirteen do align significantly to something carrying the activity
name**, at median 30.8% identity (table below). Those names are almost all automatic
annotations on uncultured or environmental genomes, with one Swiss-Prot entry among them
and that one labelled "Putative". This is the part the projection understated, and it
matters: a reviewer who repeats the search sees "beta-glucuronidase" in the hit list, so
the manuscript has to say what those hits are rather than leave the impression that nr
returns nothing of the kind.

| Protein | Activity searched | Best E to an unambiguously named enzyme | Identity | The hit |
| --- | --- | ---: | ---: | --- |
| MGYP000016927400 | β-glucuronidase | 5.6e-12 | 29.3% | Swiss-Prot, putative β-glucuronidase |
| MGYP000051719573 | β-glucuronidase | 6.5e-23 | 26.5% | *Novipirellula herctigrandis* |
| MGYP000245278657 | β-glucuronidase | 5.0e-12 | 25.1% | *Maribellus comscasis* |
| MGYP000497510748 | β-glucuronidase | 8.0e-29 | 24.7% | *Phycisphaerales* bacterium, putative |
| MGYP000574654335 | β-glucuronidase | 0.0 | 45.0% | *Anaerohalosphaera lusitana* |
| MGYP000574843407 | cytidine deaminase | 1.1e-97 | 53.1% | *Candidatus* Woesearchaeota archaeon |
| MGYP000729808394 | β-glucuronidase | 1.4e-12 | 32.2% | *Aestuariibaculum lutulentum* |
| MGYP000751766417 | β-glucuronidase | 1.6e-14 | 26.3% | *Verrucomicrobiales* bacterium |
| MGYP000968177834 | β-glucuronidase | 5.0e-14 | 33.0% | *Labilibaculum* sp. |
| MGYP001050749229 | β-glucuronidase | none | - | only "β-galactosidase/β-glucuronidase", at 3.8e-141 |
| MGYP001178215973 | cytidine deaminase | 6.2e-08 | 34.1% | *Meiothermus* sp. |
| MGYP001319153113 | cytidine deaminase | 3.7e-28 | 36.2% | *Patescibacteria* group bacterium |
| MGYP001533818318 | β-glucuronidase | 3.0e-14 | 25.3% | MAG, *Opitutia* bacterium |

A hit description in nr joins the titles of every identical sequence with ">", so each
title is classified on its own. A dual "β-galactosidase/β-glucuronidase" name is counted
separately from an unambiguous one, because it names the ambiguity rather than resolving
it, and it is the only kind of glucuronidase label MGYP001050749229 reaches.

## Where this disagrees with the projection

Section 3.8 projected that only MGYP000751766417 would be linked to its query at a
realistic database size, at E = 5.6e-5, and that MGYP001319153113 would produce no
alignment to any query at all. Against the real nr:

- MGYP000751766417 has no significant alignment to the *E. coli* enzyme, so the projection
  was optimistic for the one protein it picked out.
- MGYP001319153113 aligns to a named cytidine deaminase at E = 3.7e-28, so the projection
  was pessimistic for the protein it said would align to nothing. That protein carries no
  Pfam domain either, which is what made its case interesting; the real search shows the
  sequence is not unreachable, only unannotated in the profile libraries.
- The qualitative conclusion the section draws, that a pairwise search hands over a family
  or a loose activity label and not the specific drug-metabolising enzyme, is unchanged.

## An analysis that was run and discarded

The direct way to ask the activity question looked like restricting the search to proteins
named as the activity with an Entrez query. Every one of the thirteen came back with no
hits, which read as a clean negative result. It is an artefact: *E. coli* β-glucuronidase
searched against nr restricted to β-glucuronidase also returns no hits, and so do the
`[Title]` and free-text spellings. The restriction does not work on this endpoint. The
searches were rerun unrestricted with deep hit lists instead, which is what the table above
reports. `code/nr_activity_probe.py` is kept only so the discarded attempt is on the
record.

## Files

| File | What it holds |
| --- | --- |
| `code/run_nr_blast.py` | Submission and collection, one protein at a time |
| `code/submit_all.py` | Submits everything outstanding and writes the request ids down |
| `code/collect_nr.py` | Fetches whatever has finished; safe to rerun |
| `code/nr_deep.py` | The same searches with 5,000-hit lists |
| `code/classify_nr_hits.py` | Title-by-title classification of what the hits are named |
| `code/analyse_nr.py` | The best-hit summary, which is where the 74.1% comes from |
| `results/MGYP*.xml`, `deep_MGYP*.xml` | The searches themselves |
| `results/nr_shortlist_summary.tsv`, `nr_activity_classified.tsv` | The two tables above |
