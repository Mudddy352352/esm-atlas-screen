#!/usr/bin/env python3
"""Build the step-6 taxonomy table for the tier-1 shortlist.

What this does and does not claim
---------------------------------
MGnify assigns no taxonomy to MGYP accessions: the 2022_05 protein-database
release ships no taxonomy file, and Foldseek reports taxonomy=False for the
mgnify_esm30 database. The organism each candidate came from is therefore not
recoverable from the Atlas or MGnify metadata.

What is computed here is the taxonomy of each candidate's closest
*structurally characterised* relatives, from Foldseek searches against AFDB50,
AlphaFold/Swiss-Prot and PDB100, all of which do carry taxonomy. That describes
what the protein resembles. It is not the source organism, and a consensus over
many distant relatives is weaker still.

Inputs : raw_foldseek/<MGYP>.json[.gz]          Foldseek results (run_taxonomy_search.py)
         ../results/tier1_candidates.tsv   the shortlist
         ../results/bodysite_candidates.tsv  habitat, to put beside the taxonomy
Outputs: ../results/taxonomy_candidates.tsv   one row per candidate
         ../results/taxonomy_top_hits.tsv     per candidate x database, top hits
         lineage_cache.json                   taxId -> NCBI lineage, from UniProt
"""
import csv
import gzip
import json
import os
import sys
import time
import urllib.error
import urllib.request
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw_foldseek")
RES = os.path.join(os.path.dirname(HERE), "results")
CACHE = os.path.join(HERE, "lineage_cache.json")

TOP_N = 10          # hits written out per candidate per database
CONSENSUS_N = 20    # AFDB50 hits used for the consensus call
# UniProt's taxonomy lineage labels the top rank "domain" (older NCBI dumps call
# it "superkingdom"); accept either and normalise to "domain".
RANKS = ("domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus")


def load_cache():
    if os.path.exists(CACHE):
        return json.load(open(CACHE))
    return {}


def fetch_lineage(tax_id, cache):
    """NCBI lineage for a taxId, via the UniProt taxonomy endpoint."""
    key = str(tax_id)
    if key in cache:
        return cache[key]
    url = f"https://rest.uniprot.org/taxonomy/{key}.json"
    out = {"name": None, "ranks": {}}
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            d = json.load(r)
        out["name"] = d.get("scientificName")
        for anc in d.get("lineage", []):
            rank = (anc.get("rank") or "").lower()
            if rank in RANKS:
                out["ranks"][rank] = anc.get("scientificName")
        rank = (d.get("rank") or "").lower()
        if rank in RANKS:
            out["ranks"].setdefault(rank, d.get("scientificName"))
        if "domain" not in out["ranks"] and "superkingdom" in out["ranks"]:
            out["ranks"]["domain"] = out["ranks"]["superkingdom"]
    except urllib.error.HTTPError as exc:
        out["error"] = f"HTTP {exc.code}"
    except Exception as exc:
        out["error"] = str(exc)
    cache[key] = out
    time.sleep(0.2)
    return out


def hits_for(doc, db):
    for b in doc.get("results", []):
        if b.get("db") == db:
            al = b.get("alignments") or []
            if al and isinstance(al[0], list):
                return al[0]
            return al
    return []


def by_evalue(hits):
    return sorted(hits, key=lambda h: (float(h.get("eval", 1e9)), -float(h.get("seqId") or 0)))


def clean_target(t):
    """Trim the description Foldseek appends to the accession."""
    return (t or "").split()[0]


def describe(t):
    parts = (t or "").split(None, 1)
    return parts[1].strip() if len(parts) > 1 else ""


def main():
    shortlist = list(csv.DictReader(open(os.path.join(RES, "tier1_candidates.tsv")), delimiter="\t"))
    body = {}
    bpath = os.path.join(RES, "bodysite_candidates.tsv")
    if os.path.exists(bpath):
        body = {r["mgyp"]: r for r in csv.DictReader(open(bpath), delimiter="\t")}

    cache = load_cache()
    rows, hit_rows = [], []

    for rec in shortlist:
        mgyp = rec["target"]
        # results are stored gzipped (they carry target Ca coordinates and are bulky)
        path = os.path.join(RAW, mgyp + ".json")
        if os.path.exists(path + ".gz"):
            with gzip.open(path + ".gz", "rt") as fh:
                doc = json.load(fh)
        elif os.path.exists(path):
            doc = json.load(open(path))
        else:
            print(f"warning: no Foldseek result for {mgyp}", file=sys.stderr)
            continue

        row = {
            "mgyp": mgyp,
            "query": rec["query"],
            "site_class": rec["site_class"],
            "tm_score_qnorm": rec["tm_score_qnorm"],
            "habitat_class": (body.get(mgyp) or {}).get("habitat_classes", ""),
            "human_body_sites": (body.get(mgyp) or {}).get("human_body_sites", ""),
        }

        for db, tag in (("afdb50", "afdb50"), ("afdb-swissprot", "sprot"), ("pdb100", "pdb")):
            hits = by_evalue(hits_for(doc, db))
            row[f"n_hits_{tag}"] = len(hits)
            if not hits:
                continue
            best = hits[0]
            tax_id = best.get("taxId")
            lin = fetch_lineage(tax_id, cache) if tax_id else {"ranks": {}}
            row[f"best_{tag}_target"] = clean_target(best.get("target"))
            row[f"best_{tag}_desc"] = describe(best.get("target"))
            row[f"best_{tag}_eval"] = best.get("eval")
            row[f"best_{tag}_seqid"] = best.get("seqId")
            row[f"best_{tag}_taxname"] = best.get("taxName")
            row[f"best_{tag}_domain"] = lin.get("ranks", {}).get("domain", "")
            row[f"best_{tag}_phylum"] = lin.get("ranks", {}).get("phylum", "")

            for h in hits[:TOP_N]:
                hit_rows.append({
                    "mgyp": mgyp,
                    "database": db,
                    "target": clean_target(h.get("target")),
                    "description": describe(h.get("target")),
                    "evalue": h.get("eval"),
                    "seq_id_pct": h.get("seqId"),
                    "prob": h.get("prob"),
                    "aln_length": h.get("alnLength"),
                    "tax_id": h.get("taxId"),
                    "tax_name": h.get("taxName"),
                })

        # consensus across the closest AFDB50 relatives
        top = by_evalue(hits_for(doc, "afdb50"))[:CONSENSUS_N]
        doms, phyla = Counter(), Counter()
        for h in top:
            if not h.get("taxId"):
                continue
            lin = fetch_lineage(h["taxId"], cache).get("ranks", {})
            if lin.get("domain"):
                doms[lin["domain"]] += 1
            if lin.get("phylum"):
                phyla[lin["phylum"]] += 1
        row["consensus_n"] = len(top)
        row["consensus_domain"] = "; ".join(f"{k}:{v}" for k, v in doms.most_common(3))
        row["consensus_phylum"] = "; ".join(f"{k}:{v}" for k, v in phyla.most_common(5))
        row["consensus_phylum_distinct"] = len(phyla)
        rows.append(row)
        json.dump(cache, open(CACHE, "w"))
        print(f"  {mgyp} best afdb50 = {row.get('best_afdb50_taxname')}", flush=True)

    json.dump(cache, open(CACHE, "w"))

    fields = list(rows[0].keys())
    for r in rows:
        for k in fields:
            r.setdefault(k, "")
    with open(os.path.join(RES, "taxonomy_candidates.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(RES, "taxonomy_top_hits.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(hit_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(hit_rows)

    print(f"\n{len(rows)} candidates, {len(hit_rows)} hit rows")


if __name__ == "__main__":
    main()
