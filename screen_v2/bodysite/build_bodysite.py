#!/usr/bin/env python3
"""Build the body-site mapping table for the pass-2 tier-1 shortlist.

Inputs
  results/tier1_candidates.tsv          the 13 shortlisted MGYPs
  bodysite/raw/<MGYP>.json              MGnify Proteins API records (current release)
  bodysite/raw_samples/<SRS|ERS>.json   MGnify API sample records
  bodysite/mgy_biomes_2022_05.tsv       rows pulled from the 2022_05 release FTP table
                                        (format: MGYP | count | biome lineage)

Output
  results/bodysite_candidates.tsv       one row per MGYP, collapsed call
  results/bodysite_biome_rows.tsv       one row per MGYP x biome, with observation counts
"""
import csv
import json
import os
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BS = os.path.join(ROOT, "bodysite")
RES = os.path.join(ROOT, "results")


def host_class(lineage):
    """Collapse an MGnify biome lineage into a coarse habitat class."""
    if lineage.startswith("root:Host-associated:Human"):
        return "human"
    if lineage.startswith("root:Host-associated"):
        return "other_host"
    if lineage.startswith("root:Engineered"):
        return "engineered"
    if lineage.startswith("root:Environmental"):
        return "environmental"
    return "unclassified"


def body_site(lineage):
    """Human body site, where the lineage is a human-associated one."""
    if not lineage.startswith("root:Host-associated:Human"):
        return ""
    parts = lineage.split(":")[3:]
    return ":".join(parts) if parts else "unspecified"


def main():
    shortlist = []
    with open(os.path.join(RES, "tier1_candidates.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            shortlist.append(row)

    # 2022_05 biome rows, the release ESM Atlas was built from
    biomes_05 = defaultdict(list)
    path = os.path.join(BS, "mgy_biomes_2022_05.tsv")
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3:
                    biomes_05[f[0]].append((int(f[1]), f[2]))
    else:
        print("warning: no 2022_05 biome table found", file=sys.stderr)

    rows_out, per_biome = [], []
    for row in shortlist:
        mgyp = row["target"]

        api = {}
        p = os.path.join(BS, "raw", mgyp + ".json")
        if os.path.exists(p):
            d = json.load(open(p))
            if "mgyp" in d:
                api = d

        obs = sorted(biomes_05.get(mgyp, []), reverse=True)
        source = "2022_05_ftp"
        if not obs and api:
            obs = sorted(((b["count"], b["name"]) for b in api.get("biomes", [])), reverse=True)
            source = "current_release_api"
        if not obs:
            source = "none"

        total = sum(c for c, _ in obs)
        classes = {host_class(b) for _, b in obs}
        human_obs = sum(c for c, b in obs if host_class(b) == "human")
        sites = sorted({body_site(b) for _, b in obs if body_site(b)})

        for c, b in obs:
            per_biome.append({
                "mgyp": mgyp,
                "query": row["query"],
                "observations": c,
                "biome_lineage": b,
                "habitat_class": host_class(b),
                "human_body_site": body_site(b),
                "biome_source": source,
            })

        assemblies = ";".join(
            s["assembly_accession"] for s in api.get("study_assembly_contigs", [])
        )
        studies = ";".join(
            sorted({s["study_accession"] for s in api.get("study_assembly_contigs", [])})
        )

        rows_out.append({
            "mgyp": mgyp,
            "query": row["query"],
            "tm_score_qnorm": row["tm_score_qnorm"],
            "seq_id_pct": row["seq_id_pct"],
            "site_class": row["site_class"],
            "catalytic_class": row["catalytic_class"],
            "n_biomes": len(obs),
            "total_observations": total,
            "habitat_classes": ";".join(sorted(classes)),
            # "unknown", not "no", where the protein has no biome in any
            # reachable release: absence of metadata is not evidence of habitat
            "human_associated": ("yes" if human_obs else "no") if obs else "unknown",
            "human_observations": human_obs,
            "human_body_sites": ";".join(sites),
            "top_biome_lineage": obs[0][1] if obs else "",
            "all_biome_lineages": ";".join(b for _, b in obs),
            "biome_source": source,
            "in_current_release": "yes" if api else "no",
            "cluster_size": api.get("cluster_size", ""),
            "full_length": api.get("full_length", ""),
            "studies": studies,
            "assemblies": assemblies,
        })

    with open(os.path.join(RES, "bodysite_candidates.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)

    with open(os.path.join(RES, "bodysite_biome_rows.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_biome[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(per_biome)

    n_human = sum(1 for r in rows_out if r["human_associated"] == "yes")
    n_mapped = sum(1 for r in rows_out if r["biome_source"] != "none")
    print(f"{len(rows_out)} candidates, {n_mapped} with a biome, {n_human} human-associated")
    for r in rows_out:
        print(f"  {r['mgyp']}  {r['habitat_classes']:<28} {r['top_biome_lineage']}")


if __name__ == "__main__":
    main()
