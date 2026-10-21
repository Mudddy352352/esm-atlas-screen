#!/usr/bin/env python3
"""Emit the manuscript's tables as Markdown, built from the screen's own files.

    python3 make_tables.py > tables.md

Paste the output into the manuscript, or regenerate after a screen rerun. Same
sources as make_figures.py, so the tables and figures cannot drift apart.
"""
import csv
import gzip
import os

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
V1 = os.path.join(ROOT, "screen_v1", "results")
V2 = os.path.join(ROOT, "screen_v2", "results")
PFAM = os.path.join(ROOT, "screen_v2", "pfam")
DOCK = os.path.join(ROOT, "screen_v2", "results", "docking")


def read_tsv(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def md(rows, headers, aligns=None):
    aligns = aligns or ["---"] * len(headers)
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join(aligns) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


PRETTY = {
    "P05804_uidA_gus_w1-390": "β-glucuronidase, w1-390",
    "P05804_uidA_gus_w214-603": "β-glucuronidase, w214-603",
    "P07650_deoA_tp_w1-390": "thymidine phosphorylase, w1-390",
    "P07650_deoA_tp_w51-440": "thymidine phosphorylase, w51-440",
    "P0ABF6_cdd_ecoli": "cytidine deaminase, E. coli",
    "P0DTQ4_tyrDC_control_w1-390": "tyrosine decarboxylase, w1-390",
    "P0DTQ4_tyrDC_control_w231-620": "tyrosine decarboxylase, w231-620",
    "P15888_nfsB_nitroreductase": "nitroreductase NfsB",
    "P17117_nfsA_nitroreductase": "nitroreductase NfsA",
    "P19079_cdd_tetrameric": "cytidine deaminase, tetrameric",
    "Q7A782_azo1_azoreductase": "azoreductase, flavodoxin-like",
    "Q8KU07_azoB_azoreductase": "azoreductase, Rossmann",
}
SHORT = {
    "P05804_uidA_gus_w214-603": "β-glucuronidase",
    "P05804_uidA_gus_w1-390": "β-glucuronidase",
    "P19079_cdd_tetrameric": "cytidine deaminase (tetrameric)",
    "P0ABF6_cdd_ecoli": "cytidine deaminase (E. coli)",
}


def table2():
    rows = []
    for r in read_tsv(os.path.join(V1, "summary_per_query.tsv")):
        rows.append([
            PRETTY.get(r["query"], r["query"]),
            f'{int(r["atlas_hits_returned"]):,}',
            f'{int(r["same_fold"]):,}',
            f'{int(r["seq_detectable"]):,}',
            f'{int(r["structure_only"]):,}',
            "%.1f" % num(r["min_seq_id"]),
            "%.3f" % num(r["max_tm"]),
        ])
    tot = read_tsv(os.path.join(V1, "summary_per_query.tsv"))
    rows.append([
        "**total**",
        f'{sum(int(r["atlas_hits_returned"]) for r in tot):,}',
        f'{sum(int(r["same_fold"]) for r in tot):,}',
        f'{sum(int(r["seq_detectable"]) for r in tot):,}',
        f'{sum(int(r["structure_only"]) for r in tot):,}',
        "", "",
    ])
    return md(rows, ["Query model", "Atlas hits returned", "Same fold",
                     "Found by phmmer", "Structure-only", "Min. identity (%)",
                     "Max. TM"],
              ["---", "---:", "---:", "---:", "---:", "---:", "---:"])


def table3():
    cands = read_tsv(os.path.join(V2, "active_site_candidates.tsv"))
    n = len(cands)
    site = [("All site residues identical, within geometric thresholds", "site_complete"),
            ("70% or more identical", "site_partial"),
            ("Below 70% identical", "site_lost"),
            ("No site residue aligned", "not_checkable")]
    cat = [("All catalytic residues kept", "catalytic_kept"),
           ("Some catalytic residues kept", "catalytic_partial"),
           ("Catalytic residues lost", "catalytic_lost"),
           ("No catalytic residue annotated for the family", "no_catalytic_annotation")]
    rows = []
    for label, key in site:
        c = sum(1 for r in cands if r["site_class"] == key)
        rows.append(["Whole annotated site", label, f"{c:,}", "%.1f" % (100.0 * c / n)])
    for label, key in cat:
        c = sum(1 for r in cands if r["catalytic_class"] == key)
        rows.append(["Catalytic residues only", label, f"{c:,}", "%.1f" % (100.0 * c / n)])
    return md(rows, ["Scored on", "Class", "n", "% of 1,643"],
              ["---", "---", "---:", "---:"])


def table4():
    t1 = read_tsv(os.path.join(V2, "tier1_candidates.tsv"))
    recov = {}
    for r in read_tsv(os.path.join(PFAM, "pfam_recovery.tsv")):
        recov[r["candidate"]] = recov.get(r["candidate"], False) or \
            r["pfam_recovers_link"] == "yes"
    rows = []
    for r in t1:
        rows.append([
            "`%s`" % r["target"],
            SHORT.get(r["query"], r["query"]),
            "%.3f" % num(r["tm_score_qnorm"]),
            "%.1f" % num(r["seq_id_pct"]),
            "%s/%s" % (r["n_catalytic_identical"], r["n_catalytic"]),
            "%s/%s" % (r["n_binding_identical"], r["n_binding"]),
            "%.2f" % num(r["site_rmsd"]),
            "%.2f" % num(r["site_dist_mad"]),
            r["site_class"].replace("_", " "),
            "yes" if recov.get(r["target"]) else "no",
        ])
    return md(rows, ["Atlas protein", "Query", "TM", "Identity (%)",
                     "Catalytic kept", "Binding kept", "Site RMSD (Å)",
                     "Geometry MAD (Å)", "Overall site class", "Pfam recovers"],
              ["---", "---", "---:", "---:", "---:", "---:", "---:", "---:", "---", "---"])


def table_provenance():
    bs = read_tsv(os.path.join(V2, "bodysite_candidates.tsv"))
    rows = []
    for r in bs:
        lineage = r["top_biome_lineage"].replace("root:", "") or "no metadata"
        rows.append([
            "`%s`" % r["mgyp"],
            r["habitat_classes"] or "unplaced",
            lineage,
            r["studies"] or "—",
            r["assemblies"] or "—",
            r["biome_source"].replace("_", " ") if r["biome_source"] != "none" else "none",
        ])
    return md(rows, ["Atlas protein", "Habitat class", "Biome lineage",
                     "Study", "Assembly", "Metadata source"])


def table_relatives():
    tax = read_tsv(os.path.join(V2, "taxonomy_candidates.tsv"))
    t1 = {r["target"]: r for r in read_tsv(os.path.join(V2, "tier1_candidates.tsv"))}
    rows = []
    for r in tax:
        rows.append([
            "`%s`" % r["mgyp"],
            "%.1f" % num(t1[r["mgyp"]]["seq_id_pct"]),
            "%.1f" % num(r["best_afdb50_seqid"]),
            r["best_afdb50_desc"],
            r["best_afdb50_taxname"],
            r["best_sprot_desc"],
            r["consensus_phylum_distinct"],
        ])
    return md(rows, ["Atlas protein", "Identity to query (%)",
                     "Identity to nearest AFDB50 (%)", "Nearest AFDB50 annotation",
                     "Nearest AFDB50 organism", "Best Swiss-Prot hit",
                     "Phyla in top 20"],
              ["---", "---:", "---:", "---", "---", "---", "---:"])


def table7():
    per = {}
    for r in read_tsv(os.path.join(PFAM, "pfam_recovery.tsv")):
        d = per.setdefault(r["query_model"], [0, 0, r["query_pfam"]])
        d[0] += 1
        d[1] += r["pfam_recovers_link"] == "yes"
    qdom = {}
    for r in read_tsv(os.path.join(PFAM, "pfam_query_domains.tsv")):
        qdom.setdefault(r["sequence"], []).append(r["pfam_name"])
    rows = []
    for q in sorted(per, key=lambda k: -(100.0 * per[k][1] / per[k][0])):
        n, rec, _ = per[q]
        acc = q.split("_")[0]
        fams = sorted({f for s, fs in qdom.items() if s.split("_")[0] == acc for f in fs})
        rows.append([
            PRETTY.get(q, q), f"{n:,}", f"{rec:,}",
            "%.1f" % (100.0 * rec / n),
            ", ".join("`%s`" % f for f in fams) or "none",
        ])
    tot_n = sum(d[0] for d in per.values())
    tot_r = sum(d[1] for d in per.values())
    rows.append(["**total**", f"{tot_n:,}", f"{tot_r:,}",
                 "**%.1f**" % (100.0 * tot_r / tot_n), ""])
    return md(rows, ["Query model", "Candidates", "Recovered by Pfam",
                     "%", "Query's Pfam families"],
              ["---", "---:", "---:", "---:", "---"])


DOCK_NAMES = {
    "3K4D_A": "*E. coli* β-glucuronidase, crystal (3K4D)",
    "6D4O_A": "*E. eligens* β-glucuronidase, crystal (6D4O)",
    "6CXS_A": "*C. perfringens* β-glucuronidase, crystal (6CXS)",
    "P05804_uidA_gus_w214-603": "*E. coli* β-glucuronidase, ESMFold model",
    "193L_A": "hen egg lysozyme, crystal (193L)",
    "1TIM_A": "triose-phosphate isomerase, crystal (1TIM)",
    "P0ABF6_cdd_ecoli": "*E. coli* cytidine deaminase, ESMFold model",
    "Q8KU07_azoB_azoreductase": "NAD(P)H azoreductase, ESMFold model",
    "P0DTQ4_tyrDC_control_w231-620": "L-tyrosine decarboxylase, ESMFold model",
    "3K4D_A_offsite": "β-glucuronidase, box moved off-site (3K4D)",
}
DOCK_KIND = {
    "pos_xtal": "positive control", "pos_pred": "positive control",
    "neg_xtal": "negative control", "neg_pred": "negative control",
    "neg_offsite": "negative control", "candidate": "candidate",
}
DOCK_ORDER = {"pos_xtal": 0, "pos_pred": 0, "candidate": 1,
              "neg_xtal": 2, "neg_pred": 2, "neg_offsite": 2}


def dock_name(receptor):
    return DOCK_NAMES.get(receptor, "`%s`" % receptor)


def table8():
    """Every docked receptor: score, where the ligand sat, and what pass 2 said."""
    pose = {r["receptor"]: r for r in read_tsv(os.path.join(DOCK, "pose_scan_sn38g.tsv"))}
    hab_of = {r["mgyp"]: (r["top_biome_lineage"] or "").replace("root:", "")
              for r in read_tsv(os.path.join(V2, "bodysite_candidates.tsv"))}
    recs = read_tsv(os.path.join(DOCK, "dock_summary_sn38g.tsv"))
    recs.sort(key=lambda r: (DOCK_ORDER[r["kind"]],
                             -num(pose[r["receptor"]]["frac_poses_catalytic"])
                             if r["receptor"] in pose else 1e9))
    rows = []
    for r in recs:
        p = pose.get(r["receptor"])
        if r["kind"] == "candidate":
            hab = hab_of.get(r["receptor"], "") or "no metadata"
        else:
            hab = "—"
        if p is None:
            frac = closest = pen = "not mapped"
        else:
            frac = "%.0f" % (100.0 * num(p["frac_poses_catalytic"]))
            closest = "%.2f" % num(p["min_sugar_dist_any_pose"])
            pen = p["penalty_vs_best"]
            pen = ("%.2f" % num(pen)) if pen not in ("", "nan") else "—"
        rows.append([
            dock_name(r["receptor"]), DOCK_KIND[r["kind"]],
            "%.2f" % num(r["best_mean"]), closest, frac, pen,
            (r["site_class"] or "—").replace("_", " "), hab,
        ])
    return md(rows, ["Receptor", "Class", "Vina affinity (kcal/mol)",
                     "Closest sugar approach (Å)", "Catalytic poses (%)",
                     "Score given up (kcal/mol)", "Pass-2 site class", "Habitat"],
              ["---", "---", "---:", "---:", "---:", "---:", "---", "---"])


LIG_CLASS = [
    ("Positive controls, crystal", lambda r: r["kind"] == "pos_xtal"),
    ("Positive control, ESMFold model", lambda r: r["kind"] == "pos_pred"),
    ("Shortlisted candidates", lambda r: r["kind"] == "candidate"),
    ("Unrelated proteins, ESMFold models", lambda r: r["kind"] == "neg_pred"),
    ("Unrelated proteins, crystal", lambda r: r["kind"] == "neg_xtal"),
    ("β-glucuronidase, box moved off-site", lambda r: r["kind"] == "neg_offsite"),
]


def table9():
    """The three ligands side by side: the score tracks size, not recognition."""
    recs = read_tsv(os.path.join(DOCK, "dock_summary_sn38g.tsv"))
    rows = []
    for label, test, in LIG_CLASS:
        grp = [r for r in recs if test(r)]
        def mean(key):
            vals = [num(r[key]) for r in grp if r[key] not in ("", "nan")]
            return sum(vals) / len(vals)
        g, a = mean("best_mean"), mean("sn38_aglycone")
        rows.append([label, "%d" % len(grp), "%.2f" % g, "%.2f" % a,
                     "%.2f" % mean("glucuronic_acid"), "%+.2f" % (g - a)])
    return md(rows, ["Receptor class", "n", "SN-38 glucuronide", "SN-38 (aglycone)",
                     "β-D-glucuronic acid", "Gain from the glucuronide"],
              ["---", "---:", "---:", "---:", "---:", "---:"])


def table10():
    """The receptor set the pose criterion was built on, so it can be audited."""
    boxes = read_tsv(os.path.join(DOCK, "boxes_final.tsv"))
    kindname = {"candidate": "candidate", "pos_xtal": "positive control, crystal",
                "pos_pred": "positive control, model",
                "neg_xtal": "unrelated protein, crystal",
                "neg_pred": "unrelated protein, model",
                "neg_offsite": "off-site box"}
    boxes.sort(key=lambda r: (DOCK_ORDER[r["kind"]], r["receptor"]))
    rows = []
    for r in boxes:
        src = r["site_source"]
        if src.startswith("centroid of bound"):
            how = "centroid of the bound ligand"
        elif src.startswith("surface patch"):
            how = "surface patch 55 Å from the site"
        elif src.startswith("annotated site"):
            how = "annotated site residues (%s)" % src.split("residues ")[1].replace(";", ", ")
        else:
            res = src.split("residues ")[1].replace(";", ", ")
            how = ("catalytic residues mapped through the alignment (%s)" % res
                   if r["kind"] == "candidate"
                   else "known catalytic residues (%s)" % res)
        note = r["note"].split(";")[0].replace("beta-", "β-")
        rows.append([dock_name(r["receptor"]), kindname[r["kind"]], note,
                     "%s Å" % r["box"], how])
    return md(rows, ["Receptor", "Class", "What it is", "Box", "How the box was centred"])


if __name__ == "__main__":
    for name, fn in [("Table 2", table2), ("Table 3", table3), ("Table 4", table4),
                     ("Table 5", table_relatives), ("Table 6", table_provenance),
                     ("Table 7", table7), ("Table 8", table8),
                     ("Table 9", table9), ("Table 10", table10)]:
        print("\n\n## %s\n" % name)
        print(fn())
