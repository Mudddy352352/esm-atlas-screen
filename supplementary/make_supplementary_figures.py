#!/usr/bin/env python3
"""Supplementary Figures S1 and S4, from the deposited tables.

S1  Habitat of the 13 shortlisted proteins, with the unplaced ones shown rather
    than dropped, since the usable denominator is 7 and not 13 (section 3.7).
S4  Where SN-38 glucuronide sits, pose by pose: the distance from the glucuronide's
    sugar to the receptor's catalytic carboxylates for every pose of every seed,
    with the two proteins sampled from human faeces drawn against the positive and
    negative controls (sections 3.9 and 2.11).

Both run offline. Usage:  python3 supplementary/make_supplementary_figures.py
"""
import argparse, collections, csv, io, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))

INK, INK2 = "#1b1b1b", "#666666"
CUT = 4.0          # the catalytic-placement threshold used throughout
GUT = ("MGYP001533818318", "MGYP000245278657")   # the two human-faecal proteins
CLASS_COLOUR = {"pos_xtal": "#1f6f4a", "pos_pred": "#5aa17f",
                "candidate": "#2b5d8a", "neg_pred": "#b06a2c",
                "neg_xtal": "#9c3b3b", "neg_offsite": "#6d4c7d"}
CLASS_LABEL = {"pos_xtal": "positive control (crystal)",
               "pos_pred": "positive control (predicted)",
               "candidate": "shortlisted candidate",
               "neg_pred": "unrelated (predicted)",
               "neg_xtal": "unrelated (crystal)",
               "neg_offsite": "off-site box"}


def read(rel):
    with io.open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
        return [r for r in csv.DictReader(
            (ln for ln in fh if not ln.startswith("#")), delimiter="\t")]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def tidy(ax, grid_axis="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#cccccc")
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis=grid_axis, color="#eeeeee", lw=0.8)
    ax.set_axisbelow(True)


def figure_s1(out):
    rows = read("screen_v2/results/bodysite_candidates.tsv")
    counts = collections.Counter()
    for r in rows:
        h = (r["habitat_classes"] or "").strip()
        counts[h if h else "no biome assignment"] += 1
    order = sorted(counts, key=lambda k: (k == "no biome assignment", -counts[k], k))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.9))

    ax1.barh(range(len(order)), [counts[k] for k in order],
             color=["#bdbdbd" if k == "no biome assignment" else "#2b5d8a" for k in order],
             height=0.62)
    ax1.set_yticks(range(len(order)))
    ax1.set_yticklabels([k.replace(";", ", ") for k in order], fontsize=8)
    ax1.invert_yaxis()
    ax1.set_xlabel("shortlisted proteins", fontsize=8, color=INK2)
    ax1.set_xticks(range(0, max(counts.values()) + 2))
    for i, k in enumerate(order):
        ax1.text(counts[k] + 0.08, i, str(counts[k]), va="center", fontsize=8, color=INK2)
    ax1.set_title("a   Habitat of the 13 shortlisted proteins",
                  loc="left", fontsize=9, color=INK)
    tidy(ax1, grid_axis="x")

    placed = sum(c for k, c in counts.items() if k != "no biome assignment")
    human = sum(1 for r in rows if (r["human_associated"] or "").strip().lower()
                in ("yes", "true", "1"))
    segs = [("human-associated", human, "#2b5d8a"),
            ("other habitat", placed - human, "#8fb3cf"),
            ("unplaced", counts.get("no biome assignment", 0), "#bdbdbd")]
    left = 0
    for label, n, colour in segs:
        if not n:
            continue
        ax2.barh([0], [n], left=left, color=colour, height=0.45, label="%s (%d)" % (label, n))
        left += n
    ax2.set_yticks([])
    ax2.set_xlim(0, len(rows))
    ax2.set_xlabel("shortlisted proteins", fontsize=8, color=INK2)
    ax2.set_title("b   The denominator is 7, not 13", loc="left", fontsize=9, color=INK)
    ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=3,
               fontsize=7.5, frameon=False, labelcolor=INK2,
               handlelength=1.4, columnspacing=1.2, borderpad=0.0)
    tidy(ax2, grid_axis="x")

    fig.tight_layout(rect=(0, 0.085, 1, 1))
    fig.text(0.53, 0.035,
             "every placed protein was seen in exactly one assembly,\n"
             "so nothing here supports a claim about prevalence",
             fontsize=7, color=INK2, ha="left", va="bottom")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, "Figure_S1_habitat.%s" % ext), dpi=300)
    plt.close(fig)
    print("  Figure_S1_habitat          habitats %s" % dict(counts))


def figure_s4(out):
    rows = read("screen_v2/results/docking/pose_scan_all_poses.tsv")
    by = collections.defaultdict(list)
    kind = {}
    for r in rows:
        d = num(r["sugar_to_catalytic"])
        if d is None:
            continue
        by[r["receptor"]].append(d)
        kind[r["receptor"]] = r["kind"]

    rank = {"pos_xtal": 0, "pos_pred": 1, "candidate": 2,
            "neg_pred": 3, "neg_xtal": 4, "neg_offsite": 5}
    order = sorted(by, key=lambda rc: (rank.get(kind[rc], 9), min(by[rc])))

    fig, ax = plt.subplots(figsize=(9.6, 7.0))
    for i, rc in enumerate(order):
        vals = by[rc]
        colour = CLASS_COLOUR.get(kind[rc], INK2)
        is_gut = rc in GUT
        ax.scatter(vals, [i] * len(vals), s=14,
                   facecolors="none" if is_gut else colour,
                   edgecolors=colour, linewidths=1.0, alpha=0.85, zorder=3)
        ax.plot([min(vals), max(vals)], [i, i], color=colour, lw=0.8, alpha=0.35, zorder=2)
        frac = sum(1 for v in vals if v <= CUT) / len(vals)
        ax.text(1.005, i, "%3.0f%% of %d poses" % (100 * frac, len(vals)),
                transform=ax.get_yaxis_transform(), clip_on=False,
                va="center", ha="left", fontsize=7,
                color=INK if is_gut else INK2,
                fontweight="bold" if is_gut else "normal")
    ax.axvline(CUT, color="#9c3b3b", ls="--", lw=1.0, zorder=1)
    ax.text(CUT + 0.12, -0.9, "4 Å catalytic-placement criterion",
            fontsize=7.5, color="#9c3b3b")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([("%s  (human faeces)" % rc) if rc in GUT else rc
                        for rc in order], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(0, max(16.0, max(max(v) for v in by.values()) + 0.5))
    ax.set_xlabel("distance from the glucuronide's sugar to the catalytic carboxylates (Å)",
                  fontsize=8.5, color=INK2)
    ax.set_title("Supplementary Figure S4   Where the ligand sits, pose by pose",
                 loc="left", fontsize=10, color=INK)
    seen, handles = set(), []
    for rc in order:
        k = kind[rc]
        if k in seen:
            continue
        seen.add(k)
        handles.append(plt.Line2D([], [], marker="o", ls="none",
                                  markerfacecolor=CLASS_COLOUR.get(k, INK2),
                                  markeredgecolor=CLASS_COLOUR.get(k, INK2),
                                  markersize=5, label=CLASS_LABEL.get(k, k)))
    handles.append(plt.Line2D([], [], marker="o", ls="none", markerfacecolor="none",
                              markeredgecolor=INK, markersize=5,
                              label="open symbols: sampled from human faeces"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.08),
              ncol=3, fontsize=7.5, frameon=False, labelcolor=INK2)
    tidy(ax, grid_axis="x")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, "Figure_S4_pose_placement.%s" % ext), dpi=300)
    plt.close(fig)
    for rc in GUT:
        print("  Figure_S4  %-18s closest approach %.2f A over %d poses"
              % (rc, min(by[rc]), len(by[rc])))
    for rc in order:
        if kind[rc] == "pos_xtal":
            frac = sum(1 for v in by[rc] if v <= CUT) / len(by[rc])
            print("  Figure_S4  %-18s %.0f%% of poses placed catalytically"
                  % (rc, 100 * frac))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "supplementary", "figures"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    print("Writing supplementary figures to %s" % args.out)
    figure_s1(args.out)
    figure_s4(args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
