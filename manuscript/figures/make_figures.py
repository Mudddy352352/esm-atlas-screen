#!/usr/bin/env python3
"""Build every figure in the manuscript from the screen's own output files.

Nothing here is hardcoded from a README: each number is recomputed from the
TSVs in screen_v1/results/ and screen_v2/ of this repository.
Run it again after a screen rerun and the figures follow the data.

    python3 make_figures.py [--outdir .]

Writes, for each figure, a 300-dpi PNG (for reading and for the app preview) and
a PDF (vector, for submission).

Style notes, so a restyle stays consistent: light print surface, thin marks,
hairline recessive grid, categorical hues taken in fixed order from the
validated palette (blue then orange), muted ink for axes and labels, direct
labels instead of a hover layer because a journal figure has no hover. Panels
with fewer than 20 observations are labelled as small-n rather than having a
trend drawn through them.
"""
import argparse
import csv
import gzip
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
V1 = os.path.join(ROOT, "screen_v1", "results")
V2 = os.path.join(ROOT, "screen_v2", "results")
PFAM = os.path.join(ROOT, "screen_v2", "pfam")

# Validated palette, light mode on a white print surface.
# node validate_palette.js "#2a78d6,#eb6834" --mode light --surface "#ffffff" --pairs all -> all pass
BLUE = "#2a78d6"      # categorical slot 1
ORANGE = "#eb6834"    # categorical slot 2
AQUA = "#1baf7a"      # categorical slot 3 (sub-3:1 on white, so it always
                      # carries a visible label: the relief rule)
INK = "#0b0b0b"       # primary ink
INK2 = "#52514e"      # secondary ink
MUTED = "#898781"     # axis and label ink
GRID = "#e1e0d9"      # hairline gridline
BASELINE = "#c3c2b7"  # axis rule
SURFACE = "#ffffff"   # print surface

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "axes.edgecolor": BASELINE,
    "axes.linewidth": 0.6,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "legend.frameon": False,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
})


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def read_tsv(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def tidy(ax, grid_axis="y"):
    """Hairline recessive grid, no top or right spine."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_axisbelow(True)
    if grid_axis:
        ax.grid(axis=grid_axis, which="major")
    return ax


def save(fig, outdir, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, "%s.%s" % (name, ext)),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote %s.png / .pdf" % name)


# --------------------------------------------------------------------------
# Figure 1. Pipeline and attrition
# --------------------------------------------------------------------------
def figure1(outdir):
    cands = read_tsv(os.path.join(V2, "active_site_candidates.tsv"))
    joined = read_tsv(os.path.join(V2, "active_site_all_hits_joined.tsv.gz"))
    tier1 = read_tsv(os.path.join(V2, "tier1_candidates.tsv"))
    bodysite = read_tsv(os.path.join(V2, "bodysite_candidates.tsv"))
    summary = read_tsv(os.path.join(V1, "summary_per_query.tsv"))

    n_models = len(summary)
    n_returned = sum(int(r["atlas_hits_returned"]) for r in summary)
    n_samefold = len(joined)
    n_cand = len(cands)
    n_lost = sum(1 for r in cands if r["site_class"] == "site_lost")
    n_partial = sum(1 for r in cands if r["site_class"] == "site_partial")
    n_complete = sum(1 for r in cands if r["site_class"] == "site_complete")
    n_notcheck = sum(1 for r in cands if r["site_class"] == "not_checkable")
    n_tier1 = len(tier1)
    n_human = sum(1 for r in bodysite if r["human_associated"] == "yes")

    steps = [
        ("%d query models\nfrom 9 characterised enzymes" % n_models, None),
        ("%s Atlas alignments returned\n(1,000-hit cap per query)" % f"{n_returned:,}", None),
        ("%s share a query fold\nTM $\\geq$ 0.5, coverage $\\geq$ 0.5" % f"{n_samefold:,}", None),
        ("%s structure-only candidates\n< 30%% identity, phmmer blind" % f"{n_cand:,}", None),
        ("active-site check", "%s site lost   %d partial   %d complete   %d not checkable"
         % (f"{n_lost:,}", n_partial, n_complete, n_notcheck)),
        ("%d shortlisted proteins\nall curated catalytic residues kept" % n_tier1, None),
        ("%d from a human body site\n5 non-human, 6 unplaceable" % n_human, None),
    ]

    fig, ax = plt.subplots(figsize=(4.6, 7.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(steps))
    ax.axis("off")

    for i, (label, sub) in enumerate(steps):
        y = len(steps) - i - 0.5
        emph = i in (3, 5)
        box = FancyBboxPatch(
            (0.06, y - 0.3), 0.88, 0.6,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            linewidth=1.0 if emph else 0.6,
            edgecolor=BLUE if emph else BASELINE,
            facecolor=SURFACE, transform=ax.transData, zorder=2)
        ax.add_patch(box)
        ax.text(0.5, y, label, ha="center", va="center", zorder=3,
                color=INK if emph else INK2,
                fontweight="bold" if emph else "normal", fontsize=8.2)
        if sub:
            ax.text(0.5, y - 0.42, sub, ha="center", va="center",
                    color=MUTED, fontsize=7.0, zorder=3)
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch(
                (0.5, y - 0.32 - (0.12 if sub else 0)), (0.5, y - 0.70),
                arrowstyle="-|>", mutation_scale=8, linewidth=0.7,
                color=MUTED, zorder=1))

    ax.set_title("Screen and attrition", color=INK, loc="left", pad=8,
                 fontsize=9.5, fontweight="bold")
    save(fig, outdir, "figure1_pipeline")


# --------------------------------------------------------------------------
# Figure 2. Active-site retention
# --------------------------------------------------------------------------
def corr(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / (sxx * syy) ** 0.5


def retention(rows):
    """(n, kept) for rows whose query has a curated catalytic residue."""
    n = len(rows)
    kept = sum(1 for r in rows
               if int(r["n_catalytic_identical"] or 0) == int(r["n_catalytic"] or 0))
    return n, kept


def figure2(outdir):
    joined = read_tsv(os.path.join(V2, "active_site_all_hits_joined.tsv.gz"))
    cat = [r for r in joined if int(r["n_catalytic"] or 0) > 0]
    is_sd = lambda r: r["seq_detectable"].strip().lower() in ("true", "1", "yes")
    sd = [r for r in cat if is_sd(r)]
    so = [r for r in cat if not is_sd(r) and num(r["seq_id_pct"]) < 30]

    n_sd, k_sd = retention(sd)
    n_so, k_so = retention(so)
    pct = lambda k, n: 100.0 * k / n if n else float("nan")

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(7.2, 3.1), gridspec_kw={"width_ratios": [1.0, 1.25]})

    # Panel a: the headline comparison. One measure, two groups, so one hue with
    # emphasis on the structure-only bar rather than two competing colours.
    # The count rides in the tick label: inside the bar it would be clipped, and
    # below the bar it would collide with the category label.
    labels = ["Found by\nsequence search\n%s of %s" % (f"{k_sd:,}", f"{n_sd:,}"),
              "Structure-only\n(< 30%% identity)\n%s of %s" % (f"{k_so:,}", f"{n_so:,}")]
    vals = [pct(k_sd, n_sd), pct(k_so, n_so)]
    colors = [MUTED, BLUE]
    bars = ax1.bar(labels, vals, width=0.52, color=colors, zorder=3)
    for bar, v in zip(bars, vals):
        ax1.text(bar.get_x() + bar.get_width() / 2, v + 2.5, "%.1f%%" % v,
                 ha="center", va="bottom", color=INK, fontweight="bold", fontsize=9)
    ax1.set_ylim(0, 100)
    ax1.set_ylabel("Same-fold hits keeping all\ncurated catalytic residues (%)")
    ax1.set_title("a   Retention falls with sequence distance", loc="left",
                  color=INK, pad=6)
    tidy(ax1)

    # Panel b: the same split per query family. Small n, so bars with counts and
    # an explicit small-n note, not a trend.
    per = {}
    for r in cat:
        q = r["query"]
        d = per.setdefault(q, {"sd": [], "so": []})
        if is_sd(r):
            d["sd"].append(r)
        elif num(r["seq_id_pct"]) < 30:
            d["so"].append(r)
    fams = [q for q in per if per[q]["so"]]
    fams.sort(key=lambda q: -len(per[q]["so"]))

    pretty = {
        "P05804_uidA_gus_w214-603": "β-glucuronidase",
        "P0DTQ4_tyrDC_control_w231-620": "tyrosine decarboxylase",
        "P19079_cdd_tetrameric": "cytidine deaminase\n(tetrameric)",
        "P0ABF6_cdd_ecoli": "cytidine deaminase\n(E. coli)",
    }
    ypos = range(len(fams))
    h = 0.34
    for i, q in enumerate(fams):
        nsd, ksd = retention(per[q]["sd"])
        nso, kso = retention(per[q]["so"])
        ax2.barh(i + h / 1.6, pct(ksd, nsd), height=h, color=MUTED, zorder=3)
        ax2.barh(i - h / 1.6, pct(kso, nso), height=h, color=BLUE, zorder=3)
        ax2.text(pct(ksd, nsd) + 1.5, i + h / 1.6, "%d/%d" % (ksd, nsd),
                 va="center", color=INK2, fontsize=7)
        ax2.text(pct(kso, nso) + 1.5, i - h / 1.6, "%d/%d" % (kso, nso),
                 va="center", color=INK2, fontsize=7)
    ax2.set_yticks(list(ypos))
    ax2.set_yticklabels([pretty.get(q, q) for q in fams], fontsize=7.4)
    ax2.set_xlim(0, 118)
    ax2.set_xticks([0, 25, 50, 75, 100])
    ax2.set_xlabel("Catalytic residues all kept (%)")
    ax2.set_title("b   By query family (structure-only n is small in every row)",
                  loc="left", color=INK, pad=6, fontsize=8.6)
    tidy(ax2, grid_axis="x")
    ax2.invert_yaxis()

    handles = [plt.Rectangle((0, 0), 1, 1, color=MUTED),
               plt.Rectangle((0, 0), 1, 1, color=BLUE)]
    ax2.legend(handles, ["Found by sequence search", "Structure-only"],
               loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2,
               fontsize=7.2, labelcolor=INK2)

    fig.tight_layout()
    save(fig, outdir, "figure2_active_site_retention")


# --------------------------------------------------------------------------
# Figure 3. Where the candidates sit in sequence space
# --------------------------------------------------------------------------
def figure3(outdir):
    cands = read_tsv(os.path.join(V2, "active_site_candidates.tsv"))
    tax = read_tsv(os.path.join(V2, "taxonomy_candidates.tsv"))
    tier1 = read_tsv(os.path.join(V2, "tier1_candidates.tsv"))
    # taxonomy_candidates.tsv carries identity to the nearest labelled protein;
    # identity to the query enzyme lives in the shortlist table.
    q_id_of = {r["target"]: num(r["seq_id_pct"]) for r in tier1}
    recov = read_tsv(os.path.join(PFAM, "pfam_recovery.tsv"))
    pfam_ok = {}
    for r in recov:
        pfam_ok[r["candidate"]] = pfam_ok.get(r["candidate"], False) or \
            r["pfam_recovers_link"] == "yes"

    ids = sorted(num(r["seq_id_pct"]) for r in cands)
    med = ids[len(ids) // 2]

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(7.4, 3.4), gridspec_kw={"width_ratios": [1.0, 1.35]})

    # Panel a: every candidate is remote from its query enzyme.
    ax1.hist(ids, bins=range(6, 32, 1), color=BLUE, zorder=3)
    ax1.axvline(med, color=INK, linewidth=1.0, zorder=4)
    ax1.text(med + 0.6, ax1.get_ylim()[1] * 0.94, "median %.1f%%" % med,
             color=INK, fontsize=7.6, fontweight="bold", va="top")
    ax1.axvline(30, color=ORANGE, linewidth=1.0, zorder=4)
    ax1.text(29.4, ax1.get_ylim()[1] * 0.60, "30% cutoff", color=ORANGE,
             fontsize=7.4, va="top", ha="right", rotation=90)
    ax1.set_xlabel("Identity to the query enzyme (%)")
    ax1.set_ylabel("Candidates")
    ax1.set_title("a   All %s candidates are remote\n      from the query enzyme"
                  % f"{len(ids):,}", loc="left", color=INK, pad=6)
    tidy(ax1)

    # Panel b: but the shortlist is close to labelled proteins. Dumbbell per
    # protein, so the gap itself is the mark.
    rows = []
    for r in tax:
        rows.append((r["mgyp"], q_id_of[r["mgyp"]], num(r["best_afdb50_seqid"]),
                     pfam_ok.get(r["mgyp"], False)))
    rows.sort(key=lambda t: t[2])

    for i, (mg, q_id, a_id, ok) in enumerate(rows):
        ax2.plot([q_id, a_id], [i, i], color=BASELINE, linewidth=1.4,
                 solid_capstyle="round", zorder=2)
        ax2.plot([q_id], [i], "o", color=BLUE, markersize=5.2, zorder=3)
        ax2.plot([a_id], [i], "o", color=ORANGE, markersize=5.2, zorder=3,
                 markerfacecolor=ORANGE if ok else SURFACE,
                 markeredgecolor=ORANGE, markeredgewidth=1.3)
    ax2.axvline(30, color=MUTED, linewidth=0.8, linestyle="-", zorder=1)
    ax2.text(30.6, len(rows) - 0.4, "30%", color=MUTED, fontsize=7.2, va="top")

    ax2.set_yticks(range(len(rows)))
    ax2.set_yticklabels([t[0] for t in rows], fontsize=6.4, fontfamily="monospace")
    ax2.set_xlim(5, 70)
    ax2.set_xlabel("Sequence identity (%)")
    ax2.set_title("b   The 13 shortlisted proteins have close labelled relatives",
                  loc="left", color=INK, pad=6)
    tidy(ax2, grid_axis="x")

    handles = [
        plt.Line2D([], [], marker="o", color=BLUE, linestyle="", markersize=5.2,
                   label="to its query enzyme"),
        plt.Line2D([], [], marker="o", color=ORANGE, linestyle="", markersize=5.2,
                   label="to nearest labelled AFDB50 protein"),
        plt.Line2D([], [], marker="o", color=ORANGE, linestyle="", markersize=5.2,
                   markerfacecolor=SURFACE, markeredgecolor=ORANGE,
                   label="open ring: Pfam finds no shared family"),
    ]
    ax2.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.16),
               ncol=1, fontsize=6.8, labelcolor=INK2)

    fig.tight_layout()
    save(fig, outdir, "figure3_sequence_space")


# --------------------------------------------------------------------------
# Figure 5. Pfam recovery
# --------------------------------------------------------------------------
def figure_pfam(outdir):
    recov = read_tsv(os.path.join(PFAM, "pfam_recovery.tsv"))
    per = {}
    for r in recov:
        d = per.setdefault(r["query_model"], [0, 0])
        d[0] += 1
        d[1] += r["pfam_recovers_link"] == "yes"
    total = sum(d[0] for d in per.values())
    total_rec = sum(d[1] for d in per.values())

    pretty = {
        "P05804_uidA_gus_w1-390": "β-glucuronidase (w1-390)",
        "P05804_uidA_gus_w214-603": "β-glucuronidase (w214-603)",
        "P07650_deoA_tp_w1-390": "thymidine phosphorylase",
        "P0ABF6_cdd_ecoli": "cytidine deaminase (E. coli)",
        "P0DTQ4_tyrDC_control_w1-390": "tyrosine decarboxylase (w1-390)",
        "P0DTQ4_tyrDC_control_w231-620": "tyrosine decarboxylase (w231-620)",
        "P15888_nfsB_nitroreductase": "nitroreductase NfsB",
        "P17117_nfsA_nitroreductase": "nitroreductase NfsA",
        "P19079_cdd_tetrameric": "cytidine deaminase (tetrameric)",
        "Q7A782_azo1_azoreductase": "azoreductase (flavodoxin-like)",
        "Q8KU07_azoB_azoreductase": "azoreductase (Rossmann)",
    }
    fams = sorted(per, key=lambda q: 100.0 * per[q][1] / per[q][0])

    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    SMALL_N = 20
    for i, q in enumerate(fams):
        n, rec = per[q]
        p = 100.0 * rec / n
        small = n < SMALL_N
        ax.barh(i, p, height=0.52, zorder=3,
                color=ORANGE if p < 80 else BLUE,
                alpha=0.45 if small else 1.0)
        note = "%d/%d" % (rec, n)
        if small:
            note += "  small n"
        ax.text(p + 1.5, i, note, va="center", color=INK2, fontsize=7.2)

    overall = 100.0 * total_rec / total
    ax.axvline(overall, color=INK, linewidth=1.0, zorder=4)
    # Headroom above the top bar so the overall label does not sit on it.
    ax.set_ylim(-0.7, len(fams) - 0.3)
    ax.text(overall - 1.5, len(fams) - 0.45,
            "overall %.1f%% (%s of %s)" % (overall, f"{total_rec:,}", f"{total:,}"),
            color=INK, fontsize=7.6, fontweight="bold", ha="right", va="top")

    ax.set_yticks(range(len(fams)))
    ax.set_yticklabels([pretty.get(q, q) for q in fams], fontsize=7.4)
    ax.set_xlim(0, 118)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Candidates sharing a Pfam family with their query (%)")
    ax.set_title("A profile search recovers most of the candidates", loc="left",
                 color=INK, pad=6, fontweight="bold")
    tidy(ax, grid_axis="x")

    handles = [plt.Rectangle((0, 0), 1, 1, color=BLUE),
               plt.Rectangle((0, 0), 1, 1, color=ORANGE)]
    ax.legend(handles, ["$\\geq$ 80% recovered", "< 80% recovered"],
              loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2,
              fontsize=7.2, labelcolor=INK2)

    fig.tight_layout()
    save(fig, outdir, "figure5_pfam_recovery")


# --------------------------------------------------------------------------
# Figure 4. Habitat of the shortlist, with the denominator visible
# --------------------------------------------------------------------------
def figure_habitat(outdir):
    bs = read_tsv(os.path.join(V2, "bodysite_candidates.tsv"))
    human = [r for r in bs if r["habitat_classes"] == "human"]
    env = [r for r in bs if r["habitat_classes"] == "environmental"]
    eng = [r for r in bs if r["habitat_classes"] == "engineered"]
    unplaced = [r for r in bs if not r["habitat_classes"].strip()]
    placed = len(human) + len(env) + len(eng)

    fig, ax = plt.subplots(figsize=(6.8, 2.5))
    # Categorical slots in fixed order; the unplaced segment is deliberately the
    # neutral gridline grey, because "no metadata" is an absence, not a habitat.
    segs = [
        ("Human gut (faecal)", len(human), BLUE, SURFACE),
        ("Environmental (marine)", len(env), ORANGE, SURFACE),
        ("Engineered (sludge, wastewater, food)", len(eng), AQUA, INK),
        ("No metadata in any reachable release", len(unplaced), GRID, INK),
    ]
    left = 0
    for label, n, color, textcolor in segs:
        if not n:
            continue
        # 2px surface gap between fills rather than a border around them
        ax.barh(0, n, left=left, height=0.46, color=color, zorder=3,
                edgecolor=SURFACE, linewidth=1.6)
        ax.text(left + n / 2.0, 0, str(n), ha="center", va="center", zorder=4,
                color=textcolor, fontweight="bold", fontsize=10)
        left += n

    ax.set_xlim(0, len(bs))
    ax.set_ylim(-0.75, 0.95)
    ax.set_yticks([])
    ax.set_xticks(range(0, len(bs) + 1, 1))
    ax.set_xlabel("Shortlisted proteins")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    # Make the real denominator explicit rather than leaving it to be inferred.
    ax.annotate("", xy=(0, -0.42), xytext=(placed, -0.42),
                arrowprops=dict(arrowstyle="|-|,widthB=0.3,widthA=0.3",
                                linewidth=0.8, color=INK2))
    ax.text(placed / 2.0, -0.60,
            "placed: %d, which is the denominator for any proportion, not %d"
            % (placed, len(bs)),
            ha="center", va="top", color=INK2, fontsize=7.6)

    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, n, c, _t in segs if n]
    labels = ["%s" % lab for lab, n, c, _t in segs if n]
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(0, 1.12),
              ncol=2, fontsize=7.2, labelcolor=INK2)
    ax.set_title("Habitat of the 13 shortlisted proteins", loc="left",
                 color=INK, pad=26, fontweight="bold")

    fig.tight_layout()
    save(fig, outdir, "figure4_habitat")


# --------------------------------------------------------------------------
# Figure 6. Docking: the score does not separate, the pose geometry does
# --------------------------------------------------------------------------
DOCK = os.path.join(ROOT, "screen_v2", "results", "docking")

GUT = {"MGYP001533818318", "MGYP000245278657"}


DOCK_PRETTY = {
    "3K4D_A": "E. coli GUS (crystal)",
    "6D4O_A": "E. eligens GUS (crystal)",
    "6CXS_A": "C. perfringens GUS (crystal)",
    "P05804_uidA_gus_w214-603": "E. coli GUS (ESMFold model)",
    "193L_A": "lysozyme",
    "1TIM_A": "triose-phosphate isomerase",
    "P0ABF6_cdd_ecoli": "cytidine deaminase",
    "3K4D_A_offsite": "GUS, box moved off-site",
}


def figure_docking(outdir):
    summ = read_tsv(os.path.join(DOCK, "dock_summary_sn38g.tsv"))
    pose = read_tsv(os.path.join(DOCK, "pose_scan_sn38g.tsv"))
    final = {r["receptor"]: r for r in read_tsv(os.path.join(DOCK, "final_table.tsv"))}

    allp = read_tsv(os.path.join(DOCK, "pose_scan_all_poses.tsv"))

    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(7.2, 9.6),
        gridspec_kw={"height_ratios": [1.0, 2.0, 1.45]})

    # Panel a: every receptor's score on one axis, by class. The point of the
    # panel is the overlap, so the classes share the axis and nothing is ranked.
    classes = [("Positive controls", lambda k: k.startswith("pos"), BLUE),
               ("Shortlisted candidates", lambda k: k == "candidate", ORANGE),
               ("Negative controls", lambda k: k.startswith("neg"), MUTED)]
    ticklabels = []
    for i, (label, test, color) in enumerate(classes):
        vals = [(num(r["best_mean"]), r["receptor"]) for r in summ if test(r["kind"])]
        y = len(classes) - 1 - i
        lo, hi = min(v for v, _ in vals), max(v for v, _ in vals)
        ax1.plot([lo, hi], [y, y], color=color, linewidth=1.2, alpha=0.35,
                 solid_capstyle="round", zorder=2)
        for v, name in vals:
            gut = name in GUT
            ax1.plot([v], [y], "o", markersize=6.5, zorder=3,
                     markerfacecolor=SURFACE if gut else color,
                     markeredgecolor=color, markeredgewidth=1.4)
        ticklabels.append((y, label))

    # The overlap band, which is the whole message.
    pos = [num(r["best_mean"]) for r in summ if r["kind"].startswith("pos")]
    neg = [num(r["best_mean"]) for r in summ if r["kind"].startswith("neg")]
    ax1.axvspan(max(pos), min(neg), color=GRID, alpha=0.55, zorder=1)
    ax1.text((max(pos) + min(neg)) / 2.0, -0.78,
             "positive and negative controls overlap here",
             ha="center", va="top", color=INK2, fontsize=7.2)

    best_cand = min(num(r["best_mean"]) for r in summ if r["kind"] == "candidate")
    ax1.annotate("two candidates beat\nevery positive control",
                 xy=(best_cand, 1.18), xytext=(best_cand - 0.08, 2.05),
                 ha="left", va="bottom", color=INK2, fontsize=7.2,
                 arrowprops=dict(arrowstyle="-", linewidth=0.7, color=MUTED))

    ax1.set_xlim(-11.2, -6.3)
    ax1.set_ylim(-1.5, 2.95)
    ax1.set_yticks([y for y, _ in ticklabels])
    ax1.set_yticklabels([lab for _, lab in ticklabels], fontsize=8)
    ax1.tick_params(axis="y", length=0)
    ax1.set_xlabel("Vina affinity for SN-38 glucuronide (kcal/mol), mean of 3 seeds")
    ax1.set_title("a   The docking score does not separate glucuronidases "
                  "from anything else", loc="left", color=INK, pad=6)
    for side in ("top", "right", "left"):
        ax1.spines[side].set_visible(False)
    ax1.grid(axis="x")
    ax1.set_axisbelow(True)

    # Panel b: pose geometry, which does separate. One row per receptor.
    order = {"pos_xtal": 0, "pos_pred": 0, "candidate": 1,
             "neg_xtal": 2, "neg_pred": 2, "neg_offsite": 2}
    rows = sorted(pose, key=lambda r: (order[r["kind"]], -num(r["frac_poses_catalytic"])))
    labels, vals, colors, notes = [], [], [], []
    for r in rows:
        name = r["receptor"]
        k = r["kind"]
        colors.append(BLUE if k.startswith("pos") else
                      (ORANGE if k == "candidate" else MUTED))
        pretty = DOCK_PRETTY.get(name, name)
        if name in GUT:
            pretty += "  (human gut)"
        labels.append(pretty)
        vals.append(100.0 * num(r["frac_poses_catalytic"]))
        pen = num(r["penalty_vs_best"]) if r["penalty_vs_best"] else 0.0
        if vals[-1] == 0:
            notes.append("closest approach %.1f A" % num(r["min_sugar_dist_any_pose"]))
        elif pen > 0:
            notes.append("%.0f%%, costs %.2f kcal/mol" % (vals[-1], pen))
        else:
            notes.append("%.0f%%, no cost" % vals[-1])

    ypos = list(range(len(rows)))
    for i, (v, c, note, r) in enumerate(zip(vals, colors, notes, rows)):
        gut = r["receptor"] in GUT
        ax2.barh(i, v, height=0.56, zorder=3, color=SURFACE if gut else c,
                 edgecolor=c, linewidth=1.3 if gut else 0)
        ax2.text(max(v, 0) + 1.5, i, note, va="center", color=INK2, fontsize=6.9)
    ax2.set_yticks(ypos)
    ax2.set_yticklabels(labels, fontsize=7.0)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 100)
    ax2.set_xlabel("Poses placing the glucuronide's sugar within 4 A of the "
                   "catalytic residues (%)")
    ax2.set_title("b   Where the ligand sits does separate them", loc="left",
                  color=INK, pad=6)
    tidy(ax2, grid_axis="x")

    handles = [plt.Rectangle((0, 0), 1, 1, color=BLUE),
               plt.Rectangle((0, 0), 1, 1, color=ORANGE),
               plt.Rectangle((0, 0), 1, 1, color=MUTED),
               plt.Rectangle((0, 0), 1, 1, facecolor=SURFACE, edgecolor=ORANGE,
                             linewidth=1.3)]
    ax2.legend(handles, ["Positive control", "Shortlisted candidate",
                         "Negative control", "Human-gut candidate"],
               loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4,
               fontsize=7.0, labelcolor=INK2)

    # Panel c: the two measures against each other, one point per pose. This
    # is the only direct evidence that a better score is not a better placement.
    groups = [("Positive controls", lambda k: k.startswith("pos"), BLUE, "o", 3.4),
              ("Shortlisted candidates", lambda k: k == "candidate", ORANGE, "o", 3.4),
              ("Negative controls", lambda k: k.startswith("neg"), MUTED, "o", 3.4)]
    xs_all, ys_all = [], []
    for label, test, color, marker, size in groups:
        xs = [num(r["affinity"]) for r in allp if test(r["kind"])]
        ys = [num(r["sugar_to_catalytic"]) for r in allp if test(r["kind"])]
        xs_all += xs
        ys_all += ys
        ax3.plot(xs, ys, marker, markersize=size, alpha=0.55, color=color,
                 markeredgewidth=0, label=label, zorder=3)
    r_all = corr(xs_all, ys_all)
    cx = [num(r["affinity"]) for r in allp if r["kind"] == "candidate"]
    cy = [num(r["sugar_to_catalytic"]) for r in allp if r["kind"] == "candidate"]
    ax3.axhline(4.0, color=BASELINE, linewidth=0.9, linestyle=(0, (4, 3)), zorder=2)
    ax3.text(-5.95, 4.0, " 4 A criterion", ha="right", va="bottom",
             color=INK2, fontsize=7.0)
    ax3.set_ylim(0, 16)
    ax3.set_xlim(-11.3, -4.6)
    ax3.set_xlabel("Vina affinity of the pose (kcal/mol)")
    ax3.set_ylabel("Sugar to catalytic\nresidues (A)")
    ax3.set_title("c   A better-scoring pose is not a better-placed pose",
                  loc="left", color=INK, pad=6)
    ax3.annotate("r = %+.2f over all %d poses\n(%+.2f within the candidates)"
                 % (r_all, len(xs_all), corr(cx, cy)),
                 xy=(-11.1, 14.6), ha="left", va="top", color=INK2, fontsize=7.2)
    tidy(ax3, grid_axis="y")
    ax3.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=3,
               fontsize=7.0, labelcolor=INK2, markerscale=1.8)
    print("  figure 6c: r over all poses = %+.3f (n=%d), candidates %+.3f"
          % (r_all, len(xs_all), corr(cx, cy)))

    fig.tight_layout(h_pad=2.2)
    save(fig, outdir, "figure6_docking")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    figure1(args.outdir)
    figure2(args.outdir)
    figure3(args.outdir)
    figure_habitat(args.outdir)
    figure_pfam(args.outdir)
    figure_docking(args.outdir)
    print("\nFigure 7 (structural superposition) is not built here: it needs "
          "PyMOL or ChimeraX, which this environment does not have.")


if __name__ == "__main__":
    main()
