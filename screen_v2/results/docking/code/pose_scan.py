#!/usr/bin/env python3
"""Scan every docked pose, not just the top-scoring one.

The affinity score turned out not to separate the characterised
beta-glucuronidases from unrelated proteins, so the question becomes a
geometric one: does the pocket hold the glucuronide the way a glucuronidase
has to hold it, with the sugar's anomeric centre against the catalytic
carboxylates?

For each receptor this reports, over all nine output poses of all three seeds:
  - the smallest sugar-to-catalytic-oxygen distance found in any pose
  - the best-scoring pose that is catalytically placed (sugar within 4 A)
  - how that score compares with the receptor's own best pose overall

Usage: pose_scan.py <workdir> <out_dir>
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import SUGAR_SMARTS, sidechain_o, sugar_serials  # noqa: E402

CUTOFF = 4.0


def all_models(path):
    """-> list of (mode, affinity, serials, xyz) for the heavy atoms"""
    out = []
    mode = aff = None
    ser, xyz = [], []
    for l in open(path):
        if l.startswith('MODEL'):
            mode = int(l.split()[1])
            ser, xyz = [], []
        elif l.startswith('REMARK VINA RESULT'):
            aff = float(l.split()[3])
        elif l[:6] in ('ATOM  ', 'HETATM'):
            if l[77:79].strip().upper() in ('H', 'HD'):
                continue
            ser.append(int(l[6:11]))
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        elif l.startswith('ENDMDL'):
            out.append((mode, aff, np.array(ser), np.array(xyz)))
    return out


def main():
    work, out_dir = sys.argv[1:3]
    os.makedirs(out_dir, exist_ok=True)
    sites = pd.read_csv(os.path.join(work, 'work', 'site_boxes.tsv'), sep='\t')
    boxes = pd.read_csv(os.path.join(work, 'boxes_final.tsv'), sep='\t')
    raw = pd.read_csv(os.path.join(work, 'dock_raw.tsv'), sep='\t')

    cat = {r.target: [int(x) for x in str(r.cat_resnums).split(';') if x]
           for r in sites.itertuples()}
    cat['3K4D_A'] = [413, 504]
    cat['P05804_uidA_gus_w214-603'] = [413 - 213, 504 - 213]
    # the same test on the unrelated enzymes, so the geometric criterion has
    # negative controls of its own: their own annotated catalytic or
    # cofactor-contact residues, in each model's numbering
    # the other two gut-bacterial GUS crystals: the residues equivalent to
    # E. coli Glu413/Glu504, located by global pairwise alignment to P05804
    # (6D4O 47.3% identity, 6CXS 48.8%), both confirmed to be glutamate
    cat['6D4O_A'] = [425, 516]
    cat['6CXS_A'] = [412, 505]
    cat['P0ABF6_cdd_ecoli'] = [104]                       # cytidine deaminase
    cat['P0DTQ4_tyrDC_control_w231-620'] = [392 - 230, 420 - 230]
    cat['193L_A'] = [35, 52]                              # lysozyme
    cat['1TIM_A'] = [95, 165]                             # TIM
    cat['3K4D_A_offsite'] = [413, 504]

    sug = sugar_serials(os.path.join(work, 'ligands', 'SN38G.pdbqt'),
                        os.path.join(work, 'ligands', 'SN38G.sdf'))
    if not sug:
        sys.exit('could not identify the glucuronide atoms')

    kind = dict(zip(boxes.receptor, boxes.kind))
    rows, per_pose = [], []
    for rec, cnums in cat.items():
        # the off-site box is docked against the same prepared 3K4D structure
        pdb = os.path.join(work, 'receptors',
                           rec.replace('_offsite', '') + '.pdb')
        ox = sidechain_o(pdb, cnums)
        if not len(ox):
            continue
        seen = []
        for seed in sorted(raw[(raw.receptor == rec) &
                               (raw.ligand == 'SN38G')].seed.unique()):
            p = os.path.join(work, 'poses', f'{rec}__SN38G__s{seed}.pdbqt')
            if not os.path.exists(p):
                continue
            for mode, aff, ser, xyz in all_models(p):
                mask = np.isin(ser, sug)
                if not mask.any():
                    continue
                d = float(min(np.linalg.norm(xyz[mask] - o, axis=1).min()
                              for o in ox))
                seen.append((aff, d, seed, mode))
                per_pose.append(dict(receptor=rec, kind=kind.get(rec, ''),
                                     seed=seed, mode=mode, affinity=aff,
                                     sugar_to_catalytic=round(d, 2)))
        if not seen:
            continue
        seen.sort()
        best_any = seen[0]
        placed = [s for s in seen if s[1] <= CUTOFF]
        rows.append(dict(
            receptor=rec, kind=kind.get(rec, ''),
            n_poses=len(seen),
            best_affinity=round(best_any[0], 2),
            sugar_dist_of_best=round(best_any[1], 2),
            min_sugar_dist_any_pose=round(min(s[1] for s in seen), 2),
            n_poses_catalytic=len(placed),
            frac_poses_catalytic=round(len(placed) / len(seen), 2),
            best_catalytic_affinity=round(placed[0][0], 2) if placed else np.nan,
            penalty_vs_best=round(placed[0][0] - best_any[0], 2)
            if placed else np.nan,
            catalytic_pose=bool(placed)))

    df = pd.DataFrame(rows)
    order = {'pos_xtal': 0, 'pos_pred': 1, 'candidate': 2}
    df['_o'] = df.kind.map(lambda k: order.get(k, 3))
    df = df.sort_values(['_o', 'best_catalytic_affinity']).drop(columns='_o')
    df.to_csv(os.path.join(out_dir, 'pose_scan_sn38g.tsv'), sep='\t',
              index=False)
    pd.DataFrame(per_pose).to_csv(
        os.path.join(out_dir, 'pose_scan_all_poses.tsv'), sep='\t', index=False)
    pd.set_option('display.width', 220)
    print(df.to_string(index=False))


if __name__ == '__main__':
    main()
