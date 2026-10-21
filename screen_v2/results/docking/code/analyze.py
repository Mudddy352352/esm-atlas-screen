#!/usr/bin/env python3
"""Turn the raw Vina output into something that can be judged.

Three things are computed:

1. Seed spread. SN-38G was docked three times per receptor with different
   random seeds, so the spread across seeds says how much of a difference
   between receptors is search noise.

2. Control separation. Unless the characterised gut beta-glucuronidases score
   better than the unrelated proteins and the off-site box, the run says
   nothing and no candidate ranking should be read from it.

3. Where the pose actually sits. A score only means the ligand found somewhere
   comfortable inside a 26 A box. For each best pose this measures how far the
   glucuronide sugar ring is from the residues the screen identified as
   catalytic. A pose that parks the drug end in the pocket and leaves the sugar
   in solvent is not a glucuronidase pose, whatever it scores.

Usage: analyze.py <workdir> <out_dir>
"""
import os
import sys

import numpy as np
import pandas as pd

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "..", "..", ".."))
V2R = os.path.join(ROOT, 'screen_v2', 'results')

# the glucuronide sugar: pyranose ring plus its carboxylate
SUGAR_SMARTS = 'OC1OC(C(=O)O)C(O)C(O)C1O'


def pose_atoms(path, model=1):
    """heavy atoms of one MODEL of a Vina output pdbqt -> (serials, xyz)"""
    ser, xyz = [], []
    cur = 0
    for l in open(path):
        if l.startswith('MODEL'):
            cur = int(l.split()[1])
        elif l[:6] in ('ATOM  ', 'HETATM') and cur in (0, model):
            if l[77:79].strip().upper() in ('H', 'HD'):
                continue
            ser.append(int(l[6:11]))
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(ser), np.array(xyz)


def sugar_serials(lig_pdbqt, lig_sdf):
    """Which pdbqt atom serials belong to the glucuronide ring.

    Matched by coordinate, because Open Babel may reorder atoms when it builds
    the torsion tree, so the SDF index order cannot be assumed.
    """
    from rdkit import Chem
    m = Chem.MolFromMolFile(lig_sdf, removeHs=True)
    if m is None:
        return None
    patt = Chem.MolFromSmarts(SUGAR_SMARTS)
    match = m.GetSubstructMatch(patt)
    if not match:
        return None
    conf = m.GetConformer()
    want = np.array([list(conf.GetAtomPosition(i)) for i in match])
    ser, xyz = pose_atoms(lig_pdbqt)
    out = []
    for w in want:
        d = np.linalg.norm(xyz - w, axis=1)
        if d.min() < 0.15:
            out.append(int(ser[int(np.argmin(d))]))
    return out if len(out) >= len(match) - 1 else None


def sidechain_o(pdb_path, resnums):
    """Carboxylate / side-chain oxygens of the given residues."""
    want = set(resnums)
    pts = []
    for l in open(pdb_path):
        if l[:6] != 'ATOM  ':
            continue
        if int(l[22:26]) not in want:
            continue
        name = l[12:16].strip()
        if name in ('OE1', 'OE2', 'OD1', 'OD2'):
            pts.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(pts)


def rmsd_to_crystal(pose_pdbqt, xtal_pdb):
    """Symmetry-blind RMSD of the redocked pose to the crystal pose.

    Nearest-atom matching per crystal atom, so it is a lower bound rather than
    a true RMSD; good enough to tell a recovered pose from a wrong one.
    """
    _, pose = pose_atoms(pose_pdbqt, model=1)
    xtal = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(xtal_pdb)
                     if l[:6] == 'HETATM' and l[76:78].strip() != 'H'])
    if not len(pose) or not len(xtal):
        return np.nan
    d = [np.linalg.norm(pose - x, axis=1).min() for x in xtal]
    return float(np.sqrt(np.mean(np.square(d))))


def main():
    work, out_dir = sys.argv[1:3]
    os.makedirs(out_dir, exist_ok=True)
    raw = pd.read_csv(os.path.join(work, 'dock_raw.tsv'), sep='\t')
    boxes = pd.read_csv(os.path.join(work, 'boxes_final.tsv'), sep='\t')
    sites = pd.read_csv(os.path.join(work, 'work', 'site_boxes.tsv'), sep='\t')
    tier1 = pd.read_csv(
        os.path.join(V2R, 'tier1_candidates.tsv'), sep='\t')
    body = pd.read_csv(
        os.path.join(V2R, 'bodysite_candidates.tsv'), sep='\t')

    # ---- 1. per receptor x ligand, across seeds ----
    g = (raw[raw.ok].groupby(['receptor', 'kind', 'family', 'ligand'])
         .agg(n_seeds=('seed', 'nunique'),
              best_mean=('best', 'mean'),
              best_sd=('best', 'std'),
              best_min=('best', 'min'),
              best_max=('best', 'max'),
              per_heavy=('best_per_heavy', 'mean'))
         .reset_index().round(3))
    g.to_csv(os.path.join(out_dir, 'dock_by_receptor_ligand.tsv'),
             sep='\t', index=False)

    # ---- 2. pose geometry for SN38G ----
    cat = {r.target: [int(x) for x in str(r.cat_resnums).split(';') if x]
           for r in sites.itertuples()}
    # the same measurement on the controls whose catalytic numbering is known:
    # E. coli GUS Glu413 (proton donor) and Glu504 (nucleophile), and the same
    # two residues in the ESMFold model, which starts at residue 214
    cat['3K4D_A'] = [413, 504]
    cat['P05804_uidA_gus_w214-603'] = [413 - 213, 504 - 213]
    sug = sugar_serials(os.path.join(work, 'ligands', 'SN38G.pdbqt'),
                        os.path.join(work, 'ligands', 'SN38G.sdf'))
    geo = []
    for r in boxes.itertuples():
        if r.receptor not in cat:
            continue
        rec_pdb = os.path.join(work, 'receptors', r.receptor + '.pdb')
        ox = sidechain_o(rec_pdb, cat[r.receptor])
        for seed in sorted(raw[(raw.receptor == r.receptor) &
                               (raw.ligand == 'SN38G')].seed.unique()):
            p = os.path.join(work, 'poses',
                             f'{r.receptor}__SN38G__s{seed}.pdbqt')
            if not os.path.exists(p) or not len(ox):
                continue
            ser, xyz = pose_atoms(p, model=1)
            all_d = min(np.linalg.norm(xyz - o, axis=1).min() for o in ox)
            sd = np.nan
            if sug:
                mask = np.isin(ser, sug)
                if mask.any():
                    sd = float(min(np.linalg.norm(xyz[mask] - o, axis=1).min()
                                   for o in ox))
            geo.append(dict(receptor=r.receptor, kind=r.kind, seed=seed,
                            cat_resnums=';'.join(map(str, cat[r.receptor])),
                            min_dist_any_atom=round(float(all_d), 2),
                            min_dist_sugar=round(sd, 2) if sd == sd else np.nan))
    geo = pd.DataFrame(geo)
    if len(geo):
        geo.to_csv(os.path.join(out_dir, 'pose_geometry_sn38g.tsv'),
                   sep='\t', index=False)

    # ---- 3. redocking check ----
    red = []
    for seed in raw[raw.ligand == 'EVA_xtal'].seed.unique():
        p = os.path.join(work, 'poses', f'3K4D_A__EVA_xtal__s{seed}.pdbqt')
        if os.path.exists(p):
            red.append(dict(seed=seed, rmsd_lb=round(rmsd_to_crystal(
                p, os.path.join(work, 'ligands', 'EVA_xtal.pdb')), 2)))
    red = pd.DataFrame(red)
    if len(red):
        red.to_csv(os.path.join(out_dir, 'redock_check.tsv'),
                   sep='\t', index=False)

    # ---- 4. the headline table ----
    sn = g[g.ligand == 'SN38G'].copy()
    sna = g[g.ligand == 'SN38'].set_index('receptor').best_mean
    gcu = g[g.ligand == 'GCU'].set_index('receptor').best_mean
    sn['sn38_aglycone'] = sn.receptor.map(sna)
    sn['glucuronic_acid'] = sn.receptor.map(gcu)
    if len(geo):
        gm = geo.groupby('receptor')[['min_dist_sugar',
                                      'min_dist_any_atom']].min().round(2)
        sn = sn.merge(gm, on='receptor', how='left')
    t1 = tier1.set_index('target')
    bcol = body.columns[0]
    bs = body.set_index(bcol)
    bcoln = [c for c in body.columns if 'biome' in c.lower()
             or 'body' in c.lower() or 'lineage' in c.lower()]
    sn['site_class'] = sn.receptor.map(
        lambda x: t1.at[x, 'site_class'] if x in t1.index else '')
    sn['site_class_3d'] = sn.receptor.map(
        lambda x: t1.at[x, 'site_class_3d'] if x in t1.index else '')
    sn['tm_qnorm'] = sn.receptor.map(
        lambda x: t1.at[x, 'tm_score_qnorm'] if x in t1.index else np.nan)
    if bcoln:
        sn['habitat'] = sn.receptor.map(
            lambda x: bs.at[x, bcoln[0]] if x in bs.index else '')
    order = {'pos_xtal': 0, 'pos_pred': 1, 'candidate': 2,
             'neg_pred': 3, 'neg_xtal': 4, 'neg_offsite': 5}
    sn['_o'] = sn.kind.map(order)
    sn = sn.sort_values(['_o', 'best_mean']).drop(columns=['_o', 'ligand'])
    sn.to_csv(os.path.join(out_dir, 'dock_summary_sn38g.tsv'),
              sep='\t', index=False)

    # ---- 5. the verdict on the controls ----
    pos = sn[sn.kind.isin(['pos_xtal', 'pos_pred'])].best_mean
    neg = sn[sn.kind.isin(['neg_pred', 'neg_xtal', 'neg_offsite'])].best_mean
    cand = sn[sn.kind == 'candidate'].best_mean
    noise = float(raw[(raw.ligand == 'SN38G') & raw.ok]
                  .groupby('receptor').best.std().mean())
    verdict = dict(
        n_runs=int(raw.ok.sum()), n_failed=int((~raw.ok).sum()),
        mean_seed_sd=round(noise, 3),
        pos_mean=round(float(pos.mean()), 2), pos_best=round(float(pos.min()), 2),
        neg_mean=round(float(neg.mean()), 2), neg_best=round(float(neg.min()), 2),
        gap_pos_minus_neg=round(float(pos.mean() - neg.mean()), 2),
        pos_worst_vs_neg_best=round(float(pos.max() - neg.min()), 2),
        cand_mean=round(float(cand.mean()), 2),
        cand_best=round(float(cand.min()), 2),
        separation=('controls separate from negatives'
                    if pos.max() < neg.min() else
                    'controls overlap the negatives'))
    pd.DataFrame([verdict]).to_csv(os.path.join(out_dir, 'control_check.tsv'),
                                   sep='\t', index=False)

    pd.set_option('display.width', 250)
    print('=== control check ===')
    for k, v in verdict.items():
        print(f'  {k}: {v}')
    print('\n=== SN-38G, best affinity averaged over 3 seeds (kcal/mol) ===')
    cols = ['receptor', 'kind', 'best_mean', 'best_sd', 'sn38_aglycone',
            'glucuronic_acid', 'min_dist_sugar', 'site_class_3d']
    print(sn[[c for c in cols if c in sn.columns]].to_string(index=False))
    if len(red):
        print('\n=== redocking the crystal ligand into its own site ===')
        print(red.to_string(index=False))


if __name__ == '__main__':
    main()
