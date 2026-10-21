#!/usr/bin/env python3
"""Run AutoDock Vina over every receptor x ligand pair and collect the scores.

SN-38G is run with three different random seeds on every receptor, so that the
spread between seeds gives an honest idea of how much of any difference between
receptors is just search noise. The other ligands are run once.

Usage: run_dock.py <workdir>
"""
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

SEEDS = [20261021, 7, 999331]
EXH = 24
MODES = 9


def add_offsite(df, work):
    """A same-protein, wrong-place control: the E. coli GUS crystal with the
    box moved to a surface patch far from the catalytic site. If a score only
    reflects 'large protein with crevices', this box will reproduce it."""
    row = df[df.receptor == '3K4D_A'].iloc[0]
    site = np.array([row.cx, row.cy, row.cz])
    ca = []
    for l in open(os.path.join(work, 'receptors', '3K4D_A.pdb')):
        if l[:6] == 'ATOM  ' and l[12:16].strip() == 'CA':
            ca.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    ca = np.array(ca)
    d = np.linalg.norm(ca - site, axis=1)
    cen = ca[int(np.argmax(d))]
    # pull the box 6 A back towards the protein centre so it is not mostly solvent
    cen = cen + 6.0 * (ca.mean(0) - cen) / np.linalg.norm(ca.mean(0) - cen)
    new = row.copy()
    new.receptor = '3K4D_A_offsite'
    new.kind = 'neg_offsite'
    new.cx, new.cy, new.cz = [round(float(v), 3) for v in cen]
    new.site_source = f'surface patch {d.max():.0f} A from the catalytic site'
    new.note = 'same E. coli GUS crystal, box deliberately off the active site'
    return pd.concat([df, pd.DataFrame([new])], ignore_index=True)


def parse(out_txt):
    """Vina stdout -> list of affinities, best first."""
    aff = []
    started = False
    for l in out_txt.splitlines():
        if re.match(r'^\s*-----\+', l):
            started = True
            continue
        if started:
            m = re.match(r'^\s*(\d+)\s+(-?\d+\.?\d*)', l)
            if m:
                aff.append(float(m.group(2)))
            elif aff:
                break
    return aff


def main():
    work = sys.argv[1]
    vina = os.path.join(work, 'bin', 'vina')
    boxes = pd.read_csv(os.path.join(work, 'boxes.tsv'), sep='\t')
    boxes = add_offsite(boxes, work)
    boxes.to_csv(os.path.join(work, 'boxes_final.tsv'), sep='\t', index=False)
    ligs = pd.read_csv(os.path.join(work, 'ligands.tsv'), sep='\t')
    heavy = dict(zip(ligs.ligand, ligs.heavy_atoms))
    heavy['EVA_xtal'] = 13

    jobs = []
    for r in boxes.itertuples():
        for lig, seeds in [('SN38G', SEEDS), ('SN38', SEEDS[:1]),
                           ('GCU', SEEDS[:1])]:
            for s in seeds:
                jobs.append((r, lig, s))
        if r.receptor == '3K4D_A':
            for s in SEEDS:
                jobs.append((r, 'EVA_xtal', s))

    # controls first: if they do not separate, the candidate runs are pointless,
    # and that is worth knowing in ten minutes rather than at the end
    jobs.sort(key=lambda j: 0 if j[0].kind != 'candidate' else 1)

    pose_dir = os.path.join(work, 'poses')
    os.makedirs(pose_dir, exist_ok=True)
    log_dir = os.path.join(work, 'logs')
    os.makedirs(log_dir, exist_ok=True)

    # resume: the sandbox can be restarted mid-run, so keep the runs that
    # already finished and redo only what is missing
    raw_path = os.path.join(work, 'dock_raw.tsv')
    rows, done = [], set()
    if os.path.exists(raw_path):
        prev = pd.read_csv(raw_path, sep='\t')
        prev = prev[prev.ok.astype(str).str.lower().isin(['true', '1'])]
        rows = prev.to_dict('records')
        done = {(r['receptor'], r['ligand'], int(r['seed'])) for r in rows}
        print(f'resuming: {len(done)} runs already done', flush=True)

    for i, (r, lig, seed) in enumerate(jobs, 1):
        tag = f'{r.receptor}__{lig}__s{seed}'
        if (r.receptor, lig, seed) in done:
            continue
        pose = os.path.join(pose_dir, tag + '.pdbqt')
        cmd = [vina, '--receptor', r.receptor_pdbqt,
               '--ligand', os.path.join(work, 'ligands', lig + '.pdbqt'),
               '--center_x', str(r.cx), '--center_y', str(r.cy),
               '--center_z', str(r.cz),
               '--size_x', str(r.box), '--size_y', str(r.box),
               '--size_z', str(r.box),
               '--exhaustiveness', str(EXH), '--num_modes', str(MODES),
               '--seed', str(seed), '--cpu', '4', '--out', pose]
        p = subprocess.run(cmd, capture_output=True, text=True)
        open(os.path.join(log_dir, tag + '.log'), 'w').write(
            p.stdout + '\n--- stderr ---\n' + p.stderr)
        aff = parse(p.stdout)
        ok = p.returncode == 0 and bool(aff)
        rows.append(dict(
            receptor=r.receptor, kind=r.kind, family=r.family, ligand=lig,
            seed=seed, ok=ok,
            best=aff[0] if aff else np.nan,
            second=aff[1] if len(aff) > 1 else np.nan,
            mean_top3=float(np.mean(aff[:3])) if aff else np.nan,
            n_modes=len(aff),
            best_per_heavy=round(aff[0] / heavy[lig], 4) if aff else np.nan,
            cx=r.cx, cy=r.cy, cz=r.cz, box=r.box,
            exhaustiveness=EXH, vina='1.2.5'))
        print(f'[{i}/{len(jobs)}] {tag}  best='
              f'{aff[0] if aff else "FAIL"}', flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(work, 'dock_raw.tsv'),
                                  sep='\t', index=False)
    print('done', len(rows), 'runs')


if __name__ == '__main__':
    main()
