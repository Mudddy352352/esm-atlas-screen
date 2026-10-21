#!/usr/bin/env python3
"""Prepare receptors, ligands and docking boxes for the SN-38G docking run.

Receptor classes
  candidate  13 tier-1 Atlas proteins, ESMFold predictions, box on the
             catalytic residues located by map_sites.py
  pos_pred   ESMFold model of the E. coli beta-glucuronidase query itself:
             a positive control with the same provenance as the candidates
  pos_xtal   crystal structures of characterised gut bacterial
             beta-glucuronidases, each with a ligand bound in the catalytic
             site, box centred on that ligand
  neg_xtal   crystal structures of unrelated enzymes, box on their own
             catalytic residues
  neg_pred   ESMFold models of the screen's other query families, which are
             real enzymes but not glucuronidases

Everything is written under <out>/receptors, <out>/ligands, <out>/boxes.tsv.
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

AA3 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V', 'MSE': 'M'}

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", ".."))
V1 = os.path.join(ROOT, 'screen_v1')
V2R = os.path.join(ROOT, 'screen_v2', 'results')

# crystal positives: (pdb, chain, het code of the active-site ligand, note)
POS_XTAL = [
    ('3K4D', 'A', 'EVA', 'E. coli beta-glucuronidase + glucaro-1,5-lactam'),
    ('6D4O', 'A', 'FUV', 'Eubacterium eligens beta-glucuronidase + ligand'),
    ('6CXS', 'A', 'FJV', 'Clostridium perfringens beta-glucuronidase + ligand'),
]
# crystal negatives: (pdb, chain, [(resnum, expected one-letter)], note)
NEG_XTAL = [
    ('193L', 'A', [(35, 'E'), (52, 'D')], 'hen egg lysozyme, catalytic cleft'),
    ('1TIM', 'A', [(95, 'H'), (165, 'E')], 'chicken triose-phosphate isomerase'),
]
# ESMFold query models used as controls: (model name, uniprot acc, class, note)
PRED_CTRL = [
    ('P05804_uidA_gus_w214-603', 'P05804', 'pos_pred',
     'ESMFold model of the E. coli GUS catalytic domain (the screen query)'),
    ('P0ABF6_cdd_ecoli', 'P0ABF6', 'neg_pred',
     'ESMFold model of E. coli cytidine deaminase'),
    ('Q8KU07_azoB_azoreductase', 'Q8KU07', 'neg_pred',
     'ESMFold model of an NAD(P)H azoreductase'),
    ('P0DTQ4_tyrDC_control_w231-620', 'P0DTQ4', 'neg_pred',
     'ESMFold model of L-tyrosine decarboxylase'),
]

LIGANDS = {
    # ChEBI:8990, formula C28H28N2O11 -- the substrate of interest
    'SN38G': 'CCc1c2c(nc3ccc(O[C@@H]4O[C@H](C(=O)O)[C@@H](O)[C@H](O)[C@H]4O)cc13)'
             '-c1cc3c(c(=O)n1C2)COC(=O)[C@]3(O)CC',
    # ChEBI:8988, the aglycone product
    'SN38': 'CCc1c2c(nc3ccc(O)cc13)-c1cc3c(c(=O)n1C2)COC(=O)[C@]3(O)CC',
    # ChEBI:28860, the sugar alone
    'GCU': 'O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O',
}


def read_pdb(path):
    return [fix_atom_name(l) for l in open(path)]


def fix_atom_name(l):
    """Put the atom name back in its standard PDB columns.

    The ESM Atlas and ESMFold files left-justify the atom name from column 13,
    so a serine OG reads as the element Og and Open Babel then writes an atom
    type Vina rejects. Single-letter elements belong in column 14.
    """
    if l[:6] not in ('ATOM  ', 'HETATM'):
        return l
    name = l[12:16].strip()
    elem = l[76:78].strip()
    if len(elem) == 1 and len(name) < 4 and not l[12:13].isspace():
        name = ' ' + name.ljust(3)
    else:
        name = l[12:16]
    return l[:12] + name[:4] + l[16:]


def clean_crystal(lines, chain, keep_het=()):
    """One chain, waters and ligands dropped, MSE rewritten as MET."""
    out, seen = [], set()
    for l in lines:
        if l[:6] not in ('ATOM  ', 'HETATM'):
            continue
        if l[21] != chain:
            continue
        res = l[17:20].strip()
        if l[:6] == 'HETATM' and res not in ('MSE',) and res not in keep_het:
            continue
        if res == 'MSE':
            l = 'ATOM  ' + l[6:17] + 'MET' + l[20:]
            if l[12:16].strip() == 'SE':
                l = l[:12] + ' SD ' + l[16:76] + ' S' + l[78:]
        if l[16] not in (' ', 'A'):          # keep altloc A only
            continue
        l = l[:16] + ' ' + l[17:]
        key = (l[17:27], l[12:16])
        if key in seen:
            continue
        seen.add(key)
        out.append(l)
    return out


def ca_table(lines):
    """-> dict resnum -> (one-letter, xyz), in file order"""
    d = {}
    for l in lines:
        if l[:6] in ('ATOM  ', 'HETATM') and l[12:16].strip() == 'CA':
            d[int(l[22:26])] = (AA3.get(l[17:20].strip(), 'X'),
                                np.array([float(l[30:38]), float(l[38:46]),
                                          float(l[46:54])]))
    return d


def het_centroid(lines, chain, code):
    xyz = [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
           for l in lines if l[:6] == 'HETATM' and l[21] == chain
           and l[17:20].strip() == code and l[76:78].strip() != 'H']
    if not xyz:
        raise SystemExit(f'no {code} in chain {chain}')
    return np.mean(xyz, 0), len(xyz)


def to_pdbqt_receptor(pdb_in, pdbqt_out):
    subprocess.run(['obabel', pdb_in, '-O', pdbqt_out, '-xr', '-p', '7.4'],
                   check=True, capture_output=True)


def build_ligand(name, smiles, out_dir):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise SystemExit('bad SMILES for ' + name)
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 20261021
    cids = AllChem.EmbedMultipleConfs(m, numConfs=20, params=ps)
    if not cids:
        raise SystemExit('embedding failed for ' + name)
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    best = int(np.argmin([e for _, e in res]))
    sdf = os.path.join(out_dir, name + '.sdf')
    w = Chem.SDWriter(sdf)
    w.write(m, confId=cids[best])
    w.close()
    pdbqt = os.path.join(out_dir, name + '.pdbqt')
    # -p 7.4 puts the glucuronide carboxylate in its physiological ionised form
    subprocess.run(['obabel', sdf, '-O', pdbqt, '-p', '7.4'],
                   check=True, capture_output=True)
    heavy = m.GetNumHeavyAtoms()
    tors = sum(1 for l in open(pdbqt) if l.startswith('BRANCH'))
    return dict(ligand=name, smiles=smiles, heavy_atoms=heavy,
                formula=Chem.rdMolDescriptors.CalcMolFormula(
                    Chem.RemoveHs(m)), n_torsions=tors, sdf=sdf, pdbqt=pdbqt)


def main():
    out = sys.argv[1]
    raw_pdb = os.path.join(out, 'pdb_raw')
    rec_dir = os.path.join(out, 'receptors')
    lig_dir = os.path.join(out, 'ligands')
    for d in (rec_dir, lig_dir):
        os.makedirs(d, exist_ok=True)

    # ---------- ligands ----------
    linfo = [build_ligand(n, s, lig_dir) for n, s in LIGANDS.items()]
    pd.DataFrame(linfo).to_csv(os.path.join(out, 'ligands.tsv'),
                               sep='\t', index=False)
    print(pd.DataFrame(linfo)[['ligand', 'formula', 'heavy_atoms',
                               'n_torsions']].to_string(index=False))

    rows = []

    # ---------- candidates ----------
    boxes = pd.read_csv(os.path.join(out, 'work', 'site_boxes.tsv'), sep='\t')
    tier1 = pd.read_csv(os.path.join(V2R, 'tier1_candidates.tsv'), sep='\t')
    body = pd.read_csv(os.path.join(V2R, 'bodysite_candidates.tsv'), sep='\t')
    t1 = tier1.set_index('target')
    bs = body.set_index(body.columns[0])
    for r in boxes.itertuples():
        lines = [l for l in read_pdb(r.pdb) if l[:6] == 'ATOM  ']
        pdb_out = os.path.join(rec_dir, r.target + '.pdb')
        open(pdb_out, 'w').writelines(lines + ['END\n'])
        pq = os.path.join(rec_dir, r.target + '.pdbqt')
        to_pdbqt_receptor(pdb_out, pq)
        fam = ('GUS' if 'P05804' in r.query else 'CDD')
        rows.append(dict(
            receptor=r.target, kind='candidate', family=fam,
            receptor_pdbqt=pq, cx=r.cx, cy=r.cy, cz=r.cz,
            box=26.0, site_source=f'catalytic residues {r.cat_resnums}',
            note=f'{fam}-fold Atlas candidate; site_class='
                 f'{t1.at[r.target, "site_class"]}'))

    # ---------- predicted controls ----------
    sites = pd.read_csv(os.path.join(V2R, 'query_sites.tsv'), sep='\t')
    for name, acc, kind, note in PRED_CTRL:
        src = os.path.join(V1, 'structures', name + '.pdb')
        lines = [l for l in read_pdb(src) if l[:6] == 'ATOM  ']
        pdb_out = os.path.join(rec_dir, name + '.pdb')
        open(pdb_out, 'w').writelines(lines + ['END\n'])
        pq = os.path.join(rec_dir, name + '.pdbqt')
        to_pdbqt_receptor(pdb_out, pq)
        import re
        m = re.search(r'_w(\d+)-(\d+)$', name)
        off = int(m.group(1)) if m else 1
        cas = ca_table(lines)
        s = sites[sites.acc == acc]
        cat = s[s.kind.isin(['Active site', 'Modified residue'])]
        use = cat if len(cat) else s
        pts, tags = [], []
        for rr in use.itertuples():
            wp = int(rr.pos) - off + 1
            if wp in cas and cas[wp][0] == rr.res:
                pts.append(cas[wp][1])
                tags.append(f'{rr.res}{int(rr.pos)}')
        if not pts:                      # fall back to all annotated sites
            for rr in s.itertuples():
                wp = int(rr.pos) - off + 1
                if wp in cas and cas[wp][0] == rr.res:
                    pts.append(cas[wp][1])
                    tags.append(f'{rr.res}{int(rr.pos)}')
        c = np.mean(pts, 0)
        rows.append(dict(receptor=name, kind=kind, family=acc,
                         receptor_pdbqt=pq, cx=round(float(c[0]), 3),
                         cy=round(float(c[1]), 3), cz=round(float(c[2]), 3),
                         box=26.0,
                         site_source='annotated site residues ' + ';'.join(tags),
                         note=note))

    # ---------- crystal controls ----------
    for pdb, ch, het, note in POS_XTAL:
        lines = read_pdb(os.path.join(raw_pdb, pdb + '.pdb'))
        c, n = het_centroid(lines, ch, het)
        clean = clean_crystal(lines, ch)
        name = f'{pdb}_{ch}'
        pdb_out = os.path.join(rec_dir, name + '.pdb')
        open(pdb_out, 'w').writelines(clean + ['END\n'])
        pq = os.path.join(rec_dir, name + '.pdbqt')
        to_pdbqt_receptor(pdb_out, pq)
        rows.append(dict(receptor=name, kind='pos_xtal', family='GUS',
                         receptor_pdbqt=pq, cx=round(float(c[0]), 3),
                         cy=round(float(c[1]), 3), cz=round(float(c[2]), 3),
                         box=26.0,
                         site_source=f'centroid of bound {het} ({n} heavy atoms)',
                         note=note))
        # keep the crystal ligand for the pose-recovery check
        if pdb == '3K4D':
            lp = os.path.join(lig_dir, f'{het}_xtal.pdb')
            open(lp, 'w').writelines(
                [l for l in lines if l[:6] == 'HETATM' and l[21] == ch
                 and l[17:20].strip() == het] + ['END\n'])
            subprocess.run(['obabel', lp, '-O',
                            os.path.join(lig_dir, f'{het}_xtal.pdbqt'),
                            '-p', '7.4'], check=True, capture_output=True)

    for pdb, ch, resl, note in NEG_XTAL:
        lines = read_pdb(os.path.join(raw_pdb, pdb + '.pdb'))
        clean = clean_crystal(lines, ch)
        cas = ca_table(clean)
        pts, tags = [], []
        for num, aa in resl:
            if num not in cas or cas[num][0] != aa:
                raise SystemExit(f'{pdb}: residue {num} is not {aa} '
                                 f'({cas.get(num)})')
            pts.append(cas[num][1])
            tags.append(f'{aa}{num}')
        c = np.mean(pts, 0)
        name = f'{pdb}_{ch}'
        pdb_out = os.path.join(rec_dir, name + '.pdb')
        open(pdb_out, 'w').writelines(clean + ['END\n'])
        pq = os.path.join(rec_dir, name + '.pdbqt')
        to_pdbqt_receptor(pdb_out, pq)
        rows.append(dict(receptor=name, kind='neg_xtal', family='unrelated',
                         receptor_pdbqt=pq, cx=round(float(c[0]), 3),
                         cy=round(float(c[1]), 3), cz=round(float(c[2]), 3),
                         box=26.0, site_source='catalytic residues ' + ';'.join(tags),
                         note=note))

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, 'boxes.tsv'), sep='\t', index=False)
    print()
    print(df[['receptor', 'kind', 'family', 'cx', 'cy', 'cz',
              'site_source']].to_string(index=False))


if __name__ == '__main__':
    main()
