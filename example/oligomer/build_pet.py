#!/usr/bin/env python3
"""
build_pet.py

Builds the PET molecules used on the oligomer page as mol2 files (SYBYL
types, all hydrogens, extended all-trans geometry), with every atom assigned
to a residue in the substructure columns of the mol2:

  PTM   repeat unit  -CH2-O-C(=O)-C6H4-C(=O)-O-CH2-   (cut at the CH2-CH2 bonds)
  PTB   first unit   CH3-CH2-O-C(=O)-C6H4-C(=O)-O-CH2-   (ethyl ester end)
  PTE   last unit    -CH2-O-C(=O)-C6H4-C(=O)-O-CH2-CH3   (ethyl ester end)
  CAP   capping atoms of a fragment, removed when the residues are assembled

Fragments (the residue of interest with caps that stand in for its
neighbours) and the whole oligomer:

  python3 build_pet.py --fragment mid    -o PET_mid.mol2     # CAP-PTM-CAP, acetoxyethyl caps
  python3 build_pet.py --fragment end    -o PET_end.mol2     # PTB-CAP
  python3 build_pet.py --fragment endrev -o PET_endrev.mol2  # same molecule, named CAP-PTE
  python3 build_pet.py --fragment mid --caps methyl -o PET_mid_me.mol2   # CH3-PTM-CH3
  python3 build_pet.py --oligomer 5      -o PET5.mol2        # PTB-PTM-PTM-PTM-PTE
  python3 build_pet.py --oligomer 5 --occo 70 -o PET5_gauche.mol2   # O-C-C-O gauche

The acetoxyethyl caps (-CH2-O-C(=O)-CH3) repeat the first four atoms of the
neighbouring unit (CH2, ester O, carbonyl C and O) and replace its ring by a
methyl group; the methyl caps keep only the neighbouring CH2.
"""

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
try:
    import abcg2_conformer_charges as base
except ImportError:
    sys.exit("ERROR: abcg2_conformer_charges.py not found two directories up.")

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDistGeom

S, D, A = Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.AROMATIC

# repeat unit: (name, element), bonds; C1 is the head CH2, C10 the tail CH2
UNIT_ATOMS = [("C1", "C"), ("H11", "H"), ("H12", "H"), ("O1", "O"), ("C2", "C"),
              ("O2", "O"), ("C3", "C"), ("C4", "C"), ("H4", "H"), ("C5", "C"),
              ("H5", "H"), ("C6", "C"), ("C7", "C"), ("H7", "H"), ("C8", "C"),
              ("H8", "H"), ("C9", "C"), ("O4", "O"), ("O3", "O"), ("C10", "C"),
              ("H101", "H"), ("H102", "H")]
UNIT_BONDS = [("C1", "H11", S), ("C1", "H12", S), ("C1", "O1", S), ("O1", "C2", S),
              ("C2", "O2", D), ("C2", "C3", S), ("C3", "C4", A), ("C4", "H4", S),
              ("C4", "C5", A), ("C5", "H5", S), ("C5", "C6", A), ("C6", "C7", A),
              ("C7", "H7", S), ("C7", "C8", A), ("C8", "H8", S), ("C8", "C3", A),
              ("C6", "C9", S), ("C9", "O4", D), ("C9", "O3", S), ("O3", "C10", S),
              ("C10", "H101", S), ("C10", "H102", S)]
METHYL = [("CT", "C"), ("HT1", "H"), ("HT2", "H"), ("HT3", "H")]
METHYL_BONDS = [("CT", "HT1", S), ("CT", "HT2", S), ("CT", "HT3", S)]


def cap_atoms(side, kind):
    """Cap on the head side (P, bonded to C1) or tail side (N, bonded to C10)."""
    p = "P" if side == "head" else "N"
    if kind == "methyl":
        atoms = [(f"C{p}1", "C"), (f"H{p}1", "H"), (f"H{p}2", "H"), (f"H{p}3", "H")]
        bonds = [(f"C{p}1", f"H{p}{k}", S) for k in (1, 2, 3)]
        return atoms, bonds, f"C{p}1"
    atoms = [(f"C{p}1", "C"), (f"H{p}1", "H"), (f"H{p}2", "H"), (f"O{p}1", "O"),
             (f"C{p}2", "C"), (f"O{p}2", "O"), (f"C{p}3", "C"), (f"H{p}3", "H"),
             (f"H{p}4", "H"), (f"H{p}5", "H")]
    bonds = [(f"C{p}1", f"H{p}1", S), (f"C{p}1", f"H{p}2", S), (f"C{p}1", f"O{p}1", S),
             (f"O{p}1", f"C{p}2", S), (f"C{p}2", f"O{p}2", D), (f"C{p}2", f"C{p}3", S),
             (f"C{p}3", f"H{p}3", S), (f"C{p}3", f"H{p}4", S), (f"C{p}3", f"H{p}5", S)]
    return atoms, bonds, f"C{p}1"


class Builder:
    def __init__(self):
        self.mol = Chem.RWMol()
        self.info = []          # (name, resname, resid)

    def add(self, atoms, bonds, resname, resid):
        idx = {}
        for name, el in atoms:
            a = Chem.Atom(el)
            a.SetNoImplicit(True)
            idx[name] = self.mol.AddAtom(a)
            self.info.append((name, resname, resid))
        for i, j, t in bonds:
            self.mol.AddBond(idx[i], idx[j], t)
            if t == A:
                self.mol.GetBondBetweenAtoms(idx[i], idx[j]).SetIsAromatic(True)
                self.mol.GetAtomWithIdx(idx[i]).SetIsAromatic(True)
                self.mol.GetAtomWithIdx(idx[j]).SetIsAromatic(True)
        return idx

    def bond(self, i, j):
        self.mol.AddBond(i, j, S)


def unit(b, resname, resid, methyl_on=None):
    idx = b.add(UNIT_ATOMS, UNIT_BONDS, resname, resid)
    if methyl_on:
        m = b.add(METHYL, METHYL_BONDS, resname, resid)
        idx.update(m)
        b.bond(idx[methyl_on], m["CT"])
    return idx


def build(kind, caps, n):
    b = Builder()
    if kind == "oligomer":
        prev = None
        for k in range(1, n + 1):
            if k == 1:
                idx = unit(b, "PTB", k, methyl_on="C1")
            elif k == n:
                idx = unit(b, "PTE", k, methyl_on="C10")
            else:
                idx = unit(b, "PTM", k)
            if prev is not None:
                b.bond(prev["C10"], idx["C1"])
            prev = idx
        return b
    if kind == "mid":
        h_atoms, h_bonds, h_link = cap_atoms("head", caps)
        hc = b.add(h_atoms, h_bonds, "CAP", 1)
        idx = unit(b, "PTM", 2)
        t_atoms, t_bonds, t_link = cap_atoms("tail", caps)
        tc = b.add(t_atoms, t_bonds, "CAP", 3)
        b.bond(hc[h_link], idx["C1"])
        b.bond(idx["C10"], tc[t_link])
    elif kind == "end":
        idx = unit(b, "PTB", 1, methyl_on="C1")
        t_atoms, t_bonds, t_link = cap_atoms("tail", caps)
        tc = b.add(t_atoms, t_bonds, "CAP", 2)
        b.bond(idx["C10"], tc[t_link])
    elif kind == "endrev":
        h_atoms, h_bonds, h_link = cap_atoms("head", caps)
        hc = b.add(h_atoms, h_bonds, "CAP", 1)
        idx = unit(b, "PTE", 2, methyl_on="C10")
        b.bond(hc[h_link], idx["C1"])
    return b


def sybyl(atom):
    el = atom.GetSymbol()
    if el == "H":
        return "H"
    if el == "C":
        if atom.GetIsAromatic():
            return "C.ar"
        return "C.2" if any(b.GetBondType() == D for b in atom.GetBonds()) else "C.3"
    if el == "O":
        return "O.2" if atom.GetDegree() == 1 else "O.3"
    return el


def write_mol2(path, mol, X, info, title):
    bond_code = {S: "1", D: "2", A: "ar"}
    res = []
    for name, rn, ri in info:
        if (ri, rn) not in res:
            res.append((ri, rn))
    L = ["@<TRIPOS>MOLECULE", title,
         f"{mol.GetNumAtoms():5d} {mol.GetNumBonds():5d} {len(res):5d}     0     0",
         "SMALL", "NO_CHARGES", "", "", "@<TRIPOS>ATOM"]
    for a in mol.GetAtoms():
        i = a.GetIdx()
        name, rn, ri = info[i]
        x = X[i]
        L.append(f"{i + 1:7d} {name:<8s} {x[0]:10.4f} {x[1]:10.4f} {x[2]:10.4f} "
                 f"{sybyl(a):<6s} {ri:4d} {rn:<8s} {0.0:10.6f}")
    L.append("@<TRIPOS>BOND")
    for k, bd in enumerate(mol.GetBonds(), 1):
        L.append(f"{k:6d} {bd.GetBeginAtomIdx() + 1:5d} {bd.GetEndAtomIdx() + 1:5d} "
                 f"{bond_code[bd.GetBondType()]}")
    L.append("@<TRIPOS>SUBSTRUCTURE")
    for ri, rn in res:
        first = next(i for i, (_, r, k) in enumerate(info) if k == ri) + 1
        L.append(f"{ri:6d} {rn:<8s} {first:5d} RESIDUE           4 A     {rn:<4s}  1 ROOT")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


def set_occo(mol, X, deg):
    """O-C-C-O torsions set to deg and held there in an MMFF94s minimisation."""
    from rdkit.Chem import rdMolTransforms
    from rdkit.Geometry import Point3D
    m = Chem.Mol(mol)
    conf = m.GetConformer()
    for i, x in enumerate(X):
        conf.SetAtomPosition(i, Point3D(*x))
    tors = m.GetSubstructMatches(Chem.MolFromSmarts("[OX2]-[CH2]-[CH2]-[OX2]"))
    for t in tors:
        rdMolTransforms.SetDihedralDeg(conf, *t, deg)
    props = AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94s")
    props.SetMMFFDielectricConstant(4.0)
    ff = AllChem.MMFFGetMoleculeForceField(m, props)
    for t in tors:
        ff.MMFFAddTorsionConstraint(*t, False, deg, deg, 500.0)
    ff.Minimize(maxIts=5000)
    return m.GetConformer().GetPositions()


def main():
    p = argparse.ArgumentParser(description="Build PET fragments and oligomers as "
                                            "mol2 files with residue assignments.",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--fragment", choices=["mid", "end", "endrev"],
                   help="capped fragment: repeat unit (mid), first unit (end), or the "
                        "end fragment named as the last unit (endrev)")
    g.add_argument("--oligomer", type=int, metavar="N",
                   help="oligomer of N >= 2 units with ethyl ester ends")
    p.add_argument("--caps", choices=["acetoxyethyl", "methyl"], default="acetoxyethyl",
                   help="caps of the fragments")
    p.add_argument("-o", "--out", required=True, help="output mol2")
    p.add_argument("--seed", type=int, default=2026, help="embedding seed")
    p.add_argument("--occo", type=float, default=None, metavar="DEG",
                   help="set every O-CH2-CH2-O torsion to this value (e.g. 70 for "
                        "gauche) instead of 180, the other torsions staying extended")
    a = p.parse_args()
    if a.oligomer is not None and a.oligomer < 2:
        p.error("--oligomer needs N >= 2")
    RDLogger.DisableLog("rdApp.*")

    b = build("oligomer" if a.oligomer else a.fragment, a.caps, a.oligomer or 0)
    mol = b.mol.GetMol()
    Chem.SanitizeMol(mol)
    ep = rdDistGeom.ETKDGv3()
    ep.randomSeed = a.seed
    ep.useRandomCoords = True
    if rdDistGeom.EmbedMolecule(mol, ep) != 0:
        sys.exit("ERROR: embedding failed")
    # extended geometry: every acyclic heavy-atom torsion at 180 deg, MMFF94s
    opts = argparse.Namespace(ff="mmff94s", dielectric=4.0, dist_dielectric=False,
                              max_iters=5000)
    X = base.idealize_all_trans(mol, opts)
    if a.occo is not None:
        X = set_occo(mol, X, a.occo)
    title = os.path.splitext(os.path.basename(a.out))[0]
    write_mol2(a.out, mol, X, b.info, title)
    res = sorted({(ri, rn) for _, rn, ri in b.info})
    print(f"wrote {a.out}: {mol.GetNumAtoms()} atoms, residues "
          f"{' '.join(rn for _, rn in res)}, {Chem.MolToSmiles(Chem.RemoveHs(mol))}")


if __name__ == "__main__":
    main()
