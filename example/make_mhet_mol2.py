#!/usr/bin/env python3
"""Build an MHET (mono(2-hydroxyethyl) terephthalate, carboxylate form) mol2 with
unique atom names, SYBYL atom types and all hydrogens, as input for
abcg2_conformer_charges.py."""
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("-o", "--output", default="MHET.mol2", help="output mol2")
ap.add_argument("-rn", "--resname", default="MHT", help="residue name")
ap.add_argument("--seed", type=int, default=7, help="embedding seed")
args = ap.parse_args()

# heavy atoms in a fixed order so that names are meaningful
smi = "[O:1]([H])[CH2:2][CH2:3][O:4][C:5](=[O:6])[c:7]1[cH:8][cH:9][c:10]([cH:11][cH:12]1)[C:13](=[O:14])[O-:15]"
heavy_names = {1: "O1", 2: "C1", 3: "C2", 4: "O2", 5: "C3", 6: "O3", 7: "C4", 8: "C5",
               9: "C6", 10: "C7", 11: "C8", 12: "C9", 13: "C10", 14: "O4", 15: "O5"}
mol = Chem.AddHs(Chem.MolFromSmiles(smi))
AllChem.EmbedMolecule(mol, randomSeed=args.seed)
AllChem.MMFFOptimizeMolecule(mol)
# extended start: O1-C1-C2-O2 and C1-C2-O2-C3 anti
conf = mol.GetConformer()
idx = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
from rdkit.Chem import rdMolTransforms as T
T.SetDihedralDeg(conf, idx[1], idx[2], idx[3], idx[4], 180.0)
T.SetDihedralDeg(conf, idx[2], idx[3], idx[4], idx[5], 180.0)
AllChem.MMFFOptimizeMolecule(mol, maxIters=5000)

names, types = {}, {}
for a in mol.GetAtoms():
    m = a.GetAtomMapNum()
    if m:
        names[a.GetIdx()] = heavy_names[m]
for a in mol.GetAtoms():                       # hydrogens named after their heavy atom
    if a.GetAtomicNum() == 1:
        p = a.GetNeighbors()[0]
        hs = sorted(n.GetIdx() for n in p.GetNeighbors() if n.GetAtomicNum() == 1)
        stem = names[p.GetIdx()][1:]
        if p.GetSymbol() == "O":
            names[a.GetIdx()] = "HO" + stem
        else:
            names[a.GetIdx()] = "H" + stem + (str(hs.index(a.GetIdx()) + 1) if len(hs) > 1 else "")
for a in mol.GetAtoms():
    s = a.GetSymbol()
    if s == "H":
        t = "H"
    elif s == "C":
        t = "C.ar" if a.GetIsAromatic() else ("C.3" if a.GetHybridization() == Chem.HybridizationType.SP3 else "C.2")
    else:
        t = "O.co2" if any(n.GetAtomMapNum() == 13 for n in a.GetNeighbors()) else (
            "O.2" if any(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds()) else "O.3")
    types[a.GetIdx()] = t
assert len(set(names.values())) == mol.GetNumAtoms()

order = sorted(range(mol.GetNumAtoms()), key=lambda i: (mol.GetAtomWithIdx(i).GetAtomicNum() == 1, i))
num = {i: k + 1 for k, i in enumerate(order)}
X = conf.GetPositions()
L = ["@<TRIPOS>MOLECULE", "MHET", f"{mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0", "SMALL", "NO_CHARGES", "",
     "@<TRIPOS>ATOM"]
for i in order:
    L.append(f"{num[i]:7d} {names[i]:<6s} {X[i][0]:10.4f} {X[i][1]:10.4f} {X[i][2]:10.4f} "
             f"{types[i]:<6s} {1:4d} {args.resname:<6s} {0.0:10.6f}")
L.append("@<TRIPOS>BOND")
for k, b in enumerate(mol.GetBonds(), 1):
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    bt = "ar" if b.GetIsAromatic() else {1.0: "1", 2.0: "2"}[b.GetBondTypeAsDouble()]
    if types[i] == "O.co2" or types[j] == "O.co2":
        bt = "ar"
    L.append(f"{k:6d} {min(num[i], num[j]):5d} {max(num[i], num[j]):5d} {bt}")
L.append("@<TRIPOS>SUBSTRUCTURE")
L.append(f"     1 {args.resname:<6s}     1 RESIDUE           4 A     {args.resname}     0 ROOT")
open(args.output, "w").write("\n".join(L) + "\n")
print("wrote", args.output, mol.GetNumAtoms(), "atoms;", ", ".join(names[i] for i in order))
