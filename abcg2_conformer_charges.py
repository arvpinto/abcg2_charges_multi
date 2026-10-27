#!/usr/bin/env python3
"""
abcg2_conformer_charges.py

Conformer-averaged ABCG2 (or AM1-BCC) partial charges for a flexible small
molecule, added conformer by conformer until the average is statistically
converged. Only a mol2 file is needed:

    python3 abcg2_conformer_charges.py molecule.mol2

Protocol
  1. Read the mol2 (atom names, order and residue name are kept), perceive
     bond orders, formal charges and stereochemistry from the connectivity
     and the net charge (taken from the mol2 charges, or 0, or -nc).
  2. Conformer pool: ETKDGv3 embedding, MMFF94s minimisation with a
     screened dielectric (eps = 4), removal of stereo-inverted conformers.
  3. Filters (applied to the pool and again after the AM1 minimisation):
       - no intramolecular H-bond (H...A <= 2.5 A, D-H...A >= 120 deg,
         donor/acceptor >= 4 bonds apart);
       - no contact between formally charged atoms (<= 4 A);
       - N/O atoms >= 5 bonds apart stay >= 5 A apart (no dipolar stacking
         of distant amides, esters, ...);
       - every secondary amide trans;
       - MMFF energy within 15 kcal/mol of the lowest conformer.
     Near-duplicates (symmetry- and mirror-aware heavy-atom RMSD < 0.5 A)
     are removed and the pool is shuffled with a fixed seed, so conformers
     are used in an unbiased random order.
  4. antechamber -c abcg2 -at gaff2 (standard sqm AM1 minimisation) on the
     conformers, several in parallel; results are consumed in draw order.
  5. After every accepted conformer the per-atom running mean, SD and SEM
     are updated. The run stops when N >= --min-confs and
       max_atoms SEM                           <= --sem-tol
       max_atoms |mean(N) - mean(N - window)|  <= --drift-tol
     or when --max-confs conformers have been accepted.
  6. A single reference conformer is also charged: the input if it passes
     all filters and is fully extended, otherwise an idealised all-trans
     structure built from it. Its AM1-relaxed geometry carries the final
     charges.
  7. Final charges: ensemble mean, averaged over symmetry-equivalent atoms,
     exact net charge, rounded to 1e-6 e. Written as mol2 (GAFF2 types),
     with a parmchk2 frcmod and, if tleap is available, a checked .lib.

Outputs go to <input basename>_<method>/ (change with -o); the CSV/JSON
files there are read by abcg2_analysis.py.

Requirements: AmberTools with ABCG2 support (antechamber, sqm, parmchk2;
tleap optional), Python 3.8+, RDKit, NumPy.
"""

import argparse
import concurrent.futures as cf
import csv
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import time

import numpy as np

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, rdDetermineBonds, rdDistGeom, rdMolTransforms
    from rdkit.Geometry import Point3D
except ImportError:
    sys.exit("ERROR: RDKit is required (conda install -c conda-forge rdkit).")

PT = Chem.GetPeriodicTable()
AMIDE_SMARTS = Chem.MolFromSmarts("[OX1]=[CX3]-[NX3;H1]-[#1]")
GAUCHE_SMARTS = Chem.MolFromSmarts("[!#1]-[CH2]-[CH2]-[!#1]")
GAFF_PREFIX = [("cl", "Cl"), ("br", "Br"), ("c", "C"), ("n", "N"), ("o", "O"),
               ("s", "S"), ("p", "P"), ("h", "H"), ("f", "F"), ("i", "I")]
MIRROR = np.array([-1.0, 1.0, 1.0])


# ==========================================================================
# mol2 I/O
# ==========================================================================
class Mol2:
    """Minimal Tripos mol2 reader/writer that keeps the original lines."""

    def __init__(self, path):
        self.path = path
        with open(path) as fh:
            self.lines = fh.read().splitlines()
        self.atoms = []
        self.bonds = []
        self.mol_name = "MOL"
        self._parse()

    def _parse(self):
        section = None
        id2idx = {}
        n_mol = 0
        for ln, line in enumerate(self.lines):
            s = line.strip()
            if s.startswith("@<TRIPOS>"):
                section = s[9:].upper()
                if section == "MOLECULE":
                    n_mol += 1
                    if n_mol > 1:
                        print("WARNING: more than one MOLECULE block; only the "
                              "first is used.", file=sys.stderr)
                        self.lines = self.lines[:ln]
                        break
                    if ln + 1 < len(self.lines):
                        self.mol_name = self.lines[ln + 1].strip() or "MOL"
                continue
            if not s or s.startswith("#"):
                continue
            if section == "ATOM":
                f = s.split()
                if len(f) < 6:
                    raise ValueError(f"malformed ATOM line {ln + 1}: {line}")
                id2idx[int(f[0])] = len(self.atoms)
                self.atoms.append({
                    "line": ln, "id": int(f[0]), "name": f[1],
                    "xyz": np.array([float(f[2]), float(f[3]), float(f[4])]),
                    "type": f[5],
                    "subst_id": f[6] if len(f) > 6 else "1",
                    "subst_name": f[7] if len(f) > 7 else "MOL",
                    "charge": float(f[8]) if len(f) > 8 else None,
                    "extra": f[9:] if len(f) > 9 else []})
            elif section == "BOND":
                f = s.split()
                if len(f) < 4:
                    raise ValueError(f"malformed BOND line {ln + 1}: {line}")
                self.bonds.append((id2idx[int(f[1])], id2idx[int(f[2])], f[3]))
        if not self.atoms:
            raise ValueError(f"no atoms in {self.path}")
        if not self.bonds and len(self.atoms) > 1:
            raise ValueError(f"no BOND section in {self.path}")

    @property
    def resname(self):
        return self.atoms[0]["subst_name"]

    @property
    def names(self):
        return [a["name"] for a in self.atoms]

    def net_charge_from_file(self):
        q = [a["charge"] for a in self.atoms]
        if any(c is None for c in q) or all(abs(c) < 1e-8 for c in q):
            return None
        return int(round(sum(q)))

    def has_charges(self):
        return self.net_charge_from_file() is not None or any(
            a["charge"] not in (None, 0.0) for a in self.atoms)

    def write(self, path, coords=None, charges=None, title=None):
        out = list(self.lines)
        for i, a in enumerate(self.atoms):
            xyz = coords[i] if coords is not None else a["xyz"]
            q = charges[i] if charges is not None else (
                a["charge"] if a["charge"] is not None else 0.0)
            line = (f"{a['id']:7d} {a['name']:<8s} {xyz[0]:10.4f} "
                    f"{xyz[1]:10.4f} {xyz[2]:10.4f} {a['type']:<6s} "
                    f"{a['subst_id']:>4s} {a['subst_name']:<8s} {q:10.6f}")
            if a["extra"]:
                line += " " + " ".join(a["extra"])
            out[a["line"]] = line
        if title is not None:
            for ln, line in enumerate(out):
                if line.strip().upper() == "@<TRIPOS>MOLECULE" and ln + 1 < len(out):
                    out[ln + 1] = title
                    break
        with open(path, "w") as fh:
            fh.write("\n".join(out) + "\n")


# ==========================================================================
# RDKit molecule from the mol2
# ==========================================================================
def _valid_element(sym):
    try:
        return PT.GetAtomicNumber(sym) > 0
    except Exception:
        return False


def element_of(atom, mode):
    t, name = atom["type"], atom["name"]

    def from_type():
        base = t.split(".")[0]
        if base[:1].isupper() and _valid_element(base.capitalize()):
            return base.capitalize()            # SYBYL types
        tl = t.lower()
        for pre, el in GAFF_PREFIX:            # GAFF/GAFF2 types
            if tl.startswith(pre):
                return el
        return None                            # e.g. numeric AM1-BCC types

    def from_name():
        letters = "".join(ch for ch in name if ch.isalpha())
        if len(letters) >= 2 and letters[1].islower() and _valid_element(letters[:2]):
            return letters[:2]
        if letters and _valid_element(letters[0].upper()):
            return letters[0].upper()
        return None

    el = {"type": from_type, "name": from_name}.get(mode, lambda: None)()
    if mode == "auto":
        el = from_type() or from_name()
    if el is None:
        raise ValueError(f"cannot determine the element of atom {name} "
                         f"(type {t}); try --element-from name or type")
    return el


def build_rdkit_mol(m2, net_charge, element_mode):
    rw = Chem.RWMol()
    for a in m2.atoms:
        at = Chem.Atom(element_of(a, element_mode))
        at.SetNoImplicit(True)
        rw.AddAtom(at)
    for i, j, _ in m2.bonds:
        rw.AddBond(i, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, a in enumerate(m2.atoms):
        conf.SetAtomPosition(i, Point3D(*a["xyz"]))
    conf.Set3D(True)
    mol.AddConformer(conf, assignId=True)
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()) and \
            any(a.GetAtomicNum() == 6 for a in mol.GetAtoms()):
        raise ValueError("the mol2 has no hydrogens; add them with the "
                         "intended protonation state first")
    try:
        rdDetermineBonds.DetermineBondOrders(mol, charge=net_charge)
    except Exception as exc:
        raise RuntimeError(f"bond-order perception failed for net charge "
                           f"{net_charge} ({exc}); check the protonation "
                           f"state and -nc")
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    if Chem.GetFormalCharge(mol) != net_charge:
        raise RuntimeError("perceived formal charge differs from the net charge")
    return mol


# ==========================================================================
# geometry helpers
# ==========================================================================
def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def kabsch_rmsd(P, Q):
    Pc, Qc = P - P.mean(axis=0), Q - Q.mean(axis=0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    diff = (R @ Pc.T).T - Qc
    return math.sqrt((diff ** 2).sum() / len(P))


def kabsch_transform(P_fit, Q_fit, P_all):
    pc, qc = P_fit.mean(axis=0), Q_fit.mean(axis=0)
    U, S, Vt = np.linalg.svd((P_fit - pc).T @ (Q_fit - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return (R @ (P_all - pc).T).T + qc


def heavy_atom_permutations(mol, heavy, max_matches):
    hmol = Chem.RemoveHs(mol)
    if hmol.GetNumAtoms() != len(heavy):
        return [np.arange(len(heavy))]
    matches = hmol.GetSubstructMatches(hmol, uniquify=False, useChirality=False,
                                       maxMatches=max_matches)
    return [np.array(m) for m in matches] or [np.arange(len(heavy))]


def sym_rmsd(P, Q, perms, allow_mirror):
    best = min(kabsch_rmsd(P[p], Q) for p in perms)
    if allow_mirror:
        Pm = P * MIRROR
        best = min(best, min(kabsch_rmsd(Pm[p], Q) for p in perms))
    return best


def gauche_count(mol, X):
    seen, n = set(), 0
    for t in mol.GetSubstructMatches(GAUCHE_SMARTS):
        key = tuple(sorted((t[1], t[2])))
        if key in seen:
            continue
        seen.add(key)
        if abs(dihedral(*(X[i] for i in t))) < 120.0:
            n += 1
    return n


def hbond_sites(mol):
    donors_h, acceptors = [], []
    for at in mol.GetAtoms():
        z = at.GetAtomicNum()
        if z == 1:
            nb = at.GetNeighbors()
            if nb and nb[0].GetAtomicNum() in (7, 8, 16):
                donors_h.append((at.GetIdx(), nb[0].GetIdx()))
        elif z in (7, 8):
            if z == 7 and at.GetFormalCharge() > 0 and at.GetTotalDegree() >= 4:
                continue
            if z == 7 and at.GetIsAromatic() and at.GetTotalNumHs() > 0:
                continue
            acceptors.append(at.GetIdx())
    return donors_h, acceptors


class Checker:
    """Intramolecular-contact and amide checks for one molecule."""

    def __init__(self, mol, args):
        self.args = args
        self.heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        self.topo = Chem.GetDistanceMatrix(mol)
        self.donors_h, self.acceptors = hbond_sites(mol)
        self.charged = [a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge()]
        polar = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in (7, 8)]
        self.polar_pairs = [(i, j) for k, i in enumerate(polar) for j in polar[k + 1:]
                            if self.topo[i, j] >= args.polar_topo]
        self.amides = list(mol.GetSubstructMatches(AMIDE_SMARTS))
        self.amide_ref = self.amide_states(mol.GetConformer().GetPositions())

    def amide_states(self, X):
        return ["trans" if abs(dihedral(X[o], X[c], X[n], X[h])) > 90 else "cis"
                for o, c, n, h in self.amides]

    def amide_report(self, X, names):
        out = []
        for o, c, n, h in self.amides:
            w = dihedral(X[o], X[c], X[n], X[h])
            out.append(f"{names[o]}={names[c]}-{names[n]}-{names[h]} {w:7.1f} deg "
                       f"({'trans' if abs(w) > 90 else 'cis'})")
        return out

    def evaluate(self, X):
        a = self.args
        Xh = X[self.heavy]
        rg = math.sqrt(((Xh - Xh.mean(axis=0)) ** 2).sum(axis=1).mean())
        hb = 0
        for h, d in self.donors_h:
            for acc in self.acceptors:
                if acc == d or self.topo[d, acc] < a.min_topo:
                    continue
                if np.linalg.norm(X[h] - X[acc]) > a.hb_dist:
                    continue
                v1, v2 = X[d] - X[h], X[acc] - X[h]
                cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                if math.degrees(math.acos(max(-1.0, min(1.0, cosang)))) >= a.hb_angle:
                    hb += 1
        ion = sum(1 for k, i in enumerate(self.charged) for j in self.charged[k + 1:]
                  if self.topo[i, j] >= a.min_topo
                  and np.linalg.norm(X[i] - X[j]) <= a.ion_dist)
        dpol = min((float(np.linalg.norm(X[i] - X[j])) for i, j in self.polar_pairs),
                   default=float("inf"))
        states = self.amide_states(X)
        amide_ok = (states == self.amide_ref) if a.allow_cis_amide else \
            all(s == "trans" for s in states)
        reasons = []
        if hb:
            reasons.append(f"{hb} H-bond")
        if ion:
            reasons.append(f"{ion} ionic contact")
        if dpol < a.polar_dist:
            reasons.append(f"polar contact {dpol:.2f} A")
        if not amide_ok:
            reasons.append("cis amide" if not a.allow_cis_amide else
                           "amide cis/trans changed")
        return {"rg": rg, "n_hbond": hb, "n_ionic": ion, "min_polar_dist": dpol,
                "amide_ok": amide_ok, "ok": not reasons, "reason": "; ".join(reasons)}


def mmff_forcefield(m, args, conf_id=-1):
    variant = "MMFF94" if args.ff == "mmff94" else "MMFF94s"
    props = None if args.ff == "uff" else \
        AllChem.MMFFGetMoleculeProperties(m, mmffVariant=variant)
    if props is None:
        return AllChem.UFFGetMoleculeForceField(m, confId=conf_id), "uff"
    props.SetMMFFDielectricConstant(args.dielectric)
    props.SetMMFFDielectricModel(2 if args.dist_dielectric else 1)
    return AllChem.MMFFGetMoleculeForceField(m, props, confId=conf_id), args.ff


def idealize_all_trans(mol, args):
    """All acyclic heavy-atom torsions set to 180 deg, then minimised."""
    m = Chem.Mol(mol)
    conf = m.GetConformer()
    hdeg = {a.GetIdx(): sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1)
            for a in m.GetAtoms()}

    def pick(a, exclude):
        c = [nb.GetIdx() for nb in m.GetAtomWithIdx(a).GetNeighbors()
             if nb.GetIdx() != exclude and nb.GetAtomicNum() > 1]
        return max(c, key=lambda k: (hdeg[k], -k)) if c else None

    for b in m.GetBonds():
        if b.IsInRing() or b.GetBondType() != Chem.BondType.SINGLE:
            continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if 1 in (m.GetAtomWithIdx(i).GetAtomicNum(), m.GetAtomWithIdx(j).GetAtomicNum()):
            continue
        a, d = pick(i, j), pick(j, i)
        if a is not None and d is not None:
            rdMolTransforms.SetDihedralDeg(conf, a, i, j, d, 180.0)
    ff, _ = mmff_forcefield(m, args)
    ff.Minimize(maxIts=max(args.max_iters, 5000))
    return m.GetConformer().GetPositions()


# ==========================================================================
# conformer pool
# ==========================================================================
def build_pool(mol, checker, args, say):
    molc = Chem.Mol(mol)
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = args.seed
    p.numThreads = args.threads
    p.pruneRmsThresh = args.prune_rms
    p.useRandomCoords = True
    p.enforceChirality = True
    cids = list(rdDistGeom.EmbedMultipleConfs(molc, numConfs=args.nconf, params=p))
    if not cids:
        raise RuntimeError("ETKDG could not embed any conformer")
    energies = {}
    for cid in cids:
        ff, used = mmff_forcefield(molc, args, cid)
        ff.Minimize(maxIts=args.max_iters)
        energies[cid] = ff.CalcEnergy()
    if used != args.ff:
        say("WARNING: MMFF parameters missing, UFF used for the pool.")
    emin = min(energies.values())
    ref_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))

    rows = []
    for cid in cids:
        X = molc.GetConformer(cid).GetPositions()
        m = Chem.Mol(molc)
        m.RemoveAllConformers()
        m.AddConformer(Chem.Conformer(molc.GetConformer(cid)), assignId=True)
        Chem.AssignStereochemistryFrom3D(m)
        stereo_ok = Chem.MolToSmiles(Chem.RemoveHs(m)) == ref_smiles
        row = {"pool_id": cid, "mmff_dE": energies[cid] - emin,
               "stereo_ok": stereo_ok, "n_gauche": gauche_count(mol, X), "X": X}
        row.update(checker.evaluate(X))
        extra = []
        if not stereo_ok:
            extra.append("stereo inverted")
        if row["mmff_dE"] > args.ewin:
            extra.append("outside energy window")
        if extra:
            row["ok"] = False
            row["reason"] = "; ".join([r for r in [row["reason"]] if r] + extra)
        rows.append(row)

    heavy = checker.heavy
    perms = heavy_atom_permutations(mol, heavy, args.max_sym)
    achiral = not Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                            useLegacyImplementation=False)
    kept = []
    for r in sorted((r for r in rows if r["ok"]), key=lambda r: r["mmff_dE"]):
        if all(sym_rmsd(r["X"][heavy], k["X"][heavy], perms, achiral) >= args.div_rmsd
               for k in kept):
            kept.append(r)
        else:
            r["ok"], r["reason"] = False, "duplicate"
    random.Random(args.seed).shuffle(kept)
    return rows, kept


# ==========================================================================
# antechamber
# ==========================================================================
def read_sqm_pdb(path, names):
    xyz = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                xyz[line[12:16].strip()] = [float(line[30:38]), float(line[38:46]),
                                            float(line[46:54])]
    if set(xyz) != set(names):
        raise ValueError(f"atom names in {path} do not match the input")
    return np.array([xyz[n] for n in names])


def read_heat_of_formation(path):
    hof = None
    if os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                m = re.search(r"Heat of formation\s*=\s*(-?\d+\.\d+)\s*kcal/mol", line)
                if m:
                    hof = float(m.group(1))
    return hof


def run_antechamber(tag, mol2_path, wdir, net_charge, resname, args):
    os.makedirs(wdir, exist_ok=True)
    shutil.copy(mol2_path, os.path.join(wdir, "input.mol2"))
    outf = f"{tag}_{args.method}.mol2"
    cmd = [args.ac_exe, "-i", "input.mol2", "-fi", "mol2", "-o", outf, "-fo", "mol2",
           "-c", args.method, "-at", args.atomtype, "-nc", str(net_charge),
           "-rn", resname, "-pf", "n"]
    if args.ac_ek:
        cmd += ["-ek", args.ac_ek]
    t0 = time.time()
    with open(os.path.join(wdir, "antechamber.log"), "w") as log:
        rc = subprocess.run(cmd, cwd=wdir, stdout=log, stderr=subprocess.STDOUT).returncode
    dt = time.time() - t0
    out_path = os.path.join(wdir, outf)
    pdb = os.path.join(wdir, "sqm.pdb")
    if rc != 0 or not os.path.isfile(out_path) or not os.path.isfile(pdb):
        return {"tag": tag, "ok": False, "error": "antechamber/sqm failed", "time_s": dt}
    out = Mol2(out_path)
    res = {"tag": tag, "ok": True, "time_s": dt, "out_mol2": out_path, "sqm_pdb": pdb,
           "charges": {a["name"]: a["charge"] for a in out.atoms},
           "types": {a["name"]: a["type"] for a in out.atoms},
           "hof": read_heat_of_formation(os.path.join(wdir, "sqm.out"))}
    if not args.keep_scratch:
        for f in os.listdir(wdir):
            if f.startswith("ANTECHAMBER") or f in ("ATOMTYPE.INF", "sqm.in"):
                os.remove(os.path.join(wdir, f))
    return res


def check_amber(args):
    exes = {}
    for key, exe in (("ac_exe", args.ac_exe), ("parmchk2_exe", args.parmchk2_exe)):
        path = exe if os.path.isfile(exe) else shutil.which(exe)
        if path is None:
            sys.exit(f"ERROR: {exe} not found; load AmberTools (source amber.sh) "
                     f"or pass its path.")
        exes[key] = os.path.abspath(path)
    help_txt = subprocess.run([exes["ac_exe"], "-h"], capture_output=True,
                              text=True).stdout
    if args.method not in help_txt:
        sys.exit(f"ERROR: this antechamber does not list charge method "
                 f"'{args.method}'; ABCG2 needs a recent AmberTools.")
    tleap = None
    if not args.no_tleap:
        tleap = args.tleap_exe if os.path.isfile(args.tleap_exe) else shutil.which(args.tleap_exe)
    return exes["ac_exe"], exes["parmchk2_exe"], tleap


# ==========================================================================
# statistics
# ==========================================================================
def stats(Q):
    n = len(Q)
    mean = Q.mean(axis=0)
    if n < 2:
        return mean, np.zeros(Q.shape[1]), np.full(Q.shape[1], np.inf)
    std = Q.std(axis=0, ddof=1)
    return mean, std, std / math.sqrt(n)


def bootstrap_ci(Q, nboot, seed, level=0.95):
    rng = np.random.default_rng(seed)
    n = len(Q)
    means = np.array([Q[rng.integers(0, n, n)].mean(axis=0) for _ in range(nboot)])
    return (np.quantile(means, (1 - level) / 2, axis=0),
            np.quantile(means, 1 - (1 - level) / 2, axis=0))


def finalize_charges(mean, net_charge, classes):
    q = mean.copy()
    for cls in set(classes):
        idx = [i for i, c in enumerate(classes) if c == cls]
        q[idx] = q[idx].mean()
    q += (net_charge - q.sum()) / len(q)
    q = np.round(q, 6)
    resid = round(net_charge - q.sum(), 6)
    if abs(resid) > 0:
        counts = {c: classes.count(c) for c in set(classes)}
        unique = [i for i, c in enumerate(classes) if counts[c] == 1]
        k = unique[0] if unique else int(np.argmax(np.abs(q)))
        q[k] = round(q[k] + resid, 6)
    return q


# ==========================================================================
# command line
# ==========================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Conformer-averaged ABCG2 charges, added until statistically "
                    "converged. Usage: abcg2_conformer_charges.py molecule.mol2",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("mol2", help="input mol2 (one molecule, all H present)")

    io = p.add_argument_group("input/output")
    io.add_argument("-o", "--outdir", default=None,
                    help="output directory (None = <input basename>_<method>)")
    io.add_argument("--prefix", default=None,
                    help="file prefix (None = input basename)")
    io.add_argument("-nc", "--net-charge", type=int, default=None,
                    help="net charge (None = rounded sum of the mol2 charges, or 0)")
    io.add_argument("-rn", "--resname", default=None,
                    help="residue name (None = from the mol2)")
    io.add_argument("--element-from", choices=["auto", "type", "name"],
                    default="auto", help="how elements are read from the mol2")
    io.add_argument("--keep-scratch", action="store_true",
                    help="keep antechamber scratch files")
    io.add_argument("--overwrite", action="store_true",
                    help="allow writing into an existing output directory")

    g = p.add_argument_group("conformer pool")
    g.add_argument("-n", "--nconf", type=int, default=500,
                   help="conformers embedded with ETKDGv3")
    g.add_argument("--seed", type=int, default=2026, help="random seed")
    g.add_argument("--threads", type=int, default=0,
                   help="threads for embedding (0 = all)")
    g.add_argument("--prune-rms", type=float, default=0.3,
                   help="ETKDG pruning RMSD (A)")
    g.add_argument("--ff", choices=["mmff94s", "mmff94", "uff"], default="mmff94s",
                   help="force field for pool minimisation")
    g.add_argument("--dielectric", type=float, default=4.0,
                   help="MMFF dielectric constant")
    g.add_argument("--dist-dielectric", action="store_true",
                   help="distance-dependent MMFF dielectric")
    g.add_argument("--max-iters", type=int, default=2000,
                   help="MMFF minimisation iterations")
    g.add_argument("--ewin", type=float, default=15.0,
                   help="MMFF energy window (kcal/mol)")
    g.add_argument("--div-rmsd", type=float, default=0.5,
                   help="duplicate threshold, heavy-atom RMSD (A)")
    g.add_argument("--max-sym", type=int, default=1000,
                   help="max automorphisms in the symmetry-aware RMSD")

    f = p.add_argument_group("contact filters (pool and after sqm)")
    f.add_argument("--hb-dist", type=float, default=2.5, help="H...A cutoff (A)")
    f.add_argument("--hb-angle", type=float, default=120.0,
                   help="minimum D-H...A angle (deg)")
    f.add_argument("--min-topo", type=int, default=4,
                   help="min bonds between H-bond or ionic partners")
    f.add_argument("--ion-dist", type=float, default=4.0,
                   help="contact cutoff between charged atoms (A)")
    f.add_argument("--polar-dist", type=float, default=5.0,
                   help="minimum distance between distant N/O atoms (A)")
    f.add_argument("--polar-topo", type=int, default=5,
                   help="N/O pairs at least this many bonds apart are checked")
    f.add_argument("--allow-cis-amide", action="store_true",
                   help="keep the amide cis/trans states of the input instead "
                        "of requiring trans")

    c = p.add_argument_group("convergence")
    c.add_argument("--min-confs", type=int, default=20, help="minimum conformers")
    c.add_argument("--max-confs", type=int, default=80, help="maximum conformers")
    c.add_argument("--sem-tol", type=float, default=0.001,
                   help="max per-atom SEM of the mean (e)")
    c.add_argument("--drift-tol", type=float, default=0.0005,
                   help="max per-atom change of the running mean over the last "
                        "--drift-window conformers (e)")
    c.add_argument("--drift-window", type=int, default=5,
                   help="window for the drift criterion")
    c.add_argument("--nboot", type=int, default=5000,
                   help="bootstrap resamples for the 95%% CI")

    a = p.add_argument_group("AmberTools")
    a.add_argument("--method", default="abcg2", help="antechamber charge method")
    a.add_argument("--atomtype", default="gaff2", help="antechamber atom types")
    a.add_argument("--ac-ek", default=None, help="sqm keywords passed with -ek")
    a.add_argument("-j", "--workers", type=int, default=os.cpu_count() or 1,
                   help="antechamber runs in parallel")
    a.add_argument("--ac-exe", default="antechamber", help="antechamber")
    a.add_argument("--parmchk2-exe", default="parmchk2", help="parmchk2")
    a.add_argument("--tleap-exe", default="tleap", help="tleap (for the .lib)")
    a.add_argument("--no-tleap", action="store_true", help="skip the tleap check")
    a.add_argument("--reference", choices=["auto", "input", "ideal"], default="auto",
                   help="reference conformer: input, idealised all-trans, or auto "
                        "(input if it passes the filters with no gauche torsion)")
    args = p.parse_args()
    if args.workers < 1 or args.min_confs < 2 or args.max_confs < args.min_confs:
        p.error("need --workers >= 1 and 2 <= --min-confs <= --max-confs")
    return args


# ==========================================================================
# main
# ==========================================================================
def main():
    args = parse_args()
    RDLogger.DisableLog("rdApp.*")
    args.ac_exe, args.parmchk2_exe, tleap = check_amber(args)

    m2 = Mol2(args.mol2)
    names = m2.names
    if len(set(names)) != len(names):
        sys.exit("ERROR: atom names in the mol2 must be unique.")
    prefix = args.prefix or os.path.splitext(os.path.basename(args.mol2))[0]
    out = os.path.abspath(args.outdir or f"{prefix}_{args.method}")
    if os.path.isdir(out) and os.listdir(out) and not args.overwrite:
        sys.exit(f"ERROR: {out} exists and is not empty (use --overwrite or -o).")
    for sub in ("conformers", "runs"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    shutil.copy(args.mol2, os.path.join(out, f"{prefix}_input.mol2"))
    resname = args.resname or m2.resname
    nc = args.net_charge
    if nc is None:
        nc = m2.net_charge_from_file() or 0

    log = open(os.path.join(out, f"{prefix}_run.log"), "w")

    def say(msg):
        print(msg, flush=True)
        log.write(msg + "\n")
        log.flush()

    try:
        mol = build_rdkit_mol(m2, nc, args.element_from)
    except (ValueError, RuntimeError) as exc:
        say(f"ERROR: {exc}")
        sys.exit(1)
    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
    classes = [int(c) for c in Chem.CanonicalRankAtoms(mol, breakTies=False)]
    checker = Checker(mol, args)
    X_in = mol.GetConformer().GetPositions()
    ev_in = checker.evaluate(X_in)
    g_in = gauche_count(mol, X_in)

    say(f"Molecule {smiles}: {mol.GetNumAtoms()} atoms, net charge {nc}, "
        f"residue {resname}")
    say(f"Output: {out}")
    say(f"Input geometry: Rg {ev_in['rg']:.2f} A, {g_in} gauche, issues: "
        f"{ev_in['reason'] or 'none'}")
    for line in checker.amide_report(X_in, names):
        say(f"  input amide {line}")
    if any(s == "cis" for s in checker.amide_states(X_in)) and not args.allow_cis_amide:
        say("  NOTE: the input has a cis secondary amide; it is not used as "
            "reference, and all conformers are required to be trans.")

    use_input = args.reference == "input" or (
        args.reference == "auto" and ev_in["ok"] and g_in == 0)
    X_ref0 = X_in if use_input else idealize_all_trans(mol, args)
    ev0 = checker.evaluate(X_ref0)
    say(f"Reference conformer: {'input' if use_input else 'idealised all-trans'} "
        f"({gauche_count(mol, X_ref0)} gauche, Rg {ev0['rg']:.2f} A, issues: "
        f"{ev0['reason'] or 'none'})")
    ref_start = os.path.join(out, f"{prefix}_reference_start.mol2")
    m2.write(ref_start, coords=X_ref0, title=f"{m2.mol_name} reference start")

    # ---------------- pool ----------------
    say(f"Building the pool: {args.nconf} ETKDGv3 conformers, {args.ff} "
        f"(eps = {args.dielectric}) ...")
    rows, pool = build_pool(mol, checker, args, say)
    with open(os.path.join(out, f"{prefix}_pool.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pool_id", "mmff_dE_kcal", "rg_A", "n_gauche", "n_hbond",
                    "n_ionic", "min_polar_dist_A", "amide_ok", "stereo_ok",
                    "accepted", "reason"])
        for r in sorted(rows, key=lambda r: r["pool_id"]):
            w.writerow([r["pool_id"], f"{r['mmff_dE']:.3f}", f"{r['rg']:.3f}",
                        r["n_gauche"], r["n_hbond"], r["n_ionic"],
                        f"{r['min_polar_dist']:.3f}", int(r["amide_ok"]),
                        int(r["stereo_ok"]), int(r["ok"]), r["reason"]])
    say(f"Pool: {len(rows)} embedded, {len(pool)} pass the filters after "
        f"duplicate removal")
    if len(pool) < args.min_confs:
        say("WARNING: fewer pool conformers than --min-confs; raise -n or relax "
            "the filters.")
    if not pool:
        sys.exit(1)

    order = []
    for k, r in enumerate(pool, 1):
        tag = f"c{k:03d}"
        fn = os.path.join(out, "conformers", f"{prefix}_{tag}.mol2")
        X = kabsch_transform(r["X"][checker.heavy], X_in[checker.heavy], r["X"])
        m2.write(fn, coords=X, title=f"{m2.mol_name} pool conformer {r['pool_id']}")
        order.append((tag, fn, r))

    # ---------------- antechamber runs ----------------
    say(f"Running {args.method} with {args.workers} worker(s); stop when N >= "
        f"{args.min_confs}, max SEM <= {args.sem_tol} e and max drift over "
        f"{args.drift_window} <= {args.drift_tol} e")
    accepted, results, Q, history, sem_hist, drift_hist = [], [], [], [], [], []
    n_rejected = 0
    converged = stop = False
    ref = None
    pending, done_by_index = {}, {}
    next_i = next_use = 0

    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        fut_ref = ex.submit(run_antechamber, "reference", ref_start,
                            os.path.join(out, "runs", "reference"), nc, resname, args)
        pending[fut_ref] = "ref"

        def fill():
            nonlocal next_i
            while not stop and next_i < len(order) and len(pending) < args.workers:
                tag, fn, _ = order[next_i]
                pending[ex.submit(run_antechamber, tag, fn,
                                  os.path.join(out, "runs", tag), nc, resname,
                                  args)] = next_i
                next_i += 1

        fill()
        while pending:
            done, _ = cf.wait(list(pending), return_when=cf.FIRST_COMPLETED)
            for fut in done:
                key = pending.pop(fut)
                if key == "ref":
                    ref = fut.result()
                    if ref["ok"]:
                        Xr = read_sqm_pdb(ref["sqm_pdb"], names)
                        ev = checker.evaluate(Xr)
                        ref["X"] = kabsch_transform(Xr[checker.heavy],
                                                    X_ref0[checker.heavy], Xr)
                        say(f"  reference: done in {ref['time_s']:.0f} s; after "
                            f"sqm {gauche_count(mol, Xr)} gauche, Rg "
                            f"{ev['rg']:.2f} A, issues: {ev['reason'] or 'none'}")
                    else:
                        say("  WARNING: the reference run failed.")
                else:
                    done_by_index[key] = fut.result()
            while next_use in done_by_index:
                res = done_by_index.pop(next_use)
                tag, fn, r = order[next_use]
                next_use += 1
                e = {"tag": tag, "pool_id": r["pool_id"], "mmff_dE": r["mmff_dE"],
                     "pre_rg": r["rg"], "pre_gauche": r["n_gauche"],
                     "time_s": res["time_s"], "hof": res.get("hof")}
                results.append(e)
                if not res["ok"]:
                    e.update(accepted=False, reason=res["error"])
                    n_rejected += 1
                    say(f"  {tag}: FAILED ({res['error']})")
                    continue
                Xp = read_sqm_pdb(res["sqm_pdb"], names)
                ev = checker.evaluate(Xp)
                e.update(post_rg=ev["rg"], post_gauche=gauche_count(mol, Xp),
                         post_min_polar=ev["min_polar_dist"])
                if not ev["ok"]:
                    e.update(accepted=False, reason=ev["reason"])
                    n_rejected += 1
                    say(f"  {tag}: rejected after sqm ({ev['reason']})")
                    continue
                if converged or len(accepted) >= args.max_confs:
                    e.update(accepted=False, reason="finished before use")
                    continue
                e.update(accepted=True, reason="", types=res["types"],
                         out_mol2=res["out_mol2"])
                accepted.append(e)
                Q.append([res["charges"][n] for n in names])
                mean, std, sem = stats(np.array(Q))
                history.append(mean)
                n = len(Q)
                msem = float(sem.max()) if n > 1 else float("nan")
                drift = float(np.abs(history[-1] - history[-1 - args.drift_window]).max()) \
                    if n > args.drift_window else float("nan")
                sem_hist.append(msem)
                drift_hist.append(drift)
                say(f"  {tag}: N={n:3d}  {res['time_s']:5.0f} s  max SEM {msem:.4f}  "
                    f"max drift {drift:.4f}  ({e['post_gauche']} gauche, Rg "
                    f"{ev['rg']:.2f} A)")
                if n >= args.min_confs and msem <= args.sem_tol and \
                        drift <= args.drift_tol:
                    converged = stop = True
                    say(f"Converged at N = {n}.")
                elif n >= args.max_confs:
                    stop = True
                    say(f"Reached --max-confs ({args.max_confs}) without converging.")
            fill()

    if not accepted:
        say("ERROR: no conformer was accepted.")
        sys.exit(1)
    if not converged:
        say("WARNING: not converged; the SEM is reported. Raise -n or --max-confs.")

    # ---------------- final charges ----------------
    Qa = np.array(Q)
    mean, std, sem = stats(Qa)
    lo, hi = bootstrap_ci(Qa, args.nboot, args.seed)
    q_final = finalize_charges(mean, nc, classes)
    ref_ok = ref is not None and ref["ok"]
    ref_q = np.array([ref["charges"][n] for n in names]) if ref_ok else None
    X_final = ref["X"] if ref_ok else X_ref0
    tmpl = Mol2(accepted[0]["out_mol2"])
    if tmpl.names != names:
        sys.exit("ERROR: antechamber changed the atom order.")
    types = accepted[0]["types"]
    final_mol2 = os.path.join(out, f"{prefix}_{args.method}_ensemble.mol2")
    tmpl.write(final_mol2, coords=X_final, charges=q_final,
               title=f"{m2.mol_name} {args.method} averaged over {len(Q)} conformers")
    if ref_ok:
        tmpl.write(os.path.join(out, f"{prefix}_{args.method}_reference_only.mol2"),
                   coords=X_final, charges=ref_q,
                   title=f"{m2.mol_name} {args.method} on the reference conformer")

    frcmod = os.path.join(out, f"{prefix}.frcmod")
    subprocess.run([args.parmchk2_exe, "-i", final_mol2, "-f", "mol2", "-o", frcmod,
                    "-s", args.atomtype], stdout=subprocess.DEVNULL,
                   stderr=subprocess.STDOUT)

    leap_note = "skipped"
    if tleap:
        ldir = os.path.join(out, "tleap")
        os.makedirs(ldir, exist_ok=True)
        leaprc = "leaprc.gaff2" if args.atomtype == "gaff2" else "leaprc.gaff"
        with open(os.path.join(ldir, "leap.in"), "w") as fh:
            fh.write(f"source {leaprc}\nloadamberparams ../{prefix}.frcmod\n"
                     f"{resname} = loadmol2 ../{os.path.basename(final_mol2)}\n"
                     f"charge {resname}\ncheck {resname}\n"
                     f"saveoff {resname} {prefix}.lib\nquit\n")
        with open(os.path.join(ldir, "leap.log"), "w") as lg:
            subprocess.run([tleap, "-f", "leap.in"], cwd=ldir, stdout=lg,
                           stderr=subprocess.STDOUT)
        txt = open(os.path.join(ldir, "leap.log")).read()
        mq = re.search(r"Total unperturbed charge:\s*(-?\d+\.\d+)", txt)
        me = re.search(r"Errors = (\d+); Warnings = (\d+)", txt)
        leap_note = (f"charge {mq.group(1) if mq else '?'}, "
                     f"errors {me.group(1) if me else '?'}, "
                     f"warnings {me.group(2) if me else '?'}")

    # ---------------- tables ----------------
    with open(os.path.join(out, f"{prefix}_charges_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["atom", "type", "sym_class", "q_final", "mean", "std", "sem",
                    "ci95_lo", "ci95_hi", "min", "max", "q_reference",
                    "final_minus_reference", "q_input_file"])
        for i, n in enumerate(names):
            qin = m2.atoms[i]["charge"]
            w.writerow([n, types[n], classes[i], f"{q_final[i]:.6f}", f"{mean[i]:.6f}",
                        f"{std[i]:.6f}", f"{sem[i]:.6f}", f"{lo[i]:.6f}",
                        f"{hi[i]:.6f}", f"{Qa[:, i].min():.6f}", f"{Qa[:, i].max():.6f}",
                        f"{ref_q[i]:.6f}" if ref_ok else "",
                        f"{q_final[i] - ref_q[i]:.6f}" if ref_ok else "",
                        f"{qin:.6f}" if (qin is not None and m2.has_charges()) else ""])
    with open(os.path.join(out, f"{prefix}_charges_per_conformer.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["atom"] + [e["tag"] for e in accepted])
        for i, n in enumerate(names):
            w.writerow([n] + [f"{q:.6f}" for q in Qa[:, i]])
    with open(os.path.join(out, f"{prefix}_running_mean.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["N", "max_sem", "max_drift"] + names)
        for k, h in enumerate(history):
            w.writerow([k + 1, f"{sem_hist[k]:.6f}", f"{drift_hist[k]:.6f}"] +
                       [f"{x:.6f}" for x in h])
    with open(os.path.join(out, f"{prefix}_conformer_runs.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tag", "pool_id", "mmff_dE_kcal", "pre_rg_A", "pre_gauche",
                    "post_rg_A", "post_gauche", "post_min_polar_dist_A",
                    "am1_hof_kcal", "time_s", "accepted", "reason"])
        for e in results:
            w.writerow([e["tag"], e["pool_id"], f"{e['mmff_dE']:.3f}",
                        f"{e['pre_rg']:.3f}", e["pre_gauche"],
                        f"{e['post_rg']:.3f}" if "post_rg" in e else "",
                        e.get("post_gauche", ""),
                        f"{e['post_min_polar']:.3f}" if "post_min_polar" in e else "",
                        f"{e['hof']:.4f}" if e.get("hof") is not None else "",
                        f"{e['time_s']:.1f}", int(e["accepted"]), e["reason"]])
    with open(os.path.join(out, f"{prefix}_accepted_sqm.pdb"), "w") as fh:
        for k, e in enumerate(accepted, 1):
            fh.write(f"MODEL {k:8d}\n")
            with open(os.path.join(out, "runs", e["tag"], "sqm.pdb")) as src:
                fh.writelines(l for l in src if l.startswith(("ATOM", "HETATM")))
            fh.write("ENDMDL\n")

    summary = {
        "input": os.path.abspath(args.mol2), "prefix": prefix, "smiles": smiles,
        "resname": resname, "net_charge": nc, "method": args.method,
        "atomtype": args.atomtype,
        "input_issues": ev_in["reason"],
        "input_amides": checker.amide_report(X_in, names),
        "reference_geometry": "input" if use_input else "idealised all-trans",
        "pool_embedded": len(rows), "pool_accepted": len(pool),
        "runs": len(results), "accepted": len(accepted), "rejected": n_rejected,
        "converged": converged, "final_max_sem": float(sem.max()),
        "final_max_drift": drift_hist[-1], "max_std": float(std.max()),
        "max_std_atom": names[int(std.argmax())],
        "max_abs_final_minus_reference":
            float(np.abs(q_final - ref_q).max()) if ref_ok else None,
        "tleap": leap_note, "files": {"mol2": os.path.basename(final_mol2),
                                      "frcmod": os.path.basename(frcmod)},
        "settings": {k: v for k, v in vars(args).items()},
    }
    with open(os.path.join(out, f"{prefix}_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    say(f"\nAccepted {len(accepted)} of {len(results)} runs ({n_rejected} rejected); "
        f"converged: {converged}")
    say(f"Final max SEM {sem.max():.4f} e; largest conformational SD "
        f"{std.max():.4f} e ({names[int(std.argmax())]})")
    if ref_ok:
        d = np.abs(q_final - ref_q)
        say(f"Largest |ensemble - reference conformer| {d.max():.4f} e "
            f"({names[int(d.argmax())]})")
    say(f"tleap check: {leap_note}")
    say(f"Charges: {final_mol2}\nfrcmod:  {frcmod}")
    say(f"Analyse with: python3 abcg2_analysis.py {out}")
    log.close()


if __name__ == "__main__":
    main()
