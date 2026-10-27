#!/usr/bin/env python3
"""
abcg2_covalent_charges.py

Conformer-averaged ABCG2 (or AM1-BCC) charges for a ligand covalently bound
to an amino-acid residue. The charges are computed on the capped residue
(ACE-residue-NME, with the ligand attached to the side chain); the charges of
the capping groups are then set to fixed values that mimic the peptide bonds
to the neighbouring residues, and the remaining charges are corrected so that
the residue carries an exact integer charge. Only a mol2 file is needed:

    python3 abcg2_covalent_charges.py capped_residue.mol2 -nc -1

Protocol
  1. Read the mol2 and perceive bond orders and formal charges (as in
     abcg2_conformer_charges.py). Find the ACE and NME caps and the residue
     backbone (N, H, CA, HA, C, O) with a SMARTS pattern, or take the cap
     atoms from --ace / --nme.
  2. Conformer pool: ETKDGv3 embedding with the caps and the backbone held at
     their input coordinates (--backbone fixed, the default), so that phi,
     psi and the two peptide bonds keep the geometry of the protein, and
     MMFF94s minimisation (eps = 4) with the same atoms frozen. Only the side
     chain and the ligand are sampled. With --backbone free the whole capped
     residue is sampled.
  3. Filters as in abcg2_conformer_charges.py (intramolecular H-bonds, ionic
     and N/O contacts, cis secondary amides, MMFF energy window), except for
     contacts between two frozen atoms that are already present in the input:
     they belong to the protein geometry and are not counted. The N/O
     distance filter, meant for distant groups of the ligand folding onto
     each other, only checks pairs of free atoms (--polar-frozen for all).
  4. AM1 geometry: sqm cannot restrain atoms, so with the backbone fixed the
     script minimises the free atoms itself (L-BFGS-B, AM1 energies and
     gradients from sqm single points, caps and backbone at the protein
     geometry) and antechamber -c abcg2 -at gaff2 then computes the charges
     at that restrained AM1 minimum (sqm maxcyc=0). An unrestrained AM1
     minimisation in vacuum folds the capped residue onto itself (C7 turns,
     N-H...O contacts); it is still available with --am1-opt full, and is
     the default with --backbone free. After the AM1 step the same filters
     are applied again, and a conformer is rejected if phi or psi moved by
     more than --bb-tol from the input.
  5. Running mean, SD and SEM per atom and the same stopping rule as in
     abcg2_conformer_charges.py. The input itself (the bound pose) is the
     single reference conformer.
  6. Final charges: ensemble mean, averaged over symmetry-equivalent atoms;
     the cap atoms (and, with --fix-backbone, N, H, C and O) are set to the
     Amber ff14SB/ff19SB values (ACE and NME of aminont12.lib/aminoct12.lib,
     both summing to 0), or to values from --fixed-charges; the difference
     between the ABCG2 charge of the fixed atoms and their new values is
     spread evenly over the free atoms (--redistribute even) or put back on
     the residue atom bonded to each cap (--redistribute junction). Rounded
     to 1e-6 e with an exact net charge.
  7. Outputs: the capped model (mol2, GAFF2 types), the residue without the
     caps (mol2, net charge exactly that of the residue), a parmchk2 frcmod
     and, if tleap is available, an Amber library of the residue with head
     and tail atoms set, checked inside an ACE-residue-NME peptide built from
     the protein force field. N, H, CA, HA, C and O get the protein atom
     types (--backbone-types amber); the bonded terms between these and the
     GAFF2 types that tleap cannot find are written to a junction frcmod
     (GAFF2 analogues from parmchk2, only the missing terms).

Outputs go to <input basename>_<method>/ (change with -o); the CSV/JSON files
there are read by abcg2_analysis.py, as for abcg2_conformer_charges.py.

Requirements: abcg2_conformer_charges.py in the same directory, AmberTools
with ABCG2 support (antechamber, sqm, parmchk2; tleap optional), Python 3.8+,
RDKit, NumPy and SciPy (installed with AmberTools from conda-forge; without it
a plain NumPy L-BFGS is used).
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import abcg2_conformer_charges as base
except ImportError:
    sys.exit("ERROR: abcg2_conformer_charges.py must be in the same directory.")

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom
from rdkit.Geometry import Point3D

# ACE-X-NME backbone: CH3(ACE) C(ACE) O(ACE) N CA C O N(NME) CH3(NME)
CAPPED_SMARTS = Chem.MolFromSmarts(
    "[CX4;H3]-[CX3](=[OX1])-[NX3]-[CX4]-[CX3](=[OX1])-[NX3;H1]-[CX4;H3]")

# Amber ff14SB / ff19SB charges (aminont12.lib ACE, aminoct12.lib NME and the
# backbone N, H, C, O shared by most residues of amino12.lib / amino19.lib).
ACE_Q = {"CH3": -0.3662, "HC": 0.1123, "C": 0.5972, "O": -0.5679}
NME_Q = {"N": -0.4157, "H": 0.2719, "CH3": -0.1490, "HC": 0.0976}
BB_Q = {"N": -0.4157, "H": 0.2719, "C": 0.5973, "O": -0.5679}

# Amber atom types of the backbone (for --backbone-types amber).
BB_TYPES = {"ff14SB": {"N": "N", "H": "H", "CA": "CX", "HA": "H1", "C": "C", "O": "O"},
            "ff19SB": {"N": "N", "H": "H", "CA": "XC", "HA": "H1", "C": "C", "O": "O"}}


# ==========================================================================
# caps and backbone
# ==========================================================================
def hydrogens_on(mol, idx):
    return [nb.GetIdx() for nb in mol.GetAtomWithIdx(idx).GetNeighbors()
            if nb.GetAtomicNum() == 1]


def heavy_neighbors(mol, idx):
    return [nb.GetIdx() for nb in mol.GetAtomWithIdx(idx).GetNeighbors()
            if nb.GetAtomicNum() > 1]


def find_backbone(mol, names, ace_names, nme_names):
    """Atom indices of the caps and of the residue backbone, with roles."""
    idx = {n: i for i, n in enumerate(names)}
    if ace_names or nme_names:
        if not (ace_names and nme_names):
            raise ValueError("give both --ace and --nme, or neither")
        try:
            ace = [idx[n] for n in ace_names]
            nme = [idx[n] for n in nme_names]
        except KeyError as exc:
            raise ValueError(f"atom {exc} of --ace/--nme is not in the mol2")
        ace_c = [i for i in ace if mol.GetAtomWithIdx(i).GetAtomicNum() == 6
                 and any(mol.GetAtomWithIdx(j).GetAtomicNum() == 8
                         for j in heavy_neighbors(mol, i))]
        nme_n = [i for i in nme if mol.GetAtomWithIdx(i).GetAtomicNum() == 7]
        if len(ace_c) != 1 or len(nme_n) != 1:
            raise ValueError("--ace needs one carbonyl C and --nme one N")
        n_res = [j for j in heavy_neighbors(mol, ace_c[0]) if j not in ace]
        c_res = [j for j in heavy_neighbors(mol, nme_n[0]) if j not in nme]
        if len(n_res) != 1 or len(c_res) != 1:
            raise ValueError("each cap must be bonded to exactly one residue atom")
        n, c = n_res[0], c_res[0]
        ca = [j for j in heavy_neighbors(mol, n) if j in heavy_neighbors(mol, c)]
        if len(ca) != 1:
            raise ValueError("could not find CA between the two caps")
        ca = ca[0]
        o = [j for j in heavy_neighbors(mol, c) if mol.GetAtomWithIdx(j).GetAtomicNum() == 8]
        match = (next(i for i in ace if i != ace_c[0] and mol.GetAtomWithIdx(i).GetAtomicNum() == 6),
                 ace_c[0],
                 next(i for i in ace if mol.GetAtomWithIdx(i).GetAtomicNum() == 8),
                 n, ca, c, o[0], nme_n[0],
                 next(i for i in nme if mol.GetAtomWithIdx(i).GetAtomicNum() == 6))
    else:
        matches = mol.GetSubstructMatches(CAPPED_SMARTS)
        if len(matches) != 1:
            raise ValueError(f"{len(matches)} ACE-X-NME patterns found; give the cap "
                             f"atoms with --ace and --nme")
        match = matches[0]
    ch3a, ca_, oa, n, ca, c, o, nn, ch3n = match
    for cap_c in (ch3a, ch3n):
        if len(heavy_neighbors(mol, cap_c)) != 1:
            raise ValueError(f"cap methyl {names[cap_c]} is not terminal")
    roles = {}
    roles[ch3a] = ("ACE", "CH3")
    for h in hydrogens_on(mol, ch3a):
        roles[h] = ("ACE", "HC")
    roles[ca_] = ("ACE", "C")
    roles[oa] = ("ACE", "O")
    roles[nn] = ("NME", "N")
    for h in hydrogens_on(mol, nn):
        roles[h] = ("NME", "H")
    roles[ch3n] = ("NME", "CH3")
    for h in hydrogens_on(mol, ch3n):
        roles[h] = ("NME", "HC")
    ace = sorted(i for i, r in roles.items() if r[0] == "ACE")
    nme = sorted(i for i, r in roles.items() if r[0] == "NME")
    if ace_names and set(ace) != {idx[x] for x in ace_names}:
        raise ValueError("--ace does not describe an acetyl group (CH3, 3 H, C, O)")
    if nme_names and set(nme) != {idx[x] for x in nme_names}:
        raise ValueError("--nme does not describe an N-methyl amide (N, H, CH3, 3 H)")
    bb = {"N": [n], "H": hydrogens_on(mol, n), "CA": [ca],
          "HA": hydrogens_on(mol, ca), "C": [c], "O": [o]}
    return {"ace": ace, "nme": nme, "roles": roles, "bb": bb,
            "phi": (ca_, n, ca, c), "psi": (n, ca, c, nn),
            "omega_ace": (ch3a, ca_, n, ca), "omega_nme": (ca, c, nn, ch3n)}


def wrap(d):
    return (d + 180.0) % 360.0 - 180.0


def backbone_torsions(bk, X):
    return {k: base.dihedral(*(X[i] for i in bk[k]))
            for k in ("phi", "psi", "omega_ace", "omega_nme")}


# ==========================================================================
# filters that ignore contacts already present in the frozen part
# ==========================================================================
class CovalentChecker(base.Checker):
    """Contacts between two frozen atoms that exist in the input are part of
    the protein geometry and are not counted (they are counted if they only
    appear later, e.g. after the unrestrained AM1 minimisation)."""

    def __init__(self, mol, args, frozen, X_in):
        super().__init__(mol, args)
        fz = set(frozen)
        a = args
        self.skip_hb, self.skip_ion, self.skip_polar = set(), set(), set()
        for h, d in self.donors_h:
            for acc in self.acceptors:
                if h in fz and acc in fz and self._is_hb(X_in, h, d, acc):
                    self.skip_hb.add((h, acc))
        for k, i in enumerate(self.charged):
            for j in self.charged[k + 1:]:
                if i in fz and j in fz and np.linalg.norm(X_in[i] - X_in[j]) <= a.ion_dist:
                    self.skip_ion.add((i, j))
        for i, j in self.polar_pairs:
            if i in fz and j in fz and np.linalg.norm(X_in[i] - X_in[j]) < a.polar_dist:
                self.skip_polar.add((i, j))
        # The N/O distance filter is meant for distant groups of the ligand
        # folding onto each other; pairs with a frozen (cap or backbone) atom
        # are set by the protein geometry and are left out unless asked for.
        if fz and not a.polar_frozen:
            self.polar_pairs = [(i, j) for i, j in self.polar_pairs
                                if i not in fz and j not in fz]

    def _is_hb(self, X, h, d, acc):
        a = self.args
        if acc == d or self.topo[d, acc] < a.min_topo:
            return False
        if np.linalg.norm(X[h] - X[acc]) > a.hb_dist:
            return False
        v1, v2 = X[d] - X[h], X[acc] - X[h]
        cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return math.degrees(math.acos(max(-1.0, min(1.0, cosang)))) >= a.hb_angle

    def evaluate(self, X):
        a = self.args
        Xh = X[self.heavy]
        rg = math.sqrt(((Xh - Xh.mean(axis=0)) ** 2).sum(axis=1).mean())
        hb = sum(1 for h, d in self.donors_h for acc in self.acceptors
                 if (h, acc) not in self.skip_hb and self._is_hb(X, h, d, acc))
        ion = sum(1 for k, i in enumerate(self.charged) for j in self.charged[k + 1:]
                  if (i, j) not in self.skip_ion and self.topo[i, j] >= a.min_topo
                  and np.linalg.norm(X[i] - X[j]) <= a.ion_dist)
        dpol = min((float(np.linalg.norm(X[i] - X[j])) for i, j in self.polar_pairs
                    if (i, j) not in self.skip_polar), default=float("inf"))
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

    def contacts(self, X, names):
        """Readable list of the contacts flagged in X (for the input report)."""
        a, out = self.args, []
        for h, d in self.donors_h:
            for acc in self.acceptors:
                if (h, acc) not in self.skip_hb and self._is_hb(X, h, d, acc):
                    out.append(f"H-bond {names[d]}-{names[h]}...{names[acc]} "
                               f"{np.linalg.norm(X[h] - X[acc]):.2f} A")
        for i, j in self.polar_pairs:
            dd = float(np.linalg.norm(X[i] - X[j]))
            if (i, j) not in self.skip_polar and dd < a.polar_dist:
                out.append(f"polar contact {names[i]}...{names[j]} {dd:.2f} A")
        return out


# ==========================================================================
# conformer pool with the backbone held fixed
# ==========================================================================
def build_pool(mol, checker, frozen, X_in, args, say):
    molc = Chem.Mol(mol)
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = args.seed
    p.numThreads = args.threads
    p.pruneRmsThresh = args.prune_rms
    p.useRandomCoords = True
    p.enforceChirality = True
    if frozen:
        p.SetCoordMap({i: Point3D(*X_in[i]) for i in frozen})
    cids = list(rdDistGeom.EmbedMultipleConfs(molc, numConfs=args.nconf, params=p))
    if not cids:
        raise RuntimeError("ETKDG could not embed any conformer")
    energies, used = {}, args.ff
    for cid in cids:
        conf = molc.GetConformer(cid)
        if frozen:
            # put the frozen atoms exactly at the input coordinates
            X = conf.GetPositions()
            X = base.kabsch_transform(X[frozen], X_in[frozen], X)
            X[frozen] = X_in[frozen]
            for i, xyz in enumerate(X):
                conf.SetAtomPosition(i, Point3D(*xyz))
        ff, used = base.mmff_forcefield(molc, args, cid)
        for i in frozen:
            ff.AddFixedPoint(int(i))
        ff.Minimize(maxIts=args.max_iters)
        energies[cid] = ff.CalcEnergy()
    if used != args.ff:
        say("WARNING: MMFF parameters missing, UFF used for the pool.")
    emin = min(energies.values())
    ref_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))

    rows = []
    for cid in cids:
        X = molc.GetConformer(cid).GetPositions()
        if not frozen:
            X = base.kabsch_transform(X[checker.heavy], X_in[checker.heavy], X)
        m = Chem.Mol(molc)
        m.RemoveAllConformers()
        m.AddConformer(Chem.Conformer(molc.GetConformer(cid)), assignId=True)
        Chem.AssignStereochemistryFrom3D(m)
        stereo_ok = Chem.MolToSmiles(Chem.RemoveHs(m)) == ref_smiles
        row = {"pool_id": cid, "mmff_dE": energies[cid] - emin,
               "stereo_ok": stereo_ok, "n_gauche": base.gauche_count(mol, X), "X": X}
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
    perms = base.heavy_atom_permutations(mol, heavy, args.max_sym)
    if frozen:
        # all conformers share the frozen frame: compare the free heavy atoms
        # directly (the frozen ones would only dilute the RMSD)
        fz = set(frozen)
        sel = np.array([k for k, i in enumerate(heavy) if i not in fz])

        def rmsd(P, Q):
            return min(math.sqrt(((P[p][sel] - Q[sel]) ** 2).sum(axis=1).mean())
                       for p in perms)
    else:
        def rmsd(P, Q):
            return base.sym_rmsd(P, Q, perms, False)
    kept = []
    for r in sorted((r for r in rows if r["ok"]), key=lambda r: r["mmff_dE"]):
        if all(rmsd(r["X"][heavy], k["X"][heavy]) >= args.div_rmsd for k in kept):
            kept.append(r)
        else:
            r["ok"], r["reason"] = False, "duplicate"
    random.Random(args.seed).shuffle(kept)
    return rows, kept


# ==========================================================================
# AM1 minimisation with frozen atoms (sqm has no restraints of its own)
# ==========================================================================
FLOAT_RE = re.compile(r"-?\d+\.\d+(?:[DdEe][-+]?\d+)?")


class SqmGradient:
    """AM1 heat of formation (kcal/mol) and gradient (kcal/mol/A) from sqm
    single points (maxcyc=0, verbosity=4 prints dE/dx)."""

    def __init__(self, exe, atomic_numbers, names, charge, wdir, keywords):
        self.exe, self.z, self.names = exe, atomic_numbers, names
        self.charge, self.wdir, self.keywords = charge, wdir, keywords
        self.calls = 0

    def __call__(self, X):
        self.calls += 1
        inp = os.path.join(self.wdir, "sqm_opt.in")
        outp = os.path.join(self.wdir, "sqm_opt.out")
        with open(inp, "w") as fh:
            fh.write("AM1 single point for the restrained minimisation\n &qmmm\n"
                     f"  {self.keywords}, qmcharge={self.charge}, maxcyc=0, "
                     f"verbosity=4,\n /\n")
            for z, n, x in zip(self.z, self.names, X):
                fh.write(f"{z:4d} {n:>5s} {x[0]:14.8f} {x[1]:14.8f} {x[2]:14.8f}\n")
        subprocess.run([self.exe, "-O", "-i", inp, "-o", outp],
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        txt = open(outp).read() if os.path.isfile(outp) else ""
        m = re.search(r"Heat of formation\s*=\s*(-?\d+\.\d+)", txt)
        if not m or "Forces on QM atoms" not in txt:
            raise RuntimeError("sqm single point failed (see sqm_opt.out)")
        g = []
        for line in txt.split("Forces on QM atoms", 1)[1].splitlines():
            if "Atm" in line:
                g.append([float(v.replace("D", "E")) for v in
                          FLOAT_RE.findall(line.split(":")[-1])])
        try:
            g = np.array(g, dtype=float)
        except ValueError:
            raise RuntimeError("could not read the sqm gradient")
        if g.shape != X.shape or not np.all(np.isfinite(g)):
            raise RuntimeError("could not read the sqm gradient")
        return float(m.group(1)), g


def minimise(fun, X0, free, gtol_rms, gtol_max, maxiter, log=None):
    """Minimise over the coordinates of the free atoms (the others stay
    fixed): SciPy L-BFGS-B when available (it comes with AmberTools from
    conda-forge), otherwise the plain L-BFGS below."""
    try:
        from scipy.optimize import minimize
    except ImportError:
        return lbfgs_minimise(fun, X0, free, gtol_rms, gtol_max, maxiter, log=log)
    free = np.asarray(free)
    best = {"e": None, "X": X0.copy(), "g": None}
    step = [0]

    def f(v):
        X = X0.copy()
        X[free] = v.reshape(-1, 3)
        # a trial step of the line search that brings atoms together: refuse
        # it with a high energy instead of asking sqm
        d = np.linalg.norm(X[free][:, None, :] - X[None, :, :], axis=2)
        d[np.arange(len(free)), free] = np.inf
        if best["e"] is not None and d.min() < 0.7:
            return best["e"] + 1.0e4, -1.0e2 * best["g"].ravel()
        e, g = fun(X)
        gf = g[free]
        if best["e"] is None or e < best["e"]:
            best.update(e=e, X=X, g=gf)
        return e, gf.ravel()

    def cb(v):
        step[0] += 1
        if log and best["g"] is not None:
            gf = best["g"]
            log.write(f"{step[0]:5d} {best['e']:16.8f} "
                      f"{math.sqrt((gf ** 2).mean()):10.5f} {np.abs(gf).max():10.5f}\n")
        if best["g"] is not None and math.sqrt((best["g"] ** 2).mean()) <= gtol_rms \
                and np.abs(best["g"]).max() <= gtol_max:
            raise StopIteration

    try:
        minimize(f, X0[free].ravel(), jac=True, method="L-BFGS-B", callback=cb,
                 options={"maxiter": maxiter, "maxcor": 30, "gtol": 1e-4,
                          "ftol": 1e-14})
    except StopIteration:
        pass
    gf = best["g"]
    grms = math.sqrt((gf ** 2).mean())
    return best["X"], best["e"], step[0], float(grms), \
        bool(grms <= gtol_rms and np.abs(gf).max() <= gtol_max)


def lbfgs_minimise(fun, X0, free, gtol_rms, gtol_max, maxiter, maxstep=0.2, mem=10,
                   log=None):
    """L-BFGS over the coordinates of the free atoms; the others stay fixed."""
    mask = np.zeros(X0.shape, dtype=bool)
    mask[free] = True
    X = X0.copy()
    e, g = fun(X)
    g[~mask] = 0.0
    S, Y = [], []
    it = 0
    for it in range(1, maxiter + 1):
        gf = g[mask]
        grms, gmax = math.sqrt((gf ** 2).mean()), float(np.abs(gf).max())
        if log:
            log.write(f"{it - 1:5d} {e:16.8f} {grms:10.5f} {gmax:10.5f}\n")
        if grms <= gtol_rms and gmax <= gtol_max:
            return X, e, it - 1, float(grms), True
        q = g.ravel().copy()
        alphas = []
        for s, y in reversed(list(zip(S, Y))):
            rho = 1.0 / (y @ s)
            a = rho * (s @ q)
            q -= a * y
            alphas.append((rho, a, s, y))
        if S:
            q *= (S[-1] @ Y[-1]) / (Y[-1] @ Y[-1])
        else:
            q *= 0.05 / max(gmax, 1e-8)      # first step: at most 0.05 A
        for rho, a, s, y in reversed(alphas):
            b = rho * (y @ q)
            q += (a - b) * s
        d = -q.reshape(X.shape)
        d[~mask] = 0.0
        slope = float((g * d).sum())
        if slope >= 0:                       # not a descent direction: restart
            S, Y = [], []
            d = -g * (0.05 / max(gmax, 1e-8))
            slope = float((g * d).sum())
        dmax = float(np.linalg.norm(d, axis=1).max())
        if dmax > maxstep:
            d *= maxstep / dmax
            slope *= maxstep / dmax
        t, accepted = 1.0, False
        for _ in range(12):
            Xn = X + t * d
            en, gn = fun(Xn)
            gn[~mask] = 0.0
            if en <= e + 1e-4 * t * slope:
                accepted = True
                break
            t *= 0.5
        if not accepted:
            if not S:
                return X, e, it, grms, False
            S, Y = [], []
            continue
        s, y = (Xn - X).ravel(), (gn - g).ravel()
        if y @ s > 1e-10:
            S.append(s)
            Y.append(y)
            if len(S) > mem:
                S.pop(0)
                Y.pop(0)
        X, e, g = Xn, en, gn
    gf = g[mask]
    return X, e, it, math.sqrt((gf ** 2).mean()), False


# ==========================================================================
# antechamber (with optional reuse of an earlier run)
# ==========================================================================
def same_coordinates(f1, f2, tol=1e-4):
    try:
        a, b = base.Mol2(f1), base.Mol2(f2)
    except (OSError, ValueError, KeyError):
        return False
    return a.names == b.names and all(
        np.allclose(x["xyz"], y["xyz"], atol=tol) for x, y in zip(a.atoms, b.atoms))


def collect(tag, wdir, args, dt, extra):
    outf = os.path.join(wdir, f"{tag}_{args.method}.mol2")
    out = base.Mol2(outf)
    res = {"tag": tag, "ok": True, "time_s": dt, "out_mol2": outf,
           "sqm_pdb": os.path.join(wdir, "sqm.pdb"),
           "charges": {a["name"]: a["charge"] for a in out.atoms},
           "types": {a["name"]: a["type"] for a in out.atoms},
           "hof": base.read_heat_of_formation(os.path.join(wdir, "sqm.out"))}
    res.update(extra)
    return res


def run_charges(tag, mol2_path, wdir, net_charge, resname, args, ctx):
    """Charges of one conformer: AM1 minimisation (restrained, full or none)
    followed by antechamber; or the results of an earlier run (--reuse)."""
    os.makedirs(wdir, exist_ok=True)
    if args.reuse:
        src = os.path.join(args.reuse, "runs", tag)
        try:
            with open(os.path.join(src, "run.json")) as fh:
                meta = json.load(fh)
        except (OSError, ValueError):
            meta = {}
        if meta.get("am1_opt") == args.am1_opt and meta.get("ok") and \
                same_coordinates(os.path.join(src, "start.mol2"), mol2_path):
            for f in os.listdir(src):
                if os.path.isfile(os.path.join(src, f)):
                    shutil.copy(os.path.join(src, f), os.path.join(wdir, f))
            return collect(tag, wdir, args, 0.0,
                           {"reused": True, "am1_steps": meta.get("am1_steps"),
                            "am1_grms": meta.get("am1_grms")})
    shutil.copy(mol2_path, os.path.join(wdir, "start.mol2"))
    t0 = time.time()
    meta = {"am1_opt": args.am1_opt, "ok": False}
    extra = {}
    if args.am1_opt == "full":
        res = base.run_antechamber(tag, mol2_path, wdir, net_charge, resname, args)
    else:
        start = mol2_path
        if args.am1_opt == "restrained":
            m2 = base.Mol2(mol2_path)
            X0 = np.array([a["xyz"] for a in m2.atoms])
            fun = SqmGradient(ctx["sqm"], ctx["z"], m2.names, net_charge, wdir,
                              args.sqm_keywords)
            try:
                with open(os.path.join(wdir, "am1_restrained.log"), "w") as lg:
                    lg.write("# step  HOF (kcal/mol)  grms  gmax (kcal/mol/A), "
                             "free atoms only\n")
                    X, e, steps, grms, conv = minimise(
                        fun, X0, ctx["free_opt"], args.am1_grms, args.am1_gmax,
                        args.am1_maxiter, log=lg)
            except RuntimeError as exc:
                return {"tag": tag, "ok": False, "error": str(exc),
                        "time_s": time.time() - t0}
            if not args.keep_scratch:
                for f in ("sqm_opt.in", "sqm_opt.out"):
                    if os.path.isfile(os.path.join(wdir, f)):
                        os.remove(os.path.join(wdir, f))
            start = os.path.join(wdir, "am1_restrained.mol2")
            m2.write(start, coords=X, title=f"{tag} AM1 minimum, frozen atoms fixed")
            extra = {"am1_steps": steps, "am1_grms": grms, "am1_converged": conv,
                     "am1_calls": fun.calls}
            meta.update(am1_steps=int(steps), am1_grms=float(grms),
                        am1_converged=bool(conv))
        res = base.run_antechamber(tag, start, wdir, net_charge, resname, args)
    res["time_s"] = time.time() - t0
    res.update(extra)
    meta["ok"] = res["ok"]
    with open(os.path.join(wdir, "run.json"), "w") as fh:
        json.dump(meta, fh)
    return res


# ==========================================================================
# final charges with fixed atoms
# ==========================================================================
def fixed_values(bk, names, args):
    """{atom index: fixed charge} for the caps, the optional backbone atoms
    and the atoms of --fixed-charges (which override the others)."""
    fixed = {}
    for i, (cap, role) in bk["roles"].items():
        fixed[i] = (ACE_Q if cap == "ACE" else NME_Q)[role]
    if args.fix_backbone:
        for role in ("N", "H", "C", "O"):
            for i in bk["bb"][role]:
                fixed[i] = BB_Q[role]
    if args.fixed_charges:
        idx = {n: i for i, n in enumerate(names)}
        with open(args.fixed_charges) as fh:
            for ln, line in enumerate(fh, 1):
                s = line.split("#")[0].split()
                if not s:
                    continue
                if len(s) != 2 or s[0] not in idx:
                    raise ValueError(f"{args.fixed_charges} line {ln}: expected "
                                     f"'atom_name charge' with a name from the mol2")
                fixed[idx[s[0]]] = float(s[1])
    return fixed


def finalize_with_fixed(raw, net_charge, classes, fixed, bk, mode):
    """Symmetry-average the raw charges, set the fixed atoms and correct the
    free atoms so that the total is exactly net_charge."""
    n = len(raw)
    q = raw.copy()
    for cls in set(classes):
        idx = [i for i, c in enumerate(classes) if c == cls]
        q[idx] = q[idx].mean()
    sym = q.copy()
    free = [i for i in range(n) if i not in fixed]
    if not free:
        raise ValueError("all atoms are fixed")
    info = {}
    for cap in ("ace", "nme"):
        info[f"{cap}_abcg2_sum"] = round(float(sym[bk[cap]].sum()), 6) + 0.0
        info[f"{cap}_fixed_sum"] = round(float(sum(fixed[i] for i in bk[cap])), 6) + 0.0
    for i, v in fixed.items():
        q[i] = v
    if mode == "junction":
        n_at, c_at = bk["bb"]["N"][0], bk["bb"]["C"][0]
        if n_at in fixed or c_at in fixed:
            raise ValueError("--redistribute junction needs N and C to be free "
                             "(not with --fix-backbone)")
        q[n_at] += info["ace_abcg2_sum"] - info["ace_fixed_sum"]
        q[c_at] += info["nme_abcg2_sum"] - info["nme_fixed_sum"]
    resid = net_charge - q.sum()
    q[free] += resid / len(free)
    info["shift_per_free_atom"] = float(resid / len(free))
    q = np.round(q, 6)
    resid = round(net_charge - q.sum(), 6)
    if abs(resid) > 0:
        counts = {c: classes.count(c) for c in set(classes)}
        unique = [i for i in free if counts[classes[i]] == 1]
        k = unique[0] if unique else free[int(np.argmax(np.abs(q[free])))]
        q[k] = round(q[k] + resid, 6)
    return q, sym, info


# ==========================================================================
# mol2 of the residue without the caps
# ==========================================================================
def write_residue_mol2(path, m2, keep, coords, charges, types, resname, title):
    new = {old: k + 1 for k, old in enumerate(keep)}
    bonds = [(new[i], new[j], t) for i, j, t in m2.bonds if i in new and j in new]
    L = ["@<TRIPOS>MOLECULE", title,
         f"{len(keep):5d} {len(bonds):5d}     1     0     0", "SMALL",
         "USER_CHARGES", "", "", "@<TRIPOS>ATOM"]
    for old in keep:
        a = m2.atoms[old]
        x = coords[old]
        L.append(f"{new[old]:7d} {a['name']:<8s} {x[0]:10.4f} {x[1]:10.4f} "
                 f"{x[2]:10.4f} {types[old]:<6s} {1:>4d} {resname:<8s} "
                 f"{charges[old]:10.6f}")
    L.append("@<TRIPOS>BOND")
    for k, (i, j, t) in enumerate(bonds, 1):
        L.append(f"{k:6d} {i:5d} {j:5d} {t}")
    L += ["@<TRIPOS>SUBSTRUCTURE",
          f"     1 {resname:<8s}    1 RESIDUE           4 A     {resname:<4s}  2 ROOT", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(L))


# ==========================================================================
# tleap: library of the residue and a capped test peptide
# ==========================================================================
def run_tleap(tleap, ldir, script, stem):
    with open(os.path.join(ldir, f"{stem}.in"), "w") as fh:
        fh.write(script)
    with open(os.path.join(ldir, f"{stem}.log"), "w") as lg:
        subprocess.run([tleap, "-f", f"{stem}.in"], cwd=ldir, stdout=lg,
                       stderr=subprocess.STDOUT)
    with open(os.path.join(ldir, f"{stem}.log")) as fh:
        return fh.read()


MISSING_RE = [("bond", re.compile(r"Could not find bond parameter for atom types:\s*(.+)")),
              ("angle", re.compile(r"Could not find angle parameter for atom types:\s*(.+)")),
              ("dihe", re.compile(r"No torsion terms for atom types:\s*(.+)"))]


def missing_terms(txt):
    """Missing bonded terms reported by tleap, as (kind, (types...)) keys
    that do not depend on the direction in which the term is written."""
    out = set()
    for kind, rx in MISSING_RE:
        for m in rx.finditer(txt):
            t = tuple(x.strip() for x in m.group(1).strip().split("-"))
            out.add((kind, min(t, t[::-1])))
    return out


def leap_summary(txt):
    mq = re.findall(r"Total unperturbed charge:\s*(-?\d+\.\d+)", txt)
    me = re.findall(r"Errors = (\d+); Warnings = (\d+)", txt)
    return {"charge": mq[0] if mq else "?", "errors": me[0][0] if me else "?",
            "warnings": me[0][1] if me else "?", "missing": len(missing_terms(txt))}


def frcmod_terms(path):
    """{(kind, types): [lines]} for the BOND, ANGLE and DIHE sections."""
    width = {"BOND": 5, "ANGL": 8, "DIHE": 11}
    kind_of = {"BOND": "bond", "ANGL": "angle", "DIHE": "dihe"}
    terms, sec = {}, None
    with open(path) as fh:
        for line in fh:
            head = line.strip()[:4]
            if head in ("MASS", "BOND", "ANGL", "DIHE", "IMPR", "NONB"):
                sec = head
                continue
            if sec not in width or not line.strip():
                continue
            t = tuple(x.strip() for x in line[:width[sec]].split("-"))
            key = (kind_of[sec], min(t, t[::-1]))
            terms.setdefault(key, []).append(line.rstrip("\n"))
    return terms


def write_junction(path, wanted, terms, resname):
    sec = {"bond": "BOND", "angle": "ANGLE", "dihe": "DIHE"}
    L = [f"Missing protein/GAFF2 junction terms of {resname}, GAFF2 analogues "
         f"from parmchk2", "MASS", ""]
    found = set()
    for kind in ("bond", "angle", "dihe"):
        L.append(sec[kind])
        for key in sorted(k for k in wanted if k[0] == kind):
            if key in terms:
                L += terms[key]
                found.add(key)
        L.append("")
    L += ["IMPROPER", "", "NONBON", "", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(L))
    return wanted - found


def build_library(args, tleap, out, prefix, resname, res_mol2, frcmod, ldir, say):
    """Amber library of the residue (head N, tail C) and a check inside an
    ACE-residue-NME peptide built with the protein force field. Terms missing
    at the junction between the protein and the GAFF2 atom types are taken
    from parmchk2 (GAFF2 analogues of the protein types, PARMCHK.DAT); only
    the terms that tleap reports as missing are kept, so that no protein
    parameter is replaced."""
    os.makedirs(ldir, exist_ok=True)
    head = (f"source leaprc.protein.{args.protein_ff}\nsource leaprc.gaff2\n"
            f"loadamberparams ../{os.path.basename(frcmod)}\n")
    unit = (f"{resname} = loadmol2 ../{os.path.basename(res_mol2)}\n"
            f"set {resname} head {resname}.1.N\n"
            f"set {resname} tail {resname}.1.C\n"
            f"set {resname}.1 connect0 {resname}.1.N\n"
            f"set {resname}.1 connect1 {resname}.1.C\n"
            f"set {resname}.1 restype protein\n"
            f"set {resname}.1 name \"{resname}\"\n")
    pep = (f"pep = sequence {{ ACE {resname} NME }}\n"
           f"charge {resname}\ncheck {resname}\n"
           f"saveoff {resname} {prefix}.lib\n"
           f"savemol2 pep {prefix}_pep_types.mol2 1\n"
           f"charge pep\ncheck pep\n"
           f"saveamberparm pep {prefix}_pep.prmtop {prefix}_pep.inpcrd\n"
           f"savepdb pep {prefix}_pep.pdb\nquit\n")
    txt = run_tleap(tleap, ldir, head + unit + pep, "leap_first")
    first = leap_summary(txt)
    wanted = missing_terms(txt)
    junction, unresolved = None, set()
    if wanted and args.junction_params:
        allf = os.path.join(ldir, "parmchk2_all_terms.frcmod")
        # PARMCHK.DAT knows CX (ff14SB CA) as c3 but not XC (ff19SB CA)
        atc = os.path.join(ldir, "parmchk2_extra_types.dat")
        with open(atc, "w") as fh:
            fh.write("PARM\tXC\t0\t\t0\t12.01\t0\t6\nEQUA\tc3\n#\n")
        subprocess.run([args.parmchk2_exe, "-i", os.path.join(ldir, f"{prefix}_pep_types.mol2"),
                        "-f", "mol2", "-o", allf, "-s", args.atomtype, "-atc", atc],
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        if os.path.isfile(allf):
            junction = os.path.join(out, f"{prefix}_junction.frcmod")
            unresolved = write_junction(junction, wanted, frcmod_terms(allf), resname)
            head += f"loadamberparams ../{os.path.basename(junction)}\n"
            if unresolved:
                say("WARNING: no parmchk2 term for " + ", ".join(
                    "-".join(k[1]) for k in sorted(unresolved)))
    for f in os.listdir(ldir):
        if f.startswith(f"{prefix}_pep"):
            os.remove(os.path.join(ldir, f))
    txt = run_tleap(tleap, ldir, head + unit + pep, "leap")
    final = leap_summary(txt)
    first["terms"] = sorted("-".join(k[1]) + f" ({k[0]})" for k in wanted)
    return first, final, junction


# ==========================================================================
# command line
# ==========================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Conformer-averaged ABCG2 charges for a ligand covalently "
                    "bound to an amino-acid residue, computed on the capped "
                    "residue (ACE-X-NME) with fixed cap charges. Usage: "
                    "abcg2_covalent_charges.py capped_residue.mol2 -nc <charge>",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("mol2", help="input mol2: the capped residue with the ligand, "
                                "all H present, unique atom names")

    io = p.add_argument_group("input/output")
    io.add_argument("-o", "--outdir", default=None,
                    help="output directory (None = <input basename>_<method>)")
    io.add_argument("--prefix", default=None, help="file prefix (None = input basename)")
    io.add_argument("-nc", "--net-charge", type=int, default=None,
                    help="net charge of the capped model, equal to that of the "
                         "residue since the caps are neutral (None = rounded sum "
                         "of the mol2 charges, or 0)")
    io.add_argument("-rn", "--resname", default=None,
                    help="residue name (None = from the mol2)")
    io.add_argument("--element-from", choices=["auto", "type", "name"],
                    default="auto", help="how elements are read from the mol2")
    io.add_argument("--keep-scratch", action="store_true",
                    help="keep antechamber scratch files")
    io.add_argument("--overwrite", action="store_true",
                    help="allow writing into an existing output directory")
    io.add_argument("--reuse", default=None, metavar="DIR",
                    help="output directory of an earlier run: conformers with "
                         "identical coordinates reuse its antechamber results "
                         "(to change only how the final charges are built)")

    cv = p.add_argument_group("caps and backbone")
    cv.add_argument("--ace", default=None,
                    help="comma-separated atom names of the ACE cap (CH3, 3 H, "
                         "C, O); None = found with a SMARTS pattern")
    cv.add_argument("--nme", default=None,
                    help="comma-separated atom names of the NME cap (N, H, CH3, "
                         "3 H); None = found with a SMARTS pattern")
    cv.add_argument("--backbone", choices=["fixed", "free"], default="fixed",
                    help="keep the caps and N, H, CA, HA, C, O at the input "
                         "coordinates while the pool is built (fixed), or sample "
                         "the whole capped residue (free)")
    cv.add_argument("--keep-atoms", default=None,
                    help="comma-separated names of further atoms kept at the "
                         "input coordinates while the pool is built")
    cv.add_argument("--bb-tol", type=float, default=60.0,
                    help="largest change of phi or psi in the AM1 minimisation "
                         "(deg) for a conformer to be accepted (--backbone fixed); "
                         "0 disables the check")
    cv.add_argument("--fix-backbone", action="store_true",
                    help="also fix N, H, C and O to the Amber protein values "
                         "(-0.4157, 0.2719, 0.5973, -0.5679)")
    cv.add_argument("--fixed-charges", default=None, metavar="FILE",
                    help="file with 'atom_name charge' lines: atoms fixed to "
                         "these values (they override the cap and backbone values)")
    cv.add_argument("--redistribute", choices=["even", "junction"], default="even",
                    help="where the difference between the ABCG2 and the fixed "
                         "charges goes: spread evenly over all free atoms, or "
                         "onto N (from ACE) and C (from NME)")

    g = p.add_argument_group("conformer pool")
    g.add_argument("-n", "--nconf", type=int, default=500,
                   help="conformers embedded with ETKDGv3")
    g.add_argument("--seed", type=int, default=2026, help="random seed")
    g.add_argument("--threads", type=int, default=0,
                   help="threads for embedding (0 = all)")
    g.add_argument("--prune-rms", type=float, default=0.3, help="ETKDG pruning RMSD (A)")
    g.add_argument("--ff", choices=["mmff94s", "mmff94", "uff"], default="mmff94s",
                   help="force field for pool minimisation")
    g.add_argument("--dielectric", type=float, default=4.0, help="MMFF dielectric constant")
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
    f.add_argument("--polar-frozen", action="store_true",
                   help="also apply the N/O distance filter to pairs that "
                        "include a frozen atom (by default only pairs of free "
                        "atoms are checked)")
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
    a.add_argument("--am1-opt", choices=["auto", "restrained", "full", "none"],
                   default="auto",
                   help="geometry used for the charges: AM1 minimum with the "
                        "frozen atoms fixed (restrained), the standard sqm "
                        "minimisation of all atoms (full), or the MMFF "
                        "geometry as it is (none); auto = restrained with "
                        "--backbone fixed, full with --backbone free")
    a.add_argument("--am1-grms", type=float, default=0.02,
                   help="restrained AM1: RMS gradient of the free atoms "
                        "(kcal/mol/A) for convergence")
    a.add_argument("--am1-gmax", type=float, default=0.1,
                   help="restrained AM1: largest gradient component (kcal/mol/A)")
    a.add_argument("--am1-maxiter", type=int, default=1000,
                   help="restrained AM1: maximum L-BFGS steps")
    a.add_argument("--sqm-keywords", default="qm_theory='AM1', grms_tol=0.0005, "
                   "scfconv=1.d-10, ndiis_attempts=700",
                   help="sqm keywords (those used by antechamber by default)")
    a.add_argument("--sqm-exe", default="sqm", help="sqm (for the restrained AM1)")
    a.add_argument("-j", "--workers", type=int, default=os.cpu_count() or 1,
                   help="antechamber runs in parallel")
    a.add_argument("--ac-exe", default="antechamber", help="antechamber")
    a.add_argument("--parmchk2-exe", default="parmchk2", help="parmchk2")
    a.add_argument("--tleap-exe", default="tleap", help="tleap (for the .lib)")
    a.add_argument("--no-tleap", action="store_true", help="skip the library and test peptide")
    a.add_argument("--protein-ff", choices=["ff14SB", "ff19SB"], default="ff14SB",
                   help="protein force field of the tleap test peptide")
    a.add_argument("--backbone-types", choices=["amber", "gaff2"], default="amber",
                   help="atom types of N, H, CA, HA, C, O in the residue library: "
                        "those of the protein force field, or GAFF2")
    a.add_argument("--no-junction-params", dest="junction_params",
                   action="store_false",
                   help="do not generate the missing protein/GAFF2 junction "
                        "parameters with parmchk2")
    args = p.parse_args()
    if args.workers < 1 or args.min_confs < 2 or args.max_confs < args.min_confs:
        p.error("need --workers >= 1 and 2 <= --min-confs <= --max-confs")
    if args.reuse and not os.path.isdir(os.path.join(args.reuse, "runs")):
        p.error(f"--reuse {args.reuse}: no runs/ directory there")
    args.reference = "input"
    if args.am1_opt == "auto":
        args.am1_opt = "restrained" if args.backbone == "fixed" else "full"
    # antechamber: its own minimisation (full) or a single point on the
    # restrained AM1 minimum / the MMFF geometry
    args.ac_ek = args.sqm_keywords if args.am1_opt == "full" else \
        f"{args.sqm_keywords}, maxcyc=0,"
    return args


def split_names(s):
    return [x for x in re.split(r"[,\s]+", s.strip()) if x] if s else []


# ==========================================================================
# main
# ==========================================================================
def main():
    args = parse_args()
    RDLogger.DisableLog("rdApp.*")
    args.ac_exe, args.parmchk2_exe, tleap = base.check_amber(args)

    m2 = base.Mol2(args.mol2)
    names = m2.names
    if len(set(names)) != len(names):
        sys.exit("ERROR: atom names in the mol2 must be unique.")
    prefix = args.prefix or os.path.splitext(os.path.basename(args.mol2))[0]
    out = os.path.abspath(args.outdir or f"{prefix}_{args.method}")
    if args.reuse and os.path.abspath(args.reuse) == out:
        sys.exit("ERROR: --reuse must point to a different directory than -o.")
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
        mol = base.build_rdkit_mol(m2, nc, args.element_from)
        bk = find_backbone(mol, names, split_names(args.ace), split_names(args.nme))
        fixed = fixed_values(bk, names, args)
    except (ValueError, RuntimeError, OSError) as exc:
        say(f"ERROR: {exc}")
        sys.exit(1)
    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
    classes = [int(c) for c in Chem.CanonicalRankAtoms(mol, breakTies=False)]
    X_in = mol.GetConformer().GetPositions()

    frozen = []
    if args.backbone == "fixed":
        frozen = bk["ace"] + bk["nme"] + [i for role in ("N", "H", "CA", "HA", "C", "O")
                                          for i in bk["bb"][role]]
    keep = split_names(args.keep_atoms)
    unknown = [k for k in keep if k not in names]
    if unknown:
        say(f"ERROR: --keep-atoms {', '.join(unknown)} not in the mol2")
        sys.exit(1)
    frozen = sorted(set(frozen) | {names.index(k) for k in keep})
    checker = CovalentChecker(mol, args, frozen, X_in)
    ev_in = checker.evaluate(X_in)
    tors_in = backbone_torsions(bk, X_in)
    residue_atoms = [i for i in range(len(names)) if i not in bk["ace"] + bk["nme"]]
    free = [i for i in range(len(names)) if i not in fixed]

    say(f"Molecule {smiles}: {mol.GetNumAtoms()} atoms, net charge {nc}, "
        f"residue {resname}")
    say(f"Output: {os.path.relpath(out)}")
    say(f"Caps: ACE {' '.join(names[i] for i in bk['ace'])}; "
        f"NME {' '.join(names[i] for i in bk['nme'])}")
    say(f"Backbone: " + ", ".join(f"{r} {' '.join(names[i] for i in bk['bb'][r])}"
                                  for r in ("N", "H", "CA", "HA", "C", "O")
                                  if bk["bb"][r]))
    say(f"Input backbone: phi {tors_in['phi']:.1f}, psi {tors_in['psi']:.1f}, "
        f"omega(ACE) {tors_in['omega_ace']:.1f}, omega(NME) "
        f"{tors_in['omega_nme']:.1f} deg")
    say(f"Fixed charges: {len(fixed)} atoms ({'caps' + (' + N, H, C, O' if args.fix_backbone else '') + (' + ' + args.fixed_charges if args.fixed_charges else '')}); "
        f"{len(free)} free atoms; redistribution {args.redistribute}")
    say(f"Pool: {'caps and backbone fixed' if args.backbone == 'fixed' else 'all atoms free'}"
        f"{' + ' + ', '.join(keep) if keep else ''} ({len(frozen)} frozen atoms)")
    if checker.skip_hb or checker.skip_polar or checker.skip_ion:
        say(f"Contacts inside the frozen part (not counted): "
            f"{len(checker.skip_hb)} H-bond, {len(checker.skip_polar)} N/O, "
            f"{len(checker.skip_ion)} ionic")
    flagged = checker.contacts(X_in, names)
    say(f"Input geometry: Rg {ev_in['rg']:.2f} A, issues: "
        f"{'; '.join(flagged) if flagged else (ev_in['reason'] or 'none')}")
    if not ev_in["ok"]:
        say("  NOTE: the input (bound pose) is still used as the reference conformer.")

    X_ref0 = X_in
    ref_start = os.path.join(out, f"{prefix}_reference_start.mol2")
    m2.write(ref_start, coords=X_ref0, title=f"{m2.mol_name} reference start (input)")

    # ---------------- pool ----------------
    say(f"Building the pool: {args.nconf} ETKDGv3 conformers, {args.ff} "
        f"(eps = {args.dielectric}) ...")
    try:
        rows, pool = build_pool(mol, checker, frozen, X_in, args, say)
    except RuntimeError as exc:
        say(f"ERROR: {exc}")
        sys.exit(1)
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
        m2.write(fn, coords=r["X"], title=f"{m2.mol_name} pool conformer {r['pool_id']}")
        order.append((tag, fn, r))

    # ---------------- antechamber runs ----------------
    say(f"Running {args.method} with {args.workers} worker(s); stop when N >= "
        f"{args.min_confs}, max SEM <= {args.sem_tol} e and max drift over "
        f"{args.drift_window} <= {args.drift_tol} e"
        + (f"; reusing {os.path.relpath(args.reuse)}" if args.reuse else ""))
    accepted, results, Q, history, sem_hist, drift_hist = [], [], [], [], [], []
    n_rejected = 0
    converged = stop = False
    ref = None
    pending, done_by_index = {}, {}
    next_i = next_use = 0
    check_bb = args.backbone == "fixed" and args.bb_tol > 0

    sqm = args.sqm_exe if os.path.isfile(args.sqm_exe) else shutil.which(args.sqm_exe)
    if args.am1_opt == "restrained" and sqm is None:
        say(f"ERROR: {args.sqm_exe} not found (needed for --am1-opt restrained).")
        sys.exit(1)
    ctx = {"sqm": sqm, "z": [a.GetAtomicNum() for a in mol.GetAtoms()],
           "free_opt": [i for i in range(len(names)) if i not in frozen]}
    say(f"AM1 geometry: {args.am1_opt}"
        + (f" ({len(frozen)} atoms fixed, L-BFGS to RMS gradient <= "
           f"{args.am1_grms} and max <= {args.am1_gmax} kcal/mol/A)"
           if args.am1_opt == "restrained" else ""))
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        fut_ref = ex.submit(run_charges, "reference", ref_start,
                            os.path.join(out, "runs", "reference"), nc, resname, args, ctx)
        pending[fut_ref] = "ref"

        def fill():
            nonlocal next_i
            while not stop and next_i < len(order) and len(pending) < args.workers:
                tag, fn, _ = order[next_i]
                pending[ex.submit(run_charges, tag, fn,
                                  os.path.join(out, "runs", tag), nc, resname,
                                  args, ctx)] = next_i
                next_i += 1

        fill()
        while pending:
            done, _ = cf.wait(list(pending), return_when=cf.FIRST_COMPLETED)
            for fut in done:
                key = pending.pop(fut)
                if key == "ref":
                    ref = fut.result()
                    if ref["ok"]:
                        Xr = base.read_sqm_pdb(ref["sqm_pdb"], names)
                        ev = checker.evaluate(Xr)
                        t = backbone_torsions(bk, Xr)
                        # keep the protein frame: superpose on the frozen
                        # atoms (identity after a restrained minimisation)
                        fit = frozen if frozen else checker.heavy
                        ref["X"] = base.kabsch_transform(Xr[fit], X_ref0[fit], Xr)
                        steps = (f", {ref['am1_steps']} AM1 steps"
                                 if ref.get("am1_steps") is not None else "")
                        say(f"  reference: done in {ref['time_s']:.0f} s{steps}; "
                            f"after AM1 phi {t['phi']:.1f}, psi {t['psi']:.1f}, Rg "
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
                     "time_s": res["time_s"], "hof": res.get("hof"),
                     "am1_steps": res.get("am1_steps"), "am1_grms": res.get("am1_grms")}
                results.append(e)
                if not res["ok"]:
                    e.update(accepted=False, reason=res["error"])
                    n_rejected += 1
                    say(f"  {tag}: FAILED ({res['error']})")
                    continue
                Xp = base.read_sqm_pdb(res["sqm_pdb"], names)
                ev = checker.evaluate(Xp)
                t = backbone_torsions(bk, Xp)
                dphi = abs(wrap(t["phi"] - tors_in["phi"]))
                dpsi = abs(wrap(t["psi"] - tors_in["psi"]))
                e.update(post_rg=ev["rg"], post_gauche=base.gauche_count(mol, Xp),
                         post_min_polar=ev["min_polar_dist"], post_phi=t["phi"],
                         post_psi=t["psi"], dphi=dphi, dpsi=dpsi)
                reasons = [ev["reason"]] if ev["reason"] else []
                if check_bb and max(dphi, dpsi) > args.bb_tol:
                    reasons.append(f"backbone moved (dphi {dphi:.0f}, dpsi {dpsi:.0f} deg)")
                note = ("; AM1 minimisation not converged (--am1-maxiter)"
                        if res.get("am1_converged") is False else "")
                if reasons:
                    e.update(accepted=False, reason="; ".join(reasons))
                    n_rejected += 1
                    say(f"  {tag}: rejected after AM1 ({'; '.join(reasons)}){note}")
                    continue
                if converged or len(accepted) >= args.max_confs:
                    e.update(accepted=False, reason="finished before use")
                    continue
                e.update(accepted=True, reason="", types=res["types"],
                         out_mol2=res["out_mol2"])
                accepted.append(e)
                Q.append([res["charges"][n] for n in names])
                mean, std, sem = base.stats(np.array(Q))
                history.append(mean)
                n = len(Q)
                msem = float(sem[free].max()) if n > 1 else float("nan")
                drift = float(np.abs(history[-1] - history[-1 - args.drift_window])[free].max()) \
                    if n > args.drift_window else float("nan")
                sem_hist.append(msem)
                drift_hist.append(drift)
                steps = (f"{res['am1_steps']:4d} AM1 steps, "
                         if res.get("am1_steps") is not None else "")
                say(f"  {tag}: N={n:3d}  {res['time_s']:5.0f} s  max SEM {msem:.4f}  "
                    f"max drift {drift:.4f}  ({steps}Rg {ev['rg']:.2f} A){note}")
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
    mean, std, sem = base.stats(Qa)
    lo, hi = base.bootstrap_ci(Qa, args.nboot, args.seed)
    try:
        q_final, q_sym, info = finalize_with_fixed(mean, nc, classes, fixed, bk,
                                                   args.redistribute)
    except ValueError as exc:
        say(f"ERROR: {exc}")
        sys.exit(1)
    ref_ok = ref is not None and ref["ok"]
    ref_raw = np.array([ref["charges"][n] for n in names]) if ref_ok else None
    ref_q = finalize_with_fixed(ref_raw, nc, classes, fixed, bk, args.redistribute)[0] \
        if ref_ok else None
    X_final = ref["X"] if ref_ok else X_ref0
    tmpl = base.Mol2(accepted[0]["out_mol2"])
    if tmpl.names != names:
        sys.exit("ERROR: antechamber changed the atom order.")
    types = accepted[0]["types"]
    final_mol2 = os.path.join(out, f"{prefix}_{args.method}_ensemble.mol2")
    tmpl.write(final_mol2, coords=X_final, charges=q_final,
               title=f"{m2.mol_name} {args.method} averaged over {len(Q)} conformers, "
                     f"caps fixed")
    if ref_ok:
        tmpl.write(os.path.join(out, f"{prefix}_{args.method}_reference_only.mol2"),
                   coords=X_final, charges=ref_q,
                   title=f"{m2.mol_name} {args.method} on the reference conformer, "
                         f"caps fixed")

    frcmod = os.path.join(out, f"{prefix}.frcmod")
    subprocess.run([args.parmchk2_exe, "-i", final_mol2, "-f", "mol2", "-o", frcmod,
                    "-s", args.atomtype], stdout=subprocess.DEVNULL,
                   stderr=subprocess.STDOUT)

    res_types = [types[n] for n in names]
    if args.backbone_types == "amber":
        for role, t in BB_TYPES[args.protein_ff].items():
            for i in bk["bb"][role]:
                res_types[i] = t
    res_mol2 = os.path.join(out, f"{prefix}_{args.method}_residue.mol2")
    write_residue_mol2(res_mol2, m2, residue_atoms, X_final, q_final, res_types,
                       resname, f"{resname} {args.method}, caps removed")
    res_charge = float(q_final[residue_atoms].sum())

    leap_note, leap_first, junction = "skipped", None, None
    if tleap:
        ldir = os.path.join(out, "tleap")
        leap_first, leap_final, junction = build_library(
            args, tleap, out, prefix, resname, res_mol2, frcmod, ldir, say)
        leap_note = (f"residue charge {leap_final['charge']}, errors "
                     f"{leap_final['errors']}, warnings {leap_final['warnings']}, "
                     f"missing parameters {leap_final['missing']}")

    # ---------------- tables ----------------
    role_of = {}
    for i in range(len(names)):
        if i in bk["roles"]:
            role_of[i] = bk["roles"][i][0].lower()
        elif any(i in v for v in bk["bb"].values()):
            role_of[i] = "backbone"
        else:
            role_of[i] = "residue"
    with open(os.path.join(out, f"{prefix}_charges_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["atom", "type", "sym_class", "q_final", "mean", "std", "sem",
                    "ci95_lo", "ci95_hi", "min", "max", "q_reference",
                    "final_minus_reference", "q_input_file", "role", "fixed",
                    "residue_type"])
        for i, n in enumerate(names):
            qin = m2.atoms[i]["charge"]
            w.writerow([n, types[n], classes[i], f"{q_final[i]:.6f}", f"{mean[i]:.6f}",
                        f"{std[i]:.6f}", f"{sem[i]:.6f}", f"{lo[i]:.6f}",
                        f"{hi[i]:.6f}", f"{Qa[:, i].min():.6f}", f"{Qa[:, i].max():.6f}",
                        f"{ref_q[i]:.6f}" if ref_ok else "",
                        f"{q_final[i] - ref_q[i]:.6f}" if ref_ok else "",
                        f"{qin:.6f}" if (qin is not None and m2.has_charges()) else "",
                        role_of[i], int(i in fixed),
                        res_types[i] if i in residue_atoms else ""])
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
                    "post_phi_deg", "post_psi_deg", "dphi_deg", "dpsi_deg",
                    "am1_steps", "am1_grms", "am1_hof_kcal", "time_s",
                    "accepted", "reason"])
        for e in results:
            def fmt(k, f="{:.3f}"):
                return f.format(e[k]) if e.get(k) is not None and k in e else ""
            w.writerow([e["tag"], e["pool_id"], f"{e['mmff_dE']:.3f}",
                        f"{e['pre_rg']:.3f}", e["pre_gauche"], fmt("post_rg"),
                        e.get("post_gauche", ""), fmt("post_min_polar"),
                        fmt("post_phi", "{:.1f}"), fmt("post_psi", "{:.1f}"),
                        fmt("dphi", "{:.1f}"), fmt("dpsi", "{:.1f}"),
                        e.get("am1_steps") if e.get("am1_steps") is not None else "",
                        fmt("am1_grms", "{:.4f}"),
                        fmt("hof", "{:.4f}"), f"{e['time_s']:.1f}",
                        int(e["accepted"]), e["reason"]])
    with open(os.path.join(out, f"{prefix}_accepted_sqm.pdb"), "w") as fh:
        for k, e in enumerate(accepted, 1):
            fh.write(f"MODEL {k:8d}\n")
            with open(os.path.join(out, "runs", e["tag"], "sqm.pdb")) as src:
                fh.writelines(l for l in src if l.startswith(("ATOM", "HETATM")))
            fh.write("ENDMDL\n")

    sem_free = sem[free]
    std_free = std[free]
    summary = {
        "input": args.mol2, "prefix": prefix, "smiles": smiles,
        "resname": resname, "net_charge": nc, "method": args.method,
        "atomtype": args.atomtype,
        "input_issues": "; ".join(flagged) if flagged else ev_in["reason"],
        "input_amides": checker.amide_report(X_in, names),
        "reference_geometry": "input, bound pose",
        "pool_embedded": len(rows), "pool_accepted": len(pool),
        "runs": len(results), "accepted": len(accepted), "rejected": n_rejected,
        "converged": converged, "final_max_sem": float(sem_free.max()),
        "final_max_drift": drift_hist[-1], "max_std": float(std_free.max()),
        "max_std_atom": names[free[int(std_free.argmax())]],
        "max_abs_final_minus_reference":
            float(np.abs(q_final - ref_q)[free].max()) if ref_ok else None,
        "tleap": leap_note,
        "covalent": {
            "ace": [names[i] for i in bk["ace"]], "nme": [names[i] for i in bk["nme"]],
            "backbone": {r: [names[i] for i in v] for r, v in bk["bb"].items()},
            "fixed": {names[i]: v for i, v in sorted(fixed.items())},
            "redistribute": args.redistribute, "backbone_sampling": args.backbone,
            "frozen_in_pool": [names[i] for i in frozen],
            "input_torsions": tors_in, **info,
            "residue_atoms": [names[i] for i in residue_atoms],
            "residue_charge": round(res_charge, 6),
            "residue_mol2": os.path.basename(res_mol2),
            "junction_frcmod": os.path.basename(junction) if junction else None,
            "tleap_first_pass": leap_first},
        "files": {"mol2": os.path.basename(final_mol2),
                  "residue_mol2": os.path.basename(res_mol2),
                  "frcmod": os.path.basename(frcmod)},
        "settings": {k: v for k, v in vars(args).items()},
    }
    with open(os.path.join(out, f"{prefix}_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    say(f"\nAccepted {len(accepted)} of {len(results)} runs ({n_rejected} rejected); "
        f"converged: {converged}")
    say(f"Final max SEM {sem_free.max():.4f} e; largest conformational SD "
        f"{std_free.max():.4f} e ({names[free[int(std_free.argmax())]]})")
    say(f"ABCG2 charge of the caps: ACE {info['ace_abcg2_sum']:+.4f} e, NME "
        f"{info['nme_abcg2_sum']:+.4f} e (fixed to {info['ace_fixed_sum']:+.4f} and "
        f"{info['nme_fixed_sum']:+.4f} e); "
        + (f"shift per free atom {info['shift_per_free_atom']:+.5f} e"
           if args.redistribute == "even" else "moved onto N and C"))
    if ref_ok:
        d = np.abs(q_final - ref_q)[free]
        say(f"Largest |ensemble - reference conformer| {d.max():.4f} e "
            f"({names[free[int(d.argmax())]]})")
    say(f"Residue {resname}: {len(residue_atoms)} atoms, charge {res_charge:+.6f} e")
    if leap_first is not None:
        say(f"tleap, ACE-{resname}-NME: {leap_first['missing']} missing parameters "
            f"before the junction frcmod; final check: {leap_note}")
    rel = os.path.relpath
    say(f"Charges: {rel(final_mol2)}\nResidue: {rel(res_mol2)}\nfrcmod:  {rel(frcmod)}")
    if junction:
        say(f"Junction: {rel(junction)}")
    if tleap:
        say(f"Library: {rel(os.path.join(out, 'tleap', prefix + '.lib'))}")
    say(f"Analyse with: python3 abcg2_analysis.py {rel(out)}")
    log.close()


if __name__ == "__main__":
    main()
