#!/usr/bin/env python3
"""
abcg2_fragment_assembly.py

Residue charges for an oligomer or polymer from capped fragments whose
conformer-averaged ABCG2 charges were computed with abcg2_conformer_charges.py.
Each fragment is one residue of the chain with caps that stand in for its
neighbours; the caps are removed, every residue is given an exact integer
charge, and the residues are joined into the chain with tleap.

    python3 abcg2_fragment_assembly.py PET_residues.json -o PET_residues \
        --whole PET5.mol2

The definition file (JSON; paths relative to it):

    {
      "chain": ["PTB", "PTM", "PTM", "PTM", "PTE"],
      "residues": {
        "PTM": {"run": "PET_mid_abcg2", "head": "C1", "tail": "C10", "charge": 0,
                "mimic": {"CP1": "C10", "CN1": "C1"}},
        "PTB": {"run": "PET_end_abcg2", "tail": "C10", "charge": 0},
        "PTE": {"run": "PET_end_abcg2", "template": "PET_endrev.mol2",
                "head": "C1", "charge": 0}
      }
    }

  run       output directory of abcg2_conformer_charges.py for the fragment;
            the residue is made of the atoms whose residue name in the input
            mol2 of the run (the 4th-last column) is the key, the other atoms
            are caps
  template  optional: the same molecule as the run with other atom and residue
            names (matched atom by atom through the bond graph), e.g. to take
            the last unit of a chain from the fragment of the first one
  head/tail atoms bonded to the previous / next residue
  charge    integer charge of the residue
  mimic     optional: cap atom -> residue atom it stands for, compared in the
            report (a cap that reproduces the chain has the same charges)

Steps
  1. Charges of each fragment: the final ensemble charges of the run
     (q_final: mean over conformers, averaged over symmetry-equivalent atoms)
     and, for comparison, those of its single reference conformer.
  2. The caps are removed. The charge they carried (the difference between
     the residue charge and its integer value) is spread evenly over the
     residue (--redistribute even) or returned to the residue atom bonded to
     each cap (--redistribute junction). Rounded to 1e-6 e, exact integer.
  3. Residue mol2 files (GAFF2 types from the run), an Amber library with
     head and tail atoms, and the chain built with tleap (sequence), with
     its total charge checked.
  4. With --whole, ABCG2 charges of the whole chain in one extended
     conformer (antechamber with the standard AM1 minimisation; reused if
     already present), compared atom by atom with the assembled charges:
     a figure (600 dpi) and a report, including the net charge that each
     residue carries in the whole molecule. --whole-dir reuses an existing
     whole-chain calculation, and --subset anti|gauche builds the residues
     from only those fragment conformers whose CH2-CH2 torsions are all anti
     or all gauche.

Requirements: abcg2_conformer_charges.py in the same directory, AmberTools
(antechamber, sqm, parmchk2, tleap), Python 3.8+, RDKit, NumPy, Matplotlib.
"""

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import abcg2_conformer_charges as base
except ImportError:
    sys.exit("ERROR: abcg2_conformer_charges.py must be in the same directory.")

from rdkit import Chem, RDLogger

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]
TEXT_1, TEXT_2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


# ==========================================================================
# fragment runs
# ==========================================================================
def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def load_run(d):
    js = [f for f in os.listdir(d) if f.endswith("_summary.json")]
    if not js:
        sys.exit(f"ERROR: {d} has no *_summary.json (not an output directory of "
                 f"abcg2_conformer_charges.py?)")
    prefix = js[0][:-len("_summary.json")]
    with open(os.path.join(d, js[0])) as fh:
        summary = json.load(fh)
    method = summary["method"]
    rows = read_csv(os.path.join(d, f"{prefix}_charges_summary.csv"))
    run = {"dir": d, "prefix": prefix, "summary": summary,
           "input": base.Mol2(os.path.join(d, f"{prefix}_input.mol2")),
           "final": base.Mol2(os.path.join(d, f"{prefix}_{method}_ensemble.mol2")),
           "frcmod": os.path.join(d, f"{prefix}.frcmod"),
           "rows": {r["atom"]: r for r in rows}}
    if run["input"].names != run["final"].names:
        sys.exit(f"ERROR: atom order differs between the input and output of {d}")
    return run


def subset_charges(run, subset):
    """Mean charges over the accepted conformers whose O-C-C-O type torsions
    ([!#1]-CH2-CH2-[!#1], as counted by abcg2_conformer_charges.py) are all
    anti (subset 'anti') or all gauche ('gauche'), averaged over
    symmetry-equivalent atoms and corrected to the net charge of the fragment."""
    d, prefix = run["dir"], run["prefix"]
    per = read_csv(os.path.join(d, f"{prefix}_charges_per_conformer.csv"))
    tags = [k for k in per[0].keys() if k != "atom"]
    gauche = {r["tag"]: int(r["post_gauche"]) for r in
              read_csv(os.path.join(d, f"{prefix}_conformer_runs.csv"))
              if r["accepted"] == "1"}
    mol = rdkit_of(run["input"], int(run["summary"]["net_charge"]))
    n_tors = len({tuple(sorted(t[1:3])) for t in mol.GetSubstructMatches(base.GAUCHE_SMARTS)})
    want = 0 if subset == "anti" else n_tors
    sel = [t for t in tags if gauche.get(t) == want]
    if not sel:
        sys.exit(f"ERROR: no conformer of {d} has {want} gauche torsions (--subset {subset})")
    qa = {r["atom"]: np.array([float(r[t]) for t in sel]) for r in per}
    names = run["input"].names
    q = np.array([qa[n].mean() for n in names])
    cls = [run["rows"][n]["sym_class"] for n in names]
    for c in set(cls):
        idx = [i for i, x in enumerate(cls) if x == c]
        q[idx] = q[idx].mean()
    q += (int(run["summary"]["net_charge"]) - q.sum()) / len(q)
    return q, len(sel), n_tors


def rdkit_of(m2, charge):
    return base.build_rdkit_mol(m2, charge, "auto")


def per_atom(run, key):
    return np.array([float(run["rows"][n][key]) if run["rows"][n].get(key) not in
                     (None, "") else float("nan") for n in run["input"].names])


# ==========================================================================
# residues
# ==========================================================================
def extract(name, spec, root, run, mode, subset="all"):
    """Charges, types and coordinates of one residue from its fragment."""
    nc = int(run["summary"]["net_charge"])
    q_ens = per_atom(run, "q_final")
    n_used, n_tors = int(run["summary"]["accepted"]), None
    if subset != "all":
        q_ens, n_used, n_tors = subset_charges(run, subset)
    q_ref = per_atom(run, "q_reference")
    sd = per_atom(run, "std")
    src = run["input"]
    if spec.get("template"):
        tmpl = base.Mol2(os.path.join(root, spec["template"]))
        m_run, m_tmp = rdkit_of(src, nc), rdkit_of(tmpl, nc)
        match = m_run.GetSubstructMatch(m_tmp)
        if len(match) != m_tmp.GetNumAtoms() or m_run.GetNumAtoms() != m_tmp.GetNumAtoms():
            sys.exit(f"ERROR: template {spec['template']} is not the molecule of "
                     f"{run['dir']}")
        order = list(match)                      # template atom -> run atom
        frag = tmpl
    else:
        order = list(range(len(src.atoms)))
        frag = src
    resname_of = [a["subst_name"] for a in frag.atoms]
    names = frag.names
    res_idx = [k for k, r in enumerate(resname_of) if r == name]
    if not res_idx:
        sys.exit(f"ERROR: no atom of residue {name} in "
                 f"{spec.get('template') or run['dir']} (residue column of the mol2)")
    cap_idx = [k for k in range(len(names)) if k not in res_idx]
    m_frag = rdkit_of(frag, nc)
    res_set = set(res_idx)
    # caps: connected groups of non-residue atoms and the residue atom each joins
    caps, seen = [], set()
    for k in cap_idx:
        if k in seen:
            continue
        comp, stack = [], [k]
        while stack:
            i = stack.pop()
            if i in seen:
                continue
            seen.add(i)
            comp.append(i)
            stack += [nb.GetIdx() for nb in m_frag.GetAtomWithIdx(i).GetNeighbors()
                      if nb.GetIdx() not in res_set and nb.GetIdx() not in seen]
        joint = sorted({nb.GetIdx() for i in comp for nb in m_frag.GetAtomWithIdx(i).GetNeighbors()
                        if nb.GetIdx() in res_set})
        caps.append({"atoms": sorted(comp), "joint": joint})

    def build(q_all):
        q = np.array([q_all[order[k]] for k in range(len(names))])
        target = int(spec.get("charge", 0))
        qr = q[res_idx].copy()
        delta = target - qr.sum()
        info = {"delta": float(delta),
                "caps": [{"atoms": [names[i] for i in c["atoms"]],
                          "joint": [names[i] for i in c["joint"]],
                          "charge": float(q[c["atoms"]].sum())} for c in caps]}
        if mode == "junction":
            pos = {i: p for p, i in enumerate(res_idx)}
            for c in caps:
                if len(c["joint"]) != 1:
                    sys.exit(f"ERROR: a cap of {name} is bonded to {len(c['joint'])} "
                             f"residue atoms; use --redistribute even")
                qr[pos[c["joint"][0]]] += q[c["atoms"]].sum()
            qr += (target - qr.sum()) / len(qr)
        else:
            qr += delta / len(qr)
        qr = np.round(qr, 6)
        resid = round(target - qr.sum(), 6)
        if abs(resid) > 0:
            qr[int(np.argmax(np.abs(qr)))] += resid
            qr = np.round(qr, 6)
        return qr, q, info

    q_e, full_e, info_e = build(q_ens)
    q_r, full_r, info_r = build(q_ref)
    mimic = []
    for capname, resatom in (spec.get("mimic") or {}).items():
        if capname not in names or resatom not in names:
            sys.exit(f"ERROR: mimic {capname} -> {resatom}: atom not in the fragment")
        mimic.append((capname, resatom, float(full_e[names.index(capname)]),
                      float(full_e[names.index(resatom)])))
    types = [run["final"].atoms[order[k]]["type"] for k in res_idx]
    xyz = [run["final"].atoms[order[k]]["xyz"] for k in res_idx]
    bonds = [(i, j) for i, j, _ in frag.bonds if i in res_set and j in res_set]
    return {"name": name, "spec": spec, "run": run, "names": [names[k] for k in res_idx],
            "types": types, "xyz": xyz, "q": q_e, "q_ref": q_r,
            "sd": np.array([sd[order[k]] for k in res_idx]),
            "raw": np.array([q_ens[order[k]] for k in res_idx]),
            "info": info_e, "info_ref": info_r, "mimic": mimic,
            "n_used": n_used, "n_tors": n_tors,
            "idx": {k: p for p, k in enumerate(res_idx)}, "bonds": bonds,
            "frag_bonds": frag.bonds}


def write_residue_mol2(path, R, charges):
    new = {k: p + 1 for k, p in R["idx"].items()}
    bonds = [(new[i], new[j]) for i, j in R["bonds"]]
    btype = {}
    for i, j, t in R["frag_bonds"]:
        if i in new and j in new:
            btype[(new[i], new[j])] = t
    L = ["@<TRIPOS>MOLECULE", R["name"],
         f"{len(R['names']):5d} {len(bonds):5d}     1     0     0", "SMALL",
         "USER_CHARGES", "", "", "@<TRIPOS>ATOM"]
    for p, (n, t, x, q) in enumerate(zip(R["names"], R["types"], R["xyz"], charges), 1):
        L.append(f"{p:7d} {n:<8s} {x[0]:10.4f} {x[1]:10.4f} {x[2]:10.4f} {t:<6s} "
                 f"{1:4d} {R['name']:<8s} {q:10.6f}")
    L.append("@<TRIPOS>BOND")
    for k, (i, j) in enumerate(bonds, 1):
        L.append(f"{k:6d} {i:5d} {j:5d} {btype.get((i, j), '1')}")
    L += ["@<TRIPOS>SUBSTRUCTURE",
          f"     1 {R['name']:<8s}    1 RESIDUE           4 A     {R['name']:<4s}  1 ROOT", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(L))


# ==========================================================================
# tleap
# ==========================================================================
def tleap_build(args, out, residues, chain, frcmods, lib):
    ldir = os.path.join(out, "tleap")
    os.makedirs(ldir, exist_ok=True)
    L = [f"source {args.leaprc}"]
    for f in frcmods:
        L.append(f"loadamberparams ../{os.path.basename(f)}")
    for name, R in residues.items():
        L.append(f"{name} = loadmol2 ../{name}.mol2")
        if R["spec"].get("head"):
            L += [f"set {name} head {name}.1.{R['spec']['head']}",
                  f"set {name}.1 connect0 {name}.1.{R['spec']['head']}"]
        if R["spec"].get("tail"):
            L += [f"set {name} tail {name}.1.{R['spec']['tail']}",
                  f"set {name}.1 connect1 {name}.1.{R['spec']['tail']}"]
        L += [f"charge {name}", f"check {name}"]
    for name in residues:
        L.append(f"saveoff {name} {lib}")
    L += [f"chain = sequence {{ {' '.join(chain)} }}", "charge chain", "check chain",
          "saveamberparm chain chain.prmtop chain.inpcrd", "savepdb chain chain.pdb",
          "savemol2 chain chain.mol2 1", "quit", ""]
    with open(os.path.join(ldir, "leap.in"), "w") as fh:
        fh.write("\n".join(l for l in L if l is not None))
    with open(os.path.join(ldir, "leap.log"), "w") as lg:
        subprocess.run([args.tleap_exe, "-f", "leap.in"], cwd=ldir, stdout=lg,
                       stderr=subprocess.STDOUT)
    txt = open(os.path.join(ldir, "leap.log")).read()
    charges = re.findall(r"Total unperturbed charge:\s*(-?\d+\.\d+)", txt)
    me = re.findall(r"Errors = (\d+); Warnings = (\d+)", txt)
    missing = len(re.findall(r"Could not find (?:bond|angle) parameter|No torsion terms", txt))
    return {"chain_charge": charges[-1] if charges else "?",
            "errors": me[0][0] if me else "?", "warnings": me[0][1] if me else "?",
            "missing": missing}


# ==========================================================================
# whole molecule
# ==========================================================================
def whole_charges(args, wdir, whole_mol2):
    """ABCG2 of the whole molecule in one conformer, reused if present."""
    m2 = base.Mol2(whole_mol2)
    nc = (m2.net_charge_from_file() or 0) if args.whole_charge is None else args.whole_charge
    outf = os.path.join(wdir, f"whole_{args.method}.mol2")
    reuse = os.path.isfile(outf) and os.path.isfile(os.path.join(wdir, "sqm.pdb")) and \
        os.path.isfile(os.path.join(wdir, "input.mol2"))
    if reuse:
        old = base.Mol2(os.path.join(wdir, "input.mol2"))
        reuse = old.names == m2.names and all(
            np.allclose(a["xyz"], b["xyz"], atol=1e-4) for a, b in zip(old.atoms, m2.atoms))
    if reuse:
        print(f"Whole molecule: reusing {os.path.relpath(outf)}")
    else:
        print(f"Whole molecule: antechamber -c {args.method} on {len(m2.atoms)} atoms "
              f"(full AM1 minimisation, this can take a while) ...", flush=True)
        ns = argparse.Namespace(ac_exe=args.ac_exe, method=args.method, atomtype="gaff2",
                                ac_ek=None, keep_scratch=False)
        res = base.run_antechamber("whole", whole_mol2, wdir, nc, "MOL", ns)
        if not res["ok"]:
            sys.exit(f"ERROR: antechamber failed for the whole molecule (see {wdir})")
        print(f"  done in {res['time_s']:.0f} s")
    q = np.array([a["charge"] for a in base.Mol2(outf).atoms])
    X = base.read_sqm_pdb(os.path.join(wdir, "sqm.pdb"), m2.names) \
        if len(set(m2.names)) == len(m2.names) else None
    if X is None:
        X = read_pdb_by_order(os.path.join(wdir, "sqm.pdb"))
    return m2, q, X


def read_pdb_by_order(path):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(path) if l.startswith(("ATOM", "HETATM"))])


def gauche_states(m2, X, charge):
    mol = rdkit_of(m2, charge)
    out = []
    for t in mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2]-[CH2]-[CH2]-[OX2]")):
        out.append(base.dihedral(*(X[i] for i in t)))
    return out


# ==========================================================================
# figure and report
# ==========================================================================
def style(mpl):
    mpl.rcParams.update({
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "axes.edgecolor": TEXT_2, "axes.labelcolor": TEXT_1,
        "xtick.color": TEXT_2, "ytick.color": TEXT_2,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "legend.frameon": False, "legend.fontsize": 7,
        "lines.linewidth": 1.6, "figure.dpi": 110})


def figure(cmp, chain, plt, label_min=0.01):
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.6), constrained_layout=True,
                             gridspec_kw={"height_ratios": [2.2, 1]})
    ax = axes[0]
    x = np.arange(len(cmp["atoms"]))
    ax.axhline(0, color=TEXT_2, lw=0.8)
    ax.plot(x, 1e3 * cmp["d_ens"], "o", ms=2.8, color=SERIES[0],
            label="ensemble residue charges - whole chain")
    ax.plot(x, 1e3 * cmp["d_ref"], "s", ms=2.6, color=SERIES[1],
            label="reference-conformer residue charges - whole chain")
    starts = [k for k in range(len(x)) if k == 0 or cmp["resid"][k] != cmp["resid"][k - 1]]
    for k in starts[1:]:
        ax.axvline(k - 0.5, color=GRID, lw=1.0)
    for r, k in enumerate(starts):
        end = starts[r + 1] if r + 1 < len(starts) else len(x)
        ax.text((k + end - 1) / 2, 0.02, f"{r + 1} {chain[r]}", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=7, color=TEXT_2)
    groups = {}
    for k in x:
        d = max((cmp["d_ens"][k], cmp["d_ref"][k]), key=abs)
        if abs(d) >= label_min:
            key = (cmp["resid"][k], round(1e3 * d, 1))
            groups.setdefault(key, []).append(k)
    for (r, dy), ks in groups.items():
        # equivalent atoms with the same difference share one label
        ax.annotate(", ".join(cmp["atoms"][k] for k in ks), (ks[-1], dy), xytext=(3, 0),
                    va="center", textcoords="offset points", fontsize=6, color=TEXT_1)
    ax.set_xlim(-1, len(x) + 2)
    ax.set_xticks([])
    ax.set_ylabel("q(assembled) - q(whole) (me)")
    ax.set_xlabel("Atoms along the chain (labels: |difference| >= "
                  f"{1e3 * label_min:g} me)")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside upper center", ncol=2)
    ax2 = axes[1]
    pos = np.arange(1, len(chain) + 1)
    ax2.bar(pos, 1e3 * np.array(cmp["res_net_whole"]), width=0.55, color=SERIES[2])
    ax2.axhline(0, color=TEXT_2, lw=0.8)
    ax2.set_xticks(pos)
    ax2.set_xticklabels([f"{k} {n}" for k, n in zip(pos, chain)])
    ax2.set_ylabel("Net charge (me)")
    ax2.set_title("Net charge of each residue in the whole chain (assembled residues: 0)",
                  loc="left")
    return fig


def report(path, args, residues, chain, leap, cmp, gauche):
    """Markdown report of the residues, the chain and the comparison."""
    L = ["# Fragment-based ABCG2 residue charges", ""]
    L.append(f"Redistribution of the cap charge: {args.redistribute}; fragment "
             f"conformers: {args.subset}")
    L.append("")
    for name, R in residues.items():
        s = R["run"]["summary"]
        L += [f"## {name}", "",
              f"- Fragment: {os.path.relpath(R['run']['dir'], os.path.dirname(path))} "
              f"(`{s['smiles']}`), {R['n_used']} conformers"
              + (f" ({args.subset} subset)" if args.subset != "all" else "")
              + ", converged: "
              f"{'yes' if s['converged'] else 'no'}"
              + (f"; atoms named after {R['spec']['template']}" if R['spec'].get('template') else ""),
              f"- Residue: {len(R['names'])} atoms, charge {int(R['spec'].get('charge', 0))}"
              f", head {R['spec'].get('head', '-')}, tail {R['spec'].get('tail', '-')}"]
        for c in R["info"]["caps"]:
            L.append(f"- Cap {' '.join(c['atoms'])} (bonded to {', '.join(c['joint'])}): "
                     f"{c['charge']:+.4f} e")
        L.append(f"- Charge moved onto the residue: {R['info']['delta']:+.4f} e "
                 f"(ensemble), {R['info_ref']['delta']:+.4f} e (reference conformer)")
        if R["mimic"]:
            L.append("- Cap atoms against the residue atoms they stand for (ensemble, e): "
                     + "; ".join(f"{c} {qc:+.4f} / {r} {qr:+.4f}" for c, r, qc, qr in R["mimic"]))
        L += ["", "| atom | type | q (e) | fragment SD (e) | reference conformer (e) |",
              "|---|---|---|---|---|"]
        for n, t, q, sd, qr in zip(R["names"], R["types"], R["q"], R["sd"], R["q_ref"]):
            L.append(f"| {n} | {t} | {q:+.4f} | {sd:.4f} | {qr:+.4f} |")
        L.append("")
    L += ["## Chain", "", f"- Sequence: {' '.join(chain)}",
          f"- tleap: total charge {leap['chain_charge']}, errors {leap['errors']}, "
          f"warnings {leap['warnings']}, missing parameters {leap['missing']}", ""]
    if cmp:
        L += ["## Comparison with the whole chain", "",
              f"- Whole chain: {cmp['n']} atoms, ABCG2 on one conformer after the "
              f"standard AM1 minimisation",
              f"- O-C-C-O torsions of the whole chain after AM1: "
              + ", ".join(f"{t:.0f}" for t in gauche["whole"]) + " deg",
              f"- Residue charges from the reference conformers: RMS "
              f"{cmp['rms_ref']:.4f} e, max {cmp['max_ref']:.4f} e ({cmp['max_ref_atom']})",
              f"- Residue charges from the ensembles: RMS {cmp['rms_ens']:.4f} e, max "
              f"{cmp['max_ens']:.4f} e ({cmp['max_ens_atom']})",
              "- Net charge of each residue in the whole chain: "
              + ", ".join(f"{k + 1} {n} {q:+.4f}" for k, (n, q) in
                          enumerate(zip(chain, cmp["res_net_whole"]))), "",
              "| residue | atoms | RMS ref (e) | max ref (e) | RMS ensemble (e) | max ensemble (e) |",
              "|---|---|---|---|---|---|"]
        for k, row in enumerate(cmp["per_res"]):
            L.append(f"| {k + 1} {chain[k]} | {row['n']} | {row['rms_ref']:.4f} | "
                     f"{row['max_ref']:.4f} ({row['max_ref_atom']}) | {row['rms_ens']:.4f} | "
                     f"{row['max_ens']:.4f} ({row['max_ens_atom']}) |")
        L.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(L))


# ==========================================================================
# main
# ==========================================================================
def main():
    p = argparse.ArgumentParser(
        description="Residue charges of an oligomer from conformer-averaged ABCG2 "
                    "charges of capped fragments (abcg2_conformer_charges.py runs).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("definition", help="JSON file with the residues and the chain")
    p.add_argument("-o", "--outdir", default=None,
                   help="output directory (None = definition basename)")
    p.add_argument("--redistribute", choices=["even", "junction"], default="even",
                   help="where the charge carried by the caps goes: spread over the "
                        "residue, or onto the residue atom bonded to each cap")
    p.add_argument("--subset", choices=["all", "anti", "gauche"], default="all",
                   help="fragment conformers used for the residue charges: all "
                        "(the final ensemble charges), or only those whose "
                        "[!#1]-CH2-CH2-[!#1] torsions are all anti or all gauche "
                        "after AM1 (to compare with a whole chain in that state)")
    p.add_argument("--whole", default=None, metavar="MOL2",
                   help="mol2 of the whole chain, with the residue names and numbers "
                        "of the chain: ABCG2 on it and comparison")
    p.add_argument("--whole-dir", default=None,
                   help="directory of the whole-chain calculation (None = "
                        "<outdir>/whole); reused when it holds the same input")
    p.add_argument("--whole-charge", type=int, default=None,
                   help="net charge of the whole chain (None = from the mol2, or 0)")
    p.add_argument("--method", default="abcg2", help="antechamber charge method")
    p.add_argument("--ac-exe", default="antechamber", help="antechamber")
    p.add_argument("--tleap-exe", default="tleap", help="tleap")
    p.add_argument("--leaprc", default="leaprc.gaff2", help="leaprc for the residues")
    p.add_argument("--lib", default=None, help="library name (None = <outdir>.lib)")
    p.add_argument("--label-min", type=float, default=0.01,
                   help="label atoms whose charge differs from the whole chain by at "
                        "least this much (e)")
    p.add_argument("--dpi", type=int, default=600, help="figure resolution")
    p.add_argument("--no-show", action="store_true", help="do not open the figure")
    p.add_argument("--overwrite", action="store_true",
                   help="allow writing into an existing output directory")
    args = p.parse_args()
    RDLogger.DisableLog("rdApp.*")

    root = os.path.dirname(os.path.abspath(args.definition))
    with open(args.definition) as fh:
        D = json.load(fh)
    chain = D.get("chain", [])
    if not D.get("residues"):
        sys.exit("ERROR: the definition has no residues")
    for r in chain:
        if r not in D["residues"]:
            sys.exit(f"ERROR: residue {r} of the chain is not defined")
    out = os.path.abspath(args.outdir or os.path.splitext(os.path.basename(args.definition))[0])
    if os.path.isdir(out) and [f for f in os.listdir(out) if f != "whole"] and not args.overwrite:
        sys.exit(f"ERROR: {out} exists and is not empty (use --overwrite or -o).")
    os.makedirs(out, exist_ok=True)
    tleap = args.tleap_exe if os.path.isfile(args.tleap_exe) else shutil.which(args.tleap_exe)

    runs, residues = {}, {}
    for name, spec in D["residues"].items():
        rd = os.path.normpath(os.path.join(root, spec["run"]))
        if rd not in runs:
            runs[rd] = load_run(rd)
        residues[name] = extract(name, spec, root, runs[rd], args.redistribute, args.subset)
        R = residues[name]
        caps = "; ".join(f"{'+'.join(c['joint'])} cap {c['charge']:+.4f}"
                         for c in R["info"]["caps"])
        used = f"{R['n_used']} conformers" + (
            f", {args.subset} subset, {R['n_tors']} O-C-C-O torsions" if args.subset != "all" else "")
        print(f"{name}: {len(R['names'])} atoms from {os.path.relpath(rd)} "
              f"({used}); {caps}; "
              f"{R['info']['delta']:+.4f} e moved onto the residue")
        write_residue_mol2(os.path.join(out, f"{name}.mol2"), R, R["q"])

    frcmods = []
    for rd, run in runs.items():
        if os.path.isfile(run["frcmod"]):
            dst = os.path.join(out, f"{run['prefix']}.frcmod")
            shutil.copy(run["frcmod"], dst)
            frcmods.append(dst)
    with open(os.path.join(out, "residue_charges.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["residue", "atom", "type", "q", "q_reference_conformer",
                    "fragment_std", "fragment_q_final"])
        for name, R in residues.items():
            for n, t, q, qr, sd, raw in zip(R["names"], R["types"], R["q"], R["q_ref"],
                                            R["sd"], R["raw"]):
                w.writerow([name, n, t, f"{q:.6f}", f"{qr:.6f}", f"{sd:.6f}", f"{raw:.6f}"])

    lib = args.lib or f"{os.path.basename(out)}.lib"
    leap = {"chain_charge": "skipped", "errors": "-", "warnings": "-", "missing": "-"}
    if tleap and chain:
        leap = tleap_build(args, out, residues, chain, frcmods, lib)
        print(f"tleap: {' '.join(chain)}, total charge {leap['chain_charge']}, errors "
              f"{leap['errors']}, missing parameters {leap['missing']}")

    cmp, gauche = None, {}
    if args.whole:
        wdir = os.path.abspath(args.whole_dir or os.path.join(out, "whole"))
        os.makedirs(wdir, exist_ok=True)
        m2, qw, Xw = whole_charges(args, wdir, args.whole)
        nc_w = (m2.net_charge_from_file() or 0) if args.whole_charge is None \
            else args.whole_charge
        gauche["whole"] = gauche_states(m2, Xw, nc_w)
        resids = []
        for a in m2.atoms:
            if (a["subst_id"], a["subst_name"]) not in resids:
                resids.append((a["subst_id"], a["subst_name"]))
        if [rn for _, rn in resids] != chain:
            sys.exit(f"ERROR: residues of {args.whole} ({' '.join(rn for _, rn in resids)}) "
                     f"differ from the chain")
        atoms, rid, d_e, d_r, qws = [], [], [], [], []
        res_net = []
        per_res = []
        for k, (sid, rn) in enumerate(resids):
            R = residues[rn]
            pos = {n: p for p, n in enumerate(R["names"])}
            idx = [i for i, a in enumerate(m2.atoms) if a["subst_id"] == sid]
            if sorted(m2.atoms[i]["name"] for i in idx) != sorted(R["names"]):
                sys.exit(f"ERROR: atoms of residue {k + 1} {rn} in {args.whole} differ "
                         f"from the residue")
            res_net.append(float(qw[idx].sum()))
            de, dr, an = [], [], []
            for i in idx:
                n = m2.atoms[i]["name"]
                de.append(R["q"][pos[n]] - qw[i])
                dr.append(R["q_ref"][pos[n]] - qw[i])
                an.append(n)
            de, dr = np.array(de), np.array(dr)
            ke, kr = int(np.abs(de).argmax()), int(np.abs(dr).argmax())
            per_res.append({"n": len(idx), "rms_ens": math.sqrt((de ** 2).mean()),
                            "max_ens": abs(de[ke]), "max_ens_atom": an[ke],
                            "rms_ref": math.sqrt((dr ** 2).mean()),
                            "max_ref": abs(dr[kr]), "max_ref_atom": an[kr]})
            atoms += an
            rid += [k] * len(idx)
            d_e += list(de)
            d_r += list(dr)
            qws += [qw[i] for i in idx]
        d_e, d_r = np.array(d_e), np.array(d_r)
        ke, kr = int(np.abs(d_e).argmax()), int(np.abs(d_r).argmax())
        cmp = {"n": len(m2.atoms), "atoms": atoms, "resid": rid, "d_ens": d_e, "d_ref": d_r,
               "rms_ens": math.sqrt((d_e ** 2).mean()), "max_ens": abs(d_e[ke]),
               "max_ens_atom": f"{atoms[ke]} of residue {rid[ke] + 1}",
               "rms_ref": math.sqrt((d_r ** 2).mean()), "max_ref": abs(d_r[kr]),
               "max_ref_atom": f"{atoms[kr]} of residue {rid[kr] + 1}",
               "res_net_whole": res_net, "per_res": per_res}
        with open(os.path.join(out, "comparison_whole.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["residue_number", "residue", "atom", "q_whole", "q_assembled",
                        "q_assembled_reference", "diff", "diff_reference"])
            for k, n, q, de, dr in zip(rid, atoms, qws, d_e, d_r):
                w.writerow([k + 1, chain[k], n, f"{q:.6f}", f"{q + de:.6f}",
                            f"{q + dr:.6f}", f"{de:.6f}", f"{dr:.6f}"])
        print(f"Whole chain vs assembled residues: RMS {cmp['rms_ens']:.4f} e "
              f"(reference conformers {cmp['rms_ref']:.4f} e); residue net charges in "
              f"the whole chain " + " ".join(f"{q:+.4f}" for q in res_net))

        import matplotlib
        headless = sys.platform.startswith("linux") and not (
            os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
        show = not args.no_show and not headless
        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        style(matplotlib)
        fig = figure(cmp, chain, plt, args.label_min)
        fpath = os.path.join(out, "fragment_vs_whole.png")
        fig.savefig(fpath, dpi=args.dpi, facecolor="white")
        print(f"wrote {os.path.relpath(fpath)}")
        if show:
            plt.show()

    report(os.path.join(out, "fragment_report.md"), args, residues, chain, leap, cmp, gauche)
    print(f"Residues: {', '.join(os.path.relpath(os.path.join(out, n + '.mol2')) for n in residues)}")
    if tleap and chain:
        print(f"Library:  {os.path.relpath(os.path.join(out, 'tleap', lib))}")
    print(f"Report:   {os.path.relpath(os.path.join(out, 'fragment_report.md'))}")


if __name__ == "__main__":
    main()
