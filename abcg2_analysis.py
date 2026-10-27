#!/usr/bin/env python3
"""
abcg2_analysis.py

Analysis of one or more output directories of abcg2_conformer_charges.py.

    python3 abcg2_analysis.py EZ60_abcg2 EZ61_abcg2

Figures (saved at --dpi, 600 by default, and shown on screen unless
--no-show or no display is available):
  convergence   worst-atom SEM and running-mean drift vs N, with tolerances
  running_mean  running mean minus final mean for the most variable atoms
  charges       per symmetry-unique atom: final charge, conformational SD,
                SEM, the single reference conformer and the input-file
                charges (or other mol2 files given with --compare),
                all relative to the ensemble mean
  descriptors   charge of the most variable atoms vs conformer shape
                (radius of gyration, closest distant N/O contact) with
                Pearson r

Report (<out>/abcg2_report.md): settings, pool and run statistics,
convergence, per-atom table with 95% CI, comparison with the reference
conformer and the input charges, and charge-shape correlations.
"""

import argparse
import csv
import json
import math
import os
import re
import sys

import numpy as np
import matplotlib

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]
TEXT_1, TEXT_2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------
def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def load_run(d):
    js = [f for f in os.listdir(d) if f.endswith("_summary.json")]
    if not js:
        sys.exit(f"ERROR: {d} has no *_summary.json (not an "
                 f"abcg2_conformer_charges.py output directory?)")
    prefix = js[0][:-len("_summary.json")]
    with open(os.path.join(d, js[0])) as fh:
        summary = json.load(fh)
    run = {"dir": d, "prefix": prefix, "summary": summary,
           "atoms": read_csv(os.path.join(d, f"{prefix}_charges_summary.csv")),
           "running": read_csv(os.path.join(d, f"{prefix}_running_mean.csv")),
           "runs": read_csv(os.path.join(d, f"{prefix}_conformer_runs.csv")),
           "pool": read_csv(os.path.join(d, f"{prefix}_pool.csv"))}
    per = read_csv(os.path.join(d, f"{prefix}_charges_per_conformer.csv"))
    run["tags"] = [k for k in per[0].keys() if k != "atom"]
    run["Q"] = {r["atom"]: np.array([float(r[t]) for t in run["tags"]]) for r in per}
    return run


def read_mol2_charges(path):
    q, sec = {}, None
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("@<TRIPOS>"):
                sec = s
            elif sec == "@<TRIPOS>ATOM" and s:
                f = s.split()
                if len(f) > 8:
                    q[f[1]] = float(f[8])
    return q


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def element(name):
    letters = "".join(ch for ch in name if ch.isalpha())
    return letters[:1].upper()


def unique_atoms(run, hydrogens):
    seen, out = set(), []
    for r in run["atoms"]:
        if r["sym_class"] in seen:
            continue
        seen.add(r["sym_class"])
        if element(r["atom"]) == "H":
            if hydrogens == "none":
                continue
            if hydrogens == "polar" and not r["type"].lower().startswith(("hn", "ho", "hs")):
                continue
        out.append(r)
    return out


def most_variable(run, k):
    rows = unique_atoms(run, "all")
    rows.sort(key=lambda r: -fnum(r["std"]))
    return [r["atom"] for r in rows[:k]]


def comparison_charges(run, spec):
    """Return ({atom: q} or None, label) for the comparison series."""
    if spec is None or spec.lower() == "none":
        return None, None
    if spec.lower() == "input":
        q = {r["atom"]: fnum(r.get("q_input_file")) for r in run["atoms"]}
        if all(math.isnan(v) for v in q.values()):
            return None, None
        return q, "input mol2 charges"
    return read_mol2_charges(spec), os.path.basename(spec)


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return float("nan")
    return float(np.corrcoef(x[m], y[m])[0, 1])


def descriptors(run):
    by_tag = {r["tag"]: r for r in run["runs"]}
    d = {"Rg (A)": [], "gauche torsions": [], "closest distant N/O (A)": []}
    for t in run["tags"]:
        r = by_tag.get(t, {})
        d["Rg (A)"].append(fnum(r.get("post_rg_A")))
        d["gauche torsions"].append(fnum(r.get("post_gauche")))
        d["closest distant N/O (A)"].append(fnum(r.get("post_min_polar_dist_A")))
    return {k: np.array(v) for k, v in d.items()}


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------
def style(mpl):
    mpl.rcParams.update({
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "axes.edgecolor": TEXT_2, "axes.labelcolor": TEXT_1,
        "xtick.color": TEXT_2, "ytick.color": TEXT_2,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "legend.frameon": False, "legend.fontsize": 7,
        "lines.linewidth": 1.6, "figure.dpi": 110})


def fig_convergence(runs, labels, plt):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7), constrained_layout=True)
    s0 = runs[0]["summary"]["settings"]
    specs = [("max_sem", "Max per-atom SEM (e)", s0["sem_tol"],
              "Standard error of the mean, worst atom"),
             ("max_drift", "Max drift of the running mean (e)", s0["drift_tol"],
              f"|mean(N) - mean(N-{s0['drift_window']})|, worst atom")]
    for ax, (key, ylab, tol, title) in zip(axes, specs):
        for k, (run, lab) in enumerate(zip(runs, labels)):
            N = np.array([int(r["N"]) for r in run["running"]])
            y = np.array([fnum(r[key]) for r in run["running"]])
            m = np.isfinite(y) & (y > 0)
            if not m.any():
                continue
            ax.plot(N[m], y[m], color=SERIES[k % 4], label=lab)
            ax.annotate(lab, (N[m][-1], y[m][-1]), xytext=(4, 0),
                        textcoords="offset points", va="center", fontsize=7,
                        color=TEXT_1)
        ax.axhline(tol, color=TEXT_2, lw=0.9, ls="--")
        ax.annotate(f"tolerance {tol:g} e", (0.02, tol),
                    xycoords=("axes fraction", "data"), xytext=(0, -3),
                    textcoords="offset points", ha="left", va="top",
                    fontsize=7, color=TEXT_2)
        ax.axvline(s0["min_confs"], color=GRID, lw=1.0)
        ax.set_yscale("log")
        ax.set_xlabel("Accepted conformers, N")
        ax.set_ylabel(ylab)
        ax.set_title(title, loc="left")
    if len(runs) > 1:
        axes[0].legend(loc="upper right")
    return fig


def fig_running(runs, labels, atoms_opt, plt):
    n = len(runs)
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 2.8), sharey=True,
                             constrained_layout=True, squeeze=False)
    for ax, run, lab in zip(axes[0], runs, labels):
        atoms = (atoms_opt or most_variable(run, 4))[:4]
        N = np.array([int(r["N"]) for r in run["running"]])
        for k, a in enumerate(atoms):
            if a not in run["running"][0]:
                continue
            y = np.array([fnum(r[a]) for r in run["running"]])
            ax.plot(N, 1e3 * (y - y[-1]), color=SERIES[k], label=a, lw=1.4)
        ax.axhline(0, color=TEXT_2, lw=0.8)
        ax.set_xlabel("Accepted conformers, N")
        ax.set_title(f"{lab}: running mean - final mean", loc="left")
        ax.legend(loc="upper right", ncol=2)
    axes[0][0].set_ylabel("Charge deviation (me)")
    return fig


def fig_charges(runs, labels, comps, hydrogens, plt):
    n = len(runs)
    fig, axes = plt.subplots(n, 1, figsize=(7.0, 2.8 * n),
                             constrained_layout=True, squeeze=False)
    for ax, run, lab, (cq, clab) in zip(axes[:, 0], runs, labels, comps):
        rows = unique_atoms(run, hydrogens)
        x = np.arange(len(rows))
        std = 1e3 * np.array([fnum(r["std"]) for r in rows])
        sem = 1e3 * np.array([fnum(r["sem"]) for r in rows])
        final = np.array([fnum(r["q_final"]) for r in rows])
        ax.bar(x, 2 * sem, bottom=-sem, width=0.5, color=SERIES[0], alpha=0.25,
               lw=0, label="ensemble mean ± SEM")
        ax.errorbar(x, np.zeros(len(x)), yerr=std, fmt="o", ms=4.0, lw=1.2,
                    capsize=2.5, color=SERIES[0],
                    label=f"ensemble ± conformational SD "
                          f"(N = {run['summary']['accepted']})")
        ref = np.array([fnum(r.get("q_reference")) for r in rows])
        if np.isfinite(ref).any():
            ax.plot(x - 0.22, 1e3 * (ref - final), "s", ms=4.0, color=SERIES[1],
                    label=f"single reference conformer "
                          f"({run['summary']['reference_geometry']})")
        if cq:
            c = np.array([cq.get(r["atom"], np.nan) for r in rows])
            ax.plot(x + 0.22, 1e3 * (c - final), "D", ms=3.8, color=SERIES[2],
                    label=clab)
        ax.axhline(0, color=TEXT_2, lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{r['atom']}\n{fnum(r['q_final']):+.3f}" for r in rows],
                           fontsize=6.5)
        ax.set_ylabel("q - q(ensemble) (me)")
        ax.set_title(f"{lab}: charges relative to the ensemble mean "
                     f"(tick labels: final charge, e)", loc="left")
        ax.legend(loc="upper left", ncol=2, bbox_to_anchor=(0, -0.22))
        ax.grid(axis="x", visible=False)
    return fig


def fig_descriptors(runs, labels, atoms_opt, plt):
    keys = ["Rg (A)", "closest distant N/O (A)"]
    n = len(runs)
    fig, axes = plt.subplots(n, 2, figsize=(7.0, 2.6 * n),
                             constrained_layout=True, squeeze=False)
    for i, (run, lab) in enumerate(zip(runs, labels)):
        atoms = (atoms_opt or most_variable(run, 3))[:3]
        d = descriptors(run)
        for j, key in enumerate(keys):
            ax = axes[i, j]
            x = d[key]
            for k, a in enumerate(atoms):
                if a not in run["Q"]:
                    continue
                q = run["Q"][a]
                y = 1e3 * (q - q.mean())
                r = pearson(x, q)
                ax.plot(x, y, "o", ms=4.0, color=SERIES[k], alpha=0.85,
                        markeredgecolor="white", markeredgewidth=0.6,
                        label=f"{a} (r = {r:+.2f})")
            ax.axhline(0, color=TEXT_2, lw=0.8)
            ax.set_xlabel(f"{key}, after sqm")
            if j == 0:
                ax.set_ylabel("q - mean (me)")
            ax.set_title(f"{lab}: charge vs {key.split(' (')[0]}", loc="left")
            ax.legend(loc="best")
    return fig


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------
def report(runs, labels, comps, fig_files, path):
    L = ["# ABCG2 conformer-averaged charges: analysis", ""]
    for run, lab, (cq, clab) in zip(runs, labels, comps):
        s, st = run["summary"], run["summary"]["settings"]
        L += [f"## {lab}", "",
              f"- SMILES `{s['smiles']}`, net charge {s['net_charge']}"
              + (f", residue {s['resname']}" if s.get("resname") else "")
              + f", {s['method']} / {s['atomtype']}",
              f"- Input issues: {s.get('input_issues') or 'none'}"]
        for a in s.get("input_amides", []):
            L.append(f"  - input amide {a}")
        L += [f"- Reference conformer: {s['reference_geometry']}",
              f"- Pool: {s['pool_embedded']} embedded, {s['pool_accepted']} pass "
              f"the filters ({st['nconf']} requested, seed {st['seed']}, "
              f"{st['ff']} eps = {st['dielectric']}, window {st['ewin']} kcal/mol)",
              f"- Runs: {s['runs']}, accepted {s['accepted']}, rejected "
              f"{s.get('rejected', s.get('rejected_after_sqm', 0))}"]
        reasons = {}
        for r in run["runs"]:
            if r["accepted"] != "1" and r["reason"]:
                key = re.sub(r"\s-?\d+\.\d+ A", "", r["reason"])
                reasons[key] = reasons.get(key, 0) + 1
        for k, v in sorted(reasons.items(), key=lambda kv: -kv[1]):
            L.append(f"  - {v} x {k}")
        acc = [r for r in run["runs"] if r["accepted"] == "1"]
        rg = [fnum(r["post_rg_A"]) for r in acc]
        ga = [fnum(r["post_gauche"]) for r in acc]
        tm = [fnum(r["time_s"]) for r in acc]
        if acc:
            L.append(f"- Accepted conformers: Rg {min(rg):.2f} to {max(rg):.2f} A, "
                     f"gauche torsions {min(ga):.0f} to {max(ga):.0f} (median "
                     f"{np.median(ga):.0f}), {np.mean(tm):.0f} s per run on average")
        L += [f"- Convergence: {'yes' if s['converged'] else 'NO'}; criteria N >= "
              f"{st['min_confs']}, max SEM <= {st['sem_tol']} e, max drift over "
              f"{st['drift_window']} <= {st['drift_tol']} e",
              f"- Final: max SEM {s['final_max_sem']:.4f} e, max drift "
              f"{s['final_max_drift']:.5f} e, largest conformational SD "
              f"{s['max_std']:.4f} e ({s['max_std_atom']})"]
        if s.get("tleap"):
            L.append(f"- tleap check: {s['tleap']}")
        L.append("")

        rows = unique_atoms(run, "all")
        final = {r["atom"]: fnum(r["q_final"]) for r in run["atoms"]}
        ref = {r["atom"]: fnum(r.get("q_reference")) for r in run["atoms"]}
        allq = np.array([final[a] for a in final])
        refq = np.array([ref[a] for a in final])
        if np.isfinite(refq).all():
            d = allq - refq
            k = int(np.abs(d).argmax())
            L.append(f"Ensemble vs reference conformer: max |diff| "
                     f"{abs(d[k]):.4f} e ({list(final)[k]}), RMS "
                     f"{math.sqrt((d ** 2).mean()):.4f} e")
        if cq:
            c = np.array([cq.get(a, np.nan) for a in final])
            if np.isfinite(c).all():
                d = allq - c
                k = int(np.abs(d).argmax())
                L.append(f"Ensemble vs {clab}: max |diff| {abs(d[k]):.4f} e "
                         f"({list(final)[k]}), RMS {math.sqrt((d ** 2).mean()):.4f} e")
        L += ["", "| atom (equivalents) | type | final q (e) | SD | SEM | 95% CI | "
              "reference | " + (f"{clab} |" if cq else "") + " r(q, Rg) |",
              "|---|---|---|---|---|---|---|" + ("---|" if cq else "") + "---|"]
        desc = descriptors(run)
        counts = {}
        for r in run["atoms"]:
            counts[r["sym_class"]] = counts.get(r["sym_class"], 0) + 1
        for r in rows:
            a = r["atom"]
            line = (f"| {a} ({counts[r['sym_class']]}) | {r['type']} | "
                    f"{fnum(r['q_final']):+.4f} | {fnum(r['std']):.4f} | "
                    f"{fnum(r['sem']):.4f} | {fnum(r['ci95_lo']):+.4f} to "
                    f"{fnum(r['ci95_hi']):+.4f} | {fnum(r.get('q_reference')):+.4f} | ")
            if cq:
                line += f"{cq.get(a, float('nan')):+.4f} | "
            line += f"{pearson(desc['Rg (A)'], run['Q'][a]):+.2f} |"
            L.append(line)
        strong = []
        for r in rows:
            for key, x in desc.items():
                rr = pearson(x, run["Q"][r["atom"]])
                if np.isfinite(rr) and abs(rr) >= 0.6:
                    strong.append(f"{r['atom']} vs {key}: r = {rr:+.2f}")
        L += ["", "Charge-shape correlations with |r| >= 0.6: " +
              ("; ".join(strong) if strong else "none") + ".", ""]
    L += ["## Figures", ""] + [f"- `{os.path.basename(f)}`" for f in fig_files] + [""]
    with open(path, "w") as fh:
        fh.write("\n".join(L))


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="Figures and report for abcg2_conformer_charges.py output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("runs", nargs="+", help="output directories to analyse")
    p.add_argument("-o", "--outdir", default="abcg2_analysis",
                   help="directory for figures and report")
    p.add_argument("--labels", nargs="+", default=None,
                   help="labels for the runs (default: file prefixes)")
    p.add_argument("--compare", nargs="+", default=None,
                   help="per run: 'input' (charges in the input mol2), 'none', "
                        "or a mol2 file; one value applies to all runs")
    p.add_argument("--atoms", nargs="+", default=None,
                   help="atoms for the running-mean and descriptor plots "
                        "(default: the most variable symmetry-unique atoms)")
    p.add_argument("--hydrogens", choices=["none", "polar", "all"], default="polar",
                   help="hydrogens shown in the charge plot")
    p.add_argument("--dpi", type=int, default=600, help="figure resolution")
    p.add_argument("--format", default="png", help="figure format (png, pdf, svg)")
    p.add_argument("--no-show", action="store_true",
                   help="do not open the figures on screen")
    args = p.parse_args()

    headless = sys.platform.startswith("linux") and not (
        os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
    show = not args.no_show and not headless
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    style(matplotlib)

    runs = [load_run(d) for d in args.runs]
    labels = args.labels or [r["prefix"] for r in runs]
    if len(labels) != len(runs):
        sys.exit("ERROR: give one label per run")
    specs = args.compare or ["input"]
    if len(specs) == 1:
        specs = specs * len(runs)
    if len(specs) != len(runs):
        sys.exit("ERROR: --compare needs one value, or one per run")
    comps = [comparison_charges(r, s) for r, s in zip(runs, specs)]

    os.makedirs(args.outdir, exist_ok=True)
    figs = {"convergence": fig_convergence(runs, labels, plt),
            "running_mean": fig_running(runs, labels, args.atoms, plt),
            "charges": fig_charges(runs, labels, comps, args.hydrogens, plt),
            "descriptors": fig_descriptors(runs, labels, args.atoms, plt)}
    files = []
    for name, fig in figs.items():
        f = os.path.join(args.outdir, f"abcg2_{name}.{args.format}")
        fig.savefig(f, dpi=args.dpi, facecolor="white")
        files.append(f)
        print(f"wrote {f}")
    rep = os.path.join(args.outdir, "abcg2_report.md")
    report(runs, labels, comps, files, rep)
    print(f"wrote {rep}")
    if show:
        plt.show()
    elif not args.no_show:
        print("(no display found: figures saved but not shown)")


if __name__ == "__main__":
    main()
