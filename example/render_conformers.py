#!/usr/bin/env python3
"""Header image: accepted AM1 conformers overlaid, atoms coloured by final charge.

    python3 render_conformers.py --pdb MHET_abcg2/MHET_accepted_sqm.pdb \
        --summary MHET_abcg2/MHET_charges_summary.csv \
        --mol2 MHET_abcg2/MHET_abcg2_ensemble.mol2 -o mhet_conformers.png

For a covalent run (abcg2_covalent_charges.py) the conformers are aligned on
the frozen backbone and the caps are drawn in grey by default (--fixed grey).
"""
import argparse, csv, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--pdb", required=True, help="multi-model PDB of the accepted conformers")
ap.add_argument("--summary", required=True, help="*_charges_summary.csv")
ap.add_argument("--mol2", required=True, help="for bonds")
ap.add_argument("-o", "--out", required=True, help="image file (png, webp needs Pillow)")
ap.add_argument("--dpi", type=int, default=600)
ap.add_argument("--elev", type=float, default=18)
ap.add_argument("--azim", type=float, default=-62)
ap.add_argument("--align", nargs="+", default=None,
                help="atoms used to superpose the conformers (default: the ring and "
                     "carboxylate of MHET, or the frozen backbone of a covalent run)")
ap.add_argument("--fixed", choices=["grey", "charge", "hide"], default="grey",
                help="how atoms with fixed charges (caps) are drawn")
ap.add_argument("--qmax", type=float, default=0.75, help="charge at the ends of the colour scale")
ap.add_argument("--size", nargs=2, type=float, default=[7.2, 4.0], help="figure size (in)")
ap.add_argument("--transparent", action="store_true", help="transparent background instead of black")
ap.add_argument("--bond-color", default="#9aa3b2", help="bond colour")
a = ap.parse_args()

models, cur = [], []
for l in open(a.pdb):
    if l.startswith(("ATOM", "HETATM")):
        cur.append((l[12:16].strip(), float(l[30:38]), float(l[38:46]), float(l[46:54])))
    elif l.startswith("ENDMDL"):
        models.append(cur); cur = []
if cur:
    models.append(cur)
names = [n for n, *_ in models[0]]
X = [np.array([[x, y, z] for _, x, y, z in m]) for m in models]
rows = list(csv.DictReader(open(a.summary)))
q = {r["atom"]: float(r["q_final"]) for r in rows}
fixed = {r["atom"] for r in rows if r.get("fixed") == "1"}
role = {r["atom"]: r.get("role", "") for r in rows}
bonds, sec, ids = [], None, {}
for l in open(a.mol2):
    s = l.strip()
    if s.startswith("@<TRIPOS>"): sec = s; continue
    f = s.split()
    if sec == "@<TRIPOS>ATOM" and f: ids[int(f[0])] = f[1]
    if sec == "@<TRIPOS>BOND" and len(f) >= 3: bonds.append((names.index(ids[int(f[1])]), names.index(ids[int(f[2])])))

if a.align:
    fit = [names.index(n) for n in a.align]
elif fixed:
    fit = [i for i, n in enumerate(names) if role[n] in ("ace", "nme", "backbone")]
else:
    fit = [names.index(n) for n in ("C4", "C5", "C6", "C7", "C8", "C9", "C10", "O4", "O5")]

def kabsch(P, Q):
    pc, qc = P[fit].mean(0), Q[fit].mean(0)
    U, S, Vt = np.linalg.svd((P[fit] - pc).T @ (Q[fit] - qc)); d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return (R @ (P - pc).T).T + qc
ref = X[0] - X[0][fit].mean(0)
# orient: long axis of the whole molecule along x
u, s, vt = np.linalg.svd(ref - ref.mean(0)); ref = (ref - ref.mean(0)) @ vt.T
X = [kabsch(x, ref) for x in X]

hide = {i for i, n in enumerate(names) if n in fixed} if a.fixed == "hide" else set()
cmap = LinearSegmentedColormap.from_list("q", ["#ff3b3b", "#ff9a8a", "#f2f2f2", "#8fb4ff", "#2f6bff"])
norm = TwoSlopeNorm(vmin=-a.qmax, vcenter=0.0, vmax=a.qmax)
bg = "none" if a.transparent else "black"
fig = plt.figure(figsize=tuple(a.size), facecolor=bg)
ax = fig.add_axes([0, 0, 1, 1], projection="3d", facecolor=bg, proj_type="ortho")
ax.set_axis_off(); ax.view_init(elev=a.elev, azim=a.azim)
for k, x in enumerate(X):
    alpha = 1.0 if k == 0 else 0.7
    for i, j in bonds:
        if i in hide or j in hide:
            continue
        ax.plot(*zip(x[i], x[j]), color=a.bond_color, lw=2.2 if k == 0 else 1.4, alpha=alpha * 0.8, solid_capstyle="round")
    keep = [i for i in range(len(names)) if i not in hide]
    sizes = np.array([90 if names[i].startswith("H") else 260 for i in keep])
    c = [(0.42, 0.44, 0.48, 1.0) if (names[i] in fixed and a.fixed == "grey") else cmap(norm(q[names[i]])) for i in keep]
    ax.scatter(x[keep, 0], x[keep, 1], x[keep, 2], s=sizes, c=c, alpha=alpha, edgecolors="black", linewidths=0.4, depthshade=False)
allx = np.vstack([x[[i for i in range(len(names)) if i not in hide]] for x in X]); lo, hi = allx.min(0) - 0.9, allx.max(0) + 0.9
ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
ax.set_box_aspect(tuple(hi - lo))
fig.savefig(a.out, dpi=a.dpi, facecolor=bg, transparent=a.transparent)
print("wrote", a.out, len(X), "conformers")
