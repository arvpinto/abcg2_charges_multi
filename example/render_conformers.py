#!/usr/bin/env python3
"""Header image: accepted AM1 conformers overlaid, atoms coloured by final charge."""
import argparse, csv, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--pdb", required=True); ap.add_argument("--summary", required=True)
ap.add_argument("--mol2", required=True, help="for bonds"); ap.add_argument("-o", "--out", required=True)
ap.add_argument("--dpi", type=int, default=600); ap.add_argument("--elev", type=float, default=18)
ap.add_argument("--azim", type=float, default=-62)
ap.add_argument("--transparent", action="store_true", help="transparent background instead of black")
ap.add_argument("--bond-color", default="#9aa3b2", help="bond colour")
a = ap.parse_args()

models, cur = [], []
for l in open(a.pdb):
    if l.startswith(("ATOM", "HETATM")):
        cur.append((l[12:16].strip(), float(l[30:38]), float(l[38:46]), float(l[46:54])))
    elif l.startswith("ENDMDL"):
        models.append(cur); cur = []
names = [n for n, *_ in models[0]]
X = [np.array([[x, y, z] for _, x, y, z in m]) for m in models]
q = {r["atom"]: float(r["q_final"]) for r in csv.DictReader(open(a.summary))}
bonds, sec, ids = [], None, {}
for l in open(a.mol2):
    s = l.strip()
    if s.startswith("@<TRIPOS>"): sec = s; continue
    f = s.split()
    if sec == "@<TRIPOS>ATOM" and f: ids[int(f[0])] = f[1]
    if sec == "@<TRIPOS>BOND" and len(f) >= 3: bonds.append((names.index(ids[int(f[1])]), names.index(ids[int(f[2])])))

ring = [names.index(n) for n in ("C4", "C5", "C6", "C7", "C8", "C9", "C10", "O4", "O5")]
def kabsch(P, Q):
    pc, qc = P[ring].mean(0), Q[ring].mean(0)
    U, S, Vt = np.linalg.svd((P[ring] - pc).T @ (Q[ring] - qc)); d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return (R @ (P - pc).T).T + qc
ref = X[0] - X[0][ring].mean(0)
# orient: long axis along x
u, s, vt = np.linalg.svd(ref - ref.mean(0)); ref = (ref - ref.mean(0)) @ vt.T
X = [kabsch(x, ref) for x in X]

cmap = LinearSegmentedColormap.from_list("q", ["#ff3b3b", "#ff9a8a", "#f2f2f2", "#8fb4ff", "#2f6bff"])
norm = TwoSlopeNorm(vmin=-0.75, vcenter=0.0, vmax=0.75)
heavy = [i for i, n in enumerate(names) if not n.startswith("H")]
bg = "none" if a.transparent else "black"
fig = plt.figure(figsize=(7.2, 4.0), facecolor=bg)
ax = fig.add_axes([0, 0, 1, 1], projection="3d", facecolor=bg, proj_type="ortho")
ax.set_axis_off(); ax.view_init(elev=a.elev, azim=a.azim)
for k, x in enumerate(X):
    alpha = 1.0 if k == 0 else 0.7
    for i, j in bonds:
        ax.plot(*zip(x[i], x[j]), color=a.bond_color, lw=2.2 if k == 0 else 1.4, alpha=alpha * 0.8, solid_capstyle="round")
    sizes = np.array([90 if n.startswith("H") else 260 for n in names])
    c = [cmap(norm(q[n])) for n in names]
    ax.scatter(x[:, 0], x[:, 1], x[:, 2], s=sizes, c=c, alpha=alpha, edgecolors="black", linewidths=0.4, depthshade=False)
allx = np.vstack(X); lo, hi = allx.min(0) - 0.9, allx.max(0) + 0.9
ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
ax.set_box_aspect(tuple(hi - lo))
fig.savefig(a.out, dpi=a.dpi, facecolor=bg, transparent=a.transparent)
print("wrote", a.out, len(X), "conformers")
