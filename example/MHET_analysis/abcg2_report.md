# ABCG2 conformer-averaged charges: analysis

## MHET

- SMILES `O=C([O-])c1ccc(C(=O)OCCO)cc1`, net charge -1, residue MHT, abcg2 / gaff2
- Input issues: none
- Reference conformer: input
- Pool: 9 embedded, 4 pass the filters (500 requested, seed 2026, mmff94s eps = 4.0, window 15.0 kcal/mol)
- Runs: 4, accepted 4, rejected 0
- Accepted conformers: Rg 3.26 to 3.37 A, gauche torsions 0 to 1 (median 1), 18 s per run on average
- Convergence: NO; criteria N >= 20, max SEM <= 0.001 e, max drift over 5 <= 0.0005 e
- Final: max SEM 0.0107 e, max drift nan e, largest conformational SD 0.0214 e (C2)
- tleap check: charge -1.000000, errors 0, warnings 1

Ensemble vs reference conformer: max |diff| 0.0157 e (C1), RMS 0.0045 e

| atom (equivalents) | type | final q (e) | SD | SEM | 95% CI | reference |  r(q, Rg) |
|---|---|---|---|---|---|---|---|
| O1 (1) | oh | -0.6712 | 0.0013 | 0.0006 | -0.6725 to -0.6703 | -0.6710 | +0.86 |
| C1 (1) | c3 | +0.1379 | 0.0186 | 0.0093 | +0.1221 to +0.1536 | +0.1536 | -0.49 |
| C2 (1) | c3 | +0.1336 | 0.0214 | 0.0107 | +0.1221 to +0.1551 | +0.1236 | +0.97 |
| O2 (1) | os | -0.4778 | 0.0085 | 0.0043 | -0.4841 to -0.4706 | -0.4841 | +0.76 |
| C3 (1) | c | +0.6132 | 0.0029 | 0.0014 | +0.6107 to +0.6157 | +0.6157 | -0.30 |
| O3 (1) | o | -0.5292 | 0.0029 | 0.0014 | -0.5315 to -0.5270 | -0.5270 | -0.68 |
| C4 (1) | ca | -0.1791 | 0.0010 | 0.0005 | -0.1801 to -0.1786 | -0.1786 | +0.60 |
| C5 (2) | ca | -0.0773 | 0.0003 | 0.0001 | -0.0775 to -0.0771 | -0.0775 | -0.60 |
| C6 (2) | ca | -0.1137 | 0.0005 | 0.0003 | -0.1143 to -0.1135 | -0.1135 | -0.95 |
| C7 (1) | ca | -0.0808 | 0.0005 | 0.0003 | -0.0814 to -0.0806 | -0.0806 | -0.95 |
| C10 (1) | c | +0.6936 | 0.0000 | 0.0000 | +0.6936 to +0.6936 | +0.6936 | +nan |
| O4 (1) | o | -0.7147 | 0.0005 | 0.0003 | -0.7153 to -0.7145 | -0.7145 | -0.95 |
| O5 (1) | o | -0.7147 | 0.0005 | 0.0003 | -0.7153 to -0.7145 | -0.7145 | -0.95 |
| HO1 (1) | ho | +0.4508 | 0.0100 | 0.0050 | +0.4412 to +0.4570 | +0.4570 | -0.88 |
| H11 (2) | h1 | +0.0472 | 0.0049 | 0.0024 | +0.0427 to +0.0517 | +0.0472 | -0.95 |
| H21 (2) | h1 | +0.0609 | 0.0059 | 0.0030 | +0.0557 to +0.0659 | +0.0557 | +0.27 |
| H5 (2) | ha | +0.1105 | 0.0012 | 0.0006 | +0.1095 to +0.1115 | +0.1095 | +0.60 |
| H6 (2) | ha | +0.1417 | 0.0003 | 0.0001 | +0.1415 to +0.1419 | +0.1415 | -0.60 |

Charge-shape correlations with |r| >= 0.6: O1 vs Rg (A): r = +0.86; O1 vs gauche torsions: r = -0.66; O1 vs closest distant N/O (A): r = +0.62; C1 vs gauche torsions: r = +0.73; C1 vs closest distant N/O (A): r = +0.80; C2 vs Rg (A): r = +0.97; C2 vs gauche torsions: r = -1.00; O2 vs Rg (A): r = +0.76; O2 vs gauche torsions: r = -0.92; C3 vs closest distant N/O (A): r = +0.90; O3 vs Rg (A): r = -0.68; O3 vs gauche torsions: r = +0.87; O3 vs closest distant N/O (A): r = +0.64; C4 vs Rg (A): r = +0.60; C4 vs closest distant N/O (A): r = +0.87; C5 vs Rg (A): r = -0.60; C5 vs closest distant N/O (A): r = -0.87; C6 vs Rg (A): r = -0.95; C6 vs gauche torsions: r = +1.00; C7 vs Rg (A): r = -0.95; C7 vs gauche torsions: r = +1.00; O4 vs Rg (A): r = -0.95; O4 vs gauche torsions: r = +1.00; O5 vs Rg (A): r = -0.95; O5 vs gauche torsions: r = +1.00; HO1 vs Rg (A): r = -0.88; HO1 vs gauche torsions: r = +0.98; H11 vs Rg (A): r = -0.95; H11 vs gauche torsions: r = +0.82; H21 vs closest distant N/O (A): r = -0.92; H5 vs Rg (A): r = +0.60; H5 vs gauche torsions: r = -0.82; H5 vs closest distant N/O (A): r = -0.71; H6 vs Rg (A): r = -0.60; H6 vs closest distant N/O (A): r = -0.87.

## Figures

- `abcg2_convergence.png`
- `abcg2_running_mean.png`
- `abcg2_charges.png`
- `abcg2_descriptors.png`
