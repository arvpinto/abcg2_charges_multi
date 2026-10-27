# ABCG2 conformer-averaged charges: analysis

## SXM

- SMILES `CNC(=O)[C@H](COC(=O)c1ccc(C(=O)[O-])cc1)NC(C)=O`, net charge -1, residue SXM, abcg2 / gaff2
- Input issues: H-bond NT-HNT...OG 2.15 A
  - input amide OY=CY-N-H  -168.9 deg (trans)
  - input amide O=C-NT-HNT  -170.5 deg (trans)
- Reference conformer: input, bound pose
- Pool: 15 embedded, 5 pass the filters (500 requested, seed 2026, mmff94s eps = 4.0, window 15.0 kcal/mol)
- Runs: 5, accepted 3, rejected 2
  - 2 x 1 H-bond
- Accepted conformers: Rg 3.95 to 4.08 A, gauche torsions 0 to 0 (median 0), 104 s per run on average
- Convergence: NO; criteria N >= 20, max SEM <= 0.001 e, max drift over 5 <= 0.0005 e
- Final: max SEM 0.0077 e, max drift nan e, largest conformational SD 0.0133 e (HA)
- tleap check: residue charge -1.000000, errors 0, warnings 3, missing parameters 0
- Covalent residue: caps ACE (CAY, HY1, HY2, HY3, CY, OY) and NME (NT, HNT, CAT, HT1, HT2, HT3); backbone fixed in the pool; AM1 geometry restrained
- ABCG2 charge of the caps: ACE +0.0053 e, NME +0.0413 e; fixed to +0.0000 and +0.0000 e; 12 fixed atoms, redistribution even (+0.00185 e per free atom)
- Residue SXM without caps: 25 atoms, charge -1.000000 e (SXM_abcg2_residue.mol2)
- tleap ACE-SXM-NME before the junction parameters: 12 missing terms; junction frcmod: SXM_junction.frcmod

Ensemble vs reference conformer: max |diff| 0.0279 e (OG), RMS 0.0077 e

For fixed atoms the final charge was set, not computed; their SD, SEM and 95% CI describe the ABCG2 charges of the conformers.

| atom (equivalents) | type | final q (e) | SD | SEM | 95% CI | reference |  r(q, Rg) |
|---|---|---|---|---|---|---|---|
| CAY (1, fixed) | c3 | -0.3662 | 0.0000 | 0.0000 | -0.1801 to -0.1801 | -0.3662 | +nan |
| HY1 (3, fixed) | hc | +0.1123 | 0.0010 | 0.0006 | +0.0570 to +0.0587 | +0.1123 | -1.00 |
| CY (1, fixed) | c | +0.5972 | 0.0023 | 0.0013 | +0.6011 to +0.6051 | +0.5972 | -1.00 |
| OY (1, fixed) | o | -0.5679 | 0.0248 | 0.0143 | -0.6071 to -0.5641 | -0.5679 | +1.00 |
| N (1) | ns | -0.4987 | 0.0081 | 0.0047 | -0.5099 to -0.4959 | -0.5030 | -1.00 |
| H (1) | hn | +0.3384 | 0.0069 | 0.0040 | +0.3285 to +0.3405 | +0.3384 | -1.00 |
| CA (1) | c3 | +0.0272 | 0.0012 | 0.0007 | +0.0247 to +0.0267 | +0.0216 | +1.00 |
| HA (1) | h1 | +0.0842 | 0.0133 | 0.0077 | +0.0747 to +0.0977 | +0.0956 | +1.00 |
| CB (1) | c3 | +0.1578 | 0.0006 | 0.0003 | +0.1556 to +0.1566 | +0.1615 | -0.49 |
| HB2 (2) | h1 | +0.0604 | 0.0051 | 0.0029 | +0.0527 to +0.0617 | +0.0696 | -1.00 |
| OG (1) | os | -0.4642 | 0.0121 | 0.0070 | -0.4731 to -0.4521 | -0.4922 | +1.00 |
| C3 (1) | c | +0.6159 | 0.0006 | 0.0003 | +0.6137 to +0.6147 | +0.6206 | +1.00 |
| O3 (1) | o | -0.5308 | 0.0064 | 0.0037 | -0.5400 to -0.5290 | -0.5121 | -1.00 |
| C4 (1) | ca | -0.1794 | 0.0006 | 0.0003 | -0.1816 to -0.1806 | -0.1787 | -0.49 |
| C5 (2) | ca | -0.0748 | 0.0006 | 0.0003 | -0.0770 to -0.0760 | -0.0766 | +1.00 |
| H5 (2) | ha | +0.1182 | 0.0040 | 0.0023 | +0.1140 to +0.1210 | +0.1119 | +1.00 |
| C6 (2) | ca | -0.1131 | 0.0009 | 0.0005 | -0.1160 to -0.1145 | -0.1091 | -1.00 |
| H6 (2) | ha | +0.1437 | 0.0008 | 0.0004 | +0.1410 to +0.1425 | +0.1464 | -0.94 |
| C7 (1) | ca | -0.0791 | 0.0006 | 0.0003 | -0.0816 to -0.0806 | -0.0777 | -1.00 |
| C10 (1) | c | +0.6955 | 0.0000 | 0.0000 | +0.6936 to +0.6936 | +0.6975 | +nan |
| O4 (1) | o | -0.7138 | 0.0012 | 0.0007 | -0.7170 to -0.7150 | -0.7096 | -1.00 |
| O5 (1) | o | -0.7138 | 0.0012 | 0.0007 | -0.7170 to -0.7150 | -0.7096 | -1.00 |
| C (1) | c | +0.6000 | 0.0104 | 0.0060 | +0.5861 to +0.6041 | +0.5900 | -1.00 |
| O (1) | o | -0.6076 | 0.0040 | 0.0023 | -0.6141 to -0.6071 | -0.6272 | -1.00 |
| NT (1, fixed) | ns | -0.4157 | 0.0029 | 0.0017 | -0.5019 to -0.4969 | -0.4157 | -1.00 |
| HNT (1, fixed) | hn | +0.2719 | 0.0052 | 0.0030 | +0.3165 to +0.3255 | +0.2719 | +1.00 |
| CAT (1, fixed) | c3 | -0.1490 | 0.0000 | 0.0000 | +0.0953 to +0.0953 | -0.1490 | +nan |
| HT1 (3, fixed) | h1 | +0.0976 | 0.0017 | 0.0010 | +0.0397 to +0.0427 | +0.0976 | -1.00 |

Charge-shape correlations with |r| >= 0.6: N vs Rg (A): r = -1.00; N vs closest distant N/O (A): r = +1.00; H vs Rg (A): r = -1.00; H vs closest distant N/O (A): r = +1.00; CA vs Rg (A): r = +1.00; CA vs closest distant N/O (A): r = -1.00; HA vs Rg (A): r = +1.00; HA vs closest distant N/O (A): r = -1.00; HB2 vs Rg (A): r = -1.00; HB2 vs closest distant N/O (A): r = +1.00; OG vs Rg (A): r = +1.00; OG vs closest distant N/O (A): r = -1.00; C3 vs Rg (A): r = +1.00; C3 vs closest distant N/O (A): r = -1.00; O3 vs Rg (A): r = -1.00; O3 vs closest distant N/O (A): r = +1.00; C5 vs Rg (A): r = +1.00; C5 vs closest distant N/O (A): r = -1.00; H5 vs Rg (A): r = +1.00; H5 vs closest distant N/O (A): r = -1.00; C6 vs Rg (A): r = -1.00; C6 vs closest distant N/O (A): r = +1.00; H6 vs Rg (A): r = -0.94; H6 vs closest distant N/O (A): r = +0.94; C7 vs Rg (A): r = -1.00; C7 vs closest distant N/O (A): r = +1.00; O4 vs Rg (A): r = -1.00; O4 vs closest distant N/O (A): r = +1.00; O5 vs Rg (A): r = -1.00; O5 vs closest distant N/O (A): r = +1.00; C vs Rg (A): r = -1.00; C vs closest distant N/O (A): r = +1.00; O vs Rg (A): r = -1.00; O vs closest distant N/O (A): r = +1.00.

## Figures

- `abcg2_convergence.png`
- `abcg2_running_mean.png`
- `abcg2_charges.png`
- `abcg2_descriptors.png`
