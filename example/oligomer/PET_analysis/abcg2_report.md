# ABCG2 conformer-averaged charges: analysis

## MID

- SMILES `CC(=O)OCCOC(=O)c1ccc(C(=O)OCCOC(C)=O)cc1`, net charge 0, residue MID, abcg2 / gaff2
- Input issues: none
- Reference conformer: input
- Pool: 295 embedded, 172 pass the filters (500 requested, seed 2026, mmff94s eps = 4.0, window 15.0 kcal/mol)
- Runs: 67, accepted 66, rejected 0
  - 1 x finished before use
- Accepted conformers: Rg 4.18 to 5.45 A, gauche torsions 0 to 2 (median 2), 93 s per run on average
- Convergence: yes; criteria N >= 20, max SEM <= 0.001 e, max drift over 5 <= 0.0005 e
- Final: max SEM 0.0010 e, max drift 0.00022 e, largest conformational SD 0.0081 e (CP1)
- tleap check: charge 0.000000, errors 0, warnings 0

Ensemble vs reference conformer: max |diff| 0.0624 e (OP2), RMS 0.0198 e

| atom (equivalents) | type | final q (e) | SD | SEM | 95% CI | reference |  r(q, Rg) |
|---|---|---|---|---|---|---|---|
| CP1 (2) | c3 | +0.1464 | 0.0081 | 0.0010 | +0.1444 to +0.1484 | +0.1516 | -0.19 |
| HP1 (4) | h1 | +0.0762 | 0.0035 | 0.0004 | +0.0754 to +0.0770 | +0.0477 | -0.30 |
| OP1 (2) | os | -0.4648 | 0.0030 | 0.0004 | -0.4655 to -0.4641 | -0.4451 | +0.34 |
| CP2 (2) | c | +0.5854 | 0.0024 | 0.0003 | +0.5849 to +0.5860 | +0.5801 | -0.74 |
| OP2 (2) | o | -0.4924 | 0.0033 | 0.0004 | -0.4932 to -0.4916 | -0.4300 | -0.01 |
| CP3 (2) | c3 | -0.1502 | 0.0003 | 0.0000 | -0.1503 to -0.1501 | -0.1901 | +0.46 |
| HP3 (6) | hc | +0.0776 | 0.0004 | 0.0001 | +0.0775 to +0.0777 | +0.0757 | +0.43 |
| C1 (2) | c3 | +0.1445 | 0.0078 | 0.0010 | +0.1426 to +0.1463 | +0.1616 | +0.02 |
| H11 (4) | h1 | +0.0745 | 0.0046 | 0.0006 | +0.0733 to +0.0756 | +0.0747 | +0.46 |
| O1 (2) | os | -0.4584 | 0.0028 | 0.0003 | -0.4591 to -0.4578 | -0.4711 | -0.26 |
| C2 (2) | c | +0.5959 | 0.0021 | 0.0003 | +0.5954 to +0.5964 | +0.5977 | -0.48 |
| O2 (2) | o | -0.4874 | 0.0034 | 0.0004 | -0.4882 to -0.4866 | -0.4740 | +0.75 |
| C3 (2) | ca | -0.1033 | 0.0014 | 0.0002 | -0.1036 to -0.1029 | -0.1016 | +0.40 |
| C4 (4) | ca | -0.0719 | 0.0006 | 0.0001 | -0.0720 to -0.0718 | -0.0710 | +0.46 |
| H4 (4) | ha | +0.1470 | 0.0020 | 0.0002 | +0.1465 to +0.1475 | +0.1450 | -0.63 |

Charge-shape correlations with |r| >= 0.6: CP2 vs Rg (A): r = -0.74; CP2 vs gauche torsions: r = +0.65; CP2 vs closest distant N/O (A): r = -0.67; CP3 vs closest distant N/O (A): r = +0.66; H11 vs closest distant N/O (A): r = +0.61; O2 vs Rg (A): r = +0.75; O2 vs gauche torsions: r = -0.79; O2 vs closest distant N/O (A): r = +0.66; C3 vs closest distant N/O (A): r = +0.60; H4 vs Rg (A): r = -0.63; H4 vs closest distant N/O (A): r = -0.83.

## END

- SMILES `CCOC(=O)c1ccc(C(=O)OCCOC(C)=O)cc1`, net charge 0, residue END, abcg2 / gaff2
- Input issues: none
- Reference conformer: input
- Pool: 126 embedded, 28 pass the filters (500 requested, seed 2026, mmff94s eps = 4.0, window 15.0 kcal/mol)
- Runs: 28, accepted 28, rejected 0
- Accepted conformers: Rg 3.89 to 4.64 A, gauche torsions 0 to 1 (median 1), 47 s per run on average
- Convergence: NO; criteria N >= 20, max SEM <= 0.001 e, max drift over 5 <= 0.0005 e
- Final: max SEM 0.0023 e, max drift 0.00137 e, largest conformational SD 0.0122 e (CN1)
- tleap check: charge 0.000000, errors 0, warnings 0

Ensemble vs reference conformer: max |diff| 0.0602 e (ON2), RMS 0.0149 e

| atom (equivalents) | type | final q (e) | SD | SEM | 95% CI | reference |  r(q, Rg) |
|---|---|---|---|---|---|---|---|
| C1 (1) | c3 | +0.1636 | 0.0000 | 0.0000 | +0.1636 to +0.1636 | +0.1636 | +0.00 |
| H11 (2) | h1 | +0.0549 | 0.0027 | 0.0005 | +0.0541 to +0.0560 | +0.0547 | +0.20 |
| O1 (1) | os | -0.4632 | 0.0003 | 0.0001 | -0.4633 to -0.4631 | -0.4641 | +0.35 |
| C2 (1) | c | +0.5938 | 0.0003 | 0.0001 | +0.5937 to +0.5939 | +0.5927 | -0.19 |
| O2 (1) | o | -0.4895 | 0.0007 | 0.0001 | -0.4898 to -0.4892 | -0.4870 | +0.63 |
| C3 (1) | ca | -0.1017 | 0.0007 | 0.0001 | -0.1020 to -0.1015 | -0.0996 | +0.65 |
| C4 (2) | ca | -0.0728 | 0.0005 | 0.0001 | -0.0729 to -0.0726 | -0.0720 | +0.59 |
| H4 (2) | ha | +0.1456 | 0.0004 | 0.0001 | +0.1454 to +0.1457 | +0.1465 | +0.73 |
| C5 (2) | ca | -0.0715 | 0.0007 | 0.0001 | -0.0717 to -0.0712 | -0.0720 | -0.20 |
| H5 (2) | ha | +0.1462 | 0.0026 | 0.0005 | +0.1453 to +0.1473 | +0.1440 | -0.52 |
| C6 (1) | ca | -0.1040 | 0.0020 | 0.0004 | -0.1047 to -0.1033 | -0.1046 | -0.12 |
| C9 (1) | c | +0.5963 | 0.0031 | 0.0006 | +0.5952 to +0.5975 | +0.5987 | -0.45 |
| O4 (1) | o | -0.4874 | 0.0040 | 0.0007 | -0.4888 to -0.4860 | -0.4760 | +0.77 |
| O3 (1) | os | -0.4583 | 0.0038 | 0.0007 | -0.4596 to -0.4569 | -0.4711 | -0.23 |
| C10 (1) | c3 | +0.1439 | 0.0119 | 0.0022 | +0.1399 to +0.1483 | +0.1616 | +0.21 |
| H101 (2) | h1 | +0.0756 | 0.0062 | 0.0012 | +0.0734 to +0.0778 | +0.0747 | +0.09 |
| CT (1) | c3 | -0.1031 | 0.0068 | 0.0013 | -0.1059 to -0.1012 | -0.1011 | -0.15 |
| HT1 (3) | hc | +0.0491 | 0.0005 | 0.0001 | +0.0490 to +0.0493 | +0.0494 | +0.25 |
| CN1 (1) | c3 | +0.1431 | 0.0122 | 0.0023 | +0.1388 to +0.1475 | +0.1516 | +0.05 |
| HN1 (2) | h1 | +0.0760 | 0.0057 | 0.0011 | +0.0739 to +0.0780 | +0.0477 | -0.47 |
| ON1 (1) | os | -0.4649 | 0.0044 | 0.0008 | -0.4664 to -0.4633 | -0.4451 | +0.18 |
| CN2 (1) | c | +0.5857 | 0.0034 | 0.0007 | +0.5845 to +0.5870 | +0.5801 | -0.61 |
| ON2 (1) | o | -0.4922 | 0.0044 | 0.0008 | -0.4938 to -0.4906 | -0.4320 | +0.21 |
| CN3 (1) | c3 | -0.1502 | 0.0004 | 0.0001 | -0.1503 to -0.1500 | -0.1901 | +0.39 |
| HN3 (3) | hc | +0.0776 | 0.0006 | 0.0001 | +0.0774 to +0.0778 | +0.0757 | +0.41 |

Charge-shape correlations with |r| >= 0.6: O2 vs Rg (A): r = +0.63; O2 vs closest distant N/O (A): r = +0.66; C3 vs Rg (A): r = +0.65; C3 vs gauche torsions: r = -0.68; C3 vs closest distant N/O (A): r = +0.67; C4 vs closest distant N/O (A): r = +0.65; H4 vs Rg (A): r = +0.73; H4 vs gauche torsions: r = -0.61; H4 vs closest distant N/O (A): r = +0.68; H5 vs closest distant N/O (A): r = -0.64; O4 vs Rg (A): r = +0.77; O4 vs gauche torsions: r = -0.74; O4 vs closest distant N/O (A): r = +0.72; CN2 vs Rg (A): r = -0.61.

## MME

- SMILES `CCOC(=O)c1ccc(C(=O)OCC)cc1`, net charge 0, residue MME, abcg2 / gaff2
- Input issues: none
- Reference conformer: input
- Pool: 13 embedded, 6 pass the filters (500 requested, seed 2026, mmff94s eps = 4.0, window 15.0 kcal/mol)
- Runs: 6, accepted 6, rejected 0
- Accepted conformers: Rg 3.51 to 3.57 A, gauche torsions 0 to 0 (median 0), 10 s per run on average
- Convergence: NO; criteria N >= 20, max SEM <= 0.001 e, max drift over 5 <= 0.0005 e
- Final: max SEM 0.0055 e, max drift 0.00867 e, largest conformational SD 0.0134 e (CP1)
- tleap check: charge -0.000000, errors 0, warnings 0

Ensemble vs reference conformer: max |diff| 0.0174 e (CP1), RMS 0.0052 e

| atom (equivalents) | type | final q (e) | SD | SEM | 95% CI | reference |  r(q, Rg) |
|---|---|---|---|---|---|---|---|
| CP1 (2) | c3 | -0.1185 | 0.0134 | 0.0055 | -0.1271 to -0.1098 | -0.1011 | +0.84 |
| HP1 (6) | hc | +0.0501 | 0.0009 | 0.0004 | +0.0496 to +0.0507 | +0.0490 | -0.84 |
| C1 (2) | c3 | +0.1636 | 0.0000 | 0.0000 | +0.1636 to +0.1636 | +0.1636 | +nan |
| H11 (4) | h1 | +0.0607 | 0.0054 | 0.0022 | +0.0572 to +0.0642 | +0.0537 | -0.84 |
| O1 (2) | os | -0.4630 | 0.0008 | 0.0003 | -0.4634 to -0.4624 | -0.4641 | -0.23 |
| C2 (2) | c | +0.5943 | 0.0005 | 0.0002 | +0.5940 to +0.5947 | +0.5937 | -0.84 |
| O2 (2) | o | -0.4895 | 0.0005 | 0.0002 | -0.4898 to -0.4892 | -0.4890 | -0.48 |
| C3 (2) | ca | -0.1033 | 0.0005 | 0.0002 | -0.1036 to -0.1029 | -0.1026 | +0.84 |
| C4 (4) | ca | -0.0725 | 0.0000 | 0.0000 | -0.0725 to -0.0725 | -0.0725 | +nan |
| H4 (4) | ha | +0.1450 | 0.0000 | 0.0000 | +0.1450 to +0.1450 | +0.1450 | +nan |

Charge-shape correlations with |r| >= 0.6: CP1 vs Rg (A): r = +0.84; HP1 vs Rg (A): r = -0.84; H11 vs Rg (A): r = -0.84; C2 vs Rg (A): r = -0.84; C3 vs Rg (A): r = +0.84.

## Figures

- `abcg2_convergence.png`
- `abcg2_running_mean.png`
- `abcg2_charges.png`
- `abcg2_descriptors.png`
