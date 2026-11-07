# Torsion fallbacks: parmchk2 substitution and GAFF2 general terms (2026-10-27)

Goal (user): make the on-the-fly QM torsion fit as rare as possible, by rules
for torsions missing from OPLS/2020 that follow AmberTools parmchk2/GAFF.

Data (local, gitignored):
- `validation/torsion_fallbacks/data/`: census, mapping and comparison JSONs,
  the AmberTools 26 `parmchk2.c` the port follows, the GAFF2 frcmods of the
  parmchk2 port check.
- `validation/torsion_fallbacks/oneoff_scripts/`: census and analysis scripts
  of this session (not part of the package).
- `benchmarks/torsion_fallbacks_2026_10_27/`: the three 100-molecule batches
  (`batch_strict` = `--torsion-substitution none`, all missing torsions
  scanned and fitted; `batch_parmchk` = parmchk2 substitution only;
  `batch_gaff2general` = current default) and the 554 benchmark scans
  (`benchmark_table/`, xTB geometries kept for a DFT rescoring).

## 1. Census of missing torsions (100 test molecules, no QM)

46 molecules (48 with any missing term) miss 243 torsion quartets:

| Unique by | Count | In >= 2 molecules |
|---|---:|---:|
| class quartet A-B-C-D | 95 | 31 |
| context key (classes + SMARTS) | 142 | 17 |
| central class pair B-C | 32 | 15 |

64 of the 95 class quartets occur once: the set does not saturate, so a
quartet database would keep needing QM.

## 2. parmchk2 port (`src/opls_ff_gen/parmchk.py`)

- Follows `chk_torsion` of AmberTools 26 `parmchk2.c` with its
  `PARMCHK.DAT` (copied unchanged to `src/opls_ff_gen/data/parmchk/`): steps 2
  (EQUA), 4 (EQUA, general), 5 (CORR, penalty tor + 10 ctor + 10 ctor + tor
  + group/equtype terms), 6 (CORR, general); defaults and ctor mixing as in
  `read_parmchk_parm`.
- Check: with every type of the file and gaff2.dat it reproduces all 272
  substituted torsions of the 100 molecules (source quartet and penalty) as
  `parmchk2 -s 2` writes them.
- User decision: only the uppercase (AMBER) rows, whose names are the OPLS
  classes. Uppercase rows link to lowercase GAFF types (67 EQUA, 214 CORR),
  lowercase rows never to uppercase; CT and C have no uppercase
  correspondence (and c3 none at all), so sp3/carbonyl carbon is never
  substituted, in GAFF either: GAFF covers it with general X-b-c-X terms.
- Result alone: 43 of 95 quartets substituted (114 terms); scans 113 -> 69.

## 3. GAFF2 general terms (`src/opls_ff_gen/gaff2_torsions.py`)

OPLS/2020 has 53 general X-B-C-X terms (aromatic/double bonds), GAFF2 587.
`src/scripts/gaff2_general_torsions.py` maps them to OPLS central pairs:
4,760 molecules (NCI first_5K + test set; 4,075 before the 2026-10-28
precedence entries) typed with antechamber GAFF2 and the OPLS/2020 typer on
the same atoms (69,878 bonds); GAFF2 term as tleap/
parmchk2 would use it; AMBER PK/IDIVF [1 + cos(n phi - gamma)] -> OPLS
V_n = +-2 PK/IDIVF per quartet. A general term keeps only the periodicity of
the bond's local symmetry (n = 3 sp3-sp3, n = 2 conjugated; n <= 3 cancels on
sp2-sp3), so it does not depend on the outer atoms.

User decisions: majority term by bond count for conflicting GAFF2 pairs; the
26 OPLS/2020 general terms are kept; lookup before the type substitution
(parmchk2's order); output `src/opls_ff_gen/data/gaff2_general_torsions.par`
(322 lines since the 2026-10-28 rebuild, 309 before; BOSS format, separate
from the SI); flagged `qm-check` (term used, ITP line and warning say so): 40 terms GAFF2 itself only has by type
substitution, and the OPLS/2020 aryl-N/aryl-carbonyl general terms X-CB-N-X,
X-CB-NA-X, X-CQ-N-X, X-CA-NC-X, X-C-CB-X (GAFF2 differs strongly).

Conversion check on the 26 OPLS/2020 general terms: CA-CA, CB-CB identical
(V2 7.25), CA-CB, CM-CM, CR-NA, CW-NA close; aryl-N and aryl-carbonyl differ.

Batch (xTB, OpenMM): 100/100 pass grompp and zero-step mdrun; 243 torsions
from the general terms, no torsion scans left (24 molecules still fit bonds
or angles); batch 1.5 min (was 14.5).

## 4. QM comparison (relaxed xTB scans)

`src/scripts/compare_general_torsions.py`: the 113 scans of `batch_strict`,
QM-fit topology with the formerly missing torsions set to zero, to the
general batch's values, or fitted; same relaxation protocol and RMSE window
(reproduces the fitted RMSE of `report.json` exactly).

`src/scripts/benchmark_table_torsions.py`: reference, 52 molecules whose bonded
terms all come from the table, 554 symmetry-unique bonds scanned with the
fit's protocol.

| Torsions | Scans | Mean RMSE vs xTB (kJ/mol) | Median |
|---|---:|---:|---:|
| OPLS/2020 table (benchmark) | 554 | 4.48 | 3.51 |
| missing set to zero | 113 | 8.84 | 6.96 |
| GAFF2 general terms | 113 | 9.56 | 7.86 |
| QM fit (fitted to xTB) | 113 | 2.45 | 1.45 |

By central pair (table: benchmark scans; others: the 113 scans):

| B-C | table | zero | general | fit |
|---|---:|---:|---:|---:|
| CT-CT | 3.03 (272) | 9.57 | 10.42 | 3.24 |
| C-CT | 2.85 (19) | 9.31 | 9.73 | 3.51 |
| CA-CT | 2.43 (7) | 5.66 | 6.04 | 2.37 |
| CT-OS | 6.56 (83) | 12.46 | 13.94 | 2.08 |
| C-OS (esters) | 13.66 (14) | 5.07 | 1.69 | 1.97 |
| C-N (amides) | 11.12 (6) | 7.93 | 9.24 | 1.46 |
| C-CA | 7.71 (11) | 0.97 | 3.71 | 0.78 |

Reading:
- sp3 bonds: the table agrees well with xTB (2.4-3.0); bonds whose quartets
  are missing stay at 9-10 with zero or the general term: substituent-specific
  physics an X-B-C-X term cannot carry. Real gap.
- Conjugated polar bonds: OPLS's own ester/amide torsions miss xTB by 11-14
  (worst benchmark scans: aspirin, ethyl benzoate, salicylates, benzocaine
  C-OS); the general terms are within that range or closer to xTB. xTB
  cannot judge these bonds; a DFT reference is needed.
- General term worse than zero on: conjugated C-CM (cinnamaldehyde 3.3 ->
  21.6, coumarin), fused-ring C-N of xanthines (0.2-1.0 -> 6-8), CA-SY
  (dapsone), N-SY (majority term from alkyl/amine sulfonamide N; sulfamethoxazole
  8.4 -> 19.7), N-NT (isoniazid 1.9 -> 10.3). Better on: carbamazepine urea
  C-N (50 -> 21), trimethoprim CQ-NT, acetazolamide azole bonds, chlorzoxazone.
- User decision: do not narrow the general terms (QM fits as rare as possible).

## 5. parmchk2 substitution and a refitted general term (2026-10-28)

**parmchk2 end-atom substitution** (`compare_general_torsions.py
--general-batch batch_parmchk --tag parmchk`, `data/compare_parmchk.json`;
new field `kept_fit` = missing quartets about the bond still fitted): 57 of
the 113 scans substituted, 44 with every missing quartet of the bond.

| On the 44 scans | zero | parmchk2 | GAFF2 general | QM fit |
|---|---:|---:|---:|---:|
| mean / median RMSE | 9.00 / 6.78 | 12.18 / 9.63 | 10.63 / 9.21 | 2.47 / 1.47 |

parmchk2 is better than the general term on 23 scans, worse on 15, worse on
average (CT-CT 14.4 vs 11.0): a specific entry with one end swapped does not
carry the substituent physics either.

**Ceiling of a refitted X-B-C-X term** (`oneoff_scripts/general_ceiling.py`,
`data/general_ceiling.json`): one shared V1-V3 per central pair (32 pairs, 96
parameters) for every missing quartet of the 46 molecules, fitted jointly to
the 113 xTB scans with the fit's protocol (relaxed MM, damped Gauss-Newton,
ridge 0.1). In-sample = best any general term can do here; leave-one-molecule-
out (LOMO, refit linearized at the joint solution, pair seen only in the
held-out molecule -> zero) = what an unseen molecule would get.

| | zero | GAFF2 general | joint (in-sample) | LOMO | QM fit |
|---|---:|---:|---:|---:|---:|
| all 113 scans | 8.84 | 9.56 | 6.11 | 8.99 | 2.45 |
| CT-CT (21 scans, 16 mols) | 9.57 | 10.42 | 8.49 | 9.00 | 3.24 |
| C-CT (11, 11) | 9.31 | 9.73 | 8.07 | 9.12 | 3.51 |
| CA-CT (13, 10) | 5.66 | 6.04 | 5.73 | 7.45 | 2.37 |
| C-N (10, 6) | 7.93 | 9.24 | 5.33 | 10.07 | 1.46 |
| CT-OS (7, 7) | 12.46 | 13.94 | 5.42 | 6.36 | 2.07 |
| C-CM (5, 5) | 2.68 | 10.44 | 2.02 | 4.01 | 0.79 |

Reading:
- Even in-sample, one term per central pair reaches 6.1 (QM fit 2.45); on
  the sp3 gap bonds CT-CT and C-CT it barely moves (9.6 -> 8.5, 9.3 -> 8.1).
  The missing physics is quartet-specific, as section 4 suspected.
- Transferred (LOMO) the refit is no better than zero (8.99 vs 8.84) and
  only slightly better than GAFF2 (9.56). The fitted values absorb molecule-
  specific errors (large V1 per quartet: CT-OS -20.9, C-CT +13.0 kJ/mol).
- Only CT-OS transfers clearly (LOMO 6.4 vs zero 12.5, GAFF2 13.9) and
  consistently across 7 molecules.
- Caveats: xTB reference (wrong for the conjugated polar bonds, section 4);
  small samples per pair; LOMO linearized, not re-iterated.
- Conclusion: re-parameterizing the ~300 general terms against xTB does not
  pay; a targeted CT-OS term (and DFT for the conjugated pairs) is the most a
  refit could add.

## 6. GAFF2 in its own model against the same xTB scans (2026-10-28)

`oneoff_scripts/gaff2_vs_xtb.py`, `data/gaff2_vs_xtb.json`; GAFF2 builds in
`benchmarks/torsion_fallbacks_2026_10_27/gaff2/`. The 123 scans with central
pair CT-OS, C-OS or C-N (table benchmark and gap scans, 50 molecules); GAFF2
topology from the same atoms and coordinates (antechamber -at gaff2 -c bcc,
parmchk2 -s 2, tleap, ParmEd -> GROMACS; Amber 1-4 scaling); same relaxation
and RMSE protocol. Each bond labelled by GAFF2's own torsion source about it
(`general` = all quartets X-b-c-X, `specific`, `mixed`, `parmchk2`).

| Set, pair (GAFF2 source) | scans | GAFF2 | OPLS table | OPLS + GAFF2 general | OPLS zero | QM fit |
|---|---:|---:|---:|---:|---:|---:|
| gap CT-OS (general) | 6 | 11.51 | | 14.57 | 13.17 | 1.56 |
| gap C-N (general) | 5 | 7.83 | | 9.55 | 10.34 | 2.15 |
| gap C-N (all) | 10 | 6.69 | | 9.24 | 7.93 | 1.46 |
| gap C-OS (mixed) | 3 | 4.44 | | 1.69 | 5.06 | 1.97 |
| table CT-OS (all) | 83 | 7.18 | 6.56 | | | |
| table CT-OS (general) | 9 | 6.13 | 4.53 | | | |
| table CT-OS (mixed: ester O-CH2) | 9 | 11.90 | 9.41 | | | |
| table C-OS (specific: esters) | 13 | 6.20 | 13.38 | | | |
| table C-N (amides) | 6 | 13.43 | 11.12 | | | |

Reading:
- GAFF2's general terms miss xTB in GAFF2's own model too: 11.5 on CT-OS,
  7.8 on C-N; carried into OPLS they miss by 14.6 and 9.6. So most of the
  deviation belongs to the general term itself; carrying it into OPLS (charges,
  1-4 scaling) adds about 2-3 kJ/mol.
- Ester alkyl O-CH2 (CT-OS; benzyl, ethyl, butyl, hexyl esters, triacetin):
  xTB barrier 12-25 kJ/mol, GAFF2 35-67, OPLS (table or general) 35-57. The
  two force fields, derived independently, agree with each other and give
  2-3 times xTB's barrier. That points to xTB as the outlier on this bond; DFT
  on a few of these scans would settle it.
- Ester C(=O)-O (C-OS): GAFF2's specific terms follow xTB (6.2; barriers
  within a few kJ/mol, e.g. ethyl salicylate 54.6 vs 54.7), while OPLS/2020
  is stiffer (13.4; aspirin 80 vs 46, benzocaine 73 vs 55). Here GAFF2 backs
  xTB: the OPLS table ester terms (with CM1A-LBCC charges) look too stiff,
  which changes the reading of section 4 (observation only: the SI values
  are used as they are).
- Amides C-N: no common pattern; both force fields miss xTB by 11-13 on the
  table amides (GAFF2 overshoots salicylamide and lauramide MEA, OPLS
  overshoots acetaminophen and lauramide DEA); fused xanthine C-N bonds:
  xTB ~10, GAFF2 14-22, OPLS 19-25.

**DFT check of the ester O-CH2 reading** (`oneoff_scripts/ester_dft_check.py`,
`data/ester_dft.json`): B3LYP-D3(BJ)/def2-SVP single points on the 72 xTB
scan geometries of three scans (about 4 min on 10 threads).

| Scan | Barrier DFT | xTB | GAFF2 | OPLS | RMSE vs DFT: xTB | GAFF2 | OPLS |
|---|---:|---:|---:|---:|---:|---:|---:|
| n-butyl acetate O-CH2 (table) | 29.5 | 18.6 | 42.8 | 34.4 | 4.9 | 5.9 | 3.8 |
| ethyl benzoate O-CH2 (table) | 27.7 | 20.4 | 43.2 | 38.0 | 3.4 | 7.5 | 7.0 |
| benzyl acetate O-CH2 (gap, GAFF2 general) | 35.8 | 16.0 | 51.1 | 49.0 | 9.1 | 6.3 | 8.4 |

The reading is only partly confirmed. xTB's barrier is too low (by 7-20
kJ/mol; 45-75 % of DFT), but both force fields overshoot as well (GAFF2 by
13-15, OPLS by 5-13): DFT lies between them, and the 2-3x gap was about half
xTB and half the force fields. GAFF2 also has a spurious bump near +-125 deg
(10-13 kJ/mol, DFT 3-7). By RMSE no reference wins clearly: xTB is closest on
the two table scans (3.4-4.9), worst on benzyl acetate (9.1). Caveats: single
points on xTB geometries, not DFT-relaxed (usually raises the barrier a
little); small basis. OPLS/2020 table values stay as they are (user).

## 7. Fragment refit of the general terms, step 1: targets and fragments (2026-10-28)

User decision: refit OPLS-consistent X-B-C-X terms on small fragments, as
GAFF did (Wang et al. 2004: one generic term per model molecule, Parmscan
fit to MP4//MP2 profiles), with **xTB** as the reference (no DFT: compute
budget, DFT not clearly better).

`oneoff_scripts/fragment_targets.py` and `select_fragments.py`; output
`data/fragments/` (`nci_missing.json`, `targets.json`, `selected.json/.csv`).
NCI 5K molecules of `gaff2_types.json` minus the 100 test molecules (4,670;
703 skipped: OPLS typing ambiguities, embedding), strict lookup (no
fallbacks): 13,264 bonds with missing quartets, 87 pairs of OPLS classes
and many with a GAFF2-typed central atom (6,583 bonds; typing gaps, not
targeted). Fragment per non-ring missing bond: heavy atoms within two bonds
of the central atoms, whole aromatic ring systems, atoms on non-single
bonds, H caps, at most 18 heavy atoms; kept only if, retyped alone, all
quartets about the bond keep their classes and missing status.

Targets: the 23 OPLS-class pairs with >= 40 missing bonds and a non-ring
fragment (CT-CT 1,354 bonds, C-CT 709, CT-OS 536, C-N 499, CA-CT 420, CT-OH,
C-OS, C-CM, CT-NT, CA-NT, CT-S, C-CA, C-OH, CM-CT, CA-CM, CO-OS, CA-S, CA-SY,
CT-N, OS-P, CA-OS, CO-CT, N-NT); ring-only pairs (CW-SA, CB-CT, CW-OA) and
CZ-S (linear) left out. 150 fragments (up to 8 per pair, neutral, most
frequent first, each adding a new class quartet). Fragments with 8 picks
represent 72-98 % of the fragment bonds of most pairs, but only 38-46 % for
CT-CT, C-CT and C-N (113, 71 and 74 class quartets). Ring bonds are 41 % of
the CT-CT gaps, 45 % of CT-OS: the chain-fitted term would also be applied to
them.

## 8. Fragment refit, steps 2-3: xTB scans, fit, validation (2026-10-28)

User decisions: 24 fragments for CT-CT, C-CT and C-N (8 for the others);
ring bonds in the training set (ring fragments in proportion to each pair's
ring share). 233 fragments (26 pairs, CW-SA, CB-CT and CW-OA now included
through ring fragments), `data/fragments/fragments_{0,1,2}.csv`; strict xTB
batches `benchmarks/.../fragments_{0,1,2}` (about 3-40 s per fragment):
221 OK, 12 failed (xTB constrained optimization diverged at one scan point,
|grad| > 500, e.g. a scan turning a C=N double bond; the batch drops the whole
fragment). 3 more fragments left out because their own on-the-fly fit missed
a scan by 230-2,200 kJ/mol (broken scans; 98.6 % of the objective).

Fit (`oneoff_scripts/fit_general_fragments.py`, `data/refit/terms.json`):
218 fragments, 653 scans, one V1-V3 per central pair (97 pairs, incl.
incidental ones), chain and ring bonds together, fit protocol as in the
ceiling test; objective 889k -> 442k in 8 iterations. Several terms carry
large, poorly determined V1/V2 (C-CA V1 43.8, CT-c3 V2 -74, CB-CT V2 -38).

Validation (`data/refit/validation.json`): the 113 xTB scans of the test
molecules (none in training); pairs fitted on >= 2 fragments take the refit,
other missing quartets the current default.

| Scans | zero | GAFF2 general (current) | refit | QM fit |
|---|---:|---:|---:|---:|
| all 113 | 8.84 / 6.96 | 9.56 / 7.86 | 8.19 / 6.11 | 2.45 / 1.45 |
| 94 with a refitted pair | 8.43 / 6.36 | 9.34 / 7.74 | 7.74 / 5.52 | 2.26 / 1.49 |

(mean / median RMSE, kJ/mol). Against GAFF2 general: better on 53 scans,
worse on 26 (> 0.5 kJ/mol). By pair (mean; zero / GAFF2 / refit):
CT-OS 12.5 / 13.9 / 5.9 (7 scans), C-N 7.9 / 9.2 / 6.9 (10), CT-NT 15.9 /
17.6 / 13.6, CT-OH 12.5 / 12.2 / 10.1, CW-SA 19.7 / 13.3 / 5.6, CW-nb 19.5 /
10.3 / 5.8, C-CM 2.7 / 10.4 / 4.8, C-CT 9.3 / 9.7 / 9.1 (11), CA-CT 5.7 /
6.0 / 5.9 (13), **CT-CT 9.6 / 10.4 / 11.6 (21, worse)**.

Reading: the refit beats the current GAFF2 terms by 1.4 kJ/mol overall
(1.6 where used) and zero by 0.7; it falls short of the 2 kJ/mol bar set
before the run. Clear gains on CT-OS, the azole ring bonds and C-N; no gain
on C-CT/CA-CT; a loss on CT-CT, the most frequent pair.

Stronger ridge (`--ridge`, against zero; `data/refit_ridge{100,1000}/`), same
fragments and validation, mean RMSE:

| ridge | all 113 | 94 refitted | better / worse vs GAFF2 | CT-OS | C-N | CT-CT | largest V |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 8.19 | 7.74 | 53 / 26 | 5.9 | 6.9 | 11.6 | 84 |
| 100 | 8.74 | 8.32 | 54 / 23 | 9.3 | 7.7 | 10.5 | < 15 |
| 1000 | 8.81 | 8.46 | 57 / 11 | 11.8 | 8.2 | 9.8 | < 2 |

Reading: ridge 100 already shrinks every term below 15 kJ/mol (CT-OS V1 -17.9
-> -6.1) and the gain over zero (8.84) disappears; ridge 1000 is zero in all
but name. The gains at ridge 0.1 come from the large terms, so taming them
costs the gains; only CT-CT improves with shrinkage. Across all three,
setting missing terms to zero beats the current GAFF2 general terms.

**Adopted (user, 2026-10-28):** the ridge-0.1 terms for the six pairs with a
clear gain and sane terms: CT-OS, C-N, CT-NT, CT-OH, CW-SA, CW-nb
(`xtb-refit` lines in `src/opls_ff_gen/data/gaff2_general_torsions.par`, from
`XTB_REFIT` in `src/scripts/gaff2_general_torsions.py`; CW-SA loses its
`qm-check` flag). All other pairs keep the GAFF2 term. Check
(`data/compare_refit6.json`: the test molecules rebuilt with the new file,
`batch_refit`, scored by `src/scripts/compare_general_torsions.py --tag refit6`):

| 113 test scans | zero | GAFF2 general (before) | adopted | QM fit |
|---|---:|---:|---:|---:|
| mean / median | 8.84 / 6.96 | 9.56 / 7.86 | 8.43 / 6.34 | 2.45 / 1.45 |

30 scans better than before, 4 worse (> 0.5 kJ/mol). By pair (before ->
after): CT-OS 13.9 -> 5.9, C-N 9.2 -> 6.4, CT-NT 17.6 -> 14.4, CT-OH 12.2 ->
9.7, CW-SA 13.3 -> 5.6, CW-nb 10.3 -> 5.8. Caveat: the six pairs were chosen
on these test scans (selection on the test set); the CT-OH, CW-SA and CW-nb
numbers rest on 2 scans each.

## 9. Do bonded terms still need QM? 100 new NCI molecules (2026-10-28)

`oneoff_scripts/sample_qm_need.py`, `data/sample_qm_need.json`: 100 molecules
drawn at random (seed 20261028) from the NCI 5K set, excluding the 100 test
molecules and the 621 parents of the refit fragments; typed as the workflow
does and run through the default bonded lookup (no QM run).

- 19 of 100 cannot be typed at all: 18 "Found multiple types for atom N"
  (conflicting OPLS SMARTS matches; 694 of 4,670 NCI molecules in the census
  of section 7 fail the same way), 1 RDKit embedding failure.
- Torsions: none of the 81 typed molecules needs QM.
- Bonds and angles: 42 of 81 need QM (220 bonds, 500 angles; 72 bond and 206
  angle keys). 11 terms came from the recorded QM table.
- GAFF2-typed atoms cause most of it: 686 of the 720 terms involve a GAFF2
  class (no OPLS/2020 bond/angle for it); 27 molecules need QM only for such
  terms, 7 for both, 8 only for OPLS-class gaps (3 bonds, 31 angles; CB-NT,
  CB-SA bonds, e.g. CM-CM-CM, CM-C-OH angles).

Follow-up (user: GAFF2 bonded terms from an extended SMARTS typer, no
antechamber at runtime; xTB for the OPLS-class gaps):
- Full GAFF2 types (`src/opls_ff_gen/gaff2_bonded.py`, `src/scripts/build_gaff2_bonded.py`)
  vs antechamber on NCI 5K: gaff2.dat term identical for 99.6 % of 147,983
  bonds and 99.5 % of 253,889 angles; differences in 37 of 4,770 molecules,
  mostly cationic dyes (antechamber re-perceives their resonance).
- Same sample (`data/sample_qm_need_gaff2.json`): gaff2.dat supplies 686
  terms; 16 of 81 molecules still need QM (34 OPLS-class terms, 2 angles
  gaff2.dat lacks: CM-c-Cl).
- xTB fits of these 16 (`benchmarks/.../batch_sample_gaps`) recorded in
  `qm_bonded_terms` (now 40 molecules, 89 terms): `data/sample_qm_need_final.json`,
  0 of 81 need QM. The 100 test molecules still need none (`batch_qmtable2`,
  100/100 grompp + mdrun; acetazolamide now takes 9 terms from gaff2.dat).

## 10. Fallback for missing bonds and angles, QM fit optional (2026-10-28)

User decisions: class substitution with a default substitute as the final
fallback for bonds and angles; the default keyed by the atom (element,
neighbour count, aromaticity, bonding); an empirical formula as the last
resort; the recorded QM table stays before the substitution; the run-time
QM fit becomes an option (`--qm`, default `none`).

- Coverage first: parmchk2's own substitution (uppercase PARMCHK.DAT types)
  finds a substitute for only 46 % of the open angle keys of the 559 NCI
  molecules that still needed QM (classes CF, SZ, NO, SY, SA, CS, OA, OY, CO
  have no PARMCHK.DAT entry; C and CT no correspondences).
- Against the xTB (Seminario) fits, same terms (328 molecules first, then
  882): for angles every class substitution is further from xTB than the
  empirical rule (parmchk2-corr 2.85 vs 1.90 degrees, default 2.67 vs 1.99,
  same-element 4.80 vs 3.48; k likewise); for bonds the substitute's b0 is
  closer (0.016-0.020 vs 0.031-0.050 A) but the empirical k is (3-24 % vs
  33-106 %). Adopted: bonds b0 from the substitute with the empirical k,
  angles empirical (`src/opls_ff_gen/bonded_substitution.py`). Final, 882
  molecules: angles 2.3 degrees / k 26 %; substituted bonds 0.011-0.019 A.
- The 559 ran with xTB (`benchmarks/.../batch_nci_qm`, 554 OK; 5 fail in the
  xTB geometry step, bonds formed) and are recorded: `qm_bonded_terms` 917
  molecules, 597 terms.
- Check with the defaults (no QM): all 4,761 embeddable molecules of the
  in-scope set (NCI 5K typable by antechamber + 100 test molecules) get every
  bonded term; 18 bond/angle terms come from the fallback (the 5 failed
  molecules). 9 molecules fail before typing (RDKit conformer embedding).
  grompp + zero-step mdrun: 20 random molecules OK (`batch_fallback_noqm`);
  11 of the 559 with the recorded table off, so the fallback supplies every
  gap (1-10 terms each), all OK (`batch_fallback_only`).

## Next steps (not started)

1. Done 2026-10-28, section 5 (also the general-term refit ceiling).
2. DFT spot check: B3LYP-D3(BJ)/def2-SVP single points on the kept xTB scan
   geometries for some ester C-OS / amide C-N (table and general) and CT-CT
   gap bonds, to decide whether xTB or OPLS is off on conjugated bonds.
3. The census on the larger set (NCI 5K) for the saturation curve of
   missing central pairs.

## 11. xTB fits for the 2026-10-29 borrowed types (2026-10-29)

User decision: the bonded terms of the 2026-10-29 borrowed OPLS types
(pyrazole, triazine, epoxide, R4N+, ... ; archive `borrowed` rules dated
2026-10-29) that the OPLS/2020 table lacks are fitted to xTB and recorded:
bonds and angles in `qm_bonded_terms`, torsions as `xtb-refit` general
terms (the route of section 8; the substitution stays after the general
terms, section 5). Data: `benchmarks/borrowed_fit_2026_10_29/`.

- Need (NCI 5K + 100 new + Si/halogen set, default lookup): 255 molecules
  with a fallback term on such an atom; 72 bond/angle and 42 torsion class
  patterns. Set cover: 34 molecules (`cover34.csv`).
- `cover34_xtb` (`--qm xtb`): 33 OK, first_5K:2685 failed (xTB constrained
  optimization at one scan point). Recorded with `qm_bonded_terms.py`: 946
  molecules, 694 terms (was 917 / 597; 98 new keys, 30 existing keys moved by
  at most 4.0 deg / 0.8 pm). Bond/angle fallbacks on the borrowed atoms:
  255 molecules -> 7 (8 angles: CA-CT-N3 and CM-CT-N3 in 5-rings, recorded
  only open-chain; first_5K:2685; Br-CT-CB).
- Fix found on the way: the recorded table is keyed by the fitted topology's
  classes (an analog atom's OPLS class), so `bonded.py` now looks it up with
  the table classes first, then with the GAFF2 classes of analog atoms.
- Torsions: of the 42 patterns, only the central pairs CA-CU (141 ring, 2
  chain bonds; fused pyrazolo[3,4-d]pyrimidines) and CA-CQ (8 chain bonds;
  aryl-1,3,5-triazines) have fragments outside the validation molecules
  (`fragment_targets.py`, `data/fragments_1029/`); CS-CZ bonds are to a linear
  nitrile C. 9 fragments (7 CA-CU, 2 CA-CQ; `fragments/`), all scans OK; joint
  fit as in section 8 (`data/refit_1029/`): objective 11,098 -> 753.
- Held-out scans (`cover34_xtb` vs `cover34_default`, RMSE kJ/mol, zero /
  current default / refit / per-molecule fit):

  | Pair | Scan | zero | default | refit | fit |
  |---|---|---:|---:|---:|---:|
  | CA-CQ | first_5K:3243 (phenyl-triazine) | 9.21 | 52.61 | 5.68 | 0.11 |
  | CA-CU ring | first_5K:4875 | 13.26 | 11.63 | 1.07 | 3.60 |
  | CA-CU ring | first_5K:1643 | 11.03 | 9.18 | 0.81 | 4.20 |
  | CA-CU chain | first_5K:4090 (3-pyridylpyrazole) | 10.07 | 5.08 | 94.66 | 1.15 |

  The CA-CQ default was the in-ring aromatic CA-CA-CA-N? by parmchk2
  substitution (penalty 290.5). Adopted: X-CA-CQ-X 0.133 / 1.685 / 0.183
  kcal/mol (`gaff2_general_torsions.par` line 9323, `XTB_REFIT_NEW` in
  `src/scripts/gaff2_general_torsions.py`). CA-CU (V1-V3 -0.908 / 14.059 /
  -6.145 kcal/mol, ring-trained) is wrong for the biaryl CA-CU bond, so one
  X-CA-CU-X term cannot serve both. User decision: a ring-only general term
  (line 9324, comment `ring bonds only`; `gaff2_torsions.py` uses it only for
  torsions about a ring bond, before any term for every bond of the pair);
  chain CA-CU bonds keep the parmchk2 substitution ??-CM-CU-??. Not adopted:
  the incidental CA-NA (fitted on fused ring bonds; 40-46 vs 10-15 on N-aryl
  chain bonds).
