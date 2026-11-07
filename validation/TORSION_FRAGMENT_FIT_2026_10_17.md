# Fragment torsion fit (2026-10-17)

Archive generation for every torsion the 50-molecule batch filled with a
boss-like fallback (class-alias, terminal-ct-alias, closest-mismatch,
zero-torsion) or with a whole-molecule archive entry.

- Fragments: 37 molecules, 46 bonds (31 rotatable, 15 ring), from
  `validation/missing_torsion_fragments.json` (see `MISSING_TORSION_FRAGMENTS.md`).
- Scans: GFN2-xTB constrained relaxed scans, 1000 points, all within 0.7° of
  target. Rotatable 24 × 15°; planar rings ±30° (13 points); puckered rings
  through inversion (16–34 points).
- QM: B3LYP-D3(BJ)/def2-SVP (density fitting, PySCF), 8 threads, ~100 min.
- Fit: `python -m opls_ff_gen.fragment_fit ... --stage fit` (λ = 0.1).
  977 points (≤ 50 kJ/mol), 293 parameters for 81 fitted keys; 27 keys are
  zero by the heavy-atom convention. Archive: 108 entries (81 fitted, 27 zero),
  replacing the 4 isoniazid whole-molecule entries.
- Work files: `validation/torsion_fits/` (`plan.json`, `F##/scan_*/profile.csv`,
  `fit_report.json`, `run.log`).

## Fit problems found and fixed

1. **Pseudoinverse blow-up.** The planned minimum-norm correction
   (`x = x0 + pinv(A)(r − A x0)`, rcond 1e-3) gave |V| up to 735 kJ/mol.
   Ring torsions and rigid fragments barely move in a scan, so their columns
   are near-flat and collinear. The linear model looked fine, but relaxed MM
   with those terms distorted the geometry (F32 RMSE 1049 kJ/mol). Replaced
   by a ridge toward the prior, min |Ax − r|² + λ|x − x0|²:

   | λ | mean RMSE (relaxed MM) | scans worse than fallback | max \|V\| |
   |---|---|---|---|
   | 0.1 | 2.59 | 4 | 44.1 |
   | 1 | 3.14 | 4 | 39.8 |
   | 3 | 3.94 | 5 | 34.9 |

   Fallback mean: 6.86 kJ/mol. λ = 0.1 chosen.
2. **Key orientation.** Classes were put in canonical order but the SMARTS was
   built in quartet order, so one torsion could get two keys (split
   unknowns, duplicate archive entries). Now `archive.torsion_key` picks one
   orientation for both.
3. **Nitro O.** Formal charge in the context gave the two nitro O different
   keys, and the neutral one also matched the charged one (ambiguous lookup,
   metronidazole failed). Formal charge dropped from the context.
4. **Aromaticity depends on input route.** Caffeine's pyrimidinedione ring is
   non-aromatic from SMILES but aromatic from SDF, so the SDF-built fragment
   keys missed the SMILES-built molecule. Contexts are now built and matched
   on a Kekulized, re-sanitized copy.

## Results per scan (RMSE, kJ/mol, relaxed MM vs QM)

"Fallback" is the workflow topology before the fit; "Zeroed" has the missing
torsions set to zero. ⚠ = fitted worse than fallback.

| Fragment | SMILES | Bond | Type | Points (fit/all) | QM range | Fallback | Zeroed | Fitted |
|---|---|---|---|---|---|---|---|---|
| F01 | `C=C(C)C1CC=CC(=O)C1` | 4-5 | ring | 33/34 | 58.6 | 6.94 | 8.58 | 5.74 |
| F02 | `C=C(C)C1CC=CCC1` | 3-4 | ring | 29/29 | 40.9 | 7.25 | 6.39 | 12.76 ⚠ |
| F03 | `C=CC(C)(O)CC` | 3-6 | rotatable | 24/24 | 18.8 | 5.55 | 5.55 | 1.25 |
| F04 | `C=CC1C(C(=O)OC)C1(C)C` | 4-5 | rotatable | 24/24 | 23.1 | 11.01 | 11.00 | 1.03 |
| F05 | `C=CC1C(C=O)C1(C)C` | 2-9 | rotatable | 24/24 | 13.1 | 1.98 | 3.33 | 1.69 |
| F06 | `C=CC1CC=C(C)C(=O)C1` | 6-7 | ring | 24/24 | 38.0 | 11.33 | 9.07 | 8.83 |
| F07 | `C=CC1CC=C(C)CC1` | 3-4 | ring | 16/16 | 26.6 | 5.12 | 4.15 | 1.22 |
| F08 | `C=CCCC(=C)C` | 3-4 | rotatable | 24/24 | 25.5 | 4.86 | 1.74 | 1.48 |
| F09 | `C=CCc1ccccc1` | 4-5 | rotatable | 24/24 | 12.6 | 2.80 | 4.33 | 0.17 |
| F10 | `CC(=O)OCc1ccccc1` | 4-5 | rotatable | 24/24 | 35.2 | 5.06 | 7.64 | 0.67 |
| F11 | `CC(C(=O)O)c1ccccc1` | 5-7 | rotatable | 24/24 | 11.2 | 1.47 | 2.12 | 1.63 ⚠ |
| F12 | `CC(O)(Cn1cncn1)c1ccccc1` | 3-4 | rotatable | 24/24 | 30.3 | 14.87 | 14.49 | 12.12 |
| F13 | `CC1(C)C(C=O)C1C=C(Cl)Cl` | 3-7 | rotatable | 24/24 | 39.7 | 10.33 | 10.33 | 1.96 |
| F14 | `CC1=CCCCC1=O` | 1-2 | rotatable | 24/24 | 5.7 | 0.22 | 0.18 | 0.18 |
| F14 | `CC1=CCCCC1=O` | 2-7 | ring | 13/13 | 14.6 | 2.06 | 2.80 | 1.46 |
| F15 | `CC1COC(=O)O1` | 5-7 | ring | 13/13 | 14.1 | 2.75 | 2.75 | 1.23 |
| F16 | `CCCC=C(C)C` | 4-5 | rotatable | 24/24 | 19.2 | 2.87 | 2.87 | 0.37 |
| F17 | `CCC[N+](C)(C)[O-]` | 3-4 | rotatable | 24/24 | 34.6 | 6.75 | 6.75 | 7.95 ⚠ |
| F18 | `CCN(CC)CC(=O)Nc1ccccc1` | 6-7 | rotatable | 24/24 | 38.5 | 19.23 | 17.75 | 2.07 |
| F19 | `CC[N+](C)(C)[O-]` | 3-4 | rotatable | 24/24 | 18.8 | 4.29 | 4.29 | 0.27 |
| F20 | `Cc1cc(NS(=O)(=O)c2ccccc2)no1` | 4-5 | rotatable | 24/24 | 33.4 | 13.93 | 11.25 | 6.10 |
| F20 | `Cc1cc(NS(=O)(=O)c2ccccc2)no1` | 5-6 | rotatable | 24/24 | 29.8 | 3.58 | 7.86 | 4.31 ⚠ |
| F21 | `Cc1ccccc1Nc1ccccc1` | 7-8 | rotatable | 24/24 | 32.2 | 6.56 | 8.27 | 1.61 |
| F22 | `Clc1cccc(Cl)c1Nc1ccccc1` | 7-8 | rotatable | 24/24 | 21.9 | 2.85 | 9.96 | 1.51 |
| F23 | `Cn1c(=O)c2c(ncn2C)n(C)c1=O` | 11-12 | rotatable | 24/24 | 1.6 | 0.90 | 0.91 | 0.85 |
| F23 | `Cn1c(=O)c2c(ncn2C)n(C)c1=O` | 2-3 | ring | 13/13 | 13.2 | 3.05 | 3.30 | 0.76 |
| F23 | `Cn1c(=O)c2c(ncn2C)n(C)c1=O` | 2-13 | ring | 13/13 | 11.6 | 2.95 | 3.65 | 0.48 |
| F23 | `Cn1c(=O)c2c(ncn2C)n(C)c1=O` | 11-13 | ring | 13/13 | 15.0 | 1.59 | 3.54 | 0.93 |
| F24 | `Cn1cncc1[N+](=O)[O-]` | 4-5 | rotatable | 24/24 | 41.4 | 15.80 | 20.48 | 2.51 |
| F25 | `NC(=O)N1c2ccccc2C=Cc2ccccc21` | 2-4 | rotatable | 22/24 | 61.3 | 8.46 | 8.46 | 4.62 |
| F26 | `NNC(=O)c1ccncc1` | 2-3 | rotatable | 15/24 | 108.7 | 14.45 | 14.45 | 4.75 |
| F26 | `NNC(=O)c1ccncc1` | 1-2 | rotatable | 24/24 | 34.9 | 8.55 | 8.55 | 2.04 |
| F27 | `O=C(O)Cc1ccccc1` | 2-4 | rotatable | 24/24 | 5.7 | 0.82 | 2.51 | 0.29 |
| F28 | `O=C(OCc1ccccc1)C1CC1` | 4-6 | rotatable | 21/24 | 61.7 | 8.99 | 14.59 | 2.07 |
| F29 | `O=C(c1ccccc1)c1ccccc1` | 5-6 | rotatable | 24/24 | 12.4 | 3.60 | 6.49 | 0.46 |
| F30 | `O=C1OCCO1` | 3-4 | ring | 13/13 | 10.5 | 1.15 | 1.15 | 0.35 |
| F31 | `O=CC=Cc1ccccc1` | 4-5 | rotatable | 24/24 | 27.6 | 9.89 | 9.89 | 0.48 |
| F31 | `O=CC=Cc1ccccc1` | 2-3 | rotatable | 24/24 | 47.7 | 9.88 | 26.29 | 0.66 |
| F32 | `O=CN1c2ccccc2C=Cc2ccccc21` | 9-10 | ring | 17/23 | 111.4 | 18.92 | 18.92 | 5.65 |
| F33 | `O=Cc1ccccc1` | 4-5 | rotatable | 24/24 | 42.6 | 14.13 | 24.50 | 0.34 |
| F34 | `O=S1(=O)CCCC1` | 2-7 | ring | 13/13 | 6.7 | 3.36 | 4.48 | 0.57 |
| F35 | `O=c1ccc2ccccc2o1` | 2-3 | ring | 13/13 | 14.0 | 4.19 | 1.31 | 0.96 |
| F35 | `O=c1ccc2ccccc2o1` | 2-11 | ring | 13/13 | 10.6 | 1.63 | 0.87 | 0.58 |
| F35 | `O=c1ccc2ccccc2o1` | 4-5 | ring | 13/13 | 26.0 | 3.13 | 1.78 | 0.33 |
| F36 | `OC(Cn1cncn1)(Cn1cncn1)c1ccccc1` | 2-3 | rotatable | 22/24 | 52.7 | 12.24 | 15.29 | 5.00 |
| F37 | `OCCn1ccnc1` | 6-7 | rotatable | 24/24 | 28.0 | 14.41 | 14.41 | 6.12 |

Means: rotatable 7.75 → 2.47 (zeroed 9.57); ring 5.03 → 2.79 (zeroed 4.85).
Determined fraction (resolution diagonal) of the 293 parameters: 91 > 0.9,
121 between 0.1 and 0.9, 81 < 0.1 (left near the prior).

Still poor or worse than fallback:
- F02 (ring, 7.25 → 12.76): its keys (CT-CT-CM-HC etc.) are shared with F01,
  F06, F07; one parameter set cannot fit all four. Transferability question
  (deferred).
- F17 amine oxide (6.75 → 7.95), while its homologue F19 fits to 0.40.
- F12 (14.9 → 12.1), F36 and F37 (≈ 6): azole–CH2–C(OH) bonds; likely an
  OH···N hydrogen bond in the nonbonded terms that a torsion cannot fix.
- F11, F20 bond 5-6: < 1 kJ/mol worse.

## Batch rerun with the new archive

`benchmarks/batch_2026_10_17_fragfit/`: 50/50 pass grompp (0 warnings) and
zero-step mdrun. No proper torsion uses a fallback any more (zero-torsion,
closest-mismatch, terminal-ct-alias, class-alias all 0; 119 archive
torsions). Strict mode: 39 molecules (was 24). The 11 still boss-like need it
for atom typing (6) or for a missing angle/bond (5), not for torsions.

## Ring MD check

Gas phase, 300 K, 200 ps, Langevin; mean and SD (degrees) of the archive
ring torsions, previous batch (`batch_2026_10_17_archive`) vs new:

| Molecule | Torsion | Old mean / SD | New mean / SD |
|---|---|---|---|
| coumarin | 3-2-11-10 | 0.2 / 9.4 | −0.3 / 11.1 |
| coumarin | 4-3-2-11 | −0.2 / 7.4 | 0.4 / 8.8 |
| caffeine | 2-13-11-6 | 0.1 / 9.1 | 0.0 / 7.8 |
| caffeine | 5-3-2-13 | 0.1 / 11.2 | 0.1 / 8.8 |
| propylene carbonate | 2-7-5-4 | −2.0 / 8.9 | −3.0 / 10.5 |
| propylene carbonate | 3-4-5-7 | −0.2 / 9.4 | −0.2 / 10.4 |
| sulfolane | 4-2-7-6 | 0.4 / 17.2 | −0.6 / 23.4 |
| sulfolane | 5-4-2-7 | 0.3 / 16.9 | −0.2 / 23.2 |

Planar rings stay planar; sulfolane (puckered, pseudorotating) is softer.

## Update: self-consistent fit (same day, supersedes the fit above)

**F02 was not a transferability conflict.** In the linear model, F02 fitted
jointly with its key-sharing fragments (F01, F06, F07, F14) reached 4.86
kJ/mol, and fitted alone 4.81. The joint fit was already a simultaneous
least-squares fit over every molecule sharing a key. The loss came from the
fitting model: rows were built at QM geometries, while validation relaxes the
MM geometry with the new torsions.

**Fix.** Each cycle relaxes every scan point in MM with the current
parameters, linearizes at the relaxed geometries, and solves the same ridge
problem. The cycle with the lowest relaxed-MM objective is kept. The
objective drops from 40979 (prior) to 6545 by cycle 3, then wobbles by about
±5% without converging, so the best cycle is kept (cycle 7, 6410). The run is
deterministic.

Mean scan RMSE: fallback 6.86 → one-shot fit 2.57 → self-consistent 1.93
kJ/mol. Largest |V| is 44 kJ/mol.

| Fragment | Bond | Type | Fallback | Zeroed | One-shot fit | Self-consistent fit |
|---|---|---|---|---|---|---|
| F01 `C=C(C)C1CC=CC(=O)C1` | 4-5 | ring | 6.94 | 8.58 | 5.74 | 1.56 |
| F02 `C=C(C)C1CC=CCC1` | 3-4 | ring | 7.25 | 6.39 | 12.76 | 2.89 |
| F03 `C=CC(C)(O)CC` | 3-6 | rot | 5.55 | 5.55 | 1.25 | 0.87 |
| F04 `C=CC1C(C(=O)OC)C1(C)C` | 4-5 | rot | 11.01 | 11.00 | 1.03 | 0.79 |
| F05 `C=CC1C(C=O)C1(C)C` | 2-9 | rot | 1.98 | 3.33 | 1.69 | 1.22 |
| F06 `C=CC1CC=C(C)C(=O)C1` | 6-7 | ring | 11.33 | 9.07 | 8.83 | 5.41 |
| F07 `C=CC1CC=C(C)CC1` | 3-4 | ring | 5.12 | 4.15 | 1.22 | 1.00 |
| F08 `C=CCCC(=C)C` | 3-4 | rot | 4.86 | 1.74 | 1.48 | 1.48 |
| F09 `C=CCc1ccccc1` | 4-5 | rot | 2.80 | 4.33 | 0.17 | 0.17 |
| F10 `CC(=O)OCc1ccccc1` | 4-5 | rot | 5.06 | 7.64 | 0.67 | 0.68 |
| F11 `CC(C(=O)O)c1ccccc1` | 5-7 | rot | 1.47 | 2.12 | 1.63 | 1.63 ⚠ |
| F12 `CC(O)(Cn1cncn1)c1ccccc1` | 3-4 | rot | 14.87 | 14.49 | 12.12 | 10.81 |
| F13 `CC1(C)C(C=O)C1C=C(Cl)Cl` | 3-7 | rot | 10.33 | 10.33 | 1.96 | 1.58 |
| F14 `CC1=CCCCC1=O` | 1-2 | rot | 0.22 | 0.18 | 0.18 | 0.24 |
| F14 `CC1=CCCCC1=O` | 2-7 | ring | 2.06 | 2.80 | 1.46 | 1.98 |
| F15 `CC1COC(=O)O1` | 5-7 | ring | 2.75 | 2.75 | 1.23 | 1.22 |
| F16 `CCCC=C(C)C` | 4-5 | rot | 2.87 | 2.87 | 0.37 | 0.35 |
| F17 `CCC[N+](C)(C)[O-]` | 3-4 | rot | 6.75 | 6.75 | 7.95 | 7.93 ⚠ |
| F18 `CCN(CC)CC(=O)Nc1ccccc1` | 6-7 | rot | 19.23 | 17.75 | 2.07 | 2.07 |
| F19 `CC[N+](C)(C)[O-]` | 3-4 | rot | 4.29 | 4.29 | 0.27 | 0.21 |
| F20 `Cc1cc(NS(=O)(=O)c2ccccc2)no1` | 4-5 | rot | 13.93 | 11.25 | 6.10 | 4.67 |
| F20 `Cc1cc(NS(=O)(=O)c2ccccc2)no1` | 5-6 | rot | 3.58 | 7.86 | 4.31 | 5.35 ⚠ |
| F21 `Cc1ccccc1Nc1ccccc1` | 7-8 | rot | 6.56 | 8.27 | 2.29 | 1.91 |
| F22 `Clc1cccc(Cl)c1Nc1ccccc1` | 7-8 | rot | 2.85 | 9.96 | 1.24 | 0.74 |
| F23 `Cn1c(=O)c2c(ncn2C)n(C)c1=O` | 11-12 | rot | 0.90 | 0.91 | 0.85 | 0.93 |
| F23 `Cn1c(=O)c2c(ncn2C)n(C)c1=O` | 2-3 | ring | 3.05 | 3.30 | 0.77 | 0.69 |
| F23 `Cn1c(=O)c2c(ncn2C)n(C)c1=O` | 2-13 | ring | 2.95 | 3.65 | 0.54 | 0.95 |
| F23 `Cn1c(=O)c2c(ncn2C)n(C)c1=O` | 11-13 | ring | 1.59 | 3.54 | 0.86 | 0.51 |
| F24 `Cn1cncc1[N+](=O)[O-]` | 4-5 | rot | 15.80 | 20.48 | 2.26 | 0.34 |
| F25 `NC(=O)N1c2ccccc2C=Cc2ccccc21` | 2-4 | rot | 8.46 | 8.46 | 4.62 | 2.40 |
| F26 `NNC(=O)c1ccncc1` | 2-3 | rot | 14.45 | 14.45 | 4.75 | 1.51 |
| F26 `NNC(=O)c1ccncc1` | 1-2 | rot | 8.55 | 8.55 | 2.04 | 1.67 |
| F27 `O=C(O)Cc1ccccc1` | 2-4 | rot | 0.82 | 2.51 | 0.29 | 0.29 |
| F28 `O=C(OCc1ccccc1)C1CC1` | 4-6 | rot | 8.99 | 14.59 | 2.07 | 2.00 |
| F29 `O=C(c1ccccc1)c1ccccc1` | 5-6 | rot | 3.60 | 6.49 | 0.81 | 0.28 |
| F30 `O=C1OCCO1` | 3-4 | ring | 1.15 | 1.15 | 0.35 | 0.27 |
| F31 `O=CC=Cc1ccccc1` | 4-5 | rot | 9.89 | 9.89 | 0.48 | 0.50 |
| F31 `O=CC=Cc1ccccc1` | 2-3 | rot | 9.88 | 26.29 | 0.66 | 0.66 |
| F32 `O=CN1c2ccccc2C=Cc2ccccc21` | 9-10 | ring | 18.92 | 18.92 | 5.65 | 4.77 |
| F33 `O=Cc1ccccc1` | 4-5 | rot | 14.13 | 24.50 | 0.31 | 0.21 |
| F34 `O=S1(=O)CCCC1` | 2-7 | ring | 3.36 | 4.48 | 0.83 | 0.80 |
| F35 `O=c1ccc2ccccc2o1` | 2-3 | ring | 4.19 | 1.31 | 0.96 | 0.57 |
| F35 `O=c1ccc2ccccc2o1` | 2-11 | ring | 1.63 | 0.87 | 0.58 | 0.81 |
| F35 `O=c1ccc2ccccc2o1` | 4-5 | ring | 3.13 | 1.78 | 0.33 | 0.37 |
| F36 `OC(Cn1cncn1)(Cn1cncn1)c1ccccc1` | 2-3 | rot | 12.24 | 15.29 | 5.00 | 5.21 |
| F37 `OCCn1ccnc1` | 6-7 | rot | 14.41 | 14.41 | 6.12 | 6.12 |

rotatable: fallback 7.75, fitted 2.09
ring: fallback 5.03, fitted 1.59
determined: 87 > 0.9, 127 0.1-0.9, 79 < 0.1 of 293

Still worse than fallback: F11 (+0.2), F17 amine oxide (+1.2), F20 bond 5-6
(+1.8). F12, F36 and F37 (azole–CH2 bonds next to an OH) remain the largest
residuals (10.8, 5.2, 6.1).

Archive and batch: the self-consistent values replace the archive's 108
torsion entries (`benchmarks/batch_2026_10_17_iterfit/`). All 50 molecules
pass grompp and mdrun **in strict mode**: the atom-type archive covers the
borrowed types, and the bonded archive also holds the estimated bonds and
angles. 157 archive terms; no fallback rule is used.

## OH···N hydrogen-bond check (F12, F36, F37)

Scratch scripts: `hbond_diag.py`, `oh_lock.py`, `hb_control.py` (session
scratchpad). Archive not changed.

1. **Geometry along the scans.**
   - F37 (`OCCn1ccnc1`) has no intramolecular hydrogen bond: in QM, H···N is
     4.7–6.3 Å at every point.
   - F12 is hydrogen-bonded from −80° to 40° (H···N 2.0–2.4 Å); F36 is
     hydrogen-bonded at every point (1.8–2.1 Å).
   - In MM relaxation the hydroxyl H was free (locks only covered
     heavy-atom dihedrals). MM then turned the OH toward the ring where QM
     does not: F12 from −170° to −95°, F37 at −58°, −43°, 47° and 62°. Those
     F37 points carry the largest residuals (−11.6, −8.9 kJ/mol).
2. **Lock polar-H dihedrals (H–O–C–X, H–N–C–X) at their QM values and refit**
   (fit stage only, no new QM):

   | Scan | Fitted, OH free | Fitted, OH locked |
   |---|---|---|
   | F37 | 6.12 | **0.77** |
   | F36 | 5.21 | 3.88 |
   | F12 | 10.81 | 10.09 |
   | Mean, 46 scans | 1.93 | 1.72 |

   F26 and F20 bond 5-6 (both with polar H) also improve.
3. **Control (48 new B3LYP-D3BJ/def2-SVP points).** At each drive angle the
   OH is turned anti to the driven CH2 (xTB, both dihedrals constrained),
   and ΔE = E(OH away) − E(original) is compared in QM and MM at the same
   drive angle.

   | | ΔE error (MM − QM), fallback torsions | ΔE error (MM − QM), fitted torsions |
   |---|---|---|
   | F12 | mean +6.1, RMS 14.3 | mean +1.1, RMS 12.7 |
   | F36 | mean −4.9, RMS 8.5 | mean **+19.1**, RMS 21.7 |

   - F12, hydrogen-bonded window (−65° to 25°): MM's hydrogen bond is
     0–7 kJ/mol too weak. The larger errors (±13–15 kJ/mol, and +51 at
     −80°, where the OH-away structure stayed at 2.46 Å) are outside the
     hydrogen bond, in how the OH orientation interacts with the rest of the
     molecule.
   - F36: turning the OH away moves the hydrogen bond to the second triazole
     (H···N 2.2–2.7 Å), so the control switches acceptors rather than
     removing the bond. With fallback torsions MM is about 5 kJ/mol too
     weak; with the fitted torsions the error becomes +19 kJ/mol.
   - The fitted parameters also fit the OH-away conformers worse than the
     scans they were fitted to: F12 12.2 vs 10.9, F36 13.7 vs 7.8 kJ/mol.

**Conclusion.**
- The MM hydrogen bond is only somewhat too weak (about 3–5 kJ/mol), so it is
  not the main error.
- F37's error was a relaxation artifact (free OH rotor).
- In F12/F36 the fit pushes the error into unusually large V1 terms that
  only reproduce the scanned OH arrangement:
  - NA-CT-CT-OH V1 = 37 kJ/mol; the same class quartet in F37's context
    (no hydrogen bond) fits V1 = 0.5;
  - CA-CT-CT-NA V1 = 30 kJ/mol;
  - CT-CT-NA-NB V1 = −25 kJ/mol.

  These terms do not transfer to other OH orientations.

**Decision (user, 2026-10-17): accepted at the current level of accuracy.**
The archive keeps the self-consistent fit as written (OH free during MM
relaxation). Two possible improvements are left open:
1. Lock polar-H dihedrals during relaxation (mean 1.93 → 1.72 kJ/mol; F37
   6.1 → 0.8).
2. Fit the F12/F36 keys over several OH arrangements, so the large V1 terms
   (NA-CT-CT-OH 37, CA-CT-CT-NA 30, CT-CT-NA-NB −25 kJ/mol) stop absorbing
   OH energetics.

Known risk until then: these azole–CH2–C(OH) torsions are tuned to the
scanned OH arrangement and may transfer poorly to other conformers or
molecules.
