# Benchmark against LigParGen 2.1: DEET, icaridin, transfluthrin, C12E6

Workflow outputs are in `benchmarks/NAME/`; the LigParGen references are
`molecules_to_test/NAME.itp`. Both use the atom order of
`molecules_to_test/NAME.sdf`. Regenerate the comparison with:

```bash
for m in deet icaridin hexaoxyethylene_dodecyl_ether; do
  conda run --no-capture-output -n openMM python -m opls_ff_gen \
    --sdf molecules_to_test/$m.sdf --name $m --output-dir benchmarks/$m
done
conda run --no-capture-output -n openMM python -m opls_ff_gen \
  --sdf molecules_to_test/transfluthrin.sdf --name transfluthrin \
  --output-dir benchmarks/transfluthrin --bonded-fallback boss-like
conda run -n openMM python src/scripts/benchmark_ligpargen.py \
  deet icaridin transfluthrin hexaoxyethylene_dodecyl_ether
```

Each benchmark directory holds `benchmark_atoms.csv` (per atom),
`benchmark_terms.csv` (per interaction), `benchmark_term_families.csv`
(differences grouped by class pattern), and `benchmark_summary.json`.
Parameters are compared at LigParGen's printed precision (one unit of the
last printed digit).

## Summary

| Molecule | Atoms | Class match | LJ match | Charge MAE (e) | Max \|Δq\| (e) | Bonds | Angles | Propers | Impropers | 1-4 pairs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DEET | 31 | 31 | 30 | 0.0007 | 0.0071 | 31/31 | 54/54 | 62/68 | 1/8 | 65/65 |
| Icaridin | 39 | 39 | 37 | 0.0001 | 0.0007 | 39/39 | 74/74 | 88/105 | 1/2 | 102/102 |
| Transfluthrin† | 35 | 35 | 22 | 0.0008 | 0.0032 | 33/36 | 58/64 | 68/81 | 1/9 | 78/78 |
| C12E6 | 81 | 81 | 81 | 0.0001 | 0.0004 | 80/80 | 151/151 | 181/192 | 0/0 | 192/192 |

Term counts are `matching parameters / LigParGen terms`. For all four
molecules, the *sets* of bonds, angles, proper torsions, impropers, and pairs
are identical to LigParGen; every count difference is a parameter-value
difference. Ester carbonyl `opls_465`/`opls_466` use classes `C`/`O`, as in
LigParGen (changed 2026-10-17; before that the two transfluthrin ester atoms
were `C_2`/`O_2` class mismatches and propers were 67/81, see
`validation/BATCH_2026_10_17.md`). †Transfluthrin is generated with
`--bonded-fallback boss-like`; the other three are fully tabulated and give
identical parameters in `strict` and `boss-like` modes. All four topologies
pass `gmx grompp` (GROMACS 2026.3) with zero warnings.

## Workflow changes needed to run these molecules

Foyer 1.2.0 could not type icaridin or transfluthrin. `src/opls_ff_gen/assignment.py`
now patches three gaps, in the same way as the existing `FOYER_TYPE_PRECEDENCE`:

| Problem | Molecule, atom | Fix |
|---|---|---|
| Carbamate carbonyl matches both amide `opls_235` and ester `opls_465` (and O: `opls_236`/`opls_466`) | icaridin C6, O7 | Precedence: amide types win (LigParGen also uses classes C/O). Alkoxy O stays ester `opls_467` |
| `opls_227` (=CCl2) SMARTS only matches Cl2C=CH2; Cl atoms then fall back to chloride ion `opls_401` | transfluthrin C20, Cl21, Cl22 | Broadened `opls_227` SMARTS to any `C(Cl)(Cl)=C` |
| Tertiary-amide N-CHR2 carbon `opls_246` has no SMARTS and no class in Foyer | icaridin C13 | Added SMARTS; class `CT` (as Foyer does for `opls_245` and LigParGen does) |

The unit tests (29) pass.

## Transfluthrin: boss-like bonded fallback

GROMACS `ffbonded.itp` lacks 19 bonded class patterns for transfluthrin.
The cyclopropane class `CY` is sparsely parameterized next to ester and alkene
groups, and `Cl-CM-Cl` is absent. In the default `strict` mode the run fails
and lists all 19 patterns with their atoms. With `--bonded-fallback boss-like`
these rules are tried, in order, only after the table lookup fails:

1. `class-alias`: retry with generic classes (`C_2`/`C_3`→`C`, `O_2`/`O_3`→`O`,
   `CT_2`/`CT_3`/`CY`→`CT`), fewest substitutions first.
2. `central-class-average`: for an angle, the mean θ0 and k of all angle
   entries with the same central class. This is the rule documented in the
   BOSS/MCPRO manuals.
3. `zero-torsion`: no torsional term. A warning is raised if every torsion
   about a bond ends up zero. None does in transfluthrin.

There is no estimate for a missing bond; a bond that no alias resolves still
fails. BOSS's "closest mismatch" torsion matching is deliberately not
reproduced. Every fallback term is marked in the ITP (`; <rule>`) and the
topology JSON.

The 29 terms produced by these rules compare with LigParGen as follows:

| Rule | Terms | Match LigParGen | Differ |
|---|---:|---:|---:|
| `class-alias` | 25 | 17 | 8 |
| `central-class-average` | 1 | 0 | 1 |
| `zero-torsion` | 3 | 2 | 1 |

The ten differences:

| Classes | Rule | Workflow | LigParGen |
|---|---|---|---|
| `Cl-CM-Cl` angle | central-class-average | 119.78°, 542.3 | 118.77°, 482.7 |
| `CT-CY-CT` angle | class-alias (`CT-CT-CT`) | 112.7°, 488.3 | 114.3°, 292.9 |
| `CY-CM-HC` angle | class-alias (`CT-CM-HC`) | 117.0°, 292.9 | 135.0°, 292.9 |
| `CY-CY-C-O`, `CY-CY-CM-CM`, `CM-CM-CY-HC` torsions (4) | class-alias | non-zero CT analogues | 0 |
| `HC-CY-CM-HC` torsion | class-alias (`HC-CT-CM-HC`) | 0.665 1.996 0 -2.661 | 58.576 0 -58.576 |
| `CA-CT-OS-C` torsion | zero-torsion | 0 | -2.197 5.201 0.527 -3.531 (the `C-OS-CT-CT` values) |

LigParGen's zeros for the `CY` torsions suggest BOSS has cyclopropane-specific
entries or a zero default for them. Its `HC-CY-CM-HC` value has the shape of a
C=C double-bond torsion, consistent with the closest-mismatch matching that is
not reproduced here.

Differences in transfluthrin that do not come from the fallback rules:

- The `CY-CY-CY` ring angle is 60° in GROMACS and 83° in LigParGen.
- The `CY-CY` bond is 0.1509 nm in GROMACS and 0.152 nm in LigParGen.
- The CM impropers (`improper_Z_CM_X_Y`) are 62.76 in GROMACS and 10.46 in
  LigParGen.
- The aromatic impropers differ as in section 4.

**Enumeration fix.** Before this change the workflow also wrote 12 "proper
torsions" inside the cyclopropane ring, whose outer atom is bonded to both
central atoms, so the torsion spans a 1-3 pair. These are now omitted, and
the torsion set matches LigParGen (81/81).

## Differences

### 1. Partial charges: symmetry averaging (DEET, transfluthrin)

C12E6 and icaridin reproduce LigParGen within 7e-4 e. The larger
deviations are equal and opposite pairs. RDKit graph symmetry averages atoms
that LigParGen keeps distinct:

- DEET: the two N-CH2 carbons C11/C13 (LigParGen 0.1400 / 0.1542, workflow
  0.1471 for both) and their H atoms. These are cis and trans to the carbonyl
  O across the amide bond.
- Transfluthrin: ring C10/C14 and C11/C13 (±0.003 e), their F atoms
  (±0.0015 e), and the diastereotopic gem-dimethyl carbons C1/C23 (±0.002 e).

The pair means are preserved. LigParGen's averaging is narrower than graph
symmetry when E/Z or diastereotopic environments are involved.

### 2. Lennard-Jones parameters (Foyer/GROMACS vs LigParGen)

| Atom | Foyer type | Workflow σ (nm) / ε (kJ/mol) | LigParGen σ / ε | Molecules |
|---|---|---|---|---|
| Amide / carbamate / ester C=O carbon | `opls_235`, `opls_465` | 0.375 / 0.43932 | 0.355 / 0.29288 | DEET, icaridin, transfluthrin |
| Ester alkoxy O | `opls_467` | 0.300 / 0.71128 | 0.290 / 0.58576 | icaridin, transfluthrin |
| Cyclopropane C | `opls_712`, `opls_713` | 0.350 / 0.276144 | 0.345 / 0.347272 | transfluthrin |
| Aromatic F (polyfluorobenzene) | `opls_728` | 0.285 / 0.255224 | 0.290 / 0.25104 | transfluthrin |
| H on ester CH2 / Cα of carbonyl | `opls_469`, `opls_282` | 0.242 / 0.06276 | 0.250 / 0.12552 | transfluthrin |
| Alkene H on CH=CCl2 | `opls_144` | 0.242 / 0.12552 | 0.250 / 0.12552 | transfluthrin |

The bonded classes agree. LigParGen applies generic LJ for these types
instead of the specialized OPLS-AA values in Foyer/GROMACS.

### 3. Proper torsions (GROMACS ffbonded.itp vs BOSS)

RB coefficients C0..C5 in kJ/mol:

| Classes | Count (DEET / icaridin / C12E6) | Workflow (GROMACS) | LigParGen |
|---|---|---|---|
| `C-N-CT-CT` | 2 / 3 / – | -4.707 2.920 1.787 0 0 0 | -0.176 3.548 -0.544 -2.828 0 0 |
| `C-N-CT-HC` | 4 / 3 / – | -0.291 -0.872 0 1.163 0 0 | 0 0 0 0 0 0 |
| `CT-CT-CT-CT` | – / 5 / 9 | 2.929 -1.464 0.209 -1.674 0 0 | 1.360 -0.523 0.837 -1.674 0 0 |
| `CT-CT-CT-N` | – / 3 / – | 5.487 0.027 0 -5.515 0 0 | 4.372 -2.033 -0.665 -1.674 0 0 |
| `HC-CT-OH-HO` | – / 2 / 2 | 0.941 2.824 0 -3.766 0 0 | 0.736 2.209 0 -2.946 0 0 |
| `CT-OS-C-O` | – / 1 / – | 20.920 0 -20.920 0 0 0 | 21.439 0 -21.439 0 0 0 |

`HC-CT-OH-HO` is the same difference already seen in ethanol and
phenoxyethanol.

**`CT-CT-N-CT` (resolved).** The generic GROMACS entry for this pattern
(ffbonded.itp line 1427) is commented "from Pro CG-CD-N-CA". It is a
proline-ring fit that favours φ = 0° by 35 kJ/mol, which is unphysical for
open-chain amides. `opls_ff_gen.bonded.FIVE_RING_ONLY_TORSIONS` now applies
it only when all four atoms lie in one five-membered ring. Otherwise the
tertiary-amide entry "Pro CD-N-CA-CB" (line 1517, `CT_3 N CT_2 CT`;
6.872 -9.943 3.071 kJ/mol) is used. That entry prefers trans and penalizes
cis by 20 kJ/mol. Both values come from the GROMACS file. All five
`CT-CT-N-CT` torsions in DEET (2) and icaridin (3, in and next to the
piperidine ring) now match LigParGen. `CA-CT-N-CT` (line 1428) is the same
proline fit, but has no tertiary-amide counterpart in the file. It is not
used by these molecules and is left unchanged.

### 4. Improper torsions

| Center | Workflow k (kJ/mol) | LigParGen k | Count |
|---|---|---|---|
| Aromatic C (`improper_Z_CA_X_Y`) | 4.6024 | 10.460 | DEET 6 |
| Amide N (`improper_Z_N_X_Y`) | 4.184 | 10.460 | DEET 1, icaridin 1 |
| Carbonyl C (`improper_O_C_X_Y`) | 43.932 | 43.932 | match |

Phase (180°) and multiplicity (2) agree everywhere.

### 5. Exact agreement

All bond and angle parameters (150 bonds, 279 angles across the three ITPs),
all 1-4 pair lists, and the interaction topology match exactly.
