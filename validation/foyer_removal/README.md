# Foyer removal: equivalence evidence (2026-10-18)

Atom typing moved from Foyer 1.2.0 to the project's own engine
(`src/opls_ff_gen/opls_typer.py`), which runs the same rules
(`src/opls_ff_gen/data/opls_rules.json`, exported from Foyer's `oplsaa.xml`) plus
the atom-type archive.

- Branch: `foyer-removal`.
- Reference: tag `foyer-reference`, the last commit typing with Foyer.
- Checks: `src/scripts/check_foyer_equivalence.py`. It needs Foyer 1.2.0 and
  ParmEd, i.e. `pip install .[reference]` with the locally built wheels;
  after the removal it carries its own copy of the Foyer path.

## 1. Pattern graphs

`python src/scripts/check_foyer_equivalence.py patterns`

All 258 SMARTS definitions (Foyer's 228 plus the archive's) give the same
pattern graph as Foyer's `SMARTSGraph`: identical nodes and expression trees,
edges and bond orders, and rule-filter element. **258/258 identical.**

## 2. Molecules

`python src/scripts/check_foyer_equivalence.py molecules ...`

For every molecule, three things must be identical:
1. every atom's whitelist and blacklist after Foyer's iterative typing;
2. the final result of `assignment._assign_atom_types` with OPLS/2020
   parameters: types, classes, LJ and type rules, or the identical error;
3. the same with the legacy GROMACS parameters.

| Set | Molecules | Atoms | Identical (1, 2, 3) |
|---|---:|---:|---|
| RDKit-bundled: NCI first_5K (5000), ChEMBL2321810 series (1017), PubChem 200, BZR (163), CDK2 (47), EGFR (365), NCI first_200 (200) | 6,983 of 6,991 (8 unreadable by RDKit) | 242,914 | **6,983 / 6,983** for each |
| Project: 100 test molecules, rare-type set (100), project SDFs (17) | 217 | 4,583 | **217 / 217** for each |

Result files:
- `equivalence_large.json` (RDKit sets, run 10:51-11:13 on 9 cores);
- `equivalence_project_sets.json` (project sets, rerun after Foyer was
  removed from the package).

## 3. Package without Foyer

- The test suite (64 tests) passes with `foyer`, `parmed`, `gmso` and `lark`
  blocked from import.
- `pyproject.toml` no longer depends on Foyer or ParmEd.

## 4. End-to-end: the 100-molecule xTB workflow

The 100-molecule xTB workflow (typing, AMSOL charges, OPLS/2020 lookup,
on-the-fly xTB fits, ITP writing) was rerun with the Foyer-free code into
`benchmarks/opls2020_xtb_foyerfree`. It was compared file by file with the
Foyer-made batch `benchmarks/opls2020_xtb_2026_10_18`:

```bash
python src/scripts/compare_batches.py benchmarks/opls2020_xtb_2026_10_18 \
    benchmarks/opls2020_xtb_foyerfree --allow acetazolamide \
    --ignore-line "^; OPLS/2020 atom types, classes and Lennard-Jones parameters \((Foyer )?SMARTS typing;"
```

- **99/100 molecules identical:** ITP, GRO and TOP byte for byte, apart from
  the ignored provenance comment in the ITP header ("Foyer SMARTS typing" ->
  "SMARTS typing"). The topology and assignment JSON and the QM-fit reports
  (fitted values, RMSEs; timings excluded) are also identical.
- **Acetazolamide:** identical apart from the class label of its two
  thiadiazole N, `nd` -> `nb` (same LJ). This comes from the earlier GAFF2
  lookup commit, not from the Foyer removal.
- All 100 topologies pass `gmx grompp` (0 warnings) and a zero-step `mdrun`.
