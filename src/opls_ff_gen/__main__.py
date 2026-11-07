"""Command-line interface for complete OPLS-AA GROMACS topology generation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rdkit import Chem

from .assignment import DEFAULT_ATOMTYPE_ARCHIVE
from .mm import DEFAULT_ENGINE, ENGINES
from .qm import DEFAULT_DFT, DFTLevel
from .structure import prepare_structure
from .workflow import generate_gromacs_topology


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        epilog="Other commands: 'opls-ff-gen doctor' (status of every dependency), "
        "'opls-ff-gen config' (paths of programs you already have), "
        "'opls-ff-gen install-si' (install the OPLS/2020 SI files).",
        description=(
            "Generate a GROMACS-ready OPLS/2020 topology: OPLS/2020 atom types and "
            "parameters (SMARTS typing, GAFF2 LJ for atoms OPLS cannot type), "
            "1.14*CM1A-LBCC charges (MOPAC AM1 + CM1A), and on-the-fly QM fitting of bonded "
            "terms missing from the OPLS/2020 tables."
        )
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--structure", type=Path,
                        help="3D structure file: .sdf/.mol (one record), .mol2 or .pdb")
    source.add_argument("--sdf", type=Path, help="one-record 3D SDF (same as --structure)")
    source.add_argument("--smiles", help="SMILES; a 3D conformer is generated and optimized")
    parser.add_argument(
        "--optimize",
        choices=("none", "mmff", "xtb"),
        help="geometry optimization before the charges and fits: none = use the input "
        "coordinates (default for structure files), xtb = GFN2-xTB (default for SMILES), "
        "mmff = RDKit MMFF94. SMILES are always optimized.",
    )
    parser.add_argument("--template-smiles",
                        help="PDB or MOL2 input: assign bond orders and formal charges "
                        "from this SMILES")
    parser.add_argument("--net-charge", type=int,
                        help="net charge: used for PDB bond-order perception and checked "
                        "against the formal charges of any structure file")
    parser.add_argument("--name", help="Molecule/output name")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--charge-model",
        choices=("auto", "cm1a", "cm1a-lbcc"),
        default="auto",
        help="auto uses LBCC for neutrals and unscaled CM1A for ions",
    )
    parser.add_argument(
        "--qm",
        choices=("xtb", "dft", "none"),
        default="none",
        help="optional run-time QM fit of every bonded term missing from the OPLS/2020 "
        "table (no other lookup is then used): none = no QM (default): missing terms take "
        "the other tables and the fallbacks (bonds/angles substitution or empirical, "
        "torsions substitution, else zero); xtb = GFN2-xTB scans and Hessian; "
        "dft = DFT geometry and Hessian and DFT single points on the xTB scans (level: "
        "--dft-functional/--dft-basis/--dft-dispersion)",
    )
    parser.add_argument(
        "--torsion-substitution",
        choices=("parmchk2", "none"),
        default="parmchk2",
        help="without --qm, torsions missing from OPLS/2020: the GAFF2-derived general "
        "X-B-C-X term, substitution as AmberTools parmchk2 does (PARMCHK.DAT scores) and "
        "the torsion fallback, else zero (default); or none (strict table lookup)",
    )
    parser.add_argument(
        "--qm-bonded-table",
        choices=("use", "none"),
        default="use",
        help="without --qm, bonds and angles missing from OPLS/2020: the value recorded "
        "from earlier QM fits (src/opls_ff_gen/data/qm_bonded_terms.par; default) before the "
        "fallback; or none",
    )
    parser.add_argument("--qm-threads", type=int, default=8)
    parser.add_argument("--mm-engine", choices=ENGINES, default=DEFAULT_ENGINE,
                        help="MM engine for the restrained relaxations of the torsion fit: "
                        f"openmm or gromacs (default {DEFAULT_ENGINE}; env OPLS_FF_GEN_MM_ENGINE)")
    parser.add_argument("--dft-functional", default=DEFAULT_DFT.functional,
                        help="PySCF/libxc functional for --qm dft (default b3lyp; e.g. pbe0, "
                        "m06-2x, wb97x-v)")
    parser.add_argument("--dft-basis", default=DEFAULT_DFT.basis,
                        help="basis set for --qm dft (default def2-svp; e.g. def2-tzvp)")
    parser.add_argument("--dft-dispersion", default=DEFAULT_DFT.dispersion,
                        help="dispersion correction for --qm dft: d3bj (default), d3zero, d4 "
                        "or none (use none for functionals with built-in dispersion)")
    parser.add_argument(
        "--no-amber-fallback",
        action="store_true",
        help="fail on atoms OPLS cannot type instead of using GAFF2 types and LJ",
    )
    parser.add_argument(
        "--atomtype-archive",
        type=Path,
        default=DEFAULT_ATOMTYPE_ARCHIVE,
        help="atom-type archive: extra SMARTS, precedence, OPLS/2020 refinements and "
        "borrowed OPLS-AA types (N-oxide, hydrazide, triazole, xanthine, ...)",
    )
    parser.add_argument(
        "--no-atomtype-archive",
        action="store_true",
        help="type with the OPLS-AA rule set only (data/opls_rules.json)",
    )
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument(
        "--box-padding-nm",
        type=float,
        default=1.0,
        help="padding around the molecule in the generated GRO box",
    )
    return parser


CONFIG_HELP = """\
opls-ff-gen config: where opls-ff-gen finds programs you already have.

  opls-ff-gen config show                 list the settings and the config file
  opls-ff-gen config set KEY PATH         use this program (checked to be executable)
  opls-ff-gen config unset KEY            forget it (fall back to PATH)

KEY: mopac, xtb, gmx (programs), xtbpath (xtb parameter directory).
Environment variables (OPLS_FF_GEN_MOPAC, OPLS_FF_GEN_XTB, OPLS_FF_GEN_GMX,
OPLS_FF_GEN_XTBPATH) take precedence over the file. Check with: opls-ff-gen doctor
"""


def _config_command(argv: list[str]) -> int:
    from . import tools

    if not argv or argv[0] in ("-h", "--help", "help"):
        print(CONFIG_HELP, end="")
        return 0
    action, rest = argv[0], argv[1:]
    try:
        if action == "show" and not rest:
            values = tools.read_config()
            print(f"config file: {tools.config_path()}"
                  + ("" if tools.config_path().is_file() else " (none)"))
            for key in [*tools.PROGRAMS, *tools.SETTINGS]:
                print(f"  {key:<8} {values.get(key, '(not set)')}")
            return 0
        if action == "set" and len(rest) == 2:
            path = tools.set_config(rest[0], rest[1])
            print(f"{rest[0]} = {tools.read_config()[rest[0]]}  ({path})")
            return 0
        if action == "unset" and len(rest) == 1:
            path = tools.set_config(rest[0], None)
            print(f"{rest[0]} removed ({path})")
            return 0
    except (ValueError, RuntimeError) as exc:
        print(f"opls-ff-gen config: {exc}", file=sys.stderr)
        return 2
    print(CONFIG_HELP, end="", file=sys.stderr)
    return 2


def _doctor_command(argv: list[str]) -> int:
    from . import tools

    if argv:
        print("usage: opls-ff-gen doctor  (status of every dependency)", file=sys.stderr)
        return 2
    report, required_ok = tools.doctor()
    print(report, end="")
    return 0 if required_ok else 1


def _install_si_command(argv: list[str]) -> int:
    from . import si

    parser = argparse.ArgumentParser(
        prog="opls-ff-gen install-si",
        description="Install the OPLS/2020 Supporting Information files "
        f"(doi:{si.DOI}; not redistributable, so download them yourself from the "
        "article page). Both files are checked against the release opls-ff-gen was "
        "validated with and copied to a per-user directory that survives upgrades "
        f"(default {si.user_si_dir()}).")
    parser.add_argument("files", nargs="+", type=Path,
                        help="jp3c06602_si_002.txt and jp3c06602_si_003.txt (any order)")
    parser.add_argument("--dest", type=Path,
                        help="install here instead (then set OPLS_FF_GEN_SI_DIR to it)")
    parser.add_argument("--force", action="store_true",
                        help="install files that differ from the validated release")
    args = parser.parse_args(argv)
    try:
        report = si.install(args.files, dest=args.dest, force=args.force)
    except si.SIInstallError as exc:
        print(f"opls-ff-gen install-si: {exc}", file=sys.stderr)
        return 2
    print("\n".join(report))
    if args.dest is not None:
        print(f"Set {si.SI_DIR_ENV}={args.dest.expanduser().resolve()} so that opls-ff-gen "
              "reads them from there.")
    directory, source = si.si_dir()
    print(f"opls-ff-gen reads the SI from: {directory} [{source}]")
    return 0


COMMANDS = {"config": _config_command, "doctor": _doctor_command,
            "install-si": _install_si_command}


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] in COMMANDS:
        raise SystemExit(COMMANDS[argv[0]](argv[1:]))
    parser = build_parser()
    args = parser.parse_args(argv)
    from .tools import preflight

    errors, warnings = preflight(
        qm=args.qm, mm_engine=args.mm_engine,
        optimize_xtb=args.optimize == "xtb" or (args.optimize is None and args.smiles is not None),
    )
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    if not errors:
        from .tools import run_version_notes

        xtb_used = args.qm != "none" or args.optimize == "xtb" or (
            args.optimize is None and args.smiles is not None)
        for note in run_version_notes(xtb_used=xtb_used):
            print(f"NOTE: {note}", file=sys.stderr)
    if errors:
        print("opls-ff-gen: cannot start:\n  "
              + "\n  ".join(errors) + "\nStatus of every dependency: opls-ff-gen doctor",
              file=sys.stderr)
        raise SystemExit(2)
    structure_path = args.structure or args.sdf
    name = args.name or (structure_path.stem if structure_path else "molecule")
    molecule = prepare_structure(
        smiles=args.smiles,
        path=structure_path,
        name=name,
        optimize=args.optimize,
        template_smiles=args.template_smiles,
        net_charge=args.net_charge,
        workdir=args.output_dir / "geometry",
        threads=args.qm_threads,
    )
    result = generate_gromacs_topology(
        molecule,
        args.output_dir,
        name=name,
        charge_model=args.charge_model,
        atomtype_archive=None if args.no_atomtype_archive else args.atomtype_archive,
        amber_fallback=not args.no_amber_fallback,
        qm_method=None if args.qm == "none" else args.qm,
        qm_threads=args.qm_threads,
        mm_engine=args.mm_engine,
        torsion_substitution=args.torsion_substitution == "parmchk2",
        qm_bonded_table=args.qm_bonded_table == "use",
        dft_level=DFTLevel(args.dft_functional, args.dft_basis,
                           None if str(args.dft_dispersion).lower() == "none"
                           else args.dft_dispersion),
        timeout_seconds=args.timeout,
        box_padding_nm=args.box_padding_nm,
    )
    print(
        f"{result.assignment.name}: {result.assignment.atom_count} atoms; "
        f"{result.assignment.charge_model}; "
        f"final charge={result.assignment.final_charge_sum_e:+.10f} e"
    )
    rule_counts = ", ".join(
        f"{count} {rule}" for rule, count in result.bonded.rule_counts().items()
    )
    print(f"bonded terms: {rule_counts}")
    for warning in result.bonded.warnings:
        print(f"WARNING: {warning}")
    print(f"GROMACS ITP:     {result.itp_path}")
    print(f"GROMACS GRO:     {result.gro_path}")
    print(f"GROMACS TOP:     {result.top_path}")
    print(f"topology JSON:   {result.topology_json_path}")
    print(f"assignment JSON: {result.assignment.json_path}")
    print(f"assigned SDF:    {result.assignment.sdf_path}")
    print(f"structure read:  {Chem.MolToSmiles(Chem.RemoveHs(molecule))} "
          f"(net charge {Chem.GetFormalCharge(molecule):+d})")
    print(f"geometry:        {molecule.GetProp('Geometry_Optimization')}")
    if molecule.HasProp("Geometry_Fallback"):
        print(f"WARNING: geometry fallback: {molecule.GetProp('Geometry_Fallback')}")


if __name__ == "__main__":
    main()
