"""Chemical-context keys for bonded terms.

A key is the interaction's bonded classes plus a SMARTS chemical context
whose atoms mapped :1-:n are the interaction atoms, so one class tuple can
carry different parameters in different chemistry (for example a ring N that
is a triazole N1 in one molecule and an imidazole N1 in another). The on-the-fly
torsion fit (``qm_fit``) shares one parameter set among the missing quartets
of a molecule with the same key.
"""

from __future__ import annotations

from rdkit import Chem


def _atom_primitive(atom: Chem.Atom) -> str:
    # No formal charge: resonance-equivalent atoms (the two O of a nitro
    # group, written O= and O-) must share one context.
    aromatic = "a" if atom.GetIsAromatic() else "A"
    ring = "R" if atom.IsInRing() else "!R"
    return (f"#{atom.GetAtomicNum()};{aromatic};X{atom.GetDegree()};"
            f"H{atom.GetTotalNumHs(includeNeighbors=True)};{ring}")


def normalized(molecule: Chem.Mol) -> Chem.Mol:
    """Copy with aromaticity perceived from a Kekule structure.

    RDKit's aromatic flags depend on the input route (caffeine's
    pyrimidinedione ring is non-aromatic from its SMILES but aromatic when
    read back from an SDF), so contexts are built and matched on this copy.
    Atom indices are unchanged.
    """
    copy = Chem.RWMol(molecule)
    Chem.Kekulize(copy, clearAromaticFlags=True)
    Chem.SanitizeMol(copy)
    return copy.GetMol()


def context_smarts(molecule: Chem.Mol, atoms: tuple[int, ...]) -> str:
    """SMARTS context for a bond, angle or torsion (0-based atoms, in order).

    Every atom is described by its own properties; the central atoms (both
    atoms of a bond, the apex of an angle, the two middle atoms of a torsion)
    additionally list every neighbor (element and aromaticity), which keeps
    e.g. a triazole N1 apart from an imidazole N1.
    """
    molecule = normalized(molecule)
    central = range(1, len(atoms) + 1) if len(atoms) == 2 else range(2, len(atoms))
    parts = []
    for position, index in enumerate(atoms, start=1):
        atom = molecule.GetAtomWithIdx(index)
        primitive = _atom_primitive(atom)
        if position in central:
            neighbors = sorted(
                f"[#{n.GetAtomicNum()};{'a' if n.GetIsAromatic() else 'A'}]"
                for n in atom.GetNeighbors()
            )
            primitive += ";$([#%d]%s)" % (
                atom.GetAtomicNum(), "".join(f"(~{n})" for n in neighbors)
            )
        parts.append(f"[{primitive}:{position}]")
    return "~".join(parts)


def context_key(
    molecule: Chem.Mol, atoms: tuple[int, ...], classes: tuple[str, ...]
) -> tuple[tuple[str, ...], str]:
    """Key (classes, context SMARTS) of an interaction, in canonical orientation.

    An interaction read from either end must give the same key, and the
    classes and the SMARTS map numbers must run in the same direction.
    """
    forward = (tuple(classes), context_smarts(molecule, tuple(atoms)))
    reverse = (tuple(classes)[::-1], context_smarts(molecule, tuple(atoms)[::-1]))
    return min(forward, reverse)
