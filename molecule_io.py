from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def load_molecule(path: Path) -> Chem.Mol:
    print(f"MOL: {path}")
    print(f"MOL size: {path.stat().st_size} bytes")
    mol = Chem.MolFromMolFile(str(path), removeHs=False)
    if mol is None:
        raise RuntimeError(f"Could not read molecule from {path}")
    print(f"Loaded molecule: {mol.GetNumAtoms()} atoms")
    print(f"SMILES: {Chem.MolToSmiles(mol, canonical=False, isomericSmiles=True)}")
    return mol


def add_explicit_hydrogens(mol: Chem.Mol) -> Chem.Mol:
    """Return a copy with implicit hydrogens made explicit.

    RXNMapper can only assign atom maps to atoms present in the reaction
    SMILES.  Adding H here therefore makes the hydrogen atoms part of the
    same identity-tracking and coordinate-reordering flow as heavy atoms.
    Existing explicit hydrogens are preserved by RDKit.
    """
    with_hydrogens = Chem.AddHs(Chem.Mol(mol), addCoords=True)
    Chem.SanitizeMol(with_hydrogens)
    return with_hydrogens


def _optimize_fragment(
    fragment: Chem.Mol, random_seed: int
) -> tuple[Chem.Mol, str, int]:
    """Generate and optimize one connected fragment."""
    optimized = Chem.Mol(fragment)
    optimized.RemoveAllConformers()

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.clearConfs = True
    if AllChem.EmbedMolecule(optimized, params) != 0:
        raise RuntimeError("ETKDG failed to generate a 3D conformer for a fragment")

    mmff = AllChem.MMFFGetMoleculeProperties(optimized, mmffVariant="MMFF94s")
    if mmff is not None:
        status = AllChem.MMFFOptimizeMolecule(
            optimized, mmffVariant="MMFF94s", maxIters=500
        )
        return optimized, "MMFF94s", status

    if not AllChem.UFFHasAllMoleculeParams(optimized):
        raise RuntimeError("Neither MMFF94s nor UFF has parameters for a fragment")
    status = AllChem.UFFOptimizeMolecule(optimized, maxIters=500)
    return optimized, "UFF", status


def optimize_3d_fragments(
    mol: Chem.Mol, random_seed: int = 20260911, fragment_gap: float = 5.0
) -> tuple[Chem.Mol, list[tuple[str, int]]]:
    """Optimize each disconnected fragment independently.

    MMFF94s is preferred for the product.  Some elements/charge states are
    outside MMFF's parameter set, so UFF is used as a deterministic fallback.
    The optimized fragments are recombined in their original order and atom
    order, which keeps RXNMapper identity tracking valid.
    """
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    optimized_fragments = []
    methods = []
    for fragment_index, fragment in enumerate(fragments):
        optimized_fragment, method, status = _optimize_fragment(
            fragment, random_seed + fragment_index
        )
        optimized_fragments.append(optimized_fragment)
        methods.append((method, status))

    optimized = optimized_fragments[0]
    for fragment in optimized_fragments[1:]:
        optimized = Chem.CombineMols(optimized, fragment)

    # Each independently embedded fragment starts near the origin.  Place
    # them side by side so that atoms from different fragments cannot appear
    # at a bonding distance in the output XYZ.  The gap is measured between
    # the fragments' X bounding boxes, so it is independent of their shapes.
    if len(optimized_fragments) > 1:
        conf = optimized.GetConformer()
        next_min_x = None
        atom_offset = 0
        for fragment in optimized_fragments:
            atom_indices = range(atom_offset, atom_offset + fragment.GetNumAtoms())
            min_x = min(conf.GetAtomPosition(i).x for i in atom_indices)
            max_x = max(conf.GetAtomPosition(i).x for i in atom_indices)
            shift_x = 0.0 if next_min_x is None else next_min_x - min_x
            for atom_index in atom_indices:
                position = conf.GetAtomPosition(atom_index)
                conf.SetAtomPosition(
                    atom_index,
                    (position.x + shift_x, position.y, position.z),
                )
            next_min_x = max_x + shift_x + fragment_gap
            atom_offset += fragment.GetNumAtoms()
    return optimized, methods


def write_xyz(mol: Chem.Mol, atom_order, output_path: Path, comment=""):
    if not mol.GetNumConformers():
        raise RuntimeError(f"Molecule has no coordinates: {mol}")
    conf = mol.GetConformer()
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"{len(atom_order)}\n{comment}\n")
        for idx in atom_order:
            atom = mol.GetAtomWithIdx(idx)
            pos = conf.GetAtomPosition(idx)
            f.write(f"{atom.GetSymbol():<3s} {pos.x: .10f} {pos.y: .10f} {pos.z: .10f}\n")


def write_pdb(mol: Chem.Mol, atom_order, output_path: Path):
    """Write a PDB with atoms in the requested order and preserve bonds."""
    if not mol.GetNumConformers():
        raise RuntimeError(f"Molecule has no coordinates: {mol}")
    reordered = Chem.RenumberAtoms(mol, list(atom_order))
    Chem.MolToPDBFile(reordered, str(output_path), flavor=0)
