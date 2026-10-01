import argparse
from pathlib import Path

from atom_mapping import make_mapping, validate_mapping
from molecule_io import (
    add_explicit_hydrogens,
    load_molecule,
    optimize_3d_fragments,
    write_pdb,
)


def main():
    parser = argparse.ArgumentParser(description="Map atoms between MOL files and write reordered PDB files.")
    parser.add_argument("reactant", type=Path)
    parser.add_argument("product", type=Path)
    parser.add_argument("--reactant-out", type=Path, default=None)
    parser.add_argument("--product-out", type=Path, default=None)
    args = parser.parse_args()

    reactant_path = args.reactant if args.reactant.exists() else Path("data") / args.reactant
    product_path = args.product if args.product.exists() else Path("data") / args.product

    if not reactant_path.exists():
        raise FileNotFoundError(f"Reactant file not found: {args.reactant}")
    if not product_path.exists():
        raise FileNotFoundError(f"Product file not found: {args.product}")

    if args.reactant_out is None:
        args.reactant_out = Path("data") / f"{reactant_path.stem}_mapped.pdb"
    if args.product_out is None:
        args.product_out = Path("data") / f"{product_path.stem}_mapped.pdb"

    print("Loading reactant...")
    reactant_mol = load_molecule(reactant_path)
    reactant_mol = add_explicit_hydrogens(reactant_mol)
    reactant_mol, reactant_optimization = optimize_3d_fragments(reactant_mol)
    print(f"Reactant atoms: {reactant_mol.GetNumAtoms()}\n")
    print(f"Reactant fragment optimization: {reactant_optimization}")
    print("Loading product...")
    product_mol = load_molecule(product_path)
    product_mol = add_explicit_hydrogens(product_mol)
    product_mol, product_optimization = optimize_3d_fragments(product_mol)
    print(f"Product atoms: {product_mol.GetNumAtoms()}")
    print(f"Product fragment optimization: {product_optimization}")

    mapping, confidence, _ = make_mapping(reactant_mol, product_mol)
    print(f"Mapped atoms: {len(mapping)}")
    validate_mapping(reactant_mol, product_mol, mapping)

    reactant_order = list(range(reactant_mol.GetNumAtoms()))
    product_order = []
    unmatched_reactant = []
    unmatched_product = set(range(product_mol.GetNumAtoms()))
    for r_idx in reactant_order:
        if r_idx in mapping:
            product_order.append(mapping[r_idx])
            unmatched_product.discard(mapping[r_idx])
        else:
            unmatched_reactant.append(r_idx)
    product_order.extend(sorted(unmatched_product))

    print(f"Reactant order: {reactant_order}")
    print(f"Product order: {product_order}")
    print(f"Unmatched reactant atoms: {unmatched_reactant}")
    print(f"Unmatched product atoms: {sorted(unmatched_product)}")
    if reactant_mol.GetNumAtoms() != product_mol.GetNumAtoms():
        print("\nWARNING:\nReactant and product have different atom counts.")
        print("This must be resolved before using these PDB files as a PyDMF path.")

    args.reactant_out.parent.mkdir(parents=True, exist_ok=True)
    args.product_out.parent.mkdir(parents=True, exist_ok=True)
    write_pdb(reactant_mol, reactant_order, args.reactant_out)
    write_pdb(product_mol, product_order, args.product_out)
    print("\nDone.")
    print(f"Reactant PDB: {args.reactant_out}")
    print(f"Product PDB:  {args.product_out}")


if __name__ == "__main__":
    main()
