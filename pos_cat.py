#!/usr/bin/env python3

from pathlib import Path
import argparse

import numpy as np
from ase import Atoms
from ase.io import read, write


def center_of_mass(atoms: Atoms):
    """Geometric center of atoms."""
    return atoms.get_positions().mean(axis=0)


def translate_catalyst_away(
    substrate: Atoms,
    catalyst: Atoms,
    min_distance: float = 8.0,
):
    """
    Move catalyst away from substrate.

    The catalyst is translated without changing its internal geometry.

    min_distance:
        Minimum allowed inter-component atom distance [Å].
    """

    substrate_pos = substrate.get_positions()
    catalyst_pos = catalyst.get_positions()

    substrate_center = substrate_pos.mean(axis=0)
    catalyst_center = catalyst_pos.mean(axis=0)

    # Initial direction: catalyst center -> away from substrate center
    direction = catalyst_center - substrate_center

    norm = np.linalg.norm(direction)

    if norm < 1e-8:
        # If the centers overlap, use +X direction
        direction = np.array([1.0, 0.0, 0.0])
    else:
        direction = direction / norm

    def min_interatomic_distance():
        diff = (
            substrate_pos[:, None, :]
            - catalyst.get_positions()[None, :, :]
        )

        distances = np.linalg.norm(diff, axis=-1)

        return distances.min()

    # Move catalyst until every substrate-catalyst atom pair
    # is sufficiently separated.
    step = 2.0

    while min_interatomic_distance() < min_distance:
        catalyst.translate(direction * step)

    return catalyst


def make_combined_structure(
    substrate: Atoms,
    catalyst: Atoms,
    min_distance: float = 8.0,
):
    """
    Combine substrate and catalyst while keeping them separated.
    """

    catalyst = catalyst.copy()

    catalyst = translate_catalyst_away(
        substrate,
        catalyst,
        min_distance=min_distance,
    )

    combined = substrate + catalyst

    return combined


def validate_reactant_product(
    reactant: Atoms,
    product: Atoms,
):
    """
    Check atom correspondence between reactant and product.
    """

    if len(reactant) != len(product):
        raise ValueError(
            "Reactant and product have different numbers of atoms: "
            f"{len(reactant)} vs {len(product)}"
        )

    reactant_symbols = reactant.get_chemical_symbols()
    product_symbols = product.get_chemical_symbols()

    for i, (r, p) in enumerate(
        zip(reactant_symbols, product_symbols)
    ):
        if r != p:
            raise ValueError(
                "Reactant/product atom correspondence is inconsistent "
                f"at atom {i + 1}: {r} != {p}"
            )


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "reactant",
        help="Reactant (.pdb or .xyz)",
    )

    parser.add_argument(
        "product",
        help="Product (.pdb or .xyz)",
    )

    parser.add_argument(
        "catalyst",
        help="Catalyst (.pdb or .xyz)",
    )

    parser.add_argument(
        "--min-distance",
        type=float,
        default=8.0,
        help="Minimum substrate-catalyst distance in Å (default: 8.0)",
    )

    parser.add_argument(
        "--outdir",
        default="data",
        help="Output directory",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Read
    # --------------------------------------------------------

    def resolve_path(arg_val: str) -> Path:
        p = Path(arg_val)
        if p.exists():
            return p
        dp = Path("data") / p
        return dp if dp.exists() else p

    reactant_path = resolve_path(args.reactant)
    product_path = resolve_path(args.product)
    catalyst_path = resolve_path(args.catalyst)

    if not reactant_path.exists():
        raise FileNotFoundError(f"Reactant file not found: {args.reactant}")
    if not product_path.exists():
        raise FileNotFoundError(f"Product file not found: {args.product}")
    if not catalyst_path.exists():
        raise FileNotFoundError(f"Catalyst file not found: {args.catalyst}")

    reactant = read(reactant_path)
    product = read(product_path)
    catalyst = read(catalyst_path)

    # --------------------------------------------------------
    # Validate atom correspondence
    # --------------------------------------------------------

    validate_reactant_product(
        reactant,
        product,
    )

    # --------------------------------------------------------
    # Generate structures
    # --------------------------------------------------------

    reactant_cat = make_combined_structure(
        reactant,
        catalyst,
        min_distance=args.min_distance,
    )

    product_cat = make_combined_structure(
        product,
        catalyst,
        min_distance=args.min_distance,
    )

    # --------------------------------------------------------
    # Write
    # --------------------------------------------------------

    
    reactant_out = outdir / f"{reactant_path.stem}_cat.pdb"
    product_out = outdir / f"{product_path.stem}_cat.pdb"

    write(
        reactant_out,
        reactant_cat,
        format="proteindatabank",
    )

    write(
        product_out,
        product_cat,
        format="proteindatabank",
    )

    print(f"Reactant : {len(reactant)} atoms")
    print(f"Product  : {len(product)} atoms")
    print(f"Catalyst : {len(catalyst)} atoms")
    print()
    print(f"Output:")
    print(f"  {reactant_out}")
    print(f"  {product_out}")
    print()
    print(
        f"Minimum substrate-catalyst distance: "
        f"{args.min_distance:.1f} Å"
    )


if __name__ == "__main__":
    main()
