from pathlib import Path

import modal


app = modal.App("dmf-ts")


image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04",
        add_python="3.13",
    )
    .uv_pip_install(
        "cyipopt-wheels",
        "pydmf",
        "fairchem-core",
        "ase",
        "numpy",
    )
    .add_local_file(
        "gaussian.py",
        "/root/gaussian.py",
    )
)


volume = modal.Volume.from_name(
    "models",
    create_if_missing=True,
)


@app.function(
    image=image,
    gpu="T4",
    volumes={
        "/models": volume,
    },
    timeout=60 * 60,
)
def run_ts(
    reactions: list[tuple[str, str, bytes, bytes, int, int]],
):
    import os

    os.environ["TORCHINDUCTOR_CACHE_DIR"] = "/models"

    from io import StringIO

    import numpy as np
    from ase.io import read, write
    from fairchem.core import FAIRChemCalculator
    from fairchem.core.units.mlip_unit import load_predict_unit
    from dmf import DirectMaxFlux, interpolate_fbenm

    from gaussian import write_gaussian_log

    # -----------------------------
    # Settings
    # -----------------------------

    charge = 0
    spin = 1

    nmove = 5
    convergence = "middle"

    # -----------------------------
    # Load UMA
    # -----------------------------

    print("Loading UMA...")

    predictor = load_predict_unit(
        "/models/models/uma-s-1p2p1.pt",
        device="cuda",
    )

    print("UMA loaded.")

    results = []

    # -----------------------------
    # Run reactions
    # -----------------------------

    for reaction_index, (
        reactant_name,
        product_name,
        reactant_data,
        product_data,
        charge,
        spin
    ) in enumerate(reactions, start=1):

        print()
        print("=" * 70)
        print(
            f"Reaction {reaction_index}/{len(reactions)}"
        )
        print(f"Reactant: {reactant_name}")
        print(f"Product : {product_name}")
        print("=" * 70)

        # -----------------------------
        # Read structures
        # -----------------------------

        supported_formats = {
            "xyz": "xyz",
            "pdb": "proteindatabank"
        }

        reactant_suffix = (
            Path(reactant_name)
            .suffix
            .lower()
            .lstrip(".")
        )

        product_suffix = (
            Path(product_name)
            .suffix
            .lower()
            .lstrip(".")
        )

        if reactant_suffix not in {"xyz", "pdb"}:
            raise ValueError(
                f"Unsupported reactant format: "
                f"{reactant_name}"
            )

        if product_suffix not in {"xyz", "pdb"}:
            raise ValueError(
                f"Unsupported product format: "
                f"{product_name}"
            )

        reactant = read(
            StringIO(
                reactant_data.decode("utf-8")
            ),
            format=supported_formats[reactant_suffix],
        )

        product = read(
            StringIO(
                product_data.decode("utf-8")
            ),
            format=supported_formats[product_suffix],
        )

        ref_images = [
            reactant,
            product,
        ]

        print(
            f"Reactant: {len(reactant)} atoms"
        )
        print(
            f"Product : {len(product)} atoms"
        )

        # -----------------------------
        # Initial interpolation
        # -----------------------------

        print(
            "Generating initial DMF path..."
        )

        initial = interpolate_fbenm(
            ref_images,
            correlated=True,
        )

        print(
            f"Initial images: "
            f"{len(initial.images)}"
        )

        # -----------------------------
        # Create DMF
        # -----------------------------

        mxflx = DirectMaxFlux(
            ref_images,
            coefs=initial.coefs.copy(),
            nmove=nmove,
            update_teval=True,
        )

        print(
            f"DMF images: {len(mxflx.images)}"
        )

        # -----------------------------
        # Set calculator
        # -----------------------------

        for image in mxflx.images:
            image.info["charge"] = charge
            image.info["spin"] = spin

            image.calc = FAIRChemCalculator(
                predictor,
                task_name="omol",
            )

        # -----------------------------
        # DMF optimization
        # -----------------------------

        print(
            "Starting DMF optimization..."
        )

        mxflx.solve(
            tol=convergence,
        )

        print(
            "DMF optimization finished."
        )

        # -----------------------------
        # Energies
        # -----------------------------

        energies = np.array(
            [
                image.get_potential_energy()
                for image in mxflx.images
            ]
        )

        ts_index = int(
            np.argmax(energies)
        )

        ts_image = mxflx.images[ts_index]

        relative_energy = (
            energies[ts_index]
            - energies[0]
        )

        print(
            f"TS image       : {ts_index}"
        )
        print(
            f"TS energy      : "
            f"{energies[ts_index]:.6f} eV"
        )
        print(
            f"Relative energy: "
            f"{relative_energy:.6f} eV"
        )

        # -----------------------------
        # Serialize TS structure
        # -----------------------------

        output = StringIO()

        # DMFの全imageをGaussian log形式に変換
        write_gaussian_log(
            mxflx.images,
            output,
        )

        results.append(
            {
                "reactant": reactant_name,
                "product": product_name,
                "data": (
                    output
                    .getvalue()
                    .encode("utf-8")
                ),
                "ts_index": ts_index,
                "ts_energy": float(
                    energies[ts_index]
                ),
                "relative_energy": float(
                    relative_energy
                ),
                "energies": energies.tolist(),
                "charge": charge,
                "spin": spin,
            }
        )

    return results

@app.local_entrypoint()
def main(*args: str):

    if not args:
        raise ValueError(
            "No reaction files specified."
        )

    reactions = []

    i = 0

    while i < len(args):

        # -----------------------------
        # Reactant
        # -----------------------------

        reactant_arg = Path(args[i])
        reactant = reactant_arg if reactant_arg.exists() else Path("data") / reactant_arg
        i += 1

        if i >= len(args):
            raise ValueError(
                f"{reactant_arg}: product file is missing."
            )

        # -----------------------------
        # Product
        # -----------------------------

        product_arg = Path(args[i])
        product = product_arg if product_arg.exists() else Path("data") / product_arg
        i += 1

        # -----------------------------
        # Default settings
        # -----------------------------

        charge = 0
        spin = 1

        # -----------------------------
        # Optional --charge
        # -----------------------------

        if i < len(args) and args[i] == "--charge":
            i += 1

            if i >= len(args):
                raise ValueError(
                    f"{reactant} {product}: "
                    "missing charge value."
                )

            try:
                charge = int(args[i])
            except ValueError:
                raise ValueError(
                    f"{reactant} {product}: "
                    "charge must be an integer."
                )

            i += 1

        # -----------------------------
        # Optional --spin
        # -----------------------------

        if i < len(args) and args[i] == "--spin":
            i += 1

            if i >= len(args):
                raise ValueError(
                    f"{reactant} {product}: "
                    "missing spin value."
                )

            try:
                spin = int(args[i])
            except ValueError:
                raise ValueError(
                    f"{reactant} {product}: "
                    "spin must be an integer."
                )

            i += 1

        # -----------------------------
        # Check files
        # -----------------------------

        if not reactant.exists():
            raise FileNotFoundError(reactant)

        if not product.exists():
            raise FileNotFoundError(product)

        if reactant.suffix.lower() not in {
            ".xyz",
            ".pdb",
        }:
            raise ValueError(
                f"Unsupported reactant format: "
                f"{reactant}"
            )

        if product.suffix.lower() not in {
            ".xyz",
            ".pdb",
        }:
            raise ValueError(
                f"Unsupported product format: "
                f"{product}"
            )

        # -----------------------------
        # Add reaction
        # -----------------------------

        reactions.append(
            (
                reactant.name,
                product.name,
                reactant.read_bytes(),
                product.read_bytes(),
                charge,
                spin,
            )
        )

    # -----------------------------
    # Run Modal
    # -----------------------------

    results = run_ts.remote(reactions)

    # -----------------------------
    # Save results
    # -----------------------------

    output_dir = Path("data/TS_candidates")
    output_dir.mkdir(parents=True, exist_ok=True)

    print()
    print("=" * 70)
    print("TS Summary")
    print("=" * 70)

    for result in results:

        reactant_path = Path(
            result["reactant"]
        )

        output_path = (
            output_dir
            / f"{reactant_path.stem}_TS.log"
        )

        output_path.write_bytes(
            result["data"]
        )

        print()
        print(
            f"Reactant: {result['reactant']}"
        )
        print(
            f"Product : {result['product']}"
        )
        print(
            f"Charge  : {result['charge']}"
        )
        print(
            f"Spin    : {result['spin']}"
        )
        print(
            f"TS image       : {result['ts_index']}"
        )
        print(
            f"TS energy      : "
            f"{result['ts_energy']:.6f} eV"
        )
        print(
            f"Relative energy: "
            f"{result['relative_energy']:.6f} eV"
        )
        print(
            f"Output         : {output_path}"
        )
