from pathlib import Path
from typing import TypedDict

import modal


# ============================================================
# Types
# ============================================================

class OptimizationResult(TypedDict):
    filename: str
    data: bytes
    charge: int
    spin: int
    energy: float
    max_force: float


# ============================================================
# Modal
# ============================================================

app = modal.App("uma-optimization")


image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04",
        add_python="3.13",
    )
    .uv_pip_install(
        "fairchem-core",
        "ase",
        "numpy",
    )
)


# Modal Volume "models"
#
# models/
# └── models/
#     └── uma-s-1p2p1.pt
#
volume = modal.Volume.from_name(
    "models",
    create_if_missing=True,
)


# ============================================================
# Remote optimization
# ============================================================

@app.function(
    image=image,
    gpu="T4",
    volumes={
        "/models": volume,
    },
)
def optimize(
    files: list[tuple[str, bytes, int, int]],
) -> list[OptimizationResult]:

    import os

    os.environ["TORCHINDUCTOR_CACHE_DIR"] = "/models"

    from io import StringIO

    import numpy as np
    from ase.io import read, write
    from ase.optimize import BFGS
    from fairchem.core import FAIRChemCalculator
    from fairchem.core.units.mlip_unit import load_predict_unit


    FMAX = 0.05
    MAX_STEPS = 300

    # --------------------------------------------------------
    # Load UMA
    # --------------------------------------------------------

    predictor = load_predict_unit(
        "/models/models/uma-s-1p2p1.pt",
        device="cuda",
    )

    calculator = FAIRChemCalculator(
        predictor,
        task_name="omol",
    )

    results: list[OptimizationResult] = []

    # --------------------------------------------------------
    # Optimize each structure
    # --------------------------------------------------------

    for filename, input_data, charge, spin in files:

        print("=" * 60)
        print(f"Optimizing: {filename}")

        suffix = Path(filename).suffix.lower()

        if suffix not in {".xyz", ".pdb"}:
            print(f"Skipping unsupported file: {filename}")
            continue

        # ----------------------------------------------------
        # Load structure
        # ----------------------------------------------------

        text = input_data.decode("utf-8")

        supported_formats = {
            "xyz": "xyz",
            "pdb": "proteindatabank"
        }

        atoms = read(
            StringIO(text),
            format=supported_formats[suffix.lstrip(".")],
        )

        atoms.info["charge"] = charge
        atoms.info["spin"] = spin

        # ----------------------------------------------------
        # Set calculator
        # ----------------------------------------------------

        atoms.calc = calculator

        # ----------------------------------------------------
        # Optimization
        # ----------------------------------------------------

        optimizer = BFGS(atoms)

        optimizer.run(
            fmax=FMAX,
            steps=MAX_STEPS,
        )

        # ----------------------------------------------------
        # Results
        # ----------------------------------------------------

        energy = atoms.get_potential_energy()
        max_force = np.abs(atoms.get_forces()).max()

        # ----------------------------------------------------
        # Serialize optimized structure
        # ----------------------------------------------------

        output = StringIO()

        write(
            output,
            atoms,
            format="proteindatabank",
        )

        results.append(
            {
                "filename": filename,
                "data": output.getvalue().encode("utf-8"),
                "charge": charge,
                "spin": spin,
                "energy": float(energy),
                "max_force": float(max_force),
            }
        )

        print(f"Energy    : {energy:.8f} eV")
        print(f"Max force : {max_force:.8f} eV/Å")

    return results


# ============================================================
# Local entrypoint
# ============================================================

@app.local_entrypoint()
def main(*args: str):

    jobs = []

    i = 0

    while i < len(args):

        # filename
        filename = args[i]
        i += 1

        # デフォルト値
        charge = 0
        spin = 1

        # --charge
        if i < len(args) and args[i] == "--charge":
            i += 1

            if i >= len(args):
                raise ValueError(
                    f"{filename}: missing charge value"
                )

            try:
                charge = int(args[i])
            except ValueError:
                raise ValueError(
                    f"{filename}: charge must be an integer"
                )

            i += 1

        # --spin
        if i < len(args) and args[i] == "--spin":
            i += 1

            if i >= len(args):
                raise ValueError(
                    f"{filename}: missing spin value"
                )

            try:
                spin = int(args[i])
            except ValueError:
                raise ValueError(
                    f"{filename}: spin must be an integer"
                )

            i += 1

        jobs.append(
            {
                "filename": filename,
                "charge": charge,
                "spin": spin,
            }
        )

    if not jobs:
        raise ValueError(
            "No input files specified."
        )

    # 入力ファイルを読み込む
    files = []

    for job in jobs:
        raw_path = Path(job["filename"])
        input_path = raw_path if raw_path.exists() else Path("data") / raw_path

        if not input_path.exists():
            raise FileNotFoundError(
                f"File not found: {input_path}"
            )

        suffix = input_path.suffix.lower()

        if suffix not in {".xyz", ".pdb"}:
            raise ValueError(
                f"Unsupported file format: {input_path}"
            )

        files.append(
            (
                input_path.name,
                input_path.read_bytes(),
                job["charge"],
                job["spin"],
            )
        )

    results = optimize.remote(files)

    output_dir = Path("./data/optimized")
    output_dir.mkdir(parents=True, exist_ok=True)

    for result in results:

        input_path = Path(result["filename"])

        output_path = (
            output_dir
            / f"{input_path.stem}_opt.pdb"
        )

        output_path.write_bytes(
            result["data"]
        )

        print(
            f"{result['filename']}: "
            f"charge={result['charge']}, "
            f"spin={result['spin']}, "
            f"energy={result['energy']:.6f} eV, "
            f"max_force={result['max_force']:.6f} eV/Å"
        )

        print(f"  -> {output_path}")
