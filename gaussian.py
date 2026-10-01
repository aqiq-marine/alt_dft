EV_TO_HARTREE = 1.0 / 27.2114

IRC_LINE = (
    "IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC"
)

HEADER = """\
-----------------------------------------------------------------------
# DirectMaxFlux/UMA
-----------------------------------------------------------------------

IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC

Copyright (c) 2025
Computational Biology Laboratory＠the University of Tokyo
"""


def format_structure(atoms) -> str:
    lines = [
        IRC_LINE,
        "                         Input orientation:",
        "---------------------------------------------------------------------",
        "Center     Atomic      Atomic             Coordinates (Angstroms)",
        "Number     Number       Type             X           Y           Z",
        "---------------------------------------------------------------------",
    ]

    for i, atom in enumerate(atoms):
        x, y, z = atom.position
        lines.append(
            f"{i + 1:5d}\t"
            f"{atom.number:<2d}\t"
            f"0\t"
            f"{x: .10f}\t"
            f"{y: .10f}\t"
            f"{z: .10f}"
        )

    lines.append(
        "---------------------------------------------------------------------"
    )

    return "\n".join(lines)


def write_gaussian_log(images, stream) -> None:
    """
    DMF trajectoryをGaussian log形式でstreamに書き込む。
    """

    if not images:
        raise ValueError("No frames found in trajectory.")

    last_energy = 0.0

    stream.write(HEADER + "\n\n")

    for i, atoms in enumerate(images):

        try:
            energy = (
                float(atoms.get_potential_energy())
                * EV_TO_HARTREE
            )
            last_energy = energy
        except Exception:
            energy = last_energy

        if i:
            point = i - 1

            stream.write(
                f"{IRC_LINE}\n"
                f"Pt {point} Step number   1 out of a maximum of  1\n"
            )

            stream.write(
                f"NET REACTION COORDINATE UP TO THIS POINT = "
                f"{float(point):20.10f}\n\n"
            )

        stream.write(
            format_structure(atoms) + "\n"
        )

        stream.write(
            f"SCF Done:  E(scf) =  {energy: .10f}     A.U.\n\n"
        )

    point = len(images) - 1

    stream.write(
        f"{IRC_LINE}\n"
        f"Pt {point} Step number   1 out of a maximum of  1\n"
    )

    stream.write(
        f"NET REACTION COORDINATE UP TO THIS POINT = "
        f"{float(point):20.10f}\n\n"
    )

    stream.write(
        f"{IRC_LINE}\n"
        "Normal termination of Gaussian\n"
    )
