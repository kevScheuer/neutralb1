"""Project vector-pseudoscalar PWA fit results to moments

This script uses the conversion from partial wave complex values to moments to "project"
the potentially ambiguous PWA fit results into unique moments. It takes as input the
best fit parameters, obtained using the "-s" flag in the fit command. The moments are
computed and saved to an output csv file. Multiple files can be given as input, with
optional sorting based on the file name or path.

NOTE: This script assumes that the amplitudes are written in the vec-ps eJPmL format.
For example, the positive reflectivity, JP=1+, m=0, S-wave amplitude would be written
in the cfg file as [reaction]::RealNegSign::p1p0S.

TODO: handle free floating parameters like the D/S ratio
"""

import argparse
import itertools
import os
from dataclasses import dataclass
from typing import List, Set, TextIO

import numba
import numpy as np
import pandas as pd
import pwa_tools
import utils
from sympy.physics.quantum.cg import CG


@dataclass(frozen=True)
class Wave:
    name: str
    reflectivity: int
    spin: int
    parity: int
    m: int
    l: int
    real: float
    imaginary: float


def main(args: dict) -> None:

    # Error / default value handling
    for input_file in args["input"]:
        if not os.path.isfile(input_file):
            raise ValueError(f"Input file {input_file} does not exist")
    if args["output"] and not args["output"].endswith(".csv"):
        args["output"] += ".csv"
    else:
        args["output"] = "moments.csv"

    input_files = (
        utils.sort_input_files(args["input"]) if args["sorted"] else args["input"]
    )

    if args["preview"]:
        print("Files that will be processed:")
        for file in input_files:
            print(f"\t{file}")
        return

    # create a dataframe to store the moments, with file names as the index
    df = pd.DataFrame(index=input_files)

    # Loop to calculate moments for each input file
    for file in input_files:
        waves = get_waves(file)  # obtain waves from the input file

        # Prepare quantum numbers of moments
        Jv_array = np.array([0, 2])  # CG coefficient always makes Jv=1 be 0
        # (-) Lambda values are directly proportional to (+) ones, no need to calculate
        Lambda_array = np.arange(0, 3)

        min_J = min(wave.spin for wave in waves)
        max_J = max(wave.spin for wave in waves)
        J_array = np.arange(min_J, max_J + 1)

        min_m = max(min(wave.m for wave in waves), 0)  # don't need min below 0
        max_m = max(wave.m for wave in waves)
        M_array = np.arange(min_m, max_m + 1)  # like Lambda, -m ∝ +m moments

        # calculate each moment and add it to the dataframe
        for alpha in range(3):
            for Jv, Lambda, J, M in itertools.product(
                Jv_array, Lambda_array, J_array, M_array
            ):
                moment_str = f"H{alpha}({Jv},{Lambda},{J},{M})"
                moment_val = calculate_moment(alpha, Jv, Lambda, J, M, waves)
                df.at[file, moment_str] = moment_val

    pass


def get_waves(file: TextIO) -> Set[Wave]:
    """Obtain the set of waves, with their real and imaginary parts, from a fit result

    Args:
        file (TextIO): best fit parameters file, obtained using the "-s" flag on the fit
            command, that contains the real and imaginary parts of the best fit result
    Returns:
        Set[Wave]: all the waves and their information found in the file
    """
    waves = set()
    with open(file, "r") as f:
        for line in f:
            if (
                line.startswith("#")
                or not line.strip()
                or "isotropic" in line
                or "Bkgd" in line
            ):
                continue
            parts = line.split()

            if "::" not in parts[1]:
                raise ValueError(f"Unexpected line format in file {file}")

            # parse amplitude into its quantum number values
            amplitude = parts[1].split("::")[-1]
            reflectivity = 1 if amplitude[0] == "p" else -1
            spin = int(amplitude[1])
            parity = 1 if amplitude[2] == "p" else -1
            m = pwa_tools.char_to_int(amplitude[3])
            l = pwa_tools.char_to_int(amplitude[4])

            # obtain real and imaginary parts
            if parts[2] == "cartesian":
                re = float(parts[3])
                im = float(parts[4])
            elif parts[2] == "polar":
                magnitude = float(parts[3])
                phase = float(parts[4])
                re = magnitude * np.cos(phase)
                im = magnitude * np.sin(phase)
            else:
                raise ValueError(f"Unexpected complex number format in file {file}")

            waves.add(
                Wave(
                    name=amplitude,
                    reflectivity=reflectivity,
                    spin=spin,
                    parity=parity,
                    m=m,
                    l=l,
                    real=re,
                    imaginary=im,
                )
            )

    return waves


@numba.njit(cache=True)
def calculate_moment(
    alpha: int, Jv: int, Lambda: int, J: int, M: int, waves: List[Wave]
) -> float:
    """Calculate the moment for a given set of quantum numbers

    Vector-pseudoscalar moments are indexed by 4 quantum numbers that arise from
    combining the Wigner D functions of the resonance and vector decays.

    Args:
        alpha (int): indexes which term of the intensity the moment is associated with
            0: unpolarized
            1: polarized cos(2*Phi)
            2: polarized sin(2*Phi)
        Jv (int): Total spin resulting from the two Wigner D functions capturing the
            vector decay
        Lambda (int): Total helicity resulting from the vector decay
        J (int): Total spin resulting from the two Wigner D functions capturing the
            resonance decay
        M (int): Total m-projection resulting from the two Wigner D functions capturing
            the resonance decay
        waves (List[Wave]): List of all waves found in input file
    Returns:
        float: value of the moment
    """

    max_J = max(wave.spin for wave in waves)

    moment = 0
    # for loops are done to best match the mathematical notation
    for Ji in range(max_J + 1):
        for li in range(Ji + 1):
            for Jj in range(max_J + 1):
                factor = 1 / ((2 * Jj + 1) * 3)

                for lj in range(Jj + 1):
                    for mi in range(-Ji, Ji + 1):
                        if not any(  # skip calculation if wave doesn't exist
                            wave.spin == Ji and wave.l == li and wave.m == mi
                            for wave in waves
                        ):
                            continue
                        for mj in range(-Jj, Jj + 1):
                            if not any(  # skip calculation if wave doesn't exist
                                wave.spin == Jj and wave.l == lj and wave.m == mj
                                for wave in waves
                            ):
                                continue

                            sdme = calculate_SDME(alpha, Ji, li, mi, Jj, lj, mj, waves)
                            # if 0 then no need to calculate CG coefficients below
                            if sdme == 0.0:
                                continue

                            for lambda_i in range(-1, 2):
                                for lambda_j in range(-1, 2):
                                    cgs = (
                                        CG(li, 0, 1, lambda_i, Ji, lambda_i)
                                        * CG(lj, 0, 1, lambda_j, Jj, lambda_j)
                                        * CG(1, lambda_i, Jv, Lambda, 1, lambda_j)
                                        * CG(1, 0, Jv, 0, 1, 0)
                                        * CG(Ji, mi, J, M, Jj, mj)
                                        * CG(Ji, lambda_i, J, Lambda, Jj, lambda_j)
                                    )
                                    moment += factor * cgs * sdme

    return moment


@numba.njit
def sign(i):
    # replaces calculating costly (-1)^x powers in the SDMEs
    return 1 if i % 2 == 0 else -1


@numba.njit(cache=True)
def calculate_SDME(
    alpha: int, Ji: int, li: int, mi: int, Jj: int, lj: int, mj: int, waves: List[Wave]
) -> float:
    """Calculate the spin density matrix element using complex production coefficients

    The SDMEs are separated into 3 cases depending on the value of alpha (which indexes
    the intensity component). They each contain a sum over the two reflectivity values.
    The calculation is done by storing the 4 complex values in each sum and calculating
    the result.

    Args:
        alpha (int):indexes which term of the intensity the moment is associated with
            0: unpolarized
            1: polarized cos(2*Phi)
            2: polarized sin(2*Phi)
        Ji (int): 1st wave spin
        li (int): 1st wave angular momenta
        mi (int): 1st wave m-projection
        Jj (int): 2nd wave spin
        lj (int): 2nd wave angular momenta
        mj (int): 2nd wave spin
        waves (List[Wave]): List of all waves found in input file

    Raises:
        ValueError: if alpha value is not 0, 1, or 2

    Returns:
        float: SDME value for the corresponding alpha value
    """

    reflectivities = [-1, 1]
    result = 0

    match alpha:
        case 0:
            for e in reflectivities:
                for wave in waves:
                    conditions = {
                        "c1": wave.spin == Ji and wave.l == li and wave.m == mi,
                        "c2": wave.spin == Jj and wave.l == lj and wave.m == mj,
                        "c3": wave.spin == Ji and wave.l == li and wave.m == -mi,
                        "c4": wave.spin == Jj and wave.l == lj and wave.m == -mj,
                    }
                    if wave.reflectivity != e:
                        continue
                    if conditions["c1"]:
                        c1 = complex(wave.real, wave.imaginary)
                    if conditions["c2"]:
                        c2 = complex(wave.real, -wave.imaginary)
                    if conditions["c3"]:
                        c3 = complex(wave.real, wave.imaginary)
                    if conditions["c4"]:
                        c4 = complex(wave.real, -wave.imaginary)

                result += c1 * c2 + sign(mi + mj + li + lj + Ji + Jj) * c3 * c4
        case 1:
            for e in reflectivities:
                for wave in waves:
                    conditions = {
                        "c1": wave.spin == Ji and wave.l == li and wave.m == -mi,
                        "c2": wave.spin == Jj and wave.l == lj and wave.m == mj,
                        "c3": wave.spin == Ji and wave.l == li and wave.m == mi,
                        "c4": wave.spin == Jj and wave.l == lj and wave.m == -mj,
                    }
                    if wave.reflectivity != e:
                        continue
                    if conditions["c1"]:
                        c1 = complex(wave.real, wave.imaginary)
                    if conditions["c2"]:
                        c2 = complex(wave.real, -wave.imaginary)
                    if conditions["c3"]:
                        c3 = complex(wave.real, wave.imaginary)
                    if conditions["c4"]:
                        c4 = complex(wave.real, -wave.imaginary)

                result += e * (
                    sign(1 + mi + li + Ji) * c1 * c2 + sign(1 + mj + lj + Jj) * c3 * c4
                )
        case 2:
            for e in reflectivities:
                for wave in waves:
                    conditions = {
                        "c1": wave.spin == Ji and wave.l == li and wave.m == -mi,
                        "c2": wave.spin == Jj and wave.l == lj and wave.m == mj,
                        "c3": wave.spin == Ji and wave.l == li and wave.m == mi,
                        "c4": wave.spin == Jj and wave.l == lj and wave.m == -mj,
                    }
                    if wave.reflectivity != e:
                        continue
                    if conditions["c1"]:
                        c1 = complex(wave.real, wave.imaginary)
                    if conditions["c2"]:
                        c2 = complex(wave.real, -wave.imaginary)
                    if conditions["c3"]:
                        c3 = complex(wave.real, wave.imaginary)
                    if conditions["c4"]:
                        c4 = complex(wave.real, -wave.imaginary)

                result += e * (
                    sign(mi + li + Ji) * c1 * c2 - sign(mj + lj + Jj) * c3 * c4
                )
            result *= complex(0, 1)  # H2 is the purely imaginary moment
        case _:
            raise ValueError(f"Invalid alpha value {alpha}")

    return result


def parse_args() -> dict:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        nargs="+",
        help="Path to the best fit parameters input file(s)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="",
        help="Path to the output file",
    )
    parser.add_argument(
        "-s",
        "--sorted",
        type=bool,
        default=True,
        help=(
            "Sort the input files by last number in the file name or path. Defaults"
            " to True, so that the index of each csv row matches the ordering of the"
            " input files"
        ),
    )
    parser.add_argument(
        "-p",
        "--preview",
        action="store_true",
        help=("When passed, print out the files that will be processed and exit."),
    )

    return vars(parser.parse_args())


if __name__ == "__main__":
    main(parse_args())
