"""Project vector-pseudoscalar PWA fit results to moments

This script uses the conversion from partial wave complex values to moments to "project"
the potentially ambiguous PWA fit results into unique moments. It takes as input the
best fit parameters, obtained using the "-s" flag in the fit command. The moments are
computed and saved to an output csv file. Multiple files can be given as input, with
optional sorting based on the file name or path.

NOTE: This script assumes that the amplitudes are written in the vec-ps eJPmL format.
For example, the positive reflectivity, JP=1+, m=0, S-wave amplitude would be written
in the cfg file as [reaction]::RealNegSign::p1p0S.

If a custom scale parameter is used to multiply the complex values, the script assumes
that this parameter contains the word "scale" in the name.

TODO: handle free floating parameters like the D/S ratio
"""

BREIT_WIGNERS = {
    "1p": {"mass": 1.235, "width": 0.142},
    "1m": {"mass": 1.465, "width": 0.4},
}

import argparse
import cmath
import itertools
import os
from typing import List, Set, TextIO  # type hinting

import numba  # speed up for loop calculations
import numpy as np
import pandas as pd
import pwa_tools
import spherical
import utils
from sympy.physics.quantum.cg import CG  # clebsch-gordan coefficients

numba.config.DISABLE_JIT = True  # uncomment for debugging

# define the structure of the Wave class for numba to compile
spec = [
    ("name", numba.types.unicode_type),
    ("reflectivity", numba.int32),
    ("spin", numba.int32),
    ("parity", numba.int32),
    ("l", numba.int32),
    ("m", numba.int32),
    ("real", numba.float64),
    ("imaginary", numba.float64),
]


@numba.experimental.jitclass(spec)
class Wave:
    def __init__(self, name, reflectivity, spin, parity, m, l, real, imaginary):
        self.name = name
        self.reflectivity = reflectivity
        self.spin = spin
        self.parity = parity
        self.m = m
        self.l = l
        self.real = real
        self.imaginary = imaginary

    def __eq__(self, other):
        if isinstance(other, Wave):
            return self.name == other.name
        return False

    def __hash__(self):
        return hash(self.name)


def main(args: dict) -> None:

    # Error / default value handling
    for input_file in args["input"]:
        if not os.path.isfile(input_file):
            raise ValueError(f"Input file {input_file} does not exist")
    if args["output"] and not args["output"].endswith(".csv"):
        args["output"] += ".csv"
    else:
        args["output"] = "moments.csv"
    if not all(input_file.endswith(".txt") for input_file in args["input"]):
        raise ValueError("Input file(s) must be .txt files")

    input_files = (
        utils.sort_input_files(args["input"]) if args["sorted"] else args["input"]
    )

    if args["preview"]:
        print("Files that will be processed:")
        for file in input_files:
            print(f"\t{file}")
        return

    # create an empty dictionary to store the moments
    data_dict = {}

    # Loop to calculate moments for each input file
    for file in input_files:
        mass = get_mass(file)  # obtain center of mass bin for Breit Wigner values
        waves = get_waves(file, mass, args["breit_wigner"])  # obtain waves from file

        # Prepare quantum numbers of moments
        Jv_array = np.array([0, 2])  # CG coefficient always makes Jv=1 be 0
        # (-) Lambda values are directly proportional to (+) ones, no need to calculate
        Lambda_array = np.arange(0, 3)

        max_J = max(wave.spin for wave in waves)
        J_array = np.arange(0, 2 * max_J + 1)

        max_m = max(wave.m for wave in waves)
        M_array = np.arange(0, max_m + 1)  # like Lambda, -m ∝ +m moments

        # calculate each moment and add it to the dataframe
        for alpha in range(3):
            for Jv, Lambda, J, M in itertools.product(
                Jv_array, Lambda_array, J_array, M_array
            ):
                moment_str = f"H{alpha}({Jv},{Lambda},{J},{M})"
                moment_val = calculate_moment(alpha, Jv, Lambda, J, M, list(waves))

                # if moment string not in dictionary, add it with an empty list
                if moment_str not in data_dict.keys():
                    data_dict[moment_str] = []
                data_dict[moment_str].append(moment_val)

    # check that every data file has the same number of moments
    if not all(
        len(data) == len(data_dict["H0(0,0,0,0)"]) for data in data_dict.values()
    ):
        raise ValueError("Number of moments calculated is not consistent across files")

    # save moments to a csv file, with the input files as the index
    df = pd.DataFrame.from_dict(data_dict)
    df.index = input_files
    df.to_csv(args["output"], index_label="file")

    pass


def get_mass(file: str) -> float:
    """Obtain the mass bin's center from a subdirectory within the file path

    TODO: this is a very setup-dependent way to get this info, but only other
        option would be to implement this all in c++ to use the FitResults class

    Args:
        file (str): full file path assumed to be setup like "path/to/mass_0-100/fit.txt"
            or "path/to/mass_0-100/subdirectory/fit.txt"

    Raises:
        ValueError: since only 2 cases are handled above, this error is raised if
            the mass range cannot be obtained from the file path

    Returns:
        float: average of the two values found in the mass range
    """

    if "mass" in file.split("/")[-3]:
        mass_range = file.split("/")[-3].split("_")[-1]
    elif "mass" in file.split("/")[-2]:
        mass_range = file.split("/")[-2].split("_")[-1]
    else:
        raise ValueError("Mass range could not be obtained from file path")
    mass = (float(mass_range.split("-")[0]) + float(mass_range.split("-")[1])) / 2.0

    return mass


def get_waves(file: TextIO, mass: float, use_breit_wigner: bool) -> Set[Wave]:
    """Obtain the set of waves, with their real and imaginary parts, from a fit result

    Args:
        file (TextIO): best fit parameters file, obtained using the "-s" flag on the fit
            command, that contains the real and imaginary parts of the best fit result
        mass (float): center of mass bin for the Breit Wigner values
        use_breit_wigner (bool): whether to modify the production coefficients by the
            appropriate Breit-Wigner values
    Returns:
        Set[Wave]: all the waves and their information found in the file
    """

    scale = 1.0
    # scan the file for an optional scale parameter that multiplies the complex values
    with open(file, "r") as f:
        for line in f:
            if "parameter" in line and "scale" in line:
                scale = float(line.split()[-1])
                break

    waves = set()
    with open(file, "r") as f:
        for line in f:
            if (
                line.startswith("#")  # skip comments
                or not line.strip()  # skip empty lines
                or "isotropic" in line  # skip background wave
                or "Bkgd" in line
                or "parameter" in line  # skip scale parameter
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

            # get breit wigner if requested
            breit_wigner = 1.0
            if use_breit_wigner and amplitude[1:3] in BREIT_WIGNERS.keys():
                bw_mass = BREIT_WIGNERS[amplitude[1:3]]["mass"]
                bw_width = BREIT_WIGNERS[amplitude[1:3]]["width"]
                breit_wigner = pwa_tools.breit_wigner(mass, bw_mass, bw_width, l)

            # obtain real and imaginary parts
            match parts[2]:
                case "cartesian":
                    prod_coeff = complex(
                        scale * float(parts[3]), scale * float(parts[4])
                    )
                    prod_coeff *= breit_wigner
                case "polar":
                    # TODO: check if scale parameter multiplies the polar form, or if
                    # AmpTools scales the cartesian form
                    prod_coeff = cmath.rect(
                        scale * float(parts[3]), scale * float(parts[4])
                    )
                    prod_coeff *= breit_wigner
                case _:
                    raise ValueError(f"Unexpected complex number format in file {file}")
            re, im = prod_coeff.real, prod_coeff.imag
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


@numba.njit
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

    max_J = max([wave.spin for wave in waves])

    moment = 0
    # for loops are done to best match the mathematical notation
    for Ji in range(max_J + 1):
        for li in range(Ji + 1):
            for Jj in range(max_J + 1):
                factor = 1 / ((2 * Jj + 1) * 3)
                for lj in range(Jj + 1):
                    for mi in range(-Ji, Ji + 1):
                        for mj in range(-Jj, Jj + 1):
                            sdme = calculate_SDME(alpha, Ji, li, mi, Jj, lj, mj, waves)
                            # if 0 then no need to calculate CG coefficients below
                            if sdme.real == 0.0 and sdme.imag == 0.0:
                                continue

                            for lambda_i in range(-1, 2):
                                for lambda_j in range(-1, 2):
                                    cgs = calculate_clebsch_gordans(
                                        Jv,
                                        Lambda,
                                        J,
                                        M,
                                        Ji,
                                        li,
                                        mi,
                                        lambda_i,
                                        Jj,
                                        lj,
                                        mj,
                                        lambda_j,
                                    )
                                    moment += factor * cgs * sdme

    return moment


@numba.njit(cache=True)
def sign(i):
    # replaces calculating costly (-1)^x powers in the SDMEs
    return 1 if i % 2 == 0 else -1


@numba.njit
def calculate_SDME(
    alpha: int, Ji: int, li: int, mi: int, Jj: int, lj: int, mj: int, waves: List[Wave]
) -> complex:
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
    result = complex(0.0, 0.0)

    match alpha:
        case 0:
            for e in reflectivities:
                c1, c2, c3, c4 = 0.0, 0.0, 0.0, 0.0
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
                c1, c2, c3, c4 = 0.0, 0.0, 0.0, 0.0
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
                c1, c2, c3, c4 = 0.0, 0.0, 0.0, 0.0
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


@numba.njit(cache=True)
def calculate_clebsch_gordans(
    Jv: int,
    Lambda: int,
    J: int,
    M: int,
    Ji: int,
    li: int,
    mi: int,
    lambda_i: int,
    Jj: int,
    lj: int,
    mj: int,
    lambda_j: int,
) -> float:
    """Calculate the clebsch gordan coefficients for a generic moment

    Args:
        Jv (_type_): _description_
        Lambda (_type_): _description_
        J (_type_): _description_
        M (_type_): _description_
        Ji (_type_): _description_
        li (_type_): _description_
        mi (_type_): _description_
        lambda_i (_type_): _description_
        Jj (_type_): _description_
        lj (_type_): _description_
        mj (_type_): _description_
        lambda_j (_type_): _description_

    Returns:
        float: _description_
    """
    cgs = (
        spherical.clebsch_gordan(li, 0, 1, lambda_i, Ji, lambda_i)
        * spherical.clebsch_gordan(lj, 0, 1, lambda_j, Jj, lambda_j)
        * spherical.clebsch_gordan(1, lambda_i, Jv, Lambda, 1, lambda_j)
        * spherical.clebsch_gordan(1, 0, Jv, 0, 1, 0)
        * spherical.clebsch_gordan(Ji, mi, J, M, Jj, mj)
        * spherical.clebsch_gordan(Ji, lambda_i, J, Lambda, Jj, lambda_j)
    )

    return cgs


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
    parser.add_argument(
        "-b",
        "--breit-wigner",
        action="store_true",
        help=(
            "When passed, modify the production coefficients by the appropriate"
            " Breit-Wigner values. Note these are currently hard-coded in the script."
        ),
    )

    return vars(parser.parse_args())


if __name__ == "__main__":
    main(parse_args())
