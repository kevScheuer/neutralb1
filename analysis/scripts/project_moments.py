"""Project vector-pseudoscalar PWA fit results to moments

This script uses the conversion from partial wave complex values to moments to "project"
the potentially ambiguous PWA fit results into unique moments. The input is .fit
file(s), which are the output of an AmpTools fit The moments are computed and saved to
an output csv file. Multiple files can be given as input, with optional sorting based on
the file name or path.

NOTE: This script assumes that the amplitudes are written in the vec-ps eJPmL format.
For example, the positive reflectivity, JP=1+, m=0, S-wave amplitude would be written
in the cfg file as [reaction]::RealNegSign::p1p0S.

If a custom scale parameter is used to multiply the complex values, the script assumes
that this parameter contains the word "scale" in the name.

TODO: The D waves aren't being added to H0(0,0,0,0) moment, potentially any other moment
TODO: handle free floating parameters like the D/S ratio
TODO: Parse Breit-Wigners instead of hard-coding them
TODO: The moment matrix is calculated for every file, but only needs to be done once
"""

import argparse
import itertools
import os
from typing import Dict, List, Set, TextIO, Tuple  # type hinting

import numba  # speed up for-loop calculations
import numpy as np
import pandas as pd
import pwa_tools
import spherical
import utils

BREIT_WIGNERS = {
    "1p": {"mass": 1.235, "width": 0.142},
    "1m": {"mass": 1.465, "width": 0.4},
}

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
    ("scale", numba.float64),
]


@numba.experimental.jitclass(spec)
class Wave:
    def __init__(self, name, reflectivity, spin, parity, m, l, real, imaginary, scale):
        self.name = name
        self.reflectivity = reflectivity
        self.spin = spin
        self.parity = parity
        self.m = m
        self.l = l
        self.real = real
        self.imaginary = imaginary
        self.scale = scale

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
    elif not args["output"]:
        args["output"] = "moments.csv"
    if not all(input_file.endswith(".fit") for input_file in args["input"]):
        raise ValueError("Input file(s) must be .fit files")

    # sort the input files if requested
    input_files = (
        utils.sort_input_files(args["input"]) if args["sorted"] else args["input"]
    )

    # only print out the files that will be processed if preview flag is passed
    if args["preview"]:
        print("Files that will be processed:")
        for file in input_files:
            print(f"\t{file}")
        return

    # create empty dictionaries for
    data_dict = {}  # The moments calculated for each file
    matrix_dict = {}  # The production coefficient pairs that contribute to each moment
    # since this next dict is passed to numba functions, it must be explicitly typed
    coefficient_dict = numba.typed.Dict.empty(
        key_type=numba.types.UniTuple(numba.types.unicode_type, 2),
        value_type=numba.types.float64,
    )

    # Loop to calculate moments for each input file
    for file in input_files:
        mass = get_mass(file)  # obtain center of mass bin for Breit Wigner calculation
        waves = get_waves(file, mass, args["breit_wigner"])  # obtain waves from file

        # PREPARE QUANTUM NUMBERS OF MOMENTS
        Jv_array = np.array([0, 2])  # CG coefficient always makes Jv=1 be 0
        # (-) Lambda values are directly proportional to (+) ones, no need to calculate
        Lambda_array = np.arange(0, 3)

        max_J = max(wave.spin for wave in waves)
        J_array = np.arange(0, 2 * max_J + 1)

        max_m = max(wave.m for wave in waves)
        M_array = np.arange(0, max_m + 1)  # like Lambda, -m ∝ +m moments

        # calculate each moment and what production coefficients contribute to it
        # and add to the dictionaries
        for alpha in range(3):
            for Jv, Lambda, J, M in itertools.product(
                Jv_array, Lambda_array, J_array, M_array
            ):
                moment_str = f"H{alpha}({Jv},{Lambda},{J},{M})"

                # this is a new moment, so create list to store results for each file
                if moment_str not in data_dict.keys():
                    data_dict[moment_str] = []

                coefficient_dict.clear()
                moment_val = calculate_moment(
                    alpha, Jv, Lambda, J, M, list(waves), coefficient_dict
                )
                # save the results for this moment
                matrix_dict[moment_str] = coefficient_dict.copy()
                data_dict[moment_str].append(moment_val)

    # check that every data file has the same number of moments
    if not all(
        len(data) == len(data_dict["H0(0,0,0,0)"]) for data in data_dict.values()
    ):
        raise ValueError("Number of moments calculated is not consistent across files")

    # save moments and matrix to csv files
    df = pd.DataFrame.from_dict(data_dict)
    df.index = input_files
    df.to_csv(args["output"], index_label="file")

    matrix_df = pd.DataFrame.from_dict(matrix_dict, orient="index").T
    matrix_df.fillna(0, inplace=True)
    matrix_df.to_csv(args["output"].replace(".csv", "_matrix.csv"))

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
    mass_range = ""
    for subdir in file.split("/"):
        if subdir.startswith("mass"):
            mass_range = subdir.split("_")[-1]
    if not mass_range:
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
    waves = set()
    searching_for_scale = False
    searching_for_amplitudes = False
    with open(file, "r") as f:
        for line in f:
            # general line filtering
            if (
                line.startswith("#")  # skip comments
                or not line.strip()  # skip empty lines
                or "isotropic" in line  # skip background wave
                or "Bkgd" in line
            ):
                continue

            # this section of .fit file contains the amplitude scale values
            if "Reactions, Amplitudes, and Scale Parameters" in line:
                searching_for_scale = True
                continue
            # this section of .fit file contains the amplitude re/im parts
            if "Parameter Values and Errors" in line:
                searching_for_scale = False
                searching_for_amplitudes = True
                continue
            # stop file scanning if we've gotten to the normalization integrals
            if "Normalization Integrals" in line:
                break

            if "::" not in line:  # skip lines without amplitudes
                continue

            parts = line.split()  # split line into its components

            if searching_for_scale:
                # if we've hit this section header, stop searching for the scales
                if "Likelihood Total and Partial Sums" in line:
                    searching_for_scale = False
                    continue
                if "scale" not in line:
                    continue

                amplitude = parts[0].split("::")[-1]
                # parse amplitude into its quantum numbers
                parsed_amp = pwa_tools.parse_amplitude(amplitude.split("_")[0])
                reflectivity = pwa_tools.char_to_int(parsed_amp["e"])
                spin = int(parsed_amp["j"])
                parity = pwa_tools.char_to_int(parsed_amp["p"])
                m = pwa_tools.char_to_int(parsed_amp["m"])
                l = pwa_tools.char_to_int(parsed_amp["l"])

                # add scale parameter and amplitude info to wave set
                waves.add(
                    Wave(
                        name=amplitude,
                        reflectivity=reflectivity,
                        spin=spin,
                        parity=parity,
                        m=m,
                        l=l,
                        real=np.nan,
                        imaginary=np.nan,
                        scale=float(parts[-1]),
                    )
                )

            # find the real and imaginary parts of the amplitudes
            if searching_for_amplitudes:
                amplitude_re_im = parts[0].split("::")[-1]  # form eJPmL_re or eJPmL_im
                amplitude = amplitude_re_im.split("_")[0]  # eJPmL
                re_im_flag = amplitude_re_im.split("_")[-1]  # re or im

                # obtain real and imaginary parts and add to appropriate wave
                wave = [w for w in waves if w.name == amplitude]
                assert len(wave) != 0, f"Amplitude {amplitude} not found in file {file}"
                assert len(wave) == 1, f"File {file} contains duplicates of {amplitude}"
                wave = wave[0]
                scaled_part = wave.scale * float(parts[-1])
                match re_im_flag:
                    case "re":
                        wave.real = scaled_part
                    case "im":
                        wave.imaginary = scaled_part
                    case _:
                        raise ValueError(
                            f"Unexpected amplitude format {amplitude} in file {file}"
                        )

    # If breit wigners are used, we need to multiply the re/im parts for each wave
    # get breit wigner if requested
    if use_breit_wigner:
        for wave in waves:
            try:
                bw_mass = BREIT_WIGNERS[wave.name[1:3]]["mass"]
                bw_width = BREIT_WIGNERS[wave.name[1:3]]["width"]
            except KeyError:
                print(
                    f"Breit-Wigner parameters not found for wave {wave.name}, skipping."
                )
                continue
            breit_wigner = pwa_tools.breit_wigner(mass, bw_mass, bw_width, wave.l)
            c = complex(wave.real, wave.imaginary) * breit_wigner
            wave.real = c.real
            wave.imaginary = c.imag

    # Check that real and imaginary parts were found for all waves
    if not all(w.real and w.imaginary for w in waves):
        raise ValueError(
            f"Real and imaginary parts not found for all waves in file {file}"
        )

    return waves


@numba.njit
def calculate_moment(
    alpha: int,
    Jv: int,
    Lambda: int,
    J: int,
    M: int,
    waves: List[Wave],
    coefficient_dict: Dict[Tuple[str, str], float],
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
        coefficient_dict (Dict[Tuple[str, str], float]): dictionary to store the
            clebsch gordan coefficient values for the moment and the corresponding
            production coefficient pairs. Has the form {(wave1, wave2): value}, where
            wave2 is implicitly complex conjugated, so (wave1, wave2) != (wave2, wave1)
    Returns:
        float: value of the moment
    """

    max_J = max([wave.spin for wave in waves])

    moment = 0
    # for loops are done here to best match the mathematical notation
    for Ji in range(max_J + 1):
        for li in range(Ji + 1):
            for Jj in range(max_J + 1):
                factor = 1 / ((2 * Jj + 1) * 3)
                for lj in range(Jj + 1):
                    for mi in range(-Ji, Ji + 1):
                        for mj in range(-Jj, Jj + 1):
                            # calculate the sdme and save the production coefficient
                            # pairs that contribute to it in the dictionary
                            sdme, pairs = calculate_SDME(
                                alpha,
                                Ji,
                                li,
                                mi,
                                Jj,
                                lj,
                                mj,
                                waves,
                            )
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
                                    # add the CG coefficients to every pair of
                                    # production coefficients that contribute to SDME
                                    if cgs != 0:
                                        for pair in pairs:
                                            pair_value = coefficient_dict.get(pair, 0.0)
                                            coefficient_dict[pair] = pair_value + cgs

                                    # finally, calculate the moment
                                    moment += factor * cgs * sdme
    return moment


@numba.njit(cache=True)
def sign(i):
    # replaces calculating costly (-1)^x powers in the SDMEs
    return 1 if i % 2 == 0 else -1


@numba.njit
def process_waves(
    waves: List[Wave],
    Ji: int,
    li: int,
    mi: int,
    Jj: int,
    lj: int,
    mj: int,
    e: int,
    alpha: int,
) -> Tuple[int, int, int, int, List[Tuple[str, str]]]:
    """Return the 4 complex values for the SDME calculation and the pairs of waves

    This function iterates over all waves to find the 4 complex values that are used in
    the SDME calculation. It also stores the pairs of production coefficients that
    contribute to the SDME in a list to be used later in the clebsch gordan matrix.
    This function has common code for the 3 cases of alpha, so it is extracted for
    readability and maintainability.

    Args:
        waves (List[Wave]): List of all waves found in input file
        Ji (int): 1st wave spin
        li (int): 1st wave angular momenta
        mi (int): 1st wave m-projection
        Jj (int): 2nd wave spin
        lj (int): 2nd wave angular momenta
        mj (int): 2nd wave spin
        e (int): reflectivity of the sum being performed
        alpha (int):indexes which term of the intensity the moment is associated with
            0: unpolarized
            1: polarized cos(2*Phi)
            2: polarized sin(2*Phi)

    Returns:
        int, int, int, int, List[Tuple[str, str]]: the 4 complex values for the SDME
            calculation and the pairs of production coefficients that contribute to the
            SDME.
    """

    c1, c2, c3, c4 = 0.0, 0.0, 0.0, 0.0
    pairs = []
    c1_name, c2_name, c3_name, c4_name = "", "", "", ""
    for wave in waves:
        conditions = {
            "c1": (
                wave.spin == Ji
                and wave.l == li
                and wave.m == (mi if alpha == 0 else -mi)
            ),
            "c2": (wave.spin == Jj and wave.l == lj and wave.m == mj),
            "c3": (
                wave.spin == Ji
                and wave.l == li
                and wave.m == (-mi if alpha == 0 else mi)
            ),
            "c4": (wave.spin == Jj and wave.l == lj and wave.m == -mj),
        }
        if wave.reflectivity != e:
            continue
        if conditions["c1"]:
            c1 = complex(wave.real, wave.imaginary)
            c1_name = wave.name
        if conditions["c2"]:
            c2 = complex(wave.real, -wave.imaginary)
            c2_name = wave.name
        if conditions["c3"]:
            c3 = complex(wave.real, wave.imaginary)
            c3_name = wave.name
        if conditions["c4"]:
            c4 = complex(wave.real, -wave.imaginary)
            c4_name = wave.name

        # if the pair of waves is used in the SDME, add it to the coefficient dict, but
        # only if it hasn't been added before
        if c1_name and c2_name:
            pairs.append((c1_name, c2_name))
        if c3_name and c4_name:
            pairs.append((c3_name, c4_name))
    return c1, c2, c3, c4, pairs


@numba.njit
def calculate_SDME(
    alpha: int,
    Ji: int,
    li: int,
    mi: int,
    Jj: int,
    lj: int,
    mj: int,
    waves: List[Wave],
) -> Tuple[complex, List[Tuple[str, str]]]:
    """Calculate the spin density matrix element using complex production coefficients

    The SDMEs are separated into 3 cases depending on the value of alpha (which indexes
    the intensity component). They each contain a sum over the two reflectivity values.
    The calculation is done by storing the 4 complex values in each sum and calculating
    the result. The pairs of production coefficients that contribute to this SDME are
    also stored in a list to be used later in the clebsch gordan matrix.

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
        coefficient_dict (Dict[Tuple[str, str], float]): see calculate_moment docstring

    Raises:
        ValueError: if alpha value is not 0, 1, or 2

    Returns:
        complex: SDME value for the corresponding alpha value
        List[Tuple[str, str]]: list of pairs of production coefficients that contribute
            to the SDME
    """

    reflectivities = [-1, 1]
    result = complex(0.0, 0.0)
    pairs = []

    # calculate the sdme according to the alpha value
    match alpha:
        case 0:
            for e in reflectivities:
                c1, c2, c3, c4, new_pairs = process_waves(
                    waves, Ji, li, mi, Jj, lj, mj, e, alpha
                )
                pairs.extend(new_pairs)
                result += c1 * c2 + sign(mi + mj + li + lj + Ji + Jj) * c3 * c4
        case 1:
            for e in reflectivities:
                c1, c2, c3, c4, new_pairs = process_waves(
                    waves, Ji, li, mi, Jj, lj, mj, e, alpha
                )
                pairs.extend(new_pairs)
                result += e * (
                    sign(1 + mi + li + Ji) * c1 * c2 + sign(1 + mj + lj + Jj) * c3 * c4
                )
        case 2:
            for e in reflectivities:
                c1, c2, c3, c4, new_pairs = process_waves(
                    waves, Ji, li, mi, Jj, lj, mj, e, alpha
                )
                pairs.extend(new_pairs)
                result += e * (
                    sign(mi + li + Ji) * c1 * c2 - sign(mj + lj + Jj) * c3 * c4
                )
            result *= complex(0, 1)  # H2 is the purely imaginary moment
        case _:
            raise ValueError(f"Invalid alpha value {alpha}")

    return result, pairs


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
    """Calculate the clebsch gordan coefficients for a generic moment"""
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
