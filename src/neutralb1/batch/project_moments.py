"""DEPRECATED Project vector-pseudoscalar PWA fit results to moments

This script uses the conversion from partial wave complex values to "project" PWA fit
results into unique moments. The input is .fit file(s), which are the output of an
AmpTools fit. The moments are computed and saved to an output csv file. Multiple files
can be given as input, with optional sorting based on the file name or path.

NOTE: This script assumes that the amplitudes are written in the vec-ps eJPmL format.
For example, the positive reflectivity, JP=1+, m=0, S-wave amplitude would be written
in the cfg file as [reaction]::RealNegSign::p1p0S. If you have a different format, then
you'll have to account for it when parsing the amplitude name in the get_waves function.

TODO: handle free floating parameters like the D/S ratio
TODO: Parse Breit-Wigners instead of hard-coding them
TODO: Coefficient table is calculated for every file, but only the last one is saved
    This makes the assumption that every file has the same set of waves, and
    unnecessarily increases the time taken to process the files. Saving the values to
    build the table is deep in the calculation though, all the way down to the SDMEs,
    so a bool flag is needed through every function to avoid this redundancy. Or I could
    simply save the table for every file, but this could create a massive set of files
"""

import argparse
import concurrent.futures
import functools
import itertools
import os
import re
import sys
import timeit
from typing import Dict, List, Tuple

import numba  # speed up for-loop calculations
import numpy as np
import pandas as pd
import spherical

from src.neutralb1 import utils

# imports below need the path to work
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from analysis.scripts import pwa_tools  # type: ignore[import]

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


@numba.experimental.jitclass(spec)  # type: ignore
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


def main(args: dict) -> None:
    start_time = timeit.default_timer()
    # ===CHECK arguments===
    if args["output"] and not args["output"].endswith(".csv"):
        args["output"] += ".csv"
    elif not args["output"]:
        args["output"] = "moments.csv"

    # limit the number of workers to a max of 10 to prevent hogging resources
    args["workers"] = min(args["workers"], 10)

    # Check if args["input"] is a file containing a list of result files
    input_files = []
    if (
        len(args["input"]) == 1
        and os.path.isfile(args["input"][0])
        and not args["input"][0].endswith(".fit")
    ):
        with open(args["input"][0], "r") as file:
            input_files = [line.strip() for line in file if line.strip()]
    else:
        input_files = args["input"]

    # Ensure all input files exist and are .fit files
    if input_files == []:
        raise ValueError("No input files provided.")
    for f in input_files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"File {f} does not exist")
    if not all(f.endswith(".fit") for f in input_files):
        raise ValueError("Input file(s) must be .fit files")

    # sort the input files if requested
    input_files = (
        utils.sort_input_files(input_files, args["sort_index"])
        if args["sorted"]
        else input_files
    )

    # only print out the files that will be processed if preview flag is passed
    if args["preview"]:
        print("Files that will be processed:")
        for file in input_files:
            print(f"\t{file}")
        return

    # ===PROCESS each file in parallel===
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args["workers"]
    ) as executor:
        results = list(
            executor.map(process_file, input_files, [args] * len(input_files))
        )
    if results == []:
        raise ValueError("No results returned. Check input files.")

    # ===COLLECT each files results into dictionaries===
    moment_dict = {}
    table_dict = {}
    all_keys = set()

    # first lets collect all the moment strings
    for file_moments, _ in results:
        all_keys.update(file_moments.keys())

    # initialize the dictionary with empty lists for each moment
    for key in all_keys:
        moment_dict[key] = []

    # now we'll populate the moment and table dictionaries
    for file_moments, file_table in results:
        for key in all_keys:
            if key in file_moments:
                moment_dict[key].append(file_moments[key])
            else:
                moment_dict[key].append(complex(0, 0))
        table_dict = file_table

    if table_dict == {}:
        raise ValueError("No production coefficient pairs found. Check input files.")

    # save moments to csv file
    df = pd.DataFrame.from_dict(moment_dict)
    df = df.set_index(pd.Index(input_files))
    df.to_csv(args["output"], index_label="file")
    print("Moments saved to", args["output"])

    # save a table of the production coefficient pairs that contribute to each moment,
    # with the clebsch-gordan values that multiply them in the cells
    table_df = pd.DataFrame.from_dict(table_dict, orient="index").T
    table_df.fillna(0, inplace=True)
    table_df.reset_index(names=["wave1", "wave2"]).to_csv(
        args["output"].replace(".csv", "_table.csv"), index=False
    )
    print("Moment table saved to", args["output"].replace(".csv", "_table.csv"))

    end_time = timeit.default_timer()
    if args["verbose"]:
        print(f"Total time taken: {end_time - start_time:.4f} seconds")
    pass


def process_file(
    file: str, args: dict
) -> Tuple[Dict[str, float], Dict[str, Dict[Tuple[str, str], float]]]:
    """Process a single file to calculate the moments and production coefficient pairs

    Args:
        file (str): path to the .fit file to process
        args (dict): dictionary of arguments passed to the script

    Returns:
        Tuple[Dict[str, float], Dict[str, Dict[Tuple[str,str], float]]]:
            moments and their values
            the clebsch-gordan values that multiply the production coefficient pairs
                who contribute to each moment
    """

    # initialize the dictionary to store the production coefficient pairs and their CGs
    # we have to init it here outside the function because jit doesn't support
    # type-expression in the function
    coefficient_dict = numba.typed.Dict.empty(  # type:ignore
        key_type=numba.types.UniTuple(numba.types.unicode_type, 2),
        value_type=numba.types.float64,
    )

    if args["breit_wigner"] and not args["mass"]:
        # obtain center of mass bin for BW calculation. Split off the file name to make
        # use of function cache for repeated mass values
        mass = get_mass(file.rsplit("/", 1)[0])
    else:
        mass = args["mass"]

    # find all the partial waves in the file
    waves = get_waves(file, mass, args["breit_wigner"])

    start_time = timeit.default_timer()
    # ===PREPARE quantum numbers of moments===
    Jv_array = np.array([0, 2])  # CG coefficient always makes Jv=1 be 0
    # (-) Lambda values are directly proportional to (+) ones, no need to calculate
    Lambda_array = np.arange(0, 3)

    max_J = max(wave.spin for wave in waves)
    J_array = np.arange(0, 2 * max_J + 1)

    max_m = max(wave.m for wave in waves)
    M_array = np.arange(0, 2 * max_m + 1)  # like Lambda, -m ∝ +m moments

    # ===CALCULATE each moment and what production coefficients contribute to it===
    moment_dict = {}
    table_dict = {}
    for alpha in range(3):
        for Jv, Lambda, J, M in itertools.product(
            Jv_array, Lambda_array, J_array, M_array
        ):
            # moments with these quantum #s will always be 0, because their Wigner D
            # functions are 0
            if M > J or Lambda > J or Lambda > Jv:
                continue
            # Wishart section 5.10.2 proves H2(Jv,0,J,0) = 0 for any Jv, J
            if alpha == 2 and Lambda == 0 and M == 0:
                continue

            moment_str = f"H{alpha}_{Jv}{Lambda}{J}{M}"

            coefficient_dict.clear()
            moment_val = calculate_moment(
                alpha, Jv, Lambda, J, M, numba.typed.List(waves), coefficient_dict  # type: ignore
            )
            # save the results for this moment. The dictionary is copied to avoid
            # reference issues, and converted back to a python dict so it is pickleable
            # for the parallel processing
            table_dict[moment_str] = dict(coefficient_dict.copy())
            moment_dict[moment_str] = moment_val

    elapsed = timeit.default_timer() - start_time
    if args["verbose"]:
        print(f"Time taken to process {file}: {elapsed:.4f} seconds")

    return moment_dict, table_dict


@functools.cache
def get_mass(file: str) -> float:
    """Obtain the mass bin's center from a subdirectory within the file path

    TODO: this is a very setup-dependent way to get this info, but only other
        option would be to load the ROOT data file in the .fit file and plot its mass
        histogram and get the bin center.

    Args:
        file (str): full file path of .fit file, assumes that mass range is somewhere
            in the path as "subdir/mass_X-Y/"

    Raises:
        ValueError: raised if the mass range could not be obtained from the file path

    Returns:
        float: average of the two values found in the mass range
    """
    match = re.search(r"/mass_([\d.]+)-([\d.]+)", file)
    if not match:
        raise ValueError("Mass range could not be obtained from file path")

    mass = (float(match.group(1)) + float(match.group(2))) / 2.0
    return mass


def get_waves(file: str, mass: float, use_breit_wigner: bool) -> List[Wave]:
    """Obtain the set of waves, with their real and imaginary parts, from a fit result

    Args:
        file (str): best fit parameters file, obtained using the "-s" flag on the fit
            command, that contains the real and imaginary parts of the best fit result
        mass (float): center of mass bin, only used for the Breit Wigner calculations
        use_breit_wigner (bool): whether to modify the production coefficients by the
            appropriate Breit-Wigner values
    Returns:
        Set[Wave]: all the waves and their information found in the file
    """
    waves = {}
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

                amplitude = parts[0].split("::")[-1]
                # parse amplitude into its quantum numbers
                parsed_amp = pwa_tools.parse_amplitude(amplitude.split("_")[0])
                reflectivity = pwa_tools.char_to_int(parsed_amp["e"])
                spin = int(parsed_amp["j"])
                parity = pwa_tools.char_to_int(parsed_amp["p"])
                m = pwa_tools.char_to_int(parsed_amp["m"])
                l = pwa_tools.char_to_int(parsed_amp["l"])

                # add scale parameter and amplitude info to wave set
                waves[amplitude] = Wave(
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

            # find the real and imaginary parts of the amplitudes
            if searching_for_amplitudes:
                amplitude_re_im = parts[0].split("::")[-1]  # form eJPmL_re or eJPmL_im
                amplitude = amplitude_re_im.split("_")[0]  # eJPmL
                re_im_flag = amplitude_re_im.split("_")[-1]  # re or im

                # obtain real and imaginary parts and add to appropriate wave
                wave = waves.get(amplitude)
                if wave is None:
                    raise ValueError(f"Amplitude {amplitude} not found in file {file}")
                scaled_part = wave.scale * float(parts[-1])

                # modify component by the barrier factor
                barrier = pwa_tools.barrier_factor(mass, 0.1349768, 0.78266, wave.l)
                scaled_part *= barrier

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
    if use_breit_wigner:
        for wave in waves.values():
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
    if not all(w.real and w.imaginary for w in waves.values()):
        raise ValueError(
            f"Real and imaginary parts not found for all waves in file {file}"
        )

    return list(waves.values())


@numba.njit(cache=True)
def calculate_moment(
    alpha: int,
    Jv: int,
    Lambda: int,
    J: int,
    M: int,
    waves: List[Wave],
    coefficient_dict: Dict[Tuple[str, str], float],
) -> complex:
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
        for li in range(Ji + 2):
            for Jj in range(max_J + 1):
                norm = 1 / ((2 * Jj + 1) * 3)
                for lj in range(Jj + 2):
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
                                    moment += norm * cgs * sdme
    return moment


@numba.njit(cache=True)
def sign(i):
    # replaces calculating costly (-1)^x powers in the SDMEs
    return 1 if i % 2 == 0 else -1


@numba.njit(cache=True)
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
    also stored in a list to be used later in the clebsch gordan table.

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
        complex: SDME value for the corresponding alpha value
        List[Tuple[str, str]]: list of pairs of production coefficients that contribute
            to the SDME
    """
    reflectivities = [-1, 1]
    result = complex(0.0, 0.0)
    pairs = []

    # calculate the sdme according to the alpha value
    for e in reflectivities:
        c1, c2, c3, c4, new_pairs = process_waves(
            waves, Ji, li, mi, Jj, lj, mj, e, alpha
        )
        pairs.extend(new_pairs)
        match alpha:
            case 0:
                result += c1 * c2 + sign(mi + mj + li + lj + Ji + Jj) * c3 * c4
            case 1:
                result += e * (
                    sign(1 + mi + li + Ji) * c1 * c2 + sign(1 + mj + lj + Jj) * c3 * c4
                )
            case 2:
                result += e * (
                    sign(mi + li + Ji) * c1 * c2 - sign(mj + lj + Jj) * c3 * c4
                )
            case _:
                raise ValueError(f"Invalid alpha value {alpha}")
    if alpha == 2:
        result *= complex(0, 1)
    return result, pairs


@numba.njit(cache=True)
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
) -> Tuple[complex, complex, complex, complex, List[Tuple[str, str]]]:
    """Return the 4 complex values for the SDME calculation and the pairs of waves

    This function iterates over all waves to find the 4 complex values that are used in
    the SDME calculation. It also stores the pairs of production coefficients that
    contribute to the SDME in a list to be used later in the clebsch gordan table.
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

    c1, c2, c3, c4 = complex(0, 0), complex(0, 0), complex(0, 0), complex(0, 0)
    pairs = []
    c1_name, c2_name, c3_name, c4_name = "", "", "", ""
    for wave in waves:
        if wave.reflectivity != e:
            continue

        wave_J = wave.spin
        wave_l = wave.l
        wave_m = wave.m

        if wave_J == Ji and wave_l == li and wave_m == (mi if alpha == 0 else -mi):
            c1 = complex(wave.real, wave.imaginary)
            c1_name = wave.name
        if wave_J == Jj and wave_l == lj and wave_m == mj:
            c2 = complex(wave.real, -wave.imaginary)
            c2_name = wave.name
        if wave_J == Ji and wave_l == li and wave_m == (-mi if alpha == 0 else mi):
            c3 = complex(wave.real, wave.imaginary)
            c3_name = wave.name
        if wave_J == Jj and wave_l == lj and wave_m == -mj:
            c4 = complex(wave.real, -wave.imaginary)
            c4_name = wave.name

        # if both waves are found, add the pair to the list
        if c1_name and c2_name:
            pairs.append((c1_name, c2_name))
        if c3_name and c4_name:
            pairs.append((c3_name, c4_name))
    return c1, c2, c3, c4, pairs


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
        "-m",
        "--mass",
        type=float,
        default=None,
        help=(
            "Center of mass bin to use for the Breit-Wigner calculations. Defaults to"
            " None, meaning the mass will be obtained from the path"
        ),
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
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print out additional information during the script execution",
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=1,
        help=(
            "Number of workers to use for parallel processing. Defaults to 1, which"
            " means no parallel processing"
        ),
    )
    parser.add_argument(
        "--sort-index",
        type=int,
        default=-1,
        help=(
            "Determines what number in the file path is used for sorting. Defaults to"
            " -1, so that the last number in the path is used. See "
            " utils.sort_input_files for details."
        ),
    )

    return vars(parser.parse_args())


if __name__ == "__main__":
    main(parse_args())
