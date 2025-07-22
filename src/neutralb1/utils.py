"""A collection of some utility functions"""

# TODO: Many of these functions can be combined into a pythonized version of the
#   AmplitudeParser class

# TODO: watch this video https://www.youtube.com/watch?v=zgbUk90aQ6A to reconsider
#   how to handle pandas dataframes like wrap_phases effectively. My main concern is
#   duplicating a bunch of temporary dataframes leading to memory issues.

import itertools
import os
import pathlib
import re
import sys
from typing import Dict, List

import numpy as np
import pandas as pd


def sort_input_files(input_files: list, position: int = -1) -> list:
    """Sort the input files based off the last number in the file name or path

    Args:
        input_files (list): input files to be sorted
        position (int, optional): Index position of the number to be sorted on in the
            full path. Defaults to -1, meaning the last number is used for sorting. Be
            careful using this, as it will assume all path names have the same amount of
            distinct numbers, and thus the same indices.

    Returns:
        list: sorted list of files
    """

    def extract_last_number(full_path: str) -> float:
        numbers = re.findall(r"(?:\d*\.*\d+)", full_path)
        return float(numbers[position]) if numbers else float("inf")

    return sorted(input_files, key=extract_last_number)


def load_environment() -> None:
    """Load the shell environment variables from the .env file

    This function reads the `.env` file located in the `config` directory of the
    workspace directory, which is expected to be named `neutralb1`. It sets the
    environment variables defined in that file into the current Python environment.

    Raises:
        FileNotFoundError: If the `.env` file does not exist in the expected location.
    """

    workspace_dir = get_workspace_dir()

    if not os.path.exists(f"{workspace_dir}/config/.env"):
        raise FileNotFoundError(
            f"Environment file not found at {workspace_dir}/config/.env. "
            "Please execute `make update-env` to generate it."
        )

    env_vars = {}
    with open(f"{workspace_dir}/config/.env", "r") as env_file:
        for line in env_file:
            # Skip comments and empty lines
            if line.strip() and not line.startswith("#"):
                key, value = line.strip().split("=", 1)
                env_vars[key] = value

    os.environ.update(env_vars)
    sys.path.insert(0, workspace_dir)

    return


def get_workspace_dir() -> str:
    """Get the root directory of the workspace.

    Assumes the workspace directory is named "neutralb1" and the file calling this
    function is located somewhere within the workspace directory structure.

    Returns:
        str: The path of the workspace directory.
        Raises FileNotFoundError if the "neutralb1" directory is not found in the
    Raises:
        FileNotFoundError: If the "neutralb1" directory is not found in the path.
    """

    current_dir = pathlib.Path.cwd()
    workspace_dir = None
    while current_dir != current_dir.parent:
        if (current_dir / "neutralb1").is_dir():
            workspace_dir = str(current_dir / "neutralb1")

        current_dir = current_dir.parent
    if workspace_dir is None:
        raise FileNotFoundError("Could not find 'neutralb1' directory in the path.")

    return workspace_dir


def parse_amplitude(amp: str) -> Dict[str, str]:
    """
    Parse 'eJPmL' style amplitude string into its individual quantum numbers.

    Extracts reflectivity (e), total angular momentum (J), parity (P), m-projection (m),
    and orbital angular momentum (L) from amplitude strings like 'p1p0S'. Can also
    handle coherent sum variations like "m1m" or "p".

    Args:
        amp (str): Amplitude string, e.g. 'p1p0S', 'm1mpP', etc.

    Returns:
        dict: Dictionary with keys 'e', 'J', 'P', 'm', 'L'. Values are extracted from
        the input string, or "" if not found.
    """

    result_dict = {"e": "", "J": "", "P": "", "m": "", "L": ""}
    regex_dict = {
        "e": r"(?<![0-9]{1}[pm]{1})(?<!\d)[pm]",
        "JP": r"[0-9]{1}[pm]{1}",  # j and p always appear together
        "m": r"([0-9]{1}[pm]{1})+[lnm0pqr]",  # assumes "m" is always proceeded by "JP"
        "L": r"[SPDF]",
    }

    for quantum_number, expression in regex_dict.items():
        search = re.search(expression, amp)
        if search:
            if quantum_number == "m":
                # the search actually returns "JPm", so grab the last char if found
                result_dict["m"] = search.group()[-1]
            elif quantum_number == "JP":
                result_dict["J"] = search.group()[0]
                result_dict["P"] = search.group()[1]
            else:
                result_dict[quantum_number] = search.group()

    return result_dict


# TODO: this function modifies the dataframe in place, consider performance-efficient
#   alternative for returning a new dataframe
def wrap_phases(obj: pd.DataFrame | pd.Series) -> None:
    """Wrap phase differences to be from (-pi, pi] & convert from radians to degrees

    Accepts either a pandas DataFrame or Series. For DataFrames, only phase difference
    columns and their errors are wrapped. For Series, all values are wrapped.

    Args:
        obj (pd.DataFrame | pd.Series): DataFrame of fit results or Series of phase
            difference values.
    Raises:
        ValueError: If obj is not a DataFrame or Series.

    Returns:
        None: Edits the df or series in place
    """

    # wraps phase (in radians) to -pi < x <= pi and convert to degrees
    def wrap(phase):
        return np.rad2deg(np.angle(np.exp(1j * phase)))

    if isinstance(obj, pd.Series):
        obj.update(obj.apply(wrap))
    elif isinstance(obj, pd.DataFrame):
        phase_diffs = get_phase_differences(obj)
        for col in set(phase_diffs.values()):
            obj[col] = obj[col].apply(wrap)
            obj[f"{col}_err"] = obj[f"{col}_err"].apply(wrap)
    else:
        raise ValueError("Input must be a pandas DataFrame or Series.")

    return


def get_coherent_sums(df: pd.DataFrame) -> Dict[str, List[str]]:
    """Returns a dict of coherent sums from a fit results dataframe

    Args:
        df (pd.DataFrame): dataframe of fit results loaded from csv
    Returns:
        dict: Is of the form {Coherent sum string: set(amplitudes)}
        e.g. {"eJPmL": ["p1p0S", "m1mpP", ...]}
    """
    # create an empty list for every type of coherent sum
    sum_types = ["eJPmL", "JPmL", "eJPL", "JPL", "eJP", "JP", "e"]
    coherent_sums = {d: set() for d in sum_types}

    # grab all eJPml columns
    for column in df.columns:
        # skip the phase difference columns and any status columns
        if "_" in column or len(column) > 5:
            continue

        result_dict = parse_amplitude(column)
        # only add to key if all elements of key are in the column
        for key in coherent_sums.keys():
            split_key = list(key.lower())
            if any(result_dict[char] == "" for char in split_key):
                continue
            coh_sum = "".join([result_dict[char] for char in split_key])
            coherent_sums[key].add(coh_sum)

    return {k: sorted(v) for k, v in coherent_sums.items()}  # sort each set


def get_phase_differences(df: pd.DataFrame) -> Dict[tuple, str]:
    """Returns dict of all the phase difference columns in the dataframe

    The keys are specifically like (eJPmL_1, eJPmL_2), and there is no way to
    know how those will be ordered in the dataframe a priori. We avoid this by
    creating keys for either ordering and setting both of their values to the order
    found in the dataframe.

    Args:
        df (pd.DataFrame): dataframe of fit results loaded from csv
    Returns:
        dict: key = tuple of both possible combinations of every amplitude, val = the
            phase difference combination found in the dataframe
    """
    phase_differences = {}

    # get all possible combinations of eJPmL columns, and add their phase difference
    # column name if it exists (handles reverse ordering i.e. p1ppS_p1pmS & p1pmS_p1ppS)
    all_combos = list(itertools.combinations(get_coherent_sums(df)["eJPmL"], 2))
    columns = df.columns.values.tolist()
    for combo in all_combos:
        name = "_".join(combo)
        reverse_name = "_".join(reversed(combo))

        if name in columns:
            phase_differences[combo] = name
            phase_differences[tuple(reversed(combo))] = name
        if reverse_name in columns:
            phase_differences[combo] = reverse_name
            phase_differences[tuple(reversed(combo))] = reverse_name

    return phase_differences


def convert_amp_name(input_string: str) -> str:
    """Converts amplitude string to J^P L_m^(e) LaTeX style string

    Function can handle both amplitudes and phase differences. If input_string is not
    of eJPmL format, or subset of it i.e. eJPL, then the output will be undefined.

    Args:
        input_string (str): string in eJPmL format
    Returns:
        str: Prettier LaTeX style string in J^P L_m^(e) format. If it's a phase
            difference it's then J^P L_m^(e) - J^P L_m^(e)
    """

    pm_dict = {"m": "-", "p": "+", "": ""}
    m_proj_dict = {"n": -2, "m": -1, "0": 0, "p": +1, "q": +2, "": ""}

    # CASE 1: phase difference, always of form 'eJPmL_eJPmL'
    if "_" in input_string:
        amps = input_string.split("_")
        if len(amps) != 2:
            raise ValueError("Phase difference must be in 'eJPmL_eJPmL' format!")
        e1, j1, p1, m1, l1 = (
            pm_dict[amps[0][0]],
            amps[0][1],
            pm_dict[amps[0][2]],
            m_proj_dict[amps[0][3]],
            amps[0][-1],
        )
        e2, j2, p2, m2, l2 = (
            pm_dict[amps[1][0]],
            amps[1][1],
            pm_dict[amps[1][2]],
            m_proj_dict[amps[1][3]],
            amps[1][-1],
        )

        return (
            rf"${j1}^{{{p1}}}{l1}_{{{m1}}}^{{({e1})}}"
            rf" - {j2}^{{{p2}}}{l2}_{{{m2}}}^{{({e2})}}$"
        )

    # CASE 2: typical amplitude coherent sum
    amp_dict = parse_amplitude(input_string)

    # set each quantum number to its found value. If not found, denote it with a sum
    e = r"\Sigma\varepsilon" if amp_dict["e"] == "" else pm_dict[amp_dict["e"]]
    j = r"\Sigma J" if amp_dict["j"] == "" else amp_dict["j"]
    p = "" if amp_dict["p"] == "" else pm_dict[amp_dict["p"]]
    m = r"\Sigma m" if amp_dict["m"] == "" else m_proj_dict[amp_dict["m"]]
    l = r"\Sigma \ell" if amp_dict["l"] == "" else amp_dict["l"]

    return rf"${j}^{{{p}}}{l}_{{{m}}}^{{({e})}}$"


def human_format(num: float) -> str:
    """Converts orders of magnitude to letters i.e. 12369 -> 12.4K

    Args:
        num (float): positive floating point number

    Returns:
        str: first 3 significant digits with a character denoting its magnitude

    Raises:
        ValueError: If the number is negative
    """

    if num < 0:
        raise ValueError("Number must be positive")

    num = float(f"{num:.3g}")
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    orders = ["", "K", "M", "B", "T"]
    return f"{str(num).rstrip('0').rstrip('.')}{orders[magnitude]}"


def char_to_int(char: str) -> int:
    """Convert m-projection or L angular momenta character to integer value

    Args:
        char (str): single character. lowercase assumes m-projection and uppercase
            assumes it's an L value
    Returns:
        int: integer that char is representing
    """
    d = {
        "n": -2,
        "m": -1,
        "0": 0,
        "p": +1,
        "q": +2,
        "S": 0,
        "P": 1,
        "D": 2,
        "F": 3,
    }
    return d[char]
