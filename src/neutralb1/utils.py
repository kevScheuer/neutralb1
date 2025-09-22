"""A collection of various helper functions for partial wave analysis"""

# TODO: Many of these functions can be combined into a pythonized version of the
#   AmplitudeParser class


import itertools
import os
import pathlib
import re
import sys
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from IPython.display import HTML
from IPython.display import Image as IPyImage
from IPython.display import display
from wand.image import Image as WandImage


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


def display_pdf(file: str, page: int, resolution: int = 100) -> None:
    """Display a PDF file in a Jupyter notebook using WandImage

    Converts a pdf to png for cropping of whitespace and proper display.

    Args:
        file (str): Path to the PDF file to be displayed.
        resolution (int, optional): Resolution for rendering the PDF. Defaults to 100.
    """
    with WandImage(filename=f"{file}[{page}]", resolution=resolution) as img:
        img.format = "png"
        img.trim()  # Crop whitespace
        png_bytes = img.make_blob()
        display(IPyImage(data=png_bytes))

    return


def big_print(input_string: str, size: float) -> None:
    """Print a string in large font size in a Jupyter notebook

    Args:
        input_string (str): The string to be displayed.
        size (float): The font size multiplier (e.g., 2.0 for double size).
    Returns:
        None, just displays the string
    """
    display(HTML(f"<span style='font-size:{size}em'>{input_string}</span>"))


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


def shared_ancestor_info(
    path1: str, path2: str
) -> Tuple[Optional[pathlib.Path], Optional[int], Optional[int]]:
    """Find the deepest shared ancestor of two paths and their respective depths.

    Args:
        path1 (str): First file or directory path.
        path2 (str): Second file or directory path.

    Returns:
        Tuple[Optional[pathlib.Path], Optional[int], Optional[int]]:
        - shared_dir: the deepest shared ancestor Path (or None if none)
        - depth1: how many parent steps from path1 to shared_dir
        - depth2: how many parent steps from path2 to shared_dir
    """
    p1 = pathlib.Path(path1).resolve()
    p2 = pathlib.Path(path2).resolve()
    p1_parents = [p1.parent] + list(p1.parents)
    p2_parents = [p2.parent] + list(p2.parents)
    for i1, ancestor1 in enumerate(p1_parents):
        if ancestor1 in p2_parents:
            i2 = p2_parents.index(ancestor1)
            return ancestor1, i1, i2
    return None, None, None


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
            split_key = list(key)
            if any(result_dict[char] == "" for char in split_key):
                continue
            coh_sum = "".join([result_dict[char] for char in split_key])
            coherent_sums[key].add(coh_sum)

    return {k: sorted(v) for k, v in coherent_sums.items()}  # sort each set


def get_phase_difference_dict(df: pd.DataFrame) -> Dict[tuple, str]:
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
        elif reverse_name in columns:
            phase_differences[combo] = reverse_name
            phase_differences[tuple(reversed(combo))] = reverse_name

    return phase_differences


def get_phase_differences(df: pd.DataFrame) -> Set[str]:
    """Returns the set of all phase difference columns in the dataframe

    Args:
        df (pd.DataFrame): dataframe of fit results loaded from csv
    Returns:
        set: set of all phase difference column names in the dataframe
    """
    return set(get_phase_difference_dict(df).values())


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
    j = r"\Sigma J" if amp_dict["J"] == "" else amp_dict["J"]
    p = "" if amp_dict["P"] == "" else pm_dict[amp_dict["P"]]
    if amp_dict["m"] == "":
        m = r"\Sigma m"
    else:
        m_val = m_proj_dict[amp_dict["m"]]
        # Always show sign for m (except for 0)
        if isinstance(m_val, int) and m_val != 0:
            m = f"{m_val:+d}"
        else:
            m = f"{m_val}"
    l = r"\Sigma \ell" if amp_dict["L"] == "" else amp_dict["L"]

    return rf"${j}^{{{p}}}{l}_{{{m}}}^{{({e})}}$"


def convert_moment_name(input_string: str) -> str:
    """Converts moment string to H^alpha(Jv, Lambda, J, M) style LaTeX string

    Args:
        input_string (str): string in H<alpha>_<Jv><Lambda><J><M> format
    Returns:
        str: Prettier LaTeX style string in H^alpha(Jv, Lambda, J, M) format

    Example:
        H0_0000 -> :math:`H^0(0, 0, 0, 0)`
        H1_1011 -> $H^1(1, 0, 1, 1)$
    """

    if not input_string.startswith("H") or "_" not in input_string:
        raise ValueError("Moment must be in 'H<alpha>_<Jv><Lambda><J><M>' format!")

    try:
        alpha = int(input_string[1])
        jv = int(input_string.split("_")[1][0])
        lam = int(input_string.split("_")[1][1])
        j = int(input_string.split("_")[1][2])
        m = int(input_string.split("_")[1][3])
    except (IndexError, ValueError):
        raise ValueError("Moment must be in 'H<alpha>_<Jv><Lambda><J><M>' format!")

    return rf"$H^{alpha}({jv}, {lam}, {j}, {m})$"


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


def propagate_product_error(
    var1: pd.Series,
    var1_err: pd.Series,
    var2: pd.Series,
    var2_err: pd.Series,
    covariance: Optional[pd.Series] = None,
) -> pd.Series:
    """Calculate error for division or multiplication of two series with errors

    Uses error propagation formula:
    σ(A/B) = (A/B) * sqrt((σA/A)² + (σB/B)²)

    Args:
        numerator (pd.Series): Numerator values.
        numerator_err (pd.Series): Uncertainties in numerator.
        denominator (pd.Series): Denominator values.
        denominator_err (pd.Series): Uncertainties in denominator.
        covariance (pd.Series, optional): Covariance between numerator and
            denominator. Defaults to None, meaning no covariance is considered.

    Returns:
        pd.Series: Propagated relative errors.

    TODO: replace with uncertainties package
    """
    # Avoid division by zero
    safe_var1 = np.where(var1 == 0, np.finfo(float).eps, var1)
    safe_var2 = np.where(var2 == 0, np.finfo(float).eps, var2)

    rel_cov = (
        2 * (covariance / (safe_var1 * safe_var2)) if covariance is not None else 0
    )

    relative_error = (var1 / var2) * np.sqrt(
        np.square(var1_err / safe_var1) + np.square(var2_err / safe_var2) - rel_cov
    )

    return relative_error


def circular_residual(
    angle1: float, angle2: float, in_degrees: bool = False, low=-np.pi, high=np.pi
) -> float:
    """Calculate the circular residual between two angles

    A residual in circular data needs to account for the periodicity of the data.
    This function calculates the smallest difference between the two angles.

    Args:
        angle1 (float): First angle
        angle2 (float): Second angle
        in_degrees (bool, optional): When True, indicates that the input angles are in
            degrees. The calculation and output will also be in degrees. Defaults to
            False, in which case the input angles are in radians and the output will be
            in radians.
        low (float, optional): Lower bound of the circular range. Defaults to -π.
        high (float, optional): Upper bound of the circular range. Defaults to π.

    Returns:
        float: Residual between the two angles, in the same units as the input
    """
    if in_degrees:
        angle1 = np.deg2rad(angle1)
        angle2 = np.deg2rad(angle2)
        low = np.deg2rad(low) if low != -np.pi else -np.pi
        high = np.deg2rad(high) if high != np.pi else np.pi

    diff = angle1 - angle2
    period = high - low
    # Wrap to [low, high)
    wrapped = (diff - low) % period + low

    if in_degrees:
        wrapped = np.rad2deg(wrapped)
    return wrapped
