"""Collection of tools useful for analyzing PWA fit results
"""

import itertools
import re

import numpy as np
import pandas as pd


class Plotter:
    # TODO: finish importing all the class methods
    def __init__(self, fit_file: str, data_file: str) -> None:
        """Initialize object with paths to csv files

        Args:
            fit_file (str): path to csv that holds all FitResults, made by fitsToCsv.C
            data_file (str): path to csv describing the raw data that is being fit to
        """
        df = pd.read_csv(fit_file, index_col="index")
        data_df = pd.read_csv(data_file)

        wrap_phases(df)
        # commonly used variables in all methods
        _bin_centers = data_df["bin_centers"]
        _fit_result = df["detected_events"]
        _fit_result_err = df["detected_events_err"]

        _coherent_sums = get_coherent_sums(df)
        phase_differences = get_phase_differences(df, _coherent_sums["eJPmL"])
        pass


def wrap_phases(
    df: pd.DataFrame = pd.DataFrame([]), series: pd.Series = pd.Series([], dtype=float)
) -> None | pd.Series:
    """Wrap phase differences to be from -pi to pi

    Two options of passing either a pandas dataframe, or series. The dataframe
    case handles avoiding editing any non phase difference columns. The series
    case is much simpler, and just applies the wrapping to each value

    Args:
        df (pd.DataFrame, optional): dataframe of fit results loaded from csv
            Defaults to pd.DataFrame([]).
        series (pd.Series, optional): pandas series of phase difference values.
            Defaults to pd.Series([], dtype=float).

    Returns:
        None | pd.Series: If df given, edits the df itself. Otherwise returns the phase
        wrapped Series
    """

    if not series.empty:
        if not df.empty:
            print("Only dataframe or series should be passed, NOT both. Exiting")
            return
        new_phases = []
        for phase in series:
            if phase < np.pi and phase > -np.pi:
                new_phases.append(phase)
                continue
            phase = (phase + np.pi) % (2 * np.pi) - np.pi
            new_phases.append(phase)

        return pd.Series(new_phases)

    if df.empty:
        print("Parameters are both empty. Exiting")
        return

    jp_list = []
    # find what JP's are present in dataframe
    for col in df.columns:
        if len(col) == 2:
            jp_list.append(col)

    for col in df.columns:
        # if column isn't two of the same or different jp's, skip it
        is_same_jp = False
        num_jps = 0

        for jp in jp_list:
            if col.count(jp) == 2:
                is_same_jp = True
            elif jp in col:
                num_jps += 1

        if not is_same_jp and num_jps != 2:
            continue

        # get phase, wrap if outside of [-pi,pi], and overwrite old phase
        for i in range(df.index.max() + 1):
            phase = df.iloc[i][col]
            if phase < np.pi and phase > -np.pi:
                continue
            phase = (phase + np.pi) % (2 * np.pi) - np.pi
            df.at[i, col] = phase

    return


def get_coherent_sums(df: pd.DataFrame) -> dict:
    """Returns a dict of coherent sums from a fit results dataframe

    TODO: Add handling of more than just eJPml type sums
        use regex flagging like done in convert_amp_name

    Args:
        df (pd.DataFrame): dataframe of fit results loaded from csv
    Returns:
        dict: Is of the form
        {Coherent sum : [amplitudes]} e.g. {"eJPmL", ["p1p0S", "m1mpP",...]}
    """
    # create an empty list for every type of coherent sum
    jp_values = []
    sum_types = ["eJPmL", "JPmL", "JPm", "JPL", "eJPm", "eJPL", "eJP", "JP"]
    coherent_sums = {d: [] for d in sum_types}

    # find what JP's are present in dataframe
    for column in df.columns:
        if len(column) == 2:
            jp_values.append(column)

    # grab all eJPml columns
    for column in df.columns:
        if not any([x in column for x in jp_values]):
            continue
        if "err" in column:
            continue
        if len(column) == 5:
            coherent_sums["eJPmL"].append(column)
        if len(column) == 4 and column[0] == "1":
            coherent_sums["JPmL"].append(column)

    return coherent_sums


def convert_amp_name(amp: str) -> str:
    """Converts amplitude type string to J^P L_m^(e) LaTeX style string

    Args:
        amp (str): amplitude string in eJPmL format
    Returns:
        str: Prettier LaTeX style amplitude in J^P L_m^(e) format. If phase difference
            is then J^P L_m^(e) - J^P L_m^(e)
    """

    pm_dict = {"m": "-", "p": "+", "": ""}
    m_proj_dict = {"m": -1, "0": 0, "p": +1, "": ""}

    # CASE 1: phase difference, always of form 'eJPmL_eJPmL'
    if "_" in amp:
        amps = amp.split("_")
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
    amp_dict = {
        "e": r"(?<![0-9]{1}[pm]{1})(?<!\d)[pm]",
        "jp": r"[0-9]{1}[pm]{1}",
        "m": r"([0-9]{1}[pm]{1})+[pm0]",  # assumes always form of 'JPm'
        "l": r"[SPDF]",
    }

    for quantum_number, expression in amp_dict.items():
        search = re.search(expression, amp)
        if search is not None:
            if quantum_number == "m":  # the search actually returns JPm
                amp_dict[quantum_number] = search.group()[-1]
            else:
                amp_dict[quantum_number] = search.group()
        else:
            amp_dict[quantum_number] = ""

    # set each quantum number to its found value. If not found, denote it with a sum
    e = r"\Sigma\varepsilon" if amp_dict["e"] == "" else pm_dict[amp_dict["e"]]
    j = r"\Sigma J" if amp_dict["jp"] == "" else amp_dict["jp"][0]
    p = "" if amp_dict["jp"] == "" else pm_dict[amp_dict["jp"][1]]
    m = r"\Sigma m" if amp_dict["m"] == "" else m_proj_dict[amp_dict["m"]]
    l = r"\Sigma \ell" if amp_dict["l"] == "" else amp_dict["l"]

    return rf"${j}^{{{p}}}{l}_{{{m}}}^{{({e})}}$"


def get_phase_differences(df: pd.DataFrame, eJPmL_cols: list) -> dict:
    """Returns dict of all the phase difference columns in the dataframe

    The keys are specifically like (eJPmL_1, eJPmL_2), and there is no way to
    know how those will be ordered in the dataframe a priori. We avoid this by
    using sorted tuple keys, so that when accessing the values of those keys
    there is no concern about the ordering of the tuple so long as its accessed
    like "dict.get(tuple(sorted(key)))".

    This function could be sped up by ignoring phase differences from separate
    reflectivities, but its kept generalized for potential future cases of using
    intensities that mix reflectivities

    TODO: Implement the needed regex matching to find the eJPmL columns within this func

    Args:
        df (pd.DataFrame): dataframe of fit results loaded from csv
        eJPmL_cols (list): all columns of coherent sum type eJPmL
    Returns:
        dict: dictionary where value is the phase difference column name in the df, and
        corresponding key is the sorted list of the 2 amplitudes in the phase difference
    """
    phase_differences = {}

    # get all possible combinations of eJPmL columns, and remove those that
    #   don't exist in the dataframe
    all_combos = list(itertools.combinations(eJPmL_cols, 2))
    phase_dif_columns = df.columns.values.tolist()
    for combo in all_combos:
        name = "_".join(combo)
        reverse_name = "_".join(reversed(combo))
        if name in phase_dif_columns:
            phase_differences[tuple(sorted(combo))] = name
        if reverse_name in phase_dif_columns:
            phase_differences[tuple(sorted(combo))] = reverse_name

    return phase_differences
