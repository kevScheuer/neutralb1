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

        _bin_centers = data_df["mean"]
        _bin_centers_err = data_df["rms"]
        _fit_result = df["detected_events"]
        _fit_result_err = df["detected_events_err"]

        _coherent_sums = get_coherent_sums(df)
        _phase_differences = get_phase_differences(df)
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
            raise ValueError(
                "Only dataframe or series should be passed, NOT both. Exiting"
            )
        new_phases = []
        for phase in series:
            if phase < np.pi and phase > -np.pi:
                new_phases.append(phase)
                continue
            phase = (phase + np.pi) % (2 * np.pi) - np.pi
            new_phases.append(phase)

        return pd.Series(new_phases)

    if df.empty:
        raise ValueError("Parameters are both empty. Exiting")

    jp_list = get_coherent_sums(df)["JP"]

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


def parse_amplitude(amp: str) -> dict:
    """parse 'eJPmL' style amplitude string into its individual quantum numbers

    Args:
        amp (str): string like eJPmL, JPmL, JPm, JPL, eJPm, eJPL, eJP, JP

    Returns:
        dict: keys = quantum numbers (e, J, P, m, l). values = found values from string,
        or "" if not found
    """
    re_dict = {
        "e": r"(?<![0-9]{1}[pm]{1})(?<!\d)[pm]",
        "jp": r"[0-9]{1}[pm]{1}",  # j and p always appear together
        "m": r"([0-9]{1}[pm]{1})+[pm0]",  # assumes always form of 'JPm'
        "l": r"[SPDF]",
    }

    result_dict = {"e": "", "j": "", "p": "", "m": "", "l": ""}

    for quantum_number, expression in re_dict.items():
        search = re.search(expression, amp)
        if search is not None:
            # the search actually returns "JPm", so grab the last char if found
            if quantum_number == "m":
                result_dict["m"] = search.group()[-1]
            elif quantum_number == "jp":
                result_dict["j"] = search.group()[0]
                result_dict["p"] = search.group()[1]
            else:
                result_dict[quantum_number] = search.group()

    return result_dict


def get_coherent_sums(df: pd.DataFrame) -> dict:
    """Returns a dict of coherent sums from a fit results dataframe

    Args:
        df (pd.DataFrame): dataframe of fit results loaded from csv
    Returns:
        dict: Is of the form
        {Coherent sum : set(amplitudes)} e.g. {"eJPmL", ["p1p0S", "m1mpP",...]}
    """
    # create an empty list for every type of coherent sum
    sum_types = ["eJPmL", "JPmL", "JPm", "JPL", "eJPm", "eJPL", "eJP", "JP"]
    coherent_sums = {d: set() for d in sum_types}

    # grab all eJPml columns
    for column in df.columns:
        # skip the phase difference columns and any status columns
        if "_" in column or len(column) > 5:
            continue

        res = parse_amplitude(column)
        # only add to key if all elements of key are in the column
        for key in coherent_sums.keys():
            split_key = list(key.lower())
            if True in [res[char] == "" for char in split_key]:
                continue
            coh_sum = "".join([res[char] for char in split_key])
            coherent_sums[key].add(coh_sum)

    return coherent_sums


def convert_amp_name(amp: str) -> str:
    """Converts amplitude type string to J^P L_m^(e) LaTeX style string

    Args:
        amp (str): amplitude string in eJPmL format
    Returns:
        str: Prettier LaTeX style amplitude in J^P L_m^(e) format. If it's a phase
            difference it's then J^P L_m^(e) - J^P L_m^(e)
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
    amp_dict = parse_amplitude(amp)

    # set each quantum number to its found value. If not found, denote it with a sum
    e = r"\Sigma\varepsilon" if amp_dict["e"] == "" else pm_dict[amp_dict["e"]]
    j = r"\Sigma J" if amp_dict["j"] == "" else amp_dict["j"]
    p = "" if amp_dict["p"] == "" else pm_dict[amp_dict["p"]]
    m = r"\Sigma m" if amp_dict["m"] == "" else m_proj_dict[amp_dict["m"]]
    l = r"\Sigma \ell" if amp_dict["l"] == "" else amp_dict["l"]

    return rf"${j}^{{{p}}}{l}_{{{m}}}^{{({e})}}$"


def get_phase_differences(df: pd.DataFrame) -> dict:
    """Returns dict of all the phase difference columns in the dataframe

    The keys are specifically like (eJPmL_1, eJPmL_2), and there is no way to
    know how those will be ordered in the dataframe a priori. We avoid this by
    using sorted tuple keys, so that when accessing the values of those keys
    there is no concern about the ordering of the tuple so long as its accessed
    like "dict.get(tuple(sorted(key)))".

    This function could be sped up by ignoring phase differences from separate
    reflectivities, but its kept generalized for potential future cases of using
    intensities that mix reflectivities

    Args:
        df (pd.DataFrame): dataframe of fit results loaded from csv
    Returns:
        dict: dictionary where value is the phase difference column name in the df, and
        corresponding key is the sorted list of the 2 amplitudes in the phase difference
    """
    phase_differences = {}

    # get all possible combinations of eJPmL columns, and remove those that
    #   don't exist in the dataframe
    all_combos = list(itertools.combinations(get_coherent_sums(df)["eJPmL"], 2))
    columns = df.columns.values.tolist()
    for combo in all_combos:
        name = "_".join(combo)
        reverse_name = "_".join(reversed(combo))
        if name in columns:
            phase_differences[tuple(sorted(combo))] = name
        if reverse_name in columns:
            phase_differences[tuple(sorted(combo))] = reverse_name

    return phase_differences
