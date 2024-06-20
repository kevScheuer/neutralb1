"""Collection of tools useful for analyzing PWA fit results

NOTE: Remember to use red=pos refl, blue=neg refl standard now
"""

import itertools
import re

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Plotter:
    # TODO: finish importing all the class methods
    # TODO: add kwargs in init for figsize, dpi, and fontsizes (legend, axes, etc)
    def __init__(self, fit_file: str, data_file: str) -> None:
        """Initialize object with paths to csv files

        Args:
            fit_file (str): path to csv that holds all FitResults, made by fitsToCsv.C
            data_file (str): path to csv describing the raw data that is being fit to
        """
        self.df = pd.read_csv(fit_file, index_col="index")
        self.data_df = pd.read_csv(data_file)

        wrap_phases(self.df)

        self._mass_bins = self.data_df["mean"]
        self._mass_bins_err = self.data_df["rms"]
        self._bin_width = (self.data_df["high_edge"] - self.data_df["low_edge"])[0]

        self._coherent_sums = get_coherent_sums(self.df)
        self._phase_differences = get_phase_differences(self.df)
        pass

    def jp(self) -> None:
        """Plot each JP contribution to the total intensity as a function of mass"""

        # a map ensures consistency no matter the # of jps
        colors = matplotlib.colormaps["Dark2"].colors
        jp_map = {
            "0m": {"color": colors[0], "marker": "."},
            "1p": {"color": colors[1], "marker": "o"},
            "1m": {"color": colors[2], "marker": "s"},
            "2p": {"color": colors[3], "marker": "p"},
            "2m": {"color": colors[4], "marker": "h"},
            "3p": {"color": colors[5], "marker": "x"},
            "3m": {"color": colors[6], "marker": "d"},
        }

        fig, ax = plt.subplots()
        # plot data
        ax.errorbar(
            self._mass_bins,
            self.data_df["bin_contents"],
            self.data_df["bin_error"],
            self._mass_bins_err,
            "k.",
            label="Total (GlueX Phase-I)",
            markersize=4,
        )

        # Plot Fit Result
        ax.bar(
            self._mass_bins,
            self.df["detected_events"],
            width=self._bin_width,
            color="0.1",
            alpha=0.15,
            label="Fit Result",
        )
        ax.errorbar(
            self._mass_bins,
            self.df["detected_events"],
            self.df["detected_events_err"],
            fmt=",",
            color="0.1",
            alpha=0.2,
            markersize=0,
        )
        # plot each jp contribution
        for jp in self._coherent_sums["JP"]:
            ax.errorbar(
                self._mass_bins,
                self.df[jp],
                self.df[f"{jp}_err"],
                self._bin_width / 2,
                label=convert_amp_name(jp),
                markersize=4,
                linestyle="",
                **jp_map[jp],
            )

        ax.set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right", fontsize=18)
        ax.set_ylabel(f"Events / {self._bin_width:.3f} GeV", loc="top", fontsize=18)
        ax.set_ylim(bottom=0.0)

        ax.legend(fontsize=10)

        ax.grid(True, alpha=0.8)

        plt.show()
        pass

    def intensities(self) -> None:

        char_to_int = {"m": -1, "0": 0, "p": +1, "S": 0, "P": 1, "D": 2, "F": 3}
        int_to_char = {-1: "m", 0: "0", +1: "p"}
        pm_dict = {"m": "-", "p": "+"}

        # need to sort on the integer versions of the m-projections
        m_ints = sorted({char_to_int[JPm[-1]] for JPm in self._coherent_sums["JPm"]})

        fig, axs = plt.subplots(
            len(m_ints),
            len(self._coherent_sums["JPL"]),
            sharex=True,
            sharey=True,
            figsize=(15, 10),
            dpi=100,
        )

        # iterate through JPL (sorted like S, P, D, F wave) and sorted m-projections
        for row, jpl in enumerate(
            sorted(self._coherent_sums["JPL"], key=lambda JPL: char_to_int[JPL[-1]])
        ):
            for col, m in enumerate(m_ints):
                JPmL = f"{jpl[0:2]}{int_to_char[m]}{jpl[-1]}"

                if row == 0:
                    axs[row, col].set_title(f"m={char_to_int[JPmL[-2]]}", fontsize=16)
                if col == 0:
                    axs[row, col].set_ylabel(
                        rf"${JPmL[0]}^{{{pm_dict[JPmL[1]]}}}{JPmL[-1]}$",
                        fontsize=16,
                    )

                if "m" + JPmL in self._coherent_sums["eJPmL"]:
                    neg_plot = axs[row, col].errorbar(
                        self._mass_bins,
                        self.df["m" + JPmL],
                        self.df[f"m{JPmL}_err"],
                        self._bin_width / 2,
                        "o",
                        color="blue",
                        markersize=2,
                        label=r"$\epsilon=-1$",
                    )
                if "p" + JPmL in self._coherent_sums["eJPmL"]:
                    pos_plot = axs[row, col].errorbar(
                        self._mass_bins,
                        self.df["p" + JPmL],
                        self.df[f"p{JPmL}_err"],
                        self._bin_width / 2,
                        "o",
                        color="red",
                        markersize=2,
                        label=r"$\epsilon=+1$",
                    )

        # plot grids
        for ax in axs.reshape(-1):
            ax.grid(True, alpha=0.8)
            ax.set_ylim(bottom=0)

            # figure cosmetics
        fig.text(0.5, 0.04, r"$\omega\pi^0$ inv. mass (GeV)", ha="center", fontsize=18)
        fig.text(
            0.04,
            0.5,
            f"Events / {self._bin_width:.3f} GeV",
            ha="center",
            fontsize=18,
            rotation="vertical",
        )
        fig.legend(handles=[pos_plot, neg_plot], fontsize=18, loc="upper right")
        plt.show()
        pass

    def phase(self, amp1: str, amp2: str) -> None:
        """Plots the phase difference as a function of mass

        Order of amplitudes does not matter

        Args:
            amp1 (str): amplitude string in eJPmL format
            amp2 (str): amplitude string in eJPmL format

        Raises:
            ValueError: amp1 is not in the dataframe
            ValueError: amp2 is not in the dataframe
        """

        # first check that the amplitudes actually exist in the dataframe
        if amp1 not in self._coherent_sums["eJPmL"]:
            raise ValueError(f"Amplitude {amp1} not found in dataset")
        if amp2 not in self._coherent_sums["eJPmL"]:
            raise ValueError(f"Amplitude {amp2} not found in dataset")

        phase_dif = self._phase_differences[(amp1, amp2)]
        color = "red" if amp1[0] == "p" else "blue"

        fig, ax = plt.subplots()
        ax.errorbar(
            self._mass_bins,
            self.df[phase_dif].apply(np.rad2deg),
            self.df[phase_dif + "_err"].abs().apply(np.rad2deg),
            self._bin_width / 2,
            label=convert_amp_name(phase_dif),
            linestyle="",
            marker=".",
            color=color,
        )
        # plot the negative as well (natural sign ambiguity in the model)
        ax.errorbar(
            self._mass_bins,
            -self.df[phase_dif].apply(np.rad2deg),
            self.df[phase_dif + "_err"].abs().apply(np.rad2deg),
            self._bin_width / 2,
            linestyle="",
            marker=".",
            color=color,
        )

        ax.legend(fontsize=10)

        ax.set_yticks(np.linspace(-180, 180, 5))  # force to be in pi/2 intervals
        ax.set_ylim([-180, 180])
        ax.set_ylabel(r"Phase Diff. ($^{\circ}$)", loc="top", fontsize=18)
        ax.set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right", fontsize=18)
        ax.grid(True, alpha=0.8)

        plt.show()
        pass

    # TODO: add plot intensities method, but think about Jo's comment about the total
    #   data looking weird with it. Maybe remove data points all together?


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

    return {k: sorted(v) for k, v in coherent_sums.items()}  # sort each set


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
    creating keys for either ordering and setting both of their values to the order
    found in the dataframe.

    This function could be sped up by ignoring phase differences from separate
    reflectivities, but its kept generalized for potential future cases of using
    intensities that mix reflectivities

    Args:
        df (pd.DataFrame): dataframe of fit results loaded from csv
    Returns:
        dict: key = tuple of both possible combinations of each amplitude, val = the
            phase difference combination found in the dataframe
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
            phase_differences[combo] = name
            phase_differences[tuple(reversed(combo))] = name
        if reverse_name in columns:
            phase_differences[combo] = reverse_name
            phase_differences[tuple(reversed(combo))] = reverse_name

    return phase_differences
