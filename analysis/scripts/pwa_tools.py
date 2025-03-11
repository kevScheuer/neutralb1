"""Collection of tools useful for analyzing PWA fit results

Primary class is the Plotter, which ideally handles most standard plots of interest when
analyzing fit results
"""

import cmath
import itertools
import re
import warnings
from typing import Dict, List

import joypy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from cycler import cycler


class Plotter:
    def __init__(
        self,
        df: pd.DataFrame,
        data_df: pd.DataFrame,
        bootstrap_df: pd.DataFrame = None,
        truth_df: pd.DataFrame = None,
    ) -> None:
        """Initialize object with pandas dataframes

        TODO: The truth phase differences between generated and non-generated waves
            is non-zero, but rather should be a flat line at the generated waves phase
            value. This will need to be calculated "manually" using the re/im parts.
            This can be handled in reset_truth_phases by checking if the phase is flat
        Args:
            df (pd.DataFrame): FitResults from AmpTools
            data_df (pd.DataFrame): raw data points that AmpTools is fitting to
            bootstrap_df (pd.DataFrame, optional): all bootstrapped fit results,
                labelled by bin. Primary use is for plotting the distribution of
                each parameter and using its width to estimate the error of the nominal
                fit result in df.
            truth_df (pd.DataFrame, optional): only used when df is a from fit to
                generated signal MC. Contains the "true" parameters that were used to
                generate the MC sample. When non-empty, many plots will display a
                solid line for the truth information
        """
        self.df = df
        self.data_df = data_df
        self.bootstrap_df = bootstrap_df
        self.truth_df = truth_df

        plt.style.use(
            "/w/halld-scshelf2101/kscheuer/neutralb1/analysis/scripts/"
            "pwa_plotter.mplstyle"
        )  # load in default matplotlib style

        # Error handling
        if self.df.empty or self.data_df.empty:
            raise ValueError("dataframes are empty.")
        if self.df.shape[0] != self.data_df.shape[0]:
            raise ValueError(
                "df and data_df must have same number of rows!\n"
                f"df: {self.df.shape[0]},"
                f" data_df: {self.data_df.shape[0]}"
            )

        if any(status != 3 for status in self.df["eMatrixStatus"]):
            warnings.warn(
                (
                    "the following indices contain fit results whose"
                    " covariance matrix is not full and accurate."
                )
            )
            print(self.df[self.df["eMatrixStatus"] != 3]["eMatrixStatus"])

        # public attributes
        self.coherent_sums = get_coherent_sums(self.df)
        self.phase_differences = get_phase_differences(self.df)

        # private attributes
        self._mass_bins = self.data_df["m_center"]
        self._bin_width = (self.data_df["m_high"] - self.data_df["m_low"])[0]

        # --DATA PREPARATION--

        # wrap phases from -pi to pi range
        wrap_phases(self.df)
        if self.bootstrap_df is not None:
            wrap_phases(self.bootstrap_df)

        # Truth df needs special handling to allow for 1-1 comparison with fit results
        self.is_truth = True if self.truth_df is not None else False
        if self.is_truth:
            if self.df.shape[0] != self.truth_df.shape[0]:
                raise ValueError(
                    "Fit and truth dataframes must have same number of rows!\n"
                    f"df: {self.df.shape[0]},"
                    f" truth_df: {self.truth_df.shape[0]}"
                )
            # phases in truth fits are wrong, so must be "manually" set to the correct
            # values. Then we can wrap them
            self._reset_truth_phases(self.truth_df, self._mass_bins)
            wrap_phases(self.truth_df)

            # set non-existent columns to zero. Needed since the fit dataframe may use
            # a different waveset from the truth fit
            cols_to_add = []
            for sum_type in self.coherent_sums:
                for sum in self.coherent_sums[sum_type]:
                    if sum not in self.truth_df.columns:
                        cols_to_add.append(
                            pd.DataFrame(
                                np.zeros_like(self.df[sum]),
                                index=self.truth_df.index,
                                columns=[sum],
                            )
                        )
                        cols_to_add.append(
                            pd.DataFrame(
                                np.zeros_like(self.df[sum]),
                                index=self.truth_df.index,
                                columns=[f"{sum}_err"],
                            )
                        )
            for phase_dif in set(self.phase_differences.values()):
                # check if the reverse ordering of the phase difference is used, and
                # rename it to match the fit dataframe's ordering if so
                reverse_phase_dif = "".join(phase_dif.partition("_")[::-1])
                if reverse_phase_dif in self.truth_df.columns:
                    self.truth_df.rename(
                        columns={
                            reverse_phase_dif: phase_dif,
                            f"{reverse_phase_dif}_err": f"{phase_dif}_err",
                        },
                        inplace=True,
                    )
                # if phase difference not found, add it as columns of 0's
                if phase_dif not in self.truth_df.columns:
                    cols_to_add.append(
                        pd.DataFrame(
                            np.zeros_like(self.df[phase_dif]),
                            index=self.truth_df.index,
                            columns=[phase_dif],
                        )
                    )
                    cols_to_add.append(
                        pd.DataFrame(
                            np.zeros_like(self.df[phase_dif]),
                            index=self.truth_df.index,
                            columns=[f"{phase_dif}_err"],
                        )
                    )

            concat_list = [self.truth_df]
            concat_list.extend(cols_to_add)
            self.truth_df = pd.concat(concat_list, axis=1)
        pass

    def jp(self, data_label: str = "Total (GlueX Phase-I)") -> None:
        """Plot each JP contribution to the total intensity as a function of mass"""

        # a map ensures consistency no matter the # of jps
        colors = matplotlib.colormaps["Dark2"].colors
        jp_map = {
            "Bkgd": {"color": colors[0], "marker": "."},
            "0m": {"color": colors[1], "marker": "1"},
            "1p": {"color": colors[2], "marker": "o"},
            "1m": {"color": colors[3], "marker": "s"},
            "2p": {"color": colors[4], "marker": "p"},
            "2m": {"color": colors[5], "marker": "h"},
            "3p": {"color": colors[6], "marker": "x"},
            "3m": {"color": colors[7], "marker": "d"},
        }
        for d in jp_map.values():
            d.update({"markersize": 6})

        fig, ax = plt.subplots()
        # plot data
        ax.errorbar(
            self._mass_bins,
            self.data_df["events"],
            self.data_df["events_err"],
            self._bin_width / 2,
            fmt="k.",
            label=data_label,
        )

        # Plot Fit Result as a grey histogram with its error bars
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
        for jp, props in jp_map.items():
            if jp in self.df.columns:
                l = convert_amp_name(jp) if jp != "Bkgd" else "Iso. Background"
                ax.errorbar(
                    self._mass_bins,
                    self.df[jp],
                    self.df[f"{jp}_err"],
                    self._bin_width / 2,
                    label=l,
                    linestyle="",
                    **props,
                )

        # plot jp contributions for truth df
        if self.is_truth:
            for jp, props in jp_map.items():
                if jp in self.truth_df.columns:
                    ax.plot(
                        self._mass_bins,
                        self.truth_df[jp],
                        linestyle="-",
                        marker="",
                        color=props["color"],
                    )

        ax.set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right")
        ax.set_ylabel(f"Events / {self._bin_width:.3f} GeV", loc="top")
        ax.set_ylim(bottom=0.0)
        ax.legend()
        plt.minorticks_on()
        plt.show()
        pass

    def intensities(self, is_fit_fraction: bool = False, sharey: bool = False) -> None:
        """Plot all the amplitude intensities in a grid format

        Since the matrix plot is generally too small to see bin to bin features, this
        method plots every amplitude's intensity contribution in a grid format.
        Columns = m-projections, rows = JPL combinations. Reflectivities are plotted
        together on each subplot

        Args:
            is_fit_fraction (bool, optional): Scales all values by dividing them by the
                total intensity in each bin. Defaults to False.
            sharey (bool, optional): Sets each row to share a y-axis scale. Defaults to
                False
        """

        int_to_char = {-1: "m", 0: "0", +1: "p"}
        pm_dict = {"m": "-", "p": "+"}

        # sort the JPL values by the order of S, P, D, F waves, and the m-projections
        jpl_values = sorted(
            self.coherent_sums["JPL"], key=lambda JPL: char_to_int(JPL[-1])
        )
        m_ints = sorted({char_to_int(JPmL[-2]) for JPmL in self.coherent_sums["JPmL"]})

        fig, axs = plt.subplots(
            len(jpl_values),
            len(m_ints),
            sharex=True,
            sharey=sharey,
            figsize=(15, 10),
        )

        total = self.df["detected_events"] if is_fit_fraction else 1
        total_err = self.df["detected_events_err"] if is_fit_fraction else 0

        if self.is_truth:
            truth_total = self.truth_df["detected_events"] if is_fit_fraction else 1

        # iterate through JPL (sorted like S, P, D, F wave) and sorted m-projections
        for row, jpl in enumerate(jpl_values):
            for col, m in enumerate(m_ints):
                JPmL = f"{jpl[0:2]}{int_to_char[m]}{jpl[-1]}"

                # force sci notation so large ticklabels don't overlap with neighboring
                # plots if not plotting fit fraction
                if not is_fit_fraction:
                    axs[row, col].ticklabel_format(
                        axis="y", style="sci", scilimits=(0, 0)
                    )

                # set the title and ylabel for the first row and column
                if row == 0:
                    axs[row, col].set_title(f"m={char_to_int(JPmL[-2])}", fontsize=18)
                if col == 0:
                    axs[row, col].set_ylabel(
                        rf"${JPmL[0]}^{{{pm_dict[JPmL[1]]}}}{JPmL[-1]}$", fontsize=18
                    )

                # plot the negative reflectivity contribution
                if "m" + JPmL in self.coherent_sums["eJPmL"]:
                    neg_refl = self.df["m" + JPmL] / total
                    neg_refl_err = neg_refl * np.sqrt(
                        np.square(self.df[f"m{JPmL}_err"] / self.df["m" + JPmL])
                        + np.square(total_err / total)
                    )

                    neg_plot = axs[row, col].errorbar(
                        self._mass_bins,
                        neg_refl,
                        neg_refl_err,
                        self._bin_width / 2,
                        marker="v",
                        linestyle="",
                        markersize=8,
                        color="blue",
                        alpha=0.5,
                        label=r"$\varepsilon=-1$",
                    )
                    if self.is_truth:
                        axs[row, col].plot(
                            self._mass_bins,
                            self.truth_df["m" + JPmL] / truth_total,
                            linestyle="-",
                            marker="",
                            color="blue",
                        )
                # plot the positive reflectivity contribution
                if "p" + JPmL in self.coherent_sums["eJPmL"]:
                    pos_refl = self.df["p" + JPmL] / total
                    pos_refl_err = pos_refl * np.sqrt(
                        np.square(self.df[f"p{JPmL}_err"] / self.df["p" + JPmL])
                        + np.square(total_err / total)
                    )

                    pos_plot = axs[row, col].errorbar(
                        self._mass_bins,
                        pos_refl,
                        pos_refl_err,
                        self._bin_width / 2,
                        marker="^",
                        linestyle="",
                        markersize=8,
                        color="red",
                        alpha=0.5,
                        label=r"$\varepsilon=+1$",
                    )
                    if self.is_truth:
                        axs[row, col].plot(
                            self._mass_bins,
                            self.truth_df["p" + JPmL] / truth_total,
                            linestyle="-",
                            marker="",
                            color="red",
                        )

        # reset limits to 0 for all plots
        for ax in axs.reshape(-1):
            ax.set_ylim(bottom=0)

        # figure cosmetics
        fig.text(0.5, 0.04, r"$\omega\pi^0$ inv. mass (GeV)", ha="center", fontsize=20)
        fig.text(
            0.04,
            0.5,
            f"Events / {self._bin_width:.3f} GeV",
            ha="center",
            va="center",
            rotation="vertical",
            rotation_mode="anchor",
            fontsize=20,
        )
        fig.legend(handles=[pos_plot, neg_plot], loc="upper right")
        plt.show()
        pass

    def phase(self, amp1: str, amp2: str) -> None:
        """Plots the phase difference as a function of mass

        Order of amplitudes does not matter

        Args:
            amp1 (str): amplitude string in eJPmL format
            amp2 (str): amplitude string in eJPmL format

        Raises:
            ValueError: amp1 or amp2 is not in the dataframe
            ValueError: amplitudes are not from the same reflectivity
        """

        # first check that the amplitudes actually exist in the dataframe
        if amp1 not in self.coherent_sums["eJPmL"]:
            raise ValueError(f"Amplitude {amp1} not found in dataset")
        if amp2 not in self.coherent_sums["eJPmL"]:
            raise ValueError(f"Amplitude {amp2} not found in dataset")
        if amp1[0] != amp2[0]:
            raise ValueError(f"Amplitudes must be from same reflectivity")

        phase_dif = self.phase_differences[(amp1, amp2)]
        color = "red" if amp1[0] == "p" else "blue"

        fig, ax = plt.subplots()
        ax.errorbar(
            self._mass_bins,
            self.df[phase_dif],
            self.df[f"{phase_dif}_err"].abs(),
            self._bin_width / 2,
            label=convert_amp_name(phase_dif),
            linestyle="",
            marker=".",
            color=color,
        )
        # plot the negative as well (natural sign ambiguity in the model)
        ax.errorbar(
            self._mass_bins,
            -self.df[phase_dif],
            self.df[f"{phase_dif}_err"].abs(),
            self._bin_width / 2,
            linestyle="",
            marker=".",
            color=color,
        )

        if self.is_truth:
            ax.plot(
                self._mass_bins,
                self.truth_df[phase_dif],
                linestyle="-",
                marker="",
                color=color,
            )

        ax.legend()

        ax.set_yticks(np.linspace(-180, 180, 5))  # force to be in pi intervals
        ax.set_ylim([-180, 180])
        ax.set_ylabel(r"Phase Diff. ($^{\circ}$)", loc="top")
        ax.set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right")

        plt.show()
        pass

    def mass_phase(self, amp1: str, amp2: str, color: str = "black") -> None:
        """Plot the amplitude intensities with their phase difference together

        Args:
            amp1 (str): amplitude string in eJPmL format
            amp2 (str): amplitude string in eJPmL format
            color (str, optional): color to plot. Defaults to 'black'

        Raises:
            ValueError: amp1 or amp2 is not in the dataframe
            ValueError: amplitudes are not from the same reflectivity
        """

        # first check that the amplitudes actually exist in the dataframe
        if amp1 not in self.coherent_sums["eJPmL"]:
            raise ValueError(f"Amplitude {amp1} not found in dataset")
        if amp2 not in self.coherent_sums["eJPmL"]:
            raise ValueError(f"Amplitude {amp2} not found in dataset")
        if amp1[0] != amp2[0]:
            raise ValueError(f"Amplitudes must be from same reflectivity")

        fig, axs = plt.subplots(
            2,
            1,
            sharex=True,
            gridspec_kw={"wspace": 0.0, "hspace": 0.07},
            height_ratios=[3, 1],
        )

        # plot the two amplitudes
        axs[0].errorbar(
            self._mass_bins,
            self.df[amp1],
            self.df[f"{amp1}_err"],
            self._bin_width / 2,
            "o",
            color=color,
            label=convert_amp_name(amp1),
        )
        axs[0].errorbar(
            self._mass_bins,
            self.df[amp2],
            self.df[f"{amp2}_err"],
            self._bin_width / 2,
            "s",
            color=color,
            label=convert_amp_name(amp2),
        )

        if self.is_truth:
            axs[0].plot(
                self._mass_bins,
                self.truth_df[amp1],
                linestyle="-",
                marker="",
                color=color,
            )
            axs[0].plot(
                self._mass_bins,
                self.truth_df[amp2],
                linestyle="-",
                marker="",
                color=color,
            )

        # plot the phase difference
        phase_dif = self.phase_differences[(amp1, amp2)]
        axs[1].errorbar(
            self._mass_bins,
            self.df[phase_dif],
            self.df[f"{phase_dif}_err"].abs(),
            self._bin_width / 2,
            linestyle="",
            marker=".",
            color=color,
        )
        axs[1].errorbar(
            self._mass_bins,
            -self.df[phase_dif],
            self.df[f"{phase_dif}_err"].abs(),
            self._bin_width / 2,
            linestyle="",
            marker=".",
            color=color,
        )
        if self.is_truth:
            axs[1].plot(
                self._mass_bins,
                self.truth_df[phase_dif],
                linestyle="-",
                marker="",
                color=color,
            )

        axs[0].set_ylim(bottom=0.0)
        axs[0].set_ylabel(f"Events / {self._bin_width:.3f} GeV", loc="top")

        axs[1].set_yticks(np.linspace(-180, 180, 5))  # force to be in pi intervals
        axs[1].set_ylim([-180, 180])
        axs[1].set_ylabel(r"Phase Diff. ($^{\circ}$)", loc="center")
        axs[1].set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right")

        axs[0].legend(loc="upper right")

        plt.show()
        pass

    def matrix(self) -> None:
        """Scatter plot-like matrix of subplots, for intensities and phase differences

        The diagonal of the plot contains the intensity of each wave, in both
        reflectivities. The off diagonals are the phase difference between each wave,
        with positive reflectivity on the upper triangle and negative on the lower.
        Amplitude names are labelled on the first row and column.

        NOTE: This plot will become cramped for large models, may need dpi adjusted
        """

        fig, axs = plt.subplots(
            len(self.coherent_sums["JPmL"]),
            len(self.coherent_sums["JPmL"]),
            sharex=True,
            figsize=(13, 8),
            dpi=500,
        )
        # only the top left corner plot will have ticks for the data scale
        max_diag = max([self.df[x].max() for x in self.coherent_sums["eJPmL"]])
        data_y_ticks = np.linspace(0, max_diag, 4)
        axs[0, 0].set_yticks(data_y_ticks, [human_format(num) for num in data_y_ticks])
        for i, JPmL in enumerate(self.coherent_sums["JPmL"]):
            for j, JPmL_dif in enumerate(self.coherent_sums["JPmL"]):
                # change tick label sizes for all plots and the max
                axs[i, j].tick_params("both", labelsize=6)
                axs[i, j].set_ylim(top=max_diag)

                # write mass ticks for only the bottom row
                if i == len(self.coherent_sums["JPmL"]) - 1:
                    axs[i, j].xaxis.set_tick_params(labelbottom=True)

                # PLOT DIAGONALS
                if i == j:
                    # make outline of the plot thick to set them apart
                    spines = ["left", "right", "top", "bottom"]
                    for sp in spines:
                        axs[i, j].spines[sp].set_linewidth(2)
                        axs[i, j].spines[sp].set_color("black")

                    # top left corner plot needs both title and ylabel amplitude labels
                    if i == 0:
                        axs[i, j].set_title(
                            convert_amp_name(JPmL_dif).replace(
                                r"^{(\Sigma\varepsilon)}", ""
                            ),
                            fontsize=10,
                        )
                        axs[i, j].set_ylabel(
                            convert_amp_name(JPmL).replace(
                                r"^{(\Sigma\varepsilon)}", ""
                            ),
                            fontsize=10,
                        )
                    # otherwise remove the y ticks on the diagonal, not enough room
                    else:
                        axs[i, j].set_yticks(data_y_ticks, [""] * 4)

                    # plot the reflectivity contributions
                    neg_plot = axs[i, j].errorbar(
                        self._mass_bins,
                        self.df[f"m{JPmL}"],
                        self.df[f"m{JPmL}_err"],
                        self._bin_width / 2,
                        "s",
                        color="blue",
                        markersize=2,
                        label=r"$\varepsilon=-1$",
                    )
                    pos_plot = axs[i, j].errorbar(
                        self._mass_bins,
                        self.df[f"p{JPmL}"],
                        self.df[f"p{JPmL}_err"],
                        self._bin_width / 2,
                        "o",
                        color="red",
                        markersize=2,
                        label=r"$\varepsilon=+1$",
                    )
                    if self.is_truth:
                        axs[i, j].plot(
                            self._mass_bins,
                            self.truth_df[f"m{JPmL}"],
                            linestyle="-",
                            marker="",
                            color="blue",
                        )
                        axs[i, j].plot(
                            self._mass_bins,
                            self.truth_df[f"p{JPmL}"],
                            linestyle="-",
                            marker="",
                            color="red",
                        )

                # PLOT POSITIVE REFLECTIVITY PHASE DIFFERENCES (UPPER TRIANGLE)
                elif j > i:
                    # give same yticks as lower triangle, but with no label
                    axs[i, j].set_yticks(np.linspace(-180, 180, 5), [""] * 5)
                    axs[i, j].set_ylim([-180, 180])

                    # write amplitude titles for the top row
                    if i == 0:
                        axs[i, j].set_title(
                            convert_amp_name(JPmL_dif).replace(
                                r"^{(\Sigma\varepsilon)}", ""
                            ),
                            fontsize=10,
                        )

                    phase_dif = self.phase_differences[(f"p{JPmL}", f"p{JPmL_dif}")]

                    # plot phase and its sign ambiguity, with error bands
                    axs[i, j].plot(
                        self._mass_bins,
                        self.df[phase_dif].abs(),
                        linestyle="",
                        color="red",
                        marker=".",
                        ms=1,
                    )
                    axs[i, j].fill_between(
                        self._mass_bins,
                        self.df[phase_dif].abs() - self.df[f"{phase_dif}_err"].abs(),
                        self.df[phase_dif].abs() + self.df[f"{phase_dif}_err"].abs(),
                        color="red",
                        alpha=0.2,
                    )
                    axs[i, j].plot(
                        self._mass_bins,
                        -self.df[phase_dif].abs(),
                        linestyle="",
                        color="red",
                        marker=".",
                        ms=1,
                    )
                    axs[i, j].fill_between(
                        self._mass_bins,
                        -self.df[phase_dif].abs() - self.df[f"{phase_dif}_err"].abs(),
                        -self.df[phase_dif].abs() + self.df[f"{phase_dif}_err"].abs(),
                        color="red",
                        alpha=0.2,
                    )

                    if self.is_truth:
                        axs[i, j].plot(
                            self._mass_bins,
                            self.truth_df[phase_dif],
                            linestyle="-",
                            marker="",
                            color="red",
                        )

                # PLOT NEGATIVE REFLECTIVITY PHASE DIFFERENCES
                elif j < i:
                    axs[i, j].set_ylim([-180, 180])
                    # write amplitude titles and y-tick labels on 1st column
                    if j == 0:
                        axs[i, j].set_ylabel(
                            convert_amp_name(JPmL).replace(
                                r"^{(\Sigma\varepsilon)}", ""
                            ),
                            fontsize=10,
                        )
                        axs[i, j].set_yticks(np.linspace(-180, 180, 5))
                    else:
                        axs[i, j].set_yticks(np.linspace(-180, 180, 5), [""] * 5)

                    phase_dif = self.phase_differences[(f"m{JPmL}", f"m{JPmL_dif}")]

                    # plot phase and its sign ambiguity, with error bands
                    axs[i, j].plot(
                        self._mass_bins,
                        self.df[phase_dif].abs(),
                        linestyle="",
                        color="blue",
                        marker=".",
                        ms=1,
                    )
                    axs[i, j].fill_between(
                        self._mass_bins,
                        self.df[phase_dif].abs() - self.df[f"{phase_dif}_err"].abs(),
                        self.df[phase_dif].abs() + self.df[f"{phase_dif}_err"].abs(),
                        color="blue",
                        alpha=0.2,
                    )
                    axs[i, j].plot(
                        self._mass_bins,
                        -self.df[phase_dif].abs(),
                        linestyle="",
                        color="blue",
                        marker=".",
                        ms=1,
                    )
                    axs[i, j].fill_between(
                        self._mass_bins,
                        -self.df[phase_dif].abs() - self.df[f"{phase_dif}_err"].abs(),
                        -self.df[phase_dif].abs() + self.df[f"{phase_dif}_err"].abs(),
                        color="blue",
                        alpha=0.2,
                    )
                    if self.is_truth:
                        axs[i, j].plot(
                            self._mass_bins,
                            self.truth_df[phase_dif],
                            linestyle="-",
                            marker="",
                            color="blue",
                        )

        # Last edits to figure before printing
        fig.text(
            0.5, 0.06, r"$\omega\pi^0$ inv. mass $(GeV)$", ha="center", fontsize=15
        )
        fig.text(
            0.06,
            0.5,
            r"Phase Differences $(^{\circ})$",
            ha="center",
            fontsize=15,
            rotation="vertical",
        )
        fig.legend(handles=[pos_plot, neg_plot], fontsize=12, loc="upper right")

        plt.show()
        pass

    def ds_ratio(self) -> None:
        """Plot the ratio and phase between the D and S waves

        The plots chang depending on whether the model had a defined D/S ratio. If it is
        defined then the covariance of the ratio and phase parameters are plotted, but
        if not then its the correlation. Due to plotting style, plots have fixed sizes.

        TODO: Add a raise if D or S wave are not present, and the reflectivities
        """
        #
        if "dsratio" in self.df.columns:
            fig, axs = plt.subplots(
                3, 1, sharex=True, gridspec_kw={"wspace": 0.0, "hspace": 0.12}
            )
            ax_ratio, ax_phase, ax_corr = axs
            # RATIO PLOT
            ax_ratio.plot(  # plot the nominal ratio
                self._mass_bins,
                np.full_like(self._mass_bins, 0.27),
                "k-",
                label="E852 ratio (0.27)",
            )
            ax_ratio.errorbar(
                self._mass_bins,
                self.df["dsratio"],
                self.df["dsratio_err"],
                self._bin_width / 2,
                marker=".",
                linestyle="",
                color="gray",
            )
            # PHASE PLOT
            ax_phase.plot(  # plot the nominal phase (in degrees)
                self._mass_bins,
                np.full_like(self._mass_bins, 10.54),
                "k--",
                label=r"E852 phase $(10.54^\circ)$",
            )
            ax_phase.plot(self._mass_bins, np.full_like(self._mass_bins, -10.54), "k--")
            ax_phase.errorbar(
                self._mass_bins,
                self.df["dphase"],
                self.df["dphase_err"].abs(),
                self._bin_width / 2,
                marker=".",
                linestyle="",
                color="gray",
            )
            ax_phase.errorbar(
                self._mass_bins,
                -self.df["dphase"],
                self.df["dphase_err"].abs(),
                self._bin_width / 2,
                marker=".",
                linestyle="",
                color="gray",
            )
            # CORRELATION PLOT
            ax_corr.plot(
                self._mass_bins,
                (
                    self.df["cov_dsratio_dphase"]
                    / (self.df["dsratio_err"] * self.df["dphase_err"])
                ).fillna(0),
                marker="8",
                linestyle="",
                color="gray",
            )

            # COSMETICS
            ax_ratio.set_ylabel("D/S ratio", loc="top")
            ax_ratio.set_ylim(0, 1)

            ax_phase.set_ylabel(r"D-S (${{}}^\circ$)", loc="top")
            ax_phase.set_yticks(np.linspace(-180.0, 180.0, 5))
            ax_phase.set_ylim(-180.0, 180.0)

            ax_corr.set_xlabel(r"$\omega\pi^0$ mass inv. mass $(GeV)$", loc="right")
            ax_corr.set_ylabel(r"$\rho$ (D/S, D-S)", loc="top")
            ax_corr.set_ylim(-1, 1)

        else:
            marker_map = {"q": 6, "p": "^", "0": ".", "m": "v", "n": 7}
            color_map = {"p": "red", "m": "blue"}
            fig = plt.figure(figsize=(11, 4), dpi=300)
            subfigs = fig.subfigures(1, 2)

            ax_ratio, ax_phase = subfigs[0].subplots(2, 1, sharex=True)
            ax_corr = subfigs[1].subplots()

            # plot the nominal E852 values first
            ax_ratio.plot(self._mass_bins, np.full_like(self._mass_bins, 0.27), "k-")
            ax_phase.plot(self._mass_bins, np.full_like(self._mass_bins, 10.54), "k--")
            ax_phase.plot(self._mass_bins, np.full_like(self._mass_bins, -10.54), "k--")
            ax_corr.plot(
                np.full_like(np.linspace(-180.0, 180.0, 5), 0.27),
                np.linspace(-180.0, 180.0, 5),
                "k-",
                label="E852 ratio (0.27)",
            )
            ax_corr.plot(
                np.linspace(0, 1, 2),
                np.full_like(np.linspace(0, 1, 2), 10.54),
                "k--",
                label=r"E852 phase $(10.54^\circ)$",
            )
            for eJPmL in self.coherent_sums["eJPmL"]:
                e = eJPmL[0]
                m = eJPmL[-2]
                L = eJPmL[-1]
                if L != "D":
                    continue
                D_wave = eJPmL
                S_wave = eJPmL[:-1] + "S"

                # get ratio and phase difference
                ratio = self.df[D_wave].apply(np.sqrt) / self.df[S_wave].apply(np.sqrt)
                ratio_err = ratio * np.sqrt(
                    np.square(
                        self.df[f"{D_wave}_err"].apply(np.sqrt)
                        / self.df[D_wave].apply(np.sqrt)
                    )
                    + np.square(
                        self.df[f"{S_wave}_err"].apply(np.sqrt)
                        / self.df[S_wave].apply(np.sqrt)
                    )
                )
                phase = self.df[self.phase_differences[(D_wave, S_wave)]]
                phase_err = self.df[
                    f"{self.phase_differences[(D_wave, S_wave)]}_err"
                ].abs()

                label = convert_amp_name(eJPmL).replace("D", "(D/S)")

                ax_ratio.errorbar(
                    self._mass_bins,
                    ratio,
                    ratio_err,
                    self._bin_width / 2,
                    marker=marker_map[m],
                    color=color_map[e],
                    linestyle="",
                    markersize=6,
                )
                ax_phase.errorbar(
                    self._mass_bins,
                    phase,
                    phase_err,
                    self._bin_width / 2,
                    marker=marker_map[m],
                    color=color_map[e],
                    linestyle="",
                    markersize=6,
                )
                ax_phase.errorbar(
                    self._mass_bins,
                    -phase,
                    phase_err,
                    self._bin_width / 2,
                    marker=marker_map[m],
                    color=color_map[e],
                    linestyle="",
                    markersize=6,
                )
                ax_corr.errorbar(
                    ratio,
                    phase,
                    phase_err,
                    ratio_err,
                    marker=marker_map[m],
                    color=color_map[e],
                    linestyle="",
                    label=label,
                    markersize=6,
                )
                ax_corr.errorbar(
                    ratio,
                    -phase,
                    abs(phase_err),
                    ratio_err,
                    marker=marker_map[m],
                    color=color_map[e],
                    linestyle="",
                    label=label,
                    markersize=6,
                )
            ax_ratio.set_ylabel("D/S ratio", loc="top", fontsize=12)
            ax_ratio.set_yscale("log")

            ax_phase.set_xlabel(
                r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right", fontsize=12
            )
            ax_phase.set_ylabel(
                r"D-S phase difference $(^\circ)$", loc="top", fontsize=10
            )
            ax_phase.set_yticks(np.linspace(-180.0, 180.0, 5))
            ax_phase.set_ylim(-180.0, 180.0)

            ax_corr.set_xlabel("D/S ratio", loc="right", fontsize=12)
            ax_corr.set_xscale("log")
            ax_corr.set_ylabel(
                r"D-S phase difference $(^\circ)$", loc="top", fontsize=12
            )
            ax_corr.set_ylim(-180.0, 180)
            ax_corr.set_yticks(np.linspace(-180.0, 180.0, 5))

            subfigs[1].subplots_adjust(right=0.65)

            handles, labels = ax_corr.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            fig.legend(
                by_label.values(),
                by_label.keys(),
                loc="center right",
                bbox_to_anchor=(1, 0.5),
                fontsize=8,
            )
            plt.show()
        pass

    def mean_correlation(
        self,
        column_groups: List[List[str]] = None,
        labels: List[str] = None,
        method: str = "pearson",
    ) -> None:
        """Plots mean correlation magnitude of each column group as a function of mass

        A correlation matrix is constructed for each group of dataframe columns, and
        the mean correlation magnitude is the averaged sum of the absolute value of the
        correlation matrix. This mean value is then plotted as a function of mass.
        This essentially boils an entire corr matrix down to one parameter, defined
        from [0,1], that reports how correlated all the columns are with each other.

        Generally the more columns in each group, the less sensitive this parameter
        becomes to individual strong correlations.

        Args:
            column_groups (List[List[str]], optional): The mean correlation magnitude is
                calculated for all the columns in a sub-list (or group). This value for
                each group is then plotted together. Defaults to None, and instead uses
                all real and imaginary production parameters.
            labels (list, optional): Legend label for each group of columns. Defaults to
                None, and sets the label to be the real and imaginary production
                parameters.
            method (str, optional): Method to calculate the correlation matrix. Defaults
                to pearson, but can also be kendall or spearman.

        Raises:
            ValueError: if bootstrap df was not defined, since its optional upon init
            ValueError: need to have same number of labels as column groups
        """
        if self.bootstrap_df is None:
            raise ValueError("Bootstrap df was not defined on instantiation")

        # when no columns requested, build a list from real and imag production params
        if not column_groups:
            column_groups = [
                list(
                    itertools.chain.from_iterable(
                        (x + "_re", x + "_im") for x in self.coherent_sums["eJPmL"]
                    )
                )
            ]
            labels = [r"$[c^i]_m^{(\varepsilon)}$"]

        if len(column_groups) != len(labels):
            raise ValueError(
                "Number and position of labels must match those of each sub-list of"
                " columns"
            )

        fig, ax = plt.subplots()
        for group, label in zip(column_groups, labels):
            mcm = []
            # get mean correlation magnitude for each bin using the bootstrap df to
            # calculate the correlation matrix
            for i in self.df.index:
                # Take absolute value to make it a magnitude matrix
                abs_corr_matrix = (
                    self.bootstrap_df[self.bootstrap_df["bin"] == i][group]
                    .corr(method=method)
                    .apply(lambda x: abs(x))
                )
                mcm.append(
                    abs_corr_matrix.sum().sum()
                    / np.square(len(abs_corr_matrix.columns))
                )
            ax.plot(self._mass_bins, mcm, linestyle="-", label=label)

        # cosmetics
        ax.set_ylim(bottom=0.0, top=1.0)
        ax.set_ylabel("Mean Correlation Magnitude", loc="top")
        ax.set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right")

        ax.legend()
        plt.show()
        pass

    def correlation_matrix(
        self, columns: List[str], mass_bin: int | str, method: str = "pearson"
    ) -> None:
        """Plot the correlation matrix in a single mass bin

        The correlation matrix is calculated for the columns in the bootstrap dataframe
        for a single mass bin. The matrix is then plotted as a heatmap, with its values.

        Args:
            columns (List[str]): bootstrap df columns to calculate the correlation
                matrix for
            mass_bin (int | str): mass bin, either as integer (bin #) or string
                ("bin_low-bin_high")
            method (str, optional): Method of correlation. Defaults to "pearson", but
                can also be kendall or spearman.

        Raises:
            ValueError: if bootstrap df was not defined, since its optional upon init
            ValueError: requested columns not in the bootstrap dataframe
            ValueError: mass bin is not in dataframe, or is not string or integer
        """

        if self.bootstrap_df is None:
            raise ValueError("Bootstrap df was not defined on instantiation")

        # check if requested columns are in the dataframe
        missing_columns = [
            col for col in columns if col not in self.bootstrap_df.columns
        ]
        if missing_columns:
            raise ValueError(
                f"The following columns are not in bootstrap_df: {missing_columns}"
            )

        # calculate the correlation matrix
        if isinstance(mass_bin, str):
            corr_matrix = self.bootstrap_df[self.bootstrap_df["bin_range"] == mass_bin][
                columns
            ].corr(method=method)
        elif isinstance(mass_bin, int):
            corr_matrix = self.bootstrap_df[self.bootstrap_df["bin"] == mass_bin][
                columns
            ].corr(method=method)
        else:
            raise ValueError("Mass bin must be a string or integer")

        # check that the mass bin was found in the dataframe
        corr_matrix.dropna(inplace=True)
        if corr_matrix.empty:
            raise ValueError(f"Mass bin {mass_bin} not found in bootstrap df")

        # round to 2 decimal places
        corr_matrix = corr_matrix.round(2)

        # rename the columns to prettier LaTeX names
        renamed_titles = {}
        for col in columns:
            if any(
                col in sublist for sublist in self.coherent_sums.values()
            ) or col in set(self.phase_differences.values()):
                renamed_titles[col] = convert_amp_name(col)
        corr_matrix.rename(index=renamed_titles, columns=renamed_titles, inplace=True)

        # plot the matrix as a heatmap
        plt.figure(figsize=(10, 8))  # slightly longer figure to accommodate color bar
        sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", vmin=-1, vmax=1)

        # rotate titles for better viewing
        plt.yticks(rotation=0)
        plt.xticks(rotation=45, ha="right")

        plt.show()

        pass

    def pull_distribution(
        self,
        columns: str | List[str],
        y_limits: List[int] = None,
        scale: bool = False,
        scale_factor: int = 500,
        **legend_kwargs,
    ) -> None:
        """Plot the pull distribution of column(s)

        A pull distribution is the difference between the nominal fit value and the true
        value, divided by the error on the nominal fit value. This is then plotted as a
        function of mass. The truth information is required for this plot. The pull
        distribution is a good way to check if the nominal fit is systematically biased
        in one direction.

        Args:
            columns (str | List[str]): column(s) to plot the pull distribution for
            y_limits (List[int], optional): y-axis limits. Defaults to None, and is
                automatically set to the max pull value.
            scale (bool, optional): scale the marker size by the nominal fit values of
                the column. Works well for intensities to show relative sizes. Defaults
                to False.
            scale_factor (int, optional): factor to scale the marker size by if scale is
                True. Often needs tuning to get the right size. Defaults to 500.
            **legend_kwargs: additional keyword arguments to pass to the legend.
        Raises:
            ValueError: Truth information is required for pull distribution
        """
        if not self.is_truth:
            raise ValueError("Truth information is required for pull distribution")

        if isinstance(columns, str):
            columns = [columns]

        # Create a color cycle using the Dark2 colormap
        color_cycle = cycler(color=matplotlib.cm.Accent.colors)
        fig, ax = plt.subplots()
        ax.set_prop_cycle(color_cycle)

        if scale:
            # get normalized max values in the columns to determine marker size
            norm = max([self.df[col].max() for col in columns])

        max_pull = 0
        for column in columns:
            # calculate pull and store the max value
            pull = (self.df[column] - self.truth_df[column]) / self.df[f"{column}_err"]
            max_pull = max(max_pull, pull.abs().max())

            marker_sizes = scale_factor * (self.df[column] / norm) if scale else 1

            ax.scatter(
                self._mass_bins,
                pull,
                linestyle="",
                s=marker_sizes,
                label=convert_amp_name(column),
            )

        # central line at 0
        ax.axhline(y=0, linestyle="-", color="black")

        ax.set_ylim(y_limits) if y_limits else ax.set_ylim(-max_pull, max_pull)
        ax.set_ylabel(r"Pull $\frac{fit - truth}{\sigma_{fit}}$", loc="top")
        ax.set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right")
        ax.legend(**legend_kwargs)

        plt.minorticks_on()
        plt.show()
        pass

    def diagnose_bootstraps(
        self,
        columns: list = None,
        mean_flag: float = 0.0,
        chi2_flag: float = 1.0,
        width_flag: float = 4.0,
    ) -> None:
        """TODO: define the function so that it checks all the columns
            (default maybe re/im params?) and reports if the mean and width are
            significantly poor, and/or a gaussian fit to is poor. The flag parameters
            adjusts what 'significant' is defined to be. Might want to also pass in the
            truth information too.
            Honestly this should be in itw own "diagnostics" class

        When using stdev to flag results, separate into cases where if the stdev is
        close to 0, then instead use a percentile to determine if the nominal/truth fit
        is "far" from the distributed samples (i.e. outside the 3sigma or X percentile)
        """
        pass

    def bootstrap_matrix(
        self, bins: List[int | str], columns: List[str], **kwargs
    ) -> None:
        """Plot matrix of the bootstrap fit results with nominal fit overlaid

        View the relations between fit parameters (selected with "columns") from the
        bootstrap fit results in a single bin (chosen with bin_index). The matrix is of
        the form: histogram + its kde on the diagonal, scatter plots on the lower
        triangle, and 2D kde's on the upper triangle. Every plot has a thick green line
        overlaid to indicate the nominal fit result values. Plots framed by a solid
        black line indicate a strong correlation between those variables (>0.7). All
        amplitudes are plotted in fit fractions.

        Args:
            bins (List[int|str]): mass bins to plot
            columns (List[str]): fit parameters to plot
            **kwargs: specifically for the seaborn pairplot object

        Raises:
            ValueError: if bootstrap df was not defined, since its optional upon init
        """

        if self.bootstrap_df is None:
            raise ValueError("Bootstrap df was not defined on instantiation")
        if not bins:
            raise ValueError("No bins selected to plot")
        if not columns:
            raise ValueError("No columns selected to plot")

        # convert any string bin ranges to ints
        for i, bin in enumerate(bins):
            if isinstance(bin, str):
                bins[i] = self.bootstrap_df.loc[self.bootstrap_df["bin_range"] == bin][
                    "bin"
                ][0]

        # make selection to reduce dataset size as we work with it
        cut_df = self.df.loc[bins].copy()
        cut_bootstrap_df = self.bootstrap_df[self.bootstrap_df["bin"].isin(bins)].copy()
        if self.is_truth:
            cut_truth_df = self.truth_df.loc[bins].copy()

        # we want plotted amplitudes to be their fit fraction, so use separate df
        plotter_df = cut_bootstrap_df[columns].copy()
        for col in plotter_df:
            if any(col in sublist for sublist in self.coherent_sums.values()):
                plotter_df.loc[:, col] = plotter_df[col].div(
                    cut_bootstrap_df["detected_events"], axis="index"
                )

        # scale truth phases to be the same sign as the nominal fit phases. Due to phase
        # ambiguity in model, the nominal fit can obtain +/-|truth_phase|.
        for phase in set(self.phase_differences.values()):
            cut_truth_df[phase] = cut_truth_df[phase].where(
                np.sign(cut_truth_df[phase]) == np.sign(cut_df[phase]),
                -1.0 * cut_truth_df[phase],
            )

        plotter_df["bin_range"] = cut_bootstrap_df["bin_range"]

        # get default color palette for hue plotting
        n_colors = plotter_df["bin_range"].nunique() if "bin_range" in plotter_df else 1
        palette = sns.color_palette(n_colors=n_colors)

        pg = sns.PairGrid(plotter_df, hue="bin_range", palette=palette, **kwargs)
        # overlay a kde on the hist to compare nominal fit line to kde peak
        # if many bins are plotted, remove the histogram
        if len(bins) < 3:
            pg.map_diag(sns.histplot, kde=True)
        else:
            pg.map_diag(sns.kdeplot)
        # scatter plots on both diagonals is redundant, so plot kde on upper. Levels of
        # kde are based off 3, 2, 1 sigma widths
        pg.map_upper(sns.kdeplot, levels=[0.003, 0.05, 0.32])
        pg.map_lower(sns.scatterplot)

        num_plots = len(columns)
        for row in range(num_plots):
            for col in range(num_plots):
                # Labels from last row (for x labels) or first column (for y labels)
                col_label = pg.axes[num_plots - 1, col].xaxis.get_label().get_text()
                row_label = pg.axes[row, 0].yaxis.get_label().get_text()

                # if a fit fraction is on the main plot, then make the nominal and
                # truth info a fit fraction too
                x_scaling, y_scaling = 1.0, 1.0  # default values
                x_truth_scale, y_truth_scale = 1.0, 1.0
                if any(col_label in sublist for sublist in self.coherent_sums.values()):
                    x_scaling = cut_df["detected_events"]
                    if self.is_truth:
                        x_truth_scale = cut_truth_df["detected_events"]
                if any(row_label in sublist for sublist in self.coherent_sums.values()):
                    y_scaling = cut_df["detected_events"]
                    if self.is_truth:
                        y_truth_scale = cut_truth_df["detected_events"]

                # plot band of MINUIT uncertainty centered around nominal fit
                for i in range(len(cut_df[col_label])):
                    pg.axes[row, col].axvspan(
                        xmin=(cut_df[col_label] - cut_df[f"{col_label}_err"])
                        .div(x_scaling)
                        .iloc[i],
                        xmax=(cut_df[col_label] + cut_df[f"{col_label}_err"])
                        .div(x_scaling)
                        .iloc[i],
                        color=palette[i],
                        alpha=0.2,
                    )

                # plot marker for truth (if applicable)
                if self.is_truth and row == col:
                    for i in range(len(cut_truth_df[col_label])):
                        pg.axes[row, col].axline(
                            xy1=(
                                cut_truth_df[col_label].div(x_truth_scale).iloc[i],
                                pg.axes[row, col].get_ylim()[0],
                            ),
                            xy2=(
                                cut_truth_df[col_label].div(x_truth_scale).iloc[i],
                                pg.axes[row, col].get_ylim()[1],
                            ),
                            color=palette[i],
                            linestyle="-",
                            alpha=0.6,
                            linewidth=2.0,
                        )

                if row != col:  # off-diagonals need y-axis span
                    for i in range(len(cut_df[row_label])):
                        pg.axes[row, col].axhspan(
                            ymin=(cut_df[row_label] - cut_df[f"{row_label}_err"])
                            .div(y_scaling)
                            .iloc[i],
                            ymax=(cut_df[row_label] + cut_df[f"{row_label}_err"])
                            .div(y_scaling)
                            .iloc[i],
                            color=palette[i],
                            alpha=0.2,
                        )
                    # plot "X" marker for truth value location
                    if self.is_truth:
                        for i in range(len(cut_truth_df[row_label])):
                            pg.axes[row, col].scatter(
                                cut_truth_df[col_label].div(x_truth_scale).iloc[i],
                                cut_truth_df[row_label].div(y_truth_scale).iloc[i],
                                color=palette[i],
                                marker="X",
                                linestyle="",
                                edgecolors="black",
                            )

                    # draw black box around plot if its correlation is above |0.7|
                    # uses an average if multiple bins are being plotted
                    corrs = []
                    for index in cut_df.index:
                        corr = (
                            cut_bootstrap_df[cut_bootstrap_df["bin"] == index][
                                [col_label, row_label]
                            ]
                            .corr()
                            .iat[0, 1]
                        )
                        corrs.append(corr)
                    if np.abs(np.average(corrs)) > 0.7:
                        pg.axes[row, col].patch.set_edgecolor("black")
                        pg.axes[row, col].patch.set_linewidth(3)

                # annotate plot with ratio between bootstrap stdev and MINUIT error
                # uses an average if multiple bins are being plotted
                if row == col:
                    ratios = []
                    for index, error in zip(
                        cut_df[f"{col_label}_err"].index, cut_df[f"{col_label}_err"]
                    ):
                        series = cut_bootstrap_df[cut_bootstrap_df["bin"] == index][
                            col_label
                        ]
                        # if distribution is being split at +/- pi boundaries, then
                        # wrap it to [0, 2pi] range and use that stdev
                        # NOTE: this is sensitive to outliers and assumes distributions
                        # near the boundaries do NOT have many points crossing zero
                        if (
                            np.min(series) + 180.0 < 4.0
                            and 180.0 - np.max(series) < 4.0
                        ):
                            series = series.apply(lambda x: x % 360)
                        ratios.append(np.std(series) / error)
                    pg.axes[row, col].text(
                        0.6,
                        0.9,
                        (
                            rf"$\frac{{\sigma_{{bootstrap}}}}{{\sigma_{{MINUIT}}}}$"
                            rf" = {np.average(ratios):.2f}"
                        ),
                        transform=pg.axes[row, col].transAxes,
                    )

        # change axes labels to prettier version. Must be done last since axes labels
        # are the primary way of attaining plot information
        for row in range(num_plots):
            for col in range(num_plots):
                col_label = pg.axes[row, col].xaxis.get_label().get_text()
                row_label = pg.axes[row, col].yaxis.get_label().get_text()
                fontsize = pg.axes[row, col].yaxis.get_label().get_fontsize()

                if any(
                    col_label in sublist for sublist in self.coherent_sums.values()
                ) or col_label in set(self.phase_differences.values()):
                    pg.axes[row, col].set_xlabel(
                        convert_amp_name(col_label), fontsize=fontsize
                    )
                if any(
                    row_label in sublist for sublist in self.coherent_sums.values()
                ) or row_label in set(self.phase_differences.values()):
                    pg.axes[row, col].set_ylabel(
                        convert_amp_name(row_label), fontsize=fontsize
                    )

        # adjust legend title
        pg.add_legend()
        for ax in pg.axes.flat:
            leg = ax.get_legend()
            if leg is not None:
                break
        if leg is None:
            leg = pg.legend

        leg.set_title("mass range")

        plt.show()
        pass

    def joyplot(
        self,
        columns: List[str],
        bins: List[int | str] = None,
        colors: List[tuple] = None,
        is_sparse_labels=True,
        overlap: float = 2.0,
        **kwargs,
    ) -> None:
        """Plot a joyplot (or ridgeplot) of the bootstrap fit results

        The joyplot is a series of kde's that are stacked on top of each other, with
        each kde representing the distribution of a fit parameter in a mass bin. Note
        that the x-axis is shared, so its recommended to only plot like-parameters and
        not to mix columns like phases and intensities. If the truth values are present,
        overlay them as a vertical line on each axis.

        Args:
            columns (List[str]): bootstrap df columns to plot
            bins (List[int|str], optional): mass bins to plot on the y-axis. Can use bin
                range strings, or bin # integers, or a mix of both. Defaults to None,
                and uses all bins.
            colors (List[tuple], optional): matplotlib colormap as a list. This means
                only qualitative maps can be used. Defaults to None, and uses
                list(matplotlib.colormaps["Accent"].colors).
            is_sparse_labels (bool, optional): when True, only bins that start at values
                divisible by 10 are labelled. Defaults to True.
            overlap (float, optional): amount of overlap between the kde's. Usually the
                most important parameter for making plots look good. Defaults to 2.0
            **kwargs: additional arguments to pass to the joyplot function
        Raises:
            ValueError: if bootstrap df was not defined, since its optional upon init
            ValueError: if requested columns are not in the bootstrap dataframe
        """

        if self.bootstrap_df is None:
            raise ValueError("Bootstrap df was not defined on instantiation")
        # check if requested columns are in the dataframe
        missing_columns = [
            col for col in columns if col not in self.bootstrap_df.columns
        ]
        if missing_columns:
            raise ValueError(
                f"The following columns are not in bootstrap_df: {missing_columns}"
            )

        # use all bins if none are specified
        if not bins:
            bins = self.bootstrap_df["bin"].unique()
        else:
            # convert any string bin ranges to ints
            for i, bin in enumerate(bins):
                if isinstance(bin, str):
                    bins[i] = self.bootstrap_df.loc[
                        self.bootstrap_df["bin_range"] == bin
                    ]["bin"].iloc[0]

        # default colors if none are specified
        if colors is None:
            colors = list(matplotlib.colormaps["Accent"].colors)

        # Ensure colormap matches the length of columns
        if len(colors) < len(columns):  # repeat colormap if too short
            colors = (colors * (len(columns) // len(colors) + 1))[: len(columns)]
        else:
            colors = colors[: len(columns)]

        # select bins from dataframes and determine min/max x values
        cut_bootstrap_df = self.bootstrap_df[self.bootstrap_df["bin"].isin(bins)].copy()
        x_max = cut_bootstrap_df[columns].max().max()
        x_min = cut_bootstrap_df[columns].min().min()
        if self.is_truth:
            cut_truth_df = self.truth_df.loc[bins].copy()
            x_max = max(x_max, cut_truth_df[columns].max().max())
            x_min = min(x_min, cut_truth_df[columns].min().min())

        # force x_min to be 0 if all values are coherent sums
        if all(
            any(col in sublist for sublist in self.coherent_sums.values())
            for col in columns
        ):
            x_min = 0.0
        # force x_min/x_max to be -180/180 if all values are phase differences
        elif all(col in set(self.phase_differences.values()) for col in columns):
            x_min, x_max = -180.0, 180.0

        # get ordered set of bin_ranges and make their low edge the y-axis labels
        bin_ranges = sorted(
            list(set(cut_bootstrap_df["bin_range"])),
            key=lambda x: float(x.split("-")[1]),
        )
        y_axis_labels = []
        for bin in bin_ranges:
            low_edge, high_edge = float(bin.split("-")[0]), float(bin.split("-")[1])
            # if not sparse labels, then add all labels
            if not is_sparse_labels:
                y_axis_labels.append(f"{low_edge:.2f}")
                continue

            tolerance = 1e-5
            # if first or last bin then add label
            if (
                bin_ranges.index(bin) == 0
                or bin_ranges.index(bin) == len(bin_ranges) - 1
            ):
                y_axis_labels.append(f"{low_edge:.2f}")
            # add label if bin start is a multiple of 10 within tolerance
            elif (
                abs(low_edge * 100 % 10) < tolerance
                or abs(low_edge * 100 % 10 - 10) < tolerance
            ):
                y_axis_labels.append(f"{low_edge:.2f}")
            # else add empty string to keep the axis clean
            else:
                y_axis_labels.append("")

        # setup joyplot default args
        default_args = {
            "linewidth": 1,
            "alpha": 0.7,
            "legend": True,
            "grid": "y",
            "yrot": 0,
            "loc": "lower right",
        }
        default_args.update(kwargs)  # overwrite defaults with user input

        # avoid a bug within joyplot that occurs when colors is a single element list
        joy_color = colors[0] if len(colors) == 1 else colors

        fig, axes = joypy.joyplot(
            cut_bootstrap_df,
            by="bin_range",
            column=columns,
            labels=y_axis_labels,
            range_style="own",
            x_range=[x_min, x_max],
            color=joy_color,
            overlap=overlap,
            **default_args,
        )

        # plot the truth value as a vertical line on each axis
        if self.is_truth:
            for bin, ax in zip(bins, axes):
                # get the max y value in each bin by getting the kde max (kde's are the
                # odd values in ax.lines)
                y_maxes = [
                    line.get_ydata().max()
                    for line in ax.lines
                    if ax.lines.index(line) % 2 != 0
                ]
                for col, line_color in zip(columns, colors):
                    ax.plot(
                        [cut_truth_df.loc[bin, col]] * 2,
                        [0.0, y_maxes[columns.index(col)]],
                        color=line_color,
                        alpha=1.0,
                        linestyle="-",
                        linewidth=2,
                        zorder=200,
                    )

        # Access and modify the legend text entries to pretty LaTeX names
        for ax in axes:
            legend = ax.get_legend()
            if legend:
                for text in legend.get_texts():
                    label = text.get_text()
                    if any(
                        label in sublist for sublist in self.coherent_sums.values()
                    ) or label in set(self.phase_differences.values()):
                        text.set_text(convert_amp_name(label))

        plt.show()

        pass

    def _reset_truth_phases(self, df: pd.DataFrame, mass_bins: List[int]) -> None:
        """Reset the phase differences of mass dependent truth csv

        The phase differences in the truth csv are not natively comparable. Because the
        truth info is from a mass dependent fit, we need to obtain the true phase
        differences by modifying each phase by its associated breit wigner in each mass
        bin. This function modifies the dataframe in place.

        Args:
            df (pd.DataFrame): truth dataframe
            mass_bins (List[int]): the mass bins the index column is associated with
        """

        def new_phase_dif(
            mass: float,
            amp1_re: float,
            amp1_im: float,
            bw_mass1: float,
            bw_width1: float,
            bw_l1: int,
            amp2_re: float,
            amp2_im: float,
            bw_mass2: float,
            bw_width2: float,
            bw_l2: int,
        ):
            complex_val1 = complex(amp1_re, amp1_im)
            complex_val2 = complex(amp2_re, amp2_im)
            bw1 = breit_wigner(mass, bw_mass1, bw_width1, bw_l1)
            bw2 = breit_wigner(mass, bw_mass2, bw_width2, bw_l2)
            phase1 = cmath.phase(bw1 * complex_val1)
            phase2 = cmath.phase(bw2 * complex_val2)

            return phase1 - phase2

        for pd in set(get_phase_differences(df).values()):
            l_to_int = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}
            amp1, amp2 = pd.split("_")

            jp1 = amp1[1:3]
            l1 = l_to_int[amp1[-1]]
            jp2 = amp2[1:3]
            l2 = l_to_int[amp2[-1]]

            df[pd] = np.vectorize(new_phase_dif)(
                mass_bins,
                df[f"{amp1}_re"],
                df[f"{amp1}_im"],
                df[f"{jp1}_mass"],
                df[f"{jp1}_width"],
                l1,
                df[f"{amp2}_re"],
                df[f"{amp2}_im"],
                df[f"{jp2}_mass"],
                df[f"{jp2}_width"],
                l2,
            )
        pass


def wrap_phases(df: pd.DataFrame = None, series: pd.Series = None) -> None:
    """Wrap phase differences to be from (-pi, pi] & convert from radians to degrees

    Two options of passing either a pandas dataframe, or series. The dataframe case
    handles avoiding editing any non phase difference columns. The series case is much
    simpler, and just applies the wrapping to each value

    Args:
        df (pd.DataFrame, optional): dataframe of fit results loaded from csv
            Defaults to pd.DataFrame([]).
        series (pd.Series, optional): pandas series of phase difference values.
            Defaults to pd.Series([], dtype=float).

    Returns:
        None: Edits the df or series itself
    """

    if df is None and series is None:
        raise ValueError(
            "Both parameters are None. Provide either a dataframe or a series."
        )

    if df is not None and series is not None:
        raise ValueError("Only dataframe or series should be passed, NOT both.")

    # wraps phase (in radians) to -pi < x <= pi and convert to degrees
    def wrap(phase):
        return np.rad2deg(np.angle(np.exp(1j * phase)))

    if series is not None:
        series.apply(wrap, inplace=True)
        return

    phase_diffs = get_phase_differences(df)
    for col in set(phase_diffs.values()):
        df[col] = df[col].apply(wrap)
        df[f"{col}_err"] = df[f"{col}_err"].apply(wrap)

    return


def parse_amplitude(amp: str) -> Dict[str, str]:
    """parse 'eJPmL' style amplitude string into its individual quantum numbers

    Args:
        amp (str): string like eJPmL, JPmL, JPL, eJPL, eJP, JP

    Returns:
        dict: keys = quantum numbers (e, J, P, m, l). values = found values from string,
        or "" if not found
    """
    re_dict = {
        "e": r"(?<![0-9]{1}[pm]{1})(?<!\d)[pm]",
        "jp": r"[0-9]{1}[pm]{1}",  # j and p always appear together
        "m": r"([0-9]{1}[pm]{1})+[qp0mn]",  # assumes always form of 'JPm'
        "l": r"[SPDF]",
    }

    result_dict = {"e": "", "j": "", "p": "", "m": "", "l": ""}

    for quantum_number, expression in re_dict.items():
        search = re.search(expression, amp)
        if search:
            # the search actually returns "JPm", so grab the last char if found
            if quantum_number == "m":
                result_dict["m"] = search.group()[-1]
            elif quantum_number == "jp":
                result_dict["j"] = search.group()[0]
                result_dict["p"] = search.group()[1]
            else:
                result_dict[quantum_number] = search.group()

    return result_dict


def get_coherent_sums(df: pd.DataFrame) -> Dict[str, set]:
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

        res = parse_amplitude(column)
        # only add to key if all elements of key are in the column
        for key in coherent_sums.keys():
            split_key = list(key.lower())
            if any(res[char] == "" for char in split_key):
                continue
            coh_sum = "".join([res[char] for char in split_key])
            coherent_sums[key].add(coh_sum)

    return {k: sorted(v) for k, v in coherent_sums.items()}  # sort each set


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


def human_format(num: float) -> str:
    """Converts orders of magnitude to letters i.e. 12369 -> 12.4K

    Args:
        num (float): positive floating point number

    Returns:
        str: first 3 significant digits with a character denoting its magnitude
    """
    num = float(f"{num:.3g}")
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    orders = ["", "K", "M", "B", "T"]
    return f"{str(num).rstrip('0').rstrip('.')}{orders[magnitude]}"


def breakup_momentum(
    parent_mass: float, daughter1_mass: float, daughter2_mass: float
) -> float:
    """Breakup momentum of a parent particle into two daughter particles

    Pythonized version of https://github.com/JeffersonLab/halld_sim/blob/master/src/
    libraries/AMPTOOLS_AMPS/breakupMomentum.cc. Returned value is independent of which
    particle is considered daughter1 or daughter2. Take care to units are consistent
    between masses.

    Args:
        parent_mass (float): mass of parent particle
        daughter1_mass (float): mass of one of the daughter particles
        daughter2_mass (float): mass of the other daughter particle

    Returns:
        float: Breakup momentum of the parent particle.
    """
    return np.sqrt(
        np.abs(
            np.power(parent_mass, 4)
            + np.power(daughter1_mass, 4)
            + np.power(daughter2_mass, 4)
            - 2.0 * np.square(parent_mass) * np.square(daughter1_mass)
            - 2.0 * np.square(parent_mass) * np.square(daughter2_mass)
            - 2.0 * np.square(daughter1_mass) * np.square(daughter2_mass)
        )
    ) / (2.0 * parent_mass)


def barrier_factor(
    parent_mass: float, daughter1_mass: float, daughter2_mass: float, l: int
) -> float:
    """Blatt-Weisskopf barrier factor for a given parent-daughter system

    Pythonized version of https://github.com/JeffersonLab/halld_sim/blob/master/src/
    libraries/AMPTOOLS_AMPS/barrierFactor.cc. This parameterization often suppresses
    values near threshold, and amplifies values at higher masses, especially for large
    l values.

    Args:
        parent_mass (float): mass of the parent particle
        daughter1_mass (float): mass of the first daughter particle
        daughter2_mass (float): mass of the second daughter particle
        l (int): orbital angular momentum quantum number of the parent particle

    Returns:
        float: barrier factor
    """

    q = np.abs(breakup_momentum(parent_mass, daughter1_mass, daughter2_mass))
    z = np.square(q) / np.square(0.1973)

    match l:
        case 0:
            barrier = 1.0
        case 1:
            barrier = (2.0 * z) / (z + 1.0)
        case 2:
            barrier = (13.0 * np.square(z)) / (np.square(z - 3.0) + 9.0 * z)
        case 3:
            barrier = (277.0 * np.power(z, 3)) / (
                z * np.square(z - 15.0) + 9.0 * np.square(2.0 * z - 5.0)
            )
        case 4:
            barrier = (12746.0 * np.power(z, 4)) / (
                np.square((np.square(z) - 45.0 * z + 105.0))
                + 25.0 * z * np.square(2.0 * z - 21.0)
            )
        case _:
            barrier = 0.0

    return np.sqrt(barrier)


def breit_wigner(
    mass: float,
    bw_mass: float,
    bw_width: float,
    bw_l: int,
    daughter1_mass: float = 0.1349768,
    daughter2_mass: float = 0.78266,
) -> complex:
    """Halld_sim parameterization of the breit wigner function

    To avoid discrepancies, this function copies the parameterization of the breit
    wigner function found in https://github.com/JeffersonLab/halld_sim/blob/master/src/
    libraries/AMPTOOLS_AMPS/BreitWigner.cc.

    NOTE: all masses / widths must be in GeV. Daughter particle masses are typically
    obtained from event 4-vectors, but since they're not available post-fit, we
    approximate them by constraining their mass to the pdg value

    Args:
        mass (float): mass value to evaluate breit wigner function at. If using a mass
            bin, approximate by passing the center of the bin
        bw_mass (float): breit wigner central mass
        bw_width (float): breit wigner width
        bw_l (int): orbital angular momentum of the breit wigner i.e. rho->(omega,pi0)
            is in the P-wave, so bw_l = 1
        daughter1_mass (float, optional): mass of 1st daughter particle.
            Defaults to 0.1349768 for the pi0 mass
        daughter2_mass (float, optional): mass of 2nd daughter particle.
            Defaults to 0.78266 for the omega mass

    Returns:
        complex: value of the breit wigner at the mass value
    """

    q0 = np.abs(breakup_momentum(bw_mass, daughter1_mass, daughter2_mass))
    q = np.abs(breakup_momentum(mass, daughter1_mass, daughter2_mass))

    F0 = barrier_factor(bw_mass, daughter1_mass, daughter2_mass, bw_l)
    F = barrier_factor(mass, daughter1_mass, daughter2_mass, bw_l)

    width = bw_width * (bw_mass / mass) * (q / q0) * np.square(F / F0)

    numerator = complex(np.sqrt((bw_mass * bw_width) / np.pi), 0.0)
    denominator = complex(np.square(bw_mass) - np.square(mass), -1.0 * bw_mass * width)

    return F * numerator / denominator


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
