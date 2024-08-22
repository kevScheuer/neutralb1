"""Collection of tools useful for analyzing PWA fit results

Primary class is the Plotter, which ideally handles most standard plots of interest when 
analyzing fit results
"""

import cmath
import itertools
import re

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


class Plotter:
    def __init__(
        self,
        df: pd.DataFrame,
        data_df: pd.DataFrame,
        bootstrap_df: pd.DataFrame = None,
        truth_df: pd.DataFrame = None,
    ) -> None:
        """Initialize object with pandas dataframe

        Args:
            df (pd.DataFrame): FitResults from AmpTools, made by fitsToCsv.C
            data_df (pd.DataFrame): raw data points that AmpTools is fitting to
            bootstrap_df (pd.DataFrame, optional): all bootstrapped fit results,
                labelled by bin number. Primary use is for plotting the distribution of
                each parameter and using its width to estimate the error of the nominal
                fit result in df.
            truth_df (pd.DataFrame, optional): only used when df is a from fit to
                generated signal MC. Contains the "true" parameters that were used to
                generate the MC sample. When non-empty, every plot will display a
                solid line for the truth information
        """
        self.df = df
        self.data_df = data_df
        self.bootstrap_df = bootstrap_df
        self.truth_df = truth_df

        self.is_truth = False
        if self.truth_df is not None:
            self.is_truth = True

        self._mass_bins = self.data_df["mass_mean"]

        wrap_phases(self.df)
        if self.bootstrap_df is not None:
            wrap_phases(self.bootstrap_df)
        if self.is_truth:
            self._truth_coherent_sums = get_coherent_sums(self.truth_df)
            self._truth_phase_differences = get_phase_differences(self.truth_df)
            self._reset_truth_phases(self.truth_df, self._mass_bins)
            wrap_phases(self.truth_df)

        self._bin_width = (
            self.data_df["mass_high_edge"] - self.data_df["mass_low_edge"]
        )[0]

        self._coherent_sums = get_coherent_sums(self.df)
        self._phase_differences = get_phase_differences(self.df)
        pass

    def jp(self) -> None:
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

        fig, ax = plt.subplots()
        # plot data
        ax.errorbar(
            self._mass_bins,
            self.data_df["bin_contents"],
            self.data_df["bin_error"],
            self._bin_width / 2,
            "k.",
            label="Total (GlueX Phase-I)",
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
        # plot isotropic background if there
        if "Bkgd" in self.df.columns:
            ax.errorbar(
                self._mass_bins,
                self.df["Bkgd"],
                self.df[f"Bkgd_err"],
                self._bin_width / 2,
                label="Iso. Background",
                linestyle="",
                **jp_map["Bkgd"],
            )
        # plot each jp contribution
        for jp in self._coherent_sums["JP"]:
            ax.errorbar(
                self._mass_bins,
                self.df[jp],
                self.df[f"{jp}_err"],
                self._bin_width / 2,
                label=convert_amp_name(jp),
                linestyle="",
                **jp_map[jp],
            )

        # plot bkgd and jp contributions for truth df
        if self.is_truth:
            if "Bkgd" in self.truth_df.columns:
                ax.plot(
                    self._mass_bins,
                    self.truth_df["Bkgd"],
                    linestyle="-",
                    marker="",
                    color=jp_map["Bkgd"]["color"],
                )
            for jp in self._truth_coherent_sums["JP"]:
                ax.plot(
                    self._mass_bins,
                    self.truth_df[jp],
                    linestyle="-",
                    marker="",
                    color=jp_map[jp]["color"],
                )

        ax.set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right")
        ax.set_ylabel(f"Events / {self._bin_width:.3f} GeV", loc="top")
        ax.set_ylim(bottom=0.0)

        ax.legend()

        ax.grid(True, alpha=0.8)
        plt.show()
        pass

    def intensities(self, is_fit_fraction: bool = False, sharey=True) -> None:
        """Plot all the amplitude intensities in a grid format

        Since the matrix plot is generally too small to see bin to bin features, this
        method plots every amplitude's intensity contribution in a grid format.
        Columns = m-projections, rows = JPL combinations. Reflectivities are plotted
        together on each subplot

        Args:
            is_fit_fraction (bool, optional): Scales all values by dividing them by the
                total intensity in each bin. Defaults to False.
            sharey (bool, optional): Sets each row to share a y-axis scale. Defaults to
            True
        """

        char_to_int = {
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
        int_to_char = {-1: "m", 0: "0", +1: "p"}
        pm_dict = {"m": "-", "p": "+"}

        # need to sort on the integer versions of the m-projections
        m_ints = sorted({char_to_int[JPm[-1]] for JPm in self._coherent_sums["JPm"]})

        fig, axs = plt.subplots(
            len(self._coherent_sums["JPL"]),
            len(m_ints),
            sharex=True,
            sharey=sharey,
            figsize=(15, 10),
            dpi=100,
        )

        total = self.df["detected_events"] if is_fit_fraction else 1
        total_err = self.df["detected_events_err"] if is_fit_fraction else 0

        if self.is_truth:
            truth_total = self.truth_df["detected_events"] if is_fit_fraction else 1

        # iterate through JPL (sorted like S, P, D, F wave) and sorted m-projections
        for row, jpl in enumerate(
            sorted(self._coherent_sums["JPL"], key=lambda JPL: char_to_int[JPL[-1]])
        ):
            for col, m in enumerate(m_ints):
                JPmL = f"{jpl[0:2]}{int_to_char[m]}{jpl[-1]}"

                if row == 0:
                    axs[row, col].set_title(f"m={char_to_int[JPmL[-2]]}")
                if col == 0:
                    axs[row, col].set_ylabel(
                        rf"${JPmL[0]}^{{{pm_dict[JPmL[1]]}}}{JPmL[-1]}$",
                    )

                # plot the negative reflectivity contribution
                if "m" + JPmL in self._coherent_sums["eJPmL"]:
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
                        "o",
                        color="blue",
                        label=r"$\epsilon=-1$",
                    )
                    if (
                        self.is_truth
                        and "m" + JPmL in self._truth_coherent_sums["eJPmL"]
                    ):
                        axs[row, col].plot(
                            self._mass_bins,
                            self.truth_df["m" + JPmL] / truth_total,
                            linestyle="-",
                            marker="",
                            color="blue",
                        )
                # plot the negative reflectivity contribution
                if "p" + JPmL in self._coherent_sums["eJPmL"]:
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
                        "o",
                        color="red",
                        label=r"$\epsilon=+1$",
                    )
                    if (
                        self.is_truth
                        and "p" + JPmL in self._truth_coherent_sums["eJPmL"]
                    ):
                        axs[row, col].plot(
                            self._mass_bins,
                            self.truth_df["p" + JPmL] / truth_total,
                            linestyle="-",
                            marker="",
                            color="red",
                        )

        # plot grid lines
        for ax in axs.reshape(-1):
            ax.grid(True, alpha=0.8)
            ax.set_ylim(bottom=0)

        # figure cosmetics
        fig.text(0.5, 0.04, r"$\omega\pi^0$ inv. mass (GeV)", ha="center", fontsize=20)
        fig.text(
            0.04,
            0.5,
            f"Events / {self._bin_width:.3f} GeV",
            ha="center",
            rotation="vertical",
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
        if amp1 not in self._coherent_sums["eJPmL"]:
            raise ValueError(f"Amplitude {amp1} not found in dataset")
        if amp2 not in self._coherent_sums["eJPmL"]:
            raise ValueError(f"Amplitude {amp2} not found in dataset")
        if amp1[0] != amp2[0]:
            raise ValueError(f"Amplitudes must be from same reflectivity")

        phase_dif = self._phase_differences[(amp1, amp2)]
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

        if self.is_truth and (amp1, amp2) in self._truth_phase_differences:
            ax.plot(
                self._mass_bins,
                self.truth_df[self._truth_phase_differences[(amp1, amp2)]],
                linestyle="-",
                marker="",
                color=color,
            )

        ax.legend()

        ax.set_yticks(np.linspace(-180, 180, 5))  # force to be in pi/2 intervals
        ax.set_ylim([-180, 180])
        ax.set_ylabel(r"Phase Diff. ($^{\circ}$)", loc="top")
        ax.set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right")
        ax.grid(True, alpha=0.8)

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
        if amp1 not in self._coherent_sums["eJPmL"]:
            raise ValueError(f"Amplitude {amp1} not found in dataset")
        if amp2 not in self._coherent_sums["eJPmL"]:
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

        if (
            self.is_truth
            and amp1 in self._truth_coherent_sums["eJPmL"]
            and amp2 in self._truth_coherent_sums["eJPmL"]
        ):
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
        phase_dif = self._phase_differences[(amp1, amp2)]
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
        if self.is_truth and (amp1, amp2) in self._truth_phase_differences:
            axs[1].plot(
                self._mass_bins,
                self.truth_df[self._truth_phase_differences[(amp1, amp2)]],
                linestyle="-",
                marker="",
                color=color,
            )

        # cosmetics
        for ax in axs.reshape(-1):
            ax.grid(True, alpha=0.8)

        axs[0].set_ylim(bottom=0.0)
        axs[0].set_ylabel(f"Events / {self._bin_width:.3f} GeV", loc="top")

        axs[1].set_yticks(np.linspace(-180, 180, 5))  # force to be in pi/2 intervals
        axs[1].set_ylim([-180, 180])
        axs[1].set_ylabel(r"Phase Diff. ($^{\circ}$)", loc="center")
        axs[1].set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right")

        axs[0].legend(loc="upper right")

        plt.show()
        pass

    def matrix(self) -> None:
        """Scatterplot-like matrix of subplots, for intensities and phase differences

        The diagonal of the plot contains the intensity of each wave, in both
        reflectivities. The off diagonals are the phase difference between each wave,
        with positive reflectivity on the upper triangle and negative on the lower.
        Amplitude names are labelled on the first row and column.

        NOTE: This plot will become cramped for large models, may need dpi adjusted
        """

        fig, axs = plt.subplots(
            len(self._coherent_sums["JPmL"]),
            len(self._coherent_sums["JPmL"]),
            sharex=True,
            figsize=(13, 8),
            dpi=500,
        )
        # only the top left corner plot will have ticks for the data scale
        max_diag = max([self.df[x].max() for x in self._coherent_sums["eJPmL"]])
        data_y_ticks = np.linspace(0, max_diag, 4)
        axs[0, 0].set_yticks(data_y_ticks, [human_format(num) for num in data_y_ticks])
        for i, JPmL in enumerate(self._coherent_sums["JPmL"]):
            for j, JPmL_dif in enumerate(self._coherent_sums["JPmL"]):
                # change tick label sizes for all plots and the max
                axs[i, j].tick_params("both", labelsize=6)
                axs[i, j].set_ylim(top=max_diag)

                # write mass ticks for only the bottom row
                if i == len(self._coherent_sums["JPmL"]) - 1:
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
                        label=r"$\epsilon=-1$",
                    )
                    pos_plot = axs[i, j].errorbar(
                        self._mass_bins,
                        self.df[f"p{JPmL}"],
                        self.df[f"p{JPmL}_err"],
                        self._bin_width / 2,
                        "o",
                        color="red",
                        markersize=2,
                        label=r"$\epsilon=+1$",
                    )
                    if self.is_truth and JPmL in self._truth_coherent_sums["JPmL"]:
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

                    phase_dif = self._phase_differences[(f"p{JPmL}", f"p{JPmL_dif}")]

                    axs[i, j].errorbar(
                        self._mass_bins,
                        self.df[phase_dif],
                        abs(self.df[f"{phase_dif}_err"]),
                        self._bin_width / 2,
                        linestyle="",
                        color="red",
                        marker=".",
                        elinewidth=1,
                        ms=1,
                    )
                    # plot flipped sign due to sign ambiguity
                    axs[i, j].errorbar(
                        self._mass_bins,
                        -self.df[phase_dif],
                        abs(self.df[f"{phase_dif}_err"]),
                        self._bin_width / 2,
                        linestyle="",
                        color="red",
                        marker=".",
                        elinewidth=1,
                        ms=1,
                    )
                    if (
                        self.is_truth
                        and (f"p{JPmL}", f"p{JPmL_dif}")
                        in self._truth_phase_differences
                    ):
                        axs[i, j].plot(
                            self._mass_bins,
                            self.truth_df[
                                self._truth_phase_differences[
                                    (f"p{JPmL}", f"p{JPmL_dif}")
                                ]
                            ],
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

                    phase_dif = self._phase_differences[(f"m{JPmL}", f"m{JPmL_dif}")]

                    axs[i, j].errorbar(
                        self._mass_bins,
                        self.df[phase_dif],
                        self.df[f"{phase_dif}_err"].abs(),
                        self._bin_width / 2,
                        linestyle="",
                        color="blue",
                        marker=".",
                        elinewidth=1,
                        ms=1,
                    )
                    # plot flipped sign due to sign ambiguity
                    axs[i, j].errorbar(
                        self._mass_bins,
                        -self.df[phase_dif],
                        self.df[f"{phase_dif}_err"].abs(),
                        self._bin_width / 2,
                        linestyle="",
                        color="blue",
                        marker=".",
                        elinewidth=1,
                        ms=1,
                    )
                    if (
                        self.is_truth
                        and (f"m{JPmL}", f"m{JPmL_dif}")
                        in self._truth_phase_differences
                    ):
                        axs[i, j].plot(
                            self._mass_bins,
                            self.truth_df[
                                self._truth_phase_differences[
                                    (f"m{JPmL}", f"m{JPmL_dif}")
                                ]
                            ],
                            linestyle="-",
                            marker="",
                            color="blue",
                        )

        for ax in axs.reshape(-1):
            ax.grid(True, axis="both", alpha=0.8)

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
            for eJPmL in self._coherent_sums["eJPmL"]:
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
                phase = self.df[self._phase_differences[(D_wave, S_wave)]]
                phase_err = self.df[
                    f"{self._phase_differences[(D_wave, S_wave)]}_err"
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

    def bootstrap_matrix(self, bins: list, columns: list, **kwargs) -> None:
        """Plot matrix of the bootstrap fit results with nominal fit overlaid

        View the relations between fit parameters (selected with "columns") from the
        bootstrap fit results in a single bin (chosen with bin_index). The matrix is of
        the form: histogram + its kde on the diagonal, scatter plots on the lower
        triangle, and 2D kde's on the upper triangle. Every plot has a thick green line
        overlaid to indicate the nominal fit result values. Plots framed by a solid
        black line indicate a strong correlation between those variables (>0.7). All
        amplitudes are plotted in fit fractions.

        TODO: add a ratio for comparing percentile to stdev predictions. Might use
            this later for case delineation between gauss and non-gauss (near 0) results

        Args:
            bins (list): bins that bootstrapped fits are associated with i.e. the
                1st bin is index=0 in the data dataframe, so all bootstrap fits
                associated with that bin have a value of 0 in the "bin" column of the
                bootstrap dataframe
            columns (list): fit parameters to plot
            **kwargs: specifically for the seaborn pairplot object

        Raises:
            ValueError: if bootstrap df was not defined, since its optional upon init
        """

        if self.bootstrap_df is None:
            raise ValueError("Bootstrap df was not defined on instantiation")

        # make selection to reduce dataset size as we work with it
        cut_df = self.df.loc[bins].copy()
        cut_bootstrap_df = self.bootstrap_df[self.bootstrap_df["bin"].isin(bins)].copy()
        if self.is_truth:
            cut_truth_df = self.truth_df.loc[bins].copy()

        # we want plotted amplitudes to be their fit fraction, so use separate df
        plotter_df = cut_bootstrap_df[columns].copy()
        for col in plotter_df:
            if any(col in sublist for sublist in self._coherent_sums.values()):
                plotter_df.loc[:, col] = plotter_df[col].div(
                    cut_bootstrap_df["detected_events"], axis="index"
                )

        # scale truth phases to be the same sign as the nominal fit phases. Due to phase
        # ambiguity in model, the nominal fit can obtain +/-|truth_phase|.
        for truth_phase in set(self._truth_phase_differences.values()):
            amp1, amp2 = truth_phase.split("_")
            df_phase = self._phase_differences[(amp1, amp2)]
            cut_truth_df[truth_phase] = cut_truth_df[truth_phase].where(
                np.sign(cut_truth_df[truth_phase]) == np.sign(cut_df[df_phase]),
                -1.0 * cut_truth_df[truth_phase],
            )

        plotter_df["bin"] = cut_bootstrap_df["bin"]
        plotter_df.astype({"bin": "int8"})

        # get default color palette for hue plotting
        n_colors = plotter_df["bin"].nunique() if "bin" in plotter_df else 1
        palette = sns.color_palette(n_colors=n_colors)

        pg = sns.PairGrid(plotter_df, hue="bin", palette=palette, **kwargs)
        # overlay a kde on the hist to compare nominal fit line to kde peak
        pg.map_diag(sns.histplot, kde=True)
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
                if any(
                    col_label in sublist for sublist in self._coherent_sums.values()
                ):
                    x_scaling = cut_df["detected_events"]
                    if self.is_truth:
                        x_truth_scale = cut_truth_df["detected_events"]
                if any(
                    row_label in sublist for sublist in self._coherent_sums.values()
                ):
                    y_scaling = cut_df["detected_events"]
                    if self.is_truth:
                        y_truth_scale = cut_truth_df["detected_events"]

                # plot band of MINUIT uncertainty centered around nominal fit
                for i in range(len(cut_df[col_label])):
                    pg.axes[row, col].axvspan(
                        xmin=(cut_df[col_label] - cut_df[f"{col_label}_err"] / 2)
                        .div(x_scaling)
                        .iloc[i],
                        xmax=(cut_df[col_label] + cut_df[f"{col_label}_err"] / 2)
                        .div(x_scaling)
                        .iloc[i],
                        color=palette[i],
                        alpha=0.2,
                    )

                # plot marker for truth (if applicable). The "partition" line
                # will flip a string around the character "_". Needs to be done since
                # the truth df's phases (amp1_amp2) are not necessarily in the same
                # order in the normal df (could be amp2_amp1)
                if (
                    self.is_truth
                    and row == col
                    and (
                        col_label in cut_truth_df
                        or "".join(col_label.partition("_")[::-1]) in cut_truth_df
                    )
                ):
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
                            linestyle="--",
                            alpha=0.6,
                            linewidth=2.0,
                        )

                if row != col:  # off-diagonals need y-axis span
                    for i in range(len(cut_df[row_label])):
                        pg.axes[row, col].axhspan(
                            ymin=(cut_df[row_label] - cut_df[f"{row_label}_err"] / 2)
                            .div(y_scaling)
                            .iloc[i],
                            ymax=(cut_df[row_label] + cut_df[f"{row_label}_err"] / 2)
                            .div(y_scaling)
                            .iloc[i],
                            color=palette[i],
                            alpha=0.2,
                        )
                    # plot "X" marker for truth value location
                    if self.is_truth and (
                        row_label in cut_truth_df
                        or "".join(row_label.partition("_")[::-1]) in cut_truth_df
                    ):
                        for i in range(len(cut_truth_df[row_label])):
                            pg.axes[row, col].scatter(
                                cut_truth_df[col_label].div(x_scaling).iloc[i],
                                cut_truth_df[row_label].div(y_scaling).iloc[i],
                                color=palette[i],
                                marker="P",
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
                    # take absolute value to avoid stdev being biased by possible
                    # bi-modal nature of phase diff results (due to sign ambiguity)
                    for index, error in zip(
                        cut_df[f"{col_label}_err"].index, cut_df[f"{col_label}_err"]
                    ):
                        stdev = np.std(
                            np.abs(
                                cut_bootstrap_df[cut_bootstrap_df["bin"] == index][
                                    col_label
                                ]
                            )
                        )
                        ratios.append(stdev / error)
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
                    col_label in sublist for sublist in self._coherent_sums.values()
                ) or col_label in set(self._phase_differences.values()):
                    pg.axes[row, col].set_xlabel(
                        convert_amp_name(col_label), fontsize=fontsize
                    )
                if any(
                    row_label in sublist for sublist in self._coherent_sums.values()
                ) or row_label in set(self._phase_differences.values()):
                    pg.axes[row, col].set_ylabel(
                        convert_amp_name(row_label), fontsize=fontsize
                    )

        plt.show()
        pass

    def _reset_truth_phases(self, df: pd.DataFrame, mass_bins: list) -> None:
        """Reset the phase differences of mass dependent truth csv

        The phase differences in the truth csv are not natively comparable. Because the
        truth info is from a mass dependent fit, we need to obtain the true phase
        differences by modifying each phase by its associated breit wigner in each mass
        bin. This function modifies the dataframe in place.

        Args:
            df (pd.DataFrame): truth dataframe
            mass_bins (list): the mass bins the index column is associated with
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
            cval1 = complex(amp1_re, amp1_im)
            cval2 = complex(amp2_re, amp2_im)
            bw1 = breit_wigner(mass, bw_mass1, bw_width1, bw_l1)
            bw2 = breit_wigner(mass, bw_mass2, bw_width2, bw_l2)
            phase1 = cmath.phase(bw1 * cval1)
            phase2 = cmath.phase(bw2 * cval2)

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


def wrap_phases(
    df: pd.DataFrame = pd.DataFrame([]), series: pd.Series = pd.Series([], dtype=float)
) -> None:
    """Wrap phase differences to be from -pi/2 to pi/2 & convert from radians to degrees

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

    # wraps phase (in radians) to -pi/2 < x < pi/2 and convert to degrees
    def wrap(phase):
        if phase > np.pi / 2 or phase < -np.pi / 2:
            phase = (phase + np.pi / 2) % (np.pi) - np.pi / 2
        return np.rad2deg(phase)

    if not series.empty:
        if not df.empty:
            raise ValueError(
                "Only dataframe or series should be passed, NOT both. Exiting"
            )
        series = series.apply(wrap)

    if df.empty:
        raise ValueError("Parameters are both empty. Exiting")

    phase_diffs = get_phase_differences(df)
    for col in set(phase_diffs.values()):
        df[col] = df[col].apply(wrap)
        df[f"{col}_err"] = df[f"{col}_err"].apply(wrap)

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
    m_proj_dict = {"n": -2, "m": -1, "0": 0, "p": +1, "q": +2, "": ""}

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
    obtained from event 4-vectors, but since they're not avaialable post-fit, we
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

    def breakup_momentum(m0: float, m1: float, m2: float):
        # breakup momenta of parent (m0) -> daughter particles (m1,m2) in center of
        # momenta frame
        return np.sqrt(
            np.abs(
                np.power(m0, 4)
                + np.power(m1, 4)
                + np.power(m2, 4)
                - 2.0 * np.square(m0) * np.square(m1)
                - 2.0 * np.square(m0) * np.square(m2)
                - 2.0 * np.square(m1) * np.square(m2)
            )
        ) / (2.0 * m0)

    def barrier_factor(q: float, l: int):
        # barrier factor suppresion based on angular momenta
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

    q0 = np.abs(breakup_momentum(bw_mass, daughter1_mass, daughter2_mass))
    q = np.abs(breakup_momentum(mass, daughter1_mass, daughter2_mass))

    F0 = barrier_factor(q0, bw_l)
    F = barrier_factor(q, bw_l)

    width = bw_width * (bw_mass / mass) * (q / q0) * np.square(F / F0)

    numerator = complex(np.sqrt((bw_mass * bw_width) / np.pi), 0.0)
    denominator = complex(np.square(bw_mass) - np.square(mass), -1.0 * bw_mass * width)

    return F * numerator / denominator
