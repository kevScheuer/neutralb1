"""Collection of tools useful for analyzing PWA fit results

NOTE: Remember to use red=pos refl, blue=neg refl standard now

TODO: As cross section scripts get developed, maybe make a second plot class for them
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
    def __init__(self, df: pd.DataFrame, data_df: pd.DataFrame) -> None:
        """Initialize object with paths to csv files

        Args:
            df (pd.DataFrame): dataframe that holds all FitResults, made by fitsToCsv.C
            data_df (pd.DataFrame): dataframe of the raw data that is being fit to
        """
        self.df = df
        self.data_df = data_df

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

    def intensities(self, is_fit_fraction: bool = False) -> None:
        """Plot all the amplitude intensities in a grid format

        Since the matrix plot is generally too small to see bin to bin features, this
        method plots every amplitude's intensity contribution in a grid format.
        Columns = m-projections, rows = JPL combinations. Reflectivities are plotted
        together on each subplot

        Args:
            is_fit_fraction (bool, optional): Scales all values by dividing them by the
                total intensity in each bin. Defaults to False.
        """

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

        total = self.data_df["bin_contents"] if is_fit_fraction else 1
        total_err = self.data_df["bin_error"] if is_fit_fraction else 0

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
                        markersize=2,
                        label=r"$\epsilon=-1$",
                    )
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

    def mass_phase(self, amp1: str, amp2: str) -> None:
        """Plot the amplitude intensities with their phase difference together

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

        # setup colors for phase and intensity plot
        if amp1[0] == "p":
            phase_color = "red"
            color1 = "lightcoral"
            color2 = "orangered"
        else:
            phase_color = "blue"
            color1 = "darkturquoise"
            color2 = "darkblue"

        fig, axs = plt.subplots(
            2,
            1,
            sharex=True,
            gridspec_kw={"wspace": 0.0, "hspace": 0.07},
            height_ratios=[3, 1],
        )

        axs[0].errorbar(
            self._mass_bins,
            self.df[amp1],
            self.df[f"{amp1}_err"],
            self._bin_width / 2,
            "o",
            color=color1,
            label=convert_amp_name(amp1),
        )
        axs[0].errorbar(
            self._mass_bins,
            self.df[amp2],
            self.df[f"{amp2}_err"],
            self._bin_width / 2,
            "s",
            color=color2,
            label=convert_amp_name(amp2),
        )

        phase_dif = self._phase_differences[(amp1, amp2)]
        axs[1].errorbar(
            self._mass_bins,
            self.df[phase_dif].apply(np.rad2deg),
            self.df[phase_dif + "_err"].abs().apply(np.rad2deg),
            self._bin_width / 2,
            linestyle="",
            marker=".",
            color=phase_color,
        )
        axs[1].errorbar(
            self._mass_bins,
            -self.df[phase_dif].apply(np.rad2deg),
            self.df[phase_dif + "_err"].abs().apply(np.rad2deg),
            self._bin_width / 2,
            linestyle="",
            marker=".",
            color=phase_color,
        )

        # cosmetics
        for ax in axs.reshape(-1):
            ax.grid(True, alpha=0.8)

        axs[0].set_ylim(bottom=0.0)
        axs[0].set_ylabel(f"Events / {self._bin_width:.3f} GeV", loc="top", fontsize=12)

        axs[1].set_yticks(np.linspace(-180, 180, 5))  # force to be in pi/2 intervals
        axs[1].set_ylim([-180, 180])
        axs[1].set_ylabel(r"Phase Diff. ($^{\circ}$)", loc="center", fontsize=12)
        axs[1].set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right", fontsize=18)

        fig.legend(fontsize=14, loc="upper right")

        plt.show()
        pass

    def matrix(self) -> None:
        """Scatterplot-like matrix of subplots, for intensities and phase differences

        The diagonal of the plot contains the intensity of each wave, in both
        reflectivities, with the total data. The off diagonals are the phase difference
        between each wave, with positive reflectivity on the upper triangle and negative
        on the lower. Amplitude names are labelled on the first row and column.

        NOTE: This plot will become cramped for large models, may need dpi adjusted
        """

        fig, axs = plt.subplots(
            len(self._coherent_sums["JPmL"]),
            len(self._coherent_sums["JPmL"]),
            sharex=True,
            figsize=(13, 8),
            dpi=100,
        )
        # only the top left corner plot will have ticks for the data scale
        data_y_ticks = np.linspace(0, self.data_df["bin_contents"].max(), 4)
        axs[0, 0].set_yticks(data_y_ticks, [human_format(num) for num in data_y_ticks])
        for i, JPmL in enumerate(self._coherent_sums["JPmL"]):
            for j, JPmL_dif in enumerate(self._coherent_sums["JPmL"]):
                # change tick label sizes for all plots
                axs[i, j].tick_params("both", labelsize=6)

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

                    # first plot the data with its error
                    data_plot = axs[i, j].errorbar(
                        self._mass_bins,
                        self.data_df["bin_contents"],
                        self.data_df["bin_error"],
                        self._bin_width / 2,
                        "k.",
                        label="Data",
                        markersize=4,
                    )

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
                        self.df[phase_dif].apply(np.rad2deg),
                        abs(self.df[f"{phase_dif}_err"]).apply(np.rad2deg),
                        self._bin_width / 2,
                        linestyle="",
                        color="red",
                        marker="o",
                        ms=3,
                    )
                    # plot flipped sign due to sign ambiguity
                    axs[i, j].errorbar(
                        self._mass_bins,
                        -self.df[phase_dif].apply(np.rad2deg),
                        abs(self.df[f"{phase_dif}_err"]).apply(np.rad2deg),
                        self._bin_width / 2,
                        linestyle="",
                        color="red",
                        marker="o",
                        ms=3,
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
                        self.df[phase_dif].apply(np.rad2deg),
                        self.df[phase_dif + "_err"].abs().apply(np.rad2deg),
                        self._bin_width / 2,
                        linestyle="",
                        color="blue",
                        marker="s",
                        ms=3,
                    )
                    # plot flipped sign due to sign ambiguity
                    axs[i, j].errorbar(
                        self._mass_bins,
                        -self.df[phase_dif].apply(np.rad2deg),
                        self.df[phase_dif + "_err"].abs().apply(np.rad2deg),
                        self._bin_width / 2,
                        linestyle="",
                        color="blue",
                        marker="s",
                        ms=3,
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
        fig.legend(
            handles=[data_plot, pos_plot, neg_plot], fontsize=12, loc="upper right"
        )

        plt.show()
        pass

    def ds_ratio(self) -> None:
        """_summary_
        Whats plotted changes depending on whether the model had a defined D/S ratio.
        """
        #
        if "dsratioTEMP" in self.df.columns:
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
                self.df["dphase"].apply(np.rad2deg),
                self.df["dphase_err"].abs().apply(np.rad2deg),
                self._bin_width / 2,
                marker=".",
                linestyle="",
                color="gray",
            )
            ax_phase.errorbar(
                self._mass_bins,
                -self.df["dphase"].apply(np.rad2deg),
                self.df["dphase_err"].abs().apply(np.rad2deg),
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
            marker_map = {"p": "^", "0": ".", "m": "v"}
            color_map = {"p": "fuchsia", "0": "hotpink", "m": "darkmagenta"}
            fig = plt.figure(figsize=(11, 4))
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
                if L != "D" or e == "m":
                    continue  # negative reflectivity is too unstable
                D_wave = eJPmL
                S_wave = eJPmL[:-1] + "S"

                # get ratio and phase difference
                ratio = self.df[D_wave] / self.df[S_wave]
                ratio_err = ratio * np.sqrt(
                    np.square(self.df[D_wave + "_err"] / self.df[D_wave])
                    + np.square(self.df[S_wave + "_err"] / self.df[S_wave])
                )
                phase = self.df[self._phase_differences[(D_wave, S_wave)]].apply(
                    np.rad2deg
                )
                phase_err = self.df[
                    f"{self._phase_differences[(D_wave, S_wave)]}_err"
                ].apply(np.rad2deg)

                # label = rf"$1^{{+}}(S/D)_{{{m}}}^{{(+)}}$"
                label = convert_amp_name(eJPmL).replace("D", "(D/S)")

                ax_ratio.errorbar(
                    self._mass_bins,
                    ratio,
                    ratio_err,
                    self._bin_width / 2,
                    marker=marker_map[m],
                    color=color_map[m],
                    linestyle="",
                )
                ax_phase.errorbar(
                    self._mass_bins,
                    phase,
                    phase_err,
                    self._bin_width / 2,
                    marker=marker_map[m],
                    color=color_map[m],
                    linestyle="",
                )
                ax_phase.errorbar(
                    self._mass_bins,
                    -phase,
                    phase_err,
                    self._bin_width / 2,
                    marker=marker_map[m],
                    color=color_map[m],
                    linestyle="",
                )
                ax_corr.errorbar(
                    ratio,
                    phase,
                    phase_err,
                    ratio_err,
                    marker=marker_map[m],
                    color=color_map[m],
                    linestyle="",
                    label=label,
                )
                ax_corr.errorbar(
                    ratio,
                    -phase,
                    abs(phase_err),
                    ratio_err,
                    marker=marker_map[m],
                    color=color_map[m],
                    linestyle="",
                    label=label,
                )
            ax_ratio.set_ylabel("D/S ratio", loc="top")
            ax_ratio.set_ylim(0, 1)

            ax_phase.set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right")
            ax_phase.set_ylabel(r"D-S phase difference $(^\circ)$", loc="top")
            ax_phase.set_yticks(np.linspace(-180.0, 180.0, 5))
            ax_phase.set_ylim(-180.0, 180.0)

            ax_corr.set_xlabel("D/S ratio", loc="right")
            ax_corr.set_ylabel(r"D-S phase difference $(^\circ)$", loc="top")
            ax_corr.set_xlim(0, 1)
            ax_corr.set_ylim(-180.0, 180)
            ax_corr.set_yticks(np.linspace(-180.0, 180.0, 5))

            subfigs[1].subplots_adjust(right=0.65)

            handles, labels = ax_corr.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            fig.legend(
                by_label.values(),
                by_label.keys(),
                fontsize=9,
                loc="center right",
                bbox_to_anchor=(1, 0.5),
            )
            plt.show()
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
