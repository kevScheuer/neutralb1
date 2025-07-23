"""Collection of tools useful for analyzing PWA fit results

Primary class is the Plotter, which ideally handles most standard plots of interest when
analyzing fit results
"""

import cmath
import warnings
from typing import Dict, List, Literal, Optional, cast

import joypy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from cycler import cycler
from matplotlib.backends.backend_pdf import PdfPages


class Plotter:
    def __init__(
        self,
        df: pd.DataFrame,
        data_df: pd.DataFrame,
        bootstrap_df: "Optional[pd.DataFrame]" = None,
        truth_df: "Optional[pd.DataFrame]" = None,
    ) -> None:
        """Initialize object with pandas dataframes

        TODO: The truth phase differences between generated and non-generated waves
            is non-zero, but rather should be a flat line at the generated waves phase
            value. This will need to be calculated "manually" using the re/im parts.
            This can be handled in reset_truth_phases by checking if the phase is flat
        TODO: Implement the normality test from moment_projection here. It should use
            some form of shapiro wilkes for amp results (maybe log transformed?) and
            well defined circular data normality test. At the end of it, it should
            simply be a flag to look at distributions that are "suspected" non-normal
        Args:
            df (pd.DataFrame): FitResults from AmpTools
            data_df (pd.DataFrame): raw data points that AmpTools is fitting to
            bootstrap_df (pd.DataFrame, optional): all bootstrapped fit results.
                Primary use is for plotting the distribution of each parameter and using
                its width to estimate the error of the nominal fit result in df.
            truth_df (pd.DataFrame, optional): only used when df is a from fit to
                generated signal MC. Contains the "true" parameters that were used to
                generate the MC sample. When non-empty, many plots will display a
                solid line for the truth information
        """
        # --ERROR HANDLING--

        # # assign dataframes to public variables, copying to avoid modifying originals
        # self.fit_df = df.copy(deep=True)
        # self.data_df = data_df.copy(deep=True)
        # self.bootstrap_df = (
        #     bootstrap_df.copy(deep=True) if bootstrap_df is not None else None
        # )
        # self.truth_df = truth_df.copy(deep=True) if truth_df is not None else None

        # # Validate data
        # if self.fit_df.empty:
        #     raise ValueError("df is empty, cannot plot without fit results.")
        # if self.data_df.empty:
        #     raise ValueError("data_df is empty, cannot plot without data points.")
        # if self.bootstrap_df is not None and self.bootstrap_df.empty:
        #     raise ValueError("bootstrap_df was passed but is empty")
        # if self.truth_df is not None and self.truth_df.empty:
        #     raise ValueError("truth_df was passed but is empty")

        # if self.fit_df.shape[0] != self.data_df.shape[0]:
        #     raise ValueError(
        #         "df and data_df must have same number of rows!\n"
        #         f"df: {self.fit_df.shape[0]},\n"
        #         f" data_df: {self.data_df.shape[0]}"
        #     )
        # if self.bootstrap_df is not None and (
        #     len(self.bootstrap_df["file"].str.rsplit("/", n=1).str[0].unique())
        #     != self.fit_df.shape[0]
        # ):
        #     unique_directories = len(
        #         self.bootstrap_df["file"].str.rsplit("/", n=1).str[0].unique()
        #     )
        #     raise ValueError(
        #         "There must be a set of bootstrap fits for each fit result entry\n"
        #         f"bootstrap_df: {unique_directories},\n"
        #         f" df: {self.fit_df.shape[0]}"
        #     )
        # if self.truth_df is not None and (
        #     self.fit_df.shape[0] != self.truth_df.shape[0]
        # ):
        #     raise ValueError(
        #         "Fit and truth dataframes must have same number of rows!\n"
        #         f"df: {self.fit_df.shape[0]},"
        #         f" truth_df: {self.truth_df.shape[0]}"
        #     )

        # # load in matplotlib style for this class
        # plt.style.use(
        #     "/w/halld-scshelf2101/kscheuer/neutralb1/analysis/scripts/"
        #     "pwa_plotter.mplstyle"
        # )

        # # warn user of incomplete fits
        # if any(status != 3 for status in self.fit_df["eMatrixStatus"]):
        #     bad_files = self.fit_df[self.fit_df["eMatrixStatus"] != 3]["file"].to_list()
        #     warnings.warn(
        #         f"The following files contain fit results whose covariance matrix is"
        #         f" not full and accurate:\n{"\n".join(bad_files)}",
        #         UserWarning,
        #     )

        # public attributes
        # self.coherent_sums = get_coherent_sums(self.fit_df)
        # self.phase_differences = get_phase_differences(self.fit_df)

        # # private attributes
        # self._mass_bins = self.data_df["m_center"]
        # self._bin_width = (self.data_df["m_high"] - self.data_df["m_low"])[0]

        # --DATA PREPARATION--

        # # wrap phases from -pi to pi range
        # wrap_phases(self.fit_df)
        # if self.bootstrap_df is not None:
        #     wrap_phases(self.bootstrap_df)

        # if self.bootstrap_df is not None:
        #     # add a column for common directories between files for the bootstrap df, to
        #     # allow for easier chunking of bootstrap samples later
        #     # Create the directory column separately to avoid fragmentation, then concat
        #     directory_col = self.bootstrap_df["file"].str.rsplit("/", n=1).str[0]
        #     self.bootstrap_df = pd.concat(
        #         [self.bootstrap_df, directory_col.rename("directory")], axis=1
        #    )

        # Truth df needs special handling to allow for 1-1 comparison with fit results
        # if self.truth_df is not None:
        #     # phases in truth fits are wrong, so must be "manually" set to the correct
        #     # values. Then we can wrap them
        #     self._reset_truth_phases(self.truth_df, self._mass_bins)
        #     wrap_phases(self.truth_df)

        # set non-existent columns to zero. Needed since the fit dataframe may use
        # a different waveset from the truth fit
        # TODO: Move this to a preprocessing method
        # cols_to_add = []
        # for sum_type in self.coherent_sums:
        #     for sum in self.coherent_sums[sum_type]:
        #         if sum not in self.truth_df.columns:
        #             cols_to_add.append(
        #                 pd.DataFrame(
        #                     np.zeros_like(self.fit_df[sum]),
        #                     index=self.truth_df.index,
        #                     columns=[sum],
        #                 )
        #             )
        #             cols_to_add.append(
        #                 pd.DataFrame(
        #                     np.zeros_like(self.fit_df[sum]),
        #                     index=self.truth_df.index,
        #                     columns=[f"{sum}_err"],
        #                 )
        #             )
        # for phase_dif in set(self.phase_differences.values()):
        #     # check if the reverse ordering of the phase difference is used, and
        #     # rename it to match the fit dataframe's ordering if so
        #     reverse_phase_dif = "".join(phase_dif.partition("_")[::-1])
        #     if reverse_phase_dif in self.truth_df.columns:
        #         self.truth_df.rename(
        #             columns={
        #                 reverse_phase_dif: phase_dif,
        #                 f"{reverse_phase_dif}_err": f"{phase_dif}_err",
        #             },
        #             inplace=True,
        #         )
        #     # if phase difference not found, add it as columns of 0's
        #     if phase_dif not in self.truth_df.columns:
        #         cols_to_add.append(
        #             pd.DataFrame(
        #                 np.zeros_like(self.fit_df[phase_dif]),
        #                 index=self.truth_df.index,
        #                 columns=[phase_dif],
        #             )
        #         )
        #         cols_to_add.append(
        #             pd.DataFrame(
        #                 np.zeros_like(self.fit_df[phase_dif]),
        #                 index=self.truth_df.index,
        #                 columns=[f"{phase_dif}_err"],
        #             )
        #         )

        # concat_list = [self.truth_df]
        # concat_list.extend(cols_to_add)
        # self.truth_df = pd.concat(concat_list, axis=1)

        # use file names in dataframes to link fit indices to them for easier indexing
        # self._link_fits_to_dataframes()

        pass

    def jp(self, data_label: str = "Total (GlueX Phase-I)") -> None:
        """Plot each JP contribution to the total intensity as a function of mass"""

        # a map ensures consistency no matter the # of jps
        colors = matplotlib.colormaps["Dark2"].colors  # type: ignore
        jp_map = {
            "Bkgd": {"color": colors[0], "marker": "."},
            "0m": {"color": colors[1], "marker": "x"},
            "1p": {"color": colors[2], "marker": "o"},
            "1m": {"color": colors[3], "marker": "s"},
            "2p": {"color": colors[4], "marker": "P"},
            "2m": {"color": colors[5], "marker": "*"},
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
            self.fit_df["detected_events"],
            width=self._bin_width,
            color="0.1",
            alpha=0.15,
            label="Fit Result",
        )
        ax.errorbar(
            self._mass_bins,
            self.fit_df["detected_events"],
            self.fit_df["detected_events_err"],
            fmt=",",
            color="0.1",
            alpha=0.2,
            markersize=0,
        )

        # plot each jp contribution
        for jp, props in jp_map.items():
            if jp in self.fit_df.columns:
                l = convert_amp_name(jp) if jp != "Bkgd" else "Iso. Background"
                ax.errorbar(
                    self._mass_bins,
                    self.fit_df[jp],
                    self.fit_df[f"{jp}_err"],
                    self._bin_width / 2,
                    label=l,
                    linestyle="",
                    **props,
                )

        # plot jp contributions for truth df
        if self.truth_df is not None:
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

    def intensities(self, fractional: bool = False, sharey: bool = False) -> None:
        """Plot all the amplitude intensities in a grid format

        Since the matrix plot is generally too small to see bin to bin features, this
        method plots every amplitude's intensity contribution in a grid format.
        Columns = m-projections, rows = JPL combinations. Reflectivities are plotted
        together on each subplot

        Args:
            fractional (bool, optional): Scales all values by dividing them by the
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

        pos_plot = None
        neg_plot = None

        total = self.fit_df["detected_events"] if fractional else 1
        total_err = self.fit_df["detected_events_err"] if fractional else 0
        truth_total = 1

        if self.truth_df is not None and fractional:
            truth_total = self.truth_df["detected_events"]

        # iterate through JPL (sorted like S, P, D, F wave) and sorted m-projections
        for row, jpl in enumerate(jpl_values):
            for col, m in enumerate(m_ints):
                JPmL = f"{jpl[0:2]}{int_to_char[m]}{jpl[-1]}"

                # force sci notation so large ticklabels don't overlap with neighboring
                # plots if not plotting fit fraction
                if not fractional:
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
                    neg_refl = self.fit_df["m" + JPmL] / total
                    neg_refl_err = neg_refl * np.sqrt(
                        np.square(self.fit_df[f"m{JPmL}_err"] / self.fit_df["m" + JPmL])
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
                    if self.truth_df is not None:
                        axs[row, col].plot(
                            self._mass_bins,
                            self.truth_df["m" + JPmL] / truth_total,
                            linestyle="-",
                            marker="",
                            color="blue",
                        )
                # plot the positive reflectivity contribution
                if "p" + JPmL in self.coherent_sums["eJPmL"]:
                    pos_refl = self.fit_df["p" + JPmL] / total
                    pos_refl_err = pos_refl * np.sqrt(
                        np.square(self.fit_df[f"p{JPmL}_err"] / self.fit_df["p" + JPmL])
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
                    if self.truth_df is not None:
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
        legend_handles = []
        if pos_plot is not None:
            legend_handles.append(pos_plot)
        if neg_plot is not None:
            legend_handles.append(neg_plot)
        if legend_handles:
            fig.legend(handles=legend_handles, loc="upper right")
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
            self.fit_df[phase_dif],
            self.fit_df[f"{phase_dif}_err"].abs(),
            self._bin_width / 2,
            label=convert_amp_name(phase_dif),
            linestyle="",
            marker=".",
            color=color,
        )
        # plot the negative as well (natural sign ambiguity in the model)
        ax.errorbar(
            self._mass_bins,
            -self.fit_df[phase_dif],
            self.fit_df[f"{phase_dif}_err"].abs(),
            self._bin_width / 2,
            linestyle="",
            marker=".",
            color=color,
        )

        if self.truth_df is not None:
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
            self.fit_df[amp1],
            self.fit_df[f"{amp1}_err"],
            self._bin_width / 2,
            "o",
            color=color,
            label=convert_amp_name(amp1),
        )
        axs[0].errorbar(
            self._mass_bins,
            self.fit_df[amp2],
            self.fit_df[f"{amp2}_err"],
            self._bin_width / 2,
            "s",
            color=color,
            label=convert_amp_name(amp2),
        )

        if self.truth_df is not None:
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
            self.fit_df[phase_dif],
            self.fit_df[f"{phase_dif}_err"].abs(),
            self._bin_width / 2,
            linestyle="",
            marker=".",
            color=color,
        )
        axs[1].errorbar(
            self._mass_bins,
            -self.fit_df[phase_dif],
            self.fit_df[f"{phase_dif}_err"].abs(),
            self._bin_width / 2,
            linestyle="",
            marker=".",
            color=color,
        )
        if self.truth_df is not None:
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
        """Scatterplot-like matrix of subplots, for intensities and phase differences

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

        pos_plot = None
        neg_plot = None

        # only the top left corner plot will have ticks for the data scale
        max_diag = max([self.fit_df[x].max() for x in self.coherent_sums["eJPmL"]])
        data_y_ticks = np.linspace(0, max_diag, 4)
        axs[0, 0].set_yticks(data_y_ticks, [human_format(num) for num in data_y_ticks])

        for i, JPmL in enumerate(self.coherent_sums["JPmL"]):
            for j, JPmL_dif in enumerate(self.coherent_sums["JPmL"]):
                # change tick label sizes for all plots and the max
                axs[i, j].tick_params("both", labelsize=6)

                # write mass ticks for only the bottom row
                if i == len(self.coherent_sums["JPmL"]) - 1:
                    axs[i, j].xaxis.set_tick_params(labelbottom=True)

                # PLOT DIAGONALS
                if i == j:
                    axs[i, j].set_ylim(top=max_diag)
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
                        self.fit_df[f"m{JPmL}"],
                        self.fit_df[f"m{JPmL}_err"],
                        self._bin_width / 2,
                        "s",
                        color="blue",
                        markersize=2,
                        label=r"$\varepsilon=-1$",
                    )
                    pos_plot = axs[i, j].errorbar(
                        self._mass_bins,
                        self.fit_df[f"p{JPmL}"],
                        self.fit_df[f"p{JPmL}_err"],
                        self._bin_width / 2,
                        "o",
                        color="red",
                        markersize=2,
                        label=r"$\varepsilon=+1$",
                    )
                    if self.truth_df is not None:
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
                        self.fit_df[phase_dif].abs(),
                        linestyle="",
                        color="red",
                        marker=".",
                        ms=1,
                    )
                    axs[i, j].fill_between(
                        self._mass_bins,
                        self.fit_df[phase_dif].abs()
                        - self.fit_df[f"{phase_dif}_err"].abs(),
                        self.fit_df[phase_dif].abs()
                        + self.fit_df[f"{phase_dif}_err"].abs(),
                        color="red",
                        alpha=0.2,
                    )
                    axs[i, j].plot(
                        self._mass_bins,
                        -self.fit_df[phase_dif].abs(),
                        linestyle="",
                        color="red",
                        marker=".",
                        ms=1,
                    )
                    axs[i, j].fill_between(
                        self._mass_bins,
                        -self.fit_df[phase_dif].abs()
                        - self.fit_df[f"{phase_dif}_err"].abs(),
                        -self.fit_df[phase_dif].abs()
                        + self.fit_df[f"{phase_dif}_err"].abs(),
                        color="red",
                        alpha=0.2,
                    )

                    if self.truth_df is not None:
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
                        self.fit_df[phase_dif].abs(),
                        linestyle="",
                        color="blue",
                        marker=".",
                        ms=1,
                    )
                    axs[i, j].fill_between(
                        self._mass_bins,
                        self.fit_df[phase_dif].abs()
                        - self.fit_df[f"{phase_dif}_err"].abs(),
                        self.fit_df[phase_dif].abs()
                        + self.fit_df[f"{phase_dif}_err"].abs(),
                        color="blue",
                        alpha=0.2,
                    )
                    axs[i, j].plot(
                        self._mass_bins,
                        -self.fit_df[phase_dif].abs(),
                        linestyle="",
                        color="blue",
                        marker=".",
                        ms=1,
                    )
                    axs[i, j].fill_between(
                        self._mass_bins,
                        -self.fit_df[phase_dif].abs()
                        - self.fit_df[f"{phase_dif}_err"].abs(),
                        -self.fit_df[phase_dif].abs()
                        + self.fit_df[f"{phase_dif}_err"].abs(),
                        color="blue",
                        alpha=0.2,
                    )
                    if self.truth_df is not None:
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
        legend_handles = []
        if pos_plot is not None:
            legend_handles.append(pos_plot)
        if neg_plot is not None:
            legend_handles.append(neg_plot)
        if legend_handles:
            fig.legend(handles=legend_handles, fontsize=12, loc="upper right")

        plt.show()
        pass

    def ds_ratio(self) -> None:
        """Plot the ratio and phase between the D and S waves

        The plots change depending on whether the model had a defined D/S ratio. If it
        is defined then the covariance of the ratio and phase parameters are plotted,
        but if not, then its the correlation. Due to plotting style, plots have fixed
        sizes.

        TODO: Add a raise if D or S wave are not present, and the reflectivities
        """

        if "dsratio" in self.fit_df.columns:
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
                self.fit_df["dsratio"],
                self.fit_df["dsratio_err"],
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
                self.fit_df["dphase"],
                self.fit_df["dphase_err"].abs(),
                self._bin_width / 2,
                marker=".",
                linestyle="",
                color="gray",
            )
            ax_phase.errorbar(
                self._mass_bins,
                -self.fit_df["dphase"],
                self.fit_df["dphase_err"].abs(),
                self._bin_width / 2,
                marker=".",
                linestyle="",
                color="gray",
            )
            # CORRELATION PLOT
            ax_corr.plot(
                self._mass_bins,
                (
                    self.fit_df["cov_dsratio_dphase"]
                    / (self.fit_df["dsratio_err"] * self.fit_df["dphase_err"])
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
                ratio = self.fit_df[D_wave].apply(np.sqrt) / self.fit_df[S_wave].apply(
                    np.sqrt
                )
                ratio_err = ratio * np.sqrt(
                    np.square(
                        self.fit_df[f"{D_wave}_err"].apply(np.sqrt)
                        / self.fit_df[D_wave].apply(np.sqrt)
                    )
                    + np.square(
                        self.fit_df[f"{S_wave}_err"].apply(np.sqrt)
                        / self.fit_df[S_wave].apply(np.sqrt)
                    )
                )
                phase = self.fit_df[self.phase_differences[(D_wave, S_wave)]]
                phase_err = self.fit_df[
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

    def correlation_matrix(
        self,
        columns: List[str],
        pdf_path: str,
        method: "Literal['pearson', 'kendall', 'spearman']" = "pearson",
        directories: Optional[List[str]] = None,
    ) -> None:
        """
        Plot correlation matrix of requested columns, for every bootstrap sample

        Args:
            columns (List[str]): bootstrap df columns to calculate the correlation
                matrix for
            pdf_path (str): Save all matrices to this PDF file.
            method (str,optional): Method of correlation. Defaults to "pearson", but can
                also be kendall or spearman.
            directories (List[str], optional): Bootstrap samples are grouped by
                directory. If provided, only the directories (groups) requested will be
                plotted. Default to None, so that all samples are plotted

        Raises:
            ValueError: if bootstrap df was not defined, since its optional upon init
            ValueError: requested columns not in the bootstrap dataframe
            ValueError: file(s) not found in dataframe
        """
        if self.bootstrap_df is None:
            raise ValueError("Bootstrap df was not defined on instantiation")

        missing_columns = [
            col for col in columns if col not in self.bootstrap_df.columns
        ]
        if missing_columns:
            raise ValueError(
                f"The following columns are not in bootstrap_df: {missing_columns}"
            )

        if directories is not None:
            if not self.bootstrap_df["directory"].isin(directories).any():
                raise ValueError(
                    "Not all requested directories were found in the bootstrap df"
                )

        # Ensure directories is always a list
        if directories is None:
            directories = self.bootstrap_df["directory"].unique().tolist()
        elif not isinstance(directories, list):
            directories = list(directories)
        directories = cast(List[str], directories)

        pdf = PdfPages(pdf_path)

        for d in directories:
            selected_df = self.bootstrap_df[self.bootstrap_df["directory"] == d]
            corr_matrix = selected_df[columns].corr(method=method).round(2)
            renamed_titles = {}
            for col in columns:
                if any(
                    col in sublist for sublist in self.coherent_sums.values()
                ) or col in set(self.phase_differences.values()):
                    renamed_titles[col] = convert_amp_name(col)
            corr_matrix.rename(
                index=renamed_titles, columns=renamed_titles, inplace=True
            )

            plt.figure(figsize=(10, 8))
            sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", vmin=-1, vmax=1)
            plt.title(f"Correlation Matrix for {d}")
            plt.yticks(rotation=0)
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            pdf.savefig()
            plt.close()

        print("Saved correlation matrices to", pdf_path)
        pdf.close()

        pass

    def pull_distribution(
        self,
        columns: str | List[str],
        y_limits: Optional[List[int]] = None,
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
        if self.truth_df is None:
            raise ValueError("Truth information is required for pull distribution")

        if isinstance(columns, str):
            columns = [columns]

        # Create a color cycle using the Accent colormap
        color_cycle = cycler(color=matplotlib.cm.Accent.colors)  # type: ignore
        fig, ax = plt.subplots()
        ax.set_prop_cycle(color_cycle)

        norm = 1
        if scale:
            # get normalized max values in the columns to determine marker size
            norm = max([self.fit_df[col].max() for col in columns])

        max_pull = 0
        for column in columns:
            # calculate pull and store the max value
            pull = (self.fit_df[column] - self.truth_df[column]) / self.fit_df[
                f"{column}_err"
            ]
            max_pull = max(max_pull, pull.abs().max())

            marker_sizes = scale_factor * (self.fit_df[column] / norm) if scale else 1

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
        ax.set_ylabel(
            r"Pull $\frac{\text{fit} - \text{truth}}{\sigma_{\text{fit}}}$", loc="top"
        )
        ax.set_xlabel(r"$\omega\pi^0$ inv. mass $(GeV)$", loc="right")
        ax.legend(**legend_kwargs)

        plt.minorticks_on()
        plt.show()
        pass

    def diagnose_bootstraps(
        self,
        columns: Optional[list] = None,
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

    # TODO: move to plotter, but consider if any components can be in statistics
    def bootstrap_matrix(
        self, fits: List[str] | str, columns: List[str], **kwargs
    ) -> None:
        """Scatterplot matrix of the bootstrap fit results with nominal fit overlaid

        Bootstrap fits can reveal biases or correlations between fit parameters, so this
        plot provides a full view of this for select fit parameters, chosen with
        "columns". Multiple files can be viewed at once by passing in a list to "fits".

        The matrix is of the form: histogram + its kde on the diagonal, scatter plots on
        the lower triangle, and 2D kde's on the upper triangle. Every plot has a thick
        green line overlaid to indicate the nominal fit result values. Plots framed by a
        solid black line indicate a strong correlation between those variables (>0.7).
        All amplitudes are plotted in fit fractions.

        TODO: use mass bin values from data_df to label plots

        Args:
            fits (List[str]|str): fit files, and their associated bootstrap files
            columns (List[str]): fit parameters to plot
            **kwargs: specifically for the seaborn pairplot object

        Raises:
            ValueError: if bootstrap df was not defined, since its optional upon init
        """

        if self.bootstrap_df is None:
            raise ValueError(
                "Object must be initialized with bootstrap DataFrame to utilize this"
                " method"
            )
        if not fits:
            raise ValueError("A fit file must be chosen to be plotted")
        if not columns:
            raise ValueError("No fit parameters were chosen to be plotted")
        for col in columns:
            if col not in self.bootstrap_df.columns:
                raise KeyError(f"Column {col} not found in the bootstrap DataFrame")

        if isinstance(fits, str):
            fits = [fits]

        # ensure every fit has an associated bootstrap sample, and exists in the df
        for f in fits:
            if f not in self.fit_df["file"].to_list():
                raise FileNotFoundError(
                    f"Requested fit file {f} could not be found in the fit result"
                    " DataFrame"
                )
            if (
                not self.fit_df.index[self.fit_df["file"] == f]
                .isin(self.bootstrap_df["fit_index"])
                .any()
            ):
                raise FileNotFoundError(
                    f"Requested fit file {f} does not have any bootstrap results"
                    " associated with it. Check the bootstrap csv to ensure its file"
                    f" paths are in a '/bootstrap/' subdirectory of the requested file"
                )

        # in case only some error columns are provided, remove them and re-add them
        # separately
        for c in columns:
            columns.remove(c) if "_err" in c else None
        err_columns = [f"{c}_err" for c in columns]
        selected_columns = columns + err_columns + ["detected_events"]

        # select the samples
        cut_df = self.fit_df[self.fit_df["file"].isin(fits)][selected_columns].copy()
        cut_bootstrap_df = self.bootstrap_df[
            self.bootstrap_df["fit_index"].isin(cut_df.index)
        ][selected_columns + ["fit_index"]].copy()
        if self.truth_df is not None:
            cut_truth_df = self.truth_df[self.truth_df["fit_index"].isin(cut_df.index)][
                selected_columns + ["fit_index"]
            ].copy()
        else:
            cut_truth_df = None

        # we want plotted amplitudes to be their fit fraction for easier visuals
        plotter_df = cut_bootstrap_df.copy()
        for col in plotter_df:
            if any(col in sublist for sublist in self.coherent_sums.values()):
                plotter_df[col] = plotter_df[col].div(
                    cut_bootstrap_df["detected_events"], axis="index"
                )

        # scale truth phases to be the same sign as the nominal fit phases. Due to phase
        # ambiguity in model, the nominal fit can obtain +/-|truth_phase|.
        if cut_truth_df is not None:
            # Only process phase differences that are present in both cut_truth_df and cut_df
            available_phases = [
                phase
                for phase in set(self.phase_differences.values())
                if phase in cut_truth_df.columns and phase in cut_df.columns
            ]
            for phase in available_phases:
                cut_truth_df[phase] = cut_truth_df[phase].where(
                    np.sign(cut_truth_df[phase]) == np.sign(cut_df[phase]),
                    -1.0 * cut_truth_df[phase],
                )

        # get default color palette for hue plotting
        palette = sns.color_palette(n_colors=plotter_df["fit_index"].nunique())

        pg = sns.PairGrid(
            plotter_df[columns + ["fit_index"]],
            hue="fit_index",
            palette=palette,
            **kwargs,
        )
        # overlay a kde on the hist to compare nominal fit line to kde peak
        # if many bins are plotted, remove the histogram
        if plotter_df["fit_index"].nunique() < 4:
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
                    if cut_truth_df is not None:
                        x_truth_scale = cut_truth_df["detected_events"]
                if any(row_label in sublist for sublist in self.coherent_sums.values()):
                    y_scaling = cut_df["detected_events"]
                    if cut_truth_df is not None:
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
                if cut_truth_df is not None and row == col:
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
                    if cut_truth_df is not None:
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
                    for fit_index in cut_df.index:
                        corr = (
                            cut_bootstrap_df[
                                cut_bootstrap_df["fit_index"] == fit_index
                            ][[col_label, row_label]]
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
                    for fit_index, error in zip(
                        cut_df[f"{col_label}_err"].index, cut_df[f"{col_label}_err"]
                    ):
                        series = cut_bootstrap_df[
                            cut_bootstrap_df["fit_index"] == fit_index
                        ][col_label]

                        # phases differences (circular data) must have standard dev
                        # specially calculated
                        if col_label in set(self.phase_differences.values()):
                            # convert to rad for scipy function
                            radian_phases = np.deg2rad(series)
                            stdev = scipy.stats.circstd(
                                radian_phases, low=-np.pi, high=np.pi
                            )
                            # now back to degrees
                            stdev = np.rad2deg(stdev)
                        else:
                            stdev = np.std(series)
                        ratios.append(stdev / error)
                    pg.axes[row, col].text(
                        0.6,
                        0.9,
                        (
                            r"$\frac{{\sigma_{\text{bootstrap}}}}"
                            r"{{\sigma_{\text{MINUIT}}}}$"
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
        leg = pg.legend
        for ax in pg.axes.flat:
            if ax.get_legend() is not None:
                leg = ax.get_legend()
                break
        if leg is not None:
            leg.set_title("fit label")

        plt.show()
        pass

    def joyplot(
        self,
        columns: List[str],
        fits: Optional[List[str]] = None,
        colormap: Optional[str] = None,
        is_sparse_labels=True,
        overlap: int = 2,
        **kwargs,
    ) -> None:
        """Plot a joyplot (or ridgeplot) of the bootstrap fit results

        The joyplot is a series of kde's that are stacked on top of each other, with
        each kde representing the distribution of a fit parameter in a file. Note
        that the x-axis is shared, so its recommended to only plot like-parameters and
        not to mix columns like phases and intensities. If the truth values are present,
        overlay them as a vertical line on each axis.

        Args:
            columns (List[str]): bootstrap df columns to plot
            fits (List[str], optional): fit files (typically mass bins) to plot on the
                y-axis. Defaults to None, using all files.
            colormap (str, optional): name of matplotlib colormap. Recommend that only
                "qualitative" colormaps be used. Defaults to None, and uses the
                "Accent" map.
            is_sparse_labels (bool, optional): when True, only bins that start at values
                divisible by 10 are labelled. Defaults to True.
            overlap (float, optional): amount of overlap between the kde's. Usually the
                most important parameter for making plots look good. Defaults to 2.0
            **kwargs: additional arguments to pass to the joyplot function
        Raises:
            ValueError: if bootstrap df was not defined, since its optional upon init
            KeyError: if requested columns are not in the dataframes
        """

        if self.bootstrap_df is None:
            raise ValueError("Bootstrap df was not defined on instantiation")
        # check if requested columns are in the dataframes
        missing_columns = []
        for col in columns:
            if col not in self.bootstrap_df.columns or col not in self.fit_df.columns:
                missing_columns.append(col)
            if self.truth_df is not None and col not in self.truth_df.columns:
                missing_columns.append(col)
        if missing_columns:
            raise KeyError(
                f"The following columns were found to be missing from one or more"
                f" dataframes: {missing_columns}"
            )
        # Check if colormap is a valid matplotlib colormap name and convert to list of
        # colors, or set default value
        if isinstance(colormap, str):
            if colormap in matplotlib.colormaps:
                color_list = list(matplotlib.colormaps[colormap].colors)  # type: ignore
            else:
                raise KeyError(f"'{colormap}' is not a valid matplotlib colormap name.")
        elif colormap is not None:
            raise TypeError("Expected a string for colormap name.")
        else:
            color_list = list(matplotlib.colormaps["Accent"].colors)  # type: ignore

        # Ensure colormap matches the length of columns
        if len(color_list) < len(columns):  # repeat colormap if too short
            color_list = (color_list * (len(columns) // len(color_list) + 1))[
                : len(columns)
            ]
        else:
            color_list = color_list[: len(columns)]

        # use all fits if none are specified
        if not fits:
            fit_indices = self.fit_df.index
        else:
            fit_indices = self.fit_df[self.fit_df["file"].isin(fits)].index

        # select fits from dataframes and determine min/max x values
        cut_data_df = self.data_df[self.data_df["fit_index"].isin(fit_indices)].copy()
        cut_bootstrap_df = self.bootstrap_df[
            self.bootstrap_df["fit_index"].isin(fit_indices)
        ].copy()
        x_max = cut_bootstrap_df[columns].max().max()
        x_min = cut_bootstrap_df[columns].min().min()
        if self.truth_df is not None:
            cut_truth_df = self.truth_df[
                self.truth_df["fit_index"].isin(fit_indices)
            ].copy()
            x_max = max(x_max, cut_truth_df[columns].max().max())
            x_min = min(x_min, cut_truth_df[columns].min().min())
        else:
            cut_truth_df = None

        # force x_min to be 0 if all values are coherent sums
        if all(
            any(col in sublist for sublist in self.coherent_sums.values())
            for col in columns
        ):
            x_min = 0.0
        # force x_min/x_max to be -180/180 if all values are phase differences
        elif all(col in set(self.phase_differences.values()) for col in columns):
            x_min, x_max = -180.0, 180.0
        else:
            warnings.warn(
                "Columns chosen have no recognized common x-axis limits, which may lead"
                "to misleading plots",
                UserWarning,
            )

        # setup y-axis labels to indicate bin ranges
        y_axis_labels = []
        for idx, (low_edge, high_edge) in enumerate(
            zip(cut_data_df["m_low"].astype(float), cut_data_df["m_high"].astype(float))
        ):
            label = f"{low_edge:.2f}-{high_edge:.2f}"
            # if not sparse labels, then simply add all labels
            if not is_sparse_labels:
                y_axis_labels.append(label)
                continue

            tolerance = 1e-5
            # if first or last bin then force add label
            if idx == 0 or idx == len(cut_data_df["m_low"]):
                y_axis_labels.append(label)
            # add label if low bin edge is a multiple of 10 (within tolerance)
            elif (
                abs(low_edge * 100 % 10) < tolerance
                or abs(low_edge * 100 % 10 - 10) < tolerance
            ):
                y_axis_labels.append(f"{label}")
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

        # avoid bug within joyplot that occurs when color_list is a single element list
        joy_color = color_list[0] if len(color_list) == 1 else color_list

        fig, axes = joypy.joyplot(
            cut_bootstrap_df,
            by="fit_index",
            column=columns,
            labels=y_axis_labels,
            range_style="own",
            x_range=[x_min, x_max],
            color=joy_color,
            overlap=overlap,
            **default_args,
        )

        # plot the truth value as a vertical line on each axis
        if cut_truth_df is not None:
            for fit_index, ax in zip(cut_bootstrap_df["fit_index"].unique(), axes):
                # get the max y value by getting the kde max (kde's are the
                # odd values in the ax.lines object)
                y_maxes = [
                    line.get_ydata().max()
                    for line in ax.lines
                    if ax.lines.index(line) % 2 != 0
                ]
                for col, line_color in zip(columns, color_list):
                    truth_value = cut_truth_df[cut_truth_df["fit_index"] == fit_index][
                        col
                    ].iloc[0]
                    ax.plot(
                        [truth_value, truth_value],
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

    # TODO: move to statistics
    def normality_test(
        self,
        columns: Optional[List[str]] = None,
        method: Optional[Literal["None", "normal", "log_normal"]] = None,
    ) -> None:
        """
        Shapiro-wilkes for parameters and amplitudes. Log-transform the amps near zero.
        Use other method for circular (phase) data.

        """
        if self.bootstrap_df is None:
            raise ValueError("Bootstrap df was not defined on instantiation")

        # PLAN: loop over the fit indices of the bootstrap df and calculate the test.
        # Those that don't pass, or whose bootstrap mean is too far from the fit
        # result's are plotted in big pdf. There should be indicators for the fit mean
        # value on the QQ plot, and whether it failed due to p value or mean

        fit_indices = self.bootstrap_df["fit_index"].unique()

        # add the amplitudes and phase differences to the columns to be tested
        columns = columns if columns is not None else []
        columns.extend(self.coherent_sums["eJPmL"])
        columns.extend(set(self.phase_differences.values()))

        failed_dict = {}  # to store {fit_index: {column : [mean, p-value]}}
        for fit_index in fit_indices:
            for col in columns:
                if col in set(self.phase_differences.values()):
                    # do circular data test. Take absolute value because we have no way
                    # to distinguish positive vs negative values
                    pass

                elif any(col in sublist for sublist in self.coherent_sums.values()):
                    # do normality test on amplitudes
                    # log-transform the amps near zero
                    # type: ignore or type annotation for pylance
                    amp_values = self.bootstrap_df.loc[
                        self.bootstrap_df["fit_index"] == fit_index
                    ][col]

                    # log-transform if col is near zero
                    # use rel distance to minimum test to determine if near zero
                    # TODO: finish this

                    # perform normality test
                    stat, p_value = scipy.stats.shapiro(amp_values)
                    if p_value < 0.05:
                        failed_dict.setdefault(fit_index, {})[col] = [
                            amp_values.mean(),
                            p_value,
                        ]

                else:
                    # do normality test as usual
                    values = self.bootstrap_df.loc[
                        self.bootstrap_df["fit_index"] == fit_index
                    ][col]

                    # perform normality test
                    stat, p_value = scipy.stats.shapiro(values)
                    if p_value < 0.05:
                        failed_dict.setdefault(fit_index, {})[col] = [
                            values.mean(),
                            p_value,
                        ]

        # PLAN: here we'll loop over the failed ones and plot probdists

        pass
