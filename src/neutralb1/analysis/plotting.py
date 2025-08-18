"""A collection of classes for quick plotting of PWA fit results

The FactorPlotter interfaces with the various plotter classes to provide a
collection of methods for plotting PWA fit results, including intensity,
phase, and diagnostic plots. All plotters inherit a common base class.
"""

# TODO: use covariance from bootstrap results for error propagation

from typing import Literal, Optional

import joypy
import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
from matplotlib import colormaps
from matplotlib.backends.backend_pdf import PdfPages

from neutralb1 import utils


class FactoryPlotter:
    """Factory class to create plotters based on the type of plot needed."""

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        random_df: Optional[pd.DataFrame] = None,
        bootstrap_df: Optional[pd.DataFrame] = None,
        truth_df: Optional[pd.DataFrame] = None,
    ) -> None:
        """Initialize the factory with common data and utilities."""
        self.fit_df = fit_df
        self.data_df = data_df
        self.random_df = random_df
        self.bootstrap_df = bootstrap_df
        self.truth_df = truth_df

        # Set the matplotlib style for consistent plotting
        WORKSPACE_DIR = utils.get_workspace_dir()
        plt.style.use(f"{WORKSPACE_DIR}/config/neutralb1.mplstyle")

    def intensity(self):
        return IntensityPlotter(
            self.fit_df, self.data_df, self.bootstrap_df, self.truth_df
        )

    def phase(self):
        return PhasePlotter(self.fit_df, self.data_df, self.bootstrap_df, self.truth_df)

    def diagnostic(self):
        return DiagnosticPlotter(
            self.fit_df, self.data_df, self.bootstrap_df, self.truth_df
        )

    def random(self):
        return RandomPlotter(self.fit_df, self.data_df, self.random_df)

    def bootstrap(self):
        return BootstrapPlotter(
            self.fit_df, self.data_df, self.bootstrap_df, self.truth_df
        )


class BasePWAPlotter:
    """Base class for PWA plotting"""

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        random_df: Optional[pd.DataFrame] = None,
        bootstrap_df: Optional[pd.DataFrame] = None,
        truth_df: Optional[pd.DataFrame] = None,
        channel: Optional[str] = r"$\omega\pi^0$",
    ) -> None:
        """Initialize base plotter with common data and utilities.

        Args:
            fit_df (pd.DataFrame): Nominal fit results DataFrame.
            data_df (pd.DataFrame): Contains the data used for the fit.
            bootstrap_df (pd.DataFrame, optional): Bootstrap results for each nominal
                fit. Defaults to None.
            truth_df (pd.DataFrame, optional): Contains the ground truth values for
                the fit. Defaults to None.
            channel (str, optional): The channel label to be used in plot axes. Defaults
                to r"$\\omega\\pi^0$".
        """

        # no need to copy DataFrames here, as they will not be modified
        self.fit_df = fit_df
        self.data_df = data_df
        self.random_df = random_df
        self.bootstrap_df = bootstrap_df
        self.truth_df = truth_df
        self.channel = channel

        # Common properties
        self._masses = self.data_df["m_center"]
        self._bin_width = (self.data_df["m_high"] - self.data_df["m_low"])[0]

        # Extract coherent sums and phase differences
        self.coherent_sums = utils.get_coherent_sums(self.fit_df)
        self.phase_differences = utils.get_phase_differences(self.fit_df)


class IntensityPlotter(BasePWAPlotter):
    """Handles all strictly intensity-related plotting methods."""

    def jp(self, data_label: str = "Total (GlueX Phase-I)") -> matplotlib.axes.Axes:
        """Plot each JP contribution to the total intensity as a function of mass.

        Args:
            data_label (str): Label for the data in the plot. Defaults to "Total
                (GlueX Phase-I)".
        Returns:
            matplotlib.axes.Axes: The axes object for further customization.
        """

        # property map for consistent plotting
        colors = colormaps["Dark2"].colors  # type: ignore
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
            x=self._masses,
            xerr=self._bin_width / 2,
            y=self.data_df["events"],
            yerr=self.data_df["events_err"],
            marker=".",
            color="black",
            label=data_label,
        )

        # Plot Fit Result as a grey histogram with its error bars
        ax.bar(
            self._masses,
            self.fit_df["detected_events"],
            width=self._bin_width,
            color="0.1",
            alpha=0.15,
            label="Fit Result",
        )
        ax.errorbar(
            x=self._masses,
            y=self.fit_df["detected_events"],
            yerr=self.fit_df["detected_events_err"],
            marker=",",
            markersize=0,
            color="0.1",
            alpha=0.2,
        )

        # plot each jp contribution
        for jp, props in jp_map.items():
            if jp in self.fit_df.columns:
                l = (
                    utils.convert_amp_name(jp)
                    if "Bkgd" not in jp
                    else "Iso. Background"
                )
                ax.errorbar(
                    x=self._masses,
                    xerr=self._bin_width / 2,
                    y=self.fit_df[jp],
                    yerr=self.fit_df[f"{jp}_err"],
                    linestyle="",
                    label=l,
                    **props,
                )

        # plot jp contributions for truth df as continuous lines
        if self.truth_df is not None:
            for jp, props in jp_map.items():
                if jp in self.truth_df.columns:
                    ax.plot(
                        self._masses,
                        self.truth_df[jp],
                        linestyle="-",
                        marker="",
                        color=props["color"],
                    )

        ax.set_xlabel(rf"{self.channel} inv. mass $(GeV)$", loc="right")
        ax.set_ylabel(f"Events / {self._bin_width:.3f} GeV", loc="top")
        ax.set_ylim(bottom=0.0)
        ax.legend()
        plt.minorticks_on()

        return ax

    def waves(self, fractional: bool = False, sharey: bool = False) -> np.ndarray:
        """Plot all the wave intensities in a grid format.

        Creates a grid plot where columns represent spin-projections (m) and rows
        represent JPL combinations. Reflectivities are plotted together on each
        subplot with distinct markers and colors.

        Args:
            fractional (bool): If True, scales all values by dividing by the total
                intensity in each bin. Defaults to False.
            sharey (bool): If True, each row shares a y-axis scale. Defaults to False.

        Returns:
            np.ndarray: Array of matplotlib.axes.Axes objects for further customization.

        Raises:
            ValueError: If no wave amplitudes are found in the fit results.
        """

        # Validate that we have wave data to plot
        if not self.coherent_sums["JPmL"]:
            raise ValueError("No wave amplitudes found in fit results.")

        # Helper dictionaries for formatting TODO: use utils functions
        int_to_char = {-2: "n", -1: "m", 0: "0", 1: "p", 2: "q"}
        pm_dict = {"m": "-", "p": "+", "0": "0", "n": "-", "q": "+"}

        # Sort values in increasing order of angular momenta and spin-projection
        jpl_values = sorted(
            self.coherent_sums["JPL"], key=lambda JPL: utils.char_to_int(JPL[-1])
        )
        m_values = sorted(
            {utils.char_to_int(JPmL[-2]) for JPmL in self.coherent_sums["JPmL"]}
        )

        # Create figure
        n_rows, n_cols = len(jpl_values), len(m_values)
        fig_width = max(10, n_cols * 3)
        fig_height = max(6, n_rows * 2.5)

        fig, axs = plt.subplots(
            n_rows,
            n_cols,
            sharex=True,
            sharey=sharey,
            figsize=(fig_width, fig_height),
            squeeze=False,  # Always return 2D array even for single subplot
        )

        # Prepare normalization data
        if fractional:
            total_intensity = self.fit_df["detected_events"]
            total_intensity_err = self.fit_df["detected_events_err"]
            truth_total = (
                self.truth_df["detected_events"] if self.truth_df is not None else None
            )
        else:
            total_intensity = pd.Series(
                np.ones_like(self.fit_df["detected_events"]),
                index=self.fit_df.index,
            )
            total_intensity_err = pd.Series(
                np.zeros_like(self.fit_df["detected_events_err"]),
                index=self.fit_df.index,
            )
            truth_total = None

        # Track plot handles for legend
        pos_plot_handle = None
        neg_plot_handle = None

        # Plot each wave combination
        for row, jpl in enumerate(jpl_values):
            for col, m in enumerate(m_values):
                ax = axs[row, col]
                JPmL = f"{jpl[:2]}{int_to_char[m]}{jpl[-1]}"
                J, P, L = JPmL[0], JPmL[1], JPmL[-1]

                # Configure y-axis formatting for non-fractional plots
                if not fractional:
                    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

                # Set subplot titles and labels
                if row == 0:
                    ax.set_title(f"m = {m}", fontsize=14, pad=10)
                if col == 0:
                    wave_label = rf"${J}^{{{pm_dict[P]}}}{L}$"
                    ax.set_ylabel(wave_label, fontsize=14)

                # Plot negative reflectivity contribution
                neg_amp_name = f"m{JPmL}"
                if neg_amp_name in self.coherent_sums["eJPmL"]:
                    y_values = self.fit_df[neg_amp_name] / total_intensity
                    y_errors = utils.propagate_product_error(
                        self.fit_df[neg_amp_name],
                        self.fit_df[f"{neg_amp_name}_err"],
                        total_intensity,
                        total_intensity_err,
                    )

                    neg_plot_handle = ax.errorbar(
                        x=self._masses,
                        xerr=self._bin_width / 2,
                        y=y_values,
                        yerr=y_errors,
                        marker="v",
                        linestyle="",
                        markersize=6,
                        color="blue",
                        alpha=0.7,
                        label=r"$\varepsilon = -1$",
                        capsize=2,
                    )

                    # Plot truth data if available
                    if self.truth_df is not None:
                        truth_y = (
                            self.truth_df[neg_amp_name] / truth_total
                            if truth_total is not None
                            else self.truth_df[neg_amp_name]
                        )
                        ax.plot(
                            self._masses,
                            truth_y,
                            linestyle="-",
                            linewidth=2,
                            color="blue",
                            alpha=0.8,
                        )

                # Plot positive reflectivity contribution
                pos_amp_name = f"p{JPmL}"
                if pos_amp_name in self.coherent_sums["eJPmL"]:
                    y_values = self.fit_df[pos_amp_name] / total_intensity
                    y_errors = utils.propagate_product_error(
                        self.fit_df[pos_amp_name],
                        self.fit_df[f"{pos_amp_name}_err"],
                        total_intensity,
                        total_intensity_err,
                    )

                    pos_plot_handle = ax.errorbar(
                        x=self._masses,
                        xerr=self._bin_width / 2,
                        y=y_values,
                        yerr=y_errors,
                        marker="^",
                        linestyle="",
                        markersize=6,
                        color="red",
                        alpha=0.7,
                        label=r"$\varepsilon = +1$",
                        capsize=2,
                    )

                    # Plot truth data if available
                    if self.truth_df is not None:
                        truth_y = (
                            self.truth_df[pos_amp_name] / truth_total
                            if truth_total is not None
                            else self.truth_df[pos_amp_name]
                        )
                        ax.plot(
                            self._masses,
                            truth_y,
                            linestyle="-",
                            linewidth=2,
                            color="red",
                            alpha=0.8,
                        )

                # Set y-limit to start from zero
                ax.set_ylim(bottom=0)

        # Configure figure-wide labels and legend
        y_label = (
            "Fit Fraction" if fractional else f"Events / {self._bin_width:.3f} GeV"
        )

        fig.supxlabel(rf"{self.channel} inv. mass (GeV)", fontsize=16)
        fig.supylabel(y_label, fontsize=16)

        # Add legend if we have plot handles
        legend_handles = []
        if pos_plot_handle is not None:
            legend_handles.append(pos_plot_handle)
        if neg_plot_handle is not None:
            legend_handles.append(neg_plot_handle)

        if legend_handles:
            fig.legend(
                handles=legend_handles,
                loc="upper right",
                bbox_to_anchor=(0.98, 0.98),
                fontsize=12,
            )

        plt.tight_layout()
        plt.minorticks_on()

        return axs


class PhasePlotter(BasePWAPlotter):
    """Handles all strictly phase-related plotting methods."""

    def phase(
        self, amp1: str, amp2: str, extend_range: bool = True
    ) -> matplotlib.axes.Axes:
        """Plot the phase difference between two amplitudes as a function of mass.

        Order of amplitudes does not matter

        Args:
            amp1 (str): Name of the first amplitude in eJPmL format.
            amp2 (str): Name of the second amplitude in eJPmL format.
            extend_range (bool): Extends the y-axis from [-pi, pi) to [-2pi, 2pi) bounds
                if True, and plots mirrored transparaent values. Useful for revealing
                phase motions near the boundaries. Defaults to True.
        Returns:
            matplotlib.axes.Axes: The axes object for further customization.
        """

        phase_dif = self.phase_differences[(amp1, amp2)]
        color = "red" if amp1[0] == "p" else "blue"

        fig, ax = plt.subplots()
        ax.errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=self.fit_df[phase_dif],
            yerr=self.fit_df[f"{phase_dif}_err"].abs(),
            linestyle="",
            marker=".",
            color=color,
            label=utils.convert_amp_name(phase_dif),
        )
        # plot the negative as well (natural sign ambiguity in the model)
        ax.errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=-self.fit_df[phase_dif],
            yerr=self.fit_df[f"{phase_dif}_err"].abs(),
            linestyle="",
            marker=".",
            color=color,
        )

        # plot the truth phase if available
        if self.truth_df is not None:
            ax.plot(
                self._masses,
                self.truth_df[phase_dif],
                linestyle="-",
                marker="",
                color=color,
            )

        if extend_range:
            # add or subtract ambiguous 360 degree shifts. Positive values get shifted
            # down by 2pi, and negative values up by 2pi.
            ax.errorbar(
                x=self._masses,
                xerr=self._bin_width / 2,
                y=self.fit_df[phase_dif].abs() - 360.0,
                yerr=self.fit_df[f"{phase_dif}_err"].abs(),
                linestyle="",
                marker=".",
                color=color,
                label=utils.convert_amp_name(phase_dif),
            )
            ax.errorbar(
                x=self._masses,
                xerr=self._bin_width / 2,
                y=-(self.fit_df[phase_dif].abs()) + 360.0,
                yerr=self.fit_df[f"{phase_dif}_err"].abs(),
                linestyle="",
                marker=".",
                color=color,
            )

        ax.legend()
        if extend_range:
            ax.set_ylim(-360.0, 360.0)
            ax.set_yticks(np.linspace(-360, 360, 11))
        else:
            ax.set_ylim(-180.0, 180.0)
            ax.set_yticks(np.linspace(-180, 180, 5))  # force to be in pi intervals
        ax.set_ylabel(r"Phase Diff. ($^{\circ}$)", loc="top")
        ax.set_xlabel(rf"{self.channel} inv. mass $(GeV)$", loc="right")

        return ax

    def mass_phase(self, amp1: str, amp2: str, color: str = "black") -> np.ndarray:
        """_summary_

        Args:
            amp1 (str): _description_
            amp2 (str): _description_
            color (str, optional): _description_. Defaults to "black".

        Returns:
            matplotlib.axes.Axes: _description_
        """
        fig, axs = plt.subplots(
            2,
            1,
            sharex=True,
            gridspec_kw={"wspace": 0.0, "hspace": 0.07},
            height_ratios=[3, 1],
        )

        # plot the two amplitudes on the first subplot
        axs[0].errorbar(
            self._masses,
            self.fit_df[amp1],
            self.fit_df[f"{amp1}_err"],
            self._bin_width / 2,
            "o",
            color=color,
            label=utils.convert_amp_name(amp1),
        )
        axs[0].errorbar(
            self._masses,
            self.fit_df[amp2],
            self.fit_df[f"{amp2}_err"],
            self._bin_width / 2,
            "s",
            color=color,
            label=utils.convert_amp_name(amp2),
        )

        if self.truth_df is not None:
            axs[0].plot(
                self._masses,
                self.truth_df[amp1],
                linestyle="-",
                marker="",
                color=color,
            )
            axs[0].plot(
                self._masses,
                self.truth_df[amp2],
                linestyle="-",
                marker="",
                color=color,
            )

        # plot the phase difference on the second subplot
        phase_dif = self.phase_differences[(amp1, amp2)]
        axs[1].errorbar(
            self._masses,
            self.fit_df[phase_dif],
            self.fit_df[f"{phase_dif}_err"].abs(),
            self._bin_width / 2,
            linestyle="",
            marker=".",
            color=color,
        )
        axs[1].errorbar(
            self._masses,
            -self.fit_df[phase_dif],
            self.fit_df[f"{phase_dif}_err"].abs(),
            self._bin_width / 2,
            linestyle="",
            marker=".",
            color=color,
        )
        if self.truth_df is not None:
            axs[1].plot(
                self._masses,
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
        axs[1].set_xlabel(rf"{self.channel} inv. mass $(GeV)$", loc="right")

        axs[0].legend(loc="upper right")

        return axs


class DiagnosticPlotter(BasePWAPlotter):
    """Handles more complex diagnostic plots."""

    def ds_ratio(
        self, show_e852_reference: bool = True, figsize: tuple = (10, 8)
    ) -> matplotlib.axes.Axes | np.ndarray:
        """Plot the ratio and phase between D and S waves.

        Creates either a correlation plot (if dsratio is defined) or separate ratio/phase
        plots, comparing fit results to E852 reference values when requested.

        Args:
            show_e852_reference (bool): Whether to show E852 reference lines.
                Defaults to True.
            figsize (tuple): Figure size in inches (width, height). Defaults to (10, 8).

        Returns:
            matplotlib.axes.Axes | np.ndarray: Single axes object or array of axes
                depending on the plot type.

        Raises:
            ValueError: If no D or S waves are found in the fit results.

        Note:
            E852 reference values: D/S ratio = 0.27, D-S phase = ±10.54°
        """

        # Check if we have D and S waves
        d_waves = [amp for amp in self.coherent_sums["eJPmL"] if amp[-1] == "D"]
        s_waves = [amp for amp in self.coherent_sums["eJPmL"] if amp[-1] == "S"]

        if not d_waves or not s_waves:
            raise ValueError(
                "No D and/or S waves found in fit results. "
                "D/S ratio analysis requires both wave types."
            )

        # E852 reference values
        e852_ratio = 0.27
        e852_phase = 10.54

        # Two different plotting modes based on available data
        if "dsratio" in self.fit_df.columns:
            return self._plot_ds_ratio_with_correlation(
                show_e852_reference, e852_ratio, e852_phase, figsize
            )
        else:
            return self._plot_ds_ratio_from_amplitudes(
                d_waves, s_waves, show_e852_reference, e852_ratio, e852_phase, figsize
            )

    def _plot_ds_ratio_with_correlation(
        self,
        show_e852_reference: bool,
        e852_ratio: float,
        e852_phase: float,
        figsize: tuple,
    ) -> np.ndarray:
        """Plot D/S ratio with covariance when ratio is directly fitted."""

        fig, axs = plt.subplots(
            3,
            1,
            sharex=True,
            figsize=figsize,
            gridspec_kw={"wspace": 0.0, "hspace": 0.12},
        )

        ax_ratio, ax_phase, ax_corr = axs

        # Plot ratio
        if show_e852_reference:
            ax_ratio.axhline(
                y=e852_ratio,
                color="black",
                linestyle="-",
                alpha=0.7,
                label=f"E852 ratio ({e852_ratio})",
            )

        ax_ratio.errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=self.fit_df["dsratio"],
            yerr=self.fit_df["dsratio_err"],
            marker=".",
            linestyle="",
            color="darkblue",
            capsize=2,
            label="Fit Result",
        )

        # Plot phase
        if show_e852_reference:
            ax_phase.axhline(
                y=e852_phase,
                color="black",
                linestyle="--",
                alpha=0.7,
                label=f"E852 phase (±{e852_phase}°)",
            )
            ax_phase.axhline(y=-e852_phase, color="black", linestyle="--", alpha=0.7)

        ax_phase.errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=self.fit_df["dphase"],
            yerr=self.fit_df["dphase_err"].abs(),
            marker=".",
            linestyle="",
            color="darkred",
            capsize=2,
            label="Fit Result",
        )

        # Plot correlation
        if "cov_dsratio_dphase" in self.fit_df.columns:
            correlation = (
                self.fit_df["cov_dsratio_dphase"]
                / (self.fit_df["dsratio_err"] * self.fit_df["dphase_err"])
            ).fillna(0)

            ax_corr.plot(
                self._masses,
                correlation,
                marker="s",
                linestyle="",
                color="darkgreen",
                markersize=4,
                label="Correlation",
            )

        # Configure axes
        ax_ratio.set_ylabel("D/S Ratio", loc="top")
        ax_ratio.set_ylim(0, 1)
        ax_ratio.legend()

        ax_phase.set_ylabel("D-S Phase (°)", loc="top")
        ax_phase.set_yticks(np.linspace(-180, 180, 5))
        ax_phase.set_ylim(-180, 180)
        ax_phase.legend()

        ax_corr.set_xlabel(rf"{self.channel} inv. mass (GeV)", loc="right")
        ax_corr.set_ylabel("ρ(D/S, D-S)", loc="top")
        ax_corr.set_ylim(-1, 1)
        ax_corr.axhline(y=0, color="gray", linestyle=":", alpha=0.5)
        ax_corr.legend()

        plt.tight_layout()
        return axs

    def _plot_ds_ratio_from_amplitudes(
        self,
        d_waves: list,
        s_waves: list,
        show_e852_reference: bool,
        e852_ratio: float,
        e852_phase: float,
        figsize: tuple,
    ) -> np.ndarray:
        """Plot D/S ratio calculated from individual amplitudes."""

        # Marker and color maps for different quantum numbers
        marker_map = {"q": "h", "p": "^", "0": ".", "m": "v", "n": "s"}
        color_map = {"p": "red", "m": "blue"}

        fig = plt.figure(figsize=figsize, dpi=300)
        subfigs = fig.subfigures(1, 2, width_ratios=[2, 1])

        # Left subfigure: ratio and phase vs mass
        ax_ratio, ax_phase = subfigs[0].subplots(2, 1, sharex=True)

        # Right subfigure: correlation plot
        ax_corr = subfigs[1].subplots()

        # Plot E852 reference lines
        if show_e852_reference:
            ax_ratio.axhline(
                y=e852_ratio,
                color="black",
                linestyle="-",
                alpha=0.7,
                label=f"E852 ratio ({e852_ratio})",
            )
            ax_phase.axhline(
                y=e852_phase,
                color="black",
                linestyle="--",
                alpha=0.7,
                label=f"E852 phase (±{e852_phase}°)",
            )
            ax_phase.axhline(y=-e852_phase, color="black", linestyle="--", alpha=0.7)

            # Reference lines on correlation plot
            ax_corr.axhline(
                y=e852_phase,
                color="black",
                linestyle="--",
                alpha=0.7,
                label=f"E852 phase (±{e852_phase}°)",
            )
            ax_corr.axhline(y=-e852_phase, color="black", linestyle="--", alpha=0.7)
            ax_corr.axvline(
                x=e852_ratio,
                color="black",
                linestyle="-",
                alpha=0.7,
                label=f"E852 ratio ({e852_ratio})",
            )

        # Process each D wave and find corresponding S wave
        for d_wave in d_waves:
            # Find corresponding S wave (same quantum numbers except L)
            s_wave = d_wave[:-1] + "S"
            if s_wave not in s_waves:
                continue

            # Extract quantum numbers for styling
            epsilon = d_wave[0]  # reflectivity
            m = d_wave[-2]  # m-projection

            # Calculate ratio and errors
            ratio, ratio_err = self._calculate_ds_ratio(d_wave, s_wave)

            # Get phase difference
            if (d_wave, s_wave) in self.phase_differences:
                phase_dif_name = self.phase_differences[(d_wave, s_wave)]
            elif (s_wave, d_wave) in self.phase_differences:
                phase_dif_name = self.phase_differences[(s_wave, d_wave)]
            else:
                continue  # Skip if phase difference not available

            phase = self.fit_df[phase_dif_name]
            phase_err = self.fit_df[f"{phase_dif_name}_err"].abs()

            # Create label
            label = utils.convert_amp_name(d_wave).replace("D", "(D/S)")

            # Plot ratio vs mass
            ax_ratio.errorbar(
                x=self._masses,
                xerr=self._bin_width / 2,
                y=ratio,
                yerr=ratio_err,
                marker=marker_map[m],
                color=color_map[epsilon],
                linestyle="",
                markersize=6,
                capsize=2,
                label=(
                    label if epsilon == "p" else None
                ),  # Only label positive reflectivity
            )

            # Plot phase vs mass
            ax_phase.errorbar(
                x=self._masses,
                xerr=self._bin_width / 2,
                y=phase,
                yerr=phase_err,
                marker=marker_map[m],
                color=color_map[epsilon],
                linestyle="",
                markersize=6,
                capsize=2,
            )

            # Plot sign ambiguity
            ax_phase.errorbar(
                x=self._masses,
                xerr=self._bin_width / 2,
                y=-phase,
                yerr=phase_err,
                marker=marker_map[m],
                color=color_map[epsilon],
                linestyle="",
                markersize=6,
                alpha=0.6,
                capsize=2,
            )

            # Plot correlation (ratio vs phase)
            ax_corr.errorbar(
                x=ratio,
                xerr=ratio_err,
                y=phase,
                yerr=phase_err,
                marker=marker_map[m],
                color=color_map[epsilon],
                linestyle="",
                markersize=6,
                capsize=2,
                label=label,
            )

            # Plot sign ambiguity on correlation plot
            ax_corr.errorbar(
                x=ratio,
                xerr=ratio_err,
                y=-phase,
                yerr=phase_err,
                marker=marker_map[m],
                color=color_map[epsilon],
                linestyle="",
                markersize=6,
                alpha=0.6,
                capsize=2,
            )

        # Configure left subplot axes
        ax_ratio.set_ylabel("D/S Ratio", loc="top", fontsize=12)
        ax_ratio.set_yscale("log")
        ax_ratio.legend()

        ax_phase.set_xlabel(
            rf"{self.channel} inv. mass (GeV)", loc="right", fontsize=12
        )
        ax_phase.set_ylabel("D-S Phase (°)", loc="top", fontsize=12)
        ax_phase.set_yticks(np.linspace(-180, 180, 5))
        ax_phase.set_ylim(-180, 180)

        # Configure correlation plot
        ax_corr.set_xlabel("D/S Ratio", loc="right", fontsize=12)
        ax_corr.set_ylabel("D-S Phase (°)", loc="top", fontsize=12)
        ax_corr.set_yticks(np.linspace(-180, 180, 5))
        ax_corr.set_ylim(-180, 180)
        ax_corr.legend(fontsize=10)

        plt.tight_layout()
        return np.array([ax_ratio, ax_phase, ax_corr])

    def _calculate_ds_ratio(self, d_wave: str, s_wave: str) -> tuple:
        """Calculate D/S ratio and propagate errors."""

        # Calculate ratio using square roots (amplitude ratios)
        d_intensity = self.fit_df[d_wave]
        s_intensity = self.fit_df[s_wave]
        d_err = self.fit_df[f"{d_wave}_err"]
        s_err = self.fit_df[f"{s_wave}_err"]

        # Avoid division by zero
        safe_d = np.where(d_intensity == 0, np.finfo(float).eps, d_intensity)
        safe_s = np.where(s_intensity == 0, np.finfo(float).eps, s_intensity)

        ratio = np.sqrt(safe_d) / np.sqrt(safe_s)

        # Error propagation for sqrt(A)/sqrt(B)
        ratio_err = ratio * np.sqrt(
            0.25 * np.square(d_err / safe_d) + 0.25 * np.square(s_err / safe_s)
        )

        return ratio, ratio_err

    def matrix(
        self, dpi: int = 300, figsize: tuple = (12, 8), show_errors: bool = True
    ) -> np.ndarray:
        """Create a matrix of intensities and phase differences.

        The diagonal contains the intensity of each wave in both reflectivities.
        The off-diagonals show phase differences between waves, with positive
        reflectivity in the upper triangle and negative in the lower triangle.

        Args:
            dpi (int): Resolution of the figure. Defaults to 300. Increase for
                large models to improve readability.
            figsize (tuple): Figure size in inches (width, height). Defaults to (12, 8).
            show_errors (bool): Whether to show error bands for phase differences.
                Defaults to True.

        Returns:
            np.ndarray: 2D array of matplotlib.axes.Axes objects for further
                customization.

        Raises:
            ValueError: If no wave amplitudes are found in the fit results.

        Note:
            This plot can become cramped for large models. Consider adjusting dpi
            or figsize, or using a subset of waves for better readability.
        """

        # Validate that we have wave data to plot
        if not self.coherent_sums["JPmL"]:
            raise ValueError("No wave amplitudes found in fit results.")

        n_waves = len(self.coherent_sums["JPmL"])

        # Create figure
        fig, axs = plt.subplots(
            n_waves,
            n_waves,
            sharex=True,
            figsize=figsize,
            dpi=dpi,
        )

        # Handle single wave case (ensure axs is always 2D)
        if n_waves == 1:
            axs = np.array([[axs]])
        elif n_waves == 2:
            axs = axs.reshape(n_waves, n_waves)

        # Track plot handles for legend
        pos_plot_handle = None
        neg_plot_handle = None

        # Calculate intensity scale for diagonal plots
        max_intensity = max(
            self.fit_df[amp_name].max()
            for amp_name in self.coherent_sums["eJPmL"]
            if amp_name in self.fit_df.columns
        )

        # Create tick labels for intensity scale
        intensity_ticks = np.linspace(0, max_intensity, 4)
        intensity_labels = [utils.human_format(num) for num in intensity_ticks]

        # Set up phase tick values
        phase_ticks = np.linspace(-180, 180, 5)

        # Plot matrix elements
        for row, JPmL_row in enumerate(self.coherent_sums["JPmL"]):
            for col, JPmL_col in enumerate(self.coherent_sums["JPmL"]):
                ax = axs[row, col]

                # Configure tick parameters for all subplots
                ax.tick_params("both", labelsize=8)

                # Show x-axis labels only on bottom row
                if row == n_waves - 1:
                    ax.xaxis.set_tick_params(labelbottom=True)

                # DIAGONAL ELEMENTS: Plot wave intensities
                if row == col:
                    self._plot_diagonal_intensity(
                        ax,
                        JPmL_row,
                        max_intensity,
                        intensity_ticks,
                        intensity_labels,
                        row,
                    )

                    # Store plot handles from first diagonal element
                    if row == 0:
                        pos_plot_handle = ax.lines[-1] if ax.lines else None
                        neg_plot_handle = ax.lines[-2] if len(ax.lines) >= 2 else None

                # UPPER TRIANGLE: Positive reflectivity phase differences
                elif col > row:
                    self._plot_upper_triangle_phase(
                        ax, JPmL_row, JPmL_col, phase_ticks, row, show_errors
                    )

                # LOWER TRIANGLE: Negative reflectivity phase differences
                elif col < row:
                    self._plot_lower_triangle_phase(
                        ax, JPmL_row, JPmL_col, phase_ticks, row, col, show_errors
                    )

        # Configure figure-wide labels and styling
        fig.text(
            0.5, 0.04, rf"{self.channel} inv. mass (GeV)", ha="center", fontsize=14
        )
        fig.text(
            0.04,
            0.5,
            r"Phase Differences ($^{\circ}$)",
            ha="center",
            va="center",
            fontsize=14,
            rotation="vertical",
        )

        # Add legend if we have plot handles
        legend_handles = []
        if pos_plot_handle is not None:
            legend_handles.append(pos_plot_handle)
        if neg_plot_handle is not None:
            legend_handles.append(neg_plot_handle)

        if legend_handles:
            fig.legend(
                handles=legend_handles,
                loc="upper right",
                bbox_to_anchor=(0.98, 0.98),
                fontsize=10,
            )

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.08, left=0.08)

        return axs

    def _plot_diagonal_intensity(
        self,
        ax: matplotlib.axes.Axes,
        JPmL: str,
        max_intensity: float,
        intensity_ticks: np.ndarray,
        intensity_labels: list,
        row: int,
    ) -> None:
        """Plot wave intensities on diagonal elements of the matrix."""

        ax.set_ylim(0, max_intensity)

        # Emphasize diagonal with thick borders
        for spine in ax.spines.values():
            spine.set_linewidth(2)
            spine.set_color("black")

        # Set labels for first row/column
        if row == 0:
            title = utils.convert_amp_name(JPmL).replace(r"^{(\Sigma\varepsilon)}", "")
            ax.set_title(title, fontsize=10, pad=8)
            ax.set_ylabel(title, fontsize=10)
            ax.set_yticks(intensity_ticks, intensity_labels)
        else:
            # Remove y-tick labels for other diagonal elements to save space
            ax.set_yticks(intensity_ticks, [""] * len(intensity_ticks))

        # Plot negative reflectivity
        neg_amp = f"m{JPmL}"
        if neg_amp in self.fit_df.columns:
            ax.errorbar(
                x=self._masses,
                xerr=self._bin_width / 2,
                y=self.fit_df[neg_amp],
                yerr=self.fit_df[f"{neg_amp}_err"],
                marker="v",
                color="blue",
                markersize=3,
                linewidth=1,
                capsize=1,
                label=r"$\varepsilon = -1$",
            )

            # Plot truth data if available
            if self.truth_df is not None and neg_amp in self.truth_df.columns:
                ax.plot(
                    self._masses,
                    self.truth_df[neg_amp],
                    linestyle="-",
                    linewidth=2,
                    color="blue",
                    alpha=0.8,
                )

        # Plot positive reflectivity
        pos_amp = f"p{JPmL}"
        if pos_amp in self.fit_df.columns:
            ax.errorbar(
                x=self._masses,
                xerr=self._bin_width / 2,
                y=self.fit_df[pos_amp],
                yerr=self.fit_df[f"{pos_amp}_err"],
                marker="^",
                color="red",
                markersize=3,
                linewidth=1,
                capsize=1,
                label=r"$\varepsilon = +1$",
            )

            # Plot truth data if available
            if self.truth_df is not None and pos_amp in self.truth_df.columns:
                ax.plot(
                    self._masses,
                    self.truth_df[pos_amp],
                    linestyle="-",
                    linewidth=2,
                    color="red",
                    alpha=0.8,
                )

    def _plot_upper_triangle_phase(
        self,
        ax: matplotlib.axes.Axes,
        JPmL_row: str,
        JPmL_col: str,
        phase_ticks: np.ndarray,
        row: int,
        show_errors: bool,
    ) -> None:
        """Plot positive reflectivity phase differences in upper triangle."""

        ax.set_ylim(-180, 180)
        ax.set_yticks(phase_ticks, [""] * len(phase_ticks))  # No labels to save space

        # Set column title for top row
        if row == 0:
            title = utils.convert_amp_name(JPmL_col).replace(
                r"^{(\Sigma\varepsilon)}", ""
            )
            ax.set_title(title, fontsize=10, pad=8)

        # Get phase difference data
        amp1, amp2 = f"p{JPmL_row}", f"p{JPmL_col}"
        if (amp1, amp2) in self.phase_differences:
            phase_dif = self.phase_differences[(amp1, amp2)]
        elif (amp2, amp1) in self.phase_differences:
            phase_dif = self.phase_differences[(amp2, amp1)]
        else:
            return  # Skip if phase difference not available

        # Plot phase values and their sign ambiguity
        phase_values = self.fit_df[phase_dif]
        phase_errors = self.fit_df[f"{phase_dif}_err"].abs()

        # Plot positive values
        ax.plot(
            self._masses,
            phase_values.abs(),
            linestyle="",
            marker=".",
            markersize=2,
            color="red",
        )

        # Plot negative values (sign ambiguity)
        ax.plot(
            self._masses,
            -phase_values.abs(),
            linestyle="",
            marker=".",
            markersize=2,
            color="red",
        )

        # Add error bands if requested
        if show_errors:
            ax.fill_between(
                self._masses,
                phase_values.abs() - phase_errors,
                phase_values.abs() + phase_errors,
                color="red",
                alpha=0.15,
            )
            ax.fill_between(
                self._masses,
                -phase_values.abs() - phase_errors,
                -phase_values.abs() + phase_errors,
                color="red",
                alpha=0.15,
            )

        # Plot truth data if available
        if self.truth_df is not None and phase_dif in self.truth_df.columns:
            ax.plot(
                self._masses,
                self.truth_df[phase_dif],
                linestyle="-",
                linewidth=2,
                color="red",
                alpha=0.8,
            )

    def _plot_lower_triangle_phase(
        self,
        ax: matplotlib.axes.Axes,
        JPmL_row: str,
        JPmL_col: str,
        phase_ticks: np.ndarray,
        row: int,
        col: int,
        show_errors: bool,
    ) -> None:
        """Plot negative reflectivity phase differences in lower triangle."""

        ax.set_ylim(-180, 180)

        # Set y-axis labels only for first column
        if col == 0:
            title = utils.convert_amp_name(JPmL_row).replace(
                r"^{(\Sigma\varepsilon)}", ""
            )
            ax.set_ylabel(title, fontsize=10)
            ax.set_yticks(phase_ticks)
        else:
            ax.set_yticks(phase_ticks, [""] * len(phase_ticks))

        # Get phase difference data
        amp1, amp2 = f"m{JPmL_row}", f"m{JPmL_col}"
        if (amp1, amp2) in self.phase_differences:
            phase_dif = self.phase_differences[(amp1, amp2)]
        elif (amp2, amp1) in self.phase_differences:
            phase_dif = self.phase_differences[(amp2, amp1)]
        else:
            return  # Skip if phase difference not available

        # Plot phase values and their sign ambiguity
        phase_values = self.fit_df[phase_dif]
        phase_errors = self.fit_df[f"{phase_dif}_err"].abs()

        # Plot positive values
        ax.plot(
            self._masses,
            phase_values.abs(),
            linestyle="",
            marker=".",
            markersize=2,
            color="blue",
        )

        # Plot negative values (sign ambiguity)
        ax.plot(
            self._masses,
            -phase_values.abs(),
            linestyle="",
            marker=".",
            markersize=2,
            color="blue",
        )

        # Add error bands if requested
        if show_errors:
            ax.fill_between(
                self._masses,
                phase_values.abs() - phase_errors,
                phase_values.abs() + phase_errors,
                color="blue",
                alpha=0.15,
            )
            ax.fill_between(
                self._masses,
                -phase_values.abs() - phase_errors,
                -phase_values.abs() + phase_errors,
                color="blue",
                alpha=0.15,
            )

        # Plot truth data if available
        if self.truth_df is not None and phase_dif in self.truth_df.columns:
            ax.plot(
                self._masses,
                self.truth_df[phase_dif],
                linestyle="-",
                linewidth=2,
                color="blue",
                alpha=0.8,
            )


class RandomPlotter(BasePWAPlotter):
    """Handles plots that rely on randomized parameter fitting"""

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        random_df: Optional[pd.DataFrame] = None,
    ) -> None:
        """Ensure random DataFrame is not None and initialize the plotter."""
        super().__init__(fit_df, data_df, random_df)
        if self.random_df is None:
            raise ValueError("RandomPlotter requires a non-None random_df.")

    def chi2_likelihood(self) -> matplotlib.axes.Axes:
        fig, ax = plt.subplots()

        # calculate chi2 / ndf based off intensities and phase differences between
        #   best fit and random samples
        # plot against Delta(-2ln(L)). Idea is that best fit is in bottom left corner,
        #   and other rand fits with different chi2 should also have different
        #   likelihoods
        # color code points corresponding to eMatrixStatus

        return ax


class BootstrapPlotter(BasePWAPlotter):
    """Handles plots that rely on bootstrap results."""

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        bootstrap_df: Optional[pd.DataFrame] = None,
        truth_df: Optional[pd.DataFrame] = None,
    ) -> None:
        """Ensure bootstrap DataFrame is not None and initialize the plotter."""
        super().__init__(fit_df, data_df, bootstrap_df, truth_df)
        if self.bootstrap_df is None:
            raise ValueError("BootstrapPlotter requires a non-None bootstrap_df.")

    def correlation_matrix(
        self,
        columns: list[str],
        pdf_path: str,
        method: Literal["pearson", "kendall", "spearman"] = "pearson",
        directories: Optional[list[str]] = None,
        figsize: tuple = (10, 8),
        annot: bool = True,
        cmap: str = "RdBu_r",
    ) -> None:
        """Generate correlation matrices for bootstrap samples and save to PDF.

        Creates correlation matrix heatmaps for each bootstrap sample directory,
        with proper amplitude name formatting and statistical annotations.

        Args:
            columns (list[str]): Bootstrap DataFrame columns to calculate correlations for.
            pdf_path (str): Path to save the PDF file containing all correlation matrices.
            method (str): Correlation method. Options: 'pearson', 'kendall', 'spearman'.
                Defaults to 'pearson'.
            directories (list[str], optional): Specific bootstrap directories to plot.
                If None, plots all available directories. Defaults to None.
            figsize (tuple): Figure size in inches (width, height). Defaults to (10, 8).
            annot (bool): Whether to annotate heatmap cells with correlation values.
                Defaults to True.
            cmap (str): Colormap for the heatmap. Defaults to 'RdBu_r'.

        Raises:
            ValueError: If bootstrap DataFrame is None or contains missing columns.
            FileNotFoundError: If specified directories are not found.

        Returns:
            None: Saves plots to PDF file and prints confirmation.
        """

        # Tell type checkers that self.bootstrap_df is not None
        assert self.bootstrap_df is not None

        # Check for missing columns
        missing_columns = [
            col for col in columns if col not in self.bootstrap_df.columns
        ]
        if missing_columns:
            raise ValueError(
                f"The following columns are missing from bootstrap_df: {missing_columns}"
            )

        # Validate method
        valid_methods = ["pearson", "kendall", "spearman"]
        if method not in valid_methods:
            raise ValueError(f"Method must be one of {valid_methods}, got '{method}'")

        # Handle directories
        available_dirs = self.bootstrap_df["directory"].unique()  # type: ignore
        if directories is None:
            directories = available_dirs.tolist()
        else:
            missing_dirs = set(directories) - set(available_dirs)
            if missing_dirs:
                raise FileNotFoundError(
                    f"The following directories were not found: {missing_dirs}"
                )

        # At this point directories is guaranteed to be a list
        assert directories is not None

        # Create PDF to save all plots
        with PdfPages(pdf_path) as pdf:
            for directory in directories:
                # Filter data for current directory
                dir_data = self.bootstrap_df[
                    self.bootstrap_df["directory"] == directory
                ][columns]

                if dir_data.empty:
                    print(f"Warning: No data found for directory {directory}")
                    continue

                # Calculate correlation matrix
                corr_matrix = dir_data.corr(method=method)

                # Create prettier labels for amplitudes and phases
                label_mapping = {}
                for col in columns:
                    if any(
                        col in sublist for sublist in self.coherent_sums.values()
                    ) or col in set(self.phase_differences.values()):
                        label_mapping[col] = utils.convert_amp_name(col)
                    else:
                        label_mapping[col] = col

                # Apply label mapping
                corr_matrix = corr_matrix.rename(
                    index=label_mapping, columns=label_mapping
                )

                # Create the plot
                plt.figure(figsize=figsize)

                # Generate heatmap
                heatmap = sns.heatmap(
                    corr_matrix,
                    annot=annot,
                    cmap=cmap,
                    vmin=-1,
                    vmax=1,
                    center=0,
                    square=True,
                    fmt=".2f" if annot else "",
                    cbar_kws={"label": f"{method.capitalize()} Correlation"},
                )

                # Customize plot appearance
                plt.title(
                    f"{method.capitalize()} Correlation Matrix\n{directory}",
                    fontsize=14,
                    pad=20,
                )
                plt.xticks(rotation=45, ha="right")
                plt.yticks(rotation=0)

                # Add grid for better readability
                heatmap.set_facecolor("white")

                # Tight layout to prevent label cutoff
                plt.tight_layout()

                # Save to PDF
                pdf.savefig(bbox_inches="tight", dpi=300)
                plt.close()

        print(f"Correlation matrices saved to: {pdf_path}")
        print(f"Generated {len(directories)} correlation matrix plots")

    def pairplot(
        self,
        fit_indices: list[int],
        columns: list[str],
        show_truth: bool = True,
        show_uncertainty_bands: bool = True,
        correlation_threshold: float = 0.7,
        figsize: Optional[tuple] = None,
        **kwargs,
    ):
        """Create a comprehensive pairplot matrix of bootstrap fit results.

        Generates a scatterplot matrix showing relationships between fit parameters,
        with histograms/KDE on diagonal, scatter plots on lower triangle, and 2D KDE
        contours on upper triangle. Includes uncertainty bands and truth value overlays.

        Args:
            fit_indices (list[int]): Indices of fits to include in the analysis.
                These should correspond to indices in the linked DataFrames.
            columns (list[str]): Fit parameters to plot in the matrix.
            show_truth (bool): Whether to overlay truth values if available.
                Defaults to True.
            show_uncertainty_bands (bool): Whether to show MINUIT uncertainty bands.
                Defaults to True.
            correlation_threshold (float): Threshold for highlighting strongly
                correlated parameters with black borders. Defaults to 0.7.
            figsize (tuple, optional): Figure size. If None, calculated automatically.
            **kwargs: Additional arguments passed to seaborn.PairGrid.

        Raises:
            ValueError: If bootstrap DataFrame is None or parameters are invalid.
            IndexError: If fit_indices are not found in the DataFrames.

        Returns:
            seaborn.PairGrid: The configured PairGrid object for further customization.

        Note:
            Amplitudes are automatically converted to fit fractions for visualization.
            Phase differences use circular statistics for standard deviation
                calculations.
        """

        # Validate inputs
        if self.bootstrap_df is None:
            raise ValueError(
                "Bootstrap DataFrame required for pairplot analysis. "
                "Provide bootstrap_df during initialization."
            )

        if not fit_indices:
            raise ValueError("At least one fit index must be specified.")

        if not columns:
            raise ValueError("At least one column must be specified for plotting.")

        # Validate columns exist
        missing_columns = [
            col for col in columns if col not in self.bootstrap_df.columns
        ]
        if missing_columns:
            raise ValueError(
                "The following columns are missing from bootstrap_df:"
                f" {missing_columns}"
            )

        # Verify bootstrap samples exist for all requested fit indices
        available_indices = set(self.bootstrap_df["fit_index"].unique())
        missing_indices = [idx for idx in fit_indices if idx not in available_indices]
        if missing_indices:
            raise IndexError(
                f"No bootstrap samples found for fit indices: {missing_indices}"
            )

        # Prepare data
        data = self._prepare_pairplot_data(columns, fit_indices, show_truth)

        # Set up the plot
        if figsize is None:
            figsize = (max(10, len(columns) * 2), max(8, len(columns) * 1.5))

        # Create PairGrid
        palette = sns.color_palette(n_colors=len(fit_indices))

        grid_kwargs = {
            "hue": "fit_index",
            "palette": palette,
            "height": figsize[0] / len(columns),
            **kwargs,
        }

        pg = sns.PairGrid(data["plot_data"], **grid_kwargs)

        # Configure plot types based on number of fits
        if len(fit_indices) < 4:
            pg.map_diag(sns.histplot, kde=True, alpha=0.7)
        else:
            pg.map_diag(sns.kdeplot, alpha=0.7)

        # Map plots to grid
        pg.map_upper(sns.kdeplot, levels=[0.003, 0.05, 0.32], alpha=0.6)
        pg.map_lower(sns.scatterplot, alpha=0.7, s=30)

        # Add overlays
        self._add_pairplot_overlays(
            pg,
            data,
            columns,
            palette,
            show_truth,
            show_uncertainty_bands,
            correlation_threshold,
        )

        # Update labels to prettier versions
        self._update_pairplot_labels(pg, columns)

        # Get fit file names for title
        fit_files = self.fit_df.loc[fit_indices, "file"].tolist()
        plt.suptitle(f"Bootstrap Analysis: {', '.join(fit_files)}", fontsize=16, y=1.02)
        plt.tight_layout()

        return pg

    def _prepare_pairplot_data(
        self, columns: list, fit_indices: list, show_truth: bool
    ) -> dict:
        """Prepare and normalize data for pairplot visualization."""

        # Extract relevant data using fit_indices directly
        fit_data = self.fit_df.loc[fit_indices][
            columns + [f"{c}_err" for c in columns] + ["detected_events"]
        ].copy()

        bootstrap_data = self.bootstrap_df[  # type: ignore
            self.bootstrap_df["fit_index"].isin(fit_indices)  # type: ignore
        ][columns + ["fit_index", "detected_events"]].copy()

        truth_data = None
        if self.truth_df is not None and show_truth:
            truth_data = self.truth_df[self.truth_df.index.isin(fit_indices)][
                columns + ["detected_events"]
            ].copy()

        # Convert amplitudes to fit fractions
        plot_data = bootstrap_data.copy()
        for col in columns:
            if any(col in sublist for sublist in self.coherent_sums.values()):
                plot_data[col] = plot_data[col] / bootstrap_data["detected_events"]

        return {
            "fit_data": fit_data,
            "bootstrap_data": bootstrap_data,
            "plot_data": plot_data,
            "truth_data": truth_data,
        }

    def _add_pairplot_overlays(
        self,
        pg: sns.PairGrid,
        data: dict,
        columns: list,
        palette,
        show_truth: bool,
        show_uncertainty_bands: bool,
        correlation_threshold: float,
    ) -> None:
        """Add uncertainty bands, truth values, and correlation highlights."""

        num_plots = len(columns)
        fit_data = data["fit_data"]
        truth_data = data["truth_data"]
        bootstrap_data = data["bootstrap_data"]

        for row in range(num_plots):
            for col in range(num_plots):
                ax = pg.axes[row, col]

                # Get column labels
                col_label = columns[col]
                row_label = columns[row]

                # Calculate scaling factors for fit fractions
                x_scaling = self._get_scaling_factor(col_label, fit_data)
                y_scaling = self._get_scaling_factor(row_label, fit_data)

                # Add uncertainty bands
                if show_uncertainty_bands:
                    self._add_uncertainty_bands(
                        ax,
                        fit_data,
                        col_label,
                        row_label,
                        x_scaling,
                        y_scaling,
                        palette,
                        row == col,
                    )

                # Add truth value markers
                if show_truth and truth_data is not None:
                    self._add_truth_markers(
                        ax,
                        truth_data,
                        col_label,
                        row_label,
                        x_scaling,
                        y_scaling,
                        palette,
                        row == col,
                    )

                # Highlight strong correlations
                if row != col:
                    self._highlight_correlations(
                        ax,
                        bootstrap_data,
                        col_label,
                        row_label,
                        fit_data.index,
                        correlation_threshold,
                    )

                # Add statistical annotations on diagonal
                if row == col:
                    self._add_diagonal_annotations(
                        ax, bootstrap_data, fit_data, col_label
                    )

    def _get_scaling_factor(self, column: str, fit_data: pd.DataFrame) -> pd.Series:
        """Get appropriate scaling factor for amplitude vs phase columns.

        Returns:
            pd.Series: # of detected events for coherent sums or amplitudes, otherwise
                ones.
        """
        if any(column in sublist for sublist in self.coherent_sums.values()):
            return fit_data["detected_events"]
        return pd.Series(1.0, index=fit_data.index)

    def _add_uncertainty_bands(
        self,
        ax: matplotlib.axes.Axes,
        fit_data: pd.DataFrame,
        col_label: str,
        row_label: str,
        x_scaling: pd.Series,
        y_scaling: pd.Series,
        palette,
        is_diagonal: bool,
    ) -> None:
        """Add MINUIT uncertainty bands to the plot."""

        for idx, color in zip(fit_data.index, palette):
            # X-direction uncertainty band
            x_center = fit_data.loc[idx, col_label] / x_scaling.loc[idx]
            x_error = fit_data.loc[idx, f"{col_label}_err"] / x_scaling.loc[idx]

            ax.axvspan(x_center - x_error, x_center + x_error, color=color, alpha=0.2)

            # Y-direction uncertainty band (only for off-diagonal)
            if not is_diagonal:
                y_center = fit_data.loc[idx, row_label] / y_scaling.loc[idx]
                y_error = fit_data.loc[idx, f"{row_label}_err"] / y_scaling.loc[idx]

                ax.axhspan(
                    y_center - y_error, y_center + y_error, color=color, alpha=0.2
                )

    def _add_truth_markers(
        self,
        ax,
        truth_data: pd.DataFrame,
        col_label: str,
        row_label: str,
        x_scaling: pd.Series,
        y_scaling: pd.Series,
        palette,
        is_diagonal: bool,
    ) -> None:
        """Add truth value markers to the plot."""

        for idx, color in zip(truth_data.index, palette):
            x_truth = truth_data.loc[idx, col_label] / x_scaling.loc[idx]

            if is_diagonal:
                # Vertical line for diagonal plots
                ax.axvline(
                    x=x_truth, color=color, linestyle="-", alpha=0.8, linewidth=2
                )
            else:
                # Scatter point for off-diagonal plots
                y_truth = truth_data.loc[idx, row_label] / y_scaling.loc[idx]
                ax.scatter(
                    x_truth,
                    y_truth,
                    color=color,
                    marker="X",
                    s=100,
                    edgecolors="black",
                    linewidth=1,
                    zorder=100,
                )

    def _highlight_correlations(
        self,
        ax,
        bootstrap_data: pd.DataFrame,
        col_label: str,
        row_label: str,
        fit_indices: list,
        threshold: float,
    ) -> None:
        """Highlight strongly correlated parameters with black borders."""

        correlations = []
        for fit_idx in fit_indices:
            subset = bootstrap_data[bootstrap_data["fit_index"] == fit_idx]
            if len(subset) > 1:  # Need at least 2 points for correlation
                corr = subset[[col_label, row_label]].corr().iloc[0, 1]
                correlations.append(abs(corr))

        if correlations and np.mean(correlations) > threshold:
            ax.patch.set_edgecolor("black")
            ax.patch.set_linewidth(3)

    def _add_diagonal_annotations(
        self,
        ax,
        bootstrap_data: pd.DataFrame,
        fit_data: pd.DataFrame,
        col_label: str,
    ) -> None:
        """Add bootstrap/MINUIT ratio annotations to diagonal plots."""

        ratios = []
        for fit_idx, error in zip(fit_data.index, fit_data[f"{col_label}_err"]):
            subset = bootstrap_data[bootstrap_data["fit_index"] == fit_idx][col_label]

            if len(subset) > 1:
                # Use circular statistics for phase differences
                if col_label in set(self.phase_differences.values()):
                    radian_phases = np.deg2rad(subset)
                    stdev = scipy.stats.circstd(
                        radian_phases, low=int(-np.pi), high=int(np.pi)
                    )
                    stdev = np.rad2deg(stdev)
                else:
                    stdev = subset.std()

                ratios.append(stdev / error if error > 0 else 0)

        if ratios:
            avg_ratio = np.mean(ratios)
            ax.text(
                0.6,
                0.9,
                rf"$\frac{{\sigma_{{bootstrap}}}}{{\sigma_{{MINUIT}}}} = {avg_ratio:.2f}$",
                transform=ax.transAxes,
                fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7),
            )

    def _update_pairplot_labels(self, pg, columns: list) -> None:
        """Update axis labels to prettier LaTeX formatting."""

        num_plots = len(columns)
        for i, col in enumerate(columns):
            # Update x-axis label (bottom row)
            if hasattr(pg.axes[num_plots - 1, i], "xaxis"):
                pretty_label = self._get_pretty_label(col)
                pg.axes[num_plots - 1, i].set_xlabel(pretty_label)

            # Update y-axis label (first column)
            if hasattr(pg.axes[i, 0], "yaxis"):
                pretty_label = self._get_pretty_label(col)
                pg.axes[i, 0].set_ylabel(pretty_label)

    def _get_pretty_label(self, column: str) -> str:
        """Convert column name to pretty LaTeX formatting."""
        if any(
            column in sublist for sublist in self.coherent_sums.values()
        ) or column in set(self.phase_differences.values()):
            return utils.convert_amp_name(column)
        return column

    def joyplot(
        self,
        columns: list[str],
        fits: Optional[list[str]] = None,
        colormap: Optional[str] = None,
        sparse_labels: bool = True,
        overlap: float = 2.0,
        show_truth: bool = True,
        figsize: Optional[tuple] = None,
        **kwargs,
    ) -> np.ndarray:
        """Create a ridge plot (joyplot) of bootstrap parameter distributions.

        Generates stacked kernel density estimates for each mass bin, providing
        an intuitive view of how the bootstrap distributions change from bin-to-bin.
        Truth values are overlaid as vertical lines when available.

        Args:
            columns (list[str]): Bootstrap DataFrame columns to plot.
            fits (list[str], optional): Specific fit files to include. If None,
                uses all available fits. Defaults to None.
            colormap (str, optional): Matplotlib colormap name. Defaults to 'Accent'.
            sparse_labels (bool): If True, only label mass bins at multiples of 0.1 GeV.
                Defaults to True.
            overlap (float): Vertical overlap between distributions. Higher values
                create more compact plots. Defaults to 2.0.
            show_truth (bool): Whether to overlay truth values as vertical lines.
                Defaults to True.
            figsize (tuple, optional): Figure size. If None, calculated automatically.
            **kwargs: Additional arguments passed to joypy.joyplot.

        Raises:
            ValueError: If bootstrap DataFrame is None or columns are missing.
            TypeError: If colormap is not a valid string.

        Returns:
            None: Displays the joyplot.

        Warning:
            Mixing different parameter types (e.g., phases and intensities) may
            result in misleading visualizations due to different scales.
        """

        # Validate bootstrap DataFrame
        if self.bootstrap_df is None:
            raise ValueError("Bootstrap DataFrame required for joyplot analysis.")

        # Validate columns
        missing_columns = []
        for col in columns:
            if col not in self.bootstrap_df.columns:
                missing_columns.append(col)
            if col not in self.fit_df.columns:
                missing_columns.append(col)
            if (
                self.truth_df is not None
                and show_truth
                and col not in self.truth_df.columns
            ):
                missing_columns.append(col)

        if missing_columns:
            raise ValueError(
                f"Missing columns in one or more DataFrames: {set(missing_columns)}"
            )

        # Handle colormap
        color_list = self._process_colormap(colormap, len(columns))

        # Prepare data
        plot_data = self._prepare_joyplot_data(fits)

        # Set up axis limits and labels
        x_limits, y_labels = self._setup_joyplot_axes(
            plot_data["bootstrap_data"], plot_data["truth_data"], columns, sparse_labels
        )

        # Configure joyplot parameters
        joy_kwargs = self._setup_joyplot_kwargs(
            color_list, overlap, x_limits, y_labels, figsize, kwargs
        )

        # Create the joyplot
        fig, axes = joypy.joyplot(
            plot_data["bootstrap_data"], by="fit_index", column=columns, **joy_kwargs
        )

        # Add truth value overlays
        if show_truth and plot_data["truth_data"] is not None:
            self._add_joyplot_truth_lines(
                axes, plot_data["truth_data"], columns, color_list
            )

        # Update legend labels to pretty format
        self._update_joyplot_legend(axes)

        plt.suptitle("Bootstrap Parameter Distributions", fontsize=16, y=1.02)
        return axes

    def _process_colormap(self, colormap: Optional[str], n_colors: int) -> list:
        """Process and validate colormap, returning list of colors."""

        if colormap is None:
            colormap = "Accent"

        if not isinstance(colormap, str):
            raise TypeError("Colormap must be a string.")

        if colormap not in matplotlib.colormaps:
            available_maps = ", ".join(list(matplotlib.colormaps.keys())[:10])
            raise ValueError(
                f"'{colormap}' is not a valid matplotlib colormap. "
                f"Available options include: {available_maps}..."
            )

        # Get colors from colormap
        color_list = list(matplotlib.colormaps[colormap].colors)  # type: ignore

        # Repeat colormap if too short
        if len(color_list) < n_colors:
            repeats = (n_colors // len(color_list)) + 1
            color_list = (color_list * repeats)[:n_colors]
        else:
            color_list = color_list[:n_colors]

        return color_list

    def _prepare_joyplot_data(self, fits: Optional[list]) -> dict:
        """Prepare data for joyplot, handling fit selection and truth values."""

        # Determine fit indices to use
        if fits is None:
            fit_indices = self.fit_df.index.tolist()
        else:
            fit_indices = self.fit_df[self.fit_df["file"].isin(fits)].index.tolist()

        if not fit_indices:
            raise ValueError("No valid fit files found.")

        # Extract bootstrap data (already validated above)
        bootstrap_data = self.bootstrap_df[  # type: ignore
            self.bootstrap_df["fit_index"].isin(fit_indices)  # type: ignore
        ].copy()

        # Extract truth data if available
        truth_data = None
        if self.truth_df is not None:
            truth_data = self.truth_df[self.truth_df.index.isin(fit_indices)].copy()

        return {
            "bootstrap_data": bootstrap_data,
            "truth_data": truth_data,
            "fit_indices": fit_indices,
        }

    def _setup_joyplot_axes(
        self,
        bootstrap_data: pd.DataFrame,
        truth_data: Optional[pd.DataFrame],
        columns: list,
        sparse_labels: bool,
    ) -> tuple:
        """Set up x-axis limits and y-axis labels for joyplot."""

        # Calculate x-axis limits
        x_min = bootstrap_data[columns].min().min()
        x_max = bootstrap_data[columns].max().max()

        # Include truth data in range calculation
        if truth_data is not None:
            x_min = min(x_min, truth_data[columns].min().min())
            x_max = max(x_max, truth_data[columns].max().max())

        # Set sensible limits based on data type
        if all(
            any(col in sublist for sublist in self.coherent_sums.values())
            for col in columns
        ):
            # All intensity columns - start from zero
            x_min = 0.0
        elif all(col in set(self.phase_differences.values()) for col in columns):
            # All phase differences - use standard phase range
            x_min, x_max = -180.0, 180.0
        else:
            # Mixed or other types - add some padding
            x_range = x_max - x_min
            x_min -= 0.05 * x_range
            x_max += 0.05 * x_range

        # Create y-axis labels (mass bin ranges)
        y_labels = self._create_mass_bin_labels(bootstrap_data, sparse_labels)

        return (x_min, x_max), y_labels

    def _create_mass_bin_labels(
        self, bootstrap_data: pd.DataFrame, sparse_labels: bool
    ) -> list:
        """Create descriptive labels for mass bins."""

        # Get unique fit indices and corresponding mass ranges
        unique_indices = bootstrap_data["fit_index"].unique()
        data_subset = self.data_df[self.data_df["fit_index"].isin(unique_indices)]

        labels = []
        for idx, (low, high) in enumerate(
            zip(data_subset["m_low"], data_subset["m_high"])
        ):
            label = f"{low:.3f}-{high:.3f}"

            if not sparse_labels:
                labels.append(label)
            else:
                # Only label bins at "nice" intervals
                tolerance = 1e-5
                if (
                    idx == 0
                    or idx == len(data_subset) - 1  # First/last bin
                    or abs(low * 100 % 10) < tolerance  # Multiple of 0.1 GeV
                    or abs(low * 100 % 10 - 10) < tolerance
                ):
                    labels.append(label)
                else:
                    labels.append("")  # Empty string for unlabeled bins

        return labels

    def _setup_joyplot_kwargs(
        self,
        color_list: list,
        overlap: float,
        x_limits: tuple,
        y_labels: list,
        figsize: Optional[tuple],
        user_kwargs: dict,
    ) -> dict:
        """Configure keyword arguments for joyplot."""

        # Default joyplot parameters
        default_kwargs = {
            "labels": y_labels,
            "range_style": "own",
            "x_range": list(x_limits),
            "overlap": overlap,
            "linewidth": 1.5,
            "alpha": 0.8,
            "legend": True,
            "grid": "y",
            "yrot": 0,  # No rotation for y-axis labels
            "loc": "lower right",  # Legend location
        }

        # Handle figure size
        if figsize is not None:
            default_kwargs["figsize"] = figsize
        else:
            # Auto-calculate based on number of bins and columns
            n_bins = len(y_labels)
            width = max(10, len(color_list) * 2)
            height = max(8, n_bins * 0.8)
            default_kwargs["figsize"] = (width, height)

        # Handle color assignment (avoid joypy bug with single colors)
        if len(color_list) == 1:
            default_kwargs["color"] = color_list[0]
        else:
            default_kwargs["color"] = color_list

        # Update with user-provided kwargs
        default_kwargs.update(user_kwargs)

        return default_kwargs

    def _add_joyplot_truth_lines(
        self,
        axes: list,
        truth_data: pd.DataFrame,
        columns: list,
        color_list: list,
    ) -> None:
        """Add vertical truth value lines to each subplot."""

        fit_indices = truth_data.index.tolist()

        for fit_idx, ax in zip(fit_indices, axes):
            # Get maximum y-values for each KDE curve
            y_maxes = []
            for line in ax.lines:
                if hasattr(line, "get_ydata"):
                    y_data = line.get_ydata()
                    if len(y_data) > 0:
                        y_maxes.append(y_data.max())

            # Ensure we have enough y_max values
            while len(y_maxes) < len(columns):
                y_maxes.append(1.0)  # Default fallback

            # Plot truth lines for each column
            for col, color, y_max in zip(columns, color_list, y_maxes):
                truth_value = truth_data.loc[fit_idx, col]

                ax.plot(
                    [truth_value, truth_value],
                    [0.0, y_max],
                    color=color,
                    linestyle="-",
                    linewidth=3,
                    alpha=0.9,
                    zorder=200,  # Ensure lines appear on top
                )

    def _update_joyplot_legend(self, axes: list) -> None:
        """Update legend labels to use pretty amplitude names."""

        for ax in axes:
            legend = ax.get_legend()
            if legend is not None:
                for text in legend.get_texts():
                    label = text.get_text()
                    if any(
                        label in sublist for sublist in self.coherent_sums.values()
                    ) or label in set(self.phase_differences.values()):
                        pretty_label = utils.convert_amp_name(label)
                        text.set_text(pretty_label)
