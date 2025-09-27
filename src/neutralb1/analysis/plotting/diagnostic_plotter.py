import warnings
from typing import Tuple

import matplotlib.axes
import matplotlib.container
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import neutralb1.utils as utils
from neutralb1.analysis.plotting.base_plotter import BasePWAPlotter


class DiagnosticPlotter(BasePWAPlotter):
    """Handles more complex diagnostic plots.

    Todo:
        - Add shapiro test for bootstrap distributions, like in february study
    """

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
            layout="constrained",
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

        yerr = (
            self.get_bootstrap_error("dsratio")
            if self.bootstrap_df is not None
            else self.fit_df["dsratio_err"]
        )

        ax_ratio.errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=self.fit_df["dsratio"],
            yerr=yerr,
            color="black",
            linestyle="",
            marker=".",
            capsize=2,
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

        yerr = (
            self.get_bootstrap_error("dphase")
            if self.bootstrap_df is not None
            else self.fit_df["dphase_err"].abs()
        )
        ax_phase.errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=self.fit_df["dphase"],
            yerr=yerr,
            color="black",
            linestyle="",
            marker=".",
            capsize=2,
        )
        ax_phase.errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=-self.fit_df["dphase"],
            yerr=yerr,
            color="black",
            linestyle="",
            marker=".",
            capsize=2,
        )

        # Plot correlation
        if "cov_dsratio_dphase" in self.fit_df.columns:
            color = "tab:blue" if self.bootstrap_df is not None else "black"
            correlation = (
                self.fit_df["cov_dsratio_dphase"]
                / (self.fit_df["dsratio_err"] * self.fit_df["dphase_err"])
            ).fillna(0)

            ax_corr.plot(
                self._masses,
                correlation,
                color=color,
                linestyle="",
                marker=".",
                markersize=4,
                label="MINUIT Correlation",
            )

        if self.bootstrap_df is not None:
            color = (
                "tab:orange" if "cov_dsratio_dphase" in self.fit_df.columns else "black"
            )
            bootstrap_correlation = (
                self.bootstrap_df.groupby("fit_index")
                .corr()
                .unstack()[("dsratio", "dphase")]
            )
            ax_corr.plot(
                self._masses,
                bootstrap_correlation,
                color=color,
                linestyle="",
                marker=".",
                markersize=4,
                label="Bootstrap Correlation",
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
        ax_corr.axhline(y=0, color="black", linestyle="-")
        ax_corr.legend()

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
        """Plot D/S ratio calculated from individual amplitudes.

        Todo:
            - Implement bootstrap errors for ratio calculation.
            - Very messy and difficult to read plot. Needs improvement, though still
                helpful as it visualizes correlation as a function of mass. Would be
                nice to compare bootstrap correlation value in each bin to visual
                correlation across bins.
        """

        # Marker and color maps for different quantum numbers
        marker_map = {"q": "h", "p": "^", "0": ".", "m": "v", "n": "s"}
        color_map = {"p": "red", "m": "blue"}

        fig = plt.figure(figsize=figsize, dpi=300, layout="constrained")
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
            if (d_wave, s_wave) in self.phase_difference_dict:
                phase_dif_name = self.phase_difference_dict[(d_wave, s_wave)]
            elif (s_wave, d_wave) in self.phase_difference_dict:
                phase_dif_name = self.phase_difference_dict[(s_wave, d_wave)]
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
            np.multiply(np.square(d_err / safe_d), 0.25)
            + np.multiply(np.square(s_err / safe_s), 0.25)
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

        Todo:
            - implement bootstrap errors if available. Copy (or extract in separate
                method) phase difference handling from likelihood_comparison
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

        neg_handle, pos_handle = None, None

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
                    neg_handle, pos_handle = self._plot_diagonal_intensity(
                        ax,
                        JPmL_row,
                        max_intensity,
                        intensity_ticks,
                        intensity_labels,
                        row,
                    )

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
            0.53, 0.04, rf"{self.channel} inv. mass (GeV)", ha="center", fontsize=14
        )
        fig.text(
            0.03,
            0.55,
            r"Phase Differences ($^{\circ}$)",
            ha="center",
            va="center",
            fontsize=14,
            rotation="vertical",
        )

        # Add legend if we have plot handles
        legend_handles = []
        if pos_handle is not None:
            legend_handles.append(pos_handle)
        if neg_handle is not None:
            legend_handles.append(neg_handle)

        if legend_handles:
            fig.legend(
                handles=legend_handles,
                loc="upper right",
                bbox_to_anchor=(1.0, 1.08),
                fontsize=12,
            )

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.1, left=0.1)

        return axs

    def _plot_diagonal_intensity(
        self,
        ax: matplotlib.axes.Axes,
        JPmL: str,
        max_intensity: float,
        intensity_ticks: np.ndarray,
        intensity_labels: list,
        row: int,
    ) -> Tuple[
        matplotlib.container.ErrorbarContainer | None,
        matplotlib.container.ErrorbarContainer | None,
    ]:
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
            ax.set_ylabel(title, fontsize=10, loc="center")
            ax.set_yticks(intensity_ticks, intensity_labels)
        else:
            # Remove y-tick labels for other diagonal elements to save space
            ax.set_yticks(intensity_ticks, [""] * len(intensity_ticks))

        neg_handle, pos_handle = None, None

        # Plot negative reflectivity
        neg_amp = f"m{JPmL}"
        if neg_amp in self.fit_df.columns:
            yerr = (
                self.get_bootstrap_error(neg_amp)
                if self.bootstrap_df is not None
                else self.fit_df[f"{neg_amp}_err"]
            )
            neg_handle = ax.errorbar(
                x=self._masses,
                xerr=self._bin_width / 2,
                y=self.fit_df[neg_amp],
                yerr=yerr,
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
            yerr = (
                self.get_bootstrap_error(pos_amp)
                if self.bootstrap_df is not None
                else self.fit_df[f"{pos_amp}_err"]
            )
            pos_handle = ax.errorbar(
                x=self._masses,
                xerr=self._bin_width / 2,
                y=self.fit_df[pos_amp],
                yerr=yerr,
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

        return neg_handle, pos_handle

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
        if (amp1, amp2) in self.phase_difference_dict:
            phase_dif = self.phase_difference_dict[(amp1, amp2)]
        else:
            return  # Skip if phase difference not available

        # Plot phase values and their sign ambiguity
        phase_values = self.fit_df[phase_dif]
        phase_errors = (
            self.get_bootstrap_error(phase_dif)
            if self.bootstrap_df is not None
            else self.fit_df[f"{phase_dif}_err"].abs()
        )

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
            ax.set_ylabel(title, fontsize=10, loc="center")
            ax.set_yticks(phase_ticks)
        else:
            ax.set_yticks(phase_ticks, [""] * len(phase_ticks))

        # Get phase difference data
        amp1, amp2 = f"m{JPmL_row}", f"m{JPmL_col}"
        if (amp1, amp2) in self.phase_difference_dict:
            phase_dif = self.phase_difference_dict[(amp1, amp2)]
        else:
            return  # Skip if phase difference not available

        # Plot phase values and their sign ambiguity
        phase_values = self.fit_df[phase_dif]
        phase_errors = (
            self.get_bootstrap_error(phase_dif)
            if self.bootstrap_df is not None
            else self.fit_df[f"{phase_dif}_err"].abs()
        )

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

    def compare_results(
        self, fit1: pd.Series, fit2: pd.Series, columns: list, figsize=(15, 6)
    ) -> matplotlib.axes.Axes:
        """Compare two fit results by plotting the weighted residuals between them.

        The weighted residual is defined as (fit1 - fit2) / sqrt(err1^2 + err2^2).
        The plot also provides an averaged sum of the squared weighted residuals across
        the columns. This is similar to a chi2 / ndf, but important to distinguish as
        we don't know the number of degrees of freedom here, and so it cannot be
        strictly interpreted as such.

        Args:
            fit1_df (pd.Series): Series containing the first fit results.
            fit2_df (pd.Series): Series containing the second fit results.
            columns (list): List of column names to compare between the two fits.
            figsize (tuple): Figure size in inches (width, height). Defaults to (15, 6).
        """

        # complain if columns are not in both dataframes
        missing_cols = []
        for col in columns:
            if col not in fit1.index or col not in fit2.index:
                missing_cols.append(col)
        if missing_cols:
            warnings.warn(
                f"The following columns are missing in fit results, only those that exist"
                f" in both will be compared: {missing_cols}",
                UserWarning,
            )
            columns.remove(missing_cols)

        weighted_residuals = []
        for col in columns:
            try:
                val1, err1 = float(fit1[col]), float(fit1[f"{col}_err"])  # type: ignore
                val2, err2 = float(fit2[col]), float(fit2[f"{col}_err"])  # type: ignore
            except (ValueError, TypeError):
                warnings.warn(
                    f"Non-numeric data found in column {col}, skipping comparison.",
                    UserWarning,
                )
                continue

            # Avoid division by zero
            total_err = np.sqrt(err1**2 + err2**2)
            total_err = np.where(total_err == 0, np.finfo(float).eps, total_err)

            if col in self.phase_differences:
                # For phase differences, account for circular data
                w_res = utils.circular_residual(val1, val2) / total_err
            else:
                w_res = (val1 - val2) / total_err

            weighted_residuals.append(w_res)

        avg_squared_weighted_residuals = np.mean(np.square(weighted_residuals))

        fig, ax = plt.subplots(figsize=figsize)

        ax.bar(x=columns, height=weighted_residuals, color="darkblue", width=0.8)

        ax.set_xlabel("Fit Parameters")
        ax.set_ylabel("Weighted Residuals")
        ax.axhline(0, color="black", linestyle="-", linewidth=1.0)

        # lines at +/- 3 to denote 3 sigma threshold
        ax.axhline(3, color="red", linestyle="--", linewidth=1.0)
        ax.axhline(-3, color="red", linestyle="--", linewidth=1.0)

        # improve label readability
        ax.set_xticks(range(len(columns)))
        pretty_labels = []
        for label in columns:
            if (
                any(label in sublist for sublist in self.coherent_sums.values())
                or label in self.phase_differences
            ):
                pretty_labels.append(utils.convert_amp_name(label))
            elif label.startswith("H"):
                pretty_labels.append(utils.convert_moment_name(label))
            else:
                pretty_labels.append(label)

        ax.set_xticklabels(
            pretty_labels, rotation=45, ha="right", rotation_mode="anchor"
        )

        # put value of avg_squared_weighted_residuals in the corner
        ax.text(
            1.0,
            0.9,
            rf"$\frac{{1}}{{N}}\sum_i^N w^2_i = {avg_squared_weighted_residuals:.2f}$",
            horizontalalignment="right",
            fontsize=12,
            transform=ax.transAxes,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=2),
        )

        plt.tight_layout()
        return ax
