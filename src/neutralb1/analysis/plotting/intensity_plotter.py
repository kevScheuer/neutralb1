import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from uncertainties import ufloat, unumpy

import neutralb1.utils as utils
from neutralb1.analysis.plotting.base_plotter import BasePWAPlotter


class IntensityPlotter(BasePWAPlotter):
    """Handles all strictly intensity-related plotting methods."""

    def jp(self, data_label: str = "Total (GlueX Phase-I)") -> matplotlib.axes.Axes:
        """Plot each JP contribution to the total intensity as a function of mass.

        Args:
            data_label (str): Label for the data in the plot. Defaults to "Total
                (GlueX Phase-I)".
        Returns:
            matplotlib.axes.Axes: The axes object for further customization.

        Todo:
            - Add option to plot with MINUIT or bootstrap errors
            - Handle very small mismatch between data widths and bar widths
        """

        # property map for consistent plotting
        colors = plt.colormaps["Dark2"].colors  # type: ignore
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

        # acceptance correct the data points using the efficiency, if needed
        if self.is_acceptance_corrected:
            efficiency = (
                self.fit_df["detected_events"] / self.fit_df["generated_events"]
            )
        else:
            efficiency = 1.0

        # plot data
        ax.errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=self.data_df["events"].div(efficiency),
            yerr=self.data_df["events_err"].div(efficiency),
            marker=".",
            linestyle="",
            color="black",
            label=data_label,
        )

        # Plot Fit Result as a grey histogram with its error bars
        num_events_col = (
            "generated_events" if self.is_acceptance_corrected else "detected_events"
        )
        ax.bar(
            self._masses,
            self.fit_df[num_events_col],
            width=self._bin_width,
            color="0.1",
            alpha=0.15,
            label="Fit Result",
        )
        ax.errorbar(
            x=self._masses,
            y=self.fit_df[num_events_col],
            yerr=self.fit_df[f"{num_events_col}_err"],
            marker=",",
            linestyle="",
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
                yerr = (
                    self.get_bootstrap_error(jp)
                    if self.bootstrap_df is not None
                    else self.fit_df[f"{jp}_err"]
                )

                ax.errorbar(
                    x=self._masses,
                    xerr=self._bin_width / 2,
                    y=self.fit_df[jp],
                    yerr=yerr,
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
        plt.tight_layout()
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
        Todo:
            - create space for legend in upper right corner
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
            layout="constrained",
        )

        # Prepare normalization data. The easiest way to handle fractional plots is to
        # convert everything to unumpy arrays and let it handle the error propagation.
        num_events_col = (
            "generated_events" if self.is_acceptance_corrected else "detected_events"
        )
        if fractional:
            total_intensity = unumpy.uarray(
                self.fit_df[num_events_col].to_numpy(),
                self.fit_df[f"{num_events_col}_err"].to_numpy(),
            )
            truth_total = (
                self.truth_df[num_events_col] if self.truth_df is not None else None
            )
        else:
            # to make plotting code simpler, define total_intensity as ones with zero
            # error when not fractional
            total_intensity = unumpy.uarray(
                np.ones_like(self.fit_df[num_events_col]),
                np.zeros_like(self.fit_df[f"{num_events_col}_err"]),
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
                    ax.set_ylabel(wave_label, fontsize=14, loc="center")

                # Plot negative reflectivity contribution
                neg_amp_name = f"m{JPmL}"
                if neg_amp_name in self.coherent_sums["eJPmL"]:
                    y_errs = (
                        self.get_bootstrap_error(neg_amp_name)
                        if self.bootstrap_df is not None
                        else self.fit_df[f"{neg_amp_name}_err"]
                    )
                    y_values = unumpy.uarray(
                        self.fit_df[neg_amp_name].to_numpy(),
                        y_errs.to_numpy(),
                    )
                    y_values = y_values / total_intensity

                    neg_plot_handle = ax.errorbar(
                        x=self._masses,
                        xerr=self._bin_width / 2,
                        y=unumpy.nominal_values(y_values),
                        yerr=unumpy.std_devs(y_values),
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
                    y_errs = (
                        self.get_bootstrap_error(pos_amp_name)
                        if self.bootstrap_df is not None
                        else self.fit_df[f"{pos_amp_name}_err"]
                    )
                    y_values = unumpy.uarray(
                        self.fit_df[pos_amp_name].to_numpy(),
                        y_errs.to_numpy(),
                    )
                    y_values = y_values / total_intensity

                    pos_plot_handle = ax.errorbar(
                        x=self._masses,
                        xerr=self._bin_width / 2,
                        y=unumpy.nominal_values(y_values),
                        yerr=unumpy.std_devs(y_values),
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

        # Set y-limit to start from zero only after all plotting is done
        # This preserves the sharey behavior when enabled
        if not sharey:
            # Only set individual limits when not sharing y-axis
            for row in range(n_rows):
                for col in range(n_cols):
                    axs[row, col].set_ylim(bottom=0)
        else:
            # For shared y-axis, set the limit on one axis per row
            # which will apply to all axes in that row
            for row in range(n_rows):
                current_bottom, current_top = axs[row, 0].get_ylim()
                axs[row, 0].set_ylim(bottom=0, top=current_top)

        # Configure figure-wide labels and legend
        y_label = (
            "Fit Fraction" if fractional else f"Events / {self._bin_width:.3f} GeV"
        )

        fig.supxlabel(rf"{self.channel} inv. mass (GeV)", fontsize=16, x=0.53)
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
                bbox_to_anchor=(1.02, 1.05),
                fontsize=12,
            )

        plt.minorticks_on()

        return axs

    def moments(self, fractional: bool = False, sharey: bool = False) -> np.ndarray:
        """Plot all the projected moments in a grid format.

        Moments cannot be organized like what is done in :meth:`waves`, so they
        are simply all plotted in a grid together.

        Args:
            fractional (bool, optional): When true, normalizes all moments to the
                H0_0000 moment. Defaults to False.
            sharey (bool, optional): When true, shares the y-axis limits across all
                moments. Defaults to False.

        Returns:
            np.ndarray: Array of axes for the moments.
        """

        assert self.proj_moments_df is not None, "No projected moments data available."

        columns = [
            col
            for col in self.proj_moments_df.columns
            if col.startswith("H") and col[1].isdigit()
        ]

        n_moments = len(columns)

        fig, axs = plt.subplots(
            nrows=int(np.ceil(n_moments / 3)),
            ncols=min(3, n_moments),
            figsize=(15, 5 * np.ceil(n_moments / 3)),
            sharex=True,
            sharey=sharey,
        )

        for ax, col in zip(axs.flatten(), columns):
            ax.set_title(col, fontsize=14)
            ax.set_xlabel(rf"{self.channel} inv. mass (GeV)", fontsize=12)

            if not fractional:
                ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

            y_values = self.proj_moments_df[col]
            y_errors = (
                self.get_bootstrap_error(col)
                if self.bootstrap_proj_moments_df is not None
                else np.zeros_like(y_values)
            )

            if fractional:
                denom = self.proj_moments_df["H0_0000"]
                denom_err = (
                    self.get_bootstrap_error("H0_0000")
                    if self.bootstrap_proj_moments_df is not None
                    else np.zeros_like(denom)
                )

                # handle error propagation
                y_values_u = unumpy.uarray(y_values, y_errors)
                denom_u = unumpy.uarray(denom, denom_err)
                frac_u = y_values_u / denom_u
                y_values = unumpy.nominal_values(frac_u)
                y_errors = unumpy.std_devs(frac_u)

            ax.errorbar(
                x=self._masses,
                xerr=self._bin_width / 2,
                y=y_values,
                yerr=y_errors,
                fmt="o",
                markersize=5,
            )

            # TODO: plot truth data if available

        # figure wide labels
        y_label = (
            f"Normalized Moment / {self._bin_width:.3f} GeV"
            if fractional
            else f"Moment / {self._bin_width:.3f} GeV"
        )

        fig.supylabel(y_label, fontsize=16)

        plt.tight_layout()
        plt.minorticks_on()

        return axs
