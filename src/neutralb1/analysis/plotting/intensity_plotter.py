import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
