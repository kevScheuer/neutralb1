from typing import Any, Dict, Optional

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy

import neutralb1.utils as utils
from neutralb1.analysis.plotting.base_plotter import BasePWAPlotter


class PhasePlotter(BasePWAPlotter):
    """Handles all strictly phase-related plotting methods."""

    def phase(
        self,
        amp1: str,
        amp2: str,
        extend_range: bool = True,
        phase_kwargs: Optional[Dict[str, Any]] = None,
        truth_kwargs: Optional[Dict[str, Any]] = None,
        ax: Optional[matplotlib.axes.Axes] = None,
    ) -> matplotlib.axes.Axes:
        """Plot the phase difference between two amplitudes as a function of mass.

        Order of amplitudes does not matter

        Args:
            amp1 (str): Name of the first amplitude in eJPmL format.
            amp2 (str): Name of the second amplitude in eJPmL format.
            extend_range (bool): Extends the y-axis from [-pi, pi) to [-2pi, 2pi) bounds
                if True, and plots mirrored transparaent values. Useful for revealing
                phase motions near the boundaries. Defaults to True.
            phase_kwargs (Optional[Dict[str, Any]]): Additional kwargs passed to the
                main phase difference errorbar plot. Defaults to None.
            truth_kwargs (Optional[Dict[str, Any]]): Additional kwargs passed to the
                truth phase plot. Defaults to None.
            ax (Optional[matplotlib.axes.Axes]): directly provide axes object to plot
                on.

        Returns:
            matplotlib.axes.Axes: The axes object for further customization.
        """

        phase_dif = self.phase_difference_dict[(amp1, amp2)]

        if ax is None:
            fig, ax = plt.subplots()

        yerr = (
            self.get_bootstrap_error(phase_dif)
            if self.bootstrap_df is not None
            else self.fit_df[f"{phase_dif}_err"].abs()
        )

        # Set default parameters for main phase difference and allow kwargs to override
        color = "red" if amp1[0] == "p" else "blue"
        phase_params = {
            "x": self._masses,
            "xerr": self._bin_width / 2,
            "y": self.fit_df[phase_dif],
            "yerr": yerr,
            "linestyle": "",
            "marker": ".",
            "markersize": 3,
            "color": color,
            "alpha": 0.5,
            "label": utils.convert_amp_name(phase_dif),
        }
        if self.truth_df is not None:
            # make data more transparent if truth is also plotted, but allow override
            phase_params["alpha"] = 0.3
        if phase_kwargs is not None:
            phase_params.update(phase_kwargs)
        handle = ax.errorbar(**phase_params)

        # plot the negative as well (natural sign ambiguity in the model)
        phase_neg_params = phase_params.copy()
        phase_neg_params["y"] = -self.fit_df[phase_dif]
        phase_neg_params.pop("label", None)
        ax.errorbar(**phase_neg_params)

        # plot the truth phase if available
        if self.truth_df is not None:
            truth_params = {
                "linestyle": "-",
                "marker": "",
                "color": color,
            }
            if truth_kwargs is not None:
                truth_params.update(truth_kwargs)
            ax.plot(
                self._masses,
                self.truth_df[phase_dif],
                linestyle=truth_params["linestyle"],
                marker=truth_params["marker"],
                color=truth_params["color"],
            )

        if extend_range:
            # add or subtract ambiguous 360 degree shifts. Positive values get shifted
            # down by 2pi, and negative values up by 2pi.
            phase_params_ext1 = phase_params.copy()
            phase_params_ext1["y"] = self.fit_df[phase_dif].abs() - 360.0
            phase_params_ext1["alpha"] = phase_params_ext1["alpha"] / 2.0
            phase_params_ext1.pop("label", None)
            ax.errorbar(**phase_params_ext1)

            phase_params_ext2 = phase_params.copy()
            phase_params_ext2["y"] = -(self.fit_df[phase_dif].abs()) + 360.0
            phase_params_ext2["alpha"] = phase_params_ext2["alpha"] / 2.0
            phase_params_ext2.pop("label", None)
            ax.errorbar(**phase_params_ext2)

        if phase_params.get("label") != "":
            ax.legend()
        if extend_range:
            ax.set_ylim(-360.0, 360.0)
            ax.set_yticks(np.linspace(-360, 360, 9))  # force to be in pi intervals
        else:
            ax.set_ylim(-180.0, 180.0)
            ax.set_yticks(np.linspace(-180, 180, 5))  # force to be in pi intervals
        ax.set_ylabel(r"Phase Diff. ($^{\circ}$)", loc="top")
        ax.set_xlabel(rf"{self.channel} inv. mass $(GeV)$", loc="right")

        return ax

    def mass_phase(
        self,
        amp1: str,
        amp2: str,
        fractional: bool = False,
        amp1_kwargs: Optional[Dict[str, Any]] = None,
        amp2_kwargs: Optional[Dict[str, Any]] = None,
        phase_kwargs: Optional[Dict[str, Any]] = None,
        amp_ax: Optional[matplotlib.axes.Axes] = None,
        phase_ax: Optional[matplotlib.axes.Axes] = None,
    ) -> np.ndarray:
        """Plot the mass distributions and phase difference of two amplitudes together

        Args:
            amp1 (str): first amplitude in eJPmL format
            amp2 (str): second amplitude in eJPmL format
            fractional (bool): If True, scales all values by dividing by the total
                intensity in each bin. Defaults to False.
            amp1_kwargs (Optional[Dict[str, Any]]): Additional kwargs passed to the
                amplitude axes for amp1 (axes 0). Defaults to None.
            amp2_kwargs (Optional[Dict[str, Any]]): Additional kwargs passed to the
                amplitude axes for amp2 (axes 0). Defaults to None.
            phase_kwargs (Optional[Dict[str, Any]]): Additional kwargs passed to the
                phase axes (axes 1). Defaults to None.
            amp_ax (Optional[matplotlib.axes.Axes]): directly provide axes object to
                plot amplitudes on. If None, a new figure and axes will be created.
                Defaults to None.
            phase_ax (Optional[matplotlib.axes.Axes]): directly provide axes object to
                plot phase on. If None, a new figure and axes will be created.
                Defaults to None.

        Returns:
            matplotlib.axes.Axes: The axes object for further customization.
        """
        if amp_ax is None and phase_ax is None:
            fig, axs = plt.subplots(
                2,
                1,
                sharex=True,
                gridspec_kw={"wspace": 0.0, "hspace": 0.07},
                height_ratios=[3, 1],
            )
        elif (amp_ax is not None and phase_ax is None) or (
            amp_ax is None and phase_ax is not None
        ):
            raise ValueError(
                "Both amp_ax and phase_ax must be provided, or neither to create new."
            )
        else:
            axs = np.array([amp_ax, phase_ax])

        # Prepare normalization data if fractional plotting is requested
        num_events_col = (
            "generated_events" if self.is_acceptance_corrected else "detected_events"
        )
        if fractional:
            total_intensity = unumpy.uarray(
                self.fit_df[num_events_col].to_numpy(),
                self.fit_df[f"{num_events_col}_err"].to_numpy(),
            )
        else:
            # to make plotting code simpler, define total_intensity as ones with zero
            # error when not fractional
            total_intensity = unumpy.uarray(
                np.ones_like(self.fit_df[num_events_col]),
                np.zeros_like(self.fit_df[f"{num_events_col}_err"]),
            )

        # get errors
        amp1_err = (
            self.get_bootstrap_error(amp1)
            if self.bootstrap_df is not None
            else self.fit_df[f"{amp1}_err"]
        )
        amp2_err = (
            self.get_bootstrap_error(amp2)
            if self.bootstrap_df is not None
            else self.fit_df[f"{amp2}_err"]
        )

        # Apply normalization if fractional
        amp1_values = unumpy.uarray(
            self.fit_df[amp1].to_numpy(),
            amp1_err.to_numpy(),
        )
        amp1_values = amp1_values / total_intensity

        amp2_values = unumpy.uarray(
            self.fit_df[amp2].to_numpy(),
            amp2_err.to_numpy(),
        )
        amp2_values = amp2_values / total_intensity

        # Set default parameters for amp1 and allow kwargs to override
        amp1_params = {
            "x": self._masses,
            "xerr": self._bin_width / 2,
            "y": unumpy.nominal_values(amp1_values),
            "yerr": unumpy.std_devs(amp1_values),
            "marker": "o",
            "linestyle": "",
            "markersize": 5,
            "color": "black",
            "alpha": 1.0,
            "label": utils.convert_amp_name(amp1),
        }
        if self.truth_df is not None:
            # make data more transparent if truth is also plotted, but allow override
            amp1_params["alpha"] = 0.5
        if amp1_kwargs is not None:  # override defaults
            amp1_params.update(amp1_kwargs)
        axs[0].errorbar(**amp1_params)

        # Set default parameters for amp2 and allow kwargs to override
        amp2_params = {
            "x": self._masses,
            "xerr": self._bin_width / 2,
            "y": unumpy.nominal_values(amp2_values),
            "yerr": unumpy.std_devs(amp2_values),
            "marker": "s",
            "linestyle": "",
            "markersize": 5,
            "color": "gray",
            "alpha": 1.0,
            "label": utils.convert_amp_name(amp2),
        }
        if self.truth_df is not None:
            # make data more transparent if truth is also plotted, but allow override
            amp2_params["alpha"] = 0.5
        if amp2_kwargs is not None:  # override defaults
            amp2_params.update(amp2_kwargs)
        axs[0].errorbar(**amp2_params)

        if self.truth_df is not None:
            truth_amp1 = self.truth_df[amp1]
            truth_amp2 = self.truth_df[amp2]
            if fractional:
                truth_total = self.truth_df[num_events_col]

                truth_amp1 = truth_amp1 / truth_total
                truth_amp2 = truth_amp2 / truth_total

            axs[0].plot(
                self._masses,
                truth_amp1,
                linestyle="-",
                marker="",
                color=amp1_params["color"],
            )
            axs[0].plot(
                self._masses,
                truth_amp2,
                linestyle="-",
                marker="",
                color=amp2_params["color"],
            )

        # plot the phase difference on the second subplot
        phase_dif = self.phase_difference_dict[(amp1, amp2)]
        phase_err = (
            self.get_bootstrap_error(phase_dif)
            if self.bootstrap_df is not None
            else self.fit_df[f"{phase_dif}_err"].abs()
        )
        # Set default parameters for phase difference and allow kwargs to override
        phase_params = {
            "x": self._masses,
            "xerr": self._bin_width / 2,
            "y": self.fit_df[phase_dif],
            "yerr": phase_err,
            "linestyle": "",
            "marker": ".",
            "color": "black",
            "alpha": 1.0,
        }
        if self.truth_df is not None:
            # make data more transparent if truth is also plotted, but allow override
            phase_params["alpha"] = 0.5
        if phase_kwargs is not None:
            phase_params.update(phase_kwargs)  # override defaults
        axs[1].errorbar(**phase_params)

        # Also plot the negative phase with the same parameters (except for y values)
        phase_params_neg = phase_params.copy()
        phase_params_neg["y"] = -self.fit_df[phase_dif]
        axs[1].errorbar(**phase_params_neg)
        if self.truth_df is not None:
            axs[1].plot(
                self._masses,
                self.truth_df[phase_dif],
                linestyle="-",
                marker="",
                color=phase_params["color"],
            )

        axs[0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        axs[0].set_ylim(bottom=0.0)
        y_label = (
            f"Fit Fraction / {self._bin_width:.3f} GeV"
            if fractional
            else f"Events / {self._bin_width:.3f} GeV"
        )
        axs[0].set_ylabel(y_label, loc="top")

        axs[1].set_yticks(np.linspace(-180, 180, 5))  # force to be in pi intervals
        axs[1].set_ylim([-180, 180])
        axs[1].set_ylabel(r"Phase Diff. ($^{\circ}$)", loc="center")
        axs[1].set_xlabel(rf"${self.channel}$ inv. mass $(GeV)$", loc="right")

        if axs[0].get_legend_handles_labels()[0]:
            axs[0].legend(loc="upper right")

        return axs
