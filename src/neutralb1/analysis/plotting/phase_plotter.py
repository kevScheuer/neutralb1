from typing import Any, Dict, Optional

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np

import neutralb1.utils as utils
from neutralb1.analysis.plotting.base_plotter import BasePWAPlotter


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

        phase_dif = self.phase_difference_dict[(amp1, amp2)]
        color = "red" if amp1[0] == "p" else "blue"

        fig, ax = plt.subplots()
        yerr = (
            self.get_bootstrap_error(phase_dif)
            if self.bootstrap_df is not None
            else self.fit_df[f"{phase_dif}_err"].abs()
        )
        ax.errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=self.fit_df[phase_dif],
            yerr=yerr,
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
            yerr=yerr,
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
                yerr=yerr,
                linestyle="",
                marker=".",
                color=color,
                label=utils.convert_amp_name(phase_dif),
            )
            ax.errorbar(
                x=self._masses,
                xerr=self._bin_width / 2,
                y=-(self.fit_df[phase_dif].abs()) + 360.0,
                yerr=yerr,
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

    def mass_phase(
        self,
        amp1: str,
        amp2: str,
        amp1_kwargs: Optional[Dict[str, Any]] = None,
        amp2_kwargs: Optional[Dict[str, Any]] = None,
        phase_kwargs: Optional[Dict[str, Any]] = None,
    ) -> np.ndarray:
        """Plot the mass distributions and phase difference of two amplitudes together

        Args:
            amp1 (str): first amplitude in eJPmL format
            amp2 (str): second amplitude in eJPmL format
            amp1_kwargs (Optional[Dict[str, Any]]): Additional kwargs passed to the
                amplitude axes for amp1 (axes 0). Defaults to None.
            amp2_kwargs (Optional[Dict[str, Any]]): Additional kwargs passed to the
                amplitude axes for amp2 (axes 0). Defaults to None.
            phase_kwargs (Optional[Dict[str, Any]]): Additional kwargs passed to the
                phase axes (axes 1). Defaults to None.

        Returns:
            matplotlib.axes.Axes: The axes object for further customization.
        """
        fig, axs = plt.subplots(
            2,
            1,
            sharex=True,
            gridspec_kw={"wspace": 0.0, "hspace": 0.07},
            height_ratios=[3, 1],
        )

        # plot the two amplitudes on the first subplot
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
        axs[0].errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=self.fit_df[amp1],
            yerr=amp1_err,
            marker="o",
            linestyle="",
            markersize=5,
            color="black",
            label=utils.convert_amp_name(amp1),
            **(amp1_kwargs if amp1_kwargs is not None else {}),
        )
        axs[0].errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=self.fit_df[amp2],
            yerr=amp2_err,
            marker="s",
            linestyle="",
            markersize=5,
            color="gray",
            label=utils.convert_amp_name(amp2),
            **(amp2_kwargs if amp2_kwargs is not None else {}),
        )

        if self.truth_df is not None:
            axs[0].plot(
                self._masses,
                self.truth_df[amp1],
                linestyle="-",
                marker="",
                color="black",
            )
            axs[0].plot(
                self._masses,
                self.truth_df[amp2],
                linestyle="-",
                marker="",
                color="gray",
            )

        # plot the phase difference on the second subplot
        phase_dif = self.phase_difference_dict[(amp1, amp2)]
        phase_err = (
            self.get_bootstrap_error(phase_dif)
            if self.bootstrap_df is not None
            else self.fit_df[f"{phase_dif}_err"].abs()
        )
        axs[1].errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=self.fit_df[phase_dif],
            yerr=phase_err,
            linestyle="",
            marker=".",
            color="black",
            **(phase_kwargs if phase_kwargs is not None else {}),
        )
        axs[1].errorbar(
            x=self._masses,
            xerr=self._bin_width / 2,
            y=-self.fit_df[phase_dif],
            yerr=phase_err,
            linestyle="",
            marker=".",
            color="black",
            **(phase_kwargs if phase_kwargs is not None else {}),
        )
        if self.truth_df is not None:
            axs[1].plot(
                self._masses,
                self.truth_df[phase_dif],
                linestyle="-",
                marker="",
                color="black",
            )

        axs[0].set_ylim(bottom=0.0)
        axs[0].set_ylabel(f"Events / {self._bin_width:.3f} GeV", loc="top")

        axs[1].set_yticks(np.linspace(-180, 180, 5))  # force to be in pi intervals
        axs[1].set_ylim([-180, 180])
        axs[1].set_ylabel(r"Phase Diff. ($^{\circ}$)", loc="center")
        axs[1].set_xlabel(rf"{self.channel} inv. mass $(GeV)$", loc="right")

        axs[0].legend(loc="upper right")

        return axs
