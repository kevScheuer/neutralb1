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

        Todo:
            - Add option to plot with MINUIT or bootstrap errors
        """

        phase_dif = self.phase_difference_dict[(amp1, amp2)]
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
        phase_dif = self.phase_difference_dict[(amp1, amp2)]
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
