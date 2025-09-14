import warnings
from typing import List, Optional

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import neutralb1.utils as utils
from neutralb1.analysis.plotting.base_plotter import BasePWAPlotter


class RandomizedPlotter(BasePWAPlotter):
    """Handles plots that rely on randomized parameter fitting"""

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        proj_moments_df: Optional[pd.DataFrame] = None,
        randomized_df: Optional[pd.DataFrame] = None,
        randomized_proj_moments_df: Optional[pd.DataFrame] = None,
    ) -> None:
        """Ensure random DataFrame is not None and initialize the plotter."""
        super().__init__(
            fit_df=fit_df,
            data_df=data_df,
            proj_moments_df=proj_moments_df,
            randomized_df=randomized_df,
            randomized_proj_moments_df=randomized_proj_moments_df,
        )
        if self.randomized_df is None:
            raise ValueError("RandomPlotter requires a non-None randomized_df.")

    def likelihood_comparison(
        self,
        fit_index: int,
        contribution_threshold=0.01,
        fit_columns=None,
        moment_columns=None,
    ) -> matplotlib.axes.Axes:

        assert self.randomized_df is not None, "randomized_df must be provided"

        # extract best fits and corresponding randomized fits for the requested index
        best_series = self.fit_df.loc[fit_index].squeeze()
        if isinstance(best_series, pd.DataFrame) and len(best_series) != 1:
            raise ValueError(
                f"fit_df.loc[{fit_index}] returned {len(best_series)} rows, expected"
                " exactly one row."
            )
        assert isinstance(best_series, pd.Series), "best_series must be a Series"

        rand_df = self.randomized_df.loc[self.randomized_df["fit_index"] == fit_index]

        if self.proj_moments_df is not None:
            proj_moments_series = self.proj_moments_df.loc[
                self.proj_moments_df["fit_index"] == fit_index
            ].squeeze()
            if (
                isinstance(proj_moments_series, pd.DataFrame)
                and len(proj_moments_series) != 1
            ):
                raise ValueError(
                    f"proj_moments_df.loc[{fit_index}] returned"
                    f" {len(proj_moments_series)}"
                    " rows, expected exactly one row."
                )
            assert isinstance(
                proj_moments_series, pd.Series
            ), "proj_moments_series must be a Series"
        else:
            proj_moments_series = None

        if self.randomized_proj_moments_df is not None:
            rand_proj_moments_df = self.randomized_proj_moments_df.loc[
                self.randomized_proj_moments_df["fit_index"] == fit_index
            ]
        else:
            rand_proj_moments_df = None

        likelihoods = rand_df["likelihood"]
        best_likelihood = best_series["likelihood"]
        delta_lnL = [ll - best_likelihood for ll in likelihoods]

        # use all coherent sums and phase differences above threshold if not given
        if fit_columns is None:
            fit_columns = self._find_significant_columns(
                best_series,
                contribution_threshold,
            )

        # use all moments above threshold if not given
        if (
            moment_columns is None
            and proj_moments_series is not None
            and rand_proj_moments_df is not None
        ):
            moment_columns = self._find_significant_columns(
                proj_moments_series,
                contribution_threshold,
            )

        elif moment_columns is not None and rand_proj_moments_df is None:
            raise ValueError(
                "moment_columns provided but dataframe of projected moments from"
                " randomized fits is None."
            )

        if not fit_columns:
            raise ValueError("No significant fit columns found above the threshold.")

        fit_avg_squared_weighted_residuals = self._average_squared_weighted_residual(
            best_series=best_series,
            rand_df=rand_df,
            columns=fit_columns,
        )
        if (
            proj_moments_series is not None
            and rand_proj_moments_df is not None
            and moment_columns is not None
        ):
            moment_avg_squared_weighted_residuals = (
                self._average_squared_weighted_residual(
                    best_series=proj_moments_series,
                    rand_df=rand_proj_moments_df,
                    columns=moment_columns,
                )
            )
        else:
            moment_avg_squared_weighted_residuals = None

        # Create plot based on whether moment residuals are available
        if moment_avg_squared_weighted_residuals is not None:
            # Create 2D plot when moment residuals are available
            fig, ax = plt.subplots()

            # Use viridis colormap (colorblind friendly) with reversed order
            # so that good values (close to 0) are yellow/bright
            scatter = ax.scatter(
                moment_avg_squared_weighted_residuals,
                fit_avg_squared_weighted_residuals,
                c=delta_lnL,
                cmap="viridis_r",
                alpha=0.7,
                s=50,
            )

            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label(r"$\Delta(-2\ln(\mathcal{L}))$")

            ax.set_xlabel("Distance from Best Result (Moment Residuals)")
            ax.set_ylabel("Distance from Best Result (Fit Residuals)")
        else:
            # Create 2D plot when only fit residuals are available
            fig, ax = plt.subplots()

            ax.set_xlabel(r"$\Delta(-2\ln(\mathscr{L}))$")
            ax.set_ylabel("Fit Avg Squared Error-Weighted Residuals")
            ax.scatter(
                delta_lnL,
                fit_avg_squared_weighted_residuals,
                color="black",
                alpha=0.6,
                label="PWA Fit",
            )

        plt.tight_layout()

        return ax

    def _find_significant_columns(
        self, series: pd.Series, threshold: float
    ) -> List[str]:
        """Finds significant amplitudes, phase differences, or moments in a series."""

        significant_columns = []
        is_moment_dataframe = any(c.startswith("H0_") for c in series.index)

        if is_moment_dataframe:
            for col in series.index:
                if not col.startswith("H"):
                    continue

                moment_fraction = abs(series[col] / series["H0_0000"])
                if moment_fraction > threshold:
                    significant_columns.append(col)
        else:
            total = (
                "generated_events"
                if self.is_acceptance_corrected
                else "detected_events"
            )

            # find amplitude fit fractions above threshold
            for amp in self.coherent_sums["eJPmL"]:
                amp_fraction = series[amp] / series[total]
                if abs(amp_fraction) > threshold:
                    significant_columns.append(amp)

            for phase in self.phase_differences:
                # extract the amplitudes that make up the phase difference
                amp1, amp2 = phase.split("_")
                amp1_fraction = series[amp1] / series[total]
                amp2_fraction = series[amp2] / series[total]

                # use the average fraction of the two amplitudes to determine if phase
                # difference is significant
                if (abs(amp1_fraction) + abs(amp2_fraction)) / 2 > threshold:
                    significant_columns.append(phase)

        return significant_columns

    def _average_squared_weighted_residual(
        self,
        best_series: pd.Series,
        rand_df: pd.DataFrame,
        columns: List[str],
        bootstrap_df=None,
    ) -> List[float]:

        # separate out the phase and moment columns
        phase_columns = []
        phase_err_columns = []
        moment_columns = []
        other_columns = []
        other_err_columns = []
        for col in columns:
            if col in self.phase_differences:
                phase_columns.append(col)
                phase_err_columns.append(f"{col}_err")
            elif col.startswith("H") and col[1].isdigit():
                moment_columns.append(col)
            else:
                other_columns.append(col)
                other_err_columns.append(f"{col}_err")

        avg_squared_weighted_residuals = []
        for _, rand_series in rand_df.iterrows():
            # compute weighted residuals for this row

            # phases need special treatment since they are circular quantities
            vectorized_circular_residual = np.vectorize(utils.circular_residual)
            if phase_columns:
                rand_values = rand_series[phase_columns].to_numpy()
                best_values = best_series[phase_columns].to_numpy()
                err_values = best_series[phase_err_columns].to_numpy()
                mask = np.abs(err_values) < 1e-10
                if mask.any():
                    warnings.warn(
                        f"Very small or zero error values found in columns: "
                        f"{
                            [col for col, is_small in zip(err_values, mask) if 
                            is_small]
                        }. "
                        "This will result in NaN weighted residuals that are ignored.",
                        UserWarning,
                    )
                phase_weighted_residuals = (
                    vectorized_circular_residual(
                        rand_values,
                        best_values,
                        in_degrees=True,
                    )
                    / err_values
                )

            else:
                phase_weighted_residuals = np.array([])

            # moments currently require bootstrap uncertainties, so handled separately
            if moment_columns and bootstrap_df is None:
                if not hasattr(self, "_moment_residuals_warned"):
                    warnings.warn(
                        "Attempting to compute moment residuals without bootstrap"
                        " uncertainties, so no errors are available. Weighting will be done"
                        " with err=1.0.",
                        UserWarning,
                    )
                    self._moment_residuals_warned = True

                rand_values = rand_series[moment_columns].to_numpy()
                best_values = best_series[moment_columns].to_numpy()
                moment_weighted_residuals = (rand_values - best_values) / 1.0
            elif moment_columns and bootstrap_df is not None:
                # TODO: implement moment residuals with bootstrap uncertainties
                rand_values = rand_series[moment_columns].to_numpy()
                best_values = best_series[moment_columns].to_numpy()
                moment_weighted_residuals = (rand_values - best_values) / 1.0
                pass
            else:
                moment_weighted_residuals = np.array([])

            # compute all other weighted residuals normally
            if other_columns:
                rand_values = rand_series[other_columns].to_numpy()
                best_values = best_series[other_columns].to_numpy()
                err_values = best_series[other_err_columns].to_numpy()

                # Check for division by zero or very small errors
                mask = np.abs(err_values) < 1e-10
                if mask.any():
                    warnings.warn(
                        f"Very small or zero error values found in columns: "
                        f"{
                            [col for col, is_small in zip(other_err_columns, mask) if 
                            is_small]
                        }. "
                        "This will result in NaN weighted residuals that are ignored.",
                        UserWarning,
                    )

                other_weighted_residuals = (rand_values - best_values) / err_values

            else:
                other_weighted_residuals = np.array([])

            combined_weighted_residuals = np.concatenate(
                (
                    phase_weighted_residuals,
                    moment_weighted_residuals,
                    other_weighted_residuals,
                )
            )

            # compute the average of the squared weighted residuals
            # Ignore nan values when computing the mean
            squared_residuals = np.square(combined_weighted_residuals)
            mean_squared = np.nanmean(squared_residuals)
            avg_squared_weighted_residuals.append(mean_squared)

        return avg_squared_weighted_residuals
