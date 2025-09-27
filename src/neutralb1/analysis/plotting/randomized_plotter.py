import warnings
from typing import Dict, Optional

import matplotlib.axes
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from matplotlib.colors import Normalize

import neutralb1.utils as utils
from neutralb1.analysis.plotting.base_plotter import BasePWAPlotter


class RandomizedPlotter(BasePWAPlotter):
    """Handles plots that rely on randomized parameter fitting

    Todo:
        - The weighted_residuals function currently assumes that either all columns
            are moments or PWA results, and using it is very clunky.
    """

    def randomized_summary(
        self,
        fit_index: int,
        likelihood_threshold: float = np.inf,
        pwa_threshold: float = 0.1,
        moment_threshold: float = 0.01,
        columns: Optional[list[str]] = None,
        return_significant_columns: bool = False,
        figsize: tuple = (16, 12),
    ) -> matplotlib.figure.Figure | tuple[matplotlib.figure.Figure, list[str]]:
        """Create a 2x2 summary plot with all randomized fit analysis plots.

        This method creates a comprehensive summary figure with four subplots:
        - Upper left: Likelihood distribution
        - Upper right: Moment weighted residuals (for columns passing threshold).
            Will only plot if all projected moment dataframes are provided, due to
            bootstrap data being required for calculating errors.
        - Bottom left: PWA weighted residuals (for columns passing threshold)
        - Bottom right: Sum of the average squared weighted residuals of the PWA and
            moments, vs. delta_lnL

        Args:
            fit_index (int): Index of the fit to analyze.
            likelihood_threshold (float, optional): Only randomized fits with
                delta(-2lnL) below this threshold will be included in the plots.
                Defaults to np.inf (all fits included).
            pwa_threshold (float, optional): Threshold for PWA fits if columns not
                provided. Amplitude fit fractions, or average of the fit fractions that
                make up a phase difference must surpass this to be included. Defaults to
                0.1.
            moment_threshold (float, optional): Threshold for moment fits if columns not
                provided. Moment values must surpass this fraction of H0_0000 to be
                included. Defaults to 0.01.
            columns (Optional[list[str]]): Explicit columns to use in the weighted
                residuals plots. If None (default), significant columns above the
                thresholds are used.
            return_significant_columns (bool, optional): If True, the significant
                columns used in the weighted residuals plots are returned. Defaults to
                False.
            figsize (tuple, optional): Size of the entire figure. Defaults to (16, 12).

        Returns:
            plt.Figure: The figure object containing all four subplots.
        """

        data = self._prepare_randomized_data(
            fit_index, columns, likelihood_threshold, pwa_threshold, moment_threshold
        )

        # Create 2x2 subplot figure
        fig, axes = plt.subplots(2, 2, figsize=figsize, layout="constrained")

        # Upper left: Likelihood distribution
        self._plot_likelihood_distribution(axes[0, 0], data)

        # Upper right: Moment weighted residuals
        if (
            data["proj_moments_series"] is None
            or data["rand_proj_moments_df"] is None
            or data["bootstrap_proj_moments_df"] is None
        ):
            warnings.warn(
                "Moment residuals require best, rand, and bootstrap projected moment"
                " dataframes. Skipping moment weighted residuals plot.",
                UserWarning,
            )
            axes[0, 1].text(
                0.5,
                0.5,
                f"Moment data not provided",
                ha="center",
                va="center",
                transform=axes[0, 1].transAxes,
            )
        else:
            self._plot_weighted_residuals(
                axes[0, 1],
                data["proj_moments_series"],
                data["rand_proj_moments_df"],
                data["delta_lnL"],
                data["bootstrap_proj_moments_df"],
            )

        # Bottom left: PWA weighted residuals
        self._plot_weighted_residuals(
            axes[1, 0],
            data["best_series"],
            data["rand_df"],
            data["delta_lnL"],
            data["bootstrap_df"],
        )

        # Bottom right: Scatterplot of averaged squared weighted residuals
        self._scatterplot(axes[1, 1], data)

        # Add colorbar to the right of the entire 2x2 grid
        cmap = plt.get_cmap("cividis")
        norm = Normalize(vmin=np.min(data["delta_lnL"]), vmax=np.max(data["delta_lnL"]))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(
            sm,
            ax=axes,
            location="right",
            shrink=0.8,
            pad=0.02,
            label=r"$\Delta(-2\ln(\mathcal{L}))$",
        )

        if return_significant_columns:
            return fig, (
                data["best_series"].index.tolist()
                + data["proj_moments_series"].index.tolist()
                if data["proj_moments_series"] is not None
                else []
            )

        return fig

    def _prepare_randomized_data(
        self,
        fit_index: int,
        columns: list[str] | None,
        likelihood_threshold: float,
        pwa_threshold: float,
        moment_threshold: float,
    ) -> dict:
        """Prepare and extract relevant data for a given fit index and set of columns"""

        assert self.randomized_df is not None, "randomized_df must be provided"

        pwa_columns = []
        moment_columns = []

        if columns is None:  # uses columns above thresholds
            pwa_columns = self._find_significant_columns(
                self._assert_series(self.fit_df.loc[fit_index]),
                threshold=pwa_threshold,
            )
            if (
                self.proj_moments_df is not None
                and self.randomized_proj_moments_df is not None
            ):
                moment_columns = self._find_significant_columns(
                    self._assert_series(
                        self.proj_moments_df.loc[
                            self.proj_moments_df["fit_index"] == fit_index
                        ]
                    ),
                    threshold=moment_threshold,
                )
        else:  # explicit columns provided, separate into pwa and moment columns
            for col in columns:
                if col not in self.fit_df.columns and (
                    self.proj_moments_df is not None
                    and col not in self.proj_moments_df.columns
                ):
                    raise KeyError(
                        f"Column '{col}' not found in either fit_df or proj_moments_df."
                    )
                if col.startswith("H") and col[1].isdigit():
                    moment_columns.append(col)
                else:
                    pwa_columns.append(col)
            pwa_columns = [
                col for col in columns if not (col.startswith("H") and col[1].isdigit())
            ]
            moment_columns = [
                col for col in columns if col.startswith("H") and col[1].isdigit()
            ]

            if moment_columns and (
                self.proj_moments_df is None or self.randomized_proj_moments_df is None
            ):
                raise ValueError(
                    "moment_columns provided but dataframes of projected moments from"
                    " fits are None."
                )

        # extract best fits and corresponding randomized fits for the requested index
        best_series = self._assert_series(self.fit_df.loc[fit_index])
        rand_df = self.randomized_df.loc[self.randomized_df["fit_index"] == fit_index]

        # before we limit the rand_df to specific columns, extract delta_lnL values
        # and error columns (projected moments don't have errors)
        delta_lnL = np.array(
            [ll - best_series["likelihood"] for ll in rand_df["likelihood"]]
        )

        # filter based on likelihood threshold
        likelihood_mask = delta_lnL <= likelihood_threshold
        delta_lnL = delta_lnL[likelihood_mask]
        rand_df = rand_df.iloc[likelihood_mask]

        pwa_error_columns = [f"{col}_err" for col in pwa_columns]
        pwa_columns += pwa_error_columns

        # now limit the dataframes to just the columns we'll be using
        best_series = best_series[pwa_columns]
        rand_df = rand_df[pwa_columns]

        # check if the other dataframes exist, and do the same selection
        if self.bootstrap_df is not None:
            bootstrap_df = self.bootstrap_df.loc[
                self.bootstrap_df["fit_index"] == fit_index
            ][pwa_columns]
            bootstrap_df = bootstrap_df.iloc[likelihood_mask]
        else:
            bootstrap_df = None

        if self.proj_moments_df is not None:
            proj_moments_series = self._assert_series(
                self.proj_moments_df.loc[self.proj_moments_df["fit_index"] == fit_index]
            )[moment_columns]
        else:
            proj_moments_series = None

        if self.randomized_proj_moments_df is not None:
            rand_proj_moments_df = self.randomized_proj_moments_df.loc[
                self.randomized_proj_moments_df["fit_index"] == fit_index
            ][moment_columns]
            rand_proj_moments_df = rand_proj_moments_df.iloc[likelihood_mask]
        else:
            rand_proj_moments_df = None

        if self.bootstrap_proj_moments_df is not None:
            bootstrap_proj_moments_df = self.bootstrap_proj_moments_df.loc[
                self.bootstrap_proj_moments_df["fit_index"] == fit_index
            ][moment_columns]
            bootstrap_proj_moments_df = bootstrap_proj_moments_df.iloc[likelihood_mask]
        else:
            bootstrap_proj_moments_df = None

        return {
            "delta_lnL": delta_lnL,
            "best_series": best_series,
            "rand_df": rand_df,
            "bootstrap_df": bootstrap_df,
            "proj_moments_series": proj_moments_series,
            "rand_proj_moments_df": rand_proj_moments_df,
            "bootstrap_proj_moments_df": bootstrap_proj_moments_df,
        }

    def _assert_series(self, df: pd.DataFrame | pd.Series) -> pd.Series:
        """Asserts that the dataframe has exactly one row and returns it as a Series."""
        if isinstance(df, pd.DataFrame) and len(df) != 1:
            raise ValueError(f"DataFrame has {len(df)} rows, expected exactly one row.")
        series = df.squeeze()
        assert isinstance(series, pd.Series), "Result must be a Series"
        return series

    def _find_significant_columns(
        self, series: pd.Series, threshold: float
    ) -> list[str]:
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

    def _plot_likelihood_distribution(
        self, ax: matplotlib.axes.Axes, data: dict[str, pd.Series | pd.DataFrame | None]
    ) -> matplotlib.axes.Axes:
        """Plot the distribution of likelihoods from randomized fits."""

        # assertions for type checkers
        assert data["delta_lnL"] is not None, "delta_lnL must be provided"
        assert data["best_series"] is not None, "best fit must be provided"
        assert data["rand_df"] is not None, "randomized fits must be provided"

        # Compute histogram
        counts, bins = np.histogram(data["delta_lnL"], bins=100)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])

        # Normalize delta_lnL for colormap
        cmap = plt.get_cmap("cividis")
        norm = Normalize(vmin=np.min(data["delta_lnL"]), vmax=np.max(data["delta_lnL"]))

        # Draw each bar with its color
        for i in range(len(counts)):
            color = cmap(norm(bin_centers[i]))
            ax.bar(
                bin_centers[i],
                counts[i],
                width=(bins[1] - bins[0]),
                align="center",
                color=color,
                edgecolor="black",
                alpha=1.0,
            )

        ax.set_xlabel(r"$\Delta(-2\ln(\mathscr{L}))$")
        ax.set_ylabel("Counts")

        # Ensure left hand side of plot is at 0 (with padding)
        x_min, x_max = ax.get_xlim()
        pad_x = max(0.05, 0.05 * (x_max - x_min))
        ax.set_xlim(left=0 - pad_x)

        return ax

    def _plot_weighted_residuals(
        self,
        ax: matplotlib.axes.Axes,
        best_series: pd.Series,
        rand_df: pd.DataFrame,
        delta_lnL: np.ndarray,
        bootstrap_df: Optional[pd.DataFrame] = None,
    ) -> matplotlib.axes.Axes:
        """Calculate and plot the weighted residuals for randomized fits.

        The weighted residuals are calculated as the difference
        between the randomized fit values and the best fit values, divided by the
        uncertainty (error) in the best fit values, calculated from the `_err` columns
        or separate bootstrap samples, if provided. Points are automatically grouped
        into four categories based on their delta_lnL values, with different line styles
        and alpha values for each group.
        """

        weighted_residuals = self._calculate_weighted_residuals(
            best_series, rand_df, bootstrap_df
        )

        # remove columns who are purely NaN
        all_nan_columns = [
            col
            for col in weighted_residuals.keys()
            if np.all(np.isnan(weighted_residuals[col]))
        ]
        for col in all_nan_columns:
            warnings.warn(
                f"All weighted residuals for column '{col}' are NaN. This may be"
                " due to missing or zero errors in the best fit. Column will not"
                " be plotted.",
                UserWarning,
            )
            del weighted_residuals[col]

        # warn user if any columns contain NaN values
        contain_nan_columns = [
            col
            for col in weighted_residuals.keys()
            if np.any(np.isnan(weighted_residuals[col]))
        ]
        for col in contain_nan_columns:
            warnings.warn(
                f"Some weighted residuals for column '{col}' are NaN. This may be"
                " due to missing or zero errors in the best fit. These points will"
                " not be plotted.",
                UserWarning,
            )

        # Set up colormap
        cmap = plt.get_cmap("cividis")
        norm = Normalize(vmin=np.min(delta_lnL), vmax=np.max(delta_lnL))

        # Cluster delta_lnL values into 4 groups using quartiles
        quartiles = np.percentile(delta_lnL, [25, 50, 75])

        # Define line styles and alpha values for each group (quartile)
        # Group 0 (best fits): solid line, highest alpha
        # Group 1: dotted line
        # Group 2: dashed line
        # Group 3 (worst fits): dashdot line, lowest alpha
        line_styles = ["-", ":", "--", "-."]
        alpha_values = [1.0, 0.8, 0.6, 0.5]

        # Assign each fit to a group based on delta_lnL quartiles
        groups = np.digitize(delta_lnL, quartiles)

        # Create x-axis positions for columns
        columns = list(weighted_residuals.keys())
        values = list(weighted_residuals.values())

        for i in range(len(values[0])):
            y = [weighted_residuals[col][i] for col in columns]
            group = groups[i]
            color = cmap(norm(delta_lnL[i]))
            ax.plot(
                columns,
                y,
                color=color,
                alpha=alpha_values[group],
                linewidth=0.8,
                linestyle=line_styles[group],
                marker="o",
                markersize=3,
            )

        # Customize the plot
        ax.set_xticks(range(len(columns)))
        ax.set_xticklabels(
            [self._get_pretty_label(col) for col in columns],
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )

        # show 0, and 3 sigma lines
        ax.axhline(y=0, color="black", linestyle="-", alpha=0.5)
        ax.axhline(y=3, color="red", linestyle="--", alpha=0.5)
        ax.axhline(y=-3, color="red", linestyle="--", alpha=0.5)

        # make symmetric y-axis
        y_max = np.nanmax([np.nanmax(res) for res in weighted_residuals.values()])
        y_min = np.nanmin([np.nanmin(res) for res in weighted_residuals.values()])
        y_limit = max(abs(y_max), abs(y_min)) * 1.1  # add some padding
        ax.set_ylim(-y_limit, y_limit)

        ax.grid(True, alpha=0.3)

        return ax

    def _calculate_weighted_residuals(
        self,
        best_series: pd.Series,
        rand_df: pd.DataFrame,
        bootstrap_df: Optional[pd.DataFrame] = None,
    ) -> Dict[str, list[float]]:
        """Compute error-weighted residuals between the best fit and randomized fits.

        The value is defined to be:
            weighted_residual = (rand_value - best_value) / best_value_error
        """
        col_to_weighted_residuals = {}

        for col in best_series.index:
            if col.endswith("_err"):
                continue  # skip error columns
            # phases are circular data, so flag them for special treatment
            is_phase = True if col in self.phase_differences else False

            # extract the best, rand, and error value for this column
            best_value = best_series[col]
            rand_array = rand_df[col].to_numpy()
            if bootstrap_df is None:
                best_err = best_series.get(f"{col}_err", np.nan)
            else:
                if is_phase:
                    best_err = np.rad2deg(  # convert back to degrees
                        scipy.stats.circstd(  # circstd for circular data
                            np.deg2rad(  # convert to radians for circstd
                                np.abs(bootstrap_df[col])  # abs due to sign ambiguity
                            ),
                            low=0,
                            high=np.pi,
                        )
                    )
                else:
                    best_err = bootstrap_df[col].std()

            # calculate result, and avoid division by zero or very small errors
            if np.isnan(best_err) or abs(best_err) < 1e-10 or best_err == 0.0:
                weighted_residuals = np.full_like(rand_array, np.nan)
            else:
                if is_phase:
                    vectorized_circular_residual = np.vectorize(utils.circular_residual)
                    # again, abs due to sign ambiguity
                    residuals = vectorized_circular_residual(
                        np.abs(rand_array),
                        np.abs(best_value),
                        in_degrees=True,
                    )
                else:
                    residuals = rand_array - best_value

                weighted_residuals = residuals / best_err

            col_to_weighted_residuals[col] = weighted_residuals

        return col_to_weighted_residuals

    def _scatterplot(
        self,
        ax: matplotlib.axes.Axes,
        data: dict,
    ) -> matplotlib.axes.Axes:
        """Create a scatterplot of the averaged squared weighted residuals.

        If no moment data is provided, only the PWA residuals are plotted against
        delta_lnL. Otherwise, the moment residuals are plotted on the x-axis and
        the PWA residuals on the y-axis. Points are always colored by delta_lnL.
        """

        pwa_residuals = self._calculate_weighted_residuals(
            data["best_series"], data["rand_df"], data["bootstrap_df"]
        )
        # user has already been warned about NaN columns in _plot_weighted_residuals,
        # just remove or ignore them here
        pwa_residuals = {
            k: v for k, v in pwa_residuals.items() if not np.all(np.isnan(v))
        }
        values = np.array(list(pwa_residuals.values())).T  # shape = (N_fits, N_columns)
        pwa_averages = [
            np.nanmean(np.abs(values[row])) for row in range(values.shape[0])
        ]

        if (
            data["proj_moments_series"] is not None
            and data["rand_proj_moments_df"] is not None
            and data["bootstrap_proj_moments_df"] is not None
        ):
            # create 2D scatterplot with moment residuals on x-axis and PWA on y-axis

            # repeat process for moment residuals
            moment_residuals = self._calculate_weighted_residuals(
                data["proj_moments_series"],
                data["rand_proj_moments_df"],
                data["bootstrap_proj_moments_df"],
            )
            moment_residuals = {
                k: v for k, v in moment_residuals.items() if not np.all(np.isnan(v))
            }
            values = np.array(list(moment_residuals.values())).T
            moment_averages = [
                np.nanmean(np.abs(values[row])) for row in range(values.shape[0])
            ]

            # Create 2D scatterplot with moments on x-axis and PWA on y-axis
            ax.scatter(
                moment_averages,
                pwa_averages,
                c=data["delta_lnL"],
                cmap="cividis",
                alpha=0.6,
                s=50,
            )

            ax.set_xlabel("Distance (Moment Residuals)")
            ax.set_ylabel("Distance (PWA Residuals)")
        else:
            # Create 1D scatterplot with PWA residuals on y-axis and delta_lnL on x-axis
            ax.scatter(
                data["delta_lnL"],
                pwa_averages,
                c=data["delta_lnL"],
                cmap="cividis",
                alpha=0.6,
                s=50,
            )

            ax.set_xlabel(r"$\Delta(-2\ln(\mathscr{L}))$")
            ax.set_ylabel("Distance (PWA Residuals)")

        # Ensure lower left corner is at (0, 0) with some padding
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()

        if np.isclose(x_max, 0.0, atol=0.1):
            x_max = 0 + 0.1  # ensure some range if all x are 0
        if np.isclose(y_max, 0.0, atol=0.1):
            y_max = 0 + 0.1  # ensure some range if all y are 0

        pad_x = max(0.05, 0.05 * (x_max - x_min))
        pad_y = max(0.05, 0.05 * (y_max - y_min))
        ax.set_xlim(0 - pad_x, x_max)
        ax.set_ylim(0 - pad_y, y_max)

        return ax
