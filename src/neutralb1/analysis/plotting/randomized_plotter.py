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
    """Handles plots that rely on randomized parameter fitting"""

    def randomized_summary(
        self,
        fit_index: int,
        columns: list[str],
        likelihood_threshold: float = np.inf,
        use_truth_residual: bool = True,
        figsize: tuple = (16, 12),
    ) -> matplotlib.figure.Figure:
        """Create a 2x2 summary plot of the randomized fit results

        This method creates a comprehensive summary figure with four subplots:
        - Upper left: Likelihood distribution
        - Upper right: Weighted residuals of moments passed to columns. Will only plot
            if all projected moment dataframes are provided, due to bootstrap data being
            required for calculating errors.
        - Bottom left: Amplitude, phase, parameter, or all other weighted residuals
            from columns.
        - Bottom right: Average absolute weighted residuals of the PWA and moment data.
            If moments are not provided, only PWA residuals are plotted against
            delta_lnL.

        Args:
            fit_index (int): Index of the fit to analyze.
            columns (list[str]): Columns to use in the weighted residuals plots
            likelihood_threshold (float, optional): Only randomized fits with
                [ -2ln(L_rand) - (-2ln(L_best)) ] > threshold will be included.
                Defaults to np.inf (all fits included).
            use_truth_residual (bool, optional): Whether to calculate residuals
                relative to truth values, if available. Defaults to True. If False, or
                truth values are not available, residuals are calculated relative
                to the best fit values.
            figsize (tuple, optional): Size of the entire figure. Defaults to (16, 12).

        Returns:
            plt.Figure: The figure object containing all four subplots.
        """

        data = self._prepare_randomized_data(
            fit_index, columns, likelihood_threshold, use_truth_residual
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
                data["truth_proj_moments_series"],
                use_truth_residual,
            )

        # Bottom left: PWA weighted residuals
        self._plot_weighted_residuals(
            axes[1, 0],
            data["best_series"],
            data["rand_df"],
            data["delta_lnL"],
            data["bootstrap_df"],
            data["truth_series"],
            use_truth_residual,
        )

        # Bottom right: Scatterplot of averaged absolute weighted residuals
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

        low_mass = self.data_df[self.data_df["fit_index"] == fit_index]["m_low"].values[
            0
        ]
        high_mass = self.data_df[self.data_df["fit_index"] == fit_index][
            "m_high"
        ].values[0]

        fig.suptitle(
            rf"{low_mass:.3f} < ${self.channel}$ inv. mass < {high_mass:.3f} $(GeV)$",
            fontsize=20,
        )

        return fig

    def _prepare_randomized_data(
        self,
        fit_index: int,
        columns: list[str],
        likelihood_threshold: float,
        use_truth_residual: bool,
    ) -> dict:
        """Prepare and extract relevant data for a given fit index and set of columns"""

        assert self.randomized_df is not None, "randomized_df must be provided"

        pwa_columns = []
        moment_columns = []

        for col in columns:
            if col.startswith("H") and col[1].isdigit():
                if (
                    self.proj_moments_df is not None
                    and col not in self.proj_moments_df.columns
                ):
                    raise KeyError(f"Column '{col}' not found in proj_moments_df.")
                moment_columns.append(col)
            else:
                if col not in self.fit_df.columns:
                    raise KeyError(f"Column '{col}' not found in fit_df.")
                pwa_columns.append(col)

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
        else:
            bootstrap_proj_moments_df = None

        # Extract truth data if requested and available
        if use_truth_residual:
            if self.truth_df is not None and pwa_columns:
                truth_series = self._assert_series(self.truth_df.loc[fit_index])[
                    pwa_columns
                ]
            else:
                truth_series = None

            if self.truth_proj_moments_df is not None and moment_columns:
                truth_proj_moments_series = self._assert_series(
                    self.truth_proj_moments_df.loc[
                        self.truth_proj_moments_df["fit_index"] == fit_index
                    ]
                )[moment_columns]
            else:
                truth_proj_moments_series = None
        else:
            truth_series = None
            truth_proj_moments_series = None

        return {
            "delta_lnL": delta_lnL,
            "best_series": best_series,
            "rand_df": rand_df,
            "bootstrap_df": bootstrap_df,
            "proj_moments_series": proj_moments_series,
            "rand_proj_moments_df": rand_proj_moments_df,
            "bootstrap_proj_moments_df": bootstrap_proj_moments_df,
            "truth_series": truth_series,
            "truth_proj_moments_series": truth_proj_moments_series,
        }

    def _assert_series(self, df: pd.DataFrame | pd.Series) -> pd.Series:
        """Asserts that the dataframe has exactly one row and returns it as a Series."""
        if isinstance(df, pd.DataFrame) and len(df) != 1:
            raise ValueError(f"DataFrame has {len(df)} rows, expected exactly one row.")
        series = df.squeeze()
        assert isinstance(series, pd.Series), "Result must be a Series"
        return series

    def _plot_likelihood_distribution(
        self,
        ax: matplotlib.axes.Axes,
        data: dict[str, np.ndarray | pd.Series | pd.DataFrame | None],
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
        truth_series: Optional[pd.Series] = None,
        use_truth_residual: bool = False,
    ) -> matplotlib.axes.Axes:
        """Calculate and plot the weighted residuals for randomized fits.

        The weighted residuals are calculated as the difference between the
        randomized fit values and either the best fit values or truth values
        (if use_truth_residual=True and truth_series is provided), divided by
        the uncertainty (error) in the best fit values, calculated from the
        `_err` columns or separate bootstrap samples, if provided. Points are
        automatically grouped into four categories based on their delta_lnL
        values, with different line styles and alpha values for each group.
        """

        weighted_residuals = self._calculate_weighted_residuals(
            best_series, rand_df, bootstrap_df, truth_series
        )

        # Determine reference value for ylabel
        reference = (
            "truth" if use_truth_residual and truth_series is not None else "best"
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
        ax.set_ylabel(
            rf"Weighted Residuals $(x_\text{{rand}} - x_\text{{{reference}}})"
            rf" / \sigma_\text{{best}}$"
        )

        # show 0, and 3 sigma lines
        ax.axhline(y=0, color="black", linestyle="-", alpha=0.5)
        ax.axhline(y=3, color="red", linestyle="--", alpha=0.5)
        ax.axhline(y=-3, color="red", linestyle="--", alpha=0.5)

        # make symmetric y-axis
        y_max = max(
            np.nanmax([np.nanmax(res) for res in weighted_residuals.values()]), 3
        )
        y_min = min(
            np.nanmin([np.nanmin(res) for res in weighted_residuals.values()]), -3
        )
        y_limit = max(abs(y_max), abs(y_min)) * 1.1  # add some padding
        ax.set_ylim(-y_limit, y_limit)

        ax.grid(True, alpha=0.3)

        return ax

    def _calculate_weighted_residuals(
        self,
        best_series: pd.Series,
        rand_df: pd.DataFrame,
        bootstrap_df: Optional[pd.DataFrame] = None,
        truth_series: Optional[pd.Series] = None,
    ) -> Dict[str, list[float]]:
        """Compute error-weighted residuals between randomized fits and a reference.

        The value is defined to be:
            weighted_residual = (rand_value - reference_value) / reference_error

        Where reference_value is either best_value (default) or truth_value (if
        truth_series is provided), and reference_error is always from the best fit.
        """
        col_to_weighted_residuals = {}

        for col in best_series.index:
            if col.endswith("_err"):
                continue  # skip error columns
            # phases are circular data, so flag them for special treatment
            is_phase = True if col in self.phase_differences else False

            # extract the best, rand, truth (if available), and error value for column
            best_value = best_series[col]
            rand_array = rand_df[col].to_numpy()
            reference_value = (
                truth_series[col] if truth_series is not None else best_value
            )
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
                        np.abs(reference_value),
                        in_degrees=True,
                    )
                else:
                    residuals = rand_array - reference_value

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

        # Determine reference labels based on whether truth data is available
        pwa_reference = "truth" if data["truth_series"] is not None else "best"
        moment_reference = (
            "truth" if data["truth_proj_moments_series"] is not None else "best"
        )

        pwa_residuals = self._calculate_weighted_residuals(
            data["best_series"],
            data["rand_df"],
            data["bootstrap_df"],
            data["truth_series"],
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
                data["truth_proj_moments_series"],
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

            ax.set_xlabel(
                rf"Moments: $\frac{{1}}{{N}} \sum_i^N |(x_\text{{rand}}"
                rf" - x_\text{{{moment_reference}}})"
                rf" / \sigma_\text{{best}}|_i$"
            )
            ax.set_ylabel(
                rf"PWA: $\frac{{1}}{{N}} \sum_i^N |(x_\text{{rand}}"
                rf" - x_\text{{{pwa_reference}}})"
                rf" / \sigma_\text{{best}}|_i$"
            )
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
