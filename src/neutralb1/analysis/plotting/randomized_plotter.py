import warnings
from typing import Dict

import matplotlib.axes
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

    def likelihood_comparison(
        self,
        fit_index: int,
        pwa_threshold: float = 0.1,
        moment_threshold: float = 0.01,
        pwa_columns: list[str] | None = None,
        moment_columns: list[str] | None = None,
        figsize: tuple = (10, 8),
    ) -> matplotlib.axes.Axes:
        """Plot likelihoods of randomized fits against their distance from the best fit.

        Args:
            fit_index (int): Index of the fit to compare against.
            pwa_threshold (float, optional): Threshold for PWA fits. Amplitude fit
                fractions must surpass this to be included. For a phase difference, the
                average of the two amplitude fit fractions must surpass this.
                Defaults to 0.1.
            moment_threshold (float, optional): Threshold for moment fits. Moment values
                must surpass this fraction of H0_0000 to be included. Defaults to 0.01.
            pwa_columns (list[str], optional): Columns to use from a PWA fit to
                to calculate the distance from the best fit. If None, all coherent sums
                and phase differences above the `pwa_threshold` are used. Defaults to
                None.
            moment_columns (list[str], optional): Columns to use from the projected
                moments to calculate the distance from the best fit. If None, all
                moments above the `moment_threshold` are used. Defaults to None.

        Raises:
            ValueError: If no significant columns are found.
            ValueError: If `moment_columns` is provided but the dataframe of projected
                moments from randomized fits is None.

        Returns:
            matplotlib.axes.Axes: The axes object containing the plot.
        """

        assert self.randomized_df is not None, "randomized_df must be provided"

        data = self._prepare_randomized_data(fit_index)

        # DETERMINE COLUMNS TO USE
        # use all coherent sums and phase differences above threshold if not given
        if pwa_columns is None:
            pwa_columns = self._find_significant_columns(
                data["best_series"],
                pwa_threshold,
            )
        # use all moments above threshold if not given
        if (
            moment_columns is None
            and data["proj_moments_series"] is not None
            and data["rand_proj_moments_df"] is not None
        ):
            moment_columns = self._find_significant_columns(
                data["proj_moments_series"],
                moment_threshold,
            )

        elif moment_columns is not None and data["rand_proj_moments_df"] is None:
            raise ValueError(
                "moment_columns provided but dataframe of projected moments from"
                " randomized fits is None."
            )
        if not pwa_columns:
            raise ValueError("No significant fit columns found above the threshold.")

        # CALCULATE THE WEIGHTED RESIDUALS
        fit_avg_squared_weighted_residuals = self._average_squared_weighted_residual(
            best_series=data["best_series"],
            rand_df=data["rand_df"],
            columns=pwa_columns,
            bootstrap_df=data["bootstrap_df"],
        )
        if (
            data["proj_moments_series"] is not None
            and data["rand_proj_moments_df"] is not None
            and moment_columns is not None
        ):
            moment_avg_squared_weighted_residuals = (
                self._average_squared_weighted_residual(
                    best_series=data["proj_moments_series"],
                    rand_df=data["rand_proj_moments_df"],
                    columns=moment_columns,
                    bootstrap_df=data["bootstrap_proj_moments_df"],
                )
            )
        else:
            moment_avg_squared_weighted_residuals = None

        likelihoods = data["rand_df"]["likelihood"]
        best_likelihood = data["best_series"]["likelihood"]
        delta_lnL = [ll - best_likelihood for ll in likelihoods]

        # Create plot based on whether moment residuals are available
        if moment_avg_squared_weighted_residuals is not None:
            # Create 2D plot when moment residuals are available
            fig, ax = plt.subplots(figsize=figsize)

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
            ax.set_ylabel("Distance from Best Result (PWA Residuals)")
        else:
            # Create 2D plot when only fit residuals are available
            fig, ax = plt.subplots(figsize=figsize)

            ax.set_xlabel(r"$\Delta(-2\ln(\mathscr{L}))$")
            ax.set_ylabel("Fit Avg Squared Error-Weighted Residuals")
            ax.scatter(
                delta_lnL,
                fit_avg_squared_weighted_residuals,
                color="black",
                alpha=0.6,
            )

        plt.tight_layout()

        return ax

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

    def _average_squared_weighted_residual(
        self,
        best_series: pd.Series,
        rand_df: pd.DataFrame,
        columns: list[str],
        bootstrap_df=None,
    ) -> list[float]:

        # Get the weighted residuals dictionary from the existing function
        weighted_residuals = self._calculate_weighted_residuals(
            best_series=best_series,
            rand_df=rand_df,
            columns=columns,
            bootstrap_df=bootstrap_df,
        )

        # Convert the dict structure to compute average squared residuals for each row
        avg_squared_weighted_residuals = []
        num_rows = len(rand_df)

        for row_idx in range(num_rows):
            # Collect all weighted residuals for this row across all columns
            row_residuals = []
            for col in columns:
                residual = weighted_residuals[col][row_idx]
                if not np.isnan(residual):  # Only include non-NaN values
                    row_residuals.append(residual)

            # Compute average of squared residuals for this row
            if row_residuals:
                squared_residuals = np.square(row_residuals)
                mean_squared = np.mean(squared_residuals)
            else:
                # If all residuals are NaN, return NaN for this row
                mean_squared = np.nan

            avg_squared_weighted_residuals.append(mean_squared)

        return avg_squared_weighted_residuals

    def _calculate_weighted_residuals(
        self,
        best_series: pd.Series,
        rand_df: pd.DataFrame,
        columns: list[str],
        bootstrap_df=None,
    ) -> Dict[str, list[float]]:
        """Compute error-weighted residuals between the best fit and randomized fits."""
        # track which column is a phase, moment, or other column (e.g. coherent sum)
        phase_indices = []
        moment_indices = []
        other_indices = []
        for i, col in enumerate(columns):
            if col in self.phase_differences:
                phase_indices.append(i)
            elif col.startswith("H") and col[1].isdigit():
                moment_indices.append(i)
            else:
                other_indices.append(i)

        if moment_indices and len(moment_indices) != len(columns):
            raise ValueError(
                "Mixture of moments and pwa results currently not handled. This"
                " function expects either all moments or all pwa results in the"
                " best_series and rand_df"
            )

        weighted_residuals = {col: [] for col in columns}

        for _, rand_series in rand_df.iterrows():
            # compute weighted residuals for this row
            all_weighted_residuals = np.full(len(columns), np.nan)

            # phases need special treatment since they are circular quantities
            if phase_indices:
                phase_columns = [columns[i] for i in phase_indices]
                phase_err_columns = [f"{columns[i]}_err" for i in phase_indices]
                rand_values = rand_series[phase_columns].to_numpy()
                best_values = best_series[phase_columns].to_numpy()

                if bootstrap_df is None:
                    err_values = best_series[phase_err_columns].to_numpy()
                else:
                    err_values = np.array(
                        [
                            np.rad2deg(
                                scipy.stats.circstd(
                                    np.deg2rad(np.abs(bootstrap_df[col])),
                                    low=0,  # type: ignore
                                    high=np.pi,
                                )
                            )
                            for col in phase_columns
                        ]
                    )

                mask = np.abs(err_values) < 1e-10
                if mask.any():
                    warnings.warn(
                        f"Very small or zero error values found in columns: "
                        f"{[col for col, is_small in zip(phase_columns, mask) if is_small]}. "
                        "This will result in NaN weighted residuals that are ignored.",
                        UserWarning,
                    )

                # Only compute residuals whose errors are not too small
                phase_weighted_residuals = np.full_like(err_values, np.nan)
                valid_indices = ~mask
                vectorized_circular_residual = np.vectorize(utils.circular_residual)
                if np.any(valid_indices):
                    phase_weighted_residuals[valid_indices] = (
                        vectorized_circular_residual(
                            np.abs(rand_values[valid_indices]),
                            np.abs(best_values[valid_indices]),
                            in_degrees=True,
                        )
                        / err_values[valid_indices]
                    )

                # Store phase residuals in the correct positions
                for i, idx in enumerate(phase_indices):
                    all_weighted_residuals[idx] = phase_weighted_residuals[i]

            # moments currently require bootstrap uncertainties, so handled separately
            if moment_indices:
                moment_columns = [columns[i] for i in moment_indices]
                rand_values = rand_series[moment_columns].to_numpy()
                best_values = best_series[moment_columns].to_numpy()

                if bootstrap_df is None:
                    if not hasattr(self, "_moment_residuals_warned"):
                        warnings.warn(
                            "Attempting to compute moment residuals without bootstrap"
                            " uncertainties, so no errors are available. Weighting will"
                            " be done with err=1.0.",
                            UserWarning,
                        )
                        self._moment_residuals_warned = True
                    err_values = np.ones_like(best_values)
                else:
                    err_values = np.array(
                        [bootstrap_df[col].std() for col in moment_columns]
                    )

                moment_weighted_residuals = (rand_values - best_values) / err_values

                # Store moment residuals in the correct positions
                for i, idx in enumerate(moment_indices):
                    all_weighted_residuals[idx] = moment_weighted_residuals[i]

            # compute all other weighted residuals normally
            if other_indices:
                other_columns = [columns[i] for i in other_indices]
                other_err_columns = [f"{columns[i]}_err" for i in other_indices]
                rand_values = rand_series[other_columns].to_numpy()
                best_values = best_series[other_columns].to_numpy()

                if bootstrap_df is None:
                    err_values = best_series[other_err_columns].to_numpy()
                else:
                    err_values = np.array(
                        [bootstrap_df[col].std() for col in other_columns]
                    )

                # Check for division by zero or very small errors
                mask = np.abs(err_values) < 1e-10
                if mask.any():
                    warnings.warn(
                        f"Very small or zero error values found in columns: "
                        f"{[col for col, is_small in zip(other_columns, mask) if is_small]}. "
                        "This will result in NaN weighted residuals that are ignored.",
                        UserWarning,
                    )

                # Only compute residuals whose errors are not too small
                other_weighted_residuals = np.full_like(err_values, np.nan)
                valid_indices = ~mask
                other_weighted_residuals[valid_indices] = (
                    rand_values[valid_indices] - best_values[valid_indices]
                ) / err_values[valid_indices]

                # Store other residuals in the correct positions
                for i, idx in enumerate(other_indices):
                    all_weighted_residuals[idx] = other_weighted_residuals[i]

            # Append the weighted residuals for this row to each column's list
            for i, col in enumerate(columns):
                weighted_residuals[col].append(all_weighted_residuals[i])

        return weighted_residuals

    def _prepare_randomized_data(self, fit_index: int) -> dict:
        """Prepare and extract relevant data for a given fit index."""

        assert self.randomized_df is not None, "randomized_df must be provided"
        # extract best fits and corresponding randomized fits for the requested index
        best_series = self._assert_series(self.fit_df.loc[fit_index])
        rand_df = self.randomized_df.loc[self.randomized_df["fit_index"] == fit_index]

        if self.bootstrap_df is not None:
            bootstrap_df = self.bootstrap_df.loc[
                self.bootstrap_df["fit_index"] == fit_index
            ]
        else:
            bootstrap_df = None

        if self.proj_moments_df is not None:
            proj_moments_series = self._assert_series(
                self.proj_moments_df.loc[self.proj_moments_df["fit_index"] == fit_index]
            )
        else:
            proj_moments_series = None

        if self.randomized_proj_moments_df is not None:
            rand_proj_moments_df = self.randomized_proj_moments_df.loc[
                self.randomized_proj_moments_df["fit_index"] == fit_index
            ]
        else:
            rand_proj_moments_df = None

        if self.bootstrap_proj_moments_df is not None:
            bootstrap_proj_moments_df = self.bootstrap_proj_moments_df.loc[
                self.bootstrap_proj_moments_df["fit_index"] == fit_index
            ]
        else:
            bootstrap_proj_moments_df = None

        return {
            "best_series": best_series,
            "rand_df": rand_df,
            "bootstrap_df": bootstrap_df,
            "proj_moments_series": proj_moments_series,
            "rand_proj_moments_df": rand_proj_moments_df,
            "bootstrap_proj_moments_df": bootstrap_proj_moments_df,
        }

    def weighted_residuals(
        self,
        fit_index: int,
        pwa_columns: list[str] | None = None,
        moment_columns: list[str] | None = None,
    ) -> matplotlib.axes.Axes:
        """Calculate and plot the weighted residuals for a given fit index and columns.

        This method is particularly useful in conjuction with
        :py:func:`likelihood_comparison` to understand the residuals that are averaged
        over in that plot. The weighted residuals are calculated as the difference
        between the randomized fit values and the best fit values, divided by the
        uncertainty (error) in the best fit values. Points are automatically grouped
        into four categories based on their delta_lnL values, with different line styles
        and alpha values for each group.

        Args:
            fit_index (int): The index of the fit to analyze.
            pwa_columns (list): The list of columns to plot from a PWA fit.
            moment_columns (list): The list of columns to plot from projected moments.

        Returns:
            matplotlib.axes.Axes: The axes object containing the plot.
        """

        assert self.randomized_df is not None, "randomized_df must be provided"

        data = self._prepare_randomized_data(fit_index)

        if pwa_columns:
            pwa_weighted_residuals = self._calculate_weighted_residuals(
                best_series=data["best_series"],
                rand_df=data["rand_df"],
                columns=pwa_columns,
                bootstrap_df=data["bootstrap_df"],
            )
        else:
            pwa_weighted_residuals = {}
            pwa_columns = []
        if moment_columns:
            if (
                data["proj_moments_series"] is not None
                and data["rand_proj_moments_df"] is not None
            ):
                moment_weighted_residuals = self._calculate_weighted_residuals(
                    best_series=data["proj_moments_series"],
                    rand_df=data["rand_proj_moments_df"],
                    columns=moment_columns,
                    bootstrap_df=None,
                )
            else:
                raise ValueError(
                    "moment_columns provided but dataframe of projected moments from"
                    " randomized fits is None."
                )
        else:
            moment_weighted_residuals = {}
            moment_columns = []

        # combine the two dictionaries
        weighted_residuals = {**pwa_weighted_residuals, **moment_weighted_residuals}
        columns = pwa_columns + moment_columns

        # Calculate delta_lnL for coloring
        likelihoods = data["rand_df"]["likelihood"]
        best_likelihood = data["best_series"]["likelihood"]
        delta_lnL = np.array([ll - best_likelihood for ll in likelihoods])

        fig, ax = plt.subplots(figsize=(12, 6))

        # Set up colormap
        cmap = plt.get_cmap("viridis_r")
        norm = Normalize(vmin=np.min(delta_lnL), vmax=np.max(delta_lnL))

        # Cluster delta_lnL values into 4 groups using quartiles
        quartiles = np.percentile(delta_lnL, [25, 50, 75])

        # Define line styles and alpha values for each group
        # Group 0 (best fits): solid line, highest alpha
        # Group 1: dotted line
        # Group 2: dashed line
        # Group 3 (worst fits): dashdot line, lowest alpha
        line_styles = ["-", ":", "--", "-."]
        alpha_values = [0.8, 0.6, 0.4, 0.2]

        # Assign each fit to a group based on delta_lnL quartiles
        groups = np.digitize(delta_lnL, quartiles)

        # Create x-axis positions for columns
        x_positions = np.arange(len(columns))

        # Plot each fit as a line
        for row_idx in range(len(data["rand_df"])):
            y_values = []
            valid_x_positions = []

            # Collect non-NaN residuals and their corresponding x positions
            for col_idx, col in enumerate(columns):
                residual = weighted_residuals[col][row_idx]
                if not np.isnan(residual):
                    y_values.append(residual)
                    valid_x_positions.append(x_positions[col_idx])

            # Plot line only if there are valid points
            if y_values:
                color = cmap(norm(delta_lnL[row_idx]))
                group = groups[row_idx]
                ax.plot(
                    valid_x_positions,
                    y_values,
                    color=color,
                    alpha=alpha_values[group],
                    linewidth=0.8,
                    linestyle=line_styles[group],
                    marker="o",
                    markersize=3,
                )

        # Customize the plot
        ax.set_xticks(x_positions)
        ax.set_xticklabels(
            [self._get_pretty_label(col) for col in columns],
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
        ax.set_ylabel("Weighted Residuals")
        ax.set_title(f"Fit Index: {fit_index}")
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color="black", linestyle="-", alpha=0.5)
        ax.axhline(y=3, color="red", linestyle="--", alpha=0.5)
        ax.axhline(y=-3, color="red", linestyle="--", alpha=0.5)

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label(r"$\Delta(-2\ln(\mathcal{L}))$")

        plt.tight_layout()

        return ax

    def likelihood_distribution(
        self, fit_index: int, figsize: tuple = (10, 8), **hist_kwargs
    ) -> matplotlib.axes.Axes:
        """Plot the distribution of likelihoods from randomized fits.

        Args:
            fit_index (int): Index of the fit to compare against.
        Returns:
            matplotlib.axes.Axes: The axes object containing the plot.
        """

        assert self.randomized_df is not None, "randomized_df must be provided"

        # extract best fits and corresponding randomized fits for the requested index
        best_series = self._assert_series(self.fit_df.loc[fit_index])
        rand_df = self.randomized_df.loc[self.randomized_df["fit_index"] == fit_index]

        likelihoods = rand_df["likelihood"]
        best_likelihood = best_series["likelihood"]
        delta_lnL = [ll - best_likelihood for ll in likelihoods]

        fig, ax = plt.subplots(figsize=figsize)

        # Compute histogram
        counts, bins = np.histogram(delta_lnL, bins=30)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])

        # Normalize delta_lnL for colormap
        cmap = plt.get_cmap("viridis_r")
        norm = Normalize(vmin=np.min(delta_lnL), vmax=np.max(delta_lnL))

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
                **hist_kwargs,
            )

        ax.set_xlabel(r"$\Delta(-2\ln(\mathscr{L}))$")
        ax.set_ylabel("Counts")

        plt.tight_layout()
        return ax

    def moment_weighted_residuals(
        self,
        fit_index: int,
        moment_threshold: float = 0.01,
        moment_columns: list[str] | None = None,
    ) -> matplotlib.axes.Axes:
        """Plot weighted residuals for moment columns that pass the threshold.

        Args:
            fit_index (int): Index of the fit to analyze.
            moment_threshold (float, optional): Threshold for moment fits. Moment values
                must surpass this fraction of H0_0000 to be included. Defaults to 0.01.
            moment_columns (list[str], optional): Columns to use from the projected
                moments. If None, all moments above the `moment_threshold` are used.
                Defaults to None.

        Returns:
            matplotlib.axes.Axes: The axes object containing the plot.
        """
        assert (
            self.randomized_proj_moments_df is not None
        ), "randomized_proj_moments_df must be provided"

        data = self._prepare_randomized_data(fit_index)

        # Determine moment columns to use
        if moment_columns is None:
            if (
                data["proj_moments_series"] is not None
                and data["rand_proj_moments_df"] is not None
            ):
                moment_columns = self._find_significant_columns(
                    data["proj_moments_series"],
                    moment_threshold,
                )
            else:
                raise ValueError(
                    "No projected moments data available and moment_columns not"
                    " provided"
                )

        return self.weighted_residuals(fit_index, moment_columns=moment_columns)

    def pwa_weighted_residuals(
        self,
        fit_index: int,
        pwa_threshold: float = 0.1,
        pwa_columns: list[str] | None = None,
    ) -> matplotlib.axes.Axes:
        """Plot weighted residuals for PWA columns that pass the threshold.

        Args:
            fit_index (int): Index of the fit to analyze.
            pwa_threshold (float, optional): Threshold for PWA fits. Amplitude fit
                fractions must surpass this to be included. For a phase difference, the
                average of the two amplitude fit fractions must surpass this.
                Defaults to 0.1.
            pwa_columns (list[str], optional): Columns to use from a PWA fit to
                calculate the weighted residuals. If None, all coherent sums
                and phase differences above the `pwa_threshold` are used. Defaults to
                None.

        Returns:
            matplotlib.axes.Axes: The axes object containing the plot.
        """
        assert self.randomized_df is not None, "randomized_df must be provided"

        data = self._prepare_randomized_data(fit_index)

        # Determine PWA columns to use
        if pwa_columns is None:
            pwa_columns = self._find_significant_columns(
                data["best_series"],
                pwa_threshold,
            )

        return self.weighted_residuals(fit_index, pwa_columns=pwa_columns)

    def randomized_summary(
        self,
        fit_index: int,
        pwa_threshold: float = 0.1,
        moment_threshold: float = 0.01,
        pwa_columns: list[str] | None = None,
        moment_columns: list[str] | None = None,
        figsize: tuple = (16, 12),
    ) -> np.ndarray:
        """Create a 2x2 summary plot with all randomized fit analysis plots.

        This method creates a comprehensive summary figure with four subplots:
        - Upper left: Likelihood distribution
        - Upper right: Moment weighted residuals (for columns passing threshold)
        - Bottom left: PWA weighted residuals (for columns passing threshold)
        - Bottom right: Likelihood comparison plot

        Args:
            fit_index (int): Index of the fit to analyze.
            pwa_threshold (float, optional): Threshold for PWA fits. Amplitude fit
                fractions must surpass this to be included. Defaults to 0.1.
            moment_threshold (float, optional): Threshold for moment fits. Moment values
                must surpass this fraction of H0_0000 to be included. Defaults to 0.01.
            pwa_columns (list[str], optional): Columns to use from a PWA fit.
                If None, all coherent sums and phase differences above the
                `pwa_threshold` are used. Defaults to None.
            moment_columns (list[str], optional): Columns to use from the projected
                moments. If None, all moments above the `moment_threshold` are used.
                Defaults to None.

        Returns:
            plt.Figure: The figure object containing all four subplots.
        """
        assert self.randomized_df is not None, "randomized_df must be provided"

        # Create 2x2 subplot figure
        fig, axes = plt.subplots(2, 2, figsize=figsize)

        # Upper left: Likelihood distribution
        plt.sca(axes[0, 0])
        self.likelihood_distribution(fit_index)
        axes[0, 0].set_title("Likelihood Distribution")

        # Upper right: Moment weighted residuals
        plt.sca(axes[0, 1])
        try:
            self.moment_weighted_residuals(fit_index, moment_threshold, moment_columns)
            axes[0, 1].set_title("Moment Weighted Residuals")
        except ValueError as e:
            axes[0, 1].text(
                0.5,
                0.5,
                f"No moment data available:\n{str(e)}",
                ha="center",
                va="center",
                transform=axes[0, 1].transAxes,
            )
            axes[0, 1].set_title("Moment Weighted Residuals (N/A)")

        # Bottom left: PWA weighted residuals
        plt.sca(axes[1, 0])
        self.pwa_weighted_residuals(fit_index, pwa_threshold, pwa_columns)
        axes[1, 0].set_title("PWA Weighted Residuals")

        # Bottom right: Likelihood comparison
        plt.sca(axes[1, 1])
        self.likelihood_comparison(
            fit_index,
            pwa_threshold,
            moment_threshold,
            pwa_columns,
            moment_columns,
        )
        axes[1, 1].set_title("Likelihood Comparison")

        plt.tight_layout()
        return axes
