import warnings
from typing import Literal, Optional, Tuple

import joypy
import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

import neutralb1.utils as utils
from neutralb1.analysis.plotting.base_plotter import BasePWAPlotter


class BootstrapPlotter(BasePWAPlotter):
    """Handles plots that rely on bootstrap results."""

    def correlation_matrix(
        self,
        columns: list[str],
        report_average: bool = True,
        fit_indices: Optional[list[int]] = None,
        method: Literal["pearson", "kendall", "spearman"] = "pearson",
        pdf_path: str = "./bootstrap_correlation_matrices.pdf",
        annot: bool = True,
        cmap: str = "RdBu_r",
    ) -> None:
        """Generate correlation matrices for bootstrap samples and save to PDF.

        Creates correlation matrix heatmaps for bootstrap samples. When
        report_average is True, averages correlations across all fit_indices and
        generates a single matrix. When False, generates a separate matrix for each
        fit_index with proper amplitude name formatting and statistical annotations.

        Args:
            columns (list[str]): Bootstrap DataFrame columns to calculate correlations
                for.
            report_average (bool): If True, average correlations across all fit_indices
                and output a single matrix. If False, output one matrix per fit_index.
                Defaults to True.
            fit_indices (list[int], optional): Specific fit indices to include. If
                None, uses all available fit indices. Defaults to None.
            method (Literal["pearson", "kendall", "spearman"]): Correlation method.
                Defaults to "pearson".
            pdf_path (str): Path to save the PDF file containing correlation
                matrix/matrices. Defaults to "./bootstrap_correlation_matrices.pdf".
            annot (bool): Whether to annotate heatmap cells with correlation values.
                Defaults to True.
            cmap (str): Colormap for the heatmap. Defaults to "RdBu_r".

        Raises:
            ValueError: If bootstrap DataFrame is None, contains missing columns, or
                specified fit_indices are not found.

        Returns:
            None: Saves plots to PDF file and prints confirmation.
        """

        assert self.bootstrap_df is not None

        # Check for missing columns
        missing_columns = [
            col for col in columns if col not in self.bootstrap_df.columns
        ]
        if missing_columns:
            raise ValueError(
                f"The following columns are missing from bootstrap_df:"
                f" {missing_columns}"
            )

        # Validate method
        valid_methods = ["pearson", "kendall", "spearman"]
        if method not in valid_methods:
            raise ValueError(f"Method must be one of {valid_methods}, got '{method}'")

        # Handle fit indices
        if fit_indices is None:
            fit_indices = self.bootstrap_df["fit_index"].unique().tolist()
        else:
            missing_indices = set(fit_indices) - set(
                self.bootstrap_df["fit_index"].unique()
            )
            if missing_indices:
                raise ValueError(
                    f"The following fit indices were not found: {missing_indices}"
                )

        # At this point fit_indices is guaranteed to be a list[int]
        assert fit_indices is not None

        # Calculate figure size based on number of columns
        # Use 0.6 inches per column/row, with minimum of 6 and maximum of 16
        n_cols = len(columns)
        base_size = max(6, min(16, n_cols * 0.6))
        figsize = (base_size, base_size)

        # Create PDF to save all plots
        with PdfPages(pdf_path) as pdf:
            if report_average:
                # Average correlations across all fit indices
                correlation_matrices = []
                for fit_index in fit_indices:
                    data = self.bootstrap_df[
                        self.bootstrap_df["fit_index"] == fit_index
                    ][columns]

                    if data.empty:
                        warnings.warn(
                            f"No data found for fit index {fit_index}", UserWarning
                        )
                        continue

                    corr_matrix = data.corr(method=method)
                    correlation_matrices.append(corr_matrix)

                if not correlation_matrices:
                    raise ValueError("No valid correlation matrices computed")

                # Average all correlation matrices using pandas
                avg_corr_matrix = (
                    pd.concat(correlation_matrices).groupby(level=0).mean()
                )

                # Create prettier labels for amplitudes and phases
                label_mapping = {}
                for col in columns:
                    if (
                        any(col in sublist for sublist in self.coherent_sums.values())
                        or col in self.phase_differences
                    ):
                        label_mapping[col] = utils.convert_amp_name(col)
                    else:
                        label_mapping[col] = col

                # Apply label mapping
                avg_corr_matrix = avg_corr_matrix.rename(
                    index=label_mapping, columns=label_mapping
                )

                # Create the plot
                plt.figure(figsize=figsize)

                # Generate heatmap
                heatmap = sns.heatmap(
                    avg_corr_matrix,
                    annot=annot,
                    cmap=cmap,
                    vmin=-1,
                    vmax=1,
                    center=0,
                    square=True,
                    fmt=".2f" if annot else "",
                    cbar_kws={"label": f"{method.capitalize()} Correlation"},
                )

                # Customize plot appearance
                plt.title(
                    (
                        f"Average {method.capitalize()} Correlation Matrix "
                        f"(n={len(correlation_matrices)} fit indices)"
                    ),
                    fontsize=14,
                    pad=20,
                )
                plt.xticks(rotation=45, ha="right")
                plt.yticks(rotation=0)

                # Add grid for better readability
                heatmap.set_facecolor("white")

                # Tight layout to prevent label cutoff
                plt.tight_layout()

                # Save to PDF
                pdf.savefig(bbox_inches="tight", dpi=300)
                plt.close()

                print(f"Correlation matrix saved to: {pdf_path}")
                print(f"Averaged across {len(correlation_matrices)} fit indices")

            else:
                # Generate separate matrix for each fit_index
                num_plots = 0
                for fit_index in fit_indices:
                    # Filter data for current directory
                    data = self.bootstrap_df[
                        self.bootstrap_df["fit_index"] == fit_index
                    ][columns]

                    if data.empty:
                        warnings.warn(
                            f"No data found for fit index {fit_index}", UserWarning
                        )
                        continue

                    # Calculate correlation matrix
                    corr_matrix = data.corr(method=method)

                    # Create prettier labels for amplitudes and phases
                    label_mapping = {}
                    for col in columns:
                        if (
                            any(
                                col in sublist
                                for sublist in self.coherent_sums.values()
                            )
                            or col in self.phase_differences
                        ):
                            label_mapping[col] = utils.convert_amp_name(col)
                        else:
                            label_mapping[col] = col

                    # Apply label mapping
                    corr_matrix = corr_matrix.rename(
                        index=label_mapping, columns=label_mapping
                    )

                    # Create the plot
                    plt.figure(figsize=figsize)

                    # Generate heatmap
                    heatmap = sns.heatmap(
                        corr_matrix,
                        annot=annot,
                        cmap=cmap,
                        vmin=-1,
                        vmax=1,
                        center=0,
                        square=True,
                        fmt=".2f" if annot else "",
                        cbar_kws={"label": f"{method.capitalize()} Correlation"},
                    )

                    # Customize plot appearance
                    plt.title(
                        (
                            f"{method.capitalize()} Correlation Matrix for fit_index"
                            f" {fit_index}"
                        ),
                        fontsize=14,
                        pad=20,
                    )
                    plt.xticks(rotation=45, ha="right")
                    plt.yticks(rotation=0)

                    # Add grid for better readability
                    heatmap.set_facecolor("white")

                    # Tight layout to prevent label cutoff
                    plt.tight_layout()

                    # Save to PDF
                    pdf.savefig(bbox_inches="tight", dpi=300)
                    plt.close()
                    num_plots += 1

                print(f"Correlation matrices saved to: {pdf_path}")
                print(f"Generated {num_plots} correlation matrix plots")

    def pairplot(
        self,
        fit_indices: list[int],
        columns: list[str],
        show_truth: bool = True,
        show_uncertainty_bands: bool = True,
        correlation_threshold: float = 0.7,
        figsize: Optional[tuple] = None,
        normalize_axis_limits: bool = False,
        **kwargs,
    ) -> sns.PairGrid:
        """Create a comprehensive pairplot matrix of bootstrap fit results.

        Generates a scatterplot matrix showing relationships between fit parameters,
        with histograms/KDE on diagonal, scatter plots on lower triangle, and 2D KDE
        contours on upper triangle. Includes uncertainty bands and truth value overlays.

        Args:
            fit_indices (list[int]): Indices of fits to include in the analysis.
                These should correspond to indices in the linked DataFrames.
            columns (list[str]): Fit parameters to plot in the matrix.
            show_truth (bool): Whether to overlay truth values if available.
                Defaults to True.
            show_uncertainty_bands (bool): Whether to show MINUIT uncertainty bands.
                Defaults to True.
            correlation_threshold (float): Threshold for highlighting strongly
                correlated parameters with black borders. Defaults to 0.7.
            figsize (tuple, optional): Figure size. If None, calculated automatically.
            **kwargs: Additional arguments passed to seaborn.PairGrid.
            normalize_axis_limits (bool): If True, all amplitude and moment subplots will
                be plotted in the range [0, 1]. Phase differences will be plotted in
                the range [-180, 180]. Defaults to False.

        Raises:
            ValueError: If bootstrap DataFrame is None or parameters are invalid.
            IndexError: If fit_indices are not found in the DataFrames.

        Returns:
            seaborn.PairGrid: The configured PairGrid object for further customization.

        Note:
            Amplitudes are automatically converted to fit fractions for visualization.
            Phase differences use circular statistics for standard deviation
                calculations.
        """

        assert self.bootstrap_df is not None

        # if we want moments, ensure we have that dataframe and flag it
        if any([c.startswith("H") and c[1].isdigit() for c in columns]):
            assert self.proj_moments_df is not None
            assert self.bootstrap_proj_moments_df is not None

        if not fit_indices:
            raise ValueError("At least one fit index must be specified.")

        if not columns:
            raise ValueError("At least one column must be specified for plotting.")

        # Validate columns exist
        missing_columns = []
        for col in columns:
            if (
                self.bootstrap_proj_moments_df is None
                and col not in self.bootstrap_df.columns
            ):
                missing_columns.append(col)
            elif (
                self.bootstrap_proj_moments_df is not None
                and col not in self.bootstrap_df.columns
                and col not in self.bootstrap_proj_moments_df.columns
            ):
                missing_columns.append(col)

        if missing_columns:
            raise ValueError(
                "The following columns are missing from the bootstrap dataframes:"
                f" {missing_columns}"
            )

        # will be normalized to always be 1, and so breaks kde
        if "H0_0000" in columns:
            columns.remove("H0_0000")

        # Verify bootstrap samples exist for all requested fit indices
        available_fit_indices = self.bootstrap_df["fit_index"].unique()
        missing_fit_indices = [
            idx for idx in fit_indices if idx not in available_fit_indices
        ]
        if self.bootstrap_proj_moments_df is not None:
            available_moment_indices = self.bootstrap_proj_moments_df[
                "fit_index"
            ].unique()
            missing_moment_indices = [
                idx for idx in fit_indices if idx not in available_moment_indices
            ]
        else:
            missing_moment_indices = []

        missing_indices = set(missing_fit_indices) | set(missing_moment_indices)
        if missing_fit_indices:
            raise IndexError(
                f"No bootstrap samples found for fit indices: {missing_indices}"
            )

        # Prepare data
        data = self._prepare_pairplot_data(columns, fit_indices, show_truth)

        # Set up the plot
        if figsize is None:
            figsize = (max(10, len(columns) * 2), max(8, len(columns) * 1.5))

        # Create PairGrid
        palette = sns.color_palette(n_colors=len(fit_indices))

        grid_kwargs = {
            "hue": "fit_index",
            "palette": palette,
            "height": figsize[0] / len(columns),
            **kwargs,
        }

        pg = sns.PairGrid(data["plot_data"], **grid_kwargs)

        # Configure plot types based on number of fits
        if len(fit_indices) < 4:
            pg.map_diag(sns.histplot, kde=True, alpha=0.7)
        else:
            pg.map_diag(sns.kdeplot, alpha=0.7)

        # Map plots to grid
        pg.map_upper(sns.kdeplot, levels=[0.003, 0.05, 0.32], alpha=0.6)
        pg.map_lower(sns.scatterplot, alpha=0.7, s=30)

        # Add overlays
        self._add_pairplot_overlays(
            pg,
            data,
            palette,
            show_truth,
            show_uncertainty_bands,
            correlation_threshold,
        )

        if normalize_axis_limits:
            self._normalize_axis_limits(pg)

        # Update labels to prettier versions
        self._update_pairplot_labels(pg)

        plt.tight_layout()

        return pg

    def _prepare_pairplot_data(
        self,
        columns: list,
        fit_indices: list,
        show_truth: bool,
    ) -> dict:
        """Prepare and normalize data for pairplot visualization.

        Many of the dataframes are type: ignore'd here because the main method has
        already verified that they're not None if needed.
        """

        total_intensity = (
            "generated_events" if self.is_acceptance_corrected else "detected_events"
        )

        # separate the fit and moment columns
        pwa_cols = []
        moment_cols = []
        for col in columns:
            if col.startswith("H") and col[1].isdigit():
                moment_cols.append(col)
            else:
                pwa_cols.append(col)

        # Extract relevant data using fit_indices directly
        fit_data = self.fit_df.loc[fit_indices][
            pwa_cols + [f"{c}_err" for c in pwa_cols] + [total_intensity]
        ].copy()
        if moment_cols:
            moment_data = self.proj_moments_df.loc[fit_indices][  # type: ignore
                moment_cols + ["H0_0000"]
            ].copy()
            fit_data = pd.concat([fit_data, moment_data], axis=1)

        # Extract relevant bootstrap data
        bootstrap_data = self.bootstrap_df[  # type: ignore
            self.bootstrap_df["fit_index"].isin(fit_indices)  # type: ignore
        ][pwa_cols + ["fit_index", total_intensity]].copy()

        if moment_cols:
            moment_bootstrap = self.bootstrap_proj_moments_df[  # type: ignore
                self.bootstrap_proj_moments_df["fit_index"].isin(  # type: ignore
                    fit_indices
                )
            ][moment_cols + ["H0_0000"]].copy()
            bootstrap_data = pd.concat([bootstrap_data, moment_bootstrap], axis=1)

        # Extract truth data if available
        truth_data = None
        if self.truth_df is not None and show_truth:
            truth_data = self.truth_df[self.truth_df.index.isin(fit_indices)][
                pwa_cols + [total_intensity]
            ].copy()
        if self.truth_proj_moments_df is not None and show_truth and moment_cols:
            truth_moments = self.truth_proj_moments_df[
                self.truth_proj_moments_df.index.isin(fit_indices)
            ][moment_cols + ["H0_0000"]].copy()
            if truth_data is not None:
                truth_data = pd.concat([truth_data, truth_moments], axis=1)
            else:
                truth_data = truth_moments

        # drop the intensity from the plot data to prevent it from being plotted
        columns_to_drop = [total_intensity]
        columns_to_drop += ["H0_0000"] if moment_cols else []
        plot_data = bootstrap_data.copy().drop(columns=columns_to_drop)
        # Convert amplitudes or moments to fit fractions
        for col in columns:
            if any(col in sublist for sublist in self.coherent_sums.values()):
                plot_data[col] = plot_data[col] / bootstrap_data[total_intensity]
            if col in moment_cols:
                plot_data[col] = plot_data[col] / bootstrap_data["H0_0000"]

        return {
            "fit_data": fit_data,
            "bootstrap_data": bootstrap_data,
            "plot_data": plot_data,
            "truth_data": truth_data,
        }

    def _add_pairplot_overlays(
        self,
        pg: sns.PairGrid,
        data: dict,
        palette,
        show_truth: bool,
        show_uncertainty_bands: bool,
        correlation_threshold: float,
    ) -> None:
        """Add uncertainty bands, truth values, and correlation highlights."""

        fit_data = data["fit_data"]
        truth_data = data["truth_data"]
        bootstrap_data = data["bootstrap_data"]

        for row in range(pg.axes.shape[0]):
            for col in range(pg.axes.shape[1]):
                ax = pg.axes[row, col]

                # Get labels. There seems to be a bug where diagonal plots labels
                # are empty, so we get them from the bottom row and first column instead
                col_label, row_label = self._get_pairgrid_xy_labels(pg, row, col)

                # Calculate scaling factors for fit fractions
                x_scaling = self._get_scaling_factor(col_label, fit_data)
                y_scaling = self._get_scaling_factor(row_label, fit_data)

                # Add uncertainty bands
                if show_uncertainty_bands:
                    self._add_uncertainty_bands(
                        ax,
                        fit_data,
                        col_label,
                        row_label,
                        x_scaling,
                        y_scaling,
                        palette,
                        row == col,
                    )

                # Add truth value markers
                if show_truth and truth_data is not None:
                    self._add_truth_markers(
                        ax,
                        truth_data,
                        bootstrap_data,
                        col_label,
                        row_label,
                        x_scaling,
                        y_scaling,
                        palette,
                        row == col,
                    )

                # Highlight strong correlations
                if row != col:
                    self._highlight_correlations(
                        ax,
                        bootstrap_data,
                        col_label,
                        row_label,
                        fit_data.index,
                        correlation_threshold,
                    )

                # Add statistical annotations on diagonal
                if row == col:
                    if row_label.startswith("H") and row_label[1].isdigit():
                        pass  # No minuit uncertainty for projected moments
                    else:
                        self._add_diagonal_annotations(ax, fit_data, col_label)

    def _get_pairgrid_xy_labels(
        self, pg: sns.PairGrid, row: int, col: int
    ) -> Tuple[str, str]:
        """Get the x and y labels for the current subplot in the PairGrid.

        There seems to be a bug where diagonal plot labels are empty, so we get them
        from the bottom row and first column instead
        """
        x_label = pg.axes[-1, col].get_xlabel()
        y_label = pg.axes[row, 0].get_ylabel()
        return x_label, y_label

    def _get_scaling_factor(self, column: str, fit_data: pd.DataFrame) -> pd.Series:
        """Get appropriate scaling factor for amplitudes and moments.

        Returns:
            pd.Series: # of events for coherent sums or amplitudes, H0_0000 for moments,
            otherwise ones.
        """
        if any(column in sublist for sublist in self.coherent_sums.values()):
            if self.is_acceptance_corrected:
                return fit_data["generated_events"]
            else:
                return fit_data["detected_events"]
        elif column.startswith("H") and column[1].isdigit():
            return fit_data["H0_0000"]
        return pd.Series(1.0, index=fit_data.index)

    def _add_uncertainty_bands(
        self,
        ax: matplotlib.axes.Axes,
        fit_data: pd.DataFrame,
        col_label: str,
        row_label: str,
        x_scaling: pd.Series,
        y_scaling: pd.Series,
        palette,
        is_diagonal: bool,
    ) -> None:
        """Add MINUIT uncertainty bands to the plot.

        Projected moments have no minuit uncertainty, so those are skipped here
        """

        for idx, color in zip(fit_data.index, palette):
            # X-direction uncertainty band
            if col_label.startswith("H") and col_label[1].isdigit():
                pass
            else:
                x_center = fit_data.loc[idx, col_label] / x_scaling.loc[idx]
                x_error = fit_data.loc[idx, f"{col_label}_err"] / x_scaling.loc[idx]

                ax.axvspan(
                    x_center - x_error, x_center + x_error, color=color, alpha=0.2
                )

            # Y-direction uncertainty band (only for off-diagonal)
            if not is_diagonal:
                if row_label.startswith("H") and row_label[1].isdigit():
                    pass
                else:
                    y_center = fit_data.loc[idx, row_label] / y_scaling.loc[idx]
                    y_error = fit_data.loc[idx, f"{row_label}_err"] / y_scaling.loc[idx]

                    ax.axhspan(
                        y_center - y_error, y_center + y_error, color=color, alpha=0.2
                    )

    def _add_truth_markers(
        self,
        ax: matplotlib.axes.Axes,
        truth_data: pd.DataFrame,
        bootstrap_data: pd.DataFrame,
        col_label: str,
        row_label: str,
        x_scaling: pd.Series,
        y_scaling: pd.Series,
        palette,
        is_diagonal: bool,
    ) -> None:
        """Add truth value markers to the plot."""

        for idx, color in zip(truth_data.index, palette):
            x_truth = truth_data.loc[idx, col_label] / x_scaling.loc[idx]

            # phases are sign ambiguous, so match the sign to the bootstrap mean
            if col_label in self.phase_differences:
                subset = np.deg2rad(
                    bootstrap_data[bootstrap_data["fit_index"] == idx][col_label]
                )
                bootstrap_mean = np.rad2deg(
                    scipy.stats.circmean(subset, low=-np.pi, high=np.pi)  # type: ignore
                )
                if np.sign(x_truth) != np.sign(bootstrap_mean):
                    x_truth = -x_truth

            if is_diagonal:
                # Vertical line for diagonal plots
                ax.axvline(
                    x=x_truth, color=color, linestyle="-", alpha=0.8, linewidth=2
                )
            else:
                # Scatter point for off-diagonal plots
                y_truth = truth_data.loc[idx, row_label] / y_scaling.loc[idx]

                # apply same sign correction for y if it's a phase
                if row_label in self.phase_differences:
                    subset = np.deg2rad(
                        bootstrap_data[bootstrap_data["fit_index"] == idx][row_label]
                    )
                    bootstrap_mean = np.rad2deg(
                        scipy.stats.circmean(
                            subset, low=-np.pi, high=np.pi  # type: ignore
                        )
                    )
                    if np.sign(y_truth) != np.sign(bootstrap_mean):
                        y_truth = -y_truth

                ax.scatter(
                    x_truth,
                    y_truth,
                    color=color,
                    marker="X",
                    s=100,
                    edgecolors="black",
                    linewidth=1,
                    zorder=100,
                )

    def _highlight_correlations(
        self,
        ax,
        bootstrap_data: pd.DataFrame,
        col_label: str,
        row_label: str,
        fit_indices: list,
        threshold: float,
    ) -> None:
        """Highlight strongly correlated parameters with black borders."""

        correlations = []
        for fit_idx in fit_indices:
            subset = bootstrap_data[bootstrap_data["fit_index"] == fit_idx]
            if len(subset) > 1:  # Need at least 2 points for correlation
                corr = subset[[col_label, row_label]].corr().iloc[0, 1]
                correlations.append(abs(corr))

        if correlations and np.mean(correlations) > threshold:
            ax.patch.set_edgecolor("black")
            ax.patch.set_linewidth(3)

    def _add_diagonal_annotations(
        self,
        ax,
        fit_data: pd.DataFrame,
        col_label: str,
    ) -> None:
        """Add bootstrap/MINUIT ratio annotations to diagonal plots."""

        ratios = []
        for fit_idx, error in zip(fit_data.index, fit_data[f"{col_label}_err"]):
            stdev = self.get_bootstrap_error(col_label)[fit_idx]
            ratios.append(stdev / error if error > 0 else 0)

        if ratios:
            avg_ratio = np.mean(ratios)
            ax.text(
                0.6,
                0.9,
                rf"$\frac{{\sigma_{{bootstrap}}}}{{\sigma_{{MINUIT}}}} = {avg_ratio:.2f}$",
                transform=ax.transAxes,
                fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7),
            )

    def _update_pairplot_labels(self, pg: sns.PairGrid) -> None:
        """Update axis labels to prettier LaTeX formatting."""

        for row in range(pg.axes.shape[0]):
            for col in range(pg.axes.shape[1]):
                ax = pg.axes[row, col]

                # Get labels. There seems to be a bug where diagonal plots labels
                # are empty, so we get them from the bottom row and first column instead
                col_label, row_label = self._get_pairgrid_xy_labels(pg, row, col)
                col_label = self._get_pretty_label(col_label)
                row_label = self._get_pretty_label(row_label)

                pg.axes[pg.axes.shape[0] - 1, col].set_xlabel(col_label, loc="center")
                pg.axes[row, 0].set_ylabel(row_label, loc="center")

    def _normalize_axis_limits(self, pg: sns.PairGrid) -> None:
        """Normalize axis limits for amplitudes, moments, and phase differences"""
        for row in range(pg.axes.shape[0]):
            for col in range(pg.axes.shape[1]):
                ax = pg.axes[row, col]
                col_label, row_label = self._get_pairgrid_xy_labels(pg, row, col)

                # normalize x axis limits
                if any(
                    col_label in sublist for sublist in self.coherent_sums.values()
                ) or (col_label.startswith("H") and col_label[1].isdigit()):
                    ax.set_xlim(0, 1)
                elif col_label in self.phase_differences:
                    ax.set_xlim(-180, 180)
                # normalize y axis limits, except on diagonal
                if row == col:
                    continue
                if any(
                    row_label in sublist for sublist in self.coherent_sums.values()
                ) or (row_label.startswith("H") and row_label[1].isdigit()):
                    ax.set_ylim(0, 1)
                elif row_label in self.phase_differences:
                    ax.set_ylim(-180, 180)

    def joyplot(
        self,
        columns: list[str],
        fit_indices: Optional[list[int]] = None,
        colormap: str | list[str] = "Accent",
        sparse_labels: bool = True,
        overlap: float = 1.0,
        show_truth: bool = True,
        truth_scaling: float = 1.0,
        figsize: Optional[tuple] = None,
        **kwargs,
    ) -> np.ndarray:
        """Create a ridge plot (joyplot) of bootstrap parameter distributions.

        Generates stacked kernel density estimates for each mass bin, providing
        an intuitive view of how the bootstrap distributions change from bin-to-bin.
        Truth values are overlaid as vertical lines when available.

        Args:
            columns (list[str]): Bootstrap DataFrame columns to plot.
            fit_indices (list[int], optional): Specific fit indices to include. If None,
                uses all available fits. Defaults to None.
            colormap (str | list[str], optional): Either a matplotlib colormap name or a
                list of color specifications (hex, RGB tuples, or named colors).
                Defaults to 'Accent'.
            sparse_labels (bool): If True, only label mass bins at multiples of 0.1 GeV.
                Defaults to True.
            overlap (float): Vertical overlap between distributions. Higher values
                create more compact plots. Defaults to 2.0.
            show_truth (bool): Whether to overlay truth values as vertical lines.
                Defaults to True.
            figsize (tuple, optional): Figure size. If None, calculated automatically.
            truth_scaling (float): Truth lines are set to match the maximum of the kde,
                but this often leads to visual clutter. This scaling factor multiplies
                the truth line height for clarity. Defaults to 1.0 (no scaling).
            **kwargs: Additional arguments passed to joypy.joyplot.

        Raises:
            ValueError: If bootstrap DataFrame is None, columns are missing, or
                colormap string is invalid.
            TypeError: If colormap is neither a string nor a list.

        Returns:
            None: Displays the joyplot.

        Warning:
            Mixing different parameter types (e.g., phases and intensities) may
            result in misleading visualizations due to different scales.
        """

        # Validate bootstrap DataFrame
        assert (
            self.bootstrap_df is not None
        ), "Bootstrap DataFrame required for joyplot analysis."

        # Validate columns
        missing_columns = []
        for col in columns:
            if col not in self.bootstrap_df.columns:
                missing_columns.append(col)
            if col not in self.fit_df.columns:
                missing_columns.append(col)
            if (
                self.truth_df is not None
                and show_truth
                and col not in self.truth_df.columns
            ):
                missing_columns.append(col)

        if missing_columns:
            raise ValueError(
                f"Missing columns in one or more DataFrames: {set(missing_columns)}"
            )

        # Handle colormap
        color_list = self._process_colormap(colormap, len(columns))

        # Determine fit indices to use
        if fit_indices is None:
            fit_indices = self.fit_df.index.tolist()
        if not fit_indices:
            raise ValueError("No valid fit files found.")

        # Extract bootstrap data
        bootstrap_data = self.bootstrap_df[
            self.bootstrap_df["fit_index"].isin(fit_indices)
        ]

        # Extract fit data
        fit_data = self.fit_df.loc[fit_indices][columns]

        # Extract truth data if available
        truth_data = None
        if self.truth_df is not None:
            truth_data = self.truth_df[self.truth_df["fit_index"].isin(fit_indices)]

        # Set up axis limits and labels
        x_limits, y_labels = self._setup_joyplot_axes(
            bootstrap_data, truth_data, columns, sparse_labels
        )

        # Configure joyplot parameters
        joy_kwargs = self._setup_joyplot_kwargs(
            color_list, overlap, x_limits, y_labels, figsize, kwargs
        )

        # Create the joyplot
        fig, axes = joypy.joyplot(
            bootstrap_data, by="fit_index", column=columns, **joy_kwargs
        )

        # Add nominal fit result as dot overlay
        self._add_joyplot_fit_point(axes, fit_data, columns, color_list)

        # Add truth value overlays
        if show_truth and truth_data is not None:
            self._add_joyplot_truth_lines(
                axes, truth_data, columns, color_list, truth_scaling
            )

        # Update legend labels to pretty format
        self._update_joyplot_legend(axes)

        plt.suptitle("Bootstrap Parameter Distributions", fontsize=16, y=1.02)
        return axes

    def _process_colormap(self, colormap: str | list[str], n_colors: int) -> list:
        """Process and validate colormap, returning list of colors.

        Args:
            colormap (str | list[str]): Either a matplotlib colormap name or a list
                of color specifications.
            n_colors (int): Number of colors needed.

        Returns:
            list: List of color specifications with length n_colors.

        Raises:
            ValueError: If colormap string is not a valid matplotlib colormap.
            TypeError: If colormap is neither a string nor a list.
        """

        # Handle list of colors
        if isinstance(colormap, list):
            if not colormap:
                raise ValueError("Color list cannot be empty.")

            # Repeat color list if too short
            if len(colormap) < n_colors:
                repeats = (n_colors // len(colormap)) + 1
                color_list = (colormap * repeats)[:n_colors]
            else:
                color_list = colormap[:n_colors]

            return color_list

        # Handle colormap string
        if isinstance(colormap, str):
            if colormap not in matplotlib.colormaps:
                available_maps = ", ".join(list(matplotlib.colormaps.keys())[:10])
                raise ValueError(
                    f"'{colormap}' is not a valid matplotlib colormap. "
                    f"Available options include: {available_maps}..."
                )

            # Get colors from colormap
            color_list = list(matplotlib.colormaps[colormap].colors)  # type: ignore

            # Repeat colormap if too short
            if len(color_list) < n_colors:
                repeats = (n_colors // len(color_list)) + 1
                color_list = (color_list * repeats)[:n_colors]
            else:
                color_list = color_list[:n_colors]

            return color_list

        # Invalid type
        raise TypeError(
            f"colormap must be either a string or a list, got {type(colormap).__name__}"
        )

    def _setup_joyplot_axes(
        self,
        bootstrap_data: pd.DataFrame,
        truth_data: Optional[pd.DataFrame],
        columns: list,
        sparse_labels: bool,
    ) -> tuple:
        """Set up x-axis limits and y-axis labels for joyplot."""

        # Calculate x-axis limits
        x_min = bootstrap_data[columns].min().min()
        x_max = bootstrap_data[columns].max().max()

        # Include truth data in range calculation
        if truth_data is not None:
            x_min = min(x_min, truth_data[columns].min().min())
            x_max = max(x_max, truth_data[columns].max().max())

        # Set sensible limits based on data type
        if all(
            any(col in sublist for sublist in self.coherent_sums.values())
            for col in columns
        ):
            # All intensity columns - start from zero
            x_min = 0.0
        elif all(col in self.phase_differences for col in columns):
            # All phase differences - use standard phase range
            x_min, x_max = -180.0, 180.0
        else:
            # Mixed or other types - add some padding
            x_range = x_max - x_min
            x_min -= 0.05 * x_range
            x_max += 0.05 * x_range

        # Create y-axis labels (mass bin ranges)
        y_labels = self._create_mass_bin_labels(bootstrap_data, sparse_labels)

        return (x_min, x_max), y_labels

    def _create_mass_bin_labels(
        self, bootstrap_data: pd.DataFrame, sparse_labels: bool
    ) -> list:
        """Create descriptive labels for mass bins."""

        # Get unique fit indices and corresponding mass ranges
        unique_indices = bootstrap_data["fit_index"].unique()
        data_subset = self.data_df[self.data_df["fit_index"].isin(unique_indices)]

        labels = []
        for idx, (low, high) in enumerate(
            zip(data_subset["m_low"], data_subset["m_high"])
        ):
            label = f"{low:.3f}-{high:.3f}"

            if not sparse_labels:
                labels.append(label)
            else:
                # Only label bins at "nice" intervals
                tolerance = 1e-1
                if (
                    idx == 0  # first bin
                    or idx == len(data_subset) - 1  # last bin
                    or abs(low * 100 % 10) < tolerance  # Multiple of 0.1 GeV
                    or abs(low * 100 % 10 - 10) < tolerance
                ):
                    labels.append(label)
                else:
                    labels.append("")  # Empty string for unlabeled bins

        return labels

    def _setup_joyplot_kwargs(
        self,
        color_list: list,
        overlap: float,
        x_limits: tuple,
        y_labels: list,
        figsize: Optional[tuple],
        user_kwargs: dict,
    ) -> dict:
        """Configure keyword arguments for joyplot."""

        # Default joyplot parameters
        default_kwargs = {
            "labels": y_labels,
            "range_style": "own",
            "x_range": list(x_limits),
            "overlap": overlap,
            "linewidth": 1.0,
            "legend": True,
            "grid": "both",
            "yrot": 0,  # No rotation for y-axis labels
            "ylim": "own",
            "loc": "lower right",  # Legend location
        }
        # reduce alpha if truth is used
        default_kwargs["alpha"] = 0.6 if self.truth_df is not None else 0.8

        # Handle figure size
        if figsize is not None:
            default_kwargs["figsize"] = figsize
        else:
            # Auto-calculate based on number of bins and columns
            n_bins = len(y_labels)
            width = max(15, len(color_list) * 2)
            height = max(10, n_bins * 0.15)
            default_kwargs["figsize"] = (width, height)

        # Handle color assignment (avoid joypy bug with single colors)
        if len(color_list) == 1:
            default_kwargs["color"] = color_list[0]
        else:
            default_kwargs["color"] = color_list

        # Update with user-provided kwargs
        default_kwargs.update(user_kwargs)

        return default_kwargs

    def _add_joyplot_fit_point(
        self,
        axes: list,
        fit_data: pd.DataFrame,
        columns: list,
        color_list: list,
    ) -> None:
        """Add nominal fit result points to each subplot."""

        fit_indices = fit_data.index.tolist()

        for fit_idx, ax in zip(fit_indices, axes):
            # Get maximum y-values for each KDE curve
            y_maxes = []
            for line in ax.lines:
                if hasattr(line, "get_ydata"):
                    y_data = line.get_ydata()
                    if len(y_data) > 0:
                        y_maxes.append(y_data.max())

            y_max = max(y_maxes) if y_maxes else 1.0

            for col, color in zip(columns, color_list):
                x1 = fit_data.loc[fit_idx, col]
                ax.plot(
                    [x1],
                    [y_max * 0.05],
                    marker="o",
                    color=color,
                    markersize=2,
                    markeredgecolor="black",
                    markeredgewidth=1,
                    zorder=150,
                )

    def _add_joyplot_truth_lines(
        self,
        axes: list,
        truth_data: pd.DataFrame,
        columns: list,
        color_list: list,
        truth_scaling: float,
    ) -> None:
        """Add vertical truth value lines to each subplot."""

        fit_indices = truth_data.index.tolist()

        for fit_idx, ax in zip(fit_indices, axes):
            # Get maximum y-values for each KDE curve
            y_maxes = []
            for line in ax.lines:
                if hasattr(line, "get_ydata"):
                    y_data = line.get_ydata()
                    if len(y_data) > 0:
                        y_maxes.append(y_data.max())

            y_max = max(y_maxes) if y_maxes else 1.0

            # Plot truth lines for each column
            zorder = 200
            for col, color in zip(columns, color_list):
                true_val = truth_data.loc[fit_idx, col]

                y_max *= truth_scaling

                ax.plot(
                    [true_val, true_val],
                    [0.0, y_max],
                    color=color,
                    linestyle="-",
                    linewidth=2,
                    alpha=0.9,
                    zorder=zorder,
                )
                ax.plot(  # black outline for visibility
                    [true_val, true_val],
                    [0.0, y_max],
                    color="black",
                    linestyle="-",
                    linewidth=4,
                    alpha=0.9,
                    zorder=zorder - 1,
                )
                zorder += 2

    def _update_joyplot_legend(self, axes: list) -> None:
        """Update legend labels to use pretty amplitude names."""

        for ax in axes:
            legend = ax.get_legend()
            if legend is not None:
                for text in legend.get_texts():
                    label = text.get_text()
                    if (
                        any(label in sublist for sublist in self.coherent_sums.values())
                        or label in self.phase_differences
                    ):
                        pretty_label = utils.convert_amp_name(label)
                        text.set_text(pretty_label)
