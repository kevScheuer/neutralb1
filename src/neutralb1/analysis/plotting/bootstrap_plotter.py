from typing import Literal, Optional

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
    """Handles plots that rely on bootstrap results.

    Todo:
        - All methods that use multiple batches of fits should use fit_index as the
            selection criteria, not file names.
    """

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        bootstrap_df: Optional[pd.DataFrame] = None,
        truth_df: Optional[pd.DataFrame] = None,
    ) -> None:
        """Ensure bootstrap DataFrame is not None and initialize the plotter."""
        super().__init__(
            fit_df=fit_df, data_df=data_df, bootstrap_df=bootstrap_df, truth_df=truth_df
        )
        if self.bootstrap_df is None:
            raise ValueError("BootstrapPlotter requires a non-None bootstrap_df.")

    def correlation_matrix(
        self,
        columns: list[str],
        pdf_path: str,
        method: Literal["pearson", "kendall", "spearman"] = "pearson",
        directories: Optional[list[str]] = None,
        figsize: tuple = (10, 8),
        annot: bool = True,
        cmap: str = "RdBu_r",
    ) -> None:
        """Generate correlation matrices for bootstrap samples and save to PDF.

        Creates correlation matrix heatmaps for each bootstrap sample directory,
        with proper amplitude name formatting and statistical annotations.

        Args:
            columns (list[str]): Bootstrap DataFrame columns to calculate correlations for.
            pdf_path (str): Path to save the PDF file containing all correlation matrices.
            method (str): Correlation method. Options: 'pearson', 'kendall', 'spearman'.
                Defaults to 'pearson'.
            directories (list[str], optional): Specific bootstrap directories to plot.
                If None, plots all available directories. Defaults to None.
            figsize (tuple): Figure size in inches (width, height). Defaults to (10, 8).
            annot (bool): Whether to annotate heatmap cells with correlation values.
                Defaults to True.
            cmap (str): Colormap for the heatmap. Defaults to 'RdBu_r'.

        Raises:
            ValueError: If bootstrap DataFrame is None or contains missing columns.
            FileNotFoundError: If specified directories are not found.

        Returns:
            None: Saves plots to PDF file and prints confirmation.
        """

        # Tell type checkers that self.bootstrap_df is not None
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

        # Handle directories
        available_dirs = self.bootstrap_df["directory"].unique()
        if directories is None:
            directories = available_dirs.tolist()
        else:
            missing_dirs = set(directories) - set(available_dirs)
            if missing_dirs:
                raise FileNotFoundError(
                    f"The following directories were not found: {missing_dirs}"
                )

        # At this point directories is guaranteed to be a list
        assert directories is not None

        # Create PDF to save all plots
        with PdfPages(pdf_path) as pdf:
            for directory in directories:
                # Filter data for current directory
                dir_data = self.bootstrap_df[
                    self.bootstrap_df["directory"] == directory
                ][columns]

                if dir_data.empty:
                    print(f"Warning: No data found for directory {directory}")
                    continue

                # Calculate correlation matrix
                corr_matrix = dir_data.corr(method=method)

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
                    f"{method.capitalize()} Correlation Matrix\n{directory}",
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

        print(f"Correlation matrices saved to: {pdf_path}")
        print(f"Generated {len(directories)} correlation matrix plots")

    def pairplot(
        self,
        fit_indices: list[int],
        columns: list[str],
        is_acceptance_corrected=False,
        show_truth: bool = True,
        show_uncertainty_bands: bool = True,
        correlation_threshold: float = 0.7,
        figsize: Optional[tuple] = None,
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
            is acceptance_corrected (bool): Whether the fits are acceptance corrected.
                Defaults to False, meaning amplitude fit fractions are divided by
                the number of detected events, rather than generated events.
            show_truth (bool): Whether to overlay truth values if available.
                Defaults to True.
            show_uncertainty_bands (bool): Whether to show MINUIT uncertainty bands.
                Defaults to True.
            correlation_threshold (float): Threshold for highlighting strongly
                correlated parameters with black borders. Defaults to 0.7.
            figsize (tuple, optional): Figure size. If None, calculated automatically.
            **kwargs: Additional arguments passed to seaborn.PairGrid.

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

        # Validate inputs
        if self.bootstrap_df is None:
            raise ValueError(
                "Bootstrap DataFrame required for pairplot analysis. "
                "Provide bootstrap_df during initialization."
            )

        if not fit_indices:
            raise ValueError("At least one fit index must be specified.")

        if not columns:
            raise ValueError("At least one column must be specified for plotting.")

        # Validate columns exist
        missing_columns = [
            col for col in columns if col not in self.bootstrap_df.columns
        ]
        if missing_columns:
            raise ValueError(
                "The following columns are missing from bootstrap_df:"
                f" {missing_columns}"
            )

        # Verify bootstrap samples exist for all requested fit indices
        available_indices = set(self.bootstrap_df["fit_index"].unique())
        missing_indices = [idx for idx in fit_indices if idx not in available_indices]
        if missing_indices:
            raise IndexError(
                f"No bootstrap samples found for fit indices: {missing_indices}"
            )

        # Prepare data
        data = self._prepare_pairplot_data(
            columns, fit_indices, show_truth, is_acceptance_corrected
        )

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
            columns,
            palette,
            show_truth,
            show_uncertainty_bands,
            correlation_threshold,
            is_acceptance_corrected,
        )

        # Update labels to prettier versions
        self._update_pairplot_labels(pg, columns)

        # Get fit file names for title
        fit_files = self.fit_df.loc[fit_indices, "file"].tolist()
        plt.tight_layout()

        return pg

    def _prepare_pairplot_data(
        self,
        columns: list,
        fit_indices: list,
        show_truth: bool,
        is_acceptance_corrected: bool,
    ) -> dict:
        """Prepare and normalize data for pairplot visualization."""

        total_intensity = (
            "generated_events" if is_acceptance_corrected else "detected_events"
        )

        # Extract relevant data using fit_indices directly
        fit_data = self.fit_df.loc[fit_indices][
            columns + [f"{c}_err" for c in columns] + [total_intensity]
        ].copy()

        bootstrap_data = self.bootstrap_df[  # type: ignore
            self.bootstrap_df["fit_index"].isin(fit_indices)  # type: ignore
        ][columns + ["fit_index", total_intensity]].copy()

        truth_data = None
        if self.truth_df is not None and show_truth:
            truth_data = self.truth_df[self.truth_df.index.isin(fit_indices)][
                columns + [total_intensity]
            ].copy()

        # Convert amplitudes to fit fractions
        plot_data = bootstrap_data.copy().drop(columns=[total_intensity])
        for col in columns:
            if any(col in sublist for sublist in self.coherent_sums.values()):
                plot_data[col] = plot_data[col] / bootstrap_data[total_intensity]

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
        columns: list,
        palette,
        show_truth: bool,
        show_uncertainty_bands: bool,
        correlation_threshold: float,
        is_acceptance_corrected: bool,
    ) -> None:
        """Add uncertainty bands, truth values, and correlation highlights."""

        num_plots = len(columns)
        fit_data = data["fit_data"]
        truth_data = data["truth_data"]
        bootstrap_data = data["bootstrap_data"]

        for row in range(num_plots):
            for col in range(num_plots):
                ax = pg.axes[row, col]

                # Get column labels
                col_label = columns[col]
                row_label = columns[row]

                # Calculate scaling factors for fit fractions
                x_scaling = self._get_scaling_factor(
                    col_label, fit_data, is_acceptance_corrected
                )
                y_scaling = self._get_scaling_factor(
                    row_label, fit_data, is_acceptance_corrected
                )

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
                    self._add_diagonal_annotations(
                        ax, bootstrap_data, fit_data, col_label
                    )

    def _get_scaling_factor(
        self, column: str, fit_data: pd.DataFrame, is_acceptance_corrected: bool
    ) -> pd.Series:
        """Get appropriate scaling factor for amplitude vs phase columns.

        Returns:
            pd.Series: # of events for coherent sums or amplitudes, otherwise
                ones.
        """
        if any(column in sublist for sublist in self.coherent_sums.values()):
            if is_acceptance_corrected:
                return fit_data["generated_events"]
            else:
                return fit_data["detected_events"]
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
        """Add MINUIT uncertainty bands to the plot."""

        for idx, color in zip(fit_data.index, palette):
            # X-direction uncertainty band
            x_center = fit_data.loc[idx, col_label] / x_scaling.loc[idx]
            x_error = fit_data.loc[idx, f"{col_label}_err"] / x_scaling.loc[idx]

            ax.axvspan(x_center - x_error, x_center + x_error, color=color, alpha=0.2)

            # Y-direction uncertainty band (only for off-diagonal)
            if not is_diagonal:
                y_center = fit_data.loc[idx, row_label] / y_scaling.loc[idx]
                y_error = fit_data.loc[idx, f"{row_label}_err"] / y_scaling.loc[idx]

                ax.axhspan(
                    y_center - y_error, y_center + y_error, color=color, alpha=0.2
                )

    def _add_truth_markers(
        self,
        ax,
        truth_data: pd.DataFrame,
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

            if is_diagonal:
                # Vertical line for diagonal plots
                ax.axvline(
                    x=x_truth, color=color, linestyle="-", alpha=0.8, linewidth=2
                )
            else:
                # Scatter point for off-diagonal plots
                y_truth = truth_data.loc[idx, row_label] / y_scaling.loc[idx]
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
        bootstrap_data: pd.DataFrame,
        fit_data: pd.DataFrame,
        col_label: str,
    ) -> None:
        """Add bootstrap/MINUIT ratio annotations to diagonal plots."""

        ratios = []
        for fit_idx, error in zip(fit_data.index, fit_data[f"{col_label}_err"]):
            subset = bootstrap_data[bootstrap_data["fit_index"] == fit_idx][col_label]

            if len(subset) > 1:
                # Use circular statistics for phase differences
                if col_label in self.phase_differences:
                    radian_phases = np.deg2rad(subset)
                    stdev = scipy.stats.circstd(
                        radian_phases, low=int(-np.pi), high=int(np.pi)
                    )
                    stdev = np.rad2deg(stdev)
                else:
                    stdev = subset.std()

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

    def _update_pairplot_labels(self, pg, columns: list) -> None:
        """Update axis labels to prettier LaTeX formatting."""

        num_plots = len(columns)
        for i, col in enumerate(columns):
            # Update x-axis label (bottom row)
            if hasattr(pg.axes[num_plots - 1, i], "xaxis"):
                pretty_label = self._get_pretty_label(col)
                pg.axes[num_plots - 1, i].set_xlabel(pretty_label)

            # Update y-axis label (first column)
            if hasattr(pg.axes[i, 0], "yaxis"):
                pretty_label = self._get_pretty_label(col)
                pg.axes[i, 0].set_ylabel(pretty_label)

    def _get_pretty_label(self, column: str) -> str:
        """Convert column name to pretty LaTeX formatting."""
        if (
            any(column in sublist for sublist in self.coherent_sums.values())
            or column in self.phase_differences
        ):
            return utils.convert_amp_name(column)
        return column

    def joyplot(
        self,
        columns: list[str],
        fits: Optional[list[str]] = None,
        colormap: Optional[str] = None,
        sparse_labels: bool = True,
        overlap: float = 2.0,
        show_truth: bool = True,
        figsize: Optional[tuple] = None,
        **kwargs,
    ) -> np.ndarray:
        """Create a ridge plot (joyplot) of bootstrap parameter distributions.

        Generates stacked kernel density estimates for each mass bin, providing
        an intuitive view of how the bootstrap distributions change from bin-to-bin.
        Truth values are overlaid as vertical lines when available.

        Args:
            columns (list[str]): Bootstrap DataFrame columns to plot.
            fits (list[str], optional): Specific fit files to include. If None,
                uses all available fits. Defaults to None.
            colormap (str, optional): Matplotlib colormap name. Defaults to 'Accent'.
            sparse_labels (bool): If True, only label mass bins at multiples of 0.1 GeV.
                Defaults to True.
            overlap (float): Vertical overlap between distributions. Higher values
                create more compact plots. Defaults to 2.0.
            show_truth (bool): Whether to overlay truth values as vertical lines.
                Defaults to True.
            figsize (tuple, optional): Figure size. If None, calculated automatically.
            **kwargs: Additional arguments passed to joypy.joyplot.

        Raises:
            ValueError: If bootstrap DataFrame is None or columns are missing.
            TypeError: If colormap is not a valid string.

        Returns:
            None: Displays the joyplot.

        Warning:
            Mixing different parameter types (e.g., phases and intensities) may
            result in misleading visualizations due to different scales.
        """

        # Validate bootstrap DataFrame
        if self.bootstrap_df is None:
            raise ValueError("Bootstrap DataFrame required for joyplot analysis.")

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

        # Prepare data
        plot_data = self._prepare_joyplot_data(fits)

        # Set up axis limits and labels
        x_limits, y_labels = self._setup_joyplot_axes(
            plot_data["bootstrap_data"], plot_data["truth_data"], columns, sparse_labels
        )

        # Configure joyplot parameters
        joy_kwargs = self._setup_joyplot_kwargs(
            color_list, overlap, x_limits, y_labels, figsize, kwargs
        )

        # Create the joyplot
        fig, axes = joypy.joyplot(
            plot_data["bootstrap_data"], by="fit_index", column=columns, **joy_kwargs
        )

        # Add truth value overlays
        if show_truth and plot_data["truth_data"] is not None:
            self._add_joyplot_truth_lines(
                axes, plot_data["truth_data"], columns, color_list
            )

        # Update legend labels to pretty format
        self._update_joyplot_legend(axes)

        plt.suptitle("Bootstrap Parameter Distributions", fontsize=16, y=1.02)
        return axes

    def _process_colormap(self, colormap: Optional[str], n_colors: int) -> list:
        """Process and validate colormap, returning list of colors."""

        if colormap is None:
            colormap = "Accent"

        if not isinstance(colormap, str):
            raise TypeError("Colormap must be a string.")

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

    def _prepare_joyplot_data(self, fits: Optional[list]) -> dict:
        """Prepare data for joyplot, handling fit selection and truth values."""

        # Determine fit indices to use
        if fits is None:
            fit_indices = self.fit_df.index.tolist()
        else:
            fit_indices = self.fit_df[self.fit_df["file"].isin(fits)].index.tolist()

        if not fit_indices:
            raise ValueError("No valid fit files found.")

        # Extract bootstrap data (already validated above)
        bootstrap_data = self.bootstrap_df[  # type: ignore
            self.bootstrap_df["fit_index"].isin(fit_indices)  # type: ignore
        ].copy()

        # Extract truth data if available
        truth_data = None
        if self.truth_df is not None:
            truth_data = self.truth_df[self.truth_df.index.isin(fit_indices)].copy()

        return {
            "bootstrap_data": bootstrap_data,
            "truth_data": truth_data,
            "fit_indices": fit_indices,
        }

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
                tolerance = 1e-5
                if (
                    idx == 0
                    or idx == len(data_subset) - 1  # First/last bin
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
            "linewidth": 1.5,
            "alpha": 0.8,
            "legend": True,
            "grid": "y",
            "yrot": 0,  # No rotation for y-axis labels
            "loc": "lower right",  # Legend location
        }

        # Handle figure size
        if figsize is not None:
            default_kwargs["figsize"] = figsize
        else:
            # Auto-calculate based on number of bins and columns
            n_bins = len(y_labels)
            width = max(10, len(color_list) * 2)
            height = max(8, n_bins * 0.8)
            default_kwargs["figsize"] = (width, height)

        # Handle color assignment (avoid joypy bug with single colors)
        if len(color_list) == 1:
            default_kwargs["color"] = color_list[0]
        else:
            default_kwargs["color"] = color_list

        # Update with user-provided kwargs
        default_kwargs.update(user_kwargs)

        return default_kwargs

    def _add_joyplot_truth_lines(
        self,
        axes: list,
        truth_data: pd.DataFrame,
        columns: list,
        color_list: list,
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

            # Ensure we have enough y_max values
            while len(y_maxes) < len(columns):
                y_maxes.append(1.0)  # Default fallback

            # Plot truth lines for each column
            for col, color, y_max in zip(columns, color_list, y_maxes):
                truth_value = truth_data.loc[fit_idx, col]

                ax.plot(
                    [truth_value, truth_value],
                    [0.0, y_max],
                    color=color,
                    linestyle="-",
                    linewidth=3,
                    alpha=0.9,
                    zorder=200,  # Ensure lines appear on top
                )

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
