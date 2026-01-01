from typing import Optional

import numpy as np
import pandas as pd
import scipy.stats

import neutralb1.utils as utils


class BasePWAPlotter:
    """Base class for PWA plotting  that all sub-plotters inherit from."""

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        proj_moments_df: Optional[pd.DataFrame] = None,
        randomized_df: Optional[pd.DataFrame] = None,
        randomized_proj_moments_df: Optional[pd.DataFrame] = None,
        bootstrap_df: Optional[pd.DataFrame] = None,
        bootstrap_proj_moments_df: Optional[pd.DataFrame] = None,
        truth_df: Optional[pd.DataFrame] = None,
        truth_proj_moments_df: Optional[pd.DataFrame] = None,
        is_acceptance_corrected: bool = False,
        channel: Optional[str] = r"$\omega\pi^0$",
    ) -> None:
        """Initialize base plotter with common data and utilities.

        Args:
            fit_df (pd.DataFrame): Nominal fit results DataFrame.
            data_df (pd.DataFrame): Contains the data used for the fit.
            proj_moments_df (pd.DataFrame, optional): Contains the projected moments
                results for each nominal fit. Defaults to None.
            randomized_df (pd.DataFrame, optional): Randomized results for each
                nominal fit. Defaults to None.
            randomized_proj_moments_df (pd.DataFrame, optional): Contains the
                projected moments results for each randomized fit. Defaults to None.
            bootstrap_df (pd.DataFrame, optional): Bootstrap results for each nominal
                fit. Defaults to None.
            bootstrap_proj_moments_df (pd.DataFrame, optional): Contains the projected
                moments results for each bootstrap fit. Defaults to None.
            truth_df (pd.DataFrame, optional): Contains the ground truth values for
                the fit. Defaults to None.
            truth_proj_moments_df (pd.DataFrame, optional): Contains the ground truth
                projected moments values for the fit. Defaults to None.
            is_acceptance_corrected (bool, optional): Whether the fit is acceptance
                corrected. Defaults to False.
            channel (str, optional): The channel label to be used in plot axes. Defaults
                to r"$\\omega\\pi^0$".
        """

        # no need to copy DataFrames here, as they will not be modified
        self.fit_df = fit_df
        self.data_df = data_df
        self.randomized_df = randomized_df
        self.randomized_proj_moments_df = randomized_proj_moments_df
        self.bootstrap_df = bootstrap_df
        self.proj_moments_df = proj_moments_df
        self.bootstrap_proj_moments_df = bootstrap_proj_moments_df
        self.truth_df = truth_df
        self.truth_proj_moments_df = truth_proj_moments_df
        self.channel = channel
        self.is_acceptance_corrected = is_acceptance_corrected

        # Common properties
        self._masses = self.data_df["m_center"].astype(float).round(3)
        self._bin_width = (self.data_df["m_high"] - self.data_df["m_low"]).mean()

        # Extract coherent sums and phase differences
        self.coherent_sums = utils.get_coherent_sums(self.fit_df)
        self.phase_difference_dict = utils.get_phase_difference_dict(self.fit_df)
        self.phase_differences = utils.get_phase_differences(self.fit_df)

    def _get_pretty_label(self, label: str) -> str:
        """Convert label name to nicer LaTeX formatting."""
        if (
            any(label in sublist for sublist in self.coherent_sums.values())
            or label in self.phase_differences
        ):
            return utils.convert_amp_name(label)
        elif label.startswith("H") and label[1].isdigit():
            return utils.convert_moment_name(label)
        else:
            return label

    def get_bootstrap_error(self, label: str) -> pd.Series:
        """
        Calculate the bootstrap error (standard deviation) for the given label,
        grouped by fit index.

        Args:
            label (str): The column name for which to calculate the error.

        Returns:
            pd.Series: Standard deviation for each fit index.
        """
        if self.bootstrap_df is None:
            raise ValueError("No bootstrap DataFrame provided.")

        # Determine if label is a moment, or belongs to the PWA fit result
        if label.startswith("H") and label[1].isdigit():
            if self.bootstrap_proj_moments_df is None:
                raise KeyError(
                    f"No bootstrap projected moments DataFrame provided for moment"
                    f" label {label}."
                )
            grouped = self.bootstrap_proj_moments_df.groupby("fit_index")[label]
        else:
            if label not in self.bootstrap_df.columns:
                raise KeyError(f"Label {label} not found in bootstrap DataFrame.")
            grouped = self.bootstrap_df.groupby("fit_index")[label]

        # Use circular standard deviation for phase differences
        if label in self.phase_differences or label == "dphase":
            return grouped.apply(self._circular_std)
        else:
            return grouped.std()

    def _circular_std(self, angles: pd.Series) -> float:
        """Calculate the circular standard deviation of a series of angles.

        Implemented here due to calculation being specific to the amplitude analysis
        model used.

        Args:
            angles (pd.Series): Series of angles in radians.

        Returns:
            float: Circular standard deviation in radians.

        Note:
            Assumes data has been preprocessed to be within [-π, π] and in degrees.
        """
        angles = angles.dropna()
        if len(angles) == 0:
            return float("nan")

        angles_rad = np.deg2rad(angles)  # convert to radians
        angles_rad = np.abs(angles_rad)  # correct for sign ambiguity

        if angles_rad.max() > np.pi:
            raise ValueError(
                "Data must be within [-π, π] range in degrees. Ensure preprocessing"
                " step has been applied."
            )
        stdev = scipy.stats.circstd(angles_rad, low=0, high=np.pi)

        # convert back to degrees
        return np.rad2deg(stdev)
