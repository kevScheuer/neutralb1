from typing import Optional

import pandas as pd

import neutralb1.utils as utils


class BasePWAPlotter:
    """Base class for PWA plotting  that all sub-plotters inherit from."""

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        randomized_df: Optional[pd.DataFrame] = None,
        randomized_proj_moments_df: Optional[pd.DataFrame] = None,
        bootstrap_df: Optional[pd.DataFrame] = None,
        proj_moments_df: Optional[pd.DataFrame] = None,
        bootstrap_proj_moments_df: Optional[pd.DataFrame] = None,
        truth_df: Optional[pd.DataFrame] = None,
        channel: Optional[str] = r"$\omega\pi^0$",
    ) -> None:
        """Initialize base plotter with common data and utilities.

        Args:
            fit_df (pd.DataFrame): Nominal fit results DataFrame.
            data_df (pd.DataFrame): Contains the data used for the fit.
            randomized_df (pd.DataFrame, optional): Randomized results for each
                nominal fit. Defaults to None.
            randomized_proj_moments_df (pd.DataFrame, optional): Contains the
                projected moments results for each randomized fit. Defaults to None.
            bootstrap_df (pd.DataFrame, optional): Bootstrap results for each nominal
                fit. Defaults to None.
            proj_moments_df (pd.DataFrame, optional): Contains the projected moments
                results for each nominal fit. Defaults to None.
            bootstrap_proj_moments_df (pd.DataFrame, optional): Contains the projected
                moments results for each bootstrap fit. Defaults to None.
            truth_df (pd.DataFrame, optional): Contains the ground truth values for
                the fit. Defaults to None.
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
        self.channel = channel

        # Common properties
        self._masses = self.data_df["m_center"]
        self._bin_width = (self.data_df["m_high"] - self.data_df["m_low"])[0]

        # Extract coherent sums and phase differences
        self.coherent_sums = utils.get_coherent_sums(self.fit_df)
        self.phase_difference_dict = utils.get_phase_difference_dict(self.fit_df)
        self.phase_differences = utils.get_phase_differences(self.fit_df)
