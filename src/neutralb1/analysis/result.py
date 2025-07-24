"""Contains the ResultManager class with methods to manage and analyze fit results."""

import warnings
from typing import Optional

import pandas as pd

from neutralb1.analysis import plotting, preprocessing


class ResultManager:

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        bootstrap_df: Optional[pd.DataFrame] = None,
        truth_df: Optional[pd.DataFrame] = None,
    ) -> None:
        """Initialize the ResultManager and run preprocessing steps.

        See [link to converter script] for more details on how the DataFrames are
        constructed.

        Args:
            fit_df (pd.DataFrame): Nominal fit results DataFrame. These are typically
                the "best" fits of many randomized ones.
            data_df (pd.DataFrame): Contains the data used for the fit. bootstrap_df
            (pd.DataFrame, optional): bootstrap results for each nominal
                fit. Defaults to None.
            truth_df (pd.DataFrame, optional): Contains the ground truth values for the
                fit. Only applicable for Monte Carlo Input-Output Studies. Defaults to
                None.
        """

        # create local copies of the DataFrames to avoid modifying the originals
        self.fit_df = fit_df.copy()
        self.data_df = data_df.copy()
        self.bootstrap_df = (
            bootstrap_df.copy() if bootstrap_df is not None else pd.DataFrame()
        )
        self.truth_df = truth_df.copy() if truth_df is not None else pd.DataFrame()

        # warn user of any bad error matrices
        if any(status != 3 for status in self.fit_df["eMatrixStatus"]):
            bad_files = self.fit_df[self.fit_df["eMatrixStatus"] != 3]["file"].to_list()
            warnings.warn(
                f"The following files contain fit results whose covariance matrix is"
                f" not full and accurate:\n{"\n".join(bad_files)}",
                UserWarning,
            )

        # --PREPROCESSING--

        # warn of missing values in the DataFrames
        if preprocessing.find_null_columns(self.fit_df):
            warnings.warn(
                "The fit DataFrame contains null values. Consider checking the "
                "following columns: "
                f"{', '.join(preprocessing.find_null_columns(self.fit_df))}",
                UserWarning,
            )

        if preprocessing.find_null_columns(self.data_df):
            warnings.warn(
                "The data DataFrame contains null values. Consider checking the "
                "following columns: "
                f"{', '.join(preprocessing.find_null_columns(self.data_df))}",
                UserWarning,
            )

        if not self.bootstrap_df.empty and preprocessing.find_null_columns(
            self.bootstrap_df
        ):
            warnings.warn(
                "The bootstrap DataFrame contains null values. Consider checking the "
                "following columns: "
                f"{', '.join(preprocessing.find_null_columns(self.bootstrap_df))}",
                UserWarning,
            )

        if not self.truth_df.empty and preprocessing.find_null_columns(self.truth_df):
            warnings.warn(
                "The truth DataFrame contains null values. Consider checking the "
                "following columns: "
                f"{', '.join(preprocessing.find_null_columns(self.truth_df))}",
                UserWarning,
            )

        # link the DataFrames together by adding a 'fit_index' column
        self.data_df = preprocessing.link_dataframes(self.fit_df, self.data_df)
        if not self.bootstrap_df.empty:
            self.bootstrap_df = preprocessing.link_dataframes(
                self.fit_df, self.bootstrap_df
            )
        if not self.truth_df.empty:
            self.truth_df = preprocessing.link_dataframes(self.fit_df, self.truth_df)

        # convert unbound phase columns in radians to degrees from -180 to 180
        self.fit_df = preprocessing.wrap_phases(self.fit_df)
        self.data_df = preprocessing.wrap_phases(self.data_df)
        if not self.bootstrap_df.empty:
            self.bootstrap_df = preprocessing.wrap_phases(self.bootstrap_df)
        if not self.truth_df.empty:
            self.truth_df = preprocessing.wrap_phases(self.truth_df)

        # standardize types across DataFrames to save memory
        self.fit_df = preprocessing.standardize_types(self.fit_df)
        self.data_df = preprocessing.standardize_types(self.data_df)
        if not self.bootstrap_df.empty:
            self.bootstrap_df = preprocessing.standardize_types(self.bootstrap_df)
        if not self.truth_df.empty:
            self.truth_df = preprocessing.standardize_types(self.truth_df)

        # align phase difference names across DataFrames
        if self.bootstrap_df is not None:
            self.bootstrap_df = preprocessing.align_phase_difference_names(
                self.fit_df, self.bootstrap_df
            )
        if self.truth_df is not None:
            self.truth_df = preprocessing.align_phase_difference_names(
                self.fit_df, self.truth_df
            )

        # add missing columns to the truth DataFrame. Useful for waveset comparisons
        # where the fit DataFrame has a different set of columns than the truth
        if not self.truth_df.empty:
            self.truth_df = preprocessing.add_missing_columns(
                self.fit_df, self.truth_df
            )

        # --END PREPROCESSING--

        # Create plotter factory instance
        self.plotter_factory = plotting.FactoryPlotter(
            fit_df=self.fit_df,
            data_df=self.data_df,
            bootstrap_df=self.bootstrap_df,
            truth_df=self.truth_df,
        )

    @property
    def plot(self) -> plotting.FactoryPlotter:
        """Return a FactoryPlotter instance for plotting fit results."""
        return self.plotter_factory
