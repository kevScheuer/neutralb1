"""move here:
- phase_wrapping
- truth phase handling
- dataframe linking
- data validation
"""

import warnings
from typing import Optional

import numpy as np
import pandas as pd


class PWAPreprocessor:

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        bootstrap_df: Optional[pd.DataFrame] = None,
        truth_df: Optional[pd.DataFrame] = None,
    ) -> None:

        self.fit_df = fit_df
        self.data_df = data_df
        self.bootstrap_df = bootstrap_df if bootstrap_df is not None else None
        self.truth_df = truth_df if truth_df is not None else None

        # warn user of fits with incomplete error matrices
        if any(status != 3 for status in self.fit_df["eMatrixStatus"]):
            bad_files = self.fit_df[self.fit_df["eMatrixStatus"] != 3]["file"].to_list()
            warnings.warn(
                f"The following files contain fit results whose covariance matrix is"
                f" not full and accurate:\n{"\n".join(bad_files)}",
                UserWarning,
            )

    def validate_dataframes(self) -> None:
        """Validate that all dataframes are pandas DataFrames and not empty.

        Raises:
            TypeError: If any input is not a pandas DataFrame
            ValueError: If any dataframe is empty
        """
        if not isinstance(self.fit_df, pd.DataFrame):
            raise TypeError("fit_df must be a pandas DataFrame")
        if not isinstance(self.data_df, pd.DataFrame):
            raise TypeError("data_df must be a pandas DataFrame")
        if self.bootstrap_df is not None and not isinstance(
            self.bootstrap_df, pd.DataFrame
        ):
            raise TypeError("bootstrap_df must be a pandas DataFrame")
        if self.truth_df is not None and not isinstance(self.truth_df, pd.DataFrame):
            raise TypeError("truth_df must be a pandas DataFrame")

        if self.fit_df.empty:
            raise ValueError("fit_df cannot be empty")
        if self.data_df.empty:
            raise ValueError("data_df cannot be empty")
        if self.bootstrap_df is not None and self.bootstrap_df.empty:
            raise ValueError("bootstrap_df was passed but is empty")
        if self.truth_df is not None and self.truth_df.empty:
            raise ValueError("truth_df was passed but is empty")

    def compare_shapes(self) -> None:

        number_of_fit_rows = self.fit_df.shape[0]
        number_of_data_rows = self.data_df.shape[0]
        number_of_truth_rows = (
            self.truth_df.shape[0] if self.truth_df is not None else None
        )

        # number of bootstrap fits can vary, but each batch should live in the same
        # directory and be associated with one fit result
        number_of_bootstrap_directories = (
            self.bootstrap_df["file"].str.rsplit("/", n=1).str[0].unique()
            if self.bootstrap_df is not None
            else None
        )

        if number_of_fit_rows != number_of_data_rows:
            raise ValueError(
                "fit_df and data_df must have the same number of rows\n"
                f"fit_df: {number_of_fit_rows}"
                f", data_df: {number_of_data_rows}"
            )
        if (
            number_of_truth_rows is not None
            and number_of_fit_rows != number_of_truth_rows
        ):
            raise ValueError(
                "fit_df and truth_df must have the same number of rows\n"
                f"fit_df: {number_of_fit_rows}"
                f", truth_df: {number_of_truth_rows}"
            )
        if (
            number_of_bootstrap_directories is not None
            and number_of_fit_rows != number_of_bootstrap_directories
        ):
            raise ValueError(
                "There must be a set of bootstrap fits for each fit result.\n"
                f"fit_df: {number_of_fit_rows}"
                f", bootstrap_df: {number_of_bootstrap_directories}"
            )

    # TODO: have phase wrapper call that wraps phases of each df, and handles truth
    # in its own way. Consider having it modify dataframes in place or return new ones.
