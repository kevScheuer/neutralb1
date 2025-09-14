"""Contains the ResultManager class with methods to manage and analyze fit results."""

import warnings
from typing import Optional

import pandas as pd

import neutralb1.utils as utils
from neutralb1.analysis import preprocessing
from neutralb1.analysis.plotting.factory_plotter import FactoryPlotter


class ResultManager:

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
    ) -> None:
        """Initialize the ResultManager.

        This class is the main interface for managing and analyzing fit results. At its
        core, a fit result consists of a fit to data, collected in the `fit_df` and
        `data_df` DataFrames. Additional DataFrames can be provided to include
        randomized fits, bootstrap fits, projected moments, and ground truth values.
        The class provides preprocessing methods to clean and standardize the
        DataFrames, and plotting methods to visualize the results.

        See :py:mod:`neutralb1.batch.convert_to_csv` for more details on how the
        DataFrames are constructed.

        Args:
            fit_df (pd.DataFrame): Nominal fit results DataFrame. These are typically
                the "best" fits of many randomized ones.
            data_df (pd.DataFrame): Contains the data used for the fit.
            proj_moments_df (pd.DataFrame, optional): Contains the projected moments
                calculated from the fit results. Defaults to None.
            randomized_df (pd.DataFrame, optional): Contains all randomized fits
                associated with each nominal fit result. Defaults to None.
            randomized_proj_moments_df (pd.DataFrame, optional): Contains the projected
                moments calculated from the randomized fits. Defaults to None.
            bootstrap_df (pd.DataFrame, optional): bootstrap results for each nominal
                fit. Defaults to None.
            bootstrap_proj_moments_df (pd.DataFrame, optional): Contains the projected
                moments calculated from the bootstrap fits. Defaults to None.
            truth_df (pd.DataFrame, optional): Contains the ground truth values for the
                fit. Only applicable for Monte Carlo Input-Output Studies. Defaults to
                None.
            truth_proj_moments_df (pd.DataFrame, optional): Contains the ground truth
                projected moments for the fit, i.e. the projected moments obtained from
                the truth_df. Only applicable for Monte Carlo Input-Output Studies.
                Defaults to None.
        """

        # create local copies of the DataFrames to avoid modifying the originals
        self.fit_df = fit_df.copy()
        self.data_df = data_df.copy()
        self.proj_moments_df = (
            proj_moments_df.copy() if proj_moments_df is not None else None
        )
        self.randomized_df = randomized_df.copy() if randomized_df is not None else None
        self.randomized_proj_moments_df = (
            randomized_proj_moments_df.copy()
            if randomized_proj_moments_df is not None
            else None
        )
        self.bootstrap_df = bootstrap_df.copy() if bootstrap_df is not None else None
        self.bootstrap_proj_moments_df = (
            bootstrap_proj_moments_df.copy()
            if bootstrap_proj_moments_df is not None
            else None
        )
        self.truth_df = truth_df.copy() if truth_df is not None else None
        self.truth_proj_moments_df = (
            truth_proj_moments_df.copy() if truth_proj_moments_df is not None else None
        )

        # warn user of any bad error matrices
        if any(status != 3 for status in self.fit_df["eMatrixStatus"]):
            bad_files = self.fit_df[self.fit_df["eMatrixStatus"] != 3]["file"].to_list()
            warnings.warn(
                f"The following files contain fit results whose covariance matrix is"
                f" not full and accurate:\n{"\n".join(bad_files)}",
                UserWarning,
            )

        self.coherent_sums = utils.get_coherent_sums(self.fit_df)
        self.phase_difference_dict = utils.get_phase_difference_dict(self.fit_df)
        self.phase_differences = utils.get_phase_differences(self.fit_df)

        self.plotter_factory = None  # initialize to None until plotter is called

        return

    def preprocess(self, linker_max_depth: int = 2) -> None:
        """Preprocess the DataFrames, and modify their class copies in place.

        This includes:
            - Checking for null values in the DataFrames and warning the user
            - Linking the DataFrames together by adding a 'fit_index' column
            - Converting unbound phase columns in radians to degrees from -180 to 180
            - Standardizing types across DataFrames to save memory
            - Aligning phase difference names across DataFrames
            - Remove projected moment columns expected to be 0
                Imag(H0) = Imag(H1) = Real(H2) = 0
            - Remove real/imaginary suffixes from projected moment columns
            - Adding missing columns to the truth DataFrame

        Args:
            linker_max_depth (int, optional): Maximum depth for linking DataFrames.
                Defaults to 1, assuming that all DataFrames are a sibling or child of
                the fit DataFrame.
                See :func:`neutralb1.analysis.utils.link_dataframes` for details.

        """
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

        if self.randomized_df is not None and preprocessing.find_null_columns(
            self.randomized_df
        ):
            warnings.warn(
                "The random DataFrame contains null values. Consider checking the "
                "following columns: "
                f"{', '.join(preprocessing.find_null_columns(self.randomized_df))}",
                UserWarning,
            )

        if (
            self.randomized_proj_moments_df is not None
            and preprocessing.find_null_columns(self.randomized_proj_moments_df)
        ):
            warnings.warn(
                "The randomized projected moments DataFrame contains null values."
                " Consider checking the following columns: "
                f"{', '.join(
                    preprocessing.find_null_columns(self.randomized_proj_moments_df)
                )}",
                UserWarning,
            )

        if self.bootstrap_df is not None and preprocessing.find_null_columns(
            self.bootstrap_df
        ):
            warnings.warn(
                "The bootstrap DataFrame contains null values. Consider checking the "
                "following columns: "
                f"{', '.join(preprocessing.find_null_columns(self.bootstrap_df))}",
                UserWarning,
            )

        if self.proj_moments_df is not None and preprocessing.find_null_columns(
            self.proj_moments_df
        ):
            warnings.warn(
                "The projected moments DataFrame contains null values. Consider"
                " checking the following columns: "
                f"{', '.join(preprocessing.find_null_columns(self.proj_moments_df))}",
                UserWarning,
            )

        if (
            self.bootstrap_proj_moments_df is not None
            and preprocessing.find_null_columns(self.bootstrap_proj_moments_df)
        ):
            warnings.warn(
                "The bootstrap projected moments DataFrame contains null values."
                " Consider checking the following columns: "
                f"{', '.join(
                    preprocessing.find_null_columns(self.bootstrap_proj_moments_df)
                )}",
                UserWarning,
            )

        if self.truth_df is not None and preprocessing.find_null_columns(self.truth_df):
            warnings.warn(
                "The truth DataFrame contains null values. Consider checking the "
                "following columns: "
                f"{', '.join(preprocessing.find_null_columns(self.truth_df))}",
                UserWarning,
            )

        if self.truth_proj_moments_df is not None and preprocessing.find_null_columns(
            self.truth_proj_moments_df
        ):
            warnings.warn(
                "The truth projected moments DataFrame contains null values. Consider"
                " checking the following columns: "
                f"{', '.join(
                    preprocessing.find_null_columns(self.truth_proj_moments_df)
                )}",
                UserWarning,
            )

        # link the DataFrames together by adding a 'fit_index' column
        self.data_df = preprocessing.link_dataframes(
            self.fit_df, self.data_df, linker_max_depth
        )
        if self.randomized_df is not None:
            self.randomized_df = preprocessing.link_dataframes(
                self.fit_df, self.randomized_df, linker_max_depth
            )
        if self.randomized_proj_moments_df is not None:
            self.randomized_proj_moments_df = preprocessing.link_dataframes(
                self.fit_df, self.randomized_proj_moments_df, linker_max_depth
            )
        if self.bootstrap_df is not None:
            self.bootstrap_df = preprocessing.link_dataframes(
                self.fit_df, self.bootstrap_df, linker_max_depth
            )
        if self.proj_moments_df is not None:
            self.proj_moments_df = preprocessing.link_dataframes(
                self.fit_df, self.proj_moments_df, linker_max_depth
            )
        if self.bootstrap_proj_moments_df is not None:
            self.bootstrap_proj_moments_df = preprocessing.link_dataframes(
                self.fit_df, self.bootstrap_proj_moments_df, linker_max_depth
            )
        if self.truth_df is not None:
            self.truth_df = preprocessing.link_dataframes(
                self.fit_df, self.truth_df, linker_max_depth
            )
        if self.truth_proj_moments_df is not None:
            self.truth_proj_moments_df = preprocessing.link_dataframes(
                self.fit_df, self.truth_proj_moments_df, linker_max_depth
            )

        # convert unbound phase columns in radians to degrees from -180 to 180
        self.fit_df = preprocessing.wrap_phases(self.fit_df)
        self.data_df = preprocessing.wrap_phases(self.data_df)
        if self.randomized_df is not None:
            self.randomized_df = preprocessing.wrap_phases(self.randomized_df)
        if self.bootstrap_df is not None:
            self.bootstrap_df = preprocessing.wrap_phases(self.bootstrap_df)
        if self.truth_df is not None:
            self.truth_df = preprocessing.wrap_phases(self.truth_df)

        # Remove projected moment columns expected to be 0
        if self.proj_moments_df is not None:
            self.proj_moments_df = preprocessing.filter_projected_moments(
                self.proj_moments_df
            )
        if self.randomized_proj_moments_df is not None:
            self.randomized_proj_moments_df = preprocessing.filter_projected_moments(
                self.randomized_proj_moments_df
            )
        if self.bootstrap_proj_moments_df is not None:
            self.bootstrap_proj_moments_df = preprocessing.filter_projected_moments(
                self.bootstrap_proj_moments_df
            )
        if self.truth_proj_moments_df is not None:
            self.truth_proj_moments_df = preprocessing.filter_projected_moments(
                self.truth_proj_moments_df
            )

        # remove real/imaginary suffixes from projected moment columns
        if self.proj_moments_df is not None:
            self.proj_moments_df = preprocessing.remove_real_imag_suffixes(
                self.proj_moments_df
            )
        if self.randomized_proj_moments_df is not None:
            self.randomized_proj_moments_df = preprocessing.remove_real_imag_suffixes(
                self.randomized_proj_moments_df
            )
        if self.bootstrap_proj_moments_df is not None:
            self.bootstrap_proj_moments_df = preprocessing.remove_real_imag_suffixes(
                self.bootstrap_proj_moments_df
            )
        if self.truth_proj_moments_df is not None:
            self.truth_proj_moments_df = preprocessing.remove_real_imag_suffixes(
                self.truth_proj_moments_df
            )

        # standardize types across DataFrames to save memory
        self.fit_df = preprocessing.standardize_fit_types(self.fit_df)
        self.data_df = preprocessing.standardize_data_types(self.data_df)
        if self.randomized_df is not None:
            self.randomized_df = preprocessing.standardize_fit_types(self.randomized_df)
        if self.randomized_proj_moments_df is not None:
            self.randomized_proj_moments_df = preprocessing.standardize_moment_types(
                self.randomized_proj_moments_df
            )
        if self.bootstrap_df is not None:
            self.bootstrap_df = preprocessing.standardize_fit_types(self.bootstrap_df)
        if self.proj_moments_df is not None:
            self.proj_moments_df = preprocessing.standardize_moment_types(
                self.proj_moments_df
            )
        if self.bootstrap_proj_moments_df is not None:
            self.bootstrap_proj_moments_df = preprocessing.standardize_moment_types(
                self.bootstrap_proj_moments_df
            )
        if self.truth_df is not None:
            self.truth_df = preprocessing.standardize_fit_types(self.truth_df)
        if self.truth_proj_moments_df is not None:
            self.truth_proj_moments_df = preprocessing.standardize_moment_types(
                self.truth_proj_moments_df
            )

        # align phase difference names across DataFrames
        if self.randomized_df is not None:
            self.randomized_df = preprocessing.align_phase_difference_names(
                self.fit_df, self.randomized_df
            )
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
        if self.truth_df is not None:
            self.truth_df = preprocessing.add_missing_columns(
                self.fit_df, self.truth_df
            )

        # truth fits always include a Breit-Wigner in the amplitude, so we need to
        # reset the phases in the truth DataFrame, as AmpTools does not include
        # this modulation in its phaseDiff method
        if self.truth_df is not None:
            # make sure that the list of mass bins is in the same order as the truth_df
            # by checking the fit_index column
            mass_bins = self.data_df.set_index("fit_index").loc[
                self.truth_df["fit_index"], "m_avg"
            ]

            self.truth_df = preprocessing.restore_breit_wigner_phases(
                self.truth_df, mass_bins
            )

    @property
    def plot(self) -> FactoryPlotter:
        """Return a FactoryPlotter instance for plotting fit results."""

        if self.plotter_factory is None:
            self.plotter_factory = FactoryPlotter(
                fit_df=self.fit_df,
                data_df=self.data_df,
                proj_moments_df=self.proj_moments_df,
                randomized_df=self.randomized_df,
                randomized_proj_moments_df=self.randomized_proj_moments_df,
                bootstrap_df=self.bootstrap_df,
                bootstrap_proj_moments_df=self.bootstrap_proj_moments_df,
                truth_df=self.truth_df,
            )
        return self.plotter_factory
