import pickle
import warnings
from typing import Optional

import pandas as pd

import neutralb1.utils as utils
from neutralb1.analysis import preprocessing
from neutralb1.analysis.plotting.factory_plotter import FactoryPlotter


class ResultManager:
    """The main interface for managing and analyzing fit results.

    At its core, a fit result consists of a fit to data, collected in the `fit_df` and
    `data_df` DataFrames. Additional DataFrames can be provided to include randomized
    fits, bootstrap fits, projected moments, and truth values. The class provides
    preprocessing methods to clean and standardize the DataFrames, and plotting methods
    to visualize the results.

    See :py:mod:`neutralb1.batch.convert_to_csv` for more details on how the DataFrames
    are constructed.

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
            fit. Only applicable for Monte Carlo Input-Output Studies. Defaults to None.
        truth_proj_moments_df (pd.DataFrame, optional): Contains the ground truth
            projected moments for the fit, i.e. the projected moments obtained from the
            truth_df. Only applicable for Monte Carlo Input-Output Studies. Defaults to
            None.
        is_preprocessed (bool, optional): Whether the DataFrames have already been
            preprocessed. It's unlikely that the user would set this to True. Main
            purpose is tracking if preprocessing has been done when loading from a
            pickle file. Defaults to False.

    Todo:
        - Add methods for calculating correlations / covariances from the fit results
            and storing them internally. Allows dsratio correlation plot to work.
    """

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
        is_preprocessed: bool = False,
    ) -> None:

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

        # store commonly used derived quantities
        self.coherent_sums = utils.get_coherent_sums(self.fit_df)
        self.phase_difference_dict = utils.get_phase_difference_dict(self.fit_df)
        self.phase_differences = utils.get_phase_differences(self.fit_df)
        if self.proj_moments_df is not None:
            self.moments = utils.get_moments(self.proj_moments_df)

        self.is_preprocessed = is_preprocessed
        self.is_acceptance_corrected = self._is_fit_acceptance_corrected()

        # private attributes
        self._plotter_factory = None  # initialize to None until plotter is called

        return

    def preprocess(self, linker_max_depth: int = 1) -> None:
        """Preprocess the DataFrames, and modify their class copies in place.

        This includes:
            - Checking for null values in the DataFrames and warning the user
            - Linking the DataFrames together by adding a 'fit_index' column
            - Converting unbound phase columns in radians to degrees from -180 to 180
            - Standardizing types across DataFrames to save memory
            - Aligning phase difference names across DataFrames
            - Remove projected moment columns expected to be 0
                :math: `\\Im(H0) = \\Im(H1) = \\Re(H2) = 0`
            - Remove real/imaginary suffixes from projected moment columns
            - Adding missing columns to the truth DataFrame

        Args:
            linker_max_depth (int, optional): Maximum depth for linking DataFrames.
                Defaults to 1, assuming that all DataFrames are a sibling or child of
                the fit DataFrame.
                See :func:`neutralb1.analysis.utils.link_dataframes` for details.

        """
        if self.is_preprocessed:
            warnings.warn(
                "The DataFrames have already been preprocessed. Skipping this step.",
                UserWarning,
            )
            return

        print("Checking for null values in DataFrames...")
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

        # link dataframes by their fit_index columns
        print("Linking DataFrames...")
        self.data_df = preprocessing.link_dataframes(
            self.fit_df, self.data_df, 1  # expected to be fit sibling
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
                self.fit_df, self.proj_moments_df, 1  # expected to be fit sibling
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
        print("Wrapping phase columns...")
        self.fit_df = preprocessing.wrap_phases(self.fit_df)
        self.data_df = preprocessing.wrap_phases(self.data_df)
        if self.randomized_df is not None:
            self.randomized_df = preprocessing.wrap_phases(self.randomized_df)
        if self.bootstrap_df is not None:
            self.bootstrap_df = preprocessing.wrap_phases(self.bootstrap_df)

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

            # wrap the phases after restoring the Breit-Wigner modulation
            self.truth_df = preprocessing.wrap_phases(self.truth_df)

        # Remove projected moment columns expected to be 0
        print("Filtering projected moment columns...")
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

        print("Removing real/imaginary suffixes from projected moment columns...")
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

        print("Standardizing DataFrame types...")
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

        print("Aligning phase difference names across DataFrames...")
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
            # BUG: some sort of highly fragmented output and failure here
            print("Adding missing columns to truth DataFrame...")
            self.truth_df = preprocessing.add_missing_columns(
                self.fit_df, self.truth_df
            )

        self.is_preprocessed = True  # set flag to True after preprocessing is done
        return

    def summary(self) -> None:
        """Print a summary of the fit results contained in the DataFrames.

        Args:
            None
        Returns:
            None, just prints to console
        """

        print("Fit DataFrame Summary:")
        print(self.fit_df.info())
        print("\nData DataFrame Summary:")
        print(self.data_df.info())
        if self.proj_moments_df is not None:
            print("\nProjected Moments DataFrame Summary:")
            print(self.proj_moments_df.info())
        if self.randomized_df is not None:
            print("\nRandomized Fit DataFrame Summary:")
            print(self.randomized_df.info())
        if self.randomized_proj_moments_df is not None:
            print("\nRandomized Projected Moments DataFrame Summary:")
            print(self.randomized_proj_moments_df.info())
        if self.bootstrap_df is not None:
            print("\nBootstrap Fit DataFrame Summary:")
            print(self.bootstrap_df.info())
        if self.bootstrap_proj_moments_df is not None:
            print("\nBootstrap Projected Moments DataFrame Summary:")
            print(self.bootstrap_proj_moments_df.info())
        if self.truth_df is not None:
            print("\nTruth Fit DataFrame Summary:")
            print(self.truth_df.info())
        if self.truth_proj_moments_df is not None:
            print("\nTruth Projected Moments DataFrame Summary:")
            print(self.truth_proj_moments_df.info())

        return

    def save(self, filepath: str, preprocess=True) -> None:
        """Save all preprocessed DataFrames to a pickle file.

        This method serializes all internal DataFrames to a single pickle file,
        allowing for easy reloading into a new ResultManager instance later.

        Args:
            filepath (str): Path to the output pickle file.
            preprocess (bool, optional): Whether to preprocess the DataFrames before
                saving. Defaults to True.

        Raises:
            IOError: If the file cannot be written.

        Returns:
            None, just saves to file.

        Example:
            >>> result_manager.save("fit_results.pkl")
            >>> # Later, to load the results:
            >>> with open("fit_results.pkl", "rb") as f:
            >>>     data = pickle.load(f)
            >>> new_result_manager = ResultManager(**data)
        """
        if preprocess:
            self.preprocess()

        data = {
            "fit_df": self.fit_df,
            "data_df": self.data_df,
            "proj_moments_df": self.proj_moments_df,
            "randomized_df": self.randomized_df,
            "randomized_proj_moments_df": self.randomized_proj_moments_df,
            "bootstrap_df": self.bootstrap_df,
            "bootstrap_proj_moments_df": self.bootstrap_proj_moments_df,
            "truth_df": self.truth_df,
            "truth_proj_moments_df": self.truth_proj_moments_df,
            "is_preprocessed": self.is_preprocessed,
        }
        try:
            with open(filepath, "wb") as f:
                pickle.dump(data, f)
                print(f"DataFrames successfully saved to {filepath}")
        except IOError as e:
            raise IOError(f"Could not write to {filepath}: {e}")

        return

    @property
    def plot(self) -> FactoryPlotter:
        """Return a FactoryPlotter instance for plotting fit results.

        Example:
            >>> import matplotlib.pyplot as plt
            >>> from neutralb1.analysis import ResultManager
            >>> # Assuming fit_df and data_df are already defined
            >>> result_manager = ResultManager(fit_df, data_df)
            >>> result_manager.preprocess()  # Preprocess the DataFrames
            >>> plotter = result_manager.plot
            >>> plotter.intensity.waves(sharey=True)
            >>> plt.show()
        """

        if self._plotter_factory is None:
            self._plotter_factory = FactoryPlotter(
                fit_df=self.fit_df,
                data_df=self.data_df,
                proj_moments_df=self.proj_moments_df,
                randomized_df=self.randomized_df,
                randomized_proj_moments_df=self.randomized_proj_moments_df,
                bootstrap_df=self.bootstrap_df,
                bootstrap_proj_moments_df=self.bootstrap_proj_moments_df,
                truth_df=self.truth_df,
                truth_proj_moments_df=self.truth_proj_moments_df,
                is_acceptance_corrected=self.is_acceptance_corrected,
            )
        return self._plotter_factory

    def get_significant_amplitudes(
        self, threshold: float = 0.05, fit_indices: Optional[list[int]] = None
    ) -> set[str]:
        """Determine single amplitudes whose average fit fraction is above a threshold

        Single amplitudes refers to an 'eJPmL' amplitude, that is not coherently summed.

        Args:
            threshold (float, optional): amplitudes with an average fit fraction above
                this value are considered significant. Defaults to 0.05.
            fit_indices (list[int], optional): Indices of the fit DataFrame to consider
                for averaging. Defaults to None (all indices used)

        Returns:
            set: A set of significant single amplitudes.
        """

        if fit_indices is None:
            fit_indices = list(self.fit_df.index)

        significant_amplitudes = set()

        for amp in self.coherent_sums["eJPmL"]:
            num_events = (
                self.fit_df["generated_events"]
                if self.is_acceptance_corrected
                else self.fit_df["detected_events"]
            )
            fit_fractions = (self.fit_df[amp] / num_events).loc[fit_indices]
            is_significant = fit_fractions.mean() > threshold
            if is_significant:
                significant_amplitudes.add(amp)

        return significant_amplitudes

    def get_significant_phases(
        self, threshold: float = 0.05, fit_indices: Optional[list[int]] = None
    ) -> set[str]:
        """Determine phase differences that are only between significant amplitudes

        Only returns phase differences where both amplitudes that makeup the phase
        difference are considered significant. See 'get_significant_amplitudes' for more
        details.

        Args:
            threshold (float, optional): amplitudes with an average fit fraction above
                this value are considered significant. Defaults to 0.05.
            fit_indices (list[int], optional): Indices of the fit DataFrame to consider
                for averaging. Defaults to None (all indices used).
        Returns:
            set: A set of significant phase differences.
        """

        if fit_indices is None:
            fit_indices = list(self.fit_df.index)

        significant_phases = set()
        significant_amplitudes = self.get_significant_amplitudes(threshold, fit_indices)

        for phase_dif in self.phase_differences:
            phase1, phase2 = phase_dif.split("_")
            if phase1 in significant_amplitudes and phase2 in significant_amplitudes:
                significant_phases.add(phase_dif)

        return significant_phases

    def get_significant_moments(
        self, threshold: float = 0.05, fit_indices: Optional[list[int]] = None
    ) -> set[str]:
        """Determine moments that have an average absolute fit fraction above threshold

        Args:
            threshold (float, optional): moments with an average absolute
                fit fraction (|moment / H0_0000|) above this value are considered
                significant. Defaults to 0.05.
            fit_indices (list[int], optional): Indices of the fit DataFrame to consider
                for averaging. Defaults to None (all indices used).

        Raises:
            ValueError: If the projected moments DataFrame is not available.

        Returns:
            set[str]: A set of significant moments.
        """

        if self.proj_moments_df is None:
            raise ValueError("Projected moments DataFrame is not available.")

        if fit_indices is None:
            fit_indices = list(self.fit_df.index)

        significant_moments = set()
        for moment in self.moments:
            h0 = self.proj_moments_df["H0_0000"]
            fit_fractions = (self.proj_moments_df[moment] / h0).loc[fit_indices]
            is_significant = fit_fractions.abs().mean() > threshold
            if is_significant:
                significant_moments.add(moment)

        return significant_moments

    def _is_fit_acceptance_corrected(self) -> bool:
        """Check if the fit is using acceptance corrected amplitudes.

        The sum of the positive and negative reflectivity coherent sums will always be
        the number of events, and if this is greater than the detected events, then the
        fit must be acceptance corrected. We only check the first row, as this should be
        true for all rows. The detected_events_err is included to account for small
        numerical differences.

        Returns:
            bool: True if the fit is acceptance corrected, False otherwise.
        """

        if self.fit_df[self.coherent_sums["e"]].iloc[0].sum() > (
            self.fit_df["detected_events"].iloc[0]
            + self.fit_df["detected_events_err"].iloc[0]
        ):
            return True
        return False
