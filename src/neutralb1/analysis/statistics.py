from typing import List, Optional

import pandas as pd

import neutralb1.utils as utils

# TODO: complete this function. Uses old variables from old Plotter class
# def normality_test(
#     self,
#     columns: Optional[List[str]] = None,
#     method: Optional[Literal["None", "normal", "log_normal"]] = None,
# ) -> None:
#     """
#     Shapiro-wilkes for parameters and amplitudes. Log-transform the amps near zero.
#     Use other method for circular (phase) data.

#     """
#     if self.bootstrap_df is None:
#         raise ValueError("Bootstrap df was not defined on instantiation")

#     # PLAN: loop over the fit indices of the bootstrap df and calculate the test.
#     # Those that don't pass, or whose bootstrap mean is too far from the fit
#     # result's are plotted in big pdf. There should be indicators for the fit mean
#     # value on the QQ plot, and whether it failed due to p value or mean

#     fit_indices = self.bootstrap_df["fit_index"].unique()

#     # add the amplitudes and phase differences to the columns to be tested
#     columns = columns if columns is not None else []
#     columns.extend(self.coherent_sums["eJPmL"])
#     columns.extend(set(self.phase_differences.values()))

#     failed_dict = {}  # to store {fit_index: {column : [mean, p-value]}}
#     for fit_index in fit_indices:
#         for col in columns:
#             if col in set(self.phase_differences.values()):
#                 # do circular data test. Take absolute value because we have no way
#                 # to distinguish positive vs negative values
#                 pass

#             elif any(col in sublist for sublist in self.coherent_sums.values()):
#                 # do normality test on amplitudes
#                 # log-transform the amps near zero
#                 # type: ignore or type annotation for pylance
#                 amp_values = self.bootstrap_df.loc[
#                     self.bootstrap_df["fit_index"] == fit_index
#                 ][col]

#                 # log-transform if col is near zero
#                 # use rel distance to minimum test to determine if near zero
#                 # TODO: finish this

#                 # perform normality test
#                 stat, p_value = scipy.stats.shapiro(amp_values)
#                 if p_value < 0.05:
#                     failed_dict.setdefault(fit_index, {})[col] = [
#                         amp_values.mean(),
#                         p_value,
#                     ]

#             else:
#                 # do normality test as usual
#                 values = self.bootstrap_df.loc[
#                     self.bootstrap_df["fit_index"] == fit_index
#                 ][col]

#                 # perform normality test
#                 stat, p_value = scipy.stats.shapiro(values)
#                 if p_value < 0.05:
#                     failed_dict.setdefault(fit_index, {})[col] = [
#                         values.mean(),
#                         p_value,
#                     ]

#     # PLAN: here we'll loop over the failed ones and plot probdists

#     pass


def report_correlations(
    df: pd.DataFrame,
    fit_indices: Optional[List[int]] = None,
    correlation_threshold: float = 0.8,
    report_average: bool = True,
    drop_columns: Optional[List[str]] = None,
) -> None:
    """Report correlations between variables in a DataFrame.

    Assumes that the DataFrame contains bootstrap samples for different fit indices.

    Args:
        df (pd.DataFrame): The input DataFrame.
        fit_indices (Optional[List[int]], optional): The fit indices to consider.
            Defaults to None.
        correlation_threshold (float, optional): The correlation threshold for
            reporting. Defaults to 0.8.
        report_average (bool, optional): Whether to report average correlations.
            Defaults to True.
        drop_columns (Optional[List[str]], optional): Columns to drop from the
            analysis. Defaults to None.

    Raises:
        ValueError: If no fit indices are found in the DataFrame.

    Returns:
        None: Prints out pairs of variables with high correlations.
    """

    if fit_indices is None:
        fit_indices = df["fit_index"].unique().tolist()

    if fit_indices == []:
        raise ValueError("No fit indices provided or found in DataFrame")

    # drop error columns and non-numeric fit result info columns
    cols_to_drop = (
        [c for c in df.columns if c.endswith("_err")]
        + [
            "likelihood",
            "detected_events",
            "generated_events",
            "eMatrixStatus",
            "lastMinuitCommandStatus",
        ]
        + (drop_columns if drop_columns is not None else [])
    )

    # get coherent sums for later filtering
    coh_sums = utils.get_coherent_sums(df)

    stacked_df = (
        df.select_dtypes(include="number")  # only numeric columns
        .assign(  # resolve sign ambiguity of phase differences
            **{
                col: df[col].abs()
                for col in df.columns
                if col in utils.get_phase_differences(df)
            }
        )
        .drop(columns=cols_to_drop)
        .groupby("fit_index")  # group bootstrap samples by fit index
        .corr()  # calculate correlation matrix for each group
        .abs()  # take absolute value of correlations for later averaging
        .loc[fit_indices]  # select only the requested fit indices
        .stack()  # reshape to long format
        .reset_index()
        .rename(  # rename columns for clarity
            columns={"level_1": "var1", "level_2": "var2", 0: "correlation"}
        )
    )

    # filter dataframe
    stacked_df = stacked_df.loc[
        (stacked_df["var1"] != stacked_df["var2"])  # exclude self-correlations
        & (stacked_df["var1"] < stacked_df["var2"])  # keep only one ordering
        & (
            stacked_df.apply(
                # remove corrs. between coherent sum and non-coherent sum vars
                lambda row: _same_coherent_sum_key(row["var1"], row["var2"], coh_sums),
                axis=1,
            )
        )
        & (
            stacked_df.apply(
                # remove corrs. between production coeffs and non-prod coeffs
                lambda row: _filter_production_coefficients(row["var1"], row["var2"]),
                axis=1,
            )
        )
    ]

    if report_average:  # average the correlations over the fit indices
        avg_corrs = stacked_df.groupby(["var1", "var2"])[
            "correlation"
        ].mean()  # select just the correlations  # average over fit indices
        for var1, var2 in avg_corrs.index:
            avg_corr = avg_corrs.loc[(var1, var2)]
            if avg_corr >= correlation_threshold:  # only report high correlations
                print(f"{var1} \t{var2}\t:{avg_corr:.2f}")
    else:  # report correlations for each fit index separately
        for _, row in stacked_df.iterrows():
            if row["correlation"] >= correlation_threshold:  # only report high corrs
                print(
                    f"{row['fit_index']}"
                    f"\t{row['var1']} \t{row['var2']}"
                    f"\t:{row['correlation']:.2f}"
                )

    return None


def _same_coherent_sum_key(var1: str, var2: str, coh_sums: dict) -> bool:
    """Check if var1 and var2 belong to the same coherent sum group.

    Args:
        var1 (str): First variable name.
        var2 (str): Second variable name.
        coh_sums (dict): Dictionary of coherent sums.

    Returns:
        True if both variables belong to the same coherent sum group or are not coherent
            sums, False otherwise.
    """

    all_values = sum(coh_sums.values(), [])  # get combined list of coh sums

    # make sure var1 and var2 are coherent sums first, since the loop only checks
    # within coherent sum groups
    if var1 not in all_values or var2 not in all_values:
        return True

    for key, vals in coh_sums.items():
        if var1 in vals and var2 in vals:
            return True
    return False


def _filter_production_coefficients(var1: str, var2: str) -> bool:
    """If one is a production coefficient, the other must be too.

    Args:
        var1 (str): First variable name.
        var2 (str): Second variable name.

    Returns:
        False if one is a production coefficient and the other is not, True otherwise.
    """
    if (var1.endswith("_re") or var1.endswith("_im")) and not (
        var2.endswith("_re") or var2.endswith("_im")
    ):
        return False

    if (var2.endswith("_re") or var2.endswith("_im")) and not (
        var1.endswith("_re") or var1.endswith("_im")
    ):
        return False

    return True
