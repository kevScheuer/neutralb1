from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from matplotlib.backends.backend_pdf import PdfPages

import neutralb1.utils as utils


def normality_test(
    fit_df: pd.DataFrame,
    bootstrap_df: pd.DataFrame,
    columns: List[str],
    alpha: float = 0.05,
    output: str = "./normality_test_results.pdf",
    ignore_small_amps: bool = True,
    transform_small_amps: bool = True,
    small_amp_threshold: float = 0.05,
    skew_threshold: float = 0.5,
    n_mc_samples: int = 1000,
) -> None:
    """Test normality of bootstrap distributions for specified columns.

    This function handles 3 cases:
    1. Circular data (phase differences): uses the von Mises distribution for testing.
    2. Amplitude data: if the distribution is highly skewed and the fit fraction
        is small, this means the amplitude is affected by the boundary at zero. Most of
        the time we do not care about estimating the error on these amplitudes, and
        consider them to be negligible (ignore_small_amps=True). However, if we do want
        to test their normality, the standard Shapiro-Wilk test will often fail due to
        this boundary effect. In this case, we log-transform the data before testing for
        normality. This is configured via the `transform_small_amps`,
        `small_amp_threshold`, and `skew_threshold` parameters. Uses the Shapiro-Wilk
        test for normality.
    3. Regular data: uses the Shapiro-Wilk test for normality.

    Args:
        fit_df (pd.DataFrame): DataFrame containing nominal fit results.
        bootstrap_df (pd.DataFrame): DataFrame containing bootstrap fit results.
        columns (List[str]): List of column names to test for normality.
        alpha (float, optional): Significance level for normality tests.
            Defaults to 0.05.
        output (str, optional): Path to output PDF file for plots of failed tests.
            Defaults to "./normality_test_results.pdf".
        ignore_small_amps (bool, optional): Whether to ignore amplitude distributions
            with small fit fractions when testing for normality. Defaults to True.
        transform_small_amps (bool, optional): Whether to log-transform "small"
            amplitude distributions before testing for normality. Defaults to True.
        small_amp_threshold (float, optional): Threshold for fit fraction below which
            amplitudes are considered small. Defaults to 0.05.
        skew_threshold (float, optional): Skewness threshold above which amplitude
            distributions that are also "small" will be log-transformed before testing.
            Defaults to 1.0.
        n_mc_samples (int, optional): Number of Monte Carlo samples for circular
            normality test. Increasing will increase accuracy but also computation time.
            Defaults to 1000.
    Raises:
        KeyError: If any specified column is not found in either DataFrame.
        ValueError: If phase difference values are not within [-180, 180] degrees.
        TypeError: If amplitude fit value or number of events are not numeric types.
    Returns:
        None: Generates a PDF file with probability plots for failed normality tests.
    """

    # check that columns exist in both dataframes
    for col in columns:
        if col not in fit_df.columns:
            raise KeyError(f"Column {col} not found in fit_df")
        if col not in bootstrap_df.columns:
            raise KeyError(f"Column {col} not found in bootstrap_df")

    # if the dataframe has amplitudes or phase differences, we'll want to adjust our
    # tests accordingly
    phases = utils.get_phase_differences(fit_df)
    coherent_sums = utils.get_coherent_sums(fit_df)

    failed_dict = {}  # to store {fit_index: {column : p_value}}
    for fit_index in bootstrap_df["fit_index"].unique():
        for col in columns:
            values = bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][col]

            if col in phases:
                # do circular data test
                if values.max() > 180.0 or values.min() < -180.0:
                    raise ValueError("phase difference values must be within [-pi, pi]")

                values_rad = np.deg2rad(values.to_numpy())
                loc = scipy.stats.circmean(values_rad, high=np.pi, low=0)
                res = scipy.stats.goodness_of_fit(
                    scipy.stats.vonmises,
                    values_rad,
                    fit_params={"loc": loc},
                    statistic="cvm",
                    n_mc_samples=n_mc_samples,
                )
                p_value = res.pvalue

            elif any(col in sublist for sublist in coherent_sums.values()):
                # if the amplitude is highly skewed and near zero, its distribution
                # will be affected by the boundary at zero. In that case, we can try
                # to log-transform the data to better test for normality.
                skew = scipy.stats.skew(values)

                fit_fraction = (
                    values.mean()
                    / bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][
                        "detected_events"
                    ].mean()
                )

                if (
                    skew > skew_threshold
                    and fit_fraction < small_amp_threshold
                    and transform_small_amps
                ):
                    if ignore_small_amps:
                        # skip normality test for small amplitudes
                        continue
                    # log-transform
                    log_values = np.log1p(values)

                    # perform normality test on log-transformed values
                    stat, p_value = scipy.stats.shapiro(log_values)
                else:
                    # do normality test as usual
                    stat, p_value = scipy.stats.shapiro(values)

            else:
                # do normality test as usual
                stat, p_value = scipy.stats.shapiro(values)

            if p_value < alpha:
                if fit_index not in failed_dict:
                    failed_dict[fit_index] = {}
                failed_dict[fit_index][col] = p_value

    # Make probability plots for the failed tests against the appropriate distribution
    if not failed_dict:
        print("All normality tests passed.")
        return

    with PdfPages(output) as pdf:
        for fit_index, col_dict in failed_dict.items():
            # setup the figure as a square grid of subplots
            num_plots = len(col_dict)
            ncols = int(np.ceil(np.sqrt(num_plots)))
            nrows = int(np.ceil(num_plots / ncols))
            fig, axes = plt.subplots(
                nrows=nrows,
                ncols=ncols,
                figsize=(5 * ncols, 5 * nrows),
                layout="constrained",
            )
            axes = axes.flatten() if num_plots > 1 else [axes]

            # loop through each failed column and make the probability plot
            for ax, (col, p_value) in zip(axes, col_dict.items()):
                values = bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][col]

                if col in phases:
                    # circular data QQ plot against von Mises distribution
                    values_rad = np.deg2rad(values.to_numpy())
                    kappa, loc, scale = scipy.stats.vonmises.fit(values_rad)

                    scipy.stats.probplot(
                        values_rad,
                        dist=scipy.stats.vonmises,  # type: ignore
                        sparams=(kappa, loc, scale),
                        plot=ax,
                    )
                    ax.set_title(
                        f"{utils.convert_amp_name(col)} (circular): p={p_value:.3e}"
                    )

                    # add inset with histogram and fitted von Mises distribution
                    inset_ax = ax.inset_axes([0.6, 0.05, 0.35, 0.35])
                    inset_ax.hist(values_rad, bins=30, density=True, alpha=0.7)
                    x_range = np.linspace(
                        values_rad.min() - 0.1 * values_rad.min(),
                        values_rad.max() + 0.1 * values_rad.max(),
                        200,
                    )
                    inset_ax.plot(
                        x_range,
                        scipy.stats.vonmises.pdf(x_range, kappa, loc, scale),
                        "r-",
                        lw=2,
                    )
                    inset_ax.tick_params(labelsize=8)

                elif any(col in sublist for sublist in coherent_sums.values()):
                    # amplitude QQ plot against normal distribution
                    skew = scipy.stats.skew(values)

                    fit_fraction = (
                        values.mean()
                        / bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][
                            "detected_events"
                        ].mean()
                    )

                    if (
                        skew > skew_threshold
                        and fit_fraction < small_amp_threshold
                        and transform_small_amps
                    ):
                        # log-transform
                        log_values = np.log1p(values)
                        scipy.stats.probplot(log_values, dist="norm", plot=ax)
                        ax.set_title(
                            f"{utils.convert_amp_name(col)}"
                            f" (log-transformed): p={p_value:.3e}"
                        )

                        # add inset with histogram and fitted lognorm distribution
                        inset_ax = ax.inset_axes([0.6, 0.05, 0.35, 0.35])
                        inset_ax.hist(values, bins=30, density=True, alpha=0.7)
                        x_range = np.linspace(
                            values.min() - 0.1 * values.min(),
                            values.max() + 0.1 * values.max(),
                            200,
                        )
                        # fit lognorm to original (non-transformed) data
                        shape, loc, scale = scipy.stats.lognorm.fit(values)
                        inset_ax.plot(
                            x_range,
                            scipy.stats.lognorm.pdf(x_range, shape, loc, scale),
                            "r-",
                            lw=2,
                        )
                        inset_ax.tick_params(labelsize=8)
                    else:
                        scipy.stats.probplot(values, dist="norm", plot=ax)
                        ax.set_title(f"{utils.convert_amp_name(col)}: p={p_value:.3e}")

                        # add inset with histogram and fitted normal distribution
                        inset_ax = ax.inset_axes([0.6, 0.05, 0.35, 0.35])
                        inset_ax.hist(values, bins=30, density=True, alpha=0.7)
                        x_range = np.linspace(
                            values.min() - 0.1 * values.min(),
                            values.max() + 0.1 * values.max(),
                            200,
                        )
                        mu, sigma = values.mean(), values.std()
                        inset_ax.plot(
                            x_range,
                            scipy.stats.norm.pdf(x_range, mu, sigma),
                            "r-",
                            lw=2,
                        )
                        inset_ax.tick_params(labelsize=8)
                else:
                    scipy.stats.probplot(values, dist="norm", plot=ax)
                    ax.set_title(f"{utils.convert_amp_name(col)}: p={p_value:.3e}")

                    # add inset with histogram and fitted normal distribution
                    inset_ax = ax.inset_axes([0.6, 0.05, 0.35, 0.35])
                    inset_ax.hist(values, bins=30, density=True, alpha=0.7)
                    x_range = np.linspace(
                        values.min() - 0.1 * values.min(),
                        values.max() + 0.1 * values.max(),
                        200,
                    )
                    mu, sigma = values.mean(), values.std()
                    inset_ax.plot(
                        x_range, scipy.stats.norm.pdf(x_range, mu, sigma), "r-", lw=2
                    )
                    inset_ax.tick_params(labelsize=8)

                # finish loop over columns

            # hide any unused subplots
            for ax in axes[num_plots:]:
                ax.set_visible(False)

            fig.suptitle(
                f"Normality Test Failures for Fit Index {fit_index}", fontsize=16
            )

            pdf.savefig(fig)
            plt.close(fig)
        print("Normality test results saved to:", output)

    return


def mean_test(
    fit_df: pd.DataFrame,
    bootstrap_df: pd.DataFrame,
    columns: List[str],
) -> None:
    """
    Todo:
        this function will test whether the mean of the bootstrap samples is
        consistent with the fit result value for each column. Those that fail will be
        plotted in a big pdf with subplot for each failed column, showing its bootstrap
        distribution and the fit result value.
    """

    pass


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
