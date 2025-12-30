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

    Todo:
        - handle acc_corrected (detected_events -> generated_events)
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
        num_events = bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][
            "detected_events"
        ].mean()

        for col in columns:
            values = bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][col]

            if col in phases:
                # do circular data test
                if values.max() > 180.0 or values.min() < -180.0:
                    raise ValueError("phase difference values must be within [-pi, pi]")

                # skip normality test for phase differences with small amplitudes if
                # configured to do so
                amp1, amp2 = col.split("_")
                amp1_ff = (
                    bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][
                        amp1
                    ].mean()
                    / num_events
                )
                amp2_ff = (
                    bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][
                        amp2
                    ].mean()
                    / num_events
                )
                if ignore_small_amps and (
                    amp1_ff < small_amp_threshold or amp2_ff < small_amp_threshold
                ):
                    continue

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

                fit_fraction = values.mean() / num_events

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


def bias_test(
    fit_df: pd.DataFrame,
    bootstrap_df: pd.DataFrame,
    columns: List[str],
    output: str = "./bias_test_results.pdf",
    threshold: float = 0.05,
    small_amp_threshold: float = 0.05,
) -> None:
    """Test whether the mean of bootstrap samples is consistent with fit results.

    For each column, calculates the bias as (bootstrap mean - fit value). If the
    absolute bias exceeds the threshold, the column is flagged and plotted. Amplitudes
    are calculated and plotted as fit fractions to provide a sense of scale.

    Args:
        fit_df (pd.DataFrame): DataFrame containing nominal fit results.
        bootstrap_df (pd.DataFrame): DataFrame containing bootstrap fit results.
        columns (List[str]): List of column names to test for bias.
        output (str, optional): Path to output PDF file for plots of failed tests.
            Defaults to "./bias_test_results.pdf".
        threshold (float, optional): Threshold for bias detection. Fits with
            | (bootstrap_mean - fit_value) / d | > threshold are flagged. For circular
            data (phase differences), d = π. All other data, d = fit_value.
            Defaults to 0.05.
        small_amp_threshold (float, optional): Threshold for fit fraction below which
            amplitudes are considered small, and thus ignored in bias tests. Phase
            differences who contain a small amplitude are also ignored. Defaults to
            0.05.

    Raises:
        KeyError: If any specified column is not found in either DataFrame.

    Returns:
        None: Generates a PDF file with histograms for failed bias tests.

    Todo:
        - handle acc_corrected (detected_events -> generated_events)
    """

    # check that columns exist in both dataframes
    for col in columns:
        if col not in fit_df.columns:
            raise KeyError(f"Column {col} not found in fit_df")
        if col not in bootstrap_df.columns:
            raise KeyError(f"Column {col} not found in bootstrap_df")

    # get phase differences for circular data handling
    phases = utils.get_phase_differences(fit_df)
    coherent_sums = utils.get_coherent_sums(fit_df)

    failed_dict = {}  # to store {fit_index: {column: bias_value}}
    for fit_index in bootstrap_df["fit_index"].unique():
        fit_value_row = fit_df.loc[fit_index]
        num_events = bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][
            "detected_events"
        ].mean()

        for col in columns:
            values = bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][col]
            fit_value = fit_value_row[col]
            # avoid division by zero in fractional bias calculation
            safe_fit_value = fit_value if fit_value != 0 else 1e-6

            if col in phases:
                # calculate circular mean and bias for phase differences
                bootstrap_mean = np.rad2deg(
                    scipy.stats.circmean(
                        np.deg2rad(values.to_numpy()),
                        high=np.pi,
                        low=-np.pi,  # type: ignore
                    )
                )

                # if either amplitude in the phase difference is small, skip bias test
                amp1, amp2 = col.split("_")
                amp1_ff = (
                    bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][
                        amp1
                    ].mean()
                    / num_events
                )
                amp2_ff = (
                    bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][
                        amp2
                    ].mean()
                    / num_events
                )
                if amp1_ff < small_amp_threshold or amp2_ff < small_amp_threshold:
                    continue

                bias = utils.circular_residual(
                    bootstrap_mean, fit_value, in_degrees=True
                )
                fractional_bias = abs(bias / 180.0)  # normalize by pi for phase diffs

            elif any(col in sublist for sublist in coherent_sums.values()):
                bootstrap_mean = values.mean()

                # skip bias test for small amplitudes
                fit_fraction = bootstrap_mean / num_events
                if fit_fraction < small_amp_threshold:
                    continue

                bias = (bootstrap_mean - fit_value) / safe_fit_value
                fractional_bias = abs(bias)

            else:
                bootstrap_mean = values.mean()
                bias = bootstrap_mean - fit_value
                fractional_bias = abs(bias / safe_fit_value)

            if fractional_bias > threshold:
                if fit_index not in failed_dict:
                    failed_dict[fit_index] = {}
                failed_dict[fit_index][col] = bias

    # Make plots for the failed tests
    if not failed_dict:
        print("All bias tests passed.")
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
                sharey=True,
            )
            axes = axes.flatten() if num_plots > 1 else [axes]

            dist_handle = None
            mean_handle = None
            fit_handle = None

            # loop through each failed column and make the histogram
            for ax, (col, bias) in zip(axes, col_dict.items()):
                values = bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][col]
                fit_value = fit_df.loc[fit_index][col]
                num_events = bootstrap_df.loc[bootstrap_df["fit_index"] == fit_index][
                    "detected_events"
                ].mean()

                if col in phases:
                    # For phase differences, match sign like in bootstrap_plotter
                    bootstrap_mean = np.rad2deg(
                        scipy.stats.circmean(
                            np.deg2rad(values.to_numpy()),
                            high=np.pi,
                            low=-np.pi,  # type: ignore
                        )
                    )

                    # Match the sign of fit_df to bootstrap mean
                    if fit_value * bootstrap_mean < 0:
                        fit_value = -fit_value

                    # Create histogram
                    dist_handle = ax.hist(
                        values,
                        bins=30,
                        alpha=0.7,
                        color="blue",
                        label="Bootstrap Samples",
                    )
                    mean_handle = ax.axvline(
                        bootstrap_mean,
                        color="blue",
                        linestyle="--",
                        linewidth=2,
                        label="Bootstrap Mean",
                    )
                    fit_handle = ax.axvline(
                        fit_value,
                        color="red",
                        linestyle="-",
                        linewidth=2,
                        label="Fit Value",
                    )
                    ax.set_title(f"{utils.convert_amp_name(col)}: bias={bias:.0f}°")
                    ax.set_xlabel("Phase Difference (degrees)")

                elif any(col in sublist for sublist in coherent_sums.values()):
                    # amplitude fit fraction
                    bootstrap_mean = values.mean()
                    fit_fraction = bootstrap_mean / num_events

                    dist_handle = ax.hist(
                        values / num_events,
                        bins=30,
                        alpha=0.7,
                        color="blue",
                        label="Bootstrap Samples",
                    )
                    mean_handle = ax.axvline(
                        bootstrap_mean / num_events,
                        color="blue",
                        linestyle="--",
                        linewidth=2,
                        label="Bootstrap Mean",
                    )
                    fit_handle = ax.axvline(
                        fit_value / num_events,
                        color="red",
                        linestyle="-",
                        linewidth=2,
                        label="Fit Value",
                    )
                    ax.set_title(f"{utils.convert_amp_name(col)}: " f"bias={bias:.1e}")
                    ax.set_xlabel("Fit Fraction")

                else:
                    # Regular histogram
                    bootstrap_mean = values.mean()
                    dist_handle = ax.hist(
                        values,
                        bins=30,
                        alpha=0.7,
                        color="blue",
                        label="Bootstrap Samples",
                    )
                    mean_handle = ax.axvline(
                        bootstrap_mean,
                        color="blue",
                        linestyle="--",
                        linewidth=2,
                        label="Bootstrap Mean",
                    )
                    fit_handle = ax.axvline(
                        fit_value,
                        color="red",
                        linestyle="-",
                        linewidth=2,
                        label="Fit Value",
                    )
                    ax.set_title(f"{utils.convert_amp_name(col)}: bias={bias:.1e}")

                ax.set_ylabel("Counts")

            # hide any unused subplots
            for ax in axes[num_plots:]:
                ax.set_visible(False)

            fig.suptitle(f"Bias Test Failures for Fit Index {fit_index}", fontsize=16)

            legend_handles = []
            legend_loc = "outside upper right" if len(axes) > 1 else "outside right"
            if dist_handle is not None and isinstance(dist_handle, (list, tuple)):
                # ax.hist returns (n, bins, patches) - we want the patches
                if len(dist_handle) == 3:
                    legend_handles.append(dist_handle[2][0])
                else:
                    legend_handles.append(dist_handle[0])
            elif dist_handle is not None and hasattr(dist_handle, "get_label"):
                legend_handles.append(dist_handle)

            if mean_handle is not None and hasattr(mean_handle, "get_label"):
                legend_handles.append(mean_handle)
            if fit_handle is not None and hasattr(fit_handle, "get_label"):
                legend_handles.append(fit_handle)

            if legend_handles:
                fig.legend(handles=legend_handles, loc=legend_loc)
            pdf.savefig(fig)
            plt.close(fig)
        print("Bias test results saved to:", output)

    return


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
