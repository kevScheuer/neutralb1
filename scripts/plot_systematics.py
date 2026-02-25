#!/usr/bin/env python3
"""Plot residuals from unumpy arrays on concentric polar rings.

Each variable is represented as a concentric ring and each mass bin is mapped to
an angular location. Residuals are shown as radial offsets from each variable's
baseline ring, with radial error bars representing residual uncertainty.

"""

from __future__ import annotations

import argparse
import pickle as pkl
from typing import Any, cast

import matplotlib
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.projections.polar import PolarAxes

import neutralb1.utils as utils
from neutralb1.analysis.result import ResultManager

matplotlib.use("Agg")  # Use non-interactive backend for script usage


def main() -> None:
    args = parse_args()

    # ensure output has a filename
    if args["output"].endswith("/"):
        args["output"] = args["output"][:-1] + "/systematics_residuals"
    # ensure it doesn't already end with .pdf
    if args["output"].endswith(".pdf"):
        args["output"] = args["output"][:-4]

    # load in our data
    with open(args["nominal"], "rb") as f:
        data = pkl.load(f)
        nominal_result = ResultManager(**data)

    systematic_results: list[ResultManager] = []
    for sys_file in args["systematics"]:
        with open(sys_file, "rb") as f:
            data = pkl.load(f)
            sys_result = ResultManager(**data)
            systematic_results.append(sys_result)

    print(f"Nominal result file: {args['nominal']}")
    print("Systematic result files and their groups:")
    for res, label in zip(args["systematics"], args["labels"]):
        print(f"\t{res} -> {label}")

    mass_bins = nominal_result.data_df["m_center"].to_numpy()
    phase = nominal_result.phase_difference_dict[("p1p0S", "p1mpP")]
    t_low, t_high = nominal_result.get_t_edges()

    # for creating the csv of systematic uncertainties. the list is in the same order
    # as the mass bins
    col_to_max_abs_residuals: dict[str, list[float]] = {}

    for col in ["p1p0S", "p1mpP", phase]:
        residuals, res_errors = residuals_with_errors(
            nominal_result, systematic_results, col
        )

        fig = plot(
            mass_bins,
            residuals,
            res_errors,
            args["labels"],
            t_low,
            t_high,
            col,
        )

        fig.savefig(
            f"{args['output']}_{col}_t_{t_low}_{t_high}.pdf",
            dpi=200,
            bbox_inches="tight",
        )
        print(f"Saved plot to: {args['output']}_{col}_t_{t_low}_{t_high}.pdf")

        if args["csv"]:
            significance_threshold = 4.0
            all_res = np.array(residuals)  # (n_variations, n_bins)
            all_err = np.array(res_errors)  # (n_variations, n_bins)
            sig_mask = np.abs(all_res / all_err) > significance_threshold
            # largest significant residual per bin; 0 where no variation is significant
            masked_abs_res = np.where(sig_mask, np.abs(all_res), 0.0)
            col_to_max_abs_residuals[col] = masked_abs_res.max(axis=0).tolist()

    if args["csv"]:
        csv_df = pd.DataFrame(
            {
                "m_center": mass_bins,
                **col_to_max_abs_residuals,
            }
        )
        csv_df.to_csv(f"{args['output']}_t_{t_low}_{t_high}.csv", index=False)
        print(
            f"Saved CSV of significant unnormalized residuals to:"
            f" {args['output']}_t_{t_low}_{t_high}.csv"
        )


def find_common_base_variables(labels: list[str]) -> dict[str, list[str]]:
    """Group labels by their base variable name.

    Args:
        labels (list[str]): List of labels, e.g. ["chi2>5", "chi2>6", "Pz>5"].

    Returns:
        dict[str, list[str]]: Mapping from base variable to systematic variations, e.g.
            {
                "chi2": ["chi2>5", "chi2>6"],
                "Pz": ["Pz>5"]
            }
    """
    groups: dict[str, list[str]] = {}
    for label in labels:
        if any(sep in label for sep in [">", "<", "="]):
            var_name = label.split(">", 1)[0].split("<", 1)[0].split("=", 1)[0]
        else:
            var_name = label
        groups.setdefault(var_name, []).append(label)
    return groups


def residuals_with_errors(
    nominal_result: ResultManager, systematic_results: list[ResultManager], column: str
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Compute residuals for each systematic result

    This function computes the residuals (systematic - nominal) for a specified column
    across all mass bins, and returns a list of arrays containing the residuals
    with propagated uncertainties. Note that only MINUIT errors are used.

    Args:
        nominal_result (ResultManager): Best fit result to compare to
        systematic_results (list[ResultManager]): list of systematic results to compare
        column (str): the column in the fit_df to compute residuals for (e.g. "p1p0S")

    Returns:
        tuple[list[np.ndarray], list[np.ndarray]]: tuple of lists of residual values and
            errors across all systematic variations, in order of labels
    """
    assert (
        column in nominal_result.fit_df.columns
    ), f"Column '{column}' not found in nominal result dataframe"

    # since the columns are amplitudes and phase differences, it is safe
    # to take the absolute value. This is because amplitudes are always positive by
    # definition, and the phase differences have a sign ambiguity in the model that we
    # must account for in the residual calculation
    nominal_values = np.abs(nominal_result.fit_df[column].to_numpy())
    nominal_errors = np.abs(nominal_result.fit_df[f"{column}_err"].to_numpy())

    all_residuals: list[np.ndarray] = []
    all_errors: list[np.ndarray] = []

    vectorized_circular_residual = np.vectorize(utils.circular_residual)
    for sys_result in systematic_results:
        assert (
            column in sys_result.fit_df.columns
        ), f"Column '{column}' not found in systematic result dataframe"

        sys_values = np.abs(sys_result.fit_df[column].to_numpy())
        sys_errors = np.abs(sys_result.fit_df[f"{column}_err"].to_numpy())

        if column in nominal_result.phase_differences:
            # for phase differences, we need to account for its circular nature
            residual = vectorized_circular_residual(
                sys_values, nominal_values, in_degrees=True, low=-180, high=180
            )
        else:
            residual = sys_values - nominal_values
        residual_err = np.sqrt(
            np.abs(np.square(sys_errors) - np.square(nominal_errors))
        )

        all_residuals.append(residual)
        all_errors.append(residual_err)

    return all_residuals, all_errors


def plot(
    mass_bins: np.ndarray,
    residuals: list[np.ndarray],
    residual_errors: list[np.ndarray],
    labels: list[str],
    t_low: float,
    t_high: float,
    column: str,
) -> matplotlib.figure.Figure:

    # convert mass bin points to angles for polar plot
    n_bins = len(mass_bins)
    theta_min = 0.0
    theta_max = 1.5 * np.pi
    theta = np.linspace(theta_min, theta_max, n_bins, endpoint=True)

    # determine spacing of rings
    n_rings = len(labels)
    min_ring_radius = 2.0  # starting radius for innermost ring
    ring_radius = min_ring_radius + np.arange(n_rings, dtype=float)

    fig, ax = plt.subplots(figsize=(20, 18), subplot_kw={"projection": "polar"})
    ax = cast(PolarAxes, ax)
    ax.set_theta_offset(np.pi / 2.0)  # start from top
    ax.set_theta_direction(-1)  # clockwise
    ax.set_thetamin(np.degrees(theta_min))
    ax.set_thetamax(np.degrees(theta_max))

    # we want to map a color to each base variable (e.g. chi2, Pz, etc.), so here we
    # find the base variable names and make list of colors in the same order as labels
    color_map = plt.get_cmap("tab10")
    common_base_vars = find_common_base_variables(labels)
    colors = []
    for label in labels:
        for base_var, group_labels in common_base_vars.items():
            if label in group_labels:
                var_idx = list(common_base_vars.keys()).index(base_var)
                colors.append(color_map(var_idx % 10))
                break
        else:
            colors.append("black")  # default color if no base variable match

    # here we create a mask for each label that identifies the largest significant
    # residuals across all variations for each mass bin, so we can mark them with stars
    # on the plot
    significance_threshold = 4.0
    largest_sig_mask_dict = largest_sig_residual_mask(
        labels, residuals, residual_errors, significance_threshold
    )

    # determine a radial scale factor to keep residuals within a reasonable distance
    # from their baseline rings
    max_ring_excursion = 0.1  # max radial distance from base ring for largest residual
    largest_barlow = 4.0
    # TODO: this allows some mass bins to bleed into other rings, but its rare, and
    #   implementing this means heavily squishing all the ring guides towards the base
    # for res, res_err, label in zip(residuals, residual_errors, labels):
    #     mask = largest_sig_mask_dict[label]
    #     if np.any(mask):
    #         max_ratio = np.max(np.abs(res[mask] / res_err[mask]))
    #         if max_ratio > largest_barlow:
    #             largest_barlow = max_ratio
    radial_scale = max_ring_excursion / largest_barlow

    ring_idx = 0
    for res, res_err, label in zip(residuals, residual_errors, labels):
        base_r = ring_radius[ring_idx]

        # draw baseline ring, denotes zero residual
        ax.plot(
            np.linspace(theta_min, theta_max, 400),
            np.full(400, base_r),
            color="0.75",
            linewidth=1.0,
        )

        # draw residual / residual_errors
        ax.plot(
            theta,
            base_r + radial_scale * (res / res_err),
            "-o",
            color=colors[ring_idx],
            linewidth=1.2,
            markersize=2,
        )

        # draw guide lines for significance threshold
        guide_delta = np.full_like(res, radial_scale * significance_threshold)
        ax.plot(
            theta,
            base_r + guide_delta,
            color="red",
            linestyle=":",
            linewidth=1.0,
            alpha=0.8,
        )
        ax.plot(
            theta,
            base_r - guide_delta,
            color="red",
            linestyle=":",
            linewidth=1.0,
            alpha=0.8,
        )

        # label the rings
        ax.annotate(
            rf"{label}",
            xy=(theta_min, base_r + 0.01),  # TODO: adjust offset when all rings done
            xytext=(-10, 0),
            textcoords="offset points",
            ha="right",
            va="center",
            fontsize=9,
            fontweight="bold",
            color=colors[ring_idx],
        )

        # add star markers for the globally largest significant residual across all
        # variations at each mass bin
        if label in largest_sig_mask_dict:
            mask = largest_sig_mask_dict[label]
            theta_masked = theta[mask]
            base_r_masked = base_r + (radial_scale * res[mask] / res_err[mask])
            ax.plot(
                theta_masked,
                base_r_masked,
                marker="*",
                linestyle="",
                color=colors[ring_idx],
                markersize=6,
                markeredgecolor="black",
                markeredgewidth=0.6,
            )

        ring_idx += 1

    # draw radial lines for each mass bin over top of all lines
    for t in theta:
        ax.plot(
            [t, t],
            [0.0, ring_radius[-1] + 0.55],
            color="0.85",
            linewidth=0.8,
            zorder=0,
        )

    # plot cosmetics
    tick_labels = [f"{mass:.3f}" for mass in mass_bins]
    ax.set_xticks(theta)
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_ylim(0.0, ring_radius[-1] + 0.6)
    ax.set_yticks(ring_radius)
    ax.set_yticklabels([""] * len(ring_radius))

    # plot title according to t bin and column
    pretty_col_name = utils.convert_amp_name(column)
    ax.set_title(
        rf"${t_low} \leq -t \leq {t_high}$ GeV$^2$: {pretty_col_name}", va="bottom"
    )
    ax.grid(False)

    return fig


def largest_sig_residual_mask(
    labels: list[str],
    residuals: list[np.ndarray],
    residual_errors: list[np.ndarray],
    significance_threshold: float,
) -> dict[str, np.ndarray]:
    """Mask each label to identify the globally largest significant residual across all
    variations for each mass bin.

    A bin is only unmasked for the variation with the largest absolute residual if that
    residual also exceeds the significance threshold. Bins where no variation is
    significant are masked for all labels.

    Args:
        labels (list[str]): List of labels, e.g. ["chi2>5", "chi2>6", "Pz>5"].
        residuals (list[np.ndarray]): List of residual arrays corresponding to each
            label. Arrays should have shape (n_bins,) and be in the same order as
            labels.
        residual_errors (list[np.ndarray]): List of residual error arrays corresponding
            to each label, same shape and order as residuals.
        significance_threshold (float): Threshold for |residual / error| above which a
            residual is considered significant.

    Returns:
        dict[str, np.ndarray]: Mapping from label to boolean mask array where True
            indicates that label has the globally largest significant residual at each
            mass bin.
    """
    all_residuals = np.array(residuals)  # shape: (n_variations, n_bins)
    all_errors = np.array(residual_errors)  # shape: (n_variations, n_bins)

    # only consider bins where the residual is significant
    sig_mask = np.abs(all_residuals / all_errors) > significance_threshold

    # zero out non-significant residuals so argmax ignores them
    masked_abs_residuals = np.where(sig_mask, np.abs(all_residuals), 0.0)
    max_indices = np.argmax(masked_abs_residuals, axis=0)  # shape: (n_bins,)

    # suppress bins where no variation is significant (argmax would pick index 0)
    any_sig = sig_mask.any(axis=0)

    return {label: (max_indices == idx) & any_sig for idx, label in enumerate(labels)}


def parse_args() -> dict[str, Any]:
    parser = argparse.ArgumentParser(
        description=(
            "Create a polar plot of the normalized residuals between nominal fit"
            " results and systematic variations across mass bins, with concentric rings"
            " for each labelled systematic variation given."
        )
    )
    parser.add_argument(
        "-n",
        "--nominal",
        type=str,
        help=(
            "Path to nominal ResultsManager pickle file containing that all systematic"
            " files will be compared to"
        ),
    )
    parser.add_argument(
        "-s",
        "--systematics",
        type=str,
        nargs="+",
        help=(
            "Paths to systematic ResultsManager pickle files to compare against the"
            " nominal. Each file should correspond to a different variable subgroup"
            " (e.g. chi2>5, chi2>6, etc.) and will be plotted as a separate ring on the"
            " polar plot."
        ),
    )
    parser.add_argument(
        "-l",
        "--labels",
        type=str,
        nargs="+",
        help=(
            "Labels for each systematic file provided in --systematics, in the same"
            " order. These will be used to annotate the corresponding rings on the"
            " plot. Any labels with the same base variable name (e.g. 'chi2') will be"
            " grouped together and have an '*' marker for the largest residuals across"
            " subgroups at each mass bin."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="systematics_residuals",
        help=(
            "Path and base name for output plot file (without extension). The plot will"
            " be saved as {output}.pdf."
        ),
    )
    parser.add_argument(
        "--csv",
        action="store_true",
        help=(
            "If set, create a CSV file of the unnormalized residuals and their"
            " uncertainties. This is so they can be read in and used for plotting."
            " The table will be saved as {output}.csv."
        ),
    )

    return vars(parser.parse_args())


if __name__ == "__main__":
    main()
