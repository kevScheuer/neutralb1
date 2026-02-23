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
import numpy.typing as npt
from matplotlib.projections.polar import PolarAxes
from uncertainties import unumpy as unp

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

    for col in ["p1p0S", "p1mpP", phase]:
        normalized_residuals = normalized_residuals_with_errors(
            nominal_result, systematic_results, col
        )

        fig = plot(
            mass_bins,
            normalized_residuals,
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


def normalized_residuals_with_errors(
    nominal_result: ResultManager, systematic_results: list[ResultManager], column: str
) -> list[np.ndarray]:
    """Compute normalized residuals for each systematic result

    This function computes the normalized residuals (nominal - systematic) / nominal for
    a specified column across all mass bins, and returns a list of unumpy arrays
    containing the residuals with propagated uncertainties.


    Args:
        nominal_result (ResultManager): Best fit result to compare to
        systematic_results (list[ResultManager]): result from a systematic variation
        column (str): the column in the fit_df to compute residuals for (e.g. "p1p0S")

    Returns:
        list[np.ndarray]: list of normalized residual unumpy arrays for each systematic
            result, in the same order as the input list.
    """

    assert (
        column in nominal_result.fit_df.columns
    ), f"Column '{column}' not found in nominal result dataframe"

    # importantly, we use the MINUIT error from the best fits
    nominal_array = unp.uarray(
        nominal_result.fit_df[column].to_numpy(),
        nominal_result.fit_df[f"{column}_err"].to_numpy(),
    )

    all_residuals: list[np.ndarray] = []

    for sys_result in systematic_results:
        assert (
            column in sys_result.fit_df.columns
        ), f"Column '{column}' not found in systematic result dataframe"

        sys_array = unp.uarray(
            sys_result.fit_df[column].to_numpy(),
            sys_result.fit_df[f"{column}_err"].to_numpy(),
        )

        # Compute normalized residuals for this systematic variation
        safe_denom = np.where(
            np.abs(unp.nominal_values(nominal_array)) > 1e-12,
            unp.nominal_values(nominal_array),
            np.nan,
        )
        residual = (nominal_array - sys_array) / safe_denom
        all_residuals.append(residual)

    return all_residuals


def plot(
    mass_bins: np.ndarray,
    normalized_residuals: list[np.ndarray],
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

    # here we create a mask for each label that identifies the largest residuals across
    # subgroups with the same base variable, for each mass bin
    largest_mask_dict = largest_residual_mask(labels, normalized_residuals)

    # determine a radial scale factor to keep residuals within a reasonable distance
    # from their baseline rings
    max_ring_excursion = 0.22  # max radial distance from base ring for largest residual
    largest_abs_residual = max(
        np.max(np.abs(unp.nominal_values(residual)))
        for residual in normalized_residuals
    )
    radial_scale = (
        max_ring_excursion / largest_abs_residual if largest_abs_residual > 0 else 1.0
    )

    ring_idx = 0
    for res, label in zip(normalized_residuals, labels):
        base_r = ring_radius[ring_idx]

        # draw baseline ring, denotes zero residual
        ax.plot(
            np.linspace(theta_min, theta_max, 400),
            np.full(400, base_r),
            color="0.75",
            linewidth=1.0,
        )

        # draw normalized residuals, with transparent error bands
        ax.plot(
            theta,
            base_r + radial_scale * unp.nominal_values(res),
            "-o",
            color=colors[ring_idx],
            linewidth=1.2,
            markersize=2,
        )
        ax.fill_between(
            theta,
            base_r + radial_scale * (unp.nominal_values(res) - unp.std_devs(res)),
            base_r + radial_scale * (unp.nominal_values(res) + unp.std_devs(res)),
            color=colors[ring_idx],
            alpha=0.15,
        )

        # draw guide lines for significance threshold (e.g. 4 sigma)
        significance_threshold = 4.0
        guide_delta = radial_scale * significance_threshold * unp.std_devs(res)
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
            xy=(theta_min, base_r + 0.05),
            xytext=(-10, 0),
            textcoords="offset points",
            ha="right",
            va="center",
            fontsize=9,
            fontweight="bold",
            color=colors[ring_idx],
        )

        # add star markers for largest residuals in subgroups with same base variable
        if label in largest_mask_dict:
            mask = largest_mask_dict[label]
            theta_masked = theta[mask]
            base_r_masked = base_r + radial_scale * unp.nominal_values(res)[mask]
            ax.plot(
                theta_masked,
                base_r_masked,
                marker="*",
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


def largest_residual_mask(
    labels: list[str], residuals: list[np.ndarray]
) -> dict[str, np.ndarray]:
    """Mask each label to identify biggest residuals within subgroup for each mass bin

    Args:
        labels (list[str]): List of labels, e.g. ["chi2>5", "chi2>6", "Pz>5"].
        residuals (list[np.ndarray]): List of residual arrays corresponding to each
            label. Arrays should have shape (n_bins,) and be in the same order as
            labels.

    Returns:
        dict[str, np.ndarray]: Mapping from label to boolean mask array where True
            indicates the largest residuals for that label's base variable group at each
            mass bin.
    """
    common_base_vars = find_common_base_variables(labels)
    largest_mask: dict[str, np.ndarray] = {}

    for base_var, group_labels in common_base_vars.items():
        group_residuals = np.array(
            [res for res, label in zip(residuals, labels) if label in group_labels]
        )
        max_indices = np.argmax(abs(group_residuals), axis=0)
        for idx, label in enumerate(group_labels):
            largest_mask[label] = max_indices == idx

    return largest_mask


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

    return vars(parser.parse_args())


if __name__ == "__main__":
    main()
