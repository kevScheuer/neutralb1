#!/usr/bin/env python3
"""Plot residuals from unumpy arrays on concentric polar rings.

Each variable is represented as a concentric ring and each mass bin is mapped to
an angular location. Residuals are shown as radial offsets from each variable's
baseline ring, with radial error bars representing residual uncertainty.

Todo:
    - Later we will have this take in ResultManager objects for each variable subgroup,
        and a ResultManager for the nominal values and extract the relevant arrays for
        plotting.
    - Add a cli arg that will take in a dataframe column to be the value we'll plot.
"""

from __future__ import annotations

import argparse
from typing import Any, cast

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from matplotlib.projections.polar import PolarAxes
from uncertainties import unumpy as unp

matplotlib.use("Agg")  # Use non-interactive backend for script usage


def parse_args() -> dict[str, Any]:
    """Parse CLI arguments.

    Args:
            None: No positional arguments are required.

    Returns:
            dict[str, Any]: Parsed CLI arguments as a dictionary.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Create a polar residual plot from example unumpy arrays where "
            "variables are concentric rings and mass bins are angular sectors."
        )
    )
    parser.add_argument(
        "--save",
        type=str,
        default="",
        help="Optional output image path. If omitted, an interactive window is shown.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=7,
        help="Random seed for reproducible synthetic example data.",
    )
    parser.add_argument(
        "--n-mass-bins",
        type=int,
        default=40,
        help="Number of mass bins for the synthetic example.",
    )
    return vars(parser.parse_args())


def make_example_data(
    variables: list[str],
    subgroups: list[str],
    mass_values: np.ndarray,
    seed: int,
) -> tuple[dict[str, dict[str, npt.NDArray]], dict[str, dict[str, npt.NDArray]]]:
    """Build synthetic nominal and comparison unumpy arrays.

    Args:
                variables (list[str]): Variable names to generate.
                mass_values (np.ndarray): Mass-bin center values.
                seed (int): Random seed for reproducibility.

    Returns:
                tuple[dict[str, dict[str, npt.NDArray]], dict[str, dict[str, npt.NDArray]]]:
                                Two nested dictionaries keyed by variable and subgroup.
    """
    rng = np.random.default_rng(seed)
    nominal: dict[str, dict[str, npt.NDArray]] = {}
    comparison: dict[str, dict[str, npt.NDArray]] = {}

    for v_idx, var in enumerate(variables):
        nominal[var] = {}
        comparison[var] = {}
        for s_idx, subgroup in enumerate(subgroups):
            subgroup_shift = 0.06 * s_idx
            base = (
                1.5
                + 0.2 * v_idx
                + subgroup_shift
                + 0.25 * np.sin(2 * np.pi * mass_values)
            )
            nominal_sigma = (
                0.025
                + 0.005 * v_idx
                + 0.003 * s_idx
                + 0.005 * rng.random(len(mass_values))
            )
            nominal[var][subgroup] = unp.uarray(base, nominal_sigma)

            shift = 0.02 * np.cos(2 * np.pi * mass_values + 0.4 * v_idx + 0.3 * s_idx)
            noise = rng.normal(
                loc=0.0,
                scale=0.01 + 0.003 * v_idx + 0.002 * s_idx,
                size=len(mass_values),
            )
            comp_vals = base + shift + noise
            comp_sigma = (
                0.03
                + 0.006 * v_idx
                + 0.002 * s_idx
                + 0.004 * rng.random(len(mass_values))
            )
            comparison[var][subgroup] = unp.uarray(comp_vals, comp_sigma)

    return nominal, comparison


def residuals_with_errors(
    nominal: dict[str, dict[str, npt.NDArray]],
    comparison: dict[str, dict[str, npt.NDArray]],
) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, dict[str, np.ndarray]]]:
    """Compute residual arrays and propagated uncertainties.

    Args:
            nominal (dict[str, npt.NDArray]): Nominal measurements by variable.
            comparison (dict[str, npt.NDArray]): Measurements to compare by variable.

    Returns:
            tuple[dict[str, np.ndarray], dict[str, np.ndarray]]: Residual means and
                    standard deviations for each variable.
    """
    residual_mean: dict[str, dict[str, np.ndarray]] = {}
    residual_sigma: dict[str, dict[str, np.ndarray]] = {}

    for var in nominal:
        residual_mean[var] = {}
        residual_sigma[var] = {}
        for subgroup in nominal[var]:
            nominal_vals = unp.nominal_values(nominal[var][subgroup])
            safe_denom = np.where(np.abs(nominal_vals) > 1e-12, nominal_vals, np.nan)

            residual = nominal[var][subgroup] - comparison[var][subgroup]
            residual_mean[var][subgroup] = unp.nominal_values(residual) / safe_denom
            residual_sigma[var][subgroup] = unp.std_devs(residual) / np.abs(safe_denom)

    return residual_mean, residual_sigma


def plot_radial_residuals(
    variables: list[str],
    subgroups: list[str],
    mass_values: np.ndarray,
    residual_mean: dict[str, dict[str, np.ndarray]],
    residual_sigma: dict[str, dict[str, np.ndarray]],
):
    """Plot residuals and errors on concentric polar rings.

    Args:
            variables (list[str]): Variables in ring order from center outward.
            mass_values (np.ndarray): Mass-bin centers.
            residual_mean (dict[str, np.ndarray]): Residual means by variable.
            residual_sigma (dict[str, np.ndarray]): Residual uncertainties by variable.

    Returns:
            plt.Figure: The generated matplotlib figure.
    """
    n_bins = len(mass_values)
    theta_min = 0.0
    theta_max = 1.5 * np.pi
    theta = np.linspace(theta_min, theta_max, n_bins, endpoint=True)
    n_series = len(variables) * len(subgroups)
    min_ring_radius = 2.0
    ring_radius = min_ring_radius + np.arange(n_series, dtype=float)

    significance_threshold = 4.0
    max_ring_excursion = 0.22

    all_vals: list[np.ndarray] = []
    for var in variables:
        for subgroup in subgroups:
            all_vals.append(
                np.abs(residual_mean[var][subgroup])
                + significance_threshold * residual_sigma[var][subgroup]
            )
    max_residual_extent = float(np.max(np.concatenate(all_vals)))
    radial_scale = (
        max_ring_excursion / max_residual_extent if max_residual_extent > 0 else 1.0
    )

    fig, ax = plt.subplots(figsize=(20, 18), subplot_kw={"projection": "polar"})
    ax = cast(PolarAxes, ax)
    ax.set_theta_offset(np.pi / 2.0)
    ax.set_theta_direction(-1)
    ax.set_thetamin(np.degrees(theta_min))
    ax.set_thetamax(np.degrees(theta_max))

    color_map = plt.get_cmap("tab10")
    ring_idx = 0

    for v_idx, var in enumerate(variables):
        color = color_map(v_idx % 10)
        group_values = np.vstack(
            [residual_mean[var][subgroup] for subgroup in subgroups]
        )
        all_nan = np.all(np.isnan(group_values), axis=0)
        safe_group_values = np.where(np.isnan(group_values), -np.inf, group_values)
        max_subgroup_idx = np.argmax(safe_group_values, axis=0)

        for s_idx, subgroup in enumerate(subgroups):
            base_r = ring_radius[ring_idx]
            ring_idx += 1

            # draw baseline ring, denotes zero residual
            ax.plot(
                np.linspace(theta_min, theta_max, 400),
                np.full(400, base_r),
                color="0.75",
                linewidth=1.0,
            )

            r_points = base_r + radial_scale * residual_mean[var][subgroup]
            r_err = radial_scale * residual_sigma[var][subgroup]
            guide_delta = radial_scale * (
                significance_threshold * residual_sigma[var][subgroup]
            )

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

            ax.fill_between(
                theta,
                r_points - r_err,
                r_points + r_err,
                color=color,
                alpha=0.15,
            )
            ax.plot(
                theta,
                r_points,
                "-o",
                color=color,
                linewidth=1.2,
                markersize=2,
                label=f"{var}{subgroup}",
            )

            max_mask = (max_subgroup_idx == s_idx) & (~all_nan)
            if np.any(max_mask):
                ax.plot(
                    theta[max_mask],
                    r_points[max_mask],
                    "*",
                    color=color,
                    markersize=6,
                    markeredgecolor="black",
                    markeredgewidth=0.6,
                )

            ax.annotate(
                f"{var}{subgroup}",
                xy=(theta_min, base_r + 0.05),
                xytext=(-10, 0),
                textcoords="offset points",
                ha="right",
                va="center",
                fontsize=9,
                fontweight="bold",
                color=color,
            )

    for t in theta:
        ax.plot(
            [t, t],
            [0.0, ring_radius[-1] + 0.55],
            color="0.85",
            linewidth=0.8,
            zorder=0,
        )

    tick_labels = [f"{mass:.3f}" for mass in mass_values]
    ax.set_xticks(theta)
    ax.set_xticklabels(tick_labels, fontsize=8)

    ax.set_ylim(0.0, ring_radius[-1] + 0.6)
    ax.set_yticks(ring_radius)
    ax.set_yticklabels([""] * len(ring_radius))
    ax.set_title(
        "Normalized residuals by variable subgroups (rings) and mass bins (angles)\n"
        "point radius = ring baseline + scaled (nominal - variation) / nominal",
        va="bottom",
    )
    ax.grid(False)

    return fig


def main() -> None:
    """Run example residual generation and plotting."""
    args = parse_args()

    variables = ["chi2", "Pz", "unusedE", "unused Tracks", "Shower Quality"]
    subgroups = [">5", ">6", ">7"]
    n_mass_bins = int(args["n_mass_bins"])
    mass_values = np.linspace(1.0, 1.8, n_mass_bins, endpoint=True)

    nominal, comparison = make_example_data(
        variables=variables,
        subgroups=subgroups,
        mass_values=mass_values,
        seed=int(args["seed"]),
    )
    residual_mean, residual_sigma = residuals_with_errors(nominal, comparison)

    fig = plot_radial_residuals(
        variables=variables,
        subgroups=subgroups,
        mass_values=mass_values,
        residual_mean=residual_mean,
        residual_sigma=residual_sigma,
    )

    if args["save"]:
        fig.savefig(args["save"], dpi=200, bbox_inches="tight")
        print(f"Saved plot to: {args['save']}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
