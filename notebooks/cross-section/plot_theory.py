"""DEPRECATED: kept for constants / to update later"""

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.integrate as integrate


def main():
    # load fit files
    df_data = pd.read_csv(f"mass_1.20-1.25.csv", index_col="index")
    df_theory = pd.read_csv(f"piExch-B1neutral.txt", sep=r"\s+", names=["t", "dsigma"])
    bin_df = pd.read_csv(f"t-bin-mean_m-1.20-1.25.csv")

    # get fraction of b1 breit wigner covered by the mass bin
    total_int = integrate.quad(breit_wigner, 0, math.inf, args=(1.235, 0.142))
    bin_int = integrate.quad(breit_wigner, 1.20, 1.25, args=(1.235, 0.142))
    frac = 1 / (bin_int[0] / total_int[0])

    # plot params
    LUM = 117858000  # ub^{-1} (old 125000000)
    BF = 0.892  # omega->3pi branching fraction (from pdg)

    t_bins = bin_df["average"]
    bin_error = bin_df["error"]
    t_width = np.array([0.1, 0.1, 0.2, 0.4])

    # change some plot parameters to look better
    plt.rcParams["figure.figsize"] = (14, 8)
    plt.rcParams["xtick.labelsize"] = 18
    plt.rcParams["ytick.labelsize"] = 18
    plt.rcParams["xtick.major.width"] = 2.0
    plt.rcParams["xtick.minor.width"] = 1.8
    plt.rcParams["ytick.major.width"] = 2.0
    plt.rcParams["ytick.minor.width"] = 1.8

    # write just the negative reflectivity results to a new csv if people want to see it
    new_df = pd.DataFrame(
        {
            "data": (df_data["m1p"]) / (BF * LUM * t_width) * frac,
            "errors": (df_data["m1p_err"]) / (BF * LUM * t_width) * frac,
        }
    )

    new_df.to_csv("GlueX_results.csv")

    # plot coherent sum of negative refl m-projections
    plt.errorbar(
        t_bins,
        new_df["data"],
        new_df["errors"],
        bin_error,
        color="#004D80",
        marker=".",
        markersize=16,
        linestyle="",
        label="Unnatural (GlueX Measurement)",
    )
    plt.plot(
        df_theory["t"],
        df_theory["dsigma"],
        color="#0076BA",
        marker="",
        linestyle="-",
        linewidth=2.0,
        label="Unnatural (JPAC Calculation)",
    )

    # Cosmetics
    plt.xlabel(r"$-t$ (GeV)$^{{2}}$", loc="right", fontsize=22)
    plt.ylabel(r"$d\sigma$ ($\mu$b) / $dt$ (GeV)$^2$", loc="top", fontsize=22)
    plt.gca().set_ylim(bottom=0)

    plt.legend(fontsize=20)
    plt.savefig(
        f"cross_section_theory_comparison.pdf", format="pdf", bbox_inches="tight"
    )

    return


def breit_wigner(x, m, width):
    gamma = math.sqrt(m**2 * (m**2 + width**2))

    k = (2 * math.sqrt(2) * m * width * gamma) / (math.pi * math.sqrt(m**2 + gamma))

    result = k / ((x**2 - m**2) ** 2 + m**2 * width**2)

    return result


if __name__ == "__main__":
    main()
