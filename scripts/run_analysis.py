"""Run analysis on fit results across all mass bins.

This script aggregates fit results in a bin of -t, preprocesses the data,
and generates standard plots for analysis.

Example:
    >>> # run analysis on a t-bin directory and save results to output directory
    >>> python scripts/run_analysis.py -i /path/to/t_bin/ -o /path/to/output
    >>> # load previously saved preprocessed results and generate plots
    >>> python scripts/run_analysis.py -i preprocessed_results.pkl -o /path/to/output
"""

import argparse
import os
import pathlib
import pickle
import subprocess
import sys
from typing import Any, Dict

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from neutralb1.analysis.result import ResultManager


def main() -> int:

    args = parse_args()

    if args["input"].endswith(".pkl"):  # load previously saved preprocessed results
        if not os.path.exists(args["input"]):
            raise FileNotFoundError(f"File {args['input']} does not exist.")
        with open(args["input"], "rb") as f:
            data = pickle.load(f)
        result = ResultManager(**data)
    else:
        # collect raw csv files
        script_dir = os.path.dirname(os.path.abspath(__file__))
        pathlib.Path.mkdir(pathlib.Path(f"{args['output']}/raw"), exist_ok=True)

        command = (
            f"{script_dir}/collect_all_fits.bash -i {args['input']}"
            f" -o {args['output']}/raw/"
        )
        command += " -p" if args["preview"] else ""

        run_result = subprocess.run(command, shell=True, check=True)
        if run_result.returncode != 0:
            print(run_result.stdout)
            print(run_result.stderr)
            return run_result.returncode

        if args["preview"]:
            print("Preview mode: no files saved.")
            return 0

        print("Collected raw csv files.")

        # read in collected csv files
        ac_str = "_acceptance_corrected" if args["acceptance_corrected"] else ""
        fit_df = pd.read_csv(f"{args['output']}/raw/best{ac_str}.csv")
        data_df = pd.read_csv(f"{args['output']}/raw/data.csv")

        # following ones are optional
        proj_moments_df = (
            pd.read_csv(f"{args['output']}/raw/moments.csv")
            if os.path.exists(f"{args['output']}/raw/moments.csv")
            else None
        )
        randomized_df = (
            pd.read_csv(f"{args['output']}/raw/rand{ac_str}.csv")
            if os.path.exists(f"{args['output']}/raw/rand{ac_str}.csv")
            else None
        )
        randomized_proj_moments_df = (
            pd.read_csv(f"{args['output']}/raw/rand_moments.csv")
            if os.path.exists(f"{args['output']}/raw/rand_moments.csv")
            else None
        )
        bootstrap_df = (
            pd.read_csv(f"{args['output']}/raw/bootstrap{ac_str}.csv")
            if os.path.exists(f"{args['output']}/raw/bootstrap{ac_str}.csv")
            else None
        )
        bootstrap_proj_moments_df = (
            pd.read_csv(f"{args['output']}/raw/bootstrap_moments.csv")
            if os.path.exists(f"{args['output']}/raw/bootstrap_moments.csv")
            else None
        )
        truth_df = (
            pd.read_csv(f"{args['output']}/raw/truth{ac_str}.csv")
            if os.path.exists(f"{args['output']}/raw/truth{ac_str}.csv")
            else None
        )
        truth_proj_moments_df = (
            pd.read_csv(f"{args['output']}/raw/truth_moments.csv")
            if os.path.exists(f"{args['output']}/raw/truth_moments.csv")
            else None
        )

        # preprocess collected csv files
        result = ResultManager(
            fit_df=fit_df,
            data_df=data_df,
            proj_moments_df=proj_moments_df,
            randomized_df=randomized_df,
            randomized_proj_moments_df=randomized_proj_moments_df,
            bootstrap_df=bootstrap_df,
            bootstrap_proj_moments_df=bootstrap_proj_moments_df,
            truth_df=truth_df,
            truth_proj_moments_df=truth_proj_moments_df,
        )

        result.preprocess()

        if not args["no_save"]:
            pickle_file = f"{args['output']}/preprocessed_results{ac_str}.pkl"
            result.save(pickle_file, preprocess=False)  # preprocess already done above

    print("Results ready, plotting...")
    save_standard_plots(result, args["output"])
    print("Plotting completed.")

    return 0


def save_standard_plots(result: ResultManager, output_dir: str) -> None:
    """Generate and save standard plots using the provided ResultManager.

    Args:
        result (ResultManager): The ResultManager instance containing fit results.
        output_dir (str): Directory to save the generated plots.

    Returns:
        None, just saves plot pdfs to output_dir.
    """

    pathlib.Path.mkdir(pathlib.Path(f"{output_dir}/plots/"), exist_ok=True)

    matplotlib.use("Agg")

    # Intensity plots
    intensity = result.plot.intensity()
    intensity.jp()
    plt.gcf().savefig(f"{output_dir}/plots/jp.pdf", bbox_inches="tight")

    intensity.waves()
    plt.gcf().savefig(f"{output_dir}/plots/waves.pdf", bbox_inches="tight")

    if result.proj_moments_df is not None:
        intensity.moments()
        plt.gcf().savefig(f"{output_dir}/plots/proj_moments.pdf", bbox_inches="tight")

    # Diagnostic plots
    diagnostic = result.plot.diagnostic()
    diagnostic.matrix()
    plt.gcf().savefig(f"{output_dir}/plots/matrix.pdf", bbox_inches="tight")

    return


def parse_args() -> Dict[str, Any]:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate, preprocess, and plot standard fit results across all mass bins."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help=(
            "Path to the t directory containing mass bin subdirectories, or a pickle"
            " file containing previously saved preprocessed results."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="./",
        help="Output directory to save preprocessed results and plots.",
    )
    parser.add_argument(
        "--no-save",
        action="store_true",
        help="If passed, do not save preprocessed results to output directory.",
    )
    parser.add_argument(
        "-a",
        "--acceptance-corrected",
        action="store_true",
        help="If set, use acceptance-corrected fit results.",
    )
    parser.add_argument(
        "--preview",
        action="store_true",
        help="If set, run in preview mode (print out what files will be combined).",
    )

    return vars(parser.parse_args())


if __name__ == "__main__":
    try:
        exit_code = main()
    except Exception as e:
        print(f"Error occurred: {e}")
        exit_code = 1

    sys.exit(exit_code)
