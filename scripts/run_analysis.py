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
import glob
import os
import pathlib
import pickle
import subprocess
import sys
import warnings
from typing import Any, Dict

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from pypdf import PdfReader, PdfWriter

import neutralb1.utils as utils
from neutralb1.analysis.result import ResultManager


def main() -> int:

    args = parse_args()

    if args["input"].endswith(".pkl"):  # load previously saved preprocessed results
        if not os.path.exists(args["input"]):
            raise FileNotFoundError(f"File {args['input']} does not exist.")
        with open(args["input"], "rb") as f:
            data = pickle.load(f)
        result = ResultManager(**data)
        ac_str = (
            "_acceptance_corrected"
            if result.plot.intensity.is_fit_acceptance_corrected
            else ""
        )
        print("Loaded preprocessed results from pickle file.")
    else:
        # collect raw csv files
        script_dir = os.path.dirname(os.path.abspath(__file__))
        pathlib.Path.mkdir(pathlib.Path(f"{args['output']}/raw"), exist_ok=True)

        command = (
            f"{script_dir}/collect_all_fits.bash -i {args['input']}"
            f" -o {args['output']}/raw/"
        )
        command += " -p" if args["preview"] else ""
        command += " -f" if args["force"] else ""

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
        data_df = pd.read_csv(f"{args['output']}/raw/best_data.csv")

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

    print("Results ready.")

    if args["no_plots"]:
        print("Skipping plotting as per --no-plots flag.")
        return 0
    print("Generating standard plots...")
    save_standard_plots(result, args["output"], ac_str)
    print("Plotting completed.")

    # If input was a t-bin directory (not a .pkl file), also stitch angle PDFs
    if not args["input"].endswith(".pkl"):
        combined_angles_path = f"{args['output']}/plots/combined_angles.pdf"
        stitch_angle_pdfs(args["input"], combined_angles_path)

    return 0


def save_standard_plots(result: ResultManager, output_dir: str, ac_str: str) -> None:
    """Generate and save standard plots using the provided ResultManager.

    Args:
        result (ResultManager): The ResultManager instance containing fit results.
        output_dir (str): Directory to save the generated plots.
        ac_str (str): Suffix indicating if acceptance-corrected results are used.

    Returns:
        None, just saves plot pdfs to output_dir.
    """

    pathlib.Path.mkdir(pathlib.Path(f"{output_dir}/plots/"), exist_ok=True)

    matplotlib.use("Agg")

    # Intensity plots
    intensity = result.plot.intensity
    intensity.jp()
    plt.gcf().savefig(f"{output_dir}/plots/jp{ac_str}.pdf", bbox_inches="tight")

    intensity.waves()
    plt.gcf().savefig(f"{output_dir}/plots/waves{ac_str}.pdf", bbox_inches="tight")

    if result.proj_moments_df is not None:
        intensity.moments()
        plt.gcf().savefig(
            f"{output_dir}/plots/proj_moments{ac_str}.pdf", bbox_inches="tight"
        )
    # Diagnostic plots
    diagnostic = result.plot.diagnostic
    diagnostic.matrix()
    plt.gcf().savefig(f"{output_dir}/plots/matrix{ac_str}.pdf", bbox_inches="tight")
    return


def stitch_angle_pdfs(t_bin_dir: str, output_path: str) -> None:
    """Stitch together angle distribution PDFs from all mass bins in a t-bin.

    This function finds all PDF files matching the pattern
    "t_bin/mass*/distributions/angles.pdf", sorts them by mass bin, and combines only
    the first page of each PDF into one large PDF file.

    Args:
        t_bin_dir (str): Path to the t-bin directory containing mass bin subdirectories.
        output_path (str): Path where the combined PDF should be saved.

    Raises:
        FileNotFoundError: If no angle PDF files are found in the t-bin directory.

    Returns:
        None: Saves the combined PDF to the specified output path.
    """
    # Find all angle PDF files matching the pattern
    pdf_pattern = os.path.join(t_bin_dir, "mass*", "distributions", "angles.pdf")
    pdf_files = glob.glob(pdf_pattern)

    if not pdf_files:
        raise FileNotFoundError(
            f"No angle PDF files found matching pattern: {pdf_pattern}"
        )

    # Sort the PDF files using the same logic as collect_csv.py
    # This sorts by the last number in the path (mass bin number)
    sorted_pdf_files = utils.sort_input_files(pdf_files, position=-1)

    print(f"Found {len(sorted_pdf_files)} angle PDF files to combine:")

    # Create a PdfWriter object for the combined PDF
    pdf_writer = PdfWriter()

    # Process each PDF file
    for i, pdf_file in enumerate(sorted_pdf_files):
        # Extract mass bin from pdf_file path
        mass_bin = pathlib.Path(pdf_file).parent.parent.name
        with open(pdf_file, "rb") as file:
            pdf_reader = PdfReader(file)

            # Add only the first page of each PDF
            if len(pdf_reader.pages) > 0:
                first_page = pdf_reader.pages[0]
                pdf_writer.add_page(first_page)
                pdf_writer.add_outline_item(mass_bin, i)
            else:
                warnings.warn(f"{pdf_file} has no pages, skipping.", UserWarning)

    # Write the combined PDF to the output path
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "wb") as output_file:
        pdf_writer.write(output_file)

    print(f"Combined PDF saved to: {output_path}")


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
        "--no-plots",
        action="store_true",
        help="If passed, do not generate and save standard plots.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="If passed, overwrite existing output raw csv files.",
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
