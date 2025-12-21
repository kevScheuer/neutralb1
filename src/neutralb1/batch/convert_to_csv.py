"""This script converts AmpTools .fit file(s) or its associated ROOT files into a csv.

This script is used for two fit result purposes:
1. To aggregate the AmpTools .fit files into a single .csv file for easier analysis.
a. The covariance and correlation matrices are also included in separate .csv files
2. To convert the ROOT files that the .fit files are based off of into a .csv file.

Behind the scenes, this runs compiled c++ scripts for either situation.

Todo:
    - Eventually, the fit converter should simply include the needed data info in it
    - return codes are overwritten and so won't communicate errors properly
"""

import argparse
import os
import subprocess
import sys
import tempfile

import neutralb1.utils as utils


def main() -> int:

    parser = create_parser()
    args = vars(parser.parse_args())

    if args["output"] and not args["output"].endswith(".csv"):
        args["output"] = args["output"] + ".csv"

    # Check if args["input"] is a file containing a list of result files, and save
    # the list of files to input_files
    input_files = []
    if (
        len(args["input"]) == 1
        and os.path.isfile(args["input"][0])
        and not args["input"][0].endswith(".fit")
        and not args["input"][0].endswith(".root")
    ):
        with open(args["input"][0], "r") as file:
            input_files = [line.strip() for line in file if line.strip()]
    else:
        input_files = args["input"]

    # Check if all input files exist, and expand its full path if just a file name
    print("Checking if all input files exist...")
    for file in input_files:
        if not os.path.exists(file):
            raise FileNotFoundError(f"The file {file} does not exist")
        if not os.path.isabs(file):
            input_files[input_files.index(file)] = os.path.abspath(file)

    if not all(file.endswith(".fit") for file in input_files):
        raise ValueError("All input files must be .fit files")

    # sort the input files
    input_files = (
        utils.sort_input_files(input_files, args["sort_index"])
        if args["sorted"]
        else input_files
    )

    if args["preview"]:
        print("Files that will be processed:")
        for file in input_files:
            print(f"\t{file}")
        return 0

    # create a tempfile that contains the list of input files
    # this seems to improve the speed of subprocess.Popen
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as temp_file:
        temp_file.write("\n".join(input_files))
        temp_file_path = temp_file.name
    print(f"Temp file created at {temp_file_path}")

    # convert flags into bool integers for the c++ scripts to interpret
    is_acceptance_corrected = 1 if args["acceptance_corrected"] else 0
    extract_data = 0 if args["no_data"] else 1

    # setup c++ command with appropriate arguments
    if args["data_only"]:
        output_file_name = "data.csv" if not args["output"] else args["output"]
        command = [
            f"extract_bin_info",
            f"{temp_file_path}",
            f"{output_file_name}",
            f"{args['mass_branch']}",
        ]
    elif args["moments"]:
        output_file_name = (
            "projected_moments.csv" if not args["output"] else args["output"]
        )
        command = [
            f"project_moments",
            f"{temp_file_path}",
            f"{output_file_name}",
        ]
    else:
        output_file_name = "fits.csv" if not args["output"] else args["output"]
        command = [
            f"extract_fit_results",
            f"{temp_file_path}",
            f"{output_file_name}",
            f"{is_acceptance_corrected}",
            f"{extract_data}",
            f"{args['mass_branch']}",
        ]
    return_code = run_process(command, args["verbose"])

    if args["correlation"]:
        corr_output_name = output_file_name.replace(".csv", "_corr.csv")
        corr_command = [
            f"extract_corr_matrix",
            f"{temp_file_path}",
            f"{corr_output_name}",
        ]
        return_code = run_process(corr_command, args["verbose"])
    if args["covariance"]:
        cov_output_name = output_file_name.replace(".csv", "_cov.csv")
        cov_command = [
            f"extract_cov_matrix",
            f"{temp_file_path}",
            f"{cov_output_name}",
        ]
        return_code = run_process(cov_command, args["verbose"])

    return return_code


def run_process(command: list, is_verbose: bool) -> int:
    """Run a command with subprocess and print the output if requested"""
    print("Running command:", command)
    proc = subprocess.Popen(
        command,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if is_verbose:
        if proc.stdout is not None:
            for line in iter(proc.stdout.readline, ""):
                print(line, end="")
    proc.wait()  # wait for the process to finish and update the return code
    if proc.stderr is not None:
        stderr_output = proc.stderr.read()
        if stderr_output:
            print("Error while running command:")
            print(stderr_output, end="")
            return 1
    print("Process completed successfully")

    return 0


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--input",
        help=(
            "Input file(s). Also accepts path(s) with a wildcard '*' and finds all"
            " matching files. Can also accept a file containing a list of files"
        ),
        nargs="+",
    )
    parser.add_argument(
        "-o", "--output", default="", help="File name of output .csv file"
    )
    parser.add_argument(
        "-s",
        "--sorted",
        type=bool,
        default=True,
        help=(
            "Sort the input files by last number in the file name or path. Defaults"
            " to True, so that the index of each csv row matches the ordering of the"
            " input files"
        ),
    )
    parser.add_argument(
        "--sort-index",
        type=int,
        default=-1,
        help=(
            "Determines what number in the file path is used for sorting. Defaults to"
            " -1, so that the last number in the path is used. See "
            " utils.sort_input_files for details."
        ),
    )
    parser.add_argument(
        "-a",
        "--acceptance-corrected",
        action="store_true",
        help=(
            "When passed, the amplitude intensities are corrected for acceptance. These"
            " are the true 'generated' values with no detector effects. Defaults to"
            " False, or the 'reconstructed' values"
        ),
    )
    parser.add_argument(
        "-m",
        "--mass-branch",
        type=str,
        default="M4Pi",
        help=(
            "Name of branch for the final invariant mass combo of interest in the"
            " Amplitude Analysis. Note this is only applicable when attempting to"
            " create data csv's. Defaults to M4Pi"
        ),
    )
    parser.add_argument(
        "-p",
        "--preview",
        action="store_true",
        help=("When passed, print out the files that will be processed and exit."),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print out more information while running the script",
    )
    parser.add_argument(
        "--correlation",
        action=argparse.BooleanOptionalAction,
        help=(
            "When passed, the correlation matrix of each fit is included in a separate"
            " csv file. Pass --no-correlation to disable this. Defaults to None."
        ),
    )
    parser.add_argument(
        "--covariance",
        action=argparse.BooleanOptionalAction,
        help=(
            "When passed, the covariance matrix of each fit is included in a separate"
            " csv file. Pass --no-covariance to disable this. Defaults to None."
        ),
    )
    parser.add_argument(
        "--moments",
        action="store_true",
        help=(
            "When passed, the script will project the moments from the partial wave fit"
            " results into a csv file. The fit file must be a partial wave fit."
        ),
    )
    parser.add_argument(
        "--data-only",
        action="store_true",
        help=(
            "When passed, the script will only extract the data bin information from"
            " the fit files into a csv file."
        ),
    )
    parser.add_argument(
        "--no-data",
        action="store_true",
        help=(
            "When passed, extract_fit_results will not extract bin information. By"
            " default, bin information is extracted alongside fit results."
        ),
    )
    return parser


if __name__ == "__main__":
    try:
        exit_code = main()
        sys.exit(exit_code)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        sys.exit(1)
