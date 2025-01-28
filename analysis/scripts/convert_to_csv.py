"""This script converts AmpTools .fit file(s) or its associated ROOT files into a csv.

This script is used for two fit result purposes:
1. To aggregate the AmpTools .fit files into a single .csv file for easier analysis.
2. To convert the ROOT files that the .fit files are based off of into a .csv file.
Behind the scenes, this script calls a ROOT macro for either situation.
"""

import argparse
import os
import subprocess

import utils


def main(args: dict) -> None:

    # ERROR / VALUE HANDLING
    if not os.environ["ROOTSYS"]:
        raise EnvironmentError(
            "ROOTSYS path is not loaded. Make sure to run 'source setup_gluex.csh'\n"
        )

    if args["output"] and not args["output"].endswith(".csv"):
        args["output"] = args["output"] + ".csv"

    # Check if args["input"] is a file containing a list of result files
    input_files = []
    if len(args["input"]) == 1 and os.path.isfile(args["input"][0]):
        with open(args["input"][0], "r") as file:
            input_files = [line.strip() for line in file if line.strip()]
    else:
        input_files = args["input"]

    for f in input_files:
        if not os.path.exists(f):
            print(f"File {f} does not exist, exiting")
            return

    if all(file.endswith(".fit") for file in input_files):
        file_type = "fit"
    elif all(file.endswith(".root") for file in input_files):
        file_type = "root"
    else:
        raise ValueError(
            "All input files must be of the same type: either .fit or .root files"
        )

    # sort the input files based off the last number in the file name or path
    input_files = (
        utils.sort_input_files(input_files, args["sort_index"])
        if args["sorted"]
        else input_files
    )

    if args["preview"]:
        print("Files that will be processed:")
        for file in input_files:
            print(f"\t{file}")
        return

    # get the script directory to properly call the script with the right path
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # hand off the files to the macro as a single space-separated string
    input_files = " ".join(input_files)

    is_acceptance_corrected = 1 if args["acceptance_corrected"] else 0
    if file_type == "fit":
        output_file_name = "fits.csv" if not args["output"] else args["output"]
        command = (
            f'.x {script_dir}/extract_fit_results.cc("{input_files}",'
            f' "{output_file_name}", {is_acceptance_corrected})\n'
        )
    elif file_type == "root":
        output_file_name = "data.csv" if not args["output"] else args["output"]
        command = (
            f'.x {script_dir}/extract_bin_info.cc("{input_files}",'
            f" \"{output_file_name}\", \"{args['mass_branch']}\")\n"
        )
    else:
        raise ValueError("Invalid type. Must be either 'fit' or 'root'")

    # call the ROOT macro
    proc = subprocess.Popen(
        ["root", "-n", "-l", "-b"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if file_type == "fit":
        proc.stdin.write(".x loadAmpTools.C\n")  # load the AmpTools libraries for fits
    proc.stdin.write(command)
    proc.stdin.flush()
    stdout, stderr = proc.communicate()
    proc.stdin.close()
    proc.wait()
    print(stdout)
    print(stderr)

    if not stderr:
        print(f"Results successfully written to: {output_file_name}")

    return


def parse_args() -> dict:
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
        "-o", "--output", default="", help="File name of output .csv file"
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
            " create csv's for ROOT data files. Defaults to M4Pi"
        ),
    )
    parser.add_argument(
        "-p",
        "--preview",
        action="store_true",
        help=("When passed, print out the files that will be processed and exit."),
    )
    return vars(parser.parse_args())


if __name__ == "__main__":
    args = parse_args()
    main(args)
