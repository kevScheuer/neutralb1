"""Combine multiple CSV files into a single CSV file.

Raises:
    FileNotFoundError: one of the input files does not exist.
    ValueError: one of the input files is not a .csv file.

Returns:
    None: Writes the combined CSV file to specified output path.
"""

import argparse
import os

import pandas as pd

import neutralb1.utils as utils


def main(args: dict) -> None:
    # Check if args["output"] ends with .csv
    if args["output"] and not args["output"].endswith(".csv"):
        args["output"] += ".csv"
    elif not args["output"]:
        args["output"] = "combined.csv"

    # Check if args["input"] is a file containing a list of CSV files
    input_files = []
    if (
        len(args["input"]) == 1
        and os.path.isfile(args["input"][0])
        and not args["input"][0].endswith(".csv")
    ):
        with open(args["input"][0], "r") as file:
            input_files = [line.strip() for line in file if line.strip()]
    else:
        input_files = args["input"]

    # Ensure all input files exist and are .csv files
    if input_files == []:
        raise ValueError("No input files provided.")
    for f in input_files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"File not found: {f}")
        if not f.endswith(".csv"):
            raise ValueError(f"Input file must be a .csv file: {f}")

    # Sort the input files if requested
    input_files = (
        utils.sort_input_files(input_files, args["sort_index"])
        if args["sorted"]
        else input_files
    )

    # Only print out the files that will be processed if preview flag is passed
    if args["preview"]:
        print("Files that will be processed:")
        for file in input_files:
            print(f"\t{file}")
        return

    # Combine the CSV files
    combined_df = pd.concat([pd.read_csv(f) for f in input_files])
    combined_df.to_csv(args["output"], index=False)
    print(f"Combined CSV saved to {args['output']}")


def parse_args() -> dict:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        nargs="+",
        help="Path to the input CSV file(s).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="",
        help="Path to the output CSV file.",
    )
    parser.add_argument(
        "-s",
        "--sorted",
        type=bool,
        default=True,
        help="Sort the input files by file name. Defaults to True.",
    )
    parser.add_argument(
        "-p",
        "--preview",
        action="store_true",
        help="When passed, print out the files that will be processed and exit.",
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
    return vars(parser.parse_args())


if __name__ == "__main__":
    main(parse_args())
