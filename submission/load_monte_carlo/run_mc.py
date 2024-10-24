"""Create AmpTools Input trees from signal and phasespace MC for a given run period."""

import argparse
import subprocess

SIGNAL_PATH = "/cache/halld/home/jrsteven/REQUESTED_MC/"
SIGNAL_VERSION = "_ver03.1"
PHASESPACE_PATH = "/cache/halld/home/jrsteven/REQUESTED_MC/"
PHASESPACE_VERSION = "_ver03"
TREE_NAME = "tree_pi0pi0pippim__B4"


def main(args: dict) -> None:
    # map the run period to the initial integer in the run number
    run_period_name_to_number = {
        "2017_01": 3,
        "2018_01": 4,
        "2018_08": 5,
    }

    # ensure that phasespace files only use generated photons
    phasespace_option = args["mc_option"]
    if "noaccidental" not in phasespace_option:
        phasespace_option += "_noaccidental"

    # execute DSelector for each run period
    for run_period in args["run_periods"]:
        run_number = run_period_name_to_number[run_period]
        signal_name = (
            f"omegapi_massDepFit_{run_period}{SIGNAL_VERSION}/{TREE_NAME}/merged/"
        )
        phasespace_name = (
            f"omegapi_phasespace_{run_period}{PHASESPACE_VERSION}/{TREE_NAME}/merged/"
        )
        reconstructed_command = [
            "root.exe",
            "-l",
            "-b",
            "-q",
            "-n",
            (
                "runSelector.C("
                f'"{run_number:02d}",'
                f' "{SIGNAL_PATH}/{signal_name}",'
                f' "{args["mc_option"]}"'
                ")"
            ),
        ]
        phasespace_command = [
            "root.exe",
            "-l",
            "-b",
            "-q",
            "-n",
            (
                "runSelector.C("
                f'"{run_number:02d}",'
                f' "{PHASESPACE_PATH}/{phasespace_name}",'
                f' "{phasespace_option}"'
                ")"
            ),
        ]

        mv_reconstructed_command = [
            "mv",
            "AmpToolsInputTree.root",
            (
                "AmpToolsInputTree_sum_PARA_0_"
                f"{run_period}{SIGNAL_VERSION}_mc{args['mc_option']}.root"
            ),
        ]
        mv_phasespace_command = [
            "mv",
            "AmpToolsInputTree.root",
            (
                "anglesOmegaPiPhaseSpaceAcc_"
                f"{run_period}{PHASESPACE_VERSION}_{phasespace_option}.root"
            ),
        ]

        # run the reconstructed signal MC
        subprocess.call(reconstructed_command)
        subprocess.call(mv_reconstructed_command)

        # run the reconstructed phasespace MC
        if args["make_phasespace_trees"]:
            subprocess.call(phasespace_command)
            subprocess.call(mv_phasespace_command)

        if not args["make_gen_trees"]:
            continue

        # replace some of the command arguments to run the generated trees
        gen_reconstructed_command = reconstructed_command[:-1]
        gen_reconstructed_command.append(
            reconstructed_command[-1]
            .replace(signal_name, signal_name.replace(TREE_NAME, "tree_thrown"))
            .replace(f"{args['mc_option']}", "thrown")
        )

        gen_mv_reconstructed_command = mv_reconstructed_command[:-1]
        gen_mv_reconstructed_command.append(
            mv_reconstructed_command[-1].replace(f"mc{args['mc_option']}", "mcthrown")
        )

        gen_phasespace_command = phasespace_command[:-1]
        gen_phasespace_command.append(
            phasespace_command[-1]
            .replace(phasespace_name, phasespace_name.replace(TREE_NAME, "tree_thrown"))
            .replace(phasespace_option, "thrown")
        )

        gen_mv_phasespace_command = mv_phasespace_command[:-1]
        gen_mv_phasespace_command.append(
            mv_phasespace_command[-1]
            .replace("Acc", "Gen")
            .replace(phasespace_option, "")
        )

        # run the generated trees for signal and phasespace MC
        subprocess.call(gen_reconstructed_command)
        subprocess.call(gen_mv_reconstructed_command)
        if args["make_phasespace_trees"]:
            subprocess.call(gen_phasespace_command)
            subprocess.call(gen_mv_phasespace_command)

    pass


def parse_args() -> dict:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-o",
        "--mc_option",
        type=str,
        default="",
        choices=[
            "",
            "accept",
            "matched",
            "noaccidental",
            "accept_noaccidental",
            "matched_noaccidental",
        ],
        help=(
            "Options for the DSelector to use when processing Monte Carlo trees."
            " Default: '', which means treat just like data."
            " accept: use thrown P4s in generated order (no mis-ID, resolution or"
            " combinatorics)."
            " matched: use reconstructed P4s for matches to thrown particles"
            " (no mis-ID or combinatorics, but resolution affects apply)."
            " noaccidental: use only the beam photons which generated the event."
        ),
    )
    parser.add_argument(
        "--make_gen_trees",
        action="store_true",
        help="Option to turn on/off generated trees creation.",
    )
    parser.add_argument(
        "--make_phasespace_trees",
        action="store_true",
        help="Option to turn on/off phasespace trees creation.",
    )
    parser.add_argument(
        "-r",
        "--run_periods",
        nargs="+",
        type=str,
        default=["2017_01", "2018_01", "2018_08"],
        choices=["2017_01", "2018_01", "2018_08"],
        help=("List of GlueX run periods you'd like to analyze."),
    )
    return vars(parser.parse_args())


if __name__ == "__main__":
    args = parse_args()
    main(args)
