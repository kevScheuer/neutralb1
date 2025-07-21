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

    acc_flag = "_noaccidental" if args["remove_accidentals"] else ""

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
                f' "{args["mc_option"]}{acc_flag}"'
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
                f' "{args["mc_option"]}_noaccidental"'  # phsp always uses gen'd photons
                ")"
            ),
        ]

        mv_reconstructed_command = [
            "mv",
            "AmpToolsInputTree.root",
            (
                "AmpToolsInputTree_sum_PARA_0_"
                f"{run_period}{SIGNAL_VERSION}_mc{args['mc_option']}{acc_flag}.root"
            ),
        ]
        # if using an option, an extra underscore is needed
        phsp_option = f"_{args['mc_option']}" if args["mc_option"] else ""
        mv_phasespace_command = [
            "mv",
            "AmpToolsInputTree.root",
            (
                "anglesOmegaPiPhaseSpaceAcc_"
                f"{run_period}{PHASESPACE_VERSION}{phsp_option}.root"
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
            .replace(args["mc_option"], "thrown")
        )

        gen_mv_phasespace_command = mv_phasespace_command[:-1]
        gen_mv_phasespace_command.append(
            mv_phasespace_command[-1].replace("Acc", "Gen").replace(phsp_option, "")
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
        ],
        help=(
            "Options for the DSelector to use when processing Monte Carlo trees."
            " '' (default): apply all effects and treat like real data"
            " accept: apply detector acceptance effects but use thrown P4s."
            " matched: match the reconstructed P4s to the thrown P4s to remove"
            " bad combinatorics effects, but still use reconstructed P4 values"
        ),
    )
    parser.add_argument(
        "--remove_accidentals",
        action="store_true",
        help=(
            "Option to add 'noaccidental' flag to the mc_option, which removes"
            " sideband subtraction of accidental beam photons"
        ),
    )
    parser.add_argument(
        "--make_gen_trees",
        action="store_true",
        help=("Option to turn on/off generated trees creation."),
    )
    parser.add_argument(
        "--make_phasespace_trees",
        action="store_true",
        help=("Option to turn on/off phasespace trees creation."),
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
