"""Edits bin_template.cfg with exact bin info, and links data/phasespace files to cwd

This file is run in the run_fit.sh file, where it receives the information
needed to overwrite the template cfg file with specific bin info
"""

import os
import argparse

# CONSTANTS
POL_DICT = {
    "PARA_0": {"angle": 0, "fraction": 0.3519},
    "PERP_45": {"angle": 45, "fraction": 0.3374},
    "PERP_90": {"angle": 90, "fraction": 0.3303},
    "PARA_135": {"angle": 135, "fraction": 0.3375},
}


def main(
    data_ver: str,
    data_dir: str,
    phasespace_ver: str,
    phasespace_dir: str,
    data_out_dir: str,
    orientations: list,
    run_period: str,
    mass_range: list,
    t_range: list,
):
    E_range = [8.2, 8.8]
    f_template = "bin_template.cfg"

    # print out some info to log
    print(
        f"\ncurrent working dir:\t\t{os.getcwd()}"
        f"\ndata version and dir:\t\t{data_ver}\t{data_dir}",
        f"\nphasespace version and dir:\t{phasespace_ver}\t\t{phasespace_dir}",
    )

    # link files locally to perform fit
    linkfiles(
        data_ver,
        data_dir,
        phasespace_ver,
        phasespace_dir,
        "./",
        orientations,
        run_period,
    )
    # link files in output dir to allow the plotter to run
    linkfiles(
        data_ver,
        data_dir,
        phasespace_ver,
        phasespace_dir,
        data_out_dir,
        orientations,
        run_period,
    )

    # get template of fit config file, overwrite with bin info, and write to new cfg
    f = open(f_template, "r")
    cfg_data = f.read()
    f.close()

    # write file at correct location
    path = os.path.join(os.path.dirname(__file__), "fit.cfg")
    with open(path, "w") as f:
        f.write(cfg_data)

    return


def linkfiles(
    data_ver: str,
    data_dir: str,
    phasespace_ver: str,
    phasespace_dir: str,
    output_dir: str,
    orientations: list,
    run_period: str,
) -> None:
    """Create symbolic links between data and phasespace files to output dir

    Args:
        data_ver (str): data version in file name e.g. '_data', '_verXY.Z_mc'
        data_dir (str): source directory for data files
        phasespace_ver (str): phasespace version in file name
        phasespace_dir (str): source directory for phasespace files
        output_dir (str): where links to data and phasespace files will be created
        orientations (list): diamond orientation settings
        run_period (str): data run period e.g. '2017', 'GlueXI'

    Raises:
        FileExistsError: data file not found at data_dir
        FileExistsError: Generated phasespace file not found at phasespace_dir
        FileExistsError: Accepted phasespace file not found at phasespace_dir
    """

    TREE = "AmpToolsInputTree_sum"
    AMP_STR = "anglesOmegaPiAmplitude"
    PH_STR = "anglesOmegaPiPhaseSpace"

    # LINK DATA
    for ont in orientations:
        data_src = f"{data_dir}/{TREE}_{ont}_{run_period}_{data_ver}.root"
        if not os.path.isfile(data_src):
            raise FileExistsError(f"Path {data_src} does not exist!\n")
        os.system(
            f"ln -sf {data_src} {output_dir}/{AMP_STR}_{POL_DICT[ont]['angle']}.root"
        )

    # LINK PHASESPACE
    # setup the gen angles file
    gen_src = f"{phasespace_dir}/{PH_STR}Gen_{run_period}_{phasespace_ver}.root"
    if not os.path.isfile(gen_src):
        raise FileExistsError(f"Path {gen_src} does not exist!\n")
    os.system(f"ln -sf {gen_src} {output_dir}/{PH_STR}.root")

    # link the accepted file (if thrown, sets Accepted phasespace as Generated)
    ver = "Gen" if "mcthrown" in data_ver else "Acc"
    acc_src = f"{phasespace_dir}/{PH_STR}{ver}_{run_period}_{phasespace_ver}.root"

    if not os.path.isfile(acc_src):
        raise FileExistsError(f"Path {acc_src} does not exist!\n")

    os.system(f"ln -sf {acc_src} {output_dir}/{PH_STR}Acc.root")

    return


def check_positive_float(val) -> float:
    fl = float(val)
    if fl < 0.0:
        raise argparse.ArgumentTypeError(f"{fl} must be >= 0")
    return fl


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Overwrites a template cfg file (produced by write_config_loop)"
            " with parameters for a binned fit to be performed"
        )
    )
    parser.add_argument(
        "--data",
        nargs=2,
        metavar=("TYPE", "DIR"),
        help=(
            "specify data version and directory."
            " Data types typically in form: thrown='verXY.Z_mcthrown,"
            " MC='verXY.Z_mc' and DATA='_data'"
        ),
    )
    parser.add_argument(
        "--phasespace",
        nargs=2,
        metavar=("TYPE", "DIR"),
        help="specify phasespace version and directory",
    )
    parser.add_argument(
        "--data_out_dir", type=str, help="specify where fit results are output to"
    )

    # fit details
    parser.add_argument(
        "-o",
        "--orientations",
        nargs="*",
        choices=["PARA_0", "PARA_135", "PERP_45", "PERP_90", "ALL"],
        help="diamond orientations",
    )
    parser.add_argument("-r", "--run", type=str, help="run period of data/MC to fit")
    parser.add_argument(
        "-m",
        "--mass_range",
        nargs=2,
        type=check_positive_float,
        help="start and end of mass bins to fit to",
    )
    parser.add_argument(
        "-t",
        "--t_range",
        nargs=2,
        type=check_positive_float,
        help="start and end of t bins to fit to",
    )

    args = parser.parse_args()
    data_ver, data_dir = args.data
    phasespace_ver, phasespace_dir = args.phasespace
    data_out_dir = args.data_out_dir
    orientations = args.orientations
    run_period = args.run
    mass_range = args.mass_range
    t_range = args.t_range

    main(
        data_ver,
        data_dir,
        phasespace_ver,
        phasespace_dir,
        data_out_dir,
        orientations,
        run_period,
        mass_range,
        t_range,
    )
