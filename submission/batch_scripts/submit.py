"""Main batch script to set parameters for amplitude analysis fits.

See README for process architecture, or run with "--help" to understand variable usage

Optional improvements to make script more generalized:
    - Ability to choose polar coordinates. Means init_imag -> init_phase
    - Ability or some handling of matching the GPU architecture
        not surefire, but can check GPU arch set in $AMPTOOLS_HOME makefile
    - Arg for the reaction line, right now hardcoded to "Beam Proton Pi01 Pi02 Pi+ Pi-"
    - Ability to pass extra user options like 'omega3pi' (which is currently hardcoded)

SLURM INFO (https://scicomp.jlab.org/scicomp/slurmJob/slurmInfo)
"""

import argparse
import math
import os
import pathlib
import pwd
import subprocess
import time
from typing import List

import write_config

# Constants
USER = pwd.getpwuid(os.getuid())[0]
VOLATILE_DIR = f"/volatile/halld/home/{USER}"
CODE_DIR = str(pathlib.Path(__file__).resolve().parent) + "/"


def main(args: dict) -> None:
    """Main function. Run file with "--help" to understand variable usage"""

    # unpack some config arguments
    energy_min, energy_max = args["energy"]
    n_gpus, gpu_type = int(args["gpu"][0]), args["gpu"][1]

    if "ALL" in args["orientations"]:
        args["orientations"] = ["PARA_0", "PERP_45", "PERP_90", "PARA_135"]

    # error checks
    if args["truth_file"] and (
        args["data_version"] not in args["truth_file"]
        or not os.path.isfile(f"{CODE_DIR}{args['truth_file']}")
    ):
        raise ValueError("MC and truth versions must match and truth file must exist!")
    if f"{args['phasespace_version']}{args['phasespace_option']}" not in "\t".join(
        os.listdir(args["phasespace_dir"])
    ):
        raise FileNotFoundError(
            f"Phasespace {args['phasespace_version']}{args['phasespace_option']} does"
            f" not exist in directory: {args['phasespace_dir']}"
        )
    if f"{args['data_version']}{args['data_option']}" not in "\t".join(
        os.listdir(args["data_dir"])
    ):
        raise FileNotFoundError(
            f"Data {args['data_version']}{args['data_option']} does not exist in"
            f" directory: {args['data_dir']}"
        )
    for ont in args["orientations"]:
        if ont not in "\t".join(os.listdir(args["data_dir"])):
            raise FileNotFoundError(
                f"Orientation {ont} not found in directory: {args["data_dir"]}"
            )
    if len(args["t_momenta"]) < 2 or len(args["masses"]) < 2:
        raise ValueError("Must specify at LEAST two arguments for t and mass bins")
    if not args["waveset"]:
        raise ValueError("Must specify a waveset to fit with")
    if args["phase_reference"]:
        if any(
            phase_ref in args["remove_waves"] for phase_ref in args["phase_reference"]
        ):
            raise ValueError(
                "Phase reference waves cannot be in the list of waves to remove"
            )
        if args["phase_reference"][0][0] == args["phase_reference"][1][0]:
            raise ValueError("Phase references must be in opposite reflectivities")

    # get t and mass bins to fit over
    low_t_edges, high_t_edges = make_bins(args["t_momenta"])
    low_mass_edges, high_mass_edges = make_bins(args["masses"])

    # round the mass bins to 3 decimal places (single MeV precision)
    low_mass_edges = [round(m, 3) for m in low_mass_edges]
    high_mass_edges = [round(m, 3) for m in high_mass_edges]

    # create ROOT data files with cuts if not yet done
    create_data_files(
        low_t_edges,
        high_t_edges,
        energy_min,
        energy_max,
        low_mass_edges,
        high_mass_edges,
        args["orientations"],
        args["run_periods"],
        args["data_dir"],
        args["data_version"],
        args["data_option"],
        args["phasespace_dir"],
        args["phasespace_version"],
        args["phasespace_option"],
        args["cut_recoil_pi_mass"],
        args["reaction"],
    )

    # Create config file template "fit.cfg" that works for any bin
    if not args["truth_file"]:
        write_config.main(args)

    # when set True in the loop, will always skip asking the user if they want to
    # overwrite files
    skip_input = False

    # These loops create a job submission for every combination of
    # run period, and mass & t bins
    for run_period in args["run_periods"]:
        for low_t, high_t in zip(low_t_edges, high_t_edges):
            for low_mass, high_mass in zip(low_mass_edges, high_mass_edges):
                truth_subdir = "truth/" if args["truth_file"] else ""

                # PREPARE DIRECTORIES
                running_dir = "/".join(
                    (
                        VOLATILE_DIR,
                        "ampToolsFits",
                        args["reaction"],
                        run_period,
                        "-".join(sorted(args["orientations"])),
                        f"{args['data_version']}{args['data_option']}",
                        f"{args['phasespace_version']}{args['phasespace_option']}",
                        "_".join(sorted(args["waveset"])),
                        f"recoil-pi-mass_{args['cut_recoil_pi_mass']}",
                        f"t_{low_t:.2f}-{high_t:.2f}",
                        f"mass_{low_mass:.3f}-{high_mass:.3f}",
                        truth_subdir,
                    )
                )
                pathlib.Path(running_dir).mkdir(parents=True, exist_ok=True)

                log_dir = running_dir + "log/"
                pathlib.Path(log_dir).mkdir(parents=True, exist_ok=True)

                # truth fits don't require rand or bootstrap fits
                if not args["truth_file"]:
                    rand_dir = running_dir + "rand/"
                    bootstrap_dir = running_dir + "bootstrap/"
                    pathlib.Path(rand_dir).mkdir(parents=True, exist_ok=True)
                    pathlib.Path(bootstrap_dir).mkdir(parents=True, exist_ok=True)

                # location of pre-selected data file
                source_file_dir = volatile_path(
                    args["reaction"],
                    args["cut_recoil_pi_mass"],
                    low_t,
                    high_t,
                    energy_min,
                    energy_max,
                    low_mass,
                    high_mass,
                )

                # if a completed fit is found in the directory, ask if the user
                # is sure they want to overwrite it
                if not skip_input:
                    if os.path.isfile(f"{running_dir}best.fit") or os.path.isfile(
                        f"{running_dir}best_truth.fit"
                    ):
                        print(
                            f"Completed fit(s) already exist in {running_dir}, are"
                            " you sure you want to submit this job and overwrite the"
                            " file? (yes/no/skip_input/exit)"
                        )
                        while True:
                            ans = str(input())
                            if ans == "yes" or ans == "y" or ans == "no" or ans == "n":
                                break
                            elif ans == "skip_input":
                                skip_input = True
                                break
                            elif ans == "exit" or ans == "exit()":
                                exit()
                            else:
                                print(
                                    "Please answer yes, no, skip_input (to submit all"
                                    " jobs without asking), or exit"
                                )
                        if ans == "no" or ans == "n":
                            continue

                # copy in needed files
                if args["truth_file"]:
                    os.system(f"cp -f {CODE_DIR}{args['truth_file']} {running_dir}")
                else:
                    os.system(f"cp -f {CODE_DIR}fit.cfg {running_dir}")

                # prepare job name to be shown on slurm webpage
                job_name = "_".join(
                    (
                        args["reaction"],
                        run_period,
                        "-".join(sorted(args["orientations"])),
                        f"{args['data_version']}{args['data_option']}",
                        f"{args['phasespace_version']}{args['phasespace_option']}",
                        "_".join(sorted(args["waveset"])),
                        f"recoil-pi-mass_{args['cut_recoil_pi_mass']}",
                        f"t_{low_t:.2f}-{high_t:.2f}",
                        f"mass_{low_mass:.3f}-{high_mass:.3f}",
                    )
                )

                # handoff arguments to bash script that will run on the ifarm
                script_command = " ".join(
                    (
                        f"{CODE_DIR}run_fit.sh",
                        f"-o {'-'.join(sorted(args['orientations']))}",
                        f"-r {run_period}",
                        f"-n {args['nrand']}",
                        f"-d {args['data_version']}",
                        f"-p {args['phasespace_version']}",
                        f"-s {source_file_dir}",
                        f"-C {CODE_DIR}",
                        f"-R {args['reaction']}",
                        f"-b {args['bootstrap']}",
                    )
                )
                if args["truth_file"]:
                    script_command += f" -t {args['truth_file']}"
                if args["data_option"]:
                    script_command += f" -D {args['data_option']}"
                if args["phasespace_option"]:
                    script_command += f" -P {args['phasespace_option']}"
                submit_slurm_job(
                    job_name,
                    script_command,
                    running_dir,
                    log_dir,
                    gpu_type,
                    n_gpus,
                    args["email"],
                    args["email_type"],
                    args["time_limit"],
                )

    return


def create_data_files(
    low_t_edges: list,
    high_t_edges: list,
    energy_min: float,
    energy_max: float,
    low_mass_edges: list,
    high_mass_edges: list,
    orientations: list,
    run_periods: list,
    data_dir: str,
    data_ver: str,
    data_option: str,
    phasespace_dir: str,
    phasespace_ver: str,
    phasespace_option: str,
    cut_recoil_pi_mass: float,
    reaction: str,
) -> None:
    """Create data files with cuts and store them in volatile

    The typical ROOTDataReader method in the .cfg files reads in data much too slowly,
    and is repetitive when the same TEM region is being selected. This function will
    copy the source data files to the volatile directory, and cut them to the desired
    TEM region for quick access. Each TEM + recoil_pi_mass bin will have its own ifarm
    job. If jobs are submitted, then the program exits to avoid pwa jobs being submitted
    without the necessary data files.

    Args:
        low_t_edges (list): values of low t bin edges
        high_t_edges (list): values of high t bin edges
        energy_min (float): minimum beam energy value
        energy_max (float): maximum beam energy value
        low_mass_edges (list): values of low omega pi mass bin edges
        high_mass_edges (list): values of high omega pi mass bin edges
        orientations (list): diamond orientation settings
        run_period (list): string in file name, e.g. 2017_01, allPeriods
        data_dir (str): directory of original data files
        data_ver (str): MC version string in file name, e.g. ver03.1
        data_option (str): MC option string in file name, e.g. _mcthrown
        phasespace_dir (str): directory of original phasespace files
        phasespace_ver (str): MC version string in file name, e.g. ver03
        phasespace_option (str): MC option string in file name,
            e.g. _accept_noaccidental
        cut_recoil_pi_mass (float): removes events below the given recoil-pion mass
    Raises:
        FileExistsError: Generated phasespace file not found
        FileExistsError: Accepted phasespace file not found
        FileExistsError: Data File not found
    """

    # track what files need to be copied to what directories on volatile
    src_files_to_copy_to_dir = {}

    # if user wants to skip being asked for every job submission, this will become True
    skip_input = False

    jobs_submitted = False  # becomes True if any jobs are submitted

    for run_period in run_periods:
        # find generated phasespace file (this is always thrown, so no
        # phasespace_option is used here)
        gen_file = (
            f"{phasespace_dir}/anglesOmegaPiPhaseSpaceGen_"
            f"{run_period}_{phasespace_ver}.root"
        )
        if not os.path.isfile(gen_file):
            raise FileExistsError(f"Path {gen_file} does not exist!\n")
        src_files_to_copy_to_dir[gen_file] = []

        # find accepted phasespace file (if needed)
        if "mcthrown" not in data_option:
            acc_file = (
                f"{phasespace_dir}/anglesOmegaPiPhaseSpaceAcc_"
                f"{run_period}_{phasespace_ver}{phasespace_option}.root"
            )
            if not os.path.isfile(acc_file):
                raise FileExistsError(f"Path {acc_file} does not exist!\n")
            src_files_to_copy_to_dir[acc_file] = []

        # find data files
        data_files = []
        for ont in orientations:
            f = (
                f"{data_dir}/AmpToolsInputTree_sum_{ont}_{run_period}"
                f"_{data_ver}{data_option}.root"
            )
            if not os.path.isfile(f):
                raise FileExistsError(f"Path {f} does not exist!\n")
            data_files.append(f)
            src_files_to_copy_to_dir[f] = []

        # loop over TEM bins to determine what bins need the cut data files
        for low_mass, high_mass in zip(low_mass_edges, high_mass_edges):
            for low_t, high_t in zip(low_t_edges, high_t_edges):
                # create directory for each TEM bin if not already done
                bin_dir = volatile_path(
                    reaction,
                    cut_recoil_pi_mass,
                    low_t,
                    high_t,
                    energy_min,
                    energy_max,
                    low_mass,
                    high_mass,
                )
                pathlib.Path(bin_dir).mkdir(parents=True, exist_ok=True)

                # if files not already in volatile, add to list to be copied
                if not os.path.isfile(f"{bin_dir}/{gen_file.split('/')[-1]}"):
                    src_files_to_copy_to_dir[gen_file].append(bin_dir)
                if "mcthrown" not in data_option and not os.path.isfile(
                    f"{bin_dir}/{acc_file.split('/')[-1]}"
                ):
                    src_files_to_copy_to_dir[acc_file].append(bin_dir)
                for data_file in data_files:
                    if not os.path.isfile(f"{bin_dir}/{data_file.split('/')[-1]}"):
                        src_files_to_copy_to_dir[data_file].append(bin_dir)

        # if any files need to be cut and copied for this run period, ask user if
        # they want to submit jobs to create them
        if bool([a for a in src_files_to_copy_to_dir.values() if a != []]):
            if not skip_input:
                print(
                    f"\nFiles need to be cut and copied to the following directories"
                    f" for run period {run_period}:"
                )
                for src, dirs in src_files_to_copy_to_dir.items():
                    if dirs:
                        print(f"\n{src} ->")
                        for dir in dirs:
                            print(f"  {dir}")
                print(
                    "Do you want to submit jobs to create these files?"
                    " (yes/no/skip_input/exit)"
                )
                while True:
                    ans = str(input())
                    if ans == "yes" or ans == "y" or ans == "no" or ans == "n":
                        break
                    elif ans == "skip_input":
                        skip_input = True
                        break
                    elif ans == "exit" or ans == "exit()":
                        exit()
                    else:
                        print(
                            "Please answer yes, no, skip_input (to submit all jobs"
                            " without asking), or exit"
                        )
                if ans == "no" or ans == "n":
                    continue

            # submit jobs to create the files
            jobs_submitted = True
            for src_file, dirs in src_files_to_copy_to_dir.items():
                for dir in dirs:
                    # create log dir
                    log_dir = dir + "/log/"
                    pathlib.Path(log_dir).mkdir(parents=True, exist_ok=True)

                    # extract the t and mass bin values from the directory path
                    low_t = float(dir.split("t_")[1].split("-")[0])
                    high_t = float(dir.split("t_")[1].split("-")[1].split("/")[0])
                    low_mass = float(dir.split("/mass_")[1].split("-")[0])
                    high_mass = float(
                        dir.split("/mass_")[1].split("-")[1].split("/")[0]
                    )

                    # create command to run the ROOT macro
                    command = (
                        f"source {CODE_DIR}setup_gluex.sh && root -l -b -q"
                        f" '{CODE_DIR}copy_tree_with_cuts.C("
                        f'"{src_file}", "{dir}",'
                        f' "{cut_recoil_pi_mass}",'
                        f' "{low_t:.2f}", "{high_t:.2f}",'
                        f' "{energy_min}", "{energy_max}",'
                        f' "{low_mass:.3f}", "{high_mass:.3f}"'
                        ")'"
                    )
                    submit_slurm_job(
                        job_name=f"copy_{src_file.split('/')[-1]}_to_{dir}",
                        script_command=command,
                        running_dir=dir,
                        log_dir=log_dir,
                        gpu_type="",
                        n_gpus=0,
                        email_address="",
                        email_type=[],
                        time_limit="00:30:00",
                        n_cpus=8,
                    )

    if jobs_submitted:
        print(
            "Jobs have been submitted to create the necessary data files."
            " Please wait for them to finish before submitting PWA fits."
            " Job progress can be monitored at"
            " https://scicomp.jlab.org/scicomp/slurmJob/activeJob, or by running"
            " 'squeue -u $USER' in a terminal"
        )
        exit()

    return


def submit_slurm_job(
    job_name: str,
    script_command: str,
    running_dir: str,
    log_dir: str,
    gpu_type: str,
    n_gpus: int,
    email_address: str,
    email_type: List[str],
    time_limit: str,
    mem_per_cpu: str = "5000M",
    n_cpus: int = 32,
) -> None:
    """Submit a slurm job to the ifarm using an mpi+gpu build

    Args:
        job_name (str): shown on the scicomp webpage
        script_command (str): bash script with its arguments
        running_dir (str): /volatile/ location
        log_dir (str): where slurm log files are stored
        gpu_type (str): card type to be used
        n_gpus (int): how many gpu cards to use (supported by mpi)
        email_address (str): send email to passed address
        email_type (str): when to send email (BEGIN, END, FAIL)
        time_limit (str, optional): Max wall-time in Hour:Min:Sec. Defaults to "1:00:00"
        mem_per_cpu (str, optional): Default of 5GB appear to be min needed for fit
            jobs, though small jobs like phasespace generation can use less
        n_cpus (int, optional): Number of mpi cpus to use (only used if n_gpus=0).
            Defaults to 32.
    """

    with open("tempSlurm.txt", "w") as slurm_out:
        slurm_out.write(
            "#!/bin/sh \n"
            "#SBATCH -A halld\n"
            f"#SBATCH --time={time_limit} \n"
            f"#SBATCH --chdir={running_dir}\n"
            f"#SBATCH --error={log_dir}log.err \n"
            f"#SBATCH --output={log_dir}log.out \n"
            f"#SBATCH --job-name={job_name} \n"
            f"#SBATCH --mem-per-cpu={mem_per_cpu} \n"
            "#SBATCH --cpus-per-task=1 \n"
            "#SBATCH --ntasks-per-core=1 \n"
            "#SBATCH --threads-per-core=1 \n"
            "#SBATCH --constraint=el9 \n"
        )
        if email_address:
            mail_type = ",".join(email_type)
            slurm_out.write(
                f"#SBATCH --mail-user={email_address} \n"
                f"#SBATCH --mail-type={mail_type} \n"
            )
        # different requirements for GPU and CPU fits
        if n_gpus > 0:
            slurm_out.write(
                "#SBATCH --partition=gpu \n"
                f"#SBATCH --gres=gpu:{gpu_type}:{n_gpus} \n"
                f"#SBATCH --ntasks={n_gpus+1} \n"  # mpigpu always needs n_gpus+1
            )
        else:
            slurm_out.write(f"#SBATCH --partition=ifarm \n#SBATCH --ntasks={n_cpus} \n")
        slurm_out.write(script_command)

    # wait half a second to avoid job skip error if too many submitted quickly
    time.sleep(0.5)
    subprocess.call(["sbatch", "tempSlurm.txt"])

    # remove temporary submission file
    os.remove("tempSlurm.txt")
    return


def make_bins(args: List[float]) -> tuple[List[float], List[float]]:
    """Makes low and high bin edges given a list of values. See cases below

    Case 1: Linearly spaced bins from min to max with a specified bin width.
        List is of form [min, max, width].
    Case 2: Custom binning as defined by user. List is of form
        [min_bin, bin_2, ... bin_k, ... max_bin]

    Args:
        args (list): bin values according to either case

    Raises:
        ValueError: custom bins are not ordered sequentially
        RuntimeError: avoids user error where cases can overlap

    Returns:
        tuple[List[float], List[float]]: Paired lists defining the low and high bin
            edges
    """
    low_edges = []
    high_edges = []
    delta = 1e-15
    diff = args[-2] - args[0]
    # case 1: list = [min, max, width]. The width MUST always be smaller than the
    #   difference between the max and min bins. The "or" statement handles floating
    #   point precision if its one bin
    if args[-1] < diff or abs(args[-1] - diff) < delta:
        if len(args) != 3:
            raise (
                RuntimeError(
                    "User gave a custom list, but the last bin value is too"
                    " close to the 2nd to last value"
                )
            )
        min, max, width = args
        n_bins = math.ceil((max - min) / width)
        if n_bins == 0:
            n_bins = 1
        for i in range(n_bins):
            low_edges.append(min + width * i)
            high_edges.append(min + width * (i + 1))
    # case 2: list =[bin_1, bin_i..., bin_n. ]. Because the last value is ALWAYS
    # greater than the 2nd to last, this will never accidentally fall into case 2
    else:
        if args != sorted(args):
            raise ValueError(
                "User gave custom list, but it was not ordered from"
                " smallest to largest bin"
            )
        for i in range(len(args) - 1):
            low_edges.append(args[i])
            high_edges.append(args[i + 1])

    return low_edges, high_edges


def check_positive_float(val) -> float:
    # custom error check for argparse to ensure float is positive
    fl = float(val)
    if fl < 0.0:
        raise argparse.ArgumentTypeError(f"{fl} must be >= 0")
    return fl


def volatile_path(
    reaction: str,
    cut_recoil_pi_mass: float,
    low_t: float,
    high_t: float,
    energy_min: float,
    energy_max: float,
    low_mass: float,
    high_mass: float,
) -> str:
    """Create consistent volatile path for pre-selected data files

    The string returned here is used in both the data file creation function and the
    main function. This function ensures that the path is consistent between the two.
    """
    return "/".join(
        (
            VOLATILE_DIR,
            "TMPDIR",
            "ampToolsFits",
            reaction,
            "data_files",
            f"recoil-pi-mass_{cut_recoil_pi_mass}",
            f"t_{low_t:.2f}-{high_t:.2f}",
            f"E_{energy_min:.2f}-{energy_max:.2f}",
            f"mass_{low_mass:.3f}-{high_mass:.3f}",
        )
    )


def parse_args() -> dict:
    parser = argparse.ArgumentParser(description="Submit PWA fits various configs")

    # waveset and any modifications
    parser.add_argument(
        "-w",
        "--waveset",
        nargs="+",
        choices=[
            "0m",
            "1p",
            "1m",
            "2p",
            "2m",
            "3m",
            "iso",
            "b1",
            "b1free",
            "rho",
            "rhofree",
            "nodalitz",
        ],
        help="Waveset to fit with",
    )
    parser.add_argument(
        "--phase-reference",
        type=str,
        nargs=2,
        metavar="pJPmL, mJPmL",
        default="",
        help=(
            "Flag the waves (in 'eJPmL' format) for each reflectivity whose phase will"
            " be constrained to 0. Empty (default) picks first JP, m, L combination"
        ),
    )
    parser.add_argument(
        "--phaselock",
        action="store_true",
        help=(
            "Option to turn on the 'phaselock' model, where phases are common across"
            " m-projections for a particular eJPL combination"
        ),
    )
    parser.add_argument(
        "-ds",
        "--ds-ratio",
        type=str,
        default="",
        choices=["free", "fixed", "split"],
        help=(
            "option to modify the ratio & phase between the D/S waves."
            " Leaving empty (default) lets them float within some bounds."
            " 'Free' removes the parameters, allowing them to float freely."
            " 'Fixed' sets to E852 nominal values"
            " (ratio=0.27 & phase=0.184 radians)."
            " 'Split' gives each reflectivity its own ratio & phase"
        ),
    )
    parser.add_argument(
        "--frame",
        type=str,
        default="",
        choices=["GJ", "Adair"],
        metavar="decay frame option",
        help="change decay frame used, empty default means helicity will be used",
    )
    parser.add_argument(
        "--force-refl",
        type=int,
        default=0,
        choices=[-1, 1],
        help="only allow a single reflectivity by setting to +1 or -1",
    )
    parser.add_argument(
        "--init-refl",
        type=int,
        default=0,
        choices=[-1, 1],
        help="initialize single reflectivity by setting to +1 or -1",
    )
    parser.add_argument(
        "--init-real",
        type=float,
        default=100.0,
        metavar=("real_part"),
        help="value to initialize real cartesian part of amplitude to",
    )
    parser.add_argument(
        "--init-imag",
        type=float,
        default=100.0,
        metavar=("imag_part"),
        help="value to initialize imaginary cartesian part of amplitude to",
    )
    parser.add_argument(
        "--remove-waves",
        nargs="+",
        default=[],
        help="waves given in 'eJPmL' format to remove from the waveset",
    )

    # details of fit type
    parser.add_argument(
        "--reaction",
        type=str,
        default="omegapi",
        help="base reaction name to be used in cfg files",
    )
    parser.add_argument(
        "-n", "--nrand", type=int, default=20, help="number of random fits"
    )
    parser.add_argument(
        "-o",
        "--orientations",
        nargs="*",
        default=["PARA_0"],
        choices=["PARA_0", "PARA_135", "PERP_45", "PERP_90", "ALL"],
        help="diamond orientations",
    )
    parser.add_argument(
        "-r",
        "--run-periods",
        nargs="*",
        default=["allPeriods"],
        choices=["allPeriods", "2017_01", "2018_01", "2018_08"],
        help="run periods of data/MC to fit",
    )
    parser.add_argument(
        "-m",
        "--masses",
        nargs="+",
        type=check_positive_float,
        default=[1.0, 1.5, 0.05],
        help=(
            "start, end, and width of mass bins to fit to. Alternatively the user may"
            " pass their own list of mass bins, ex: -m 1.0 1.1 1.5"
        ),
    )
    parser.add_argument(
        "-t",
        "--t-momenta",
        nargs="+",
        type=check_positive_float,
        default=[0.1, 0.5, 0.1],
        help=(
            "start, end, and width of t bins to fit to. Alternatively the user may"
            " pass their own list of t bins, ex: -t 0.2 0.3 0.6"
        ),
    )
    parser.add_argument(
        "-e",
        "--energy",
        nargs=2,
        type=check_positive_float,
        default=[8.2, 8.8],
        metavar=("low_E_edge", "high_E_edge"),
        help="range of beam energy to fit to. Default is coherent peak region",
    )
    parser.add_argument(
        "--truth_file",
        type=str,
        default="",
        help=(
            "cfg file used to generate signal MC, with fixed parameters. Used for"
            " input-output testing"
        ),
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="/w/halld-scshelf2101/kscheuer/neutralb1/submission/source_files/data",
        help="directory where data files are stored",
    )
    parser.add_argument(
        "--data-version",
        type=str,
        default="data",
        help=(
            "data version to use for fits. GlueX data is typically 'data', and Monte"
            " uses the actual version number i.e. 'verXY.Z'"
        ),
    )
    parser.add_argument(
        "--data-option",
        type=str,
        default="",
        help=(
            "Monte Carlo option used in DSelector. Default assumes real data is used."
            " Options are typically '_mcthrown', '_mc', etc."
        ),
    )
    parser.add_argument(
        "--phasespace_dir",
        type=str,
        default=(
            "/w/halld-scshelf2101/kscheuer/neutralb1/submission/source_files/phasespace"
        ),
        help="directory where phasespace files are stored",
    )
    parser.add_argument(
        "--phasespace-version",
        type=str,
        default="ver03",
        help="phasespace version to use for fits",
    )
    parser.add_argument(
        "--phasespace-option",
        type=str,
        default="",
        help=(
            "Monte Carlo option used in DSelector. Default assumes no special options"
            " were used. An option like '_accept_noaccidental' could be used."
        ),
    )
    parser.add_argument(
        "-c",
        "--cut-recoil-pi-mass",
        type=float,
        default=1.4,
        help=(
            "Cuts events below the given value in the recoil-pion mass spectrum, i.e."
            " a value of 1.3 selects events like 'MRecoilPi > 1.3' when the trees are"
            " copied. This flag assumes the ROOT data files has an 'MRecoilPi' leaf"
        ),
    )
    parser.add_argument(
        "-b",
        "--bootstrap",
        type=int,
        default=0,
        help=(
            "When non-zero value passed, the requested number of bootstrap fits will"
            " be performed starting from the nominal values found at the end of all the"
            " --nrand fits in each bin"
        ),
    )

    # other arguments
    parser.add_argument(
        "--template-name",
        type=str,
        default="template.cfg",
        help="template with default values to copy and overwrite",
    )
    parser.add_argument(
        "-g",
        "--gpu",
        nargs=2,
        default=["0", ""],
        choices=["1", "2", "3", "4", "T4", "TitanRTX", "A100", "A800"],
        metavar=("#GPUs", "CARD"),
        help="set # of GPUs to use for a card. Default assumes only CPU fits",
    )
    parser.add_argument(
        "--email",
        type=str,
        default="",
        help=("when email address given, mails address when a job starts/stops/fails."),
    )
    parser.add_argument(
        "--email-type",
        type=str,
        default=["BEGIN", "END", "FAIL"],
        nargs="+",
        choices=["BEGIN", "END", "FAIL"],
        help=(
            "If email flag is used, this argument handles the cases when an email is"
            " sent. Default is 'BEGIN,END,FAIL'"
        ),
    )
    parser.add_argument(
        "--time-limit",
        type=str,
        default="01:00:00",
        help=("Max walltime for each slurm job. Default assumes quick jobs (1 hr)"),
    )

    args = parser.parse_args()

    return vars(args)


if __name__ == "__main__":
    args = parse_args()
    main(args)
