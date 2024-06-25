""" Main batch script to set parameters for amplitude analysis fits.

See README for process architecture
	
TODO:        
    mails broken (says 'requested node config is not available')
	Add ability to choose polar coordinates. Means init_imag -> init_phase	
	add ability or some handling of matching the GPU architecture
        not surefire, but can check GPU arch set in $AMPTOOLS_HOME makefile
    Add arg for reaction line, right now hardcoded to "Beam Proton Pi01 Pi02 Pi+ Pi-"
    Add ability to pass extra user options like 'omega3pi' (which is hardcoded rn)

SLURM INFO (https://scicomp.jlab.org/scicomp/slurmJob/slurmInfo)
"""

import argparse
import os
import pathlib
import pwd
import subprocess
import time

import write_config

# Constants
USER = pwd.getpwuid(os.getuid())[0]
VOLATILE_DIR = f"/volatile/halld/home/{USER}"
CODE_DIR = f"/w/halld-scshelf2101/{USER}/neutralb1/submission/batch_scripts/"


def main(
    waveset: list,
    phase_reference: str,
    is_phaselock: bool,
    ds_option: str,
    frame: str,
    force_refl: int,
    init_refl: int,
    init_real: float,
    init_imag: float,
    reaction: str,
    num_rand_fits: int,
    orientations: list,
    run_periods: list,
    t_args: list,
    energy_min: float,
    energy_max: float,
    m_args: list,
    data_version: str,
    data_dir: str,
    phasespace_version: str,
    phasespace_dir: str,
    cut_recoil_pi_mass: float,
    template_name: str,
    truth_file: str,
    n_gpus: int,
    gpu_type: str,
    is_send_mail: bool,
    time_limit: str,
):
    """Main function. Run file with "--help" to understand variable usage"""

    if "ALL" in orientations:
        orientations = ["PARA_0", "PERP_45", "PERP_90", "PARA_135"]

    # error checks
    if truth_file and data_version not in truth_file:
        raise ValueError("MC and truth versions must match!")
    if phasespace_version not in "\t".join(os.listdir(phasespace_dir)):
        raise FileNotFoundError(
            f"Phasespace version {phasespace_version} does"
            f" not exist in directory: {phasespace_dir}"
        )
    if data_version not in "\t".join(os.listdir(data_dir)):
        raise FileNotFoundError(
            f"Data version {data_version} does not exist in directory: {data_dir}"
        )
    for ont in orientations:
        if ont not in "\t".join(os.listdir(data_dir)):
            raise FileNotFoundError(
                f"Orientation {ont} not found in directory: {data_dir}"
            )
    if len(t_args) < 2 or len(m_args) < 2:
        raise ValueError("Must specify at LEAST two arguments for t and mass bins")

    # get t and mass bins to fit over
    low_t_edges, high_t_edges = make_bins(t_args)
    low_mass_edges, high_mass_edges = make_bins(m_args)

    # create ROOT data files with cuts if not yet done
    create_data_files(
        low_t_edges,
        high_t_edges,
        energy_min,
        energy_max,
        low_mass_edges,
        high_mass_edges,
        orientations,
        run_periods,
        data_version,
        data_dir,
        phasespace_version,
        phasespace_dir,
        cut_recoil_pi_mass,
    )

    # make strings for directory creation
    waveset_str = "_".join(sorted(waveset))
    waveset = " ".join(waveset)  # spaced version, so it can be used in below command
    orientations_str = "-".join(sorted(orientations))

    # Create config file template "fit.cfg" that works for any bin
    write_config.main(
        waveset,
        phase_reference,
        is_phaselock,
        ds_option,
        frame,
        force_refl,
        init_refl,
        init_real,
        init_imag,
        reaction,
        template_name,
        truth_file,
        orientations,
    )

    # when set True in the loop, will always skip asking the user if they want to
    # overwrite files
    is_skip_all = False

    # These loops create a job submission for every combination of
    # run period, and mass & t bins
    for run_period in run_periods:
        for low_t, high_t in zip(low_t_edges, high_t_edges):
            for low_mass, high_mass in zip(low_mass_edges, high_mass_edges):
                # prepare dirs and create those that don't exist
                running_dir = "/".join(
                    (
                        VOLATILE_DIR,
                        "TMPDIR",
                        "ampToolsFits",
                        reaction,
                        run_period,
                        orientations_str,
                        data_version,
                        phasespace_version,
                        waveset_str,
                        f"recoil-pi-mass_{cut_recoil_pi_mass}",
                        f"t_{low_t:.2f}-{high_t:.2f}",
                        f"mass_{low_mass:.3f}-{high_mass:.3f}/",
                    )
                )

                data_out_dir = running_dir.replace("TMPDIR/", "")
                log_dir = running_dir + "log/"

                source_file_dir = "/".join(
                    (
                        f"{VOLATILE_DIR}",
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

                pathlib.Path(running_dir).mkdir(parents=True, exist_ok=True)
                pathlib.Path(data_out_dir).mkdir(parents=True, exist_ok=True)
                pathlib.Path(log_dir).mkdir(parents=True, exist_ok=True)

                # if a completed fit is found in the output directory, ask if the user
                # is sure they want to overwrite it
                if not is_skip_all:
                    if os.path.isfile(f"{data_out_dir}best.fit"):
                        print(
                            f"best.fit already exists at {data_out_dir}, are"
                            " you sure you want to submit this job and overwrite the"
                            " file? (Answer 'skip_all' to not show this prompt again)"
                        )
                        while True:
                            ans = str(input())
                            if ans == "yes" or ans == "y" or ans == "no" or ans == "n":
                                break
                            elif ans == "skip_all":
                                is_skip_all = True
                                break
                            elif ans == "exit" or ans == "exit()":
                                exit()
                            else:
                                print("Please answer yes, no, or skip_all")
                        if ans == "no" or ans == "n":
                            continue

                # copy in needed files
                os.system(f"cp -f {CODE_DIR}fit.cfg {running_dir}")
                if truth_file:
                    os.system(f"cp -f {CODE_DIR}{truth_file} {running_dir}")

                # prepare names for batch submission
                job_name = "_".join(
                    (
                        reaction,
                        run_period,
                        orientations_str,
                        data_version,
                        phasespace_version,
                        waveset_str,
                        f"recoil-pi-mass_{cut_recoil_pi_mass}",
                        f"t_{low_t:.2f}-{high_t:.2f}",
                        f"mass_{low_mass:.3f}-{high_mass:.3f}",
                    )
                )

                # prepare bash script to be run on slurm node
                script_command = " ".join(
                    (
                        f"{CODE_DIR}run_fit.sh",
                        f"-o {orientations_str}",
                        f" -r {run_period}",
                        f" -n {num_rand_fits}",
                        f" -d {data_version}",
                        f" -p {phasespace_version}",
                        f" -s {source_file_dir}",
                        f" -D {data_out_dir}",
                        f" -C {CODE_DIR}",
                        f" -R {reaction}",
                    )
                )
                submit_slurm_job(
                    job_name,
                    script_command,
                    running_dir,
                    log_dir,
                    gpu_type,
                    n_gpus,
                    is_send_mail,
                    time_limit,
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
    data_ver: str,
    data_dir: str,
    phasespace_ver: str,
    phasespace_dir: str,
    cut_recoil_pi_mass: float,
) -> None:
    """Create data files with TEM region pre-selected

    The typical ROOTDATAREADER method in the .cfg files reads in data much too slowly,
    and is repetitive when the same TEM region is being selected. This function creates
    a dir in volatile under the 'reaction' arg, and dirs for each TEM selection. The
    same dirs are given to run_fits.sh for where to source the data files from.

    Args:
        low_t_edges (list): values of low t bin edges
        high_t_edges (list): values of high t bin edges
        energy_min (float): minimum beam energy value
        energy_max (float): maximum beam energy value
        low_mass_edges (list): values of low omega pi mass bin edges
        high_mass_edges (list): values of high omega pi mass bin edges
        orientations (list): diamond orientation settings
        run_period (list): string in file name, e.g. 2017_01, allPeriods
        data_ver (str): version string in file name, e.g. data, mc_thrown
        data_dir (str): directory of original data files
        phasespace_ver (str): MC version string in file name, e.g. ver03
        phasespace_dir (str): directory of original phasespace files
        cut_recoil_pi_mass (float): removes events below the given recoil-pion mass
    Raises:
        EnvironmentError: ROOT isn't loaded, and so scripts can't run
        FileExistsError: Generated phasespace file not found
        FileExistsError: Accepted phasespace file not found
        FileExistsError: Data File not found
    """

    # first check if ROOT is loaded into the shell environment
    if not os.environ["ROOTSYS"]:
        raise EnvironmentError("ROOTSYS path is not loaded\n")

    for run_period in run_periods:
        # find generated phasespace file
        gen_file = f"anglesOmegaPiPhaseSpaceGen_{run_period}_{phasespace_ver}.root"
        gen_src = f"{phasespace_dir}/{gen_file}"
        if not os.path.isfile(gen_src):
            raise FileExistsError(f"Path {gen_src} does not exist!\n")

        # find accepted phasespace file
        acc_file = gen_file.replace("Gen", "Acc")
        acc_src = f"{phasespace_dir}/{acc_file}"
        if not os.path.isfile(acc_src):
            raise FileExistsError(f"Path {acc_src} does not exist!\n")

        # find data files
        data_files = []
        data_srcs = []
        for ont in orientations:
            f = f"AmpToolsInputTree_sum_{ont}_{run_period}_{data_ver}.root"
            if not os.path.isfile(f"{data_dir}/{f}"):
                raise FileExistsError(f"Path {data_dir}/{f} does not exist!\n")
            data_files.append(f)
            data_srcs.append(f"{data_dir}/{f}")

        # loop over TEM bins
        for low_mass, high_mass in zip(low_mass_edges, high_mass_edges):
            for low_t, high_t in zip(low_t_edges, high_t_edges):
                # create directory for each TEM bin
                dir = "/".join(
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
                pathlib.Path(dir).mkdir(parents=True, exist_ok=True)

                # this root macro will make the given TEM cuts to the 'kin' Tree in the
                # given ROOT file (arg 1) and save the output with the same file name
                # to a specified directory (arg 2)
                command = (
                    "root -l -b -q"
                    f" '{CODE_DIR}copy_tree_with_cuts.C("
                    f'"{gen_src}", "{dir}",'
                    f' "{cut_recoil_pi_mass}",'
                    f' "{low_t}", "{high_t}",'
                    f' "{energy_min}", "{energy_max}",'
                    f' "{low_mass}", "{high_mass}"'
                    ")'"
                )

                # create PhasespaceGen file if it doesn't exist
                if not os.path.isfile(f"{dir}/{gen_file}"):
                    os.system(command)
                # TODO: check if this will accidentally skip making a Acc file if thrown
                #   is selected
                if "mcthrown" not in data_ver and not os.path.isfile(
                    f"{dir}/{acc_file}"
                ):
                    os.system(command.replace(gen_src, acc_src))
                for ont, data_src, data_file in zip(
                    orientations, data_srcs, data_files
                ):
                    if not os.path.isfile(f"{dir}/{data_file}"):
                        os.system(command.replace(gen_src, data_src))

    return


def submit_slurm_job(
    job_name: str,
    script_command: str,
    running_dir: str,
    log_dir: str,
    gpu_type: str,
    n_gpus: int,
    is_send_mail: bool,
    time_limit: str,
    mem_per_cpu: str = "5000M",
    n_cpus: int = 32,
) -> None:
    """Submit a slurm job to the ifarm using an mpi+gpu build

    Args:
        job_name (str): shown on the scicomp webpage
        script_command (str): bash script with its arguments
        running_dir (str): /volatile/TMPDIR location
        log_dir (str): where slurm log files are stored
        gpu_type (str): card type to be used
        n_gpus (int): how many gpu cards to use (supported by mpi)
        is_send_mail (bool): send user mail when job submits/fails/succeeds
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
            "#SBATCH -w sciml2402\n"  # temp: issue with sciml2401 node
            "#SBATCH --constraint=el9 \n"
        )
        if is_send_mail:
            slurm_out.write(
                f"#SBATCH --mail-user={USER}@jlab.org \n"
                "#SBATCH --mail-type=begin \n"
                "#SBATCH --mail-type=end \n "
                "#SBATCH --mail-type=fail \n "
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

    time.sleep(0.5)  # jobs can be skipped if too many submitted quickly
    subprocess.call(["sbatch", "tempSlurm.txt"])

    # remove temporary submission file
    os.remove("tempSlurm.txt")
    return


def check_positive_float(val) -> float:
    fl = float(val)
    if fl < 0.0:
        raise argparse.ArgumentTypeError(f"{fl} must be >= 0")
    return fl


def make_bins(args: list) -> tuple[float, float]:
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
        tuple[float, float]: low_bin_edges and high_bin edges. A paired list defining
        each bin edge
    """
    low_edges = []
    high_edges = []
    delta = 1e-15
    diff = args[-2] - args[0]
    # case 1: list = [min, max, width]. The width MUST always be smaller than the
    #   difference between the max and min bins. The "or" statement handles floating point
    #   precision if its one bin
    if args[-1] < diff or abs(args[-1] - diff) < delta:
        if len(args) != 3:
            raise (
                RuntimeError(
                    "User gave a custom list, but the last bin value is too"
                    " close to the 2nd to last value"
                )
            )
        min, max, width = args
        n_bins = int((max - min) / width)
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


if __name__ == "__main__":
    # parse arguments passed to script at command line. Reference "help"
    #   arguments for details about each variable
    parser = argparse.ArgumentParser(description="Submit PWA fits to ifarm via slurm")

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
        "--phase_reference",
        type=str,
        metavar="JPmL",
        default="",
        help=(
            "Flag the wave (in 'JPmL' format) whose phase will be constrained to 0."
            " Empty (default) picks lowest JP, m, L combination"
        ),
    )
    parser.add_argument(
        "--phaselock",
        type=bool,
        metavar="bool",
        default=False,
        help=(
            "If True uses 'phaselock' model, where phases are common across"
            " m-projections for a particular eJPL combination"
        ),
    )
    parser.add_argument(
        "-ds",
        "--dsratio",
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
        "--force_refl",
        type=int,
        default=0,
        choices=[-1, 1],
        help="only allow a single reflectivity by setting to +1 or -1",
    )
    parser.add_argument(
        "--init_refl",
        type=int,
        default=0,
        choices=[-1, 1],
        help="initialize single reflectivity by setting to +1 or -1",
    )
    parser.add_argument(
        "--init_re_im",
        type=float,
        nargs=2,
        default=[100, 100],
        metavar=("re", "im"),
        help="value to initialize real and imaginary parts of amplitudes to",
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
        "--runs",
        nargs="*",
        default=["allPeriods"],
        choices=["allPeriods", "2017_01", "2018_01", "2018_08"],
        help="run periods of data/MC to fit",
    )
    parser.add_argument(
        "-m",
        "--mass",
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
        "--t",
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
        "--truth",
        type=str,
        default="",
        help="cfg file that contains truth info used to generate MC",
    )
    parser.add_argument(
        "-d",
        "--data",
        nargs=2,
        default=[
            "data",
            "/w/halld-scshelf2101/kscheuer/neutralb1/submission/source_files/data",
        ],
        metavar=("TYPE", "DIR"),
        help=(
            "specify data version and directory."
            " Data types typically in form: thrown='verXY.Z_mcthrown,"
            " MC='verXY.Z_mc' and DATA='data'"
        ),
    )
    parser.add_argument(
        "-p",
        "--phasespace",
        nargs=2,
        default=[
            "ver03",
            "/w/halld-scshelf2101/kscheuer/neutralb1/submission/source_files/phasespace",
        ],
        metavar=("TYPE", "DIR"),
        help="specify phasespace version and directory.",
    )
    parser.add_argument(
        "-c",
        "--cut_recoil_pi_mass",
        type=float,
        default=0.0,
        help=(
            "Cuts events below the given value in the recoil-pion mass spectrum, i.e."
            " a value of 1.3 selects events like 'MRecoilPi > 1.3' when the trees are"
            " copied. This flag assumes the ROOT data files has an 'MRecoilPi' leaf"
        ),
    )

    # other arguments
    parser.add_argument(
        "--template_name",
        type=str,
        default="template.cfg",
        help="template with default values to copy and overwrite",
    )
    parser.add_argument(
        "-g",
        "--gpu",
        nargs=2,
        default=["0", ""],
        choices=["1", "2", "3", "4", "T4", "TitanRTX", "A800"],
        metavar=("#GPUs", "CARD"),
        help="set # of GPUs to use for a card. Default assumes only CPU fits",
    )
    parser.add_argument(
        "--mail",
        type=bool,
        default=False,
        help=(
            "mails user when a job starts/stops/fails."
            " NOTE this assumes the ifarm username matches the username of"
            "USER@jlab.org"
        ),
    )
    parser.add_argument(
        "--time_limit",
        type=str,
        default="01:00:00",
        help=("Max walltime for each slurm job. Default assumes quick jobs (1 hr)"),
    )

    args = parser.parse_args()

    waveset = args.waveset
    phase_reference = args.phase_reference.lower()
    is_phaselock = args.phaselock
    ds_option = args.dsratio
    frame = args.frame
    force_refl = args.force_refl
    init_refl = args.init_refl
    init_real, init_imag = args.init_re_im

    reaction = args.reaction
    num_rand_fits = args.nrand
    orientations = args.orientations
    run_periods = args.runs
    m_args = args.mass
    t_args = args.t
    energy_min, energy_max = args.energy
    truth_file = args.truth
    data_version, data_dir = args.data
    phasespace_version, phasespace_dir = args.phasespace
    cut_recoil_pi_mass = args.cut_recoil_pi_mass

    template_name = args.template_name
    n_gpus, gpu_type = int(args.gpu[0]), args.gpu[1]
    is_send_mail = args.mail
    time_limit = args.time_limit

    main(
        waveset,
        phase_reference,
        is_phaselock,
        ds_option,
        frame,
        force_refl,
        init_refl,
        init_real,
        init_imag,
        reaction,
        num_rand_fits,
        orientations,
        run_periods,
        t_args,
        energy_min,
        energy_max,
        m_args,
        data_version,
        data_dir,
        phasespace_version,
        phasespace_dir,
        cut_recoil_pi_mass,
        template_name,
        truth_file,
        n_gpus,
        gpu_type,
        is_send_mail,
        time_limit,
    )
