"""Main submission manager for PWA fits.

This module orchestrates the entire submission process by coordinating
data preparation, configuration writing, and job submission.
"""

import math
import os
import pathlib
import pwd
import shutil
from typing import List

import neutralb1.utils as utils

from .config_manager import ConfigManager
from .config_models import PWAConfig
from .data_manager import DataManager
from .job_submitter import JobSubmitter

WORKSPACE_DIR = utils.get_workspace_dir()


class SubmissionManager:
    """Orchestrates the entire PWA submission process.

    This class coordinates data preparation, configuration file writing,
    and job submission for PWA fits across multiple kinematic bins.
    """

    def __init__(self):
        """Initialize the submission manager with component managers."""
        self.config_manager = ConfigManager()
        self.data_manager = DataManager()
        self.job_submitter = JobSubmitter()

    def submit_fits(self, config: PWAConfig) -> List[str]:
        """Submit PWA fits for the given configuration.

        Args:
            config (PWAConfig): Complete configuration object.

        Returns:
            List[str]: List of submitted job IDs.

        Raises:
            ValueError: If configuration validation fails.
            FileNotFoundError: If required files don't exist.
        """
        # Validate configuration
        errors = self.config_manager.validate_config(config)
        if errors:
            raise ValueError(f"Configuration validation failed:\n" + "\n".join(errors))

        # create bins if custom bins not requested
        if config.data.mass_bins == []:
            bins = self._make_bins(
                config.data.mass_min, config.data.mass_max, config.data.mass_width
            )
            bins = [round(b, 3) for b in bins]  # round to single MeV precision
            config.data.mass_bins = bins
        if config.data.t_bins == []:
            bins = self._make_bins(
                config.data.t_min, config.data.t_max, config.data.t_width
            )
            config.data.t_bins = bins

        # Create ROOT data files with cuts if not yet done
        data_job_ids = self.data_manager.prepare_data_files(config)
        if data_job_ids:
            # Data files are being created, exit early
            return data_job_ids

        # Create AmpTools config file template if not using a truth file
        amptools_cfg = (
            config.data.truth_file
            if config.data.truth_file
            else self.config_manager.create_amptools_config(config)
        )

        moment_cfg = self.config_manager.create_moment_config(config)

        # Submit jobs for each combination of run period and kinematic bins
        job_ids = []
        for run_period in config.data.run_periods:
            for low_t, high_t in zip(config.data.t_bins[:-1], config.data.t_bins[1:]):
                for low_mass, high_mass in zip(
                    config.data.mass_bins[:-1], config.data.mass_bins[1:]
                ):
                    job_dir = self._get_job_directory(
                        config, run_period, low_t, high_t, low_mass, high_mass
                    )

                    data_path = self.data_manager.get_volatile_path(
                        config.general.reaction,
                        config.data.cut_recoil_pi_mass,
                        low_t,
                        high_t,
                        config.data.energy[0],
                        config.data.energy[1],
                        low_mass,
                        high_mass,
                    )

                    if config.compute.test:  # provide cfg file to user if running test
                        shutil.copy(amptools_cfg, f"{os.getcwd()}/fit.cfg")
                        shutil.copy(moment_cfg, f"{os.getcwd()}/moments.cfg")
                    else:  # otherwise prepare the job directory by linking needed files
                        self._prepare_job_directory(
                            config,
                            job_dir,
                            amptools_cfg,
                            moment_cfg,
                            data_path,
                            run_period,
                        )

                    # submit and capture the slurm job id
                    try:
                        job_id = self.job_submitter.submit_fit_job(config, job_dir)
                    except Exception as e:
                        print(
                            f"Failed to submit job for {run_period} "
                            f"t=[{low_t:.2f}, {high_t:.2f}] "
                            f"mass=[{low_mass:.3f}, {high_mass:.3f}]: {e}"
                        )
                        job_id = ""

                    if job_id:
                        job_ids.append(job_id)

        return job_ids

    def _make_bins(self, min: float, max: float, width: float):

        delta = 1e-10
        bins = []

        n_bins = math.floor((max - min) / width + delta)

        if n_bins < 1:
            raise ValueError(
                "Invalid binning parameters, (max-min)/width produces 0 bins"
            )
        for i in range(n_bins + 1):  # n_bins+1 to capture last bin edge
            edge = min + width * i
            bins.append(edge)

        return bins

    def _get_job_directory(
        self,
        config: PWAConfig,
        run_period: str,
        low_t: float,
        high_t: float,
        low_mass: float,
        high_mass: float,
    ) -> str:
        """Get the standardized running directory for a specific job

        Args:
            config (PWAConfig): Configuration object.
            run_period (str): Run period string.
            low_t (float): Low t edge.
            high_t (float): High t edge.
            low_mass (float): Low mass edge.
            high_mass (float): High mass edge.

        Returns:
            str: The running directory path.
        """

        volatile_dir = f"/volatile/halld/home/{pwd.getpwuid(os.getuid())[0]}"
        truth_subdir = ""
        if config.data.truth_file:
            if "init" in config.data.truth_file:
                truth_subdir = "truth-init"
            else:
                truth_subdir = "truth"
        return "/".join(
            [
                volatile_dir,
                "ampToolsFits",
                config.general.reaction,
                run_period,
                "_".join(sorted(config.data.orientations)),
                f"{config.data.data_version}{config.data.data_option}",
                f"{config.data.phasespace_version}{config.data.phasespace_option}",
                "_".join(sorted(config.physics.waveset)),
                f"recoil_pi_mass_{config.data.cut_recoil_pi_mass}",
                f"t_{low_t:.2f}-{high_t:.2f}",
                f"mass_{low_mass:.3f}-{high_mass:.3f}",
                truth_subdir,
            ]
        )

    def _prepare_job_directory(
        self,
        config: PWAConfig,
        dir: str,
        amptools_config: str,
        moment_config: str,
        selected_data_path: str,
        run_period: str,
    ) -> None:
        """Create and link the necessary files for the job to run in the kinematic bin

        Args:
            config (PWAConfig): Configuration object.
            dir (str): Directory for the job.
            amptools_config (str): Path to the amptools configuration file.
            moment_config (str): Path to the moment configuration file.
            selected_data_path (str): Path to the data directory holding the
                pre-selected data files. See DataManager.
            run_period (str): Run period string.
        """

        if "init" in config.data.truth_file:
            self._prepare_init_directory(dir)

        # Create directories
        pathlib.Path(dir).mkdir(parents=True, exist_ok=True)

        if config.compute.bootstrap != 0:
            bootstrap_dir = f"{dir}/bootstrap/"
            pathlib.Path(bootstrap_dir).mkdir(parents=True, exist_ok=True)

        if config.compute.nrand != 0:
            rand_dir = f"{dir}/rand/"
            pathlib.Path(rand_dir).mkdir(parents=True, exist_ok=True)

        shutil.copy(amptools_config, f"{dir}/fit.cfg")
        shutil.copy(moment_config, f"{dir}/moments.cfg")
        # save the submission config to the dir so its clear what parameters were used
        self.config_manager.save_yaml_config(config, f"{dir}/config.yaml")

        # copy gluex environment scripts
        if int(config.compute.gpu[0]) > 0:
            shutil.copy(
                f"{WORKSPACE_DIR}/config/setup_gpu_gluex.sh", f"{dir}/setup_gluex.sh"
            )
            shutil.copy(f"{WORKSPACE_DIR}/config/version_gpu.xml", f"{dir}/version.xml")
        else:
            shutil.copy(
                f"{WORKSPACE_DIR}/config/setup_gluex.sh", f"{dir}/setup_gluex.sh"
            )
            shutil.copy(f"{WORKSPACE_DIR}/config/version.xml", f"{dir}/version.xml")

        # link data phasespace files
        for ont in config.data.orientations:
            # data file should be labelled with just the orientation angle to match
            # the AmpTools cfg file
            ont_number = ont.replace("PARA_", "").replace("PERP_", "")
            self._symlink_force(
                (
                    f"{selected_data_path}/AmpToolsInputTree_sum_{ont}_{run_period}"
                    f"_{config.data.data_version}{config.data.data_option}.root"
                ),
                (f"{dir}/anglesOmegaPiAmplitude_{ont_number}.root"),
            )
        # gen phasespace
        self._symlink_force(
            (
                f"{selected_data_path}/anglesOmegaPiPhaseSpaceGen_{run_period}"
                f"_{config.data.phasespace_version}{config.data.phasespace_option}.root"
            ),
            (f"{dir}/anglesOmegaPiPhaseSpace.root"),
        )

        # acc phasespace (won't apply for thrown MC studies)
        label = "Gen" if "mcthrown" in config.data.data_option else "Acc"
        self._symlink_force(
            (
                f"{selected_data_path}/anglesOmegaPiPhaseSpace{label}_{run_period}"
                f"_{config.data.phasespace_version}{config.data.phasespace_option}.root"
            ),
            (f"{dir}/anglesOmegaPiPhaseSpaceAcc.root"),
        )

    def _prepare_init_directory(
        self, running_dir: str, scale_par_name: str = "intensity_scale"
    ) -> None:
        """Create the necessary scale.txt file for truth initialized fits

        A truth initialized fit requires a scale.txt file that simply has an
        intensity_scale parameter that multiplies all amplitudes by a constant factor,
        in order to properly initialize the amplitudes to the right value. This function
        creates that file by assuming that a truth fit has already been run and the
        scale factor is stored in the best_truth.fit file. If the best_truth.fit file is
        not found, the submission stops.

        Args:
            running_dir (str): Directory where the scale.txt file will be created.
            scale_par_name (str): Name of the scale parameter written in the truth-init
                cfg file that needs to be sourced from a previous fixed truth fit

        Raises:
            FileNotFoundError: If the best_truth.fit file is not found.
            ValueError: If the scale factor is not found in the best_truth.fit file.
        """

        # running_dir (truth-init) is a sibling of the truth directory
        truth_dir = f"{running_dir.rsplit('/', 2)[0]}/truth/"

        # Find the scale factor from the best_truth.fit file
        truth_fit = f"{truth_dir}best_truth.fit"
        if not os.path.isfile(truth_fit):
            raise FileNotFoundError(
                f"File {truth_fit} not found! Truth initialized fits require this file."
                " Please run a truth fit first."
            )

        scale_factor = 0.0
        is_searching = False
        with open(truth_fit, "r") as f:
            for line in f:
                # Begin search within the parameter values and errors section
                if "Parameter Values and Errors" in line:
                    is_searching = True
                    continue
                if not is_searching:
                    continue

                if scale_par_name in line.split()[0]:
                    scale_factor = float(line.split()[1])
                    break

                # If we've hit the normalization integrals, we've gone too far
                if "Normalization Integrals" in line:
                    break

        if scale_factor == 0.0:
            raise ValueError("Scale factor was not found in the best_truth.fit file!")

        # Write the scale factor to the scale.txt file
        with open(f"{running_dir}/scale.txt", "w") as f:
            f.write(f"parameter intensity_scale {scale_factor} fixed")

    def _symlink_force(self, target: str, link_name: str) -> None:
        """Forcefully create a symlink, overwriting any existing link or file.

        Args:
            target (str): The target file to link to.
            link_name (str): The name of the symlink to create.
        """
        try:
            os.symlink(target, link_name)
        except FileExistsError:
            os.remove(link_name)
            os.symlink(target, link_name)
        except Exception as e:
            print(f"Failed to create symlink {link_name} -> {target}: {e}")
