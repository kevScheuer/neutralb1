"""Data management module for PWA submissions.

This module handles the creation and management of data files for PWA fits,
including cutting data to specific kinematic regions and managing file paths.
"""

import os
import pathlib
import pwd
import shutil

import neutralb1.utils as utils

from .config_models import PWAConfig
from .job_submitter import JobSubmitter

WORKSPACE_DIR = utils.get_workspace_dir()


class DataManager:
    """Manages data file creation and preparation for PWA fits.

    This class handles the copying and cutting of ROOT data files to the
    volatile directory with appropriate kinematic cuts applied.
    """

    def __init__(self):
        self.user = pwd.getpwuid(os.getuid())[0]
        self.volatile_dir = f"/volatile/halld/home/{self.user}"
        self.job_submitter = JobSubmitter()

    def prepare_data_files(self, config: PWAConfig) -> bool:
        """Create data files with cuts and store them in volatile directory.

        The typical ROOTDataReader method in .cfg files reads in data much too slowly,
        and is repetitive when the same TEM region is being selected. This function
        will copy the source data files to the volatile directory, and cut them to the
        desired TEM region for quick access.

        Args:
            config (PWAConfig): Complete configuration object.

        Returns:
            bool: True if files are ready, False if jobs were submitted to create them.

        Raises:
            FileExistsError: If required source files don't exist.
        """
        # Extract energy range
        energy_min, energy_max = config.data.energy

        # Track what files need to be copied to what directories on volatile
        src_files_to_copy_to_dir = {}

        # If user wants to skip being asked about every job submission,
        # this will become True
        skip_input = False

        job_ids = []

        for run_period in config.data.run_periods:
            # Find generated phasespace file
            gen_file = (
                f"{config.data.phasespace_dir}/anglesOmegaPiPhaseSpaceGen_"
                f"{run_period}_{config.data.phasespace_version}.root"
            )
            if not os.path.isfile(gen_file):
                raise FileNotFoundError(f"Path {gen_file} does not exist!")
            src_files_to_copy_to_dir[gen_file] = []

            # Find accepted phasespace file. Thrown MC accepted files will be a copy
            # of the generated phasespace file
            if "mcthrown" not in config.data.data_option:
                acc_file = (
                    f"{config.data.phasespace_dir}/anglesOmegaPiPhaseSpaceAcc_"
                    f"{run_period}_{config.data.phasespace_version}"
                    f"{config.data.phasespace_option}.root"
                )
                if not os.path.isfile(acc_file):
                    raise FileNotFoundError(f"Path {acc_file} does not exist!")
                src_files_to_copy_to_dir[acc_file] = []

            # Find data files
            for ont in config.data.orientations:
                data_file = (
                    f"{config.data.data_dir}/AmpToolsInputTree_sum_{ont}_{run_period}"
                    f"_{config.data.data_version}{config.data.data_option}.root"
                )
                if not os.path.isfile(data_file):
                    raise FileNotFoundError(f"Path {data_file} does not exist!")
                src_files_to_copy_to_dir[data_file] = []

            # Loop over TEM bins to determine what bins need the cut data files
            for low_mass, high_mass in zip(
                config.data.mass_bins[:-1], config.data.mass_bins[1:]
            ):
                for low_t, high_t in zip(
                    config.data.t_bins[:-1], config.data.t_bins[1:]
                ):
                    volatile_path = self.get_volatile_path(
                        config.general.reaction,
                        config.data.cut_recoil_pi_mass,
                        low_t,
                        high_t,
                        energy_min,
                        energy_max,
                        low_mass,
                        high_mass,
                    )

                    # Check if files exist for this bin
                    for src_file in src_files_to_copy_to_dir.keys():
                        base_name = os.path.basename(src_file)
                        cut_file = f"{volatile_path}/{base_name}"

                        if not os.path.isfile(cut_file):
                            src_files_to_copy_to_dir[src_file].append(volatile_path)

            # If any files need to be cut and copied for this run period
            if any(dirs for dirs in src_files_to_copy_to_dir.values() if dirs):
                if not skip_input and not config.compute.test:
                    response = input(
                        f"Data files for {run_period} need to be created. "
                        f"Submit jobs to create them? (y/n/a for all): "
                    )
                    if response.lower() == "a":
                        skip_input = True
                    elif response.lower() != "y":
                        continue

                # Submit jobs to create the files
                for src_file, dirs in src_files_to_copy_to_dir.items():
                    if config.compute.test:
                        print(
                            f"Test mode: file {src_file} would be copied to directories"
                            f"\n\t{'\n\t'.join(dirs)}"
                        )

                    for dir in dirs:
                        if not config.compute.test:
                            self._prepare_job_directory(dir)
                        try:
                            job_id = self.job_submitter.submit_data_job(
                                config, src_file, dir
                            )
                        except Exception as e:
                            print(
                                f"Failed to submit job to copy {src_file} to {dir}: {e}"
                            )
                            job_id = ""

                        if job_id:
                            job_ids.append(job_id)

        if job_ids and not config.compute.test:
            print(
                "Jobs have been submitted to create the necessary data files."
                " Please wait for them to finish before submitting PWA fits."
                " Job progress can be monitored at"
                " https://scicomp.jlab.org/scicomp/slurmJob/activeJob, or by running"
                " 'squeue --me' in a terminal"
            )
            return False

        return True

    def get_volatile_path(
        self,
        reaction: str,
        cut_recoil_pi_mass: float,
        low_t: float,
        high_t: float,
        energy_min: float,
        energy_max: float,
        low_mass: float,
        high_mass: float,
    ) -> str:
        """Create consistent volatile path for pre-selected data files.

        Args:
            reaction (str): Reaction name.
            cut_recoil_pi_mass (float): Recoil pion mass cut.
            low_t (float): Low t edge.
            high_t (float): High t edge.
            energy_min (float): Minimum beam energy.
            energy_max (float): Maximum beam energy.
            low_mass (float): Low mass edge.
            high_mass (float): High mass edge.

        Returns:
            str: Path to volatile directory for this kinematic bin.
        """
        return "/".join(
            [
                self.volatile_dir,
                "ampToolsFits",
                reaction,
                "data_files",
                f"recoil-pi-mass_{cut_recoil_pi_mass}",
                f"t_{low_t:.2f}-{high_t:.2f}",
                f"E_{energy_min:.2f}-{energy_max:.2f}",
                f"mass_{low_mass:.3f}-{high_mass:.3f}",
            ]
        )

    def _prepare_job_directory(self, dir: str) -> None:
        """Prepare the job directory for a data job.

        Args:
            dir (str): The directory to prepare.
        """

        pathlib.Path(dir).mkdir(parents=True, exist_ok=True)
        # copy in cpu-based gluex env
        shutil.copy(
            f"{WORKSPACE_DIR}/config/setup_gluex.sh",
            f"{dir}/setup_gluex.sh",
        )
        shutil.copy(f"{WORKSPACE_DIR}/config/version.xml", f"{dir}/version.xml")

        return
