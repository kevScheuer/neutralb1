"""Job submission module for PWA fits.

This module handles SLURM job submission and management for PWA fits,
including GPU and CPU configurations.
"""

import os
import pathlib
import re
import shutil
import subprocess
import tempfile
import time
from typing import List

from .config_models import PWAConfig


class JobSubmitter:
    """Handles SLURM job submission for PWA fits.

    This class manages the creation and submission of SLURM batch scripts
    for both GPU and CPU-based fits.
    """

    def submit_fit_job(
        self,
        config: PWAConfig,
        job_dir: str,
    ) -> str:
        """Submit a PWA fit job via SLURM to the ifarm

        Args:
            config (PWAConfig): Configuration object.
            job_dir (str): Directory where the job will run.

        Returns:
            str: Job ID of the submitted job.
        """
        job_name = job_dir.replace("/", "_")
        script_dir = os.path.dirname(os.path.abspath(__file__))

        # Build truth file argument
        truth_arg = 1 if config.data.truth_file else 0

        script_command = (
            f"{script_dir}/run_fit.sh "
            f"-n {config.compute.nrand} "
            f"-r {config.general.reaction} "
            f"-b {config.compute.bootstrap} "
            f"-s {config.compute.random_seed} "
            f"-t {truth_arg}"
        )

        log_dir = self._get_log_directory(config, job_dir)

        # Extract GPU configuration
        n_gpus = int(config.compute.gpu[0])
        gpu_type = config.compute.gpu[1] if len(config.compute.gpu) > 1 else ""

        job_id = self.submit_slurm_job(
            job_name,
            script_command,
            job_dir,
            log_dir,
            gpu_type,
            n_gpus,
            config.compute.email,
            config.compute.email_type,
            config.compute.time_limit,
            config.compute.mem_per_cpu,
            config.compute.n_cpus,
            config.compute.test,
        )

        return job_id

    def submit_data_job(self, config: PWAConfig, src_file: str, job_dir: str) -> str:
        """Submit a job to create cut data files.

        Args:
            config (PWAConfig): Configuration object.
            src_file (str): Source ROOT file path.
            job_dir (str): Directory where the job will run.
        """
        job_name = job_dir.replace("/", "_")

        # Use regex to extract t and mass bin values from the directory path
        t_match = re.search(r"t_([0-9.]+)-([0-9.]+)", job_dir)
        mass_match = re.search(r"mass_([0-9.]+)-([0-9.]+)", job_dir)
        energy_match = re.search(r"E_([0-9.]+)-([0-9.]+)", job_dir)

        if t_match and mass_match and energy_match:
            low_t = float(t_match.group(1))
            high_t = float(t_match.group(2))
            low_mass = float(mass_match.group(1))
            high_mass = float(mass_match.group(2))
            low_energy = float(energy_match.group(1))
            high_energy = float(energy_match.group(2))
        else:
            raise ValueError(f"Could not parse t/mass bins from directory: {job_dir}")

        # build shell commands directly into files since its simpler
        script_command = (
            "source setup_gluex.sh &&"
            "export PATH="  # path export hard-coded for now
            '"/w/halld-scshelf2101/kscheuer/neutralb1/build/release/bin:$PATH"'
            " && copy_tree_with_cuts"
            f" {src_file}"
            f" {job_dir}"
            f" {config.data.cut_recoil_pi_mass}"
            f" {low_t} {high_t}"
            f" {low_energy} {high_energy}"
            f" {low_mass} {high_mass}"
        )

        # Create directories
        log_dir = self._get_log_directory(config, job_dir)

        job_id = self.submit_slurm_job(
            job_name,
            script_command,
            job_dir,
            log_dir,
            "",  # GPU unnecessary for data file jobs
            0,
            config.compute.email,
            config.compute.email_type,
            "00:30:00",  # copy_tree_with_cuts should not take long
            "5000M",  # 5 G should be enough memory
            8,  # MPI cant help, so request small number of GPUs
            config.compute.test,
        )

        return job_id

    def submit_slurm_job(
        self,
        job_name: str,
        script_command: str,
        running_dir: str,
        log_dir: str,
        gpu_type: str,
        n_gpus: int,
        email_address: str,
        email_type: List[str],
        time_limit: str,
        mem_per_cpu: str,
        n_cpus: int,
        is_test: bool = False,
    ) -> str:
        """Submit a SLURM job to the ifarm

        Args:
            job_name (str): Name shown on the scicomp webpage.
            script_command (str): Bash script with its arguments.
            running_dir (str): /volatile/ location.
            log_dir (str): Where SLURM log files are stored.
            gpu_type (str): Card type to be used.
            n_gpus (int): How many GPU cards to use (supported by MPI).
            email_address (str): Send email to passed address.
            email_type (List[str]): When to send email (BEGIN, END, FAIL).
            time_limit (str): Max wall-time in Hour:Min:Sec.
            mem_per_cpu (str): Memory per CPU
            n_cpus (int): Number of MPI CPUs to use (only used if n_gpus=0).
            is_test (bool): If True, will not submit the job to the farm.

        Returns:
            str: Job ID of submitted job, or "" if in test mode.
        """
        slurm_script_content = self._generate_slurm_script(
            job_name=job_name,
            script_command=script_command,
            running_dir=running_dir,
            log_dir=log_dir,
            gpu_type=gpu_type,
            n_gpus=n_gpus,
            email_address=email_address,
            email_type=email_type,
            time_limit=time_limit,
            mem_per_cpu=mem_per_cpu,
            n_cpus=n_cpus,
        )

        temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
        temp_file.write(slurm_script_content)
        temp_file.close()

        if is_test:  # provide slurm script to user if running test
            shutil.copy(temp_file.name, f"{os.getcwd()}/slurm.txt")
            return ""

        # Submit the job
        time.sleep(0.5)  # Avoid job skip error if too many submitted quickly
        result = subprocess.run(
            ["sbatch", temp_file.name],
            capture_output=True,
            text=True,
        )

        # Extract job ID from sbatch output
        if result.returncode == 0:
            # sbatch output is typically "Submitted batch job JOBID"
            job_id = result.stdout.strip().split()[-1]
            return job_id
        else:
            raise RuntimeError(f"Failed to submit job: {result.stderr}")

    def _generate_slurm_script(
        self,
        job_name: str,
        script_command: str,
        running_dir: str,
        log_dir: str,
        gpu_type: str,
        n_gpus: int,
        email_address: str,
        email_type: List[str],
        time_limit: str,
        mem_per_cpu: str,
        n_cpus: int,
    ) -> str:
        """Create slurm script using defaults and user-provided values.

        Args:
            job_name (str): SLURM job name.
            script_command (str): Command to execute.
            running_dir (str): Working directory.
            log_dir (str): Log directory.
            gpu_type (str): GPU type.
            n_gpus (int): Number of GPUs.
            email_address (str): Email for notifications.
            email_type (List[str]): Email notification types.
            time_limit (str): Time limit.
            mem_per_cpu (str): Memory per CPU.
            n_cpus (int): Number of CPUs.

        Returns:
            str: Complete SLURM script content.
        """
        script_lines = [
            "#!/bin/sh",
            "#SBATCH -A halld",
            f"#SBATCH --time={time_limit}",
            f"#SBATCH --chdir={running_dir}",
            f"#SBATCH --error={log_dir}/log.err",
            f"#SBATCH --output={log_dir}/log.out",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --mem-per-cpu={mem_per_cpu}",
            "#SBATCH --cpus-per-task=1",
            "#SBATCH --ntasks-per-core=1",
            "#SBATCH --threads-per-core=1",
            "#SBATCH --constraint=el9",
        ]

        # Add email configuration if provided
        if email_address:
            mail_type = ",".join(email_type)
            script_lines.extend(
                [
                    f"#SBATCH --mail-user={email_address}",
                    f"#SBATCH --mail-type={mail_type}",
                ]
            )

        # Add GPU or CPU specific configuration
        if n_gpus > 0:
            script_lines.extend(
                [
                    "#SBATCH --partition=gpu",
                    f"#SBATCH --gres=gpu:{gpu_type}:{n_gpus}",
                    f"#SBATCH --ntasks={n_gpus+1}",  # mpigpu always needs n_gpus+1
                ]
            )
        else:
            script_lines.extend(
                [
                    "#SBATCH --partition=ifarm",
                    f"#SBATCH --ntasks={n_cpus}",
                ]
            )

        # Add the actual command
        script_lines.append(script_command)

        return "\n".join(script_lines)

    def _get_log_directory(self, config: PWAConfig, job_dir: str) -> str:
        """Get the log directory for a job.

        Creates the log dir if not running a test case. Log directory is in dedicated
        farm directory for I/O heavy operations.

        Args:
            config (PWAConfig): The PWA configuration object.
            job_dir (str): The job directory.

        Returns:
            str: The log directory path.
        """
        log_dir = f"{job_dir.replace('/volatile/halld/home/', '/farm_out/')}/log/"
        if not config.compute.test:
            pathlib.Path(log_dir).mkdir(parents=True, exist_ok=True)

        return log_dir
