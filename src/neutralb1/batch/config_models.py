"""Configuration models for the PWA submission system.

This module defines the data structures for storing PWA fit configurations
using dataclasses and YAML serialization.
"""

# TODO; see if there's some way to easily set defaults, and convey what options are
# available at submission

from dataclasses import dataclass, field
from typing import List

import neutralb1.utils as utils

WORKSPACE_DIR = utils.get_workspace_dir()


@dataclass
class PhysicsConfig:
    """Configuration for physics parameters of the fit.

    Args:
        waveset (List[str]): Waveset to fit with (e.g., ["1p", "1m", "2p"]).
        max_moment_J (int): Maximum moment J value for fitting with moments. When not
            set to zero, all moments up to this value will be included, not just those
            automatically determined from the waveset
        phase_reference (List[str]): Phase reference waves in eJPmL format.
        phaselock (bool): Enable phaselock model for common phases across m-projections.
        ds_ratio (str): D/S wave ratio constraint ("", "free", "fixed", "split").
        frame (str): Decay frame ("", "GJ", "Adair").
        single_refl (int): Single reflectivity constraint (-1, 0, or 1).
        init_refl (int): Initialize reflectivity (-1, 0, or 1).
        init_real (float): Real part initialization value.
        init_imag (float): Imaginary part initialization value.
        remove_waves (List[str]): Waves to remove from waveset in eJPmL format.
    """

    waveset: List[str] = field(default_factory=list)
    max_moment_J: int = 0
    phase_reference: List[str] = field(default_factory=list)
    phaselock: bool = False
    ds_ratio: str = ""
    frame: str = ""
    single_refl: int = 0
    init_refl: int = 0
    init_real: float = 100.0
    init_imag: float = 100.0
    remove_waves: List[str] = field(default_factory=list)


@dataclass
class DataConfig:
    """Configuration for data files and binning.

    Args:
        orientations (List[str]): Diamond orientations.
        run_periods (List[str]): GlueX run periods.
        mass_bins (List[float]): Mass bin edges
        mass_min (float): Minimum resonance mass value.
        mass_max (float): Maximum resonance mass value.
        mass_width (float): Width of resonance mass bins.
        t_bins (List[float]): Mandelstam -t bin edges
        t_min (float): Minimum Mandelstam -t value.
        t_max (float): Maximum Mandelstam -t value.
        t_width (float): Width of Mandelstam -t bins.
        energy (List[float]): Beam energy range [min, max].
        data_dir (str): Directory containing data files.
        data_version (str): Data version string.
        data_option (str): Data option string.
        phasespace_dir (str): Directory containing phasespace files.
        phasespace_version (str): Phasespace version string.
        phasespace_option (str): Phasespace option string.
        cut_recoil_pi_mass (float): Recoiling_proton + pion mass cut.
        truth_file (str): Optional truth file for MC studies.
    """

    orientations: List[str] = field(default_factory=lambda: ["PARA_0"])
    run_periods: List[str] = field(default_factory=lambda: ["allPeriods"])
    tree_name: str = "ntFSGlueX_100_112"
    num_final_state_particles: int = 5
    mass_bins: List[float] = field(default_factory=list)
    mass_min: float = 0
    mass_max: float = 0
    mass_width: float = 0
    t_bins: List[float] = field(default_factory=list)
    t_min: float = 0
    t_max: float = 0
    t_width: float = 0
    energy: List[float] = field(default_factory=lambda: [8.2, 8.8])
    data_dir: str = f"{WORKSPACE_DIR}/data/FSRoot/GlueX"
    data_version: str = "data"
    data_option: str = ""
    phasespace_dir: str = f"{WORKSPACE_DIR}/data/FSRoot/phasespace"
    phasespace_version: str = "ver03"
    phasespace_option: str = ""
    truth_file: str = ""
    # TODO: for systematics, could add additional cut parameters here


@dataclass
class ComputeConfig:
    """Configuration for computational resources and job settings.

    Args:
        nrand (int): Number of random fits.
        random_seed (int): Seed for randomized fits (0 = use time).
        bootstrap (int): Number of bootstrap fits.
        mem_per_cpu (str): Memory per CPU
        n_cpus (int): Number of CPUs for MPI jobs
        gpu (List[str]): GPU configuration [n_gpus, gpu_type].
        email (str): Email address for notifications.
        email_type (List[str]): When to send emails.
        time_limit (str): SLURM time limit.
        test (bool): Test mode (don't actually submit jobs).
    """

    nrand: int = 20
    random_seed: int = 0
    bootstrap: int = 0
    mem_per_cpu: str = "5000M"
    n_cpus: int = 32
    gpu: List[str] = field(default_factory=lambda: ["0", ""])
    email: str = ""
    email_type: List[str] = field(default_factory=lambda: ["BEGIN", "END", "FAIL"])
    time_limit: str = "01:00:00"
    test: bool = False


@dataclass
class GeneralConfig:
    """General configuration options.

    Args:
        reaction (str): Base reaction name.
        template_name (str): Template configuration file name.
        n_events (int): Number of events to initialize / sample moments with
    """

    reaction: str = "omegapi"
    template_name: str = "template.cfg"
    n_events: int = 1000


@dataclass
class PWAConfig:
    """Complete PWA fit configuration.

    Args:
        name (str): Configuration name for identification.
        description (str): Human-readable description.
        physics (PhysicsConfig): Physics parameters.
        data (DataConfig): Data and binning configuration.
        compute (ComputeConfig): Computational resources.
        general (GeneralConfig): General options.
    """

    name: str = ""
    description: str = ""
    physics: PhysicsConfig = field(default_factory=PhysicsConfig)
    data: DataConfig = field(default_factory=DataConfig)
    compute: ComputeConfig = field(default_factory=ComputeConfig)
    general: GeneralConfig = field(default_factory=GeneralConfig)
