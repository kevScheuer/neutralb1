"""Configuration models for the PWA submission system.

This module defines the data structures for storing PWA fit configurations
using dataclasses and YAML serialization.
"""

# TODO; see if there's some way to easily set defaults, and convey what options are
# available at submission

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import neutralb1.utils as utils

WORKSPACE_DIR = utils.get_workspace_dir()


def _format_cut_value(value: float) -> str:
    """Format cut values consistently for paths and CLI arguments.

    Args:
        value (float): Numeric cut value.

    Returns:
        str: String representation with at least one decimal place.
    """

    formatted = f"{value:.6f}".rstrip("0").rstrip(".")
    if "." not in formatted:
        formatted = f"{formatted}.0"
    return formatted


# Central definition of the nominal (default) systematic-cut ranges.
# Keys must match the variable names used by copy_tree_with_cuts.
NOMINAL_CUT_RANGES: Dict[str, Tuple[float, float]] = {
    "unusedE": (-0.1, 0.1),
    "unusedTracks": (-0.1, 1.0),
    "chi2": (-0.1, 5.0),
    "PzCMrecoilPi": (-0.1, 100.0),
    "shQualityP2a": (0.5, 1.1),
    "shQualityP2b": (0.5, 1.1),
    "shQualityP3a": (0.5, 1.1),
    "shQualityP3b": (0.5, 1.1),
}

# shQuality sub-field names affected by the shQuality group cut.
_SH_QUALITY_SUFFIXES = ("P2a", "P2b", "P3a", "P3b")
_SH_QUALITY_FIELDS = tuple(f"shQuality{s}" for s in _SH_QUALITY_SUFFIXES)


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
        ds_ratio (str): D/S wave ratio constraint ("", "free", "fixed", "split",
            "positive", "negative").
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
        shQuality (List[float]): Convenience field to set all four shQuality*
            sub-fields (P2a, P2b, P3a, P3b) simultaneously. When non-empty,
            overrides the individual shQualityP2a/b/P3a/b fields on
            initialization. Leave empty to set each sub-field individually.
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
    # Below are the nominal cut values, which we can vary for systematics

    # cuts that have lower bound of -0.1, aside from PzCMrecoilPi, cant' go negative,
    # and so this ensures that the nominal selection includes all events at or above 0.
    unusedE: List[float] = field(default_factory=lambda: [-0.1, 0.1])
    unusedTracks: List[float] = field(default_factory=lambda: [-0.1, 1.0])
    chi2: List[float] = field(default_factory=lambda: [-0.1, 5.0])
    PzCMrecoilPi: List[float] = field(
        default_factory=lambda: [-0.1, 100.0]
    )  # 100 basically ensures no upper cut
    shQualityP2a: List[float] = field(
        default_factory=lambda: [0.5, 1.1]
    )  # upper limit it can be is 1.0
    shQualityP2b: List[float] = field(default_factory=lambda: [0.5, 1.1])
    shQualityP3a: List[float] = field(default_factory=lambda: [0.5, 1.1])
    shQualityP3b: List[float] = field(default_factory=lambda: [0.5, 1.1])
    # Group convenience field — when set, propagates to all four sub-fields above.
    shQuality: List[float] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Propagate the shQuality group field to all four sub-fields.

        If shQuality is provided (non-empty), it overrides shQualityP2a,
        shQualityP2b, shQualityP3a, and shQualityP3b so that a single YAML
        entry controls all four photon shower-quality cuts simultaneously.
        """
        if self.shQuality:
            for field_name in _SH_QUALITY_FIELDS:
                setattr(self, field_name, list(self.shQuality))

    def get_cut_ranges(self) -> Dict[str, Tuple[float, float]]:
        """Return the configured cut ranges as a mapping.

        Returns:
            Dict[str, Tuple[float, float]]: Map of variable name to (min, max).
        """

        return {
            "unusedE": (float(self.unusedE[0]), float(self.unusedE[1])),
            "unusedTracks": (
                float(self.unusedTracks[0]),
                float(self.unusedTracks[1]),
            ),
            "chi2": (float(self.chi2[0]), float(self.chi2[1])),
            "PzCMrecoilPi": (
                float(self.PzCMrecoilPi[0]),
                float(self.PzCMrecoilPi[1]),
            ),
            "shQualityP2a": (
                float(self.shQualityP2a[0]),
                float(self.shQualityP2a[1]),
            ),
            "shQualityP2b": (
                float(self.shQualityP2b[0]),
                float(self.shQualityP2b[1]),
            ),
            "shQualityP3a": (
                float(self.shQualityP3a[0]),
                float(self.shQualityP3a[1]),
            ),
            "shQualityP3b": (
                float(self.shQualityP3b[0]),
                float(self.shQualityP3b[1]),
            ),
        }

    def get_systematics_tag(self) -> str:
        """Compute the systematic tag for paths based on deviations from nominal.

        This enforces that at most one cut range differs from the nominal range, so
        that the systematic is uniquely identifiable (e.g. "chi2_0.0_4.0"). If
        nothing differs, returns "nominal".

        Raises:
            ValueError: If more than one cut range differs from nominal.

        Returns:
            str: Systematics tag string for directory naming.
        """

        diffs = self._get_non_nominal_cuts()
        if not diffs:
            return "nominal"

        if len(diffs) > 1:
            diff_names = ", ".join(sorted(diffs.keys()))
            raise ValueError(
                "Multiple systematic cuts differ from nominal. Set only one at a "
                f"time. Differences: {diff_names}"
            )

        var, (low, high) = next(iter(diffs.items()))
        return f"{var}_{_format_cut_value(low)}_{_format_cut_value(high)}"

    def get_optional_cut_args(self) -> List[str]:
        """Build optional cut arguments for copy_tree_with_cuts.

        By default, this passes *all* nominal cut ranges to copy_tree_with_cuts so the
        produced per-bin ROOT files always include the nominal selection. If exactly
        one configured cut differs from nominal, that one cut range replaces the
        nominal range in the arguments.

        Raises:
            ValueError: If more than one cut range differs from nominal.

        Returns:
            List[str]: Optional cut arguments for copy_tree_with_cuts, e.g.
            ["chi2:0.0:5.0", "unusedE:0.0:0.1", ...].
        """

        diffs = self._get_non_nominal_cuts()
        if len(diffs) > 1:
            diff_names = ", ".join(sorted(diffs.keys()))
            raise ValueError(
                "Multiple systematic cuts differ from nominal. Set only one at a "
                f"time. Differences: {diff_names}"
            )

        merged: Dict[str, Tuple[float, float]] = dict(NOMINAL_CUT_RANGES)
        if diffs:
            var, (low, high) = next(iter(diffs.items()))
            if var == "shQuality":
                # Expand the group key back to the four individual CLI keys.
                for field_name in _SH_QUALITY_FIELDS:
                    merged[field_name] = (float(low), float(high))
            else:
                merged[var] = (float(low), float(high))

        return [
            f"{var}:{_format_cut_value(low)}:{_format_cut_value(high)}"
            for var, (low, high) in sorted(merged.items())
        ]

    def get_systematic_cut_range(self) -> Optional[Tuple[str, float, float]]:
        """Return the unique non-nominal systematic cut, if any.

        Raises:
            ValueError: If more than one cut range differs from nominal.

        Returns:
            Optional[Tuple[str, float, float]]: (variable, low, high) if a single cut
            differs from nominal; otherwise None.
        """

        diffs = self._get_non_nominal_cuts()
        if not diffs:
            return None

        if len(diffs) > 1:
            diff_names = ", ".join(sorted(diffs.keys()))
            raise ValueError(
                "Multiple systematic cuts differ from nominal. Set only one at a "
                f"time. Differences: {diff_names}"
            )

        var, (low, high) = next(iter(diffs.items()))
        return var, low, high

    def _get_non_nominal_cuts(self) -> Dict[str, Tuple[float, float]]:
        """Return mapping of all cut ranges that differ from nominal.

        When all four shQuality* sub-fields share the same non-nominal range,
        they are collapsed into a single ``"shQuality"`` entry so that the
        caller sees at most one logical systematic variation.
        """

        non_nominal: Dict[str, Tuple[float, float]] = {}
        current = self.get_cut_ranges()
        for var, (low, high) in current.items():
            nominal_low, nominal_high = NOMINAL_CUT_RANGES[var]
            if (low, high) != (nominal_low, nominal_high):
                non_nominal[var] = (low, high)

        # Collapse all four shQuality sub-fields into a single "shQuality" entry
        # when they all carry the same non-nominal range.
        sh_diffs = {v: r for v, r in non_nominal.items() if v in _SH_QUALITY_FIELDS}
        if sh_diffs and len(sh_diffs) == len(_SH_QUALITY_FIELDS):
            ranges = set(sh_diffs.values())
            if len(ranges) == 1:  # all four are identical
                for v in _SH_QUALITY_FIELDS:
                    del non_nominal[v]
                non_nominal["shQuality"] = ranges.pop()

        return non_nominal


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
        cpus_per_task (int): Number of CPU cores per task (only used if n_gpus=0).
        threads_per_core (int): Number of threads per core (only used if n_gpus=0).
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
    cpus_per_task: int = 1
    threads_per_core: int = 1
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
