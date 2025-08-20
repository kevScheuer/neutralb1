"""Configuration file management for PWA submissions.

This module handles AmpTools configuration file creation and YAML
configuration file parsing.
"""

import os
import pathlib
import tempfile
from typing import Any, Dict, List

import yaml

from .amptools_config_writer import AmpToolsConfigWriter
from .config_models import (
    ComputeConfig,
    DataConfig,
    GeneralConfig,
    PhysicsConfig,
    PWAConfig,
)
from .moment_config_writer import MomentConfigWriter


class ConfigManager:
    """
    Manages PWA configuration files and YAML parsing.

    This class handles the creation of AmpTools configuration files
    and the loading/saving of YAML configuration files.
    """

    def __init__(self):
        self.code_dir = str(pathlib.Path(__file__).resolve().parent) + "/"

    def load_yaml_config(self, config_path: str) -> PWAConfig:
        """Load configuration from YAML file.

        Args:
            config_path (str): Path to YAML configuration file.

        Returns:
            PWAConfig: Parsed configuration object.

        Raises:
            FileNotFoundError: If config file doesn't exist.
            yaml.YAMLError: If YAML parsing fails.
        """
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        try:
            with open(config_path, "r") as f:
                config_dict = yaml.safe_load(f)

            # Convert YAML structure to PWAConfig
            return self._yaml_dict_to_config(config_dict)

        except yaml.YAMLError as e:
            raise yaml.YAMLError(f"Error parsing YAML file {config_path}: {e}")

    def save_yaml_config(self, config: PWAConfig, config_path: str) -> None:
        """Save configuration to YAML file.

        Args:
            config (PWAConfig): Configuration object to save.
            config_path (str): Path where to save the YAML file.
        """
        config_dict = self._config_to_yaml_dict(config)

        with open(config_path, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)

    def create_amptools_config(self, config: PWAConfig) -> str:
        """Create AmpTools configuration file.

        Args:
            config (PWAConfig): Configuration object with user parameters
        Returns:
            str: Path to the created .cfg file.
        """
        cfg_writer = AmpToolsConfigWriter(config)
        amptools_cfg = tempfile.NamedTemporaryFile(delete=False, mode="w")
        cfg_writer.write_config_to_file(amptools_cfg.name)
        amptools_cfg.close()

        return amptools_cfg.name

    def create_moment_config(self, config: PWAConfig) -> str:
        """Create moment configuration file.

        Args:
            config (PWAConfig): Configuration object with user parameters

        Returns:
            str: Path to the created moment configuration file.
        """
        cfg_writer = MomentConfigWriter(config)
        moment_cfg = tempfile.NamedTemporaryFile(delete=False, mode="w")
        cfg_writer.write_config_to_file(moment_cfg.name)
        moment_cfg.close()

        return moment_cfg.name

    def create_template_config(self, output_path: str) -> None:
        """Create a template YAML configuration file.

        Args:
            output_path (str): Path where to save the template file.
        """
        template_config = PWAConfig(
            name="example_fit",
            description="Example PWA fit configuration",
            physics=PhysicsConfig(
                waveset=["1p", "1m"],
                ds_ratio="",
                phaselock=False,
                remove_waves=[],
            ),
            data=DataConfig(
                orientations=["PARA_0", "PERP_45"],
                run_periods=["allPeriods"],
                mass_min=1.0,
                mass_max=2.0,
                mass_width=0.1,
                t_bins=[0.1, 0.2, 0.4],
                energy=[8.2, 8.8],
                cut_recoil_pi_mass=1.4,
            ),
            compute=ComputeConfig(
                nrand=20,
                bootstrap=0,
                gpu=["1", "A100"],
                time_limit="02:00:00",
            ),
        )

        self.save_yaml_config(template_config, output_path)

    def validate_config(self, config: PWAConfig) -> List[str]:
        """Validate configuration and return list of errors.

        This function ensure no argument passed conflict with each other and produce
        undefined behavior. It also ensures that files passed do exist.

        Args:
            config (PWAConfig): Configuration to validate.

        Returns:
            List[str]: List of validation error messages. Empty if valid.
        """
        errors = []

        # ===PHYSICS===
        if not config.physics.waveset:
            errors.append("Physics: waveset cannot be empty")

        if config.physics.phase_reference:
            if (
                config.physics.phase_reference[0][0]
                == config.physics.phase_reference[1][0]
                and not config.physics.single_refl
            ):
                errors.append(
                    "Physics: phase_references must be in separate reflectivities"
                )
            if any(
                ref in config.physics.remove_waves
                for ref in config.physics.phase_reference
            ):
                errors.append(
                    "Physics: phase_reference waves cannot be in remove_waves"
                )

        # ===DATA===
        # ensure only one set of bin arguments are passed
        if config.data.mass_bins and config.data.mass_width != 0:
            errors.append(
                "Data: cannot specify custom mass bins and width. Please use bins, or"
                " specify min/max/width for automatic bin creation"
            )
        if config.data.t_bins and config.data.t_width != 0:
            errors.append(
                "Data: cannot specify custom t bins and width. Please use bins, or"
                " specify min/max/width for automatic bin creation"
            )

        if config.data.mass_bins:
            if len(config.data.mass_bins) < 2:
                errors.append("Data: custom mass_bins must have at least 2 values")
        elif (
            config.data.mass_min == config.data.mass_max
            or config.data.mass_width == 0
            or config.data.mass_min > config.data.mass_max
        ):
            errors.append(
                "Invalid mass bin parameters. Ensure mass_max > mass_min and mass_width"
                " is non-zero"
            )

        if config.data.t_bins:
            if len(config.data.t_bins) < 2:
                errors.append("Data: custom t_bins must have at least 2 values")
        elif (
            config.data.t_min == config.data.t_max
            or config.data.t_width == 0
            or config.data.t_min > config.data.t_max
        ):
            errors.append(
                "Invalid t bin parameters. Ensure t_max > t_min and t_width"
                " is non-zero"
            )

        # data
        if not os.path.exists(config.data.data_dir):
            errors.append(f"Data: data_dir does not exist: {config.data.data_dir}")
        else:
            data_files = os.listdir(config.data.data_dir)

            for ont in config.data.orientations:
                if all(ont not in fname for fname in data_files):
                    errors.append(
                        f"Orientation {ont} not found in directory:"
                        f" {config.data.data_dir}"
                    )

            data_pattern = f"{config.data.data_version}{config.data.data_option}"
            if all(data_pattern not in fname for fname in data_files):
                errors.append(
                    f"Data {data_pattern} does not exist in directory:"
                    f" {config.data.data_dir}"
                )

        # phasespace
        if not os.path.exists(config.data.phasespace_dir):
            errors.append(
                f"Data: phasespace_dir does not exist: {config.data.phasespace_dir}"
            )
        else:
            phasespace_files = os.listdir(config.data.phasespace_dir)
            phasespace_pattern = (
                f"{config.data.phasespace_version}{config.data.phasespace_option}"
            )
            if all(phasespace_pattern not in fname for fname in phasespace_files):
                errors.append(
                    f"Phasespace {phasespace_pattern} does not exist in directory:"
                    f" {config.data.phasespace_dir}"
                )

        # optional truth file
        if config.data.truth_file and not os.path.exists(config.data.truth_file):
            errors.append(
                f"Truth file {config.data.truth_file} does not exist,"
                " please ensure the absolute path is used"
            )

        # ===COMPUTATIONAL===
        if config.compute.nrand < 0:
            errors.append("Compute: nrand must be non-negative")

        if config.compute.bootstrap < 0:
            errors.append("Compute: bootstrap must be non-negative")

        # randomized parameters ruin the purpose of truth files,
        # and bootstrapping is unnecessary for them
        if config.data.truth_file and (
            config.compute.nrand != 0 or config.compute.bootstrap != 0
        ):
            errors.append(
                "Compute: nrand and bootstrap must be 0 when using truth files"
            )

        return errors

    def _expand_orientations(self, orientations: List[str]) -> List[str]:
        """Expand orientation shortcuts to full list.

        Args:
            orientations (List[str]): List of orientations, possibly including "ALL".

        Returns:
            List[str]: Expanded list of orientations.
        """
        if "ALL" in orientations:
            return ["PARA_0", "PERP_45", "PERP_90", "PARA_135"]
        return orientations

    def _yaml_dict_to_config(self, config_dict: Dict[str, Any]) -> PWAConfig:
        """Convert dictionary from YAML to PWAConfig object.

        Args:
            config_dict (Dict[str, Any]): Dictionary from YAML file.

        Returns:
            PWAConfig: Converted configuration object.
        """
        # Extract sections with defaults
        physics_dict = config_dict.get("physics", {})
        data_dict = config_dict.get("data", {})
        compute_dict = config_dict.get("compute", {})
        general_dict = config_dict.get("general", {})

        # Create config objects
        physics = PhysicsConfig(**physics_dict)
        data = DataConfig(**data_dict)
        compute = ComputeConfig(**compute_dict)
        general = GeneralConfig(**general_dict)

        # if "ALL" is used, expand to 4 orientation names
        data.orientations = self._expand_orientations(data.orientations)

        return PWAConfig(
            name=config_dict.get("name", ""),
            description=config_dict.get("description", ""),
            physics=physics,
            data=data,
            compute=compute,
            general=general,
        )

    def _config_to_yaml_dict(self, config: PWAConfig) -> Dict[str, Any]:
        """Convert PWAConfig object to dictionary for YAML serialization.

        Args:
            config (PWAConfig): Configuration object.

        Returns:
            Dict[str, Any]: Dictionary suitable for YAML serialization.
        """
        from dataclasses import asdict

        config_dict = {
            "name": config.name,
            "description": config.description,
            "physics": asdict(config.physics),
            "data": asdict(config.data),
            "compute": asdict(config.compute),
            "general": asdict(config.general),
        }

        return config_dict
