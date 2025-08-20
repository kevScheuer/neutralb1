"""Class for writing an AmpTools fit cfg file for binned fits to be performed

Its role is to write a fit cfg file that will be run in each bin of the fit. It uses the
user requested waveset and all its modifications to write the necessary AmpTools
recognized lines to the cfg file.

NOTE: This file (and other batch files) are tuned specifically to omegapi0
"""

from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from typing import List, TextIO

from .config_models import PWAConfig

# CONSTANTS
# orientations have to be looped over manually here (without AmpTools LOOPs),
#   so this sets up necessary variables. Polarization angles are in degrees and
#   fractions denote the fraction of the intensity that is polarized in that direction,
#   whose values are standard in halld
POL_DICT = {
    "PARA_0": {"angle": 0, "fraction": 0.3519},
    "PERP_45": {"angle": 45, "fraction": 0.3374},
    "PERP_90": {"angle": 90, "fraction": 0.3303},
    "PARA_135": {"angle": 135, "fraction": 0.3375},
}
REFLECTIVITY_DICT = {
    "ImagNegSign": -1,
    "RealPosSign": -1,
    "ImagPosSign": 1,
    "RealNegSign": 1,
}


@dataclass(frozen=True)
class Wave:
    name: str
    reflectivity: int
    spin: int
    parity: int
    m: int
    l: int


class AmpToolsConfigWriter:
    def __init__(self, config: PWAConfig):
        self.config = config

    def write_config_to_file(self, file_name: str) -> None:
        """Write the AmpTools PWA fit cfg file for a specific configuration in the class

        Args:
            file_name (str): Path to the output configuration file.
        """

        # convert the user requested wave string into lists of waves/bw's that contain all
        # the necessary info for each wave/bw
        waves, breit_wigners = self._get_waves_and_breit_wigners(
            self.config.physics.waveset,
            self.config.physics.remove_waves,
            self.config.physics.single_refl,
        )

        # parse out special options aside from the waves and bw's
        is_iso = True if "iso" in self.config.physics.waveset else False
        is_dalitz = False if "nodalitz" in self.config.physics.waveset else True

        # copy common part of all omegapi config files
        with open(self.config.general.template_name, "r") as f:
            template_data = f.read()
        template_data = template_data.replace("FITNAME", self.config.general.reaction)

        with open(file_name, "w") as cfg_file:
            # write the common part of the cfg file
            cfg_file.write(template_data)

            # Write lines to load in data files for each orientation
            self._write_data_lines(
                cfg_file,
                self.config.data.orientations,
                self.config.general.reaction,
                is_iso,
            )

            # Write all the lines for each wave in the waveset
            self._write_waves(
                cfg_file,
                waves,
                self.config.physics.init_refl,
                self.config.physics.init_real,
                self.config.physics.init_imag,
                is_dalitz,
                self.config.physics.frame,
                self.config.data.orientations,
                self.config.general.reaction,
            )

            # write lines if isotropic background is present
            if is_iso:
                self._write_isotropic_background(
                    cfg_file,
                    is_dalitz,
                    self.config.data.orientations,
                    self.config.general.reaction,
                )

            # check that the reference waves are in the waveset if requested, otherwise
            # set the reference to be the first wave in the waveset that has both
            # reflectivities
            phase_reference = []
            if self.config.physics.phase_reference:
                for wave in waves:
                    if wave.name.lower() in [
                        ref.lower() for ref in self.config.physics.phase_reference
                    ]:
                        phase_reference.append(wave)
                if len(phase_reference) != 2:
                    raise ValueError(
                        "The phase reference waves"
                        f" {self.config.physics.phase_reference} were not"
                        " found in the requested waveset"
                    )
            else:
                for wave1, wave2 in product(waves, waves):
                    if (
                        wave1.reflectivity == -1
                        and wave2.reflectivity == 1
                        and wave1.spin == wave2.spin
                        and wave1.parity == wave2.parity
                        and wave1.m == wave2.m
                        and wave1.l == wave2.l
                    ):
                        phase_reference.append(wave1)
                        phase_reference.append(wave2)
                        break
            if not phase_reference:
                raise ValueError(
                    "Phase reference wave was not specified, but could not be"
                    " automatically determined as no wave has both reflectivities"
                )

            # phaselock function needs its own reference wave method
            if not self.config.physics.phaselock:
                # write lines to remove overall unconstrained phase factor
                self._write_phase_convention(
                    cfg_file,
                    phase_reference,
                    self.config.data.orientations,
                    self.config.general.reaction,
                )

            # writes breit wigners (if present in waveset)
            if breit_wigners:
                self._write_breit_wigners(
                    cfg_file,
                    waves,
                    breit_wigners,
                    self.config.data.orientations,
                    self.config.general.reaction,
                )
            # write ds ratio if 1p amp is used, and if DS ratio isn't free
            if self.config.physics.ds_ratio != "free" and any(
                wave.spin == 1 and wave.parity == 1 for wave in waves
            ):
                self._write_ds_ratio(
                    cfg_file,
                    waves,
                    self.config.physics.ds_ratio,
                    self.config.physics.phaselock,
                    self.config.data.orientations,
                    self.config.general.reaction,
                )

            # user option to constrain phases between m projections
            if self.config.physics.phaselock:
                self._write_phaselock(
                    cfg_file,
                    waves,
                    phase_reference,
                    self.config.physics.single_refl,
                    self.config.physics.init_refl,
                    self.config.physics.init_real,
                    self.config.data.orientations,
                    self.config.general.reaction,
                )

        return

    def _get_waves_and_breit_wigners(
        self,
        requested_jp_and_bw: List[str],
        waves_to_remove: List[str],
        single_refl: int,
    ) -> tuple[List[Wave], List[dict]]:
        """Converts list of JP waves and breit wigners into code interpretable lists

        This function takes JP (spin+parity) quantum numbers and converts them into
        a list of Wave dataclasses that contain all necessary info for each wave. All the
        possible waves of each JP value are used, unless the user requests to remove them.
        If breit wigners are requested by name, they are added to a separate list of dicts
        that contain all necessary info for each breit wigner.

        Args:
            requested_jp_and_bw (List[str]): user requested jp waves and breit wigners to be
                built in cfg file e.g. ["1p", "1m", "b1"]
            waves_to_remove (List[str]): user passed waves to remove from waveset in eJPmL
                format.
            single_refl (int): when non-zero, only waves with the chosen reflectivity will
                be written
        Returns:
            Two lists. The first is a list of Wave dataclasses that contain all necessary
            info for each wave requested. The second is a list of dicts that contain all
            necessary info for each breit wigner requested.
        """

        # define wave dictionary that contains all possible quantum numbers for each wave
        wave_dict = {
            "0m": {
                "reflectivity": [-1, 1],
                "spin": 0,
                "parity": -1,
                "m": [0],
                "l": [1],
            },
            "1p": {
                "reflectivity": [-1, 1],
                "spin": 1,
                "parity": +1,
                "m": [-1, 0, 1],
                "l": [0, 2],
            },
            "1m": {
                "reflectivity": [-1, 1],
                "spin": 1,
                "parity": -1,
                "m": [-1, 0, 1],
                "l": [1],
            },
            "2p": {
                "reflectivity": [-1, 1],
                "spin": 2,
                "parity": +1,
                "m": [-2, -1, 0, 1, 2],
                "l": [2],
            },
            "2m": {
                "reflectivity": [-1, 1],
                "spin": 2,
                "parity": -1,
                "m": [-2, -1, 0, 1, 2],
                "l": [1, 3],
            },
            "3m": {
                "reflectivity": [-1, 1],
                "spin": 3,
                "parity": -1,
                "m": [-3, -2, -1, 0, 1, 2, 3],
                "l": [3],
            },
        }

        # Create a list of all possible waves using the Wave dataclass
        waveset_all = [
            Wave(
                name=(
                    f"{self._int_to_char(refl)}{spin}{self._int_to_char(parity)}"
                    f"{self._int_to_char(m)}{self._int_to_char(l, True)}"
                ),
                reflectivity=refl,
                spin=spin,
                parity=parity,
                m=m,
                l=l,
            )
            for val_dict in wave_dict.values()
            for refl, spin, parity, m, l in product(
                val_dict["reflectivity"],
                [val_dict["spin"]],
                [val_dict["parity"]],
                val_dict["m"],
                val_dict["l"],
            )
        ]
        # Check if any requested wave to remove is not in the waveset
        all_wave_names = [wave.name.lower() for wave in waveset_all]
        for rm_wave in waves_to_remove:
            if rm_wave.lower() not in all_wave_names:
                raise ValueError(
                    f"The {rm_wave} wave that was requested to be removed is not a valid"
                    " wave"
                )

        # Add requested waves, if not to be removed or do not match requested reflectivity
        waveset = [
            wave
            for wave in waveset_all
            if wave.name.lower()[1:3] in requested_jp_and_bw
            and wave.name.lower()
            not in [rm_wave.lower() for rm_wave in waves_to_remove]
            and (single_refl == 0 or wave.reflectivity == single_refl)
        ]

        breit_wigner_dict = {
            "b1": {
                "jp": "1p",
                "spin": +1,
                "mass": 1.235,
                "mass_low": 1.1,
                "mass_high": 1.35,
                "width": 0.142,
                "width_low": 0.1,
                "width_high": 0.3,
                "indices": "2 345",
                "fix": True,
            },
            "rho": {
                "jp": "1m",
                "spin": +1,
                "mass": 1.465,
                "mass_low": 1.3,
                "mass_high": 1.6,
                "width": 0.4,
                "width_low": 0.2,
                "width_high": 0.6,
                "indices": "2 345",
                "fix": True,
            },
        }
        if "b1free" in requested_jp_and_bw:
            breit_wigner_dict["b1"]["fix"] = False
        if "rhofree" in requested_jp_and_bw:
            breit_wigner_dict["rho"]["fix"] = False

        # Add requested breit wigner dictionaries to a list
        breit_wigners = []
        for key, val in breit_wigner_dict.items():
            if key in requested_jp_and_bw:
                breit_wigners.append(val)

        # check that amplitudes are present for any requested breit wigners
        if breit_wigners:
            jp_set = set()
            for wave in waveset:
                jp_set.add(wave.name[1:3])
            for bw in breit_wigners:
                if bw["jp"] not in jp_set:
                    raise ValueError(
                        f"Breit Wigner JP={bw['jp']} not found in requested waveset"
                    )

        if waveset == []:
            raise ValueError("No waves found for requested waveset")

        return waveset, breit_wigners

    def _write_data_lines(
        self, cfg_file: TextIO, orientations: List[str], reaction: str, is_iso: bool
    ) -> None:
        """Writes the AmpTools 'loop' statements to read in data files for each orientation

        This function also handles the 'reaction' and 'sum' lines

        Args:
            cfg_file (TextIO): cfg file being written to
            orientations (List[str]): diamond orientation settings
            reaction (str): user-defined reaction name
            is_iso (bool): if set, isotropic background is present
        """

        cfg_file.write(f"\n\n{'#'*8} LOOP OVER POLARIZATION STATES {'#'*8}\n")

        # create the parScale parameters for each orientation
        for i, ont in enumerate(orientations):
            pol_scales = {
                "PARA_0": "1.0",  # TODO: potential for an arg that changes these values
                "PERP_45": "1.0",
                "PERP_90": "1.0",
                "PARA_135": "1.0",
            }
            fix_string = "fixed" if i == 0 else ""
            cfg_file.write(
                f"parameter parScale_{POL_DICT[ont]['angle']}"
                f" {pol_scales[ont]} {fix_string}\n"
            )

        # create reaction, scale, and data loops for each orientation
        reaction_and_angles = ""
        loop_data_str = ""
        for ont in orientations:
            reaction_and_angles += f" {reaction}_{POL_DICT[ont]['angle']}"
            loop_data_str += f" anglesOmegaPiAmplitude_{POL_DICT[ont]['angle']}.root"

        cfg_file.write(f"\nloop LOOPREAC{reaction_and_angles}\n")

        # the phasespace files are the same for each orientation
        cfg_file.write(
            f"loop LOOPGENMC {'anglesOmegaPiPhaseSpace.root ' * len(orientations)} \n"
            f"loop LOOPACCMC {'anglesOmegaPiPhaseSpaceAcc.root ' * len(orientations)} \n"
            f"loop LOOPDATA{loop_data_str}\n"
        )

        bkgd = "Bkgd" if is_iso else ""

        # write reaction, data reader, and sum lines that utilize the loops
        cfg_file.write(
            f"\n{'#'*8} DATA, REACTIONS, AND SUMS {'#'*8}\n"
            "reaction LOOPREAC Beam Proton Pi01 Pi02 Pi+ Pi-\n"
            "genmc LOOPREAC ROOTDataReader LOOPGENMC\n"
            "accmc LOOPREAC ROOTDataReader LOOPACCMC\n"
            "data LOOPREAC ROOTDataReader LOOPDATA\n\n"
            f"sum LOOPREAC ImagNegSign RealNegSign RealPosSign ImagPosSign {bkgd}\n\n"
        )

        return

    ####################### WRITE AMPLITUDES #######################
    def _write_waves(
        self,
        cfg_file: TextIO,
        waves: List[Wave],
        init_refl: int,
        init_real: float,
        init_imag: float,
        is_dalitz: bool,
        frame: str,
        orientations: List[str],
        reaction: str,
    ) -> None:
        """Write all waves to the cfg file for vecps_refl amplitudes

        Each wave and its necessary lines are each written out individually to the cfg
        file. This method is preferred over AmpTools loops, so we can manipulate individual
        waves much easier. The amplitude lines are written in the following order:
        - Comment the wave info
        - Write the vecps_refl amplitude line
        - Write the OmegaDalitz amplitude line (if not removed)
        - Write the initialization line
        - Write the parScale line
        - Repeat all of this for each orientation
        - Repeat for the other reflectivity sum

        Args:
            cfg_file (TextIO): cfg file being written to
            waves (List[Wave]): list of Wave dataclasses that contain all necessary info
                See 'get_waves_and_breit_wigners' for details
            init_refl (int): if non-zero, the real and imaginary parts of the amplitudes
                that don't match the chosen reflectivity will be initialized to 0
            init_real (float): initialization value for real part of amplitudes.
            init_imag (float): initialization value for imaginary part of amplitudes.
            is_dalitz (bool): controls whether OmegaDalitz amplitudes are written
            frame (str): decay frame used to define the angles
            orientations (List[str]): diamond orientation settings
            reaction (str): user-defined reaction name
        """

        cfg_file.write(f"\n\n{'#'*8} AMPLITUDES {'#'*8}\n")

        # write the amplitude lines for each wave in the waveset
        for wave in waves:
            cfg_file.write(f"\n#{'-'*4} {wave.name} {'-'*4}\n")

            # each reflectivity has 2 coherent sums, so loop over them
            for sum_label, refl_sign in REFLECTIVITY_DICT.items():
                if refl_sign != wave.reflectivity:
                    continue

                real = 1 if "Real" in sum_label else -1
                sign = 1 if "Pos" in sum_label else -1

                # last loop is over the polarization orientations
                for ont in orientations:
                    cfg_file.write(
                        f"\n{'#'*2} Sum = {sum_label}, Orientation = {ont} {'#'*2}\n"
                    )

                    # common amplitude string used in many lines
                    amplitude = (
                        f"{reaction}_{POL_DICT[ont]['angle']}::{sum_label}::{wave.name}"
                    )

                    # Write amplitude line NOTE: 'omega3pi' is specific to omega->3pi decays
                    cfg_file.write(
                        f"amplitude {amplitude} Vec_ps_refl {wave.spin} {wave.m} {wave.l}"
                        f" {real} {sign}"
                        f" {POL_DICT[ont]['angle']} {POL_DICT[ont]['fraction']}"
                        f" omega3pi {frame}\n"
                    )
                    # write dalitz line ('dalitz' has defined params in template file)
                    if is_dalitz:
                        cfg_file.write(f"amplitude {amplitude} OmegaDalitz dalitz\n")

                    # initialize the amplitude
                    cfg_file.write(
                        "initialize {amp} cartesian {re} {im}\n".format(
                            amp=amplitude,
                            re=(
                                0
                                if init_refl and wave.reflectivity != init_refl
                                else init_real
                            ),
                            im=(
                                0
                                if init_refl and wave.reflectivity != init_refl
                                else init_imag
                            ),
                        )
                    )

                    # attach polarization scale parameter to amplitude
                    cfg_file.write(
                        f"scale {amplitude} [parScale_{POL_DICT[ont]['angle']}]\n"
                    )

                    # constrain coherent sums to the first orientation. This ensures
                    # that different orientations are manipulating the same amplitudes
                    cfg_file.write(
                        f"constrain"
                        f" {reaction}_{POL_DICT[orientations[0]]['angle']}::"
                        f"{sum_label}::{wave.name}"
                        f" {reaction}_{POL_DICT[ont]['angle']}::"
                        f"{sum_label}::{wave.name}\n"
                    )

            # constrain coherent sums across reflectivities. Each reflectivity has two sums
            # so they must be constrained to each other. This only constrains the first
            # orientations to each other, since all other orientations are constrained to
            # the first orientation
            cfg_file.write("\nconstrain ")
            for sum_label, refl_sign in REFLECTIVITY_DICT.items():
                if refl_sign != wave.reflectivity:
                    continue
                cfg_file.write(
                    f"{reaction}_{POL_DICT[orientations[0]]['angle']}::"
                    f"{sum_label}::{wave.name} "
                )
            cfg_file.write("\n")

        return

    def _write_isotropic_background(
        self, cfg_file: TextIO, is_dalitz: bool, orientations: List[str], reaction: str
    ) -> None:
        """Writes lines in cfg file template for isotropic background amplitude

        Args:
            cfg_file (TextIO): cfg file being written to
            is_dalitz (bool): controls whether OmegaDalitz amplitudes are written
            orientations (List[str]): diamond orientation settings
            reaction (str): user-defined reaction name
        """
        cfg_file.write(f"\n{'#'*8} isotropic background {'#'*8}")

        for ont in orientations:
            if is_dalitz:
                dalitz = (
                    f"amplitude {reaction}_{POL_DICT[ont]['angle']}"
                    "::Bkgd::isotropic OmegaDalitz dalitz"
                )
            else:
                dalitz = ""

            cfg_file.write(
                f"\namplitude {reaction}_{POL_DICT[ont]['angle']}::Bkgd::isotropic Uniform"
                f"\n{dalitz}"
                f"\ninitialize {reaction}_{POL_DICT[ont]['angle']}"
                "::Bkgd::isotropic cartesian 100 0 real"
                f"\nscale {reaction}_{POL_DICT[ont]['angle']}::Bkgd::isotropic"
                f" [parScale_{POL_DICT[ont]['angle']}]"
                f"\nconstrain"
                f" {reaction}_{POL_DICT[orientations[0]]['angle']}::"
                f"Bkgd::isotropic"
                f" {reaction}_{POL_DICT[ont]['angle']}::"
                f"Bkgd::isotropic\n"
            )

        return

    def _write_phase_convention(
        self,
        cfg_file: TextIO,
        waves: List[Wave],
        orientations: List[str],
        reaction: str,
    ) -> None:
        """Constrains 1 amplitude in each reflectivity to be purely real

        Due to the construction of the intensity formulation, each of the 2 unique
        coherent sums can have an overall phase factor that goes unaccounted for.
        To remedy this, one of the amplitudes is set to be purely real, so that the
        phases of every other amplitude are defined relative to this one, and the
        overall unconstrained phase is removed.

        Args:
            cfg_file (TextIO): cfg file being written to
            waves (List[Wave]): the two Wave objects that will be set to be purely real.
            orientations (List[str]): diamond orientation settings
            reaction (str): user-defined reaction name
        """
        cfg_file.write("\n# fix phase\n")

        for ont in orientations:
            for sum_label, refl_sign in REFLECTIVITY_DICT.items():
                for wave in waves:
                    if wave.reflectivity != refl_sign:
                        continue
                    cfg_file.write(
                        f"initialize {reaction}_{POL_DICT[ont]['angle']}"
                        f"::{sum_label}::{wave.name} cartesian 100 0 real\n"
                    )

        return

    def _write_breit_wigners(
        self,
        cfg_file: TextIO,
        waves: List[Wave],
        breit_wigners: List[dict],
        orientations: List[str],
        reaction: str,
    ) -> None:
        """Write breit wigner lines to the cfg file

        Args:
            cfg_file (TextIO): cfg file being written to
            waves (List[Wave]): list of Wave dataclasses that contain all necessary info.
                See 'get_waves_and_breit_wigners' for details
            breit_wigners (List[dict]): list of bw dicts that contain all necessary info.
                See 'get_waves_and_breit_wigners' for details
            orientations (List[str]): diamond orientation settings
            reaction (str): user-defined reaction name
        """

        cfg_file.write(f"\n{'#'*8} Mass dependent amplitudes {'#'*8}\n")

        # write mass and width parameters
        for bw in breit_wigners:

            # fix breit wigner parameters if requested, or bound their allowed range
            mass_option = (
                "fixed"
                if bw["fix"]
                else (
                    f"bounded"
                    f" {bw['mass_low']} {bw['mass_high']}"
                    f"\nparRange {bw['jp']}_mass"
                    f" {bw['mass_low']} {bw['mass_high']}"
                )
            )
            width_option = (
                "fixed"
                if bw["fix"]
                else (
                    f"bounded"
                    f" {bw['width_low']} {bw['width_high']}"
                    f"\nparRange {bw['jp']}_width"
                    f" {bw['width_low']} {bw['width_high']}"
                )
            )
            cfg_file.write(f"parameter {bw['jp']}_mass {bw['mass']} {mass_option}\n")
            cfg_file.write(f"parameter {bw['jp']}_width {bw['width']} {width_option}\n")

        for wave in waves:
            for bw in breit_wigners:
                # only attach the breit wigner to the waves with matching jp values
                if bw["jp"] != wave.name[1:3]:
                    continue
                for sum_label, refl_sign in REFLECTIVITY_DICT.items():
                    if refl_sign != wave.reflectivity:
                        continue
                    for ont in orientations:
                        cfg_file.write(
                            f"amplitude {reaction}_{POL_DICT[ont]['angle']}"
                            f"::{sum_label}::{wave.name}"
                            f" BreitWigner [{bw['jp']}_mass] [{bw['jp']}_width]"
                            f" {bw['spin']} {bw['indices']}\n"
                        )

        return

    def _write_ds_ratio(
        self,
        cfg_file: TextIO,
        waves: List[Wave],
        ds_option: str,
        is_phaselock: bool,
        orientations: List[str],
        reaction: str,
    ) -> None:
        """Writes lines in cfg file template for a D to S wave ratio to be set

        Initialized (and optionally fixed) values are E852 parameters. The ratio can be
        fixed to E852 values, split between two reflectivities, or left to float between 0
        and 1.

        Args:
            cfg_file (TextIO): cfg file being written to
            waves (List[Wave]): list of Wave dataclasses that contain all necessary info
                See 'get_waves_and_breit_wigners' for details
            ds_option (str): user option to fix the ratio to E852 values, split it between,
                reflectivities, or let float between 0 and 1
            is_phaselock (bool): user option to constrain phases between eJPl amplitudes.
                In this function it fixes the "dsphase" parameter to 0
            orientations (List[str]): diamond orientation settings
            reaction (str): user-defined reaction name
        """

        cfg_file.write(
            "\n# constrain S and D waves to same amplitude and set scale factor"
            " for D/S ratio"
        )

        # setup whether the ratio and phase will be fixed values, or bounded within a range
        ratio_opt = "fixed" if ds_option == "fixed" else "bounded 0.0 1.0"
        phase_opt = (
            "fixed"
            if ds_option == "fixed" or is_phaselock
            else "bounded -3.14159 3.14159"
        )
        phase_val = "0" if is_phaselock else "0.184"

        # setup the parameter names for AmpTools, and split between reflectivities if
        # requested
        if ds_option == "split":
            ratios = ["dsratio_p", "dsratio_m"]
            phases = ["dphase_p", "dphase_m"]
        else:
            ratios = ["dsratio"]
            phases = ["dphase"]

        for ont in orientations:
            for ratio, phase in zip(ratios, phases):
                cfg_file.write(
                    f"\nparameter {ratio} 0.27 {ratio_opt}\n"
                    f"parameter {phase} {phase_val} {phase_opt}\n"
                )
                if ds_option != "fixed":
                    cfg_file.write(f"parRange {ratio} 0.0 1.0\n")
                if not is_phaselock and ds_option != "fixed":
                    cfg_file.write(f"parRange {phase} -3.14159 3.14159\n")

                for sum_label, refl_sign in REFLECTIVITY_DICT.items():
                    for wave in waves:
                        if wave.reflectivity != refl_sign or (
                            wave.spin != 1 or wave.parity != 1 or wave.l != 2
                        ):
                            continue

                        # write out the D wave amplitude line
                        cfg_file.write(
                            f"amplitude {reaction}_{POL_DICT[ont]['angle']}::{sum_label}::"
                            f"{wave.name} ComplexCoeff [{ratio}] [{phase}] MagPhi\n"
                        )
                        # constrain the D wave to the S wave
                        cfg_file.write(
                            f"constrain {reaction}_{POL_DICT[ont]['angle']}::{sum_label}::"
                            f"{wave.name} {reaction}_{POL_DICT[ont]['angle']}::"
                            f"{sum_label}::{wave.name.replace("D", "S")}\n"
                        )

        return

    def _write_phaselock(
        self,
        cfg_file: TextIO,
        waves: List[Wave],
        reference_waves: List[str],
        single_refl: int,
        init_refl: int,
        init_real: float,
        orientations: List[str],
        reaction: str,
    ) -> None:
        """Lock m-projections to have the same phase

        Phaselocking refers to the model that the m-projections for a certain JPL
        are all given the same phase value, meaning no phase difference can occur
        within an eJPL combination

        Args:
            cfg_file (TextIO): cfg file being written to
            waves (List[Wave]): list of Wave dataclasses that contain all necessary info
                See 'get_waves_and_breit_wigners' for details
            reference_waves (List[str]): All m projections that match jpl combinations
                in the list will be fixed to be real.
            single_refl (int): if non-zero, only the chosen reflectivity will be written
            init_refl (int): if non-zero, the real and imaginary parts of the amplitudes
                that don't match the chosen reflectivity will be initialized to 0
            init_real (float): initialization value for real part of amplitudes.
            orientations (list): diamond orientation settings
            reaction (str): user-defined reaction name
        """

        cfg_file.write(f"\n{'#'*8} PHASELOCK {'#'*8}")

        # create a dictionary of all waves with the same jpl values
        jpl_sets = defaultdict(set)
        for wave in waves:
            jp = wave.name[1:3]
            l = wave.name[-1]
            jpl_sets[jp + l].add(wave)

        # write a common phase parameter for each jpl value to be shared among m projections
        for jpl, wave_group in jpl_sets.items():
            cfg_file.write(
                f"\n\n{'#'*2} set the phase for all m projections of {jpl}{'#'*2}\n"
            )

            # write the parameter for each reflectivity once
            for refl_sign in set(REFLECTIVITY_DICT.values()):
                if single_refl and refl_sign != single_refl:
                    continue
                phaselock_string = f"phase_{self._int_to_char(refl_sign)}{jpl}"
                if [r[:2] + r[-1] == jpl for r in reference_waves]:
                    cfg_file.write(f"parameter {phaselock_string} 0 fixed\n")
                else:
                    cfg_file.write(
                        f"parameter {phaselock_string} 0 bounded -3.14159 3.14159\n"
                        f"parRange {phaselock_string} -3.14159 3.14159\n"
                    )

            # write the amplitude lines for each wave in the waveset
            for coh_sum_label, refl_sign in REFLECTIVITY_DICT.items():
                phaselock_string = f"phase_{self._int_to_char(refl_sign)}{jpl}"

                for wave in wave_group:
                    if wave.reflectivity != refl_sign:
                        continue
                    for ont in orientations:
                        amplitude = (
                            f"{reaction}_{POL_DICT[ont]['angle']}"
                            f"::{coh_sum_label}::{wave.name}"
                        )
                        # re-initialize the amplitude to be purely real
                        cfg_file.write(
                            f"initialize {amplitude}"
                            " cartesian {re} 0 real\n".format(
                                re=(
                                    0
                                    if init_refl and refl_sign != init_refl
                                    else init_real
                                )
                            )
                        )
                        # attach the common phase parameter to the amplitude
                        cfg_file.write(
                            f"amplitude {amplitude} PhaseOffset [{phaselock_string}]\n"
                        )

        return

    def _int_to_char(self, x: int, is_ang_momentum: bool = False) -> str:
        """Convert quantum # integers to character amplitude notation

        i.e. m=-1 -> m, m=+1 -> p

        Args:
            x: integer to be converted
            is_ang_momentum: if true, uses a separate dictionary to perform the
                conversion. Uses the standard "S,P,D,F,G" wave notation
        Returns:
            single character
        """
        # TODO: this letter map is terrible, and should really use 2 characters to denote
        # the m-projection. This would require a pretty massive rewrite though since the
        # single character "eJPmL" format is baked into many scripts
        map = {-3: "l", -2: "n", -1: "m", 0: "0", +1: "p", +2: "q", +3: "r"}
        if is_ang_momentum:
            map = {0: "S", 1: "P", 2: "D", 3: "F", 4: "G"}

        return map[x]
