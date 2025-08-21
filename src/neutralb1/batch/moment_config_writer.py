"""Class for writing a moment cfg file for binned fits to be performed

Similar to the amptools_config_writer, this file contains the class that can write
a fit configuration file based off a user-requested config, including all necessary
information for the moment fits.
"""

import itertools
from dataclasses import dataclass
from typing import List, TextIO

from .config_models import PWAConfig


@dataclass(frozen=True)
class Moment:
    name: str
    alpha: int
    Jv: int
    Lambda: int
    J: int
    M: int


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


class MomentConfigWriter:
    def __init__(self, config: PWAConfig):
        self.config = config

    def write_config_to_file(self, file_name: str) -> None:
        """Write the AmpTools Moment fit cfg for a specific configuration

        Args:
            file_name (str): Path to the output configuration file.
        """

        with open(file_name, "w") as cfg_file:
            self._write_preamble(
                cfg_file, self.config.general.reaction, self.config.data.orientations
            )

            self._write_data_lines(
                cfg_file, self.config.data.orientations, self.config.general.reaction
            )

            moment_list = self._get_moment_list(
                self.config.physics.max_moment_J,
                self.config.physics.waveset,
            )

            self._write_parameters(cfg_file, moment_list, self.config.general.n_events)

            self._write_moments(
                cfg_file,
                moment_list,
                self.config.general.reaction,
                self.config.data.orientations,
            )

        return

    def _write_preamble(
        self, cfg_file: TextIO, reaction: str, orientations: List[str]
    ) -> None:
        """Write the first setup lines in the cfg

        Args:
            cfg_file (TextIO): cfg being written to
            reaction (str): user-defined reaction name
            orientations (List[str]): The list of polarization orientations.
        """
        cfg_file.write(f"\n\n{'#'*8} SETUP {'#'*8}\n")

        cfg_file.write(
            f"fit {reaction}_moment\n"  # avoid conflict with pwa fit names
            "keyword parRange 3 3\n"  # allows for random sampling of parameters
        )

        # create a normalization integral file for each reaction
        for ont in orientations:
            cfg_file.write(
                f"normintfile {reaction}_{POL_DICT[ont]['angle']}"
                f" {reaction}_{POL_DICT[ont]['angle']}_moment.ni\n"
            )

    def _write_data_lines(
        self, cfg_file: TextIO, orientations: List[str], reaction: str
    ) -> None:

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

        # write reaction, data reader, and sum lines that utilize the loops
        cfg_file.write(
            f"\n{'#'*8} DATA, REACTIONS, AND SUMS {'#'*8}\n"
            "reaction LOOPREAC Beam Proton Pi01 Pi02 Pi+ Pi-\n"
            "genmc LOOPREAC ROOTDataReader LOOPGENMC\n"
            "accmc LOOPREAC ROOTDataReader LOOPACCMC\n"
            "data LOOPREAC ROOTDataReader LOOPDATA\n\n"
            f"sum LOOPREAC vecPSMoment\n\n"
        )

        return

    def _get_moment_list(self, max_J: int, waveset: List[str]) -> List[Moment]:
        """Generate a set of moments based on the waveset configuration

        Assumes that the all moments below the maximum J value should be generated e.g.
        a waveset with just J=3 waves will produce moments with J = 0, 1, ... 6.

        Args:
            max_J (int): The maximum J value for the moments.
            waveset (List[str]): The list of waves requested.

        Raises:
            ValueError: If no moments are generated.

        Returns:
            List[Moment]: A list of generated moments.
        """

        # if max_J not set by user, then auto determine from waveset
        if max_J == 0:
            max_wave_J = max([int(JP[0]) for JP in waveset if JP != "iso"])
            max_J = 2 * max_wave_J  # moments go from 0 -> J_i + J_j (indexed by wave)

        alpha_array = [0, 1, 2]
        Jv_array = [0, 2]  # CGs always ensure Jv==1 moments are 0
        Lambda_array = [0, 1, 2]
        J_array = [x for x in range(max_J + 1)]
        M_array = [x for x in range(max_J + 1)]

        moment_list = []
        for alpha, Jv, Lambda, J, M in itertools.product(
            alpha_array, Jv_array, Lambda_array, J_array, M_array
        ):

            # Filter out non-physical moments
            if M > J or Lambda > J or Lambda > Jv:
                continue  # Wigner D functions are 0 in these conditions
            if alpha == 2 and Lambda == 0 and M == 0:
                continue  # Wishart section 5.10.2

            moment_list.append(
                Moment(
                    name=f"H{alpha}_{Jv}{Lambda}{J}{M}",
                    alpha=alpha,
                    Jv=Jv,
                    Lambda=Lambda,
                    J=J,
                    M=M,
                )
            )

        if not moment_list:
            raise ValueError("No moments generated, check your waveset configuration.")

        return moment_list

    def _write_parameters(
        self, cfg_file: TextIO, moment_list: List[Moment], n_events: int
    ) -> None:
        """Write the AmpTools parameters and parRange lines for each moment

        Args:
            cfg_file (TextIO): cfg being written to
            moment_list (List[Moment]): Moment parameters that will be fit with
            n_events (int): Number of events for setting moment initializations or
                random sample ranges
        """

        cfg_file.write(f"\n\n{'#'*8} MOMENTS {'#'*8}\n")

        # write parameter and parRange lines
        for mom in moment_list:
            if mom.name == "H0_0000":
                cfg_file.write(f"parameter {mom.name} {n_events}\n\n")
                continue  # avoid random sampling due to specific initialization

            # initialize to 0 if not random
            cfg_file.write(f"parameter {mom.name} 0.0\n")

            # randomly sample between +/- 1% of # of events
            cfg_file.write(f"parRange {mom.name} -{n_events*0.01} {n_events*0.01}\n\n")

        return

    def _write_moments(
        self,
        cfg_file: TextIO,
        moment_list: List[Moment],
        reaction: str,
        orientations: List[str],
    ) -> None:
        """Write the moment amplitudes and initializations

        Args:
            cfg_file (TextIO): cfg being written to
            moment_list (List[Moment]): Moment parameters that will be fit with
            reaction (str): Reaction name
            orientations (List[str]): List of polarizations
        """

        cfg_file.write(
            "# Define one moment amplitude for each reaction and all the same"
            " moment parameters\n"
            "# Fix the 'amplitude' to be real and 1, so there is no extra free"
            " parameter\n"
            "# Constrain lines are likely unnecessary, but included for extra assurance"
            " that amplitudes are constrained across orientations\n"
        )

        # structure of vector-pseudoscalar moment is:
        # <reaction>::<sum>::<name> Vec_ps_moment <angle> <fraction> [mom_1] [mom_2] ...
        for ont in orientations:
            cfg_file.write(
                f"amplitude {reaction}_{POL_DICT[ont]['angle']}::"
                "vecPSMoment::vecPSMoment Vec_ps_moment"
                f" {POL_DICT[ont]['angle']} {POL_DICT[ont]['fraction']}"
            )
            for mom in moment_list:
                cfg_file.write(f" [{mom.name}]")
            cfg_file.write("\n")

            # moment parameters do all the work, the amplitude production coefficient
            # is not useful so is fixed to 1
            cfg_file.write(
                f"initialize {reaction}_{POL_DICT[ont]['angle']}::"
                "vecPSMoment::vecPSMoment cartesian 1 0 fixed\n"
            )
            cfg_file.write(
                f"scale {reaction}_{POL_DICT[ont]['angle']}::"
                f"vecPSMoment::vecPSMoment [parScale_{POL_DICT[ont]['angle']}]\n"
            )
            cfg_file.write(
                f"constrain"
                f" {reaction}_{POL_DICT[orientations[0]]['angle']}::"
                f"vecPSMoment::vecPSMoment"
                f" {reaction}_{POL_DICT[ont]['angle']}::"
                f"vecPSMoment::vecPSMoment\n\n"
            )

        return
