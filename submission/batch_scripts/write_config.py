"""Writes a template cfg file for binned fits to be performed

This script is part of the batch submission process, called with submit.py. Its
role is to write a cfg file that acts as a template for all the binned fits
to overwrite. This file accepts a users waveset, along with any modifications to
it, and writes the necessary AmpTools cfg file lines to fit using that waveset

NOTE: This file (and other batch files) are tuned specifically to omegapi0
TODO:
    Make refl_dict a constant and handle appropriately in each function
    Have phaselock function handle dominant eJPmL arg

TODO: finish converting all the methods to instead write out each amplitude line-by-line
    instead of the AmpTools LOOP method. Also remember that the naming scheme is changed
    such that reflectivity is explicit, and no longer needs to be inferred from
    the "ImagPosSign" type strings
"""

from dataclasses import dataclass
from itertools import product
from typing import List, TextIO

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


@dataclass
class Wave:
    name: str
    reflectivity: int
    spin: int
    parity: int
    m: int
    l: int


def main(args: dict):
    """Write AmpTools config file using the parsed args from submit.py"""

    # check for conflicting arguments
    jpml_set = {wave_string[1:] for wave_string in args["remove_waves"]}
    if args["phase_reference"] in jpml_set:
        raise ValueError(
            f"User requested to remove the reference wave {args['phase_reference']} from"
            " the waveset"
        )

    # convert the user requested wave string into lists of waves/bw's that contain all
    # the necessary info for each wave/bw
    waves, breit_wigners = get_waves_and_breit_wigners(
        args["waveset"], args["remove_waves"]
    )

    # parse out special options aside from the waves and bw's
    is_iso = True if "iso" in args["waveset"] else False
    is_dalitz = False if "nodalitz" in args["waveset"] else True

    # copy common part of all omegapi config files
    with open(f"submission/batch_scripts/{args['template_name']}", "r") as f:
        template_data = f.read()
    template_data = template_data.replace("FITNAME", args["reaction"])

    with open(f"submission/batch_scripts/fit.cfg", "w") as cfg_file:
        # write the common part of the cfg file
        cfg_file.write(template_data)

        # Write lines to load in data files for each orientation
        write_data_lines(cfg_file, args["orientations"], args["reaction"], is_iso)

        # Write all the lines for each wave in the waveset
        write_waves(
            cfg_file,
            waves,
            args["force_refl"],
            args["init_refl"],
            args["init_real"],
            args["init_imag"],
            is_dalitz,
            args["frame"],
            args["orientations"],
            args["reaction"],
        )

        # write lines if isotropic background is present
        if is_iso:
            write_isotropic_background(
                cfg_file, is_dalitz, args["orientations"], args["reaction"]
            )

        # check that the reference wave is in the waveset if requested, otherwise set
        # the reference to be the first wave in the waveset
        if args["phase_reference"]:
            if not is_reference_wave_in_waveset(args["phase_reference"], waves):
                raise ValueError(f"{args['phase_reference']} is not in waveset")
            phase_reference = args["phase_reference"]
        else:
            phase_reference = waves[0].name

        # write lines to remove overall unconstrained phase factor
        write_phase_convention(
            cfg_file,
            phase_reference,
            args["force_refl"],
            args["orientations"],
            args["reaction"],
        )

        # writes breit wigners (if present in waveset)
        if breit_wigners:
            write_breit_wigners(
                cfg_file,
                waves,
                breit_wigners,
                args["force_refl"],
                args["orientations"],
                args["reaction"],
            )
        exit()
        # write ds ratio if 1p amp is used, and if DS ratio isn't free
        is_1p_present = False
        for wave in waves:
            if wave.get("spin") == 1 and wave.get("parity") == 1:
                is_1p_present = True
        if args["ds_ratio"] != "free" and is_1p_present:
            write_ds_ratio(
                cfg_file,
                args["ds_ratio"],
                args["force_refl"],
                args["phaselock"],
                args["orientations"],
                args["reaction"],
            )

        # user option to constrain phases between m projections
        if args["phaselock"]:
            write_phaselock(
                cfg_file,
                waves,
                args["force_refl"],
                args["init_refl"],
                args["init_real"],
                args["orientations"],
                args["reaction"],
            )

    return


def get_waves_and_breit_wigners(
    requested_jp_and_bw: List[str], waves_to_remove: List[str]
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
    Returns:
        Two lists. The first is a list of Wave dataclasses that contain all necessary
        info for each wave requested. The second is a list of dicts that contain all
        necessary info for each breit wigner requested.
    """

    # define wave dictionary that contains all possible quantum numbers for each wave
    wave_dict = {
        "0m": {"reflectivity": [-1, 1], "spin": 0, "parity": -1, "m": [0], "l": [1]},
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
                f"{int_to_char(refl)}{spin}{int_to_char(parity)}"
                f"{int_to_char(m)}{int_to_char(l, True)}"
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

    # Remove waves if requested, and add requested waves to final waveset
    waveset = [
        wave
        for wave in waveset_all
        if wave.name.lower()[1:3] in requested_jp_and_bw
        and wave.name.lower() not in [rm_wave.lower() for rm_wave in waves_to_remove]
    ]

    breit_wigner_dict = {
        "b1": {
            "jp": "1p",
            "spin": +1,
            "mass": 1.235,
            "mass_low": 1.15,
            "mass_high": 1.15,
            "width": 0.142,
            "width_low": 0.1,
            "width_high": 0.2,
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


def write_data_lines(
    cfg_file: TextIO, orientations: List[str], reaction: str, is_iso: bool
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
    parScale_and_angles = ""
    loop_data_str = ""
    for ont in orientations:
        reaction_and_angles += f" {reaction}_{POL_DICT[ont]['angle']}"
        parScale_and_angles += f" [parScale_{POL_DICT[ont]['angle']}]"
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
def write_waves(
    cfg_file: TextIO,
    waves: List[Wave],
    force_refl: int,
    init_refl: int,
    init_real: float,
    init_imag: float,
    is_dalitz: bool,
    frame: str,
    orientations: List[str],
    reaction: str,
) -> None:
    """Write all waves to the cfg file for vecps_refl amplitudes

    Each wave is and its necessary lines are each written out individually to the cfg
    file. This method is preferred over AmpTools loops to manipulate individual waves
    much easier. The amplitude lines are written in the following order:
    - Comment the wave name
    - Write the vecps_refl amplitude line
    - Write the OmegaDalitz amplitude line
    - Write the initialization line
    - Write the parScale line
    - Repeat all of this for each orientation
    - Repeat for the other reflectivity sum

    Args:
        cfg_file (TextIO): cfg file being written to
        waves (List[Wave]): list of Wave dataclasses that contain all necessary info
            See 'get_waves_and_breit_wigners' for details
        force_refl (int): if non-zero, only the chosen reflectivity will be written
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
        if force_refl and wave.reflectivity != force_refl:
            continue

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


def write_isotropic_background(
    cfg_file: TextIO, is_dalitz: bool, orientations: List[str], reaction: str
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


def write_phase_convention(
    cfg_file: TextIO,
    wave: str,
    force_refl: int,
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
        wave (str): wave written in jpml format that will be set to be purely real
        force_refl (int): if non-zero, only the chosen reflectivity will be written
        orientations (List[str]): diamond orientation settings
        reaction (str): user-defined reaction name
    """

    # ensure last character of wave (L value) is capitalized
    wave = wave[:-1] + wave[-1].upper()

    cfg_file.write("\n# fix phase\n")

    short_refl_dict = {"ImagNegSign": -1, "ImagPosSign": 1}
    for ont in orientations:
        for sum_label, refl_sign in short_refl_dict.items():
            if force_refl and refl_sign != force_refl:
                continue
            cfg_file.write(
                f"initialize {reaction}_{POL_DICT[ont]['angle']}"
                f"::{sum_label}::{wave} cartesian 100 0 real\n"
            )

    return


def write_breit_wigners(
    cfg_file: TextIO,
    waves: List[Wave],
    breit_wigners: List[dict],
    force_refl: int,
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
        force_refl (int): if non-zero, only the chosen reflectivity will be written
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
                if refl_sign != wave.reflectivity or (
                    force_refl and refl_sign != force_refl
                ):
                    continue
                for ont in orientations:
                    cfg_file.write(
                        f"amplitude {reaction}_{POL_DICT[ont]['angle']}"
                        f"::{sum_label}::{wave.name}"
                        f" BreitWigner [{bw['jp']}_mass] [{bw['jp']}_width]"
                        f" {bw['spin']} {bw['indices']}\n"
                    )

    return


# TODO: finish converting all the methods to instead write out each amplitude line-by-line
def write_ds_ratio(
    cfg_file: TextIO,
    ds_option: str,
    force_refl: int,
    is_phaselock: bool,
    orientations: List[str],
    reaction: str,
) -> None:
    """Writes lines in cfg file template for a D to S wave ratio to be set

    Args:
        cfg_file (TextIO): cfg file being written to
        ds_option (str): user option to fix the ratio to E852 values, split it between,
            reflectivities, or let float between 0 and 1
        force_refl (int): if non-zero, only the chosen reflectivity will be written
        is_phaselock (bool): user option to constrain phases between eJPl amplitudes.
            In this function it fixes the "dsphase" parameter to 0
        orientations (List[str]): diamond orientation settings
        reaction (str): user-defined reaction name
    """

    cfg_file.write(
        "\n# constrain S and D waves to same amplitude and set scale factor"
        " for D/S ratio\n"
    )

    # setup variables used in all cases
    if ds_option == "fixed":
        ratio_opt = "fixed"
        phase_opt = "fixed"
    elif is_phaselock:
        ratio_opt = "bounded 0.0 1.0"
        phase_opt = "fixed"
    else:
        ratio_opt = "bounded 0.0 1.0"
        phase_opt = "bounded -3.14159 3.14159"

    phase_val = "0" if is_phaselock else "0.184"

    refl_dict = {
        "ImagNegSign": -1,
        "RealPosSign": -1,
        "ImagPosSign": 1,
        "RealNegSign": 1,
    }

    if force_refl == -1:
        del refl_dict["RealNegSign"]
        del refl_dict["ImagPosSign"]
    if force_refl == +1:
        del refl_dict["RealPosSign"]
        del refl_dict["ImagNegSign"]

    params, phases = ["dsratio"], ["dphase"]
    if ds_option == "split":
        params, phases = ["dsratio_p", "dsratio_m"], ["dphase_p", "dphase_m"]

    # loop over parameter and phase options and write out their lines
    for par, phase in zip(params, phases):
        if ds_option == "fixed":
            par_range = ""
        else:
            par_range = f"parRange {par} 0.0 1.0"

        if is_phaselock or ds_option == "fixed":
            phase_range = ""
        else:
            phase_range = f"parRange {phase} -3.14159 3.14159"

        cfg_file.write(
            f"parameter {par} 0.27 {ratio_opt}\n"
            f"parameter {phase} {phase_val} {phase_opt}\n"
            f"{par_range}\n"
            f"{phase_range}\n"
        )

        # if is_DS_split then certain reflectivities are paired to certain pars
        refls = ""
        for key, val in refl_dict.items():
            if "_p" in (par or phase) and val == 1:
                refls += f"{key} "
            elif "_m" in (par or phase) and val == -1:
                refls += f"{key} "

        if refls == "":
            refls = " ".join(refl_dict.keys())

        cfg_file.write(f"\nloop LOOPSUM {refls}\n")

    # perform D/S constraints for each orientation
    for ont in orientations:
        for par, phase in zip(params, phases):
            cfg_file.write(
                f"constrain {reaction}_{POL_DICT[ont]['angle']}::LOOPSUM::1ppd"
                f" {reaction}_{POL_DICT[ont]['angle']}::LOOPSUM::1pps"
                f"\nconstrain {reaction}_{POL_DICT[ont]['angle']}::LOOPSUM::1p0d"
                f" {reaction}_{POL_DICT[ont]['angle']}::LOOPSUM::1p0s"
                f"\nconstrain {reaction}_{POL_DICT[ont]['angle']}::LOOPSUM::1pmd"
                f" {reaction}_{POL_DICT[ont]['angle']}::LOOPSUM::1pms"
                f"\namplitude {reaction}_{POL_DICT[ont]['angle']}::LOOPSUM::1ppd"
                f" ComplexCoeff [{par}] [{phase}] MagPhi"
                f"\namplitude {reaction}_{POL_DICT[ont]['angle']}::LOOPSUM::1p0d"
                f" ComplexCoeff [{par}] [{phase}] MagPhi"
                f"\namplitude {reaction}_{POL_DICT[ont]['angle']}::LOOPSUM::1pmd"
                f" ComplexCoeff [{par}] [{phase}] MagPhi\n\n"
            )

    return


def write_phaselock(
    cfg_file: TextIO,
    waves: list,
    force_refl: int,
    init_refl: int,
    init_real: float,
    orientations: List[str],
    reaction: str,
) -> None:
    """Lock m-projections to have the same phase

    Phaselocking refers to the model that the m-projections for a certain JPL
    are all given the same phase value, meaning no phase difference can occur
    within a eJPL combination

    Args:
        cfg_file (TextIO): cfg file being written to
        waves (list): list of dict entries that contains JP, spin, parity, m, and l.
            See 'get_waves_and_breit_wigners' for details
        force_refl (int): if non-zero, only the chosen reflectivity will be written
        init_refl (int): user option to init nonzero values for one reflectivity
        init_real (float): user option to init the real value for amplitudes. Default
            is 100 (or 0 if init/force_refl is set)
        orientations (list): diamond orientation settings
        reaction (str): user-defined reaction name
    """

    cfg_file.write(f"\n{'#'*8} phaselock {'#'*8}\n")

    # write parameters for each phase, one each for plus and minus reflectivity.
    #   first wave is fixed to be 0
    fix_phase = ""

    for wave in waves:
        for l in wave["l"]:
            if fix_phase == "":
                fix_phase = l

            if l == fix_phase:
                option = "fixed"
                par_range = ""
            else:
                option = "bounded -3.14159 3.14159"
                par_range = (
                    "parRange"
                    f" phase_p{int_to_char(l,True).capitalize()}"
                    " -3.14159 3.14159\n"
                )

            if force_refl != -1:
                cfg_file.write(
                    f"parameter phase_p{int_to_char(l,True).capitalize()} 0"
                    f" {option}\n"
                    f"{par_range}"
                )
            if force_refl != 1:
                cfg_file.write(
                    f"parameter phase_m{int_to_char(l,True).capitalize()} 0"
                    f" {option}\n"
                    f"{par_range}"
                )

    # write new initialization line to overwrite previous ones in "write_wave_loops",
    #   and add PhaseOffset amplitude line
    for ont in orientations:
        for wave in waves:
            j = wave["spin"]
            p = int_to_char(wave["parity"])
            refl_dict = {
                "ImagNegSign": -1,
                "RealPosSign": -1,
                "ImagPosSign": 1,
                "RealNegSign": 1,
            }

            for coh_sum_label, refl in refl_dict.items():
                if force_refl and refl != force_refl:
                    continue
                cfg_file.write("\n")
                refl_str = "p" if refl == 1 else "m"

                for m in wave["m"]:
                    for l in wave["l"]:

                        amplitude = (
                            f"{reaction}_{POL_DICT[ont]['angle']}::{coh_sum_label}::"
                            f"{j}{p}{int_to_char(m)}{int_to_char(l,True)}"
                        )
                        cfg_file.write(
                            f"initialize {amplitude}"
                            " cartesian {re} 0 real\n".format(
                                re=0 if init_refl and refl != init_refl else init_real
                            )
                        )
                        cfg_file.write(
                            f"amplitude {amplitude} PhaseOffset"
                            f" [phase_"
                            f"{refl_str}{int_to_char(l,True).capitalize()}]\n"
                        )
    return


def int_to_char(x: int, is_ang_momentum: bool = False) -> str:
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


def is_reference_wave_in_waveset(reference_wave: str, waves: List[Wave]) -> bool:
    """Check if reference wave is in waveset

    Args:
        reference_wave (str): user chosen wave that will be fixed to be real, in jpml
            format
        waves (List[Wave]): list of Wave dataclasses that contain all necessary info.
            See 'get_waves_and_breit_wigners' for details

    Returns:
        bool: True if reference wave is in waveset, False otherwise
    """

    if reference_wave == "":
        return False
    elif len(reference_wave) != 4:
        raise ValueError("Reference wave must be in jpml format")

    for wave in waves:
        jpml = wave.name[1:]
        if jpml.lower() == reference_wave.lower():
            return True
    return False
