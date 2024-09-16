""" Writes a template cfg file for binned fits to be performed

This script is part of the batch submission process, called with submit.py. Its
role is to write a cfg file that acts as a template for all the binned fits
to overwrite. This file accepts a users waveset, along with any modifications to
it, and writes the necessary AmpTools cfg file lines to fit using that waveset

NOTE: This file (and other batch files) are tuned specifically to omegapi0
TODO:     
    Make refl_dict a constant and handle appropriately in each function
    Have phaselock function handle dominant eJPmL arg
"""

import argparse
import typing
from typing import List

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

# orientations: list,


def main(args: dict):
    """Write config using parsed args from submit.py"""

    # convert the user requested wave string into lists of waves/bw's that contain all
    # the necessary info for each wave/bw
    waves, breit_wigners = get_waves_and_breit_wigners(args["waveset"])

    # check that the reference wave is in the waveset if requested, otherwise set the
    # reference to be the first wave in the waveset
    if args["phase_reference"]:
        if not is_reference_wave_in_waveset(args["phase_reference"], waves):
            raise ValueError(f"{args['phase_reference']} is not in waveset")
    else:
        args["phase_reference"] = waves[0]

    # parse out special options aside from waves and bw's
    is_iso = True if "iso" in args["waveset"] else False
    is_dalitz = False if "nodalitz" in args["waveset"] else True

    # copy common part of omegapi config files to a new file that will be the
    #   template for each bin
    ftemplate = open(f"submission/batch_scripts/{args['template_name']}", "r")
    template_data = ftemplate.read()
    ftemplate.close()
    template_data = template_data.replace("FITNAME", args["reaction"])
    fout = open("submission/batch_scripts/fit.cfg", "w")
    fout.write(template_data)

    # Write LOOP over polarization orientations and setup corresponding reaction names
    write_orientation_loop(fout, args["orientations"], args["reaction"])

    # Write AmpTools LOOP for each JP
    write_wave_loops(
        fout,
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
            fout, is_dalitz, args["orientations"], args["reaction"]
        )

    # write lines to remove overall unconstrained phase factor
    write_phase_convention(
        fout,
        args["phase_reference"],
        args["force_refl"],
        args["orientations"],
        args["reaction"],
    )

    # writes breit wigners (if present in waveset)
    if breit_wigners:
        write_breit_wigners(
            fout,
            waves,
            breit_wigners,
            args["force_refl"],
            args["orientations"],
            args["reaction"],
        )

    # write ds ratio if 1p amp is used, and if DS ratio isn't free
    is_1p_present = False
    for wave in waves:
        if wave.get("spin") == 1 and wave.get("parity") == 1:
            is_1p_present = True
    if args["ds_ratio"] != "free" and is_1p_present:
        write_ds_ratio(
            fout,
            args["ds_ratio"],
            args["force_refl"],
            args["phaselock"],
            args["orientations"],
            args["reaction"],
        )

    # user option to constrain phases between m projections
    if args["phaselock"]:
        write_phaselock(
            fout,
            waves,
            args["force_refl"],
            args["init_refl"],
            args["init_real"],
            args["orientations"],
            args["reaction"],
        )

    fout.close()

    return


def get_waves_and_breit_wigners(waveset: List[str]) -> tuple[list, list]:
    """Converts user passed waveset into detailed dicts

    Args:
        waveset: user passed waveset to be built in cfg file e.g. ["1p", "1m", "iso"]
    Returns:
        Two lists containing dicts, that convert each requested wave into the
        proper quantum numbers Ex:
            1p -> {"spin":1, "parity":-1, "m": [-1,0,1], "l":1}
    """
    # only consider m <=2 for photon beam
    wave_dict = {
        "0m": {"spin": 0, "parity": -1, "m": [0], "l": [1]},
        "1p": {"spin": 1, "parity": +1, "m": [-1, 0, 1], "l": [0, 2]},
        "1m": {"spin": 1, "parity": -1, "m": [-1, 0, 1], "l": [1]},
        "2p": {"spin": 2, "parity": +1, "m": [-2, -1, 0, 1, 2], "l": [2]},
        "2m": {"spin": 2, "parity": -1, "m": [-2, -1, 0, 1, 2], "l": [1, 3]},
        "3m": {"spin": 3, "parity": -1, "m": [-2, -1, 0, 1, 2], "l": [3]},
    }
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
    if "b1free" in waveset:
        breit_wigner_dict["b1"]["fix"] = False
    if "rhofree" in waveset:
        breit_wigner_dict["rho"]["fix"] = False

    # add requested waves and breit wigners to lists
    waves = []
    breit_wigners = []
    for key, val in wave_dict.items():
        if key in waveset:
            waves.append(val)
    for key, val in breit_wigner_dict.items():
        if key in waveset:
            breit_wigners.append(val)

    # check that amplitudes are present for any requested breit wigners
    if breit_wigners:
        jp_list = []
        for wave in waves:
            jp_list.append(str(wave["spin"]) + int_to_char(wave["parity"]))
        for bw in breit_wigners:
            if bw["jp"] not in jp_list:
                raise ValueError(
                    f"Breit Wigner JP={bw['jp']} not found in wave jp values: {jp_list}"
                )

    return waves, breit_wigners


def write_orientation_loop(
    fout: typing.TextIO, orientations: List[str], reaction: str
) -> None:
    """Writes the AmpTools 'loop' statements to handle polarization orientations

    Args:
        fout (typing.TextIO): cfg file being written to
        orientations (List[str]): diamond orientation settings
        reaction (str): user-defined reaction name
    """

    fout.write(f"\n\n{'#'*8} LOOP OVER POLARIZATION STATES {'#'*8}\n")

    fix_string = "fixed"
    for ont in orientations:
        pol_scales = {
            "PARA_0": "1.0",  # potential for arg that changes these values
            "PERP_45": "1.0",
            "PERP_90": "1.0",
            "PARA_135": "1.0",
        }
        fout.write(
            f"parameter parScale_{POL_DICT[ont]['angle']}"
            f" {pol_scales[ont]} {fix_string}\n"
        )
        if fix_string == "fixed":  # effectively "fixes" the first ont in orientations
            fix_string = ""

    # write necessary starting loops
    reaction_and_angles = ""
    parScale_and_angles = ""
    loop_data_str = ""
    for ont in orientations:
        reaction_and_angles += f" {reaction}_{POL_DICT[ont]['angle']}"
        parScale_and_angles += f" [parScale_{POL_DICT[ont]['angle']}]"
        loop_data_str += f" anglesOmegaPiAmplitude_{POL_DICT[ont]['angle']}.root"

    fout.write(f"\nloop LOOPREAC{reaction_and_angles}\n")
    fout.write(f"loop LOOPSCALE{parScale_and_angles}\n")

    fout.write(
        f"loop LOOPGENMC {'anglesOmegaPiPhaseSpace.root ' * len(orientations)} \n"
        f"loop LOOPACCMC {'anglesOmegaPiPhaseSpaceAcc.root ' * len(orientations)} \n"
        f"loop LOOPDATA{loop_data_str}\n"
    )

    # write reaction, data reader, and sum lines
    fout.write(
        f"\n{'#'*8} DATA, REACTIONS, AND SUMS {'#'*8}\n"
        "reaction LOOPREAC Beam Proton Pi01 Pi02 Pi+ Pi-\n"
        "genmc LOOPREAC ROOTDataReader LOOPGENMC\n"
        "accmc LOOPREAC ROOTDataReader LOOPACCMC\n"
        "data LOOPREAC ROOTDataReader LOOPDATA\n\n"
        "sum LOOPREAC ImagNegSign RealNegSign RealPosSign ImagPosSign\n\n"
    )

    return


####################### WRITE AMPLITUDES #######################
def write_wave_loops(
    fout: typing.TextIO,
    waves: list,
    force_refl: int,
    init_refl: int,
    init_real: float,
    init_imag: float,
    is_dalitz: bool,
    frame: str,
    orientations: List[str],
    reaction: str,
) -> None:
    """Write all AmpTools LOOP statements for the waves

    Each JP value has a finite combination of m-projections and angular momenta
    values "l". These can be looped over in the cfg file to write, initialize,
    and constrain all amplitudes. This function writes those loops, with
    edits according to the user optional arguments. Note that Orientation loops must be
    handled here "manually" (without the use of AmpTools LOOPS)

    Args:
        fout (typing.TextIO): cfg file being written to
        waves (list): list of dict entries that contains JP, spin, parity, m, and l.
            See 'get_waves_and_breit_wigners' for details
        force_refl (int): user option to only use 1 reflectivity
        init_refl (int): user option to init nonzero values for one reflectivity
        init_real (float): user option to init the real value for amplitudes. Default
            is 100 (or 0 if init/force_refl is set)
        init_imag (float): user option to init the imaginary value for amplitudes.
            Default is 100 (or 0 if init/force_refl is set)
        is_dalitz (bool): controls whether OmegaDalitz amplitudes are written
        frame (str): option for using a special decay angle frame
        orientations (List[str]): diamond orientation settings
        reaction (str): user-defined reaction name
    """

    refl_dict = {
        "ImagNegSign": -1,
        "RealPosSign": -1,
        "ImagPosSign": 1,
        "RealNegSign": 1,
    }

    # Create loops for each J^P value
    for wave in waves:
        j = wave["spin"]
        p = int_to_char(wave["parity"])

        fout.write(f"{'#'*8} spin {j} parity {wave['parity']} {'#'*8}\n")

        # Write loop header lines for ampltiudes, m projections, and l values
        amp_loop = f"loop LOOPAMP{j}{p}"
        m_loop = f"loop LOOPM{j}{p}"
        l_loop = f"loop LOOPL{j}{p}"

        for m in wave["m"]:
            for l in wave["l"]:
                amp_loop += f" {j}{p}{int_to_char(m)}{int_to_char(l,True)}"
                m_loop += f" {m}"
                l_loop += f" {l}"

        fout.write(f"{amp_loop}\n{m_loop}\n{l_loop}\n")

        # each orientation needs it own "reaction" in each amplitude
        for ont in orientations:
            fout.write(f"{'#'*4} Orientation {ont} {'#'*4}\n")

            # 4 coherent sums for reflectivity and sign of (1 +/- P_gamma)
            for coh_sum_label, refl in refl_dict.items():
                if force_refl and refl != force_refl:
                    continue

                real = 1 if "Real" in coh_sum_label else -1
                sign = 1 if "Pos" in coh_sum_label else -1

                # Write amplitude line NOTE: 'omega3pi' is specific to omega->3pi decays
                amplitude = (
                    f"{reaction}_{POL_DICT[ont]['angle']}::"
                    f"{coh_sum_label}::LOOPAMP{j}{p}"
                )
                fout.write(
                    f"amplitude {amplitude} Vec_ps_refl {j} LOOPM{j}{p}"
                    f" LOOPL{j}{p} {real} {sign}"
                    f" {POL_DICT[ont]['angle']} {POL_DICT[ont]['fraction']}"
                    f" omega3pi {frame}\n"
                )
                if is_dalitz:
                    fout.write(f"amplitude {amplitude} OmegaDalitz dalitz\n")

                # write amplitude's initialization line
                fout.write(
                    "initialize {amp} cartesian {re} {im}\n".format(
                        amp=amplitude,
                        re=0 if init_refl and refl != init_refl else init_real,
                        im=0 if init_refl and refl != init_refl else init_imag,
                    )
                )

                # attach polarization scale parameter to amplitude
                fout.write(f"scale {amplitude} [parScale_{POL_DICT[ont]['angle']}]\n\n")

            # constrain amplitudes across coherent sums
            constraint_line = (
                f"constrain"
                f" {reaction}_{POL_DICT[ont]['angle']}::ImagNegSign::LOOPAMP{j}{p}"
                f" {reaction}_{POL_DICT[ont]['angle']}::RealPosSign::LOOPAMP{j}{p}\n"
                f"constrain"
                f" {reaction}_{POL_DICT[ont]['angle']}::RealNegSign::LOOPAMP{j}{p}"
                f" {reaction}_{POL_DICT[ont]['angle']}::ImagPosSign::LOOPAMP{j}{p}"
            )
            if force_refl == -1:
                constraint_line = constraint_line.split("\n")[0]
            elif force_refl == 1:
                constraint_line = constraint_line.split("\n")[1]

            fout.write(constraint_line + "\n\n")

            # constrain amplitudes across polarizations
            #   all amplitudes are constrained to the first orientation in the list
            for coh_sum_label, refl in refl_dict.items():
                if force_refl and refl != force_refl:
                    continue
                fout.write(
                    f"constrain"
                    f" {reaction}_{POL_DICT[orientations[0]]['angle']}::"
                    f"{coh_sum_label}::LOOPAMP{j}{p}"
                    f" {reaction}_{POL_DICT[ont]['angle']}::"
                    f"{coh_sum_label}::LOOPAMP{j}{p}\n"
                )
            fout.write("\n")

    return


def write_isotropic_background(
    fout: typing.TextIO, is_dalitz: bool, orientations: List[str], reaction: str
) -> None:
    """Writes lines in cfg file template for isotropic background amplitude

    Args:
        fout (typing.TextIO): cfg file being written to
        is_dalitz (bool): controls whether OmegaDalitz amplitudes are written
        orientations (List[str]): diamond orientation settings
        reaction (str): user-defined reaction name
    """
    fout.write(f"{'#'*8} isotropic background {'#'*8}")

    for ont in orientations:
        if is_dalitz:
            dalitz = (
                f"amplitude {reaction}_{POL_DICT[ont]['angle']}"
                "::Bkgd::isotropic OmegaDalitz dalitz"
            )
        else:
            dalitz = ""

        fout.write(
            f"\nsum {reaction}_{POL_DICT[ont]['angle']} Bkgd"
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
    fout: typing.TextIO,
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
        fout (typing.TextIO): cfg file being written to
        wave (str): spin, parity, m, and l of the reference wave
        force_refl (int): user option to only use 1 reflectivity
        orientations (List[str]): diamond orientation settings
        reaction (str): user-defined reaction name
    """

    fout.write("\n# fix phase\n")

    refl_dict = {"ImagNegSign": -1, "ImagPosSign": 1}
    for ont in orientations:
        for coh_sum_label, refl in refl_dict.items():
            if force_refl and refl != force_refl:
                continue
            fout.write(
                f"initialize {reaction}_{POL_DICT[ont]['angle']}"
                f"::{coh_sum_label}::{wave} cartesian 100 0 real\n"
            )

    return


def write_breit_wigners(
    fout: typing.TextIO,
    waves: list,
    breit_wigners: list,
    force_refl: int,
    orientations: List[str],
    reaction: str,
) -> None:
    """Write cfg file lines for breit wigners. only b1 and rho supported.

    A check is put before this function to ensure that breit_wigners is a non-
    empty list.

    Args:
        fout (typing.TextIO): cfg file being written to
        waves (list): list of dict entries that contains JP, spin, parity, m, and l.
            See 'get_waves_and_breit_wigners' for details
        breit_wigners (list): list of bw dicts that contain all necessary info
        force_refl (int): user option to only use 1 reflectivity
        orientations (List[str]): diamond orientation settings
        reaction (str): user-defined reaction name
    """

    fout.write("\n# Mass dependent amplitudes\n")
    refls = "ImagNegSign RealNegSign RealPosSign ImagPosSign"
    if force_refl == 1:
        refls = refls.replace("RealNegSign", "")
        refls = refls.replace("ImagPosSign", "")
    if force_refl == -1:
        refls = refls.replace("RealPosSign", "")
        refls = refls.replace("ImagNegSign", "")

    fout.write(f"loop LOOPSUMBW {refls}\n")

    # write mass and width parameters
    for bw in breit_wigners:
        fout.write(
            f"parameter {bw['jp']}_mass {bw['mass']:0.3f}"
            " {option}\n".format(
                option=(
                    "fixed"
                    if bw["fix"]
                    else (
                        f"bounded {bw['mass_low']:.3f} {bw['mass_high']:.3f}"
                        f"\nparRange {bw['jp']}_mass"
                        f" {bw['mass_low']:.3f} {bw['mass_high']:.3f}"
                    )
                )
            )
        )
        fout.write(
            f"parameter {bw['jp']}_width {bw['width']:.3f}"
            " {option}\n".format(
                option=(
                    "fixed"
                    if bw["fix"]
                    else (
                        f"bounded {bw['width_low']:.3f} {bw['width_high']:.3f}"
                        f"\nparRange {bw['jp']}_width"
                        f" {bw['width_low']:.3f} {bw['width_high']:.3f}"
                    )
                )
            )
        )

    jpml_set = set()

    # get all the different jpml combinations present in the waveset
    for wave in waves:
        j = str(wave["spin"])
        p = int_to_char(wave["parity"])
        for m in wave["m"]:
            for l in wave["l"]:
                jpml = j + p + int_to_char(m) + int_to_char(l, True)
                jpml_set.add(jpml)

    # construct the bw lines to be written
    for jpml in jpml_set:
        for bw in breit_wigners:
            for ont in orientations:
                fout.write(
                    f"amplitude {reaction}_{POL_DICT[ont]['angle']}::LOOPSUMBW::{jpml}"
                    f" BreitWigner [{bw['jp']}_mass] [{bw['jp']}_width]"
                    f" {bw['spin']} {bw['indices']}\n"
                )

    return


def write_ds_ratio(
    fout: typing.TextIO,
    ds_option: str,
    force_refl: int,
    is_phaselock: bool,
    orientations: List[str],
    reaction: str,
) -> None:
    """Writes lines in cfg file template for a D to S wave ratio to be set

    Args:
        fout (typing.TextIO): cfg file being written to
        ds_option (str): user option to fix the ratio to E852 values, split it between,
            reflectivities, or let float between 0 and 1
        force_refl (int): user option to only use 1 reflectivity
        is_phaselock (bool): user option to constrain phases between eJPl amplitudes.
            In this function it fixes the "dsphase" parameter to 0
        orientations (List[str]): diamond orientation settings
        reaction (str): user-defined reaction name
    """

    fout.write(
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

        fout.write(
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

        fout.write(f"\nloop LOOPSUM {refls}\n")

    # perform D/S constraints for each orientation
    for ont in orientations:
        for par, phase in zip(params, phases):
            fout.write(
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
    fout: typing.TextIO,
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
        fout (typing.TextIO): cfg file being written to
        waves (list): list of dict entries that contains JP, spin, parity, m, and l.
            See 'get_waves_and_breit_wigners' for details
        force_refl (int): user option to only use 1 reflectivity
        init_refl (int): user option to init nonzero values for one reflectivity
        init_real (float): user option to init the real value for amplitudes. Default
            is 100 (or 0 if init/force_refl is set)
        orientations (list): diamond orientation settings
        reaction (str): user-defined reaction name
    """

    fout.write(f"\n{'#'*8} phaselock {'#'*8}\n")

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
                fout.write(
                    f"parameter phase_p{int_to_char(l,True).capitalize()} 0"
                    f" {option}\n"
                    f"{par_range}"
                )
            if force_refl != 1:
                fout.write(
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
                fout.write("\n")
                refl_str = "p" if refl == 1 else "m"

                for m in wave["m"]:
                    for l in wave["l"]:

                        amplitude = (
                            f"{reaction}_{POL_DICT[ont]['angle']}::{coh_sum_label}::"
                            f"{j}{p}{int_to_char(m)}{int_to_char(l,True)}"
                        )
                        fout.write(
                            f"initialize {amplitude}"
                            " cartesian {re} 0 real\n".format(
                                re=0 if init_refl and refl != init_refl else init_real
                            )
                        )
                        fout.write(
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
        is_ang_momentum: if true, uses a seperate dictionary to perform the
            conversion. Uses the standard "s,p,d,f,g" wave notation
    Returns:
        single character
    """
    map = {-2: "n", -1: "m", 0: "0", +1: "p", +2: "q"}
    if is_ang_momentum:
        map = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g"}

    return map[x]


def is_reference_wave_in_waveset(reference_wave: str, waves: list) -> bool:
    """Check if reference wave is in waveset

    Args:
        reference_wave (str): user chosen wave that will be fixed to be real
        waves (list): list of dict entries that contains spin, parity, m, and l quantum
            numbers. See 'get_waves_and_breit_wigners' for details

    Returns:
        bool: True if reference wave is in waveset, False otherwise
    """
    for wave in waves:
        j = str(wave["spin"])
        p = int_to_char(wave["parity"])
        for m in wave["m"]:
            for l in wave["l"]:
                jpml = j + p + int_to_char(m) + int_to_char(l, True)
                if reference_wave == jpml:
                    return True
    return False
