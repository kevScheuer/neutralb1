"""This module contains functions for calculating various physics-related quantities

Most of these are pythonized versions of functions from the halld_sim C++ codebase.
"""

import numpy as np


def breakup_momentum(
    parent_mass: float, daughter1_mass: float, daughter2_mass: float
) -> float:
    """Breakup momentum of a parent particle into two daughter particles

    Pythonized version of https://github.com/JeffersonLab/halld_sim/blob/master/src/
    libraries/AMPTOOLS_AMPS/breakupMomentum.cc. Returned value is independent of which
    particle is considered daughter1 or daughter2. Take care to units are consistent
    between masses.

    Args:
        parent_mass (float): mass of parent particle
        daughter1_mass (float): mass of one of the daughter particles
        daughter2_mass (float): mass of the other daughter particle

    Returns:
        float: Breakup momentum of the parent particle.
    """
    return np.sqrt(
        np.abs(
            np.power(parent_mass, 4)
            + np.power(daughter1_mass, 4)
            + np.power(daughter2_mass, 4)
            - 2.0 * np.square(parent_mass) * np.square(daughter1_mass)
            - 2.0 * np.square(parent_mass) * np.square(daughter2_mass)
            - 2.0 * np.square(daughter1_mass) * np.square(daughter2_mass)
        )
    ) / (2.0 * parent_mass)


def barrier_factor(
    parent_mass: float, daughter1_mass: float, daughter2_mass: float, l: int
) -> float:
    """Blatt-Weisskopf barrier factor for a given parent-daughter system

    Pythonized version of https://github.com/JeffersonLab/halld_sim/blob/master/src/
    libraries/AMPTOOLS_AMPS/barrierFactor.cc. This parameterization often suppresses
    values near threshold, and amplifies values at higher masses, especially for large
    l values.

    Args:
        parent_mass (float): mass of the parent particle
        daughter1_mass (float): mass of the first daughter particle
        daughter2_mass (float): mass of the second daughter particle
        l (int): orbital angular momentum quantum number of the parent particle

    Returns:
        float: barrier factor
    """

    q = np.abs(breakup_momentum(parent_mass, daughter1_mass, daughter2_mass))
    z = np.square(q) / np.square(0.1973)

    match l:
        case 0:
            barrier = 1.0
        case 1:
            barrier = (2.0 * z) / (z + 1.0)
        case 2:
            barrier = (13.0 * np.square(z)) / (np.square(z - 3.0) + 9.0 * z)
        case 3:
            barrier = (277.0 * np.power(z, 3)) / (
                z * np.square(z - 15.0) + 9.0 * np.square(2.0 * z - 5.0)
            )
        case 4:
            barrier = (12746.0 * np.power(z, 4)) / (
                np.square((np.square(z) - 45.0 * z + 105.0))
                + 25.0 * z * np.square(2.0 * z - 21.0)
            )
        case _:
            barrier = 0.0

    return np.sqrt(barrier)


def breit_wigner(
    mass: float,
    bw_mass: float,
    bw_width: float,
    bw_l: int,
    daughter1_mass: float = 0.1349768,
    daughter2_mass: float = 0.78266,
) -> complex:
    """Halld_sim parameterization of the breit wigner function

    To avoid discrepancies, this function copies the parameterization of the breit
    wigner function found in https://github.com/JeffersonLab/halld_sim/blob/master/src/
    libraries/AMPTOOLS_AMPS/BreitWigner.cc.

    NOTE: all masses / widths must be in GeV. Daughter particle masses are typically
    obtained from event 4-vectors, but since they're not available post-fit, we
    approximate them by constraining their mass to the pdg value

    Args:
        mass (float): mass value to evaluate breit wigner function at. If using a mass
            bin, approximate by passing the center of the bin
        bw_mass (float): breit wigner central mass
        bw_width (float): breit wigner width
        bw_l (int): orbital angular momentum of the breit wigner i.e. rho->(omega,pi0)
            is in the P-wave, so bw_l = 1
        daughter1_mass (float, optional): mass of 1st daughter particle.
            Defaults to 0.1349768 for the pi0 mass
        daughter2_mass (float, optional): mass of 2nd daughter particle.
            Defaults to 0.78266 for the omega mass

    Returns:
        complex: value of the breit wigner at the mass value
    """

    q0 = np.abs(breakup_momentum(bw_mass, daughter1_mass, daughter2_mass))
    q = np.abs(breakup_momentum(mass, daughter1_mass, daughter2_mass))

    F0 = barrier_factor(bw_mass, daughter1_mass, daughter2_mass, bw_l)
    F = barrier_factor(mass, daughter1_mass, daughter2_mass, bw_l)

    width = bw_width * (bw_mass / mass) * (q / q0) * np.square(F / F0)

    numerator = complex(np.sqrt((bw_mass * bw_width) / np.pi), 0.0)
    denominator = complex(np.square(bw_mass) - np.square(mass), -1.0 * bw_mass * width)

    return F * numerator / denominator
