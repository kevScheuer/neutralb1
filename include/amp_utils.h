/**
 * @file amp_utils.cc
 * @brief Utility functions for handling amplitudes in PWA fits.
 *
 */

#ifndef AMP_UTILS_H
#define AMP_UTILS_H

#include <string>
#include <tuple>

/**
 * Parses an amplitude string into its quantum number components.
 * 
 * Assumes the amplitude is the "full" amplitude, that is, 
 * <reaction>::<reflectivitySum>::<eJPmL>, where the <eJPmL> amplitude name is in the 
 * vector-pseudoscalar style format:
 *    e = reflectivity (either p [+] or m [-])
 *    J = total angular momentum (positive integers including 0)
 *    P = parity (either p [+] or m [-])
 *    m = m-projection (l [-3], n [-2], m [-1], 0 [0], p [1], q [2], r [3])
 *    L = orbital angular momentum (standard letter convention: S, P, D, F, ...).
 * @param amplitude The full amplitude string to parse
 * @return A tuple containing the components:
 *         (reflectivity, total angular momentum, parity, m-projection,
 *         orbital angular momentum)
 * @throw std::invalid_argument If the amplitude string is unrecognized
 */
std::tuple<std::string, std::string, std::string, std::string, std::string> parse_amplitude(std::string amplitude);

#endif // AMP_UTILS_H