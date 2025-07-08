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
 * Parses an amplitude string into its components. Assumes the amplitude is in the
 * vector-pseudoscalar style 'eJPmL' format, where:
 *    e = reflectivity (either p [+] or m [-])
 *    J = total angular momentum (positive integers including 0)
 *    P = parity (either p [+] or m [-])
 *    m = m-projection (either p [+], m [-], or 0)
 * @param amplitude The amplitude string to parse
 * @return A tuple containing the components:
 *         (reflectivity, total angular momentum, parity, m-projection,
 *         orbital angular momentum)
 * @throw std::invalid_argument If the format is invalid, returns empty strings for all
 *        components.
 */
std::tuple<std::string, std::string, std::string, std::string, std::string> parse_amplitude(std::string amplitude);

#endif // AMP_UTILS_H