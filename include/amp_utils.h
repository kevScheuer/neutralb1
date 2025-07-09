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
 * @struct QuantumNumbers
 * @brief Represents the quantum numbers extracted from an amplitude string.
 */
struct QuantumNumbers {
    std::string e;
    std::string J;
    std::string P;
    std::string m;
    std::string L;
    
    /**
     * @brief Convert reflectivity to integer (+1 for 'p', -1 for 'm')
     * @return Integer representation of reflectivity
     */
    int get_e_int() const;
    
    /**
     * @brief Convert total angular momentum to integer
     * @return Integer value of J
     */
    int get_J_int() const;
    
    /**
     * @brief Convert parity to integer (+1 for 'p', -1 for 'm')
     * @return Integer representation of parity
     */
    int get_P_int() const;
    
    /**
     * @brief Convert m-projection to integer (-3 to +3)
     * @return Integer value of m-projection
     */
    int get_m_int() const;
    
    /**
     * @brief Convert orbital angular momentum to integer (0=S, 1=P, 2=D, etc.)
     * @return Integer value of L
     */
    int get_L_int() const;
};

/**
 * @struct QuantumNumbersInt
 * @brief Represents quantum numbers as integers for calculations.
 */
struct QuantumNumbersInt {
    int e;  ///< Reflectivity: +1 or -1
    int J;  ///< Total angular momentum: 0, 1, 2, ...
    int P;  ///< Parity: +1 or -1
    int m;  ///< m-projection: -3 to +3
    int L;  ///< Orbital angular momentum: 0=S, 1=P, 2=D, ...
};

/**
 * @brief Parses an amplitude string into its quantum number components.
 * 
 * @details
 * Assumes the amplitude is the "full" amplitude, that is, 
 * <reaction>::<reflectivitySum>::<eJPmL>, where the <eJPmL> amplitude name is in the 
 * vector-pseudoscalar style format:
 *    e = reflectivity (either p [+] or m [-])
 *    J = total angular momentum (positive integers including 0)
 *    P = parity (either p [+] or m [-])
 *    m = m-projection (l [-3], n [-2], m [-1], 0 [0], p [1], q [2], r [3])
 *    L = orbital angular momentum (standard letter convention: S, P, D, F, ...).
 * @param amplitude The full amplitude string to parse
 * @return A QuantumNumbers struct containing the quantum numbers
 * @throw std::invalid_argument If the amplitude string is unrecognized
 */
QuantumNumbers parse_amplitude(std::string amplitude);

/**
 * @brief Convert string quantum numbers to integer representation
 * @param qn QuantumNumbers struct with string values
 * @return QuantumNumbersInt struct with integer values
 * @throw std::invalid_argument If conversion fails
 */
QuantumNumbersInt to_integers(const QuantumNumbers& qn);

#endif // AMP_UTILS_H