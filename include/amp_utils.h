/**
 * @file amp_utils.h
 * @brief Utility functions for handling amplitudes in PWA fits.
 *
 */

#ifndef AMP_UTILS_H
#define AMP_UTILS_H

#include <string>

/**
 * @class AmplitudeParser
 * @brief A class for parsing and representing quantum numbers from amplitude strings.
 *
 * This class can be constructed from either an amplitude string or individual quantum
 * numbers, and provides access to both string and integer representations of the
 * quantum numbers.
 */
class AmplitudeParser
{
private:
    std::string e_str; ///< Reflectivity as string: "p" or "m"
    std::string J_str; ///< Total angular momentum as string: "0" to "9"
    std::string P_str; ///< Parity as string: "p" or "m"
    std::string m_str; ///< m-projection as string: "l", "n", "m", "0", "p", "q", "r"
    std::string L_str; ///< Orbital angular momentum as string: "S", "P", "D", "F", "G", ...

    int e_int; ///< Reflectivity as integer: +1 or -1
    int J_int; ///< Total angular momentum as integer: 0, 1, 2, ...
    int P_int; ///< Parity as integer: +1 or -1
    int m_int; ///< m-projection as integer: -3 to +3
    int L_int; ///< Orbital angular momentum as integer: 0=S, 1=P, 2=D, ...

    /**
     * @brief Convert string representations to integer representations
     */
    void compute_integers();

    /**
     * @brief Convert integer representations to string representations
     */
    void compute_strings();

    /**
     * @brief Parse an amplitude string into quantum number components
     * @param amplitude The full amplitude string to parse
     */
    void parse_from_amplitude(const std::string &amplitude);

public:
    /**
     * @brief Construct from amplitude string
     * @param amplitude The full amplitude string in format reaction::sum::eJPmL
     * @throw std::invalid_argument If the amplitude string is unrecognized
     */
    AmplitudeParser(const std::string &amplitude);

    /**
     * @brief Construct from individual quantum numbers as strings
     * @param e Reflectivity (p or m)
     * @param J Total angular momentum (0-9)
     * @param P Parity (p or m)
     * @param m m-projection (l, n, m, 0, p, q, r)
     * @param L Orbital angular momentum (S, P, D, F, G, ...)
     * @throw std::invalid_argument If any quantum number is invalid
     */
    AmplitudeParser(const std::string &e, const std::string &J, const std::string &P,
                    const std::string &m, const std::string &L);

    /**
     * @brief Construct from individual quantum numbers as integers
     * @param e Reflectivity (+1 or -1)
     * @param J Total angular momentum (0, 1, 2, ...)
     * @param P Parity (+1 or -1)
     * @param m m-projection (-3 to +3)
     * @param L Orbital angular momentum (0=S, 1=P, 2=D, ...)
     * @throw std::invalid_argument If any quantum number is invalid
     */
    AmplitudeParser(int e, int J, int P, int m, int L);

    // String accessors
    std::string get_e_str() const { return e_str; }
    std::string get_J_str() const { return J_str; }
    std::string get_P_str() const { return P_str; }
    std::string get_m_str() const { return m_str; }
    std::string get_L_str() const { return L_str; }

    // Integer accessors
    int get_e_int() const { return e_int; }
    int get_J_int() const { return J_int; }
    int get_P_int() const { return P_int; }
    int get_m_int() const { return m_int; }
    int get_L_int() const { return L_int; }

    /**
     * @brief Get the amplitude name (eJPmL format)
     * @return String representation of the amplitude name
     */
    std::string get_amplitude_name() const;
};

#endif // AMP_UTILS_H