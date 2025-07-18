/**
 * @file AmplitudeParser.h
 * @brief Parse amplitude into its quantum number strings and integers
 *
 */

#ifndef AMPLITUDE_PARSER_H
#define AMPLITUDE_PARSER_H

#include <string>

// TODO: Later, have an amplitude style be inherited from somewhere, where the parsing
//  becomes dependent on vec-ps, two-ps, or other defined styles with associated regex
//  strings. Provide an option to override the style somewhere here too.

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
    std::string m_e_str; ///< Reflectivity as string: "p" or "m"
    std::string m_J_str; ///< Total angular momentum as string: "0" to "9"
    std::string m_P_str; ///< Parity as string: "p" or "m"
    std::string m_m_str; ///< m-projection as string: "l", "n", "m", "0", "p", "q", "r"
    std::string m_L_str; ///< Orbital angular momentum as string: "S", "P", "D", "F", "G", ...

    std::string m_reaction; ///< Reaction part of the amplitude name
    std::string m_sum;      ///< Sum part of the amplitude name

    int m_e_int; ///< Reflectivity as integer: +1 or -1
    int m_J_int; ///< Total angular momentum as integer: 0, 1, 2, ...
    int m_P_int; ///< Parity as integer: +1 or -1
    int m_m_int; ///< m-projection as integer: -3 to +3
    int m_L_int; ///< Orbital angular momentum as integer: 0=S, 1=P, 2=D, ...

    /**
     * @brief Convert string representations to integer representations
     */
    void compute_integers();

    /**
     * @brief Convert integer representations to string representations
     */
    void compute_strings();

    /**
     * @brief Parse an amplitude string into its reaction, sum, and quantum number components
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
    std::string get_e_str() const { return m_e_str; } ///< Reflectivity as string
    std::string get_J_str() const { return m_J_str; } ///< Total angular momentum as string
    std::string get_P_str() const { return m_P_str; } ///< Parity as string
    std::string get_m_str() const { return m_m_str; } ///< m-projection as string
    std::string get_L_str() const { return m_L_str; } ///< Orbital angular momentum as string

    // Integer accessors
    int get_e_int() const { return m_e_int; } ///< Reflectivity as integer
    int get_J_int() const { return m_J_int; } ///< Total angular momentum as integer
    int get_P_int() const { return m_P_int; } ///< Parity as integer
    int get_m_int() const { return m_m_int; } ///< m-projection as integer
    int get_L_int() const { return m_L_int; } ///< Orbital angular momentum as integer

    // Amplitude accessors
    std::string get_amplitude_reaction() const { return m_reaction; }; ///< Reaction part of the amplitude name
    std::string get_amplitude_sum() const { return m_sum; }; ///< Sum part of the amplitude name
    /**
     * @brief Get the amplitude name in the format eJPmL
     */
    std::string get_amplitude_name() const 
    {
        return m_e_str + m_J_str + m_P_str + m_m_str + m_L_str;
    };
    /**
     * @brief Get the full amplitude name in the format reaction::sum::eJPmL
     */
    std::string get_full_amplitude() const
    {
        return m_reaction + "::" + m_sum + "::" + get_amplitude_name();
    };

    // Setters for amplitude components. Needed since we can construct amplitudes from just the quantum numbers.
    void set_amplitude_reaction(std::string reaction) { m_reaction = reaction; }; ///< Set the reaction part of the amplitude name
    void set_amplitude_sum(std::string sum) { m_sum = sum; }; ///< Set the sum part of the amplitude name
    // no setter for name, since constructors handle it
};

#endif // AMPLITUDE_PARSER_H