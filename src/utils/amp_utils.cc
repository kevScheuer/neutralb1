/**
 * @file amp_utils.cc
 * @brief Implementation of utility functions for handling amplitudes in PWA fits.
 *
 */

#include <string>
#include <stdexcept>
#include <regex>

#include "amp_utils.h"

// AmplitudeParser implementation

AmplitudeParser::AmplitudeParser(const std::string &amplitude)
{
    parse_from_amplitude(amplitude);
    compute_integers();
}

AmplitudeParser::AmplitudeParser(const std::string &e, const std::string &J, const std::string &P,
                                 const std::string &m, const std::string &L)
    : e_str(e), J_str(J), P_str(P), m_str(m), L_str(L)
{
    compute_integers();
}

AmplitudeParser::AmplitudeParser(int e, int J, int P, int m, int L)
    : e_int(e), J_int(J), P_int(P), m_int(m), L_int(L)
{
    compute_strings();
}

void AmplitudeParser::parse_from_amplitude(const std::string &amplitude)
{
    // Extract amp_name from a format of reaction::sum::amp_name
    size_t last_colon = amplitude.rfind("::");
    if (last_colon == std::string::npos)
    {
        throw std::invalid_argument(
            "Amplitude string does not contain '::' separator. Ensure that the full"
            " amplitude is provided in the format 'reaction::sum::amp_name'.");
    }
    std::string amp_name = amplitude.substr(last_colon + 2);

    // Use regex to parse the components according to expected format:
    // e = reflectivity (p or m)
    // J = total angular momentum (0-9)
    // P = parity (p or m)
    // m = m-projection (l, n, m, 0, p, q, r for -3, -2, -1, 0, 1, 2, 3)
    // L = orbital angular momentum (S, P, D, F, G, H, I, etc.)

    std::regex amp_regex("([pm])([0-9])([pm])([lnm0pqr])([SPDFG])");
    std::smatch matches;

    if (!std::regex_match(amp_name, matches, amp_regex))
    {
        throw std::invalid_argument(
            "Amplitude name '" + amp_name + "' does not match expected format. "
                                            "Expected format: [pm][0-9][pm][lnm0pqr][SPDFG...] where "
                                            "first character is reflectivity (p/m), "
                                            "second is total angular momentum (0-9), "
                                            "third is parity (p/m), "
                                            "fourth is m-projection (l/n/m/0/p/q/r), "
                                            "and remaining characters are orbital angular momentum (S/P/D/F/...).");
    }

    e_str = matches[1].str(); // reflectivity
    J_str = matches[2].str(); // total angular momentum
    P_str = matches[3].str(); // parity
    m_str = matches[4].str(); // m-projection
    L_str = matches[5].str(); // orbital angular momentum
}

void AmplitudeParser::compute_integers()
{
    // Convert reflectivity
    if (e_str == "p")
        e_int = 1;
    else if (e_str == "m")
        e_int = -1;
    else
        throw std::invalid_argument("Invalid reflectivity: " + e_str);

    // Convert total angular momentum
    J_int = std::stoi(J_str);

    // Convert parity
    if (P_str == "p")
        P_int = 1;
    else if (P_str == "m")
        P_int = -1;
    else
        throw std::invalid_argument("Invalid parity: " + P_str);

    // Convert m-projection
    if (m_str == "l")
        m_int = -3;
    else if (m_str == "n")
        m_int = -2;
    else if (m_str == "m")
        m_int = -1;
    else if (m_str == "0")
        m_int = 0;
    else if (m_str == "p")
        m_int = 1;
    else if (m_str == "q")
        m_int = 2;
    else if (m_str == "r")
        m_int = 3;
    else
        throw std::invalid_argument("Invalid m-projection: " + m_str);

    // Convert orbital angular momentum
    if (L_str == "S")
        L_int = 0;
    else if (L_str == "P")
        L_int = 1;
    else if (L_str == "D")
        L_int = 2;
    else if (L_str == "F")
        L_int = 3;
    else if (L_str == "G")
        L_int = 4;
    // Add more as needed...
    else
        throw std::invalid_argument("Invalid orbital angular momentum: " + L_str);
}

void AmplitudeParser::compute_strings()
{
    // Convert reflectivity
    if (e_int == 1)
        e_str = "p";
    else if (e_int == -1)
        e_str = "m";
    else
        throw std::invalid_argument("Invalid reflectivity integer: " + std::to_string(e_int));

    // Convert total angular momentum
    if (J_int < 0 || J_int > 9)
    {
        throw std::invalid_argument("Invalid total angular momentum: " + std::to_string(J_int));
    }
    J_str = std::to_string(J_int);

    // Convert parity
    if (P_int == 1)
        P_str = "p";
    else if (P_int == -1)
        P_str = "m";
    else
        throw std::invalid_argument("Invalid parity integer: " + std::to_string(P_int));

    // Convert m-projection
    if (m_int == -3)
        m_str = "l";
    else if (m_int == -2)
        m_str = "n";
    else if (m_int == -1)
        m_str = "m";
    else if (m_int == 0)
        m_str = "0";
    else if (m_int == 1)
        m_str = "p";
    else if (m_int == 2)
        m_str = "q";
    else if (m_int == 3)
        m_str = "r";
    else
        throw std::invalid_argument("Invalid m-projection integer: " + std::to_string(m_int));

    // Convert orbital angular momentum
    if (L_int == 0)
        L_str = "S";
    else if (L_int == 1)
        L_str = "P";
    else if (L_int == 2)
        L_str = "D";
    else if (L_int == 3)
        L_str = "F";
    else if (L_int == 4)
        L_str = "G";
    // Add more as needed...
    else
        throw std::invalid_argument("Invalid orbital angular momentum integer: " + std::to_string(L_int));
}

std::string AmplitudeParser::get_amplitude_name() const
{
    return e_str + J_str + P_str + m_str + L_str;
}