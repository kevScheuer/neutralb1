/**
 * @file AmplitudeParser.cc
 * @brief Implementation of the AmplitudeParser class
 *
 */

#include <string>
#include <stdexcept>
#include <regex>

#include "utils/AmplitudeParser.h"

AmplitudeParser::AmplitudeParser(const std::string &amplitude)
{
    parse_from_amplitude(amplitude);
    compute_integers();
}

AmplitudeParser::AmplitudeParser(const std::string &e, const std::string &J, const std::string &P,
                                 const std::string &m, const std::string &L)
    : m_e_str(e), m_J_str(J), m_P_str(P), m_m_str(m), m_L_str(L)
{
    compute_integers();
}

AmplitudeParser::AmplitudeParser(int e, int J, int P, int m, int L)
    : m_e_int(e), m_J_int(J), m_P_int(P), m_m_int(m), m_L_int(L)
{
    compute_strings();
}

void AmplitudeParser::parse_from_amplitude(const std::string &amplitude)
{
    // Extract reaction, sum, and amp_name from a format of reaction::sum::amp_name
    size_t first_colon = amplitude.find("::");
    size_t second_colon = amplitude.find("::", first_colon + 2);
    if (first_colon == std::string::npos || second_colon == std::string::npos)
    {
        throw std::invalid_argument(
            "Amplitude string does not contain two '::' separators. Ensure that the full"
            " amplitude is provided in the format 'reaction::sum::amp_name'.");
    }
    m_reaction = amplitude.substr(0, first_colon);
    m_sum = amplitude.substr(first_colon + 2, second_colon - (first_colon + 2));
    std::string amp_name = amplitude.substr(second_colon + 2);

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

    m_e_str = matches[1].str(); // reflectivity
    m_J_str = matches[2].str(); // total angular momentum
    m_P_str = matches[3].str(); // parity
    m_m_str = matches[4].str(); // m-projection
    m_L_str = matches[5].str(); // orbital angular momentum
}

void AmplitudeParser::compute_integers()
{
    // Convert reflectivity
    if (m_e_str == "p")
        m_e_int = 1;
    else if (m_e_str == "m")
        m_e_int = -1;
    else
        throw std::invalid_argument("Invalid reflectivity: " + m_e_str);

    // Convert total angular momentum
    m_J_int = std::stoi(m_J_str);

    // Convert parity
    if (m_P_str == "p")
        m_P_int = 1;
    else if (m_P_str == "m")
        m_P_int = -1;
    else
        throw std::invalid_argument("Invalid parity: " + m_P_str);

    // Convert m-projection
    if (m_m_str == "l")
        m_m_int = -3;
    else if (m_m_str == "n")
        m_m_int = -2;
    else if (m_m_str == "m")
        m_m_int = -1;
    else if (m_m_str == "0")
        m_m_int = 0;
    else if (m_m_str == "p")
        m_m_int = 1;
    else if (m_m_str == "q")
        m_m_int = 2;
    else if (m_m_str == "r")
        m_m_int = 3;
    else
        throw std::invalid_argument("Invalid m-projection: " + m_m_str);

    // Convert orbital angular momentum
    if (m_L_str == "S")
        m_L_int = 0;
    else if (m_L_str == "P")
        m_L_int = 1;
    else if (m_L_str == "D")
        m_L_int = 2;
    else if (m_L_str == "F")
        m_L_int = 3;
    else if (m_L_str == "G")
        m_L_int = 4;
    // Add more as needed...
    else
        throw std::invalid_argument("Invalid orbital angular momentum: " + m_L_str);
}

void AmplitudeParser::compute_strings()
{
    // Convert reflectivity
    if (m_e_int == 1)
        m_e_str = "p";
    else if (m_e_int == -1)
        m_e_str = "m";
    else
        throw std::invalid_argument("Invalid reflectivity integer: " + std::to_string(m_e_int));

    // Convert total angular momentum
    if (m_J_int < 0 || m_J_int > 9)
    {
        throw std::invalid_argument("Invalid total angular momentum: " + std::to_string(m_J_int));
    }
    m_J_str = std::to_string(m_J_int);

    // Convert parity
    if (m_P_int == 1)
        m_P_str = "p";
    else if (m_P_int == -1)
        m_P_str = "m";
    else
        throw std::invalid_argument("Invalid parity integer: " + std::to_string(m_P_int));

    // Convert m-projection
    if (m_m_int == -3)
        m_m_str = "l";
    else if (m_m_int == -2)
        m_m_str = "n";
    else if (m_m_int == -1)
        m_m_str = "m";
    else if (m_m_int == 0)
        m_m_str = "0";
    else if (m_m_int == 1)
        m_m_str = "p";
    else if (m_m_int == 2)
        m_m_str = "q";
    else if (m_m_int == 3)
        m_m_str = "r";
    else
        throw std::invalid_argument("Invalid m-projection integer: " + std::to_string(m_m_int));

    // Convert orbital angular momentum
    if (m_L_int == 0)
        m_L_str = "S";
    else if (m_L_int == 1)
        m_L_str = "P";
    else if (m_L_int == 2)
        m_L_str = "D";
    else if (m_L_int == 3)
        m_L_str = "F";
    else if (m_L_int == 4)
        m_L_str = "G";
    // Add more as needed...
    else
        throw std::invalid_argument("Invalid orbital angular momentum integer: " + std::to_string(m_L_int));
}