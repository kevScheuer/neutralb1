/**
 * @file AmplitudeParser.cc
 * @brief Implementation of the AmplitudeParser class
 *
 */

#include <map>
#include <string>
#include <stdexcept>
#include <regex>

#include "neutralb1/AmplitudeParser.h"

AmplitudeParser::AmplitudeParser(const std::string &amplitude)
{
    parse_from_amplitude(amplitude);
    compute_integers();
}

AmplitudeParser::AmplitudeParser(const char e, const char J, const char P,
                                 const char m, const char L)
    : m_e_char(e), m_J_char(J), m_P_char(P), m_m_char(m), m_L_char(L)
{
    compute_integers();
}

AmplitudeParser::AmplitudeParser(const int &e, const int &J, const int &P, const int &m, const int &L)
    : m_e_int(e), m_J_int(J), m_P_int(P), m_m_int(m), m_L_int(L)
{
    compute_chars();
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

    m_e_char = matches[1].str()[0]; // reflectivity
    m_J_char = matches[2].str()[0]; // total angular momentum
    m_P_char = matches[3].str()[0]; // parity
    m_m_char = matches[4].str()[0]; // m-projection
    m_L_char = matches[5].str()[0]; // orbital angular momentum
}

void AmplitudeParser::compute_integers()
{
    // Convert reflectivity
    if (m_e_char == 'p')
        m_e_int = 1;
    else if (m_e_char == 'm')
        m_e_int = -1;
    else
        throw std::invalid_argument("Invalid reflectivity: " + std::string(1, m_e_char));

    // Convert total angular momentum    
    m_J_int = m_J_char - '0';
    if (m_J_int < 0 || m_J_int > 9)
    {
        throw std::invalid_argument("Invalid total angular momentum: " + std::to_string(m_J_int));
    }

    // Convert parity
    if (m_P_char == 'p')
        m_P_int = 1;
    else if (m_P_char == 'm')
        m_P_int = -1;
    else
        throw std::invalid_argument("Invalid parity: " + std::string(1, m_P_char));

    // Convert m-projection
    if (m_m_char == 'l')
        m_m_int = -3;
    else if (m_m_char == 'n')
        m_m_int = -2;
    else if (m_m_char == 'm')
        m_m_int = -1;
    else if (m_m_char == '0')
        m_m_int = 0;
    else if (m_m_char == 'p')
        m_m_int = 1;
    else if (m_m_char == 'q')
        m_m_int = 2;
    else if (m_m_char == 'r')
        m_m_int = 3;
    else
        throw std::invalid_argument("Invalid m-projection: " + std::string(1, m_m_char));

    // Convert orbital angular momentum
    if (m_L_char == 'S')
        m_L_int = 0;
    else if (m_L_char == 'P')
        m_L_int = 1;
    else if (m_L_char == 'D')
        m_L_int = 2;
    else if (m_L_char == 'F')
        m_L_int = 3;
    else if (m_L_char == 'G')
        m_L_int = 4;
    // Add more as needed...
    else
        throw std::invalid_argument("Invalid orbital angular momentum: " + std::string(1, m_L_char));
}

void AmplitudeParser::compute_chars()
{
    // Convert reflectivity
    if (m_e_int == 1)
        m_e_char = 'p';
    else if (m_e_int == -1)
        m_e_char = 'm';
    else
        throw std::invalid_argument("Invalid reflectivity integer: " + std::to_string(m_e_int));

    // Convert total angular momentum
    if (m_J_int < 0 || m_J_int > 9)
    {
        throw std::invalid_argument("Invalid total angular momentum: " + std::to_string(m_J_int));
    }
    m_J_char = '0' + m_J_int;

    // Convert parity
    if (m_P_int == 1)
        m_P_char = 'p';
    else if (m_P_int == -1)
        m_P_char = 'm';
    else
        throw std::invalid_argument("Invalid parity integer: " + std::to_string(m_P_int));

    // Convert m-projection
    if (m_m_int == -3)
        m_m_char = 'l';
    else if (m_m_int == -2)
        m_m_char = 'n';
    else if (m_m_int == -1)
        m_m_char = 'm';
    else if (m_m_int == 0)
        m_m_char = '0';
    else if (m_m_int == 1)
        m_m_char = 'p';
    else if (m_m_int == 2)
        m_m_char = 'q';
    else if (m_m_int == 3)
        m_m_char = 'r';
    else
        throw std::invalid_argument("Invalid m-projection integer: " + std::to_string(m_m_int));

    // Convert orbital angular momentum
    if (m_L_int == 0)
        m_L_char = 'S';
    else if (m_L_int == 1)
        m_L_char = 'P';
    else if (m_L_int == 2)
        m_L_char = 'D';
    else if (m_L_int == 3)
        m_L_char = 'F';
    else if (m_L_int == 4)
        m_L_char = 'G';
    // Add more as needed...
    else
        throw std::invalid_argument("Invalid orbital angular momentum integer: " + std::to_string(m_L_int));
}

std::string AmplitudeParser::get_latex_amplitude() const
{
    // convert reflectivity and parity to +/- characters
    char e_sign = (m_e_int == 1) ? '+' : '-';
    char P_sign = (m_P_int == 1) ? '+' : '-';    
    
    // Format: J^{P}L_m^e
    std::ostringstream oss;
    oss << m_J_char << "^{" << P_sign << "}" << m_L_char << "_{" << m_m_int << "}^{" << e_sign << "}";
    return oss.str();
}
