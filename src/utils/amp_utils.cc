#include <string>
#include <tuple>
#include <stdexcept>
#include <regex>

#include "amp_utils.h"

QuantumNumbers parse_amplitude(std::string amplitude)
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

    // where e = reflectivity, J = total angular momentum, P = parity,
    // m = m-projection, L = orbital angular momentum
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
    
    std::string e = matches[1].str(); // reflectivity
    std::string J = matches[2].str(); // total angular momentum
    std::string P = matches[3].str(); // parity
    std::string m = matches[4].str(); // m-projection
    std::string L = matches[5].str(); // orbital angular momentum

    return {e, J, P, m, L};
}

// QuantumNumbers member function implementations
int QuantumNumbers::get_e_int() const {
    if (e == "p") return 1;
    if (e == "m") return -1;
    throw std::invalid_argument("Invalid reflectivity: " + e);
}

int QuantumNumbers::get_J_int() const {
    return std::stoi(J);
}

int QuantumNumbers::get_P_int() const {
    if (P == "p") return 1;
    if (P == "m") return -1;
    throw std::invalid_argument("Invalid parity: " + P);
}

int QuantumNumbers::get_m_int() const {
    if (m == "l") return -3;
    if (m == "n") return -2;
    if (m == "m") return -1;
    if (m == "0") return 0;
    if (m == "p") return 1;
    if (m == "q") return 2;
    if (m == "r") return 3;
    throw std::invalid_argument("Invalid m-projection: " + m);
}

int QuantumNumbers::get_L_int() const {
    if (L == "S") return 0;
    if (L == "P") return 1;
    if (L == "D") return 2;
    if (L == "F") return 3;
    if (L == "G") return 4;
    // Add more as needed...

    throw std::invalid_argument("Invalid orbital angular momentum: " + L);
}

// Conversion function implementation
QuantumNumbersInt to_integers(const QuantumNumbers& qn) {
    return {
        qn.get_e_int(),
        qn.get_J_int(), 
        qn.get_P_int(),
        qn.get_m_int(),
        qn.get_L_int()
    };
}