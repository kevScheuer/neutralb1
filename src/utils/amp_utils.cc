#include <string>
#include <tuple>
#include <stdexcept>

#include "file_utils.h"

std::tuple<std::string, std::string, std::string, std::string, std::string> parse_amplitude(std::string amplitude)
{
    // Example format: eJPmL
    // where e = reflectivity, J = total angular momentum, P = parity,
    // m = m-projection, L = orbital angular momentum
    std::string e, J, P, m, L;

    // Assuming the format is always correct and well-formed
    if (amplitude.length() < 5)
    {
        throw std::invalid_argument(
            "Amplitude string is too short. Expected at least 5 characters.");
    }

    e = amplitude.substr(0, 1); // First character
    J = amplitude.substr(1, 1); // Next character
    P = amplitude.substr(2, 1); // Next character
    m = amplitude.substr(3, 1); // Next character
    L = amplitude.substr(4);    // Remaining characters

    return {e, J, P, m, L};
}