/*Project vector-pseudoscalar PWA fit results to moments

This script uses the conversion from partial wave complex values to "project" PWA fit
results into unique moments. The file take AmpTools .fit files as input, which contain
the results of a PWA fit. The moments are computed and saved to an output csv file.

NOTE: This script assumes that the amplitudes are written in the vec-ps eJPmL format.
For example, the positive reflectivity, JP=1+, m=0, S-wave amplitude would be written
in the cfg file as [reaction]::RealNegSign::p1p0S.
*/

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <complex>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <tuple>

#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/NormIntInterface.h"
#include "AMPTOOLS_AMPS/barrierFactor.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "file_utils.h"
#include "AmplitudeParser.h"

struct Moment
{
    int alpha;
    int Jv;
    int Lambda;
    int J;
    int M;

    std::string name() const
    {
        return "H" + std::to_string(alpha) + "_" + std::to_string(Jv) + std::to_string(Lambda) + std::to_string(J) + std::to_string(M);
    }
};

// Structure to represent SDME calculation parameters as a cache key
struct SDMEKey
{
    int alpha;
    int Ji;
    int li;
    int mi;
    int Jj;
    int lj;
    int mj;

    bool operator<(const SDMEKey &other) const
    {
        return std::tie(alpha, Ji, li, mi, Jj, lj, mj) <
               std::tie(other.alpha, other.Ji, other.li, other.mi, other.Jj, other.lj, other.mj);
    }
};

// Global cache for SDME values
static std::map<SDMEKey, std::complex<double>> sdme_cache;

// forward declarations
std::vector<Moment> initialize_moments(const FitResults &results);
complex<double> calculate_moment(const Moment &moment, const FitResults &results);
complex<double> calculate_SDME(
    int alpha, int Ji, int li, int mi,
    int Jj, int lj, int mj,
    const FitResults &results);
complex<double> get_production_coefficient(
    int e, int J, int m, int L,
    const FitResults &results);
int sign(int i);
int find_max_m(const FitResults &results);
int find_max_J(const FitResults &results);
void clear_sdme_cache();
void precompute_sdme_cache(const FitResults &results);

int main(int argc, char *argv[])
{
    // Check if we have the required arguments
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_csv>" << std::endl;
        std::cerr << "  input_file: Text file containing list of .fit files" << std::endl;
        std::cerr << "  output_csv: Name of the output CSV file" << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string csv_name = argv[2];

    // input file is a text file with a list of .fit results, each on a newline.
    // load this into a vector using the utility function
    std::vector<std::string> file_vector = read_file_list(input_file);

    if (file_vector.empty())
    {
        std::cerr << "Error: Could not read file list from " << input_file << std::endl;
        return 1;
    }

    // initialize the map of moment names to their values
    std::map<std::string, complex<double>> moment_results;
    std::vector<Moment> moments; // vector of all moments to be calculated    

    // Collect all rows in a stringstream to minimize I/O operations
    std::stringstream csv_data;

    // track whether the header has been written, and the order of the columns
    bool is_header_written = false;
    std::vector<std::string> column_order;

    // ==== BEGIN FILE ITERATION ====
    for (const std::string &file : file_vector)
    {
        std::cout << "Analyzing File: " << file << "\n";
        FitResults results(file);
        if (!results.valid())
        {
            std::cout << "Invalid fit results in file: " << file << "\n";
            continue;
        }

        // before getting this file's info, clear the results from the last file
        moment_results.clear();
        moments.clear();
        clear_sdme_cache(); // Clear SDME cache for new file

        // initialize the set of moments we can have from this file's waveset
        moments = initialize_moments(results);

        if (!is_header_written)
        {
            // write the header row if this is the first file
            csv_data << "file,";
            for (auto it = moments.begin(); it != moments.end(); ++it)
            {
                csv_data << (it->name()) << "_real," << (it->name()) << "_imag";
                if (std::next(it) != moments.end())
                { // ensure no extra comma at the end
                    csv_data << ",";
                }
                column_order.push_back(it->name()); // keep track of the order of columns
            }
            csv_data << "\n"; // end header row
            is_header_written = true;
        }

        // Ensure column_order matches the moments vector for every file
        if (column_order.size() != moments.size() ||
            !std::equal(
                column_order.begin(), column_order.end(), moments.begin(),
                [](const std::string &col, const Moment &mom)
                {
                    return col == mom.name();
                }))
        {
            std::cerr
                << "Error: The set/order of moments has changed between files. "
                   "All files must have the same moments in the same order."
                << "\n";
            return 1;
        }

        // Precompute all SDME values for this file
        precompute_sdme_cache(results);

        for (const Moment &moment : moments)
        {
            // calculate the value for this moment and save it to the map
            moment_results[moment.name()] = calculate_moment(moment, results);
        }

        // ==== WRITE TO CSV ====
        // write the file name as the first column
        csv_data << file << ",";

        for (auto it = moments.begin(); it != moments.end(); ++it)
        {
            const Moment &moment = *it;
            const complex<double> &val = moment_results[moment.name()];
            csv_data << std::real(val) << "," << std::imag(val);
            if (std::next(it) != moments.end())
            {
                csv_data << ",";
            }
        }
        csv_data << "\n";

    } // end of file iteration

    // Write all collected data to the CSV file at once    
    std::ofstream csv_file(csv_name, std::ios::out | std::ios::trunc);
    if (!csv_file.is_open())
    {
        std::cerr << "Error: Could not open output file " << csv_name << " for writing." << std::endl;
        return 1;
    }
    csv_file << csv_data.str();
    csv_file.close();

    std::cout << "Projected moments written to " << csv_name << "\n";

    return 0;
}

/**
 * @brief Initialize the set of moments based on the fit results.
 *
 * @param[in] results The fit results containing the necessary data.
 * @return std::vector<Moment> A vector of all possible moments constructed from the
 *  fit results.
 */
std::vector<Moment> initialize_moments(const FitResults &results)
{
    std::vector<Moment> moments;
    int max_J = find_max_J(results);
    int max_m = find_max_m(results);

    // prepare moment quantum numbers for the moments
    std::vector<int> alpha_vector = {0, 1, 2};
    std::vector<int> Jv_vector = {0, 2}; // CGs coefficient ensure Jv=1 is always 0
    // moments with negative Lambda values are proportional to positive ones
    std::vector<int> Lambda_vector = {0, 1, 2};
    std::vector<int> J_vector;
    for (int J = 0; J <= max_J + 1; ++J) // J defined to be >= 0
    {
        J_vector.push_back(J);
    }
    std::vector<int> M_vector; // similar to lambda, (-) M values ar proportional to (+)
    for (int m = 0; m <= max_m + 1; ++m)
    {
        M_vector.push_back(m);
    }
    for (int alpha : alpha_vector)
    {
        for (int Jv : Jv_vector)
        {
            for (int Lambda : Lambda_vector)
            {
                for (int J : J_vector)
                {
                    for (int M : M_vector)
                    {
                        // FILTER non-physical moments
                        // Wigner D functions for these are always 0
                        if (M > J || Lambda > J || Lambda > Jv)
                            continue;
                        // Wishart seciton 5.10.2 proves H2(Jv,0,J,0) = 0 for any Jv, J
                        if (alpha == 2 && Lambda == 0 && M == 0)
                            continue;

                        // create a moment from these quantum numbers
                        Moment moment;
                        moment.alpha = alpha;
                        moment.Jv = Jv;
                        moment.Lambda = Lambda;
                        moment.J = J;
                        moment.M = M;

                        moments.push_back(moment);
                    }
                }
            }
        }
    }
    return moments;
}

/**
 * @brief Calculate the value of a moment based on its quantum numbers and fit results.
 *
 * @param[in] moment The moment for which to calculate the value.
 * @param[in] results The fit results containing the necessary data.
 * @return complex<double> The calculated value of the moment.
 */
complex<double> calculate_moment(const Moment &moment, const FitResults &results)
{
    complex<double> moment_value = 0.0;
    int max_J = find_max_J(results);

    // for-loops below are done to best match the mathematical definition
    for (int Ji = 0; Ji <= max_J; ++Ji)
    {
        for (int li = 0; li <= Ji + 1; ++li)
        {
            for (int Jj = 0; Jj <= max_J; ++Jj)
            {
                for (int lj = 0; lj <= Jj + 1; ++lj)
                {
                    for (int mi = -Ji; mi <= Ji; ++mi)
                    {
                        for (int mj = -Jj; mj <= Jj; ++mj)
                        {
                            complex<double> sdme = calculate_SDME(moment.alpha, Ji, li, mi, Jj, lj, mj, results);
                            if (sdme == 0.0)
                                continue; // skip lambda loops since result is 0

                            double norm = (std::sqrt(2 * li + 1) * std::sqrt(2 * lj + 1)) / (2.0 * Jj + 1);

                            for (int lambda_i = -1; lambda_i <= 1; ++lambda_i)
                            {
                                for (int lambda_j = -1; lambda_j <= 1; ++lambda_j)
                                {
                                    double clebsch = 1.0;
                                    clebsch *= clebschGordan(li, 1, 0, lambda_i, Ji, lambda_i);
                                    clebsch *= clebschGordan(lj, 1, 0, lambda_j, Jj, lambda_j);
                                    clebsch *= clebschGordan(1, moment.Jv, lambda_i, moment.Lambda, 1, lambda_j);
                                    clebsch *= clebschGordan(1, moment.Jv, 0, 0, 1, 0);
                                    clebsch *= clebschGordan(Ji, moment.J, mi, moment.M, Jj, mj);
                                    clebsch *= clebschGordan(Ji, moment.J, lambda_i, moment.Lambda, Jj, lambda_j);
                                    if (clebsch == 0.0)
                                        continue; // result for this loop is 0

                                    // add the contribution to the moment value
                                    moment_value += norm * clebsch * sdme;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    moment_value *= 0.5; // multiply by 1/2 due to double counting in SDME. This should be a more formal normalization factor
    moment_value *= ((2.0 * moment.J + 1) * (2 * moment.Jv + 1));
    return moment_value;
}

/**
 * @brief Calculate the Spin Density Matrix Element (SDME) for given quantum numbers.
 *
 * @param[in] alpha indexes the intensity term
 *  0: unpolarized, 1: polarized cos(2*Phi), 2: polarized sin(2*Phi)
 * @param[in] Ji 1st wave spin
 * @param[in] li 1st wave orbital angular momentum
 * @param[in] mi 1st wave m-projection
 * @param[in] Jj 2nd wave spin
 * @param[in] lj 2nd wave orbital angular momentum
 * @param[in] mj 2nd wave m-projection
 * @param[in] results The fit results containing the necessary data.
 * @return complex<double> The calculated SDME.
 */
complex<double> calculate_SDME(
    int alpha, int Ji, int li, int mi,
    int Jj, int lj, int mj,
    const FitResults &results)
{
    // Create cache key
    SDMEKey key = {alpha, Ji, li, mi, Jj, lj, mj};

    // Check if value is already cached
    auto cache_it = sdme_cache.find(key);
    if (cache_it != sdme_cache.end())
    {
        return cache_it->second;
    }

    // Calculate SDME value if not cached
    complex<double> sdme = 0.0;
    complex<double> c1, c2, c3, c4; // the four complex production coefficients
    double s1, s2;                  // sign factors in the SDME formula
    std::vector<double> reflectivities = {-1, 1};
    for (double e : reflectivities)
    {
        switch (alpha)
        {
        case 0:
            c1 = get_production_coefficient(e, Ji, mi, li, results);
            c2 = std::conj(get_production_coefficient(e, Jj, mj, lj, results));
            c3 = get_production_coefficient(e, Ji, -mi, li, results);
            c4 = std::conj(get_production_coefficient(e, Jj, -mj, lj, results));

            s1 = 1.0;
            s2 = sign(mi + mj + li + lj + Ji + Jj);

            sdme += c1 * c2 + s2 * c3 * c4;
            break;
        case 1:
            c1 = get_production_coefficient(e, Ji, -mi, li, results);
            c2 = std::conj(get_production_coefficient(e, Jj, mj, lj, results));
            c3 = get_production_coefficient(e, Ji, mi, li, results);
            c4 = std::conj(get_production_coefficient(e, Jj, -mj, lj, results));

            s1 = sign(1 + mi + li + Ji);
            s2 = sign(1 + mj + lj + Jj);

            sdme += e * (s1 * c1 * c2 + s2 * c3 * c4);
            break;
        case 2:
            c1 = get_production_coefficient(e, Ji, -mi, li, results);
            c2 = std::conj(get_production_coefficient(e, Jj, mj, lj, results));
            c3 = get_production_coefficient(e, Ji, mi, li, results);
            c4 = std::conj(get_production_coefficient(e, Jj, -mj, lj, results));

            s1 = sign(mi + li + Ji);
            s2 = sign(mj + lj + Jj);

            sdme += e * (s1 * c1 * c2 - s2 * c3 * c4);
            break;
        }
    }

    if (alpha == 2)
        sdme *= complex<double>(0, 1); // needs a factor of i outside the loop

    // Cache the calculated value
    sdme_cache[key] = sdme;

    return sdme;
}

/**
 * @brief Get the production coefficient for the wave's quantum numbers.
 *
 * @details
 * Since we aren't explicitly looping over the parity values, we'll need to infer them
 * from the J and L values. This function is thus hard-coded for omega-pi production
 * processes for now, and is reaction and sum independent. Any fit with multiple
 * reactions and non-constrained sums could be subject to undefined behavior.
 *
 * @param[in] e reflectivity
 * @param[in] J total angular momentum
 * @param[in] m m-projection
 * @param[in] L orbital angular momentum
 * @return The production coefficient if found, or 0 if not found.
 */
complex<double> get_production_coefficient(
    int e, int J, int m, int L,
    const FitResults &results)
{
    // determine the parity of the amplitude being requested
    int P_omega = -1;                // parity of the omega
    int P_pi = -1;                   // parity of the pion
    int P_L = (L % 2 == 0) ? 1 : -1; // parity due to relative orbital angular momentum of the 2-particle system
    int P = P_omega * P_pi * P_L;    // total parity of the system

    // construct the amplitude name in the vec-ps format
    AmplitudeParser parser(e, J, P, m, L);
    std::string amp_name_to_get = parser.get_amplitude_name();

    // we'll need to find the amplitude's value by looping through the available amplitude list
    std::vector<std::string> amp_list = results.ampList();

    // filter out isotropic background amplitudes
    amp_list.erase(
        std::remove_if(
            amp_list.begin(),
            amp_list.end(),
            [](const std::string &amp)
            {
                return amp.find("isotropic") != std::string::npos ||
                       amp.find("Background") != std::string::npos;
            }),
        amp_list.end());

    // the production coefficient from constrained sums is equal, but their
    // normalization aren't, and so we need to capture both sums
    std::vector<std::string> matching_amplitudes;
    for (const std::string &i_amp : amp_list)
    {
        AmplitudeParser i_amp_parser(i_amp);
        std::string i_amp_name = i_amp_parser.get_amplitude_name();
        if (i_amp_name == amp_name_to_get)
        {
            matching_amplitudes.push_back(i_amp);
        }
    }

    // run some error checks
    if (matching_amplitudes.empty())
    {
        return 0.0; // amplitude not found, return 0
    }
    else if (matching_amplitudes.size() > 2)
    {
        throw std::runtime_error(
            "Amplitude '" + amp_name_to_get +
            "' was found in more than two coherent sums");
    }
    if (results.scaledProductionParameter(matching_amplitudes[0]) !=
        results.scaledProductionParameter(matching_amplitudes[1]))
    {
        throw std::runtime_error(
            "Amplitude '" + amp_name_to_get +
            "' has different production parameters in coherent sums,"
            " check that it is constrained in the fit");
    }

    // get the production coefficient, and the normalization integrals from the 2 sums
    complex<double> production_coefficient = results.scaledProductionParameter(matching_amplitudes[0]);
    complex<double> sum_normalization_integrals = 0.0;

    for (const std::string &amp : matching_amplitudes)
    {
        // NOTE: currently only uses normInt, so no acceptance correction done
        AmplitudeParser amp_parser(amp);
        const NormIntInterface *norm_interface = results.normInt(amp_parser.get_amplitude_reaction());
        complex<double> N = norm_interface->normInt(amp, amp);
        sum_normalization_integrals += std::real(N);
    }

    production_coefficient *= std::sqrt(sum_normalization_integrals);

    // TODO: Temporary, unsure if required, but doing fixed values to make sure we
    // replicate the python script
    // multiply the production coefficient by the barrier factor
    double mass = 1.21;
    double pion_mass = 0.1349768; // PDG value
    double omega_mass = 0.78265;  // PDG value
    // production_coefficient *= barrierFactor(mass, L, pion_mass, omega_mass);

    return production_coefficient;
}

int sign(int i)
{
    return (i % 2 == 0) ? 1 : -1;
}

int find_max_m(const FitResults &results)
{
    int max_m = 0;
    for (const auto &reaction : results.reactionList())
    {
        for (const std::string &amplitude : results.ampList(reaction))
        {
            if (amplitude.find("isotropic") != std::string::npos ||
                amplitude.find("Background") != std::string::npos)
                continue; // skip isotropic and background amplitudes

            AmplitudeParser parser(amplitude);
            // Extract the m value from the amplitude
            int m_value = parser.get_m_int();
            if (m_value > max_m)
            {
                max_m = m_value;
            }
        }
    }

    return max_m;
}

int find_max_J(const FitResults &results)
{
    int max_J = 0;
    for (const auto &reaction : results.reactionList())
    {
        for (const std::string &amplitude : results.ampList(reaction))
        {
            if (amplitude.find("isotropic") != std::string::npos ||
                amplitude.find("Background") != std::string::npos)
                continue; // skip isotropic and background amplitudes

            AmplitudeParser parser(amplitude);
            // Extract the J value from the amplitude
            int J_value = parser.get_J_int();
            if (J_value > max_J)
            {
                max_J = J_value;
            }
        }
    }

    return max_J;
}

/**
 * @brief Clear the SDME cache.
 *
 * This should be called when switching to a new fit file to ensure
 * cached values from the previous file don't interfere.
 */
void clear_sdme_cache()
{
    sdme_cache.clear();
}

/**
 * @brief Precompute all possible SDME values for the given fit results.
 *
 * This function calculates and caches all SDME values that could be needed
 * for moment calculations, providing maximum performance at the cost of
 * upfront computation time and memory usage.
 *
 * @param[in] results The fit results to precompute SDME values for.
 */
void precompute_sdme_cache(const FitResults &results)
{
    int max_J = find_max_J(results);

    // Precompute all possible SDME values
    for (int alpha = 0; alpha <= 2; ++alpha)
    {
        for (int Ji = 0; Ji <= max_J; ++Ji)
        {
            for (int li = 0; li <= Ji + 1; ++li)
            {
                for (int Jj = 0; Jj <= max_J; ++Jj)
                {
                    for (int lj = 0; lj <= Jj + 1; ++lj)
                    {
                        for (int mi = -Ji; mi <= Ji; ++mi)
                        {
                            for (int mj = -Jj; mj <= Jj; ++mj)
                            {
                                // This will calculate and cache the SDME value
                                calculate_SDME(alpha, Ji, li, mi, Jj, lj, mj, results);
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "Precomputed " << sdme_cache.size() << " SDME values\n";
}