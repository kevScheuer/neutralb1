/*Project vector-pseudoscalar PWA fit results to moments

This script uses the conversion from partial wave complex values to "project" PWA fit
results into unique moments. The file take AmpTools .fit files as input, which contain
the results of a PWA fit. The moments are computed and saved to an output csv file.

NOTE: This script assumes that the amplitudes are written in the vec-ps eJPmL format.
For example, the positive reflectivity, JP=1+, m=0, S-wave amplitude would be written
in the cfg file as [reaction]::RealNegSign::p1p0S. It also assumes that each amplitude's
reaction is associated with a polarization orientation. If multiple orientations are
fit, then a moment is the sum of moments calculated independently for each orientation.
*/

// TODO: add the ability to optionally printout csv of each orientation
// TODO: Drastically needs performance improvements

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <complex>
#include <map>

#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/NormIntInterface.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "neutralb1/file_utils.h"
#include "neutralb1/AmplitudeParser.h"

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

    std::string reaction;

    bool operator<(const SDMEKey &other) const
    {
        return std::tie(alpha, Ji, li, mi, Jj, lj, mj, reaction) <
               std::tie(other.alpha, other.Ji, other.li, other.mi, other.Jj, other.lj, other.mj, other.reaction);
    }
};

struct CGKey
{
    int Ji;
    int li;
    int mi;
    int Jj;
    int lj;
    int mj;
    int lambda_i;
    int lambda_j;
    int Jv;
    int Lambda;
    int M;
    int J;

    bool operator<(const CGKey &other) const
    {
        return std::tie(Ji, li, mi, Jj, lj, mj, lambda_i, lambda_j, Jv, Lambda, M, J) <
               std::tie(
                   other.Ji, other.li, other.mi, other.Jj, other.lj, other.mj,
                   other.lambda_i, other.lambda_j,
                   other.Jv, other.Lambda, other.M, other.J);
    }
};

// Global caches
static std::map<SDMEKey, std::complex<double>> sdme_cache;
static std::map<CGKey, double> cg_cache;

// forward declarations
std::vector<Moment> initialize_moments(const FitResults &results);
complex<double> calculate_moment(
    const Moment &moment, const std::string &reaction, const FitResults &results);
complex<double> calculate_SDME(
    int alpha, int Ji, int li, int mi,
    int Jj, int lj, int mj,
    const std::string &reaction, const FitResults &results);
double calculate_CGs(
    int Ji, int li, int mi, int Jj, int lj, int mj,
    int lambda_i, int lambda_j, const Moment &moment);
complex<double> get_production_coefficient(
    int e, int J, int m, int L,
    const std::string &reaction, const FitResults &results);
int sign(int i);
int find_max_J(const FitResults &results);
void clear_sdme_cache();
void precompute_caches(
    const std::vector<std::string> &reactions,
    const std::vector<Moment> &moments,
    const FitResults &results);
double calculate_intensity(const FitResults &results);

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

    // initialize maps
    std::map<std::string, complex<double>> total_moment_results;                           // sum of all moments across reactions
    std::map<std::string, std::map<std::string, complex<double>>> reaction_moment_results; // per-reaction moment results
    std::map<std::string, complex<double>> production_coefficients;                        // production coefficients for each reaction
    std::vector<std::string> reactions;                                                    // tracks reactions in a fit result
    std::vector<Moment> moments;                                                           // vector of all moments to be calculated

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
        total_moment_results.clear();
        reaction_moment_results.clear();
        reactions.clear();
        moments.clear();
        // Clear SDME cache for new file. The Clebsch Gordan cache does not need
        // clearing as it does not depend on the fit result values, only different
        // combinations of the quantum numbers, which will be common across files
        clear_sdme_cache();

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

        reactions = results.reactionList();

        // Precompute all SDME values for this file
        precompute_caches(reactions, moments, results);

        // initialize the total moment results to 0 so we can add to it in the next loop
        for (const Moment &moment : moments)
        {
            total_moment_results[moment.name()] = complex<double>(0.0, 0.0);
        }

        for (const Moment &moment : moments)
        {
            for (const std::string &reaction : reactions)
            {
                // calculate the moment for this reaction
                // right now we don't strictly need these in memory, but later
                // we'll want the ability to write out the reaction moments individually
                reaction_moment_results[moment.name()][reaction] = calculate_moment(moment, reaction, results);

                // add this reaction's moment to the total value of the moments
                total_moment_results[moment.name()] += reaction_moment_results[moment.name()][reaction];
            }
        }

        // ==== WRITE TO CSV ====
        // write the file name as the first column
        csv_data << file << ",";

        for (auto it = moments.begin(); it != moments.end(); ++it)
        {
            const Moment &moment = *it;
            const complex<double> &val = total_moment_results[moment.name()];
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
    int max_wave_J = find_max_J(results);

    // prepare moment quantum numbers for the moments
    std::vector<int> alpha_vector = {0, 1, 2};
    std::vector<int> Jv_vector = {0, 2}; // CGs coefficient ensure Jv=1 is always 0
    // moments with negative Lambda values are proportional to positive ones
    std::vector<int> Lambda_vector = {0, 1, 2};
    std::vector<int> J_vector;
    for (int J = 0; J <= 2 * max_wave_J; ++J) // J defined to be >= 0
    {
        J_vector.push_back(J);
    }
    std::vector<int> M_vector; // similar to lambda, (-) M values are proportional to (+)
    for (int m = 0; m <= max_wave_J; ++m)
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
 * @param[in] reaction The reaction string (for polarization orientation)
 * @param[in] results The fit results containing the necessary data.
 * @return complex<double> The calculated value of the moment.
 */
complex<double> calculate_moment(
    const Moment &moment, const std::string &reaction, const FitResults &results)
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
                            complex<double> sdme = calculate_SDME(moment.alpha, Ji, li, mi, Jj, lj, mj, reaction, results);
                            if (sdme == 0.0)
                                continue; // skip lambda loops since result is 0

                            double norm = (std::sqrt(2 * li + 1) * std::sqrt(2 * lj + 1)) / (2.0 * Jj + 1);

                            for (int lambda_i = -1; lambda_i <= 1; ++lambda_i)
                            {
                                for (int lambda_j = -1; lambda_j <= 1; ++lambda_j)
                                {
                                    double clebsch = calculate_CGs(Ji, li, mi, Jj, lj, mj, lambda_i, lambda_j, moment);
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

    if (moment.Lambda == 0 && moment.M == 0)
    {
        // cannot determine where this 1/2 factor is from...
        // Is definitely needed to get H0_0000 == # of events, and the values to align
        // with the fitted moments
        moment_value *= 0.5;
    }
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
 * @param[in] reaction The reaction string (for polarization orientation)
 * @param[in] results The fit results containing the necessary data.
 * @return complex<double> The calculated SDME.
 */
complex<double> calculate_SDME(
    int alpha, int Ji, int li, int mi,
    int Jj, int lj, int mj,
    const std::string &reaction, const FitResults &results)
{
    // Create cache key
    SDMEKey key = {alpha, Ji, li, mi, Jj, lj, mj, reaction};

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
            c1 = get_production_coefficient(e, Ji, mi, li, reaction, results);
            c2 = std::conj(get_production_coefficient(e, Jj, mj, lj, reaction, results));
            c3 = get_production_coefficient(e, Ji, -mi, li, reaction, results);
            c4 = std::conj(get_production_coefficient(e, Jj, -mj, lj, reaction, results));

            s1 = 1.0;
            s2 = sign(mi + mj + li + lj + Ji + Jj);

            sdme += s1 * c1 * c2 + s2 * c3 * c4;
            break;
        case 1:
            c1 = get_production_coefficient(e, Ji, -mi, li, reaction, results);
            c2 = std::conj(get_production_coefficient(e, Jj, mj, lj, reaction, results));
            c3 = get_production_coefficient(e, Ji, mi, li, reaction, results);
            c4 = std::conj(get_production_coefficient(e, Jj, -mj, lj, reaction, results));

            s1 = sign(1 + mi + li + Ji);
            s2 = sign(1 + mj + lj + Jj);

            sdme += e * (s1 * c1 * c2 + s2 * c3 * c4);
            break;
        case 2:
            c1 = get_production_coefficient(e, Ji, -mi, li, reaction, results);
            c2 = std::conj(get_production_coefficient(e, Jj, mj, lj, reaction, results));
            c3 = get_production_coefficient(e, Ji, mi, li, reaction, results);
            c4 = std::conj(get_production_coefficient(e, Jj, -mj, lj, reaction, results));

            s1 = sign(mi + li + Ji);
            s2 = sign(mj + lj + Jj);

            sdme += e * (s1 * c1 * c2 - s2 * c3 * c4);
            break;
        }
    }

    if (alpha == 2)
    {
        sdme *= complex<double>(0, 1); // needs a factor of i outside the loop
        sdme *= -1.0;
    }

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
 * @param[in] reaction The reaction string (for polarization orientation)
 * @return The production coefficient if found, or 0 if not found.
 */
complex<double> get_production_coefficient(
    int e, int J, int m, int L,
    const std::string &reaction, const FitResults &results)
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

        if (i_amp_parser.get_amplitude_reaction() != reaction)
        {
            continue; // skip amplitudes not in the requested reaction
        }

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

    return production_coefficient;
}

/**
 * @brief Calculate the set of clebsch gordan coefficients for a set of quantum numbers
 *
 * @param Ji spin 1
 * @param li angular momenta 1
 * @param mi spin projection 1
 * @param Jj spin 2
 * @param lj angular momenta 2
 * @param mj spin projection 2
 * @param lambda_i helicity 1
 * @param lambda_j helicity 2
 * @param moment the moment object containing additional quantum numbers
 * @return double the calculated Clebsch-Gordan coefficient
 */
double calculate_CGs(
    int Ji, int li, int mi, int Jj, int lj, int mj, int lambda_i, int lambda_j,
    const Moment &moment)
{

    // create cache key
    CGKey key = {Ji, li, mi, Jj, lj, mj, lambda_i, lambda_j, moment.Jv, moment.Lambda, moment.M, moment.J};

    // check if value already cached
    auto cache_it = cg_cache.find(key);
    if (cache_it != cg_cache.end())
    {
        return cache_it->second;
    }

    // Calculate the set of Clebsch-Gordan coefficients for the given quantum numbers
    double clebsch = 1.0;

    clebsch *= clebschGordan(li, 1, 0, lambda_i, Ji, lambda_i);
    clebsch *= clebschGordan(lj, 1, 0, lambda_j, Jj, lambda_j);
    clebsch *= clebschGordan(1, moment.Jv, lambda_i, moment.Lambda, 1, lambda_j);
    clebsch *= clebschGordan(1, moment.Jv, 0, 0, 1, 0);
    clebsch *= clebschGordan(Ji, moment.J, mi, moment.M, Jj, mj);
    clebsch *= clebschGordan(Ji, moment.J, lambda_i, moment.Lambda, Jj, lambda_j);

    cg_cache[key] = clebsch; // cache the calculated value

    return clebsch;
}

int sign(int i)
{
    return (i % 2 == 0) ? 1 : -1;
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
 * @brief Precompute all possible SDME and Clebsch Gordan values
 *
 * This function calculates and caches all SDME values for the fit result, and all
 * Clebsch-Gordan coefficients possible for the set of moments. This has a high cost
 * upfront but drastically reduces computation time, especially for larger moment sets.
 *
 * @param[in] reactions The list of reactions to precompute SDME values for.
 * @param[in] moments The list of moments to precompute CG values for.
 * @param[in] results The total fit result to precompute SDME values for.
 */
void precompute_caches(
    const std::vector<std::string> &reactions,
    const std::vector<Moment> &moments,
    const FitResults &results)
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
                                for (const std::string &reaction : reactions)
                                {
                                    // This will calculate and cache the SDME value
                                    calculate_SDME(alpha, Ji, li, mi, Jj, lj, mj, reaction, results);
                                }

                                for (int lambda_i = -1; lambda_i <= 1; ++lambda_i)
                                {
                                    for (int lambda_j = -1; lambda_j <= 1; ++lambda_j)
                                    {
                                        for (const Moment mom : moments)
                                        {
                                            // This will calculate and cache the CG value
                                            calculate_CGs(Ji, li, mi, Jj, lj, mj, lambda_i, lambda_j, mom);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "Precomputed " << sdme_cache.size() << " SDME values\n";
}

/**
 * @brief Calculate the intensity of the fit results.
 *
 * Useful for comparing to H0_0000 moment, or the AmpTools reported intensity, to check
 * that \ref get_production_coefficients "the production coefficient calculator" is
 * working properly.
 *
 * @param results The fit results to calculate the intensity for.
 * @return complex<double> The calculated intensity.
 */
double calculate_intensity(const FitResults &results)
{
    complex<double> intensity = 0.0;

    for (const std::string &reaction : results.reactionList())
    {
        for (const std::string &amplitude : results.ampList(reaction))
        {
            // skip background wave
            if (amplitude.find("iso") != std::string::npos ||
                amplitude.find("Bkgd") != std::string::npos)
            {
                continue;
            }

            AmplitudeParser parser(amplitude);
            complex<double> c = get_production_coefficient(
                parser.get_e_int(),
                parser.get_J_int(),
                parser.get_m_int(),
            parser.get_L_int(),
            reaction,
            results);
            intensity += c * std::conj(c);
        }
    }

    if (intensity.imag() != 0.0)
    {
        std::cerr << "Warning: Non-zero imaginary part in intensity calculation\n";
    }

    return intensity.real();
}