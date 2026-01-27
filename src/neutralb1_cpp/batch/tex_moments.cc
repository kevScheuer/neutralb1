/**
 * @file tex_moments.cc
 * @author Kevin Scheuer
 * @brief Generate TeX representations of projected moment expressions.
 *
 * @details
 * This program calculates determines the set of available moments from a fit results
 * file, and writes out their algebraic expressions in LaTeX format. The coefficients
 * of the moments (normalization constants, Clebsch-Gordan coefficients, etc.) are
 * calculated and combined, resulting in expressions in terms of the amplitude
 * intensities (|c_i|^2) and interference terms (Re(c_i c_j*).
 *
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/NormIntInterface.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "neutralb1/file_utils.h"
#include "neutralb1/AmplitudeParser.h"
#include "neutralb1/fit_utils.h"

/**
 * @brief Structure to uniquely identify a production coefficient pair
 *
 * This struct is used to group contributions to moments by their complex production
 * coefficients pairings, which are of the form c_i c_j*. Here we access ci
 * quantum numbers (e, J, l, m) and cj quantum numbers (e_conj, J_conj, l_conj, m_conj).
 *
 */
struct ProdCoefficientPair
{
    // Unique quantum number identifiers for the two waves
    int e, J, m, l;
    int e_conj, J_conj, m_conj, l_conj;

    // accumulated coefficient (Clebsch-Gordans, normalization, etc.)
    complex<double> coefficient;

    bool operator<(const ProdCoefficientPair &other) const
    {
        return std::tie(e, J, m, l, e_conj, J_conj, m_conj, l_conj) <
               std::tie(other.e, other.J, other.m, other.l,
                        other.e_conj, other.J_conj, other.m_conj, other.l_conj);
    }

    bool operator==(const ProdCoefficientPair &other) const
    {
        return e == other.e && J == other.J && m == other.m && l == other.l &&
               e_conj == other.e_conj && J_conj == other.J_conj &&
               m_conj == other.m_conj && l_conj == other.l_conj;
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

struct ProdCoefficientPairKey
{
    int alpha;
    int e;
    int Ji;
    int li;
    int mi;
    int Jj;
    int lj;
    int mj;
    std::string reaction;

    bool operator<(const ProdCoefficientPairKey &other) const
    {
        return std::tie(alpha, e, Ji, li, mi, Jj, lj, mj, reaction) <
               std::tie(other.alpha, other.e, other.Ji, other.li, other.mi,
                        other.Jj, other.lj, other.mj, other.reaction);
    }
};

// Global caches
static std::map<CGKey, double> cg_cache;
static std::map<
    ProdCoefficientPairKey,
    std::pair<ProdCoefficientPair, ProdCoefficientPair>>
    prod_coeff_pair_cache;

// Threshold for considering a coefficient as effectively zero
static const double THRESHOLD = 1e-10;

// forward declarations
std::map<Moment, std::vector<ProdCoefficientPair>> calculate_moment_expressions(
    const std::vector<Moment> &moments,
    const std::string &reaction,
    const FitResults &results);
std::pair<ProdCoefficientPair, ProdCoefficientPair> get_sdme_pairs(
    const int alpha, const int e,
    const int Ji, const int mi, const int li,
    const int Jj, const int mj, const int lj,
    const std::string &reaction, const FitResults &results);
std::map<Moment, std::vector<ProdCoefficientPair>> combine_like_terms(
    const std::map<Moment, std::vector<ProdCoefficientPair>> &input_expressions);
bool verify_reflectivities(
    const std::map<Moment, std::vector<ProdCoefficientPair>> &combined_moment_expressions);
void write_moment_expressions_tex(
    const std::map<Moment,
                   std::vector<ProdCoefficientPair>> &combined_moment_expressions,
    const std::string &filename);
double calculate_CGs(
    int Ji, int li, int mi, int Jj, int lj, int mj,
    int lambda_i, int lambda_j, const Moment &moment);
int sign(int i);
void clear_prod_coeff_pair_cache();

int main(int argc, char *argv[])
{
    // STEP 1: Parse command line and initialize moments
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <fit_result_file>" << std::endl;
        std::cerr << "  fit_result_file: Path to AmpTools .fit file" << std::endl;
        return 1;
    }

    std::string fit_file = argv[1];

    // Load fit results
    FitResults results(fit_file);
    if (!results.valid())
    {
        std::cerr << "Error: Could not load fit results from " << fit_file << std::endl;
        return 1;
    }

    std::cout << "Loaded fit results from: " << fit_file << std::endl;

    // Determine the set of all possible moments according to the waves present
    std::vector<Moment> moments = initialize_moments(results);
    std::cout << "Found " << moments.size() << " possible moments\n";

    // Get the first reaction (for now, assuming amplitudes are common across reactions)
    std::vector<std::string> reactions = results.reactionList();
    if (reactions.empty())
    {
        std::cerr << "Error: No reactions found in fit results\n";
        return 1;
    }
    std::string reaction = reactions[0];

    // Find all non-zero production coefficient pairs for each moment
    std::map<Moment, std::vector<ProdCoefficientPair>> moment_expressions =
        calculate_moment_expressions(
            moments,
            reaction,
            results);
    std::cout << "Moment expressions built successfully\n";

    // Combine the coefficients of all like pairs
    std::map<Moment, std::vector<ProdCoefficientPair>> combined_moment_expressions =
        combine_like_terms(moment_expressions);
    if (!verify_reflectivities(combined_moment_expressions))
    {
        std::cerr << "Error: Reflectivity terms verification failed\n";
        return 1;
    }
    // since reflectivity terms are correct, we can drop half of them since we'll write
    // it out as a sum
    for (auto &entry : combined_moment_expressions)
    {
        std::vector<ProdCoefficientPair> &pairs = entry.second;
        pairs.erase(
            std::remove_if(
                pairs.begin(), pairs.end(),
                [](const ProdCoefficientPair &pair)
                { return pair.e == -1; }),
            pairs.end());
    }

    std::cout << "Combined like terms in moment expressions\n";

    // At this point, we have combinations like c_i c_i* and
    // c_i c_j* +/- c_j c_i*, with respective coefficients. These can be combined and
    // written into a LaTeX format.
    write_moment_expressions_tex(
        combined_moment_expressions,
        "moment_expressions.tex");

    return 0;
}

std::map<Moment, std::vector<ProdCoefficientPair>> calculate_moment_expressions(
    const std::vector<Moment> &moments,
    const std::string &reaction,
    const FitResults &results)
{
    std::map<Moment, std::vector<ProdCoefficientPair>> expressions;

    int max_J = find_max_J(results);

    std::cout << "Processing " << moments.size() << " moments\n";

    // Process each moment
    for (const Moment &moment : moments)
    {
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
                                for (int lambda_i = -1; lambda_i <= 1; ++lambda_i)
                                {
                                    for (int lambda_j = -1; lambda_j <= 1; ++lambda_j)
                                    {
                                        // get the 4 production coefficient pairs from
                                        // the SDME
                                        ProdCoefficientPair pair1_neg, pair2_neg;
                                        ProdCoefficientPair pair1_pos, pair2_pos;
                                        std::tie(pair1_neg, pair2_neg) = get_sdme_pairs(
                                            moment.alpha, -1,
                                            Ji, mi, li,
                                            Jj, mj, lj,
                                            reaction, results);
                                        std::tie(pair1_pos, pair2_pos) = get_sdme_pairs(
                                            moment.alpha, 1,
                                            Ji, mi, li,
                                            Jj, mj, lj,
                                            reaction, results);

                                        // normalization factor
                                        double norm = (std::sqrt(2 * li + 1) * std::sqrt(2 * lj + 1)) / (2.0 * Jj + 1);
                                        pair1_neg.coefficient *= norm;
                                        pair2_neg.coefficient *= norm;
                                        pair1_pos.coefficient *= norm;
                                        pair2_pos.coefficient *= norm;

                                        // special case due to -M and -Lambda symmetry
                                        if (moment.Lambda == 0 && moment.M == 0)
                                        {
                                            pair1_neg.coefficient *= 0.5;
                                            pair2_neg.coefficient *= 0.5;
                                            pair1_pos.coefficient *= 0.5;
                                            pair2_pos.coefficient *= 0.5;
                                        }

                                        // multiply by clebsch gordan coefficient
                                        double clebsch = calculate_CGs(
                                            Ji, li, mi, Jj, lj, mj,
                                            lambda_i, lambda_j, moment);
                                        pair1_neg.coefficient *= clebsch;
                                        pair2_neg.coefficient *= clebsch;
                                        pair1_pos.coefficient *= clebsch;
                                        pair2_pos.coefficient *= clebsch;

                                        // accumulate contributions if above threshold
                                        if (std::abs(pair1_neg.coefficient) > THRESHOLD)
                                            expressions[moment].push_back(pair1_neg);
                                        if (std::abs(pair2_neg.coefficient) > THRESHOLD)
                                            expressions[moment].push_back(pair2_neg);
                                        if (std::abs(pair1_pos.coefficient) > THRESHOLD)
                                            expressions[moment].push_back(pair1_pos);
                                        if (std::abs(pair2_pos.coefficient) > THRESHOLD)
                                            expressions[moment].push_back(pair2_pos);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return expressions;
}

/**
 * @brief Get the two production coefficient pairs from an SDME calculation.
 *
 * @details
 * An SDME is typically of the form sum_e(s_1 c_i c_j* +/- s_2 c_k c_l*), where each c_x
 * is a complex production coefficient, characterized by quantum numbers (e, J, m, l),
 * and s1 and s2 are sign factors. This function extracts the two unique coefficient
 * pairs (c_i c_j* and c_k c_l*) involved in the SDME, for a particular set of quantum
 * numbers. In project_moments.cc, the reflectivities are summed over so are not
 * specified, but here we explicitly include them to identify each pairing.
 *
 * @param alpha intensity component (0, 1, or 2)
 * @param e reflectivity of first wave
 * @param Ji total spin of first wave
 * @param mi spin-projection of first wave
 * @param li orbital angular momentum of first wave
 * @param ej reflectivity of second wave
 * @param Jj total spin of second wave
 * @param mj spin-projection of second wave
 * @param lj orbital angular momentum of second wave
 * @param reaction The reaction string. This assumes all waves are present across all
 *  reactions, so only one reaction needs to be given here.
 * @param results The fit results containing production coefficients
 * @return std::pair<ProdCoefficientPair, ProdCoefficientPair> The two unique production coefficient pairs
 */
std::pair<ProdCoefficientPair, ProdCoefficientPair> get_sdme_pairs(
    const int alpha, const int e,
    const int Ji, const int mi, const int li,
    const int Jj, const int mj, const int lj,
    const std::string &reaction, const FitResults &results)
{
    // Create cache key
    ProdCoefficientPairKey key = {alpha, e, Ji, li, mi, Jj, lj, mj, reaction};

    // Check if value is already cached
    auto cache_it = prod_coeff_pair_cache.find(key);
    if (cache_it != prod_coeff_pair_cache.end())
    {
        return cache_it->second;
    }

    // the intensity component (alpha) determines what sign factors multiply the pair
    switch (alpha)
    {
    case 0:
    {
        std::complex<double> pair1_val = get_production_coefficient_pair(e, Ji, mi, li, e, Jj, mj, lj, reaction, results);
        ProdCoefficientPair first_pair = {e, Ji, mi, li, e, Jj, mj, lj, 0.0};
        std::complex<double> pair2_val = get_production_coefficient_pair(e, Ji, -mi, li, e, Jj, -mj, lj, reaction, results);
        ProdCoefficientPair second_pair = {e, Ji, -mi, li, e, Jj, -mj, lj, 0.0};

        double s1 = 1.0;
        double s2 = sign(mi + mj + li + lj + Ji + Jj);

        // if either amplitude was not found, the pair value will be 0. In this
        // case, we set the coefficient to zero so that it is removed later
        if (pair1_val == 0.0)
            first_pair.coefficient = 0.0;
        else
            first_pair.coefficient = s1;

        if (pair2_val == 0.0)
            second_pair.coefficient = 0.0;
        else
            second_pair.coefficient = s2;

        auto result = std::make_pair(first_pair, second_pair);
        prod_coeff_pair_cache[key] = result;
        return result;
    }
    case 1:
    {
        std::complex<double> pair1_val = get_production_coefficient_pair(e, Ji, -mi, li, e, Jj, mj, lj, reaction, results);
        ProdCoefficientPair first_pair = {e, Ji, -mi, li, e, Jj, mj, lj, 0.0};
        std::complex<double> pair2_val = get_production_coefficient_pair(e, Ji, mi, li, e, Jj, -mj, lj, reaction, results);
        ProdCoefficientPair second_pair = {e, Ji, mi, li, e, Jj, -mj, lj, 0.0};

        double s1 = sign(1 + mi + li + Ji);
        double s2 = sign(1 + mj + lj + Jj);

        if (pair1_val == 0.0)
            first_pair.coefficient = 0.0;
        else
            first_pair.coefficient = e * s1;

        if (pair2_val == 0.0)
            second_pair.coefficient = 0.0;
        else
            second_pair.coefficient = e * s2;

        auto result = std::make_pair(first_pair, second_pair);
        prod_coeff_pair_cache[key] = result;
        return result;
    }
    case 2:
    {
        std::complex<double> pair1_val = get_production_coefficient_pair(e, Ji, -mi, li, e, Jj, mj, lj, reaction, results);
        ProdCoefficientPair first_pair = {e, Ji, -mi, li, e, Jj, mj, lj, 0.0};
        std::complex<double> pair2_val = get_production_coefficient_pair(e, Ji, mi, li, e, Jj, -mj, lj, reaction, results);
        ProdCoefficientPair second_pair = {e, Ji, mi, li, e, Jj, -mj, lj, 0.0};

        double s1 = sign(mi + li + Ji);
        double s2 = sign(mj + lj + Jj);
        if (pair1_val == 0.0)
            first_pair.coefficient = 0.0;
        else
            first_pair.coefficient = complex<double>(0, e * -1.0 * s1);

        if (pair2_val == 0.0)
            second_pair.coefficient = 0.0;
        else
            second_pair.coefficient = complex<double>(0, e * s2);
        auto result = std::make_pair(first_pair, second_pair);
        prod_coeff_pair_cache[key] = result;
        return result;
    }
    default:
        // Should never reach here if alpha is valid (0, 1, or 2)
        throw std::invalid_argument("Invalid alpha value in get_sdme_pairs: " + std::to_string(alpha));
    }
}

/**
 * @brief Combine like terms in the moment expressions.
 *
 * @details
 * The moment expressions calculated in calculate_moment_expressions will contain
 * multiple instances of the same production coefficient pair, each with their own
 * accumulated coefficient. This function combines those like terms by summing their
 * coefficients. Any terms with coefficients that are effectively zero (below a defined
 * threshold) are removed from the final expressions.
 *
 * @param input_expressions moment and their associated production coefficient pairs
 * @return std::map<Moment, std::vector<ProdCoefficientPair>> combined moment
 *  expressions with like terms combined
 */
std::map<Moment, std::vector<ProdCoefficientPair>> combine_like_terms(
    const std::map<Moment, std::vector<ProdCoefficientPair>> &input_expressions)
{
    std::map<Moment, std::vector<ProdCoefficientPair>> combined_expressions;

    for (const auto &entry : input_expressions)
    {
        const Moment &moment = entry.first;
        const std::vector<ProdCoefficientPair> &pairs = entry.second;

        std::map<ProdCoefficientPair, complex<double>> temp_map;

        // Combine like terms
        for (const auto &pair : pairs)
        {
            temp_map[pair] += pair.coefficient;
        }

        // Filter out near-zero coefficients and store in combined_expressions
        for (const auto &temp_entry : temp_map)
        {
            if (std::abs(temp_entry.second) > THRESHOLD)
            {
                ProdCoefficientPair combined_pair = temp_entry.first;
                combined_pair.coefficient = temp_entry.second;
                combined_expressions[moment].push_back(combined_pair);
            }
        }
    }

    return combined_expressions;
}

/**
 * @brief Write the moment expressions to a TeX file.
 *
 * @details
 * All intensity terms |c_i|^2 are written first, followed by all interference terms
 * Re(c_i c_j*). Checks are in place to ensure that coefficients are purely real (H0, H1)
 * or purely imaginary (H2) as expected, and that every interference term has a matching
 * conjugate pair to form the real part. It also ensures these pairs have matching
 * coefficients. Any discrepancies will throw runtime errors.
 *
 * @param combined_moment_expressions map of moments to their combined production
 *  coefficient pairs
 * @param filename output TeX file path
 */
void write_moment_expressions_tex(
    const std::map<Moment,
                   std::vector<ProdCoefficientPair>> &combined_moment_expressions,
    const std::string &filename)
{
    std::map<int, char> sign_map = {
        {-1, '-'},
        {1, '+'}};
    std::map<int, char> l_int_to_char_map = {
        {0, 'S'},
        {1, 'P'},
        {2, 'D'},
        {3, 'F'},
        {4, 'G'}};

    std::ostringstream tex_content;
    tex_content << "\\documentclass{article}\n";
    tex_content << "\\usepackage{amsmath}\n";
    tex_content << "\\begin{document}\n\n\n";

    for (const auto &entry : combined_moment_expressions)
    {
        const Moment &moment = entry.first;
        // make local copy, since we'll remove items as we process it
        std::vector<ProdCoefficientPair> pairs = entry.second;

        tex_content << "\\begin{gather*}\n";
        if (moment.alpha == 0 || moment.alpha == 1)
        {
            tex_content << "H^{" << moment.alpha << "}("
                        << moment.Jv << "," << moment.Lambda << ","
                        << moment.J << "," << moment.M << ") = \\sum_\\varepsilon \n";
        }
        else if (moment.alpha == 2)
        { // alpha == 2 moments are purely imaginary, so divide by i
            tex_content << "\\frac{H^{" << moment.alpha << "}("
                        << moment.Jv << "," << moment.Lambda << ","
                        << moment.J << "," << moment.M << ")}{i}"
                        << " = \\sum_\\varepsilon \n";
        }

        // H1, H2 moments are multiplied by sign of reflectivity
        if (moment.alpha == 1 || moment.alpha == 2)
        {
            tex_content << "\\varepsilon \\bigg\\{";
        }

        // First, write the amplitude terms |c_i|^2
        int counter = 0;
        std::vector<ProdCoefficientPair> pairs_to_remove;
        for (auto it = pairs.begin(); it != pairs.end(); ++it)
        {
            const ProdCoefficientPair &pair = *it;
            if (pair.e == pair.e_conj && pair.J == pair.J_conj &&
                pair.m == pair.m_conj && pair.l == pair.l_conj)
            {
                int parity = calculate_system_parity(pair.l);

                double coeff;
                if (moment.alpha == 2)
                {
                    if (pair.coefficient.real() != 0.0)
                    {
                        throw std::runtime_error("Warning: Amplitude term with non-zero real part found for H2 moment " + moment.name() + ": " + std::to_string(pair.e) + std::to_string(pair.J) + std::to_string(pair.m) + std::to_string(pair.l) + " with coefficient " + std::to_string(pair.coefficient.real()) + " + " + std::to_string(pair.coefficient.imag()) + "i\n");
                    }
                    coeff = pair.coefficient.imag();
                }
                else
                {
                    if (pair.coefficient.imag() != 0.0)
                    {
                        throw std::runtime_error("Warning: Amplitude term with non-zero imaginary part found for moment " + moment.name() + ": " + std::to_string(pair.e) + std::to_string(pair.J) + std::to_string(pair.m) + std::to_string(pair.l) + " with coefficient " + std::to_string(pair.coefficient.real()) + " + " + std::to_string(pair.coefficient.imag()) + "i\n");
                    }
                    coeff = pair.coefficient.real();
                }

                if (it != pairs.begin() && coeff >= 0.0)
                    tex_content << "+ ";

                // Amplitude term
                if (std::abs(coeff - 1.0) < THRESHOLD)
                {
                    tex_content << "|"
                                << pair.J << "^{" << sign_map[parity] << "}"
                                << l_int_to_char_map.at(pair.l)
                                << "_{" << pair.m << "}"
                                << "^{(\\varepsilon)}"
                                << "|^2\n";
                }
                else
                {
                    tex_content << coeff << "|"
                                << pair.J << "^{" << sign_map[parity] << "}"
                                << l_int_to_char_map.at(pair.l)
                                << "_{" << pair.m << "}"
                                << "^{(\\varepsilon)}"
                                << "|^2\n";
                }

                if (counter++ % 4 == 3)
                    tex_content << "\\\\\n";

                // Mark for removal from main list
                pairs_to_remove.push_back(pair);
            }
        }

        // remove the amplitude terms we just wrote out from the main list
        for (const auto &rem_pair : pairs_to_remove)
        {
            pairs.erase(
                std::remove_if(
                    pairs.begin(), pairs.end(),
                    [&rem_pair](const ProdCoefficientPair &p)
                    {
                        return p == rem_pair;
                    }),
                pairs.end());
        }        

        // Next, write the interference terms
        std::vector<bool> used(pairs.size(), false);
        for (size_t i = 0; i < pairs.size(); ++i)
        {
            if (used[i])
                continue;

            const ProdCoefficientPair &pair = pairs[i];

            size_t conj_index = pairs.size();
            for (size_t j = i + 1; j < pairs.size(); ++j)
            {
                if (used[j])
                    continue;

                const ProdCoefficientPair &p = pairs[j];
                if (p.e == pair.e_conj && p.J == pair.J_conj &&
                    p.m == pair.m_conj && p.l == pair.l_conj &&
                    p.e_conj == pair.e && p.J_conj == pair.J &&
                    p.m_conj == pair.m && p.l_conj == pair.l)
                {
                    conj_index = j;
                    break;
                }
            }

            // all interference pairs should have a matching conjugate term, otherwise
            // our calculation failed somehow
            if (conj_index == pairs.size())
            {
                throw std::runtime_error(
                    "Could not find conjugate pair for interference term in moment " + moment.name() + ": pair " + std::to_string(pair.e) + std::to_string(pair.J) + std::to_string(pair.m) + std::to_string(pair.l) + "\n");
            }

            const ProdCoefficientPair &conj_pair = pairs[conj_index];

            // we expect conjugate pairs to have matching coefficients, otherwise
            // something has gone wrong in our calculations
            if (std::abs(pair.coefficient - conj_pair.coefficient) > THRESHOLD)
            {
                throw std::runtime_error("Warning: Interference term with unequal coefficients found for moment " + moment.name() + ": pairs " + std::to_string(pair.e) + std::to_string(pair.J) + std::to_string(pair.m) + std::to_string(pair.l) + " and " + std::to_string(conj_pair.e) + std::to_string(conj_pair.J) + std::to_string(conj_pair.m) + std::to_string(conj_pair.l) + " with coefficients " + std::to_string(pair.coefficient.real()) + " + " + std::to_string(pair.coefficient.imag()) + "i and " + std::to_string(conj_pair.coefficient.real()) + " + " + std::to_string(conj_pair.coefficient.imag()) + "i\n");
            }

            // having verified the coefficients are the same, we now have
            // C (c_i c_j* + c_j c_i*) = 2C Re(c_i c_j*)

            int parity = calculate_system_parity(pair.l);
            int parity_conj = calculate_system_parity(pair.l_conj);
            

            // H2 moments are purely imaginary
            double coeff;
            if (moment.alpha == 2)
                coeff = pair.coefficient.imag();
            else
                coeff = pair.coefficient.real();

            // negative terms already have minus sign, so write plus sign if needed.
            // First condition ensures that plus sign is written after amplitude terms
            if ((!pairs_to_remove.empty() || i > 0) && coeff >= 0.0)
                tex_content << "+ ";

            tex_content << 2.0 * coeff
                        << "\t"
                        << "\\Re\\left[ "
                        << pair.J << "^{" << sign_map[parity] << "}"
                        << l_int_to_char_map.at(pair.l)
                        << "_{" << pair.m << "}"
                        << "^{(\\varepsilon)}"
                        << pair.J_conj << "^{" << sign_map[parity_conj] << "}"
                        << l_int_to_char_map.at(pair.l_conj)
                        << "_{" << pair.m_conj << "}"
                        << "^{(\\varepsilon)\\ast}"
                        << " \\right]\n";

            if (counter++ % 2 == 1)
                tex_content << "\\\\\n";

            used[i] = true;
            used[conj_index] = true;

        } // finish loop over interference terms

        // Close big bracket for reflectivity multiplication
        if (moment.alpha == 1 || moment.alpha == 2)
        {
            tex_content << "\\bigg\\}";
        }

        tex_content << "\\end{gather*}\n\n";
    }

    tex_content << "\n\n\n\\end{document}\n";

    // write to file
    std::string output_string = tex_content.str();
    std::ofstream output_file(filename);
    if (!output_file.is_open())
    {
        std::cerr << "Error: Could not open output file " << filename
                  << " for writing." << "\n";
    }
    output_file.write(output_string.c_str(), output_string.size());
    output_file.close();

    std::cout << "Wrote moment expressions to TeX file: " << filename << "\n";

    return;
}

/**
 * @brief Check that all reflectivity pairs are present and have correct coefficients.
 *
 * @details
 * We don't need to include all reflectivity pairs in the final moment expressions,
 * since they are summed over. The get_sdme_pairs function unfortunately has to
 * explicitly include them to identify unique production coefficient pairings. So here
 * we check that every pair has its opposite reflectivity pair present, and that their
 * coefficients are correct (same for H0, opposite for H1 and H2). If any issues are
 * found, return false.
 *
 * @param combined_moment_expressions map of moments to their combined production
 *  coefficient pairs
 * @return true All reflectivity pairs verified
 * @return false Reflectivity pair verification failed
 */
bool verify_reflectivities(
    const std::map<Moment, std::vector<ProdCoefficientPair>> &combined_moment_expressions)
{
    for (const auto &entry : combined_moment_expressions)
    {
        const Moment &moment = entry.first;
        const std::vector<ProdCoefficientPair> &pairs = entry.second;

        std::vector<bool> used(pairs.size(), false);
        for (size_t i = 0; i < pairs.size(); ++i)
        {
            if (used[i])
                continue;

            const ProdCoefficientPair &pair = pairs[i];

            // find the production coefficient pair that has the opposite reflectivity
            // values
            size_t conj_index = pairs.size();
            for (size_t j = i + 1; j < pairs.size(); ++j)
            {
                if (used[j])
                    continue;

                const ProdCoefficientPair &p = pairs[j];
                if (p.e == (-1) * pair.e_conj && p.J == pair.J_conj && p.m == pair.m_conj && p.l == pair.l_conj && p.e_conj == (-1) * pair.e && p.J_conj == pair.J && p.m_conj == pair.m && p.l_conj == pair.l)
                {
                    conj_index = j;
                    break;
                }
            }

            if (conj_index == pairs.size())
            {
                std::cerr << "Could not find opposite reflectivity pair for moment "
                          << moment.name() << ": pair "
                          << std::to_string(pair.e) << std::to_string(pair.J)
                          << std::to_string(pair.m) << std::to_string(pair.l) << "\n";
                return false;
            }

            const ProdCoefficientPair &conj_pair = pairs[conj_index];

            // verify that its coefficient is the same (H0) or opposite (H1, H2) sign.
            // This also ensures that magnitudes are equal.
            if (moment.alpha == 0)
            {
                if (std::abs(pair.coefficient - conj_pair.coefficient) > THRESHOLD)
                {
                    std::cerr << "Warning: Opposite reflectivity pair with unequal coefficients found for H0 moment "
                              << moment.name() << ": pairs "
                              << std::to_string(pair.e) << std::to_string(pair.J) << std::to_string(pair.m) << std::to_string(pair.l) << " and "
                              << std::to_string(conj_pair.e) << std::to_string(conj_pair.J) << std::to_string(conj_pair.m) << std::to_string(conj_pair.l)
                              << " with coefficients "
                              << std::to_string(pair.coefficient.real()) << " + " << std::to_string(pair.coefficient.imag()) << "i and "
                              << std::to_string(conj_pair.coefficient.real()) << " + " << std::to_string(conj_pair.coefficient.imag()) << "i\n";
                    return false;
                }
            }
            else // H1 or H2
            {
                if (std::abs(pair.coefficient + conj_pair.coefficient) > THRESHOLD)
                {
                    std::cerr << "Warning: Opposite reflectivity pair with non-opposite coefficients found for H1/H2 moment "
                              << moment.name() << ": pairs "
                              << std::to_string(pair.e) << std::to_string(pair.J) << std::to_string(pair.m) << std::to_string(pair.l) << " and "
                              << std::to_string(conj_pair.e) << std::to_string(conj_pair.J) << std::to_string(conj_pair.m) << std::to_string(conj_pair.l)
                              << " with coefficients "
                              << std::to_string(pair.coefficient.real()) << " + " << std::to_string(pair.coefficient.imag()) << "i and "
                              << std::to_string(conj_pair.coefficient.real()) << " + " << std::to_string(conj_pair.coefficient.imag()) << "i\n";
                    return false;
                }
            }
            used[i] = true;
            used[conj_index] = true;
        } // finish loop over pairs
    }
    return true;
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

/**
 * @brief Clear the production coefficient pair cache.
 *
 * This should be called when switching to a new fit file to ensure
 * cached values from the previous file don't interfere.
 */
void clear_prod_coeff_pair_cache()
{
    prod_coeff_pair_cache.clear();
}
