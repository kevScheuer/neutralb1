/*
The idea behind this file is to produce a TeX file that writes out all the 
components of the projected moments for interpretation. 

Steps: 
1. It will initially function like projected_moments.cc, by taking a fit result file and
   determining the set of moments from the waves present. A fit result path should be
   provided as a command line argument, since only one file is need to determine the
   moment expressions.

2. It will then calculate all the clebsch gordan coefficients and SDMEs. The critical
   difference is that the SDMEs will primarily be to determine the complex coefficient
   pairs c_i c_j* for each moment. We do want to check if the SDME is zero to skip
   unnecessary calculations, but we need to retain the full complex coefficient pair
   rather than just the real value of the SDME. Since the SDME is the form of
   c_i c_j* +/- c_k c_l*, and the moment is like sum CG_ijkl * SDME_ijkl, we can
   distribute the CG coefficients to each coefficient pair.

3. Within a moment, we can now group by coefficient pairs. Each unique pair c_i c_j*
   will have a sum of CG coefficients multiplying it. Those that are near zero can be
   skipped.   

4. From there, it will need to algebraically combine all c_i c_j* pairs, with their
   clebsch gordan coefficients, to determine full expressions in terms of the amplitude
   terms (|c_i|^2) and interference terms ( rho_i rho_j cos(phi_i - phi_j) and 
   rho_i rho_j sin(phi_i - phi_j) ). Amplitude terms are easy since we just check 
   if i == j. Interference terms are really from Re(c_i c_j*) and Im(c_i c_j*), which
   can be expanded to the cosine and sine forms. So we have to look for combinations
   like (c_i c_hj* + c_j c_i*) and (c_i c_j* - c_j c_i*), which correspond to 
   the real and imaginary parts respectively.
   

5. Finally, it will need to write out the TeX file with all the expressions for each
   moment. This should probably be done with all intensity terms first, then all 
   interference terms, to make it easier to read. We'll need some way to auto-wrap
   these very long equations too.

*/

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <complex>
#include <map>
#include <cmath>

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

    complex<double> coefficient; // accumulated coefficient (Clebsch-Gordans, normalization, etc.)

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

// Global caches
static std::map<CGKey, double> cg_cache;

// Threshold for considering a CG coefficient as effectively zero
static const double CG_THRESHOLD = 1e-10;

// forward declarations
std::map<Moment, ProdCoefficientPair> calculate_moment_expressions(
    const std::vector<Moment> &moments,
    const std::string &reaction,
    const FitResults &results);
std::pair<ProdCoefficientPair, ProdCoefficientPair> get_sdme_pairs(
    const int alpha, const int e, const int Ji, const int mi, const int const li,
    const int ej, const int Jj, const int mj, const int lj,
    const std::string &reaction, const FitResults &results
);
double calculate_CGs(
    int Ji, int li, int mi, int Jj, int lj, int mj,
    int lambda_i, int lambda_j, const Moment &moment);
int sign(int i);

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

    // STEPS 2 & 3: Calculate SDMEs and CG coefficients, then group by coefficient pairs
    std::map<Moment, ProdCoefficientPair> moment_expressions = calculate_moment_expressions(
        moments,
        reaction,
        results);

    std::cout << "\nMoment expressions built successfully\n";

    // // Print summary statistics
    // for (size_t i = 0; i < moment_expressions.size() && i < 5; ++i)
    // {
    //     const MomentExpression &expr = moment_expressions[i];
    //     std::cout << "  " << expr.moment.name() << ": " 
    //               << expr.coefficient_pair_cgs.size() << " unique coefficient pairs\n";
    // }
    // if (moment_expressions.size() > 5)
    // {
    //     std::cout << "  ... and " << (moment_expressions.size() - 5) << " more moments\n";
    // }

    return 0;
}


std::map<Moment, ProdCoefficientPair> calculate_moment_expressions(
    const std::vector<Moment> &moments,
    const std::string &reaction,
    const FitResults &results)
{
    std::map<Moment, ProdCoefficientPair> expressions;

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
                                ProdCoefficientPair pair1_neg, pair2_neg;
                                ProdCoefficientPair pair1_pos, pair2_pos;
                                std::tie(pair1_neg, pair2_neg) = get_sdme_pairs(
                                    moment.alpha, 
                                    -1, Ji, mi, li,
                                    -1, Jj, mj, lj,
                                    reaction, results);
                                std::tie(pair1_pos, pair2_pos) = get_sdme_pairs(
                                    moment.alpha,
                                    1, Ji, mi, li,
                                    1, Jj, mj, lj,
                                    reaction, results);

                                double norm = (std::sqrt(2 * li + 1) * std::sqrt(2 * lj + 1)) / (2.0 * Jj + 1);
                                pair1_neg.coefficient *= norm;
                                pair2_neg.coefficient *= norm;
                                pair1_pos.coefficient *= norm;
                                pair2_pos.coefficient *= norm;

                                for (int lambda_i = -1; lambda_i <= 1; ++lambda_i)
                                {
                                    for (int lambda_j = -1; lambda_j <= 1; ++lambda_j)
                                    {
                                        double clebsch = calculate_CGs(
                                            Ji, li, mi, Jj, lj, mj, lambda_i, lambda_j, moment);
                                            
                                        // TODO: from here multiply coeff by CG, and determine if sum easy sum over the coefficients
                                        // can be done. Otherwise, just store the full ProdCoefficientPair with accumulated coefficient,
                                        // so long as it is above threshold (rename threshold limit to not be CG specific anymore).

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
 * @param alpha {text}
 * @param e {text}
 * @param Ji {text}
 * @param mi {text}
 * @param li {text}
 * @param ej {text}
 * @param Jj {text}
 * @param mj {text}
 * @param lj {text}
 * @param reaction {text}
 * @param results {text}
 * @return std::pair<ProdCoefficientPair, ProdCoefficientPair> {text}
 */
std::pair<ProdCoefficientPair, ProdCoefficientPair> get_sdme_pairs(
    const int alpha, const int e, const int Ji, const int mi, const int const li,
    const int ej, const int Jj, const int mj, const int lj,
    const std::string &reaction, const FitResults &results
)
{
    // TODO: caching of the pair would go here

    std::complex<double> pair1_val, pair2_val;  // the two pairs of complex production coefficients
    double s1, s2;                              // sign factors in the SDME formula    

    // the intensity component (alpha) determines what sign factors multiply the pair
    switch(alpha)
    {
        case 0:
            pair1_val = get_production_coefficient_pair(e, Ji, mi, li, e, Jj, mj, lj, reaction, results);
            ProdCoefficientPair first_pair = {e, Ji, mi, li, e, Jj, mj, lj};
            pair2_val = get_production_coefficient_pair(e, Ji, -mi, li, e, Jj, -mj, lj, reaction, results);
            ProdCoefficientPair second_pair = {e, Ji, -mi, li, e, Jj, -mj, lj};

            s1 = 1.0;
            s2 = sign(mi + mj + li + lj + Ji + Jj);

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

            return {first_pair, second_pair}; 
        case 1:
            pair1_val = get_production_coefficient_pair(e, Ji, -mi, li, e, Jj, mj, lj, reaction, results);
            ProdCoefficientPair first_pair = {e, Ji, -mi, li, e, Jj, mj, lj};
            pair2_val = get_production_coefficient_pair(e, Ji, mi, li, e, Jj, -mj, lj, reaction, results);
            ProdCoefficientPair second_pair = {e, Ji, mi, li, e, Jj, -mj, lj};            

            s1 = sign(1 + mi + li + Ji);
            s2 = sign(1 + mj + lj + Jj);

            if (pair1_val == 0.0)
                first_pair.coefficient = 0.0;
            else
                first_pair.coefficient = e * s1;                                
            
            if (pair2_val == 0.0)
                second_pair.coefficient = 0.0;
            else
                second_pair.coefficient = e * s2;

            return {first_pair, second_pair};
        case 2:
            pair1_val = get_production_coefficient_pair(e, Ji, -mi, li, e, Jj, mj, lj, reaction, results);
            ProdCoefficientPair first_pair = {e, Ji, -mi, li, e, Jj, mj, lj};
            pair2_val = get_production_coefficient_pair(e, Ji, mi, li, e, Jj, -mj, lj, reaction, results);
            ProdCoefficientPair second_pair = {e, Ji, mi, li, e, Jj, -mj, lj};

            s1 = sign(mi + li + Ji);
            s2 = sign(mj + lj + Jj);

            if (pair1_val == 0.0)
                first_pair.coefficient = 0.0;
            else
                first_pair.coefficient = complex<double>(0, e * -1.0 * s1);
            
            if (pair2_val == 0.0)
                second_pair.coefficient = 0.0;
            else
                second_pair.coefficient = complex<double>(0, e * s2);
    }
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

    cg_cache[key] = clebsch;  // cache the calculated value

    return clebsch;
}

int sign(int i)
{
    return (i % 2 == 0) ? 1 : -1;
}
