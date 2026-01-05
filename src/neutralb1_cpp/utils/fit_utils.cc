/**
 * @file fit_utils.cc
 * @author Kevin Scheuer
 * @brief Implementation of utility functions for handling AmpTools fit results
 *
 */

#include <algorithm>
#include <iostream>

#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/NormIntInterface.h"
#include "neutralb1/AmplitudeParser.h"
#include "neutralb1/fit_utils.h"

int calculate_system_parity(int L)
{
    int P_omega = -1;                // parity of the omega
    int P_pi = -1;                   // parity of the pion
    int P_L = (L % 2 == 0) ? 1 : -1; // parity due to relative orbital angular momentum of the 2-particle system
    int P = P_omega * P_pi * P_L;    // total parity of the system
    return P;
}

std::vector<std::string> find_matching_amplitudes(
    const AmplitudeParser &amp_to_find, 
    const FitResults &results,
    const bool skip_background)
{
    std::vector<std::string> matching_amplitudes;
    std::string amp_name_to_find = amp_to_find.get_amplitude_name();

    std::vector<std::string> amp_list = results.ampList();

    // filter out isotropic background amplitudes if requested
    if (skip_background)
    {
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
    }

    // if reaction or sum is not specified, match all in that category
    bool reaction_specified = !amp_to_find.get_amplitude_reaction().empty();
    bool sum_specified = !amp_to_find.get_amplitude_sum().empty();

    for (const std::string &i_amp : amp_list)
    {
        AmplitudeParser i_amp_parser(i_amp);

        if (reaction_specified && 
            i_amp_parser.get_amplitude_reaction() != amp_to_find.get_amplitude_reaction())
        {
            continue; // skip amplitudes not in the requested reaction
        }

        if (sum_specified && 
            i_amp_parser.get_amplitude_sum() != amp_to_find.get_amplitude_sum())
        {
            continue; // skip amplitudes not in the requested sum
        }

        std::string i_amp_name = i_amp_parser.get_amplitude_name();
        if (i_amp_name == amp_name_to_find)
        {
            matching_amplitudes.push_back(i_amp);
        }
    }
    return matching_amplitudes;
}

complex<double> get_production_coefficient_pair(
    int e, int J, int m, int L,
    int e_conj, int J_conj, int m_conj, int L_conj,
    const std::string &reaction, const FitResults &results, bool acceptance_corrected)
{
    // determine the parity values for each amplitude
    int P = calculate_system_parity(L);
    int P_conj = calculate_system_parity(L_conj);

    // construct the amplitude names in the vec-ps format
    AmplitudeParser parser(e, J, P, m, L);    
    AmplitudeParser parser_conj(e_conj, J_conj, P_conj, m_conj, L_conj);
    
    // the production coefficient from constrained sums is equal, but their
    // normalization integrals aren't, so we do not specify the sums here
    parser.set_amplitude_reaction(reaction);
    parser_conj.set_amplitude_reaction(reaction);

    std::vector<std::string> matching_amplitudes = find_matching_amplitudes(
        parser, results);
    std::vector<std::string> matching_amplitudes_conj = find_matching_amplitudes(
        parser_conj, results);


    // run some error checks
    if (matching_amplitudes.empty() || matching_amplitudes_conj.empty())
    {
        return 0.0; // amplitude not found, return 0
    }
    if (matching_amplitudes.size() > 2)
    {
        throw std::runtime_error(
            "Amplitude '" + parser.get_amplitude_name() +
            "' was found in more than two coherent sums");
    }
    if (matching_amplitudes_conj.size() > 2)
    {
        throw std::runtime_error(
            "Amplitude '" + parser_conj.get_amplitude_name() +
            "' was found in more than two coherent sums");
    }
    if (results.scaledProductionParameter(matching_amplitudes[0]) !=
        results.scaledProductionParameter(matching_amplitudes[1]))
    {
        throw std::runtime_error(
            "Amplitude '" + parser.get_amplitude_name() +
            "' has different production parameters in coherent sums,"
            " check that it is constrained in the fit");
    }
    if (results.scaledProductionParameter(matching_amplitudes_conj[0]) !=
        results.scaledProductionParameter(matching_amplitudes_conj[1]))
    {
        throw std::runtime_error(
            "Amplitude '" + parser_conj.get_amplitude_name() +
            "' has different production parameters in coherent sums,"
            " check that it is constrained in the fit");
    }

    // Find the normalization integrals of the 2 constrained sums for the pair
    complex<double> sum_normalization_integrals = 0.0;
    for (const std::string &amp : matching_amplitudes)
    {
        for (const std::string &amp_conj : matching_amplitudes_conj)
        {
            AmplitudeParser amp_parser(amp);
            AmplitudeParser amp_parser_conj(amp_conj);

            if (amp_parser.get_amplitude_sum() != amp_parser_conj.get_amplitude_sum())
            {
                // only want the normalization integral from the same sum. This should
                // already be handled in it, but just in case, we skip it
                continue;
            }

            const NormIntInterface *norm_interface = results.normInt(reaction);
            complex<double> N;
            if (acceptance_corrected)
            {
                N = norm_interface->ampInt(amp, amp_conj);
            }
            else
            {
                N = norm_interface->normInt(amp, amp_conj);
            }
            sum_normalization_integrals += N;
        }
    }

    // get production coefficients with their scale values
    complex<double> prod1 = results.scaledProductionParameter(matching_amplitudes[0]);
    complex<double> prod2 = std::conj(results.scaledProductionParameter(matching_amplitudes_conj[0]));

    return prod1 * prod2 * sum_normalization_integrals;
}

double calculate_intensity(const FitResults &results, bool acceptance_corrected)
{
    complex<double> intensity = 0.0;

    for (const std::string &reaction : results.reactionList())
    {
        // first calculate the |c|^2 terms
        for (const std::string &amplitude : results.ampList(reaction))
        {
            // skip background wave
            if (amplitude.find("iso") != std::string::npos ||
                amplitude.find("Bkgd") != std::string::npos)
            {
                continue;
            }

            // skip one of each of the constrained coherent sums, since we already
            // account for it in get_production_pair
            if (amplitude.find("ImagNegSign") != std::string::npos || amplitude.find("ImagPosSign") != std::string::npos)
            {
                continue;
            }

            AmplitudeParser parser(amplitude);
            complex<double> result = get_production_coefficient_pair(
                parser.get_e_int(),
                parser.get_J_int(),
                parser.get_m_int(),
                parser.get_L_int(),
                parser.get_e_int(),
                parser.get_J_int(),
                parser.get_m_int(),
                parser.get_L_int(),
                reaction,
                results,
                acceptance_corrected);
            intensity += result;
        }

        // now the interference terms
        for (size_t i = 0; i < results.ampList(reaction).size(); ++i)
        {
            for (size_t j = i + 1; j < results.ampList(reaction).size(); ++j)
            {
                std::string amp1 = results.ampList(reaction)[i];
                std::string amp2 = results.ampList(reaction)[j];

                // skip one of each of the constrained coherent sums, since we already
                // account for it in get_production_pair
                if (amp1.find("ImagNegSign") != std::string::npos || amp1.find("ImagPosSign") != std::string::npos)
                {
                    continue;
                }
                if (amp2.find("ImagNegSign") != std::string::npos || amp2.find("ImagPosSign") != std::string::npos)
                {
                    continue;
                }

                // skip background wave
                if (amp1.find("iso") != std::string::npos ||
                    amp1.find("Bkgd") != std::string::npos)
                {
                    continue;
                }
                if (amp2.find("iso") != std::string::npos ||
                    amp2.find("Bkgd") != std::string::npos)
                {
                    continue;
                }

                AmplitudeParser parser1(amp1);
                AmplitudeParser parser2(amp2);
                complex<double> result = get_production_coefficient_pair(
                    parser1.get_e_int(),
                    parser1.get_J_int(),
                    parser1.get_m_int(),
                    parser1.get_L_int(),
                    parser2.get_e_int(),
                    parser2.get_J_int(),
                    parser2.get_m_int(),
                    parser2.get_L_int(),
                    reaction,
                    results,
                    acceptance_corrected);
                intensity += 2 * std::real(result);
            }
        }
    }

    if (intensity.imag() != 0.0)
    {
        std::cerr << "Warning: Non-zero imaginary part in intensity calculation\n";
    }

    return intensity.real();
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

TString join_keys(const std::map<TString, Int_t> &m, const TString &delimiter)
{
    TString result;
    for (auto it = m.begin(); it != m.end(); ++it)
    {
        if (it != m.begin())
            result += delimiter;
        result += it->first;
    }
    return result;
}

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
    for (int m = 0; m <= 2 * max_wave_J; ++m)
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
                        // Wishart section 5.10.2 proves H2(Jv,0,J,0) = 0 for any Jv, J
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