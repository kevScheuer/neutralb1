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
    std::string amp_name_to_get = parser.get_amplitude_name();
    std::string amp_conj_name_to_get = parser_conj.get_amplitude_name();

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
    // normalization integrals aren't, and so we need to capture both sums
    std::vector<std::string> matching_amplitudes;
    std::vector<std::string> matching_amplitudes_conj;
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
        if (i_amp_name == amp_conj_name_to_get)
        {
            matching_amplitudes_conj.push_back(i_amp);
        }
    }

    // run some error checks
    if (matching_amplitudes.empty() || matching_amplitudes_conj.empty())
    {
        return 0.0; // amplitude not found, return 0
    }
    if (matching_amplitudes.size() > 2)
    {
        throw std::runtime_error(
            "Amplitude '" + amp_name_to_get +
            "' was found in more than two coherent sums");
    }
    if (matching_amplitudes_conj.size() > 2)
    {
        throw std::runtime_error(
            "Amplitude '" + amp_conj_name_to_get +
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
    if (results.scaledProductionParameter(matching_amplitudes_conj[0]) !=
        results.scaledProductionParameter(matching_amplitudes_conj[1]))
    {
        throw std::runtime_error(
            "Amplitude '" + amp_conj_name_to_get +
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
                    results);
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