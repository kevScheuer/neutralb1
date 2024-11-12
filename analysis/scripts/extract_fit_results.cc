/* Extract the fit results from a list of AmpTools output files and writes them to a csv
The csv file will have the following columns, and their errors if applicable:
    - eMatrixStatus
    - lastMinuitCommandStatus
    - likelihood
    - detected_events
    - generated_events
    - all AmpTools parameters
    - all production coefficients
    - all amplitude coherent sums
    - all phase differences

This script assumes the amplitudes are named in the vector-pseudoscalar style 'eJPmL'
format, where:
    e = reflectivity (either p [+] or m [-])
    J = total spin (positive integers including 0)
    P = parity (either p [+] or m [-])
    m = m-projection (either p [+], m [-], or 0)
    L = orbital angular momentum (standard letter convention: S, P, D, F, ...)
It also assumes that reflectivity sums do not mix and are constrained across sums,
meaning that phase differences can not be computed between negative and positive
reflectivity waves.


NOTE: this script is "reaction" independent, meaning if multiple reactions are in the
    fit, it assumes the phase differences are common across all of them. It also means
    that coherent sums are calculated over all reactions. This is so that multiple
    orientations, typically denoted using the "reaction", can be fit simultaneously and
    have their fit results extracted in one go.
*/

#include <cstring>
#include <fstream> // for writing csv
#include <iostream>
#include <map>
#include <sstream> // for std::istringstream
#include <string>
#include <vector>

#include "IUAmpTools/FitResults.h"

// forward declarations
void fill_maps(
    FitResults &results,
    std::map<std::string, double> &standard_results,
    std::map<std::string, std::complex<double>> &production_coefficients,
    std::map<std::string, std::map<std::string, std::vector<std::string>>> &coherent_sums,
    std::map<std::string, std::pair<std::string, std::string>> &phase_diffs);

std::tuple<std::string, std::string, std::string, std::string> parse_amplitude(std::string amplitude);

void extract_fit_results(std::string files, std::string csv_name, bool is_acceptance_corrected)
{
    // store space-separated list of files into a vector
    std::vector<std::string> file_vector;
    std::istringstream iss(files);
    std::string file;
    while (iss >> file)
    {
        file_vector.push_back(file);
    }

    // ==== MAP INITIALIZATION ====

    // map of the standard AmpTools outputs that will be common to any fit result
    std::map<std::string, double> standard_results;

    // map for all phase differences between amplitudes, whose keys are in
    // "eJPmL_eJPmL" format, and whose values are the pair of full AmpTools
    // amplitude names
    std::map<std::string, std::pair<std::string, std::string>> phase_diffs;

    // This next map stores all the different coherent sum types. The keys are the
    // coherent sum type in eJPmL format, and the values are the strings of each
    // amplitude in that type. These strings are mapped to a vector of amplitudes
    // which are the full AmpTools names of the amplitudes that match that coherent sum.
    // Any quantum numbers that are dropped from the key means they have been coherently
    // summed over.

    // EXAMPLES
    // An individual amplitude such as the positive reflectivity, JP=1+, S-wave with a
    // +1 m-projection is then stored like:
    // "eJPmL -> "p1p0S" -> {xx::ImagPosSign::p1p0S, xx::RealNegSign::p1p0S}
    // An example of a coherent sum such as the over all JP=1+ states would be:
    // "JP" -> "1p" -> {xx::ImagNegSign::m1p0S, xx::RealNegSign::p1ppD, ...}

    // Initialize the map for all coherent sum types
    std::map<std::string, std::map<std::string, std::vector<std::string>>> coherent_sums;
    std::vector<std::string> coherent_sum_types = {
        "eJPmL", // single amplitudes
        "JPmL",  // sum reflectivity
        "eJPL",  // sum m-projection
        "JPL",   // sum {reflectivity, m-projection}
        "eJP",   // sum {m-projection, angular momenta}
        "JP",    // sum {reflectivity, m-projection, angular momenta}
        "e"      // sum all except reflectivity
    };
    for (const auto &key : coherent_sum_types)
    {
        coherent_sums[key] = std::map<std::string, std::vector<std::string>>();
    }

    // lastly a map for the production coefficients (in eJPmL format)
    std::map<std::string, std::complex<double>> production_coefficients;

    // open csv file for writing
    std::ofstream csv_file;
    csv_file.open(csv_name);

    // ==== BEGIN FILE ITERATION ====
    // Iterate over each file, and add their results as a row in the csv
    for (std::string file : file_vector)
    {
        cout << "Analyzing File: " << file << "\n";
        FitResults results(file);
        if (!results.valid())
        {
            cout << "Invalid fit results in file: " << file << "\n";
            continue;
        }

        // before getting this file's info, clear the results from the last file
        standard_results.clear();
        production_coefficients.clear();
        for (auto &pair : coherent_sums)
        {
            pair.second.clear();
        }
        phase_diffs.clear();

        // fill all the maps for this file
        fill_maps(results, standard_results, production_coefficients, coherent_sums, phase_diffs);

        // == WRITE TO CSV ==
        // write the header row if this is the first file
        if (csv_file.tellp() == 0)
        {
            // 1. standard results (these already have _err values)
            for (const auto &pair : standard_results)
            {
                csv_file << pair.first << ",";
            }
            // 2. AmpTools parameter names
            for (const auto &par_name : results.parNameList())
            {
                // skip amplitude-based parameters
                if (par_name.find("::") != std::string::npos)
                {
                    continue;
                }
                csv_file << par_name << "," << par_name << "_err,";
            }
            // 3. production parameters in eJPmL_(re/im) format
            for (const auto &pair : production_coefficients)
            {
                csv_file << pair.first << "_re" << ",";
                csv_file << pair.first << "_im" << ",";
            }
            // 4. eJPmL based coherent sum titles
            for (const auto &pair : coherent_sums)
            {
                for (const auto &pair : coherent_sums[pair.first])
                {
                    csv_file << pair.first << "," << pair.first << "_err,";
                }
            }
            // 5. phase difference names in eJPmL_eJPmL format
            // use a different iterator method to avoid adding an extra comma at the end
            for (auto it = phase_diffs.begin(); it != phase_diffs.end(); ++it)
            {
                csv_file << it->first << "," << it->first << "_err";
                if (std::next(it) != phase_diffs.end())
                {
                    csv_file << ",";
                }
            }
            csv_file << "\n";
        } // end header row

        // now write the values in the same order of map loops
        // 1. standard results
        for (const auto &pair : standard_results)
        {
            csv_file << pair.second << ",";
        }
        // 2. AmpTools parameters
        for (const auto &par_name : results.parNameList())
        {
            // skip amplitude-based parameters
            if (par_name.find("::") != std::string::npos)
            {
                continue;
            }
            csv_file << results.parValue(par_name) << ",";
            csv_file << results.parError(par_name) << ",";
        }
        // 3. production parameters
        for (const auto &pair : production_coefficients)
        {
            csv_file << pair.second.real() << ",";
            csv_file << pair.second.imag() << ",";
        }
        // 4. coherent sums
        for (const auto &pair : coherent_sums)
        {
            for (const auto &sub_pair : coherent_sums[pair.first])
            {
                csv_file << results.intensity(sub_pair.second, is_acceptance_corrected).first << ",";
                csv_file << results.intensity(sub_pair.second, is_acceptance_corrected).second << ",";
            }
        }
        // 5. phase differences, again avoiding an extra comma at the end
        for (auto it = phase_diffs.begin(); it != phase_diffs.end(); ++it)
        {
            std::string phase1 = it->second.first;
            std::string phase2 = it->second.second;
            csv_file << results.phaseDiff(phase1, phase2).first << ","; // value
            csv_file << results.phaseDiff(phase1, phase2).second;       // error
            if (std::next(it) != phase_diffs.end())
            {
                csv_file << ",";
            }
        }
        csv_file << "\n"; // end of row, move on to next file
    }
}

// fill all the fit results maps for a single file
void fill_maps(
    FitResults &results,
    std::map<std::string, double> &standard_results,
    std::map<std::string, std::complex<double>> &production_coefficients,
    std::map<std::string, std::map<std::string, std::vector<std::string>>> &coherent_sums,
    std::map<std::string, std::pair<std::string, std::string>> &phase_diffs)
{
    // Store the standard AmpTools fit outputs
    standard_results["eMatrixStatus"] = results.eMatrixStatus();
    standard_results["lastMinuitCommandStatus"] = results.lastMinuitCommandStatus();
    standard_results["likelihood"] = results.likelihood();
    standard_results["detected_events"] = results.intensity(false).first;
    standard_results["detected_events_err"] = results.intensity(false).second;
    standard_results["generated_events"] = results.intensity(true).first;
    standard_results["generated_events_err"] = results.intensity(true).second;

    // fill the coherent sum and phase difference maps by iterating over all amps
    for (auto reaction : results.reactionList())
    {
        for (std::string amplitude : results.ampList(reaction))
        {
            // 'amplitude' is the full name stored by AmpTools in the format:
            // "reaction::reflectivitySum::eJPmL"

            // put isotropic background into the single amplitude category
            if (amplitude.find("Bkgd") != std::string::npos ||
                amplitude.find("iso") != std::string::npos)
            {
                coherent_sums["eJPmL"]["Bkgd"].push_back(amplitude);
                continue;
            }

            // split the "eJPmL" part of the amplitude into its components
            std::string e, JP, m, L;
            std::tie(e, JP, m, L) = parse_amplitude(amplitude);

            std::string eJPmL = e + JP + m + L;

            // store the production coefficients
            production_coefficients[eJPmL] = results.scaledProductionParameter(amplitude);

            // store the amplitudes in the coherent sum maps
            coherent_sums["eJPmL"][eJPmL].push_back(amplitude);
            coherent_sums["JPmL"][JP + m + L].push_back(amplitude);
            coherent_sums["eJPL"][e + JP + L].push_back(amplitude);
            coherent_sums["JPL"][JP + L].push_back(amplitude);
            coherent_sums["eJP"][e + JP].push_back(amplitude);
            coherent_sums["JP"][JP].push_back(amplitude);
            coherent_sums["e"][e].push_back(amplitude);

            // store the phase differences
            for (std::string pd_amplitude : results.ampList(reaction))
            {
                std::string pd_e, pd_JP, pd_m, pd_L;
                std::tie(pd_e, pd_JP, pd_m, pd_L) = parse_amplitude(pd_amplitude);
                std::string pd_eJPmL = pd_e + pd_JP + pd_m + pd_L;

                if (pd_eJPmL == eJPmL)
                    continue; // don't compare to itself
                // isotropic background cannot have a phase difference
                if (pd_amplitude.find("Bkgd") != std::string::npos ||
                    pd_amplitude.find("iso") != std::string::npos)
                {
                    continue;
                }

                // avoid duplicates due to reverse ordering of names
                if (phase_diffs.find(pd_eJPmL + "_" + eJPmL) != phase_diffs.end())
                {
                    continue;
                }

                // avoid phase differences between different reflectivities
                if (eJPmL[0] != pd_eJPmL[0])
                {
                    continue;
                }

                phase_diffs[eJPmL + "_" + pd_eJPmL] = std::make_pair(amplitude, pd_amplitude);
            }
        }
    }
}

// grab the "eJPmL" part of the amplitude and split into its components
std::tuple<std::string, std::string, std::string, std::string> parse_amplitude(std::string amplitude)
{
    // NOTE: this assumes the amplitude is named in the "eJPmL" format already, and will
    // need to be adjusted if the naming convention is different
    std::string eJPmL = amplitude.substr(amplitude.rfind("::") + 2);
    std::string e, JP, m, L;
    e = eJPmL.substr(0, 1);
    JP = eJPmL.substr(1, 2);
    m = eJPmL.substr(3, 1);
    L = eJPmL.substr(4);
    return std::make_tuple(e, JP, m, L);
}