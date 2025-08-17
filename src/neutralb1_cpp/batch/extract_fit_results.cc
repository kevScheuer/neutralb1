/* Extract the fit results from a list of AmpTools output files and writes them to a csv
The csv file will have the following columns, and their errors if applicable:
    - file name
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
#include "utils/file_utils.h" // for read_file_list
#include "utils/AmplitudeParser.h" // for parse_amplitude

// forward declarations
void fill_maps(
    FitResults &results,
    std::map<std::string, double> &standard_results,
    std::map<std::string, std::complex<double>> &production_coefficients,
    std::map<std::string, std::map<std::string, std::vector<std::string>>> &coherent_sums,
    std::map<std::string, std::pair<std::string, std::string>> &phase_diffs);

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <file_list.txt> <output.csv> [acceptance_corrected (0/1)]\n";
        return 1;
    }

    std::string file_path = argv[1];
    std::string csv_name = argv[2];
    bool is_acceptance_corrected = false;
    if (argc > 3)
    {
        is_acceptance_corrected = (std::string(argv[3]) == "1");
    }

    // file path is a text file with a list of AmpTools output files, each on a newline
    // load this into a vector using the utility function
    std::vector<std::string> file_vector = read_file_list(file_path);
    
    if (file_vector.empty()) {
        std::cerr << "Error: Could not read file list from " << file_path << std::endl;
        return 1;
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
    // A coherent sum such as the one over all JP=1+ states would be:
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

    // Collect all rows in a stringstream to minimize I/O operations
    std::stringstream csv_data;
    bool is_header_written = false;

    // ==== BEGIN FILE ITERATION ====
    // Iterate over each file, and add their results as a row in the csv
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
        if (!is_header_written)
        {
            // 1. standard results (these already have _err values)
            csv_data << "file" << ",";
            for (const auto &pair : standard_results)
            {
                csv_data << pair.first << ",";
            }
            // 2. AmpTools parameter names
            for (const auto &par_name : results.parNameList())
            {
                // skip amplitude-based parameters
                if (par_name.find("::") != std::string::npos)
                {
                    continue;
                }
                csv_data << par_name << "," << par_name << "_err,";
            }
            // 3. production parameters in eJPmL_(re/im) format
            for (const auto &pair : production_coefficients)
            {
                csv_data << pair.first << "_re" << ",";
                csv_data << pair.first << "_im" << ",";
            }
            // 4. eJPmL based coherent sum titles
            for (const auto &pair : coherent_sums)
            {
                for (const auto &sub_pair : coherent_sums[pair.first])
                {
                    csv_data << sub_pair.first << "," << sub_pair.first << "_err,";
                }
            }
            // 5. phase difference names in eJPmL_eJPmL format
            // use a different iterator method to avoid adding an extra comma at the end
            for (auto it = phase_diffs.begin(); it != phase_diffs.end(); ++it)
            {
                csv_data << it->first << "," << it->first << "_err";
                if (std::next(it) != phase_diffs.end())
                {
                    csv_data << ",";
                }
            }
            csv_data << "\n";
            is_header_written = true;
        } // end header row

        // now write the values in the same order of map loops
        // 1. standard results
        csv_data << file << ",";
        for (const auto &pair : standard_results)
        {
            csv_data << pair.second << ",";
        }
        // 2. AmpTools parameters
        for (const auto &par_name : results.parNameList())
        {
            // skip amplitude-based parameters
            if (par_name.find("::") != std::string::npos)
            {
                continue;
            }
            csv_data << results.parValue(par_name) << ",";
            csv_data << results.parError(par_name) << ",";
        }
        // 3. production parameters
        for (const auto &pair : production_coefficients)
        {
            csv_data << pair.second.real() << ",";
            csv_data << pair.second.imag() << ",";
        }
        // 4. coherent sums
        for (const auto &pair : coherent_sums)
        {
            for (const auto &sub_pair : coherent_sums[pair.first])
            {
                csv_data << results.intensity(sub_pair.second, is_acceptance_corrected).first << ",";
                csv_data << results.intensity(sub_pair.second, is_acceptance_corrected).second << ",";
            }
        }
        // 5. phase differences, again avoiding an extra comma at the end
        for (auto it = phase_diffs.begin(); it != phase_diffs.end(); ++it)
        {
            std::string phase1 = it->second.first;
            std::string phase2 = it->second.second;
            auto phase_diff = results.phaseDiff(phase1, phase2);
            csv_data << phase_diff.first << ","; // value
            csv_data << phase_diff.second;       // error
            if (std::next(it) != phase_diffs.end())
            {
                csv_data << ",";
            }
        }
        csv_data << "\n"; // end of row, move on to next file
    }

    // Write all collected data to the CSV file at once
    csv_file << csv_data.str();
    csv_file.close();

    return 0;
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
            std::string e, J, P, m, L;
            try {
                AmplitudeParser parser(amplitude);
                e = parser.get_e_str();
                J = parser.get_J_str();
                P = parser.get_P_str();
                m = parser.get_m_str();
                L = parser.get_L_str();
            } catch (const std::exception &e) {
                std::cerr << "Error parsing amplitude '" << amplitude << "': " << e.what() << std::endl;
                continue;
            }            
            std::string eJPmL = e + J + P + m + L;

            // store the production coefficients
            production_coefficients[eJPmL] = results.scaledProductionParameter(amplitude);

            // store the amplitudes in the coherent sum maps
            coherent_sums["eJPmL"][eJPmL].push_back(amplitude);
            coherent_sums["JPmL"][J + P + m + L].push_back(amplitude);
            coherent_sums["eJPL"][e + J + P + L].push_back(amplitude);
            coherent_sums["JPL"][J + P + L].push_back(amplitude);
            coherent_sums["eJP"][e + J + P].push_back(amplitude);
            coherent_sums["JP"][J + P].push_back(amplitude);
            coherent_sums["e"][e].push_back(amplitude);

            // store the phase differences
            for (std::string pd_amplitude : results.ampList(reaction))
            {
                // isotropic background cannot have a phase difference
                if (pd_amplitude.find("Bkgd") != std::string::npos ||
                    pd_amplitude.find("iso") != std::string::npos)
                {
                    continue;
                }

                std::string pd_e, pd_J, pd_P, pd_m, pd_L;
                AmplitudeParser pd_parser(pd_amplitude);
                pd_e = pd_parser.get_e_str();
                pd_J = pd_parser.get_J_str();
                pd_P = pd_parser.get_P_str();
                pd_m = pd_parser.get_m_str();
                pd_L = pd_parser.get_L_str();
                std::string pd_eJPmL = pd_e + pd_J + pd_P + pd_m + pd_L;

                if (pd_eJPmL == eJPmL)
                    continue; // don't compare to itself

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