/* Read .fit files and output csv containing coherent sums and parameters.

   Two main functionalities are to combine many .fit files in one
   directory, or combine .fit files with a particular name in many
   subdirectories. Execute script with "-h" or "--help" to see options

 NOTE: If adding columns, be very careful since the header line and value line are
 separate, and so they must be done in the same order in the code

 TODO: Save the real and imaginary parts of every amplitude to columns. Practically with
 the full amp name I can do name=real(imag)ProdParName(full_amp_name), then
 parValue(name), parError(name) to obtain the information. This will be nice to
 histogram for bootstrap fits to see if these actually are gaussian shaped
 */

#include <algorithm>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "IUAmpTools/FitResults.h"
#include "TMath.h"

// Forward declarations of functions. Comments in main should hopefully make clear how
// the program executes.
void GetListOfFitFiles(
    std::vector<std::string> &file_list,
    std::string parent_dir,
    std::string search_for_file);
float GetLastNumOfString(std::string str);
bool NaturalComp(std::string s1, string s2);
std::string ConvertFullAmplitudeName(std::string full_amp);
void WriteToCsv(
    ofstream &csv_file,
    FitResults &results,
    std::map<std::string, std::vector<std::string>> &mapCohSum,
    bool is_header,
    bool is_acceptance_corrected);
void WriteToCsv(
    ofstream &csv_file,
    FitResults &results,
    std::map<std::string, std::pair<std::string, std::string>> &mapCohSum,
    bool is_header,
    bool is_acceptance_corrected);

void fitsToCsv(std::string args = "")
{
    // initialize default arguments
    std::string search_for_file = "";
    std::string parent_dir = "./";
    std::string output_csv_name = "combinedFitPars.csv";
    bool is_acceptance_corrected = false;

    // ==== PARSE COMMAND LINE ARGS ====
    size_t flag_pos = 0;
    size_t arg_pos = 0;
    while (args.length() != 0)
    {
        // get positions of flag (-flag) and argument that follows the flag
        flag_pos = args.find(" ");
        arg_pos = args.find(" ", flag_pos + 1);

        if (flag_pos == std::string::npos && (args != "-h" && args != "--help"))
        {
            cout << "Invalid input. Use \"--help\" to see options" << "\n";
            exit(1);
        }

        std::string flag = args.substr(0, flag_pos);
        std::string arg = args.substr(flag_pos + 1, arg_pos - flag_pos - 1);

        // check against flags and arguments
        if (flag == "-f" || flag == "--file")
        {
            search_for_file = arg;
            if (search_for_file.substr(search_for_file.length() - 4) != ".fit")
            {
                search_for_file += ".fit";
            }
        }
        if (flag == "-o" || flag == "--output")
        {
            output_csv_name = arg;
            if (output_csv_name.substr(output_csv_name.length() - 4) != ".csv")
            {
                output_csv_name += ".csv";
            }
        }
        if (flag == "-d" || flag == "--directory")
        {
            parent_dir = arg;
            if (parent_dir.back() != '/')
                parent_dir += '/';
        }
        if (flag == "-a" || flag == "--acceptance_corrected")
        {
            // convert string into bool and assign
            std::istringstream(arg) >> std::boolalpha >> is_acceptance_corrected;
        }
        if (flag == "-h" || flag == "--help")
        {
            cout << "usage: root -l 'fitsToCsv.C(\"[args]\")' [-h] [-f FILE_NAME] "
                 << "[-o OUTPUT_FILE_NAME] [-d PARENT_DIRECTORY] [-a BOOL]"
                 << "\n\n"
                 << "Aggregate results from many .fit files into a CSV.\n"
                 << "Running without argument combines .fit files in current dir.\n"
                 << " Files will be sorted by the last integer in the full path and\n"
                 << " the csv's index will match this sorting.\n"
                 << "\n\n"
                 << "optional arguments:\n"
                 << "\t-h, --help"
                 << "\n\t\tshow this help message and exit\n"
                 << "\t-f, --file FILE_NAME\n"
                 << "\t\tcombine files with matching name in all subdirectories of "
                    "the parent directory\n"
                 << "\t-o, --output OUTPUT_FILE_NAME\n"
                 << "\t\tchange output file name (default: combinedFitPars.csv)\n"
                 << "\t-d, --directory PARENT_DIRECTORY\n"
                 << "\t\tspecify parent dir (default: ./)\n"
                 << "\t-a, --acceptance_corrected BOOL\n"
                 << "\t\tIf true, corrects amp intensities for detector effects\n"
                 << "\t\tDefaults to False i.e. uses 'detected values";
            exit(0);
        }
        args.erase(0, arg_pos);
        args.erase(0, 1);
    }
    // ==== END COMMAND LINE PARSING ====

    // make list of all .fit files
    std::vector<std::string> file_list;
    GetListOfFitFiles(file_list, parent_dir, search_for_file);

    if (file_list.size() == 0)
    {
        cout << "No files found!" << "\n";
        exit(1);
    }

    // sort file vector using last integer in element
    // This ensures the first row/index in the csv corresponds to the file with the
    // smallest last integer in it. So files can be sorted by mass bin, rand fit index,
    // t bin, etc.
    std::sort(file_list.begin(), file_list.end(), NaturalComp);

    int index = 0; // setup index and create csv file
    ofstream csv_file;
    csv_file.open(output_csv_name);

    // initialize containers
    std::vector<std::string> amp_list;
    std::vector<std::string> par_list;

    // setup map of eJPmL based header name, to a vector that stores the full
    // amplitude names with that eJPmL value. Excluded variables mean they are summed
    // ex: sum of JP=1+,l=0 waves is:
    // < 1pS , < xx::ImagNegSign::1pps, xx::RealPosSign:1pms, ... > >
    // note that all reactions (orientations) will be included in each container
    std::map<std::string, std::vector<std::string>> amp_sum_eJPmL;
    std::map<std::string, std::vector<std::string>> amp_sum_eJPmL_re;
    std::map<std::string, std::vector<std::string>> amp_sum_eJPmL_im;
    std::map<std::string, std::vector<std::string>> amp_sum_JPmL;
    std::map<std::string, std::vector<std::string>> amp_sum_JPm;
    std::map<std::string, std::vector<std::string>> amp_sum_JPL;
    std::map<std::string, std::vector<std::string>> amp_sum_eJPm;
    std::map<std::string, std::vector<std::string>> amp_sum_eJPL;
    std::map<std::string, std::vector<std::string>> amp_sum_eJP;
    std::map<std::string, std::vector<std::string>> amp_sum_JP;
    std::map<std::string, std::vector<std::string>> amp_sum_e;
    std::map<std::string, std::pair<std::string, std::string>> amp_phase_diffs; // eJPmL_eJPmL

    // ==== BEGIN FILE ITERATION ====
    // Iterate over each file, and add their results as a row in the csv
    for (std::string file : file_list)
    {
        // load .fit file
        cout << "Analyzing File: " << file << "\n";

        FitResults results(file);
        if (!results.valid())
        {
            cout << "Invalid fit results in file: " << file << "\n";
            continue;
        }

        cout << "Fit results loaded" << "\n";

        // WRITE HEADER LINE
        if (amp_list.size() == 0)
        {
            // each orientation is stored as a "reaction", so loop over to get all amps
            std::vector<std::string> reactions = results.reactionList();
            for (auto reaction : reactions)
            {
                for (std::string amp : results.ampList(reaction))
                {
                    amp_list.push_back(amp);
                }
            }
            par_list = results.parNameList();

            // write non-amplitude headers first
            csv_file << "index"
                     << ",eMatrixStatus"
                     << ",lastMinuitCommandStatus"
                     << ",likelihood"
                     << ",detected_events"
                     << ",detected_events_err"
                     << ",generated_events"
                     << ",generated_events_err";

            // write headers for AmpTools defined parameters
            for (std::string par : par_list)
            {
                // only want Parameters, not reaction-based ones
                if (par.find("::") != std::string::npos)
                {
                    continue;
                }
                csv_file << "," << par << "," << par + "_err";

                // write covariance between any D/S params if they exist
                if (par.find("dsratio") != std::string::npos)
                {
                    // get reflectivity of ratio parameter if it exists
                    std::string ratio_refl = "";
                    if (par.find("_") != std::string::npos)
                    {
                        ratio_refl = par.back();
                    }

                    if (ratio_refl != "")
                    {
                        csv_file << ","
                                 << "cov_dsratio_dphase_" << ratio_refl;
                    }
                    else
                    {
                        csv_file << ","
                                 << "cov_dsratio_dphase";
                    }
                }
            }

            // add full amplitude name to corresponding eJPmL key
            for (unsigned int i = 0; i < amp_list.size(); i++)
            {
                std::string full_amp = amp_list[i];

                // handle background amplitude separately
                if (full_amp.find("Bkgd") != std::string::npos)
                {
                    amp_sum_eJPmL["Bkgd"].push_back(full_amp);
                    continue;
                }

                std::string eJPmL = ConvertFullAmplitudeName(full_amp);
                std::string e = eJPmL.substr(0, 1), JP = eJPmL.substr(1, 2),
                            m = eJPmL.substr(3, 1), L = eJPmL.substr(4, 1);

                amp_sum_eJPmL[eJPmL].push_back(full_amp);
                amp_sum_eJPmL_re[eJPmL + "_re"].push_back(full_amp);
                amp_sum_eJPmL_im[eJPmL + "_im"].push_back(full_amp);
                amp_sum_JPmL[JP + m + L].push_back(full_amp);
                amp_sum_JPm[JP + m].push_back(full_amp);
                amp_sum_JPL[JP + L].push_back(full_amp);
                amp_sum_eJPm[e + JP + m].push_back(full_amp);
                amp_sum_eJPL[e + JP + L].push_back(full_amp);
                amp_sum_eJP[e + JP].push_back(full_amp);
                amp_sum_JP[JP].push_back(full_amp);
                amp_sum_e[e].push_back(full_amp);

                // second loop to get phase difference headers
                for (unsigned int j = i + 1; j < amp_list.size(); j++)
                {
                    std::string pd_full_amp = amp_list[j];

                    if (pd_full_amp.find("isotropic") != std::string::npos)
                    {
                        continue; // avoid trying to get phase diff with background
                    }

                    std::string pd_eJPmL = ConvertFullAmplitudeName(pd_full_amp);

                    // only make phase differences between the same coherent sum
                    if (pd_eJPmL.substr(0, 1) != eJPmL.substr(0, 1))
                    {
                        continue;
                    }
                    // avoid making phase differences between same amplitudes from
                    // different reactions (they're constrained to be the same
                    // at the config level)
                    if (pd_eJPmL == eJPmL)
                    {
                        continue;
                    }
                    // avoid writing the reverse ordering of the phase difference.
                    // Otherwise we end up getting two columns of essentially same data
                    if (amp_phase_diffs.find(pd_eJPmL + "_" + eJPmL) != amp_phase_diffs.end())
                    {
                        continue;
                    }

                    amp_phase_diffs[eJPmL + "_" + pd_eJPmL] = std::make_pair(full_amp, pd_full_amp);
                }
            } // end mapping for-loop

            // write amplitude headers
            WriteToCsv(csv_file, results, amp_sum_eJPmL, true, is_acceptance_corrected);
            WriteToCsv(csv_file, results, amp_sum_eJPmL_re, true, is_acceptance_corrected);
            WriteToCsv(csv_file, results, amp_sum_eJPmL_im, true, is_acceptance_corrected);
            WriteToCsv(csv_file, results, amp_sum_JPmL, true, is_acceptance_corrected);
            WriteToCsv(csv_file, results, amp_sum_JPm, true, is_acceptance_corrected);
            WriteToCsv(csv_file, results, amp_sum_JPL, true, is_acceptance_corrected);
            WriteToCsv(csv_file, results, amp_sum_eJPm, true, is_acceptance_corrected);
            WriteToCsv(csv_file, results, amp_sum_eJPL, true, is_acceptance_corrected);
            WriteToCsv(csv_file, results, amp_sum_eJP, true, is_acceptance_corrected);
            WriteToCsv(csv_file, results, amp_sum_JP, true, is_acceptance_corrected);
            WriteToCsv(csv_file, results, amp_sum_e, true, is_acceptance_corrected);
            WriteToCsv(csv_file, results, amp_phase_diffs, true, is_acceptance_corrected);

            csv_file
                << "\n";
        } // end header line

        // ==== WRITE VALUES TO CSV !! MUST BE IN SAME ORDER AS HEADERS !! ====
        double detected_events = results.intensity(false).first;
        double detected_events_err = results.intensity(false).second;
        double generated_events = results.intensity().first;
        double generated_events_err = results.intensity().second;

        if (TMath::IsNaN(detected_events))
        {
            detected_events = 0;
        }
        if (TMath::IsNaN(detected_events_err))
        {
            detected_events_err = 0;
        }
        if (TMath::IsNaN(generated_events))
        {
            generated_events = 0;
        }
        if (TMath::IsNaN(generated_events_err))
        {
            generated_events_err = 0;
        }

        csv_file << index;
        csv_file << "," << results.eMatrixStatus();
        csv_file << "," << results.lastMinuitCommandStatus();
        csv_file << "," << results.likelihood();
        csv_file << "," << detected_events
                 << "," << detected_events_err;
        csv_file << "," << generated_events
                 << "," << generated_events_err;

        for (std::string par : par_list)
        {
            // only want Parameters, not reaction-based ones
            if (par.find("::") != std::string::npos)
            {
                continue;
            }
            csv_file << "," << results.parValue(par)
                     << "," << results.parError(par);

            // write D/S ratio phase covariance if it exists
            if (par.find("dsratio") != std::string::npos)
            {
                // get refl of ratio param if exists
                std::string ratio_refl = "";
                if (par.find("_") != std::string::npos)
                {
                    ratio_refl = par.back();
                }

                if (ratio_refl != "")
                {
                    std::string r = "dsratio_" + ratio_refl;
                    std::string ph = "dphase_" + ratio_refl;
                    csv_file << "," << results.covariance(r, ph);
                }
                else
                {
                    csv_file << "," << results.covariance("dsratio", "dphase");
                }
            }
        }

        // finally write the sums
        // MUST be in same order as for-loops when writing headers
        WriteToCsv(csv_file, results, amp_sum_eJPmL, false, is_acceptance_corrected);
        WriteToCsv(csv_file, results, amp_sum_eJPmL_re, false, is_acceptance_corrected);
        WriteToCsv(csv_file, results, amp_sum_eJPmL_im, false, is_acceptance_corrected);
        WriteToCsv(csv_file, results, amp_sum_JPmL, false, is_acceptance_corrected);
        WriteToCsv(csv_file, results, amp_sum_JPm, false, is_acceptance_corrected);
        WriteToCsv(csv_file, results, amp_sum_JPL, false, is_acceptance_corrected);
        WriteToCsv(csv_file, results, amp_sum_eJPm, false, is_acceptance_corrected);
        WriteToCsv(csv_file, results, amp_sum_eJPL, false, is_acceptance_corrected);
        WriteToCsv(csv_file, results, amp_sum_eJP, false, is_acceptance_corrected);
        WriteToCsv(csv_file, results, amp_sum_JP, false, is_acceptance_corrected);
        WriteToCsv(csv_file, results, amp_sum_e, false, is_acceptance_corrected);
        WriteToCsv(csv_file, results, amp_phase_diffs, false, is_acceptance_corrected);

        csv_file << "\n";
        index += 1;
    }

    cout << "All Fit Results successfully written to " << output_csv_name
         << "\n";
}

// Returns list of fit files. If search_for_file is non-empty, it will search
// every subdir of the parent_dir for .fit files matching "search_for_file"
void GetListOfFitFiles(std::vector<std::string> &file_list,
                       std::string parent_dir, std::string search_for_file)
{
    DIR *dir;
    struct dirent *ent;

    // Open parent directory
    dir = opendir(parent_dir.c_str());
    if (dir == NULL)
    {
        cout << "Directory " << parent_dir << " could not be opened"
             << "\n";
        exit(1);
    }

    // read directory (directories) and add .fit files to file_list
    while ((ent = readdir(dir)))
    {
        if (ent->d_name[0] == '.')
        {
            continue;
        }
        // if searching for matching file, then recursively search the subdirs
        if (ent->d_type == DT_DIR && search_for_file != "")
        {
            std::string path = parent_dir + ent->d_name + '/';
            GetListOfFitFiles(file_list, path, search_for_file);
        }

        std::string file = ent->d_name;
        if (file.size() < 4)
        {
            continue;
        }
        if (file.substr(file.size() - 4) != ".fit")
        {
            continue;
        }
        if (search_for_file != "" && file != search_for_file)
        {
            continue;
        }
        file_list.push_back(parent_dir + file);
    }
}

// Returns last integer in string. If none is found, returns maximum int
// allowed. This ensures in "NaturalComp" that files not indexed by an
// integer at the end of the csv
float GetLastNumOfString(std::string str)
{
    float last_num = 2147483647;
    size_t begin = 0, end = 0;

    std::string numbers = ".0123456789";

    begin = str.find_first_of(numbers);

    while (begin != std::string::npos)
    {
        // avoid cases where is just single "." in substring
        if (str.at(begin) == '.' &&
            numbers.find(str.at(begin + 1)) == std::string::npos)
        {
            str = str.substr(begin + 1, str.length() - (begin + 1));
            begin = str.find_first_of(numbers);
            continue;
        }

        end = str.find_first_not_of(numbers, begin);
        std::string num = str.substr(begin, end - begin);
        last_num = std::atof(num.c_str());
        str = str.substr(end, str.length() - end);
        begin = str.find_first_of(numbers);
    }

    return last_num;
}

bool NaturalComp(std::string s1, std::string s2)
{
    return (GetLastNumOfString(s1) < GetLastNumOfString(s2));
}

// Converts full amplitude named "reaction::reflectivity::JPmL" to
// "eJPmL" format, where reaction name is dropped
std::string ConvertFullAmplitudeName(std::string full_amp)
{
    std::string eJPmL;
    std::string delim = "::";

    // extract reflectivity and eJPmL amplitude
    // assumes amplitude always in "reaction::refl::amp" format
    size_t positionRefl = full_amp.find(delim) + delim.length();
    size_t positionAmp = full_amp.find(delim, positionRefl);
    std::string refl = full_amp.substr(positionRefl, positionAmp - positionRefl);
    eJPmL = full_amp.substr(positionAmp + delim.length(),
                            full_amp.length() - positionAmp + delim.length());

    // capitalize L value to avoid ambiguity between L=p and m=p values
    size_t last = std::strlen(eJPmL.c_str()) - 1;
    eJPmL[last] = toupper(eJPmL[last]);

    // add reflectivity to string
    if (refl == "ImagNegSign" || refl == "RealPosSign")
    {
        eJPmL = "m" + eJPmL;
    }
    else
    {
        eJPmL = "p" + eJPmL;
    }

    return eJPmL;
}

// Write values from maps to csv file.
// Args:
//      csv_file: file to write results to
//      results: current FitResults file loaded
//      mapCohSum: map of headers to corresponding full amplitude names
//      is_header: writes keys of mapCohSum (headers) to file when true
//      is_acceptance_corrected: when true, scales up any intensity values to correct
//          for the detector acceptance. False corresponds to the "detected" events
// NOTE: NaNs are currently treated as 0's to avoid issues when reading into
// ROOT Tree. Real and Imaginary parameters are also assumed to be constrained across
// both reflectivities, so only one of the reflectivities is accessed to get the value
void WriteToCsv(ofstream &csv_file, FitResults &results,
                std::map<std::string, std::vector<std::string>> &mapCohSum,
                bool is_header, bool is_acceptance_corrected)
{
    double val, err;
    for (auto iter = mapCohSum.begin(); iter != mapCohSum.end(); iter++)
    {
        if (is_header) // write the keys to the csv if writing the header line
        {
            csv_file << "," << iter->first << "," << iter->first + "_err";
            continue;
        }
        // handle cases for writing the real or imaginary parameters of amps
        else if ((iter->first).find("_re") != std::string::npos)
        {
            val = results.parValue(results.realProdParName((iter->second)[0]));
            err = results.parError(results.realProdParName((iter->second)[0]));
        }
        else if ((iter->first).find("_im") != std::string::npos)
        {
            val = results.parValue(results.imagProdParName((iter->second)[0]));
            err = results.parError(results.imagProdParName((iter->second)[0]));
        }
        else
        {
            val = results.intensity(iter->second, is_acceptance_corrected).first;
            err = results.intensity(iter->second, is_acceptance_corrected).second;
        }
        if (TMath::IsNaN(val))
        {
            val = 0;
        }
        if (TMath::IsNaN(err))
        {
            err = 0;
        }

        csv_file << "," << val << "," << err;
    }
}

// overloaded function to handle mapping for phase differences
void WriteToCsv(ofstream &csv_file, FitResults &results,
                std::map<std::string, pair<std::string, std::string>> &mapCohSum,
                bool is_header, bool is_acceptance_corrected)
{
    for (auto iter = mapCohSum.begin(); iter != mapCohSum.end(); iter++)
    {
        if (is_header) // write the keys to the csv if writing the header line
        {
            csv_file << "," << iter->first << "," << iter->first + "_err";
        }
        else
        {
            pair<double, double> phaseDiff = results.phaseDiff((iter->second).first,
                                                               (iter->second).second);

            if (TMath::IsNaN(phaseDiff.first))
            {
                phaseDiff.first = 0;
            }
            if (TMath::IsNaN(phaseDiff.second))
            {
                phaseDiff.second = 0;
            }

            csv_file << "," << phaseDiff.first << "," << phaseDiff.second;
        }
    }
}
