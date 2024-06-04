/* Read .fit files and output csv containing coherent sums and parameters.

   Two main functionalities are to combine many .fit files in one
   directory, or combine .fit files with a particular name in many
   subdirectories. Execute script with "-h" or "--help" to see options

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

// Forward declarations of functions. Comments in main should make clear how
// the program executes.
void GetListOfFitFiles(std::vector<std::string> &file_list,
                       std::string parent_dir, std::string search_for_file);
float GetLastNumOfString(std::string str);
bool NaturalComp(std::string s1, string s2);
std::string ConvertFullAmplitudeName(std::string full_amp);
void WriteValueToCsv(ofstream &csv_file, FitResults &results,
                     std::map<std::string,
                              std::vector<std::string>> &mapCohSum);

void fitsToCsv(string args = "")
{
    std::string search_for_file = "";
    std::string parent_dir = "./";
    std::string output_csv_name = "combinedFitPars.csv";

    // PARSE COMMAND LINE ARGS
    size_t flag_pos = 0;
    size_t arg_pos = 0;
    while (args.length() != 0)
    {

        // get positions of flag (-flag) and argument that follows the flag
        flag_pos = args.find(" ");
        arg_pos = args.find(" ", flag_pos + 1);

        if (flag_pos == string::npos && (args != "-h" && args != "--help"))
        {
            cout << "Invalid input. Use \"--help\" to see options"
                 << "\n";
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
        if (flag == "-h" || flag == "--help")
        {
            cout << "usage: analyzeFitResults [-h] [-f FILE_NAME] "
                 << "[-o OUTPUT_FILE_NAME] [-d PARENT_DIRECTORY]"
                 << "\n\n"
                 << "Aggregate results from many .fit files into a CSV.\n"
                 << "Running without argument combines .fit files in current dir."
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
                 << "\t\tspecify parent dir (default: ./)\n";
            exit(0);
        }
        args.erase(0, arg_pos);
        args.erase(0, 1);
    }

    // GET FILES AND SORT THEM IF INDEXED
    vector<string> file_list;
    GetListOfFitFiles(file_list, parent_dir, search_for_file);

    if (file_list.size() == 0)
    {
        cout << "No files found!"
             << "\n";
        exit(1);
    }

    // sort vector using last integer in element
    std::sort(file_list.begin(), file_list.end(), NaturalComp);

    // index at 0, unless files begin indexing at 1
    int index = 0;

    ofstream csv_file;
    csv_file.open(output_csv_name);

    // INITIALIZE CONTAINERS
    vector<string> amp_list;
    vector<string> par_list;

    // setup map of eJPmL based header name, to a vector that stores the full
    // amplitude name. Excluded variables mean they are summed over
    // ex: sum of JP=1+,l=0 waves is:
    // < 1ps , < xx::ImagNegSign::1pps, xx::RealPosSign:1pms, ... > >
    std::map<string, vector<string>> amp_sum_eJPmL;
    std::map<string, vector<string>> amp_sum_JPmL;
    std::map<string, vector<string>> amp_sum_JPm;
    std::map<string, vector<string>> amp_sum_JPL;
    std::map<string, vector<string>> amp_sum_eJPm;
    std::map<string, vector<string>> amp_sum_eJPL;
    std::map<string, vector<string>> amp_sum_eJP;
    std::map<string, vector<string>> amp_sum_JP;
    std::map<string, pair<string, string>> ampPhaseDiffs; // eJPmL_eJPmL
    // TODO: add 2 coherent sums for total reflectivity 

    for (string file : file_list)
    {
        // load .fit file
        cout << "Analyzing File: " << file << "\n";

        FitResults results(file);
        if (!results.valid())
        {
            cout << "Invalid fit results in file:  " << file << "\n";
            continue;
        }

        cout << "Fit results loaded"
             << "\n";

        // WRITE HEADER LINE
        if (amp_list.size() == 0)
        {
            vector<string> reactions = results.reactionList();
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

            // write parameter headers
            for (string par : par_list)
            {
                // only want Parameters, not reaction-based ones
                if (par.find("::") != string::npos)
                {
                    continue;
                }
                csv_file << "," << par << "," << par + "_err";

                // write covariance between any D/S params if they exist
                if (par.find("dsratio") != string::npos)
                {
                    // get refl of ratio param if it exists
                    string ratio_refl = "";
                    if (par.find("_") != string::npos)
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

            // add full amplitude name (value) to corresponding eJPmL key
            for (unsigned int i = 0; i < amp_list.size(); i++)
            {
                string full_amp = amp_list[i];

                // add bkgd to "eJPmL" vector once
                if (full_amp.find("Bkgd") != string::npos)
                {
                    amp_sum_eJPmL["Bkgd"].push_back(full_amp);
                    continue;
                }

                string eJPmL = ConvertFullAmplitudeName(full_amp);
                string e = eJPmL.substr(0, 1), JP = eJPmL.substr(1, 2),
                       m = eJPmL.substr(3, 1), L = eJPmL.substr(4, 1);

                amp_sum_eJPmL[eJPmL].push_back(full_amp);
                amp_sum_JPmL[JP + m + L].push_back(full_amp);
                amp_sum_JPm[JP + m].push_back(full_amp);
                amp_sum_JPL[JP + L].push_back(full_amp);
                amp_sum_eJPm[e + JP + m].push_back(full_amp);
                amp_sum_eJPL[e + JP + L].push_back(full_amp);
                amp_sum_eJP[e + JP].push_back(full_amp);
                amp_sum_JP[JP].push_back(full_amp);

                // second loop to get phase difference headers
                for (unsigned int j = i + 1; j < amp_list.size(); j++)
                {
                    string pd_full_amp = amp_list[j];

                    if (pd_full_amp.find("isotropic") != string::npos)
                    {
                        continue;
                    }

                    string pd_eJPmL = ConvertFullAmplitudeName(pd_full_amp);

                    // keep amplitudes from same coherent sum
                    if (pd_eJPmL.substr(0, 1) != eJPmL.substr(0, 1))
                    {
                        continue;
                    }

                    ampPhaseDiffs[eJPmL + "_" + pd_eJPmL] = std::make_pair(full_amp, pd_full_amp);
                }
            } // end mapping for-loop

            // write amplitude headers !! KEEP ORDER WHEN WRITING VALUES LATER!!
            for (auto iter = amp_sum_eJPmL.begin(); iter != amp_sum_eJPmL.end(); iter++)
            {
                csv_file << "," << iter->first << "," << iter->first + "_err";
            }
            for (auto iter = amp_sum_JPmL.begin(); iter != amp_sum_JPmL.end(); iter++)
            {
                csv_file << "," << iter->first << "," << iter->first + "_err";
            }
            for (auto iter = amp_sum_JPm.begin(); iter != amp_sum_JPm.end(); iter++)
            {
                csv_file << "," << iter->first << "," << iter->first + "_err";
            }
            for (auto iter = amp_sum_JPL.begin(); iter != amp_sum_JPL.end(); iter++)
            {
                csv_file << "," << iter->first << "," << iter->first + "_err";
            }
            for (auto iter = amp_sum_eJPm.begin(); iter != amp_sum_eJPm.end(); iter++)
            {
                csv_file << "," << iter->first << "," << iter->first + "_err";
            }
            for (auto iter = amp_sum_eJPL.begin(); iter != amp_sum_eJPL.end(); iter++)
            {
                csv_file << "," << iter->first << "," << iter->first + "_err";
            }
            for (auto iter = amp_sum_eJP.begin(); iter != amp_sum_eJP.end(); iter++)
            {
                csv_file << "," << iter->first << "," << iter->first + "_err";
            }
            for (auto iter = amp_sum_JP.begin(); iter != amp_sum_JP.end(); iter++)
            {
                csv_file << "," << iter->first << "," << iter->first + "_err";
            }
            for (auto iter = ampPhaseDiffs.begin(); iter != ampPhaseDiffs.end(); iter++)
            {
                csv_file << "," << iter->first << "," << iter->first + "_err";
            }

            csv_file << "\n";
        } // end header line

        // WRITE VALUES TO CSV
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

        for (string par : par_list)
        {
            // only want Parameters, not reaction-based ones
            if (par.find("::") != string::npos)
            {
                continue;
            }
            csv_file << "," << results.parValue(par)
                     << "," << results.parError(par);

            // write D/S ratio phase covariance if it exists
            if (par.find("dsratio") != string::npos)
            {
                // get refl of ratio param if exists
                string ratio_refl = "";
                if (par.find("_") != string::npos)
                {
                    ratio_refl = par.back();
                }

                if (ratio_refl != "")
                {
                    string r = "dsratio_" + ratio_refl;
                    string ph = "dphase_" + ratio_refl;
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
        WriteValueToCsv(csv_file, results, amp_sum_eJPmL);
        WriteValueToCsv(csv_file, results, amp_sum_JPmL);
        WriteValueToCsv(csv_file, results, amp_sum_JPm);
        WriteValueToCsv(csv_file, results, amp_sum_JPL);
        WriteValueToCsv(csv_file, results, amp_sum_eJPm);
        WriteValueToCsv(csv_file, results, amp_sum_eJPL);
        WriteValueToCsv(csv_file, results, amp_sum_eJP);
        WriteValueToCsv(csv_file, results, amp_sum_JP);

        // results.rotateResults(); // TEMP
        for (auto iter = ampPhaseDiffs.begin(); iter != ampPhaseDiffs.end(); iter++)
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

        string file = ent->d_name;
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
float GetLastNumOfString(string str)
{
    float last_num = 2147483647;
    size_t begin = 0, end = 0;

    string numbers = ".0123456789";

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
        string num = str.substr(begin, end - begin);
        last_num = std::atof(num.c_str());
        str = str.substr(end, str.length() - end);
        begin = str.find_first_of(numbers);
    }

    return last_num;
}

bool NaturalComp(string s1, string s2)
{
    return (GetLastNumOfString(s1) < GetLastNumOfString(s2));
}

// Converts full amplitude named "reaction::reflectivity::JPmL" to
// "eJPmL" format, where reaction name is dropped
string ConvertFullAmplitudeName(string full_amp)
{
    string eJPmL;
    string delim = "::";

    // extract reflectivity and eJPmL amplitude
    // assumes amplitude always in "reaction::refl::amp" format
    size_t positionRefl = full_amp.find(delim) + delim.length();
    size_t positionAmp = full_amp.find(delim, positionRefl);
    string refl = full_amp.substr(positionRefl, positionAmp - positionRefl);
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
// NOTE: NaNs are currently treated as 0's to avoid issues when reading into
// ROOT Tree. Amplitudes give their "generated" intensity when "false" is passed to 
// results.intensity()
void WriteValueToCsv(ofstream &csv_file, FitResults &results,
                     std::map<string, vector<string>> &mapCohSum)
{
    for (auto iter = mapCohSum.begin(); iter != mapCohSum.end(); iter++)
    {
        double val = results.intensity(iter->second, false).first;
        double err = results.intensity(iter->second, false).second;

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
