/* Extract the bin information from a list of pre-cut ROOT data files used in fits

The csv file will have columns for:
    - The low and high bin edges for the t, beam energy, and mass histograms
    - The center, average, and RMS values for the t, beam energy, and mass histograms
    - The total number of events and the error on the total number of events
 */

#include <cmath>   // for power function
#include <filesystem> // for getting paths
#include <fstream> // for writing csv
#include <iostream>
#include <sstream> // for std::istringstream
#include <string>
#include <vector>

#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

#include "neutralb1/file_utils.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigurationInfo.h"

// forward declarations
std::pair<double, double> get_hist_edges(TH1D *h, int round_to_decimals);
void fill_map_with_variable(
    std::map<std::string, 
    double> &map, const std::string &key, 
    TH1D *hist, 
    int round_to_decimals);
std::vector<std::string> get_data_files(std::string fit_file);
std::vector<std::string> get_background_files(std::string fit_file);
std::string get_tree_name(FitResults& results);
std::string get_weight_branch(std::string file, std::string tree_name, std::string fit_file = "");
std::string get_beamE_branch(std::string root_file, std::string tree_name);
std::string get_t_branch(std::string root_file, std::string tree_name);

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] 
        << " <file_list.txt> <output.csv> <mass_branch>\n";
        
        return 1;
    }

    std::string file_path = argv[1];
    std::string csv_name = argv[2];
    std::string mass_branch = argv[3];
    
    // file path is to a text file with a list of fit files, each on a newline
    std::vector<std::string> file_vector = read_file_list(file_path);
    if (file_vector.empty())
    {
        std::cerr << "Error: No files found in " << file_path << "\n";
        return 1;
    }

    // Define header row
    std::vector<std::string> headers = {
        "file",
        "t_low",
        "t_high",
        "t_center",
        "t_avg",
        "t_rms",
        "e_low",
        "e_high",
        "e_center",
        "e_avg",
        "e_rms",
        "m_low",
        "m_high",
        "m_center",
        "m_avg",
        "m_rms",
        "events",
        "events_err",
    };
    // setup map to store values for each file
    std::vector<std::map<std::string, double>> values;

    for (const std::string &file : file_vector)
    {
        std::map<std::string, double> value_map; // initialize map for this file

        FitResults results(file);
        if (!results.valid())
        {
            throw std::runtime_error("Invalid FitResults for file: " + file);
        }

        // read in data and background files from the fit results
        std::vector<std::string> data_files = get_data_files(file);
        if (data_files.empty())
        {
            throw std::runtime_error("No data files found in fit results for file: " + file);
        }
        std::vector<std::string> background_files = get_background_files(file);
        if (background_files.empty())
        {
            std::cerr << "No background files found in fit results for file: " 
                      << file << ". Asssuming weights are stored directly in the data"
                      << " trees.\n";            
        }

        // get tree name, weight branch name, and beamE branch name
        std::string tree_name = get_tree_name(results);
        std::string weight_branch;
        if (background_files.empty())
            weight_branch = get_weight_branch(data_files[0], tree_name);
        else
            weight_branch = get_weight_branch(background_files[0], tree_name, file);

        std::string beamE_branch = get_beamE_branch(data_files[0], tree_name);
        std::string t_branch = get_t_branch(data_files[0], tree_name);

        // now add all the data (signal) file histograms together
        TH1D *h_t_total = nullptr;
        TH1D *h_e_total = nullptr;
        TH1D *h_m_total = nullptr;
        TCanvas *c = new TCanvas("c", "c", 100, 100);    
        for (const std::string &data_file : data_files)
        {
            // Open the ROOT file
            std::unique_ptr<TFile> f(TFile::Open(data_file.c_str()));
            if (!f || f->IsZombie())
            {
                throw std::runtime_error("Error opening data file: " + data_file);
            }

            TH1D *h_t = nullptr;
            TH1D *h_e = nullptr;
            TH1D *h_m = nullptr;

            // Check that tree can be opened
            TTree *tree = f->Get<TTree>(tree_name.c_str());
            if (!tree)
            {
                std::cout << "'" << tree_name << "' tree could not be opened in file: " 
                          << data_file << "\n";
                exit(1);
            }

            // fill the file histograms, and include weight if applicable
            if (background_files.empty())
            {
                tree->Draw((t_branch + ">>h_t").c_str(), weight_branch.c_str(), "goff");
                tree->Draw((beamE_branch + ">>h_e").c_str(), weight_branch.c_str(), "goff");
                tree->Draw((mass_branch + ">>h_m").c_str(), weight_branch.c_str(), "goff");            
            }
            else
            {
                tree->Draw((t_branch + ">>h_t").c_str(), "", "goff");
                tree->Draw((beamE_branch + ">>h_e").c_str(), "", "goff");
                tree->Draw((mass_branch + ">>h_m").c_str(), "", "goff");            
            }

            if (!h_t || !h_m || !h_e)
            {
                throw std::runtime_error("Error: Could not retrieve histograms from data file: " + data_file);
            }

            // Add current histograms to the total histograms
            if (!h_t_total)
                h_t_total = (TH1D *)h_t->Clone("h_t_total");
            else
                h_t_total->Add(h_t);
            if (!h_e_total)
                h_e_total = (TH1D *)h_e->Clone("h_e_total");
            else
                h_e_total->Add(h_e);
            if (!h_m_total)
                h_m_total = (TH1D *)h_m->Clone("h_m_total");
            else
                h_m_total->Add(h_m);

            // Clean up temporary histograms
            delete h_t;
            delete h_e;
            delete h_m;           
        }

        // now do the same for the background files, if applicable
        if (!background_files.empty())
        {
            for (const std::string &bg_file : background_files)
            {
                // Open the ROOT file
                std::unique_ptr<TFile> f(TFile::Open(bg_file.c_str()));
                if (!f || f->IsZombie())
                {
                    throw std::runtime_error("Error opening background file: " + bg_file);
                }
                TH1D *h_t = nullptr;
                TH1D *h_e = nullptr;
                TH1D *h_m = nullptr;
                TTree *tree = f->Get<TTree>(tree_name.c_str());
                if (!tree)
                {
                    std::cout << "'" << tree_name << "' tree could not be opened in file: " 
                              << bg_file << "\n";
                    exit(1);
                }

                tree->Draw("t>>h_t", weight_branch.c_str(), "goff");
                tree->Draw((beamE_branch + ">>h_e").c_str(), weight_branch.c_str(), "goff");
                tree->Draw((mass_branch + ">>h_m").c_str(), weight_branch.c_str(), "goff");

                if (!h_t || !h_m || !h_e)
                {
                    throw std::runtime_error("Error: Could not retrieve histograms from background file: " + bg_file);
                }

                // Subtract current histograms from the total histograms
                h_t_total->Add(h_t, -1.0);
                h_e_total->Add(h_e, -1.0);
                h_m_total->Add(h_m, -1.0);

                // Clean up temporary histograms
                delete h_t;
                delete h_e;
                delete h_m;           
            }
        }
        delete c;


        // use one of the hists to get the total number of events and its error
        if (!h_m_total || !h_t_total || !h_e_total)
        {
            throw std::runtime_error("Error: Total histograms not properly filled for file: " 
                + file + "\n");
        }
        double bin_contents_error;
        value_map["events"] = h_m_total->IntegralAndError(0, h_m_total->GetNbinsX() + 1, bin_contents_error);
        value_map["events_err"] = bin_contents_error;

        // Fill the map. Round -t and E_beam to 2nd decimal, and the mass values to
        // the third decimal (1 MeV)
        fill_map_with_variable(value_map, "t", h_t_total, 2);
        fill_map_with_variable(value_map, "e", h_e_total, 2);
        fill_map_with_variable(value_map, "m", h_m_total, 3);
        
        values.push_back(value_map);
    } // end loop over files

    // open csv file for writing
    std::ofstream csv_file;
    csv_file.open(csv_name);

    // Write the header line
    for (size_t i = 0; i < headers.size(); ++i)
    {
        csv_file << headers[i];
        if (i < headers.size() - 1)
        {
            csv_file << ",";
        }
    }
    csv_file << "\n";

    // Write values
    for (size_t i = 0; i < values.size(); ++i)
    {
        csv_file << file_vector[i] << ",";
        for (size_t j = 1; j < headers.size(); ++j)
        {
            csv_file << values[i].at(headers[j]);
            if (j < headers.size() - 1)
            {
                csv_file << ",";
            }
        }
        csv_file << "\n";
    }

    csv_file.close();
    return 0;
}

/**
 * @brief Get the min / max edges of a histogram
 * 
 * @param[in] h input histogram
 * @param[in] round_to_decimals number of decimal places to round to
 * @return std::pair<double, double> min and max edges of the histogram
 */
std::pair<double, double> get_hist_edges(TH1D *h, int round_to_decimals)
{
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();

    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
        double bin_content = h->GetBinContent(i);
        if (bin_content > 0 && h->GetXaxis()->GetBinLowEdge(i) < min)
        {
            min = h->GetXaxis()->GetBinLowEdge(i);
        }
        if (bin_content > 0 && h->GetXaxis()->GetBinUpEdge(i) > max)
        {
            max = h->GetXaxis()->GetBinUpEdge(i);
        }
    }
    // round values to requested decimal place
    min = std::round((min * std::pow(10, round_to_decimals))) / std::pow(10, round_to_decimals);
    max = std::round((max * std::pow(10, round_to_decimals))) / std::pow(10, round_to_decimals);
    return std::make_pair(min, max);
}

/**
 * @brief Fill a map with statistical information from a histogram
 * 
 * @param[out] map map to fill with statistical information
 * @param key key prefix for the map entries
 * @param hist input histogram
 * @param round_to_decimals number of decimal places to round to
 * 
 * Specifically, this function fills the map with the following entries:
 * - key_low: the low edge of the histogram
 * - key_high: the high edge of the histogram
 * - key_center: the center of the histogram (average of low and high edges)
 * - key_avg: the average x-axis value of the histogram
 * - key_rms: the x-axis RMS value of the histogram
 */
void fill_map_with_variable(
    std::map<std::string, 
    double> &map, const std::string &key, 
    TH1D *hist, 
    int round_to_decimals)
{
    if (!hist)
    {
        std::cerr << "Error: Histogram is null for key: " << key << "\n";
        return;
    }
    double low, high;
    std::tie(low, high) = get_hist_edges(hist, round_to_decimals);
    map[key + "_low"] = low;
    map[key + "_high"] = high;
    map[key + "_center"] = (low + high) / 2.0;
    map[key + "_avg"] = hist->GetMean();
    map[key + "_rms"] = hist->GetRMS();
}


/**
 * @brief Get the data files from a fit results file
 * 
 * @param fit_file path to the fit results file
 * @return std::vector<std::string> vector of data file absolute paths used in the fit
 */
std::vector<std::string> get_data_files(std::string fit_file)
{
    FitResults results(fit_file);

    // get the directory of the fit file
    std::filesystem::path dir = std::filesystem::path(fit_file).parent_path();

    std::vector<std::string> data_files;

    const std::vector<std::string> reaction_list = results.reactionList();
    for (std::string reaction : reaction_list)
    {
        const std::pair< std::string, std::vector<std::string> > data_pair = 
            results.configInfo()->reaction(reaction)->data();

        std::string data_reader_name = data_pair.first;
        std::vector<std::string> data_file_vector = data_pair.second;

        // NOTE: this assumes only one data file per reaction, and that it is the first
        // argument passed to the data reader. If this is not the case, one can perform
        // data_reader specific parsing here based off the data_reader_name variable.
        std::filesystem::path data_path = data_file_vector[0]; 

        // assumes that relative paths are relative to the fit file directory
        if (!data_path.is_absolute())
            data_path = dir / data_path;
        
        // canonical will remove any ../ or ./ from the path, while also checking
        // for existence
        data_files.push_back(std::filesystem::canonical(data_path).string());
    }
    
    return data_files;    
}

/**
 * @brief Get the background files from a fit results file
 * 
 * @param fit_file path to the fit results file
 * @return std::vector<std::string> vector of background file absolute paths used in the fit
 */
std::vector<std::string> get_background_files(std::string fit_file)
{
    FitResults results(fit_file);

    // get the directory of the fit file
    std::filesystem::path dir = std::filesystem::path(fit_file).parent_path();

    std::vector<std::string> background_files;

    const std::vector<std::string> reaction_list = results.reactionList();
    for (std::string reaction : reaction_list)
    {
        const std::pair< std::string, std::vector<std::string> > bkgd_pair = 
            results.configInfo()->reaction(reaction)->bkgnd();

        std::string data_reader_name = bkgd_pair.first;
        std::vector<std::string> bkgd_file_vector = bkgd_pair.second;

        // NOTE: this assumes only one background file per reaction, and that it is the first
        // argument passed to the data reader. If this is not the case, one can perform
        // data_reader specific parsing here based off the data_reader_name variable.
        std::filesystem::path bkgd_path = bkgd_file_vector[0]; 

        // assumes that relative paths are relative to the fit file directory
        if (!bkgd_path.is_absolute())
            bkgd_path = dir / bkgd_path;
        
        // canonical will remove any ../ or ./ from the path, while also checking
        // for existence
        background_files.push_back(std::filesystem::canonical(bkgd_path).string());
    }
    
    return background_files;    
}


/**
 * @brief Get the tree name from a FitResults object
 * 
 * Unique tree names are specific to the data reader used, so this function will extract
 * the tree name for hard-coded known data readers. If the data reader is unknown, it 
 * will default to "kin".
 * 
 * @param results FitResults object
 * @return std::string tree name
 */
std::string get_tree_name(FitResults& results)
{
    const std::vector<std::string> reaction_list = results.reactionList();

    // assume all reactions use the same tree name, so just get the first one
    std::string reaction = reaction_list[0];
    const std::pair< std::string, std::vector<std::string> > data_pair = 
        results.configInfo()->reaction(reaction)->data();

    std::string data_reader_name = data_pair.first;
    std::vector<std::string> data_file_vector = data_pair.second;

    // unfortunately, different data readers have different argument orderings and 
    // the DataReader only stores the arguments passed to it, so we have to build 
    // in knowledge of the argument orderings here.
    if (data_reader_name == "ROOTDataReader" && data_file_vector.size() == 2)
    {
        return data_file_vector[1];        
    }
    else if (data_reader_name == "ROOTDataReaderBootstrap" && data_file_vector.size() == 3)
    {
        return data_file_vector[2];        
    }
    else if (data_reader_name == "FSRootDataReader")
    {
        return data_file_vector[1];
    }
    else if (data_reader_name == "FSRootDataReaderBootstrap")
    {
        return data_file_vector[1];
    }
    else
    {
        std::cerr << "Warning: Unknown data reader or insufficient arguments to"
                  << " determine tree name. Defaulting to 'kin'\n";
        return "kin"; // default tree name
    }
}

/**
 * @brief Get the name of the weight branch
 * 
 * The optional fit file is for cases where the weight branch is stored in a friend tree,
 * such as with the FSRootDataReader. If no fit file is provided, it will simply check
 * for "weight" or "Weight" branches in the main tree. If those branches do not exist,
 * it will return an empty string (weight of 1.0).
 * 
 * @param root_file base root file (signal or background)
 * @param tree_name tree name
 * @param fit_file path to the fit results file (optional)
 * @return std::string weight branch name
 */
std::string get_weight_branch(std::string root_file, std::string tree_name, std::string fit_file)
{
    // FSRoot data readers can optionally pass a friend weight branch name as an 
    // argument, so check for that first
    if (!fit_file.empty())
    {
        FitResults results(fit_file);

        const std::vector<std::string> reaction_list = results.reactionList();

        // assume all reactions use the same weight branch, so just get the first one
        std::string reaction = reaction_list[0];
        const std::pair< std::string, std::vector<std::string> > data_pair = 
            results.configInfo()->reaction(reaction)->data();

        std::string data_reader_name = data_pair.first;
        std::vector<std::string> data_file_vector = data_pair.second;

        // FSDataReader only case where friend tree weights are used
        if (data_reader_name.find("FSRootDataReader") != std::string::npos 
            && data_file_vector.size() >= 6)
        {
            std::string friend_file_name = data_file_vector[3];
            std::string friend_tree_name = data_file_vector[4];
            std::string weight_branch_name = data_file_vector[5];

            // need to add friend tree to get the weight branch
            TFile *f = TFile::Open(root_file.c_str());
            if (!f || f->IsZombie())
            {
                throw std::runtime_error("Error opening ROOT file: " + root_file);
            }
            TTree *tree = f->Get<TTree>(tree_name.c_str());
            if (!tree)
            {
                throw std::runtime_error("Error: Tree '" + tree_name + "' not found in file: " + root_file);
            }
            tree->AddFriend(friend_tree_name.c_str(), friend_file_name.c_str());

            return friend_tree_name + "." + weight_branch_name;
        }
        else // not FSDataReader, so no friend tree weights
        {
            // proceed to check for weight branches in the main tree
        }
    }

    // otherwise, simply check if "weight" or "Weight" branch exists in the tree
    TFile *f = TFile::Open(root_file.c_str());
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Error opening ROOT file: " + root_file);
    }
    TTree *tree = f->Get<TTree>(tree_name.c_str());
    if (!tree)
    {
        throw std::runtime_error("Error: Tree '" + tree_name + "' not found in file: " + root_file);
    }

    // check if weight branch exists before returning it. Most analyses use "weight" or 
    // "Weight" branches
    if (tree->GetBranch("weight"))
    {
        f->Close();
        return "weight";
    }
    else if (tree->GetBranch("Weight"))
    {
        f->Close();
        return "Weight";
    }
    else
    {
        f->Close();
        std::cerr << "Error: Weight branch not found in tree '" << tree_name 
                  << "' in file: " << root_file << ". Setting weights to 1.0\n";
        return ""; // return empty string to indicate no weight branch found
    }
}


/**
 * @brief Get the t branch name
 * 
 * Currently supports "t" and "-t" branch names, as these are most common.
 * 
 * @param root_file data root file path
 * @param tree_name tree name
 * @return std::string t branch name
 */
std::string get_t_branch(std::string root_file, std::string tree_name)
{
    // Open the ROOT file
    TFile *f = TFile::Open(root_file.c_str());
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Error opening ROOT file: " + root_file);
    }
    TTree *tree = f->Get<TTree>(tree_name.c_str());
    if (!tree)
    {
        throw std::runtime_error("Error: Tree '" + tree_name + "' not found in file: " + root_file);
    }

    // check if t branch exists before returning it. Most analyses use "t" or 
    // "-t" branches
    if (tree->GetBranch("t"))
    {
        f->Close();
        return "t";
    }
    else if (tree->GetBranch("-t"))
    {
        f->Close();
        return "-t";
    }
    else
    {
        f->Close();
        throw std::runtime_error(
            "Error: t branch not found in tree '" + tree_name 
            + "' in file: " + root_file + ".");
    }
}

/**
 * @brief Get the beam energy branch name
 * 
 * Currently supports "E_Beam" and "EnPB" branch names, as these are most common.
 *
 * @param root_file data root file path
 * @param tree_name tree name
 * @return std::string beam energy branch name
 */
std::string get_beamE_branch(std::string root_file, std::string tree_name)
{
    TFile *f = TFile::Open(root_file.c_str());
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Error opening ROOT file: " + root_file);
    }
    TTree *tree = f->Get<TTree>(tree_name.c_str());
    if (!tree)
    {
        throw std::runtime_error("Error: Tree '" + tree_name + "' not found in file: " + root_file);
    }

    // check if beamE branch exists before returning it
    // 
    
    if (tree->GetBranch("E_Beam"))
    {
        f->Close();
        return "E_Beam";
    }
    else if (tree->GetBranch("EnPB"))
    {
        f->Close();
        return "EnPB";
    }
    else
    {
        f->Close();
        throw std::runtime_error(
            "Error: beamE branch not found in tree '" + tree_name 
            + "' in file: " + root_file + ".");
    }
}