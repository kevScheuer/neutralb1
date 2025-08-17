/* Extract the bin information from a list of pre-cut ROOT data files used in fits

The csv file will have columns for:
    - The low and high bin edges for the t, E_beam, and mass histograms
    - The center, average, and RMS values for the t, E_beam, and mass histograms
    - The total number of events and the error on the total number of events

NOTE:
This script assumes that the original Flat Tree data files have been cut to their
respective bin, and that they contain the t, E_beam, and Weight branches. The
Weight branch is used for sideband subtraction, and so if a separate "background" file
is used, then it will need to be implemented here.
 */

#include <cmath>   // for power function
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

// forward declarations
std::pair<double, double> get_hist_edges(TH1D *h, int round_to_decimals);
void fill_map_with_variable(std::map<std::string, double> &map, const std::string &key, TH1D *hist);

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        std::cerr << "Usage: " << argv[0] << " <file_list.txt> <output.csv> <mass_branch>\n";
        return 1;
    }

    std::string file_path = argv[1];
    std::string csv_name = argv[2];
    std::string mass_branch = argv[3];

    // file path is a text file with a list of ROOT files, each on a newline
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
    // setup map to store values
    std::vector<std::map<std::string, double>> values;

    for (const std::string &file : file_vector)
    {
        // determine if this is a data file or a gen_vec_ps file
        // this will determine how the trees are processed
        bool is_gen = false;
        if (file.find("gen_vec_ps") != std::string::npos)
            is_gen = true;
        
        std::map<std::string, double> value_map; // initialize map for this file

        // Open the ROOT file
        std::unique_ptr<TFile> f(TFile::Open(file.c_str()));
        if (!f || f->IsZombie())
        {
            std::cerr << "Error opening file: " << file << "\n";
            continue; // skip to the next file
        }

        TCanvas *c1 = new TCanvas("c1", "c1", 1800, 1800);
        // assign the histograms based on whether it is a gen file or not
        TH1D *h_t = nullptr;
        TH1D *h_e = nullptr;
        TH1D *h_m = nullptr;
        if (is_gen)
        {
            // NOTE: for now gen_vec_ps_diagnostic does not contain the E_beam hist.
            //      Since E_beam is not crucial to the analysis and reading the gen
            //      file is a simple operation rarely used, we can skip it
            h_t = (TH1D *)f->Get("t");
            h_m = (TH1D *)f->Get("M");            
        }
        else
        {
            TTree *tree = f->Get<TTree>("kin");
            if (!tree)
            {
                std::cout << "'kin' tree could not be opened in file: " << file << "\n";
                exit(1);
            }
            tree->Draw("t>>h_t", "Weight");
            h_t = (TH1D *)gPad->GetPrimitive("h_t");
            tree->Draw("E_Beam>>h_e", "Weight");
            h_e = (TH1D *)gPad->GetPrimitive("h_e");
            tree->Draw((mass_branch + ">>h_m").c_str(), "Weight");
            h_m = (TH1D *)gPad->GetPrimitive("h_m");
        }

        if (!h_t || !h_m)
        {
            std::cerr << "Error: Could not retrieve histograms from file: " << file << "\n";
            continue; // skip to the next file
        }

        // Fill the map. Round -t and E_beam to 2nd decimal, and the mass values to
        // the third decimal (1 MeV)
        double bin_contents_error;
        value_map["events"] = h_m->IntegralAndError(0, h_m->GetNbinsX() + 1, bin_contents_error);
        value_map["events_err"] = bin_contents_error;

        fill_map_with_variable(value_map, "t", h_t);
        if (is_gen)
            // gen files do not have E_beam histograms
            value_map["e_low"] = value_map["e_high"] = value_map["e_center"] = value_map["e_avg"] = value_map["e_rms"] = 0.0;
        else
            fill_map_with_variable(value_map, "e", h_e);
        fill_map_with_variable(value_map, "m", h_m);
        
        values.push_back(value_map);

        // Clean up
        delete h_t;
        delete h_e;
        delete h_m;
        delete c1;
    }

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

// Get the min / max non-zero bin edges of a histogram
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

void fill_map_with_variable(std::map<std::string, double> &map, const std::string &key, TH1D *hist)
{
    if (!hist)
    {
        std::cerr << "Error: Histogram is null for key: " << key << "\n";
        return;
    }
    double low, high;
    // round mass to 3 decimals for single MeV precision
    if (key == "m")
        std::tie(low, high) = get_hist_edges(hist, 3);
    else
        std::tie(low, high) = get_hist_edges(hist, 2);
    map[key + "_low"] = low;
    map[key + "_high"] = high;
    map[key + "_center"] = (low + high) / 2.0;
    map[key + "_avg"] = hist->GetMean();
    map[key + "_rms"] = hist->GetRMS();
}