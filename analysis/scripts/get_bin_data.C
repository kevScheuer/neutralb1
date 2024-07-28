/*Output bin content of data to a csv file

This circumvents the need to edit the data reader of vecps_plotter to handle
the mass binning. This file instead reads the AmpTools Tree, and makes the same TEM
cuts I used. I then bin it to match the fit bins, and weight using
the "Weight" leaf of the Tree.
*/

#include <string>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1F.h"

/*
    bin_width: in GeV
*/
void get_bin_data(TString file, double bin_width, std::string f_name = "",
                  std::string min_recoil_pi_mass = "0.0",
                  std::string t_low = "0.4", std::string t_high = "0.5",
                  std::string E_low = "8.2", std::string E_high = "8.8",
                  std::string m_low = "1.0", std::string m_high = "1.5")
{
    // get file and tree
    std::unique_ptr<TFile> f(TFile::Open(file));
    if (!f)
    {
        cout << "File DNE" << "\n";
        exit(1);
    }

    TTree *tree = f->Get<TTree>("kin");
    if (!tree)
    {
        cout << "Tree DNE" << "\n";
        exit(1);
    }
    // Create list of cuts and apply to tree
    std::string t_range[2] = {t_low, t_high};
    std::string E_range[2] = {E_low, E_high};
    std::string m_recoil_pi_cut = "MRecoilPi>" + min_recoil_pi_mass;
    std::string t_cut = "t>" + t_range[0] + " && t<" + t_range[1];
    std::string E_cut = "E_Beam>" + E_range[0] + " && E_Beam<" + E_range[1];
    std::string selection = "Weight*(" +
                            m_recoil_pi_cut + " && " + t_cut + " && " + E_cut + ")";

    // setup mass bin details
    std::string mass_range[2] = {m_low, m_high};
    int n_bins = (std::stod(mass_range[1]) -
                  std::stod(mass_range[0])) /
                 bin_width;

    std::string h_str = "h(" + std::to_string(n_bins) + "," + mass_range[0] +
                        "," + mass_range[1] + ")";

    TCanvas *c1 = new TCanvas("c1", "c1", 1800, 1800);

    // draw histogram with correct bins & cuts, then grab it from pad
    tree->Draw(("M4Pi>>" + h_str).c_str(), selection.c_str());
    TH1F *h = (TH1F *)gPad->GetPrimitive("h");

    // loop over the mass bins to get the bin average and rms
    double bin_low = std::stod(mass_range[0]);
    double bin_high = std::stod(mass_range[0]) + bin_width;
    std::vector<double> bin_average_vec;
    std::vector<double> bin_rms_vec;
    for (int i = 0; i < n_bins; i++)
    {
        // make each mass bin a hist of 100 bins to get reasonable mean and rms values
        tree->Draw(("M4Pi>>m(100," + std::to_string(bin_low) + "," + std::to_string(bin_high) + ")").c_str(), selection.c_str());
        TH1F *m = (TH1F *)gPad->GetPrimitive("m");
        bin_average_vec.push_back(m->GetMean());
        bin_rms_vec.push_back(m->GetRMS());
        bin_low += bin_width;
        bin_high += bin_width;
    }

    // create csv file
    ofstream csv;
    if (f_name == "")
    {
        std::string f_name = "bin_data_t_" + t_range[0] + "-" + t_range[1] + "_" +
                             "m_" + mass_range[0] + "-" + mass_range[1] + "-" +
                             std::to_string(bin_width) + ".csv";
    }
    csv.open(f_name);

    // write header of bins to csv
    csv << "low_edge,high_edge,bin_contents,bin_error,mean,rms\n";

    // write bin contents to csv
    for (int bin = 1; bin <= h->GetNbinsX(); bin++)
    {
        csv << h->GetBinLowEdge(bin) << "," << h->GetBinLowEdge(bin) + h->GetBinWidth(bin) << "," << h->GetBinContent(bin) << "," << h->GetBinError(bin) << "," << bin_average_vec[bin - 1] << "," << bin_rms_vec[bin - 1] << "\n";
    }

    return;
}