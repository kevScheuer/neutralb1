/**
 * @file find_t_bin_stats.cc
 * @author Kevin Scheuer
 * @brief Quickly determine statistics in each t-bin for finalized trees
 * 
 * Generally want roughly equal statistics in each bin for PWA fits, so this
 * utility helps determine the appropriate bin edges.
 */

#include <iostream>

#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TH1D.h"


void find_t_bin_stats(
    const std::string& input_signal,
    const std::string& input_background,
    const std::string& tree_name, 
    const float& low_bin1,
    const float& low_bin2,
    const float& low_bin3,
    const float& low_bin4,
    const float& high_bin4
)
{
    TFile *f_signal = TFile::Open(input_signal.c_str(), "READ");
    if (!f_signal || f_signal->IsZombie())
    {
        std::cerr << "Error: Could not open file " << input_signal << "\n";
        return;
    }
    TFile *f_background = TFile::Open(input_background.c_str(), "READ");
    if (!f_background || f_background->IsZombie())
    {
        std::cerr << "Error: Could not open file " << input_background << "\n";
        f_signal->Close();
        return;
    }
    TTree *signal_tree = dynamic_cast<TTree*>(f_signal->Get(tree_name.c_str()));
    if (!signal_tree)
    {
        std::cerr << "Error: Could not find tree " << tree_name << " in file " << input_signal << "\n";
        f_signal->Close();
        f_background->Close();
        return;
    }
    TTree *background_tree = dynamic_cast<TTree*>(f_background->Get(tree_name.c_str()));
    if (!background_tree)
    {
        std::cerr << "Error: Could not find tree " << tree_name << " in file " << input_background << "\n";
        f_signal->Close();
        f_background->Close();
        return;
    }
    
    std::vector<double> t_bin_edges = { // example edges
        low_bin1, // 0.1
        low_bin2, // 0.2
        low_bin3, // 0.3 
        low_bin4, // 0.5
        high_bin4 // 1.0
    };

    for (size_t i = 0; i < t_bin_edges.size() - 1; ++i)
    {
        double t_low = t_bin_edges[i];
        double t_high = t_bin_edges[i + 1];

        TString t_cut = TString::Format("(t>%f && t<%f && M4Pi>1.0 && M4Pi<1.8)", t_low, t_high);
        signal_tree->Draw(
            TString::Format("t>>h_sig(1, %f, %f)", t_low, t_high), 
            t_cut, 
            "goff"
        );
        TH1D* h_signal = (TH1D*)gDirectory->Get("h_sig");

        // add weight branch for background
        t_cut = TString::Format("(t>%f && t<%f && M4Pi>1.0 && M4Pi<1.8)*weight", t_low, t_high);
        background_tree->Draw(
            TString::Format("t>>h_bkg(1, %f, %f)", t_low, t_high), 
            t_cut, 
            "goff"
        );
        TH1D* h_background = (TH1D*)gDirectory->Get("h_bkg");

        h_signal->Add(h_background);
        double n_entries = h_signal->Integral(0, h_signal->GetNbinsX() + 1);

        std::cout << "t-bin [" << t_low << ", " << t_high << "]: " << static_cast<int>(n_entries) << " entries\n";
    }

    f_signal->Close();
    f_background->Close();
}