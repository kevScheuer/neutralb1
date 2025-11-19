/**
 * @file sideband_variation.cc
 * @author Kevin Scheuer    
 * @brief Plot omega pi0 mass distributions with different sideband widths
 * 
 * Utilize the norwegian method of sideband subtraction to plot the omega pi0 mass
 * distribution using different sideband widths to see how it affects the final
 * distribution.
 */

#include <iostream>
#include <tuple>
#include <vector>

#include "TH1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TArrow.h"

#include "FSBasic/FSHistogram.h"
#include "FSBasic/FSCut.h"
#include "FSBasic/FSTree.h"
#include "FSMode/FSModeHistogram.h"

#include "load_broad_cuts.cc"
#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"

// forward declarations
void plot_sideband_variation(
    TString NT,
    TString CATEGORY,
    TString cuts,
    TString input_files,
    double omega_stdev,
    double signal_half_width,
    double sideband_full_width,
    std::vector<int> sideband_starts,
    int period,
    bool mc);


void sideband_variation(int period, bool mc = false)
{
    TString input_files;
    if (mc)
    {
        input_files = TString::Format(
            "/lustre24/expphy/volatile/halld/home/kscheuer/"
            "FSRoot-skimmed-trees/best-chi2/"
            "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_ver03.1_mc.root",
            period);
    }
    else
    {
        input_files = TString::Format(
            "/lustre24/expphy/volatile/halld/home/kscheuer/"
            "FSRoot-skimmed-trees/best-chi2/"
            "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_data.root",
            period);
    }
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(true);
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();
    TString cuts = join_keys(cut_color_map);

    // our signal stdev will be the Voigtian FWHM / 2sqrt(2ln2), using the voigtian 
    // fit results from previous studies
    double omega_stdev = 0.020;
    double signal_half_width = omega_stdev * 3; // 3 sigma on each side

    // setup the sideband locations to study. Remember that each sideband is 1/2 as wide
    // as the entire signal region. So the sideband regions will cover from
    // N*sigma to (N+3)*sigma
    double sideband_full_width = signal_half_width;
    std::vector<int> sideband_starts = {3, 4, 5, 6}; // in units of omega_stdev

    plot_sideband_variation(
        NT,
        CATEGORY,
        cuts,
        input_files,
        omega_stdev,
        signal_half_width,
        sideband_full_width,
        sideband_starts,
        period,
        mc);

    return;
}


void  plot_sideband_variation(
    TString NT,
    TString CATEGORY,
    TString cuts,
    TString input_files,
    double omega_stdev,
    double signal_half_width,
    double sideband_full_width,
    std::vector<int> sideband_starts,
    int period,
    bool mc) 
{        
    const double omega_mass = 0.78266; // PDG omega mass in GeV    
    std::vector<int> sideband_colors = {
        kRed+3, kRed+2, kRed+1, kRed
    };
    TCanvas *c = new TCanvas("c", "c", 800, 600);

    // Our particles are given in the order
    //   0       1          2             3              4              5
    // [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (omega)] [pi0 (bachelor)]
    //
    // we have two possible pi0 assignments for the omega decay
    TString data_omega_mass_1 = "MASS(2,3,4)";
    TString data_omega_mass_2 = "MASS(2,3,5)";

    // ==== Draw omega histogram for both combinations ====
     TH1F *h_omega_mass_1 = FSModeHistogram::getTH1F(
        input_files,
        NT,
        CATEGORY,
        data_omega_mass_1,
        "(100, 0.6, 1.0)",
        TString::Format(
            "CUT(%s)&&CUTWT(rf)",
            cuts.Data()));
    TH1F *h_omega_mass_2 = FSModeHistogram::getTH1F(
        input_files,
        NT,
        CATEGORY,
        data_omega_mass_2,
        "(100, 0.6, 1.0)",
        TString::Format(
            "CUT(%s)&&CUTWT(rf)",
            cuts.Data()));
        
    // add both permutations together
    h_omega_mass_1->Add(h_omega_mass_2);

    // ==== Highlight the signal region ====
    TH1F* h_signal = (TH1F *)h_omega_mass_1->Clone("h_signal");
    h_signal->SetFillColorAlpha(kGreen + 1, 0.5);
    h_signal->SetFillStyle(1001);
    h_signal->SetLineWidth(0); // No line for the fill histogram

    // Zero out bins outside the selection region
    double signal_low_edge = omega_mass - signal_half_width;
    double signal_high_edge = omega_mass + signal_half_width;
    for (int j = 1; j <= h_signal->GetNbinsX(); ++j)
    {
        double bin_center = h_signal->GetBinCenter(j);
        if (bin_center < signal_low_edge || bin_center > signal_high_edge)
        {
            h_signal->SetBinContent(j, 0);
        }
    }

    // draw a horizontal line at 0
    TLine *zero_line = new TLine(h_omega_mass_1->GetXaxis()->GetXmin(), 0, h_omega_mass_1->GetXaxis()->GetXmax(), 0);
    zero_line->SetLineColor(kBlack);
    zero_line->SetLineWidth(1);    
    
    double line_spacing = h_omega_mass_1->GetMaximum() * 0.03; // 3% of max height spacing between lines    

    for (int i=0; i<sideband_starts.size(); ++i)
    {        
        double sideband_low_edge = omega_mass + sideband_starts[i] * omega_stdev;
        double sideband_high_edge = sideband_low_edge + sideband_full_width;

        // ==== Draw lines for each sideband region ====
        
        // Set histogram minimum to accommodate the horizontal lines below zero        
        double new_min = -(6 * line_spacing);
        h_omega_mass_1->SetMinimum(new_min);
        
        // Create horizontal lines for the sideband region
        double line_y = -line_spacing - (i * line_spacing);
        TArrow *sideband_line = new TArrow(
            sideband_low_edge, 
            line_y, 
            sideband_high_edge, 
            line_y, 
            0.05, 
            "<|>");
        sideband_line->SetLineColor(sideband_colors[i]);
        sideband_line->SetLineWidth(2);
        
        h_omega_mass_1->SetXTitle("#pi^{+}#pi^{-}#pi^{0}_{i} inv. mass (GeV)");
        double bin_width = get_bin_width(h_omega_mass_1);
        h_omega_mass_1->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
        
        // Draw everything
        if (i == 0) {
            h_omega_mass_1->Draw("HIST");
            h_signal->Draw("HIST SAME");
            zero_line->Draw("SAME");
        }
        sideband_line->Draw("SAME");
    }
    FSHistogram::dumpHistogramCache();


    c->SaveAs(TString::Format("sideband_variation_%d%s.pdf",
        period,
        mc ? "mc" : "data"
    ));
}