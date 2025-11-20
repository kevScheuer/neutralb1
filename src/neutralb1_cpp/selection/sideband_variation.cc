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
#include <utility>
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
std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>
plot_sideband_variation(
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
void plot_omega_pi0_sideband_variations(
    TString NT,
    TString CATEGORY,
    TString cuts,
    TString input_files,
    std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> sideband_edges,
    double signal_half_width,
    int period,
    bool mc);
TH1F *get_omega_mass(
    TString input_files,
    TString NT,
    TString CATEGORY,
    TString cuts);
TH1F *get_signal_region(
    TH1F *h_omega,
    double omega_mass,
    double signal_half_width);
std::pair<double, double> get_signal_edges(
    double omega_mass,
    double signal_half_width);

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

    // our signal stdev will use the voigtian fit results from previous studies
    double omega_stdev = 0.020;
    double signal_half_width = omega_stdev * 3; // 3 sigma on each side

    // setup the sideband locations to study. Remember that each sideband is 1/2 as wide
    // as the entire signal region. So the sideband regions will cover from
    // N*sigma to (N+3)*sigma
    double sideband_full_width = signal_half_width;
    std::vector<int> sideband_starts = {3, 4, 5, 6}; // in units of omega_stdev

    std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> sideband_edges =
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

    plot_omega_pi0_sideband_variations(
        NT,
        CATEGORY,
        cuts,
        input_files,
        sideband_edges,
        signal_half_width,
        period,
        mc);

    return;
}

std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>
plot_sideband_variation(
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
        kRed + 3, kRed + 2, kRed + 1, kRed};
    TCanvas *c = new TCanvas("c", "c", 800, 600);

    TH1F *h_omega = get_omega_mass(input_files, NT, CATEGORY, cuts);            // full omega mass histogram
    TH1F *h_signal = get_signal_region(h_omega, omega_mass, signal_half_width); // highlight signal region

    // draw a horizontal line at 0
    TLine *zero_line = new TLine(h_omega->GetXaxis()->GetXmin(), 0, h_omega->GetXaxis()->GetXmax(), 0);
    zero_line->SetLineColor(kBlack);
    zero_line->SetLineWidth(1);

    double line_spacing = h_omega->GetMaximum() * 0.03; // 3% of max height spacing between lines

    // vector to store left/right sideband pair, and their pairs of low/high edges
    std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> sideband_edges;
    for (int i = 0; i < sideband_starts.size(); ++i)
    {
        double sideband_right_low_edge = omega_mass + sideband_starts[i] * omega_stdev;
        double sideband_right_high_edge = sideband_right_low_edge + sideband_full_width;

        double sideband_left_low_edge = omega_mass - sideband_starts[i] * omega_stdev - sideband_full_width;
        double sideband_left_high_edge = sideband_left_low_edge + sideband_full_width;

        sideband_edges.push_back(
            {{sideband_left_low_edge, sideband_left_high_edge},
             {sideband_right_low_edge, sideband_right_high_edge}});

        // ==== Draw lines for each sideband region ====

        // Set histogram minimum to accommodate the horizontal lines below zero
        double new_min = -(6 * line_spacing);
        h_omega->SetMinimum(new_min);

        // Create horizontal lines for the sideband region
        double line_y = -line_spacing - (i * line_spacing);
        TArrow *sideband_right_line = new TArrow(
            sideband_right_low_edge,
            line_y,
            sideband_right_high_edge,
            line_y,
            0.01,
            "<|>");
        sideband_right_line->SetLineColor(sideband_colors[i]);
        sideband_right_line->SetLineWidth(2);
        sideband_right_line->SetFillStyle(1001);
        sideband_right_line->SetFillColor(sideband_colors[i]);

        TArrow *sideband_left_line = new TArrow(
            sideband_left_low_edge,
            line_y,
            sideband_left_high_edge,
            line_y,
            0.01,
            "<|>");
        sideband_left_line->SetLineColor(sideband_colors[i]);
        sideband_left_line->SetLineWidth(2);
        sideband_left_line->SetFillStyle(1001);
        sideband_left_line->SetFillColor(sideband_colors[i]);

        // Draw everything
        if (i == 0)
        {
            h_omega->SetXTitle("#pi^{+}#pi^{-}#pi^{0}_{i} inv. mass (GeV)");
            double bin_width = get_bin_width(h_omega);
            h_omega->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
            h_omega->SetTitle("");
            h_omega->Draw("HIST");
            h_signal->Draw("HIST SAME");
            zero_line->Draw("SAME");
        }
        sideband_right_line->Draw();
        sideband_left_line->Draw();
    }
    FSHistogram::dumpHistogramCache();

    c->SaveAs(TString::Format("sideband_variation_%d%s.pdf",
                              period,
                              mc ? "_mc" : "_data"));

    return sideband_edges;
}

TH1F *get_omega_mass(TString input_files, TString NT, TString CATEGORY, TString cuts)
{
    // Our particles are given in the order
    //   0       1          2             3              4              5
    // [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (omega)] [pi0 (bachelor)]
    //
    // we have two possible pi0 assignments for the omega decay
    TString data_omega_mass_1 = "MASS(2,3,4)";
    TString data_omega_mass_2 = "MASS(2,3,5)";

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

    return h_omega_mass_1;
}

TH1F *get_signal_region(
    TH1F *h_omega,
    double omega_mass,
    double signal_half_width)
{
    TH1F *h_signal = (TH1F *)h_omega->Clone("h_signal");
    h_signal->SetFillColorAlpha(kGreen + 1, 0.5);
    h_signal->SetFillStyle(1001);
    h_signal->SetLineWidth(0); // No line for the fill histogram

    // Zero out bins outside the selection region
    double signal_low_edge, signal_high_edge;
    std::tie(signal_low_edge, signal_high_edge) = get_signal_edges(omega_mass, signal_half_width);
    for (int j = 1; j <= h_signal->GetNbinsX(); ++j)
    {
        double bin_center = h_signal->GetBinCenter(j);
        if (bin_center < signal_low_edge || bin_center > signal_high_edge)
        {
            h_signal->SetBinContent(j, 0);
        }
    }
    return h_signal;
}

std::pair<double, double> get_signal_edges(
    double omega_mass,
    double signal_half_width)
{
    double signal_low_edge = omega_mass - signal_half_width;
    double signal_high_edge = omega_mass + signal_half_width;
    return {signal_low_edge, signal_high_edge};
}

void plot_omega_pi0_sideband_variations(
    TString NT,
    TString CATEGORY,
    TString cuts,
    TString input_files,
    std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> sideband_edges,
    double signal_half_width,
    int period,
    bool mc)
{
    const double omega_mass = 0.78266; // PDG omega mass in GeV
    TCanvas *c = new TCanvas("c_omega_pi0", "c_omega_pi0", 800, 600);    
    std::vector<int> sideband_colors = {
        kRed + 3, kRed + 2, kRed + 1, kRed};

    // first lets get the omega pi0 mass distribution after all cuts, but no omega selections
    TH1F *h_omega_pi0 = FSModeHistogram::getTH1F(
        input_files,
        NT,
        CATEGORY,
        "MASS(2,3,4,5)",
        "(100, 1.0, 2.0)",
        TString::Format(
            "CUT(%s)*CUTWT(rf)",
            cuts.Data()));

    h_omega_pi0->SetXTitle("#omega#pi^{0} inv. mass (GeV)");
    h_omega_pi0->SetLineColor(kBlack);
    double bin_width = get_bin_width(h_omega_pi0);
    h_omega_pi0->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    h_omega_pi0->SetTitle("");

    double signal_low_edge, signal_high_edge;
    std::tie(signal_low_edge, signal_high_edge) = get_signal_edges(omega_mass, signal_half_width);

    // next we'll obtain just the signal region histogram for omega pi0
    TString signal_cut = TString::Format(
        "((MASS(2,3,4)>%f&&MASS(2,3,4)<%f)||(MASS(2,3,5)>%f&&MASS(2,3,5)<%f))",
        signal_low_edge,
        signal_high_edge,
        signal_low_edge,
        signal_high_edge);
    TH1F *h_signal_highlight = FSModeHistogram::getTH1F(
        input_files,
        NT,
        CATEGORY,
        "MASS(2,3,4,5)",
        "(100, 1.0, 2.0)",
        TString::Format(
            "%s*CUT(%s)*CUTWT(rf)",
            signal_cut.Data(),
            cuts.Data()));
    h_signal_highlight->SetLineWidth(1);
    h_signal_highlight->SetLineColor(kGreen + 1);    

    // now we can loop over the different sideband selections
    std::vector<TH1F *> sideband_subtracted_histograms;
    std::vector<TH1F *> sideband_highlight_histograms;
    for(const auto& edges : sideband_edges)
    {
        double sideband_left_low_edge = edges.first.first;
        double sideband_left_high_edge = edges.first.second;
        double sideband_right_low_edge = edges.second.first;
        double sideband_right_high_edge = edges.second.second;

        // since the sidebands are defined by the arrows from before, its easier just 
        // to store those values and use them directly here instead of the simpler 
        // abs() method.
        TString low_sb_cut = TString::Format(
            "(MASS(2,3,4)>%f&&MASS(2,3,4)<%f)||(MASS(2,3,5)>%f&&MASS(2,3,5)<%f)",
            sideband_left_low_edge,
            sideband_left_high_edge,
            sideband_left_low_edge,
            sideband_left_high_edge);
        TString high_sb_cut = TString::Format(
            "(MASS(2,3,4)>%f&&MASS(2,3,4)<%f)||(MASS(2,3,5)>%f&&MASS(2,3,5)<%f)",
            sideband_right_low_edge,
            sideband_right_high_edge,
            sideband_right_low_edge,
            sideband_right_high_edge);

        // this reads: "If either 3pi permutation is in the low or high sideband, and 
        // neither is in the signal region"
        TString sideband_cut = TString::Format(
            "(%s||%s)&&(!%s)",
            low_sb_cut.Data(),
            high_sb_cut.Data(),
            signal_cut.Data());

        FSCut::defineCut("selection", signal_cut, sideband_cut, 1.0);

        TH1F *h_sideband_subtracted = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            "MASS(2,3,4,5)",
            "(100, 1.0, 2.0)",
            TString::Format(
                "CUT(%s)*CUTWT(rf)*CUTWT(selection)",
                cuts.Data()));        
        h_sideband_subtracted->SetLineWidth(1);
        h_sideband_subtracted->SetLineColor(sideband_colors[&edges - &sideband_edges[0]]);
        sideband_subtracted_histograms.push_back(h_sideband_subtracted);

        TH1F *h_sideband_highlight = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            "MASS(2,3,4,5)",
            "(100, 1.0, 2.0)",
            TString::Format(
                "CUT(%s)*CUTWT(rf)*CUTSB(selection)*(-1)",
                cuts.Data()));
        h_sideband_highlight->SetLineWidth(1);
        h_sideband_highlight->SetLineColor(sideband_colors[&edges - &sideband_edges[0]]);
        sideband_highlight_histograms.push_back(h_sideband_highlight);
    }

    // draw everything
    h_omega_pi0->Draw("HIST");
    h_signal_highlight->Draw("HIST SAME");
    for (int i = 0; i < sideband_subtracted_histograms.size(); ++i)
    {
        sideband_subtracted_histograms[i]->Draw("HIST SAME");
    }
    FSHistogram::dumpHistogramCache();
    c->SaveAs(TString::Format("omega_pi0_sideband_variation_%d%s.pdf",
                              period,
                              mc ? "_mc" : "_data"));
    // NOTE: for now the negative sideband highlights are not drawn to reduce clutter    
}