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

#include "../load_broad_cuts.cc"
#include "neutralb1/fit_utils.h"
#include "../fsroot_setup.cc"

// forward declarations
std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>
plot_sideband_variation(
    TString NT,
    TString CATEGORY,
    TString cuts,
    TString input_files,
    double signal_stdev,
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
    std::tie(NT, CATEGORY) = setup(false);
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();
    TString cuts = join_keys(cut_color_map);

    // to define signal and sideband points, we'll use the signal width from the 
    // fits to data in omega_fits.cc. Then to establish sideband regions, we'll
    // use multiples of the signal region standard deviation (assuming 6 sigma
    // covers the full width)
    double signal_full_width = 0.056;
    double signal_half_width = signal_full_width / 2;
    double signal_stdev = signal_full_width / 6;

    // setup the sideband locations to study. Remember that each sideband is 1/2 as wide
    // as the entire signal region. So the sideband regions will cover from
    // N*sigma to (N+3)*sigma
    double sideband_full_width = signal_half_width;
    std::vector<int> sideband_starts = {3, 4, 5, 6}; // in units of signal_stdev

    std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> sideband_edges =
        plot_sideband_variation(
            NT,
            CATEGORY,
            cuts,
            input_files,
            signal_stdev,
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
    double signal_stdev,
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
        double sideband_right_low_edge = omega_mass + sideband_starts[i] * signal_stdev;
        double sideband_right_high_edge = sideband_right_low_edge + sideband_full_width;

        double sideband_left_low_edge = omega_mass - sideband_starts[i] * signal_stdev - sideband_full_width;
        double sideband_left_high_edge = sideband_left_low_edge + sideband_full_width;

        std::cout << "Sideband " << i << ": Left [" << sideband_left_low_edge << ", " << sideband_left_high_edge
                  << "], Right [" << sideband_right_low_edge << ", " << sideband_right_high_edge << "]" << std::endl;

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
        "(100, 0.4, 1.200)",
        TString::Format(
            "CUT(%s)*CUTWT(rf)",
            cuts.Data()));
    TH1F *h_omega_mass_2 = FSModeHistogram::getTH1F(
        input_files,
        NT,
        CATEGORY,
        data_omega_mass_2,
        "(100, 0.4, 1.200)",
        TString::Format(
            "CUT(%s)*CUTWT(rf)",
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

    // For accurate comparison, we need to plot both permutations together on the same
    // omega pi0 plot. We'll copy the method in friend_reorder to have a 1:1 methodology
    std::map<int, std::vector<TString>> fs_permutation_to_order_map = {
        {1, {"B", "1", "2", "3", "4", "5"}}, // particle 4 = pi0 from omega
        {2, {"B", "1", "2", "3", "5", "4"}}  // particle 5 = pi0 from omega
    };

    TH1F *h_omega_pi0[2], *h_signal_highlight[2];
    std::vector<TH1F *> sideband_subtracted_histograms;    

    for (std::map<int, std::vector<TString>>::const_iterator perm_it = fs_permutation_to_order_map.begin();
         perm_it != fs_permutation_to_order_map.end(); ++perm_it)
    {
        int perm = perm_it->first;
        std::vector<TString> perm_particles = perm_it->second;

        int beam = 0;
        int proton = 1;
        int pi_plus = 2;
        int pi_minus = 3;
        int pi0_omega = 4;
        int pi0_bachelor = 5;

        TString data_omega_mass = TString::Format(
            "MASS(%s,%s,%s)",
            perm_particles[pi_plus].Data(),
            perm_particles[pi_minus].Data(),
            perm_particles[pi0_omega].Data());

        TString data_omega_pi0_mass = TString::Format(
            "MASS(%s,%s,%s,%s)",            
            perm_particles[pi_plus].Data(),
            perm_particles[pi_minus].Data(),
            perm_particles[pi0_omega].Data(),
            perm_particles[pi0_bachelor].Data());

        // get the omega pi0 mass distribution for this permutation, with broad cuts applied
        h_omega_pi0[perm-1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            data_omega_pi0_mass,
            "(100, 1.0, 2.0)",
            TString::Format(
                "CUT(%s)*CUTWT(rf)",
                cuts.Data()));
        
        // define the constant signal region
        TString signal_region = TString::Format(
            "abs(%s-%f)<%f",
            data_omega_mass.Data(),
            omega_mass,
            signal_half_width
        );

        // now loop over each sideband variation
        for(int sb_iter=0; sb_iter < sideband_edges.size(); ++sb_iter)
        {
            double sideband_left_low_edge = sideband_edges[sb_iter].first.first;
            double sideband_left_high_edge = sideband_edges[sb_iter].first.second;
            double sideband_right_low_edge = sideband_edges[sb_iter].second.first;
            double sideband_right_high_edge = sideband_edges[sb_iter].second.second;

            double sideband_gap = sideband_right_low_edge - omega_mass;
            double sideband_width = sideband_right_high_edge - sideband_right_low_edge;        

            // check that sidebands are symmetric about the omega mass
            if (std::abs((omega_mass - sideband_left_high_edge) - sideband_gap) > 1e-6)
            {
                std::cerr << "Warning: Sideband gap inconsistency detected!" << std::endl;
            }
            if (std::abs((sideband_left_high_edge - sideband_left_low_edge) - sideband_width) > 1e-6)
            {
                std::cerr << "Warning: Sideband width inconsistency detected!" << std::endl;
            }

            TString sideband_region_cut = TString::Format(
                "abs(%s-%f)>(%f)&&abs(%s-%f)<(%f+%f)",
                data_omega_mass.Data(),
                omega_mass,
                sideband_gap,
                data_omega_mass.Data(),
                omega_mass,
                sideband_gap,
                sideband_width);

            std::cout << "Sideband variation " << sb_iter << ": "
                      << sideband_region_cut
                      << std::endl;

            FSCut::defineCut("selection", signal_region, sideband_region_cut, 1.0);

            TH1F *h_sideband_subtracted = FSModeHistogram::getTH1F(
                input_files,
                NT,
                CATEGORY,
                data_omega_pi0_mass,
                "(100, 1.0, 2.0)",
                TString::Format(
                    "CUT(%s)*CUTWT(rf, selection)",
                    cuts.Data()));        
            
            // add the permutations together
            if(perm == 1)
            {
                sideband_subtracted_histograms.push_back(h_sideband_subtracted);                
            }
            else
            {
                sideband_subtracted_histograms[sb_iter]->Add(h_sideband_subtracted);
            }            
        }
        // this can be done after FSCut defined because it selects a constant signal 
        // region
        h_signal_highlight[perm-1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            data_omega_pi0_mass,
            "(100, 1.0, 2.0)",
            TString::Format(
                "CUT(%s,selection)*CUTWT(rf)",                
                cuts.Data()));
    } // end loop over permutations

    // add the total and signal histograms together
    h_omega_pi0[0]->Add(h_omega_pi0[1]);
    h_signal_highlight[0]->Add(h_signal_highlight[1]);

    // plot
    h_omega_pi0[0]->SetXTitle("#omega#pi^{0} inv. mass (GeV)");
    h_omega_pi0[0]->SetMinimum(0);
    h_omega_pi0[0]->SetLineColor(kBlack);
    double bin_width = get_bin_width(h_omega_pi0[0]);
    h_omega_pi0[0]->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    h_omega_pi0[0]->SetTitle("");
    h_omega_pi0[0]->Draw("HIST");

    h_signal_highlight[0]->SetLineWidth(1);
    h_signal_highlight[0]->SetLineColor(kGreen + 1);    
    h_signal_highlight[0]->Draw("HIST SAME");

    for (int i = 0; i < sideband_subtracted_histograms.size(); ++i)
    {
        sideband_subtracted_histograms[i]->SetLineWidth(1);
        sideband_subtracted_histograms[i]->SetLineColor(sideband_colors[i]);
        sideband_subtracted_histograms[i]->Draw("HIST SAME");
    }
    FSHistogram::dumpHistogramCache();
    c->SaveAs(TString::Format("omega_pi0_sideband_variation_%d%s.pdf",
                              period,
                              mc ? "_mc" : "_data"));
    return; 
}