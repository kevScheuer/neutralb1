/**
 * @file compare_sideband_methods.cc
 * @author Kevin Scheuer, Amy Schertz
 * @brief Compare different sideband subtraction methods for omega pi0
 *
 * A priori, there is no way to know which pi0 came from the omega, and which is the
 * "bachelor" pi0 from the resonance decay. Normally a sideband subtraction to pick the
 * omega and remove background involves a simple 1D subtraction by taking event from
 * sideband regions and using them to remove the background beneath the signal region.
 * However, in this case we have 2 pi0 candidates for the omega, complicating the
 * subtraction. We can either sideband subtract using a complicated 2D weighting scheme,
 * or we can simply sideband subtract each pi0 pi+ pi- combo independently. This file
 * compares these 2 methods to see how they affect the final omega pi0 mass
 * distribution.
 *
 * Authorship is also credited to Amy Schertz, as the 2d sideband code is directly
 * taken from her FSRoot tutorial for this channel:
 * https://github.com/JeffersonLab/gluex_workshops/blob/master/tutorial_2025/session2d/plots.C
 * 
 */

#include <iostream>
#include <map>
#include <vector>

#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLine.h"

#include "FSBasic/FSHistogram.h"
#include "FSBasic/FSCut.h"
#include "FSBasic/FSTree.h"
#include "FSMode/FSModeCollection.h"
#include "FSMode/FSModeHistogram.h"

#include "load_broad_cuts.cc"
#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"

TString NT("ntFSGlueX_MODECODE");
TString CATEGORY("pi0pi0pippim");

Int_t NORWEGIAN_BLUE = TColor::GetColor("#00205B");
Int_t NORWEGIAN_RED = TColor::GetColor("#BA0C2F");

// we make these TStrings for the FSCut definitions
double OMEGA_MASS = 0.78266;             // PDG omega mass in GeV
double SIGNAL_FULL_WIDTH = 0.168;      // full width of signal region in GeV (6 sigma)
double SIGNAL_HALF_WIDTH = SIGNAL_FULL_WIDTH / 2; // half width of signal region in GeV (3 sigma)
double SIGNAL_STDEV = SIGNAL_FULL_WIDTH / 6;      // "standard deviation" of signal region in GeV (1 sigma)

double SIGNAL_HIGH = OMEGA_MASS + SIGNAL_HALF_WIDTH;
double SIGNAL_LOW = OMEGA_MASS - SIGNAL_HALF_WIDTH;

// Forward declarations
std::vector<TH1F *> sideband_individual(
    std::map<TString, Int_t> &cut_color_map,
    TString input_files);
std::vector<TH1F *> sideband_2d(
    std::map<TString, Int_t> &cut_color_map,
    TString input_files,
    bool create_friends = false);

void compare_sideband_methods(bool mc=false, bool create_friend_trees=false)
{
    TString input_files;
    if (mc)
    {
        input_files = "/lustre24/expphy/volatile/halld/home/kscheuer/FSRoot-skimmed-trees/best-chi2/tree_pi0pi0pippim__B4_bestChi2_SKIM_03_ver03.1_mc.root";
    }
    else
    {
        input_files = "/lustre24/expphy/volatile/halld/home/kscheuer/FSRoot-skimmed-trees/best-chi2/tree_pi0pi0pippim__B4_bestChi2_SKIM_03_data.root";
    }        
    setup(false);
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();

    std::vector<TH1F *> h_2d_vector;
    std::vector<TH1F *> h_individual_vector;
    h_2d_vector = sideband_2d(cut_color_map, input_files, create_friend_trees);
    h_individual_vector = sideband_individual(cut_color_map, input_files);

    // Plotting
    TCanvas *c = new TCanvas("c", "Sideband Subtraction Methods", 800, 600);

    // plot omega mass, highlighting signal and sideband regions
    h_individual_vector[0]->SetLineColor(kGray);
    h_individual_vector[0]->SetLineWidth(2);
    h_individual_vector[0]->SetXTitle("#pi^{+}#pi^{-}#pi^{0}_{i} Inv. Mass (GeV)");
    double bin_width_omega = get_bin_width(h_individual_vector[0]);
    h_individual_vector[0]->SetYTitle(TString::Format("Events / %.3f GeV", bin_width_omega));
    h_individual_vector[0]->SetTitle("");

    h_individual_vector[1]->SetLineColor(NORWEGIAN_BLUE);
    h_individual_vector[1]->SetLineWidth(0);
    h_individual_vector[1]->SetFillColorAlpha(NORWEGIAN_BLUE, 0.35);

    h_individual_vector[2]->SetLineColor(NORWEGIAN_RED);
    h_individual_vector[2]->SetLineWidth(0);
    h_individual_vector[2]->SetFillColorAlpha(NORWEGIAN_RED, 0.35);
    h_individual_vector[0]->Draw("HIST");
    h_individual_vector[1]->Draw("HIST SAME");
    h_individual_vector[2]->Draw("HIST SAME");

    // draw lines for historical 2d window used in Amy's tutorial
    double historical_signal_low = 0.760;
    double historical_signal_high = 0.805;
    double historical_sideband_low_low = 0.690;
    double historical_sideband_low_high = 0.735;
    double historical_sideband_high_low = 0.830;
    double historical_sideband_high_high = 0.875;

    double max = h_individual_vector[0]->GetMaximum();
    double short_max = max * 0.8;
    TLine *line_historical_signal_low = new TLine(historical_signal_low, 0, historical_signal_low, max);
    TLine *line_historical_signal_high = new TLine(historical_signal_high, 0, historical_signal_high, max);
    TLine *line_historical_sideband_low_low = new TLine(historical_sideband_low_low, 0, historical_sideband_low_low, max);
    TLine *line_historical_sideband_low_high = new TLine(historical_sideband_low_high, 0, historical_sideband_low_high, max);
    TLine *line_historical_sideband_high_low = new TLine(historical_sideband_high_low, 0, historical_sideband_high_low, max);
    TLine *line_historical_sideband_high_high = new TLine(historical_sideband_high_high, 0, historical_sideband_high_high, max);
    line_historical_signal_low->SetLineColor(kBlack);
    line_historical_signal_low->SetLineWidth(1);
    line_historical_signal_high->SetLineColor(kBlack);
    line_historical_signal_high->SetLineWidth(1);
    line_historical_sideband_low_low->SetLineColor(kBlack);
    line_historical_sideband_low_low->SetLineWidth(1);
    line_historical_sideband_low_low->SetLineStyle(9);
    line_historical_sideband_low_high->SetLineColor(kBlack);
    line_historical_sideband_low_high->SetLineWidth(1);
    line_historical_sideband_low_high->SetLineStyle(9);
    line_historical_sideband_high_low->SetLineColor(kBlack);
    line_historical_sideband_high_low->SetLineWidth(1);
    line_historical_sideband_high_low->SetLineStyle(9);
    line_historical_sideband_high_high->SetLineColor(kBlack);
    line_historical_sideband_high_high->SetLineWidth(1);
    line_historical_sideband_high_high->SetLineStyle(9);
    line_historical_signal_low->Draw("SAME");
    line_historical_signal_high->Draw("SAME");
    line_historical_sideband_low_low->Draw("SAME");
    line_historical_sideband_low_high->Draw("SAME");
    line_historical_sideband_high_low->Draw("SAME");
    line_historical_sideband_high_high->Draw("SAME");

    // draw lines for adjusted 2d window to be more comparable to our regions    
    double sideband_low_low = SIGNAL_LOW - SIGNAL_HALF_WIDTH;
    double sideband_low_high = SIGNAL_LOW;    
    double sideband_high_low = SIGNAL_HIGH;
    double sideband_high_high = SIGNAL_HIGH + SIGNAL_HALF_WIDTH;

    TLine *line_signal_low = new TLine(SIGNAL_LOW, 0, SIGNAL_LOW, short_max);
    TLine *line_signal_high = new TLine(SIGNAL_HIGH, 0, SIGNAL_HIGH, short_max);
    TLine *line_sideband_low_low = new TLine(sideband_low_low, 0, sideband_low_low, short_max);
    TLine *line_sideband_low_high = new TLine(sideband_low_high, 0, sideband_low_high, short_max);
    TLine *line_sideband_high_low = new TLine(sideband_high_low, 0, sideband_high_low, short_max);
    TLine *line_sideband_high_high = new TLine(sideband_high_high, 0, sideband_high_high, short_max);
    line_signal_low->SetLineColor(kYellow - 1);
    line_signal_low->SetLineWidth(1);
    line_signal_high->SetLineColor(kYellow - 1);
    line_signal_high->SetLineWidth(1);
    line_sideband_low_low->SetLineColor(kYellow - 1);
    line_sideband_low_low->SetLineWidth(1);
    line_sideband_low_low->SetLineStyle(9);
    line_sideband_low_high->SetLineColor(kYellow - 1);
    line_sideband_low_high->SetLineWidth(1);
    line_sideband_low_high->SetLineStyle(9);
    line_sideband_high_low->SetLineColor(kYellow - 1);
    line_sideband_high_low->SetLineWidth(1);
    line_sideband_high_low->SetLineStyle(9);
    line_sideband_high_high->SetLineColor(kYellow - 1);
    line_sideband_high_high->SetLineWidth(1);
    line_sideband_high_high->SetLineStyle(9);
    line_signal_low->Draw("SAME");
    line_signal_high->Draw("SAME");
    line_sideband_low_low->Draw("SAME");
    line_sideband_low_high->Draw("SAME");
    line_sideband_high_low->Draw("SAME");
    line_sideband_high_high->Draw("SAME");

    // Finally, the legend
    TLegend *legend_omega = new TLegend(0.65, 0.75, 0.95, 0.95);
    legend_omega->SetBorderSize(1);
    legend_omega->AddEntry(h_individual_vector[0], "Data", "l");
    legend_omega->AddEntry(h_individual_vector[1], "Signal", "f");
    legend_omega->AddEntry(h_individual_vector[2], "Sidebands", "f");
    legend_omega->AddEntry(line_signal_low, "Adjusted 2D Signal", "l");
    legend_omega->AddEntry(line_sideband_low_low, "Adjusted 2D Sidebands", "l");
    legend_omega->AddEntry(line_historical_signal_low, "Historical 2D Signal", "l");
    legend_omega->AddEntry(line_historical_sideband_low_low, "Historical 2D Sidebands", "l");
    legend_omega->Draw();

    TString output_suffix = mc ? "_mc" : "_data";
    c->SaveAs("sideband_omega_mass" + output_suffix + ".pdf");
    c->Clear();

    // plot historical 2d sideband method result
    h_2d_vector[0]->SetLineColor(kBlack);
    h_2d_vector[0]->SetLineWidth(2);
    h_2d_vector[0]->SetXTitle("#pi^{+}#pi^{-}#pi^{0}_{1}#pi^{0}_{2} Inv. Mass (GeV)");
    double bin_width_2d = get_bin_width(h_2d_vector[0]);
    h_2d_vector[0]->SetYTitle(TString::Format("Events / %.3f GeV", bin_width_2d));
    h_2d_vector[0]->SetTitle("");

    // plot 2d sideband with adjusted windows to match my selection regions
    h_2d_vector[1]->SetLineColor(kYellow - 1);
    h_2d_vector[1]->SetLineWidth(1);

    // plot data with cosmetic changes
    h_individual_vector[3]->SetLineColor(kGray);
    h_individual_vector[3]->SetLineWidth(2);

    // plot norwegian result
    h_individual_vector[4]->SetLineColor(kMagenta);
    h_individual_vector[4]->SetLineWidth(1);

    h_2d_vector[0]->Draw("HIST");
    h_2d_vector[1]->Draw("HIST SAME");
    h_individual_vector[3]->Draw("HIST SAME");
    h_individual_vector[4]->Draw("HIST SAME");
    h_2d_vector[0]->SetMinimum(h_individual_vector[4]->GetMinimum() * 1.1);
    h_2d_vector[0]->SetMaximum(h_individual_vector[3]->GetMaximum() * 1.1);

    TLegend *legend = new TLegend(0.65, 0.75, 0.95, 0.95);
    legend->SetBorderSize(0);
    legend->AddEntry(h_individual_vector[3], "Data (no SB)", "l");
    legend->AddEntry(h_2d_vector[0], "Historical 2D Result", "l");
    legend->AddEntry(h_2d_vector[1], "Adjusted 2D Result", "l");
    legend->AddEntry(h_individual_vector[4], "Norwegian Result", "l");
    legend->Draw();

    c->SaveAs("sideband" + output_suffix + ".pdf");
}

/**
 * @brief Perform the historical 2D sideband subtraction method, and adjusted version
 * 
 * Adjusted refers to doing the method of the historical version, but adjusting the 
 * signal and sideband widths/regions to better match the norwegian one
 * 
 * @param cut_color_map map of broad cuts to colors
 * @param input_files files to plot from
 * @param create_friends whether to create friend trees
 * @return std::vector<TH1F *> vector of histograms
 */
std::vector<TH1F *> sideband_2d(
    std::map<TString, Int_t> &cut_color_map,
    TString input_files,
    bool create_friends = false)
{
    TString cuts = join_keys(cut_color_map);

    if (create_friends)
    {
        vector<pair<TString, TString>> friendTreeContents;
        friendTreeContents.push_back(pair<TString, TString>("M234", "MASS(2,3,4)"));
        friendTreeContents.push_back(pair<TString, TString>("M235", "MASS(2,3,5)"));
        FSTree::createFriendTree(input_files, "ntFSGlueX_100_112", "M3PI", friendTreeContents);
    }
    FSTree::addFriendTree("M3PI");

    // ____Historical 2D sideband subtraction method____
    // copied from Amy Schertz's FSRoot tutorial

    // cuts for the signal region
    TString SIG4("(M234>0.760&&M234<0.805)");
    TString SIG5("(M235>0.760&&M235<0.805)");
    TString SIG4orSIG5("((" + SIG4 + ")||(" + SIG5 + "))");

    // cuts for the sideband regions
    TString SB4("(((M234>0.690&&M234<0.735)||"
                "(M234>0.830&&M234<0.875))&&!(" +
                SIG4 + ")&&!(" + SIG5 + "))");
    TString SB5("(((M235>0.690&&M235<0.735)||"
                "(M235>0.830&&M235<0.875))&&!(" +
                SIG4 + ")&&!(" + SIG5 + "))");
    TString SB4andSB5("((" + SB4 + ")&&(" + SB5 + "))");
    TString SB4orSB5("((" + SB4 + ")||(" + SB5 + "))");
    TString SB4xorSB5("(((" + SB4 + ")||(" + SB5 + "))&&!(" + SB4andSB5 + "))");

    // event weights
    TString OMEGAWT("(1.00*(" + SIG4orSIG5 + ")-0.50*(" + SB4xorSB5 + ")-1.25*(" + SB4andSB5 + "))");

    TH1F *h_omega_pi0_mass = FSModeHistogram::getTH1F(
        input_files,
        NT,
        CATEGORY,
        "MASS(2,3,4,5)",
        "(100,1.0,2.0)",
        TString::Format("%s*CUT(%s)*CUTWT(rf)", OMEGAWT.Data(), cuts.Data()));

    // ____the same process, but now adjustd to my analysis windows____
    SIG4 = TString::Format("(M234>%f&&M234<%f)", SIGNAL_LOW, SIGNAL_HIGH);
    SIG5 = TString::Format("(M235>%f&&M235<%f)", SIGNAL_LOW, SIGNAL_HIGH);
    SIG4orSIG5 = "((" + SIG4 + ")||(" + SIG5 + "))";

    // cuts for the sideband regions

    // my analysis uses sidebands that start at the end of the signal region i.e.
    // SIGNAL_HIGH to SIGNAL_HIGH + SIGNAL_HALF_WIDTH on the high side
    // SIGNAL_LOW - SIGNAL_HALF_WIDTH to SIGNAL_LOW on the low side. To better compare
    // to the 2D method of weighting, I'll double the width of the sideband regions to
    // match the full width of the signal region

    SB4 = TString::Format(
        "(((M234>%f&&M234<%f)||(M234>%f&&M234<%f))&&!(%s)&&!(%s))",
        SIGNAL_LOW - SIGNAL_FULL_WIDTH, SIGNAL_LOW,
        SIGNAL_HIGH, SIGNAL_HIGH + SIGNAL_FULL_WIDTH,
        SIG4.Data(), SIG5.Data()
    );
    SB5 = TString::Format(
        "(((M235>%f&&M235<%f)||(M235>%f&&M235<%f))&&!(%s)&&!(%s))",
        SIGNAL_LOW - SIGNAL_FULL_WIDTH, SIGNAL_LOW,
        SIGNAL_HIGH, SIGNAL_HIGH + SIGNAL_FULL_WIDTH,
        SIG4.Data(), SIG5.Data()
    );
    SB4andSB5 = "((" + SB4 + ")&&(" + SB5 + "))";
    SB4orSB5 = "((" + SB4 + ")||(" + SB5 + "))";
    SB4xorSB5 = "(((" + SB4 + ")||(" + SB5 + "))&&!(" + SB4andSB5 + "))";

    // event weights
    OMEGAWT = "(1.00*(" + SIG4orSIG5 + ")-0.50*(" + SB4xorSB5 + ")-1.25*(" + SB4andSB5 + "))";

    TH1F *h_omega_pi0_mass_adjusted = FSModeHistogram::getTH1F(
        input_files,
        NT,
        CATEGORY,
        "MASS(2,3,4,5)",
        "(100,1.0,2.0)",
        TString::Format("%s*CUT(%s)*CUTWT(rf)", OMEGAWT.Data(), cuts.Data()));

    return std::vector<TH1F *>{h_omega_pi0_mass, h_omega_pi0_mass_adjusted};
}

std::vector<TH1F *> sideband_individual(
    std::map<TString, Int_t> &cut_color_map,
    TString input_files)
{
    TString cuts = join_keys(cut_color_map);

    // Our particles are given in the order
    //   0       1          2             3              4              5
    // [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (omega)] [pi0 (bachelor)]
    //
    // We need to define each permutation of the pi0
    std::map<int, std::vector<TString>> fs_permutation_to_order_map = {
        {1, {"B", "1", "2", "3", "4", "5"}}, // particle 4 = pi0 from omega
        {2, {"B", "1", "2", "3", "5", "4"}}  // particle 5 = pi0 from omega
    };

    // setup histograms for each permutation
    TH1F *h_omega_mass[2], *h_omega_mass_sig[2], *h_omega_mass_sb[2];
    TH1F *h_omega_pi0_mass[2], *h_omega_pi0_mass_result[2];

    // Begin loop over permutations
    for (std::map<int, std::vector<TString>>::const_iterator perm_it = fs_permutation_to_order_map.begin();
         perm_it != fs_permutation_to_order_map.end(); ++perm_it)
    {
        int perm_number = perm_it->first;
        std::vector<TString> particles = perm_it->second;
        int beam = 0;
        int proton = 1;
        int pi_plus = 2;
        int pi_minus = 3;
        int pi0_omega = 4;
        int pi0_bachelor = 5;

        // ==== omega sideband subtraction cut definition ====
        TString data_omega_mass = TString::Format(
            "MASS(%s,%s,%s)",
            particles[pi_plus].Data(),
            particles[pi_minus].Data(),
            particles[pi0_omega].Data());
        TString signal_region_cut = TString::Format(
            "abs(%s-%f)<%f",
            data_omega_mass.Data(),
            OMEGA_MASS,
            SIGNAL_HALF_WIDTH);        
        
        // We'll use two sidebands in the region
        // (OMEGA_MASS +/- SIGNAL_HALF_WIDTH) to 
        // +/- (OMEGA_MASS +/- SIGNAL_HALF_WIDTH + SIGNAL_HALF_WIDTH)
        // This way, the sideband regions together have the same width as the signal
        // region and we don't have to handle any weighting factors
        TString sideband_region_cut = TString::Format(
            "abs(%s-%f)>(%f)&&abs(%s-%f)<(%f+%f)",
            data_omega_mass.Data(),
            OMEGA_MASS,
            SIGNAL_HALF_WIDTH,
            data_omega_mass.Data(),
            OMEGA_MASS,
            SIGNAL_HALF_WIDTH,
            SIGNAL_HALF_WIDTH);
        FSCut::defineCut(
            "omega_sb_subtraction",
            signal_region_cut,
            sideband_region_cut,
            1.0);

        // ==== omega mass ====
        // only broad cuts, since we don't want to see effect of SB subtraction here
        h_omega_mass[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            data_omega_mass,
            "(100,0.6,1.0)",
            TString::Format("CUT(%s)*CUTWT(rf)", cuts.Data()));
        h_omega_mass_sig[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            data_omega_mass,
            "(100,0.6,1.0)",
            TString::Format("CUT(%s,omega_sb_subtraction)*CUTWT(rf)", cuts.Data()));
        h_omega_mass_sb[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            data_omega_mass,
            "(100,0.6,1.0)",
            TString::Format("CUT(%s)*CUTSB(omega_sb_subtraction)*CUTWT(rf)", cuts.Data()));

        // ==== omega pi0 mass ====
        TString omega_pi0_mass = TString::Format(
            "MASS(%s,%s,%s,%s)",
            particles[pi_plus].Data(),
            particles[pi_minus].Data(),
            particles[pi0_omega].Data(),
            particles[pi0_bachelor].Data());
        TString omega_pi0_mass_cut = TString::Format(
            "%s>%f&&%s<%f",
            omega_pi0_mass.Data(),
            1.0, // min omega pi0 mass
            omega_pi0_mass.Data(),
            2.0 // max omega pi0 mass
        );

        FSCut::defineCut(
            "omega_pi0_mass_cut",
            omega_pi0_mass_cut);

        // full omega pi0 mass distribution with historical cuts and the sideband
        // subtraction applied
        h_omega_pi0_mass[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            omega_pi0_mass,
            "(100,1.0,2.0)",
            TString::Format("CUT(omega_pi0_mass_cut,%s)*CUTWT(rf)", cuts.Data()));
        h_omega_pi0_mass_result[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            omega_pi0_mass,
            "(100,1.0,2.0)",
            TString::Format("CUT(omega_pi0_mass_cut,%s)*CUTWT(omega_sb_subtraction,rf)", cuts.Data()));
    } // end loop over permutations

    // Now combine the two permutations' histograms
    TH1F *h_omega_mass_total = (TH1F *)h_omega_mass[0]->Clone("h_omega_mass_total");
    h_omega_mass_total->Add(h_omega_mass[1]);
    TH1F *h_omega_mass_sig_total = (TH1F *)h_omega_mass_sig[0]->Clone("h_omega_mass_sig_total");
    h_omega_mass_sig_total->Add(h_omega_mass_sig[1]);
    TH1F *h_omega_mass_sb_total = (TH1F *)h_omega_mass_sb[0]->Clone("h_omega_mass_sb_total");
    h_omega_mass_sb_total->Add(h_omega_mass_sb[1]);
    TH1F *h_omega_pi0_mass_total = (TH1F *)h_omega_pi0_mass[0]->Clone("h_omega_pi0_mass_total");
    h_omega_pi0_mass_total->Add(h_omega_pi0_mass[1]);
    TH1F *h_omega_pi0_mass_result_total = (TH1F *)h_omega_pi0_mass_result[0]->Clone("h_omega_pi0_mass_result_total");
    h_omega_pi0_mass_result_total->Add(h_omega_pi0_mass_result[1]);

    return std::vector<TH1F *>{
        h_omega_mass_total,
        h_omega_mass_sig_total,
        h_omega_mass_sb_total,
        h_omega_pi0_mass_total,
        h_omega_pi0_mass_result_total};
}