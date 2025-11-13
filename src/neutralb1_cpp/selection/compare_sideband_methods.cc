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

// we make these TStrings for the FSCut definitions
float OMEGA_MASS = 0.7826;             // PDG omega mass in GeV
float SIGNAL_WIDTH = 0.00868 * 3;      // exactly 3 sigma of PDG omega width
float SIDEBAND_GAP = SIGNAL_WIDTH * 3; // 9 sigma distance from signal region

// Forward declarations
TString join_keys(const std::map<TString, Int_t> &m, const TString &delimiter = ",");
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

    // plot omega mass, highglighting signal and sideband regions
    h_individual_vector[0]->SetLineColor(kGray);
    h_individual_vector[0]->SetLineWidth(2);
    h_individual_vector[0]->SetXTitle("#pi^{+}#pi^{-}#pi^{0}_{i} Inv. Mass (GeV)");
    double bin_width_omega = get_bin_width(h_individual_vector[0]);
    h_individual_vector[0]->SetYTitle(TString::Format("Events / %.3f GeV", bin_width_omega));
    h_individual_vector[0]->SetTitle("");

    h_individual_vector[1]->SetLineColor(kBlue);
    h_individual_vector[1]->SetLineWidth(0);
    h_individual_vector[1]->SetFillColorAlpha(kBlue, 0.35);

    h_individual_vector[2]->SetLineColor(kRed - 2);
    h_individual_vector[2]->SetLineWidth(0);
    h_individual_vector[2]->SetFillColorAlpha(kRed - 2, 0.35);
    h_individual_vector[0]->Draw("HIST");
    h_individual_vector[1]->Draw("HIST SAME");
    h_individual_vector[2]->Draw("HIST SAME");

    // draw lines for standard 2d window used in Amy's tutorial
    double standard_signal_low = 0.760;
    double standard_signal_high = 0.805;
    double standard_sideband_low_low = 0.690;
    double standard_sideband_low_high = 0.735;
    double standard_sideband_high_low = 0.830;
    double standard_sideband_high_high = 0.875;

    double max = h_individual_vector[0]->GetMaximum();
    double short_max = max * 0.8;
    TLine *line_standard_signal_low = new TLine(standard_signal_low, 0, standard_signal_low, max);
    TLine *line_standard_signal_high = new TLine(standard_signal_high, 0, standard_signal_high, max);
    TLine *line_standard_sideband_low_low = new TLine(standard_sideband_low_low, 0, standard_sideband_low_low, max);
    TLine *line_standard_sideband_low_high = new TLine(standard_sideband_low_high, 0, standard_sideband_low_high, max);
    TLine *line_standard_sideband_high_low = new TLine(standard_sideband_high_low, 0, standard_sideband_high_low, max);
    TLine *line_standard_sideband_high_high = new TLine(standard_sideband_high_high, 0, standard_sideband_high_high, max);
    line_standard_signal_low->SetLineColor(kBlack);
    line_standard_signal_low->SetLineWidth(1);
    line_standard_signal_high->SetLineColor(kBlack);
    line_standard_signal_high->SetLineWidth(1);
    line_standard_sideband_low_low->SetLineColor(kBlack);
    line_standard_sideband_low_low->SetLineWidth(1);
    line_standard_sideband_low_low->SetLineStyle(9);
    line_standard_sideband_low_high->SetLineColor(kBlack);
    line_standard_sideband_low_high->SetLineWidth(1);
    line_standard_sideband_low_high->SetLineStyle(9);
    line_standard_sideband_high_low->SetLineColor(kBlack);
    line_standard_sideband_high_low->SetLineWidth(1);
    line_standard_sideband_high_low->SetLineStyle(9);
    line_standard_sideband_high_high->SetLineColor(kBlack);
    line_standard_sideband_high_high->SetLineWidth(1);
    line_standard_sideband_high_high->SetLineStyle(9);
    line_standard_signal_low->Draw("SAME");
    line_standard_signal_high->Draw("SAME");
    line_standard_sideband_low_low->Draw("SAME");
    line_standard_sideband_low_high->Draw("SAME");
    line_standard_sideband_high_low->Draw("SAME");
    line_standard_sideband_high_high->Draw("SAME");

    // draw lines for adjusted 2d window to be more comparable to our regions
    double signal_low = OMEGA_MASS - SIGNAL_WIDTH;
    double signal_high = OMEGA_MASS + SIGNAL_WIDTH;
    double sideband_low_center = OMEGA_MASS - SIDEBAND_GAP - SIGNAL_WIDTH * 0.5;
    double sideband_low_low = sideband_low_center - SIGNAL_WIDTH;
    double sideband_low_high = sideband_low_center + SIGNAL_WIDTH;
    double sideband_high_center = OMEGA_MASS + SIDEBAND_GAP + SIGNAL_WIDTH * 0.5;
    double sideband_high_low = sideband_high_center - SIGNAL_WIDTH;
    double sideband_high_high = sideband_high_center + SIGNAL_WIDTH;

    TLine *line_signal_low = new TLine(signal_low, 0, signal_low, short_max);
    TLine *line_signal_high = new TLine(signal_high, 0, signal_high, short_max);
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
    legend_omega->AddEntry(line_standard_signal_low, "Standard 2D Signal", "l");
    legend_omega->AddEntry(line_standard_sideband_low_low, "Standard 2D Sidebands", "l");
    legend_omega->Draw();

    TString output_suffix = mc ? "_mc" : "_data";
    c->SaveAs("sideband_omega_mass" + output_suffix + ".pdf");
    c->Clear();

    // plot standard 2d sideband method result
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

    h_individual_vector[4]->SetLineColor(kBlue);
    h_individual_vector[4]->SetLineWidth(1);

    // add signal and sideband histograms
    h_individual_vector[5]->SetLineColor(kViolet);
    h_individual_vector[5]->SetLineWidth(1);
    h_individual_vector[6]->SetLineColor(kRed - 2);
    h_individual_vector[6]->SetLineWidth(1);
    h_individual_vector[6]->SetFillColorAlpha(kRed - 2, 0.35);

    h_2d_vector[0]->Draw("HIST");
    h_2d_vector[1]->Draw("HIST SAME");
    h_individual_vector[3]->Draw("HIST SAME");
    h_individual_vector[4]->Draw("HIST SAME");
    h_individual_vector[5]->Draw("HIST SAME");
    h_individual_vector[6]->Draw("HIST SAME");
    h_2d_vector[0]->SetMinimum(h_individual_vector[6]->GetMinimum() * 1.2);
    h_2d_vector[0]->SetMaximum(h_individual_vector[3]->GetMaximum() * 1.2);

    TLegend *legend = new TLegend(0.65, 0.75, 0.95, 0.95);
    legend->SetBorderSize(0);
    legend->AddEntry(h_individual_vector[3], "Data (no SB)", "l");
    legend->AddEntry(h_2d_vector[0], "Standard 2D Result", "l");
    legend->AddEntry(h_2d_vector[1], "Adjusted 2D Result", "l");
    legend->AddEntry(h_individual_vector[4], "Simple Signal", "l");
    legend->AddEntry(h_individual_vector[6], "Simple Sidebands", "f");
    legend->AddEntry(h_individual_vector[5], "Simple Result", "l");
    legend->Draw();

    c->SaveAs("sideband" + output_suffix + ".pdf");
}

std::vector<TH1F *> sideband_2d(
    std::map<TString, Int_t> &cut_color_map,
    TString input_files,
    bool create_friends = false)
{
    TString cuts = join_keys(cut_color_map);
    cuts += ",rf"; // will add RF cut to all plots too

    if (create_friends)
    {
        vector<pair<TString, TString>> friendTreeContents;
        friendTreeContents.push_back(pair<TString, TString>("M234", "MASS(2,3,4)"));
        friendTreeContents.push_back(pair<TString, TString>("M235", "MASS(2,3,5)"));
        FSTree::createFriendTree(input_files, "ntFSGlueX_100_112", "M3PI", friendTreeContents);
    }
    FSTree::addFriendTree("M3PI");

    // ____Amy's 2D sideband subtraction method____
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
    SIG4 = "(M234>0.757&&M234<0.809)";
    SIG5 = "(M235>0.757&&M235<0.809)";
    SIG4orSIG5 = "((" + SIG4 + ")||(" + SIG5 + "))";

    // cuts for the sideband regions

    // for this next part, I'll refer to the full width of the omega as being
    // +/- 3 sigma from the peak (so full width = 6 sigma). In terms of the constants,
    // this means OMEGA_MASS +/- SIGNAL_WIDTH

    // my analysis covers the 9 sigma - 12 sigma sidebands (each sideband is 1/2 full
    // omega width) and the 2D method sidebands are 6 sigma wide (each one covers full
    // omega width). I'll make the new regions 10.5 sigma +/- 3 sigma. That way the
    // way the widths are the same and I can keep the weighting scheme, while
    // covering a comparable region
    SB4 = "(((M234>0.666&&M234<0.718)||(M234>0.848&&M234<0.900))&&!(" + SIG4 + ")&&!(" + SIG5 + "))";
    SB5 = "(((M235>0.666&&M235<0.718)||(M235>0.848&&M235<0.900))&&!(" + SIG4 + ")&&!(" + SIG5 + "))";
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
    TH1F *h_omega_pi0_mass[2], *h_omega_pi0_mass_result[2], *h_omega_pi0_mass_sig[2], *h_omega_pi0_mass_sb[2];

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
            SIGNAL_WIDTH);

        // Instead of using two sidebands centered at +/- SIDEBAND_GAP from the signal
        // region, we'll use two sidebands in the region
        // +/- SIDEBAND_GAP to +/- (SIDEBAND_GAP + SIGNAL_WIDTH)
        // This way, the sideband regions together have the same width as the signal
        // region and we don't have to handle any weighting factors
        TString sideband_region_cut = TString::Format(
            "abs(%s-%f)>(%f)&&abs(%s-%f)<(%f+%f)",
            data_omega_mass.Data(),
            OMEGA_MASS,
            SIDEBAND_GAP,
            data_omega_mass.Data(),
            OMEGA_MASS,
            SIDEBAND_GAP,
            SIGNAL_WIDTH);

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

        // full omega pi0 mass distribution with standard cuts and the sideband
        // subtraction applied
        h_omega_pi0_mass[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            omega_pi0_mass,
            "(100,1.0,2.0)",
            TString::Format("CUT(omega_pi0_mass_cut,%s)*CUTWT(rf)", cuts.Data()));
        h_omega_pi0_mass_sig[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            omega_pi0_mass,
            "(100,1.0,2.0)",
            TString::Format("CUT(omega_pi0_mass_cut,omega_sb_subtraction,%s)*CUTWT(rf)", cuts.Data()));
        h_omega_pi0_mass_result[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            omega_pi0_mass,
            "(100,1.0,2.0)",
            TString::Format("CUT(omega_pi0_mass_cut,%s)*CUTWT(omega_sb_subtraction,rf)", cuts.Data()));
        h_omega_pi0_mass_sb[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            omega_pi0_mass,
            "(100,1.0,2.0)",
            TString::Format("CUT(omega_pi0_mass_cut,%s)*CUTSB(omega_sb_subtraction)*CUTWT(rf)*(-1.0)", cuts.Data()));
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
    TH1F *h_omega_pi0_mass_sig_total = (TH1F *)h_omega_pi0_mass_sig[0]->Clone("h_omega_pi0_mass_sig_total");
    h_omega_pi0_mass_sig_total->Add(h_omega_pi0_mass_sig[1]);
    TH1F *h_omega_pi0_mass_sb_total = (TH1F *)h_omega_pi0_mass_sb[0]->Clone("h_omega_pi0_mass_sb_total");
    h_omega_pi0_mass_sb_total->Add(h_omega_pi0_mass_sb[1]);

    return std::vector<TH1F *>{
        h_omega_mass_total,
        h_omega_mass_sig_total,
        h_omega_mass_sb_total,
        h_omega_pi0_mass_total,
        h_omega_pi0_mass_sig_total,
        h_omega_pi0_mass_result_total,
        h_omega_pi0_mass_sb_total};
}

/**
 * @brief Join the keys of a color map into a single TString separated by a delimiter
 *
 * @param m map of cut TStrings to color integers
 * @param delimiter delimiter to separate keys, default is ","
 * @return TString joined keys
 */
TString join_keys(const std::map<TString, Int_t> &m, const TString &delimiter = ",")
{
    TString result;
    for (auto it = m.begin(); it != m.end(); ++it)
    {
        if (it != m.begin())
            result += delimiter;
        result += it->first;
    }
    return result;
}