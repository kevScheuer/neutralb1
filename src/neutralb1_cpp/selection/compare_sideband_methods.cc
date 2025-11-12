/**
 * @file compare_sideband_methods.cc
 * @author Kevin Scheuer
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
 */


#include <iostream>
#include <map>
#include <vector>

#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TColor.h"

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
float OMEGA_MASS = 0.783; // PDG omega mass in GeV
float SIGNAL_WIDTH = 0.03; // exactly 3 sigma of PDG omega width
float SIDEBAND_GAP = 0.08;  // 9 sigma distance from signal region

// Forward declarations
TString join_keys(const std::map<TString, Int_t> &m, const TString &delimiter = ",");
std::vector<TH1F *> sideband_individual(
    std::map<TString, Int_t> &cut_color_map, 
    TString input_files);
// std::vector<TH1F *> sideband_2d(
//     std::map<TString, Int_t> &cut_color_map, 
//     TString input_files);

// TODO: should have a function for each method which takes in input files (data or MC)
// and returns the omega pi0 mass histogram after sideband subtraction


void compare_sideband_methods()
{
    TString input_data_files = "/lustre24/expphy/volatile/halld/home/kscheuer/FSRoot-skimmed-trees/best-chi2/tree_pi0pi0pippim__B4_bestChi2_SKIM_03_data.root";
    TString input_mc_files = "/lustre24/expphy/volatile/halld/home/kscheuer/FSRoot-skimmed-trees/best-chi2/tree_pi0pi0pippim__B4_bestChi2_SKIM_03_ver03.1_mc.root";
    setup(false);    
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();

    std::vector<TH1F *> h_2d_vector;
    std::vector<TH1F *> h_individual_vector;
    // h_2d_vector = sideband_2d(cut_color_map, input_data_files);
    h_individual_vector = sideband_individual(cut_color_map, input_data_files);

    // Plotting
    TCanvas *c = new TCanvas("c", "Sideband Subtraction Methods", 800, 600);

    // TODO: plot both together, including the omega mass, but for now just see if 
    // individual method works

    // plot omega mass, highglighting signal and sideband regions
    h_individual_vector[0]->SetLineColor(kBlack);
    h_individual_vector[0]->SetLineWidth(2);
    h_individual_vector[0]->SetXTitle("#omega Inv. Mass (GeV)");
    double bin_width_omega = get_bin_width(h_individual_vector[0]);
    h_individual_vector[0]->SetYTitle(TString::Format("Events / %.3f GeV", bin_width_omega));
    h_individual_vector[0]->SetTitle("");

    h_individual_vector[1]->SetLineColor(kGreen);
    h_individual_vector[1]->SetLineWidth(0);
    h_individual_vector[1]->SetFillColorAlpha(kGreen, 0.35);

    h_individual_vector[2]->SetLineColor(kRed);
    h_individual_vector[2]->SetLineWidth(0);
    h_individual_vector[2]->SetFillColorAlpha(kRed, 0.35);

    h_individual_vector[0]->Draw("HIST");
    h_individual_vector[1]->Draw("HIST SAME");
    h_individual_vector[2]->Draw("HIST SAME");

    c->SaveAs("sideband_individual_omega_mass.pdf");
    c->Clear();

    // plot data with cosmetic changes
    h_individual_vector[3]->SetLineColor(kGray);
    h_individual_vector[3]->SetLineWidth(2);
    h_individual_vector[3]->SetXTitle("#omega#pi^{0} Inv. Mass (GeV)");
    double bin_width = get_bin_width(h_individual_vector[3]);
    h_individual_vector[3]->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    h_individual_vector[3]->SetTitle("");

    h_individual_vector[4]->SetLineColor(kBlack);
    h_individual_vector[4]->SetLineWidth(2);    

    // add signal and sideband histograms
    h_individual_vector[5]->SetLineColor(kGreen);
    h_individual_vector[5]->SetLineWidth(1);
    h_individual_vector[5]->SetFillColorAlpha(kGreen, 0.35);
    h_individual_vector[6]->SetLineColor(kRed);
    h_individual_vector[6]->SetLineWidth(1);
    h_individual_vector[6]->SetFillColorAlpha(kRed, 0.35);

    h_individual_vector[3]->Draw("HIST");
    h_individual_vector[4]->Draw("HIST SAME");
    h_individual_vector[5]->Draw("HIST SAME");
    h_individual_vector[6]->Draw("HIST SAME");
    h_individual_vector[3]->SetMinimum(h_individual_vector[6]->GetMinimum() * 1.2);    

    c->SaveAs("sideband_individual_method.pdf");
    
}


// std::vector<TH1F *> sideband_2d(
//     std::map<TString, Int_t> &cut_color_map, 
//     TString input_files)
// {
//     TString cuts = join_keys(cut_color_map);


// }

std::vector<TH1F *> sideband_individual(
    std::map<TString, Int_t> &cut_color_map, 
    TString input_files)
{
    TString cuts = join_keys(cut_color_map);
    cuts += ",rf"; // will add RF cut to all plots too

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
    for(std::map<int, std::vector<TString>>::const_iterator perm_it = fs_permutation_to_order_map.begin();
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
            particles[pi0_omega].Data()
        );
        TString signal_region_cut = TString::Format(
            "abs(%s-%f)<%f",
            data_omega_mass.Data(),
            OMEGA_MASS,
            SIGNAL_WIDTH
        );
        
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
            SIGNAL_WIDTH
        );

        FSCut::defineCut(
            "omega_sb_subtraction",
            signal_region_cut,
            sideband_region_cut,
            1.0
        );

        // ==== omega mass ====
        // only broad cuts, since we don't want to see effect of SB subtraction here
        h_omega_mass[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            data_omega_mass,
            "(100,0.6,1.0)",
            TString::Format("CUT(%s)", cuts.Data()) 
        );
        h_omega_mass_sig[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            data_omega_mass,
            "(100,0.6,1.0)",
            TString::Format("CUT(%s)*CUTWT(omega_sb_subtraction)", cuts.Data()) 
        );
        h_omega_mass_sb[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            data_omega_mass,
            "(100,0.6,1.0)",
            TString::Format("CUT(%s)*CUTSBWT(omega_sb_subtraction)", cuts.Data()) 
        );

        // ==== omega pi0 mass ====
        TString omega_pi0_mass = TString::Format(
            "MASS(%s,%s,%s,%s)",
            particles[pi_plus].Data(),
            particles[pi_minus].Data(),
            particles[pi0_omega].Data(),
            particles[pi0_bachelor].Data()
        );
        TString omega_pi0_mass_cut = TString::Format(
            "%s>%f&&%s<%f",
            omega_pi0_mass.Data(),
            1.0, // min omega pi0 mass
            omega_pi0_mass.Data(),
            2.0  // max omega pi0 mass
        );

        FSCut::defineCut(
            "omega_pi0_mass_cut",
            omega_pi0_mass_cut
        );        

        // full omega pi0 mass distribution with standard cuts and the sideband 
        // subtraction applied
        h_omega_pi0_mass[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            omega_pi0_mass,
            "(100,1.0,2.0)",
            TString::Format("CUT(omega_pi0_mass_cut,%s)", cuts.Data())
        );
        h_omega_pi0_mass_result[perm_number -1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            omega_pi0_mass,
            "(100,1.0,2.0)",
            TString::Format("CUT(omega_pi0_mass_cut,omega_sb_subtraction,%s)", cuts.Data())
        );
        h_omega_pi0_mass_sig[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            omega_pi0_mass,
            "(100,1.0,2.0)",
            TString::Format("CUT(omega_pi0_mass_cut,%s)*CUTWT(omega_sb_subtraction)", cuts.Data())
        );
        h_omega_pi0_mass_sb[perm_number - 1] = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            omega_pi0_mass,
            "(100,1.0,2.0)",
            TString::Format("CUT(omega_pi0_mass_cut,%s)*CUTSBWT(omega_sb_subtraction)*(-1.0)", cuts.Data())
        );
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
        h_omega_pi0_mass_result_total,
        h_omega_pi0_mass_sig_total,
        h_omega_pi0_mass_sb_total
    };
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