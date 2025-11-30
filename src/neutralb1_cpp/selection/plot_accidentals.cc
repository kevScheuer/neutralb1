/**
 * @file plot_accidentals.cc
 * @author Kevin Scheuer
 * @brief Plot RF accidentals and their weights
 * 
 * Has to be separate from plot_cuts.cc due to the "og" histograms already having their
 * RF prompt peak cut applied.
 */

#include <tuple>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TH1.h"
#include "TLegend.h"
#include "TString.h"

#include "FSMode/FSModeHistogram.h"
#include "FSBasic/FSCut.h"

#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"

/**
 * @brief Plot the rf accidentals and their weights
 * 
 * @param input_data_files data files to plot from
 * @param NT name of the tree in the ROOT file
 * @param CATEGORY category name in the ROOT file
 * @param period data taking period
 */
void plot_accidentals()
{
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(false);

    TString input_data3 = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/best-chi2/"
        "tree_pi0pi0pippim__B4_bestChi2_SKIM_03_data.root";
    TString input_data4 = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/best-chi2/"
        "tree_pi0pi0pippim__B4_bestChi2_SKIM_04_data.root";
    TString input_data5 = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/best-chi2/"
        "tree_pi0pi0pippim__B4_bestChi2_SKIM_05_data.root";
        
    std::vector<TString> data_files = {
        input_data3,
        input_data4,
        input_data5
    };
    std::vector <TH1F*> og_hists, prompt_hists, sideband_hists;

    FSCut::defineCut("rf", "abs(RFDeltaT)<2.0", "abs(RFDeltaT)>2.0", 0.125);

    // loop over run periods and get respective histograms
    for (const auto& input_data : data_files)
    {
        // no cut applied
        TH1F *h_og = FSModeHistogram::getTH1F(
        input_data,
        NT,
        CATEGORY,
        "RFDeltaT",
        "(400,-20.0,20.0)",
        "");
        // sideband region with weight
        TH1F *h_sideband = FSModeHistogram::getTH1F(
            input_data,
            NT,
            CATEGORY,
            "RFDeltaT",
            "(400,-20.0,20.0)",
            "CUTSBWT(rf)");
        // prompt peak region
        TH1F *h_prompt = FSModeHistogram::getTH1F(
            input_data,
            NT,
            CATEGORY,
            "RFDeltaT",
            "(40,-2.0,2.0)",
            "CUT(rf)");

        og_hists.push_back(h_og);
        prompt_hists.push_back(h_prompt);
        sideband_hists.push_back(h_sideband);
    }

    // combine histograms from all run periods
    for (size_t i = 1; i < og_hists.size(); ++i)
    {
        og_hists[0]->Add(og_hists[i]);
        prompt_hists[0]->Add(prompt_hists[i]);
        sideband_hists[0]->Add(sideband_hists[i]);
    }

    TCanvas *c_rf = new TCanvas("c_rf", "Accidental Subtraction", 800, 600);
    TH1F *h_data = og_hists[0];
    TH1F *h_prompt = prompt_hists[0];
    TH1F *h_sideband = sideband_hists[0];

    h_data->SetLineColor(kGray);
    h_data->SetLineWidth(2);

    h_prompt->SetLineWidth(0);
    h_prompt->SetFillColorAlpha(kGreen + 1, 0.5);
    h_prompt->SetFillStyle(1001);

    h_sideband->SetLineWidth(0);
    h_sideband->SetFillColorAlpha(TColor::GetColor("#c849a9"), 0.5);
    h_sideband->SetFillStyle(1001);

    h_data->Draw("HIST");
    h_prompt->Draw("HIST SAME");
    h_sideband->Draw("HIST SAME");


    TLegend *legend_acc = new TLegend(0.65, 0.65, 0.88, 0.88);
    legend_acc->AddEntry(h_data, "Original Data", "l");
    legend_acc->AddEntry(h_prompt, "In-Time Signal Region", "f");
    legend_acc->AddEntry(h_sideband, "Weighted Out-of-Time Combos", "f");
    legend_acc->Draw();

    h_data->SetTitle("");
    h_data->SetXTitle("RF #DeltaT (ns)");
    double bin_width = get_bin_width(h_data);
    h_data->SetYTitle(TString::Format("Combos / %.3f ns", bin_width));
    c_rf->Update();
    c_rf->SaveAs("GlueXI_Accidental_Subtraction.pdf");        
}