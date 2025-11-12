/**
 * @file plot_broad_cuts.cc
 * @author Kevin Scheuer
 * @brief First step to select broad cuts in pi0pi0pi+pi- final state
 *
 * This file creates plots of several variables with broad cuts to select the physical
 * region of interest for a pi0pi0pi+pi- final state. Selection of omega signal
 * events is not yet done here, and so no combinatorics are studied yet.
 *
 */

#include <iostream>
#include <map>

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
std::map<TString, TString> CUT_TO_LEGEND = {
    {"unusedE", "Unused Shower Energy"},
    {"unusedTracks", "Unused Charged Tracks"},
    {"z", "Production Vertex"},
    {"MM2", "Missing Mass^{2}"},
    {"eBeam", "Beam Energy"},
    {"chi2", "#chi^{2} / NDF"},
    {"t", "-t"},
    {"shQuality", "Shower Quality"}};

// Forward declarations
TCanvas *setup_canvas(bool logy = false);
TH1F *plot_variable(
    TCanvas *c,
    const TString input_data_files,
    const TString cut_variable,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const double cut_lower_bound,
    const double cut_upper_bound,
    const std::map<TString, Int_t> &cut_color_map);

void plot_broad_cuts(
    bool mc = true,
    bool bggen = false,
    bool read_cache = false,
    bool dump_cache = false)
{
    TString input_data_files;
    TString input_mc_files;

    input_data_files = "/lustre24/expphy/volatile/halld/home/kscheuer/FSRoot-skimmed-trees/best-chi2/tree_pi0pi0pippim__B4_bestChi2_SKIM_03_data.root";
    input_mc_files = "/lustre24/expphy/volatile/halld/home/kscheuer/FSRoot-skimmed-trees/best-chi2/tree_pi0pi0pippim__B4_bestChi2_SKIM_03_ver03.1_mc.root";

    // gROOT->ProcessLine(".L glueXstyle.C");
    // gROOT->ProcessLine("gluex_style();");
    // gStyle->SetOptStat(0); // force stats box off

    setup(read_cache);

    // load our broad cuts
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();
    double bin_width;

    // =================== Accidental Subtraction ===================
    // Before we look at any other cuts let's look at the accidental subtraction itself
    // to see how well it performs. All other cuts will have it applied
    TCanvas *c_rf = new TCanvas("c_rf", "Accidental Subtraction", 800, 600);
    TH1F *h_acc_data = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        "RFDeltaT",
        "(400,-20.0,20.0)",
        "");
    TH1F *h_acc_subtracted = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        "RFDeltaT",
        "(400,-20.0,20.0)",
        "CUTSBWT(rf)");
    TH1F *h_acc_signal = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        "RFDeltaT",
        "(40,-2.0,2.0)",
        "CUTWT(rf)");
    h_acc_data->SetLineColor(kGray);
    h_acc_data->SetLineWidth(2);

    h_acc_signal->SetLineWidth(0);
    h_acc_signal->SetFillColorAlpha(kGreen + 1, 0.5);
    h_acc_signal->SetFillStyle(1001);

    h_acc_subtracted->SetLineWidth(0);
    h_acc_subtracted->SetFillColorAlpha(TColor::GetColor("#c849a9"), 0.5);
    h_acc_subtracted->SetFillStyle(1001);

    h_acc_data->Draw("HIST");
    h_acc_subtracted->Draw("HIST SAME");
    h_acc_signal->Draw("HIST SAME");

    TLegend *legend_acc = new TLegend(0.65, 0.65, 0.88, 0.88);
    legend_acc->AddEntry(h_acc_data, "Original Data", "l");
    legend_acc->AddEntry(h_acc_signal, "In-Time Signal Region", "f");
    legend_acc->AddEntry(h_acc_subtracted, "Weighted Out-of-Time Combos", "f");
    legend_acc->Draw();

    h_acc_data->SetXTitle("RF #DeltaT (ns)");    
    bin_width = get_bin_width(h_acc_data);
    h_acc_data->SetYTitle(TString::Format("Combos / %.3f ns", bin_width));
    c_rf->Update();
    c_rf->SaveAs("Accidental_Subtraction.pdf");
    delete c_rf;

    TCanvas *c;
    // =================== Unused Shower Energy ===================
    c = setup_canvas(true);
    TH1F *h;
    h = plot_variable(
        c,
        input_data_files,
        "unusedE",
        "EnUnusedSh",
        "100",
        "0.0",
        "1.0",
        0.0,
        0.1,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("Unused Shower Energy (GeV)");
    h->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    c->Update();
    c->SaveAs("Unused_Shower_Energy.pdf");
    delete c;

    c = setup_canvas(true);    
    h = plot_variable(
        c,
        input_mc_files,
        "unusedE",
        "EnUnusedSh",
        "100",
        "0.0",
        "1.0",
        0.0,
        0.1,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("Unused Shower Energy (GeV)");
    h->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    c->Update();
    c->SaveAs("Unused_Shower_Energy_MC.pdf");
    delete c;

    // =================== Number Unused Tracks ===================
    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_data_files,
        "unusedTracks",
        "NumUnusedTracks",
        "6",
        "0.0",
        "6.0",
        0.0,
        1.0,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("Number of Unused Tracks");
    h->SetYTitle("Events / Track");
    c->Update();
    c->SaveAs("Number_Unused_Tracks.pdf");
    delete c;

    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_mc_files,
        "unusedTracks",
        "NumUnusedTracks",
        "6",
        "0.0",
        "6.0",
        0.0,
        1.0,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("Number of Unused Tracks");
    h->SetYTitle("Events / Track");
    c->Update();
    c->SaveAs("Number_Unused_Tracks_MC.pdf");
    delete c;

    // =================== Production Vertex ===================
    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_data_files,
        "z",
        "ProdVz",
        "100",
        "0.0",
        "100.0",
        51.2,
        78.8,
        cut_color_map);

    bin_width = get_bin_width(h);
    h->SetXTitle("Production Vertex Z (cm)");
    h->SetYTitle(TString::Format("Events / %.3f cm", bin_width));
    c->Update();
    c->SaveAs("Production_Vertex_Z.pdf");
    delete c;

    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_mc_files,
        "z",
        "ProdVz",
        "100",
        "0.0",
        "100.0",
        51.2,
        78.8,
        cut_color_map);

    bin_width = get_bin_width(h);
    h->SetXTitle("Production Vertex Z (cm)");
    h->SetYTitle(TString::Format("Events / %.3f cm", bin_width));
    c->Update();
    c->SaveAs("Production_Vertex_Z_MC.pdf");
    delete c;

    // =================== Missing Mass^2 ===================
    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_data_files,
        "MM2",
        "RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5)",
        "100",
        "-0.1",
        " 0.1",
        -0.05,
        0.05,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("Missing Mass^{2} (GeV^{2})");
    h->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    c->Update();
    c->SaveAs("MM2.pdf");
    delete c;

    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_mc_files,
        "MM2",
        "RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5)",
        "100",
        "-0.1",
        " 0.1",
        -0.05,
        0.05,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("Missing Mass^{2} (GeV^{2})");
    h->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    c->Update();
    c->SaveAs("MM2_MC.pdf");
    delete c;

    // =================== Beam Energy ===================
    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_data_files,
        "eBeam",
        "EnPB",
        "100",
        "8.0",
        "9.0",
        8.2,
        8.8,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("Beam Energy (GeV)");
    h->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    c->Update();
    c->SaveAs("Beam_Energy.pdf");
    delete c;

    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_mc_files,
        "eBeam",
        "EnPB",
        "100",
        "8.0",
        "9.0",
        8.2,
        8.8,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("Beam Energy (GeV)");
    h->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    c->Update();
    c->SaveAs("Beam_Energy_MC.pdf");
    delete c;

    //=================== Chi2 / NDF ===================
    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_data_files,
        "chi2",
        "Chi2DOF",
        "40",
        "0.0",
        "8.0",
        0.0,
        5.0,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("#chi^{2} / NDF");
    h->SetYTitle(TString::Format("Events / %.3f", bin_width));
    c->Update();
    c->SaveAs("Chi2.pdf");
    delete c;

    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_mc_files,
        "chi2",
        "Chi2DOF",
        "40",
        "0.0",
        "8.0",
        0.0,
        5.0,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("#chi^{2} / NDF");
    h->SetYTitle(TString::Format("Events / %.3f", bin_width));
    c->Update();
    c->SaveAs("Chi2_MC.pdf");
    delete c;

    // =================== four momentum transfer -t ===================
    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_data_files,
        "t",
        "abs(-1*MASS2([proton],-GLUEXTARGET))",
        "100",
        "0.0",
        "2.5",
        0.0,
        1.0,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("-t (GeV^{2})");
    h->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    c->Update();
    c->SaveAs("t.pdf");
    delete c;

    c = setup_canvas(true);
    h = plot_variable(
        c,
        input_mc_files,
        "t",
        "abs(-1*MASS2([proton],-GLUEXTARGET))",
        "100",
        "0.0",
        "2.5",
        0.0,
        1.0,
        cut_color_map);
    bin_width = get_bin_width(h);
    h->SetXTitle("-t (GeV^{2})");
    h->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    c->Update();
    c->SaveAs("t_MC.pdf");
    delete c;

    // =================== Shower Quality  ===================
    // Create a larger canvas for a 2x2 grid, to see all 4 photons at once
    TCanvas *c_shower = new TCanvas("c_shower", "Shower Quality", 1200, 900);
    c_shower->Divide(1, 2, 0, 0);

    // Set up the legend pad at the top
    c_shower->cd(1);
    gPad->SetPad(0, 0.93, 0.9, 0.98);

    // Set up the main plotting area
    c_shower->cd(2);
    gPad->SetPad(0, 0, 1.0, 0.92);
    gPad->Divide(2, 2, 0.01, 0.01); // Divide into 2x2 grid with small margins

    TLegend *common_legend = new TLegend(0.05, 0.1, 0.95, 0.9);
    common_legend->SetNColumns(5);
    common_legend->SetBorderSize(0);
    common_legend->SetFillStyle(0);

    // Array of shower quality variables and their descriptions
    TString shower_vars[4] = {"ShQualityP4a", "ShQualityP4b", "ShQualityP5a", "ShQualityP5b"};
    TString photon_names[4] = {"Photon 1 (#pi^{0}_{1})", "Photon 2 (#pi^{0}_{1})",
                               "Photon 3 (#pi^{0}_{2})", "Photon 4 (#pi^{0}_{2})"};

    for (int i = 0; i < 4; i++)
    {
        c_shower->cd(2);
        gPad->cd(i + 1);
        gPad->SetLogy(1);

        TH1F *h_og = FSModeHistogram::getTH1F(
            input_data_files,
            NT,
            CATEGORY,
            shower_vars[i],
            "(100,0.0,1.0)",
            "CUT(rfSignal)");

        h_og->SetMinimum(1);
        h_og->SetLineColor(kGray);
        h_og->SetLineWidth(1);
        h_og->Draw("HIST");

        if (i == 0)
            common_legend->AddEntry(h_og, "Original Data", "l");

        // apply the special rf cut first
        TH1F *h_rf = FSModeHistogram::getTH1F(
            input_data_files,
            NT,
            CATEGORY,
            shower_vars[i],
            "(100,0.0,1.0)",
            "CUTWT(rf)");
        h_rf->SetMinimum(1);
        h_rf->SetLineColor(TColor::GetColor("#c849a9"));
        h_rf->SetLineWidth(1);
        h_rf->Draw("HIST SAME");

        if (i == 0)
            common_legend->AddEntry(h_rf, "RF", "l");

        // now iterate through our cuts and apply them one by one, except the one we are plotting
        TString cuts_so_far = "";
        TH1F *h_next = nullptr;
        for (map<TString, Int_t>::const_iterator it = cut_color_map.begin();
             it != cut_color_map.end(); ++it)
        {
            TString cut_name = it->first;
            if (cut_name == "shQuality")
                continue; // skip the cut we are plotting
            cuts_so_far += (cuts_so_far.IsNull() ? "" : ",") + cut_name;

            // apply this cut
            h_next = FSModeHistogram::getTH1F(
                input_data_files,
                NT,
                CATEGORY,
                shower_vars[i],
                "(100,0.0,1.0)",
                TString::Format("CUT(%s)*CUTWT(rf)", cuts_so_far.Data()));
            h_next->SetMinimum(1); // set minimum again after fixing zeros
            h_next->SetLineColor(it->second);
            h_next->SetLineWidth(1);

            h_next->Draw("HIST SAME");
            if (i == 0)
                common_legend->AddEntry(h_next, CUT_TO_LEGEND[cut_name], "l");
        }

        // Create a copy of the final histogram for the selection region fill
        TH1F *h_selection = (TH1F *)h_next->Clone("h_selection");
        h_selection->SetFillColorAlpha(kGreen + 1, 0.5);
        h_selection->SetFillStyle(1001);
        h_selection->SetLineWidth(0); // No line for the fill histogram

        // Zero out bins outside the selection region
        for (int j = 1; j <= h_selection->GetNbinsX(); ++j)
        {
            double bin_center = h_selection->GetBinCenter(j);
            if (bin_center < 0.5 || bin_center > 1.0)
            {
                h_selection->SetBinContent(j, 1); // small value for log scale
            }
        }

        // Draw the selection region fill
        h_selection->Draw("HIST SAME");
        if (i == 0)
            common_legend->AddEntry(h_selection, "Selection", "f");

        h_og->SetMaximum(h_og->GetMaximum() * 1.1); // add some headroom
        h_og->SetTitle("");
        bin_width = get_bin_width(h_og);
        h_og->SetXTitle(TString::Format("%s Quality", photon_names[i].Data()));
        h_og->SetYTitle(TString::Format("Events / %.3f", bin_width));

        // Adjust margins for the grid layout
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.1);
        gPad->SetRightMargin(0.05);
    }
    c_shower->cd(1);
    common_legend->Draw();

    c_shower->Update();
    c_shower->SaveAs("Shower_Quality.pdf");
    c_shower->Clear();

    // ==== MC Shower Quality ====
    delete c_shower;
    delete common_legend;
    c_shower = new TCanvas("c_shower", "Shower Quality", 1200, 900);
    c_shower->Divide(1, 2, 0, 0);

    // Set up the legend pad at the top
    c_shower->cd(1);
    gPad->SetPad(0, 0.93, 0.9, 0.98);

    // Set up the main plotting area
    c_shower->cd(2);
    gPad->SetPad(0, 0, 1.0, 0.92);
    gPad->Divide(2, 2, 0.01, 0.01); // Divide into 2x2 grid with small margins

    common_legend = new TLegend(0.05, 0.1, 0.95, 0.9);
    common_legend->SetNColumns(5);
    common_legend->SetBorderSize(0);
    common_legend->SetFillStyle(0);

    for (int i = 0; i < 4; i++)
    {
        c_shower->cd(2);
        gPad->cd(i + 1);
        gPad->SetLogy(1);

        TH1F *h_og = FSModeHistogram::getTH1F(
            input_mc_files,
            NT,
            CATEGORY,
            shower_vars[i],
            "(100,0.0,1.0)",
            "CUT(rfSignal)");

        h_og->SetMinimum(1);
        h_og->SetLineColor(kGray);
        h_og->SetLineWidth(1);
        h_og->Draw("HIST");

        if (i == 0)
            common_legend->AddEntry(h_og, "Original Data", "l");

        // apply the special rf cut first
        TH1F *h_rf = FSModeHistogram::getTH1F(
            input_mc_files,
            NT,
            CATEGORY,
            shower_vars[i],
            "(100,0.0,1.0)",
            "CUTWT(rf)");
        h_rf->SetMinimum(1);
        h_rf->SetLineColor(TColor::GetColor("#c849a9"));
        h_rf->SetLineWidth(1);
        h_rf->Draw("HIST SAME");

        if (i == 0)
            common_legend->AddEntry(h_rf, "RF", "l");

        // now iterate through our cuts and apply them one by one, except the one we are plotting
        TString cuts_so_far = "";
        TH1F *h_next = nullptr;
        for (map<TString, Int_t>::const_iterator it = cut_color_map.begin();
             it != cut_color_map.end(); ++it)
        {
            TString cut_name = it->first;
            if (cut_name == "shQuality")
                continue; // skip the cut we are plotting
            cuts_so_far += (cuts_so_far.IsNull() ? "" : ",") + cut_name;

            // apply this cut
            h_next = FSModeHistogram::getTH1F(
                input_mc_files,
                NT,
                CATEGORY,
                shower_vars[i],
                "(100,0.0,1.0)",
                TString::Format("CUT(%s)*CUTWT(rf)", cuts_so_far.Data()));
            h_next->SetMinimum(1); // set minimum again after fixing zeros
            h_next->SetLineColor(it->second);
            h_next->SetLineWidth(1);

            h_next->Draw("HIST SAME");
            if (i == 0)
                common_legend->AddEntry(h_next, CUT_TO_LEGEND[cut_name], "l");
        }

        // Create a copy of the final histogram for the selection region fill
        TH1F *h_selection = (TH1F *)h_next->Clone("h_selection");
        h_selection->SetFillColorAlpha(kGreen + 1, 0.5);
        h_selection->SetFillStyle(1001);
        h_selection->SetLineWidth(0); // No line for the fill histogram

        // Zero out bins outside the selection region
        for (int j = 1; j <= h_selection->GetNbinsX(); ++j)
        {
            double bin_center = h_selection->GetBinCenter(j);
            if (bin_center < 0.5 || bin_center > 1.0)
            {
                h_selection->SetBinContent(j, 1); // small value for log scale
            }
        }

        // Draw the selection region fill
        h_selection->Draw("HIST SAME");
        if (i == 0)
            common_legend->AddEntry(h_selection, "Selection", "f");

        h_og->SetMaximum(h_og->GetMaximum() * 1.1); // add some headroom
        h_og->SetTitle("");
        bin_width = get_bin_width(h_og);
        h_og->SetXTitle(TString::Format("%s Quality", photon_names[i].Data()));
        h_og->SetYTitle(TString::Format("Events / %.3f", bin_width));

        // Adjust margins for the grid layout
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.1);
        gPad->SetRightMargin(0.05);
    }
    c_shower->cd(1);
    common_legend->Draw();

    c_shower->Update();
    c_shower->SaveAs("Shower_Quality_MC.pdf");
    c_shower->Clear();

    if (dump_cache)
        FSHistogram::dumpHistogramCache();

    return;
}

/**
 * @brief Plot a variable with all other broad cuts applied
 *
 * This function creates a TCanvas containing a histogram of the specified variable,
 * applying the given broad cuts, except the one being plotted. If mc is true,
 * a Monte Carlo histogram is also drawn on top for comparison. This function must
 * be run AFTER load_broad_cuts() to ensure cuts are defined.
 *
 * @param[in] c pointer to TCanvas to draw on
 * @param[in] input_data_files data files to plot from
 * @param[in] cut_variable name of variable to exclude from cuts string
 * @param[in] tree_variable tree variable name in the ROOT file
 * @param[in] bins histogram bins
 * @param[in] lower_bound histogram lower bound
 * @param[in] upper_bound histogram upper bound
 * @param[in] cut_lower_bound lower boundary for the selection region visualization
 * @param[in] cut_upper_bound upper boundary for the selection region visualization
 * @param[in] cut_color_map map of cut names to colors
 * @return TH1F* pointer to the created histogram
 */
TH1F *plot_variable(
    TCanvas *c,
    const TString input_data_files,
    const TString cut_variable,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const double cut_lower_bound,
    const double cut_upper_bound,
    const std::map<TString, Int_t> &cut_color_map)
{
    // create histograms in format: file name, tree name, category, variable, bins and
    // bounds, then cuts. Our cuts are already loaded, so just their names are enough
    // for FSCut to find them

    c->cd(2); // draw operations in the bottom pad

    // first, plot the data without any cuts applied
    TH1F *h_og = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%s,%s,%s)", bins.Data(), lower_bound.Data(), upper_bound.Data()),
        "CUT(rfSignal)");

    bool isLogY = gPad->GetLogy();
    h_og->SetMinimum(isLogY ? 1 : 0);
    h_og->SetLineColor(kGray);
    h_og->SetLineWidth(1);
    h_og->Draw("HIST");

    TLegend *legend = new TLegend(0.15, 0.02, 1.0, 1.0);
    legend->SetNColumns(5);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(h_og, "Original Data", "l");

    // apply the special rf cut first
    TH1F *h_rf = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%s,%s,%s)", bins.Data(), lower_bound.Data(), upper_bound.Data()),
        "CUTWT(rf)");
    h_rf->SetMinimum(isLogY ? 1 : 0);
    h_rf->SetLineColor(TColor::GetColor("#c849a9"));
    h_rf->SetLineWidth(1);
    h_rf->Draw("HIST SAME");
    legend->AddEntry(h_rf, "RF Accidentals", "l");

    // now iterate through our cuts and apply them one by one, except the one we are plotting
    TString cuts_so_far = "";
    TH1F *h_next = nullptr;
    for (map<TString, Int_t>::const_iterator it = cut_color_map.begin();
         it != cut_color_map.end(); ++it)
    {
        TString cut_name = it->first;
        if (cut_name == cut_variable)
            continue; // skip the cut we are plotting
        cuts_so_far += (cuts_so_far.IsNull() ? "" : ",") + cut_name;

        // apply this cut
        h_next = FSModeHistogram::getTH1F(
            input_data_files,
            NT,
            CATEGORY,
            tree_variable,
            TString::Format("(%s,%s,%s)", bins.Data(), lower_bound.Data(), upper_bound.Data()),
            TString::Format("CUT(%s)*CUTWT(rf)", cuts_so_far.Data()));
        h_next->SetMinimum(isLogY ? 1 : 0); // set minimum again after fixing zeros
        h_next->SetLineColor(it->second);
        h_next->SetLineWidth(1);

        h_next->Draw("HIST SAME");
        legend->AddEntry(h_next, CUT_TO_LEGEND[cut_name], "l");
    }

    // Create a copy of the final histogram for the selection region fill
    TH1F *h_selection = (TH1F *)h_next->Clone("h_selection");
    h_selection->SetFillColorAlpha(kGreen + 1, 0.5);
    h_selection->SetFillStyle(1001);
    h_selection->SetLineWidth(0); // No line for the fill histogram

    // Zero out bins outside the selection region
    for (int i = 1; i <= h_selection->GetNbinsX(); ++i)
    {
        double bin_center = h_selection->GetBinCenter(i);
        if (bin_center < cut_lower_bound || bin_center > cut_upper_bound)
        {
            if (isLogY)
                h_selection->SetBinContent(i, 1); // small value for log scale
            else
                h_selection->SetBinContent(i, 0);
        }
    }

    // Draw the selection region fill
    h_selection->Draw("HIST SAME");
    legend->AddEntry(h_selection, "Selection", "f");

    // finalize some cosmetics
    h_og->SetMaximum(h_og->GetMaximum() * 1.1); // add some headroom
    h_og->SetTitle("");

    // finally, draw legend on first pad
    c->cd(1);
    legend->Draw();

    return h_og;
}

TCanvas *setup_canvas(bool logy = false)
{
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->Divide(1, 2, 0, 0);
    // Make top pad slim for legend, bottom for plot
    c->cd(1);
    gPad->SetPad(0, 0.93, 0.9, 0.98);
    gPad->SetTopMargin(0.05);
    c->cd(2);
    gPad->SetPad(0, 0, 0.98, 0.93);
    if (logy)
        gPad->SetLogy(1);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    c->SetRightMargin(0.02);
    c->SetTopMargin(0.02);
    return c;
}