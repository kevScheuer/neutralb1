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

// Forward declarations
TString join_keys(const std::map<TString, Int_t> &m, const TString &delimiter = ",");
TH1F *plot_variable(
    const TString input_data_files,
    const TString cut_variable,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const double cut_lower_bound,
    const double cut_upper_bound,
    const std::map<TString, Int_t> &cut_color_map,
    const TString input_mc_files,
    bool bggen);

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

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    double bin_width;

    // =================== Accidental Subtraction ===================
    // Before we look at any other cuts let's look at the accidental subtraction itself
    // to see how well it performs. All other cuts will have it applied
    c->cd();
    gPad->SetLogy(0);
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
        "CUTSBWT(rf)*(-1.0)");
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
    h_acc_signal->SetFillColor(kGreen+1);
    h_acc_signal->SetFillStyle(3001);

    h_acc_subtracted->SetLineWidth(0);
    h_acc_subtracted->SetFillColor(kRed+2);
    h_acc_subtracted->SetFillStyle(3002);

    h_acc_data->Draw("HIST");
    h_acc_subtracted->Draw("HIST SAME");
    h_acc_signal->Draw("HIST SAME");        

    TLegend *legend_acc = new TLegend(0.65, 0.65, 0.88, 0.88);
    legend_acc->AddEntry(h_acc_data, "Original Data", "l");
    legend_acc->AddEntry(h_acc_signal, "In-Time Signal Region", "f");
    legend_acc->AddEntry(h_acc_subtracted, "Weighted Out-of-Time Combos", "f");
    legend_acc->Draw();

    h_acc_data->SetXTitle("RF #DeltaT (ns)");
    h_acc_data->SetMinimum(h_acc_subtracted->GetMinimum()*1.1);
    bin_width = get_bin_width(h_acc_data);
    h_acc_data->SetYTitle(TString::Format("Combos / %.3f ns", bin_width));
    c->Update();
    c->SaveAs("Accidental_Subtraction.pdf");
    c->Clear();
    
    // =================== Unused Shower Energy ===================
    c->cd();
    gPad->SetLogy(1);
    TH1F *h;
    h = plot_variable(
        input_data_files,
        "unusedE",
        "EnUnusedSh",
        "100",
        "0.0",
        "1.0",
        0.0,
        0.1,
        cut_color_map,
        input_mc_files,
        bggen);
    bin_width = get_bin_width(h);
    h->SetMinimum(1);
    h->SetXTitle("Unused Shower Energy (GeV)");
    h->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    c->Update();
    c->SaveAs("Unused_Shower_Energy.pdf");
    c->Clear();

    // =================== Number Unused Tracks ===================
    c->cd();
    gPad->SetLogy(0);
    h = plot_variable(
        input_data_files,
        "unusedTracks",
        "NumUnusedTracks",
        "10",
        "0.0",
        "10.0",
        0.0,
        1.0,
        cut_color_map,
        input_mc_files,
        bggen);
    bin_width = get_bin_width(h);
    h->SetXTitle("Number of Unused Tracks");
    h->SetYTitle("Events / Track");
    c->Update();
    c->SaveAs("Number_Unused_Tracks.pdf");
    c->Clear();

    // =================== Production Vertex ===================
    c->cd();
    gPad->SetLogy(0);
    h = plot_variable(
        input_data_files,
        "z",
        "ProdVz",
        "100",
        "0.0",
        "100.0",
        51.2,
        78.8,
        cut_color_map,
        input_mc_files,
        bggen);

    bin_width = get_bin_width(h);
    h->SetXTitle("Production Vertex Z (cm)");
    h->SetYTitle(TString::Format("Events / %.3f cm", bin_width));
    c->Update();
    c->SaveAs("Production_Vertex_Z.pdf");
    c->Clear();

    // =================== Missing Mass^2 ===================
    c->cd();
    gPad->SetLogy(1);
    h = plot_variable(
        input_data_files,
        "MM2",
        "RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5)",
        "100",
        "-0.1",
        " 0.1",
        -0.05,
        0.05,
        cut_color_map,
        input_mc_files,
        bggen);
    bin_width = get_bin_width(h);
    h->SetMinimum(1);
    h->SetXTitle("Missing Mass^{2} (GeV^{2})");
    h->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    c->Update();
    c->SaveAs("MM2.pdf");
    c->Clear();

    // =================== Beam Energy ===================
    c->cd();
    gPad->SetLogy(0);
    h = plot_variable(
        input_data_files,
        "eBeam",
        "EnPB",
        "100",
        "8.0",
        "9.0",
        8.2,
        8.8,
        cut_color_map,
        input_mc_files,
        bggen);
    bin_width = get_bin_width(h);
    h->SetXTitle("Beam Energy (GeV)");
    h->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    c->Update();
    c->SaveAs("Beam_Energy.pdf");
    c->Clear();

    //=================== Chi2 / NDF ===================
    c->cd();
    gPad->SetLogy(0);
    h = plot_variable(
        input_data_files,
        "chi2",
        "Chi2DOF",
        "40",
        "0.0",
        "8.0",
        0.0,
        5.0,
        cut_color_map,
        input_mc_files,
        bggen);
    bin_width = get_bin_width(h);
    h->SetXTitle("#chi^{2} / NDF");
    h->SetYTitle(TString::Format("Events / %.3f", bin_width));
    c->Update();
    c->SaveAs("Chi2.pdf");
    c->Clear();

    // =================== four momentum transfer -t ===================
    c->cd();
    gPad->SetLogy(0);
    h = plot_variable(
        input_data_files,
        "t",
        "abs(-1*MASS2([proton],-GLUEXTARGET))",
        "100",
        "0.0",
        "2.5",
        0.0,
        1.0,
        cut_color_map,
        input_mc_files,
        bggen);
    bin_width = get_bin_width(h);
    h->SetXTitle("-t (GeV^{2})");
    h->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    c->Update();
    c->SaveAs("t.pdf");
    c->Clear();

    // TODO: Consider adding the recoil proton momentum cut
    // https://halldweb.jlab.org/DocDB/0049/004924/002/TaskForce_ProtonFiducialVertex.pdf

    // // =================== Shower Quality  ===================
    // // Create a larger canvas for a 2x2 grid, to see all 4 photons at once
    // TCanvas *c_shower = new TCanvas("c_shower", "Shower Quality", 1200, 900);
    // c_shower->Divide(2, 2);

    // // Array of shower quality variables and their descriptions
    // TString shower_vars[4] = {"ShQualityP4a", "ShQualityP4b", "ShQualityP5a", "ShQualityP5b"};
    // TString photon_names[4] = {"Photon 1 (#pi^{0}_{1})", "Photon 2 (#pi^{0}_{1})",
    //                            "Photon 3 (#pi^{0}_{2})", "Photon 4 (#pi^{0}_{2})"};

    // for (int i = 0; i < 4; i++)
    // {
    //     c_shower->cd(i + 1);
    //     gPad->SetLogy();

    //     // For shower quality, we exclude the shQuality cut to see the full range
    //     std::map<TString, Int_t> filtered_cuts = cut_color_map;
    //     filtered_cuts.erase("shQuality");
    //     // Join all keys of filtered_cuts into a comma separated TString
    //     TString cut_string = join_keys(filtered_cuts);

    //     // Create histogram for this shower quality variable
    //     TH1F *h_sq_data = FSModeHistogram::getTH1F(
    //         input_data_files,
    //         NT,
    //         CATEGORY,
    //         shower_vars[i],
    //         "(100,0.0,1.0)",
    //         TString::Format("CUT(%s)", cut_string.Data()));
    //     h_sq_data->SetLineColor(kBlack);
    //     h_sq_data->SetLineWidth(2);
    //     h_sq_data->SetMinimum(1); // avoid log(0)
    //     h_sq_data->Draw("HIST");

    //     // Create a copy for the selection region fill
    //     TH1F *h_sq_selection = (TH1F *)h_sq_data->Clone(TString::Format("h_sq_selection_%d", i));
    //     h_sq_selection->SetFillColor(kGreen);
    //     h_sq_selection->SetFillStyle(3002);
    //     h_sq_selection->SetLineWidth(0);

    //     // Zero out bins outside the selection region (>0.5)
    //     for (int j = 1; j <= h_sq_selection->GetNbinsX(); ++j)
    //     {
    //         double bin_center = h_sq_selection->GetBinCenter(j);
    //         if (bin_center < 0.5)
    //         {
    //             h_sq_selection->SetBinContent(j, 1); // small value for log scale
    //         }
    //     }

    //     // Draw the selection region fill
    //     h_sq_selection->Draw("HIST SAME");

    //     // Redraw the data histogram on top
    //     h_sq_data->Draw("HIST SAME");

    //     double max_val = h_sq_data->GetMaximum();

    //     if (input_mc_files != "")
    //     {
    //         TH1F *h_sq_mc = FSModeHistogram::getTH1F(
    //             input_mc_files,
    //             NT,
    //             CATEGORY,
    //             shower_vars[i],
    //             "(50,0.0,1.0)",
    //             TString::Format("CUT(%s)", filtered_cuts.Data()));
    //         // scale MC to data using their integrals
    //         double scale = h_sq_data->Integral() / h_sq_mc->Integral();
    //         h_sq_mc->Scale(scale);

    //         TLegend *legend = new TLegend(0.15, 0.15, 0.38, 0.38);
    //         legend->AddEntry(h_sq_data, "Data", "l");
    //         legend->AddEntry(h_sq_mc, TString::Format("MC (scale=%.3f)", scale), "p");
    //         legend->AddEntry(h_sq_selection, "Selection", "f");
    //         legend->Draw();

    //         // draw MC on top as blue points with error bars
    //         h_sq_mc->SetMarkerStyle(24);
    //         h_sq_mc->SetMarkerColor(kBlue);
    //         h_sq_mc->SetLineColor(kBlue);
    //         h_sq_mc->Draw("E1 SAME");

    //         max_val = std::max(h_sq_data->GetMaximum(), h_sq_mc->GetMaximum());
    //     }
    //     else
    //     {
    //         // Create legend even when no MC data
    //         TLegend *legend = new TLegend(0.15, 0.15, 0.38, 0.38);
    //         legend->AddEntry(h_sq_data, "Data", "l");
    //         legend->AddEntry(h_sq_selection, "Selection", "f");
    //         legend->Draw();
    //     }

    //     h_sq_data->SetMaximum(max_val * 1.1); // add some headroom

    //     bin_width = get_bin_width(h_sq_data);
    //     h_sq_data->SetXTitle(TString::Format("%s Quality", photon_names[i].Data()));
    //     h_sq_data->SetYTitle(TString::Format("Combos / %.3f", bin_width));

    //     h_sq_data->SetMinimum(1);

    //     // Adjust margins for the grid layout
    //     gPad->SetLeftMargin(0.15);
    //     gPad->SetBottomMargin(0.15);
    //     gPad->SetTopMargin(0.1);
    //     gPad->SetRightMargin(0.05);
    // }
    // c_shower->Update();
    // c_shower->SaveAs("Shower_Quality.pdf");
    // c_shower->Clear();

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
 * @param[in] input_data_files data files to plot from
 * @param[in] cut_variable name of variable to exclude from cuts string
 * @param[in] tree_variable tree variable name in the ROOT file
 * @param[in] bins histogram bins
 * @param[in] lower_bound histogram lower bound
 * @param[in] upper_bound histogram upper bound
 * @param[in] cut_lower_bound lower boundary for the selection region visualization
 * @param[in] cut_upper_bound upper boundary for the selection region visualization
 * @param[in] cut_color_map map of cut names to colors
 * @param[in] input_mc_files Monte Carlo data files
 * @param[in] bggen whether to include background generation
 * @return TH1F* pointer to the created histogram
 */
TH1F *plot_variable(
    const TString input_data_files,
    const TString cut_variable,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const double cut_lower_bound,
    const double cut_upper_bound,
    const std::map<TString, Int_t> &cut_color_map,
    const TString input_mc_files = "",
    bool bggen = false)
{    
    // create histograms in format: file name, tree name, category, variable, bins and
    // bounds, then cuts. Our cuts are already loaded, so just their names are enough 
    // for FSCut to find them

    // first, plot the data without any cuts applied
    TH1F* h_og = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%s,%s,%s)", bins.Data(), lower_bound.Data(), upper_bound.Data()),
        "");
    
    bool isLogY = gPad->GetLogy();
    h_og->SetMinimum(isLogY ? 1 : 0); 
    h_og->SetLineColor(kGray);
    h_og->SetLineWidth(1);
    h_og->Draw("HIST");

    // TODO: Legend will be huge. Probably best to create a thin canvas above the plot
    // that has a long legend only. Need to pass canvas to function.
    // TLegend *legend = new TLegend(l_x1, l_y1, l_x2, l_y2);    
    // legend->AddEntry(h_og, "Original Data", "l");
    // legend->AddEntry(h_mc, TString::Format("MC (scale=%.3f)", scale), "p");
    // legend->AddEntry(h_selection, "Selection", "f");

    // apply the special rf cut first
    TH1F* h_rf = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%s,%s,%s)", bins.Data(), lower_bound.Data(), upper_bound.Data()),
        "CUTWT(rf)");
    h_rf->SetMinimum(isLogY ? 1 : 0);
    h_rf->SetLineColor(kRed+2);
    h_rf->SetLineWidth(1);
    h_rf->Draw("HIST SAME");

    // now iterate through our cuts and apply them one by one, except the one we are plotting    
    TString cuts_so_far = "";
    TH1F* h_next = nullptr;
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
            TString::Format("CUT(%s)", cuts_so_far.Data())
            );
        h_next->SetMinimum(isLogY ? 1 : 0); // set minimum again after fixing zeros
        h_next->SetLineColor(it->second);
        h_next->SetLineWidth(1);     

        h_next->Draw("HIST SAME");      
        // TODO: add to legend
    }

    // Create a copy of the final histogram for the selection region fill
    TH1F *h_selection = (TH1F *)h_next->Clone("h_selection");
    h_selection->SetFillColor(kGreen+1);
    h_selection->SetFillStyle(3001);
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

    // for the legend, upper right hand corner works well EXCEPT for production vertex
    float l_x1, l_x2, l_y1, l_y2;
    if (cut_variable == "z")
    {
        // move legend to upper left corner
        l_x1 = 0.15;
        l_x2 = 0.38;
        l_y1 = 0.65;
        l_y2 = 0.88;
    }
    else
    {
        // default upper right corner
        l_x1 = 0.65;
        l_x2 = 0.88;
        l_y1 = 0.65;
        l_y2 = 0.88;
    }

    // if (input_mc_files != "")
    // {
    //     TH1F *h_mc = FSModeHistogram::getTH1F(
    //         input_mc_files,
    //         NT,
    //         CATEGORY,
    //         tree_variable,
    //         TString::Format("(%s,%s,%s)", bins.Data(), lower_bound.Data(), upper_bound.Data()),
    //         TString::Format("CUT(%s)", filtered_cuts.Data()));
    //     // scale MC to data using their integrals
    //     double scale = h_data->Integral() / h_mc->Integral();
    //     h_mc->Scale(scale);

        
    //     legend->Draw();

    //     // draw MC on top as blue points with error bars
    //     h_mc->SetMarkerStyle(24);
    //     h_mc->SetMarkerColor(kBlue);
    //     h_mc->SetLineColor(kBlue);
    //     h_mc->Draw("E1 SAME");

    //     max_val = std::max(h_data->GetMaximum(), h_mc->GetMaximum());
    // }
    // else
    // {
    //     // Create legend even when no MC data
    //     TLegend *legend = new TLegend(l_x1, l_y1, l_x2, l_y2);
    //     legend->AddEntry(h_data, "Data", "l");
    //     legend->AddEntry(h_selection, "Selection", "f");
    //     legend->Draw();
    // }

    h_og->SetMaximum(h_og->GetMaximum() * 1.1); // add some headroom

    return h_og;
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