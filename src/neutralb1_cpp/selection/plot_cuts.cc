/**
 * @file plot_broad_cuts.cc
 * @author Kevin Scheuer
 * @brief First step to select broad cuts in pi0pi0pi+pi- final state
 *
 * This file creates plots of several variables with broad cuts to select the physical
 * region of interest for a pi0pi0pi+pi- final state. Selection of omega signal
 * events is not yet done here, and so no combinatorics are studied yet.
 *
 * TODO: Have the main plotter function add the 2 permutation hists properly.
 * Then return the hist and the legend, that way the canvas can be set to log AFTER the 
 * hist is made and the minimum is set, to avoid dumb log issues.
 * 
 * reincorprate the direct MC comparison, by adding the sideband subtraction
 * to all the cuts. Have it be done last in the list of cuts (make sure to use * not &&)
 * This way the signal MC is directly comparable.
 */

#include <iostream>
#include <map>
#include <tuple>

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
#include "load_friend_cuts.cc"
#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"


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
void plot_accidentals(TString input_data_files, const TString NT, const TString CATEGORY);
std::pair<std::vector<TH1F*>, TLegend*> get_iterated_histograms(
    const TString input_data_files,
    const TString NT,
    const TString CATEGORY,
    const TString cut_variable,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const double cut_lower_bound,
    const double cut_upper_bound,
    const std::map<TString, Int_t> &cut_color_map
);
TH1F* get_og_histogram(
    const TString input_data_files,
    const TString NT,
    const TString CATEGORY,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound
);
TH1F* get_rf_histogram(
    const TString input_data_files,
    const TString NT,
    const TString CATEGORY,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound
);
std::vector<TH1F*> get_histograms_with_cuts(
    const TString input_data_files,
    const TString NT,
    const TString CATEGORY,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const TString cut_variable,
    const std::map<TString, Int_t> &cut_color_map,
    TLegend *legend
);
TH1F* get_selection_histogram(
    TH1F* h_source,
    const double cut_lower_bound,
    const double cut_upper_bound
);
TH1F* get_mc_histogram(
    const TString input_mc_files,
    const TString NT,
    const TString CATEGORY,
    const TString cut_variable,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const std::map<TString, Int_t> &cut_color_map,
    TH1F* h_data,
    const TString scale_choice,
    TLegend *legend
);

void plot_cuts(
    const int period,    
    bool read_cache = false,
    bool dump_cache = false)
{
    TString input_data_files;
    TString input_mc_files = "";

    TString input_files;
    
    // input_data_files = TString::Format(
    //     "/lustre24/expphy/volatile/halld/home/kscheuer/"
    //     "FSRoot-skimmed-trees/best-chi2/"
    //     "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_data.root",
    //     period);
    input_mc_files = TString::Format(
        "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/best-chi2/"
        "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_ver03.1_mc.root",
        period);
    input_data_files = input_mc_files; // TODO: for testing purposes
    
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(read_cache);

    // load our cuts
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();
    int n_permutations = load_friend_cuts();

    double bin_width;

    // =================== Accidental Subtraction ===================
    // Before we look at any other cuts let's look at the accidental subtraction itself
    // to see how well it performs. All other cuts will have it applied
    plot_accidentals(input_data_files, NT, CATEGORY);

    // create pointers for all upcoming plots
    TCanvas *c;
    std::vector<TH1F*> histograms;
    TH1F* mc_hist;
    TLegend *legend;
    // =================== Unused Shower Energy ===================
    c = setup_canvas(true);
    
    // TODO: all of this stuff below is almost the same, just not the x and y naming
    std::tie(histograms, legend) = get_iterated_histograms(
        input_data_files,
        NT,
        CATEGORY,
        "unusedE",
        "EnUnusedSh",
        "100",
        "0.0",
        "1.0",
        0.0,
        0.1,
        cut_color_map
    );
    c->cd(2);
    for(int i=0; i<histograms.size(); i++) {
        if (gPad->GetLogy()) {
            histograms[i]->SetMinimum(1);
        } 
        else {
            histograms[i]->SetMinimum(0);
        }
        if(i==0) {            
            histograms[i]->SetXTitle("Unused Shower Energy (GeV)");
            bin_width = get_bin_width(histograms[i]);
            histograms[i]->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
            histograms[i]->Draw("HIST");
        } else {
            histograms[i]->Draw("HIST SAME");
        }
    }
    mc_hist = get_mc_histogram(
        input_mc_files,
        NT,
        CATEGORY,
        "unusedE",
        "EnUnusedSh",
        "100",
        "0.0",
        "1.0",
        cut_color_map,
        histograms.back(),
        "max",
        legend
    );
    if (gPad->GetLogy()) {
        mc_hist->SetMinimum(1);
    } 
    else {
        mc_hist->SetMinimum(0);
    }
    mc_hist->Draw("E0 SAME");    

    c->cd(1);
    legend->Draw();        
    c->Update();
    c->SaveAs("Unused_Shower_Energy.pdf");
    delete c;

    // // =================== Number Unused Tracks ===================
    // c = setup_canvas(true);
    // h = plot_variable(
    //     c,
    //     input_data_files,
    //     "unusedTracks",
    //     "NumUnusedTracks",
    //     "6",
    //     "0.0",
    //     "6.0",
    //     0.0,
    //     1.0,
    //     cut_color_map,
    //     "max",
    //     input_mc_files);
    // bin_width = get_bin_width(h);
    // h->SetXTitle("Number of Unused Tracks");
    // h->SetYTitle("Events / Track");
    // c->Update();
    // c->SaveAs("Number_Unused_Tracks.pdf");
    // delete c;

    // // =================== Production Vertex ===================
    // c = setup_canvas(true);
    // h = plot_variable(
    //     c,
    //     input_data_files,
    //     "z",
    //     "ProdVz",
    //     "100",
    //     "0.0",
    //     "100.0",
    //     51.2,
    //     78.8,
    //     cut_color_map,
    //     "integral",
    //     input_mc_files);

    // bin_width = get_bin_width(h);
    // h->SetXTitle("Production Vertex Z (cm)");
    // h->SetYTitle(TString::Format("Events / %.3f cm", bin_width));
    // c->Update();
    // c->SaveAs("Production_Vertex_Z.pdf");
    // delete c;

    // // =================== Missing Mass^2 ===================
    // c = setup_canvas(true);
    // h = plot_variable(
    //     c,
    //     input_data_files,
    //     "MM2",
    //     "RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5)",
    //     "100",
    //     "-0.1",
    //     " 0.1",
    //     -0.05,
    //     0.05,
    //     cut_color_map,
    //     "max",
    //     input_mc_files);
    // bin_width = get_bin_width(h);
    // h->SetXTitle("Missing Mass^{2} (GeV^{2})");
    // h->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    // c->Update();
    // c->SaveAs("MM2.pdf");
    // delete c;

    // // =================== Beam Energy ===================
    // c = setup_canvas(true);
    // h = plot_variable(
    //     c,
    //     input_data_files,
    //     "eBeam",
    //     "EnPB",
    //     "100",
    //     "8.0",
    //     "9.0",
    //     8.2,
    //     8.8,
    //     cut_color_map,
    //     "integral",
    //     input_mc_files);
    // bin_width = get_bin_width(h);
    // h->SetXTitle("Beam Energy (GeV)");
    // h->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    // c->Update();
    // c->SaveAs("Beam_Energy.pdf");
    // delete c;

    // //=================== Chi2 / NDF ===================
    // c = setup_canvas(true);
    // h = plot_variable(
    //     c,
    //     input_data_files,
    //     "chi2",
    //     "Chi2DOF",
    //     "40",
    //     "0.0",
    //     "8.0",
    //     0.0,
    //     5.0,
    //     cut_color_map,
    //     "max",
    //     input_mc_files);
    // bin_width = get_bin_width(h);
    // h->SetXTitle("#chi^{2} / NDF");
    // h->SetYTitle(TString::Format("Events / %.3f", bin_width));
    // c->Update();
    // c->SaveAs("Chi2.pdf");
    // delete c;

    // // =================== four momentum transfer -t ===================
    // c = setup_canvas(true);
    // h = plot_variable(
    //     c,
    //     input_data_files,
    //     "t",
    //     "abs(-1*MASS2([proton],-GLUEXTARGET))",
    //     "100",
    //     "0.0",
    //     "2.5",
    //     0.0,
    //     1.0,
    //     cut_color_map,
    //     "max",
    //     input_mc_files);
    // bin_width = get_bin_width(h);
    // h->SetXTitle("-t (GeV^{2})");
    // h->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    // c->Update();
    // c->SaveAs("t.pdf");
    // delete c;

    // // =================== Missing Energy  ===================
    // c = setup_canvas(true);
    // h = plot_variable(
    //     c,
    //     input_data_files,
    //     "",
    //     "ENERGY(GLUEXTARGET, B) - ENERGY(1, 2, 3, 4, 5)",
    //     "100",
    //     "-3.0",
    //     "2.0",
    //     -3.0,
    //     3.0,
    //     cut_color_map);
    // bin_width = get_bin_width(h);
    // h->SetXTitle("Missing Energy (GeV)");
    // h->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    // c->Update();
    // c->SaveAs("Missing_Energy.pdf");
    // delete c;

    // // =================== Omega Mass ===================
    // TODO: add this

    // // =================== Omega Pi0 Mass ===================
    // c = setup_canvas(true);
    // h = plot_variable(
    //     c,
    //     input_data_files,
    //     "",
    //     "MASS(2,3,4,5)",
    //     "100",
    //     "1.0",
    //     "2.0",
    //     1.0,
    //     2.0,
    //     cut_color_map);
    // bin_width = get_bin_width(h);
    // h->SetXTitle("#omega#pi^{0} inv. mass (GeV)");
    // h->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    // c->Update();
    // c->SaveAs("Omega_Pi0_Mass.pdf");
    // delete c;

    // // =================== Proton Bachelor Mass ===================
    // TODO: fill in

    // // =================== pi01 momenta ===================
    // // TODO: once we figure out how to get these momenta properly, plot them

    // // =================== pi02 momenta ===================

    // // =================== Shower Quality  ===================
    // // Create a larger canvas for a 2x2 grid, to see all 4 photons at once
    // TCanvas *c_shower = new TCanvas("c_shower", "Shower Quality", 1200, 900);
    // c_shower->Divide(1, 2, 0, 0);

    // // Set up the legend pad at the top
    // c_shower->cd(1);
    // gPad->SetPad(0, 0.93, 0.9, 0.98);

    // // Set up the main plotting area
    // c_shower->cd(2);
    // gPad->SetPad(0, 0, 1.0, 0.92);
    // gPad->Divide(2, 2, 0.01, 0.01); // Divide into 2x2 grid with small margins

    // TLegend *common_legend = new TLegend(0.05, 0.1, 0.95, 0.9);
    // common_legend->SetNColumns(5);
    // common_legend->SetBorderSize(0);
    // common_legend->SetFillStyle(0);

    // // Array of shower quality variables and their descriptions
    // TString shower_vars[4] = {"ShQualityP4a", "ShQualityP4b", "ShQualityP5a", "ShQualityP5b"};
    // TString photon_names[4] = {"Photon 1 (#pi^{0}_{1})", "Photon 2 (#pi^{0}_{1})",
    //                            "Photon 3 (#pi^{0}_{2})", "Photon 4 (#pi^{0}_{2})"};

    // for (int i = 0; i < 4; i++)
    // {
    //     c_shower->cd(2);
    //     gPad->cd(i + 1);
    //     gPad->SetLogy(1);

    //     TH1F *h_og = FSModeHistogram::getTH1F(
    //         input_data_files,
    //         NT,
    //         CATEGORY,
    //         shower_vars[i],
    //         "(100,0.0,1.0)",
    //         "CUT(rfSignal)");

    //     h_og->SetMinimum(1);
    //     h_og->SetLineColor(kGray);
    //     h_og->SetLineWidth(1);
    //     h_og->Draw("HIST");

    //     if (i == 0)
    //         common_legend->AddEntry(h_og, "Original Data", "l");

    //     // apply the special rf cut first
    //     TH1F *h_rf = FSModeHistogram::getTH1F(
    //         input_data_files,
    //         NT,
    //         CATEGORY,
    //         shower_vars[i],
    //         "(100,0.0,1.0)",
    //         "CUTWT(rf)");
    //     h_rf->SetMinimum(1);
    //     h_rf->SetLineColor(TColor::GetColor("#c849a9"));
    //     h_rf->SetLineWidth(1);
    //     h_rf->Draw("HIST SAME");

    //     if (i == 0)
    //         common_legend->AddEntry(h_rf, "RF", "l");

    //     // now iterate through our cuts and apply them one by one, except the one we are plotting
    //     TString cuts_so_far = "";
    //     TH1F *h_next = nullptr;
    //     for (map<TString, Int_t>::const_iterator it = cut_color_map.begin();
    //          it != cut_color_map.end(); ++it)
    //     {
    //         TString cut_name = it->first;
    //         if (cut_name == "shQuality")
    //             continue; // skip the cut we are plotting
    //         cuts_so_far += (cuts_so_far.IsNull() ? "" : ",") + cut_name;

    //         // apply this cut
    //         h_next = FSModeHistogram::getTH1F(
    //             input_data_files,
    //             NT,
    //             CATEGORY,
    //             shower_vars[i],
    //             "(100,0.0,1.0)",
    //             TString::Format("CUT(%s)*CUTWT(rf)", cuts_so_far.Data()));
    //         h_next->SetMinimum(1); // set minimum again after fixing zeros
    //         h_next->SetLineColor(it->second);
    //         h_next->SetLineWidth(1);

    //         h_next->Draw("HIST SAME");
    //         if (i == 0)
    //             common_legend->AddEntry(h_next, CUT_TO_LEGEND[cut_name], "l");
    //     }

    //     // Create a copy of the final histogram for the selection region fill
    //     TH1F *h_selection = (TH1F *)h_next->Clone("h_selection");
    //     h_selection->SetFillColorAlpha(kGreen + 1, 0.5);
    //     h_selection->SetFillStyle(1001);
    //     h_selection->SetLineWidth(0); // No line for the fill histogram

    //     // Zero out bins outside the selection region
    //     for (int j = 1; j <= h_selection->GetNbinsX(); ++j)
    //     {
    //         double bin_center = h_selection->GetBinCenter(j);
    //         if (bin_center < 0.5 || bin_center > 1.0)
    //         {
    //             h_selection->SetBinContent(j, 1); // small value for log scale
    //         }
    //     }

    //     // Draw the selection region fill
    //     h_selection->Draw("HIST SAME");
    //     if (i == 0)
    //         common_legend->AddEntry(h_selection, "Selection", "f");

    //     if (input_mc_files != "")
    //     {
    //         TH1F *h_mc = FSModeHistogram::getTH1F(
    //             input_mc_files,
    //             NT,
    //             CATEGORY,
    //             shower_vars[i],
    //             "(100,0.0,1.0)",
    //             TString::Format("CUT(%s)*CUTWT(rf)", cuts_so_far.Data()));

    //         double scale;
    //         // scale MC to data using their maximum bin content
    //         scale = h_next->GetMaximum() / h_mc->GetMaximum();
    //         h_mc->Scale(scale);
    //         TString legend_scale = TString::Format("MC (max scaling=%.3f)", scale);
    //         if (i == 0)
    //             common_legend->AddEntry(h_mc, legend_scale, "p");

    //         // draw MC on top as blue points with error bars
    //         h_mc->SetMarkerStyle(24);
    //         h_mc->SetMarkerColor(kBlack);
    //         h_mc->SetLineColor(kBlack);
    //         h_mc->Draw("E0 SAME");
    //     }

    //     h_og->SetMaximum(h_og->GetMaximum() * 1.1); // add some headroom
    //     h_og->SetTitle("");
    //     bin_width = get_bin_width(h_og);
    //     h_og->SetXTitle(TString::Format("%s Quality", photon_names[i].Data()));
    //     h_og->SetYTitle(TString::Format("Events / %.3f", bin_width));

    //     // Adjust margins for the grid layout
    //     gPad->SetLeftMargin(0.15);
    //     gPad->SetBottomMargin(0.15);
    //     gPad->SetTopMargin(0.1);
    //     gPad->SetRightMargin(0.05);
    // }
    // c_shower->cd(1);
    // common_legend->Draw();

    // c_shower->Update();
    // c_shower->SaveAs("Shower_Quality.pdf");
    // c_shower->Clear();

    if (dump_cache)
        FSHistogram::dumpHistogramCache();

    return;
}

void plot_accidentals(TString input_data_files, const TString NT, const TString CATEGORY)
{
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
        "CUTSB(rf)");
    TH1F *h_acc_signal = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        "RFDeltaT",
        "(40,-2.0,2.0)",
        "CUT(rf)");
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
    double bin_width = get_bin_width(h_acc_data);
    h_acc_data->SetYTitle(TString::Format("Combos / %.3f ns", bin_width));
    c_rf->Update();
    c_rf->SaveAs("Accidental_Subtraction.pdf");    
}


/**
 * @brief Get histogram of one variable, with all other cuts iteravely applied
 *
 * @param[in] input_data_files data files to plot from
 * @param[in] NT name of the tree in the ROOT file
 * @param[in] CATEGORY category name in the ROOT file
 * @param[in] cut_variable name of variable to exclude from cuts string
 * @param[in] tree_variable tree variable name in the ROOT file
 * @param[in] bins histogram bins
 * @param[in] lower_bound histogram lower bound
 * @param[in] upper_bound histogram upper bound
 * @param[in] cut_lower_bound lower boundary for the selection region visualization
 * @param[in] cut_upper_bound upper boundary for the selection region visualization
 * @param[in] cut_color_map map of cut names to colors
 * @return 
 */
std::pair<std::vector<TH1F*>, TLegend*> get_iterated_histograms(
    const TString input_data_files,
    const TString NT,
    const TString CATEGORY,
    const TString cut_variable,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const double cut_lower_bound,
    const double cut_upper_bound,
    const std::map<TString, Int_t> &cut_color_map
)
{
    std::vector<TH1F*> histograms;

    // setup legend for all histograms and fill as we go
    TLegend *legend = new TLegend(0.15, 0.02, 1.0, 1.0);
    legend->SetNColumns(5);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    // first, plot the data without any cuts applied
    TH1F *h_og = get_og_histogram(
        input_data_files,
        NT,
        CATEGORY,
        tree_variable,
        bins,
        lower_bound,
        upper_bound
    );
    legend->AddEntry(h_og, "Original Data", "l");
    histograms.push_back(h_og);

    // apply the special rf cut as our first cut
    TH1F *h_rf = get_rf_histogram(
        input_data_files,
        NT,
        CATEGORY,
        tree_variable,
        bins,
        lower_bound,
        upper_bound
    );
    legend->AddEntry(h_rf, "RF Accidentals", "l");
    histograms.push_back(h_rf);

    // now iterate through our cuts and apply them one by one, excluding the one we 
    // are plotting. Note the legend is added to within the function
    std::vector<TH1F*> h_iter_cuts = get_histograms_with_cuts(
        input_data_files,
        NT,
        CATEGORY,
        tree_variable,
        bins,
        lower_bound,
        upper_bound,
        cut_variable,
        cut_color_map,
        legend
    );
    histograms.insert(
        histograms.end(),
        h_iter_cuts.begin(),
        h_iter_cuts.end()
    );
    
    // TODO: here is where sideband cut will be specially applied by adding the sideband
    // histogram to the signal histogram. It doesn't need the cut map since the friend
    // tree cuts is already setup to do that.
    
    // Create a copy of the final histogram for the selection region fill
    TH1F* h_selection = get_selection_histogram(
        h_iter_cuts.back(), // TODO: change to sideband result later
        cut_lower_bound,
        cut_upper_bound
    );
    legend->AddEntry(h_selection, "Selection", "f");
    histograms.push_back(h_selection);

    return std::make_pair(histograms, legend);
}

TCanvas *setup_canvas(bool logy = false)
{
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->Divide(1, 2, 0, 0);
    // Make top pad slim for legend, bottom for plot
    c->cd(1);
    gPad->SetPad(0, 0.93, 0.9, 0.95);
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

/**
 * @brief Get the histogram for data with no cuts applied
 * 
 * One cut is made, but it simply selects the prompt RF peak, where signal lies.
 * 
 * @param[in] input_data_files data files to plot from
 * @param[in] cut_variable name of variable to exclude from cuts string
 * @param[in] tree_variable tree variable name in the ROOT file
 * @param[in] bins histogram bins
 * @param[in] lower_bound histogram lower bound
 * @param[in] upper_bound histogram upper bound
 * @return TH1F* histogram with no cuts applied
 */
TH1F* get_og_histogram(
    const TString input_data_files,
    const TString NT,
    const TString CATEGORY,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound
)
{
    TH1F *h_og = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%s,%s,%s)", bins.Data(), lower_bound.Data(), upper_bound.Data()),
        "CUT(rfSignal)");

    h_og->SetLineColor(kGray);
    h_og->SetLineWidth(1);
    h_og->SetTitle(""); // this will be the first plotted histogram, so remove title

    return h_og;
}

/**
 * @brief Get the histogram for data with only the RF accidentals cut applied
 * 
 * @param[in] input_data_files data files to plot from
 * @param[in] NT tree name
 * @param[in] CATEGORY category name
 * @param[in] tree_variable tree variable name in the ROOT file
 * @param[in] bins histogram bins
 * @param[in] lower_bound histogram lower bound
 * @param[in] upper_bound histogram upper bound
 * @return TH1F* histogram with only RF accidentals cut applied
 */
TH1F* get_rf_histogram(
    const TString input_data_files,
    const TString NT,
    const TString CATEGORY,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound
)
{
    TH1F *h_rf = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%s,%s,%s)", bins.Data(), lower_bound.Data(), upper_bound.Data()),
        "CUTWT(rf)");
    h_rf->SetLineColor(TColor::GetColor("#c849a9"));
    h_rf->SetLineWidth(1);
    return h_rf;    
}


/**
 * @brief Get the set of histograms with each cut iteratively applied
 * 
 * @param[in] input_data_files data files to plot from
 * @param[in] NT tree name
 * @param[in] CATEGORY category name
 * @param[in] tree_variable tree variable name in the ROOT file
 * @param[in] bins histogram bins
 * @param[in] lower_bound histogram lower bound
 * @param[in] upper_bound histogram upper bound
 * @param[in] cut_variable cut variable name to exclude from cuts string
 * @param[in] cut_color_map map of cut names to colors
 * @param[in] legend legend object to add entries to
 * @return std::vector<TH1F*> vector of histograms with cuts applied
 */
std::vector<TH1F*> get_histograms_with_cuts(
    const TString input_data_files,
    const TString NT,
    const TString CATEGORY,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const TString cut_variable,
    const std::map<TString, Int_t> &cut_color_map,
    TLegend *legend
)
{
    std::vector<TH1F*> histograms;
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
        h_next->SetLineColor(it->second);
        h_next->SetLineWidth(1);        
        legend->AddEntry(h_next, CUT_TO_LEGEND[cut_name], "l");
        histograms.push_back(h_next);
    }

    return histograms;
}

TH1F* get_selection_histogram(
    TH1F* h_source,
    const double cut_lower_bound,
    const double cut_upper_bound)
{
    TH1F *h_selection = (TH1F *)h_source->Clone("h_selection");
    h_selection->SetFillColorAlpha(kGreen + 1, 0.5);
    h_selection->SetFillStyle(1001);
    h_selection->SetLineWidth(0); // No line for the fill histogram

    // Zero out bins outside the selection region
    for (int i = 1; i <= h_selection->GetNbinsX(); ++i)
    {
        double bin_center = h_selection->GetBinCenter(i);
        if (bin_center < cut_lower_bound || bin_center > cut_upper_bound)
        {
            h_selection->SetBinContent(i, 1); // small value for log scale
        }
    }

    return h_selection;
}

/**
 * @brief Get the signal mc histogram object
 * 
 * @param[in] input_mc_files files to get MC from
 * @param[in] NT tree name
 * @param[in] CATEGORY category name
 * @param[in] cut_variable cut variable name to exclude from cuts string
 * @param[in] tree_variable tree variable name in the ROOT file
 * @param[in] bins histogram bins
 * @param[in] lower_bound histogram lower bound
 * @param[in] upper_bound histogram upper bound
 * @param[in] cut_color_map map of cut names to colors
 * @param[in] h_data final data histogram to scale MC to
 * @param[in] scale_choice scaling method ("integral" or "max")
 * @param legend legend object to add entries to
 * @return TH1F* scaled MC histogram
 */
TH1F* get_mc_histogram(
    const TString input_mc_files,
    const TString NT,
    const TString CATEGORY,
    const TString cut_variable,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const std::map<TString, Int_t> &cut_color_map,
    TH1F* h_data,
    const TString scale_choice,
    TLegend *legend
)
{
    // remove the cut we are plotting from the cuts to apply to MC
    TString cuts = join_keys(cut_color_map);
    cuts.ReplaceAll(cut_variable, ""); // remove the cut we are plotting
    cuts.ReplaceAll(",,", ","); // clean up any double commas
    cuts = cuts.Strip(TString::kBoth, ','); // remove leading/trailing commas

    TH1F *h_mc = FSModeHistogram::getTH1F(
        input_mc_files,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%s,%s,%s)", bins.Data(), lower_bound.Data(), upper_bound.Data()),
        TString::Format("CUT(%s)", cuts.Data()));

    double scale;
    TString legend_scale;
    if (scale_choice == "integral")
    {
        // scale MC to data using their integrals
        scale = h_data->Integral() / h_mc->Integral();
        h_mc->Scale(scale);
        legend_scale = TString::Format("MC (integral scaling=%.3f)", scale);
    }
    else if (scale_choice == "max")
    {
        // scale MC to data using their maximum bin content
        scale = h_data->GetMaximum() / h_mc->GetMaximum();
        h_mc->Scale(scale);
        legend_scale = TString::Format("MC (max scaling=%.3f)", scale);
    }
    else
    {
        std::cerr << "Unknown scale choice: " << scale_choice << ". No scaling applied." << std::endl;
    }
    legend->AddEntry(h_mc, legend_scale, "p");

    // draw MC on top as black empty points with error bars
    h_mc->SetMarkerStyle(24);
    h_mc->SetMarkerColor(kBlack);
    h_mc->SetLineColor(kBlack);
    h_mc->Draw("E0 SAME");
    return h_mc;
}