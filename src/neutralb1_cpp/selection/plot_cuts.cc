/**
 * @file plot_broad_cuts.cc
 * @author Kevin Scheuer
 * @brief View effects of all cuts on tree variables
 *
 * This file creates plots of several variables to select the physical
 * region of interest for a pi0pi0pi+pi- final state.
 *
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

// CONSTANT MAPS
std::map<TString, TString> CUT_TO_LEGEND = {
    {"unusedE", "Unused Shower Energy"},
    {"unusedTracks", "Unused Charged Tracks"},
    {"z", "Production Vertex"},
    {"MM2", "Missing Mass^{2}"},
    {"eBeam", "Beam Energy"},
    {"chi2", "#chi^{2} / NDF"},
    {"t", "-t"},
    {"shQuality", "Shower Quality"}};
std::map<const int, TString> PERIOD_TO_LABEL = {
    {3, "S2017"},
    {4, "S2018"},
    {5, "F2018"}};

// Forward declarations
TCanvas *setup_canvas(bool logy = false);
void plot_accidentals(
    TString input_data_files, 
    const TString NT, 
    const TString CATEGORY, 
    int period);
void plot(
    TCanvas *c,
    const TString input_data_files,
    const TString input_mc_files,
    const TString NT,
    const TString CATEGORY,
    const TString cut_variable,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const double cut_lower_bound,
    const double cut_upper_bound,
    std::map<TString, Int_t> &cut_color_map,
    const TString x_title,
    const TString y_units,
    const TString scale_choice
);
std::vector<TH1F*> get_iterated_histograms(
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
    const std::map<TString, Int_t> &cut_color_map,
    TLegend *legend
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
    TLegend *legend,
    bool fill_legend = true
);
TH1F* get_selection_histogram(
    TH1F* h_source,
    const double cut_lower_bound,
    const double cut_upper_bound
);
TH1F* get_sideband_histogram(
    const TString input_data_files,
    const TString NT,    
    const TString CATEGORY,
    const TString cut_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound    
);
TH1F* get_mc_histogram(
    const TString input_mc_files,
    const TString NT,
    const TString CATEGORY,
    const TString cut_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
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
    
    input_data_files = TString::Format(
        "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/best-chi2/"
        "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_data.root",
        period);
    input_mc_files = TString::Format(
        "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/best-chi2/"
        "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_ver03.1_mc.root",
        period);    
    
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(read_cache);

    // load our broad cuts for the chi2 trees
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();
    double bin_width;

    // =================== Accidental Subtraction ===================
    // Before we look at any other cuts let's look at the accidental subtraction itself
    // to see how well it performs. All other cuts will have it applied
    plot_accidentals(input_data_files, NT, CATEGORY, period);

    // create pointers for all upcoming plots
    TCanvas *c;    
    // =================== Unused Shower Energy ===================
    c = setup_canvas(true);
    plot(
        c,
        input_data_files,
        input_mc_files,
        NT,
        CATEGORY,
        "unusedE",
        "EnUnusedSh",
        "100",
        "0.0",
        "1.0",
        0.0,
        0.1,
        cut_color_map,
        "Unused Shower Energy (GeV)",
        "GeV",
        "max"
    );
    c->SaveAs(TString::Format(
        "%s_Unused_Shower_Energy.pdf", 
        PERIOD_TO_LABEL.at(period).Data()));
    delete c;

    // =================== Number Unused Tracks ===================
    c = setup_canvas(true);
    plot(
        c,
        input_data_files,
        input_mc_files,
        NT,
        CATEGORY,
        "unusedTracks",
        "NumUnusedTracks",
        "6",
        "0.0",
        "6.0",
        0.0,
        1.0,
        cut_color_map,
        "Number of Unused Tracks",
        "Track",
        "max"
    );
    c->SaveAs(TString::Format(
        "%s_Number_Unused_Tracks.pdf", 
        PERIOD_TO_LABEL.at(period).Data()));
    delete c;

    // =================== Production Vertex ===================
    c = setup_canvas(true);
    plot(
        c,
        input_data_files,
        input_mc_files,
        NT,
        CATEGORY,
        "z",
        "ProdVz",
        "100",
        "0.0",
        "100.0",
        51.2,
        78.8,
        cut_color_map,
        "Production Vertex Z (cm)",
        "cm",
        "integral"
    );    
    c->SaveAs(TString::Format(
        "%s_Production_Vertex_Z.pdf", 
        PERIOD_TO_LABEL.at(period).Data()));
    delete c;

    // =================== Missing Mass^2 ===================
    c = setup_canvas(true);
    plot(
        c,
        input_data_files,
        input_mc_files,
        NT,
        CATEGORY,
        "MM2",
        "RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5)",
        "100",
        "-0.1",
        " 0.1",
        -0.05,
        0.05,
        cut_color_map,
        "Missing Mass^{2} (GeV^{2})",
        "GeV^{2}",
        "max"
    );
    c->SaveAs(TString::Format("%s_MM2.pdf", PERIOD_TO_LABEL.at(period).Data()));
    delete c;

    // =================== Beam Energy ===================
    c = setup_canvas(true);
    plot(
        c,
        input_data_files,
        input_mc_files,
        NT,
        CATEGORY,
        "eBeam",
        "EnPB",
        "100",
        "8.0",
        "9.0",
        8.2,
        8.8,
        cut_color_map,
        "Beam Energy (GeV)",
        "GeV",
        "integral"
    );
    c->SaveAs(TString::Format("%s_Beam_Energy.pdf", PERIOD_TO_LABEL.at(period).Data()));
    delete c;

    //=================== Chi2 / NDF ===================
    c = setup_canvas(true);
    plot(
        c,
        input_data_files,
        input_mc_files,
        NT,
        CATEGORY,
        "chi2",
        "Chi2DOF",
        "40",
        "0.0",
        "8.0",
        0.0,
        5.0,
        cut_color_map,
        "#chi^{2} / NDF",
        "",
        "max"
    );
    c->SaveAs(TString::Format("%s_Chi2.pdf", PERIOD_TO_LABEL.at(period).Data()));
    delete c;

    // =================== four momentum transfer -t ===================
    c = setup_canvas(true);
    plot(
        c,
        input_data_files,
        input_mc_files,
        NT,
        CATEGORY,
        "t",
        "abs(-1*MASS2([proton],-GLUEXTARGET))",
        "100",
        "0.0",
        "2.5",
        0.0,
        1.0,
        cut_color_map,
        "-t (GeV^{2})",
        "GeV^{2}",
        "max"
    );
    c->SaveAs(TString::Format("%s_t.pdf", PERIOD_TO_LABEL.at(period).Data()));
    delete c;

    // =================== Missing Energy  ===================
    c = setup_canvas(true);
    plot(
        c,
        input_data_files,
        input_mc_files,
        NT,
        CATEGORY,
        "missingEnergy",
        "ENERGY(GLUEXTARGET, B) - ENERGY(1, 2, 3, 4, 5)",
        "100",
        "-3.0",
        "3.0",
        -3.0,
        3.0,
        cut_color_map,
        "Missing Energy (GeV)",
        "GeV",
        "max"
    );
    c->SaveAs(TString::Format(
        "%s_Missing_Energy.pdf", 
        PERIOD_TO_LABEL.at(period).Data()));
    delete c;

    // Some of the next plots have an && in their tree name, which is specially handled
    // within the functions to plot both permutations
    // // =================== Omega Mass ===================
    c = setup_canvas(true);    
    plot(
        c,
        input_data_files,
        input_mc_files,
        NT,
        CATEGORY,
        "omegaMass",
        "MASS(2,3,4)&&MASS(2,3,5)",
        "100",
        "0.3",
        "1.8",
        0.6986,
        0.8666,
        cut_color_map,
        "#omega inv. mass (GeV)",
        "GeV",
        "max"
    );
    c->SaveAs(TString::Format(
        "%s_Omega_Mass.pdf", PERIOD_TO_LABEL.at(period).Data()));
    delete c;

    // =================== Omega Pi0 Mass ===================
    c = setup_canvas(true);
    plot(
        c,
        input_data_files,
        input_mc_files,
        NT,
        CATEGORY,
        "omegaPi0Mass",
        "MASS(2,3,4,5)",
        "100",
        "1.0",
        "2.0",
        1.0,
        2.0,
        cut_color_map,
        "#omega#pi^{0} inv. mass (GeV)",
        "GeV",
        "integral"
    );
    c->SaveAs(TString::Format(
        "%s_Omega_Pi0_Mass.pdf", PERIOD_TO_LABEL.at(period).Data()));
    delete c;

    // // =================== Proton Bachelor Mass ===================
    c = setup_canvas(true);
    plot(
        c,
        input_data_files,
        input_mc_files,
        NT,
        CATEGORY,
        "protonBachelorMass",
        "MASS(1,5)&&MASS(1,4)",
        "400",
        "1.0",
        "4.0",
        1.0,
        4.0,
        cut_color_map,
        "p'#pi^{0}_{bachelor} inv. mass (GeV)",
        "GeV",
        "integral"
    );
    c->SaveAs(TString::Format(
        "%s_Proton_Bachelor_Mass.pdf", PERIOD_TO_LABEL.at(period).Data()));
    delete c;

    // // =================== pi01 momenta ===================
    // // TODO: once we figure out how to get these momenta properly, plot them

    // // =================== pi02 momenta ===================

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
    TString friend_shower_vars[4] = {"shQuality_P3a", "shQuality_P3b", "shQuality_P2a", "shQuality_P2b"};
    TString photon_names[4] = {"Photon 1 (#pi^{0}_{1})", "Photon 2 (#pi^{0}_{1})",
                               "Photon 3 (#pi^{0}_{2})", "Photon 4 (#pi^{0}_{2})"};

    for (int i = 0; i < 4; i++)
    {
        c_shower->cd(2);
        gPad->cd(i + 1);
        gPad->SetLogy(1);

        TH1F *h_og = get_og_histogram(
            input_data_files,
            NT,
            CATEGORY,
            shower_vars[i],
            "100",
            "0.0",
            "1.0");

        if (i == 0)
            common_legend->AddEntry(h_og, "Original Data", "l");

        // apply the special rf cut first
        TH1F *h_rf = get_rf_histogram(
            input_data_files,
            NT,
            CATEGORY,
            shower_vars[i],
            "100",
            "0.0",
            "1.0"
        );
        if (i == 0)
            common_legend->AddEntry(h_rf, "RF", "l");

        // now iterate through our cuts and apply them one by one, except the one we are plotting
        bool fill_legend = (i == 0); // only fill legend on first plot
        std::vector<TH1F*> iterated_hists = get_histograms_with_cuts(
            input_data_files,
            NT,
            CATEGORY,
            shower_vars[i],
            "100",
            "0.0",
            "1.0",
            "shQuality",
            cut_color_map,
            common_legend,
            fill_legend
        );        

        // Create a copy of the final histogram for the selection region fill
        TH1F *h_selection = get_selection_histogram(
            iterated_hists.back(),
            0.5,
            1.0
        );
        if (i == 0)
            common_legend->AddEntry(h_selection, "Selection", "f");
        
        auto sh_map = cut_color_map;
        sh_map.erase("shQuality");
        TString sh_cuts = join_keys(sh_map);
        TH1F *h_mc = FSModeHistogram::getTH1F(
            input_mc_files,
            NT,
            CATEGORY,
            shower_vars[i],
            "(100,0.0,1.0)",
            TString::Format("CUT(%s)*CUTWT(rf)", sh_cuts.Data()));

        double scale;
        // scale MC to data using the maximum bin content
        scale = iterated_hists.back()->GetMaximum() / h_mc->GetMaximum();
        h_mc->Scale(scale);
        TString legend_scale = TString::Format("MC (max scaling=%.3f)", scale);
        if (i == 0)
            common_legend->AddEntry(h_mc, legend_scale, "p");

        // draw MC on top as blue points with error bars
        h_mc->SetMarkerStyle(24);
        h_mc->SetMarkerColor(kBlack);
        h_mc->SetLineColor(kBlack);        

        // set minimums for log scale
        h_og->SetMinimum(1); 
        h_rf->SetMinimum(1);        
        for (auto h_it : iterated_hists)
        {
            h_it->SetMinimum(1);
        }
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
        
        bin_width = get_bin_width(h_og);
        h_og->SetXTitle(TString::Format("%s Quality", photon_names[i].Data()));
        h_og->SetYTitle(TString::Format("Events / %.3f", bin_width));

        // draw all hists
        h_og->Draw("HIST");
        h_rf->Draw("HIST SAME");
        for (auto h_it : iterated_hists)
        {
            h_it->Draw("HIST SAME");
        }        
        h_selection->Draw("HIST SAME");
        h_mc->Draw("E0 SAME");

        // Adjust margins for the grid layout
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.1);
        gPad->SetRightMargin(0.05);
    }
    c_shower->cd(1);
    common_legend->Draw();

    c_shower->Update();
    c_shower->SaveAs(TString::Format("%s_Shower_Quality.pdf", PERIOD_TO_LABEL.at(period).Data()));
    c_shower->Clear();

    if (dump_cache)
        FSHistogram::dumpHistogramCache();

    return;
}
/**
 * @brief Plot the rf accidentals and their weights
 * 
 * @param input_data_files data files to plot from
 * @param NT name of the tree in the ROOT file
 * @param CATEGORY category name in the ROOT file
 * @param period data taking period
 */
void plot_accidentals(
    TString input_data_files, 
    const TString NT, 
    const TString CATEGORY, 
    int period)
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

    h_acc_data->SetTitle("");
    h_acc_data->SetXTitle("RF #DeltaT (ns)");
    double bin_width = get_bin_width(h_acc_data);
    h_acc_data->SetYTitle(TString::Format("Combos / %.3f ns", bin_width));
    c_rf->Update();
    c_rf->SaveAs(TString::Format(
        "%s_Accidental_Subtraction.pdf", 
        PERIOD_TO_LABEL.at(period).Data()));    
}


/**
 * @brief Plot the variable with iterative cuts applied
 * 
 * @param c TCanvas to plot on
 * @param input_data_files data files to plot from
 * @param input_mc_files Monte Carlo files to plot from
 * @param NT name of the tree in the ROOT file
 * @param CATEGORY category name in the ROOT file
 * @param cut_variable name of variable to exclude from cuts string
 * @param tree_variable tree variable name in the ROOT file
 * @param bins binning specification for the histogram
 * @param lower_bound lower bound of the histogram range
 * @param upper_bound upper bound of the histogram range
 * @param cut_lower_bound lower bound for the cut on the variable
 * @param cut_upper_bound upper bound for the cut on the variable
 * @param cut_color_map map of cut names to colors
 * @param x_title x-axis title for the histogram
 * @param y_units units for the y-axis
 * @param scale_choice MC scaling option for the histogram
 */
void plot(
    TCanvas *c,
    const TString input_data_files,
    const TString input_mc_files,
    const TString NT,
    const TString CATEGORY,
    const TString cut_variable,
    const TString tree_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    const double cut_lower_bound,
    const double cut_upper_bound,
    std::map<TString, Int_t> &cut_color_map,
    const TString x_title,
    const TString y_units,
    const TString scale_choice
)
{    
    std::vector<TH1F*> histograms; // will contain all data histograms to be plotted
    TH1F* mc_hist; // for the monte carlo histogram with all cuts applied
    double bin_width;

    TLegend *legend = new TLegend(0.15, 0.0, 0.93, 1.0);
    legend->SetNColumns(5);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    histograms = get_iterated_histograms(
        input_data_files,
        NT,
        CATEGORY,
        cut_variable,
        tree_variable,
        bins,
        lower_bound,
        upper_bound,
        cut_lower_bound,
        cut_upper_bound,
        cut_color_map,
        legend
    );

    c->cd(2); // histograms are on second pad
    for(int i=0; i<histograms.size(); i++) { // first plot data histograms with iterative cuts
        if (gPad->GetLogy()) {
            histograms[i]->SetMinimum(1);
        } 
        else {
            histograms[i]->SetMinimum(0);
        }
        if(i==0) { // first histogram (og) gets all the cosmetic changes
            histograms[i]->SetXTitle(x_title);
            bin_width = get_bin_width(histograms[i]);
            histograms[i]->SetYTitle(TString::Format("Events / %.3f %s", bin_width, y_units.Data()));
            histograms[i]->Draw("HIST");
        } else {
            histograms[i]->Draw("HIST SAME");
        }
    }
    // next plot the monte carlo hist
    mc_hist = get_mc_histogram(
        input_mc_files,
        NT,
        CATEGORY,
        cut_variable,
        bins,
        lower_bound,
        upper_bound,
        histograms.back(),
        scale_choice,
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
    return;
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
 * @param[in] cut_color_map map of cut names to colors\
 * @param[out] legend legend for all histograms
 * @return vector of histograms, each with a new cut applied
 */
std::vector<TH1F*> get_iterated_histograms(
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
    const std::map<TString, Int_t> &cut_color_map,
    TLegend *legend
)
{
    std::vector<TH1F*> histograms;

    // tree variables with "&&" are for cases which need permutations of the same variable
    // plotted together
    std::vector<TString> tree_variables_to_process;
    if (tree_variable.Contains("&&"))
    {
        Ssiz_t pos = 0;
        TString token;
        TString str = tree_variable;
        while ((pos = str.Index("&&")) != kNPOS)
        {
            token = str(0, pos);
            tree_variables_to_process.push_back(token);
            str.Remove(0, pos + 2);
        }
        tree_variables_to_process.push_back(str); // last token
    }
    else
    {
        tree_variables_to_process.push_back(tree_variable);
    }

    // ==== Get Data with no cuts applied ====    
    for (TString var : tree_variables_to_process)
    {
        TH1F *h_og_part = get_og_histogram(
            input_data_files,
            NT,
            CATEGORY,
            var,
            bins,
            lower_bound,
            upper_bound
        );
        if (histograms.size() == 0)
        {
            histograms.push_back(h_og_part);
        }
        else
        {
            histograms[0]->Add(h_og_part);
        }
    }  
    legend->AddEntry(histograms[0], "Original Data", "l");    

    // apply the special rf cut as our first cut
    for (TString var : tree_variables_to_process)
    {
        TH1F *h_rf_part = get_rf_histogram(
            input_data_files,
            NT,
            CATEGORY,
            var,
            bins,
            lower_bound,
            upper_bound
        );
        if (histograms.size() == 1)
        {
            histograms.push_back(h_rf_part);
        }
        else
        {
            histograms[1]->Add(h_rf_part);
        }
    }
    legend->AddEntry(histograms[1], "RF Accidentals", "l");

    // now iterate through our cuts and apply them one by one, excluding the one we 
    // are plotting. Note the legend is added to within the function
    for (size_t tree_iter = 0; tree_iter < tree_variables_to_process.size(); tree_iter++)    
    {   
        // this bool will avoid double filling of the legend for multiple permutations
        bool fill_legend = true ? tree_iter == 0 : false;        
        std::vector<TH1F*> h_iter_cuts_part = get_histograms_with_cuts(
            input_data_files,
            NT,
            CATEGORY,
            tree_variables_to_process[tree_iter],
            bins,
            lower_bound,
            upper_bound,
            cut_variable,
            cut_color_map,
            legend,
            fill_legend
        );
        // merge the parts together
        for (size_t i = 0; i < h_iter_cuts_part.size(); i++)
        {
            if (histograms.size() <= i + 2)
            {
                histograms.push_back(h_iter_cuts_part[i]);
            }
            else
            {
                histograms[i + 2]->Add(h_iter_cuts_part[i]);
            }
        }
    }

    // sideband subtract. These don't require multiple permutations since it is built
    // into the cut
    TH1F* h_sideband = get_sideband_histogram(
        input_data_files,
        NT,
        CATEGORY,
        cut_variable,        
        bins,
        lower_bound,
        upper_bound
    );
    legend->AddEntry(h_sideband, "Sideband Subtracted", "l");
    histograms.push_back(h_sideband);
    
    // Create a copy of the final histogram for the selection region fill
    TH1F* h_selection = get_selection_histogram(
        h_sideband,
        cut_lower_bound,
        cut_upper_bound
    );
    legend->AddEntry(h_selection, "Selection", "f");
    histograms.push_back(h_selection);

    return histograms;
}

TCanvas *setup_canvas(bool logy = false)
{
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->Divide(1, 2, 0, 0);
    // Make top pad slim for legend, bottom for plot
    c->cd(1);
    gPad->SetPad(0, 0.93, 0.98, 1.0);
    gPad->SetTopMargin(0.2);
    c->cd(2);
    gPad->SetPad(0, 0, 0.98, 0.93);
    if (logy)
        gPad->SetLogy(1);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.01);
    c->SetRightMargin(0.02);
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
    TLegend *legend,
    bool fill_legend = true
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
        if (fill_legend)   
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
 * @brief Get the sideband subtracted histogram for our cut variable
 * 
 * The friend trees require a lot of special handling. We've made them such that the 
 * cut variables are named the same as the branches. We also have cuts on those 
 * branches, which match the cut names for the broad cuts but with a "P" i.e. "chi2"
 * cut is "Pchi2" for the friend trees. This function will add the two friend tree
 * permutations and subtract their sidebands, to properly perform this sideband
 * subtraction.
 * 
 * @param input_files base data file name
 * @param NT tree name
 * @param CATEGORY category name
 * @param cut_variable variable we want to observe
 * @param bins histogram bins
 * @param lower_bound histogram lower bound
 * @param upper_bound histogram upper bound 
 * @return TH1F* sideband subtracted histogram
 */
TH1F* get_sideband_histogram(
    const TString input_files,
    const TString NT,    
    const TString CATEGORY,
    const TString cut_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound
)
{
    // load in the friend tree cuts, which uses our specialized branches we created
    // in each file. Remove the cut we want to observe
    std::map<TString, Int_t> friend_cut_map = load_friend_cuts();    
    friend_cut_map.erase(TString::Format("P%s", cut_variable.Data()));
    TString cuts = join_keys(friend_cut_map);

    // For those cuts with an "&&" in their tree variable name, we can specially handle
    // them here because the sideband subtraction collapses to one permutation only.
    TString new_cut_var = cut_variable;
    std::map<TString,TString> cut_mapping = {
        {"omegaMass", "MASS(3,4,5)"},
        {"omegaPi0Mass", "MASS(2,3,4,5)"},
        {"protonBachelorMass", "MRecoilPi"},
        {"missingEnergy", "ENERGY(GLUEXTARGET, B) - ENERGY(1, 2, 3, 4, 5)"}    
    };
    if(cut_mapping.find(cut_variable) != cut_mapping.end()) {
        new_cut_var = cut_mapping[cut_variable];
    }
    
    // now fill a histogram for each pi0 permutations signal and sidebands
    TH1F *h_signal[2];
    TH1F *h_sideband[2];

    for(int i=1; i<3; i++) {
        TString input_file = TString::Format("%s.permutation_%d", input_files.Data(), i);
        TString tree_name = TString::Format("%s_permutation_%d", NT.Data(), i);        
        h_signal[i-1] = FSModeHistogram::getTH1F(
            input_file,
            tree_name,
            CATEGORY,
            new_cut_var, // the friend tree branch names are conveniently named the same as the cut variable
            TString::Format("(%s,%s,%s)", bins.Data(), lower_bound.Data(), upper_bound.Data()),
            TString::Format("CUT(%s,signal)*CUTWT(Prf)", cuts.Data())
        );

        h_sideband[i-1] = FSModeHistogram::getTH1F(
            input_file,
            tree_name,
            CATEGORY,
            new_cut_var,
            TString::Format("(%s,%s,%s)", bins.Data(), lower_bound.Data(), upper_bound.Data()),
            TString::Format("CUT(%s,sideband)*CUTWT(Prf)", cuts.Data())
        );
    }
    h_signal[0]->Add(h_signal[1]); // add the two permutations together
    h_signal[0]->Add(h_sideband[0], -1); // subtract sideband permutation 1
    h_signal[0]->Add(h_sideband[1], -1); // subtract sideband permutation 2

    h_signal[0]->SetLineColor(TColor::GetColor("#bd1f01"));
    h_signal[0]->SetLineWidth(1);    
    
    return h_signal[0];
}

/**
 * @brief Get the signal mc histogram object
 * 
 * @param[in] input_mc_files files to get MC from
 * @param[in] NT tree name
 * @param[in] CATEGORY category name
 * @param[in] cut_variable cut variable name to exclude from cuts string 
 * @param[in] bins histogram bins
 * @param[in] lower_bound histogram lower bound
 * @param[in] upper_bound histogram upper bound 
 * @param[in] h_data final data histogram to scale MC to
 * @param[in] scale_choice scaling method ("integral" or "max") 
 * @return TH1F* scaled MC histogram
 */
TH1F* get_mc_histogram(
    const TString input_mc_files,
    const TString NT,
    const TString CATEGORY,
    const TString cut_variable,
    const TString bins,
    const TString lower_bound,
    const TString upper_bound,
    TH1F* h_data,
    const TString scale_choice,
    TLegend *legend
)
{
    TH1F *h_mc = get_sideband_histogram(
        input_mc_files,
        NT,
        CATEGORY,
        cut_variable,
        bins,
        lower_bound,
        upper_bound
    );

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
    return h_mc;
}