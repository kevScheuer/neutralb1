/**
 * @file plot_broad_cuts.cc
 * @author Kevin Scheuer
 * @brief View effects of all cuts on tree variables
 *
 * NOTE: that "og" histograms use the particle ordering of FSRoot 
 *    0       1          2             3              4              5
 * [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (omega)] [pi0 (bachelor)]
 * 
 * while the "cut" and "mc" histograms use the AmpTools ordering 
 *    0       1            2               3            4             5
 * [beam] [proton] [pi0 (bachelor)] [pi0 (omega)] [pi+ (omega)] [pi- (omega)]
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

#include "load_broad_branch_cuts.cc"
#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"

// Forward declarations
TString remove_key_from_cuts(
    std::map<TString, Int_t> &cut_color_map, 
    const TString key_to_remove);
TH1F* get_og_histogram(
    const TString input_data_files,
    const TString NT,
    const TString CATEGORY,
    const TString tree_variable,
    const int bins,
    const double lower_bound,
    const double upper_bound
);
TH1F* get_cut_histogram(
    const TString input_data_signal,
    const TString input_data_sideband,
    const TString NT,
    const TString CATEGORY,
    const TString tree_variable,
    const TString cuts,
    const int bins,
    const double lower_bound,
    const double upper_bound
);
TH1F* get_selection_histogram(
    TH1F* h_source,
    const double cut_lower_bound,
    const double cut_upper_bound
);
TH1F* get_mc_histogram(
    TH1F* h_data,
    const TString input_mc_signal,
    const TString input_mc_sideband,
    const TString NT,
    const TString CATEGORY,
    const TString tree_variable,
    const TString cuts,
    const int bins,
    const double lower_bound,
    const double upper_bound,
    const TString scale_choice,
    TLegend *legend
);

void plot_cuts(    
    bool log_scale = false,
    bool read_cache = false,
    bool dump_cache = false
)
{    
    TString input_data_og = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/best-chi2/"
        "tree_pi0pi0pippim__B4_bestChi2_SKIM_allPeriods_data_rfPrompt.root";

    TString input_data_signal = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/reordered/"
        "tree_pi0pi0pippim__B4_reordered_SKIM_allPeriods_data_signal.root";
    TString input_data_sideband = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/reordered/"
        "tree_pi0pi0pippim__B4_reordered_SKIM_allPeriods_data_sideband.root";

    TString input_mc_signal = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/reordered/"
        "tree_pi0pi0pippim__B4_reordered_SKIM_allPeriods_ver03.1_mc_signal.root";
    TString input_mc_sideband = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/reordered/"
        "tree_pi0pi0pippim__B4_reordered_SKIM_allPeriods_ver03.1_mc_sideband.root";

    
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(read_cache);

    // load our broad cuts for the chi2 trees
    std::map<TString, Int_t> cut_color_map = load_broad_branch_cuts();    

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    if (log_scale)
        c->SetLogy();
    TLegend *legend = new TLegend(0.63, 0.63, 0.88, 0.88);
    TString cuts;
    TH1F *h_og, *h_cut, *h_selection, *h_mc;
    double bin_width;
    legend->SetBorderSize(1);    

    // NOTE: that for all original data files, the tree variable is slightly different
    // as it uses the original branch name. We've renamed the branches to match the
    // easier to use cut definitions in the signal and sideband files
    // =================== Unused Shower Energy ===================
    
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "EnUnusedSh",
        100,
        0.0,
        1.0
    );
    h_og->SetXTitle("Unused Shower Energy (GeV)");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    legend->AddEntry(h_og, "Original Data", "l");

    cuts = remove_key_from_cuts(
        cut_color_map,
        "unusedE"
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "unusedE",
        cuts,
        100,
        0.0,
        1.0
    );
    legend->AddEntry(h_cut, "After Cuts", "l");

    h_selection = get_selection_histogram(
        h_cut,
        0.0,
        0.1
    );
    legend->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "unusedE",
        cuts,
        100,
        0.0,
        1.0,
        "max",
        legend
    );
    
    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    legend->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_Unused_Shower_Energy%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    legend->Clear();
    

    // =================== Number Unused Tracks ===================
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "NumUnusedTracks",
        6,
        0.0,
        6.0
    );
    h_og->SetXTitle("Number of Unused Tracks");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.1f Track", bin_width));
    legend->AddEntry(h_og, "Original Data", "l");

    cuts = remove_key_from_cuts(
        cut_color_map,
        "unusedTracks"
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "unusedTracks",
        cuts,
        6,
        0.0,
        6.0
    );
    legend->AddEntry(h_cut, "After Cuts", "l");

    h_selection = get_selection_histogram(
        h_cut,
        0.0,
        1.0
    );
    legend->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "unusedTracks",
        cuts,
        6,
        0.0,
        6.0,
        "max",
        legend
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    legend->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_Number_Unused_Tracks%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    legend->Clear();

    // =================== Production Vertex ===================
    TLegend *leg_z = new TLegend(0.2, 0.63, 0.45, 0.88); // custom legend location
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "ProdVz",
        100,
        0.0,
        100.0
    );
    h_og->SetXTitle("Production Vertex Z (cm)");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f cm", bin_width));
    leg_z->AddEntry(h_og, "Original Data", "l");

    cuts = remove_key_from_cuts(
        cut_color_map,
        "z"
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "z",
        cuts,
        100,
        0.0,
        100.0
    );
    leg_z->AddEntry(h_cut, "After Cuts", "l");

    h_selection = get_selection_histogram(
        h_cut,
        51.2,
        78.8
    );
    leg_z->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "z",
        cuts,
        100,
        0.0,
        100.0,
        "integral",
        leg_z
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    leg_z->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_Production_Vertex_Z%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    delete leg_z;

    // =================== Missing Mass^2 ===================
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5)",
        100,
        -0.1,
        0.1
    );
    h_og->SetXTitle("Missing Mass^{2} (GeV^{2})");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    legend->AddEntry(h_og, "Original Data", "l");

    cuts = remove_key_from_cuts(
        cut_color_map,
        "MM2"
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "MM2",
        cuts,
        100,
        -0.1,
        0.1
    );
    legend->AddEntry(h_cut, "After Cuts", "l");

    h_selection = get_selection_histogram(
        h_cut,
        -0.05,
        0.05
    );
    legend->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "MM2",
        cuts,
        100,
        -0.1,
        0.1,
        "max",
        legend
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    legend->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_MM2%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    legend->Clear();
    
    // =================== Beam Energy ===================
    TLegend *leg_beam = new TLegend(0.12, 0.8, 0.88, 0.88); // custom legend location
    leg_beam->SetNColumns(4); // spread out across top of hist
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "EnPB",
        50,
        8.0,
        9.0
    );
    h_og->SetXTitle("Beam Energy (GeV)");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    leg_beam->AddEntry(h_og, "Original Data", "l");

    cuts = remove_key_from_cuts(
        cut_color_map,
        "EnPB"
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "EnPB",
        cuts,
        50,
        8.0,
        9.0
    );
    leg_beam->AddEntry(h_cut, "After Cuts", "l");

    h_selection = get_selection_histogram(
        h_cut,
        8.2,
        8.8
    );
    leg_beam->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "EnPB",
        cuts,
        50,
        8.0,
        9.0,
        "integral",
        leg_beam
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    leg_beam->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_Beam_Energy%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    delete leg_beam;

    //=================== Chi2 / NDF ===================
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "Chi2DOF",
        40,
        0.0,
        8.0
    );
    h_og->SetXTitle("#chi^{2} / NDF");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f", bin_width));
    legend->AddEntry(h_og, "Original Data", "l");

    cuts = remove_key_from_cuts(
        cut_color_map,
        "chi2"
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "chi2",
        cuts,
        40,
        0.0,
        8.0
    );
    legend->AddEntry(h_cut, "After Cuts", "l");

    h_selection = get_selection_histogram(
        h_cut,
        0.0,
        5.0
    );
    legend->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "chi2",
        cuts,
        40,
        0.0,
        8.0,
        "max",
        legend
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    legend->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_Chi2%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    legend->Clear();

    // =================== four momentum transfer -t ===================
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "abs(-1*MASS2([proton],-GLUEXTARGET))",
        100,
        0.0,
        2.5
    );
    h_og->SetXTitle("-t (GeV^{2})");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    legend->AddEntry(h_og, "Original Data", "l");

    cuts = remove_key_from_cuts(
        cut_color_map,
        "t"
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "t",
        cuts,
        100,
        0.0,
        2.5
    );
    legend->AddEntry(h_cut, "After Cuts", "l");

    h_selection = get_selection_histogram(
        h_cut,
        0.0,
        1.0
    );
    legend->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "t",
        cuts,
        100,
        0.0,
        2.5,
        "max",
        legend
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    legend->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_t%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    legend->Clear();

    // =================== Missing Energy  ===================
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "RENERGY(GLUEXTARGET, B) - RENERGY(1, 2, 3, 4, 5)",
        100,
        -3.0,
        3.0
    );
    h_og->SetXTitle("Missing Energy (GeV)");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    legend->AddEntry(h_og, "Original Data", "l");

    // missing Energy cut is not included, so no need to remove it from cuts
    cuts = remove_key_from_cuts(
        cut_color_map,
        ""
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "missingE",
        cuts,
        100,
        -3.0,
        3.0
    );
    legend->AddEntry(h_cut, "After Cuts", "l");

    // selection same as analysis launch window
    h_selection = get_selection_histogram(
        h_cut,
        -3.0,
        3.0
    );
    legend->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "missingE",
        cuts,
        100,
        -3.0,
        3.0,
        "max",
        legend
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    legend->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_Missing_Energy%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    legend->Clear();

    // =================== Omega Mass ===================    
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "MASS(2,3,4)",
        100,
        0.3,
        1.8
    );
    TH1F* h_og_perm2 = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "MASS(2,3,5)",
        100,
        0.3,
        1.8
    );
    h_og->Add(h_og_perm2);
    h_og->SetXTitle("#pi^{+}#pi^{-}#pi^{0} inv. mass (GeV)");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    legend->AddEntry(h_og, "Original Data", "l");

    // no cut on omega mass to exclude in broad cuts
    cuts = remove_key_from_cuts(
        cut_color_map,
        ""
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "MASS(3,4,5)",
        cuts,
        100,
        0.3,
        1.8
    );
    legend->AddEntry(h_cut, "After Cuts", "l");

    h_selection = get_selection_histogram(
        h_cut,
        0.755,
        0.811
    );
    legend->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "MASS(3,4,5)",
        cuts,
        100,
        0.3,
        1.8,
        "max",
        legend
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    legend->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_Omega_Mass%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    legend->Clear();

    // =================== Omega Pi0 Mass ===================
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "MASS(2,3,4,5)",
        100,
        1.0,
        2.0
    );    
    h_og->SetXTitle("#omega#pi^{0} inv. mass (GeV)");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    legend->AddEntry(h_og, "Original Data", "l");

    // no cut on omega mass to exclude in broad cuts
    cuts = remove_key_from_cuts(
        cut_color_map,
        ""
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "MASS(2,3,4,5)",
        cuts,
        100,
        1.0,
        2.0
    );
    legend->AddEntry(h_cut, "After Cuts", "l");

    h_selection = get_selection_histogram(
        h_cut,
        1.0,
        2.0
    );
    legend->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "MASS(2,3,4,5)",
        cuts,
        100,
        1.0,
        2.0,
        "max",
        legend
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    legend->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_Omega_Pi0_Mass%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    legend->Clear();

    // =================== Proton Bachelor Mass ===================
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "MASS(1,5)",
        400,
        1.0,
        4.0
    );
    h_og_perm2 = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "MASS(1,4)",
        400,
        1.0,
        4.0
    );
    h_og->Add(h_og_perm2);
    h_og->SetXTitle("p'#pi^{0}_{bachelor} inv. mass (GeV)");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    legend->AddEntry(h_og, "Original Data", "l");

    // no cut on proton pi0 mass to exclude in broad cuts
    cuts = remove_key_from_cuts(
        cut_color_map,
        ""
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "MRecoilPi",
        cuts,
        400,
        1.0,
        4.0
    );
    legend->AddEntry(h_cut, "After Cuts", "l");

    // no selection region for this variable

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "MRecoilPi",
        cuts,
        400,
        1.0,
        4.0,
        "integral",
        legend
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);        
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");    
    h_mc->Draw("SAME");
    legend->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_Proton_Bachelor_Mass%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    legend->Clear();

    // =================== pi0 (omega) momenta ===================
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "MOMENTUMZBOOST(4;B,GLUEXTARGET)",
        200,
        -0.5,
        1.5
    );
    h_og_perm2 = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "MOMENTUMZBOOST(5;B,GLUEXTARGET)",
        200,
        -0.5,
        1.5
    );
    h_og->Add(h_og_perm2);
    h_og->SetXTitle("P_{z}^{CM}(#pi^{0}_{omega}) (GeV)");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    legend->AddEntry(h_og, "Original Data", "l");

    // no cut on omega pi0 momenta to exclude in broad cuts
    cuts = remove_key_from_cuts(
        cut_color_map,
        ""
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "MOMENTUMZBOOST(3;B,GLUEXTARGET)",
        cuts,
        200,
        -0.5,
        1.5
    );
    legend->AddEntry(h_cut, "After Cuts", "l");

    h_selection = get_selection_histogram(
        h_cut,
        -0.1,
        1.5
    );
    legend->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "MOMENTUMZBOOST(3;B,GLUEXTARGET)",
        cuts,
        200,
        -0.5,
        1.5,
        "integral",
        legend
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    legend->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_pz_omega%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    legend->Clear();

    // =================== pi0 (bachelor) momenta ===================
    h_og = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "MOMENTUMZBOOST(4;B,GLUEXTARGET)",
        200,
        -0.5,
        1.5
    );
    h_og_perm2 = get_og_histogram(
        input_data_og,
        NT,
        CATEGORY,
        "MOMENTUMZBOOST(5;B,GLUEXTARGET)",
        200,
        -0.5,
        1.5
    );
    h_og->Add(h_og_perm2);
    h_og->SetXTitle("P_{z}^{CM}(#pi^{0}_{bachelor}) (GeV)");
    bin_width = get_bin_width(h_og);
    h_og->SetYTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
    legend->AddEntry(h_og, "Original Data", "l");

    cuts = remove_key_from_cuts(
        cut_color_map,
        "pzPi0"
    );
    h_cut = get_cut_histogram(
        input_data_signal,
        input_data_sideband,
        NT,
        CATEGORY,
        "MOMENTUMZBOOST(2;B,GLUEXTARGET)",
        cuts,
        200,
        -0.5,
        1.5
    );
    legend->AddEntry(h_cut, "After Cuts", "l");

    h_selection = get_selection_histogram(
        h_cut,
        -0.1,
        1.5
    );
    legend->AddEntry(
        h_selection, 
        "Selection",
        "f"
    );

    h_mc = get_mc_histogram(
        h_cut,
        input_mc_signal,
        input_mc_sideband,
        NT,
        CATEGORY,
        "MOMENTUMZBOOST(2;B,GLUEXTARGET)",
        cuts,
        200,
        -0.5,
        1.5,
        "integral",
        legend
    );

    if (log_scale)
    {
        h_og->SetMinimum(1);
        h_cut->SetMinimum(1);
        h_selection->SetMinimum(1);
        h_mc->SetMinimum(1);
    }

    h_og->Draw("HIST");
    h_cut->Draw("HIST SAME");
    h_selection->Draw("HIST SAME");
    h_mc->Draw("SAME");
    legend->Draw("SAME");
    c->SaveAs(TString::Format("GlueXI_pz_bachelor%s.pdf", 
        log_scale ? "_log" : ""));
    c->Clear();
    legend->Clear();

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
    // TString friend_shower_vars[4] = {"shQuality_P3a", "shQuality_P3b", "shQuality_P2a", "shQuality_P2b"};
    // TString photon_names[4] = {"Photon 1 (#pi^{0}_{1})", "Photon 2 (#pi^{0}_{1})",
    //                            "Photon 3 (#pi^{0}_{2})", "Photon 4 (#pi^{0}_{2})"};

    // for (int i = 0; i < 4; i++)
    // {
    //     c_shower->cd(2);
    //     gPad->cd(i + 1);
    //     gPad->SetLogy(1);

    //     TH1F *h_og = get_og_histogram(
    //         input_data_files,
    //         NT,
    //         CATEGORY,
    //         shower_vars[i],
    //         "100",
    //         "0.0",
    //         "1.0");

    //     if (i == 0)
    //         common_legend->AddEntry(h_og, "Original Data", "l");

    //     // apply the special rf cut first
    //     TH1F *h_rf = get_rf_histogram(
    //         input_data_files,
    //         NT,
    //         CATEGORY,
    //         shower_vars[i],
    //         "100",
    //         "0.0",
    //         "1.0"
    //     );
    //     if (i == 0)
    //         common_legend->AddEntry(h_rf, "RF", "l");

    //     // now iterate through our cuts and apply them one by one, except the one we are plotting
    //     bool fill_legend = (i == 0); // only fill legend on first plot
    //     std::vector<TH1F*> iterated_hists = get_histograms_with_cuts(
    //         input_data_files,
    //         NT,
    //         CATEGORY,
    //         shower_vars[i],
    //         "100",
    //         "0.0",
    //         "1.0",
    //         "shQuality",
    //         cut_color_map,
    //         common_legend,
    //         fill_legend
    //     );        

    //     // Create a copy of the final histogram for the selection region fill
    //     TH1F *h_selection = get_selection_histogram(
    //         iterated_hists.back(),
    //         0.5,
    //         1.0
    //     );
    //     if (i == 0)
    //         common_legend->AddEntry(h_selection, "Selection", "f");
        
    //     auto sh_map = cut_color_map;
    //     sh_map.erase("shQuality");
    //     TString sh_cuts = join_keys(sh_map);
    //     TH1F *h_mc = FSModeHistogram::getTH1F(
    //         input_mc_files,
    //         NT,
    //         CATEGORY,
    //         shower_vars[i],
    //         "(100,0.0,1.0)",
    //         TString::Format("CUT(%s)*CUTWT(rf)", sh_cuts.Data()));

    //     double scale;
    //     // scale MC to data using the maximum bin content
    //     scale = iterated_hists.back()->GetMaximum() / h_mc->GetMaximum();
    //     h_mc->Scale(scale);
    //     TString legend_scale = TString::Format("MC (max scaling=%.3f)", scale);
    //     if (i == 0)
    //         common_legend->AddEntry(h_mc, legend_scale, "p");

    //     // draw MC on top as blue points with error bars
    //     h_mc->SetMarkerStyle(24);
    //     h_mc->SetMarkerColor(kBlack);
    //     h_mc->SetLineColor(kBlack);        

    //     // set minimums for log scale
    //     h_og->SetMinimum(1); 
    //     h_rf->SetMinimum(1);        
    //     for (auto h_it : iterated_hists)
    //     {
    //         h_it->SetMinimum(1);
    //     }
    //     h_selection->SetMinimum(1);
    //     h_mc->SetMinimum(1);
        
    //     bin_width = get_bin_width(h_og);
    //     h_og->SetXTitle(TString::Format("%s Quality", photon_names[i].Data()));
    //     h_og->SetYTitle(TString::Format("Events / %.3f", bin_width));

    //     // draw all hists
    //     h_og->Draw("HIST");
    //     h_rf->Draw("HIST SAME");
    //     for (auto h_it : iterated_hists)
    //     {
    //         h_it->Draw("HIST SAME");
    //     }        
    //     h_selection->Draw("HIST SAME");
    //     h_mc->Draw("E0 SAME");

    //     // Adjust margins for the grid layout
    //     gPad->SetLeftMargin(0.15);
    //     gPad->SetBottomMargin(0.15);
    //     gPad->SetTopMargin(0.1);
    //     gPad->SetRightMargin(0.05);
    // }
    // c_shower->cd(1);
    // common_legend->Draw();

    // c_shower->Update();
    // c_shower->SaveAs(TString::Format("%s_Shower_Quality.pdf", PERIOD_TO_LABEL.at(period).Data()));
    // c_shower->Clear();

    if (dump_cache)
        FSHistogram::dumpHistogramCache();

    return;
}


/**
 * @brief Remove a key from the cut color map and return the joined keys
 * 
 * @param cut_color_map Map of cut names to colors
 * @param key_to_remove Key to remove from the map
 * @return TString Joined keys after removal
 */
TString remove_key_from_cuts(
    std::map<TString, Int_t> &cut_color_map, 
    const TString key_to_remove)
{
    std::map<TString, Int_t> map_copy = cut_color_map;
    std::map<TString, Int_t>::iterator it = map_copy.find(key_to_remove);
    if (it != map_copy.end()) {
        map_copy.erase(it);
    }
    return join_keys(map_copy);
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
    const int bins,
    const double lower_bound,
    const double upper_bound
)
{
    TH1F *h_og = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%d,%f,%f)", bins, lower_bound, upper_bound),
        "");

    h_og->SetLineColor(kGray);
    h_og->SetLineWidth(1);
    h_og->SetTitle(""); // this will be the first plotted histogram, so remove title

    return h_og;
}

/**
 * @brief Get the histogram with all cuts and sideband subtraction applied
 * 
 * @param input_data_signal data consisting of only omega,rf signal events
 * @param input_data_sideband data consisting of only omega,rf sideband events
 * @param NT name of the tree in the ROOT file
 * @param CATEGORY category name in the ROOT file
 * @param tree_variable tree variable name in the ROOT file
 * @param cuts cuts string to apply to the data
 * @param bins number of bins for the histogram
 * @param lower_bound lower bound of the histogram range
 * @param upper_bound upper bound of the histogram range
 * @return TH1F* histogram with all cuts and sideband subtraction applied
 */
TH1F* get_cut_histogram(
    const TString input_data_signal,
    const TString input_data_sideband,
    const TString NT,
    const TString CATEGORY,
    const TString tree_variable,
    const TString cuts,
    const int bins,
    const double lower_bound,
    const double upper_bound
)
{
    TH1F *h_signal = FSModeHistogram::getTH1F(
        input_data_signal,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%d,%f,%f)", bins, lower_bound, upper_bound),
        TString::Format("CUT(%s)", cuts.Data())
    );
    TH1F *h_sideband = FSModeHistogram::getTH1F(
        input_data_sideband,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%d,%f,%f)", bins, lower_bound, upper_bound),
        TString::Format("CUT(%s)", cuts.Data())
    );
    h_signal->Add(h_sideband, -1.0); // subtract sideband from signal
    h_signal->SetLineColor(kBlue);
    h_signal->SetLineWidth(1);
    return h_signal;
}

/**
 * @brief Get the selection histogram object
 * 
 * @param[in] h_source histogram to base the selection histogram on
 * @param[in] cut_lower_bound lower bound for the cut on the variable
 * @param[in] cut_upper_bound upper bound for the cut on the variable
 * @return TH1F* histogram with selection region highlighted
 */
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
 * @param h_data data to scale MC to
 * @param input_mc_signal MC consisting of only omega,rf signal events
 * @param input_mc_sideband MC consisting of only omega,rf sideband events
 * @param NT name of the tree in the ROOT file
 * @param CATEGORY category name in the ROOT file
 * @param tree_variable tree variable name in the ROOT file
 * @param cuts cuts string to apply to the data
 * @param bins number of bins for the histogram
 * @param lower_bound lower bound of the histogram range
 * @param upper_bound upper bound of the histogram range
 * @param scale_choice choice of scaling method ("integral" or "max")
 * @param legend legend to add the MC histogram entry to
 * @return TH1F* scaled MC histogram
 */
TH1F* get_mc_histogram(
    TH1F* h_data,
    const TString input_mc_signal,
    const TString input_mc_sideband,
    const TString NT,
    const TString CATEGORY,
    const TString tree_variable,
    const TString cuts,
    const int bins,
    const double lower_bound,
    const double upper_bound,
    const TString scale_choice,
    TLegend *legend
)
{
    TH1F *h_signal = FSModeHistogram::getTH1F(
        input_mc_signal,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%d,%f,%f)", bins, lower_bound, upper_bound),
        TString::Format("CUT(%s)", cuts.Data())
    );
    TH1F *h_sideband = FSModeHistogram::getTH1F(
        input_mc_sideband,
        NT,
        CATEGORY,
        tree_variable,
        TString::Format("(%d,%f,%f)", bins, lower_bound, upper_bound),
        TString::Format("CUT(%s)", cuts.Data())
    );
    h_signal->Add(h_sideband, -1.0); // subtract sideband from signal

    double scale;
    TString legend_scale;
    if (scale_choice == "integral") {
        // scale MC to data using the integrals
        scale = h_data->Integral() / h_signal->Integral();
        h_signal->Scale(scale);
        legend_scale = TString::Format("MC (integral scaling=%.3f)", scale);
    } else if (scale_choice == "max") {
        // scale MC to data using the maximum bin content
        scale = h_data->GetMaximum() / h_signal->GetMaximum();
        h_signal->Scale(scale);
        legend_scale = TString::Format("MC (max scaling=%.3f)", scale);
    } else {
        throw std::invalid_argument("Invalid scale_choice");
    }
    legend->AddEntry(h_signal, legend_scale, "p");
    h_signal->SetMarkerStyle(24);
    h_signal->SetMarkerColor(kBlack);    
    return h_signal;
}