/**
 * @file plot_final_distributions.cc
 * @author Kevin Scheuer
 * @brief Make plots from finalized amptools trees to show what goes into PWA
 * 
 * These plots do not distinguish between different polarizations or run periods,
 * they are just the final distributions after all cuts and sideband subtractions
 * have been applied.
 */

#include <iostream>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "TString.h"
#include "TH1.h"
#include "TLegend.h"

#include "FSBasic/FSTree.h"
#include "FSMode/FSModeTree.h"
#include "FSMode/FSModeHistogram.h"

#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"

// forward declarations
void plot_mass_spectra(
    TString NT,
    TString CATEGORY,
    TString data_signal,
    TString data_sideband,
    TString mc_signal,
    TString mc_sideband
);

void plot_final_distributions()
{
    TString tree_dir = "/lustre24/expphy/volatile/halld/home/kscheuer/"
    "FSRoot-skimmed-trees/final-amptools-trees/";

    // input files 
    TString data_signal = tree_dir + "allPeriods_data_signal.root";
    TString data_sideband = tree_dir + "allPeriods_data_sideband.root";
    TString mc_signal = tree_dir + "allPeriods_ver03.1_mc_signal.root";
    TString mc_sideband = tree_dir + "allPeriods_ver03.1_mc_sideband.root";
    TString acc_phasespace = tree_dir + "allPeriods_ver03_phasespace.root";
    TString gen_phasespace = tree_dir + "allPeriods_ver03_gen_phasespace.root";

    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(false);

    plot_mass_spectra(NT, CATEGORY, data_signal, data_sideband, mc_signal, mc_sideband);

    plot_1d_angles(
        NT,
        CATEGORY,
        data_signal,
        data_sideband,
        mc_signal,
        mc_sideband,
        acc_phasespace,
        gen_phasespace
    );
    /* TODO: implement
        - 2D histograms of angle correlations (data, signal MC)
        - 2D acceptance histograms of angle correlations
        - omega pi0 and proton pi0 mass spectrum
        - lambda distribution
    
        angles will need to import, or at least copy the vec ps angles calculations

        for 2D, use col0 and plot data and signal MC separately
    */
}


/**
 * @brief Plot the omega pi0 and proton pi0 mass spectra for data and MC
 * 
 * @param NT tree name
 * @param CATEGORY category name
 * @param data_signal data events from signal file
 * @param data_sideband data events from sideband (background) file
 * @param mc_signal MC events from signal file
 * @param mc_sideband MC events from sideband (background) file
 */
void plot_mass_spectra(
    TString NT,
    TString CATEGORY,
    TString data_signal,
    TString data_sideband,
    TString mc_signal,
    TString mc_sideband
)
{
    // get histograms
    TH1F* h_omega_pi0_data_signal = FSModeHistogram::getTH1F(
        data_signal,
        NT,
        CATEGORY,
        "M4Pi",
        "(100,1.0,2.0)",
        ""        
    );
    TH1F *h_omega_pi0_data_sideband = FSModeHistogram::getTH1F(
        data_sideband,
        NT,
        CATEGORY,
        "M4Pi",
        "(100,1.0,2.0)",
        ""        
    );
    TH1F* h_omega_pi0_mc_signal = FSModeHistogram::getTH1F(
        mc_signal,
        NT,
        CATEGORY,
        "M4Pi",
        "(100,1.0,2.0)",
        ""        
    );
    TH1F *h_omega_pi0_mc_sideband = FSModeHistogram::getTH1F(
        mc_sideband,
        NT,
        CATEGORY,
        "M4Pi",
        "(100,1.0,2.0)",
        ""        
    );
    TH1F *h_omega_pi0_data_total = 
        (TH1F*)h_omega_pi0_data_signal->Clone("h_omega_pi0_data_total");
    h_omega_pi0_data_total->Add(h_omega_pi0_data_sideband);

    TH1F* h_proton_pi0_data_signal = FSModeHistogram::getTH1F(
        data_signal,
        NT,
        CATEGORY,
        "MRecoilPi",
        "(200,1.0,3.0)",
        ""        
    );
    TH1F *h_proton_pi0_data_sideband = FSModeHistogram::getTH1F(
        data_sideband,
        NT,
        CATEGORY,
        "MRecoilPi",
        "(200,1.0,3.0)",
        ""        
    );
    TH1F* h_proton_pi0_mc_signal = FSModeHistogram::getTH1F(
        mc_signal,
        NT,
        CATEGORY,
        "MRecoilPi",
        "(200,1.0,3.0)",
        ""        
    );
    TH1F *h_proton_pi0_mc_sideband = FSModeHistogram::getTH1F(
        mc_sideband,
        NT,
        CATEGORY,
        "MRecoilPi",
        "(200,1.0,3.0)",
        ""        
    );
    TH1F *h_proton_pi0_data_total = 
        (TH1F*)h_proton_pi0_data_signal->Clone("h_proton_pi0_data_total");
    h_proton_pi0_data_total->Add(h_proton_pi0_data_sideband);

    // customize the histograms
    double bin_width = get_bin_width( h_omega_pi0_data_total );
    h_omega_pi0_data_total->SetTitle("");
    h_omega_pi0_data_total->SetXTitle("#omega#pi^{0} inv. mass (GeV)");
    h_omega_pi0_data_total->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    h_omega_pi0_data_total->SetMarkerStyle(1);
    h_omega_pi0_data_total->SetMarkerColor(kBlack);

    h_omega_pi0_data_signal->SetLineColor(kBlue);
    h_omega_pi0_data_signal->SetLineWidth(2);

    h_omega_pi0_mc_signal->SetLineColor(kCyan+2);
    h_omega_pi0_mc_signal->SetLineWidth(1);
    h_omega_pi0_mc_signal->SetLineStyle(2);

    h_omega_pi0_data_sideband->SetLineColor(kRed+1);
    h_omega_pi0_data_sideband->SetLineWidth(2);

    h_omega_pi0_mc_sideband->SetLineColor(kRed-7);
    h_omega_pi0_mc_sideband->SetLineWidth(1);
    h_omega_pi0_mc_sideband->SetLineStyle(2);


    bin_width = get_bin_width( h_proton_pi0_data_total );
    h_proton_pi0_data_total->SetTitle("");
    h_proton_pi0_data_total->SetXTitle("p'#pi^{0} inv. mass (GeV)");
    h_proton_pi0_data_total->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    h_proton_pi0_data_total->SetMarkerStyle(1);
    h_proton_pi0_data_total->SetMarkerColor(kBlack);

    h_proton_pi0_data_signal->SetLineColor(kBlue);
    h_proton_pi0_data_signal->SetLineWidth(2);

    h_proton_pi0_mc_signal->SetLineColor(kCyan+2);
    h_proton_pi0_mc_signal->SetLineWidth(1);
    h_proton_pi0_mc_signal->SetLineStyle(2);

    h_proton_pi0_data_sideband->SetLineColor(kRed+1);
    h_proton_pi0_data_sideband->SetLineWidth(2);

    h_proton_pi0_mc_sideband->SetLineColor(kRed-7);
    h_proton_pi0_mc_sideband->SetLineWidth(1);
    h_proton_pi0_mc_sideband->SetLineStyle(2);

    // create legend
    TLegend* legend_omega_pi0 = new TLegend(0.6,0.6,0.89,0.89);
    legend_omega_pi0->AddEntry(h_omega_pi0_data_total, "Data", "lp");
    legend_omega_pi0->AddEntry(h_omega_pi0_data_signal, "Data Signal", "l");
    legend_omega_pi0->AddEntry(h_omega_pi0_mc_signal, "MC Signal", "l");
    legend_omega_pi0->AddEntry(h_omega_pi0_data_sideband, "Data Sideband", "l");
    legend_omega_pi0->AddEntry(h_omega_pi0_mc_sideband, "MC Sideband", "l");

    TLegend* legend_proton_pi0 = new TLegend(0.6,0.6,0.89,0.89);
    legend_proton_pi0->AddEntry(h_proton_pi0_data_total, "Data", "lp");
    legend_proton_pi0->AddEntry(h_proton_pi0_data_signal, "Data Signal", "l");
    legend_proton_pi0->AddEntry(h_proton_pi0_mc_signal, "MC Signal", "l");
    legend_proton_pi0->AddEntry(h_proton_pi0_data_sideband, "Data Sideband", "l");
    legend_proton_pi0->AddEntry(h_proton_pi0_mc_sideband, "MC Sideband", "l");

    // draw the histograms
    h_omega_pi0_data_total->Draw("E");
    h_omega_pi0_data_signal->Draw("HIST SAME");
    h_omega_pi0_mc_signal->Draw("HIST SAME");
    h_omega_pi0_data_sideband->Draw("HIST SAME");
    h_omega_pi0_mc_sideband->Draw("HIST SAME");
    legend_omega_pi0->Draw("SAME");

    h_proton_pi0_data_total->Draw("E");
    h_proton_pi0_data_signal->Draw("HIST SAME");
    h_proton_pi0_mc_signal->Draw("HIST SAME");
    h_proton_pi0_data_sideband->Draw("HIST SAME");
    h_proton_pi0_mc_sideband->Draw("HIST SAME");
    legend_proton_pi0->Draw("SAME");

    return;
}

void plot_1d_angles(
    TString NT,
    TString CATEGORY,
    TString data_signal,
    TString data_sideband,
    TString mc_signal,
    TString mc_sideband,
    TString acc_phasespace,
    TString gen_phasespace
)
{  
    // - 1D histograms of the 5 angles (data, signal MC), with acceptance line
}