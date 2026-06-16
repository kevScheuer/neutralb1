/**
 * @file mass_spectrum.cc
 * @author Kevin Scheuer
 * @brief Simple script to plot the final omega pi0 mass spectrum for data
 * 
 */

#include <iostream>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TPad.h"
#include "TStyle.h"

#include "neutralb1/fit_utils.h"


void mass_spectrum()
{
    gStyle->SetOptStat(0);
    TString tree_dir = "/home/kscheuer/work/neutralb1/data/FSRoot/GlueX/";

    // input files 
    TString data_signal = tree_dir + "allPeriods_data_signal.root";
    TString data_sideband = tree_dir + "allPeriods_data_background.root";

    TString NT = "ntFSGlueX_100_112";
    TString CATEGORY = "pi0pi0pippim";    

    // get files and trees
    TFile *f_data_signal = TFile::Open(data_signal);
    TFile *f_data_sideband = TFile::Open(data_sideband);
    TTree *tree_data_signal = (TTree*)f_data_signal->Get(NT);
    TTree *tree_data_sideband = (TTree*)f_data_sideband->Get(NT);

    // setup histograms for omega pi0 and proton pi0
    TH1F *h_omega_pi0_data_signal = new TH1F(
        "h_omega_pi0_data_signal", 
        "", 
        50, 1.0, 2.0);
    TH1F *h_omega_pi0_data_sideband = new TH1F(
        "h_omega_pi0_data_sideband", 
        "", 
        50, 1.0, 2.0);

    // fill histograms
    tree_data_signal->Draw("M4Pi>>h_omega_pi0_data_signal", "", "goff");
    tree_data_sideband->Draw("M4Pi>>h_omega_pi0_data_sideband", "weight", "goff");
    
    
    // get final data as signal minus sideband
    TH1F *h_omega_pi0_data_total = 
        (TH1F*)h_omega_pi0_data_signal->Clone("h_omega_pi0_data_total");
    h_omega_pi0_data_total->Add(h_omega_pi0_data_sideband, -1.0);

    // customize the histograms
    double bin_width = get_bin_width( h_omega_pi0_data_total );
    h_omega_pi0_data_total->SetTitle("");
    h_omega_pi0_data_total->SetXTitle("#omega#pi^{0} inv. mass (GeV)");
    h_omega_pi0_data_total->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    h_omega_pi0_data_total->GetXaxis()->SetTitleSize(0.05);
    h_omega_pi0_data_total->GetYaxis()->SetTitleSize(0.05);
    h_omega_pi0_data_total->GetXaxis()->SetTitleOffset(0.8);
    h_omega_pi0_data_total->GetYaxis()->SetTitleOffset(0.8);
    h_omega_pi0_data_total->SetMarkerStyle(8);
    h_omega_pi0_data_total->SetMarkerSize(1.02);
    h_omega_pi0_data_total->SetMarkerColor(kBlack);
    h_omega_pi0_data_total->SetLineColor(kBlack);
    h_omega_pi0_data_total->SetMinimum(0);

    // create legend
    TLegend* legend_omega_pi0 = new TLegend(0.55,0.78,0.85,0.88);
    legend_omega_pi0->SetFillColorAlpha(kWhite, 0.8);
    legend_omega_pi0->AddEntry(h_omega_pi0_data_total, "GlueX-I Data", "lp");
    legend_omega_pi0->SetTextSize(0.04);

    TCanvas *c = new TCanvas("c", "c", 800, 600);

    // draw the histograms
    h_omega_pi0_data_total->Draw("E");
    legend_omega_pi0->Draw("SAME");

    c->SaveAs("simple_mass_spectrum.pdf");
    c->Clear();

    return;
}