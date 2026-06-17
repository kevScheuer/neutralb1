/**
 * @file fit_t_distribution.cc
 * @author Kevin Scheuer
 * @brief fit the -t distribution of the data for MC generation later
  * @date 2026-06-17
  * 
 */

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TString.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TText.h"

#include "neutralb1/fit_utils.h"

void fit_t_distribution() {
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

    if (!f_data_signal || f_data_signal->IsZombie()) {
        std::cerr << "Failed to open signal file: " << data_signal << std::endl;
        return;
    }

    if (!f_data_sideband || f_data_sideband->IsZombie()) {
        std::cerr << "Failed to open sideband file: " << data_sideband << std::endl;
        return;
    }

    TTree *tree_data_signal = (TTree*)f_data_signal->Get(NT);
    TTree *tree_data_sideband = (TTree*)f_data_sideband->Get(NT);

    if (!tree_data_signal) {
        std::cerr << "Failed to find tree " << NT << " in " << data_signal << std::endl;
        return;
    }

    if (!tree_data_sideband) {
        std::cerr << "Failed to find tree " << NT << " in " << data_sideband << std::endl;
        return;
    }

    // setup histograms
    TH1F *h_t_signal = new TH1F(
        "h_t_signal", 
        "", 
        50, 0.0, 1.0);
    TH1F *h_t_sideband = new TH1F(
        "h_t_sideband", 
        "", 
        50, 0.0, 1.0);

    // fill histograms
    tree_data_signal->Draw("t>>h_t_signal", "", "goff");
    tree_data_sideband->Draw("t>>h_t_sideband", "weight", "goff");
    
    
    // get final data as signal minus sideband
    TH1F *h_t_total = 
        (TH1F*)h_t_signal->Clone("h_t_total");
    h_t_total->Add(h_t_sideband, -1.0);

    // customize the histograms
    double bin_width = get_bin_width( h_t_total );
    h_t_total->SetTitle("");
    h_t_total->SetXTitle("-t (GeV)^{2}");
    h_t_total->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    h_t_total->GetXaxis()->SetTitleSize(0.05);
    h_t_total->GetYaxis()->SetTitleSize(0.05);
    h_t_total->GetXaxis()->SetTitleOffset(0.8);
    h_t_total->GetYaxis()->SetTitleOffset(0.8);
    h_t_total->SetMarkerStyle(1);    
    h_t_total->SetMarkerColor(kBlack);
    h_t_total->SetLineColor(kBlack);
    h_t_total->SetMinimum(0);


    // fit with exponential
    TF1 *f_expo_15 = new TF1("f_expo_15", "[0] * exp(-[1] * x)", 0.16, 1.0);
    TF1 *f_expo_2 = new TF1("f_expo_2", "[0] * exp(-[1] * x)", 0.2, 1.0);
    TF1 *f_expo_25 = new TF1("f_expo_25", "[0] * exp(-[1] * x)", 0.26, 1.0);

    double max_val = h_t_total->GetMaximum();

    // initial guess for scale
    f_expo_15->SetParameters(max_val, 1.0);
    f_expo_2->SetParameters(max_val, 1.0);
    f_expo_25->SetParameters(max_val, 1.0);

    h_t_total->Fit(f_expo_15, "QR0");
    h_t_total->Fit(f_expo_2, "QR0");
    h_t_total->Fit(f_expo_25, "QR0");
    
    double scale_15 = f_expo_15->GetParameter(0);
    double scale_2 = f_expo_2->GetParameter(0);
    double scale_25 = f_expo_25->GetParameter(0);

    double slope_15 = f_expo_15->GetParameter(1);
    double slope_2 = f_expo_2->GetParameter(1);
    double slope_25 = f_expo_25->GetParameter(1);

    double scale_15_err = f_expo_15->GetParError(0);
    int slope_15_err = std::round(f_expo_15->GetParError(1)*1000);
    double scale_2_err = f_expo_2->GetParError(0);
    int slope_2_err = std::round(f_expo_2->GetParError(1)*1000);
    double scale_25_err = f_expo_25->GetParError(0);
    int slope_25_err = std::round(f_expo_25->GetParError(1)*1000);

    std::cout << "0.16 < -t < 1.0 - Scale: " << scale_15 << " +/- " << scale_15_err << std::endl;
    std::cout << "\tSlope: " << slope_15 << " +/- " << slope_15_err << std::endl;
    std::cout << "0.2 < -t < 1.0 - Scale: " << scale_2 << " +/- " << scale_2_err << std::endl;
    std::cout << "\tSlope: " << slope_2 << " +/- " << slope_2_err << std::endl;
    std::cout << "0.26 < -t < 1.0 - Scale: " << scale_25 << " +/- " << scale_25_err << std::endl;
    std::cout << "\tSlope: " << slope_25 << " +/- " << slope_25_err << std::endl;

    TCanvas *c = new TCanvas("c", "c", 800, 600);

    // draw the histograms
    h_t_total->Draw("E");
    
    f_expo_15->SetLineColorAlpha(kRed, 0.7);
    f_expo_2->SetLineColorAlpha(kBlue, 0.7);
    f_expo_25->SetLineColorAlpha(kViolet, 0.7);

    f_expo_15->SetLineStyle(1);
    f_expo_2->SetLineStyle(2);
    f_expo_25->SetLineStyle(9);

    f_expo_15->Draw("same");
    f_expo_2->Draw("same");
    f_expo_25->Draw("same");

    // add fit parameters to the plot
    TPaveText *pt = new TPaveText(0.50, 0.50, 0.88, 0.89, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.028);
    pt->AddText("Exponential Fits");
    TText *t15 = pt->AddText(TString::Format("0.16 < -t < 1.0: A*exp[-%.3f(%i)*t]", slope_15, slope_15_err));
    TText *t2 = pt->AddText(TString::Format("0.2 < -t < 1.0: A*exp[-%.3f(%i)*t]", slope_2, slope_2_err));
    TText *t25 = pt->AddText(TString::Format("0.26 < -t < 1.0: A*exp[-%.3f(%i)*t]", slope_25, slope_25_err));

    t15->SetTextColor(kRed);
    t2->SetTextColor(kBlue);
    t25->SetTextColor(kViolet);

    pt->Draw();

    c->SaveAs("t_fit.pdf");
    return;
}