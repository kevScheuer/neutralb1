/**
 * @file t_resolution.cc
 * @author Kevin Scheuer
 * @brief Quick script to determine -t resolution from signal MC
 * 
 * Our signal MC files can be used to determine the -t resolution by subtracting 
 * the generated and reconstructed 4-momenta of the photon and proton difference,
 * and fitting the resulting distribution to a Gaussian.
 */

#include <iostream>

#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TPaveText.h"

#include "FSMode/FSModeHistogram.h"

#include "../fsroot_setup.cc"

TString combine_signal_mc_files(TString mc_dir, bool make_file);

void t_resolution(bool make_file = true)
{
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(false);
    
    // first, combine the non-AmpTools signal MC files into one big signal file. We
    // don't care about permutations here, since MASS(2,3,4,5) and MASS(2,3,5,4) are the
    // same quantity
    TString mc_dir = "/lustre24/expphy/volatile/halld/home/kscheuer/"
    "FSRoot-skimmed-trees/final-amptools-trees/";
    TString combined_file = combine_signal_mc_files(mc_dir, make_file);

    // get the t difference histogram
    TString t_diff_expr = "(-1*MASS2(1,-GLUEXTARGET) - -1*MCMASS2(1,-GLUEXTARGET))";
    TH1F *h_t = FSModeHistogram::getTH1F(
        combined_file,
        NT,
        CATEGORY,
        t_diff_expr,
        "(100, -0.1, 0.1)",
        ""
    );
    
    // fit a double gaussian and report the average sigma 
    TF1 *double_gaus = new TF1(
        "double_gaus",
        "gaus(0) + gaus(3)",
        -0.1,
        0.1
    );
    
    // Set initial parameter estimates for better convergence
    double max_val = h_t->GetMaximum();
    double mean_val = h_t->GetMean();
    double rms_val = h_t->GetRMS();
    
    // First Gaussian: narrower core
    double_gaus->SetParameter(0, 0.7 * max_val);  // amplitude
    double_gaus->SetParameter(1, mean_val);       // mean
    double_gaus->SetParameter(2, 0.5 * rms_val);  // sigma
    
    // Second Gaussian: wider tail
    double_gaus->SetParameter(3, 0.3 * max_val);  // amplitude
    double_gaus->SetParameter(4, mean_val);       // mean
    double_gaus->SetParameter(5, 1.5 * rms_val);  // sigma
    
    h_t->Fit("double_gaus", "R");
    
    TF1 *fit = h_t->GetFunction("double_gaus");
    double amp1 = fit->GetParameter(0);
    double mean1 = fit->GetParameter(1);
    double sigma1 = fit->GetParameter(2);
    double amp2 = fit->GetParameter(3);
    double mean2 = fit->GetParameter(4);
    double sigma2 = fit->GetParameter(5);
    
    double amp1_err = fit->GetParError(0);
    double mean1_err = fit->GetParError(1);
    double sigma1_err = fit->GetParError(2);
    double amp2_err = fit->GetParError(3);
    double mean2_err = fit->GetParError(4);
    double sigma2_err = fit->GetParError(5);
    
    double sigma_avg = (sigma1 + sigma2) / 2.0;
    double sigma_avg_err = sqrt(sigma1_err*sigma1_err + sigma2_err*sigma2_err) / 2.0;
    std::cout << "Mass resolution (sigma): " << sigma_avg << " +/- " << sigma_avg_err << " GeV\n";
    
    // draw histogram with fit and save as PDF
    TCanvas *c = new TCanvas("c_t_diff", "t Resolution", 800, 600);
    h_t->SetTitle(
        ";-t_{kinfit} - -t_{generated} [GeV^{2}];Entries"
    );
    h_t->Draw();
    fit->Draw("same");
    
    // Add fit parameters to plot
    TPaveText *pt = new TPaveText(0.58, 0.50, 0.88, 0.89, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.028);
    pt->AddText("Double Gaussian Fit");
    pt->AddText(TString::Format("Gaus 1: A = %.1f #pm %.1f", amp1, amp1_err));
    pt->AddText(TString::Format("#mu = %.4f #pm %.4f GeV^{2}", mean1, mean1_err));
    pt->AddText(TString::Format("#sigma = %.4f #pm %.4f GeV^{2}", sigma1, sigma1_err));
    pt->AddText(TString::Format("Gaus 2: A = %.1f #pm %.1f", amp2, amp2_err));
    pt->AddText(TString::Format("#mu = %.4f #pm %.4f GeV^{2}", mean2, mean2_err));
    pt->AddText(TString::Format("#sigma = %.4f #pm %.4f GeV^{2}", sigma2, sigma2_err));
    pt->AddText(TString::Format("Avg #sigma = %.4f #pm %.4f GeV^{2}", sigma_avg, sigma_avg_err));
    pt->Draw();
    
    c->SaveAs("t_resolution.pdf");
    delete c;

    return;
}


TString combine_signal_mc_files(TString mc_dir, bool make_file)
{
    TString input_pattern = mc_dir + 
    "tree_pi0pi0pippim__B4_finalAmptools_SKIM_0*_ver03.1_mc_pol0_perm*_signal.root";

    TString output_file = mc_dir +
    "tree_pi0pi0pippim__B4_finalAmptools_SKIM_allPeriods_ver03.1_mc_pol0_signal.root";


    TString hadd_command = TString::Format(
        "hadd -f %s %s",
        output_file.Data(),
        input_pattern.Data()
    );

    if (make_file)
        system(hadd_command.Data());
        
    return output_file;
}