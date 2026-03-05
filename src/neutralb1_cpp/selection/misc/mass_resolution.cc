/**
 * @file mass_resolution.cc
 * @author Kevin Scheuer
 * @brief Quick script to determine mass resolution from signal MC
 * 
 * Our signal MC files can be used to determine the mass resolution by subtracting 
 * the masses of the generated and reconstructed 4-momenta of the omega pi0 system, 
 * and fitting the resulting distribution to a Gaussian.
 */

#include <iostream>

#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TText.h"
#include "TPaveText.h"

#include "FSMode/FSModeHistogram.h"

#include "../fsroot_setup.cc"

TString combine_signal_mc_files(TString mc_dir, bool make_file);

void mass_resolution(bool make_file = true)
{
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(false);
    
    // first, combine the non-AmpTools signal MC files into one big signal file. We
    // don't care about permutations here, since MASS(2,3,4,5) and MASS(2,3,5,4) are the
    // same quantity
    TString mc_dir = "/lustre24/expphy/volatile/halld/home/kscheuer/"
    "FSRoot-skimmed-trees/final-amptools-trees/";
    TString combined_file = combine_signal_mc_files(mc_dir, make_file);

    // get the mass difference histogram
    TH1F *h_mass_diff = FSModeHistogram::getTH1F(
        combined_file,
        NT,
        CATEGORY,
        "MASS(2,3,4,5)-MCMASS(2,3,4,5)",
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
    double max_val = h_mass_diff->GetMaximum();
    double mean_val = h_mass_diff->GetMean();
    double rms_val = h_mass_diff->GetRMS();
    
    // First Gaussian: narrower core
    double_gaus->SetParameter(0, 0.7 * max_val);  // amplitude
    double_gaus->SetParameter(1, mean_val);       // mean
    double_gaus->SetParameter(2, 0.5 * rms_val);  // sigma
    
    // Second Gaussian: wider tail
    double_gaus->SetParameter(3, 0.3 * max_val);  // amplitude
    double_gaus->SetParameter(4, mean_val);       // mean
    double_gaus->SetParameter(5, 1.5 * rms_val);  // sigma
    
    h_mass_diff->Fit("double_gaus", "R");
    
    TF1 *fit = h_mass_diff->GetFunction("double_gaus");
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

    // Style the combined fit as solid red
    fit->SetLineColor(kRed);
    fit->SetLineStyle(1);
    fit->SetLineWidth(2);

    // Create individual Gaussians from fitted parameters
    TF1 *gaus1_func = new TF1("gaus1_func", "gaus", -0.1, 0.1);
    gaus1_func->SetParameters(amp1, mean1, sigma1);
    gaus1_func->SetLineStyle(3);       // dotted
    gaus1_func->SetLineColor(kAzure+7); // light blue
    gaus1_func->SetLineWidth(2);

    TF1 *gaus2_func = new TF1("gaus2_func", "gaus", -0.1, 0.1);
    gaus2_func->SetParameters(amp2, mean2, sigma2);
    gaus2_func->SetLineStyle(2);       // dashed
    gaus2_func->SetLineColor(kBlue+2); // dark blue
    gaus2_func->SetLineWidth(2);

    // draw histogram with fit and save as PDF
    TCanvas *c = new TCanvas("c_mass_diff", "Mass Resolution", 800, 600);
    h_mass_diff->SetTitle(
        ";M_{kinfit}(#omega#pi^{0}) - M_{generated}(#omega#pi^{0}) [GeV];Entries"
    );
    h_mass_diff->Draw();
    fit->Draw("same");
    gaus1_func->Draw("same");
    gaus2_func->Draw("same");
    
    // Add fit parameters to plot, colored to match each Gaussian
    TPaveText *pt = new TPaveText(0.58, 0.50, 0.88, 0.89, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.028);
    ((TText*)pt->AddText("Double Gaussian Fit"))->SetTextColor(kRed);
    ((TText*)pt->AddText(TString::Format("Gaus 1: A = %.0f", amp1)))->SetTextColor(kAzure+7);
    ((TText*)pt->AddText(TString::Format("#mu = %.4f GeV", mean1)))->SetTextColor(kAzure+7);
    ((TText*)pt->AddText(TString::Format("#sigma = %.2f GeV", sigma1)))->SetTextColor(kAzure+7);
    ((TText*)pt->AddText(TString::Format("Gaus 2: A = %.0f", amp2)))->SetTextColor(kBlue+2);
    ((TText*)pt->AddText(TString::Format("#mu = %.4f GeV", mean2)))->SetTextColor(kBlue+2);
    ((TText*)pt->AddText(TString::Format("#sigma = %.2f GeV", sigma2)))->SetTextColor(kBlue+2);
    ((TText*)pt->AddText(TString::Format("Avg #sigma = %.2f GeV", sigma_avg)))->SetTextColor(kRed);
    pt->Draw();

    c->SaveAs("mass_resolution.pdf");
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