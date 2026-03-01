/**
 * @file omega_resolution.cc
 * @author Kevin Scheuer
 * @brief Quick script to plot measured vs kinfit omega mass distribution
 * 
 * Reminder that this uses the best chi2 trees, so we plot the both permutations
 * together of the pi0s.
 */


#include <iostream>

#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TLegend.h"

#include "FSBasic/FSHistogram.h"
#include "FSBasic/FSCut.h"
#include "FSBasic/FSTree.h"
#include "FSMode/FSModeHistogram.h"

#include "neutralb1/fit_utils.h"
#include "../fsroot_setup.cc"

void omega_resolution()
{
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(false);

    TString input_files = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/best-chi2/"
        "tree_pi0pi0pippim__B4_bestChi2_SKIM_03_data.root";

    TH1F *h_measured_mass_perm1 = FSModeHistogram::getTH1F(
        input_files,
        NT,
        CATEGORY,
        "RMASS(2,3,4)",
        "(100, 0.5, 1.1)",
        ""
    );
    TH1F *h_measured_mass_perm2 = FSModeHistogram::getTH1F(
        input_files,
        NT,
        CATEGORY,
        "RMASS(2,3,5)",
        "(100, 0.5, 1.1)",
        ""
    );

    TH1F *h_kinfit_mass_perm1 = FSModeHistogram::getTH1F(
        input_files,
        NT,
        CATEGORY,
        "MASS(2,3,4)",
        "(100, 0.5, 1.1)",
        ""
    );
     TH1F *h_kinfit_mass_perm2 = FSModeHistogram::getTH1F(
        input_files,
        NT,
        CATEGORY,
        "MASS(2,3,5)",
        "(100, 0.5, 1.1)",
        ""
    );

    h_measured_mass_perm1->Add(h_measured_mass_perm2);
    h_kinfit_mass_perm1->Add(h_kinfit_mass_perm2);

    TCanvas *c = new TCanvas("c_omega_mass", "c_omega_mass", 800, 600);
    h_kinfit_mass_perm1->SetLineColor(kBlue);
    h_kinfit_mass_perm1->SetLineWidth(2);
    h_kinfit_mass_perm1->Draw("HIST");

    h_kinfit_mass_perm1->SetXTitle("#pi^{#plus}#pi^{#minus}#pi^{0}_{i} inv. mass (GeV)");
    h_kinfit_mass_perm1->SetYTitle("Events / 6 MeV");
    h_kinfit_mass_perm1->SetTitle("");

    h_measured_mass_perm1->SetLineColor(kRed);
    h_measured_mass_perm1->SetLineWidth(2);
    h_measured_mass_perm1->Draw("HIST SAME");
    
    TLegend *legend = new TLegend(0.55, 0.55, 0.88, 0.88);
    legend->SetBorderSize(1);
    legend->SetFillStyle(0);
    legend->AddEntry(h_measured_mass_perm1, "Measured", "l");
    legend->AddEntry(h_kinfit_mass_perm1, "Kinematic Fit", "l");
    legend->Draw();
    
    c->SaveAs("omega_measured_mass_comparison.pdf");
    
    return;
}