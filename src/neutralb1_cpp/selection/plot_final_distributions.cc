/**
 * @file plot_final_distributions.cc
 * @author Kevin Scheuer
 * @brief Make plots from finalized amptools trees to show what goes into PWA
 * 
 * These plots do not distinguish between different polarizations or run periods,
 * they are just the final distributions after all cuts and sideband subtractions
 * have been applied.
 * 
 * TODO: have angular calculation be applied to 2D case as well. Don't copy/paste but 
 * actually separate the function
 */

#include <iostream>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TPad.h"

#include "FSBasic/FSTree.h"
#include "FSMode/FSModeTree.h"
#include "FSMode/FSModeHistogram.h"

#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"
#include "/work/halld/kscheuer/my_build/cpu/halld_sim/src/libraries/AMPTOOLS_AMPS/vecPsAngles.h"
#include "/work/halld/kscheuer/my_build/cpu/halld_sim/src/libraries/AMPTOOLS_AMPS/vecPsAngles.cc"

// forward declarations
void plot_mass_spectra(
    TString NT,
    TString CATEGORY,
    TString data_signal,
    TString data_sideband,
    TString mc_signal,
    TString mc_sideband
);
void plot_1d_angles(
    TString NT,
    TString CATEGORY,
    TString data_signal,
    TString data_sideband,
    TString mc_signal,
    TString mc_sideband,
    TString acc_phasespace,
    TString gen_phasespace
);

void plot_final_distributions()
{
    TString tree_dir = "/lustre24/expphy/volatile/halld/home/kscheuer/"
    "FSRoot-skimmed-trees/final-amptools-trees/";

    // input files 
    TString data_signal = tree_dir + "allPeriods_data_signal.root";
    TString data_sideband = tree_dir + "allPeriods_data_background.root";
    TString mc_signal = tree_dir + "PARA_0_allPeriods_ver03.1_mc_signal.root";
    TString mc_sideband = tree_dir + "PARA_0_allPeriods_ver03.1_mc_background.root";
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

    TCanvas *c = new TCanvas("c", "c", 800, 600);

    // draw the histograms
    h_omega_pi0_data_total->Draw("E");
    h_omega_pi0_data_signal->Draw("HIST SAME");
    h_omega_pi0_mc_signal->Draw("HIST SAME");
    h_omega_pi0_data_sideband->Draw("HIST SAME");
    h_omega_pi0_mc_sideband->Draw("HIST SAME");
    legend_omega_pi0->Draw("SAME");

    c->SaveAs("final_omega_pi0_mass_spectrum.pdf");
    c->Clear();

    h_proton_pi0_data_total->Draw("E");
    h_proton_pi0_data_signal->Draw("HIST SAME");
    h_proton_pi0_mc_signal->Draw("HIST SAME");
    h_proton_pi0_data_sideband->Draw("HIST SAME");
    h_proton_pi0_mc_sideband->Draw("HIST SAME");
    legend_proton_pi0->Draw("SAME");

    c->SaveAs("final_proton_pi0_mass_spectrum.pdf");
    c->Clear();

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
    // Open the trees to read 4-momenta and calculate angles
    TFile* f_data_signal = TFile::Open(data_signal);
    TFile* f_data_sideband = TFile::Open(data_sideband);
    TFile* f_mc_signal = TFile::Open(mc_signal);
    TFile* f_mc_sideband = TFile::Open(mc_sideband);
    TFile* f_acc_phasespace = TFile::Open(acc_phasespace);
    TFile* f_gen_phasespace = TFile::Open(gen_phasespace);
    
    TTree* tree_data_signal = (TTree*)f_data_signal->Get(NT);
    TTree* tree_data_sideband = (TTree*)f_data_sideband->Get(NT);
    TTree* tree_mc_signal = (TTree*)f_mc_signal->Get(NT);
    TTree* tree_mc_sideband = (TTree*)f_mc_sideband->Get(NT);
    TTree* tree_acc_phasespace = (TTree*)f_acc_phasespace->Get(NT);
    TTree* tree_gen_phasespace = (TTree*)f_gen_phasespace->Get(NT);
    
    // Create histograms for angles
    // cos(theta) and phi for X resonance decay e.g. b1 -> omega pi0
    TH1F* h_theta_data_signal = new TH1F("h_theta_data_signal", "", 100, -1, 1);
    TH1F* h_theta_data_sideband = new TH1F("h_theta_data_sideband", "", 100, -1, 1);
    TH1F* h_theta_mc_signal = new TH1F("h_theta_mc_signal", "", 100, -1, 1);
    TH1F* h_theta_mc_sideband = new TH1F("h_theta_mc_sideband", "", 100, -1, 1);
    TH1F* h_theta_acc = new TH1F("h_theta_acc", "", 100, -1, 1);
    TH1F* h_theta_gen = new TH1F("h_theta_gen", "", 100, -1, 1);
    
    TH1F* h_phi_data_signal = new TH1F("h_phi_data_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_data_sideband = new TH1F("h_phi_data_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_mc_signal = new TH1F("h_phi_mc_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_mc_sideband = new TH1F("h_phi_mc_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_acc = new TH1F("h_phi_acc", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_gen = new TH1F("h_phi_gen", "", 100, -TMath::Pi(), TMath::Pi());
    
    TH1F* h_Phi_data_signal = new TH1F("h_Phi_data_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_Phi_data_sideband = new TH1F("h_Phi_data_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_Phi_mc_signal = new TH1F("h_Phi_mc_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_Phi_mc_sideband = new TH1F("h_Phi_mc_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_Phi_acc = new TH1F("h_Phi_acc", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_Phi_gen = new TH1F("h_Phi_gen", "", 100, -TMath::Pi(), TMath::Pi());
    
    // cos(theta_h) and phi_h for omega decay (normal to pi+pi- plane in omega rest frame)
    TH1F* h_theta_h_data_signal = new TH1F("h_theta_h_data_signal", "", 100, -1, 1);
    TH1F* h_theta_h_data_sideband = new TH1F("h_theta_h_data_sideband", "", 100, -1, 1);
    TH1F* h_theta_h_mc_signal = new TH1F("h_theta_h_mc_signal", "", 100, -1, 1);
    TH1F* h_theta_h_mc_sideband = new TH1F("h_theta_h_mc_sideband", "", 100, -1, 1);
    TH1F* h_theta_h_acc = new TH1F("h_theta_h_acc", "", 100, -1, 1);
    TH1F* h_theta_h_gen = new TH1F("h_theta_h_gen", "", 100, -1, 1);
    
    TH1F* h_phi_h_data_signal = new TH1F("h_phi_h_data_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_h_data_sideband = new TH1F("h_phi_h_data_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_h_mc_signal = new TH1F("h_phi_h_mc_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_h_mc_sideband = new TH1F("h_phi_h_mc_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_h_acc = new TH1F("h_phi_h_acc", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_h_gen = new TH1F("h_phi_h_gen", "", 100, -TMath::Pi(), TMath::Pi());
    
    // lambda for omega decay
    TH1F* h_lambda_data_signal = new TH1F("h_lambda_data_signal", "", 100, 0, 1);
    TH1F* h_lambda_data_sideband = new TH1F("h_lambda_data_sideband", "", 100, 0, 1);
    TH1F* h_lambda_mc_signal = new TH1F("h_lambda_mc_signal", "", 100, 0, 1);
    TH1F* h_lambda_mc_sideband = new TH1F("h_lambda_mc_sideband", "", 100, 0, 1);
    TH1F* h_lambda_acc = new TH1F("h_lambda_acc", "", 100, 0, 1);
    TH1F* h_lambda_gen = new TH1F("h_lambda_gen", "", 100, 0, 1);
    
    // Helper function to fill angle histograms from a tree
    auto fillAngles = [](
        TTree* tree, 
        TH1F* h_theta, TH1F* h_phi, TH1F* h_Phi, 
        TH1F* h_theta_h, TH1F* h_phi_h, TH1F* h_lambda,
        double weight_scale = 1.0) {
        // Set up branch addresses for 4-momenta
        // AmpTools ordering: B(beam), 1(proton), 2(pi0 bachelor), 3(pi0 omega), 4(pi+ omega), 5(pi- omega)
        double EnPB, PxPB, PyPB, PzPB;
        double EnP1, PxP1, PyP1, PzP1;
        double EnP2, PxP2, PyP2, PzP2;
        double EnP3, PxP3, PyP3, PzP3;
        double EnP4, PxP4, PyP4, PzP4;
        double EnP5, PxP5, PyP5, PzP5;
        double weight = 1.0;
        double polAngle = 0.0;
        
        tree->SetBranchAddress("EnPB", &EnPB);
        tree->SetBranchAddress("PxPB", &PxPB);
        tree->SetBranchAddress("PyPB", &PyPB);
        tree->SetBranchAddress("PzPB", &PzPB);
        
        tree->SetBranchAddress("EnP1", &EnP1);
        tree->SetBranchAddress("PxP1", &PxP1);
        tree->SetBranchAddress("PyP1", &PyP1);
        tree->SetBranchAddress("PzP1", &PzP1);
        
        tree->SetBranchAddress("EnP2", &EnP2);
        tree->SetBranchAddress("PxP2", &PxP2);
        tree->SetBranchAddress("PyP2", &PyP2);
        tree->SetBranchAddress("PzP2", &PzP2);
        
        tree->SetBranchAddress("EnP3", &EnP3);
        tree->SetBranchAddress("PxP3", &PxP3);
        tree->SetBranchAddress("PyP3", &PyP3);
        tree->SetBranchAddress("PzP3", &PzP3);
        
        tree->SetBranchAddress("EnP4", &EnP4);
        tree->SetBranchAddress("PxP4", &PxP4);
        tree->SetBranchAddress("PyP4", &PyP4);
        tree->SetBranchAddress("PzP4", &PzP4);
        
        tree->SetBranchAddress("EnP5", &EnP5);
        tree->SetBranchAddress("PxP5", &PxP5);
        tree->SetBranchAddress("PyP5", &PyP5);
        tree->SetBranchAddress("PzP5", &PzP5);
        
        // Try to get weight and polAngle if they exist
        if (tree->GetBranch("Weight")) {
            tree->SetBranchAddress("Weight", &weight);
        }
        if (tree->GetBranch("PolAngle")) {
            tree->SetBranchAddress("PolAngle", &polAngle);
        }
        
        Long64_t nentries = tree->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree->GetEntry(i);
            
            // Create TLorentzVectors
            TLorentzVector beam(PxPB, PyPB, PzPB, EnPB);
            TLorentzVector proton(PxP1, PyP1, PzP1, EnP1);
            TLorentzVector pi0_bachelor(PxP2, PyP2, PzP2, EnP2);
            TLorentzVector pi0_omega(PxP3, PyP3, PzP3, EnP3);
            TLorentzVector pi_plus(PxP4, PyP4, PzP4, EnP4);
            TLorentzVector pi_minus(PxP5, PyP5, PzP5, EnP5);
            
            // Calculate combined 4-vectors
            TLorentzVector beamP = beam + proton;  // center of mass
            TLorentzVector X = pi0_bachelor + pi0_omega + pi_plus + pi_minus;
            TLorentzVector omega = pi0_omega + pi_plus + pi_minus;
            
            // Get b1 decay angles (theta, phi, Phi)
            vector<double> X_angles = getXDecayAngles(polAngle, beam, beamP, X, omega);
            double theta = X_angles[0];
            double phi = X_angles[1];
            double Phi = X_angles[2];
            
            // Get omega decay angles (theta_h, phi_h, lambda)
            vector<double> omega_angles = getVectorDecayAngles(beamP, X, omega, pi_plus, pi_minus);
            double theta_h = omega_angles[0];
            double phi_h = omega_angles[1];
            double lambda = omega_angles[2];
            
            // Fill histograms with weight (use cos for theta angles)
            double final_weight = weight * weight_scale;
            h_theta->Fill(TMath::Cos(theta), final_weight);
            h_phi->Fill(phi, final_weight);
            h_Phi->Fill(Phi, final_weight);
            h_theta_h->Fill(TMath::Cos(theta_h), final_weight);
            h_phi_h->Fill(phi_h, final_weight);
            h_lambda->Fill(lambda, final_weight);
        }
    };
    
    // Fill histograms
    std::cout << "filling data signal" << "\n";
    fillAngles(tree_data_signal, h_theta_data_signal, h_phi_data_signal, h_Phi_data_signal,
               h_theta_h_data_signal, h_phi_h_data_signal, h_lambda_data_signal);
    std::cout << "filling data sideband" << "\n";
    fillAngles(tree_data_sideband, h_theta_data_sideband, h_phi_data_sideband, h_Phi_data_sideband,
               h_theta_h_data_sideband, h_phi_h_data_sideband, h_lambda_data_sideband);
    std::cout << "filling mc signal" << "\n";
    fillAngles(tree_mc_signal, h_theta_mc_signal, h_phi_mc_signal, h_Phi_mc_signal,
               h_theta_h_mc_signal, h_phi_h_mc_signal, h_lambda_mc_signal);
    std::cout << "filling mc sideband" << "\n";
    fillAngles(tree_mc_sideband, h_theta_mc_sideband, h_phi_mc_sideband, h_Phi_mc_sideband,
               h_theta_h_mc_sideband, h_phi_h_mc_sideband, h_lambda_mc_sideband);
    std::cout << "filling acceptance phase space" << "\n";
    fillAngles(tree_acc_phasespace, h_theta_acc, h_phi_acc, h_Phi_acc,
               h_theta_h_acc, h_phi_h_acc, h_lambda_acc);
    std::cout << "filling generated phase space" << "\n";
    fillAngles(tree_gen_phasespace, h_theta_gen, h_phi_gen, h_Phi_gen,
               h_theta_h_gen, h_phi_h_gen, h_lambda_gen);
    
    // Combine data signal and sideband
    TH1F* h_theta_data_total = (TH1F*)h_theta_data_signal->Clone("h_theta_data_total");
    h_theta_data_total->Add(h_theta_data_sideband);
    
    TH1F* h_phi_data_total = (TH1F*)h_phi_data_signal->Clone("h_phi_data_total");
    h_phi_data_total->Add(h_phi_data_sideband);
    
    TH1F* h_Phi_data_total = (TH1F*)h_Phi_data_signal->Clone("h_Phi_data_total");
    h_Phi_data_total->Add(h_Phi_data_sideband);
    
    TH1F* h_theta_h_data_total = (TH1F*)h_theta_h_data_signal->Clone("h_theta_h_data_total");
    h_theta_h_data_total->Add(h_theta_h_data_sideband);
    
    TH1F* h_phi_h_data_total = (TH1F*)h_phi_h_data_signal->Clone("h_phi_h_data_total");
    h_phi_h_data_total->Add(h_phi_h_data_sideband);
    
    TH1F* h_lambda_data_total = (TH1F*)h_lambda_data_signal->Clone("h_lambda_data_total");
    h_lambda_data_total->Add(h_lambda_data_sideband);
    
    // Calculate acceptance
    TH1F* h_theta_acceptance = (TH1F*)h_theta_acc->Clone("h_theta_acceptance");
    h_theta_acceptance->Divide(h_theta_gen);
    
    TH1F* h_phi_acceptance = (TH1F*)h_phi_acc->Clone("h_phi_acceptance");
    h_phi_acceptance->Divide(h_phi_gen);
    
    TH1F* h_Phi_acceptance = (TH1F*)h_Phi_acc->Clone("h_Phi_acceptance");
    h_Phi_acceptance->Divide(h_Phi_gen);
    
    TH1F* h_theta_h_acceptance = (TH1F*)h_theta_h_acc->Clone("h_theta_h_acceptance");
    h_theta_h_acceptance->Divide(h_theta_h_gen);
    
    TH1F* h_phi_h_acceptance = (TH1F*)h_phi_h_acc->Clone("h_phi_h_acceptance");
    h_phi_h_acceptance->Divide(h_phi_h_gen);
    
    TH1F* h_lambda_acceptance = (TH1F*)h_lambda_acc->Clone("h_lambda_acceptance");
    h_lambda_acceptance->Divide(h_lambda_gen);
    
    // Style histograms - cos(theta)
    double bin_width = get_bin_width(h_theta_data_total);
    h_theta_data_total->SetTitle("");
    h_theta_data_total->SetXTitle("cos(#theta)");
    h_theta_data_total->SetYTitle(TString::Format("Events / %.3f", bin_width));
    h_theta_data_total->SetMarkerStyle(1);
    h_theta_data_total->SetMarkerColor(kBlack);
    
    h_theta_data_signal->SetLineColor(kBlue);
    h_theta_data_signal->SetLineWidth(2);
    
    h_theta_mc_signal->SetLineColor(kCyan+2);
    h_theta_mc_signal->SetLineWidth(1);
    h_theta_mc_signal->SetLineStyle(2);
    
    h_theta_data_sideband->SetLineColor(kRed+1);
    h_theta_data_sideband->SetLineWidth(2);
    
    h_theta_mc_sideband->SetLineColor(kRed-7);
    h_theta_mc_sideband->SetLineWidth(1);
    h_theta_mc_sideband->SetLineStyle(2);
    
    // Style histograms - phi
    bin_width = get_bin_width(h_phi_data_total);
    h_phi_data_total->SetTitle("");
    h_phi_data_total->SetXTitle("#phi (rad)");
    h_phi_data_total->SetYTitle(TString::Format("Events / %.3f rad", bin_width));
    h_phi_data_total->SetMarkerStyle(1);
    h_phi_data_total->SetMarkerColor(kBlack);
    
    h_phi_data_signal->SetLineColor(kBlue);
    h_phi_data_signal->SetLineWidth(2);
    
    h_phi_mc_signal->SetLineColor(kCyan+2);
    h_phi_mc_signal->SetLineWidth(1);
    h_phi_mc_signal->SetLineStyle(2);
    
    h_phi_data_sideband->SetLineColor(kRed+1);
    h_phi_data_sideband->SetLineWidth(2);
    
    h_phi_mc_sideband->SetLineColor(kRed-7);
    h_phi_mc_sideband->SetLineWidth(1);
    h_phi_mc_sideband->SetLineStyle(2);
    
    // Style histograms - Phi
    bin_width = get_bin_width(h_Phi_data_total);
    h_Phi_data_total->SetTitle("");
    h_Phi_data_total->SetXTitle("#Phi (rad)");
    h_Phi_data_total->SetYTitle(TString::Format("Events / %.3f rad", bin_width));
    h_Phi_data_total->SetMarkerStyle(1);
    h_Phi_data_total->SetMarkerColor(kBlack);
    
    h_Phi_data_signal->SetLineColor(kBlue);
    h_Phi_data_signal->SetLineWidth(2);
    
    h_Phi_mc_signal->SetLineColor(kCyan+2);
    h_Phi_mc_signal->SetLineWidth(1);
    h_Phi_mc_signal->SetLineStyle(2);
    
    h_Phi_data_sideband->SetLineColor(kRed+1);
    h_Phi_data_sideband->SetLineWidth(2);
    
    h_Phi_mc_sideband->SetLineColor(kRed-7);
    h_Phi_mc_sideband->SetLineWidth(1);
    h_Phi_mc_sideband->SetLineStyle(2);
    
    // Style histograms - cos(theta_h)
    bin_width = get_bin_width(h_theta_h_data_total);
    h_theta_h_data_total->SetTitle("");
    h_theta_h_data_total->SetXTitle("cos(#theta_{H})");
    h_theta_h_data_total->SetYTitle(TString::Format("Events / %.3f", bin_width));
    h_theta_h_data_total->SetMarkerStyle(1);
    h_theta_h_data_total->SetMarkerColor(kBlack);
    
    h_theta_h_data_signal->SetLineColor(kBlue);
    h_theta_h_data_signal->SetLineWidth(2);
    
    h_theta_h_mc_signal->SetLineColor(kCyan+2);
    h_theta_h_mc_signal->SetLineWidth(1);
    h_theta_h_mc_signal->SetLineStyle(2);
    
    h_theta_h_data_sideband->SetLineColor(kRed+1);
    h_theta_h_data_sideband->SetLineWidth(2);
    
    h_theta_h_mc_sideband->SetLineColor(kRed-7);
    h_theta_h_mc_sideband->SetLineWidth(1);
    h_theta_h_mc_sideband->SetLineStyle(2);
    
    // Style histograms - phi_h
    bin_width = get_bin_width(h_phi_h_data_total);
    h_phi_h_data_total->SetTitle("");
    h_phi_h_data_total->SetXTitle("#phi_{H} (rad)");
    h_phi_h_data_total->SetYTitle(TString::Format("Events / %.3f rad", bin_width));
    h_phi_h_data_total->SetMarkerStyle(1);
    h_phi_h_data_total->SetMarkerColor(kBlack);
    
    h_phi_h_data_signal->SetLineColor(kBlue);
    h_phi_h_data_signal->SetLineWidth(2);
    
    h_phi_h_mc_signal->SetLineColor(kCyan+2);
    h_phi_h_mc_signal->SetLineWidth(1);
    h_phi_h_mc_signal->SetLineStyle(2);
    
    h_phi_h_data_sideband->SetLineColor(kRed+1);
    h_phi_h_data_sideband->SetLineWidth(2);
    
    h_phi_h_mc_sideband->SetLineColor(kRed-7);
    h_phi_h_mc_sideband->SetLineWidth(1);
    h_phi_h_mc_sideband->SetLineStyle(2);
    
    // Style histograms - lambda
    bin_width = get_bin_width(h_lambda_data_total);
    h_lambda_data_total->SetTitle("");
    h_lambda_data_total->SetXTitle("#lambda");
    h_lambda_data_total->SetYTitle(TString::Format("Events / %.3f", bin_width));
    h_lambda_data_total->SetMarkerStyle(1);
    h_lambda_data_total->SetMarkerColor(kBlack);
    
    h_lambda_data_signal->SetLineColor(kBlue);
    h_lambda_data_signal->SetLineWidth(2);
    
    h_lambda_mc_signal->SetLineColor(kCyan+2);
    h_lambda_mc_signal->SetLineWidth(1);
    h_lambda_mc_signal->SetLineStyle(2);
    
    h_lambda_data_sideband->SetLineColor(kRed+1);
    h_lambda_data_sideband->SetLineWidth(2);
    h_lambda_data_sideband->SetLineStyle(2);
    
    // Style acceptance histograms
    h_theta_acceptance->SetLineColor(kGray+2);
    h_theta_acceptance->SetLineWidth(2);
    h_theta_acceptance->SetLineStyle(1);
    
    h_phi_acceptance->SetLineColor(kGray+2);
    h_phi_acceptance->SetLineWidth(2);
    h_phi_acceptance->SetLineStyle(1);
    
    h_Phi_acceptance->SetLineColor(kGray+2);
    h_Phi_acceptance->SetLineWidth(2);
    h_Phi_acceptance->SetLineStyle(1);
    
    h_theta_h_acceptance->SetLineColor(kGray+2);
    h_theta_h_acceptance->SetLineWidth(2);
    h_theta_h_acceptance->SetLineStyle(1);
    
    h_phi_h_acceptance->SetLineColor(kGray+2);
    h_phi_h_acceptance->SetLineWidth(2);
    h_phi_h_acceptance->SetLineStyle(1);

    h_lambda_acceptance->SetLineColor(kGray+2);
    h_lambda_acceptance->SetLineWidth(2);
    h_lambda_acceptance->SetLineStyle(1);
    
    // Create canvas for saving plots
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    
    // Draw theta
    TLegend* legend_theta = new TLegend(0.6, 0.55, 0.89, 0.89);
    legend_theta->AddEntry(h_theta_data_total, "Data", "lp");
    legend_theta->AddEntry(h_theta_data_signal, "Data Signal", "l");
    legend_theta->AddEntry(h_theta_mc_signal, "MC Signal", "l");
    legend_theta->AddEntry(h_theta_data_sideband, "Data Sideband", "l");
    legend_theta->AddEntry(h_theta_mc_sideband, "MC Sideband", "l");
    legend_theta->AddEntry(h_theta_acceptance, "Acceptance", "l");
    h_theta_data_total->Draw("E");
    h_theta_data_signal->Draw("HIST SAME");
    h_theta_mc_signal->Draw("HIST SAME");
    h_theta_data_sideband->Draw("HIST SAME");
    h_theta_mc_sideband->Draw("HIST SAME");
    
    // Draw acceptance on secondary y-axis
    gPad->Update();
    double y_max = gPad->GetUymax();
    double acc_max = h_theta_acceptance->GetMaximum();
    double scale_factor = y_max / (acc_max * 1.1);  // 1.1 for some headroom
    TH1F* h_theta_acc_scaled = (TH1F*)h_theta_acceptance->Clone("h_theta_acc_scaled");
    h_theta_acc_scaled->Scale(scale_factor);
    h_theta_acc_scaled->Draw("HIST SAME");
    
    // Create right axis for acceptance
    TGaxis* axis_theta = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                                     gPad->GetUxmax(), gPad->GetUymax(),
                                     0, acc_max * 1.1, 510, "+L");
    axis_theta->SetLineColor(kGray+2);
    axis_theta->SetLabelColor(kGray+2);
    axis_theta->SetTitleColor(kGray+2);
    axis_theta->SetTitle("Acceptance");
    axis_theta->Draw();
    
    legend_theta->Draw("SAME");
    c1->SaveAs("final_cos_theta.pdf");
    c1->Clear();
    delete legend_theta;
    
    // Draw phi
    TLegend* legend_phi = new TLegend(0.6, 0.55, 0.89, 0.89);
    legend_phi->AddEntry(h_phi_data_total, "Data", "lp");
    legend_phi->AddEntry(h_phi_data_signal, "Data Signal", "l");
    legend_phi->AddEntry(h_phi_mc_signal, "MC Signal", "l");
    legend_phi->AddEntry(h_phi_data_sideband, "Data Sideband", "l");
    legend_phi->AddEntry(h_phi_mc_sideband, "MC Sideband", "l");
    legend_phi->AddEntry(h_phi_acceptance, "Acceptance", "l");
    h_phi_data_total->Draw("E");
    h_phi_data_signal->Draw("HIST SAME");
    h_phi_mc_signal->Draw("HIST SAME");
    h_phi_data_sideband->Draw("HIST SAME");
    h_phi_mc_sideband->Draw("HIST SAME");
    
    // Draw acceptance on secondary y-axis
    gPad->Update();
    y_max = gPad->GetUymax();
    acc_max = h_phi_acceptance->GetMaximum();
    scale_factor = y_max / (acc_max * 1.1);
    TH1F* h_phi_acc_scaled = (TH1F*)h_phi_acceptance->Clone("h_phi_acc_scaled");
    h_phi_acc_scaled->Scale(scale_factor);
    h_phi_acc_scaled->Draw("HIST SAME");
    
    TGaxis* axis_phi = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                                   gPad->GetUxmax(), gPad->GetUymax(),
                                   0, acc_max * 1.1, 510, "+L");
    axis_phi->SetLineColor(kGray+2);
    axis_phi->SetLabelColor(kGray+2);
    axis_phi->SetTitleColor(kGray+2);
    axis_phi->SetTitle("Acceptance");
    axis_phi->Draw();
    
    legend_phi->Draw("SAME");
    c1->SaveAs("final_phi.pdf");
    c1->Clear();
    delete legend_phi;
    
    // Draw Phi
    TLegend* legend_Phi = new TLegend(0.6, 0.55, 0.89, 0.89);
    legend_Phi->AddEntry(h_Phi_data_total, "Data", "lp");
    legend_Phi->AddEntry(h_Phi_data_signal, "Data Signal", "l");
    legend_Phi->AddEntry(h_Phi_mc_signal, "MC Signal", "l");
    legend_Phi->AddEntry(h_Phi_data_sideband, "Data Sideband", "l");
    legend_Phi->AddEntry(h_Phi_mc_sideband, "MC Sideband", "l");
    legend_Phi->AddEntry(h_Phi_acceptance, "Acceptance", "l");
    
    h_Phi_data_total->Draw("E");
    h_Phi_data_signal->Draw("HIST SAME");
    h_Phi_mc_signal->Draw("HIST SAME");
    h_Phi_data_sideband->Draw("HIST SAME");
    h_Phi_mc_sideband->Draw("HIST SAME");
    
    // Draw acceptance on secondary y-axis
    gPad->Update();
    y_max = gPad->GetUymax();
    acc_max = h_Phi_acceptance->GetMaximum();
    scale_factor = y_max / (acc_max * 1.1);
    TH1F* h_Phi_acc_scaled = (TH1F*)h_Phi_acceptance->Clone("h_Phi_acc_scaled");
    h_Phi_acc_scaled->Scale(scale_factor);
    h_Phi_acc_scaled->Draw("HIST SAME");
    
    TGaxis* axis_Phi = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                                   gPad->GetUxmax(), gPad->GetUymax(),
                                   0, acc_max * 1.1, 510, "+L");
    axis_Phi->SetLineColor(kGray+2);
    axis_Phi->SetLabelColor(kGray+2);
    axis_Phi->SetTitleColor(kGray+2);
    axis_Phi->SetTitle("Acceptance");
    axis_Phi->Draw();
    
    legend_Phi->Draw("SAME");
    c1->SaveAs("final_Phi.pdf");
    c1->Clear();
    delete legend_Phi;
    
    // Draw theta_h
    TLegend* legend_theta_h = new TLegend(0.6, 0.55, 0.89, 0.89);
    legend_theta_h->AddEntry(h_theta_h_data_total, "Data", "lp");
    legend_theta_h->AddEntry(h_theta_h_data_signal, "Data Signal", "l");
    legend_theta_h->AddEntry(h_theta_h_mc_signal, "MC Signal", "l");
    legend_theta_h->AddEntry(h_theta_h_data_sideband, "Data Sideband", "l");
    legend_theta_h->AddEntry(h_theta_h_mc_sideband, "MC Sideband", "l");
    legend_theta_h->AddEntry(h_theta_h_acceptance, "Acceptance", "l");
    
    h_theta_h_data_total->Draw("E");
    h_theta_h_data_signal->Draw("HIST SAME");
    h_theta_h_mc_signal->Draw("HIST SAME");
    h_theta_h_data_sideband->Draw("HIST SAME");
    h_theta_h_mc_sideband->Draw("HIST SAME");
    
    // Draw acceptance on secondary y-axis
    gPad->Update();
    y_max = gPad->GetUymax();
    acc_max = h_theta_h_acceptance->GetMaximum();
    scale_factor = y_max / (acc_max * 1.1);
    TH1F* h_theta_h_acc_scaled = (TH1F*)h_theta_h_acceptance->Clone("h_theta_h_acc_scaled");
    h_theta_h_acc_scaled->Scale(scale_factor);
    h_theta_h_acc_scaled->Draw("HIST SAME");
    
    TGaxis* axis_theta_h = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                                       gPad->GetUxmax(), gPad->GetUymax(),
                                       0, acc_max * 1.1, 510, "+L");
    axis_theta_h->SetLineColor(kGray+2);
    axis_theta_h->SetLabelColor(kGray+2);
    axis_theta_h->SetTitleColor(kGray+2);
    axis_theta_h->SetTitle("Acceptance");
    axis_theta_h->Draw();
    
    legend_theta_h->Draw("SAME");
    c1->SaveAs("final_cos_theta_h.pdf");
    c1->Clear();
    delete legend_theta_h;
    
    // Draw phi_h
    TLegend* legend_phi_h = new TLegend(0.6, 0.55, 0.89, 0.89);
    legend_phi_h->AddEntry(h_phi_h_data_total, "Data", "lp");
    legend_phi_h->AddEntry(h_phi_h_data_signal, "Data Signal", "l");
    legend_phi_h->AddEntry(h_phi_h_mc_signal, "MC Signal", "l");
    legend_phi_h->AddEntry(h_phi_h_data_sideband, "Data Sideband", "l");
    legend_phi_h->AddEntry(h_phi_h_mc_sideband, "MC Sideband", "l");
    legend_phi_h->AddEntry(h_phi_h_acceptance, "Acceptance", "l");
    
    h_phi_h_data_total->Draw("E");
    h_phi_h_data_signal->Draw("HIST SAME");
    h_phi_h_mc_signal->Draw("HIST SAME");
    h_phi_h_data_sideband->Draw("HIST SAME");
    h_phi_h_mc_sideband->Draw("HIST SAME");
    
    // Draw acceptance on secondary y-axis
    gPad->Update();
    y_max = gPad->GetUymax();
    acc_max = h_phi_h_acceptance->GetMaximum();
    scale_factor = y_max / (acc_max * 1.1);
    TH1F* h_phi_h_acc_scaled = (TH1F*)h_phi_h_acceptance->Clone("h_phi_h_acc_scaled");
    h_phi_h_acc_scaled->Scale(scale_factor);
    h_phi_h_acc_scaled->Draw("HIST SAME");
    
    TGaxis* axis_phi_h = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                                     gPad->GetUxmax(), gPad->GetUymax(),
                                     0, acc_max * 1.1, 510, "+L");
    axis_phi_h->SetLineColor(kGray+2);
    axis_phi_h->SetLabelColor(kGray+2);
    axis_phi_h->SetTitleColor(kGray+2);
    axis_phi_h->SetTitle("Acceptance");
    axis_phi_h->Draw();
    
    legend_phi_h->Draw("SAME");
    c1->SaveAs("final_phi_h.pdf");
    c1->Clear();
    delete legend_phi_h;
    
    // Draw lambda
    TLegend* legend_lambda = new TLegend(0.6, 0.55, 0.89, 0.89);
    legend_lambda->AddEntry(h_lambda_data_total, "Data", "lp");
    legend_lambda->AddEntry(h_lambda_data_signal, "Data Signal", "l");
    legend_lambda->AddEntry(h_lambda_mc_signal, "MC Signal", "l");
    legend_lambda->AddEntry(h_lambda_data_sideband, "Data Sideband", "l");
    legend_lambda->AddEntry(h_lambda_mc_sideband, "MC Sideband", "l");
    legend_lambda->AddEntry(h_lambda_acceptance, "Acceptance", "l");
    
    h_lambda_data_total->Draw("E");
    h_lambda_data_signal->Draw("HIST SAME");
    h_lambda_mc_signal->Draw("HIST SAME");
    h_lambda_data_sideband->Draw("HIST SAME");
    h_lambda_mc_sideband->Draw("HIST SAME");
    
    // Draw acceptance on secondary y-axis
    gPad->Update();
    y_max = gPad->GetUymax();
    acc_max = h_lambda_acceptance->GetMaximum();
    scale_factor = y_max / (acc_max * 1.1);
    TH1F* h_lambda_acc_scaled = (TH1F*)h_lambda_acceptance->Clone("h_lambda_acc_scaled");
    h_lambda_acc_scaled->Scale(scale_factor);
    h_lambda_acc_scaled->Draw("HIST SAME");
    
    TGaxis* axis_lambda = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                                      gPad->GetUxmax(), gPad->GetUymax(),
                                      0, acc_max * 1.1, 510, "+L");
    axis_lambda->SetLineColor(kGray+2);
    axis_lambda->SetLabelColor(kGray+2);
    axis_lambda->SetTitleColor(kGray+2);
    axis_lambda->SetTitle("Acceptance");
    axis_lambda->Draw();
    
    legend_lambda->Draw("SAME");
    c1->SaveAs("final_lambda.pdf");
    c1->Clear();
    delete legend_lambda;
    delete c1;
    
    // Clean up
    f_data_signal->Close();
    f_data_sideband->Close();
    f_mc_signal->Close();
    f_mc_sideband->Close();
    f_acc_phasespace->Close();
    f_gen_phasespace->Close();
    
    return;
}