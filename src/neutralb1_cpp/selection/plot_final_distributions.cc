/**
 * @file plot_final_distributions.cc
 * @author Kevin Scheuer
 * @brief Make plots from finalized amptools trees to show what goes into PWA
 * 
 * These plots do not distinguish between different polarizations or run periods,
 * they are just the final distributions after all cuts and sideband subtractions
 * have been applied.
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
#include "/work/halld/kscheuer/my_build/cpu/halld_sim/src/libraries/AMPTOOLS_AMPS/vecPsAngles.h"
#include "/work/halld/kscheuer/my_build/cpu/halld_sim/src/libraries/AMPTOOLS_AMPS/vecPsAngles.cc"

// Structure to hold calculated angles
struct AngleData {
    std::vector<double> theta;
    std::vector<double> phi;
    std::vector<double> Phi;
    std::vector<double> theta_h;
    std::vector<double> phi_h;
    std::vector<double> lambda;
    std::vector<double> weights;
};

// forward declarations
void plot_mass_spectra(
    TString NT,
    TString CATEGORY,
    TString data_signal,
    TString data_sideband,
    TString mc_signal,
    TString mc_sideband
);
AngleData calculate_angles(
    TTree* tree,
    double weight_scale = 1.0
);
void plot_1d_angles(
    TString NT,
    TString CATEGORY,
    TString data_signal,
    TString data_sideband,
    TString mc_signal,
    TString mc_sideband,
    TString acc_phasespace,
    TString gen_phasespace,
    AngleData& angles_data_signal,
    AngleData& angles_data_sideband,
    AngleData& angles_mc_signal,
    AngleData& angles_mc_sideband,
    AngleData& angles_acc,
    AngleData& angles_gen
);
void plot_2d_angles(
    TString NT,
    TString CATEGORY,
    const AngleData& angles_signal,
    const AngleData& angles_sideband,
    TString output_prefix
);
void plot_2d_acceptance(
    TString NT,
    TString CATEGORY,
    const AngleData& angles_acc,
    const AngleData& angles_gen
);

void plot_final_distributions()
{
    gStyle->SetOptStat(0);
    TString tree_dir = "/lustre24/expphy/volatile/halld/home/kscheuer/"
    "FSRoot-skimmed-trees/final-amptools-trees/";

    // input files 
    TString data_signal = tree_dir + "allPeriods_data_signal.root";
    TString data_sideband = tree_dir + "allPeriods_data_background.root";
    TString mc_signal = tree_dir + "PARA_0_allPeriods_ver03.1_mc_signal.root";
    TString mc_sideband = tree_dir + "PARA_0_allPeriods_ver03.1_mc_background.root";
    TString acc_phasespace = tree_dir + "allPeriods_ver03_phasespace.root";
    TString gen_phasespace = tree_dir + "allPeriods_ver03_gen_phasespace.root";

    TString NT = "ntFSGlueX_100_112";
    TString CATEGORY = "pi0pi0pippim";    

    plot_mass_spectra(NT, CATEGORY, data_signal, data_sideband, mc_signal, mc_sideband);

    // Calculate angles once for all datasets
    std::cout << "Calculating angles for all datasets..." << std::endl;
    
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
    
    std::cout << "Calculating data signal angles..." << std::endl;
    AngleData angles_data_signal = calculate_angles(tree_data_signal);
    std::cout << "Calculating data sideband angles..." << std::endl;
    AngleData angles_data_sideband = calculate_angles(tree_data_sideband);
    std::cout << "Calculating MC signal angles..." << std::endl;
    AngleData angles_mc_signal = calculate_angles(tree_mc_signal);
    std::cout << "Calculating MC sideband angles..." << std::endl;
    AngleData angles_mc_sideband = calculate_angles(tree_mc_sideband);
    std::cout << "Calculating acceptance phase space angles..." << std::endl;
    AngleData angles_acc = calculate_angles(tree_acc_phasespace);
    std::cout << "Calculating generated phase space angles..." << std::endl;
    //AngleData angles_gen = calculate_angles(tree_gen_phasespace);
    AngleData angles_gen;

    // Close files
    f_data_signal->Close();
    f_data_sideband->Close();
    f_mc_signal->Close();
    f_mc_sideband->Close();
    f_acc_phasespace->Close();
    f_gen_phasespace->Close();
    
    // Now create plots using pre-calculated angles
    plot_1d_angles(
        NT,
        CATEGORY,
        data_signal,
        data_sideband,
        mc_signal,
        mc_sideband,
        acc_phasespace,
        gen_phasespace,
        angles_data_signal,
        angles_data_sideband,
        angles_mc_signal,
        angles_mc_sideband,
        angles_acc,
        angles_gen
    );
    
    plot_2d_angles(
        NT,
        CATEGORY,
        angles_data_signal,
        angles_data_sideband,
        "data"
    );
    
    plot_2d_angles(
        NT,
        CATEGORY,
        angles_mc_signal,
        angles_mc_sideband,
        "mc"
    );
    
    // plot_2d_acceptance(
    //     NT,
    //     CATEGORY,
    //     angles_acc,
    //     angles_gen
    // );
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
    // get files and trees
    TFile *f_data_signal = TFile::Open(data_signal);
    TFile *f_data_sideband = TFile::Open(data_sideband);
    TFile *f_mc_signal = TFile::Open(mc_signal);
    TFile *f_mc_sideband = TFile::Open(mc_sideband);

    TTree *tree_data_signal = (TTree*)f_data_signal->Get(NT);
    TTree *tree_data_sideband = (TTree*)f_data_sideband->Get(NT);
    TTree *tree_mc_signal = (TTree*)f_mc_signal->Get(NT);
    TTree *tree_mc_sideband = (TTree*)f_mc_sideband->Get(NT);

    // setup histograms for omega pi0 and proton pi0
    TH1F *h_omega_pi0_data_signal = new TH1F(
        "h_omega_pi0_data_signal", 
        "", 
        100, 1.0, 2.0);
    TH1F *h_omega_pi0_data_sideband = new TH1F(
        "h_omega_pi0_data_sideband", 
        "", 
        100, 1.0, 2.0);
    TH1F *h_omega_pi0_mc_signal = new TH1F(
        "h_omega_pi0_mc_signal", 
        "", 
        100, 1.0, 2.0);
    TH1F *h_omega_pi0_mc_sideband = new TH1F(
        "h_omega_pi0_mc_sideband", 
        "", 
        100, 1.0, 2.0);

    TH1F *h_proton_pi0_data_signal = new TH1F(
        "h_proton_pi0_data_signal", 
        "", 
        200, 1.0, 3.0);
    TH1F *h_proton_pi0_data_sideband = new TH1F(
        "h_proton_pi0_data_sideband", 
        "", 
        200, 1.0, 3.0);
    TH1F *h_proton_pi0_mc_signal = new TH1F(
        "h_proton_pi0_mc_signal", 
        "", 
        200, 1.0, 3.0);
    TH1F *h_proton_pi0_mc_sideband = new TH1F(
        "h_proton_pi0_mc_sideband", 
        "", 
        200, 1.0, 3.0);

    // fill histograms
    tree_data_signal->Draw("M4Pi>>h_omega_pi0_data_signal", "", "goff");
    tree_data_sideband->Draw("M4Pi>>h_omega_pi0_data_sideband", "weight", "goff");
    tree_mc_signal->Draw("M4Pi>>h_omega_pi0_mc_signal", "", "goff");
    tree_mc_sideband->Draw("M4Pi>>h_omega_pi0_mc_sideband", "weight", "goff");
    
    tree_data_signal->Draw("MRecoilPi>>h_proton_pi0_data_signal", "", "goff");
    tree_data_sideband->Draw("MRecoilPi>>h_proton_pi0_data_sideband", "weight", "goff");
    tree_mc_signal->Draw("MRecoilPi>>h_proton_pi0_mc_signal", "", "goff");
    tree_mc_sideband->Draw("MRecoilPi>>h_proton_pi0_mc_sideband", "weight", "goff");
    
    // get final data as signal minus sideband
    TH1F *h_omega_pi0_data_total = 
        (TH1F*)h_omega_pi0_data_signal->Clone("h_omega_pi0_data_total");
    h_omega_pi0_data_total->Add(h_omega_pi0_data_sideband, -1.0);
    TH1F *h_proton_pi0_data_total = 
        (TH1F*)h_proton_pi0_data_signal->Clone("h_proton_pi0_data_total");
    h_proton_pi0_data_total->Add(h_proton_pi0_data_sideband, -1.0);

    // customize the histograms
    double bin_width = get_bin_width( h_omega_pi0_data_total );
    h_omega_pi0_data_total->SetTitle("");
    h_omega_pi0_data_total->SetXTitle("#omega#pi^{0} inv. mass (GeV)");
    h_omega_pi0_data_total->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    h_omega_pi0_data_total->SetMarkerStyle(1);
    h_omega_pi0_data_total->SetMarkerColor(kBlack);
    h_omega_pi0_data_total->SetLineColor(kBlack);
    h_omega_pi0_data_total->SetMinimum(0);
    h_omega_pi0_data_total->SetMaximum(h_omega_pi0_data_signal->GetMaximum() * 1.1);

    h_omega_pi0_data_signal->SetLineColor(kBlue);
    h_omega_pi0_data_signal->SetLineWidth(2);

    h_omega_pi0_mc_signal->SetMarkerColor(kCyan+2);
    h_omega_pi0_mc_signal->SetMarkerStyle(24);
    h_omega_pi0_mc_signal->SetMarkerSize(0.8);
    h_omega_pi0_mc_signal->SetLineColor(kCyan+2);

    h_omega_pi0_data_sideband->SetLineColor(kRed+1);
    h_omega_pi0_data_sideband->SetLineWidth(2);

    h_omega_pi0_mc_sideband->SetMarkerColor(kRed-7);
    h_omega_pi0_mc_sideband->SetMarkerStyle(25);
    h_omega_pi0_mc_sideband->SetMarkerSize(0.8);
    h_omega_pi0_mc_sideband->SetLineColor(kRed-7);


    bin_width = get_bin_width( h_proton_pi0_data_total );
    h_proton_pi0_data_total->SetTitle("");
    h_proton_pi0_data_total->SetXTitle("p'#pi^{0} inv. mass (GeV)");
    h_proton_pi0_data_total->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
    h_proton_pi0_data_total->SetMarkerStyle(1);
    h_proton_pi0_data_total->SetMarkerColor(kBlack);
    h_proton_pi0_data_total->SetLineColor(kBlack);
    h_proton_pi0_data_total->SetMinimum(0);
    h_proton_pi0_data_total->SetMaximum(h_proton_pi0_data_signal->GetMaximum() * 1.1);

    h_proton_pi0_data_signal->SetLineColor(kBlue);
    h_proton_pi0_data_signal->SetLineWidth(2);

    h_proton_pi0_mc_signal->SetMarkerColor(kCyan+2);
    h_proton_pi0_mc_signal->SetMarkerStyle(24);
    h_proton_pi0_mc_signal->SetMarkerSize(0.8);
    h_proton_pi0_mc_signal->SetLineColor(kCyan+2);

    h_proton_pi0_data_sideband->SetLineColor(kRed+1);
    h_proton_pi0_data_sideband->SetLineWidth(2);

    h_proton_pi0_mc_sideband->SetMarkerColor(kRed-7);
    h_proton_pi0_mc_sideband->SetMarkerStyle(25);
    h_proton_pi0_mc_sideband->SetMarkerSize(0.8);
    h_proton_pi0_mc_sideband->SetLineColor(kRed-7);

    // Scale MC to match data maximum
    double scale_omega_signal = h_omega_pi0_data_signal->GetMaximum() / h_omega_pi0_mc_signal->GetMaximum();
    h_omega_pi0_mc_signal->Scale(scale_omega_signal);
    
    double scale_omega_sideband = h_omega_pi0_data_sideband->GetMaximum() / h_omega_pi0_mc_sideband->GetMaximum();
    h_omega_pi0_mc_sideband->Scale(scale_omega_sideband);
    
    double scale_proton_signal = h_proton_pi0_data_signal->GetMaximum() / h_proton_pi0_mc_signal->GetMaximum();
    h_proton_pi0_mc_signal->Scale(scale_proton_signal);
    
    double scale_proton_sideband = h_proton_pi0_data_sideband->GetMaximum() / h_proton_pi0_mc_sideband->GetMaximum();
    h_proton_pi0_mc_sideband->Scale(scale_proton_sideband);

    // create legend
    TLegend* legend_omega_pi0 = new TLegend(0.7,0.7,0.88,0.88);
    legend_omega_pi0->AddEntry(h_omega_pi0_data_total, "Final Data", "lp");
    legend_omega_pi0->AddEntry(h_omega_pi0_data_signal, "Data Signal", "l");
    legend_omega_pi0->AddEntry(h_omega_pi0_mc_signal, TString::Format("MC Signal (scale=%.3f)", scale_omega_signal), "l");
    legend_omega_pi0->AddEntry(h_omega_pi0_data_sideband, "Data Sideband", "l");
    legend_omega_pi0->AddEntry(h_omega_pi0_mc_sideband, TString::Format("MC Sideband (scale=%.3f)", scale_omega_sideband), "l");

    TLegend* legend_proton_pi0 = new TLegend(0.15,0.7,0.33,0.88);
    legend_proton_pi0->AddEntry(h_proton_pi0_data_total, "Final Data", "lp");
    legend_proton_pi0->AddEntry(h_proton_pi0_data_signal, "Data Signal", "l");
    legend_proton_pi0->AddEntry(h_proton_pi0_mc_signal, TString::Format("MC Signal (scale=%.3f)", scale_proton_signal), "l");
    legend_proton_pi0->AddEntry(h_proton_pi0_data_sideband, "Data Sideband", "l");
    legend_proton_pi0->AddEntry(h_proton_pi0_mc_sideband, TString::Format("MC Sideband (scale=%.3f)", scale_proton_sideband), "l");

    TCanvas *c = new TCanvas("c", "c", 800, 600);

    // draw the histograms
    h_omega_pi0_data_total->Draw("E");
    h_omega_pi0_data_signal->Draw("HIST SAME");
    h_omega_pi0_mc_signal->Draw("E SAME");
    h_omega_pi0_data_sideband->Draw("HIST SAME");
    h_omega_pi0_mc_sideband->Draw("E SAME");
    legend_omega_pi0->Draw("SAME");

    c->SaveAs("final_omega_pi0_mass_spectrum.pdf");
    c->Clear();

    h_proton_pi0_data_total->Draw("E");
    h_proton_pi0_data_signal->Draw("HIST SAME");
    h_proton_pi0_mc_signal->Draw("E SAME");
    h_proton_pi0_data_sideband->Draw("HIST SAME");
    h_proton_pi0_mc_sideband->Draw("E SAME");
    legend_proton_pi0->Draw("SAME");

    c->SaveAs("final_proton_pi0_mass_spectrum.pdf");
    c->Clear();

    return;
}

/**
 * @brief Calculate all angles from a tree and store them in an AngleData structure
 * 
 * @param tree TTree containing 4-momenta data
 * @param weight_scale Optional weight scaling factor
 * @return AngleData structure containing all calculated angles and weights
 */
AngleData calculate_angles(
    TTree* tree,
    double weight_scale
)
{
    AngleData angles;
    
    // Set up branch addresses for 4-momenta
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
    if (tree->GetBranch("weight")) {
        tree->SetBranchAddress("weight", &weight);
    }
    if (tree->GetBranch("PolAngle")) {
        tree->SetBranchAddress("PolAngle", &polAngle);
    }
    
    Long64_t nentries = tree->GetEntries();
    angles.theta.reserve(nentries);
    angles.phi.reserve(nentries);
    angles.Phi.reserve(nentries);
    angles.theta_h.reserve(nentries);
    angles.phi_h.reserve(nentries);
    angles.lambda.reserve(nentries);
    angles.weights.reserve(nentries);
    
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
        angles.theta.push_back(X_angles[0]);
        angles.phi.push_back(X_angles[1]);
        angles.Phi.push_back(X_angles[2]);
        
        // Get omega decay angles (theta_h, phi_h, lambda)
        vector<double> omega_angles = getVectorDecayAngles(beamP, X, omega, pi_plus, pi_minus);
        angles.theta_h.push_back(omega_angles[0]);
        angles.phi_h.push_back(omega_angles[1]);
        angles.lambda.push_back(omega_angles[2]);
        
        angles.weights.push_back(weight * weight_scale);
    }
    
    return angles;
}

void plot_1d_angles(
    TString NT,
    TString CATEGORY,
    TString data_signal,
    TString data_sideband,
    TString mc_signal,
    TString mc_sideband,
    TString acc_phasespace,
    TString gen_phasespace,
    AngleData& angles_data_signal,
    AngleData& angles_data_sideband,
    AngleData& angles_mc_signal,
    AngleData& angles_mc_sideband,
    AngleData& angles_acc,
    AngleData& angles_gen
)
{
    // Create histograms for angles
    // cos(theta) and phi for X resonance decay e.g. b1 -> omega pi0
    TH1F* h_theta_data_signal = new TH1F("h_theta_data_signal", "", 100, -1, 1);
    TH1F* h_theta_data_sideband = new TH1F("h_theta_data_sideband", "", 100, -1, 1);
    TH1F* h_theta_mc_signal = new TH1F("h_theta_mc_signal", "", 100, -1, 1);
    TH1F* h_theta_mc_sideband = new TH1F("h_theta_mc_sideband", "", 100, -1, 1);
    TH1F* h_theta_acc = new TH1F("h_theta_acc", "", 100, -1, 1);
    // TH1F* h_theta_gen = new TH1F("h_theta_gen", "", 100, -1, 1);
    
    TH1F* h_phi_data_signal = new TH1F("h_phi_data_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_data_sideband = new TH1F("h_phi_data_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_mc_signal = new TH1F("h_phi_mc_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_mc_sideband = new TH1F("h_phi_mc_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_acc = new TH1F("h_phi_acc", "", 100, -TMath::Pi(), TMath::Pi());
    // TH1F* h_phi_gen = new TH1F("h_phi_gen", "", 100, -TMath::Pi(), TMath::Pi());
    
    TH1F* h_Phi_data_signal = new TH1F("h_Phi_data_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_Phi_data_sideband = new TH1F("h_Phi_data_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_Phi_mc_signal = new TH1F("h_Phi_mc_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_Phi_mc_sideband = new TH1F("h_Phi_mc_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_Phi_acc = new TH1F("h_Phi_acc", "", 100, -TMath::Pi(), TMath::Pi());
    // TH1F* h_Phi_gen = new TH1F("h_Phi_gen", "", 100, -TMath::Pi(), TMath::Pi());
    
    // cos(theta_h) and phi_h for omega decay (normal to pi+pi- plane in omega rest frame)
    TH1F* h_theta_h_data_signal = new TH1F("h_theta_h_data_signal", "", 100, -1, 1);
    TH1F* h_theta_h_data_sideband = new TH1F("h_theta_h_data_sideband", "", 100, -1, 1);
    TH1F* h_theta_h_mc_signal = new TH1F("h_theta_h_mc_signal", "", 100, -1, 1);
    TH1F* h_theta_h_mc_sideband = new TH1F("h_theta_h_mc_sideband", "", 100, -1, 1);
    TH1F* h_theta_h_acc = new TH1F("h_theta_h_acc", "", 100, -1, 1);
    // TH1F* h_theta_h_gen = new TH1F("h_theta_h_gen", "", 100, -1, 1);
    
    TH1F* h_phi_h_data_signal = new TH1F("h_phi_h_data_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_h_data_sideband = new TH1F("h_phi_h_data_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_h_mc_signal = new TH1F("h_phi_h_mc_signal", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_h_mc_sideband = new TH1F("h_phi_h_mc_sideband", "", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_phi_h_acc = new TH1F("h_phi_h_acc", "", 100, -TMath::Pi(), TMath::Pi());
    // TH1F* h_phi_h_gen = new TH1F("h_phi_h_gen", "", 100, -TMath::Pi(), TMath::Pi());
    
    // lambda for omega decay
    TH1F* h_lambda_data_signal = new TH1F("h_lambda_data_signal", "", 100, 0, 1);
    TH1F* h_lambda_data_sideband = new TH1F("h_lambda_data_sideband", "", 100, 0, 1);
    TH1F* h_lambda_mc_signal = new TH1F("h_lambda_mc_signal", "", 100, 0, 1);
    TH1F* h_lambda_mc_sideband = new TH1F("h_lambda_mc_sideband", "", 100, 0, 1);
    TH1F* h_lambda_acc = new TH1F("h_lambda_acc", "", 100, 0, 1);
    // TH1F* h_lambda_gen = new TH1F("h_lambda_gen", "", 100, 0, 1);
    
    // Helper function to fill histograms from pre-calculated AngleData
    auto fillHistograms = [](
        const AngleData& angles,
        TH1F* h_theta, TH1F* h_phi, TH1F* h_Phi, 
        TH1F* h_theta_h, TH1F* h_phi_h, TH1F* h_lambda) {
        
        size_t nentries = angles.theta.size();
        for (size_t i = 0; i < nentries; i++) {
            double weight = angles.weights[i];
            // Fill histograms with weight (use cos for theta angles)
            h_theta->Fill(TMath::Cos(angles.theta[i]), weight);
            h_phi->Fill(angles.phi[i], weight);
            h_Phi->Fill(angles.Phi[i], weight);
            h_theta_h->Fill(TMath::Cos(angles.theta_h[i]), weight);
            h_phi_h->Fill(angles.phi_h[i], weight);
            h_lambda->Fill(angles.lambda[i], weight);
        }
    };
    
    // Fill histograms from pre-calculated angles
    std::cout << "Filling 1D histograms from pre-calculated angles..." << std::endl;
    fillHistograms(angles_data_signal, h_theta_data_signal, h_phi_data_signal, h_Phi_data_signal,
                   h_theta_h_data_signal, h_phi_h_data_signal, h_lambda_data_signal);
    fillHistograms(angles_data_sideband, h_theta_data_sideband, h_phi_data_sideband, h_Phi_data_sideband,
                   h_theta_h_data_sideband, h_phi_h_data_sideband, h_lambda_data_sideband);
    fillHistograms(angles_mc_signal, h_theta_mc_signal, h_phi_mc_signal, h_Phi_mc_signal,
                   h_theta_h_mc_signal, h_phi_h_mc_signal, h_lambda_mc_signal);
    fillHistograms(angles_mc_sideband, h_theta_mc_sideband, h_phi_mc_sideband, h_Phi_mc_sideband,
                   h_theta_h_mc_sideband, h_phi_h_mc_sideband, h_lambda_mc_sideband);
    fillHistograms(angles_acc, h_theta_acc, h_phi_acc, h_Phi_acc,
                   h_theta_h_acc, h_phi_h_acc, h_lambda_acc);
    // fillHistograms(angles_gen, h_theta_gen, h_phi_gen, h_Phi_gen,
    //                h_theta_h_gen, h_phi_h_gen, h_lambda_gen);
    
    // Combine data signal and sideband (signal - sideband)
    TH1F* h_theta_data_total = (TH1F*)h_theta_data_signal->Clone("h_theta_data_total");
    h_theta_data_total->Add(h_theta_data_sideband, -1.0);
    
    TH1F* h_phi_data_total = (TH1F*)h_phi_data_signal->Clone("h_phi_data_total");
    h_phi_data_total->Add(h_phi_data_sideband, -1.0);
    
    TH1F* h_Phi_data_total = (TH1F*)h_Phi_data_signal->Clone("h_Phi_data_total");
    h_Phi_data_total->Add(h_Phi_data_sideband, -1.0);
    
    TH1F* h_theta_h_data_total = (TH1F*)h_theta_h_data_signal->Clone("h_theta_h_data_total");
    h_theta_h_data_total->Add(h_theta_h_data_sideband, -1.0);
    
    TH1F* h_phi_h_data_total = (TH1F*)h_phi_h_data_signal->Clone("h_phi_h_data_total");
    h_phi_h_data_total->Add(h_phi_h_data_sideband, -1.0);
    
    TH1F* h_lambda_data_total = (TH1F*)h_lambda_data_signal->Clone("h_lambda_data_total");
    h_lambda_data_total->Add(h_lambda_data_sideband, -1.0);
    
    // // Calculate acceptance
    // TH1F* h_theta_acceptance = (TH1F*)h_theta_acc->Clone("h_theta_acceptance");
    // h_theta_acceptance->Divide(h_theta_gen);
    
    // TH1F* h_phi_acceptance = (TH1F*)h_phi_acc->Clone("h_phi_acceptance");
    // h_phi_acceptance->Divide(h_phi_gen);
    
    // TH1F* h_Phi_acceptance = (TH1F*)h_Phi_acc->Clone("h_Phi_acceptance");
    // h_Phi_acceptance->Divide(h_Phi_gen);
    
    // TH1F* h_theta_h_acceptance = (TH1F*)h_theta_h_acc->Clone("h_theta_h_acceptance");
    // h_theta_h_acceptance->Divide(h_theta_h_gen);
    
    // TH1F* h_phi_h_acceptance = (TH1F*)h_phi_h_acc->Clone("h_phi_h_acceptance");
    // h_phi_h_acceptance->Divide(h_phi_h_gen);
    
    // TH1F* h_lambda_acceptance = (TH1F*)h_lambda_acc->Clone("h_lambda_acceptance");
    // h_lambda_acceptance->Divide(h_lambda_gen);
    
    // Style histograms - cos(theta)
    double bin_width = get_bin_width(h_theta_data_total);
    h_theta_data_total->SetTitle("");
    h_theta_data_total->SetXTitle("cos(#theta)");
    h_theta_data_total->SetYTitle(TString::Format("Events / %.3f", bin_width));
    h_theta_data_total->SetMarkerStyle(1);
    h_theta_data_total->SetMarkerColor(kBlack);
    h_theta_data_total->SetLineColor(kBlack);
    h_theta_data_total->SetMinimum(0);
    h_theta_data_total->SetMaximum(h_theta_data_signal->GetMaximum() * 1.1);
    
    h_theta_data_signal->SetLineColor(kBlue);
    h_theta_data_signal->SetLineWidth(2);
    
    h_theta_mc_signal->SetMarkerColor(kCyan+2);
    h_theta_mc_signal->SetMarkerStyle(24);
    h_theta_mc_signal->SetMarkerSize(0.8);
    h_theta_mc_signal->SetLineColor(kCyan+2);
    
    h_theta_data_sideband->SetLineColor(kRed+1);
    h_theta_data_sideband->SetLineWidth(2);
    
    h_theta_mc_sideband->SetMarkerColor(kRed-7);
    h_theta_mc_sideband->SetMarkerStyle(25);
    h_theta_mc_sideband->SetMarkerSize(0.8);
    h_theta_mc_sideband->SetLineColor(kRed-7);
    
    // Scale MC to match data maximum
    double scale_theta_signal = h_theta_data_signal->GetMaximum() / h_theta_mc_signal->GetMaximum();
    h_theta_mc_signal->Scale(scale_theta_signal);
    
    double scale_theta_sideband = h_theta_data_sideband->GetMaximum() / h_theta_mc_sideband->GetMaximum();
    h_theta_mc_sideband->Scale(scale_theta_sideband);
    
    // Style histograms - phi
    bin_width = get_bin_width(h_phi_data_total);
    h_phi_data_total->SetTitle("");
    h_phi_data_total->SetXTitle("#phi (rad)");
    h_phi_data_total->SetYTitle(TString::Format("Events / %.3f rad", bin_width));
    h_phi_data_total->SetMarkerStyle(1);
    h_phi_data_total->SetMarkerColor(kBlack);
    h_phi_data_total->SetLineColor(kBlack);
    h_phi_data_total->SetMinimum(0);
    h_phi_data_total->SetMaximum(h_phi_data_signal->GetMaximum() * 1.1);
    
    h_phi_data_signal->SetLineColor(kBlue);
    h_phi_data_signal->SetLineWidth(2);
    
    h_phi_mc_signal->SetMarkerColor(kCyan+2);
    h_phi_mc_signal->SetMarkerStyle(24);
    h_phi_mc_signal->SetMarkerSize(0.8);
    h_phi_mc_signal->SetLineColor(kCyan+2);
    
    h_phi_data_sideband->SetLineColor(kRed+1);
    h_phi_data_sideband->SetLineWidth(2);
    
    h_phi_mc_sideband->SetMarkerColor(kRed-7);
    h_phi_mc_sideband->SetMarkerStyle(25);
    h_phi_mc_sideband->SetMarkerSize(0.8);
    h_phi_mc_sideband->SetLineColor(kRed-7);
    
    // Scale MC to match data maximum
    double scale_phi_signal = h_phi_data_signal->GetMaximum() / h_phi_mc_signal->GetMaximum();
    h_phi_mc_signal->Scale(scale_phi_signal);
    
    double scale_phi_sideband = h_phi_data_sideband->GetMaximum() / h_phi_mc_sideband->GetMaximum();
    h_phi_mc_sideband->Scale(scale_phi_sideband);
    
    // Style histograms - Phi
    bin_width = get_bin_width(h_Phi_data_total);
    h_Phi_data_total->SetTitle("");
    h_Phi_data_total->SetXTitle("#Phi (rad)");
    h_Phi_data_total->SetYTitle(TString::Format("Events / %.3f rad", bin_width));
    h_Phi_data_total->SetMarkerStyle(1);
    h_Phi_data_total->SetMarkerColor(kBlack);
    h_Phi_data_total->SetLineColor(kBlack);
    h_Phi_data_total->SetMinimum(0);
    h_Phi_data_total->SetMaximum(h_Phi_data_signal->GetMaximum() * 1.1);
    
    h_Phi_data_signal->SetLineColor(kBlue);
    h_Phi_data_signal->SetLineWidth(2);
    
    h_Phi_mc_signal->SetMarkerColor(kCyan+2);
    h_Phi_mc_signal->SetMarkerStyle(24);
    h_Phi_mc_signal->SetMarkerSize(0.8);
    h_Phi_mc_signal->SetLineColor(kCyan+2);
    
    h_Phi_data_sideband->SetLineColor(kRed+1);
    h_Phi_data_sideband->SetLineWidth(2);
    
    h_Phi_mc_sideband->SetMarkerColor(kRed-7);
    h_Phi_mc_sideband->SetMarkerStyle(25);
    h_Phi_mc_sideband->SetMarkerSize(0.8);
    h_Phi_mc_sideband->SetLineColor(kRed-7);
    
    // Scale MC to match data maximum
    double scale_Phi_signal = h_Phi_data_signal->GetMaximum() / h_Phi_mc_signal->GetMaximum();
    h_Phi_mc_signal->Scale(scale_Phi_signal);
    
    double scale_Phi_sideband = h_Phi_data_sideband->GetMaximum() / h_Phi_mc_sideband->GetMaximum();
    h_Phi_mc_sideband->Scale(scale_Phi_sideband);
    
    // Style histograms - cos(theta_h)
    bin_width = get_bin_width(h_theta_h_data_total);
    h_theta_h_data_total->SetTitle("");
    h_theta_h_data_total->SetXTitle("cos(#theta_{H})");
    h_theta_h_data_total->SetYTitle(TString::Format("Events / %.3f", bin_width));
    h_theta_h_data_total->SetMarkerStyle(1);
    h_theta_h_data_total->SetMarkerColor(kBlack);
    h_theta_h_data_total->SetLineColor(kBlack);
    h_theta_h_data_total->SetMinimum(0);
    h_theta_h_data_total->SetMaximum(h_theta_h_data_signal->GetMaximum() * 1.1);
    
    h_theta_h_data_signal->SetLineColor(kBlue);
    h_theta_h_data_signal->SetLineWidth(2);
    
    h_theta_h_mc_signal->SetMarkerColor(kCyan+2);
    h_theta_h_mc_signal->SetMarkerStyle(24);
    h_theta_h_mc_signal->SetMarkerSize(0.8);
    h_theta_h_mc_signal->SetLineColor(kCyan+2);
    
    h_theta_h_data_sideband->SetLineColor(kRed+1);
    h_theta_h_data_sideband->SetLineWidth(2);
    
    h_theta_h_mc_sideband->SetMarkerColor(kRed-7);
    h_theta_h_mc_sideband->SetMarkerStyle(25);
    h_theta_h_mc_sideband->SetMarkerSize(0.8);
    h_theta_h_mc_sideband->SetLineColor(kRed-7);
    
    // Scale MC to match data maximum
    double scale_theta_h_signal = h_theta_h_data_signal->GetMaximum() / h_theta_h_mc_signal->GetMaximum();
    h_theta_h_mc_signal->Scale(scale_theta_h_signal);
    
    double scale_theta_h_sideband = h_theta_h_data_sideband->GetMaximum() / h_theta_h_mc_sideband->GetMaximum();
    h_theta_h_mc_sideband->Scale(scale_theta_h_sideband);
    
    // Style histograms - phi_h
    bin_width = get_bin_width(h_phi_h_data_total);
    h_phi_h_data_total->SetTitle("");
    h_phi_h_data_total->SetXTitle("#phi_{H} (rad)");
    h_phi_h_data_total->SetYTitle(TString::Format("Events / %.3f rad", bin_width));
    h_phi_h_data_total->SetMarkerStyle(1);
    h_phi_h_data_total->SetMarkerColor(kBlack);
    h_phi_h_data_total->SetLineColor(kBlack);
    h_phi_h_data_total->SetMinimum(0);
    h_phi_h_data_total->SetMaximum(h_phi_h_data_signal->GetMaximum() * 1.1);
    
    h_phi_h_data_signal->SetLineColor(kBlue);
    h_phi_h_data_signal->SetLineWidth(2);
    
    h_phi_h_mc_signal->SetMarkerColor(kCyan+2);
    h_phi_h_mc_signal->SetMarkerStyle(24);
    h_phi_h_mc_signal->SetMarkerSize(0.8);
    h_phi_h_mc_signal->SetLineColor(kCyan+2);
    
    h_phi_h_data_sideband->SetLineColor(kRed+1);
    h_phi_h_data_sideband->SetLineWidth(2);
    
    h_phi_h_mc_sideband->SetMarkerColor(kRed-7);
    h_phi_h_mc_sideband->SetMarkerStyle(25);
    h_phi_h_mc_sideband->SetMarkerSize(0.8);
    h_phi_h_mc_sideband->SetLineColor(kRed-7);
    
    // Scale MC to match data maximum
    double scale_phi_h_signal = h_phi_h_data_signal->GetMaximum() / h_phi_h_mc_signal->GetMaximum();
    h_phi_h_mc_signal->Scale(scale_phi_h_signal);
    
    double scale_phi_h_sideband = h_phi_h_data_sideband->GetMaximum() / h_phi_h_mc_sideband->GetMaximum();
    h_phi_h_mc_sideband->Scale(scale_phi_h_sideband);
    
    // Style histograms - lambda
    bin_width = get_bin_width(h_lambda_data_total);
    h_lambda_data_total->SetTitle("");
    h_lambda_data_total->SetXTitle("#lambda");
    h_lambda_data_total->SetYTitle(TString::Format("Events / %.3f", bin_width));
    h_lambda_data_total->SetMarkerStyle(1);
    h_lambda_data_total->SetMarkerColor(kBlack);
    h_lambda_data_total->SetLineColor(kBlack);
    h_lambda_data_total->SetMinimum(0);
    h_lambda_data_total->SetMaximum(h_lambda_data_signal->GetMaximum() * 1.1);
    
    h_lambda_data_signal->SetLineColor(kBlue);
    h_lambda_data_signal->SetLineWidth(2);
    
    h_lambda_mc_signal->SetMarkerColor(kCyan+2);
    h_lambda_mc_signal->SetMarkerStyle(24);
    h_lambda_mc_signal->SetMarkerSize(0.8);
    h_lambda_mc_signal->SetLineColor(kCyan+2);
    
    h_lambda_data_sideband->SetLineColor(kRed+1);
    h_lambda_data_sideband->SetLineWidth(2);
    h_lambda_data_sideband->SetLineStyle(2);
    
    h_lambda_mc_sideband->SetMarkerColor(kRed-7);
    h_lambda_mc_sideband->SetMarkerStyle(25);
    h_lambda_mc_sideband->SetMarkerSize(0.8);
    h_lambda_mc_sideband->SetLineColor(kRed-7);
    
    // Scale MC to match data maximum
    double scale_lambda_signal = h_lambda_data_signal->GetMaximum() / h_lambda_mc_signal->GetMaximum();
    h_lambda_mc_signal->Scale(scale_lambda_signal);
    
    double scale_lambda_sideband = h_lambda_data_sideband->GetMaximum() / h_lambda_mc_sideband->GetMaximum();
    h_lambda_mc_sideband->Scale(scale_lambda_sideband);
    
    // // Style acceptance histograms
    // h_theta_acceptance->SetLineColor(kGray+2);
    // h_theta_acceptance->SetLineWidth(2);
    // h_theta_acceptance->SetLineStyle(1);
    
    // h_phi_acceptance->SetLineColor(kGray+2);
    // h_phi_acceptance->SetLineWidth(2);
    // h_phi_acceptance->SetLineStyle(1);
    
    // h_Phi_acceptance->SetLineColor(kGray+2);
    // h_Phi_acceptance->SetLineWidth(2);
    // h_Phi_acceptance->SetLineStyle(1);
    
    // h_theta_h_acceptance->SetLineColor(kGray+2);
    // h_theta_h_acceptance->SetLineWidth(2);
    // h_theta_h_acceptance->SetLineStyle(1);
    
    // h_phi_h_acceptance->SetLineColor(kGray+2);
    // h_phi_h_acceptance->SetLineWidth(2);
    // h_phi_h_acceptance->SetLineStyle(1);

    // h_lambda_acceptance->SetLineColor(kGray+2);
    // h_lambda_acceptance->SetLineWidth(2);
    // h_lambda_acceptance->SetLineStyle(1);
    
    // Create canvas for saving plots
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    
    // Draw theta
    TLegend* legend_theta = new TLegend(0.15, 0.7, 0.33, 0.88);
    legend_theta->AddEntry(h_theta_data_total, "Final Data", "lp");
    legend_theta->AddEntry(h_theta_data_signal, "Data Signal", "l");
    legend_theta->AddEntry(h_theta_mc_signal, TString::Format("MC Signal (scale=%.3f)", scale_theta_signal), "l");
    legend_theta->AddEntry(h_theta_data_sideband, "Data Sideband", "l");
    legend_theta->AddEntry(h_theta_mc_sideband, TString::Format("MC Sideband (scale=%.3f)", scale_theta_sideband), "l");
    // legend_theta->AddEntry(h_theta_acceptance, "Acceptance", "l");
    h_theta_data_total->Draw("E");
    h_theta_data_signal->Draw("HIST SAME");
    h_theta_mc_signal->Draw("E SAME");
    h_theta_data_sideband->Draw("HIST SAME");
    h_theta_mc_sideband->Draw("E SAME");
    
    // Draw acceptance on secondary y-axis
    // gPad->Update();
    // double y_max = gPad->GetUymax();
    // double acc_max = h_theta_acceptance->GetMaximum();
    // double scale_factor = y_max / (acc_max * 1.1);  // 1.1 for some headroom
    // TH1F* h_theta_acc_scaled = (TH1F*)h_theta_acceptance->Clone("h_theta_acc_scaled");
    // h_theta_acc_scaled->Scale(scale_factor);
    // h_theta_acc_scaled->Draw("HIST SAME");
    
    // // Create right axis for acceptance
    // TGaxis* axis_theta = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
    //                                  gPad->GetUxmax(), gPad->GetUymax(),
    //                                  0, acc_max * 1.1, 510, "+L");
    // axis_theta->SetLineColor(kGray+2);
    // axis_theta->SetLabelColor(kGray+2);
    // axis_theta->SetTitleColor(kGray+2);
    // axis_theta->SetTitle("Acceptance");
    // axis_theta->Draw();
    
    legend_theta->Draw("SAME");
    c1->SaveAs("final_cos_theta.pdf");
    c1->Clear();
    delete legend_theta;
    
    // Draw phi
    TLegend* legend_phi = new TLegend(0.41, 0.7, 0.59, 0.88);
    legend_phi->AddEntry(h_phi_data_total, "Final Data", "lp");
    legend_phi->AddEntry(h_phi_data_signal, "Data Signal", "l");
    legend_phi->AddEntry(h_phi_mc_signal, TString::Format("MC Signal (scale=%.3f)", scale_phi_signal), "l");
    legend_phi->AddEntry(h_phi_data_sideband, "Data Sideband", "l");
    legend_phi->AddEntry(h_phi_mc_sideband, TString::Format("MC Sideband (scale=%.3f)", scale_phi_sideband), "l");
    // legend_phi->AddEntry(h_phi_acceptance, "Acceptance", "l");
    h_phi_data_total->Draw("E");
    h_phi_data_signal->Draw("HIST SAME");
    h_phi_mc_signal->Draw("E SAME");
    h_phi_data_sideband->Draw("HIST SAME");
    h_phi_mc_sideband->Draw("E SAME");
    
    // Draw acceptance on secondary y-axis
    // gPad->Update();
    // y_max = gPad->GetUymax();
    // acc_max = h_phi_acceptance->GetMaximum();
    // scale_factor = y_max / (acc_max * 1.1);
    // TH1F* h_phi_acc_scaled = (TH1F*)h_phi_acceptance->Clone("h_phi_acc_scaled");
    // h_phi_acc_scaled->Scale(scale_factor);
    // h_phi_acc_scaled->Draw("HIST SAME");
    
    // TGaxis* axis_phi = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
    //                                gPad->GetUxmax(), gPad->GetUymax(),
    //                                0, acc_max * 1.1, 510, "+L");
    // axis_phi->SetLineColor(kGray+2);
    // axis_phi->SetLabelColor(kGray+2);
    // axis_phi->SetTitleColor(kGray+2);
    // axis_phi->SetTitle("Acceptance");
    // axis_phi->Draw();
    
    legend_phi->Draw("SAME");
    c1->SaveAs("final_phi.pdf");
    c1->Clear();
    delete legend_phi;
    
    // Draw Phi
    TLegend* legend_Phi =  new TLegend(0.7, 0.7, 0.88, 0.88);
    legend_Phi->AddEntry(h_Phi_data_total, "Final Data", "lp");
    legend_Phi->AddEntry(h_Phi_data_signal, "Data Signal", "l");
    legend_Phi->AddEntry(h_Phi_mc_signal, TString::Format("MC Signal (scale=%.3f)", scale_Phi_signal), "l");
    legend_Phi->AddEntry(h_Phi_data_sideband, "Data Sideband", "l");
    legend_Phi->AddEntry(h_Phi_mc_sideband, TString::Format("MC Sideband (scale=%.3f)", scale_Phi_sideband), "l");
    // legend_Phi->AddEntry(h_Phi_acceptance, "Acceptance", "l");
    
    h_Phi_data_total->Draw("E");
    h_Phi_data_signal->Draw("HIST SAME");
    h_Phi_mc_signal->Draw("E SAME");
    h_Phi_data_sideband->Draw("HIST SAME");
    h_Phi_mc_sideband->Draw("E SAME");
    
    // // Draw acceptance on secondary y-axis
    // gPad->Update();
    // y_max = gPad->GetUymax();
    // acc_max = h_Phi_acceptance->GetMaximum();
    // scale_factor = y_max / (acc_max * 1.1);
    // TH1F* h_Phi_acc_scaled = (TH1F*)h_Phi_acceptance->Clone("h_Phi_acc_scaled");
    // h_Phi_acc_scaled->Scale(scale_factor);
    // h_Phi_acc_scaled->Draw("HIST SAME");
    
    // TGaxis* axis_Phi = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
    //                                gPad->GetUxmax(), gPad->GetUymax(),
    //                                0, acc_max * 1.1, 510, "+L");
    // axis_Phi->SetLineColor(kGray+2);
    // axis_Phi->SetLabelColor(kGray+2);
    // axis_Phi->SetTitleColor(kGray+2);
    // axis_Phi->SetTitle("Acceptance");
    // axis_Phi->Draw();
    
    legend_Phi->Draw("SAME");
    c1->SaveAs("final_Phi.pdf");
    c1->Clear();
    delete legend_Phi;
    
    // Draw theta_h
    TLegend* legend_theta_h =  new TLegend(0.7, 0.7, 0.88, 0.88);
    legend_theta_h->AddEntry(h_theta_h_data_total, "Final Data", "lp");
    legend_theta_h->AddEntry(h_theta_h_data_signal, "Data Signal", "l");
    legend_theta_h->AddEntry(h_theta_h_mc_signal, TString::Format("MC Signal (scale=%.3f)", scale_theta_h_signal), "l");
    legend_theta_h->AddEntry(h_theta_h_data_sideband, "Data Sideband", "l");
    legend_theta_h->AddEntry(h_theta_h_mc_sideband, TString::Format("MC Sideband (scale=%.3f)", scale_theta_h_sideband), "l");
    // legend_theta_h->AddEntry(h_theta_h_acceptance, "Acceptance", "l");
    
    h_theta_h_data_total->Draw("E");
    h_theta_h_data_signal->Draw("HIST SAME");
    h_theta_h_mc_signal->Draw("E SAME");
    h_theta_h_data_sideband->Draw("HIST SAME");
    h_theta_h_mc_sideband->Draw("E SAME");
    
    // Draw acceptance on secondary y-axis
    // gPad->Update();
    // y_max = gPad->GetUymax();
    // acc_max = h_theta_h_acceptance->GetMaximum();
    // scale_factor = y_max / (acc_max * 1.1);
    // TH1F* h_theta_h_acc_scaled = (TH1F*)h_theta_h_acceptance->Clone("h_theta_h_acc_scaled");
    // h_theta_h_acc_scaled->Scale(scale_factor);
    // h_theta_h_acc_scaled->Draw("HIST SAME");
    
    // TGaxis* axis_theta_h = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
    //                                    gPad->GetUxmax(), gPad->GetUymax(),
    //                                    0, acc_max * 1.1, 510, "+L");
    // axis_theta_h->SetLineColor(kGray+2);
    // axis_theta_h->SetLabelColor(kGray+2);
    // axis_theta_h->SetTitleColor(kGray+2);
    // axis_theta_h->SetTitle("Acceptance");
    // axis_theta_h->Draw();
    
    legend_theta_h->Draw("SAME");
    c1->SaveAs("final_cos_theta_h.pdf");
    c1->Clear();
    delete legend_theta_h;
    
    // Draw phi_h
    TLegend* legend_phi_h =  new TLegend(0.7, 0.7, 0.88, 0.88);
    legend_phi_h->AddEntry(h_phi_h_data_total, "Final Data", "lp");
    legend_phi_h->AddEntry(h_phi_h_data_signal, "Data Signal", "l");
    legend_phi_h->AddEntry(h_phi_h_mc_signal, TString::Format("MC Signal (scale=%.3f)", scale_phi_h_signal), "l");
    legend_phi_h->AddEntry(h_phi_h_data_sideband, "Data Sideband", "l");
    legend_phi_h->AddEntry(h_phi_h_mc_sideband, TString::Format("MC Sideband (scale=%.3f)", scale_phi_h_sideband), "l");
    // legend_phi_h->AddEntry(h_phi_h_acceptance, "Acceptance", "l");
    
    h_phi_h_data_total->Draw("E");
    h_phi_h_data_signal->Draw("HIST SAME");
    h_phi_h_mc_signal->Draw("E SAME");
    h_phi_h_data_sideband->Draw("HIST SAME");
    h_phi_h_mc_sideband->Draw("E SAME");
    
    // Draw acceptance on secondary y-axis
    // gPad->Update();
    // y_max = gPad->GetUymax();
    // acc_max = h_phi_h_acceptance->GetMaximum();
    // scale_factor = y_max / (acc_max * 1.1);
    // TH1F* h_phi_h_acc_scaled = (TH1F*)h_phi_h_acceptance->Clone("h_phi_h_acc_scaled");
    // h_phi_h_acc_scaled->Scale(scale_factor);
    // h_phi_h_acc_scaled->Draw("HIST SAME");
    
    // TGaxis* axis_phi_h = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
    //                                  gPad->GetUxmax(), gPad->GetUymax(),
    //                                  0, acc_max * 1.1, 510, "+L");
    // axis_phi_h->SetLineColor(kGray+2);
    // axis_phi_h->SetLabelColor(kGray+2);
    // axis_phi_h->SetTitleColor(kGray+2);
    // axis_phi_h->SetTitle("Acceptance");
    // axis_phi_h->Draw();
    
    legend_phi_h->Draw("SAME");
    c1->SaveAs("final_phi_h.pdf");
    c1->Clear();
    delete legend_phi_h;
    
    // Draw lambda
    TLegend* legend_lambda = new TLegend(0.15, 0.7, 0.33, 0.88);
    legend_lambda->AddEntry(h_lambda_data_total, "Final Data", "lp");
    legend_lambda->AddEntry(h_lambda_data_signal, "Data Signal", "l");
    legend_lambda->AddEntry(h_lambda_mc_signal, TString::Format("MC Signal (scale=%.3f)", scale_lambda_signal), "l");
    legend_lambda->AddEntry(h_lambda_data_sideband, "Data Sideband", "l");
    legend_lambda->AddEntry(h_lambda_mc_sideband, TString::Format("MC Sideband (scale=%.3f)", scale_lambda_sideband), "l");
    // legend_lambda->AddEntry(h_lambda_acceptance, "Acceptance", "l");
    
    h_lambda_data_total->Draw("E");
    h_lambda_data_signal->Draw("HIST SAME");
    h_lambda_mc_signal->Draw("E SAME");
    h_lambda_data_sideband->Draw("HIST SAME");
    h_lambda_mc_sideband->Draw("E SAME");
    
    // Draw acceptance on secondary y-axis
    // gPad->Update();
    // y_max = gPad->GetUymax();
    // acc_max = h_lambda_acceptance->GetMaximum();
    // scale_factor = y_max / (acc_max * 1.1);
    // TH1F* h_lambda_acc_scaled = (TH1F*)h_lambda_acceptance->Clone("h_lambda_acc_scaled");
    // h_lambda_acc_scaled->Scale(scale_factor);
    // h_lambda_acc_scaled->Draw("HIST SAME");
    
    // TGaxis* axis_lambda = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
    //                                   gPad->GetUxmax(), gPad->GetUymax(),
    //                                   0, acc_max * 1.1, 510, "+L");
    // axis_lambda->SetLineColor(kGray+2);
    // axis_lambda->SetLabelColor(kGray+2);
    // axis_lambda->SetTitleColor(kGray+2);
    // axis_lambda->SetTitle("Acceptance");
    // axis_lambda->Draw();
    
    legend_lambda->Draw("SAME");
    c1->SaveAs("final_lambda.pdf");
    c1->Clear();
    delete legend_lambda;
    delete c1;
    
    return;
}

/**
 * @brief Plot 2D angle distributions for signal and sideband data
 * 
 * @param NT tree name (unused but kept for API consistency)
 * @param CATEGORY category name (unused but kept for API consistency)
 * @param angles_signal pre-calculated angles for signal events
 * @param angles_sideband pre-calculated angles for sideband (background) events
 * @param output_prefix prefix for output filenames ("data" or "mc")
 */
void plot_2d_angles(
    TString NT,
    TString CATEGORY,
    const AngleData& angles_signal,
    const AngleData& angles_sideband,
    TString output_prefix
)
{
    std::cout << "Creating 2D angle plots for " << output_prefix << "..." << std::endl;
    
    // Create 2D histograms
    TH2F* h_theta_phi_signal = new TH2F("h_theta_phi_signal", "", 50, -1, 1, 50, -TMath::Pi(), TMath::Pi());
    TH2F* h_theta_phi_sideband = new TH2F("h_theta_phi_sideband", "", 50, -1, 1, 50, -TMath::Pi(), TMath::Pi());
    
    TH2F* h_theta_h_phi_h_signal = new TH2F("h_theta_h_phi_h_signal", "", 50, -1, 1, 50, -TMath::Pi(), TMath::Pi());
    TH2F* h_theta_h_phi_h_sideband = new TH2F("h_theta_h_phi_h_sideband", "", 50, -1, 1, 50, -TMath::Pi(), TMath::Pi());
    
    TH2F* h_Phi_phi_signal = new TH2F("h_Phi_phi_signal", "", 50, -TMath::Pi(), TMath::Pi(), 50, -TMath::Pi(), TMath::Pi());
    TH2F* h_Phi_phi_sideband = new TH2F("h_Phi_phi_sideband", "", 50, -TMath::Pi(), TMath::Pi(), 50, -TMath::Pi(), TMath::Pi());
    
    // Fill signal histograms
    for (size_t i = 0; i < angles_signal.theta.size(); i++) {
        double weight = angles_signal.weights[i];
        h_theta_phi_signal->Fill(TMath::Cos(angles_signal.theta[i]), angles_signal.phi[i], weight);
        h_theta_h_phi_h_signal->Fill(TMath::Cos(angles_signal.theta_h[i]), angles_signal.phi_h[i], weight);
        h_Phi_phi_signal->Fill(angles_signal.Phi[i], angles_signal.phi[i], weight);
    }
    
    // Fill sideband histograms
    for (size_t i = 0; i < angles_sideband.theta.size(); i++) {
        double weight = angles_sideband.weights[i];
        h_theta_phi_sideband->Fill(TMath::Cos(angles_sideband.theta[i]), angles_sideband.phi[i], weight);
        h_theta_h_phi_h_sideband->Fill(TMath::Cos(angles_sideband.theta_h[i]), angles_sideband.phi_h[i], weight);
        h_Phi_phi_sideband->Fill(angles_sideband.Phi[i], angles_sideband.phi[i], weight);
    }
    
    // Create total histograms (signal - sideband)
    TH2F* h_theta_phi_total = (TH2F*)h_theta_phi_signal->Clone("h_theta_phi_total");
    h_theta_phi_total->Add(h_theta_phi_sideband, -1.0);
    
    TH2F* h_theta_h_phi_h_total = (TH2F*)h_theta_h_phi_h_signal->Clone("h_theta_h_phi_h_total");
    h_theta_h_phi_h_total->Add(h_theta_h_phi_h_sideband, -1.0);
    
    TH2F* h_Phi_phi_total = (TH2F*)h_Phi_phi_signal->Clone("h_Phi_phi_total");
    h_Phi_phi_total->Add(h_Phi_phi_sideband, -1.0);
    
    // Style histograms
    h_theta_phi_total->SetTitle("");
    h_theta_phi_total->SetXTitle("cos(#theta)");
    h_theta_phi_total->SetYTitle("#phi (rad)");
    h_theta_phi_total->SetZTitle("Events");
    
    h_theta_h_phi_h_total->SetTitle("");
    h_theta_h_phi_h_total->SetXTitle("cos(#theta_{H})");
    h_theta_h_phi_h_total->SetYTitle("#phi_{H} (rad)");
    h_theta_h_phi_h_total->SetZTitle("Events");
    
    h_Phi_phi_total->SetTitle("");
    h_Phi_phi_total->SetXTitle("#Phi (rad)");
    h_Phi_phi_total->SetYTitle("#phi (rad)");
    h_Phi_phi_total->SetZTitle("Events");
    
    // Create canvas and draw
    TCanvas* c2d = new TCanvas("c2d", "c2d", 800, 600);
    c2d->SetRightMargin(0.15);
    
    // Draw cos(theta) vs phi
    h_theta_phi_total->Draw("COL0");
    c2d->SaveAs("final_2d_" + output_prefix + "_cos_theta_vs_phi.pdf");
    c2d->Clear();
    
    // Draw cos(theta_h) vs phi_h
    h_theta_h_phi_h_total->Draw("COL0");
    c2d->SaveAs("final_2d_" + output_prefix + "_cos_theta_h_vs_phi_h.pdf");
    c2d->Clear();
    
    // Draw Phi vs phi
    h_Phi_phi_total->Draw("COL0");
    c2d->SaveAs("final_2d_" + output_prefix + "_Phi_vs_phi.pdf");
    c2d->Clear();
    
    delete c2d;
    delete h_theta_phi_signal;
    delete h_theta_phi_sideband;
    delete h_theta_phi_total;
    delete h_theta_h_phi_h_signal;
    delete h_theta_h_phi_h_sideband;
    delete h_theta_h_phi_h_total;
    delete h_Phi_phi_signal;
    delete h_Phi_phi_sideband;
    delete h_Phi_phi_total;
    
    return;
}

/**
 * @brief Plot 2D acceptance distributions from phase space Monte Carlo
 * 
 * @param NT tree name (unused but kept for API consistency)
 * @param CATEGORY category name (unused but kept for API consistency)
 * @param angles_acc pre-calculated angles for accepted phase space events
 * @param angles_gen pre-calculated angles for generated phase space events
 */
void plot_2d_acceptance(
    TString NT,
    TString CATEGORY,
    const AngleData& angles_acc,
    const AngleData& angles_gen
)
{
    std::cout << "Creating 2D acceptance plots..." << std::endl;
    
    // Create 2D histograms for accepted and generated events
    TH2F* h_theta_phi_acc = new TH2F("h_theta_phi_acc", "", 50, -1, 1, 50, -TMath::Pi(), TMath::Pi());
    TH2F* h_theta_phi_gen = new TH2F("h_theta_phi_gen", "", 50, -1, 1, 50, -TMath::Pi(), TMath::Pi());
    
    TH2F* h_theta_h_phi_h_acc = new TH2F("h_theta_h_phi_h_acc", "", 50, -1, 1, 50, -TMath::Pi(), TMath::Pi());
    TH2F* h_theta_h_phi_h_gen = new TH2F("h_theta_h_phi_h_gen", "", 50, -1, 1, 50, -TMath::Pi(), TMath::Pi());
    
    TH2F* h_Phi_phi_acc = new TH2F("h_Phi_phi_acc", "", 50, -TMath::Pi(), TMath::Pi(), 50, -TMath::Pi(), TMath::Pi());
    TH2F* h_Phi_phi_gen = new TH2F("h_Phi_phi_gen", "", 50, -TMath::Pi(), TMath::Pi(), 50, -TMath::Pi(), TMath::Pi());
    
    // Fill accepted events histograms
    for (size_t i = 0; i < angles_acc.theta.size(); i++) {
        double weight = angles_acc.weights[i];
        h_theta_phi_acc->Fill(TMath::Cos(angles_acc.theta[i]), angles_acc.phi[i], weight);
        h_theta_h_phi_h_acc->Fill(TMath::Cos(angles_acc.theta_h[i]), angles_acc.phi_h[i], weight);
        h_Phi_phi_acc->Fill(angles_acc.Phi[i], angles_acc.phi[i], weight);
    }
    
    // Fill generated events histograms
    for (size_t i = 0; i < angles_gen.theta.size(); i++) {
        double weight = angles_gen.weights[i];
        h_theta_phi_gen->Fill(TMath::Cos(angles_gen.theta[i]), angles_gen.phi[i], weight);
        h_theta_h_phi_h_gen->Fill(TMath::Cos(angles_gen.theta_h[i]), angles_gen.phi_h[i], weight);
        h_Phi_phi_gen->Fill(angles_gen.Phi[i], angles_gen.phi[i], weight);
    }
    
    // Calculate acceptance (acc/gen)
    TH2F* h_theta_phi_acceptance = (TH2F*)h_theta_phi_acc->Clone("h_theta_phi_acceptance");
    h_theta_phi_acceptance->Divide(h_theta_phi_gen);
    
    TH2F* h_theta_h_phi_h_acceptance = (TH2F*)h_theta_h_phi_h_acc->Clone("h_theta_h_phi_h_acceptance");
    h_theta_h_phi_h_acceptance->Divide(h_theta_h_phi_h_gen);
    
    TH2F* h_Phi_phi_acceptance = (TH2F*)h_Phi_phi_acc->Clone("h_Phi_phi_acceptance");
    h_Phi_phi_acceptance->Divide(h_Phi_phi_gen);
    
    // Style histograms
    h_theta_phi_acceptance->SetTitle("");
    h_theta_phi_acceptance->SetXTitle("cos(#theta)");
    h_theta_phi_acceptance->SetYTitle("#phi (rad)");
    h_theta_phi_acceptance->SetZTitle("Acceptance");
    h_theta_phi_acceptance->SetMinimum(0);
    h_theta_phi_acceptance->SetMaximum(1);
    
    h_theta_h_phi_h_acceptance->SetTitle("");
    h_theta_h_phi_h_acceptance->SetXTitle("cos(#theta_{H})");
    h_theta_h_phi_h_acceptance->SetYTitle("#phi_{H} (rad)");
    h_theta_h_phi_h_acceptance->SetZTitle("Acceptance");
    h_theta_h_phi_h_acceptance->SetMinimum(0);
    h_theta_h_phi_h_acceptance->SetMaximum(1);
    
    h_Phi_phi_acceptance->SetTitle("");
    h_Phi_phi_acceptance->SetXTitle("#Phi (rad)");
    h_Phi_phi_acceptance->SetYTitle("#phi (rad)");
    h_Phi_phi_acceptance->SetZTitle("Acceptance");
    h_Phi_phi_acceptance->SetMinimum(0);
    h_Phi_phi_acceptance->SetMaximum(1);
    
    // Create canvas and draw
    TCanvas* c2d_acc = new TCanvas("c2d_acc", "c2d_acc", 800, 600);
    c2d_acc->SetRightMargin(0.15);
    
    // Draw cos(theta) vs phi acceptance
    h_theta_phi_acceptance->Draw("COL0");
    c2d_acc->SaveAs("final_2d_acceptance_cos_theta_vs_phi.pdf");
    c2d_acc->Clear();
    
    // Draw cos(theta_h) vs phi_h acceptance
    h_theta_h_phi_h_acceptance->Draw("COL0");
    c2d_acc->SaveAs("final_2d_acceptance_cos_theta_h_vs_phi_h.pdf");
    c2d_acc->Clear();
    
    // Draw Phi vs phi acceptance
    h_Phi_phi_acceptance->Draw("COL0");
    c2d_acc->SaveAs("final_2d_acceptance_Phi_vs_phi.pdf");
    c2d_acc->Clear();
    
    delete c2d_acc;
    delete h_theta_phi_acc;
    delete h_theta_phi_gen;
    delete h_theta_phi_acceptance;
    delete h_theta_h_phi_h_acc;
    delete h_theta_h_phi_h_gen;
    delete h_theta_h_phi_h_acceptance;
    delete h_Phi_phi_acc;
    delete h_Phi_phi_gen;
    delete h_Phi_phi_acceptance;
    
    return;
}