/**
 * @file plot_dalitz_xy.cc
 * @author Kevin Scheuer
 * @brief Plot Dalitz plot in XY coordinates for signal MC and data events
 *
 * As a reminder, signal MC is only in PARA_0 orientation, so we only compare to
 * data in PARA_0. Also, since the dalitz functions are intended to be used in the
 * AmpTools PWA framework and subsequent plotting, we'll need to use the data files
 * that contain the relabelled Amptools parameters.
 */

#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH2.h"

struct DalitzXY
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> weights;
};

DalitzXY calculate_dalitz_xy(TTree *tree);
void plot_dalitz(
    const DalitzXY &dalitz_signal,
    const DalitzXY &dalitz_sideband,
    const TString &output_file);

void plot_dalitz_xy()
{
    gStyle->SetOptStat(0);
    TString tree_dir = "/lustre24/expphy/volatile/halld/home/kscheuer/"
                       "FSRoot-skimmed-trees/final-amptools-trees/";

    // input files
    TString data_signal = tree_dir + "PARA_0_allPeriods_data_signal.root";
    TString data_sideband = tree_dir + "PARA_0_allPeriods_data_background.root";
    TString mc_signal = tree_dir + "PARA_0_allPeriods_ver03.1_mc_signal.root";
    TString mc_sideband = tree_dir + "PARA_0_allPeriods_ver03.1_mc_background.root";

    TString NT = "ntFSGlueX_100_112";
    TString CATEGORY = "pi0pi0pippim";

    // Open files once
    TFile *f_data_signal = TFile::Open(data_signal);
    TFile *f_data_sideband = TFile::Open(data_sideband);
    TFile *f_mc_signal = TFile::Open(mc_signal);
    TFile *f_mc_sideband = TFile::Open(mc_sideband);

    TTree *tree_data_signal = (TTree *)f_data_signal->Get(NT);
    TTree *tree_data_sideband = (TTree *)f_data_sideband->Get(NT);
    TTree *tree_mc_signal = (TTree *)f_mc_signal->Get(NT);
    TTree *tree_mc_sideband = (TTree *)f_mc_sideband->Get(NT);

    // calculate dalitz XY parameters for all events
    std::cout << "Calculating Dalitz XY for data signal..." << std::endl;
    DalitzXY dalitz_data_signal = calculate_dalitz_xy(tree_data_signal);
    std::cout << "Calculating Dalitz XY for data sideband..." << std::endl;
    DalitzXY dalitz_data_sideband = calculate_dalitz_xy(tree_data_sideband);
    std::cout << "Calculating Dalitz XY for MC signal..." << std::endl;
    DalitzXY dalitz_mc_signal = calculate_dalitz_xy(tree_mc_signal);
    std::cout << "Calculating Dalitz XY for MC sideband..." << std::endl;
    DalitzXY dalitz_mc_sideband = calculate_dalitz_xy(tree_mc_sideband);

    plot_dalitz(dalitz_data_signal, dalitz_data_sideband, "dalitz_xy_data.pdf");
    plot_dalitz(dalitz_mc_signal, dalitz_mc_sideband, "dalitz_xy_mc.pdf");

    // close files
    f_data_signal->Close();
    f_data_sideband->Close();
    f_mc_signal->Close();
    f_mc_sideband->Close();

    return;
}

DalitzXY calculate_dalitz_xy(TTree *tree)
{
    DalitzXY dalitz;

    // Set up branch addresses for 4-momenta
    double EnP3, PxP3, PyP3, PzP3;
    double EnP4, PxP4, PyP4, PzP4;
    double EnP5, PxP5, PyP5, PzP5;
    double weight = 1.0;
    double M4Pi = 0.0;

    tree->SetBranchAddress("M4Pi", &M4Pi);

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

    // Try to get weight if it exists
    if (tree->GetBranch("weight"))
    {
        tree->SetBranchAddress("weight", &weight);
    }

    Long64_t nentries = tree->GetEntries();
    dalitz.x.reserve(nentries);
    dalitz.y.reserve(nentries);
    dalitz.weights.reserve(nentries);

    for (Long64_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        // Create TLorentzVectors
        TLorentzVector pi0_omega(PxP3, PyP3, PzP3, EnP3);
        TLorentzVector pi_plus(PxP4, PyP4, PzP4, EnP4);
        TLorentzVector pi_minus(PxP5, PyP5, PzP5, EnP5);
        TLorentzVector omega = pi0_omega + pi_plus + pi_minus;

        double dalitz_s = (pi_plus + pi_minus).M2();   // s=M2(pip pim)
        double dalitz_t = (pi_plus + pi0_omega).M2();  // t=M2(pip pi0)
        double dalitz_u = (pi_minus + pi0_omega).M2(); // u=M2(pim pi0)
        double m3pi = (2 * pi_plus.M()) + pi0_omega.M();
        double dalitz_d = 2 * omega.M() * (omega.M() - m3pi);
        double dalitz_sc = (1 / 3.) * (omega.M2() + pi_plus.M2() + pi_minus.M2() + pi0_omega.M2());
        double dalitzx = sqrt(3) * (dalitz_t - dalitz_u) / dalitz_d;
        double dalitzy = 3 * (dalitz_sc - dalitz_s) / dalitz_d;

        dalitz.x.push_back(dalitzx);
        dalitz.y.push_back(dalitzy);
        dalitz.weights.push_back(weight);
    }

    return dalitz;
}

void plot_dalitz(
    const DalitzXY &dalitz_signal,
    const DalitzXY &dalitz_sideband,
    const TString &output_file)
{
    // create histograms
    TH2F *h_dalitz_xy_signal = new TH2F(
        "h_dalitz_xy_signal",
        "; X; Y",
        50, -1.0, 1.0,
        50, -1.0, 1.0);
    TH2F *h_dalitz_xy_sideband = new TH2F(
        "h_dalitz_xy_sideband",
        "; X; Y",
        50, -1.0, 1.0,
        50, -1.0, 1.0);

    // fill histogram with dalitz x,y data
    for (size_t i = 0; i < dalitz_signal.x.size(); i++)
    {
        h_dalitz_xy_signal->Fill(
            dalitz_signal.x[i],
            dalitz_signal.y[i],
            dalitz_signal.weights[i]);
    }

    for (size_t i = 0; i < dalitz_sideband.x.size(); i++)
    {
        h_dalitz_xy_sideband->Fill(
            dalitz_sideband.x[i],
            dalitz_sideband.y[i],
            dalitz_sideband.weights[i]);
    }

    // Create total histogram (signal - sideband)
    TH2F *h_dalitz_xy_total =
        (TH2F *)h_dalitz_xy_signal->Clone("h_dalitz_xy_total");
    h_dalitz_xy_total->Add(h_dalitz_xy_sideband, -1.0);

    // Set any negative bin contents to zero for better representation of holes
    for (int binx = 1; binx <= h_dalitz_xy_total->GetNbinsX(); ++binx)
    {
        for (int biny = 1; biny <= h_dalitz_xy_total->GetNbinsY(); ++biny)
        {
            if (h_dalitz_xy_total->GetBinContent(binx, biny) < 0)
                h_dalitz_xy_total->SetBinContent(binx, biny, 0);
        }
    }

    h_dalitz_xy_total->SetZTitle("Events");

    // create canvas and draw
    TCanvas *c1 = new TCanvas("c1", "Dalitz XY Plot", 800, 600);
    h_dalitz_xy_total->Draw("COLZ");
    c1->SaveAs(output_file);
    delete c1;
    delete h_dalitz_xy_signal;
    delete h_dalitz_xy_sideband;
    delete h_dalitz_xy_total;

    return;
}