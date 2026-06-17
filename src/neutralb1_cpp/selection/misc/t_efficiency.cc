/**
 * @file t_efficiency.cc
 * @author Kevin Scheuer
 * @brief Quick check of t efficiency / acceptance to know how much MC will be modified.
 * 
 */


#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TString.h"
#include "TH1.h"
#include "TCanvas.h"

#include "neutralb1/fit_utils.h"

void t_efficiency() {
    gStyle->SetOptStat(0);
    TString tree_dir = "/home/kscheuer/work/neutralb1/data/FSRoot/phasespace/";

    // input files 
    TString acc_file = tree_dir + "allPeriods_ver03_phasespace.root";
    TString gen_file = tree_dir + "allPeriods_ver03_gen_phasespace.root";

    TString NT = "ntFSGlueX_100_112";
    TString CATEGORY = "pi0pi0pippim";

    // get files and trees
    TFile *f_acc = TFile::Open(acc_file);
    TFile *f_gen = TFile::Open(gen_file);

    if (!f_acc || f_acc->IsZombie()) {
        std::cerr << "Failed to open acceptance file: " << acc_file << std::endl;
        return;
    }

    if (!f_gen || f_gen->IsZombie()) {
        std::cerr << "Failed to open generation file: " << gen_file << std::endl;
        return;
    }

    TTree *tree_acc = (TTree*)f_acc->Get(NT);
    TTree *tree_gen = (TTree*)f_gen->Get(NT);

    if (!tree_acc) {
        std::cerr << "Failed to find tree " << NT << " in " << acc_file << std::endl;
        return;
    }

    if (!tree_gen) {
        std::cerr << "Failed to find tree " << NT << " in " << gen_file << std::endl;
        return;
    }

    TH1F *h_t_acc = new TH1F(
        "h_t_acc", 
        "", 
        100, 0.0, 1.0);
    TH1F *h_t_gen = new TH1F(
        "h_t_gen", 
        "", 
        100, 0.0, 1.0);

    // Fill the histograms
    tree_acc->Draw("t>>h_t_acc", "weight", "goff");
    tree_gen->Draw("t>>h_t_gen", "", "goff");

    TH1F *h_t_eff = (TH1F*)h_t_acc->Clone("h_t_eff");
    h_t_eff->Divide(h_t_gen);

    double bin_width = get_bin_width( h_t_eff );

    h_t_eff->SetTitle("Acceptance Efficiency");
    h_t_eff->GetXaxis()->SetTitle("-t (GeV^{2})");
    h_t_eff->SetYTitle(TString::Format("Efficiency / %.3f GeV", bin_width));

    TCanvas *c1 = new TCanvas("c1", "Efficiency", 800, 600);
    h_t_eff->Draw();
    c1->SaveAs("t_efficiency.pdf");
}