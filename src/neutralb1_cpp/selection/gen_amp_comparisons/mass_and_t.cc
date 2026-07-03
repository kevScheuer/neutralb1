/**
 * @file mass_and_t.cc
 * @author Kevin Scheuer
 * @brief Plot the acceptance-corrected mass and -t distributions for data against MC
 * @date 2026-06-19
 * 
 * This script is to check if some gen_amp_V2 Monte Carlo data that we generated and
 * applied the acceptance to is consistent with the data. The acceptance is, of course,
 * not perfect, but provides some quick estimate for how our MC will be modified in the
 * full detector simulation. Thus if our MC is very different from data, we save time by
 * not running the full detector simulation on it, and we can generate new MC with
 * different parameters instead.
 * 
 * It also plots the original phasespace distributions from gen_amp_V2. Even though we
 * do not want the Breit-Wigners in our final acc/gen phasespace files, as we want the
 * distributions to be flat to allow our PWA model to weight them, we can use the
 * weighted events to get the phasespace distributions. gen_amp_V2 will track the
 * weights of the events, such that when using them it will "undo" the Breit-Wigner
 * weighting and give us the original phasespace distributions. This provides increased
 * precision in the areas where the Breit-Wigners are large, improving statistics over a
 * flat phasespace generation. This script is to plot the original phasespace
 * distributions 
 */

#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TH1.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

#include "FSBasic/FSHistogram.h"

#include "neutralb1/fit_utils.h"



void mass_and_t(TString gen_amp_V2_file)
{   
    gStyle->SetOptStat(0);

	const int n_bins = 50;
	const double t_min = 0.0;
	const double t_max = 1.0;
    const double m_min = 1.0;
    const double m_max = 2.0;

    std::cout << "Opening data and phasespace files..." << "\n";

	// Data files
	TString data_dir = "/home/kscheuer/work/neutralb1/data/FSRoot/GlueX/";
	TString data_signal = data_dir + "allPeriods_data_signal.root";
	TString data_sideband = data_dir + "allPeriods_data_background.root";

	// Phase-space files
	TString phase_space_dir = "/home/kscheuer/work/neutralb1/data/FSRoot/phasespace/";
	TString acc_file = phase_space_dir + "allPeriods_ver03_phasespace.root";
	TString gen_file = phase_space_dir + "allPeriods_ver03_gen_phasespace.root";

	TString NT = "ntFSGlueX_100_112";

	TFile *f_data_signal = TFile::Open(data_signal);
	TFile *f_data_sideband = TFile::Open(data_sideband);
	TFile *f_acc = TFile::Open(acc_file);
	TFile *f_gen = TFile::Open(gen_file);

	if (!f_data_signal || f_data_signal->IsZombie()) {
		std::cerr << "Failed to open signal file: " << data_signal << std::endl;
		return;
	}
	if (!f_data_sideband || f_data_sideband->IsZombie()) {
		std::cerr << "Failed to open sideband file: " << data_sideband << std::endl;
		return;
	}
	if (!f_acc || f_acc->IsZombie()) {
		std::cerr << "Failed to open acceptance file: " << acc_file << std::endl;
		return;
	}
	if (!f_gen || f_gen->IsZombie()) {
		std::cerr << "Failed to open generated phase-space file: " << gen_file << std::endl;
		return;
	}

	TTree *tree_data_signal = (TTree*)f_data_signal->Get(NT);
	TTree *tree_data_sideband = (TTree*)f_data_sideband->Get(NT);
	TTree *tree_acc = (TTree*)f_acc->Get(NT);
	TTree *tree_gen = (TTree*)f_gen->Get(NT);

	if (!tree_data_signal || !tree_data_sideband || !tree_acc || !tree_gen) {
		std::cerr << "Failed to find tree " << NT << " in one or more input files." << std::endl;
		return;
	}

    std::cout << "Building data histograms..." << "\n";

	// Build observed data histogram (signal - sideband)
	TH1F *h_t_signal = new TH1F("h_t_signal", "", n_bins, t_min, t_max);
	TH1F *h_t_sideband = new TH1F("h_t_sideband", "", n_bins, t_min, t_max);
    TH1F *h_m_signal = new TH1F("h_m_signal", "", n_bins, m_min, m_max);
    TH1F *h_m_sideband = new TH1F("h_m_sideband", "", n_bins, m_min, m_max);

	tree_data_signal->Draw("t>>h_t_signal", "", "goff");
	tree_data_sideband->Draw("t>>h_t_sideband", "weight", "goff");
    tree_data_signal->Draw("M4Pi>>h_m_signal", "", "goff");
    tree_data_sideband->Draw("M4Pi>>h_m_sideband", "weight", "goff");

	TH1F *h_t_data = (TH1F*)h_t_signal->Clone("h_t_data");
	h_t_data->Add(h_t_sideband, -1.0);
    TH1F *h_m_data = (TH1F*)h_m_signal->Clone("h_m_data");
    h_m_data->Add(h_m_sideband, -1.0);

    std::cout << "Calculating Acceptance..." << "\n";

	// Build acceptance efficiency histogram from phase-space MC
	TH1F *h_t_acc = new TH1F("h_t_acc", "", n_bins, t_min, t_max);
	TH1F *h_t_gen_ps = new TH1F("h_t_gen_ps", "", n_bins, t_min, t_max);
    TH1F *h_m_acc = new TH1F("h_m_acc", "", n_bins, m_min, m_max);
    TH1F *h_m_gen_ps = new TH1F("h_m_gen_ps", "", n_bins, m_min, m_max);

	tree_acc->Draw("t>>h_t_acc", "weight", "goff");
	tree_gen->Draw("t>>h_t_gen_ps", "", "goff");
    tree_acc->Draw("M4Pi>>h_m_acc", "weight", "goff");
    tree_gen->Draw("M4Pi>>h_m_gen_ps", "", "goff");

	TH1F *h_t_eff = (TH1F*)h_t_acc->Clone("h_t_eff");
	h_t_eff->SetDirectory(0);
	h_t_eff->Divide(h_t_gen_ps);

    TH1F *h_m_eff = (TH1F*)h_m_acc->Clone("h_m_eff");
    h_m_eff->SetDirectory(0);
    h_m_eff->Divide(h_m_gen_ps);

    // Get the mass and t distributions from the Monte Carlo file
    // Need to use FSROOT and a different tree due to default gen_amp_V2 output
    std::cout << "Getting MC histograms..." << "\n";

    TString t_expr = "-1.*MASS2(GLUEXTARGET,-1)";

    TH1F* h_t_mc = FSHistogram::getTH1F(
        gen_amp_V2_file, 
        "nt", 
        t_expr,
        TString::Format("(%i, %f, %f)", n_bins, t_min, t_max),
        ""
    );
    TH1F* h_t_mc_gen = FSHistogram::getTH1F(
        gen_amp_V2_file, 
        "nt", 
        t_expr,
        TString::Format("(%i, %f, %f)", n_bins, t_min, t_max),
        "weight"
    );

    TString m_expr = "MASS(2,3,4,5)";

    TH1F* h_m_mc = FSHistogram::getTH1F(
        gen_amp_V2_file, 
        "nt", 
        m_expr,
        TString::Format("(%i, %f, %f)", n_bins, m_min, m_max),
        ""
    );
    TH1F* h_m_mc_gen = FSHistogram::getTH1F(
        gen_amp_V2_file, 
        "nt", 
        m_expr,
        TString::Format("(%i, %f, %f)", n_bins, m_min, m_max),
        "weight"        
    );

    // modify the MC histograms to account for acceptance
    h_t_mc->Multiply(h_t_eff);
    h_m_mc->Multiply(h_m_eff);

    // multiply the MC histograms by the ratio of integrals or max to scale to data
    double t_scale = h_t_data->Integral() / h_t_mc->Integral();
    double m_scale = h_m_data->GetMaximum() / h_m_mc->GetMaximum();

    h_t_mc->Scale(t_scale);
    h_m_mc->Scale(m_scale);

     // plot the phasespace distributions for data and MC
    TCanvas *c_t_ps = new TCanvas("c_t_ps", "t distribution", 800, 600);
    const double t_bin_width_gen = get_bin_width(h_t_mc_gen);
    
    h_t_mc_gen->SetTitle("Generated t distribution");
    h_t_mc_gen->GetXaxis()->SetTitle("t (GeV^{2})");
    h_t_mc_gen->GetYaxis()->SetTitle(Form("Events / %.3f GeV^{2}", t_bin_width_gen));
    h_t_mc_gen->GetXaxis()->SetTitleSize(0.05);
    h_t_mc_gen->GetYaxis()->SetTitleSize(0.05);
    h_t_mc_gen->GetXaxis()->SetTitleOffset(0.8);
    h_t_mc_gen->GetYaxis()->SetTitleOffset(0.8);
    h_t_mc_gen->SetLineColor(kBlack);
    h_t_mc_gen->SetMarkerStyle(1);
    h_t_mc_gen->SetMarkerColor(kBlack);
    h_t_mc_gen->SetMinimum(0.0);

    h_t_mc_gen->Draw("E");
    c_t_ps->SaveAs("t_gen_amp.pdf");
    delete c_t_ps;

    TCanvas *c_m_ps = new TCanvas("c_m_ps", "M4Pi distribution", 800, 600);
    const double m_bin_width_gen = get_bin_width(h_m_mc_gen);

    h_m_mc_gen->SetTitle("Generated #pi^{+}#pi^{-}#pi^{0}#pi^{0} distribution");
    h_m_mc_gen->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}#pi^{0}#pi^{0}} (GeV)");
    h_m_mc_gen->GetYaxis()->SetTitle(Form("Events / %.3f GeV", m_bin_width_gen));
    h_m_mc_gen->GetXaxis()->SetTitleSize(0.05);
    h_m_mc_gen->GetYaxis()->SetTitleSize(0.05);
    h_m_mc_gen->GetXaxis()->SetTitleOffset(0.8);
    h_m_mc_gen->GetYaxis()->SetTitleOffset(0.8);
    h_m_mc_gen->SetLineColor(kBlack);
    h_m_mc_gen->SetMarkerStyle(1);
    h_m_mc_gen->SetMarkerColor(kBlack);
    h_m_mc_gen->SetMinimum(0.0);

    h_m_mc_gen->Draw("E");
    c_m_ps->SaveAs("m_gen_amp.pdf");
    delete c_m_ps;

    // plot the mass and t distributions for data and MC
    TCanvas *c_t = new TCanvas("c_t", "t distribution", 800, 600);
    
    const double t_bin_width = get_bin_width(h_t_data);   
    h_t_data->SetLineColor(kBlack); 
    h_t_data->SetTitle("");
    h_t_data->GetXaxis()->SetTitle("t (GeV^{2})");
    h_t_data->GetYaxis()->SetTitle(Form("Events / %.3f GeV^{2}", t_bin_width));
    h_t_data->GetXaxis()->SetTitleSize(0.05);
    h_t_data->GetYaxis()->SetTitleSize(0.05);
    h_t_data->GetXaxis()->SetTitleOffset(0.8);
    h_t_data->GetYaxis()->SetTitleOffset(0.8);
    h_t_mc->SetMarkerStyle(1);
    h_t_mc->SetMarkerColor(kRed);
    h_t_mc->SetLineColor(kRed);
    h_t_data->SetMinimum(0.0);
    h_t_data->SetMaximum(std::max(h_t_data->GetMaximum(), h_t_mc->GetMaximum()) * 1.1);\

    TLegend *t_legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    t_legend->AddEntry(h_t_data, "Data", "lep");
    t_legend->AddEntry(h_t_mc, TString::Format("MC (integral scale: %.2f)", t_scale), "lep");

    h_t_data->Draw("E");
    h_t_mc->Draw("SAME");
    t_legend->Draw("SAME");
    c_t->SaveAs("t_data_vs_gen_amp.pdf");
    delete c_t;

    TCanvas *c_m = new TCanvas("c_m", "M4Pi distribution", 800, 600);
    
    const double m_bin_width = get_bin_width(h_m_data);   
    h_m_data->SetLineColor(kBlack);
    h_m_data->SetTitle("");
    h_m_data->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}#pi^{0}#pi^{0}} (GeV)");
    h_m_data->GetYaxis()->SetTitle(Form("Events / %.3f GeV", m_bin_width));
    h_m_data->GetXaxis()->SetTitleSize(0.05);
    h_m_data->GetYaxis()->SetTitleSize(0.05);
    h_m_data->GetXaxis()->SetTitleOffset(0.8);
    h_m_data->GetYaxis()->SetTitleOffset(0.8);
    h_m_mc->SetMarkerStyle(1);
    h_m_mc->SetLineColor(kRed);
    h_m_mc->SetMarkerColor(kRed);
    h_m_data->SetMinimum(0.0);
    h_m_data->SetMaximum(std::max(h_m_data->GetMaximum(), h_m_mc->GetMaximum()) * 1.1);

    TLegend *m_legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    m_legend->AddEntry(h_m_data, "Data", "lep");
    m_legend->AddEntry(h_m_mc, TString::Format("MC (max scale: %.2f)", m_scale), "lep");

    h_m_data->Draw("E");
    h_m_mc->Draw("SAME");
    m_legend->Draw("SAME");
    c_m->SaveAs("m_data_vs_gen_amp.pdf");
    delete c_m;

    return;
}