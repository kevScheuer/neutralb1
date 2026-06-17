/**
 * @file generate_t_with_acceptance.cc
 * @author Kevin Scheuer
 * @brief Generate t distribution for some slope and multiply by acceptance
 * 
 * This is to get a brief understanding of the t distribution we want to generate for
 * Monte Carlo events will be modified by the acceptance.
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

#include "neutralb1/fit_utils.h"


void generate_t_with_acceptance(double slope, std::string plot_log_with_generated = "true") {

	if (slope <= 0.0) {
		std::cerr << "Slope must be positive for exp(-slope * t). Received: "
				  << slope << std::endl;
		return;
	}

	gStyle->SetOptStat(0);

	const int n_bins = 50;
	const double t_min = 0.0;
	const double t_max = 1.0;

	// Data files (for observed t spectrum)
	TString data_dir = "/home/kscheuer/work/neutralb1/data/FSRoot/GlueX/";
	TString data_signal = data_dir + "allPeriods_data_signal.root";
	TString data_sideband = data_dir + "allPeriods_data_background.root";

	// Phase-space files (for acceptance efficiency)
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

	// Build observed data histogram (signal - sideband)
	TH1F *h_t_signal = new TH1F("h_t_signal", "", n_bins, t_min, t_max);
	TH1F *h_t_sideband = new TH1F("h_t_sideband", "", n_bins, t_min, t_max);
	tree_data_signal->Draw("t>>h_t_signal", "", "goff");
	tree_data_sideband->Draw("t>>h_t_sideband", "weight", "goff");

	TH1F *h_t_data = (TH1F*)h_t_signal->Clone("h_t_data");
	h_t_data->Add(h_t_sideband, -1.0);

	// Build acceptance efficiency histogram from phase-space MC
	TH1F *h_t_acc = new TH1F("h_t_acc", "", n_bins, t_min, t_max);
	TH1F *h_t_gen_ps = new TH1F("h_t_gen_ps", "", n_bins, t_min, t_max);
	tree_acc->Draw("t>>h_t_acc", "weight", "goff");
	tree_gen->Draw("t>>h_t_gen_ps", "", "goff");

	TH1F *h_t_eff = (TH1F*)h_t_acc->Clone("h_t_eff");
	h_t_eff->SetDirectory(0);
	h_t_eff->Divide(h_t_gen_ps);

	// Build generated exponential template and acceptance-modified template.
	TH1F *h_t_generated = new TH1F("h_t_generated", "", n_bins, t_min, t_max);
	TH1F *h_t_accepted = new TH1F("h_t_accepted", "", n_bins, t_min, t_max);

	TF1 f_expo_model("f_expo_model", "exp(-[0]*x)", t_min, t_max);
	f_expo_model.SetParameter(0, slope);

	for (int i = 1; i <= n_bins; ++i) {
		const double t_center = h_t_generated->GetBinCenter(i);
		const double gen_val = f_expo_model.Eval(t_center);
		const double eff_val = h_t_eff->GetBinContent(i);

		h_t_generated->SetBinContent(i, gen_val);
		h_t_accepted->SetBinContent(i, gen_val * eff_val);
	}

	// Scale model histograms so accepted integral matches observed data integral.
	const double data_integral = h_t_data->Integral();
	const double accepted_integral = h_t_accepted->Integral();
	if (accepted_integral > 0.0 && data_integral > 0.0) {
		const double scale = data_integral / accepted_integral;
		h_t_generated->Scale(scale);
		h_t_accepted->Scale(scale);
	}

	const double bin_width = get_bin_width(h_t_data);

	h_t_data->SetTitle("");
	h_t_data->GetXaxis()->SetTitle("-t (GeV)^{2}");
	h_t_data->GetYaxis()->SetTitle(TString::Format("Events / %.3f GeV^{2}", bin_width));
	h_t_data->GetXaxis()->SetTitleSize(0.05);
	h_t_data->GetYaxis()->SetTitleSize(0.05);
	h_t_data->GetXaxis()->SetTitleOffset(0.8);
	h_t_data->GetYaxis()->SetTitleOffset(0.8);
	h_t_data->SetMarkerStyle(8);
	h_t_data->SetMarkerSize(1.0);
	h_t_data->SetMarkerColor(kBlack);
	h_t_data->SetLineColor(kBlack);
	h_t_data->SetMinimum(0.0);

	h_t_generated->SetLineColor(kBlue + 1);
	h_t_generated->SetLineWidth(3);

	h_t_accepted->SetLineColor(kRed + 1);
	h_t_accepted->SetLineWidth(3);

	double y_max;
	double y_min;
    if (plot_log_with_generated == "true") {
        y_min = 0.01;
        y_max = std::max(
            h_t_data->GetMaximum(),
            std::max(h_t_generated->GetMaximum(), h_t_accepted->GetMaximum())
        );
    }
    else {
        y_min = 0.0;
        y_max = h_t_data->GetMaximum();
    }
    h_t_data->SetMaximum(1.20 * y_max);    
    h_t_data->SetMinimum(y_min);

	TCanvas *c = new TCanvas("c_t_gen_acc", "Generated vs accepted vs data", 900, 700);
	h_t_data->Draw("E");
	h_t_accepted->Draw("HIST SAME");
	h_t_data->Draw("E SAME");

	TLegend *leg = new TLegend(0.52, 0.68, 0.88, 0.88);
	leg->SetFillColorAlpha(kWhite, 0.85);
	leg->SetTextSize(0.03);
	leg->AddEntry(h_t_data, "GlueX-I data", "lep");
	
	leg->AddEntry(h_t_accepted, "Accepted", "l");
	leg->Draw();

	std::cout << "Input slope: " << slope << std::endl;
	std::cout << "Data integral: " << h_t_data->Integral() << std::endl;
	
	std::cout << "Accepted integral (scaled): " << h_t_accepted->Integral() << std::endl;    

	
    TString out_name;
    if (plot_log_with_generated == "true") {
        h_t_generated->Draw("HIST SAME");
        leg->AddEntry(h_t_generated, TString::Format("Generated: exp(-%.3f t)", slope), "l");
        std::cout << "Generated integral (scaled): " << h_t_generated->Integral() << std::endl;
        c->SetLogy();
        out_name = TString::Format("t_generated_with_acceptance_slope_log_%.3f.pdf", slope);	
    }
    else
    {
        out_name = TString::Format("t_generated_with_acceptance_slope_%.3f.pdf", slope);	
    }

	c->SaveAs(out_name);

}