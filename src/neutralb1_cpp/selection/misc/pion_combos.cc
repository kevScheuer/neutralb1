/**
 * @file pion_combos.cc
 * @author Kevin Scheuer
 * @brief Make some plots of the different proton+pion combos to look for backgrounds
 * 
 */

#include <iostream>
#include <vector>
#include <utility>

#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"

#include "FSMode/FSModeHistogram.h"

#include "../fsroot_setup.cc"


void pion_combos()
{
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(false);
    
    // also show sigmc scale params as usual
    TString dir = "/lustre24/expphy/volatile/halld/home/kscheuer/"
    "FSRoot-skimmed-trees/final-amptools-trees/";
    TString signal_file = dir + "PARA_0_allPeriods_data_signal.root";
    TString background_file = dir + "PARA_0_allPeriods_data_background.root";
    TString mc_signal_file = dir + "PARA_0_allPeriods_ver03.1_mc_signal.root";
    TString mc_background_file = dir + "PARA_0_allPeriods_ver03.1_mc_background.root";

    // AmpTools indices
    int PROTON = 1; 
    int PIO_BACH = 2;
    int PI0_OMEGA = 3;
    int PIPLUS = 4;
    int PIMINUS = 5;

    std::vector<std::pair<double, double>> omega_mass_windows = {
        {1.20, 1.25},
        {1.45, 1.50},
        {1.65, 1.70}
    };

    for (const auto& window : omega_mass_windows) {
        double low = window.first;
        double high = window.second;

        TString omega_mass = TString::Format(
            "MASS(%i, %i, %i, %i)",
            PIO_BACH, PI0_OMEGA, PIPLUS, PIMINUS
        );

        TString selection = TString::Format(
            "(%s > %.3f && %s < %.3f)",
            omega_mass.Data(), low,
            omega_mass.Data(), high
        );    
    
        // ==== proton + pi0 from bachelor =====
        TH1F *h_p_pi0bach_sig_data = FSModeHistogram::getTH1F(
            signal_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIO_BACH),
            "(100, 1.0, 3.0)",
            selection.Data()
        );
        TH1F *h_p_pi0bach_bkg_data = FSModeHistogram::getTH1F( 
            background_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIO_BACH),
            "(100, 1.0, 3.0)",
            (selection+"*weight").Data()
        );
        TH1F *h_p_pi0bach_sig_mc = FSModeHistogram::getTH1F(
            mc_signal_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIO_BACH),
            "(100, 1.0, 3.0)",
            selection.Data()
        );
        TH1F *h_p_pi0bach_bkg_mc = FSModeHistogram::getTH1F( 
            mc_background_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIO_BACH),
            "(100, 1.0, 3.0)",
            (selection+"*weight").Data()
        );
        TH1F *h_p_pi0bach_final_data = (TH1F*)h_p_pi0bach_sig_data->Clone("h_p_pi0bach_final_data");
        h_p_pi0bach_final_data->Add(h_p_pi0bach_bkg_data, -1.0);

        // ==== proton + pi0 from omega =====
        TH1F *h_p_pi0omega_sig_data = FSModeHistogram::getTH1F( // where we typically see Delta+
            signal_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PI0_OMEGA),
            "(100, 1.0, 3.0)",
            selection.Data()
        );
        TH1F *h_p_pi0omega_bkg_data = FSModeHistogram::getTH1F( // where we typically see Delta+
            background_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PI0_OMEGA),
            "(100, 1.0, 3.0)",
            (selection+"*weight").Data()
        );
        TH1F *h_p_pi0omega_sig_mc = FSModeHistogram::getTH1F( // where we typically see Delta+
            mc_signal_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PI0_OMEGA),
            "(100, 1.0, 3.0)",
            selection.Data()
        );
        TH1F *h_p_pi0omega_bkg_mc = FSModeHistogram::getTH1F( // where we typically see Delta+
            mc_background_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PI0_OMEGA),
            "(100, 1.0, 3.0)",
            (selection+"*weight").Data()
        );
        TH1F *h_p_pi0omega_final_data = (TH1F*)h_p_pi0omega_sig_data->Clone("h_p_pi0omega_final_data");
        h_p_pi0omega_final_data->Add(h_p_pi0omega_bkg_data, -1.0);
        
        // ==== proton + pi+ =====
         TH1F *h_p_piplus_sig_data = FSModeHistogram::getTH1F(
            signal_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIPLUS),
            "(100, 1.0, 3.0)",
            selection.Data()
        );
        TH1F *h_p_piplus_bkg_data = FSModeHistogram::getTH1F(
            background_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIPLUS),
            "(100, 1.0, 3.0)",
            (selection+"*weight").Data()
        );
        TH1F *h_p_piplus_sig_mc = FSModeHistogram::getTH1F(
            mc_signal_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIPLUS),
            "(100, 1.0, 3.0)",
            selection.Data()
        );
        TH1F *h_p_piplus_bkg_mc = FSModeHistogram::getTH1F(
            mc_background_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIPLUS),
            "(100, 1.0, 3.0)",
            (selection+"*weight").Data()
        );
        TH1F *h_p_piplus_final_data = (TH1F*)h_p_piplus_sig_data->Clone("h_p_piplus_final_data");
        h_p_piplus_final_data->Add(h_p_piplus_bkg_data, -1.0);

        // ==== proton + pi- =====
         TH1F *h_p_piminus_sig_data = FSModeHistogram::getTH1F(
            signal_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIMINUS),
            "(100, 1.0, 3.0)",
            selection.Data()
        );
        TH1F *h_p_piminus_bkg_data = FSModeHistogram::getTH1F(
            background_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIMINUS),
            "(100, 1.0, 3.0)",
            (selection+"*weight").Data()
        );
        TH1F *h_p_piminus_sig_mc = FSModeHistogram::getTH1F(
            mc_signal_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIMINUS),
            "(100, 1.0, 3.0)",
            selection.Data()
        );
        TH1F *h_p_piminus_bkg_mc = FSModeHistogram::getTH1F(
            mc_background_file,
            NT,
            CATEGORY,
            TString::Format("MASS(%i,%i)", PROTON, PIMINUS),
            "(100, 1.0, 3.0)",
            (selection+"*weight").Data()
        );
        TH1F *h_p_piminus_final_data = (TH1F*)h_p_piminus_sig_data->Clone("h_p_piminus_final_data");
        h_p_piminus_final_data->Add(h_p_piminus_bkg_data, -1.0);


        // scale mc histograms to match data maximum
        h_p_pi0bach_sig_mc->Scale(
            h_p_pi0bach_sig_data->GetMaximum() / h_p_pi0bach_sig_mc->GetMaximum()
        );
        h_p_pi0bach_bkg_mc->Scale(
            h_p_pi0bach_bkg_data->GetMaximum() / h_p_pi0bach_bkg_mc->GetMaximum()
        );
        h_p_pi0omega_sig_mc->Scale(
            h_p_pi0omega_sig_data->GetMaximum() / h_p_pi0omega_sig_mc->GetMaximum()
        );
        h_p_pi0omega_bkg_mc->Scale(
            h_p_pi0omega_bkg_data->GetMaximum() / h_p_pi0omega_bkg_mc->GetMaximum()
        );
        h_p_piplus_sig_mc->Scale(
            h_p_piplus_sig_data->GetMaximum() / h_p_piplus_sig_mc->GetMaximum()
        );
        h_p_piplus_bkg_mc->Scale(
            h_p_piplus_bkg_data->GetMaximum() / h_p_piplus_bkg_mc->GetMaximum()
        );
        h_p_piminus_sig_mc->Scale(
            h_p_piminus_sig_data->GetMaximum() / h_p_piminus_sig_mc->GetMaximum()
        );
        h_p_piminus_bkg_mc->Scale(
            h_p_piminus_bkg_data->GetMaximum() / h_p_piminus_bkg_mc->GetMaximum()
        );

                
        // draw histograms and save as PDF
        TCanvas *c = new TCanvas("c_pion_combos", "Pion Combinations", 800, 400);
        c->Divide(2, 2);        

        double bin_width = h_p_pi0bach_sig_data->GetXaxis()->GetBinWidth(1);

        // change some plot properties for the first histogram to be plotted on each canvas
        h_p_pi0bach_sig_data->GetXaxis()->SetLabelSize(0.05);
        h_p_pi0bach_sig_data->GetYaxis()->SetLabelSize(0.05);
        h_p_pi0omega_sig_data->GetXaxis()->SetLabelSize(0.05);
        h_p_pi0omega_sig_data->GetYaxis()->SetLabelSize(0.05);
        h_p_piplus_sig_data->GetXaxis()->SetLabelSize(0.05);
        h_p_piplus_sig_data->GetYaxis()->SetLabelSize(0.05);
        h_p_piminus_sig_data->GetXaxis()->SetLabelSize(0.05);
        h_p_piminus_sig_data->GetYaxis()->SetLabelSize(0.05);

        h_p_pi0bach_sig_data->GetXaxis()->SetTitleSize(0.06);
        h_p_pi0bach_sig_data->GetYaxis()->SetTitleSize(0.06);
        h_p_pi0omega_sig_data->GetXaxis()->SetTitleSize(0.06);
        h_p_pi0omega_sig_data->GetYaxis()->SetTitleSize(0.06);
        h_p_piplus_sig_data->GetXaxis()->SetTitleSize(0.06);
        h_p_piplus_sig_data->GetYaxis()->SetTitleSize(0.06);
        h_p_piminus_sig_data->GetXaxis()->SetTitleSize(0.06);
        h_p_piminus_sig_data->GetYaxis()->SetTitleSize(0.06);

        h_p_pi0bach_sig_data->GetYaxis()->SetTitleOffset(0.7);
        h_p_pi0omega_sig_data->GetYaxis()->SetTitleOffset(0.7);
        h_p_piplus_sig_data->GetYaxis()->SetTitleOffset(0.7);
        h_p_piminus_sig_data->GetYaxis()->SetTitleOffset(0.7);

        // signal mc signal is kCyan+2
        // background is kRed+1
        // signal mc background is kRed-7
        
        c->cd(1);
        h_p_pi0bach_sig_data->SetTitle(
            TString::Format(";p'#pi^{0}_{bach} inv. mass (GeV);Entries / %.3f GeV", bin_width));        
        h_p_pi0bach_sig_data->SetLineColor(kBlue);
        h_p_pi0bach_bkg_data->SetLineColor(kRed+1);
        h_p_pi0bach_sig_mc->SetLineColorAlpha(kCyan+2, 0.6);
        h_p_pi0bach_sig_mc->SetMarkerColorAlpha(kCyan+2, 0.6);
        h_p_pi0bach_sig_mc->SetMarkerStyle(24);
        h_p_pi0bach_sig_mc->SetMarkerSize(0.4);
        h_p_pi0bach_bkg_mc->SetLineColorAlpha(kRed-7, 0.6);
        h_p_pi0bach_bkg_mc->SetMarkerColorAlpha(kRed-7, 0.6);
        h_p_pi0bach_bkg_mc->SetMarkerStyle(25);
        h_p_pi0bach_bkg_mc->SetMarkerSize(0.4);
        h_p_pi0bach_final_data->SetLineColor(kBlack);
        h_p_pi0bach_final_data->SetMarkerStyle(1);
        h_p_pi0bach_sig_data->Draw("HIST");
        h_p_pi0bach_bkg_data->Draw("HIST SAME");
        h_p_pi0bach_final_data->Draw("E SAME");
        h_p_pi0bach_sig_mc->Draw("E SAME");
        h_p_pi0bach_bkg_mc->Draw("E SAME");        
        
        c->cd(2);
        h_p_pi0omega_sig_data->SetTitle(
            TString::Format(";p'#pi^{0} inv. mass (GeV);Entries / %.3f GeV", bin_width));
        h_p_pi0omega_sig_data->SetLineColor(kBlue);
        h_p_pi0omega_bkg_data->SetLineColor(kRed+1);
        h_p_pi0omega_sig_mc->SetLineColorAlpha(kCyan+2, 0.6);
        h_p_pi0omega_sig_mc->SetMarkerColorAlpha(kCyan+2, 0.6);
        h_p_pi0omega_sig_mc->SetMarkerStyle(24);
        h_p_pi0omega_sig_mc->SetMarkerSize(0.4);
        h_p_pi0omega_bkg_mc->SetLineColorAlpha(kRed-7, 0.6);
        h_p_pi0omega_bkg_mc->SetMarkerColorAlpha(kRed-7, 0.6);
        h_p_pi0omega_bkg_mc->SetMarkerStyle(25);
        h_p_pi0omega_bkg_mc->SetMarkerSize(0.4);
        h_p_pi0omega_final_data->SetLineColor(kBlack);
        h_p_pi0omega_final_data->SetMarkerStyle(1);
        h_p_pi0omega_sig_data->Draw("HIST");
        h_p_pi0omega_bkg_data->Draw("HIST SAME");
        h_p_pi0omega_final_data->Draw("E SAME");
        h_p_pi0omega_sig_mc->Draw("E SAME");
        h_p_pi0omega_bkg_mc->Draw("E SAME");

        
        c->cd(3);
        h_p_piplus_sig_data->SetTitle(
            TString::Format(";p'#pi^{+} inv. mass (GeV);Entries / %.3f GeV", bin_width));
        h_p_piplus_sig_data->SetLineColor(kBlue);
        h_p_piplus_bkg_data->SetLineColor(kRed+1);
        h_p_piplus_sig_mc->SetLineColorAlpha(kCyan+2, 0.6);
        h_p_piplus_sig_mc->SetMarkerColorAlpha(kCyan+2, 0.6);
        h_p_piplus_sig_mc->SetMarkerStyle(24);
        h_p_piplus_sig_mc->SetMarkerSize(0.4);
        h_p_piplus_bkg_mc->SetLineColorAlpha(kRed-7, 0.6);
        h_p_piplus_bkg_mc->SetMarkerColorAlpha(kRed-7, 0.6);
        h_p_piplus_bkg_mc->SetMarkerStyle(25);
        h_p_piplus_bkg_mc->SetMarkerSize(0.4);
        h_p_piplus_final_data->SetLineColor(kBlack);
        h_p_piplus_final_data->SetMarkerStyle(1);

        TLegend *leg = new TLegend(0.6, 0.6, 0.88, 0.88);
        leg->SetFillColorAlpha(kWhite, 0.8);
        leg->SetBorderSize(1);
        leg->AddEntry(h_p_pi0bach_final_data, "Data Total", "lp");
        leg->AddEntry(h_p_pi0bach_sig_data, "Data Signal", "l");
        leg->AddEntry(h_p_pi0bach_bkg_data, "Data Background", "l");
        leg->AddEntry(h_p_pi0bach_sig_mc, "MC Signal", "p");
        leg->AddEntry(h_p_pi0bach_bkg_mc, "MC Background", "p");
        leg->SetTextSize(0.04);

        h_p_piplus_sig_data->Draw("HIST");
        h_p_piplus_bkg_data->Draw("HIST SAME");
        h_p_piplus_final_data->Draw("E SAME");
        h_p_piplus_sig_mc->Draw("E SAME");
        h_p_piplus_bkg_mc->Draw("E SAME");    
        leg->Draw();    
        
        c->cd(4);
        h_p_piminus_sig_data->SetTitle(
            TString::Format(";p'#pi^{-} inv. mass (GeV);Entries / %.3f GeV", bin_width));
        h_p_piminus_sig_data->SetLineColor(kBlue);
        h_p_piminus_bkg_data->SetLineColor(kRed+1);
        h_p_piminus_sig_mc->SetLineColorAlpha(kCyan+2, 0.6);
        h_p_piminus_sig_mc->SetMarkerColorAlpha(kCyan+2, 0.6);
        h_p_piminus_sig_mc->SetMarkerStyle(24);
        h_p_piminus_sig_mc->SetMarkerSize(0.4);
        h_p_piminus_bkg_mc->SetLineColorAlpha(kRed-7, 0.6);
        h_p_piminus_bkg_mc->SetMarkerColorAlpha(kRed-7, 0.6);
        h_p_piminus_bkg_mc->SetMarkerStyle(25);
        h_p_piminus_bkg_mc->SetMarkerSize(0.4);
        h_p_piminus_final_data->SetLineColor(kBlack);
        h_p_piminus_final_data->SetMarkerStyle(1);
        h_p_piminus_sig_data->Draw("HIST");
        h_p_piminus_bkg_data->Draw("HIST SAME");
        h_p_piminus_final_data->Draw("E SAME");
        h_p_piminus_sig_mc->Draw("E SAME");
        h_p_piminus_bkg_mc->Draw("E SAME");        
        
        c->SaveAs(TString::Format("pion_combos_%.2f_%.2f.pdf", low, high));
        delete c;
        delete leg;
    }
    return;
}