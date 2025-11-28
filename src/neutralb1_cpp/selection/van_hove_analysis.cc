/**
 * @file van_hove_analysis.cc
 * @author Kevin Scheuer, Amy Schertz
 * @brief Plot van Hove distributions for data and MC
 *
 * Once again, Amy Schertz deserves credit for the original authorship of this code.
 */

#include <iostream>
#include <tuple>
#include <map>
#include <vector>

#include "TArrow.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TStyle.h"

#include "FSBasic/FSHistogram.h"
#include "FSBasic/FSCut.h"
#include "FSMode/FSModeHistogram.h"

#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"

// CONSTANTS
// AmpTools indices for particles in the final state
const TString TARGET = "GLUEXTARGET";
const TString BEAM = "B";
const int PROTON = 1;
const int PI0_BACH = 2;
const int PI0_OM = 3;
const int PIPLUS = 4;
const int PIMINUS = 5;

const TString M_OMEGA = TString::Format("MASS(%d,%d,%d)", PIPLUS, PIMINUS, PI0_OM);
const TString M_OMEGA_PI0 = TString::Format("MASS(%d,%d,%d,%d)", PIPLUS, PIMINUS, PI0_OM, PI0_BACH);

// forward declarations
void plot_relations(
    TString NT,
    TString CATEGORY,
    TString input_signal,
    TString input_sideband,
    bool mc);
void plot_van_hove(
    TString NT,
    TString CATEGORY,
    TString input_signal,
    TString input_sideband,
    bool mc);
void plot_pi_correlations(
    TString NT, 
    TString CATEGORY, 
    TString input_signal, 
    TString input_sideband, 
    bool mc,
    bool apply_delta_cut = false);
void plot_pi_cuts(
    TString NT,
    TString CATEGORY,
    TString input_signal,
    TString input_sideband,
    bool mc);

void van_hove_analysis(
    bool mc = false,
    bool plot_relations_flag = false,
    bool plot_van_hove_flag = false,
    bool plot_pi0_corr_flag = false,
    bool plot_pi0_cuts_flag = false)
{
    gStyle->SetPalette(kBird);

    TString input_signal;
    TString input_sideband;
    if (mc)
    {
        std::cout << "MC" << "\n";
        input_signal = "/lustre24/expphy/volatile/halld/home/kscheuer/"
                       "FSRoot-skimmed-trees/reordered/"
                       "tree_pi0pi0pippim__B4_reordered_SKIM_allPeriods_ver03.1_mc_signal.root";
        input_sideband = "/lustre24/expphy/volatile/halld/home/kscheuer/"
                         "FSRoot-skimmed-trees/reordered/"
                         "tree_pi0pi0pippim__B4_reordered_SKIM_allPeriods_ver03.1_mc_sideband.root";
    }
    else
    {
        input_signal = "/lustre24/expphy/volatile/halld/home/kscheuer/"
                       "FSRoot-skimmed-trees/reordered/"
                       "tree_pi0pi0pippim__B4_reordered_SKIM_allPeriods_data_signal.root";
        input_sideband = "/lustre24/expphy/volatile/halld/home/kscheuer/"
                         "FSRoot-skimmed-trees/reordered/"
                         "tree_pi0pi0pippim__B4_reordered_SKIM_allPeriods_data_sideband.root";
    }
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(false);

    FSCut::defineCut("broad", "cut==1"); // uses special branch to implement all broad cuts

    if (plot_relations_flag)
        plot_relations(NT, CATEGORY, input_signal, input_sideband, mc);
    if (plot_van_hove_flag)
        plot_van_hove(NT, CATEGORY, input_signal, input_sideband, mc);
    if (plot_pi0_corr_flag)
    {
        plot_pi_correlations(NT, CATEGORY, input_signal, input_sideband, mc, false);
        plot_pi_correlations(NT, CATEGORY, input_signal, input_sideband, mc, true);
    }
    if (plot_pi0_cuts_flag)
        plot_pi_cuts(NT, CATEGORY, input_signal, input_sideband, mc);

    return;
}

/**
 * @brief Plot helicity angle vs omega pi0 and a dalitz plot to look for baryon resonances
 *
 * @param[in] NT FSRoot tree name
 * @param[in] CATEGORY FSRoot mode category
 * @param[in] input_signal signal input file
 * @param[in] input_sideband sideband input file
 * @param[in] mc whether the input files are from Monte Carlo simulation
 */
void plot_relations(
    TString NT,
    TString CATEGORY,
    TString input_signal,
    TString input_sideband,
    bool mc)
{
    TCanvas *c_corr = new TCanvas("ccorr", "ccorr", 800, 600);

    // Under X-> omega pi0 p', we need to input into HELCOSTHETA(omega;pi0;p')
    TH2F *h_corr[2];
    h_corr[0] = FSModeHistogram::getTH2F(
        input_signal,
        NT,
        CATEGORY,
        TString::Format(
            "HELCOSTHETA(%d,%d,%d;%d;%d):%s",
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON, M_OMEGA_PI0.Data()),
        "(100, 1.0, 2.0, 200, -1.0, 1.0)",
        "CUT(broad)");
    h_corr[1] = FSModeHistogram::getTH2F(
        input_sideband,
        NT,
        CATEGORY,
        TString::Format(
            "HELCOSTHETA(%d,%d,%d;%d;%d):%s",
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON, M_OMEGA_PI0.Data()),
        "(100, 1.0, 2.0, 200, -1.0, 1.0)",
        "CUT(broad)");
    h_corr[0]->Add(h_corr[1], -1); // signal - sideband

    h_corr[0]->GetXaxis()->SetTitle("#omega#pi^{0} inv. mass (GeV)");
    h_corr[0]->GetYaxis()->SetTitle("cos #theta");
    h_corr[0]->SetTitle("");
    h_corr[0]->Draw("colz");
    c_corr->SaveAs(
        TString::Format(
            "costhetah_mass_corr%s.pdf",
            mc ? "_mc" : "_data"));

    // correlation plot of cos theta against p' pi0 mass
    TCanvas *c_corr_ppi0 = new TCanvas("ccorr_ppi0", "ccorr_ppi0", 800, 600);
    TH2F *h_corr_ppi0[2];
    h_corr_ppi0[0] = FSModeHistogram::getTH2F(
        input_signal,
        NT,
        CATEGORY,
        TString::Format(
            "HELCOSTHETA(%d,%d,%d;%d;%d):MASS(%d,%d)",
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON,
            PROTON, PI0_BACH),
        "(200, 1.0, 3.0, 200, -1.0, 1.0)",
        "CUT(broad)");
    h_corr_ppi0[1] = FSModeHistogram::getTH2F(
        input_sideband,
        NT,
        CATEGORY,
        TString::Format(
            "HELCOSTHETA(%d,%d,%d;%d;%d):MASS(%d,%d)",
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON,
            PROTON, PI0_BACH),
        "(200, 1.0, 3.0, 200, -1.0, 1.0)",
        "CUT(broad)");
    h_corr_ppi0[0]->Add(h_corr_ppi0[1], -1); // signal - sideband
    h_corr_ppi0[0]->GetXaxis()->SetTitle("p' #pi^{0} inv. mass (GeV)");
    h_corr_ppi0[0]->GetYaxis()->SetTitle("cos #theta");
    h_corr_ppi0[0]->SetTitle("");
    h_corr_ppi0[0]->Draw("colz");
    c_corr_ppi0->SaveAs(
        TString::Format(
            "costhetah_mass_ppi0_corr%s.pdf",
            mc ? "_mc" : "_data"));

    // now dalitz plot of M^2(p' pi0) vs M^2(omega pi0)
    TCanvas *c_dalitz = new TCanvas("cdalitz", "cdalitz", 800, 600);
    TH2F *h_dalitz[2];
    h_dalitz[0] = FSModeHistogram::getTH2F(
        input_signal,
        NT,
        CATEGORY,
        TString::Format(
            "MASS2(%d,%d):MASS2(%d,%d,%d, %d)",
            PI0_BACH, PROTON,
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH),
        "(300, 1.0, 4.0, 900, 1.0, 9.0)",
        "CUT(broad)");
    h_dalitz[1] = FSModeHistogram::getTH2F(
        input_sideband,
        NT,
        CATEGORY,
        TString::Format(
            "MASS2(%d,%d):MASS2(%d,%d,%d, %d)",
            PI0_BACH, PROTON,
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH),
        "(300, 1.0, 4.0, 900, 1.0, 9.0)",
        "CUT(broad)");
    h_dalitz[0]->Add(h_dalitz[1], -1); // signal - sideband
    h_dalitz[0]->GetXaxis()->SetTitle("M^{2}(#omega#pi^{0}) (GeV^{2})");
    h_dalitz[0]->GetYaxis()->SetTitle("M^{2}(p' #pi^{0}) (GeV^{2})");
    h_dalitz[0]->SetTitle("");
    h_dalitz[0]->Draw("colz");
    c_dalitz->SaveAs(
        TString::Format(
            "dalitz_omegapi0_ppi0%s.pdf",
            mc ? "_mc" : "_data"));

    // FSHistogram::dumpHistogramCache();

    return;
}

/**
 * @brief plot van Hove longitudinal phasespace with / without cut on p'pi0 mass 
 * 
 * @param[in] NT FSRoot tree name
 * @param[in] CATEGORY FSRoot mode category
 * @param[in] input_signal signal input file
 * @param[in] input_sideband sideband input file
 * @param[in] mc whether the input files are from Monte Carlo simulation
 */
void plot_van_hove(
    TString NT,
    TString CATEGORY,
    TString input_signal,
    TString input_sideband,
    bool mc)
{

    // ==== Plot the longitudinal phasespace in Van Hove XY with normal cuts ==== 
    // Three particle final state (omega, pi_bach, and recoil proton) is fed into
    // van hove function
    TH2F *h_van_hove[2];
    h_van_hove[0] = FSModeHistogram::getTH2F(
        input_signal,
        NT,
        CATEGORY,
        TString::Format(
            "VANHOVEY(%d,%d,%d;%d;%d):VANHOVEX(%d,%d,%d;%d;%d)",
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON,
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON),
        "(400, -4.0, 4.0, 400, -4.0, 4.0)",
        "CUT(broad)");
    h_van_hove[1] = FSModeHistogram::getTH2F(
        input_sideband,
        NT,
        CATEGORY,
        TString::Format(
            "VANHOVEY(%d,%d,%d;%d;%d):VANHOVEX(%d,%d,%d;%d;%d)",
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON,
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON),
        "(400, -4.0, 4.0, 400, -4.0, 4.0)",
        "CUT(broad)");

    // draw lines for van Hove regions
    TF1 *vHLine_pi0 = new TF1("vhl1", "-1.0*x*TMath::Tan( TMath::Pi() / 3.0 )", -4, 4);
    TF1 *vHLine_omega = new TF1("vhl2", "0*x", -4, 4.1);
    TF1 *vHLine_proton = new TF1("vhl3", "x*TMath::Tan( TMath::Pi() / 3.0 )", -4, 4);
    vHLine_pi0->SetLineColor(kGray+1);
    vHLine_pi0->SetLineWidth(2);
    vHLine_omega->SetLineColor(kGray+1);
    vHLine_omega->SetLineWidth(2);
    vHLine_proton->SetLineColor(kGray+1);
    vHLine_proton->SetLineWidth(2);

    // label each region
    TLatex *omega_label = new TLatex(3.7, 0.2, "#omega");
    omega_label->SetTextColor(kGray+1);
    TArrow *arrow_omega = new TArrow(3.5, 0.0, 3.5, 1.0, 0.01, "|>");
    arrow_omega->SetLineWidth(2);
    arrow_omega->SetLineColor(kGray+1);
    arrow_omega->SetFillColor(kGray+1);
    arrow_omega->SetFillStyle(1001);

    TLatex *baryon_label = new TLatex(-1.5, -3.5, "p'");
    baryon_label->SetTextColor(kGray+1);
    TArrow *arrow_baryon = new TArrow(
        -1.5, vHLine_proton->Eval(-1.5), 
        -1.1, -3.0, 
        0.01, "|>");
    arrow_baryon->SetLineWidth(2);
    arrow_baryon->SetLineColor(kGray+1);
    arrow_baryon->SetFillColor(kGray+1);
    arrow_baryon->SetFillStyle(1001);

    TLatex *bachelor_label = new TLatex(-2.5, 3.2, "#pi^{0}");
    bachelor_label->SetTextColor(kGray+1);
    TArrow *arrow_bachelor = new TArrow(
        -2.0, vHLine_pi0->Eval(-2.0), 
        -2.4, 2.9, 
        0.01, "|>");
    arrow_bachelor->SetLineWidth(2);
    arrow_bachelor->SetLineColor(kGray+1);
    arrow_bachelor->SetFillColor(kGray+1);
    arrow_bachelor->SetFillStyle(1001);


    h_van_hove[0]->Add(h_van_hove[1], -1); // signal - sideband
    h_van_hove[0]->GetXaxis()->SetTitle("X");
    h_van_hove[0]->GetYaxis()->SetTitle("Y");
    TCanvas *c_van_hove = new TCanvas("cvanhove", "cvanhove", 800, 600);
    h_van_hove[0]->SetTitle("");
    h_van_hove[0]->Draw("colz");
    vHLine_pi0->Draw("SAME");
    vHLine_omega->Draw("SAME");
    vHLine_proton->Draw("SAME");
    omega_label->Draw();
    arrow_omega->Draw();
    baryon_label->Draw();
    arrow_baryon->Draw();
    bachelor_label->Draw();
    arrow_bachelor->Draw();
    c_van_hove->SaveAs(
        TString::Format(
            "van_hove%s.pdf",
            mc ? "_mc" : "_data"));

    // ==== Plot the longitudinal phasespace in Van Hove XY, selecting Delta+ ==== 
    FSCut::defineCut("delta", TString::Format("MASS(%d,%d) < 1.4", PI0_BACH, PROTON));
    TH2F *h_van_hove_delta[2];
    h_van_hove_delta[0] = FSModeHistogram::getTH2F(
        input_signal,
        NT,
        CATEGORY,
        TString::Format(
            "VANHOVEY(%d,%d,%d;%d;%d):VANHOVEX(%d,%d,%d;%d;%d)",
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON,
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON),
        "(400, -4.0, 4.0, 400, -4.0, 4.0)",
        "CUT(delta,broad)");
    h_van_hove_delta[1] = FSModeHistogram::getTH2F(
        input_sideband,
        NT,
        CATEGORY,
        TString::Format(
            "VANHOVEY(%d,%d,%d;%d;%d):VANHOVEX(%d,%d,%d;%d;%d)",
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON,
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON),
        "(400, -4.0, 4.0, 400, -4.0, 4.0)",
        "CUT(delta,broad)");

    h_van_hove_delta[0]->Add(h_van_hove_delta[1], -1); // signal - sideband
    h_van_hove_delta[0]->GetXaxis()->SetTitle("X");
    h_van_hove_delta[0]->GetYaxis()->SetTitle("Y");
    TCanvas *c_van_hove_delta = new TCanvas("cvanhove_delta", "cvanhove_delta", 800, 600);
    h_van_hove_delta[0]->SetTitle("");
    h_van_hove_delta[0]->Draw("colz");
    vHLine_pi0->Draw("SAME");
    vHLine_omega->Draw("SAME");
    vHLine_proton->Draw("SAME");
    omega_label->Draw();
    arrow_omega->Draw();
    baryon_label->Draw();
    arrow_baryon->Draw();
    bachelor_label->Draw();
    arrow_bachelor->Draw();
    c_van_hove_delta->SaveAs(
        TString::Format(
            "van_hove_delta%s.pdf",
            mc ? "_mc" : "_data"));
    return;
}

/**
 * @brief plot correlations of bachelor pi0 momentum to other distributions
 * 
 * Other distributions include omega pi0 mass, proton pi0 mass, and cos theta
 * 
 * @param[in] NT FSRoot tree name
 * @param[in] CATEGORY FSRoot mode category
 * @param[in] input_signal signal input file
 * @param[in] input_sideband sideband input file
 * @param[in] mc whether the input files are from Monte Carlo simulation
 * @param[in] apply_delta_cut whether to apply a cut on the p' pi0 mass to select Delta+ resonance
 */
void plot_pi_correlations(
    TString NT, 
    TString CATEGORY, 
    TString input_signal, 
    TString input_sideband, 
    bool mc,
    bool apply_delta_cut = false)
{
    TCanvas *c_pi_corr = new TCanvas("cpi_corr", "cpi_corr", 800, 600);

    FSCut::defineCut("delta", TString::Format("MASS(%d,%d) < 1.4", PI0_BACH, PROTON));
    TString cut_string = "CUT(broad)";
    TString file_suffix = "";
    if (apply_delta_cut)
    {
        cut_string = "CUT(delta,broad)";
        file_suffix = "_delta";
    }

    // correlation of bachelor pi0 momentum vs omega pi0 mass
    TH2F *h_pi_corr_omegapi0[2];
    h_pi_corr_omegapi0[0] = FSModeHistogram::getTH2F(
        input_signal,
        NT,
        CATEGORY,
        TString::Format(
            "MOMENTUMZBOOST(%d;%s,%s):%s",
            PI0_BACH, BEAM.Data(), TARGET.Data(), M_OMEGA_PI0.Data()),
        "(200, 1.0, 2.0, 200, -0.5, 1.5)",
        cut_string);
    h_pi_corr_omegapi0[1] = FSModeHistogram::getTH2F(
        input_sideband,
        NT,
        CATEGORY,
        TString::Format(
            "MOMENTUMZBOOST(%d;%s,%s):%s",
            PI0_BACH, BEAM.Data(), TARGET.Data(), M_OMEGA_PI0.Data()),
        "(200, 1.0, 2.0, 200, -0.5, 1.5)",
        cut_string);
    h_pi_corr_omegapi0[0]->Add(h_pi_corr_omegapi0[1], -1); // signal - sideband
    h_pi_corr_omegapi0[0]->GetXaxis()->SetTitle("#omega#pi^{0} inv. mass (GeV)");
    h_pi_corr_omegapi0[0]->GetYaxis()->SetTitle("P_{z}^{CM}(#pi^{0}) (GeV)");
    h_pi_corr_omegapi0[0]->SetTitle("");
    h_pi_corr_omegapi0[0]->Draw("colz");
    c_pi_corr->SaveAs(
        TString::Format(
            "pi0_pz_omegapi0_corr%s%s.pdf",
            mc ? "_mc" : "_data",
            file_suffix.Data()));

    // correlation of bachelor pi0 momentum vs proton pi0 mass
    TCanvas *c_pi_corr_ppi0 = new TCanvas("cpi_corr_ppi0", "cpi_corr_ppi0", 800, 600);
    TH2F *h_pi_corr_ppi0[2];
    h_pi_corr_ppi0[0] = FSModeHistogram::getTH2F(
        input_signal,
        NT,
        CATEGORY,
        TString::Format(
            "MOMENTUMZBOOST(%d;%s,%s):MASS(%d,%d)",
            PI0_BACH, BEAM.Data(), TARGET.Data(),
            PROTON, PI0_BACH),
        "(200, 1.0, 3.0, 200, -0.5, 1.5)",
        cut_string);
    h_pi_corr_ppi0[1] = FSModeHistogram::getTH2F(
        input_sideband,
        NT,
        CATEGORY,
        TString::Format(
            "MOMENTUMZBOOST(%d;%s,%s):MASS(%d,%d)",
            PI0_BACH, BEAM.Data(), TARGET.Data(),
            PROTON, PI0_BACH),
        "(200, 1.0, 3.0, 200, -0.5, 1.5)",
        cut_string);
    h_pi_corr_ppi0[0]->Add(h_pi_corr_ppi0[1], -1); // signal - sideband
    h_pi_corr_ppi0[0]->GetXaxis()->SetTitle("p' #pi^{0} inv. mass (GeV)");
    h_pi_corr_ppi0[0]->GetYaxis()->SetTitle("P_{z}^{CM}(#pi^{0}) (GeV)");
    h_pi_corr_ppi0[0]->SetTitle("");
    h_pi_corr_ppi0[0]->Draw("colz");
    c_pi_corr_ppi0->SaveAs(
        TString::Format(
            "pi0_pz_ppi0_corr%s%s.pdf",
            mc ? "_mc" : "_data",
            file_suffix.Data()));
    
    // correlation of bachelor pi0 momentum vs cos theta
    TCanvas *c_pi_corr_costheta = new TCanvas("cpi_corr_costheta", "cpi_corr_costheta", 800, 600);
    TH2F *h_pi_corr_costheta[2];
    h_pi_corr_costheta[0] = FSModeHistogram::getTH2F(
        input_signal,
        NT,
        CATEGORY,
        TString::Format(
            "MOMENTUMZBOOST(%d;%s,%s):HELCOSTHETA(%d,%d,%d;%d;%d)",
            PI0_BACH, BEAM.Data(), TARGET.Data(),
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON),
        "(200, -1.0, 1.0, 200, -0.5, 1.5)",
        cut_string);
    h_pi_corr_costheta[1] = FSModeHistogram::getTH2F(
        input_sideband,
        NT,
        CATEGORY,
        TString::Format(
            "MOMENTUMZBOOST(%d;%s,%s):HELCOSTHETA(%d,%d,%d;%d;%d)",
            PI0_BACH, BEAM.Data(), TARGET.Data(),
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON),
        "(200, -1.0, 1.0, 200, -0.5, 1.5)",
        cut_string);
    h_pi_corr_costheta[0]->Add(h_pi_corr_costheta[1], -1); // signal - sideband
    h_pi_corr_costheta[0]->GetXaxis()->SetTitle("cos #theta");
    h_pi_corr_costheta[0]->GetYaxis()->SetTitle("P_{z}^{CM}(#pi^{0}) (GeV)");
    h_pi_corr_costheta[0]->SetTitle("");
    h_pi_corr_costheta[0]->Draw("colz");
    c_pi_corr_costheta->SaveAs(
        TString::Format(
            "pi0_pz_costheta_corr%s%s.pdf",
            mc ? "_mc" : "_data",
            file_suffix.Data()));

    return;
}


void plot_pi_cuts(
    TString NT,
    TString CATEGORY,
    TString input_signal,
    TString input_sideband,
    bool mc)
{
    // cuts on pi0 momentum to loop over, from least to most restrictive
    std::map<float, Int_t> pi0_mom_cut_map = {
        {-1.0, kGray+1}, // effectively no cut
        {-0.4, kPink-7}, // colors become lighter
        {-0.3, kPink+4},
        {-0.2, kPink+3},
        {-0.1, kPink+2},
        { 0.0, kPink+1}
    };
    
    // setup containers for each type of histogram
    std::map<float, TH1F*> h_costheta_map;
    std::map<float, TH1F*> h_omegapi0_mass_map;
    std::map<float, TH1F*> h_ppi0_mass_map;
    std::map<float, TH2F*> h_costheta_phi_map;    

    for (std::map<float, Int_t>::iterator it = pi0_mom_cut_map.begin(); it != pi0_mom_cut_map.end(); ++it)
    {
        float pi0_cut = it->first;
        Int_t color = it->second;  
        
        // define the cut and put into common string to use for histograms
        FSCut::defineCut(
            TString::Format("pi0_%1.1f", pi0_cut),
            TString::Format(
                "MOMENTUMZBOOST(%d;%s,%s) > %f",
                PI0_BACH, BEAM.Data(), TARGET.Data(), pi0_cut)
        );
        TString cut_string = TString::Format(
            "CUT(broad, %s)",
            TString::Format("pi0_%1.1f", pi0_cut).Data()
        );
        
        // cut effect on cos theta
        TH1F *h_costheta[2];
        h_costheta[0] = FSModeHistogram::getTH1F(
            input_signal,
            NT,
            CATEGORY,
            TString::Format(
                "HELCOSTHETA(%d,%d,%d;%d;%d)",
                PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON),
            "(200, -1.0, 1.0)",
            cut_string);       
        h_costheta[1] = FSModeHistogram::getTH1F(
            input_sideband,
            NT,
            CATEGORY,
            TString::Format(
                "HELCOSTHETA(%d,%d,%d;%d;%d)",
                PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON),
            "(200, -1.0, 1.0)",
            cut_string);
        h_costheta[0]->Add(h_costheta[1], -1);
        h_costheta[0]->SetLineColor(kBlack);
        h_costheta[0]->SetFillColor(color);
        h_costheta[0]->SetFillStyle(1001);
        h_costheta[0]->SetMarkerSize(1);
        h_costheta_map[pi0_cut] = h_costheta[0];

        // cut effect on omega pi0 mass
        TH1F *h_omegapi0_mass[2];
        h_omegapi0_mass[0] = FSModeHistogram::getTH1F(
            input_signal,
            NT,
            CATEGORY,
            M_OMEGA_PI0,
            "(100, 1.0, 2.0)",
            cut_string);
        h_omegapi0_mass[1] = FSModeHistogram::getTH1F(
            input_sideband,
            NT,
            CATEGORY,
            M_OMEGA_PI0,
            "(100, 1.0, 2.0)",
            cut_string);
        h_omegapi0_mass[0]->Add(h_omegapi0_mass[1], -1);
        h_omegapi0_mass[0]->SetLineColor(kBlack);
        h_omegapi0_mass[0]->SetFillColor(color);
        h_omegapi0_mass[0]->SetFillStyle(1001);
        h_omegapi0_mass[0]->SetMarkerSize(1);
        h_omegapi0_mass_map[pi0_cut] = h_omegapi0_mass[0];

        // cut effect on p' pi0 mass
        TH1F *h_ppi0_mass[2];
        h_ppi0_mass[0] = FSModeHistogram::getTH1F(
            input_signal,
            NT,
            CATEGORY,
            TString::Format("MASS(%d,%d)", PROTON, PI0_BACH),
            "(200, 1.0, 3.0)",
            cut_string);
        h_ppi0_mass[1] = FSModeHistogram::getTH1F(
            input_sideband,
            NT,
            CATEGORY,
            TString::Format("MASS(%d,%d)", PROTON, PI0_BACH),
            "(200, 1.0, 3.0)",
            cut_string);
        h_ppi0_mass[0]->Add(h_ppi0_mass[1], -1);
        h_ppi0_mass[0]->SetLineColor(kBlack);
        h_ppi0_mass[0]->SetFillColor(color);
        h_ppi0_mass[0]->SetFillStyle(1001);
        h_ppi0_mass[0]->SetMarkerSize(1);
        h_ppi0_mass_map[pi0_cut] = h_ppi0_mass[0];

        // cut effect on helcostheta vs helphi
        TH2F *h_costheta_phi[2];
        h_costheta_phi[0] = FSModeHistogram::getTH2F(
            input_signal,
            NT,
            CATEGORY,
            TString::Format( 
                "HELCOSTHETA(%d,%d,%d;%d;%d):HELPHI(%d,%d,%d;%d;%d;%s)",                
                PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON,
                PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON, BEAM.Data()),
            "(200, -3.14, 3.14, 200, 1.0, 1.0)",
            cut_string);
        h_costheta_phi[1] = FSModeHistogram::getTH2F(
            input_sideband,
            NT,
            CATEGORY,
            TString::Format( 
                "HELCOSTHETA(%d,%d,%d;%d;%d):HELPHI(%d,%d,%d;%d;%d;%s)",                
                PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON,
                PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON, BEAM.Data()),
            "(200, -3.14, 3.14, 200, 1.0, 1.0)",
            cut_string);
        h_costheta_phi[0]->Add(h_costheta_phi[1], -1);
        h_costheta_phi[0]->SetTitle("");
        h_costheta_phi_map[pi0_cut] = h_costheta_phi[0];        
    }

    // Now that all histograms are filled, plot them
    TCanvas *c_pi0_cut = new TCanvas("cpi0_cut", "cpi0_cut", 800, 600);
    double bin_width;

    // Plot cos theta
    TLegend *leg_costheta = new TLegend(0.18, 0.58, 0.51, 0.89);
    leg_costheta->SetTextSize(0.03);
    leg_costheta->AddEntry(h_costheta_map.begin()->second, "No cut", "f"); // special legend entry for no cut
    h_costheta_map.begin()->second->SetTitle("");
    h_costheta_map.begin()->second->GetXaxis()->SetTitle("cos #theta");
    bin_width = get_bin_width(h_costheta_map.begin()->second);
    h_costheta_map.begin()->second->GetYaxis()->SetTitle(TString::Format("Events / %.3f", bin_width));
    h_costheta_map.begin()->second->SetMinimum(0.0);
    h_costheta_map.begin()->second->Draw("HIST");
    // start map at second element to avoid re-drawing no cut
    std::map<float, TH1F*>::iterator it = h_costheta_map.begin();
    ++it;
    for (; it != h_costheta_map.end(); ++it)
    {
        it->second->Draw("HIST SAME");
        leg_costheta->AddEntry(
            it->second,
            TString::Format("P_{z}^{CM}(#pi^{0}) > %1.1f GeV", it->first), "f");
    }
    leg_costheta->Draw();
    c_pi0_cut->SaveAs(
        TString::Format(
            "pi0_mom_cut_costheta%s.pdf",
            mc ? "_mc" : "_data"));
    c_pi0_cut->Clear();

    // Plot omega pi0 mass
    TLegend *leg_omegapi0 = new TLegend(0.5, 0.6, 0.8, 0.88);
    leg_omegapi0->SetTextSize(0.03);
    leg_omegapi0->AddEntry(h_omegapi0_mass_map.begin()->second, "No Cut", "f");
    h_omegapi0_mass_map.begin()->second->SetTitle("");
    h_omegapi0_mass_map.begin()->second->GetXaxis()->SetTitle("#omega#pi^{0} inv. mass (GeV)");
    bin_width = get_bin_width(h_omegapi0_mass_map.begin()->second);
    h_omegapi0_mass_map.begin()->second->GetYaxis()->SetTitle(TString::Format("Events / %.3f", bin_width));
    h_omegapi0_mass_map.begin()->second->SetMinimum(0.0);
    h_omegapi0_mass_map.begin()->second->Draw("HIST");
    it = h_omegapi0_mass_map.begin();
    ++it;
    for (; it != h_omegapi0_mass_map.end(); ++it)
    {
        it->second->Draw("HIST SAME");
        leg_omegapi0->AddEntry(
            it->second,
            TString::Format("P_{z}^{CM}(#pi^{0}) > %1.1f GeV", it->first), "f");
    }
    leg_omegapi0->Draw();
    c_pi0_cut->SaveAs(
        TString::Format(
            "pi0_mom_cut_omegapi0_mass%s.pdf",
            mc ? "_mc" : "_data"));
    c_pi0_cut->Clear();

    // Plot p' pi0 mass
    TLegend *leg_ppi0 = new TLegend(0.68, 0.7, 0.98, 0.98);
    leg_ppi0->SetTextSize(0.03);
    leg_ppi0->AddEntry(h_ppi0_mass_map.begin()->second, "No Cut", "f");
    h_ppi0_mass_map.begin()->second->SetTitle("");
    h_ppi0_mass_map.begin()->second->GetXaxis()->SetTitle("p' #pi^{0} inv. mass (GeV)");
    bin_width = get_bin_width(h_ppi0_mass_map.begin()->second);
    h_ppi0_mass_map.begin()->second->GetYaxis()->SetTitle(TString::Format("Events / %.3f", bin_width));
    h_ppi0_mass_map.begin()->second->SetMinimum(0.0);
    h_ppi0_mass_map.begin()->second->Draw("HIST");
    it = h_ppi0_mass_map.begin();
    ++it;
    for (; it != h_ppi0_mass_map.end(); ++it)
    {
        it->second->Draw("HIST SAME");
        leg_ppi0->AddEntry(
            it->second,
            TString::Format("P_{z}^{CM}(#pi^{0}) > %1.1f GeV", it->first), "f");
    }
    leg_ppi0->Draw();
    c_pi0_cut->SaveAs(
        TString::Format(
            "pi0_mom_cut_ppi0_mass%s.pdf",
            mc ? "_mc" : "_data"));

    // plot helcostheta vs helphi
    TCanvas *c_costheta_phi = new TCanvas("ccostheta_phi", "ccostheta_phi", 800, 600);
    for(std::map<float, TH2F*>::iterator it = h_costheta_phi_map.begin(); it != h_costheta_phi_map.end(); ++it)
    {
        // set any negative bin contents to zero for better representation of holes
        for(int binx = 1; binx <= it->second->GetNbinsX(); ++binx)
        {
            for(int biny = 1; biny <= it->second->GetNbinsY(); ++biny)
            {
                int bin_content = it->second->GetBinContent(binx, biny);
                if (bin_content < 0)
                    it->second->SetBinContent(binx, biny, 0);
            }
        }        
        it->second->GetXaxis()->SetTitle("#phi (rad)");
        it->second->GetXaxis()->SetRangeUser(-3.14, 3.14);
        it->second->GetYaxis()->SetTitle("cos #theta");
        it->second->GetYaxis()->SetRangeUser(-1.0, 1.0);
        if (it->first < -0.9)
            it->second->SetTitle("No cut");
        else
            it->second->SetTitle(TString::Format("P_{z}^{CM}(#pi^{0}) > %1.1f GeV", it->first));
        
        it->second->Draw("colz");
        c_costheta_phi->SaveAs(
            TString::Format(
                "pi0_mom_cut_costheta_phi_%1.1f%s.pdf",
                it->first,
                mc ? "_mc" : "_data"));
        c_costheta_phi->Clear();
    }

    return;
}