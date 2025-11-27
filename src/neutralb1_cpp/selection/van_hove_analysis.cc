/**
 * @file van_hove_analysis.cc
 * @author Kevin Scheuer, Amy Schertz
 * @brief Plot van Hove distributions for data and MC
 *
 * Once again, Amy Schertz deserves credit for the original authorship of this code.
 */

#include <iostream>
#include <tuple>

#include "TArrow.h"
#include "TCanvas.h"
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

    h_corr[0]->GetXaxis()->SetTitle("#omega #pi^{0} inv. mass (GeV)");
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
    h_dalitz[0]->GetXaxis()->SetTitle("M^{2}(#omega #pi^{0}) (GeV^{2})");
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
    h_pi_corr_omegapi0[0]->GetXaxis()->SetTitle("#omega #pi^{0} inv. mass (GeV)");
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
    // TODO: plot the different pi0 momenta cuts on:
    //      cos theta
    //      omega pi0 mass
    //      p' pi0 mass
    // also plot the helcostheta vs helphi 2D distribution

    return;
}