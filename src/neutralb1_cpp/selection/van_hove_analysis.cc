/**
 * @file van_hove_analysis.cc
 * @author Kevin Scheuer, Amy Schertz
 * @brief Plot van Hove distributions for data and MC 
 * 
 * Once again, Amy Schertz deserves credit for the original authorship of this code. 
 */

#include <iostream>
#include <tuple>

#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"

#include "FSBasic/FSHistogram.h"
#include "FSBasic/FSCut.h"
#include "FSMode/FSModeHistogram.h"

#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"

// CONSTANTS
// AmpTools indices for particles in the final state
const int BEAM = 0; // four vectors named like energy = "EnPB"
const int PROTON = 1;
const int PI0_BACH = 2;
const int PI0_OM = 3;
const int PIPLUS = 4;
const int PIMINUS = 5;

// proton, pi+, pi-
const TString M_OMEGA = TString::Format("MASS(%d,%d,%d)", PIPLUS, PIMINUS, PI0_OM);
const TString M_OMEGA_PI0 = TString::Format("MASS(%d,%d,%d,%d)", PIPLUS, PIMINUS, PI0_OM, PI0_BACH);

// forward declarations
void plot_relations(
    TString NT, 
    TString CATEGORY,    
    TString input_signal,
    TString input_sideband,
    bool mc);

void van_hove_analysis(bool mc = false)
{
    gStyle->SetPalette( kBird );

    TString input_signal;
    TString input_sideband;
    if (mc)
    {
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

    plot_relations(NT, CATEGORY, input_signal, input_sideband, mc);    

    /* TODO: want to follow these steps
        1. DONE. Plot the cos(theta) vs baryon+bachelor mass for data and MC 
        2. Plot the van Hove distributions for data and MC
        3. See where Delta contribution peaks in van Hove 
        4. Apply cut on resulting bachelor pion momentum
        5. Observe cut's effect on phasespace angles and the helicity angles and
            baryon+bachelor system
    */
    return;
}

/**
 * @brief Plot helicity angle vs omega pi0 and a dalitz plot to look for baryon resonances
 * 
 * @param NT FSRoot tree name
 * @param CATEGORY FSRoot mode category
 * @param input_signal signal input file
 * @param input_sideband sideband input file
 * @param mc whether the input files are from Monte Carlo simulation
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
        "(100, 1.0, 2.0, 100, -1.0, 1.0)",
        "cut==1"
    );
    h_corr[1] = FSModeHistogram::getTH2F(
        input_sideband,
        NT,
        CATEGORY,
        TString::Format(
            "HELCOSTHETA(%d,%d,%d;%d;%d):%s", 
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON, M_OMEGA_PI0.Data()),
        "(100, 1.0, 2.0, 100, -1.0, 1.0)",
        "cut==1"
    );
    h_corr[0]->Add(h_corr[1], -1); // signal - sideband
    

    h_corr[0]->GetXaxis()->SetTitle("#omega #pi^{0} inv. mass (GeV)");
    h_corr[0]->GetYaxis()->SetTitle("cos #theta_{H}");
    h_corr[0]->SetTitle("");
    h_corr[0]->Draw("colz");
    c_corr->SaveAs(
        TString::Format(
            "costhetah_mass_corr%s.pdf",
            mc ? "_mc" : "_data"
        )
    );

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
        "(100, 1.0, 3.0, 100, 1.0, 4.0)",
        "cut==1"
    );
    h_dalitz[1] = FSModeHistogram::getTH2F(
        input_sideband,
        NT,
        CATEGORY,
        TString::Format(
            "MASS2(%d,%d):MASS2(%d,%d,%d, %d)",             
            PI0_BACH, PROTON,
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH),
        "(100, 0.0, 2.0, 100, 0.0, 3.0)",
        "cut==1"
    );
    h_dalitz[0]->Add(h_dalitz[1], -1); // signal - sideband
    h_dalitz[0]->GetXaxis()->SetTitle("M^{2}(#omega #pi^{0}) (GeV^{2})");
    h_dalitz[0]->GetYaxis()->SetTitle("M^{2}(p' #pi^{0}) (GeV^{2})");
    h_dalitz[0]->SetTitle("");
    h_dalitz[0]->Draw("colz");
    c_dalitz->SaveAs(
        TString::Format(
            "dalitz_omegapi0_ppi0%s.pdf",
            mc ? "_mc" : "_data"          
        )
    );

    // FSHistogram::dumpHistogramCache();

    return;
}

void plot_van_hove(
    TString NT, 
    TString CATEGORY,    
    TString input_signal,
    TString input_sideband,
    bool mc)
{

    // Three particle final state (omega, pi_bach, and recoil proton) is fed into 
    // van hove function
    TH2F* h_van_hove[2];
    h_van_hove[0] = FSModeHistogram::getTH2F(
        input_signal,
        NT,
        CATEGORY,
        TString::Format(
            "VANHOVEY(%d,%d,%d;%d;%d):VANHOVEX(%d,%d,%d;%d;%d)",
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON,
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON
        ),
        "(100, -4.0, 4.0, 100, -4.0, 4.0)",
        ""
    );
    h_van_hove[1] = FSModeHistogram::getTH2F(
        input_sideband,
        NT,
        CATEGORY,
        TString::Format(
            "VANHOVEY(%d,%d,%d;%d;%d):VANHOVEX(%d,%d,%d;%d;%d)",
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON,
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON
        ),
        "(100, -4.0, 4.0, 100, -4.0, 4.0)",
        ""
    );
    h_van_hove[0]->Add(h_van_hove[1], -1); // signal - sideband
    h_van_hove[0]->GetXaxis()->SetTitle("X");
    h_van_hove[0]->GetYaxis()->SetTitle("Y");
    TCanvas *c_van_hove = new TCanvas("cvanhove", "cvanhove", 800, 600);
    h_van_hove[0]->Draw("colz");
    c_van_hove->SaveAs(
        TString::Format(
            "van_hove%s.pdf",
            mc ? "_mc" : "_data"
        )
    );
    
}