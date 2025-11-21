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

#include "load_broad_cuts.cc"
#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"

// CONSTANTS
// TODO: these are FSROOT indices, but if we're reading from AmpTools reordered files
// this will be changed
const int BEAM = 0;
const int PROTON = 1;
const int PIPLUS = 2;
const int PIMINUS = 3;
const int PI0_OM = 4;
const int PI0_BACH = 5;

// proton, pi+, pi-
const TString M_OMEGA = TString::Format("MASS(%d,%d,%d)", PIPLUS, PIMINUS, PI0_OM);
const TString M_OMEGA_PI0 = TString::Format("MASS(%d,%d,%d,%d)", PIPLUS, PIMINUS, PI0_OM, PI0_BACH);

// forward declarations
void plot_correlations(
    TString NT, 
    TString CATEGORY,    
    TString input_files,
    bool mc);

void van_hove_analysis(int period, bool mc = false)
{
    gStyle->SetPalette( kBird );

    TString input_files;
    if (mc)
    {
        input_files = ""; // TODO: this should be data file with all cuts and background subtraction done
    }
    else
    {
        input_files = ""; // TODO: this should be mc file with all cuts and background subtraction done
    }
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(true);
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();
    TString cuts = join_keys(cut_color_map);

    plot_correlations(NT, CATEGORY, input_files, mc);

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
 * @brief Plot helicity angle vs omega pi0 invariant mass to look for correlations
 * 
 * @param NT FSRoot tree name
 * @param CATEGORY FSRoot mode category
 * @param input_files files to be analyzed
 * @param mc whether the input files are from Monte Carlo simulation
 */
void plot_correlations(
    TString NT, 
    TString CATEGORY,    
    TString input_files,
    bool mc)
{
    TCanvas *c_corr = new TCanvas("ccorr", "ccorr", 600, 800);

    // Under X-> omega pi0 p', we need to input into HELCOSTHETA(omega;pi0;p')
    TH2F *h_corr = FSModeHistogram::getTH2F(
        input_files,
        NT,
        CATEGORY,
        TString::Format(
            "HELCOSTHETA(%d,%d,%d;%d;%d):%s", 
            PIPLUS, PIMINUS, PI0_OM, PI0_BACH, PROTON, M_OMEGA_PI0.Data()),
        "(100, 1.0, 2.0, 100, -1.0, 1.0)",
        ""
    );
    h_corr->GetXaxis()->SetTitle("#omega #pi^{0} inv. mass (GeV)");
    h_corr->GetYaxis()->SetTitle("cos #theta_{H}");

    h_corr->Draw("colz");
    c_corr->SaveAs(
        TString::Format(
            "costhetah_mass_corr%s.pdf",
            mc ? "_mc" : "_data",
            CATEGORY.Data()
        )
    );

    FSHistogram::dumpHistogramCache();

    return;
}

