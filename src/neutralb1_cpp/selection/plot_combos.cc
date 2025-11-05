/**
 * @file plot_combos.cc
 * @author Kevin Scheuer
 * @brief FSRoot script to plot combinatorics and signal vs background selections
 * 
 * Once the broad cuts are completed, the next step is to study the different
 * available combinations of showers and tracks to select the best combo
 * representing the physical final state. Remember that we have only selected
 * pi0pi0pi+pi- final states so far, but not yet selected the omega signal.
 * This script also visualizes the signal vs background regions in the selection of 
 * omega events.
 */


/* TODO: This script needs to include
    1. Visualize the hybrid method
    2. Plot combinatorics arising from the two pi0s
    3. Look at the sideband subtraction and its effect on omega pi0 mass (red negative box for bkg)
    4. Visualize pi0 momentum cut and its affect on p pi0 system
        4a. Dalitz plot for baryonic systems
        4b. The classic omega dalitz plot. Want to show that we need weighting for 
        the OmegaDalitz amplitude 
*/

#include <iostream>

#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"

#include "FSBasic/FSHistogram.h"
#include "FSBasic/FSCut.h"
#include "FSBasic/FSTree.h"
#include "FSMode/FSModeCollection.h"
#include "FSMode/FSModeHistogram.h"

#include "load_broad_cuts.cc"
#include "neutralb1/fit_utils.h"

TString NT("ntFSGlueX_MODECODE");
TString CATEGORY("pi0pi0pippim");

// Forward declarations
void plot_rf();
void plot_combinatorics();


void plot_combos(
    bool mc = true,
    bool bggen = false,
    bool read_cache = false,
    bool dump_cache = false)
{
    TString input_data_files;
    TString input_mc_files;

    input_data_files = "/cache/halld/home/jrsteven/flattened/omegapi_gx1_pwa/ver01/tree_pi0pi0pippim__B4/data/tree_pi0pi0pippim__B4_FSROOT_030568.root";
    input_mc_files = "/cache/halld/home/jrsteven/flattened/omegapi_gx1_pwa/ver01/tree_pi0pi0pippim__B4/omegapi_phasespace_2017_01_ver03/tree_pi0pi0pippim__B4_FSROOT_030568.root";

    if (FSModeCollection::modeVector().size() != 0)
        return;

    if (read_cache)
        FSHistogram::readHistogramCache();

    // https://github.com/JeffersonLab/hd_utilities/blob/master/FlattenForFSRoot/Documentation/GlueXFSRootFormat.pdf
    // unused numbers in front dropped, so this reads as 1 proton,
    // then 1 pi+, 1 pi-, 2 pi0s
    FSModeCollection::addModeInfo("100_112")->addCategory(CATEGORY);

    // load our broad cuts
    TString cuts = load_broad_cuts();

    // RF Subtraction
    plot_rf();

    // Plot Combinatorics Histograms
    plot_combinatorics();



}

void plot_rf(
    const TString input_data_files,
    const TString input_mc_files,
    const TString cuts,
    bool bggen = false    
)
{
    TCanvas *c_rf = new TCanvas("c_rf", "RF Delta T", 800, 600);

    TH1F *h_data = FSModeHistogram::getTH1F(
        input_data_files,
        NT,
        CATEGORY,
        "RFDeltaT",
        "(400, -20, 20)",
        cuts);
    h_data->SetXTitle("RF #Delta t (ns)");
    double bin_width = get_bin_width(h_data);
    h_data->SetYTitle(TString::Format("Combos / %.1f ns", bin_width));
    h_data->SetLineColor(kBlack);
    h_data->SetLineWidth(2);
    h_data->Draw();

    c_rf->Update();
    c_rf->SaveAs("RF_Delta_T.pdf");
}

// TODO: photon combinatorics too?
void plot_combinatorics()
{
    std::cout << "Plotting combinatorics histograms..." << std::endl;

    // First lets visualize the number of combinations per particle
    TCanvas *c_particle_combos = new TCanvas("c_particle_combos", "Particle Combinatorics", 800, 600);
    
}