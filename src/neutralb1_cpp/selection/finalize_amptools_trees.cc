/**
 * @file finalize_amptools_trees.cc
 * @author Kevin Scheuer
 * @brief Produce the AmpTools signal and background trees for all files
 * 
 * This script applies all cuts and sideband subtractions to create the final 
 * flattened AmpTools trees for fitting. Cuts are appropriately applied to data, 
 * signal MC, and phasespace MC files. A lot of this is repetitive from 
 * friend_reorder.cc, but its better to start from the chi2 trees to start
 * "from scratch" and ensure all cuts are properly applied.
 * 
 */

#include <iostream>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "FSBasic/FSCut.h"
#include "FSBasic/FSTree.h"
#include "FSMode/FSModeTree.h"

#include "load_broad_cuts.cc"
#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"

// forward declarations
void load_final_cuts(
    int perm_id,
    std::vector<TString> perm_particles,
    bool nominal_cuts
);
std::vector<std::pair<TString, TString>> 
get_amptools_branch_mappings(std::vector<TString> perm_particles);

// CONSTANTS
const double OMEGA_MASS = 0.7826; // PDG omega mass in GeV
const double SIGNAL_HALF_WIDTH = 0.028; // half width of selected signal region in GeV
const double SIDEBAND_GAP = SIGNAL_HALF_WIDTH*2; // how far from omega mass the sideband starts

void finalize_amptools_trees(
    bool signal_mc = true, 
    bool phasespace_mc = true, 
    bool nominal_cuts = true)
{
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(false);    
    
    std::vector<int> run_periods = {3, 4, 5};
    for (int period : run_periods)
    {
        TString data_file = TString::Format(
            "/lustre24/expphy/volatile/halld/home/kscheuer/"
            "FSRoot-skimmed-trees/best-chi2/"
            "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_data.root", period);
        TString mc_file = TString::Format(
            "/lustre24/expphy/volatile/halld/home/kscheuer/"
            "FSRoot-skimmed-trees/best-chi2/"
            "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_ver03.1_mc.root", period);
        TString phsp_file = TString::Format(
            "/lustre24/expphy/volatile/halld/home/kscheuer/"
            "FSRoot-skimmed-trees/best-chi2/"
            "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_ver03.root", period);

        // Our particle orders determine what decay they are from
        //   0       1          2             3              4              5
        // [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (omega)] [pi0 (bachelor)]
        //
        // We need to permute the pi0's to try both possibilities
        std::map<int, std::vector<TString>> fs_permutation_to_order_map = {
            {1, {"B", "1", "2", "3", "4", "5"}}, // particle 4 = pi0 from omega
            {2, {"B", "1", "2", "3", "5", "4"}}  // particle 5 = pi0 from omega
        };

        // loop over both pi0 permutations
        for (std::map<int, std::vector<TString>>::iterator perm_it = fs_permutation_to_order_map.begin();
             perm_it != fs_permutation_to_order_map.end(); ++perm_it)
        {
            int perm_id = perm_it->first;
            std::vector<TString> perm_particles = perm_it->second;

            // get the 4-momenta branch mappings to AmpTools for this permutation
            std::vector<std::pair<TString, TString>> amptools_branch_mappings = 
                get_amptools_branch_mappings(perm_particles);

            load_final_cuts(perm_id, perm_particles, nominal_cuts);


        } // end loop over permutations

    }

    // TODO: apply cuts, set branches for everything including pz momenta
    // separate by signal and background for data and MC

    // reference Amy / Malte for signal and background files and what cuts to make 
    // on phasespace
    

    // TODO: extract polarization orientation
}


/**
 * @brief Get the branch mappings between FSRoot and AmpTools for a given permutation of particles
 * 
 * @param perm_particles The permutation of particles in FSRoot order. This will 
 * typically be either {B, 1, 2, 3, 4, 5} or {B, 1, 2, 3, 5, 4} depending on
 * which pi0 is assigned to the omega decay.
 * @return std::vector<std::pair<TString, TString>> A vector of pairs where the first 
 * element is the AmpTools branch name and the second element is the FSRoot branch name
 */
std::vector<std::pair<TString, TString>> 
get_amptools_branch_mappings(std::vector<TString> perm_particles)
{
    // FSRoot particle ordering
    //   0       1          2             3              4              5
    // [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (omega)] [pi0 (bachelor)]

    // The ordering for amptools will be different; we need to map the index of fs root
    // particles to the corresponding amptools particle that follows the order
    //   0       1            2               3            4             5
    // [beam] [proton] [pi0 (bachelor)] [pi0 (omega)] [pi+ (omega)] [pi- (omega)]
    std::map<int, TString> fs_index_to_amptools_map = {
        {0, "B"}, {1, "1"}, {5, "2"}, {4, "3"}, {2, "4"}, {3, "5"}};

    std::vector<TString> p4_components = {"EnP", "PxP", "PyP", "PzP"};

    // setup branches for our friend tree in this loop
    std::vector<std::pair<TString, TString>> branches;

    // create the 4-momenta branches. We want to pair the amptools-ordered string
    // to the fsroot-ordered string for this permutation
    for (int i = 0; i < perm_particles.size(); i++)
    {
        TString fs_particle = perm_particles[i];
        TString amptools_particle = fs_index_to_amptools_map[i];

        for (const TString comp : p4_components)
        {
            TString amptools_branch_name = comp + amptools_particle;
            TString fs_branch_name = comp + fs_particle;
            branches.push_back(std::make_pair(amptools_branch_name, fs_branch_name));
        }
    }

    return branches;
}


/**
 * @brief Define and load all the final cuts to be applied to the AmpTools trees
 * 
 * Cuts are divided into "stable", "nominal" and "loose" cuts depending on whether they
 * will be subject to systematics in the future. Stable cuts are not varied, nominal
 * cuts use the values for the main result, and loose cuts are looser to allow for
 * variations. The sideband regions are a potential systematic, but are not varied 
 * here since the subtraction is baked into the whole script, not just a simple branch.
 * 
 * @param perm_id permutation id (1 or 2)
 * @param perm_particles The permutation of particles in FSRoot order. This will 
 * typically be either {B, 1, 2, 3, 4, 5} or {B, 1, 2, 3, 5, 4} depending on
 * which pi0 is assigned to the omega decay.
 * @param nominal_cuts If true, use nominal cuts; if false, use loose cuts
 */
void load_final_cuts(
    int perm_id,
    std::vector<TString> perm_particles,
    bool nominal_cuts
)
{
    // define some FSRoot index labels for readability
    int beam = 0;
    int proton = 1;
    int pi_plus = 2;
    int pi_minus = 3;
    int pi0_omega = 4;
    int pi0_bachelor = 5;

    // ==== omega sideband subtraction cut definition ====
    TString data_omega_mass = TString::Format(
        "MASS(%s,%s,%s)",
        perm_particles[pi_plus].Data(),
        perm_particles[pi_minus].Data(),
        perm_particles[pi0_omega].Data());
    TString signal_region_cut = TString::Format(
        "abs(%s-%f)<%f",
        data_omega_mass.Data(),
        OMEGA_MASS,
        SIGNAL_HALF_WIDTH);
    TString sideband_region_cut = TString::Format(
        "abs(%s-%f)>(%f)&&abs(%s-%f)<(%f+%f)",
        data_omega_mass.Data(),
        OMEGA_MASS,
        2*SIDEBAND_GAP,
        data_omega_mass.Data(),
        OMEGA_MASS,
        2*SIDEBAND_GAP,
        SIGNAL_HALF_WIDTH);
    
    FSCut::defineCut(
        TString::Format("omega_sb_subtraction_perm%d", perm_id),
        signal_region_cut,
        sideband_region_cut,
        1.0);
    
    // ==== omega pi0 mass region ====
    TString omega_pi0_mass = TString::Format(
        "MASS(%s,%s,%s,%s)",
        perm_particles[pi_plus].Data(),
        perm_particles[pi_minus].Data(),
        perm_particles[pi0_omega].Data(),
        perm_particles[pi0_bachelor].Data());
    TString omega_pi0_mass_cut = TString::Format(
        "%s>%f&&%s<%f",
        omega_pi0_mass.Data(),
        1.0, // min omega pi0 mass
        omega_pi0_mass.Data(),
        2.0 // max omega pi0 mass
    ); 
    FSCut::defineCut(
        TString::Format("omega_pi0_mass_cut_perm%d", perm_id),
        omega_pi0_mass_cut);

    // ==== stable cuts ====
    // now define some cuts that will always remain consistent (no systematic variations)
    // regardless of nominal or systematic cut settings

    // require -t < 1.0 GeV^2
    FSCut::defineCut("t", "abs(-1*MASS2(1,-GLUEXTARGET))<1.0");
    // avoid interactions with walls of target
    FSCut::defineCut("z", "ProdVz>=51.2&&ProdVz<=78.8");
    // rf sideband subtraction
    FSCut::defineCut("rf", "abs(RFDeltaT)<0.2", "abs(RFDeltaT)>2.0", 0.125);

    // ==== systematic cuts ====
    // these cuts may be varied later for systematic checks. We'll separate them into
    // "nominal" and "loose" cuts. Nominal cuts are the ones we'll use for the main
    // result, while loose cuts are, well, looser versions of those cuts so that we can
    // perform systematic variations on them.
    if (nominal_cuts)
    {     
        FSCut::defineCut("unusedE", "EnUnusedSh<0.1"); // energy of unmatched showers in time with combo but not part of the combo
        FSCut::defineCut("unusedTracks", "NumUnusedTracks<1"); // number charged particle tracks not used in combo
        FSCut::defineCut("MM2", "abs(RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5))<0.05"); // remove events with large missing mass2
        FSCut::defineCut("chi2", "Chi2DOF<5");  // kinematic fit quality cut
        FSCut::defineCut( // classifier probability that showers are from neutral particles
            "shQuality", 
            "ShQualityP4a>0.5&&ShQualityP4b>0.5&&ShQualityP5a>0.5&&ShQualityP5b>0.5" 
        );
        FSCut::defineCut( // require bachelor pi0 moving forward in CM frame
            TString::Format("pzPi0_perm%d", perm_id),
            TString::Format(
                "MOMENTUMZBOOST(%d;B,GLUEXTARGET)>-0.1", pi0_bachelor
            )); 
    }
    else 
    {
        FSCut::defineCut("unusedE", "EnUnusedSh<0.5"); // energy of unmatched showers in time with combo but not part of the combo
        FSCut::defineCut("unusedTracks", "NumUnusedTracks<4"); // number charged particle tracks not used in combo
        FSCut::defineCut("MM2", "abs(RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5))<0.08"); // remove events with large missing mass2
        FSCut::defineCut("chi2", "Chi2DOF<10");  // kinematic fit quality cut
        FSCut::defineCut( // classifier probability that showers are from neutral particles
            "shQuality", 
            "ShQualityP4a>0.5&&ShQualityP4b>0.5&&ShQualityP5a>0.5&&ShQualityP5b>0.3" 
        );
        FSCut::defineCut( // require bachelor pi0 moving forward in CM frame
            TString::Format("pzPi0_perm%d", perm_id),
            TString::Format(
                "MOMENTUMZBOOST(%d;B,GLUEXTARGET)>-0.3", pi0_bachelor
            )); 
    }
    
    return;
}