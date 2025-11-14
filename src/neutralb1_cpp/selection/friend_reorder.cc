/**
 * @file skim_reorder.cc
 * @author Kevin Scheuer
 * @brief Apply broad cuts and tag the omega and bachelor pi0's, create friend trees
 * 
 * @note Run this code once to create friend trees, then use the friend trees
 *       in subsequent analyses to have the correct pi0 assignments.
 */

#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "TString.h"

#include "FSBasic/FSCut.h"
#include "FSBasic/FSTree.h"

#include "load_broad_cuts.cc"
#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"


// CONSTANTS
TString NT("ntFSGlueX_MODECODE");
TString CATEGORY("pi0pi0pippim");

double OMEGA_MASS = 0.7826;             // PDG omega mass in GeV
double SIGNAL_WIDTH = 0.00868 * 3;      // exactly 3 sigma of PDG omega width
double SIDEBAND_GAP = SIGNAL_WIDTH * 3; // 9 sigma distance from signal region


void friend_reorder(int period, bool mc=false)
{
    TString input_files;
    if (mc)
    {
        input_files = TString::Format(
            "/lustre24/expphy/volatile/halld/home/kscheuer/"
            "FSRoot-skimmed-trees/best-chi2/"
            "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_ver03.1_mc.root",
            period
        );
    }
    else
    {
        input_files = TString::Format(
            "/lustre24/expphy/volatile/halld/home/kscheuer/"
            "FSRoot-skimmed-trees/best-chi2/"
            "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_data.root",
            period
        );
    }   
    setup(false);
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();


    // Our particle orders determine what decay they are from
    //   0       1          2             3              4              5
    // [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (omega)] [pi0 (bachelor)]
    //
    // We need to permute the pi0's to try both possibilities
    std::map<int, std::vector<TString>> fs_permutation_to_order_map = {
        {1, {"B", "1", "2", "3", "4", "5"}}, // particle 4 = pi0 from omega
        {2, {"B", "1", "2", "3", "5", "4"}}  // particle 5 = pi0 from omega
    };

    // The ordering for amptools will be different; we need to map the index of fs root
    // particles to the corresponding amptools particle that follows the order
    //   0       1            2               3            4             5
    // [beam] [proton] [pi0 (bachelor)] [pi0 (omega)] [pi+ (omega)] [pi- (omega)]
    std::map<int, TString> fs_index_to_amptools_map = {
        {0, "B"}, {1, "1"}, {5, "2"}, {4, "3"}, {2, "4"}, {3, "5"}
    };

    std::vector<TString> p4_components = {"EnP", "PxP", "PyP", "PzP"};

    // begin loop over permutations
    for (std::map<int, std::vector<TString>>::const_iterator perm_it = fs_permutation_to_order_map.begin();
         perm_it != fs_permutation_to_order_map.end(); ++perm_it)
    {
        int p = perm_it->first;
        std::vector<TString> perm_particles = perm_it->second;

        // setup branches for our friend tree in this loop
        std::vector<std::pair<TString, TString>> branches;

        // create the 4-momenta branches. We want to pair the amptools-ordered string
        // to the fsroot-ordered string for this permutation
        for(int i=0; i<perm_particles.size(); i++)
        {
            TString fs_particle = perm_particles[i];
            TString amptools_particle = fs_index_to_amptools_map[i];

            for(const TString comp : p4_components)
            {
                TString amptools_branch_name = comp + amptools_particle;
                TString fs_branch_name = comp + fs_particle;
                branches.push_back(std::make_pair(amptools_branch_name, fs_branch_name));
            }
        }

        // with the mapping done, lets define some index labels for readability
        int beam = 0;
        int proton = 1;
        int pi_plus = 2;
        int pi_minus = 3;
        int pi0_omega = 4;
        int pi0_bachelor = 5;

        // Now we'll want to define the signal and sideband cuts for the omega
        TString data_omega_mass = TString::Format(
            "MASS(%s,%s,%s)",
            perm_particles[pi_plus].Data(),
            perm_particles[pi_minus].Data(),
            perm_particles[pi0_omega].Data()
        );
        TString signal_region_cut = TString::Format(
            "abs(%s-%f)<%f",
            data_omega_mass.Data(),
            OMEGA_MASS,
            SIGNAL_WIDTH
        );
        TString sideband_region_cut = TString::Format(
            "abs(%s-%f)>(%f)&&abs(%s-%f)<(%f+%f)",
            data_omega_mass.Data(),
            OMEGA_MASS,
            SIDEBAND_GAP,
            data_omega_mass.Data(),
            OMEGA_MASS,
            SIDEBAND_GAP,
            SIGNAL_WIDTH
        );
        FSCut::defineCut(
            "omega_sb_subtraction",
            signal_region_cut,
            sideband_region_cut,
            1.0
        );

        // define a cut for the omega pi mass region of interest
        TString omega_pi0_mass = TString::Format(
            "MASS(%s,%s,%s,%s)",
            perm_particles[pi_plus].Data(),
            perm_particles[pi_minus].Data(),
            perm_particles[pi0_omega].Data(),
            perm_particles[pi0_bachelor].Data()
        );
        TString omega_pi0_mass_cut = TString::Format(
            "%s>%f&&%s<%f",
            omega_pi0_mass.Data(),
            1.0, // min omega pi0 mass
            omega_pi0_mass.Data(),
            2.0 // max omega pi0 mass
        );
        FSCut::defineCut("omega_pi0_mass_cut", omega_pi0_mass_cut);

        // Now we can create the branches for the signal and sideband regions. We'll
        // also make a "cut" branch where RF+broad cuts are applied, along with some 
        // other useful distributions we'll want to bin in or study for our PWA
        branches.push_back(std::make_pair(
            "cut", TString::Format("CUT(%s,omega_pi0_mass_cut)*CUTWT(rf)", join_keys(cut_color_map).Data())
        ));
        branches.push_back(std::make_pair("signal", "CUT(omega_sb_subtraction)"));
        branches.push_back(std::make_pair("sideband", "CUTSB(omega_sb_subtraction)"));
        branches.push_back(std::make_pair("t", "-1*MASS2(1,-GLUEXTARGET)"));
        branches.push_back(std::make_pair("M4Pi", omega_pi0_mass));
        branches.push_back(std::make_pair("MRecoilPi", TString::Format(
            "MASS(%s,%s)",            
            perm_particles[proton].Data(),
            perm_particles[pi0_bachelor].Data()
        )));

        // finally create the friend tree for this permutation with these branches
        TString friend_tree_name = TString::Format("permutation_%d", p);
        FSTree::createFriendTree(input_files, NT, friend_tree_name, branches);
    }
    return;
}