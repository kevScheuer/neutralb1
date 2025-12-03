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
std::pair<TString, TString> load_final_cuts(
    int perm_id,
    std::vector<TString> perm_particles,
    bool nominal_cuts
);
std::vector<std::pair<TString, TString>> 
get_amptools_branch_mappings(std::vector<TString> perm_particles);
std::map<TString, TString> create_var_branches_map(
    int perm_id,
    std::vector<TString> perm_particles
);

// CONSTANTS
const double OMEGA_MASS = 0.7826; // PDG omega mass in GeV
const double SIGNAL_HALF_WIDTH = 0.028; // half width of selected signal region in GeV
const double SIDEBAND_GAP = SIGNAL_HALF_WIDTH*2; // how far from omega mass the sideband starts

// FSRoot particle ordering
//   0       1          2             3              4              5
// [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (omega)] [pi0 (bachelor)]

// The ordering for amptools will be different; we need to map the index of fs root
// particles to the corresponding amptools particle that follows the order
//   0       1            2               3            4             5
// [beam] [proton] [pi0 (bachelor)] [pi0 (omega)] [pi+ (omega)] [pi- (omega)]
const std::map<int, TString> FS_INDEX_TO_AMPTOOLS_MAP = {
    {0, "B"}, {1, "1"}, {5, "2"}, {4, "3"}, {2, "4"}, {3, "5"}};

// define some FSRoot index labels for readability
const int BEAM = 0;
const int PROTON = 1;
const int PI_PLUS = 2;
const int PI_MINUS = 3;
const int PI0_OMEGA = 4;
const int PI0_BACHELOR = 5;


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
        TString best_chi2_dir = TString::Format(
            "/lustre24/expphy/volatile/halld/home/kscheuer/"
            "FSRoot-skimmed-trees/best-chi2/");
        TString data_file = TString::Format(            
            "%stree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_data.root", 
            best_chi2_dir.Data(), period);
        TString mc_file = TString::Format(
            "%stree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_ver03.1_mc.root", 
            best_chi2_dir.Data(), period);
        TString phsp_file = TString::Format(
            "%stree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_ver03.root", 
            best_chi2_dir.Data(), period);

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
            std::vector<std::pair<TString, TString>> branches = 
                get_amptools_branch_mappings(perm_particles);

            std::pair<TString, TString> final_cuts = 
                load_final_cuts(perm_id, perm_particles, nominal_cuts);

            std::map<TString, TString> var_to_branch = 
                create_var_branches_map(perm_id, perm_particles);

            // add the variable branches to the AmpTools mapping branches
            for (std::map<TString, TString>::iterator var_it = var_to_branch.begin();
                 var_it != var_to_branch.end(); ++var_it)
            {
                branches.push_back(
                    std::make_pair(var_it->first, var_it->second)
                );
            }

            // print the branch mapping for debugging
            for (auto pair : branches)
                std::cout << "Branch: " << pair.first << " -> " << pair.second << "\n";

            // skim the trees for data, signal MC, and phasespace MC
            skim_trees(
                NT,
                CATEGORY, 
                data_file, 
                period, 
                "data", 
                perm_id, 
                final_cuts,
                branches
            );

            // TODO: same cuts for signal MC, but phasespace might want other function?

        } // end loop over permutations
    } // end loop over run periods

    // TODO: once all run periods done, then we can hadd the final trees together and 
    // rename them to the more usual AmpToolInputTree format

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
    

    std::vector<TString> p4_components = {"EnP", "PxP", "PyP", "PzP"};

    // setup branches for our friend tree in this loop
    std::vector<std::pair<TString, TString>> branches;

    // create the 4-momenta branches. We want to pair the amptools-ordered string
    // to the fsroot-ordered string for this permutation
    for (int i = 0; i < perm_particles.size(); i++)
    {
        TString fs_particle = perm_particles[i];
        TString amptools_particle = FS_INDEX_TO_AMPTOOLS_MAP.at(i);

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
 * @return std::pair<TString, TString> A pair of strings where the first element is
 * the comma separated normal cuts and the second element is the comma separated 
 * sideband cuts
 */
std::pair<TString, TString> load_final_cuts(
    int perm_id,
    std::vector<TString> perm_particles,
    bool nominal_cuts
)
{
    // we'll combine these vectors together at the end, but this makes tracking them 
    // easier
    std::vector<TString> single_vec; 
    std::vector<TString> sideband_vec;

    // ==== omega sideband subtraction cut definition ====
    TString data_omega_mass = TString::Format(
        "MASS(%s,%s,%s)",
        perm_particles[PI_PLUS].Data(),
        perm_particles[PI_MINUS].Data(),
        perm_particles[PI0_OMEGA].Data());
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
    sideband_vec.push_back(
        TString::Format("omega_sb_subtraction_perm%d", perm_id)
    );
    
    // ==== omega pi0 mass region ====
    TString omega_pi0_mass = TString::Format(
        "MASS(%s,%s,%s,%s)",
        perm_particles[PI_PLUS].Data(),
        perm_particles[PI_MINUS].Data(),
        perm_particles[PI0_OMEGA].Data(),
        perm_particles[PI0_BACHELOR].Data());
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
    single_vec.push_back(
        TString::Format("omega_pi0_mass_cut_perm%d", perm_id)
    );

    // ==== stable cuts ====
    // now define some cuts that will always remain consistent (no systematic variations)
    // regardless of nominal or systematic cut settings

    // require -t < 1.0 GeV^2
    FSCut::defineCut("t", "abs(-1*MASS2(1,-GLUEXTARGET))<1.0");
    // avoid interactions with walls of target
    FSCut::defineCut("z", "ProdVz>=51.2&&ProdVz<=78.8");
    // rf sideband subtraction
    FSCut::defineCut("rf", "abs(RFDeltaT)<0.2", "abs(RFDeltaT)>2.0", 0.125);
    single_vec.push_back("t");
    single_vec.push_back("z");
    sideband_vec.push_back("rf");

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
                "MOMENTUMZBOOST(%d;B,GLUEXTARGET)>-0.1", PI0_BACHELOR
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
                "MOMENTUMZBOOST(%d;B,GLUEXTARGET)>-0.3", PI0_BACHELOR   
            )); 
    }
    single_vec.push_back("unusedE");
    single_vec.push_back("unusedTracks");
    single_vec.push_back("MM2");
    single_vec.push_back("chi2");
    single_vec.push_back("shQuality");
    single_vec.push_back(TString::Format("pzPi0_perm%d", perm_id));

    // now combine our single and sideband cut vectors into comma-separated strings
    TString all_single_cuts;
    TString all_sideband_cuts;
    for (const TString& cut : single_vec)
    {
        if (all_single_cuts.IsNull())
        {
            all_single_cuts = cut;
        }
        else
        {
            all_single_cuts += "," + cut;
        }
    }
    for (const TString& cut : sideband_vec)
    {
        if (all_sideband_cuts.IsNull())
        {
            all_sideband_cuts = cut;
        }
        else
        {
            all_sideband_cuts += "," + cut;
        }
    }
    
    return std::make_pair(all_single_cuts, all_sideband_cuts);
}


std::map<TString, TString> create_var_branches_map(
    int perm_id,
    std::vector<TString> perm_particles
)
{
    // define some FSRoot index labels for readability
    int beam = 0;
    int proton = 1;
    int pi_plus = 2;
    int pi_minus = 3;
    int pi0_omega = 4;
    int pi0_bachelor = 5;

    std::map<TString, TString> var_to_branch;
    var_to_branch["rf"] = "RFDeltaT";
    var_to_branch["M4Pi"] = TString::Format(
        "MASS(%s,%s,%s,%s)",
        perm_particles[pi_plus].Data(),
        perm_particles[pi_minus].Data(),
        perm_particles[pi0_omega].Data(),
        perm_particles[pi0_bachelor].Data());
    var_to_branch["unusedE"] = "EnUnusedSh";
    var_to_branch["unusedTracks"] = "NumUnusedTracks";
    var_to_branch["z"] = "ProdVz";
    var_to_branch["MM2"] = "RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5)";
    var_to_branch["missingE"] = "RENERGY(GLUEXTARGET, B) - RENERGY(1, 2, 3, 4, 5)";
    var_to_branch["chi2"] = "Chi2DOF";
    var_to_branch["t"] = "-1*MASS2(1,-GLUEXTARGET)";
    var_to_branch["MRecoilPi"] = TString::Format(
        "MASS(%s,%s)",
        perm_particles[proton].Data(),
        perm_particles[pi0_bachelor].Data());
    var_to_branch["PzCMrecoilPi"] = TString::Format(
        "MOMENTUMZBOOST(%d;B,GLUEXTARGET)",
        pi0_bachelor
    );

    // shower quality variables need to be remapped to match AmpTools indexing
    var_to_branch[
        TString::Format(
            "shQualityP%sa", 
            FS_INDEX_TO_AMPTOOLS_MAP.at(PI0_OMEGA).Data()
        )] = TString::Format("ShQualityP%sa", perm_particles[PI0_OMEGA].Data());
    var_to_branch[
        TString::Format(
            "shQualityP%sb", 
            FS_INDEX_TO_AMPTOOLS_MAP.at(PI0_OMEGA).Data()
        )] = TString::Format("ShQualityP%sb", perm_particles[PI0_OMEGA].Data());
    var_to_branch[
        TString::Format(
            "shQualityP%sa", 
            FS_INDEX_TO_AMPTOOLS_MAP.at(PI0_BACHELOR).Data()
        )] = TString::Format("ShQualityP%sa", perm_particles[PI0_BACHELOR].Data());
    var_to_branch[
        TString::Format(
            "shQualityP%sb", 
            FS_INDEX_TO_AMPTOOLS_MAP.at(PI0_BACHELOR).Data()
        )] = TString::Format("ShQualityP%sb", perm_particles[PI0_BACHELOR].Data());

    return var_to_branch;
}


void skim_trees(
    TString NT, 
    TString CATEGORY,
    TString file, 
    int period,
    TString label,
    int perm_id, 
    std::pair<TString, TString> final_cuts,
    std::vector<std::pair<TString, TString>> branches
)
{
    // make copy since we have to add weight branch later
    std::vector<std::pair<TString, TString>> branches_copy = branches;

    std::cout << "Skimming file: " << file << "\n";
    std::cout << "Permutation " << perm_id << "\n";
    std::cout << "Applying cuts: " << final_cuts.first.Data() << "\n";
    std::cout << "Applying sideband cuts: " << final_cuts.second.Data() << "\n";

    TString out_dir = TString::Format(
    "/lustre24/expphy/volatile/halld/home/kscheuer/"
    "FSRoot-skimmed-trees/final-amptools-trees/");

    // define cuts for the 4 polarization orientations
    FSCut::defineCut("pol0", "PolarizationAngle>-5&&PolarizationAngle<5");
    FSCut::defineCut("pol45", "PolarizationAngle>40&&PolarizationAngle<50");
    FSCut::defineCut("pol90", "PolarizationAngle>85&&PolarizationAngle<95");
    FSCut::defineCut("pol135", "PolarizationAngle>130&&PolarizationAngle<140");

    std::vector<int> pol_angles = {0, 45, 90, 135};
    for (int angle : pol_angles)
    {
        std::cout << "Polarization angle: " << angle << "\n";
        TString pol_cut_name = TString::Format("pol%d", angle);    

        TString pol_angle_name = TString::Format("pol_%d", angle);
    
        // signal
        TString signal_out = TString::Format(
            "%stree_pi0pi0pippim__B4_finalAmptools_SKIM_0%d_%s_%s_perm%d_signal.root",
            out_dir.Data(),
            period,
            label.Data(),
            pol_angle_name.Data(),
            perm_id);
        // apply cuts, selected in signal region, for this pol orientation
        FSModeTree::skimTree(
            file,
            NT,
            CATEGORY,
            signal_out,
            TString::Format(
                "CUT(%s,%s,%s)", 
                final_cuts.first.Data(), final_cuts.second.Data(), pol_cut_name.Data())
        );
        // create friend with all of the branches
        FSTree::createFriendTree(
            signal_out,
            NT,                
            TString::Format("amptools_branches_perm%d", perm_id),
            branches_copy
        );

        // background
        TString data_out_background = TString::Format(
            "%stree_pi0pi0pippim__B4_finalAmptools_SKIM_0%d_%s_%s_perm%d_background.root",
            out_dir.Data(),
            period,
            label.Data(),
            pol_angle_name.Data(),
            perm_id);
         // apply cuts, with sideband weighting, for this pol orientation
        FSModeTree::skimTree(
            file,
            NT,
            CATEGORY,
            data_out_background,
            TString::Format(
                "CUTSBWT(%s)*CUT(%s)", 
                final_cuts.second.Data(), pol_cut_name.Data())
        );

        // create friend with all of the branches, including a "weight" one for background
        branches_copy.push_back(std::make_pair(
            "weight", 
            TString::Format(
                "CUTSBWT(%s)", final_cuts.second.Data())
        ));
        FSTree::createFriendTree(
            data_out_background,
            NT,                
            TString::Format("amptools_branches_perm%d", perm_id),
            branches_copy
        );
    } // end loop over polarization angles
}