/**
 * @file relabel_gen_phasespace.cc
 * @author Kevin Scheuer
 * @brief Relabel gen phasespace branches for AmpTools
 *
 * After the gen phasespace trees have been skimmed, we need to relabel the branches 
 * to match the AmpTools expectations, similar to what is done in 
 * finalize_amptools_trees.cc for the signal and background trees. We don't need to 
 * apply any permutation handling here, since the trees are generated knowing the
 * correct particle assignments.
 */

#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "FSBasic/FSTree.h"
#include "FSMode/FSModeTree.h"

#include "fsroot_setup.cc"

// forward declarations
std::vector<std::pair<TString, TString>> 
get_amptools_branch_mappings_gen(std::vector<TString> perm_particles);
std::map<TString, TString> create_var_branches_map_gen(
    std::vector<TString> perm_particles
);
void process_gen_phasespace_file(
    TString NT,
    TString CATEGORY,
    TString input_file,
    int period,
    std::vector<std::pair<TString, TString>> branches
);
void create_common_key(TString friend_file, TString NT, TString friend_extension);

// CONSTANTS
const TString TREE_OUTPUT_DIR = "/lustre24/expphy/volatile/halld/home/kscheuer/"
    "FSRoot-skimmed-trees/final-amptools-trees/";

// FSRoot particle ordering for gen phasespace (bachelor pi0 is P4, omega pi0 is P5)
//   0       1          2             3              4              5
// [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (bachelor)] [pi0 (omega)]

// The ordering for amptools will be different; we need to map the index of fs root
// particles to the corresponding amptools particle that follows the order
//   0       1            2               3            4             5
// [beam] [proton] [pi0 (bachelor)] [pi0 (omega)] [pi+ (omega)] [pi- (omega)]
const std::map<int, TString> FS_INDEX_TO_AMPTOOLS_MAP = {
    {0, "B"}, {1, "1"}, {4, "2"}, {5, "3"}, {2, "4"}, {3, "5"}};

// define some FSRoot index labels for readability
const int BEAM = 0;
const int PROTON = 1;
const int PI_PLUS = 2;
const int PI_MINUS = 3;
const int PI0_BACHELOR = 4;
const int PI0_OMEGA = 5;


void relabel_gen_phasespace()
{
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(false);    

    TString gen_phasespace_mc_files_3 = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/general/"
        "tree_pi0pi0pippim__B4_GENERAL_SKIM_03_ver03_gen.root";
    TString gen_phasespace_mc_files_4 = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/general/"
        "tree_pi0pi0pippim__B4_GENERAL_SKIM_04_ver03_gen.root";
    TString gen_phasespace_mc_files_5 = "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/general/"
        "tree_pi0pi0pippim__B4_GENERAL_SKIM_05_ver03_gen.root";
    
    // Particle ordering (no permutations needed for gen phasespace)
    //   0       1          2             3              4              5
    // [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (bachelor)] [pi0 (omega)]
    std::vector<TString> perm_particles = {"B", "1", "2", "3", "4", "5"};

    // Get branch mappings for AmpTools
    std::vector<std::pair<TString, TString>> branches = 
        get_amptools_branch_mappings_gen(perm_particles);
    
    // Add variable branches (M4Pi, MRecoilPi, t)
    std::map<TString, TString> var_branches = create_var_branches_map_gen(perm_particles);
    for (const auto& pair : var_branches)
    {
        branches.push_back(pair);
    }

    // Process each run period
    std::vector<int> run_periods = {3, 4, 5};
    std::vector<TString> files = {
        gen_phasespace_mc_files_3,
        gen_phasespace_mc_files_4,
        gen_phasespace_mc_files_5
    };

    for (size_t i = 0; i < run_periods.size(); i++)
    {
        process_gen_phasespace_file(
            NT,
            CATEGORY,
            files[i],
            run_periods[i],
            branches
        );
    }

    // Combine all run periods together
    std::cout << "Combining all run periods..." << std::endl;
    system(TString::Format(
        "cd %s && hadd -f allPeriods_ver03_gen_phasespace.root "
        "2017_01_ver03_gen_phasespace.root "
        "2018_01_ver03_gen_phasespace.root "
        "2018_08_ver03_gen_phasespace.root",
        TREE_OUTPUT_DIR.Data()
    ).Data());
}


/**
 * @brief Get the branch mappings between FSRoot and AmpTools for gen phasespace
 * 
 * @param perm_particles The permutation of particles in FSRoot order
 * @return std::vector<std::pair<TString, TString>> A vector of pairs where the first 
 * element is the AmpTools branch name and the second element is the FSRoot branch name
 */
std::vector<std::pair<TString, TString>> 
get_amptools_branch_mappings_gen(std::vector<TString> perm_particles)
{
    std::vector<TString> p4_components = {"EnP", "PxP", "PyP", "PzP"};

    // setup branches for our friend tree
    std::vector<std::pair<TString, TString>> branches;

    // create the 4-momenta branches. We want to pair the amptools-ordered string
    // to the fsroot-ordered string for gen phasespace (with MC prefix)
    for (int i = 0; i < perm_particles.size(); i++)
    {
        TString fs_particle = perm_particles[i];
        TString amptools_particle = FS_INDEX_TO_AMPTOOLS_MAP.at(i);

        for (const TString comp : p4_components)
        {
            TString amptools_branch_name = comp + amptools_particle;
            // For gen phasespace, use MC prefix
            TString fs_branch_name = "MC" + comp + fs_particle;
            branches.push_back(std::make_pair(amptools_branch_name, fs_branch_name));
        }
    }

    return branches;
}


/**
 * @brief Sets up a map for the variable branches (M4Pi, MRecoilPi, t)
 * 
 * @param perm_particles The permutation of particles in FSRoot order
 * @return std::map<TString, TString> A map where the key is the variable name
 * and the value is the corresponding branch name/formula in the FSRoot tree
 */
std::map<TString, TString> create_var_branches_map_gen(
    std::vector<TString> perm_particles
)
{
    std::map<TString, TString> var_to_branch;
    
    // M4Pi: mass of all four pions
    var_to_branch["M4Pi"] = TString::Format(
        "MCMASS(%s,%s,%s,%s)",
        perm_particles[PI_PLUS].Data(),
        perm_particles[PI_MINUS].Data(),
        perm_particles[PI0_OMEGA].Data(),
        perm_particles[PI0_BACHELOR].Data());
    
    // MRecoilPi: mass of proton and bachelor pi0
    var_to_branch["MRecoilPi"] = TString::Format(
        "MCMASS(%s,%s)",
        perm_particles[PROTON].Data(),
        perm_particles[PI0_BACHELOR].Data());
    
    // t: four-momentum transfer squared
    var_to_branch["t"] = "abs(-1*MCMASS2(1,-GLUEXTARGET))";

    return var_to_branch;
}


/**
 * @brief Process a gen phasespace file for a given run period
 * 
 * @param NT tree name
 * @param CATEGORY tree category
 * @param input_file input file path
 * @param period run period
 * @param branches vector of branch mappings to write to friend tree
 */
void process_gen_phasespace_file(
    TString NT,
    TString CATEGORY,
    TString input_file,
    int period,
    std::vector<std::pair<TString, TString>> branches
)
{
    std::cout << "Processing file: " << input_file << std::endl;
    
    std::map<int, TString> period_to_name_map = {
        {3, "2017_01"},
        {4, "2018_01"},
        {5, "2018_08"}
    };

    // Output file name
    TString output_file = TString::Format(
        "%s%s_ver03_gen_phasespace.root",
        TREE_OUTPUT_DIR.Data(),
        period_to_name_map.at(period).Data()
    );

    // Create friend tree with relabeled branches
    std::cout << "Creating friend tree with relabeled branches..." << std::endl;
    FSTree::createFriendTree(
        input_file,
        NT,
        "amptools_branches_gen",
        branches
    );

    // Create common key for the friend tree
    create_common_key(
        input_file,
        NT,
        "amptools_branches_gen"
    );

    // copy the friend tree to be the output file
    std::cout << "Copying friend tree to output file: " << output_file << std::endl;
    system(TString::Format(
        "cp -f %s.amptools_branches_gen %s",
        input_file.Data(),
        output_file.Data()
    ).Data());
}


/**
 * @brief Add a common key to the friend trees
 * 
 * Renames the friend tree to match the original tree name so they can be 
 * properly combined with hadd.
 * 
 * @param original_file ROOT file the friend tree is built off of
 * @param NT Original tree name
 * @param friend_extension friend tree name that extends on NT
 */
void create_common_key(TString original_file, TString NT, TString friend_extension)
{
    TString friend_file = original_file + "." + friend_extension;
    TString tree_name = NT + "_" + friend_extension;

    TFile *f = TFile::Open(friend_file, "UPDATE");    
    TTree *t = (TTree*)f->Get(tree_name);
    t->SetName(NT);
    t->Write("",TObject::kOverwrite);
    f->Close();
}