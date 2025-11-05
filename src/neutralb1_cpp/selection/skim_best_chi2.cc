/**
 * @file skim_flattened_trees.cc
 * @author Kevin Scheuer
 * @brief Skim down the trees to reduce their size in preparation for broad cuts
 *
 * Once trees have been flattened for FSRoot analysis, they can still be quite large.
 * This script selects the best chi2 combination per event and beam photon to
 * prepare for the hybrid method, and skim out unnecessary events. It will create new
 * files that have a very loose chi2 cut imposed and then create chi2 ranked trees
 * used to select the best combination, and these will be saved to new files.
 */

#include "FSBasic/FSHistogram.h"
#include "FSBasic/FSCut.h"
#include "FSBasic/FSTree.h"
#include "FSMode/FSModeTree.h"
#include "FSMode/FSModeCollection.h"
#include "FSMode/FSModeHistogram.h"

#include "fsroot_setup.cc"

TString NT("ntFSGlueX_MODECODE");
TString CATEGORY("pi0pi0pippim");

void skim_best_chi2()
{
    TString input_files = "/cache/halld/home/jrsteven/flattened/omegapi_gx1_pwa/ver01/tree_pi0pi0pippim__B4/data/tree_pi0pi0pippim*.root";
    // signal includes thrown and reconstructed
    TString signal_mc_files = "/cache/halld/home/jrsteven/flattened/omegapi_gx1_pwa/ver01/tree_pi0pi0pippim__B4/omegapi_massDepFit_201*ver03.1/tree_pi0pi0pippim*.root";
    TString phasespace_mc_files = "/cache/halld/home/jrsteven/flattened/omegapi_gx1_pwa/ver01/tree_pi0pi0pippim__B4/omegapi_phasespace*ver03/tree_pi0pi0pippim*.root";

    setup(false);

    // Following this document https://halldweb.jlab.org/DocDB/0065/006572/003/RankingMethods_in_FSRoot.pdf
    // only a few cuts commute with the chi2 ranking method. We will thus only implement
    // a loose chi2 cut here to produce a skimmed tree that can then be ranked properly.

    FSCut::defineCut("chi2", "Chi2DOF<20");
    TString output_data = "/w/halld-scshelf2101/kscheuer/neutralb1/data/FSRoot_skimmed_trees/general/tree_pi0pi0pippim__B4_GENERAL_SKIM_allPeriods_data.root";
    TString output_signal_mc = "/w/halld-scshelf2101/kscheuer/neutralb1/data/FSRoot_skimmed_trees/general/tree_pi0pi0pippim__B4_GENERAL_SKIM_allPeriods_ver03.1_mc.root";
    TString output_phasespace_mc = "/w/halld-scshelf2101/kscheuer/neutralb1/data/FSRoot_skimmed_trees/general/tree_pi0pi0pippim__B4_GENERAL_SKIM_allPeriods_ver03.root";

    // Now skim the trees and save to new files
    FSModeTree::skimTree(
        input_files,
        NT,
        CATEGORY,
        output_data,
        "CUT(chi2)");
    FSModeTree::skimTree(
        signal_mc_files,
        NT,
        CATEGORY,
        output_signal_mc,
        "CUT(chi2)");
    FSModeTree::skimTree(
        phasespace_mc_files,
        NT,
        CATEGORY,
        output_phasespace_mc,
        "CUT(chi2)");

    // Following the doc from above, we'll group by run, event, and beam energy to
    // create chi2 ranked trees
    FSModeTree::createRankingTree(
        output_data,
        NT,
        CATEGORY,
        "Chi2RankHybrid",
        "Chi2*1000",
        "1==1",
        "Run",
        "Event",
        "EnPB*1000");
    FSModeTree::createRankingTree(
        output_signal_mc,
        NT,
        CATEGORY,
        "Chi2RankHybrid",
        "Chi2*1000",
        "1==1",
        "Run",
        "Event",
        "EnPB*1000");
    FSModeTree::createRankingTree(
        output_phasespace_mc,
        NT,
        CATEGORY,
        "Chi2RankHybrid",
        "Chi2*1000",
        "1==1",
        "Run",
        "Event",
        "EnPB*1000");

    // now add the new friend tree, and select only the best chi2 combo per event per
    // beam photon. Save the result to a new file
    TString chi2_output_data = output_data;
    chi2_output_data.ReplaceAll("GENERAL", "bestChi2").ReplaceAll("general/", "best-chi2/");
    TString chi2_output_signal_mc = output_signal_mc;
    chi2_output_signal_mc.ReplaceAll("GENERAL", "bestChi2").ReplaceAll("general/", "best-chi2/");
    TString chi2_output_phasespace_mc = output_phasespace_mc;
    chi2_output_phasespace_mc.ReplaceAll("GENERAL", "bestChi2").ReplaceAll("general/", "best-chi2/");

    FSCut::defineCut("chi2rank", "Chi2RankHybrid==1");
    FSTree::addFriendTree("Chi2RankHybrid");
    FSModeTree::skimTree(
        output_data,
        NT,
        CATEGORY,
        chi2_output_data,
        "CUT(chi2rank)");
    FSModeTree::skimTree(
        output_signal_mc,
        NT,
        CATEGORY,
        chi2_output_signal_mc,
        "CUT(chi2rank)");
    FSModeTree::skimTree(
        output_phasespace_mc,
        NT,
        CATEGORY,
        chi2_output_phasespace_mc,
        "CUT(chi2rank)");
}