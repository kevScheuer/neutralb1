/*This file loads in the reaction filtered trees and applies "skimming" cuts

The skim cuts are standard loose cuts, intended to reduce file size while keeping almost
all of the signal region
*/

TString NT("ntFSGlueX_MODECODE");
TString cut_variables;

void setup()
{
    if (FSModeCollection::modeVector().size() != 0)
        return;
    FSHistogram::readHistogramCache();
    FSModeCollection::addModeInfo("100_112")->addCategory("pi0pi0pippim");

    // FIXED CUTS
    FSCut::defineCut("eBeam", "(EnPB>8.2&&EnPB<8.8)");
    FSCut::defineCut("chi2", "Chi2DOF<20");

    // Select prompt RF region
    FSCut::defineCut("rf", "abs(RFDeltaT)<2.004");

    // Select omega region
    FSCut::defineCut("omega", "MASS(2,3,4)<1.0 || MASS(2,3,5)<1.0");

    // SET SOME DEFAULT SKIM CUTS
    cut_variables = "eBeam,chi2,rf,omega";
}

void create_skim_trees(int period = 3, bool is_bggen = false)
{
    // TODO: bggen appears unused, ask Justin if the commented out bggen lines are linked to this

    setup();

    const char* parent_cache_dir = "/cache/halld/home/jrsteven/flattened/omegapi_gx1_pwa/ver01/tree_pi0pi0pippim__B4";
    const char* target_dir = "/w/halld-scshelf2101/kscheuer/neutralb1/selector/skimmed_trees";
    TString reaction = "pi0pi0pippim";
    const char* tree_name = "tree_pi0pi0pippim__B4"; // B4 indicates mass of ?? is constrained

    // Get all the data files
    TString FND_DATA = Form("%s/data/%s_FSROOT_%02d*.root", parent_cache_dir, tree_name, period);
    TString FND_BGGEN = Form("/volatile/halld/home/jrsteven/flattened/%s/*bggen*/tree_pi0pi0pippim__B4_FSROOT_0*.root", tree_name);
    TString FND_SIGMC = Form("%s/omegapi_massDepFit_*_ver03.1/%s_FSROOT_%02d*.root", parent_cache_dir, tree_name, period);
    TString FND_PSMC = Form("%s/omegapi_phasespace_*_ver03/%s_FSROOT_%02d*.root", parent_cache_dir, tree_name, period);

    // Skim the trees by applying loose cuts
    FSModeTree::skimTree(
        FND_DATA, NT, reaction,
        Form("%s/%s_GENERAL_SKIM_%02d.root", target_dir, tree_name, period),
        "CUT(" + cut_variables + ")");
    // FSModeTree::skimTree(
    //     FND_BGGEN, NT, reaction,
    //     Form("%s_BGGEN_GENERAL_SKIM.root", tree_name),
    //     "CUT(" + cut_variables + ")");
    FSModeTree::skimTree(
        FND_SIGMC, NT, reaction,
        Form("%s/%s_SIGMC_GENERAL_SKIM_%02d.root", target_dir, tree_name, period),
        "CUT(" + cut_variables + ")");
    FSModeTree::skimTree(
        FND_PSMC, NT, reaction,
        Form("%s/%s_PHASESPACE_GENERAL_SKIM_%02d.root", target_dir, tree_name, period),
        "CUT(" + cut_variables + ")");

    // Create the associated chi2 rank trees from the skimmed ones
    // TODO: why is the rf cut applied again here?
    FSModeTree::createChi2RankingTree(
        Form("%s/%s_GENERAL_SKIM_%02d.root", target_dir, tree_name, period),
        NT, reaction, "CUT(rf)");
    // FSModeTree::createChi2RankingTree(
    //     Form("./%s_BGGEN_GENERAL_SKIM.root", tree_name),
    //     NT, reaction, "CUT(rf)");
    FSModeTree::createChi2RankingTree(
        Form("%s/%s_SIGMC_GENERAL_SKIM_%02d.root", target_dir, tree_name, period),
        NT, reaction, "CUT(rf)");
    FSModeTree::createChi2RankingTree(
        Form("%s/%s_PHASESPACE_GENERAL_SKIM_%02d.root", target_dir, tree_name, period),
        NT, reaction, "CUT(rf)");

    // skim only best chi2 candidates
    FSCut::defineCut("chi2rank", "Chi2Rank==1");
    FSTree::addFriendTree("Chi2Rank");

    FSModeTree::skimTree(
        Form("%s/%s_GENERAL_SKIM_%02d.root", target_dir, tree_name, period),
        NT, reaction,
        Form("%s/%s_BestChi2_SKIM_%02d.root", tree_name, target_dir, period),
        "CUT(chi2rank)");
    FSModeTree::skimTree(
        Form("%s/%s_SIGMC_GENERAL_SKIM_%02d.root", target_dir, tree_name, period),
        NT, reaction,
        Form("%s/%s_SIGMC_BestChi2_SKIM_%02d.root", tree_name, target_dir, period),
        "CUT(chi2rank)");
    FSModeTree::skimTree(
        Form("%s/%s_PHASESPACE_GENERAL_SKIM_%02d.root", target_dir, tree_name, period),
        NT, reaction,
        Form("%s/%s_PHASESPACE_BestChi2_SKIM_%02d.root", tree_name, target_dir, period),
        "CUT(chi2rank)");

    return;
}