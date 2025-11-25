/**
 * @file skim_for_van_hove.cc
 * @author Kevin Scheuer
 * @brief Create skimmed FSRoot trees for more convenient van Hove analysis
 * 
 * This code applies broad cuts and the sideband subtraction method to create a signal
 * and background skimmed tree for each run period. This tree can then be used for van
 * Hove analysis and other studies without having to reapply cuts and sideband
 * subtraction each time.
 * 
 */

#include <tuple>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include "FSMode/FSModeTree.h"

#include "fsroot_setup.cc"


void skim_for_van_hove(int period)
{
    TString input_data_file = TString::Format(
        "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/best-chi2/"
        "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_data.root",
        period);
    TString input_mc_file = TString::Format(
        "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/best-chi2/"
        "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_ver03.1_mc.root",
        period);

    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(false);

    TString output_base_data = TString::Format(
        "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/van_hove/"
        "tree_pi0pi0pippim__B4_vanHove_SKIM_0%d_data.root", 
        period);
    TString output_base_mc = TString::Format(
        "/lustre24/expphy/volatile/halld/home/kscheuer/"
        "FSRoot-skimmed-trees/van_hove/"
        "tree_pi0pi0pippim__B4_vanHove_SKIM_0%d_ver03.1_mc.root", 
        period);

    for(int p=1; p<=2; p++) {
        // Note that cut==1 applies all broad cuts
        // signal==1 selects signal in the omega and rf regions
        // sideband==1 selects sideband in the omega and rf regions
        FSModeTree::skimTree(
            TString::Format("%s.permutation_%d", input_data_file.Data(), p),
            TString::Format("%s_permutation_%d", NT.Data(), p),
            CATEGORY,
            TString::Format("%s.permutation_%d_signal", output_base_data.Data(), p),
            "CUT(cut==1,signal==1)");
        FSModeTree::skimTree(
            TString::Format("%s.permutation_%d", input_data_file.Data(), p),
            TString::Format("%s_permutation_%d", NT.Data(), p),
            CATEGORY,
            TString::Format("%s.permutation_%d_sideband", output_base_data.Data(), p),
            "CUT(cut==1,sideband==1)");
        // rename the tree to remove the permutation tag
        TFile *f = TFile::Open(
            TString::Format("%s.permutation_%d_sideband", output_base_data.Data(), p).Data(),
            "UPDATE");
        TTree *t = (TTree*)f->Get(TString::Format("%s_permutation_%d", NT.Data(), p).Data());
        t->SetName(NT.Data());
        t->Write("",TObject::kOverwrite);
        f->Close();

        // repeat for mc
        FSModeTree::skimTree(
            TString::Format("%s.permutation_%d", input_mc_file.Data(), p),
            TString::Format("%s_permutation_%d", NT.Data(), p),
            CATEGORY,
            TString::Format("%s.permutation_%d_signal", output_base_mc.Data(), p),
            "CUT(cut==1,signal==1)");
        FSModeTree::skimTree(
            TString::Format("%s.permutation_%d", input_mc_file.Data(), p),
            TString::Format("%s_permutation_%d", NT.Data(), p),
            CATEGORY,
            TString::Format("%s.permutation_%d_sideband", output_base_mc.Data(), p),
            "CUT(cut==1,sideband==1)");
        f = TFile::Open(
            TString::Format("%s.permutation_%d_sideband", output_base_mc.Data(), p).Data(),
            "UPDATE");
        t = (TTree*)f->Get(TString::Format("%s_permutation_%d", NT.Data(), p).Data());
        t->SetName(NT.Data());
        t->Write("",TObject::kOverwrite);
        f->Close();
    }
    
    // now hadd the two permutations together for signal and sideband
    // data
    TString signal_output = output_base_data;
    signal_output.ReplaceAll(".root", "_signal.root");
    TString sideband_output = output_base_data;
    sideband_output.ReplaceAll(".root", "_sideband.root");
    system(TString::Format(
        "hadd -f %s %s.permutation_1_signal %s.permutation_2_signal",
        signal_output.Data(),
        output_base_data.Data(),
        output_base_data.Data()
    ).Data());
    system(TString::Format(
        "hadd -f %s %s.permutation_1_sideband %s.permutation_2_sideband",
        sideband_output.Data(),
        output_base_data.Data(),
        output_base_data.Data()
    ).Data());
    
    // signal mc
    signal_output = output_base_mc;
    signal_output.ReplaceAll(".root", "_signal.root");
    sideband_output = output_base_mc;
    sideband_output.ReplaceAll(".root", "_sideband.root");
    system(TString::Format(
        "hadd -f %s %s.permutation_1_signal %s.permutation_2_signal",
        signal_output.Data(),
        output_base_mc.Data(),
        output_base_mc.Data()
    ).Data());
    system(TString::Format(
        "hadd -f %s %s.permutation_1_sideband %s.permutation_2_sideband",
        sideband_output.Data(),
        output_base_mc.Data(),
        output_base_mc.Data()
    ).Data());
}
