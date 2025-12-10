/**
 * @file copy_tree_with_cuts.cc
 * @author Kevin Scheuer
 * @brief Copy ROOT Tree to new directory, with specified cuts applied
 * 
 * Primarily used for splitting the data file into TEM bins for PWA, but can be used
 * for systematics on the other kinematic variables as well. Note that it can only make
 * cuts on branches already present in the tree, though one can use the 4-momenta to
 * calculate other variables if needed.
 */

#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"


int main(int argc, char *argv[])
{
    auto print_help = []()
    {
        std::cout << "Usage: copy_tree_with_cuts [file] [output_path] [tree_name]"
                  << " [t_low] [t_high] [E_low] [E_high] [m_low] [m_high]"
                  << " [optional_cuts...]" << "\n\n"
                  << "Required arguments:\n"
                  << "  file: Path to the ROOT file containing the tree\n"
                  << "  output_path: Path to the output directory\n"
                  << "  tree_name: Name of the tree to copy\n"
                  << "  t_low: Lower cut on the Mandelstam t variable\n"
                  << "  t_high: Upper cut on the Mandelstam t variable\n"
                  << "  E_low: Lower cut on the photon beam energy\n"
                  << "  E_high: Upper cut on the photon beam energy\n"
                  << "  m_low: Lower cut on the invariant omega-pi0 mass\n"
                  << "  m_high: Upper cut on the invariant omega-pi0 mass\n\n"
                  << "Optional cuts (can specify multiple):\n"
                  << "  Each optional cut should be in the format: variable:min:max\n"
                  << "  Available variables: rf, M4Pi, unusedE, unusedTracks, z, MM2, missingE, chi2, MRecoilPi, PzCMrecoilPi\n"
                  << "  Examples:\n"
                  << "    chi2:0:5          (cut on Chi2DOF between 0 and 5)\n"
                  << "    MM2:-0.05:0.05    (cut on missing mass squared)\n"
                  << "    unusedE:0:0.1     (cut on unused shower energy)\n";
    };
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help")
        {
            print_help();
            return 0;
        }
    }

    if (argc < 10)
    {
        std::cerr << "Error: Invalid number of arguments.\n";
        print_help();
        return 1;
    }

    std::string file = argv[1];
    std::string output_path = argv[2];
    std::string tree_name = argv[3];
    TString t_low = argv[4];
    TString t_high = argv[5];
    TString E_low = argv[6];
    TString E_high = argv[7];
    TString m_low = argv[8];
    TString m_high = argv[9];
    
    // Parse optional cuts from remaining arguments
    std::vector<std::string> optional_cuts;
    for (int i = 10; i < argc; i++)
    {
        optional_cuts.push_back(argv[i]);
    }

    // open original file and get the tree which holds all the 4 vectors we need
    TFile *f_input = TFile::Open(file.c_str());
    if (!f_input || f_input->IsZombie())
    {
        std::cerr << "Error: Could not open file " << file << "\n";
        return 1;
    }

    TTree *input_tree = (TTree *)f_input->Get(tree_name.c_str());
    if (!input_tree)
    {
        std::cerr << "Error: Could not find tree '" << tree_name << "' in file " << file << "\n";
        return 1;
    }

    // we'll use the original file basename to copy to the new file path
    size_t pos = file.find_last_of("/\\");
    std::string basename = (pos == std::string::npos) ? file : file.substr(pos + 1);

    // Ensure output_path ends with "/"
    if (!output_path.empty() && output_path.back() != '/')
        output_path += '/';

    std::string output_file = output_path + basename;

    std::cout << "Copying Tree for: " << file << " to " << output_file << "\n";

    // open new file
    TFile *f_output = new TFile(output_file.c_str(), "new");

    // setup required selection windows
    // NOTE: these assume the kin trees have the leaf names used below
    TString t_cut = "t>" + t_low + " && t<" + t_high;
    TString E_cut = "EnPB>" + E_low + " && EnPB<" + E_high;
    TString m_cut = "M4Pi>" + m_low + " && M4Pi<" + m_high;

    TString selection = t_cut + " && " + E_cut + " && " + m_cut;
    
    // Parse and add optional cuts
    for (const auto& cut_str : optional_cuts)
    {
        size_t first_colon = cut_str.find(':');
        size_t second_colon = cut_str.find(':', first_colon + 1);
        
        if (first_colon == std::string::npos || second_colon == std::string::npos)
        {
            std::cerr << "Warning: Skipping malformed optional cut: " << cut_str << "\n";
            std::cerr << "  Expected format: variable:min:max\n";
            continue;
        }
        
        std::string var_name = cut_str.substr(0, first_colon);
        std::string min_val = cut_str.substr(first_colon + 1, second_colon - first_colon - 1);
        std::string max_val = cut_str.substr(second_colon + 1);
        
        TString optional_cut = TString::Format(
            "%s>%s && %s<%s",
            var_name.c_str(),
            min_val.c_str(),
            var_name.c_str(),
            max_val.c_str()
        );
        
        selection += " && " + optional_cut;
        std::cout << "  Adding optional cut: " << optional_cut << "\n";
    }
    
    std::cout << "Final selection: " << selection << "\n";

    TTree *output_tree = input_tree->CopyTree(selection.Data());

    output_tree->Write();
    f_output->Close();

    return 0;
}