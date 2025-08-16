/* Copy ROOT Tree to new directory, with data cuts applied to TEM regions

T = mandelstam t variable
E = beam energy (typically 8.2 - 8.8)
M = invariant mass
*/

#include <iostream>

#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"

int main(int argc, char *argv[])
{
    auto print_help = []()
    {
        std::cout << "Usage copy_tree_with_cuts [file] [output_path]"
                  << " [min_recoil_pi_mass]"
                  << " [t_low] [t_high] [E_low] [E_high] [m_low] [m_high]" << "\n"
                  << "file: Path to the ROOT file containing the tree\n"
                  << "output_path: Path to the output directory\n"
                  << "min_recoil_pi_mass: remove events below the recoil-proton + bachelor pion mass cut\n"
                  << "t_low: Lower cut on the Mandelstam t variable\n"
                  << "t_high: Upper cut on the Mandelstam t variable\n"
                  << "E_low: Lower cut on the photon beam energy\n"
                  << "E_high: Upper cut on the photon beam energy\n"
                  << "m_low: Lower cut on the invariant omega-pi0 mass\n"
                  << "m_high: Upper cut on the invariant omega-pi0 mass\n";
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
    TString min_recoil_pi_mass = argv[3];
    TString t_low = argv[4];
    TString t_high = argv[5];
    TString E_low = argv[6];
    TString E_high = argv[7];
    TString m_low = argv[8];
    TString m_high = argv[9];

    // open original file and get the kin tree which holds all the 4 vectors we need
    TFile *f_input = TFile::Open(file.c_str());
    if (!f_input || f_input->IsZombie())
    {
        std::cerr << "Error: Could not open file " << file << "\n";
        return 1;
    }

    TTree *input_tree = (TTree *)f_input->Get("kin");
    if (!input_tree)
    {
        std::cerr << "Error: Could not find tree 'kin' in file " << file << "\n";
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

    // setup selection windows
    // NOTE: these assume the kin trees have the leaf names used below
    TString cut_recoil_pi_mass = "MRecoilPi>" + min_recoil_pi_mass;
    TString t_cut = "t>" + t_low + " && t<" + t_high;
    TString E_cut = "E_Beam>" + E_low + " && E_Beam<" + E_high;
    TString m_cut = "M4Pi>" + m_low + " && M4Pi<" + m_high;

    TString selection = cut_recoil_pi_mass + " && " + t_cut + " && " + E_cut + " && " + m_cut;

    TTree *output_tree = input_tree->CopyTree(selection.Data());

    output_tree->Write();
    f_output->Close();

    return 0;
}