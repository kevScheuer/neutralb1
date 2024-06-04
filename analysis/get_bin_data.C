/*Output bin content of data to a csv file

This circumvents the need to edit the data reader of vecps_plotter to handle
the mass binning. This file insteads reads the AmpTools Tree, and makes the same TEM 
cuts I used (set in the script). I then bin it to match the fit bins, and weight using 
the "Weight" leaf of the Tree. 
*/

#include <string>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1F.h"

void get_bin_data() {
    // setup details of selection
    std::string t_range[2]    = {"0.4", "0.5"};
    std::string E_range[2]    = {"8.2", "8.8"};
    std::string mass_range[2] = {"1.0", "1.5"};
    double bin_width          = 0.01;

    int n_bins = (std::stod(mass_range[1])-
                  std::stod(mass_range[0])) / bin_width;
    
    // get file and tree
    std::unique_ptr<TFile> f(TFile::Open("All.root"));    
    if(!f) {cout << "File DNE" << "\n"; exit(1);}

    TTree *tree = f->Get<TTree>("kin");
    if(!tree) {cout << "Tree DNE" << "\n"; exit(1);}

    TCanvas *c1 = new TCanvas("c1", "c1", 1800, 1800);

    // Create list of cuts        
    std::string t_cut = "t>"+t_range[0] + " && t<"+t_range[1];
    std::string E_cut = "E_Beam>"+E_range[0] + " && E_Beam<"+E_range[1];

    std::string selection = "Weight*("+t_cut+" && "+E_cut+")";
    
    std::string h_str = "h(" + std::to_string(n_bins) + "," + mass_range[0] + 
                        "," + mass_range[1] + ")";    
    
    // draw histogram with correct bins & cuts, then grab it from pad
    tree->Draw(("M4Pi>>"+h_str).c_str(), selection.c_str());
    TH1F *h = (TH1F*)gPad->GetPrimitive("h");

    // create csv file
    ofstream csv;
    std::string f_name = "bin_data_t_"+t_range[0]+"-"+t_range[1]+"_"+
                         "m_"+mass_range[0]+"-"+mass_range[1]+"-"+
                         std::to_string(bin_width)+".csv";
    csv.open(f_name);
    
    // write header of bins to csv
    csv << "bin_contents,bin_error\n";

    // write bin contents to csv
    for(int bin=1; bin<=h->GetNbinsX(); bin++) {
        csv << h->GetBinContent(bin) << "," << h->GetBinError(bin) << "\n";        
    }
    
    return;
}