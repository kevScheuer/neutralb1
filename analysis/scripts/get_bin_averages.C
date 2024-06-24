/*Get averages and errors for each t bin in one mass bin, and output to csv file

There's no built in way to do this in ROOT, so here a histogram is made for each bin
from the tree and the mean and RMS is obtained
*/

#include <string>
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1F.h"

void get_bin_averages(TString file_name) {
    // setup details of selection
    std::vector<std::string> t_bins = {"0.1", "0.2", "0.3", "0.5", "0.9"};        
    std::string E_range[2]    = {"8.2", "8.8"};
    std::string mass_range[2] = {"1.22", "1.24"};    

    std::string E_cut = "E_Beam>"+E_range[0] + " && E_Beam<"+E_range[1];
    std::string m_cut = "M4Pi>"+mass_range[0] + " && M4Pi<"+mass_range[1];
    std::string selection = "Weight*("+m_cut+" && "+E_cut+")";
    
    // get file and tree
    std::unique_ptr<TFile> f(TFile::Open(file_name));
    if(!f) {std::cout << "File DNE" << "\n"; exit(1);}

    TTree *tree = f->Get<TTree>("kin");
    if(!tree) {std::cout << "Tree DNE" << "\n"; exit(1);}

    // setup csv file
    ofstream csv;
    std::string csv_name = "t-bin-mean_m-"+mass_range[0]+"-"+mass_range[1]+".csv"; 
    csv.open(csv_name);
    csv << "low_edge,high_edge,average,error\n"; // header line

    TCanvas *c1 = new TCanvas("c1", "c1", 1800, 1800);

    for(int i=0; i<(t_bins.size()-1); i++) {
        std::string bin_low = t_bins[i];
        std::string bin_high = t_bins[i+1];
        std::cout << bin_low << "\t" << bin_high << "\n";        
        tree->Draw(("t>>h(100,"+bin_low+","+bin_high+")").c_str(), selection.c_str());
        c1->SaveAs((bin_low+"-"+bin_high+".pdf").c_str());
        TH1F *h = (TH1F*)gPad->GetPrimitive("h");
        float mean = h->GetMean();        
        float rms = h->GetRMS();
        csv << bin_low <<","<< bin_high <<","<< mean <<"," << rms <<"\n";
    }    
    
    return;
}