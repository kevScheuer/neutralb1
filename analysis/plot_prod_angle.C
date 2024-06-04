/* Plots the reflectivity contribution to the production angle Phi

The Phi production angle holds most of the shape for which reflectivitiy contributes
to it. This script plots Phi for multiple t bins (hard coded) split by data, fit result,
and +/- reflectivity using the .root files output by vecps_plotter
*/

#include <iostream>
#include <string>

#include "glueXstyle.C"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TMath.h"

void plot_prod_angle() {
    gluex_style();
    gStyle->SetOptStat(0);    

    TString parent_dir = "/lustre19/expphy/volatile/halld/home/kscheuer/ampToolsFits/"
    "omegapi/GlueXI/PARA_0-PARA_135-PERP_45-PERP_90/data/ver03/1m_1p_iso/";
    
    TFile *t1 = TFile::Open(parent_dir+"t_0.10-0.20/mass_1.200-1.250/vecps_plot.root");
    TFile *t2 = TFile::Open(parent_dir+"t_0.20-0.30/mass_1.200-1.250/vecps_plot.root");
    TFile *t3 = TFile::Open(parent_dir+"t_0.30-0.50/mass_1.200-1.250/vecps_plot.root");
    TFile *t4 = TFile::Open(parent_dir+"t_0.50-0.90/mass_1.200-1.250/vecps_plot.root");

    std::vector<TString> bin_names = {"0.1-0.2", "0.2-0.3", "0.3-0.5", "0.5-0.9"};
    std::vector<TFile*> file_vec;
    file_vec.push_back(t1);
    file_vec.push_back(t2);
    file_vec.push_back(t3);
    file_vec.push_back(t4);

    if(!t1 || !t2 || !t3 || !t4) {
        cout << "A data file does not exist!" << "\n"; exit(1);
    }

    // these trees already have TEM cuts applied, all thats needed is to plot 
    std::map<TString, std::map<TString, TString>> plot_map;

    plot_map["data"]["plot_name"] = "Prod_Angdat";    
    plot_map["data"]["color"] = "1";
    plot_map["data"]["marker_style"] = "6";
    plot_map["data"]["legend_title"] = "GlueX-I Data"; 

    plot_map["fit"]["plot_name"] = "Prod_Angacc";
    plot_map["fit"]["draw_option"] = "h";
    plot_map["fit"]["color"] = "16";
    plot_map["fit"]["legend_title"] = "Fit Result";

    plot_map["pos"]["plot_name"] = "Prod_Angacc_PosRefl_1p";
    plot_map["pos"]["draw_option"] = "ep";
    plot_map["pos"]["color"] = "#B51700";
    plot_map["pos"]["marker_style"] = "20";
    plot_map["pos"]["legend_title"] = "b_{1} (#varepsilon = +1)";

    plot_map["neg"]["plot_name"] = "Prod_Angacc_NegRefl_1p";
    plot_map["neg"]["draw_option"] = "ep";
    plot_map["neg"]["color"] = "#004D80";
    plot_map["neg"]["marker_style"] = "21";
    plot_map["neg"]["legend_title"] = "b_{1} (#varepsilon = -1)";

    TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
    TLegend *leg = new TLegend(0.65, 0.9, 0.92, 1.0);
    leg->SetEntrySeparation(0.01);
    leg->SetNColumns(2);
    leg->SetColumnSeparation(0.2);
    leg->SetMargin(0.2);
    leg->SetFillColor(0);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetLineColor(kWhite);

    int i=0;
    for(auto f : file_vec) {
        leg->Clear();
        cc->Clear();        

        // // draw data histogram and make aesthetic edits
        TH1F *hdat = (TH1F*)f->Get((plot_map.at("data").at("plot_name")));
        if(!hdat) {cout << "hdat doesn't exist!" << "\n"; exit(1);}
        hdat->SetLineColor((plot_map.at("data").at("color").Atoi()));
        hdat->SetLineWidth(3);
        // hdat->SetLabelSize(0.06, "xy");
        // hdat->SetTitleSize(0.08, "xy");
        // hdat->SetTitleOffset(0.88, "x");
        // hdat->SetTitleOffset(0.9, "y");
        hdat->SetMinimum(0);        
        hdat->SetMarkerStyle((plot_map.at("data").at("marker_style").Atoi()));
        hdat->SetMarkerSize(0.5);
        hdat->GetXaxis()->SetTitle("Prod Angle #Phi (rad)");
        hdat->GetYaxis()->SetTitle("Candidates / 0.13 rad");
        hdat->Draw();
        leg->AddEntry(hdat, (plot_map.at("data").at("legend_title")), "ep");

        // draw fit result
        TH1F *hacc = (TH1F*)f->Get((plot_map.at("fit").at("plot_name")));
        if(!hacc) {cout << "hacc doesn't exist!" << "\n"; exit(1);}
        hacc->SetFillColorAlpha((plot_map.at("fit").at("color").Atoi()), 0.5);
        hacc->SetLineColorAlpha((plot_map.at("fit").at("color").Atoi()), 0.5);
        hacc->Draw("same"+(plot_map.at("fit").at("draw_option")));
        leg->AddEntry(hacc, (plot_map.at("fit").at("legend_title")), "f");

        // draw positive reflectivity
        TH1F *hpos = (TH1F*)f->Get((plot_map.at("pos").at("plot_name")));
        if(!hpos) {cout << "hpos doesn't exist!" << "\n"; exit(1);}
        hpos->SetLineColor(TColor::GetColor((plot_map.at("pos").at("color"))));
        hpos->SetMarkerColor(TColor::GetColor((plot_map.at("pos").at("color"))));
        hpos->SetMarkerSize(0.6);
        hpos->SetMarkerStyle((plot_map.at("pos").at("marker_style").Atoi()));
        hpos->Draw("same"+(plot_map.at("pos").at("draw_option")));
        leg->AddEntry(hpos, (plot_map.at("pos").at("legend_title")));

        // draw negative reflectivity
        TH1F *hneg = (TH1F*)f->Get((plot_map.at("neg").at("plot_name")));
        if(!hneg) {cout << "hneg doesn't exist!" << "\n"; exit(1);}
        hneg->SetLineColor(TColor::GetColor((plot_map.at("neg").at("color"))));
        hneg->SetMarkerColor(TColor::GetColor((plot_map.at("neg").at("color"))));
        hneg->SetMarkerSize(0.6);
        hneg->SetMarkerStyle((plot_map.at("neg").at("marker_style").Atoi()));
        hneg->Draw("same"+(plot_map.at("neg").at("draw_option")));
        leg->AddEntry(hneg, (plot_map.at("neg").at("legend_title")));


        leg->Draw("same");
        cc->Print("prod_angles_" + bin_names[i] + ".pdf");
        i+=1;
    }


    return;
}