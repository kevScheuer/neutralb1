/*Plot fit results for angular distributions and other variables of interest

There's no simple way to place these results with the .csv of all the other amplitude
fit results, so this script handles their plotting. The primary use of these is for
diagnostic purposes, to check that the grey "Fit Result" reasonably matches the black
data points. Each total J^P contribution is also plotted, to observe any interference
effects those waves create.

*/

#include <iostream>

#include "glueXstyle.C"
#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TLegend.h"

void plot1D(TFile *f, TString dir, TString data_title, TString reac);
void plot2D(TFile *f, TString dir, TString reac);

void angle_plotter(TString file_name = "vecps_plot.root",
                   TString dir = "./", TString data_title = "GlueX Data",
                   TString reac = "")
{
    gluex_style();
    gStyle->SetOptStat(0);

    TFile *f = TFile::Open(dir + file_name);
    if (!f)
    {
        cout << "File" << dir + file_name << "doesn't exist! Exiting" << "\n";
        exit(1);
    }

    plot1D(f, dir, data_title, reac);
    plot2D(f, dir, reac);

    return;
}

void plot1D(TFile *f, TString dir, TString data_title, TString reac)
{
    /* Make the 1D histograms of interest from vecps_plotter
     */
    std::vector<TString> plot_names = {"CosTheta", "Phi", "CosTheta_H", "Phi_H",
                                       "Prod_Ang", "MVecPs", "MProtonPs"};
    std::map<TString, std::map<TString, TString>> amp_map;

    // TODO: Iterate through all plots and extract all the JP values present in the
    // waveset, and add them to the map below. Probably best to iterate through eJPml
    // strings, that way I can list all the L waves that contribute to a JP

    amp_map["fit"]["plot_name"] = "";
    amp_map["fit"]["title"] = "Fit Result";
    amp_map["fit"]["draw_option"] = "f";
    amp_map["fit"]["color"] = "16";

    amp_map["1p"]["plot_name"] = "_0m";
    amp_map["1p"]["title"] = "#[]{0^{#minus}}^{(#pm)} (P)";
    amp_map["1p"]["draw_option"] = "ep";
    amp_map["1p"]["color"] = "6";

    amp_map["1p"]["plot_name"] = "_1p";
    amp_map["1p"]["title"] = "#[]{1^{#plus}}^{(#pm)} (S+D)";
    amp_map["1p"]["draw_option"] = "ep";
    amp_map["1p"]["color"] = "2";

    amp_map["1m"]["plot_name"] = "_1m";
    amp_map["1m"]["title"] = "#[]{1^{#minus}}^{(#pm)} P";
    amp_map["1m"]["draw_option"] = "ep";
    amp_map["1m"]["color"] = "4";

    // create canvas and setup some plot parameters
    TCanvas *cc = new TCanvas("cc", "cc", 1800, 1000);
    cc->Divide(3, 3);

    double textSize = 0.10;
    TLegend *leg1 = new TLegend(0.1, 0.1, 0.5, 0.9);
    leg1->SetEntrySeparation(0.01);
    leg1->SetNColumns(2);
    leg1->SetColumnSeparation(1.0);
    leg1->SetMargin(0.2);
    leg1->SetFillColor(0);
    leg1->SetTextSize(textSize);
    leg1->SetBorderSize(0);
    leg1->SetLineColor(kWhite);

    int plot_count = 0;
    for (auto plot_name : plot_names)
    {
        cc->cd(plot_count + 2);
        // first plot the data for this variable
        TH1F *hdat = (TH1F *)f->Get(reac + plot_name + "dat");
        if (!hdat)
        {
            cout << Form("hdat Plot %s doesn't exist! exiting", (reac + plot_name + "dat").Data()) << "\n";
            exit(1);
        }

        // since its the first in the subplot, setup all the aesthetics
        hdat->SetLineColor(kBlack);
        hdat->SetLabelSize(0.06, "xy");
        hdat->SetTitleSize(0.08, "xy");
        hdat->SetTitleOffset(0.88, "x");
        hdat->SetTitleOffset(0.9, "y");
        hdat->SetMinimum(0);
        hdat->SetMarkerStyle(20);
        hdat->SetMarkerSize(0.5);
        hdat->Draw();
        if (plot_count == 0)
            leg1->AddEntry(hdat, data_title, "ep");

        // now for each variable, we iterate through each JP's contribution that is
        // modified by the detector acceptance (acc)
        for (auto const &map_itr : amp_map)
        {
            auto key = map_itr.first;
            auto val = map_itr.second;
            TH1F *hacc = (TH1F *)f->Get(reac + plot_name + "acc" + val.at("plot_name"));

            if (!hacc)
            {
                cout << Form("hacc Plot %s/%s doesn't exist! exiting", plot_name.Data(), val.at("plot_name").Data()) << "\n";
                exit(1);
            }

            if (val.at("plot_name") == "")
            { // fit result plotted with gray fill style
                hacc->SetFillColorAlpha((val.at("color")).Atoi(), 0.2);
                hacc->Draw("same HIST");
            }
            else
            { // all other plots get markers
                hacc->SetMarkerColor((val.at("color")).Atoi());
                hacc->SetMarkerSize(0.6);
                hacc->SetMarkerStyle(20);
                hacc->Draw("same" + val.at("draw_option"));
            }
            hacc->SetLineColor((val.at("color")).Atoi());

            if (plot_count == 0)
                leg1->AddEntry(hacc, val.at("title"), val.at("draw_option"));
        }
        plot_count += 1;
        // hdat->Draw("same");
    }

    // finally draw the legend on the 1st subplot
    cc->cd(1);
    leg1->Draw();

    // save file
    cc->Print(dir + "fit.pdf");

    return;
}

void plot2D(TFile *f, TString dir, TString reac)
{
    /* Similar to plot1D, but the individual amplitude contributions cannot be seen, and
    so no legend is required here
    */

    // Note Psi = phi - Prod_Ang (Phi)
    std::vector<std::pair<TString, TString>> plot_pairs = {
        {"Psi", "CosTheta"},
        {"Psi", "CosTheta_H"},
        {"Psi", "Phi_H"},
        {"Prod_Ang", "Phi"},
        {"Phi", "CosTheta"},
        {"Phi_H", "CosTheta_H"},
        {"MProtonPs", "CosTheta"},
        {"MVecPs", "CosTheta"}};

    std::vector<TString> plot_titles = {
        ";#Psi [rad.];cos#theta",
        ";#Psi [rad.];cos#theta_{H}",
        ";#Psi [rad.];#phi_{H} [rad.]",
        ";Prod_Ang [rad.];#phi [rad.]",
        ";#phi [rad.];cos#theta",
        ";#phi_{H} [rad.];cos#theta_{H}",
        ";p+bachelor Ps [GeV];cos#theta",
        ";Vec+Ps [GeV];cos#theta",
    };

    // make 2d hists for just data first
    TCanvas *cdat = new TCanvas("cdat", "cdat", 1800, 1000);

    cdat->Divide(4, 2);

    int pair_count = 0;
    for (auto pair : plot_pairs)
    {
        cdat->cd(pair_count + 1);

        TString plot_name = pair.first + "Vs" + pair.second;

        TH2F *hdat = (TH2F *)f->Get(reac + plot_name + "dat");
        if (!hdat)
        {
            cout << Form("hdat Plot %s doesn't exist! exiting", (reac + plot_name + "dat").Data()) << "\n";
            exit(1);
        }
        hdat->SetLabelSize(0.06, "xy");
        hdat->SetTitleSize(0.08, "xy");
        hdat->SetTitleOffset(0.88, "x");
        hdat->SetTitleOffset(0.9, "y");
        hdat->SetMinimum(0);
        hdat->SetTitle(plot_titles[pair_count]);
        hdat->Draw("colz");

        pair_count += 1;
    }
    cdat->Print(dir + "fit2D_dat.pdf");

    // now do the same for the accepted fit result. Note that these have to be done in
    // a separate for-loop to void the canvases interfering with each other
    TCanvas *cacc = new TCanvas("cacc", "cacc", 1800, 1000);
    cacc->Divide(4, 2);

    pair_count = 0;
    for (auto pair : plot_pairs)
    {
        TString plot_name = pair.first + "Vs" + pair.second;
        cacc->cd(pair_count + 1);

        TH2F *hacc = (TH2F *)f->Get(reac + plot_name + "acc");

        if (!hacc)
        {
            cout << Form("hacc Plot %s doesn't exist! exiting", plot_name.Data()) << "\n";
            exit(1);
        }
        hacc->SetLabelSize(0.06, "xy");
        hacc->SetTitleSize(0.08, "xy");
        hacc->SetTitleOffset(0.88, "x");
        hacc->SetTitleOffset(0.9, "y");
        hacc->SetMinimum(0);
        hacc->SetTitle(plot_titles[pair_count]);
        hacc->Draw("colz");

        pair_count += 1;
    }

    cacc->Print(dir + "fit2D_acc.pdf");

    return;
}
