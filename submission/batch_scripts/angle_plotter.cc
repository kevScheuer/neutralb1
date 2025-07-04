/*Plot fit results for angular distributions and other variables of interest

There's no simple way to place these results with the .csv of all the other amplitude
fit results, so this script handles their plotting. The primary use of these is for
diagnostic purposes, to check that the grey "Fit Result" reasonably matches the black
data points. Each total J^P contribution is also plotted, to observe any interference
effects those waves create.

Usage: angle_plotter [file_path] [data_title] [reaction] [output_dir] [--gluex-style]
  file_path:  Full path to the ROOT file (default: "./vecps_plot.root")
  data_title: Title for data in legend (default: "GlueX Data")
  reaction:   Reaction prefix for histogram names (default: "")
  output_dir: Directory to save PDF files (default: current directory)
  --gluex-style: Apply GlueX collaboration style (optional)

Output PDFs will be saved in the current directory or specified output directory.

*/

#include <iostream>
#include <regex>
#include <set>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TList.h"
#include "TString.h"
#include "TStyle.h"


// struct to hold the properties of each JP contribution
struct JP_props
{
    TString legend;
    TColor *color;
    int marker;
};

// forward declarations
void plot1D(TFile *f, TString dir, TString data_title, TString reaction);
void plot2D(TFile *f, TString dir, TString reaction);
std::vector<TColor *> create_custom_colors();
const std::map<TString, JP_props> create_jp_map();

int main(int argc, char *argv[])
{
    std::string file_path = "./vecps_plot.root";
    std::string data_title = "GlueX Data";
    std::string reaction = "";
    std::string output_dir = "./";
    bool use_gluex_style = false;

    // Help message
    auto print_help = []() {
        std::cout << "Usage: angle_plotter [file_path] [data_title] [reaction] [output_dir] [--gluex-style]\n"
                  << "  file_path:     Full path to the ROOT file (default: ./vecps_plot.root)\n"
                  << "  data_title:    Title for data in legend (default: GlueX Data)\n"
                  << "  reaction:      Reaction prefix for histogram names (default: \"\")\n"
                  << "  output_dir:    Directory to save PDF files (default: current directory)\n"
                  << "  --gluex-style: Apply GlueX collaboration style (optional)\n";
    };

    // Check for help flag
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_help();
            return 0;
        }
    }

    if (argc > 6 || argc < 2) {
        std::cerr << "Error: Invalid number of arguments.\n";
        print_help();
        return 1;
    }

    if (argc > 1)
        file_path = argv[1];
    if (argc > 2)
        data_title = argv[2];
    if (argc > 3)
        reaction = argv[3];
    if (argc > 4)
        output_dir = argv[4];    
    
    // Check for --gluex-style flag in any position
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--gluex-style") {
            use_gluex_style = true;
            break;
        }
    }

    // Ensure output directory ends with a slash
    if (!output_dir.empty() && output_dir.back() != '/' && output_dir.back() != '\\')
    {
        output_dir += "/";
    }

    // Conditionally load and apply GlueX style
    if (use_gluex_style) {
        gROOT->ProcessLine(".L glueXstyle.C");
        gROOT->ProcessLine("gluex_style();");
    }
    gStyle->SetOptStat(0);

    TFile *f = TFile::Open(file_path.c_str());
    if (!f || f->IsZombie())
    {
        throw std::runtime_error(
            Form("File %s doesn't exist or is corrupted!", file_path.c_str()));
    }

    plot1D(f, output_dir, data_title, reaction);
    plot2D(f, output_dir, reaction);

    return 0;
}

// Make the 1D histograms of interest from vecps_plotter
void plot1D(TFile *f, TString dir, TString data_title, TString reaction)
{
    std::vector<TString> distributions = {"CosTheta", "Phi", "CosTheta_H", "Phi_H",
                                          "Prod_Ang", "MVecPs", "MProtonPs"};
    std::map<TString, TString> x_axis_titles = {
        {"CosTheta", "cos#theta"},
        {"Phi", "#phi [rad.]"},
        {"CosTheta_H", "cos#theta_{H}"},
        {"Phi_H", "#phi_{H} [rad.]"},
        {"Prod_Ang", "#Phi [rad.]"},
        {"MVecPs", "Vec+Ps Inv. Mass [GeV]"},
        {"MProtonPs", "p+bachelor Ps Inv. Mass [GeV]"}
    };
    const std::map<TString, JP_props> jp_properties = create_jp_map();

    // iterate through the histograms in the file and find the JP contributions
    std::set<TString> jp_keys_found;
    TList *keys = f->GetListOfKeys();
    for (int i = 0; i < keys->GetSize(); i++)
    {
        TString key_name = keys->At(i)->GetName();

        // regex pattern for capturing a JP key (number + letter after an "_")
        std::regex regex_pattern(".*_(\\d[a-zA-Z]).*");
        std::smatch match;
        std::string key_name_str = key_name.Data();
        if (std::regex_match(key_name_str, match, regex_pattern))
        {
            TString jp_key = match[1].str();
            jp_keys_found.insert(jp_key);
        }
    }

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
    for (auto distribution : distributions)
    {
        // avoid plotting on the first subplot, since it will be used for the legend
        cc->cd(plot_count + 2);
        
        // Set proper margins to prevent axis title cutoff
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.17);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);

        // first plot the data for this distribution
        TH1F *hdat = (TH1F *)f->Get(reaction + distribution + "dat");
        if (!hdat)
        {
            throw std::runtime_error(
                Form("hdat Plot %s doesn't exist!", (reaction + distribution + "dat").Data()));
        }

        // since its the first in the subplot, setup all the aesthetics
        hdat->SetTitle("");
        hdat->SetXTitle(x_axis_titles[distribution]);
        hdat->SetLineColor(kBlack);
        hdat->SetLabelSize(0.06, "xy");
        hdat->SetTitleSize(0.08, "xy");
        hdat->SetTitleOffset(0.88, "x");
        hdat->SetTitleOffset(1.0, "y");
        hdat->SetMinimum(0);
        hdat->SetMarkerStyle(20);
        hdat->SetMarkerSize(0.5);
        hdat->Draw();
        if (plot_count == 0)
            leg1->AddEntry(hdat, data_title, "ep");

        // now plot the fit result, which is a special case
        TH1F *hfit = (TH1F *)f->Get(reaction + distribution + "acc");
        if (!hfit)
        {
            throw std::runtime_error(
                Form("hfit Plot %s doesn't exist!", (reaction + distribution + "acc").Data()));
        }
        hfit->SetFillStyle(1001);
        hfit->SetFillColorAlpha(16, 0.2);
        hfit->SetLineColor(16);
        hfit->Draw("same HIST");
        if (plot_count == 0)
            leg1->AddEntry(hfit, "Fit Result", "f");

        // now we can loop through the JP contributions and plot them
        for (const TString &jp : jp_keys_found)
        {
            TString hist_name = reaction + distribution + "acc" + "_" + jp;
            TH1F *hjp = (TH1F *)f->Get(hist_name);

            if (!hjp)
            {
                throw std::runtime_error(
                    Form("Plot %s doesn't exist!", hist_name.Data()));
            }
            if (hjp->GetEntries() == 0 || hjp->Integral() == 0)
                continue;

            hjp->SetMarkerColor(jp_properties.at(jp).color->GetNumber());
            hjp->SetLineColor(jp_properties.at(jp).color->GetNumber());
            hjp->SetMarkerSize(0.6);
            hjp->SetMarkerStyle(jp_properties.at(jp).marker);
            hjp->Draw("same ep");

            if (plot_count == 0)
                leg1->AddEntry(hjp, jp_properties.at(jp).legend, "ep");
        }
        plot_count += 1;
    }

    // finally draw the legend on the 1st subplot
    cc->cd(1);
    leg1->Draw();

    // save file
    cc->Print(dir + "fit.pdf");

    return;
}

/* Similar to plot1D, but the individual amplitude contributions cannot be seen, and
so no legend is required here
*/
void plot2D(TFile *f, TString dir, TString reaction)
{
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
        
        // Set proper margins to prevent axis title cutoff
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.15);  // Extra space for colorbar
        gPad->SetTopMargin(0.05);

        TString plot_name = pair.first + "Vs" + pair.second;

        TH2F *hdat = (TH2F *)f->Get(reaction + plot_name + "dat");
        if (!hdat)
        {
            throw std::runtime_error(
                Form("hdat Plot %s doesn't exist!", (reaction + plot_name + "dat").Data()));
        }
        hdat->SetLabelSize(0.06, "xy");
        hdat->SetTitleSize(0.08, "xy");
        hdat->SetTitleOffset(1.2, "x");
        hdat->SetTitleOffset(1.3, "y");
        hdat->SetMinimum(0);
        hdat->SetTitle(plot_titles[pair_count]);
        hdat->Draw("colz");

        pair_count += 1;
    }
    cdat->Print(dir + "fit2D_dat.pdf");

    // now do the same for the accepted fit result. Note that these have to be done in
    // a separate for-loop to avoid the canvases interfering with each other
    TCanvas *cacc = new TCanvas("cacc", "cacc", 1800, 1000);
    cacc->Divide(4, 2);

    pair_count = 0;
    for (auto pair : plot_pairs)
    {
        TString plot_name = pair.first + "Vs" + pair.second;
        cacc->cd(pair_count + 1);
        
        // Set proper margins to prevent axis title cutoff
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.15);  // Extra space for colorbar
        gPad->SetTopMargin(0.05);

        TH2F *hacc = (TH2F *)f->Get(reaction + plot_name + "acc");

        if (!hacc)
        {
            throw std::runtime_error(
                Form("hacc Plot %s doesn't exist!", plot_name.Data()));
        }
        hacc->SetLabelSize(0.06, "xy");
        hacc->SetTitleSize(0.08, "xy");
        hacc->SetTitleOffset(1.2, "x");
        hacc->SetTitleOffset(1.3, "y");
        hacc->SetMinimum(0);
        hacc->SetTitle(plot_titles[pair_count]);
        hacc->Draw("colz");

        pair_count += 1;
    }

    cacc->Print(dir + "fit2D_acc.pdf");

    return;
}

// Mimic the colors of matplotlib's Dark2 colormap for plotting consistency
std::vector<TColor *> create_custom_colors()
{
    // rgb values from matplotlib's Dark2 colormap
    std::vector<std::vector<float>> mpl_Dark2_rgb_values = {
        {0.10588235294117647, 0.6196078431372549, 0.4666666666666667},
        {0.8509803921568627, 0.37254901960784315, 0.00784313725490196},
        {0.4588235294117647, 0.4392156862745098, 0.7019607843137254},
        {0.9058823529411765, 0.1607843137254902, 0.5411764705882353},
        {0.4, 0.6509803921568628, 0.11764705882352941},
        {0.9019607843137255, 0.6705882352941176, 0.00784313725490196},
        {0.6509803921568628, 0.4627450980392157, 0.11372549019607843},
        {0.4, 0.4, 0.4}};

    // Create TColor objects for each RGB value
    // and store them in a vector
    std::vector<TColor *> custom_colors;
    for (const auto &rgb : mpl_Dark2_rgb_values)
    {
        int color_index = TColor::GetFreeColorIndex();
        TColor *color = new TColor(color_index, rgb[0], rgb[1], rgb[2]);
        custom_colors.push_back(color);
    }
    return custom_colors;
}

/* Create a map of the JP contributions to their properties. Matches the jp_map
    used in pwa_tools.py for consistency.
*/
const std::map<TString, JP_props> create_jp_map()
{
    std::vector<TColor *> custom_colors = create_custom_colors();

    const std::map<TString, JP_props> jp_map = {
        {"0m", JP_props{"#[]{0^{#minus}}^{(#pm)} (P)", custom_colors[1], 5}},
        {"1p", JP_props{"#[]{1^{#plus}}^{(#pm)} (S+D)", custom_colors[2], 20}},
        {"1m", JP_props{"#[]{1^{#minus}}^{(#pm)} (P)", custom_colors[3], 21}},
        {"2m", JP_props{"#[]{2^{#minus}}^{(#pm)} (P+F)", custom_colors[4], 29}},
        {"2p", JP_props{"#[]{2^{#plus}}^{(#pm)} (D)", custom_colors[5], 34}},
        {"3m", JP_props{"#[]{3^{#minus}}^{(#pm)} (F)", custom_colors[6], 33}},
    };

    return jp_map;
}