/**
 * @file angle_plotter.cc
 * @author Kevin Scheuer
 * @brief Plot data + fit angular distributions and mass variables from PWA
 * 
 * There's no simple way to place these results with the .csv of all the other amplitude
 * fit results, so this script handles their plotting. The primary use of these is for
 * diagnostic purposes, to check that the grey "Fit Result" reasonably matches the black
 * data points. Each total J^P contribution is also plotted, to observe any interference
 * effects those waves create.

 * Usage: angle_plotter [file_path] [data_title] [output_dir]
 * file_path:  Full path to the ROOT file (default: "./vecps_plot.root")
 * data_title: Title for data in legend (default: "GlueX Data")
 * output_dir: Directory to save PDF files (default: current directory) 
 *
 * Output PDFs will be saved in the current directory or specified output directory.
 * 
 */

#include <iostream>
#include <regex>
#include <set>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TList.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include "neutralb1/AmplitudeParser.h"

// struct to hold the properties of each JP contribution
struct JP_props
{
    TString legend;
    TColor *color;
    int marker;
};

// forward declarations
TCanvas *summary_plot(TFile *f, TString data_title);
std::map<std::string, TCanvas *> amplitude_plots(TFile *f);
void plot2D(TFile *f, TString dir);
std::vector<TColor *> create_custom_colors();
const std::map<TString, JP_props> create_jp_map();
const std::set<TString> get_unique_jp_keys(TFile *f);
bool is_background_present(TFile *f);
const std::set<std::string> get_unique_amp_keys(TFile *f);

int main(int argc, char *argv[])
{
    std::string file_path = "./vecps_plot.root";
    std::string data_title = "GlueX Data";
    std::string output_dir = "./";    

    // Help message
    auto print_help = []()
    {
        std::cout << "Usage: angle_plotter [-f file_path] [-l legend_title] [-o output_dir]\n"
                  << "  -f file_path:     Full path to the ROOT file (default: ./vecps_plot.root)\n"
                  << "  -l legend_title:  Title for data in legend (default: GlueX Data)\n"
                  << "  -o output_dir:    Directory to save PDF files (default: current directory)\n"
                  << "  -h, --help:       Show this help message\n";
    };

    // Parse command-line arguments using flag style
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help")
        {
            print_help();
            return 0;
        }
        else if ((arg == "-f" || arg == "--file-path") && i + 1 < argc)
            file_path = argv[++i];
        else if ((arg == "-l" || arg == "--legend-title") && i + 1 < argc)
            data_title = argv[++i];
        else if ((arg == "-o" || arg == "--output-dir") && i + 1 < argc)
            output_dir = argv[++i];
        else
        {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
            print_help();
            return 1;
        }
    }

    // Ensure output directory ends with a slash
    if (!output_dir.empty() && output_dir.back() != '/' && output_dir.back() != '\\')
    {
        output_dir += "/";
    }

    gStyle->SetOptStat(0); // force stats box off

    TFile *f = TFile::Open(file_path.c_str());
    if (!f || f->IsZombie())
    {
        throw std::runtime_error(
            Form("File %s doesn't exist or is corrupted!", file_path.c_str()));
    }

    TCanvas *summary_canvas = summary_plot(f, data_title);
    std::map<std::string, TCanvas *> amp_canvases = amplitude_plots(f);
    TCanvas *test_canvas = summary_plot(f, data_title);

    TString pdf_path = output_dir + "angles.pdf";

    summary_canvas->Print(pdf_path + "(", "Title:Summary");
    for (std::map<std::string, TCanvas *>::iterator it = amp_canvases.begin(); it != amp_canvases.end(); ++it)
    {
        const std::string &amp = it->first;
        TCanvas *canvas = it->second;
        canvas->Print(pdf_path, Form("Title:%s", amp.c_str()));
    }
    test_canvas->Print(pdf_path + ")", "Title:Last Plot");

    plot2D(f, output_dir);

    return 0;
}

/**
 * @brief Create a set of summary plots from the vecps_plotter output
 *
 * This function plots the fit result, data, and JP contributions for the 5 angular
 * distributions we fit to, the resonance mass, proton + recoil pi0 mass, and lambda
 * distribution.
 *
 * @param f ROOT file ouptut by vecps_plotter
 * @param data_title title for data markers to use in the legend
 * @return TCanvas* ROOT canvas with subdivisions for each plot
 */
TCanvas *summary_plot(TFile *f, TString data_title)
{
    std::vector<TString> distributions = {
        "CosTheta",
        "Phi",
        "CosTheta_H",
        "Phi_H",
        "Prod_Ang",
        "MVecPs",
        "MProtonPs",
        "Lambda"};
    std::map<TString, TString> x_axis_titles = {
        {"CosTheta", "cos#theta"},
        {"Phi", "#phi [rad.]"},
        {"CosTheta_H", "cos#theta_{H}"},
        {"Phi_H", "#phi_{H} [rad.]"},
        {"Prod_Ang", "#Phi [rad.]"},
        {"MVecPs", "Vec+Ps Inv. Mass [GeV]"},
        {"MProtonPs", "p+bachelor Ps Inv. Mass [GeV]"},
        {"Lambda", "#lambda_{#omega}"}};
    const std::map<TString, JP_props> jp_properties = create_jp_map();

    // iterate through the histograms in the file and find the JP contributions
    std::set<TString> jp_keys = get_unique_jp_keys(f);

    // if background is present, add it to the set of jp_keys
    if (is_background_present(f))
        jp_keys.insert("Bkgd");

    // create a unique canvas and setup some plot parameters
    static int canvas_count = 0;
    TCanvas *cc = new TCanvas(Form("cc_%d", canvas_count++), "summary plots", 1800, 1000);
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
        TH1F *hdat = (TH1F *)f->Get(distribution + "dat");
        if (!hdat)
        {
            throw std::runtime_error(
                Form("hdat Plot %s doesn't exist!", (distribution + "dat").Data()));
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
        TH1F *hfit = (TH1F *)f->Get(distribution + "acc");
        if (!hfit)
        {
            throw std::runtime_error(
                Form("hfit Plot %s doesn't exist!", (distribution + "acc").Data()));
        }
        hfit->SetFillStyle(1001);
        hfit->SetFillColorAlpha(16, 0.2);
        hfit->SetLineColor(16);
        hfit->Draw("same HIST");
        if (plot_count == 0)
            leg1->AddEntry(hfit, "Fit Result", "f");

        // next, if a background file is present, plot it
        TH1F *hbkgd = (TH1F *)f->Get(distribution + "bkgd");
        if (hbkgd)
        {
            hbkgd->SetLineColor(kRed);
            hbkgd->SetLineWidth(1);
            hbkgd->SetMarkerStyle(0);
            hbkgd->SetFillStyle(1001);   
            hbkgd->SetFillColorAlpha(kRed, 0.3);         
            hbkgd->Draw("same HIST");
            if (plot_count == 0)
                leg1->AddEntry(hbkgd, "Background", "l");
        }

        // now we can loop through the JP contributions and plot them
        for (const TString &jp : jp_keys)
        {
            TString hist_name = distribution + "acc" + "_" + jp;
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

    return cc;
}

/**
 * @brief Create plots of 2D angular distributions for each amplitude
 *
 * @param f ROOT file output by vecps_plotter
 * @return std::map<std::string, TCanvas *> eJPmL amplitude names, corresponding canvas pointers
 */
std::map<std::string, TCanvas *> amplitude_plots(TFile *f)
{
    std::map<std::string, TCanvas *> amplitude_canvases;

    std::set<std::string> amp_keys = get_unique_amp_keys(f);

    // determine the set of 2D plots we want to see
    // Note Psi = phi - Prod_Ang (Phi)
    std::vector<std::pair<TString, TString>> plot_pairs = {
        {"Psi", "CosTheta"},
        {"Phi", "CosTheta"},
        {"Psi", "Phi_H"},
        {"Prod_Ang", "Phi"},
        {"Psi", "CosTheta_H"},
        {"Phi_H", "CosTheta_H"},
        {"MProtonPs", "CosTheta"},
    }; // TODO: Make this last plot the Dalitz plot of omegapi0 vs p recoil pi0

    std::vector<TString> plot_titles = {
        ";#Psi [rad.];cos#theta",
        ";#phi [rad.];cos#theta",
        ";#Psi [rad.];#phi_{H} [rad.]",
        ";#Phi [rad.];#phi [rad.]",
        ";#Psi [rad.];cos#theta_{H}",
        ";#phi_{H} [rad.];cos#theta_{H}",
        ";p+bachelor Ps [GeV];cos#theta",
    };

    static int plot_count = 0;
    for (const std::string &amp : amp_keys)
    {
        // create a unique canvas for each amplitude
        TCanvas *cc = new TCanvas(Form("AmpPlots_%d", plot_count++), "AmpPlots", 2400, 1000);

        cc->Divide(4, 2);

        int pair_count = 0;
        for (auto pair : plot_pairs)
        {
            cc->cd(pair_count + 1);

            // Set proper margins to prevent axis title cutoff
            gPad->SetLeftMargin(0.2);
            gPad->SetBottomMargin(0.22);
            gPad->SetRightMargin(0.2); // Extra space for colorbar
            gPad->SetTopMargin(0.15);

            TString plot_name = pair.first + "Vs" + pair.second;
            TH2F *hamp = (TH2F *)f->Get(plot_name + "acc_" + amp);
            if (!hamp)
            {
                throw std::runtime_error(
                    Form("hamp Plot %s doesn't exist!", (plot_name + "acc_" + amp).Data()));
            }
            hamp->SetLabelSize(0.06, "xy");
            hamp->SetTitleSize(0.08, "xy");
            hamp->SetTitleOffset(1.2, "x");
            hamp->SetTitleOffset(1.3, "y");
            hamp->SetMinimum(0);
            hamp->SetTitle(plot_titles[pair_count]);
            hamp->Draw("colz");
            pair_count += 1;
        }

        // save amplitude name in eJPmL format. Currently in "<sign>Refl_JPmD" format
        char refl;
        if (amp.find("PosRefl") != std::string::npos)
            refl = 'p';
        else if (amp.find("NegRefl") != std::string::npos)
            refl = 'm';
        else
        {
            throw std::runtime_error(
                Form("Amplitude %s does not contain PosRefl or NegRefl!", amp.c_str()));
        }

        // we want to use the eJPmL string since that's what will be saved as the title
        // in the pdf ToC, which allows easy access later and follows other analysis
        // conventions
        std::string eJPmL = refl + amp.substr(amp.find('_') + 1);

        cc->cd();
        // draw the eJPmL name on the canvas itself as a title
        TPad *title_pad = new TPad("title_pad", "title_pad", 0, 0, 1, 1);
        title_pad->SetFillStyle(4000);
        title_pad->Draw();
        title_pad->cd();
        TLatex *lat = new TLatex();
        AmplitudeParser parser(eJPmL[0], eJPmL[1], eJPmL[2], eJPmL[3], eJPmL[4]);
        lat->DrawLatexNDC(.45, .95, parser.get_latex_amplitude().c_str());

        amplitude_canvases[eJPmL] = cc;
    }

    return amplitude_canvases;
}

// TODO: change this to plot either total data, acc, or the subtraction of the two
// Also have it return canvas or set of canvases
// Also reorder according to amplitude_plots

/* Similar to plot1D, but the individual amplitude contributions cannot be seen, and
so no legend is required here
*/
void plot2D(TFile *f, TString dir)
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
        ";#Phi [rad.];#phi [rad.]",
        ";#phi [rad.];cos#theta",
        ";#phi_{H} [rad.];cos#theta_{H}",
        ";p+bachelor Ps [GeV];cos#theta",
        ";Vec+Ps [GeV];cos#theta",
    };

    // make 2d hists for just data first
    TCanvas *cdat = new TCanvas("cdat", "cdat", 2400, 1000);

    cdat->Divide(4, 2);

    int pair_count = 0;
    for (auto pair : plot_pairs)
    {
        cdat->cd(pair_count + 1);

        // Set proper margins to prevent axis title cutoff
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
        gPad->SetRightMargin(0.2); // Extra space for colorbar
        gPad->SetTopMargin(0.05);

        TString plot_name = pair.first + "Vs" + pair.second;

        TH2F *hdat = (TH2F *)f->Get(plot_name + "dat");
        if (!hdat)
        {
            throw std::runtime_error(
                Form("hdat Plot %s doesn't exist!", (plot_name + "dat").Data()));
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
    TCanvas *cacc = new TCanvas("cacc", "cacc", 2400, 1000);
    cacc->Divide(4, 2);

    pair_count = 0;
    for (auto pair : plot_pairs)
    {
        TString plot_name = pair.first + "Vs" + pair.second;
        cacc->cd(pair_count + 1);

        // Set proper margins to prevent axis title cutoff
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
        gPad->SetRightMargin(0.2); // Extra space for colorbar
        gPad->SetTopMargin(0.05);

        TH2F *hacc = (TH2F *)f->Get(plot_name + "acc");

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
        {"Bkgd", JP_props{"Iso. Bkgd", custom_colors[0], 1}},
        {"0m", JP_props{"#[]{0^{#minus}}^{(#pm)} (P)", custom_colors[1], 5}},
        {"1p", JP_props{"#[]{1^{#plus}}^{(#pm)} (S+D)", custom_colors[2], 20}},
        {"1m", JP_props{"#[]{1^{#minus}}^{(#pm)} (P)", custom_colors[3], 21}},
        {"2m", JP_props{"#[]{2^{#minus}}^{(#pm)} (P+F)", custom_colors[4], 29}},
        {"2p", JP_props{"#[]{2^{#plus}}^{(#pm)} (D)", custom_colors[5], 34}},
        {"3m", JP_props{"#[]{3^{#minus}}^{(#pm)} (F)", custom_colors[6], 33}},
    };

    return jp_map;
}

bool is_background_present(TFile *f)
{
    TString key = "isotropic";
    TList *keys = f->GetListOfKeys();
    for (int i = 0; i < keys->GetSize(); i++)
    {
        TString key_name = keys->At(i)->GetName();
        if (key_name.Contains(key))
            return true;
    }
    return false;
}

/**
 * @brief Extract the set of unique jp keys from all histograms
 *
 * @param f ROOT file output by vecps_plotter
 * @return const std::set<TString> unique set of jp strings
 * @note Assumes "JP" values can be found as a number+letter following an "_" character
 */
const std::set<TString> get_unique_jp_keys(TFile *f)
{
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

    return jp_keys_found;
}

/**
 * @brief Extract the set of unique amplitude strings from all histograms
 *
 * @param f ROOT file output by vecps_plotter
 * @return const std::set<TString> unique set of amplitude strings
 * @note Assumes "Amp" values can be found as "word with 'refl' in it, followed by an
 * underscore and a word that begins with a number and is at least 3 characters long "
 */
const std::set<std::string> get_unique_amp_keys(TFile *f)
{
    std::set<std::string> amp_keys_found;

    TList *keys = f->GetListOfKeys();
    for (int i = 0; i < keys->GetSize(); i++)
    {
        std::string key_name = keys->At(i)->GetName();

        // regex pattern for capturing a JP key (number + letter after an "_")
        std::regex regex_pattern("([A-Za-z]+Refl)_([0-9]\\w{2,})");
        std::smatch match;
        if (std::regex_search(key_name, match, regex_pattern))
        {
            std::string amp_key = match[1].str() + "_" + match[2].str();
            amp_keys_found.insert(amp_key);
        }
    }

    return amp_keys_found;
}