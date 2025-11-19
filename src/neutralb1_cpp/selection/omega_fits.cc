/**
 * @file omega_fits.cc
 * @author Kevin Scheuer
 * @brief Fit the omega mass in mass and t bins to determine signal widths
 *
 * Sideband subtraction needs a width for the signal region, and we can't simply assume
 * the omega width is sufficient due to detector resolution effects. This file fits
 * the omega mass peak with a voigtian in bins of omega-pi0 mass and t to determine an
 * appropriate width. The resolution and mass parameters are plotted as a function of
 * omega-pi0 mass and t to check how consistent they are.
 */

#include <iostream>
#include <tuple>

#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TColor.h"
#include "TMath.h"
#include "TGraphErrors.h"

#include "FSBasic/FSHistogram.h"
#include "FSBasic/FSCut.h"
#include "FSBasic/FSTree.h"
#include "FSMode/FSModeHistogram.h"

#include "load_broad_cuts.cc"
#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"

// forward declarations
std::map<double, TH1F *> create_omega_histograms(
    TString NT,
    TString CATEGORY,
    TString cuts,
    TString input_files,
    double omega_pi0_bin_width,
    double t_bin_low_edge,
    double t_bin_high_edge);

std::map<std::pair<float, float>, std::map<double, std::vector<double>>>
fit_and_plot_omega_histograms(
    std::map<std::pair<float, float>, std::map<double, TH1F *>> omega_hist_map,
    bool mc);
void plot_fit_parameters(
    std::map<std::pair<float, float>, std::map<double, std::vector<double>>> fit_params_map,
    bool mc);    

Double_t voigt(Double_t *x, Double_t *par);
Double_t linear(Double_t *x, Double_t *par);
Double_t cubic(Double_t *x, Double_t *par);
Double_t lin_omega_fit(Double_t *x, Double_t *par);
Double_t cubic_omega_fit(Double_t *x, Double_t *par);

void omega_fits(int period, bool mc = false)
{
    TString input_files;
    if (mc)
    {
        input_files = TString::Format(
            "/lustre24/expphy/volatile/halld/home/kscheuer/"
            "FSRoot-skimmed-trees/best-chi2/"
            "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_ver03.1_mc.root",
            period);
    }
    else
    {
        input_files = TString::Format(
            "/lustre24/expphy/volatile/halld/home/kscheuer/"
            "FSRoot-skimmed-trees/best-chi2/"
            "tree_pi0pi0pippim__B4_bestChi2_SKIM_0%d_data.root",
            period);
    }
    TString NT, CATEGORY;
    std::tie(NT, CATEGORY) = setup(true);
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();
    TString cuts = join_keys(cut_color_map);

    // create all the histograms in mass bins for each t bin
    double omega_pi0_bin_width = 0.050; // GeV. from 1.0 to 2.0 GeV we'll have 20 bins
    std::map<std::pair<float, float>, std::map<double, TH1F *>> t_to_omega_hist_map;
    std::vector<float> t_bin_edges = {0.1, 0.2, 0.3, 0.5, 1.0};

    for (size_t i = 0; i < t_bin_edges.size() - 1; ++i)
    {
        double t_bin_low_edge = t_bin_edges[i];
        double t_bin_high_edge = t_bin_edges[i + 1];
        std::map<double, TH1F *> omega_hist_map = create_omega_histograms(
            NT,
            CATEGORY,
            cuts,
            input_files,
            omega_pi0_bin_width,
            t_bin_low_edge,
            t_bin_high_edge);

        t_to_omega_hist_map[std::make_pair(t_bin_low_edge, t_bin_high_edge)] = omega_hist_map;
    }

    std::map<std::pair<float, float>, std::map<double, std::vector<double>>>
        fit_params_map = fit_and_plot_omega_histograms(t_to_omega_hist_map, mc);
    
    plot_fit_parameters(fit_params_map, mc);

    return;
}

/**
 * @brief Create omega mass histograms in bins of omega-pi0 mass for a given t bin
 *
 * @param NT FSRoot tree name
 * @param CATEGORY FSRoot category name
 * @param cuts Pre-defined selection cuts to apply
 * @param input_files Input ROOT files
 * @param omega_pi0_bin_width Bin width for omega-pi0 mass
 * @param t_bin_low_edge Lower edge of t bin
 * @param t_bin_high_edge Upper edge of t bin
 * @return std::map<double, TH1F*> center of omega-pi0 mass bin and the corresponding
 *  omega mass histogram in that bin
 */
std::map<double, TH1F *> create_omega_histograms(
    TString NT,
    TString CATEGORY,
    TString cuts,
    TString input_files,
    double omega_pi0_bin_width,
    double t_bin_low_edge,
    double t_bin_high_edge)
{
    // first we need to replace the t cut with the one for this bin
    TString t_var = "abs(-1*MASS2([proton],-GLUEXTARGET))";
    TString t_cut_str = TString::Format(
        "%s<%f&&%s>%f",
        t_var.Data(), t_bin_high_edge, 
        t_var.Data(), t_bin_low_edge);
    FSCut::defineCut("t_bin_cut", t_cut_str);
    if (cuts.Contains("t,"))
        cuts.ReplaceAll("t,", "t_bin_cut,");
    else if (cuts.Contains(",t"))
        cuts.ReplaceAll(",t", ",t_bin_cut");
    else
        throw std::runtime_error("[omega_fits::create_omega_histograms] Warning: t cut not found in cuts string!");        
    

    // Reminder, particles are ordered as:
    //   0       1          2             3              4              5
    // [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (omega)] [pi0 (bachelor)]
    // but pi0s 4 and 5 are interchangeable

    // we have two possible pi0 assignments for the omega decay
    TString data_omega_mass_1 = "MASS(2,3,4)";
    TString data_omega_mass_2 = "MASS(2,3,5)";

    std::map<double, TH1F *> omega_hist_map;
    for (float mass_bin_low = 1.0; mass_bin_low < 2.0; mass_bin_low += omega_pi0_bin_width)
    {
        float mass_bin_high = mass_bin_low + omega_pi0_bin_width;
        double mass_bin_center = mass_bin_low + omega_pi0_bin_width / 2.0;
        TH1F *h_omega_mass_1 = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            data_omega_mass_1,
            "(100, 0.6, 1.0)",
            TString::Format(
                "(MASS(2,3,4,5)<%f&&MASS(2,3,4,5)>%f)&&CUT(%s)&&CUTWT(rf)",
                mass_bin_high, mass_bin_low, cuts.Data()));
        TH1F *h_omega_mass_2 = FSModeHistogram::getTH1F(
            input_files,
            NT,
            CATEGORY,
            data_omega_mass_2,
            "(100, 0.6, 1.0)",
            TString::Format(
                "(MASS(2,3,4,5)<%f&&MASS(2,3,4,5)>%f)&&CUT(%s)&&CUTWT(rf)",
                mass_bin_high, mass_bin_low, cuts.Data()));

        // add both permutations together
        h_omega_mass_1->Add(h_omega_mass_2);

        // define some cosmetic properties easily here
        h_omega_mass_1->SetTitle(
            TString::Format("%.3f < M_{#omega#pi^{0}} < %.3f GeV", mass_bin_low, mass_bin_high));
        h_omega_mass_1->SetXTitle("#pi^{+}#pi^{-}#pi^{0}_{i} inv. mass (GeV)");
        double bin_width = get_bin_width(h_omega_mass_1);
        h_omega_mass_1->SetYTitle(TString::Format("Events / %.3f GeV", bin_width));
        h_omega_mass_1->SetStats(true);

        omega_hist_map[mass_bin_center] = h_omega_mass_1;
    }
    FSHistogram::dumpHistogramCache();
    return omega_hist_map;
}

/**
 * @brief Fit and plot omega histograms for each t bin
 *
 * @param omega_hist_map Map of t bin ranges to map of omega-pi0 mass bin centers to omega mass histograms
 * @return std::map<std::pair<float, float>, std::map<double, std::vector<double>>> Map of t bin ranges
 *  to omega pi0 mass bin centers to vector of voigtian fit results. The vector contains
 *  the fit parameters for each mass bin, and the double is a length 3 array of
 *  {mass, sigma, gamma}.
 */
std::map<std::pair<float, float>, std::map<double, std::vector<double>>> fit_and_plot_omega_histograms(
    std::map<std::pair<float, float>, std::map<double, TH1F *>> omega_hist_map,
    bool mc)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111); // show reduced chi2, errors, and values of fit params

    // for storing the voigtain fit parameters
    std::map<std::pair<float, float>, std::map<double, std::vector<double>>> fit_params_map;

    // loop over t bins
    for (std::map<std::pair<float, float>, std::map<double, TH1F *>>::const_iterator map_iter = omega_hist_map.begin();
         map_iter != omega_hist_map.end();
         ++map_iter)
    {
        // 4x5 canvas for each omega pi mass bin
        TCanvas *c1 = new TCanvas("c1", "Omega Fits", 1200, 1000);
        c1->Divide(5, 4); // 5 columns, 4 rows

        float t_bin_low_edge = map_iter->first.first;
        float t_bin_high_edge = map_iter->first.second;
        std::map<double, TH1F *> omega_hist_map = map_iter->second;

        // loop over mass bins
        int pad_index = 1; // keep track of which canvas pad we're on
        for (std::map<double, TH1F *>::const_iterator hist_iter = omega_hist_map.begin();
             hist_iter != omega_hist_map.end();
             ++hist_iter)
        {
            double bin_center = hist_iter->first;
            TH1F *h_omega_mass = hist_iter->second;

            c1->cd(pad_index);
            pad_index++;

            // Define the fit
            TF1 *omega_fit = new TF1("omega_fit", cubic_omega_fit, 0.6, 1.0, 8);
            omega_fit->SetLineColor(kMagenta);
            omega_fit->SetLineWidth(1);

            omega_fit->SetParNames(""
                                   "const",
                                   "lin",
                                   "quad",
                                   "cubic",
                                   "A",
                                   "M_{#omega}",
                                   "#sigma",
                                   "#Gamma_{#omega}");

            // initial parameter guesses
            double end_difference = h_omega_mass->GetBinContent(h_omega_mass->GetNbinsX()) - h_omega_mass->GetBinContent(1);

            omega_fit->SetParameters(           // cubic background + voigtian
                h_omega_mass->GetBinContent(1), // background constant
                end_difference,                 // background slope
                0.0,                            // background quadratic coeff
                0.0,                            // background cubic coeff
                h_omega_mass->GetMaximum(),     // voigt amplitude
                0.78266,                          // omega mass
                0.010,                          // sigma (Gaussian width)
                0.00868);                       // gamma (omega decay width)

            // parameter limits
            // widths should be positive and no wider than 100 MeV. The omega mass
            // should be reasonably near the PDG one
            omega_fit->SetParLimits(4, 0.0, h_omega_mass->GetMaximum() * 2); // amplitude
            omega_fit->FixParameter(5, 0.78266);                            // omega mass
            omega_fit->SetParLimits(6, 0.0, 0.5);                            // sigma
            omega_fit->FixParameter(7, 0.00868);                             // gamma fixed to PDG value for now

            h_omega_mass->Fit(omega_fit, "V", "e");
            h_omega_mass->Draw("E");
            omega_fit->Draw("SAME");

            // get background and signal by themselves
            TF1 *background = new TF1("background", cubic, 0.6, 1.0, 4);
            background->SetLineColor(kRed);
            background->SetLineWidth(1);
            TF1 *signal = new TF1("signal", voigt, 0.6, 1.0, 4);
            signal->SetLineColor(kBlue);
            signal->SetLineWidth(1);
            double par[8];

            omega_fit->GetParameters(par); // fill array with result parameters
            background->SetParameters(par);
            background->Draw("SAME");
            signal->SetParameters(&par[4]); // signal parameters start at index 4
            signal->Draw("SAME");

            // store voigtian fit parameters
            // (ignore the amplitude term, only want mass and widths)
            fit_params_map[std::make_pair(t_bin_low_edge, t_bin_high_edge)][bin_center] = {par[6], omega_fit->GetParError(6)}; // sigma, sigma_err
        } // end loop over mass bins

        TString mc_tag = mc ? "_mc" : "";

        c1->SetTitle(TString::Format("Omega Fits for %.1f < |t| < %.1f GeV^{2}", t_bin_low_edge, t_bin_high_edge));
        c1->SaveAs(TString::Format(
                       "omega_fits%s_%.3f_%.3f.pdf",
                       mc_tag.Data(), t_bin_low_edge, t_bin_high_edge)
                       .Data());
        delete c1;
    } // end loop over t bin to omega histograms map

    return fit_params_map;
}

/**
 * @brief Voigtian function for fitting omega mass peaks
 *
 * @param x Pointer to the array of x values (omega mass)
 * @param par Pointer to the array of parameters:
 *            par[0] = Voigtian peak amplitude
 *            par[1] = Voigtian peak position (mass)
 *            par[2] = Voigtian Gaussian width (sigma)
 *            par[3] = Voigtian Lorentzian width (gamma)
 * @return Double_t The value of the function at x
 */
Double_t voigt(Double_t *x, Double_t *par)
{
    return par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3]);
}

/**
 * @brief Linear function for fitting background
 *
 * @param x Pointer to the array of x values (omega mass)
 * @param par Pointer to the array of parameters:
 *            par[0] = background constant offset
 *            par[1] = background slope
 * @return Double_t The value of the function at x
 */
Double_t linear(Double_t *x, Double_t *par)
{
    return par[0] + par[1] * x[0];
}

/**
 * @brief Cubic function for fitting background
 *
 * @param x Pointer to the array of x values (omega mass)
 * @param par Pointer to the array of parameters:
 *            par[0] = background constant offset
 *            par[1] = background linear coefficient
 *            par[2] = background quadratic coefficient
 *            par[3] = background cubic coefficient
 * @return Double_t The value of the function at x
 */
Double_t cubic(Double_t *x, Double_t *par)
{
    return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
}

/**
 * @brief Linear background plus Voigtian function for fitting omega mass peaks
 *
 * @param x Pointer to the array of x values (omega mass)
 * @param par Pointer to the array of parameters:
 *            par[0] = background constant offset
 *            par[1] = background slope
 *            par[2] = Voigtian peak amplitude
 *            par[3] = Voigtian peak position (mass)
 *            par[4] = Voigtian Gaussian width (sigma)
 *            par[5] = Voigtian Lorentzian width (gamma)
 * @return Double_t The value of the function at x
 */
Double_t lin_omega_fit(Double_t *x, Double_t *par)
{
    return linear(x, par) + voigt(x, &par[2]);
}

/**
 * @brief Cubic background plus Voigtian function for fitting omega mass peaks
 *
 * @param x Pointer to the array of x values (omega mass)
 * @param par Pointer to the array of parameters:
 *            par[0] = background constant offset
 *            par[1] = background linear coefficient
 *            par[2] = background quadratic coefficient
 *            par[3] = background cubic coefficient
 *            par[4] = Voigtian peak amplitude
 *            par[5] = Voigtian peak position (mass)
 *            par[6] = Voigtian Gaussian width (sigma)
 *            par[7] = Voigtian Lorentzian width (gamma)
 * @return Double_t The value of the function at x
 */
Double_t cubic_omega_fit(Double_t *x, Double_t *par)
{
    return cubic(x, par) + voigt(x, &par[4]);
}


void plot_fit_parameters(std::map<std::pair<float, float>, std::map<double, std::vector<double>>> fit_params_map, bool mc)
{
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);

    // Create canvas
    TCanvas *c1 = new TCanvas("fit_params_canvas", "Omega Resolution vs Mass Bin Center", 800, 600);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);

    // Define colors and marker styles for different t-bins    
    int colors[4] = {kOrange-6, kGreen+1, kCyan+1, kViolet-8};
    int marker_styles[4] = {20, 21, 22, 23}; // circle, square, triangle, diamond

    std::vector<TGraphErrors*> sigma_graphs;
    std::vector<TString> legend_labels;

    int t_bin_index = 0;

    // Process each t-bin
    for (std::map<std::pair<float, float>, std::map<double, std::vector<double>>>::const_iterator map_iter = fit_params_map.begin();
         map_iter != fit_params_map.end();
         ++map_iter)
    {
        float t_bin_low_edge = map_iter->first.first;
        float t_bin_high_edge = map_iter->first.second;
        std::map<double, std::vector<double>> mass_to_fit_map = map_iter->second;

        // Create graph with error bars for sigma parameter
        TGraphErrors* sigma_graph = new TGraphErrors();

        int point_index = 0;
        for (std::map<double, std::vector<double>>::const_iterator mass_bin_iter = mass_to_fit_map.begin();
             mass_bin_iter != mass_to_fit_map.end();
             ++mass_bin_iter)
        {
            double mass_bin_center = mass_bin_iter->first;
            std::vector<double> fit_params = mass_bin_iter->second; // {sigma, sigma_err}

            if (fit_params.size() >= 2) {
                sigma_graph->SetPoint(point_index, mass_bin_center, fit_params[0]); // sigma parameter
                sigma_graph->SetPointError(point_index, 0, fit_params[1]); // sigma error
                point_index++;
            }
        }

        // Set graph properties
        int color = colors[t_bin_index];
        int marker_style = marker_styles[t_bin_index];

        sigma_graph->SetMarkerColor(color);
        sigma_graph->SetMarkerStyle(marker_style);
        sigma_graph->SetMarkerSize(1.0);
        sigma_graph->SetLineColor(color);

        sigma_graphs.push_back(sigma_graph);

        // Create legend label        
        legend_labels.push_back(TString::Format(
            "%.1f-%.1f",
            t_bin_low_edge,
            t_bin_high_edge
        ));

        t_bin_index++;
    }

    // Draw sigma parameters
    if (!sigma_graphs.empty()) {
        sigma_graphs[0]->SetMinimum(0);
        sigma_graphs[0]->Draw("AP");        
        sigma_graphs[0]->GetXaxis()->SetTitle("#omega#pi^{0} inv. mass (GeV)");
        sigma_graphs[0]->GetXaxis()->SetTitleSize(0.045);
        sigma_graphs[0]->GetXaxis()->SetLabelSize(0.04);
        sigma_graphs[0]->GetYaxis()->SetTitle("Gaussian #sigma (GeV)");
        sigma_graphs[0]->GetYaxis()->SetTitleSize(0.045);
        sigma_graphs[0]->GetYaxis()->SetLabelSize(0.04);

        for (size_t i = 1; i < sigma_graphs.size(); ++i) {
            sigma_graphs[i]->Draw("P SAME");
        }

        // Create legend
        TLegend *legend = new TLegend(0.65, 0.65, 0.9, 0.9);
        legend->SetHeader("-t bins (GeV^{2})", "C");
        for (size_t i = 0; i < sigma_graphs.size(); ++i) {
            legend->AddEntry(sigma_graphs[i], legend_labels[i].Data(), "p");
        }
        legend->Draw();
    }

    // Save the plot
    TString mc_tag = mc ? "_mc" : "";
    c1->SaveAs(TString::Format("omega_fit_parameters%s.pdf", mc_tag.Data()).Data());    

    // Clean up
    delete c1;

    return;
}