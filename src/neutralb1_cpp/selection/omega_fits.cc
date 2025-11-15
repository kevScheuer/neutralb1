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

#include "FSBasic/FSHistogram.h"
#include "FSBasic/FSCut.h"
#include "FSBasic/FSTree.h"
#include "FSMode/FSModeHistogram.h"

#include "load_broad_cuts.cc"
#include "neutralb1/fit_utils.h"
#include "fsroot_setup.cc"

// forward declarations
std::vector<TH1F *> create_omega_histograms(
    TString NT,
    TString CATEGORY,
    TString cuts,
    TString input_files,
    double omega_pi0_bin_width,
    double t_bin_low_edge,
    double t_bin_high_edge);

std::map<std::pair<float, float>, std::vector<TFitResultPtr>> fit_and_plot_omega_histograms(
    std::map<std::pair<float, float>, std::vector<TH1F *>> omega_hist_map);
    
Double_t lin_voigt(Double_t *x, Double_t *par);

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
    std::tie(NT, CATEGORY) = setup(false);
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();
    TString cuts = join_keys(cut_color_map);

    // create all the histograms in mass bins for each t bin
    double omega_pi0_bin_width = 0.050; // GeV. from 1.0 to 2.0 GeV we'll have 20 bins
    std::map<std::pair<float, float>, std::vector<TH1F *>> omega_hist_map;
    std::vector<float> t_bin_edges = {0.1, 0.2, 0.3, 0.5, 1.0};

    for (size_t i = 0; i < t_bin_edges.size() - 1; ++i)
    {
        double t_bin_low_edge = t_bin_edges[i];
        double t_bin_high_edge = t_bin_edges[i + 1];
        std::vector<TH1F *> omega_hists = create_omega_histograms(
            NT,
            CATEGORY,
            cuts,
            input_files,
            omega_pi0_bin_width,
            t_bin_low_edge,
            t_bin_high_edge);
        omega_hist_map[std::make_pair(t_bin_low_edge, t_bin_high_edge)] = omega_hists;
        break; // TODO: remove break after testing
    }

    std::map<std::pair<float, float>, std::vector<TFitResultPtr>>
        fit_params_map = fit_and_plot_omega_histograms(omega_hist_map);

    // TODO: below function should use TFitResultPtr to extract fit parameters and plot them
    // plot_fit_parameters();

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
 * @return std::vector<TH1F *> omega mass histograms for each omega-pi0 mass bin
 */
std::vector<TH1F *> create_omega_histograms(
    TString NT,
    TString CATEGORY,
    TString cuts,
    TString input_files,
    double omega_pi0_bin_width,
    double t_bin_low_edge,
    double t_bin_high_edge)
{
    // Reminder, particles are ordered as:
    //   0       1          2             3              4              5
    // [beam] [proton] [pi+ (omega)] [pi- (omega)] [pi0 (omega)] [pi0 (bachelor)]
    // but pi0s 4 and 5 are interchangeable

    // we have two possible pi0 assignments for the omega decay
    TString data_omega_mass_1 = "MASS(2,3,4)";
    TString data_omega_mass_2 = "MASS(2,3,5)";

    std::vector<TH1F *> omega_histograms;
    for (float mass_bin_low = 1.0; mass_bin_low < 2.0; mass_bin_low += omega_pi0_bin_width)
    {
        float mass_bin_high = mass_bin_low + omega_pi0_bin_width;
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

        omega_histograms.push_back(h_omega_mass_1);
        break; // TODO: remove break after testing
    }
    return omega_histograms;
}

/**
 * @brief Fit and plot omega histograms for each t bin
 * 
 * @param omega_hist_map Map of t bin ranges to omega mass histograms
 * @return std::map<std::pair<float, float>, std::vector<TFitResultPtr>> Map of t bin ranges to fit results
 */
std::map<std::pair<float, float>, std::vector<TFitResultPtr>> fit_and_plot_omega_histograms(
    std::map<std::pair<float, float>, std::vector<TH1F *>> omega_hist_map)
{

    // for storing fit parameters
    std::map<std::pair<float, float>, std::vector<TFitResultPtr>> fit_params_map;

    // loop over t bins
    for (std::map<std::pair<float, float>, std::vector<TH1F *>>::const_iterator map_iter = omega_hist_map.begin();
         map_iter != omega_hist_map.end();
         ++map_iter)
    {
        // 4x5 canvas for each omega pi mass bin
        TCanvas *c1 = new TCanvas("c1", "Omega Fits", 1200, 1000);
        c1->Divide(5, 4); // 5 columns, 4 rows

        float t_bin_low_edge = map_iter->first.first;
        float t_bin_high_edge = map_iter->first.second;
        std::vector<TH1F *> omega_hists = map_iter->second;

        std::vector<TFitResultPtr> fit_params_vector;

        // loop over mass bins
        for (int i = 0; i < omega_hists.size(); ++i)
        {
            TH1F *h_omega_mass = omega_hists[i];

            int pad_index = i + 1; // pads are 1-indexed
            c1->cd(pad_index);
            h_omega_mass->Draw("E0");

            // fit the histogram with a voigtian
            TF1 *voigt_fit = new TF1("lin_voigt", lin_voigt, 0.6, 1.0);

            voigt_fit->SetParNames(""
                                   "Bkg offset",
                                   "Bkg slope",
                                   "M_{#omega}",
                                   "#sigma_{resolution}",
                                   "#Gamma_{#omega}");

            // initial parameter guesses
            voigt_fit->SetParameter(0, h_omega_mass->GetMinimum()); // background constant
            voigt_fit->SetParameter(1, 1.0);                        // positive background slope
            voigt_fit->SetParameter(2, 0.782);                      // omega mass
            voigt_fit->SetParLimits(2, 0.770, 0.790);               // limit mass to reasonable range
            voigt_fit->SetParameter(3, 0.010);                      // sigma (Gaussian width)
            voigt_fit->SetParameter(4, 0.00868);                    // gamma

            TFitResultPtr r = h_omega_mass->Fit(voigt_fit, "Q");
            fit_params_vector.push_back(r);

            h_omega_mass->GetFunction("lin_voigt")->Draw("SAME");

            // clean up
            delete voigt_fit;
        } // end loop over mass bins

        fit_params_map[std::make_pair(t_bin_low_edge, t_bin_high_edge)] = fit_params_vector;
        c1->SaveAs(TString::Format(
            "omega_fits_%.3f_%.3f.pdf", 
            t_bin_low_edge, t_bin_high_edge
        ).Data());
        delete c1;
    } // end loop over t bin to omega histograms map

    return fit_params_map;
}

/**
 * @brief Linear background plus Voigtian function for fitting omega mass peaks
 *
 * @param x Pointer to the array of x values (omega mass)
 * @param par Pointer to the array of parameters:
 *            par[0] = background constant offset
 *            par[1] = background slope
 *            par[2] = Voigtian peak position (mass)
 *            par[3] = Voigtian Gaussian width (sigma)
 *            par[4] = Voigtian Lorentzian width (gamma)
 * @return Double_t The value of the function at x
 */
Double_t lin_voigt(Double_t *x, Double_t *par)
{
    Double_t background = par[0] + par[1] * x[0];
    Double_t voigtian = TMath::Voigt(x[0], par[4], par[3]);

    return background + voigtian;
}