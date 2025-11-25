/**
 * @file integrate_voigt.cc
 * @author Kevin Scheuer
 * @brief Quick integration of omega signal (voigtian) to determine appropriate sideband subtraction widths
 * 
 * We want to integrate the function, fixed at the omega mass and width and using 
 * the optimal gaussian width determined from omega_fits.cc, to determine the fraction
 * of signal contained within a given window. We'll integrate symmetrically around the
 * omega mass, until we reach 95% of the area covered, equivalent to a 2-sigma cut
 * in a pure gaussian case. This is done because the voigtian has no obvious standard
 * deviation to use for cuts.
 */

#include <iostream>

#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"


/**
 * @brief Voigtian function for fitting omega mass peaks
 *
 * @param x Pointer to the array of x values (omega mass)
 * @param par Pointer to the array of parameters:
 *            par[0] = Voigtian peak position (mass)
 *            par[1] = Voigtian Gaussian width (sigma)
 *            par[2] = Voigtian Lorentzian width (gamma)
 * @return Double_t The value of the function at x
 */
Double_t voigt(Double_t *x, Double_t *par)
{
    return TMath::Voigt(x[0] - par[0], par[1], par[2]);
}

Double_t lorentzian(Double_t *x, Double_t *par)
{            
    double num = par[1] / (2.0 * TMath::Pi());
    double denom = (x[0] - par[0]) * (x[0] - par[0]) + (par[1] / 2.0) * (par[1] / 2.0);
    return num / denom;
}

void integrate_voigt()
{
    // omega mass and width parameters
    const double omega_mass = 0.78266; // PDG omega mass in GeV
    const double omega_gamma = 0.00868; // PDG omega width in GeV

    // detector resolution (gaussian sigma) determined from omega_fits.cc
    const double sigma = 0.018; // in GeV

    TF1 *function = new TF1("function", voigt, 0.6, 1.0, 3);
    function->SetParameters(omega_mass, sigma, omega_gamma);
    function->FixParameter(0, omega_mass);
    function->FixParameter(1, sigma);
    function->FixParameter(2, omega_gamma);

    // integrate from mass (- delta to mass) to  (+ delta) until we reach a % of area
    const double target_area = 0.9545;
    const double step_size = 0.001; // step 1 MeV at a time
    const double start = 0.010; // begin integration 10 MeV from peak
    
    double lower_result, upper_result;

    int iterations = 0;
    for (double delta = step_size; delta < 0.2; delta += step_size)
    {
        double lower_bound = omega_mass - start - delta;
        double upper_bound = omega_mass + start + delta;
        
        double area = function->Integral(lower_bound, upper_bound);
        double total_area = function->Integral(0.6, 1.0, 1.0E-6); // integrate over full range
        double fraction = area / total_area;

        if (fraction >= target_area)
        {
            std::cout << "Integrated area from " << lower_bound << " to " << upper_bound 
                      << " is " << fraction * 100 << "% of total area." << std::endl;
            std::cout << "Required window half-width: " << start + delta << " GeV." << std::endl;
            lower_result = lower_bound;
            upper_result = upper_bound;
            break;
        }
        iterations++;
        if (iterations > 500)
        {
            std::cout << "Reached maximum iterations without achieving target area." << std::endl;
            break;
        }
        if (iterations % 10 == 0)
        {
            std::cout << "Current delta: " << delta << " GeV, fraction: " << fraction * 100 << "%" << std::endl;
        }
    }

    // For plotting purposes, look at the gaussian and BW components separately
    TF1 *gauss = new TF1("gauss", "gaus", 0.6, 1.0);
    gauss->SetParameters(1, omega_mass, sigma);
    gauss->SetLineColor(kRed-7);
    gauss->SetLineWidth(1);
    gauss->SetLineStyle(2);
    gauss->SetNpx(1000);

    TF1 *bw = new TF1("bw", lorentzian, 0.6, 1.0, 2);
    bw->SetParameters(omega_mass, omega_gamma);
    bw->SetLineColor(kRed+2);
    bw->SetLineWidth(1);
    bw->SetLineStyle(3);
    bw->SetNpx(1000);

    // Draw the functions and shade the integrated area
    TCanvas *c1 = new TCanvas("c1", "Voigtian Integration", 800, 600);

    function->SetLineColor(kBlack);
    function->SetLineWidth(2);
    function->Draw();
    TH1D *h_voigt = (TH1D*)function->GetHistogram();
    h_voigt->SetTitle("");
    h_voigt->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    h_voigt->GetYaxis()->SetTitle("Arbitrary Units");

    gauss->Draw("same");
    bw->Draw("same");
    
    TH1D *h_area = (TH1D*)h_voigt->Clone("h_area");

    //zero out bins outside the integrated area
    for (int i = 1; i <= h_area->GetNbinsX(); ++i)
    {
        double bin_center = h_area->GetBinCenter(i);
        if (bin_center < lower_result || bin_center > upper_result)
        {
            h_area->SetBinContent(i, 0);;
        }
    }

    h_area->SetLineWidth(0);
    h_area->SetLineColor(kBlue);
    h_area->SetFillStyle(1001);
    h_area->SetFillColorAlpha(kBlue, 0.5);
    h_area->Draw("HIST same");

    // Add legend
    TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.85);
    leg->AddEntry(function, "Voigt Profile", "l");
    leg->AddEntry(gauss, "Gaussian Component", "l");
    leg->AddEntry(bw, "Lorentzian Component", "l");
    leg->AddEntry(h_area, "Integrated area", "f");
    leg->Draw();    

    h_voigt->SetMaximum(bw->GetMaximum() * 1.1);

    c1->Update();
    c1->SaveAs("voigt_integration.pdf");
}