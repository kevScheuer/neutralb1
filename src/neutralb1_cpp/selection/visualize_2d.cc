/**
 * @file visualize_2d.cc
 * @author Kevin Scheuer
 * @brief Simple diagrams to show weighting regions for sideband subtraction
 * 
 * To better compare the norwegian method sideband subtraction and 2D method, we can
 * make some quick diagrams to visualize the weightings for the appropriate regions.
 */

#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TColor.h"
#include "TBox.h"
#include "TGraph.h"
#include "TPolyLine.h"
#include "TLegend.h"

TH2F* draw_2d_sideband_regions(
    double omega_mass,
    double signal_width, 
    double sideband_gap
);
TH2F* draw_norwegian_sideband_regions(
    double omega_mass,
    double signal_width,
    double sideband_gap
);

void visualize_2d()
{
    const double OMEGA_MASS = 0.7826;             // PDG omega mass in GeV
    const double SIGNAL_WIDTH = 0.00868 * 3;      // exactly 3 sigma of PDG omega width
    const double SIDEBAND_GAP = SIGNAL_WIDTH * 3; // 9 sigma distance from signal region

    TCanvas *c = new TCanvas("c", "c", 800, 800);  
    c->SetMargin(0.12, 0.12, 0.12, 0.12);

    TH2F *h = draw_2d_sideband_regions(
        OMEGA_MASS,
        SIGNAL_WIDTH,
        SIDEBAND_GAP        
    );
    h->SetStats(0);
    c->SaveAs("2d_sideband_regions.pdf");
    c->Clear();
    delete h;

    TH2F *h2 = draw_norwegian_sideband_regions(
        OMEGA_MASS,
        SIGNAL_WIDTH,
        SIDEBAND_GAP
    );
    h2->SetStats(0);
    c->SaveAs("norwegian_sideband_regions.pdf");
    c->Clear();

    return;
}

TH2F* draw_2d_sideband_regions(
    double omega_mass,
    double signal_width, 
    double sideband_gap
)
{
    // draw axes    
    TH2F *h_axes = new TH2F(
        "h_axes", 
        ";M(\\pi^{+}\\pi^{-}\\pi^{0}_{1}) (GeV);M(\\pi^{+}\\pi^{-}\\pi^{0}_{2}) (GeV)", 
        100, 0.5652, 1.0, 100, 0.5652, 1.0);
    h_axes->Draw(); 

    // draw y-axis lines
    double signal_low = omega_mass - signal_width;
    double signal_high = omega_mass + signal_width;
    double sideband_low_center = omega_mass - sideband_gap - signal_width * 0.5;
    double sideband_low_low = sideband_low_center - signal_width;
    double sideband_low_high = sideband_low_center + signal_width;
    double sideband_high_center = omega_mass + sideband_gap + signal_width * 0.5;
    double sideband_high_low = sideband_high_center - signal_width;
    double sideband_high_high = sideband_high_center + signal_width;

    double min_y = h_axes->GetYaxis()->GetXmin();
    double max_y = h_axes->GetYaxis()->GetXmax();
    double min_x = h_axes->GetXaxis()->GetXmin();
    double max_x = h_axes->GetXaxis()->GetXmax();

    TLine *line_y_signal_low = new TLine(signal_low, min_y, signal_low, max_y);
    TLine *line_y_signal_high = new TLine(signal_high, min_y, signal_high, max_y);
    TLine *line_y_sideband_low_low = new TLine(sideband_low_low, min_y, sideband_low_low, max_y);
    TLine *line_y_sideband_low_high = new TLine(sideband_low_high, min_y, sideband_low_high, max_y);
    TLine *line_y_sideband_high_low = new TLine(sideband_high_low, min_y, sideband_high_low, max_y);
    TLine *line_y_sideband_high_high = new TLine(sideband_high_high, min_y, sideband_high_high, max_y);
    line_y_signal_low->SetLineColor(kBlack);
    line_y_signal_low->SetLineWidth(1);
    line_y_signal_high->SetLineColor(kBlack);
    line_y_signal_high->SetLineWidth(1);
    line_y_sideband_low_low->SetLineColor(kBlack);
    line_y_sideband_low_low->SetLineWidth(1);
    line_y_sideband_low_low->SetLineStyle(9);
    line_y_sideband_low_high->SetLineColor(kBlack);
    line_y_sideband_low_high->SetLineWidth(1);   
    line_y_sideband_low_high->SetLineStyle(9);
    line_y_sideband_high_low->SetLineColor(kBlack);
    line_y_sideband_high_low->SetLineWidth(1);
    line_y_sideband_high_low->SetLineStyle(9);
    line_y_sideband_high_high->SetLineColor(kBlack);
    line_y_sideband_high_high->SetLineWidth(1);
    line_y_sideband_high_high->SetLineStyle(9);
    line_y_signal_low->Draw("SAME");
    line_y_signal_high->Draw("SAME");
    line_y_sideband_low_low->Draw("SAME");
    line_y_sideband_low_high->Draw("SAME");
    line_y_sideband_high_low->Draw("SAME");
    line_y_sideband_high_high->Draw("SAME");

    // draw x-axis lines
    TLine *line_x_signal_low = new TLine(min_x, signal_low, max_x, signal_low);
    TLine *line_x_signal_high = new TLine(min_x, signal_high, max_x, signal_high);
    TLine *line_x_sideband_low_low = new TLine(min_x, sideband_low_low, max_x, sideband_low_low);
    TLine *line_x_sideband_low_high = new TLine(min_x, sideband_low_high, max_x, sideband_low_high);
    TLine *line_x_sideband_high_low = new TLine(min_x, sideband_high_low, max_x, sideband_high_low);
    TLine *line_x_sideband_high_high = new TLine(min_x, sideband_high_high, max_x, sideband_high_high);

    line_x_signal_low->SetLineColor(kBlack);
    line_x_signal_low->SetLineWidth(1);
    line_x_signal_high->SetLineColor(kBlack);
    line_x_signal_high->SetLineWidth(1);
    line_x_sideband_low_low->SetLineColor(kBlack);
    line_x_sideband_low_low->SetLineWidth(1);
    line_x_sideband_low_low->SetLineStyle(9);
    line_x_sideband_low_high->SetLineColor(kBlack);
    line_x_sideband_low_high->SetLineWidth(1);
    line_x_sideband_low_high->SetLineStyle(9);
    line_x_sideband_high_low->SetLineColor(kBlack);
    line_x_sideband_high_low->SetLineWidth(1);
    line_x_sideband_high_low->SetLineStyle(9);
    line_x_sideband_high_high->SetLineColor(kBlack);
    line_x_sideband_high_high->SetLineWidth(1);
    line_x_sideband_high_high->SetLineStyle(9);

    line_x_signal_low->Draw("SAME");
    line_x_signal_high->Draw("SAME");
    line_x_sideband_low_low->Draw("SAME");
    line_x_sideband_low_high->Draw("SAME");
    line_x_sideband_high_low->Draw("SAME");
    line_x_sideband_high_high->Draw("SAME");

    // Fill sideband regions with hatched red
    TBox *box_sideband_low_x = new TBox(sideband_low_low, min_y, sideband_low_high, max_y);
    box_sideband_low_x->SetFillColorAlpha(kRed - 2, 0.35);
    box_sideband_low_x->SetFillStyle(3345);
    box_sideband_low_x->SetLineColor(0);
    box_sideband_low_x->Draw("SAME");

    TBox *box_sideband_high_x = new TBox(sideband_high_low, min_y, sideband_high_high, max_y);
    box_sideband_high_x->SetFillColorAlpha(kRed - 2, 0.35);
    box_sideband_high_x->SetFillStyle(3345);
    box_sideband_high_x->SetLineColor(0);

    box_sideband_high_x->Draw("SAME");

    TBox *box_sideband_low_y = new TBox(min_x, sideband_low_low, max_x, sideband_low_high);
    box_sideband_low_y->SetFillColorAlpha(kRed - 2, 0.35);
    box_sideband_low_y->SetFillStyle(3354);
    box_sideband_low_y->SetLineColor(0);
    box_sideband_low_y->Draw("SAME");

    TBox *box_sideband_high_y = new TBox(min_x, sideband_high_low, max_x, sideband_high_high);
    box_sideband_high_y->SetFillColorAlpha(kRed - 2, 0.35);
    box_sideband_high_y->SetFillStyle(3354);
    box_sideband_high_y->SetLineColor(0);
    box_sideband_high_y->Draw("SAME");

    // Fill both signal regions with semi-transparent blue, but avoid overlap

    // Signal region along x (excluding overlap)
    TBox *box_signal_x_left = new TBox(signal_low, min_y, signal_high, signal_low);
    box_signal_x_left->SetFillColor(kBlue);
    box_signal_x_left->SetLineColor(0);
    box_signal_x_left->Draw("SAME");

    TBox *box_signal_x_right = new TBox(signal_low, signal_high, signal_high, max_y);
    box_signal_x_right->SetFillColor(kBlue);
    box_signal_x_right->SetLineColor(0);
    box_signal_x_right->Draw("SAME");

    // Signal region along y (excluding overlap)
    TBox *box_signal_y_bottom = new TBox(min_x, signal_low, signal_low, signal_high);
    box_signal_y_bottom->SetFillColor(kBlue);
    box_signal_y_bottom->SetLineColor(0);
    box_signal_y_bottom->Draw("SAME");

    TBox *box_signal_y_top = new TBox(signal_high, signal_low, max_x, signal_high);
    box_signal_y_top->SetFillColor(kBlue);
    box_signal_y_top->SetLineColor(0);
    box_signal_y_top->Draw("SAME");

    // draw a blue X in overlap region to separate regions
    TPolyLine *x1 = new TPolyLine(2);
    x1->SetPoint(0, signal_low, signal_low);
    x1->SetPoint(1, signal_high, signal_high);
    x1->SetLineColor(kBlue);
    x1->SetLineWidth(1);
    x1->Draw("SAME");

    TPolyLine *x2 = new TPolyLine(2);
    x2->SetPoint(0, signal_low, signal_high);
    x2->SetPoint(1, signal_high, signal_low);
    x2->SetLineColor(kBlue);
    x2->SetLineWidth(1);
    x2->Draw("SAME");
    
    // Fill top and bottom of X (signal region where pi01 is chosen) with horizontal    
    TPolyLine *triangle_bottom = new TPolyLine(4);
    triangle_bottom->SetPoint(0, signal_low, signal_low);
    triangle_bottom->SetPoint(1, signal_high, signal_low);
    triangle_bottom->SetPoint(2, (signal_low + signal_high) / 2.0, (signal_low + signal_high) / 2.0);
    triangle_bottom->SetPoint(3, signal_low, signal_low);
    triangle_bottom->SetFillColor(kBlue);
    triangle_bottom->SetFillStyle(3350);         
    triangle_bottom->SetLineColor(kBlue);
    triangle_bottom->SetLineWidth(1);
    triangle_bottom->Draw("f SAME"); 
    triangle_bottom->Draw("SAME");   

    TPolyLine *triangle_top = new TPolyLine(4);
    triangle_top->SetPoint(0, signal_low, signal_high);
    triangle_top->SetPoint(1, signal_high, signal_high);
    triangle_top->SetPoint(2, (signal_low + signal_high) / 2.0, (signal_low + signal_high) / 2.0);
    triangle_top->SetPoint(3, signal_low, signal_high); 
    triangle_top->SetFillColor(kBlue);
    triangle_top->SetFillStyle(3350);
    triangle_top->SetLineColor(kBlue);
    triangle_top->SetLineWidth(1);
    triangle_top->Draw("f SAME"); 
    triangle_top->Draw("SAME");   

    // Fill left and right of X (signal region where pi02 is chosen) with vertical 
    // lines pattern
    TPolyLine *triangle_left = new TPolyLine(4);
    triangle_left->SetPoint(0, signal_low, signal_low);
    triangle_left->SetPoint(1, signal_low, signal_high);
    triangle_left->SetPoint(2, (signal_low + signal_high) / 2.0, (signal_low + signal_high) / 2.0);
    triangle_left->SetPoint(3, signal_low, signal_low);
    triangle_left->SetFillColor(kBlue);
    triangle_left->SetFillStyle(3359);
    triangle_left->SetLineColor(kBlue);
    triangle_left->SetLineWidth(1);
    triangle_left->Draw("f SAME");
    triangle_left->Draw("SAME");  

    // Fill right of X with vertical lines pattern
    TPolyLine *triangle_right = new TPolyLine(4);
    triangle_right->SetPoint(0, signal_high, signal_low);
    triangle_right->SetPoint(1, signal_high, signal_high);
    triangle_right->SetPoint(2, (signal_low + signal_high) / 2.0, (signal_low + signal_high) / 2.0);
    triangle_right->SetPoint(3, signal_high, signal_low);  
    triangle_right->SetFillColor(kBlue);
    triangle_right->SetFillStyle(3359);
    triangle_right->SetLineColor(kBlue);
    triangle_right->SetLineWidth(1);
    triangle_right->Draw("f SAME"); 
    triangle_right->Draw("SAME");   

    TLegend *legend = new TLegend(0.75, 0.75, 0.95, 0.95);
    legend->AddEntry(box_sideband_low_x, "-0.5", "f");
    legend->AddEntry(box_sideband_low_y, "-0.5", "f");
    legend->AddEntry(box_signal_x_left, "+1.0", "f");
    legend->AddEntry(triangle_right, "+1.0 (choose #pi^{0}_{2})", "f");
    legend->AddEntry(triangle_bottom, "+1.0 (choose #pi^{0}_{1})", "f");
    TBox *dummy_fill = new TBox(0, 0, 0, 0);
    dummy_fill->SetFillStyle(3244);
    dummy_fill->SetFillColor(kRed-2);
    legend->AddEntry(dummy_fill, "-5/8", "f");
    legend->Draw();

    return h_axes;
}

TH2F* draw_norwegian_sideband_regions(
    double omega_mass,
    double signal_width,
    double sideband_gap
)
{
    // draw axes    
    TH2F *h_ax = new TH2F(
        "h_axes", 
        ";M(\\pi^{+}\\pi^{-}\\pi^{0}_{1}) (GeV);M(\\pi^{+}\\pi^{-}\\pi^{0}_{2}) (GeV)", 
        100, 0.5652, 1.0, 100, 0.5652, 1.0);
    h_ax->Draw(); 

    
    double signal_low = omega_mass - signal_width;
    double signal_high = omega_mass + signal_width;    
    double sideband_low_low = omega_mass - sideband_gap;
    double sideband_low_high = omega_mass - sideband_gap - signal_width;    
    double sideband_high_low = omega_mass + sideband_gap;
    double sideband_high_high = omega_mass + sideband_gap + signal_width;

    double min_y = h_ax->GetYaxis()->GetXmin();
    double max_y = h_ax->GetYaxis()->GetXmax();
    double min_x = h_ax->GetXaxis()->GetXmin();
    double max_x = h_ax->GetXaxis()->GetXmax();

    Int_t norwegian_blue = TColor::GetColor("#00205B");
    Int_t norwegian_red = TColor::GetColor("#BA0C2F");

    // fill sideband regions with solid red
    TBox *box_sideband_low_x = new TBox(sideband_low_low, min_y, sideband_low_high, max_y);
    box_sideband_low_x->SetFillColor(norwegian_red);
    box_sideband_low_x->SetLineColor(0);
    box_sideband_low_x->Draw("SAME");

    TBox *box_sideband_high_x = new TBox(sideband_high_low, min_y, sideband_high_high, max_y);
    box_sideband_high_x->SetFillColor(norwegian_red);
    box_sideband_high_x->SetLineColor(0);
    box_sideband_high_x->Draw("SAME");

    TBox *box_sideband_low_y = new TBox(min_x, sideband_low_low, max_x, sideband_low_high);
    box_sideband_low_y->SetFillColor(norwegian_red);
    box_sideband_low_y->SetLineColor(0);
    box_sideband_low_y->Draw("SAME");

    TBox *box_sideband_high_y = new TBox(min_x, sideband_high_low, max_x, sideband_high_high);
    box_sideband_high_y->SetFillColor(norwegian_red);
    box_sideband_high_y->SetLineColor(0);
    box_sideband_high_y->Draw("SAME");

    TBox *box_signal_x = new TBox(signal_low, min_y, signal_high, max_y);
    box_signal_x->SetFillColor(norwegian_blue);
    box_signal_x->SetLineColor(0);
    box_signal_x->Draw("SAME");

    TBox *box_signal_y = new TBox(min_x, signal_low, max_x, signal_high);
    box_signal_y->SetFillColor(norwegian_blue);
    box_signal_y->SetLineColor(0);
    box_signal_y->Draw("SAME");

    TLegend *legend = new TLegend(0.75, 0.75, 0.85, 0.85);
    legend->AddEntry(box_sideband_low_x, "-1.0", "f");
    legend->AddEntry(box_signal_x, "+1.0", "f");
    legend->Draw();

    return h_ax;
}