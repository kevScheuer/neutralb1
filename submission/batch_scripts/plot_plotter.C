#include "glueXstyle.C"

void plot_plotter(TString fileName = "omegapi_plot.root", 
		  TString dir = "./", TString dataTitle = "GlueX Data",
		  TString reac = "") {

  gluex_style();
  gStyle->SetOptStat(0);

  const int maxPlots = 8;
  TString plotNames[maxPlots] = {"CosTheta","Phi","CosTheta_H","Phi_H",
				 "Prod_Ang","MVecPs","MProtonPs","MRecoilPs"};
  // Left out "MRecoil", in above plot for now to allow vecps mass to be plotted
  const int maxAmps = 3; // MUST be true max, otherwise legend/plots get buggy
  TString ampNames[maxAmps] = {"","_1p","_1m"}; //,"_2m","_3m"};
  TString ampTitles[maxAmps] = {"Fit Result","#[]{1^{#plus}}^{(#pm)} (S+D)",
				"#[]{1^{#minus}}^{(#pm)} P"}; 
  //,"#[]{2^{#minus}}^{(#pm)} (P+F)", "#[]{3^{#minus}}^{(#pm)} (F)"};
  TString ampDrawOpt[maxAmps] = {"h","ep","ep"}; //,"ep","ep"};
  int ampColors[maxAmps] = {1, 4, 2}; //, 7, 8};
    
  TFile *f = TFile::Open(dir+fileName);  
  if(!f) {cout << "File doesn't exist! Exiting" << "\n"; exit(1);}

  TCanvas *cc = new TCanvas("cc","cc",1800,1000);
  cc->Divide(3,3); 

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
    
  for(int i=0; i<maxPlots; i++) {
    TH1F *hdat = (TH1F*)f->Get(reac+plotNames[i]+"dat");
    if(!hdat) {
      cout << Form("hdat Plot i=%i doesn't exist! exiting", i) << "\n";
      exit(1);
    }

    //cout<<plotNames[i].Data()<<endl;
    hdat->SetLineColor(kBlack);
    hdat->SetLabelSize(0.06, "xy");
    hdat->SetTitleSize(0.08, "xy");
    hdat->SetTitleOffset(0.88, "x");
    hdat->SetTitleOffset(0.9, "y");
    hdat->SetMinimum(0);
    cc->cd(i+2);
    hdat->SetMarkerStyle(20);
    hdat->SetMarkerSize(0.5);
    hdat->Draw();
    if(i==0) leg1->AddEntry(hdat, dataTitle, "ep");
        
    for(int j=0; j<maxAmps; j++) {
      TH1F *hacc = (TH1F*)f->Get(reac+plotNames[i]+"acc"+ampNames[j]);
      //cout<<ampNames[j].Data()<<endl;        
		    
      if(!hacc) {
        cout << Form("hacc Plot j=%i doesn't exist! exiting", j) << "\n"; 
        exit(1);
      }
		    
      if(j==0) {
        hacc->SetFillColor(16);
        hacc->SetLineColor(16);
        hacc->Draw("same" + ampDrawOpt[j]);
        if(i==0) leg1->AddEntry(hacc, ampTitles[j], "f");
      }
      else {
        hacc->SetLineColor(ampColors[j]);
        hacc->SetMarkerColor(ampColors[j]);
        hacc->SetMarkerSize(0.6);
        hacc->SetMarkerStyle(20);
        hacc->Draw("same" + ampDrawOpt[j]);
        if(i==0) leg1->AddEntry(hacc, ampTitles[j], ampDrawOpt[j]);
      }
    }
        
    hdat->Draw("same");
  }
	
  cc->cd(1);
  leg1->Draw();

  cc->Print(dir+"fit.pdf");
    
  return;
}
