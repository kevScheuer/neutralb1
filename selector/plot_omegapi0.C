TString NT("ntFSGlueX_MODECODE");
TString DEFAULT_CUTS = "";
TString REACTION = "pi0pi0pippim"

    std::map<TString, std::map<TString, TString>>
    setup()
{
    if (FSModeCollection::modeVector().size() != 0)
        return;
    FSHistogram::readHistogramCache();
    FSModeCollection::addModeInfo("100_112")->addCategory(REACTION);

    // Map of cuts to loop over.
    //      value = mathematical expression of values to cut on
    //      tree_variable = actual variable name in tree
    //      binning = histogram binning in ROOT format (# bins, start, stop)
    //      units = for plotting the axis names
    std::map<TString, std::map<TString, TString>> cut_map;

    cut_map["unusedE"]["value"] = "EnUnusedSh<0.1";
    cut_map["unusedE"]["tree_variable"] = "EnUnusedSh";
    cut_map["unusedE"]["binning"] = "(100,0.0,1.0)";
    cut_map["unusedE"]["units"] = "[GeV]";

    cut_map["z"]["value"] = "ProdVz>=51.2&&ProdVz<=78.8";
    cut_map["z"]["tree_variable"] = "ProdVz";
    cut_map["z"]["binning"] = "(100,0.,100.0)";
    cut_map["z"]["units"] = "[GeV]";

    cut_map["MM2"]["value"] = "abs(RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5))<0.05";
    cut_map["MM2"]["tree_variable"] = "RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5)";
    cut_map["MM2"]["binning"] = "(100,-0.1,0.1)";
    cut_map["MM2"]["units"] = "[GeV^{2}]";

    cut_map["eBeam"]["value"] = "(EnPB>8.2&&EnPB<8.8)";
    cut_map["eBeam"]["tree_variable"] = "EnPB";
    cut_map["eBeam"]["binning"] = "(125,5,12)";
    cut_map["eBeam"]["units"] = "[GeV]";

    cut_map["chi2"]["value"] = "Chi2DOF<5";
    cut_map["chi2"]["tree_variable"] = "Chi2DOF";
    cut_map["chi2"]["binning"] = "(40,0,20)";
    cut_map["chi2"]["units"] = "";

    cut_map["t"]["value"] = "abs(-1*MASS2([proton],-GLUEXTARGET))<1.0";
    cut_map["t"]["tree_variable"] = "abs(-1*MASS2([proton],-GLUEXTARGET))";
    cut_map["t"]["binning"] = "(100,0,5)";
    cut_map["t"]["units"] = "[GeV^{2}]";

    // cuts we won't plot
    cut_map["unusedTracks"]["value"] = "NumUnusedTracks<1";
    cut_map["chi2rank"]["value"] = "Chi2Rank==1";
    cut_map["shQuality"]["value"] = "ShQualityP4a>0.5&&ShQualityP4b>0.5&&ShQualityP5a>0.5&&ShQualityP5b>0.5";

    for (it = cut_map.begin(); it != cut_map.end(); it++)
    {
        FSCut::defineCut(it->first, cut_map[it->first]["value"]);
        DEFAULT_CUTS += it->first + ",";
    }

    return cut_map;
}

void plot_omegapi0()
{
    // Basic plots for flattened and skimmed DATA trees:
    TString FND_DATA = "./selector/skimmed_trees/tree_pi0pi0pippim__B4_GENERAL_SKIM_05.root";
    TString FND_MC = FND_DATA;

    FSTree::addFriendTree("Chi2Rank");

    auto cut_map = setup();
    system("rm -rf ./selector/plots");
    system("mkdir ./selector/plots");

    TString CUTS;
    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->Divide(3, 2); // TODO: automate this

    uint i = 1; // enumeration for canvases
    for (it = cut_map.begin(); it != cut_map.end(); it++)
    {
        TString name = it->first;

        if (name == "chi2rank" || name == "shQuality" || name == "unusedTracks")
            continue; // these variables aren't interesting to plot

        c1->cd(i);
        cuts = DEFAULT_CUTS;             // reset cuts we'll use to default
        cuts.ReplaceAll("," + name, ""); // remove cut we want to plot

        // plot the data and set the axis titles
        TH1F *h = FSModeHistogram::getTH1F(
            FND_DATA, NT, REACTION,
            cut_map[name]["variable"], cut_map[name]["binning"],
            Form("CUT(%s)", cuts.Data()));

        h->SetXtitle(Form("%s %s", cut_map[name]["Xtitle"], cut_map[name]["units"]));
        h->SetYtitle(Form("Events / %.2f %s", h->GetBinWidth(), cut_map[name]["units"]));
        h->Draw();

        // Draw the MC with it in red, scaled to match the maximum values
        TH1F *h_mc = FSModeHistogram::getTH1F(
            FND_MC, NT, REACTION,
            cut_map[name]["variable"], cut_map[name]["binning"],
            Form("CUT(%s)", cuts.Data()));
        h_mc->Scale(h->GetMaximum() / h_mc->GetMaximum());
        h_mc->MarkerColor(kRed);
        h_mc->Draw("same");

        // TODO: instead of fixed Vlines for cuts, just make a histogram with the cut
        // applied and set its fill color to show its selected region
        i += 1;
    }

    FSHistogram::dumpHistogramCache();

    // c1->cd(1);
    // CUTS = DEFAULT_CUTS;
    // CUTS.ReplaceAll(",unusedE", "");
    // TH1F *hEnUnusedSh = FSModeHistogram::getTH1F(FND_DATA, NT, "pi0pi0pippim", "EnUnusedSh", "(100,0.0,1.0)", Form("CUT(%s)", CUTS.Data()));
    // hEnUnusedSh->SetXTitle("E_{unused}  [GeV/c^{2}]");
    // hEnUnusedSh->SetYTitle("Events");
    // hEnUnusedSh->Draw();
    // TH1F *hEnUnusedShMC = FSModeHistogram::getTH1F(FND_MC, NT, "pi0pi0pippim", "EnUnusedSh", "(100,0.0,1.0)", Form("CUT(%s)", CUTS.Data()));
    // hEnUnusedShMC->Scale(hEnUnusedSh->GetMaximum() / hEnUnusedShMC->GetMaximum());
    // hEnUnusedShMC->SetMarkerColor(kRed);
    // hEnUnusedShMC->Draw("same");
    // TLine *cutUnusedE = new TLine(0.1, 0, 0.1, hEnUnusedSh->GetMaximum());
    // cutUnusedE->SetLineColor(kRed);
    // cutUnusedE->Draw("same");

    // c1->cd(2);
    // CUTS = DEFAULT_CUTS;
    // CUTS.ReplaceAll(",z", "");
    // TH1F *hProdVz = FSModeHistogram::getTH1F(FND_DATA, NT, "pi0pi0pippim", "ProdVz", "(100,0.,100.0)", Form("CUT(%s)", CUTS.Data()));
    // hProdVz->SetXTitle("ProdVz  [GeV/c^{2}]");
    // hProdVz->SetYTitle("Events");
    // hProdVz->Draw();
    // TH1F *hProdVzMC = FSModeHistogram::getTH1F(FND_MC, NT, "pi0pi0pippim", "ProdVz", "(100,0.,100.0)", Form("CUT(%s)", CUTS.Data()));
    // hProdVzMC->Scale(hProdVz->GetMaximum() / hProdVzMC->GetMaximum());
    // hProdVzMC->SetMarkerColor(kRed);
    // hProdVzMC->Draw("same");
    // TLine *cutVz_low = new TLine(52, 0, 52, hProdVz->GetMaximum());
    // cutVz_low->SetLineColor(kRed);
    // cutVz_low->Draw("same");
    // TLine *cutVz_hi = new TLine(78, 0, 78, hProdVz->GetMaximum());
    // cutVz_hi->SetLineColor(kRed);
    // cutVz_hi->Draw("same");

    // c1->cd(3);
    // CUTS = DEFAULT_CUTS;
    // CUTS.ReplaceAll(",t", "");
    // TH1F *htk = FSModeHistogram::getTH1F(FND_DATA, NT, "pi0pi0pippim", "abs(-1*MASS2([proton],-GLUEXTARGET))", "(100,0,5)", Form("CUT(%s)", CUTS.Data()));
    // htk->SetXTitle("|-t| [GeV^{2}]");
    // htk->SetYTitle("Entries");
    // htk->Draw();
    // TH1F *htkMC = FSModeHistogram::getTH1F(FND_MC, NT, "pi0pi0pippim", "abs(-1*MASS2([proton],-GLUEXTARGET))", "(100,0,5)", Form("CUT(%s)", CUTS.Data()));
    // htkMC->Scale(htk->GetMaximum() / htkMC->GetMaximum());
    // htkMC->SetMarkerColor(kRed);
    // htkMC->Draw("same");

    // c1->cd(4);
    // CUTS = DEFAULT_CUTS;
    // CUTS.ReplaceAll(",eBeam", "");
    // TH1F *hEnPB = FSModeHistogram::getTH1F(FND_DATA, NT, "pi0pi0pippim", "EnPB", "(125,5,12)", Form("CUT(%s)", CUTS.Data()));
    // hEnPB->SetXTitle("E_{beam} [GeV]");
    // hEnPB->SetYTitle("Entries");
    // hEnPB->Draw();
    // TH1F *hEnPBMC = FSModeHistogram::getTH1F(FND_MC, NT, "pi0pi0pippim", "EnPB", "(125,5,12)", Form("CUT(%s)", CUTS.Data()));
    // hEnPBMC->Scale(hEnPB->GetMaximum() / hEnPBMC->GetMaximum());
    // hEnPBMC->SetMarkerColor(kRed);
    // hEnPBMC->Draw("same");

    // c1->cd(5);
    // CUTS = DEFAULT_CUTS;
    // CUTS.ReplaceAll(",MM2", "");
    // TH1F *hMM2 = FSModeHistogram::getTH1F(FND_DATA, NT, "pi0pi0pippim", "RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5)", "(100,-0.1,0.1)", Form("CUT(%s)", CUTS.Data()));
    // hMM2->SetXTitle("Missing Mass Squared [GeV^{2}]");
    // hMM2->SetYTitle("Entries");
    // hMM2->Draw();
    // TH1F *hMM2MC = FSModeHistogram::getTH1F(FND_MC, NT, "pi0pi0pippim", "RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5)", "(100,-0.1,0.1)", Form("CUT(%s)", CUTS.Data()));
    // hMM2MC->Scale(hMM2->GetMaximum() / hMM2MC->GetMaximum());
    // hMM2MC->SetMarkerColor(kRed);
    // hMM2MC->Draw("same");

    // c1->cd(6);
    // CUTS = DEFAULT_CUTS;
    // CUTS.ReplaceAll(",chi2,", ",");
    // TH1F *hChi2DOF = FSModeHistogram::getTH1F(FND_DATA, NT, "pi0pi0pippim", "Chi2DOF", "(40,0,20)", Form("CUT(%s)", CUTS.Data()));
    // hChi2DOF->SetXTitle("#chi^{2}/dof");
    // hChi2DOF->SetYTitle("Events");
    // hChi2DOF->Draw();
    // TH1F *hChi2DOFMC = FSModeHistogram::getTH1F(FND_MC, NT, "pi0pi0pippim", "Chi2DOF", "(40,0,20)", Form("CUT(%s)", CUTS.Data()));
    // hChi2DOFMC->Scale(hChi2DOF->GetMaximum() / hChi2DOFMC->GetMaximum());
    // hChi2DOFMC->SetMarkerColor(kRed);
    // hChi2DOFMC->Draw("same");

    return;
}