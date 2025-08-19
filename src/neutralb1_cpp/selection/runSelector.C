// macro to process analysis TTree with TSelector
#include <iostream> 

#include "DSelector/DPROOFLiteManager.h"
#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TTree.h"

R__LOAD_LIBRARY(libDSelector.so)

void runSelector(TString runNumber = "30496", TString myPath = "/sciclone/gluex10/RunPeriod-2017-01/analysis/ver35/tree_pi0pi0pippim__B4/merged/", TString myOption = "") 

{
  // Load DSelector library
  gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
  int Proof_Nthreads = 8;

  // process signal 
  TString sampleDir = myPath;
  std::cout<<"running selector on files in: "<<sampleDir.Data() << "\n";
  
  TString treeName = "pi0pi0pippim__B4_Tree";
  if(myPath.Contains("tree_thrown")) treeName = "Thrown_Tree";
  TChain *chain = new TChain(treeName); 
  TSystemDirectory dir(sampleDir, sampleDir);
  TList *files = dir.GetListOfFiles();
  int ifile = 0;
  if(files) {
	  TSystemFile *file;
	  TString fileName;
	  TIter next(files);
	  
	  // loop over files
	  while ((file=(TSystemFile*)next())) {
		  fileName = file->GetName();
		  if(fileName.Contains(runNumber)) {
			  std::cout<<fileName.Data()<< "\n";
			  
			  // check if file corrupted
			  TFile f(sampleDir+fileName);
			  if(f.TestBit(TFile::kRecovered)) {
				  std::cout<<"file corrupted -> skipping"<< "\n";
				  continue;
			  }
			  if(f.IsZombie()) {
				  std::cout<<"file is a Zombie -> skipping"<< "\n";
				  continue;
			  }
			  
			  // add file to chain
			  chain->Add(sampleDir+fileName);
			  ifile++;
		  }
	  }

	  std::cout<<"total entries in TChain = "<<chain->GetEntries()<<" from "<<ifile<<" files"<<endl;
	  DPROOFLiteManager::Process_Chain(chain, "DSelector_omegapi_all.C++", Proof_Nthreads, Form("hist_pomegapi_omega3pi_%s.acc.root", runNumber.Data()), Form("tree_pomegapi_omega3pi_%s.acc.root", runNumber.Data()), Form("%s DEFAULTFLATOFF",myOption.Data()));
  }

  return;
}
