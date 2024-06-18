# Introduction

:warning: **This repository must be on a JLab ifarm node**:warning:

Within this repository lies all scripts needed to submit amplitude analysis fits to the ifarm, aggregate them into csv files, and analyze those csv's using a python framework. The scripts are tuned to $\omega \pi^0$ analysis, though are constantly being generalized so that others may easily copy its functionality. To run the scripts, you must have a working version of `conda` in order to source and start the python environment here. [Miniforge](https://github.com/conda-forge/miniforge#mambaforge) and [Miniconda](https://docs.anaconda.com/free/miniconda/index.html) are recommended. 

To create the environment and activate it, run:
```
conda env create --file environment.yml
conda activate neutralb1
```

You will also need to load the load the AmpTools libraries on ROOT startup. This is done with 2 files placed in your home directory. The first is a `.rootrc` file:
```
Proof.MaxOldSessions 500
ProofLite.Sandbox /YOUR_WORK_DIRECTORY/proof/
Unix.*.Root.MacroPath:   .:$(AMPTOOLS):$(AMPTOOLS)/ROOT:
Rint.Logon: ~/rootlogon.C
```

The second is the `rootlogon.C` script mentioned above:
```
{
  printf("Load AmpTools Classes");
  // below line needs macropath edited in .rootrc file
  gROOT->ProcessLine(".x loadAmpTools.C");

  if(getenv("ROOT_ANALYSIS_HOME"))       
  	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/programs/MakePROOFPackage/SETUP.C");
}
```