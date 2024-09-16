# Introduction

 :bangbang:**This repository must be on a JLab ifarm node**:bangbang:

Within this repository lies all scripts needed to perform a full partial wave analysis of the neutral $\omega\pi^0$ channel from start to finish, which is dominated by a neutral $b_1$ resonance for which this repo is named after. The repo is structured such that one can:
1. Begin in [selector](./selector/) to load GlueX ROOT data files and select events of interest to prepare them for a partial wave analysis.
2. Submit partial wave fits to the selected ROOT trees in [submission](./submission/) and obtain `.fit` result files in bins of mass and $-t$.
3. Analyze those scripts in [analysis](./analysis/) using jupyter notebooks and an [all-in-one pwa plotter](./analysis/scripts/pwa_tools.py). 

The scripts are tuned to $\omega \pi^0$ analysis, though are constantly being generalized so that others may easily copy its functionality.

# Setup
To run the scripts, you must have a working version of `conda` in order to source and start the python environment here. [Miniforge](https://github.com/conda-forge/miniforge#mambaforge) and [Miniconda](https://docs.anaconda.com/free/miniconda/index.html) are recommended. 

To create the environment and activate it, run:
```
conda env create --file environment.yml
conda activate neutralb1
```

:warning: If perl is installed with the base conda environment, you may get an error when attempting to source the setup scripts in the [batch_scripts folder](./submission/batch_scripts) due to the setup script trying to use the ifarm's default perl binary. Simply remove perl from conda and the error should be resolved.

## ROOT setup
You will need to load the load the AmpTools libraries on ROOT startup for many of these scripts. If you plan to use the [selector](./selector/) scripts these depend on FSROOT and will also require to be loaded on startup This is done with 2 files placed in your home directory. The first is a `.rootrc` file:
```
Unix.*.Root.DynamicPath: .:$(FSROOT):$(ROOTSYS)/lib:
Unix.*.Root.MacroPath:   .:$(FSROOT):$(AMPTOOLS):$(AMPTOOLS)/ROOT:
Rint.Logon:     $(FSROOT)/rootlogon.FSROOT.C
Rint.Logoff:    $(FSROOT)/rootlogoff.FSROOT.C

Proof.MaxOldSessions 500
ProofLite.Sandbox /work/halld/kscheuer/proof/
```

The second is the `rootlogon.C` script mentioned above:
```
{
    printf("Load Classes");
    // below line needs macropath edited in .rootrc file
    gROOT->ProcessLine(".x loadAmpTools.C");

    if(getenv("ROOT_ANALYSIS_HOME"))       
        gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/programs/MakePROOFPackage/SETUP.C");
}
```