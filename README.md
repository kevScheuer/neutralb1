# Introduction

 :bangbang:**This repository must be on a JLab ifarm node**:bangbang:

Within this repository lies all scripts needed to perform a full partial wave analysis of the neutral $\omega\pi^0$ channel from start to finish, which is dominated by a neutral $b_1$ resonance for which this repo is named after. 

# Setup
To run the scripts, you must have a working version of `conda` in order to source and start the python environment here. [Miniforge](https://github.com/conda-forge/miniforge#mambaforge) and [Miniconda](https://docs.anaconda.com/free/miniconda/index.html) are recommended. 

To create the environment and activate it, run:
```
conda env create --file environment.yml
conda activate neutralb1
```

:warning: If perl is installed with the base conda environment, you may get an error when attempting to source the setup scripts in the [batch_scripts folder](./submission/batch_scripts) due to the setup script trying to use the ifarm's default perl binary. Simply remove perl from conda and the error should be resolved.