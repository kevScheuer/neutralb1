:bangbang: The selector is not yet setup. That directory will host the source files, but for now the files are stored [here](./batch_scripts/). :bangbang:

We've now made the selections on our data and have ROOT trees ready for partial wave analysis. We will typically want to be submitting "mass-independent" fits (fits performed in individual bins of $\omega\pi^0$ mass), but writing and submitting all those fits is a pain. Instead we can use an automated [submitter](./batch_scripts/submit.py) that can handle all that for us. 

# Architecture
1. The user requests via [submit.py](./batch_scripts/submit.py) fits to MC or data with a particular waveset in bins of mass and t (or to the entire mass range if doing a mass dependent fit). This is done through the use of several command line arguments, and is the only file the user needs to interact with.
2. [write_config.py](./batch_scripts/write_config.py) then creates an AmpTools fit configuration file based off the user requested parameters. 
3. Each bin to be fit is given a temporary directory path within `/volatile/halld/home/{USER}/TMPDIR` and all the necessary data files / config files are sent or linked to the directories.
4. A job is submitted in each bin, which executes `run_fit.sh` with a set of parameters specific to that bin. 
5. The `.fit` files and diagnostic plots are copied to the same directory path but now within `/volatile/halld/home/{USER}/ampToolsFits/` where they are now ready to be aggregated for [analysis](../analysis/)

# Monte Carlo files
You'll also notice the [load_monte_carlo](./load_monte_carlo) directory. This holds the scripts necessary to load pre-made Monte Carlo off the cache directory (created by Justin) and apply the DSelector to those files. `run_mc.py` is the main interactive tool to load these files.