# Summary
This directory contains the event selection scripts that gives us the flat trees for our AmpTools PWA. Note that **the DSelector (and runSelector) is deprecated**, and to refer to the FSRoot based scripts for the cuts in use. The old selectors are kept in the [DSelector folder](./DSelector/) for reference to past analyses. All scripts are intended to be run via an interactive ROOT session rather than compiled. Note that due to the large statistics of the channel, it is best to run many of the scripts on HPC systems.

# Primary Analysis Steps
1. This analysis uses the *hybrid method*, and so we can begin by only selecting the best $\chi^2 / NDF$ for each beam photon in each event. This will massively reduce our file sizes and simplify the process for the rest of the analysis. This is done in [`skim_best_chi2.cc`](./skim_best_chi2.cc), along with a very loose cut on $\chi^2/ndf$ and the selection of the coherent peak $8.2<E_\gamma<8.8$ GeV region, to reduce the file size further.
2. Then we want to see how our general cuts on the physics of interest will affect our distributions. Run `plot_cuts.cc` to view these effects for each variable being cut, for data and signal MC separately.
    a. Use the `load_broad_cuts` function to load and apply the predefined broad cuts to your dataset. This is separated so as to remain consistent between plotting scripts, and avoid mistakes of one cut being edited but not the other.
3. Once we are happy with the broad cuts, we can move on to selecting the $\omega$ signal from our final state. To do this we will produce new friend trees of FILL HERE LATER
4. Once you are happy and tuned the cuts to your preference, the trees can be saved with `finalize_amptools_trees.cc`.

# Specialized Scripts
Almost every script utilizes to establish common tree names, categories, mode info, etc. Below are scripts for specific studies done to motivate the choices made in the primary analysis steps.

## Convenient Functions
* `fsroot_setup.cc` is used by nearly every script to load common parameters and setup the session for FSRoot
* `load_broad_cuts.cc` defines the FSCuts on RF, chi2, and all other non-sideband cuts. By defining them in one function, every script is ensured to use the same cut regions.

## Sideband Studies
* `omega_fits.cc` fits the $3\pi$ spectrum with a voigt profile for the signal, and some polynomial for the background in data and MC. The fits are done in bins of $M_{\omega\pi^0}$ and $-t$ and plotted individually, along with the resulting gaussian width parameter from the Voigt function.
* `integrate_voigt.cc` reports the boundaries that provide a 95% coverage of the area under the voigt profile used in `omega_fits.cc`. This result gives us the signal region for our sideband subtraction.
* `sideband_variation.cc` studies what effects the sideband locations has on the final $\omega\pi^0$ mass spectra.
* `compare_sideband_methods.cc` studies how the $\omega\pi^0$ distribution (after broad cuts have been applied) is affected by the 2 options of sideband subtraction below
  * Subtracting via the very complicated 2D sideband method with a variety of weights for different regions
  * Simply sideband subtracting each individual $\pi^0_i\pi^+\pi^-$ combo independently.

## Other
* `friend_reorder.cc` creates convenient friend trees that establish signal and sideband regions for every event, such that loading the friend tree directly allows one to select them by choosing `signal==1` or `sideband==1`. The friend trees are split by permutation, so both can be added to the same histogram as intended by the norwegian method.
* `report_integrals.cc` will report the number of combos leftover after each cut
* `van_hove_analysis.cc` studies where the baryon systems show up in the different momentum regions, and how cuts on the momenta can influence the angular distributions