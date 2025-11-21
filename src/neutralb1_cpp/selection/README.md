# Summary
This directory contains the event selection scripts that gives us the flat trees for our AmpTools PWA. Note that **the DSelector (and runSelector) is deprecated**, and to refer to the FSRoot based scripts for the cuts in use. The old selectors are included for reference to past analyses. All scripts are intended to be run via an interactive ROOT session rather than compiled, to keep to the intent of FSRoot being interactive 

# Analysis Steps
1. This analysis uses the *hybrid method*, and so we can begin by only selecting the best $\chi^2 / NDF$ for each beam photon in each event. This will massively reduce our file sizes and simplify the process for the rest of the analysis. This is done in [`skim_best_chi2.cc`](./skim_best_chi2.cc), along with a very loose cut on $\chi^2/ndf$ and the selection of the coherent peak $8.2<E_\gamma<8.8$ GeV region, to reduce the file size further.
2. Then we want to see how our general cuts on the physics of interest will affect our distributions. Run `plot_cuts.cc` to view these effects for each variable being cut, for data and signal MC separately.
    a. Use the `load_broad_cuts` function to load and apply the predefined broad cuts to your dataset. This is separated so as to remain consistent between plotting scripts, and avoid mistakes of one cut being edited but not the other.
    b. These broad cuts only select the final $\pi^0\pi^0\pi^+\pi^-p'$ final state.
4. Once we are happy with the broad cuts, we can move on to selecting the $\omega$ signal from our final state. To do this we will produce new friend trees of FILL HERE LATER
5. Now we can study all the combinatorics, and view the signal and background events with `plot_combos.cc`.
6. Once you are happy and tuned the cuts to your preference, the trees can be saved with `finalize_amptools_trees.cc`.

## Specialized Scripts
Almost every script utilizes `fsroot_setup.cc` to establish common tree names, categories, mode info, etc. Below are scripts for specific studies done to motivate the choices made in the analysis steps.

### Sideband Studies
* `omega_fits.cc` fits the $3\pi$ spectrum with a voigt profile for the signal, and some polynomial for the background in data and MC. The fits are done in bins of $M_{\omega\pi^0}$ and $-t$ and plotted individually, along with the resulting gaussian width parameter from the Voigt function.
* `integrate_voigt.cc` reports the boundaries that provide a 95% coverage of the area under the voigt profile used in `omega_fits.cc`. This result gives us the signal region for our sideband subtraction.
* `sideband_variation.cc` studies what effects the sideband locations has on the final $\omega\pi^0$ mass spectra.
* `compare_sideband_methods.cc` studies how the $\omega\pi^0$ distribution (after broad cuts have been applied) is affected by the 2 options of sideband subtraction below
  * Subtracting via the very complicated 2D sideband method with a variety of weights for different regions
  * Simply sideband subtracting each individual $\pi^0_i\pi^+\pi^-$ combo independently.

### Other
* `report_integrals.cc` will report the number of combos leftover after each cut
* `van_hove_analysis.cc` studies where the baryon systems show up in the different momentum regions, and how cuts on the momenta can influence the angular distributions