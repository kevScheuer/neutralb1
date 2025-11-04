# Summary
This directory contains the event selection scripts that gives us the flat trees for our AmpTools PWA. Note that **the DSelector (and runSelector) is deprecated**, and to refer to the FSRoot based scripts for the cuts in use. The old selectors are included for reference to past analyses. All scripts are intended to be run via an interactive ROOT session rather than compiled, to keep to the intent of FSRoot being interactive 

# Analysis Steps
1. We first want to see how our general cuts on the physics of interest will affect our distributions. Run `plot_broad_cuts.cc` to view these effects for each variable being cut, with comparisons to Monte Carlo.
    a. Use the `load_broad_cuts` function to load and apply the predefined broad cuts to your dataset. This is separated so as to remain consistent between plotting scripts, and avoid mistakes of one cut being edited but not the other.
2. Next we can study all the combinatorics, and view the signal and background events with `plot_combos.cc`.
3. Once you are happy and tuned the cuts to your preference, the trees can be saved with `skim.cc`.

## Other Scripts
* `print_integrals.cc` will report the number of combos leftover after each broad cut