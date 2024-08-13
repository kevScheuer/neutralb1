This directory contains all the scripts needed to convert the reaction filtered ROOT trees into "Flat" trees, denoted as AmpToolsInputTrees in [the dataa source files directory](../source_files/data/), that are ready for amplitude analysis. 

# Background
To summarize where the ROOT trees come from:
1. The detector takes raw data during a run period
2. The Reconstructed Event Storage (REST) processes the raw data into event-level entries by reconstructing the particle tracks from all the detector information
3. Users request a particular reaction, such as $\gamma p \rightarrow \omega \pi^0 p$ from the REST dataset, and a "reaction filter" is applied. Often the "exclusive" reaction is used, meaning all requested particles are present
   1. The filter analyzes all the possible combinations in the REST data to apply selection to. For $\omega\pi^0$ we need a tagged beam photon, two positively charged tracks ($p$ & $\pi^+$) and one negative track ($\pi^-$) with four neutral showers are detected. 
   2. Common reactions with certain configurations are bundled into "analysis launches", which have loose cuts applied to the data. These cuts help reduce file sizes, while still being loose enough to keep almost all true signal events
4. A kinematic fit is applied, and converged fits are stored in a ROOT tree with an associated $\chi^2$ value reported for each event. These ROOT trees are then analyzed here, and further 

# FSROOT

## Particle ID\# Mapping
For $\gamma p \rightarrow \omega \pi^0 p$:

| Particle Index | Particle (parent) |
|:--------------:|:----------------- |
| 0              | beam photon       |
| 1              | proton            |
| 2              | $\pi^+$ ($\omega$)|
| 3              | $\pi^-$ ($\omega$)|
| 4              | $\pi^0$ ($\omega$)|
| 5              | $\pi^0$ (bachelor)|