/**
 * @file load_broad_cuts.cc
 * @author Kevin Scheuer
 * @brief Defines function to load broad cuts for pi0pi0pi+pi- final state selection
 *
 * To keep cuts consistent across different selection and plotting scripts, this file
 * defines a function to load the broad cuts used in the pi0pi0pi+pi- final state
 * selection.
 */

#include <map>
#include <iostream>

#include "TString.h"
#include "FSBasic/FSCut.h"
#include "TColor.h"

std::map<TString, Int_t> load_broad_cuts()
{
    TString cuts;

    std::map<TString, TString> cut_map;

    // All values in map below must use FSRoot defined branch names in the trees
    cut_map["unusedE"] = "EnUnusedSh<0.1";                                                           // energy of unmatched showers in time with combo but not part of the combo
    cut_map["unusedTracks"] = "NumUnusedTracks<1";                                                   // number charged particle tracks not used in combo
    cut_map["z"] = "ProdVz>=51.2&&ProdVz<=78.8";                                                     //
    cut_map["MM2"] = "abs(RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5))<0.05";                               // remove events with large missing mass2
    cut_map["chi2"] = "Chi2DOF<5";                                                                   // kinematic fit quality cut
    cut_map["t"] = "abs(-1*MASS2([proton],-GLUEXTARGET))<1.0";                                       // define the momentum transfer with the proton at the lower vertex
    cut_map["shQuality"] = "ShQualityP4a>0.5&&ShQualityP4b>0.5&&ShQualityP5a>0.5&&ShQualityP5b>0.5"; // Ask Justin

    for (std::map<TString, TString>::const_iterator it = cut_map.begin(); it != cut_map.end(); ++it)
    {
        FSCut::defineCut(it->first, it->second);
    }

    // Add the accidental subtraction sideband cut. Use 1/8 weight since there are 8
    // out of time RF buckets in the data.
    // Don't add to returned cuts TString since it needs separate calls for signal and
    // sideband regions
    FSCut::defineCut("rf", "OR(abs(RFDeltaT)<2.0)", "abs(RFDeltaT)>2.0", 0.125);

    // Define a cut for just the signal region, this will be the "original data" to
    // compare successive cuts to
    FSCut::defineCut("rfSignal", "(abs(RFDeltaT)<2.0)");

    for (std::map<TString, TString>::const_iterator it = cut_map.begin(); it != cut_map.end(); ++it)
    {
        cuts += it->first;
        if (std::next(it) != cut_map.end())
        {
            cuts += ",";
        }
    }

    // setup a consistent mapping between each cut and a TColor for plotting
    // mimics the 10 color scheme from
    // https://root.cern.ch/doc/v636/classTColor.html
    std::map<TString, Int_t> cut_color_map;
    cut_color_map["unusedE"] = TColor::GetColor("#92dadd");
    cut_color_map["unusedTracks"] = TColor::GetColor("#b9ac70");
    cut_color_map["z"] = TColor::GetColor("#e76300");
    cut_color_map["MM2"] = TColor::GetColor("#a96b59");
    cut_color_map["chi2"] = TColor::GetColor("#832db6");
    cut_color_map["t"] = TColor::GetColor("#3f90da");
    cut_color_map["shQuality"] = TColor::GetColor("#ffa90e");

    return cut_color_map;
}