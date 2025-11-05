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

TString load_broad_cuts()
{
    TString cuts;

    std::map<TString, TString> cut_map;

    // All values in map below must use FSRoot defined branch names in the trees
    cut_map["unusedE"] = "EnUnusedSh<0.1";                             // energy of unmatched showers in time with combo but not part of the combo
    cut_map["unusedTracks"] = "NumUnusedTracks<1";                     // number charged particle tracks not used in combo
    cut_map["z"] = "ProdVz>=51.2&&ProdVz<=78.8";                       // 
    cut_map["MM2"] = "abs(RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5))<0.05"; // remove events with large missing mass2
    cut_map["eBeam"] = "(EnPB>8.2&&EnPB<8.8)";                         // coherent peak region for PWA
    cut_map["chi2"] = "Chi2DOF<5";                                     // kinematic fit quality cut
    // cut_map["chi2rank"]       = "Chi2Rank==1";
    cut_map["t"] = "abs(-1*MASS2([proton],-GLUEXTARGET))<1.0";                                       // define the momentum transfer with the proton at the lower vertex
    cut_map["shQuality"] = "ShQualityP4a>0.5&&ShQualityP4b>0.5&&ShQualityP5a>0.5&&ShQualityP5b>0.5"; // Ask Justin

    for (std::map<TString, TString>::const_iterator it = cut_map.begin(); it != cut_map.end(); ++it)
    {
        FSCut::defineCut(it->first, it->second);
    }

    for (std::map<TString, TString>::const_iterator it = cut_map.begin(); it != cut_map.end(); ++it)
    {
        cuts += it->first;
        if (std::next(it) != cut_map.end())
        {
            cuts += ",";
        }
    }

    return cuts;
}