/**
 * @file load_broad_branch_cuts.cc
 * @author Kevin Scheuer
 * @brief Like broad cuts, but use the branches created in the friend reordering
 *
  */

#include "TColor.h"

#include "FSBasic/FSCut.h"

std::map<TString, Int_t> load_broad_branch_cuts()
{
    std::map<TString, TString> cut_map;    

    // These are the same as the broad cuts, but using the friend tree branch names. We
    // have to do this so that we can make the same selective cuts for plotting.
    cut_map["unusedE"] = "unusedE<0.1";                                                              // energy of unmatched showers in time with combo but not part of the combo
    cut_map["unusedTracks"] = "unusedTracks<1";                                                      // number charged particle tracks not used in combo
    cut_map["z"] = "z>=51.2&&z<=78.8";                                                               // avoid interactions with walls of target
    cut_map["MM2"] = "abs(MM2)<0.05";                                                                // remove events with large missing mass2
    cut_map["chi2"] = "chi2<5";                                                                      // kinematic fit quality cut
    cut_map["t"] = "abs(t)<1.0";                                                                     // define the momentum transfer with the proton at the lower vertex
    cut_map["shQuality"] = "shQualityP2a>0.5&&shQualityP2b>0.5&&shQualityP3a>0.5&&shQualityP3b>0.5"; // classifier probability that showers are from neutral particles
    cut_map["pzPi0"] = "MOMENTUMZBOOST(2;B,GLUEXTARGET)>-0.1";                                     // require pi0 moving forward in CM frame

    TString cuts;
    for (std::map<TString, TString>::const_iterator it = cut_map.begin(); it != cut_map.end(); ++it)
    {
        FSCut::defineCut(it->first, it->second);
        cuts += (cuts.IsNull() ? "" : ",") + it->first;
    }

    FSCut::defineCut("rf", "abs(rf)<2.0", "abs(rf)>2.0", 0.125);

    FSCut::defineCut("cut", "cut==1"); // This is the equivalent of all the above cuts applied together
    FSCut::defineCut("signal", "signal==1");
    FSCut::defineCut("sideband", "sideband==1");

    // not strictly necessary, but consistent return type with broad cuts loader
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