/**
 * @file report_integrals.cc
 * @author Kevin Scheuer
 * @brief Print out the omega pi0 integrals after each cut to see # events remaining  
 */

#include <iostream>
#include <map>


#include "TString.h"

#include "load_broad_cuts.cc"
#include "fsroot_setup.cc"

void report_integrals(
    bool mc = true,
    bool bggen = false)
{
    TString input_data_files = "data_all.root";
    TString input_mc_files = "mc_all.root";

    // Load broad cuts
    std::map<TString, Int_t> cut_color_map = load_broad_cuts();

    // TODO: Iterate through ALL cuts (including RF and SB subtraction), and report 
    // the integrals after each cut is applied.
    
}