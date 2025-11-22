/**
 * @file fsroot_setup.cc
 * @author Kevin Scheuer
 * @brief Setup FSRoot modes and categories for neutral b1 analysis
 * 
 */

#include <utility>

#include "FSMode/FSModeCollection.h"
#include "FSBasic/FSHistogram.h"


/**
 * @brief FSRoot setup function to define modes and categories
 * 
 * @param read_cache whether to read histogram cache 
 */
std::pair<TString, TString> setup(bool read_cache = false)
{    
    TString CATEGORY("pi0pi0pippim");

    // make sure modes & categories aren't defined already
    if (FSModeCollection::modeVector().size() != 0)
        std::cerr << "Warning: FSModeCollection already has modes defined!" << std::endl;

    if (read_cache)
        FSHistogram::readHistogramCache();
        
    // https://github.com/JeffersonLab/hd_utilities/blob/master/FlattenForFSRoot/Documentation/GlueXFSRootFormat.pdf
    // unused numbers in front dropped, so this reads as 1 proton,
    // then 1 pi+, 1 pi-, 2 pi0s
    FSModeCollection::addModeInfo("100_112")->addCategory(CATEGORY);
    TString NT = "ntFSGlueX_100_112"; // TODO: not yet confirmed if this breaks modecode histograms

    return std::make_pair(NT, CATEGORY);
}