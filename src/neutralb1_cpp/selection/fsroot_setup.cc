/**
 * @file fsroot_setup.cc
 * @author Kevin Scheuer
 * @brief Setup FSRoot modes and categories for neutral b1 analysis
 * 
 */

#include "FSMode/FSModeCollection.h"
#include "FSBasic/FSHistogram.h"

void setup(bool read_cache = false, bool dump_cache = false)
{
    TString NT("ntFSGlueX_MODECODE");
    TString CATEGORY("pi0pi0pippim");

    // make sure modes & categories aren't defined already
    if (FSModeCollection::modeVector().size() != 0)
        return;

    if (read_cache)
        FSHistogram::readHistogramCache();
        
    // https://github.com/JeffersonLab/hd_utilities/blob/master/FlattenForFSRoot/Documentation/GlueXFSRootFormat.pdf
    // unused numbers in front dropped, so this reads as 1 proton,
    // then 1 pi+, 1 pi-, 2 pi0s
    FSModeCollection::addModeInfo("100_112")->addCategory(CATEGORY);
}