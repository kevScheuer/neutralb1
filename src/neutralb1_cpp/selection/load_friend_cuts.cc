/**
 * @file load_friend_cuts.cc
 * @author Kevin Scheuer
 * @brief Load FSRoot cuts for implementing the omega sideband subtraction
 * 
 * Uses the same branch definitions as in friend_reorder.cc to define
 * the signal and sideband cuts for omega selection.
 */

#include "FSBasic/FSCut.h"

const int load_friend_cuts()
{
    FSTree::addFriendTree("cut");
    FSTree::addFriendTree("signal");
    FSTree::addFriendTree("sideband");
    FSTree::addFriendTree("t");
    FSTree::addFriendTree("M4Pi");
    FSTree::addFriendTree("MRecoilPi");

    FSCut::defineCut("cut", "cut==1");
    FSCut::defineCut("signal", "signal==1");
    FSCut::defineCut("sideband", "sideband==1");

    const int n_permutations = 2;
    return n_permutations; // return number of pi0 permutations
}