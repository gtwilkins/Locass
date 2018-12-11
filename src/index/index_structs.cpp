/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "index_structs.h"

void CharCount::clear()
{
    counts[0] = counts[1] = counts[2] = counts[3] = 0;
    endCounts = 0;
}

int CharCount::getBranchCount()
{
    return bool(counts[0]) + bool(counts[1]) + bool(counts[2]) + bool(counts[3]);
}

int CharCount::getMaxBranch()
{
    int best = counts[0] > 0 ? 0 : ( counts[1] > 0 ? 1 : ( counts[2] > 0 ? 2 : ( counts[3] > 0 ? 3 : 4 ) ) );
    bool bad = false;
    for ( int i = best + 1; i < 4; i++ )
    {
        if ( !counts[i] ) continue;
        if ( counts[i] == counts[best] ) bad = true;
        else if ( counts[i] > counts[best] )
        {
            best = i;
            bad = false;
        }
    }
    return bad ? 4 : best;
}

void CharCount::operator -=( CharCount &rhs )
{
    for ( int i = 0; i < 4; i++ ) counts[i] -= rhs.counts[i];
    endCounts -= rhs.endCounts;
}