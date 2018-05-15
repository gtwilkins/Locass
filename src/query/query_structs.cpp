/*
 * Copyright (C) 2017 Glen T. Wilkins <glen.t.wilkins@gmail.com>
 * Written by Glen T. Wilkins
 * 
 * This file is part of the Locass software package <https://github.com/gtwilkins/Locass>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "query_structs.h"
#include <algorithm>

void MappedSeqs::setBest( string &seq )
{
    int bestLen = 0;
    int bestErrors = 1;
    ReadStruct bestRead;
    
    for ( ReadStruct &read : reads )
    {
        int errors = 1;
        int len = min( read.coords[1], (int)seq.length() ) - max( 0, read.coords[0] );
        int diff = -read.coords[0];
        for ( int i( max( 0, -diff ) ); i < len; i++ )
        {
            errors += seq[i] != read.seq[i + diff];
        }
        
        if ( (float)len / (float)errors > (float)bestLen / (float)bestErrors )
        {
            bestLen = len;
            bestErrors = errors;
            bestRead = read;
        }
    }
    
    string bestSeq = bestRead.seq.substr( max( 0, -bestRead.coords[0] ), seq.length() - max( 0, bestRead.coords[0] ) );
    int limits[2] = { max( -bestRead.coords[0], 0 ), min( (int)bestRead.seq.length(), bestRead.coords[1] ) - bestRead.coords[0] };
    int32_t coords[2] = { max( 0, bestRead.coords[0] ), max( 0, bestRead.coords[0] ) + (int)bestSeq.length() };
    vector<ReadStruct> newReads, falseReads;
    for ( ReadStruct &read : reads )
    {
        int diff = read.coords[0] - coords[0];
        int i = max( coords[0], max( 0, diff ) );
        int j = min( coords[1], (int)read.seq.length() + diff );
        while ( bestSeq[i] == read.seq[i - diff] && i < j )
        {
            i++;
        }
        if ( i == j )
        {
            newReads.push_back( read );
        }
        else
        {
            falseReads.push_back( read );
        }
    }
    
    reads = newReads;
}

void MappedSeqs::sort()
{
    std::sort( reads.begin(), reads.end(), []( ReadStruct &a, ReadStruct &b ){ 
        return a.coords[0] == b.coords[0] ? a.coords[1] > b.coords[1] : a.coords[0] < b.coords[0];
    });
}

void MappedSeqs::updateTethers( string &seq )
{
    for ( ReadStruct &read : reads )
    {
        while ( read.tether[0] > read.coords[0] 
                && seq[read.tether[0]-1] == read.seq[ read.tether[0]-1 - read.coords[0] ] )
        {
            read.tether[0]--;
        }
        while ( read.tether[1] < seq.length() 
                && read.tether[1] < read.coords[1]
                && read.seq[ read.tether[1] - read.coords[0] ] == seq[ read.tether[1] ] )
        {
            read.tether[1]++;
        }
    }
}