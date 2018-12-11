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

#ifndef QUERY_STATE_H
#define QUERY_STATE_H

#include "types.h"
#include "parameters.h"
#include "query_structs.h"

struct QueryState
{
    QueryState( uint8_t* query, int seqLength, int minOverlap ) : q( query ), seqLen( seqLength ), minOver( minOverlap ){};
    int setFirst();
    void updateSeq( string &seq, int off, bool drxn );
    uint8_t* q;
    int seqLen, minOver;
    vector<uint8_t> endOverlaps;
    vector<CharId> endRanks, endCounts;
};

struct QueryCorrectState
{
    QueryCorrectState( uint8_t* query, int corrLength, int seqLength, int minOverlap );
    uint8_t* q;
    vector<uint8_t> endOverlaps, alts, altIts;
    vector<uint32_t> ends;
    vector<CharId> endRanks, endCounts, altRanks, altCounts;
    int corrLen, seqLen, minOver;
    uint32_t endCount;
    uint8_t endCutoff, altCutoff;
    bool fresh;
};

struct QueryKmerState
{
    QueryKmerState( uint8_t* query ) : q( query ) { errors = cleans = 0; };
    
    uint8_t* q;
    int errors, cleans;
    CharId rank, edge, count;
};

struct QuerySeedState
{
    QuerySeedState( string &seq, MappedSeqs &ms, int errorRate );
    void add( ReadStruct &read, bool drxn );
    bool advance( int it, int k, int &thisErrors, bool drxn );
    void countChunks( int &count, bool drxn );
    bool doAdd( int errors, bool drxn );
    void setup( int it, bool drxn );
    
    uint8_t* query[2];
    MappedSeqs &ms;
    int i, j, iBegin, iStart, iAdd, iError, len, minLen;
    int maxErrors, errorsPer, errorsLeft;
    int leftCount, rightCount;
    bool exact = true;
    CharId rank, count;
    vector<int32_t> chunks[2];
};



#endif /* QUERY_STATE_H */

