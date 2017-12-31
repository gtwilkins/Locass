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

#include "query_state.h"
#include <cassert>

extern Parameters params;

QuerySeedState::QuerySeedState( string &seq, MappedSeqs &ms, int errorRate )
: len( seq.length() ), ms( ms )
{
    maxErrors = ( len * errorRate ) / 100;
    errorsPer = min( maxErrors, ( params.readLen * errorRate ) / 100 );
    minLen = min( len, max( 33, min( int( params.readLen / 2.5 ), len / 2 ) ) );
    
    
    int chunkLen = len / ( maxErrors + 1 );
    int remainder = ( maxErrors + 1 ) - ( len % ( maxErrors + 1 ) );
    rightCount = ( maxErrors + 1 ) / 2;
    leftCount = ( maxErrors + 1 ) - rightCount;
    chunks[0].push_back( 0 );
    
    for ( int k ( 0 ); k < ( maxErrors + 1 ); k++ )
    {
        chunks[0].push_back( chunks[0][k] + chunkLen + ( k >= remainder ) );
    }
    
    for ( int k ( chunks[0].size() ); --k >= 0; )
    {
        chunks[1].push_back( len - chunks[0][k] );
    }
    
    assert( chunks[0].back() == len );
    
    ms.chunks = chunks[0];
}

void QuerySeedState::add( ReadStruct &read, bool drxn )
{
    if ( ms.usedIds.find( read.readId ) == ms.usedIds.end() )
    {
        read.tether[0] = drxn ? len - min( i, iError ) : iBegin;
        read.tether[1] = drxn ? len - iBegin : min( i, iError );
        read.coords[0] = drxn ? len - i : i - read.seq.length();
        read.coords[1] = drxn ? len - i + read.seq.length() : i;
        
        int32_t limits[2] = { 0, len };
        int32_t tether[2] = { len, 0 };
        for ( int32_t &chunk : chunks[0] )
        {
            limits[0] = chunk < read.tether[0] ? chunk : limits[0];
            limits[1] = read.tether[1] < chunk ? chunk : limits[1];
            tether[0] = read.tether[0] <= chunk ? min( tether[0], chunk ) : tether[0];
            tether[1] = chunk <= read.tether[1] ? max( tether[1], chunk ) : tether[1];
        }
        
        if ( tether[!drxn] != read.tether[!drxn] && read.coords[!drxn] != read.tether[!drxn] )
        {
            read.tether[!drxn] = tether[!drxn];
        }
        
        ms.reads.push_back( read );
        ms.usedIds.insert( read.readId );
    }
}

bool QuerySeedState::advance( int it, int k, int &thisErrors, bool drxn )
{
    i = it + 1;
    
    if ( i + 1 >= len )
    {
        if ( !thisErrors ) iError = len;
        j = k;
        return true;
    }
    
    assert( i < len );
    
    if ( query[drxn][i] != k )
    {
        if ( !thisErrors )
        {
            iError = i;
        }
        thisErrors++;
    }
    
    if ( thisErrors == 0 )
    {
        iError = i + 1;
    }
    
    int thisChunks = 0;
    countChunks( thisChunks, drxn );
    int thisMaxErrors = max( 0, min( min( errorsPer, errorsLeft ), ( thisChunks - 1 ) * 2 ) );
    
    j = k;
    
    return thisErrors <= thisMaxErrors && i < len - 1;
}

void QuerySeedState::countChunks( int &count, bool drxn )
{
    bool inRange = false, outRange = false;
    for ( int32_t &chunk : chunks[drxn] )
    {
        inRange = inRange || iBegin <= chunk;
        outRange = inRange && i < chunk;
        count += inRange && !outRange;
    }
}

bool QuerySeedState::doAdd( int errors, bool drxn )
{
    int count = 0;
    countChunks( count, drxn );
    
    return errors <= count;
}

void QuerySeedState::setup( int it, bool drxn )
{
    iBegin = chunks[drxn][it];
    iStart = chunks[drxn][it+1];
    errorsLeft = maxErrors - it;
    i = iBegin + 1;
    iAdd = iBegin + minLen;
    iError = len;
    j = query[drxn][i];
}

