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

#include "index_reader.h"
#include "index_structs.h"
#include <cassert>
#include <string.h>
#include <iostream>

IndexReader::IndexReader( Filenames* fns )
{
    FILE* bin;
    fns->setIndex( bin, bwt, idx );
    CharId binId, bwtId, idxId;
    fseek( bin, 1, SEEK_SET );
    fread( &binId, 8, 1, bin );
    fclose( bin );
    
    fread( &beginBwt, 1, 1, bwt );
    fread( &bwtId, 8, 1, bwt );
    fread( &bwtSize, 8, 1, bwt );
    fread( &charCounts[4], 8, 1, bwt );
    fread( &charCounts, 8, 4, bwt );
    
    fread( &beginIdx, 1, 1, idx );
    fread( &idxId, 8, 1, idx );
    fread( &bwtPerIndex, 4, 1, idx );
    fread( &indexPerMark, 4, 1, idx );
    
    if ( binId != bwtId || binId != idxId )
    {
        cerr << "Error: disagreement among data files. They may be corrupted, incomplete or from different sessions." << endl;
        exit( EXIT_FAILURE );
    }
    
    charRanks[0] = charCounts[4];
    charRanks[1] = charRanks[0] + charCounts[0];
    charRanks[2] = charRanks[1] + charCounts[1];
    charRanks[3] = charRanks[2] + charCounts[2];
    
    // Load index
    sizePerIndex = 37;
    fread( &indexSize, 8, 1, idx );
    fread( &markSize, 8, 1, idx );
    index_ = new uint8_t[indexSize * sizePerIndex];
    marks_ = new ReadId[markSize];
    buff = new uint8_t[bwtPerIndex*2];
    fread( index_, 1, indexSize * sizePerIndex, idx );
    fread( marks_, 4, markSize, idx );
    
    runFlag = 1 << 7;
    runMask = ~runFlag;
    memset( &isBaseRun, false, 256 );
    isBaseRun[255] = true;
    for ( int i ( 0 ); i < 4; i++ )
    {
        isBaseRun[i * 63 + 62] = true;
        memset( &decodeBaseChar[ i * 63 ], i, 63 );
        decodeBaseChar[ 252 + i ] = 4;
        decodeBaseRun[ 252 + i ] = i + 1;
        for ( int j ( 0 ); j < 63; j++ )
        {
            decodeBaseRun[ i * 63 + j ] = j + 1;
        }
    }
    
    CharCount ranks;
    setRank( 0, 0, ranks );
    memcpy( &baseCounts[0][0], &ranks.counts, 32 );
    for ( int i ( 0 ); i < 4; i++ )
    {
        setRank( i, charCounts[i], ranks );
        memcpy( &baseCounts[i + 1][0], &ranks.counts, 32 );
        setRank( i, baseCounts[0][i], ranks );
        memcpy( &midRanks[i][0], &ranks.counts, 32 );
    }
}

IndexReader::~IndexReader()
{
    if ( buff ) delete[] buff;
    if ( index_ ) delete[] index_;
    if ( marks_ ) delete[] marks_;
}

void IndexReader::countRange( uint8_t i, CharId rank, CharId count, CharCount &ranks, CharCount &counts )
{
    setRank( i, rank, ranks );
    setRank( i, rank + count, counts );
    for ( int j ( 0 ); j < 5; j++ )
    {
        counts.counts[j] -= ranks.counts[j];
    }
}

void IndexReader::setBaseAll( uint8_t i, uint8_t j, CharId &rank, CharId &count )
{
    rank = baseCounts[i][j];
    count = baseCounts[ i + 1 ][j] - rank;
}

void IndexReader::setBaseMap( uint8_t i, uint8_t j, CharId &rank, CharId &count )
{
    rank = baseCounts[i][j];
    count = midRanks[i][j] - rank;
}

void IndexReader::setBaseOverlap( uint8_t i, uint8_t j, CharId &rank, CharId &count )
{
    rank = midRanks[i][j];
    count = baseCounts[ i + 1 ][j] - rank;
}

void IndexReader::setRank( uint8_t i, CharId rank, CharCount &ranks )
{
    rank += charRanks[i];
    ReadId rankMark = rank / indexPerMark;
    ReadId rankIndex = marks_[rankMark];
    
    CharId totalCount = setRankIndex( rankIndex, ranks );
    assert( totalCount <= rank );
    if ( rank - totalCount >= bwtPerIndex && rankIndex + 1 < indexSize )
    {
        CharCount tmpRanks;
        CharId tmpTotal = setRankIndex( rankIndex + 1, tmpRanks );
        while ( tmpTotal <= rank )
        {
            memcpy( &ranks.counts, &tmpRanks.counts, 32 );
            memcpy( &ranks.endCounts, &tmpRanks.endCounts, 4 );
            totalCount = tmpTotal;
            ++rankIndex;
            if ( rank - totalCount < bwtPerIndex || rankIndex + 1 == indexSize ) break;
            tmpTotal = setRankIndex( rankIndex + 1, tmpRanks );
        }
    }
    
    ReadId rankBwt = rankIndex * bwtPerIndex;
    ReadId rankChunk = bwtPerIndex;
    uint8_t offset = index_[(rankIndex * sizePerIndex)+36];
    rankBwt -= offset;
    rankChunk += offset;
    CharId rankLeft = rank - totalCount;
    if ( rankLeft < rankChunk ) rankChunk = rankLeft;
    
    fseek( bwt, rankBwt + beginBwt, SEEK_SET );
    fread( buff, 1, rankChunk, bwt );
    
    uint8_t c;
    ReadId p = 0, thisRun, addRun;
    
    while ( rankLeft )
    {
        c = decodeBaseChar[ buff[p] ];
        thisRun = decodeBaseRun[ buff[p] ];
        if ( isBaseRun[ buff[p++] ] )
        {
            addRun = buff[p] & runMask;
            uint8_t byteCount = 0;
            while ( buff[p++] & runFlag )
            {
                addRun ^= ( buff[p] & runMask ) << ( 7 * ++byteCount );
            }
            thisRun += addRun;
        }
        
        if ( thisRun > rankLeft ) thisRun = rankLeft;
        if ( c == 4 )
        {
            ranks.endCounts += thisRun;
        }
        else
        {
            ranks.counts[c] += thisRun;
        }
        rankLeft -= thisRun;
    }
}

CharId IndexReader::setRankIndex( ReadId rankIndex, CharCount &ranks )
{
    CharId indexBegin = rankIndex * sizePerIndex;
    memcpy( &ranks.counts, &index_[indexBegin], 32 );
    memcpy( &ranks.endCounts, &index_[indexBegin+32], 4 );
    return ( ranks[0] + ranks[1] + ranks[2] + ranks[3] + ranks.endCounts );
}

