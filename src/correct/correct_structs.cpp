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

#include "correct_structs.h"
#include <cassert>
#include <fstream>

Fastq::Fastq( string line, uint8_t qualCutoff )
: cutoff( qualCutoff )
{
    size_t it = line.find( "paired" );
    if ( !it ) paired = true;
    if ( it ) it = line.find( "single" );
    assert( paired || !it );
    it = line.find_first_of( ' ' );
    assert( it != line.npos );
    line = line.substr( it + 1 );
    fp.open( line );
    assert( fp.is_open() && fp.good() );
}

bool Fastq::getSeq( string &seq, int nCoords[2], int qCoords[2], int i )
{
    bool anyErrors = false;
    nCoords[0] = qCoords[0] = 0;
    nCoords[1] = qCoords[1] = seqs[i].length();
    seq = seqs[i];
    
    assert( seqs[i].length() == quals[i].length() );
    int bestLen = 0, currStart = 0;
    for ( int j = 0; j < seqs[i].length(); j++ )
    {
        if ( seqs[i][j] == 'N' )
        {
            anyErrors = true;
            int len = j - currStart;
            if ( len > bestLen )
            {
                nCoords[0] = currStart;
                nCoords[1] = j;
                bestLen = len;
            }
            currStart = j + 1;
        }
    }
    
    if ( seqs[i].length() - currStart > bestLen )
    {
        nCoords[0] = currStart;
        nCoords[1] = seqs[i].length();
    }
    
    bestLen = 16;
    currStart = nCoords[0];
    for ( int j = 0; j < seqs[i].length(); j++ )
    {
        if ( quals[i][j] < cutoff )
        {
            anyErrors = true;
            if ( seq[j] == 'A' ) seq[j] = 'a';
            else if ( seq[j] == 'C' ) seq[j] = 'c';
            else if ( seq[j] == 'G' ) seq[j] = 'g';
            else if ( seq[j] == 'T' ) seq[j] = 't';
            if ( j < qCoords[0] || j >= qCoords[1] ) continue;
            int len = j - currStart;
            if ( len > bestLen )
            {
                qCoords[0] = currStart;
                qCoords[1] = j;
                bestLen = len;
            }
            currStart = j + 1;
        }
    }
    
    if ( seqs[i].length() - currStart > bestLen )
    {
        qCoords[0] = currStart;
        qCoords[1] = seqs[i].length();
    }
    
    return anyErrors;
}

bool Fastq::setNext()
{
    for ( int i = 0; i < 1 + paired; i++ )
    {
        if ( !getline( fp, seqs[i] ) )
        {
            assert( !i );
            return false;
        }
        getline( fp, seqs[i] );
        getline( fp, quals[i] );
        getline( fp, quals[i] );
    }
    
    return true;
}
