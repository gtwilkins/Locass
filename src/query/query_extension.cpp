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

#include "query_extension.h"
#include <algorithm>

extern Parameters params;

void Overlap::offset( int offset, bool drxn )
{
    seq = ( drxn ? seq.substr( offset ) : seq.substr( 0, seq.length() - offset ) );
    overLen += offset;
    extLen -= offset;
}

void Overlap::truncate( bool drxn )
{
    seq = ( drxn ? seq.substr( seq.length() - extLen ) : seq.substr( 0, extLen ) );
}

Extension::Extension( Overlap &overlap, int readCount, bool drxn )
: seq( overlap.seq ), maxOverLen( overlap.overLen ), readCount( readCount ), drxn( drxn )
{
    anyPe = params.isReadPe( overlap.readId );
    valid = readCount > 1;
    overlaps.push_back( overlap );
    maxExtLen = seq.size();
}

Extension::Extension( Overlap &overlap, int readCount, bool drxn, vector<Extension> &inFwdExts, bool doAdd )
: Extension( overlap, readCount, drxn )
{
    doAdd = doAdd && inFwdExts.size() <= 5;
    for ( const Extension &ext : inFwdExts )
    {
        maxExtLen = max( maxExtLen, ext.maxExtLen );
        if ( doAdd )
        {
            fwdExts.push_back( ext );
        }
    }
}

void Extension::addOverlap( Overlap &overlap )
{
    bool isUnique = overlaps.empty() || overlaps.back().overLen != overlap.overLen;
    anyPe = anyPe || params.isReadPe( overlap.readId );
    overlaps.push_back( overlap );
    if ( !overlap.redundant )
    {
        maxOverLen = overlap.overLen;
        readCount += isUnique;
    }
}

bool Extension::canAdd( Overlap &overlap )
{
    if ( !overlap.extLen || drxn ? 
        seq.substr( 0, overlap.extLen ) == overlap.seq : 
        seq.substr( seq.length() - overlap.extLen ) == overlap.seq )
    {
        if ( overlap.overLen < maxOverLen )
        {
            overlap.redundant = true;
        }
        return true;
    }
    return false;
}

bool Extension::checkLoop( Overlap &overlap )
{
    if ( this->overlaps[0].readId == overlap.readId )
    {
        int i = 1;
        while ( overlaps[i].extLen == overlaps[0].extLen && i < overlaps.size() )
        {
            i++;
        }
        overlaps.erase( overlaps.begin(), overlaps.begin() + i );
        return overlaps.empty();
    }
    return false;
}

bool Extension::getExtendSum() const
{
    int sum = 0;
    for ( const Overlap &overlap : overlaps )
    {
        sum += overlap.extLen;
    }
    return sum;
}

bool Extension::getOverlapSum() const
{
    int sum = 0;
    for ( const Overlap &overlap : overlaps )
    {
        sum += overlap.overLen;
    }
    return sum;
}

bool Extension::isCongruent( Overlap &overlap )
{
    int len = min( overlap.extLen, (uint16_t)seq.length() );
    return ( drxn ? overlap.seq.substr( 0, len ) == seq.substr( 0, len ) 
                  : overlap.seq.substr( overlap.seq.length() - len ) == seq.substr( seq.length() - len ) );
}

bool Extension::isRepeat( string &inSeq )
{
    uint16_t maxLen = min( params.readLen, (uint16_t)seq.length() + maxOverLen );
    uint16_t maxShift = maxLen - maxOverLen;
    uint16_t bgnCoord = drxn ? inSeq.length() - maxLen : seq.length() - maxShift;
    fullSeq = drxn ? inSeq.substr( inSeq.length() - maxOverLen ) + seq 
                   : seq + inSeq.substr( 0, maxOverLen );
    for ( uint16_t i ( 1 ); i < maxShift; i++ )
    {
        uint16_t j = 0;
        while ( ( j < maxLen - i ) && ( drxn ? inSeq[bgnCoord + i + j] == fullSeq[j]
                                             : inSeq[j] == fullSeq[bgnCoord + i + j - 1] ) )
        {
            j++;
        }
        
        if ( j == maxLen - i )
        {
            return true;
        }
    }
    
    return false;
}

bool Extension::isSuperior( Extension &ext, uint16_t cutoff )
{
    if ( readCount >= 4 )
    {
        cutoff *= ext.readCount >= 3 ? 1 : float(ext.readCount + 2) / (float)6 ;
        
        // Lack of unique overlaps
        if ( ext.overlaps.size() <= 1 && ext.fwdExts.empty() && ext.maxOverLen < maxOverLen * min( (float)2, (float)readCount / (float)4 ) )
        {
            return true;
        }
        
        // Much better extension and much better overlap
        if ( maxOverLen > ext.maxOverLen + cutoff && maxExtLen > ext.maxExtLen - cutoff )
        {
            return true;
        }
    }
    return false;
}

bool Extension::isValid()
{
    if( !valid && readCount > 1 )
    {
        if ( fwdExts.empty() )
        {
            for ( int i( 1 ); i < overlaps.size(); i++ )
            {
                if ( overlaps[i].extLen != overlaps[0].extLen )
                {
                    seq = drxn ? seq.substr( 0, overlaps[i].extLen ) : seq.substr( seq.length() - overlaps[i].extLen );
                    overlaps.erase( overlaps.begin(), overlaps.begin() + i );
                    break;
                }
            }
        }
        else
        {
            for ( Extension &ext : fwdExts )
            {
                if ( !ext.isValid() )
                {
                    fwdExts.clear();
                    break;
                }
            }
        }
        valid = true;
    }
    return valid;
}

void Extension::offsetFwdExts( int offset )
{
    if ( offset )
    {
        seq = drxn ? seq.substr( offset ) : seq.substr( 0, seq.length() - offset );
        maxOverLen += offset;
        for ( Overlap &ol : overlaps )
        {
            ol.extLen -= offset;
            ol.overLen += offset;
        }
    }
    if ( !fwdExts.empty() )
    {
        offset += seq.length();
        for ( Extension &ext : fwdExts )
        {
            ext.offsetFwdExts( offset );
        }
    }
    
}

bool Extension::operator <(const Extension &rhs) const
{
    if ( maxOverLen == rhs.maxOverLen ){
        int lSum = getOverlapSum();
        int rSum = rhs.getOverlapSum();
        if ( lSum == rSum )
        {
            return getExtendSum() > rhs.getExtendSum();
        }
        return lSum > rSum;
    }
    return maxOverLen > rhs.maxOverLen;
}

void Extension::resetMaxOverLen()
{
    maxOverLen = 0;
    for ( Overlap &read : overlaps )
    {
        maxOverLen = max( maxOverLen, read.overLen ); 
    }
}

void Extension::setFullSeq( string &inSeq )
{
    fullSeq = drxn ? inSeq.substr( inSeq.length() - maxOverLen ) + seq 
                   : seq + inSeq.substr( 0, maxOverLen );
    for ( Extension &ext : fwdExts )
    {
        ext.setFullSeq( fullSeq );
    }
}

void Extension::setScores( int minOver, int overCount )
{
    score = 0;
    vector<uint16_t> lens;
    for ( Overlap &over : overlaps )
    {
        if ( over.overLen >= minOver )
        {
            lens.push_back( over.overLen - minOver );
        }
    }
    sort( lens.begin(), lens.end(), []( uint16_t a, uint16_t b ) { return a > b; } );
    int counts = min( overCount, int( lens.size() ) );
    for ( int i( 0 ); i < counts; i++ )
    {
        score += lens[i];
    }
    if ( overCount > lens.size() && readCount > lens.size() )
    {
        score = ( score * min( readCount, overCount ) ) / lens.size();
    }
}

void Extension::trimToExtLen( uint16_t extCutoff )
{
    if ( overlaps[0].extLen > extCutoff && overlaps[0].extLen != overlaps.back().extLen )
    {
        fwdExts.clear();
        auto it = overlaps.begin();
        while ( (*it).extLen > extCutoff && (*it).extLen > overlaps.back().extLen )
        {
            it++;
        }
        overlaps.erase( overlaps.begin(), it );
    }
    else
    {
        bool doClear = true;
        for ( Extension &fwd : fwdExts )
        {
            doClear = doClear && fwd.overlaps.back().extLen > extCutoff;
        }
        if ( doClear )
        {
            fwdExts.clear();
        }
        for ( Extension &fwd : fwdExts )
        {
            fwd.trimToExtLen( extCutoff );
        }
    }
    
    uint16_t maxExt = overlaps[0].extLen;
    if ( maxExt < seq.length() )
    {
        seq = drxn ? seq.substr( 0, maxExt ) : seq.substr( seq.length() - maxExt );
    }
}

bool Extension::trimToOverlap( uint16_t &overlapCutoff )
{
    uint16_t extendLen = 0;
    vector<Overlap> newOverlaps;
    for ( Overlap &overlap : overlaps )
    {
        if ( overlap.overLen > overlapCutoff )
        {
            newOverlaps.push_back( overlap );
            extendLen = max ( overlap.extLen, extendLen );
        }
    }
    overlapCutoff += extendLen;
    overlaps = newOverlaps;
    seq = drxn ? seq.substr( 0, extendLen ) : seq.substr( seq.length() - extendLen ) ;
    fwdExts.clear();
    return overlaps.empty();
}

bool Extension::trimToSeeds( vector<Overlap> &seedOverlaps, unordered_set<SeqNum> &seeds )
{
    for ( Overlap &ol : overlaps )
    {
        if ( seeds.find( ol.readId ) != seeds.end() )
        {
            return isValid();
        }
    }
    
    for ( Overlap &ol : seedOverlaps )
    {
        if ( isCongruent( ol ) )
        {
            for ( auto it = fwdExts.begin(); it != fwdExts.end(); )
            {
                if ( !(*it).trimToSeeds( seedOverlaps, seeds ) )
                {
                    it = fwdExts.erase( it );
                    continue;
                }
                it++;
            }
            if ( fwdExts.size() == 1 )
            {
                vector<Extension> tmpFwdExts = fwdExts;
                vector<Overlap> tmpOverlaps = overlaps;
                overlaps = fwdExts[0].overlaps;
                overlaps.insert( overlaps.end(), tmpOverlaps.begin(), tmpOverlaps.end() );
                maxExtLen = fwdExts[0].maxExtLen;
                seq = fwdExts[0].seq;
                fwdExts = tmpFwdExts;
            }
            return isValid();
        }
    }
    
    return false;
}

