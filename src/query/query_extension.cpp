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
#include <cassert>

extern Parameters params;

Ext::Ext( string& seq, ReadId id, int ol, bool drxn )
: ext_( drxn ? seq.substr( ol ) : seq.substr( 0, seq.size()-ol ) ), reads_{ ExtRead( id, seq.size()-ol, ol ) }
{}

Ext::Ext( Ext* base, int skim, bool drxn )
: ext_( base->ext_ ), exts_( base->exts_ ), reads_( base->reads_.begin()+skim, base->reads_.end() )
{
    int ext = base->reads_[skim-1].ext_;
    assert( skim < base->reads_.size() );
    base->reads_.erase( base->reads_.begin()+skim, base->reads_.end() );
    base->exts_.clear();
    base->ext_ = drxn ? base->ext_.substr( 0, ext ) : base->ext_.substr( base->ext_.size()-ext );
    for ( int i = 0; i < base->redundant_.size(); i++)
    {
        if ( reads_[0].ol_ <= base->redundant_[i].ol_ ) redundant_.push_back( base->redundant_[i] );
        if ( ext < base->redundant_[i].ext_ ) base->redundant_.erase( base->redundant_.begin()+i-- );
    }
    shift( ext, drxn );
}

Ext::~Ext()
{
    for ( Ext* e : exts_ ) delete e;
}


bool Ext::add( ReadId id, int ext, int ol, int ins )
{
    if ( ( ext <= reads_.back().ext_ ) && ( ( reads_.back().ext_+reads_.back().ol_ ) > ( ext+ol ) ) ) redundant_.push_back( ExtRead( id, ext, ol ) );
    else reads_.insert( reads_.begin()+ins, ExtRead( id, ext, ol ) );
    return true;
}

int Ext::count( bool inclMp, int pairDrxn )
{
    count_ = 0;
    for ( Ext* e : exts_ ) count_ = max( count_, e->count( inclMp, pairDrxn ) );
    
    Lib* lib;
    for ( ExtRead& er : reads_ ) if ( ( lib = params.getLib( er.id_ ) ) 
                                   && ( inclMp || lib->isPe ) 
                                   && ( pairDrxn == 2 || lib->getDrxn( er.id_ ) == pairDrxn ) ) count_++;
    
    return count_;
}

void Ext::sanitise( int minOl )
{
    for ( int i = 0; i < exts_.size(); i++ ) if ( exts_[i]->reads_[0].ol_ < minOl )
    {
        delete exts_[i];
        exts_.erase( exts_.begin() + i-- );
    }
    for ( Ext* ext : exts_ ) ext->sanitise( minOl );
}

void Ext::set( string& seq, bool drxn )
{
    seq_ = drxn ? seq.substr( seq.size()-reads_[0].ol_ ) + ext_ : ext_ + seq.substr( 0, reads_[0].ol_ );
    for ( Ext* e : exts_ ) e->set( seq_, drxn );
}

void Ext::shift( int ext, bool drxn )
{
    assert( ext < ext_.size() );
    ext_ = drxn ? ext_.substr( ext ) : ext_.substr( 0, ext_.size()-ext );
    for ( ExtRead& er : reads_ )
    {
        er.ext_ -= ext;
        er.ol_ += ext;
    }
    for ( ExtRead& er : redundant_ )
    {
        er.ext_ -= ext;
        er.ol_ += ext;
    }
}

Exts::Exts( string& base, int coord, bool drxn )
: seq_( drxn ? base.substr( 0, coord ) : base.substr( base.size()-coord ) ), coord_( drxn ? coord : base.size() - coord )
{
    
}

Exts::~Exts()
{
    for ( Ext* e : exts_ ) delete e;
}

bool Exts::add( vector<Ext*>& exts, string seq, ReadId id, int ol, bool drxn )
{
    int ext = seq.size()-ol;
    bool added = false;
    for ( Ext* e : exts )
    {
        int i = 0, j = 0, limit = min( (int)e->ext_.size(), ext );
        while ( i < limit && ( drxn ? seq[i+ol] == e->ext_[i] : seq.end()[-i-ol-1] == e->ext_.end()[-i-1] ) ) i++;
        
        // Perfect match
        if ( i+ol == seq.size() && ( added = true ) && e->add( id, ext, ol, e->reads_.size() ) ) continue;
        
        // Matched and surpassed the existing ext seq
        if ( i == e->ext_.size() )
        {
            // Continue to pre-existing branches
            if ( !e->exts_.empty() ) return add( e->exts_, seq, id, i+ol, drxn );
            
            // append unbranched extend
            e->ext_ = drxn ? seq.substr( ol ) : seq.substr( 0, seq.size()-ol );
            return e->add( id, ext, ol, e->reads_.size() );
        }
        
        while ( j < e->reads_.size() && e->reads_[j].ext_ <= i ) j++;
        
        // Fork after read
        if ( j && i+ol < seq.size() )
        {
            Ext* alt = new Ext( e, j, drxn );
            e->exts_ = { alt, new Ext( seq, id, ol, drxn ) };
            e->exts_.back()->shift( e->ext_.size(), drxn );
            assert( !e->exts_[0]->ext_.empty() && !e->exts_[1]->ext_.empty() );
            return true;
        }
    }
    
    // Create new branch
    if ( !added ) exts.push_back( new Ext( seq, id, ol, drxn ) );
    
    return true;
}

//bool Exts::add( vector<Ext*>& exts, string seq, ReadId id, int ol, bool drxn )
//{
//    int ext = seq.size()-ol;
//    vector<Ext*> hits;
//    for ( Ext* e : exts )
//    {
//        int i = 0, j = 0, limit = min( (int)e->ext_.size(), ext );
//        while ( i < limit && ( drxn ? seq[i+ol] == e->ext_[i] : seq.end()[-i-ol-1] == e->ext_.end()[-i-1] ) ) i++;
//        
//        // Equalled or surpassed the existing ext seq
//        if ( i == e->ext_.size() )
//        {
//            assert( hits.empty() );
//            // Continue to pre-existing branches
//            if ( !e->exts_.empty() ) return add( e->exts_, seq, id, i+ol, drxn );
//            
//            // append unbranched extend
//            e->ext_ = drxn ? seq.substr( ol ) : seq.substr( 0, seq.size()-ol );
//            return e->add( id, ext, ol, e->reads_.size() );
//        }
//        
//        while ( j < e->reads_.size() && e->reads_[j].ext_ <= i ) j++;
//        
//        // Fits within ext or splits it
//        if ( j )
//        {
//            assert( hits.empty() );
//            // Add read between others
//            if ( i+ol == seq.size() ) return e->add( id, ext, ol, j );
//            
//            // Fork after read
//            Ext* alt = new Ext( e, j, drxn );
//            e->exts_ = { alt, new Ext( seq, id, ol, drxn ) };
//            e->exts_.back()->shift( e->ext_.size(), drxn );
//            assert( !e->exts_[0]->ext_.empty() && !e->exts_[1]->ext_.empty() );
//            return true;
//        }
//        
//        // Potentially a fork
//        if ( i+ol == seq.size() ) hits.push_back( e );
//    }
//    
//    // Add read to single branch
//    if ( hits.size() == 1 ) return hits[0]->add( id, ext, ol, 0 );
//    
//    // Fork multiple branches
//    for ( Ext* e : hits )
//    {
//        e->shift( ext, drxn );
//        exts.erase( remove( exts.begin(), exts.end(), e ), exts.end() );
//    }
//    
//    // Create new branch
//    exts.push_back( new Ext( seq, id, ol, drxn ) );
//    exts.back()->exts_ = hits;
//    
//    return true;
//}

void Overlap::offset( int offset, bool drxn )
{
    seq = ( drxn ? seq.substr( offset ) : seq.substr( 0, seq.length() - offset ) );
    overLen += offset;
    extLen -= offset;
}

void Overlap::sortByExt( vector<Overlap> &ols )
{
    sort( ols.begin(), ols.end(), []( const Overlap &a, const Overlap &b ) { 
        return a.extLen == b.extLen ? a.overLen > b.overLen : a.extLen > b.extLen; 
    });
    
    for ( int i = 0; i < ols.size(); )
    {
        bool doErase = false;
        for ( int j = i + 1; j < ols.size(); j++ )
        {
            if ( ols[i].readId == ols[j].readId )
            {
                doErase = true;
                break;
            }
        }
        if ( doErase ) ols.erase( ols.begin() + i );
        else i++;
    }
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

Extension::Extension( Overlap &overlap, int readCount, bool drxn, vector<Extension> &inFwdExts )
: Extension( overlap, readCount, drxn )
{
    for ( const Extension &ext : inFwdExts )
    {
        maxExtLen = max( maxExtLen, ext.maxExtLen );
        fwdExts.push_back( ext );
    }
}

void Extension::test()
{
    for ( Extension &ext : fwdExts ) ext.test();
    bool good = false;
    for ( Overlap &ol : overlaps )
    {
        if ( ol.overLen == maxOverLen ) good = true;
    }
    assert( good );
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

void Extension::cullFwd( int cutoff )
{
    if ( maxExtLen - seq.length() < cutoff || fwdExts.size() > 4 ) fwdExts.clear();
    for ( Extension &ext : fwdExts )
    {
        if ( ext.readCount > 4 ) continue;
        fwdExts.clear();
        return;
    }
    for ( Extension &ext : fwdExts )
    {
        ext.cullFwd( cutoff );
    }
    
}

void Extension::getExtendLens( vector<int> &lens )
{
    for ( Overlap &ol : overlaps ) lens.push_back( ol.extLen );
    for ( Extension &ext : fwdExts ) ext.getExtendLens( lens );
}

//bool Extension::getExtendSum() const
//{
//    int sum = 0;
//    for ( const Overlap &overlap : overlaps )
//    {
//        sum += overlap.extLen;
//    }
//    return sum;
//}

int Extension::getExtendVolume( int minExt )
{
    int vol = 0;
    for ( int i = fwdExts.empty() && !valid; i < overlaps.size(); i++ )
    {
        if ( overlaps[i].extLen < minExt ) break;
        if ( overlaps[i].redundant ) continue;
        vol += overlaps[i].extLen - minExt;
    }
    
    for ( Extension &ext : fwdExts ) vol += ext.getExtendVolume( minExt );
    
    return vol;
}

//bool Extension::getOverlapSum() const
//{
//    int sum = 0;
//    for ( const Overlap &overlap : overlaps )
//    {
//        sum += overlap.overLen;
//    }
//    return sum;
//}

int Extension::getReadCount()
{
    int maxCounted = 0;
    for ( Extension& ext : fwdExts ) maxCounted = max( maxCounted, ext.getReadCount() );
    return maxCounted + overlaps.size();
}

string Extension::getSeq( int extCount )
{
    vector<int> lens;
    getExtendLens( lens );
    sort( lens.rbegin(), lens.rend() );
    
    int len = seq.length();
    if ( lens.size() > extCount ) len = min( len, lens[extCount] );
    else len = min( len, 1 );
    
    return seq.substr( 0, len );
}

bool Extension::isCongruent( Overlap &overlap )
{
    int len = min( overlap.extLen, (uint16_t)seq.length() );
    return ( drxn ? overlap.seq.substr( 0, len ) == seq.substr( 0, len ) 
                  : overlap.seq.substr( overlap.seq.length() - len ) == seq.substr( seq.length() - len ) );
}

bool Extension::isRepeat( string &inSeq )
{
    uint16_t maxLen = min( params.readLen, (int)min( inSeq.length(), seq.length() + maxOverLen ) );
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

bool Extension::isSuperior( Extension &ext, int minLen, float expectedPer, int cutoff )
{
    if ( !cutoff )
    {
        return readCount >= ext.readCount * 5;
    }
    
    int len = max( maxExtLen, ext.maxExtLen ) - minLen;
    if ( len < cutoff ) return false;
    
    int expected = ( len * len * expectedPer / 8 ) - len;
    int extVolumes[2] = { min( expected, getExtendVolume( minLen ) )
                        , min( expected, ext.getExtendVolume( minLen ) ) };
    if ( readCount <= 1 ) extVolumes[0] = min( extVolumes[0], maxExtLen - minLen );
    if ( ext.readCount <= 1 ) extVolumes[1] = min( extVolumes[1], ext.maxExtLen - minLen );
    
    // Not enough evidence
    if ( extVolumes[0] < cutoff * expectedPer ) return false;
    
    // Probably sequencing error
    if ( extVolumes[0] > extVolumes[1] * 4 ) return true;
    
    return maxOverLen > ext.maxOverLen + cutoff * 2;
    
//    if ( readCount >= 4 )
//    {
//        cutoff *= ext.readCount >= 3 ? 1 : float(ext.readCount + 2) / (float)6 ;
//        
//        // Lack of unique overlaps
//        if ( ext.overlaps.size() <= 1 && ext.fwdExts.empty() && ext.maxOverLen < maxOverLen * min( (float)2, (float)readCount / (float)4 ) )
//        {
//            return true;
//        }
//        
//        // Much better extension and much better overlap
//        if ( maxOverLen > ext.maxOverLen + cutoff && maxExtLen > ext.maxExtLen - cutoff )
//        {
//            return true;
//        }
//    }
//    return false;
}

bool Extension::isValid()
{
    if( !valid && readCount > 1 )
    {
        if ( fwdExts.empty() )
        {
            for ( int i( 1 ); i < overlaps.size(); i++ )
            {
                if ( overlaps[i].extLen != overlaps[0].extLen && overlaps[i].overLen > overlaps[0].overLen )
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

//bool Extension::operator <(const Extension &rhs) const
//{
//    if ( maxOverLen == rhs.maxOverLen ){
//        int lSum = getOverlapSum();
//        int rSum = rhs.getOverlapSum();
//        if ( lSum == rSum )
//        {
//            return getExtendSum() > rhs.getExtendSum();
//        }
//        return lSum > rSum;
//    }
//    return maxOverLen > rhs.maxOverLen;
//}

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

void Extension::trimToExtLen( uint16_t extCutoff )
{
    if ( overlaps[0].extLen > extCutoff && overlaps[0].extLen != overlaps.back().extLen )
    {
        fwdExts.clear();
        auto it = overlaps.begin();
        while ( (*it).redundant || ( (*it).extLen > extCutoff && (*it).extLen > overlaps.back().extLen ) )
        {
            it++;
            if ( it == overlaps.end() ) return;
        }
        if ( it != overlaps.end() ) overlaps.erase( overlaps.begin(), it );
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

