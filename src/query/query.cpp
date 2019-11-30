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

#include "query.h"
#include "parameters.h"
#include "shared_functions.h"
#include "timer.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>

extern struct Parameters params;

Querier::Querier( IndexReader* ir, QueryBinaries* qb )
: ir_( ir ), qb_( qb )
{
    minOver_ = min( 50, int( params.readLen * 0.45 ) );
    maxSeqs_ = max( 300, int( params.cover * 20 ) );
    constCutoff_ = ( 5 + ( 700 / params.cover ) );
    expectedPer_ = params.cover / params.readLen;
    
    float cover = params.isCalibrated ? params.cover : params.readLen / 2;
    olLimits_[0] = (float)params.readLen * max( (float)1, min( (float)2, cover * 2 / 100 ) ) / 4;
    olLimits_[1] = params.readLen * 8 / 10;
    for ( int i = 0; i < params.readLen; i++ ) countLimits_[i] = cover * i / ( 2 * params.readLen );
    for ( int i = 0; i < params.readLen/2; i++ ) countLimits_[i+params.readLen/2] += cover * i * 3 / params.readLen;
    countLimits_[ olLimits_[1] ] = 10 * cover;
}

Querier::~Querier()
{
    delete ir_;
    delete qb_;
}

vector<Extension> Querier::compileExtensions( vector<Overlap> &overlaps, bool drxn, bool doTrim )
{
    vector<Extension> exts;
    
    for( Overlap &ol : overlaps )
    {
        if ( ol.seq.empty() ) continue;
        vector<int> found;
        int readCount = 1;
        bool isRedundant = false;
        
        for ( int i ( 0 ); i < exts.size(); )
        {
            if ( exts[i].checkLoop( ol ) )
            {
                exts.erase( exts.begin() + i );
                continue;
            }
            i++;
        }
        
        // Try adding to existing extensions
        for ( int i = 0; i < exts.size(); i++)
        {
            if ( exts[i].canAdd( ol ) )
            {
                bool thisRedundant = ol.overLen < exts[i].maxOverLen;
                if ( isRedundant && !thisRedundant )
                {
                    continue;
                }
                else if ( thisRedundant && !isRedundant )
                {
                    found.clear();
                }
                isRedundant = isRedundant ? : thisRedundant;
                found.push_back( i );
            }
        }
        
        if ( found.size() > 1 && !isRedundant )
        {
            vector<Extension> tmpExts, fwdExts;
            
            for ( int i( 0 ); i < exts.size(); i++ )
            {
                vector<Extension>* whichPile = find( found.begin(), found.end(), i ) == found.end() ? &tmpExts : &fwdExts;
                (*whichPile).push_back( exts[i] );
            }
            
            for ( auto it1 = fwdExts.begin(); it1 != fwdExts.end(); )
            {
                readCount += (*it1).readCount;
                bool rmv = false;
                for ( auto it2 = fwdExts.begin(); doTrim && it2 != fwdExts.end() && !rmv; it2++ )
                {
                    rmv = it1 != it2 && (*it2).isSuperior( *it1, ol.extLen, expectedPer_, constCutoff_ );
                }
                if ( rmv )
                {
                    it1 = fwdExts.erase( it1 );
                    continue;
                }
                it1++;
            }
            
            if ( fwdExts.size() == 1 )
            {
                fwdExts[0].addOverlap( ol );
                tmpExts.push_back( fwdExts[0] );
            }
            else
            {
                tmpExts.push_back( Extension( ol, readCount, drxn, fwdExts ) );
            }
            
            exts = tmpExts;
            for ( Extension &ext : exts ) ext.test();
        }
        else
        {
            for ( int i : found )
            {
                exts[i].addOverlap( ol );
            }
        }
        
        // Add new extension if none compatible or multiple
        if( found.empty() )
        {
            exts.push_back( Extension( ol, readCount, drxn ) );
        }
        for ( Extension &ext : exts ) ext.test();
    }
    
    for ( Extension &ext : exts ) ext.test();
    
    return exts;
}

vector<QueryBranch> Querier::mapBranches( string& seq, vector<bool>* block, int minLen, bool drxn )
{
    vector<QueryBranch> branches;
    
    uint8_t q[seq.size()];
    setQuery( seq, q, seq.size(), drxn );
    int len[2]{ minLen, minLen }, cutoff = min( 500, max( 100, (int)params.cover ) );
    
    for ( int i = seq.size()+1 - minLen; i-- > 0; )
    {
        if ( block && ( drxn ? block->end()[ -i-len[0] ] : (*block)[ i+len[0]-1 ] ) ) continue;
        CharId rank, count;
        CharCount ranks, counts;
        for ( int ol = ir_->primeOverlap( &q[i], rank, count ); ol <= len[0]; ol++ )
        {
            bool excess = i && count > cutoff;
            ir_->countRange( q[i+ol-1], rank, count, ranks, counts );
            if ( ol >= minLen )
            {
                bool blocked = block && ( drxn ? block->end()[ -i-ol ] : (*block)[ i+ol-1 ] );
                len[1] = excess && !blocked ? ol+1 : min( len[1], ol );
                if ( counts.endCounts && !excess && !blocked ) branches.push_back( QueryBranch( ranks.endCounts, counts.endCounts, drxn ? seq.size()-i-ol: i, ol ) );
            }
            if ( ol == len[0] || !counts[ q[i+ol] ] ) break;
            rank = ranks[ q[i+ol] ];
            count = counts[ q[i+ol] ];
        }
        len[0] = len[1];
    }
    
    return branches;
}

CorrectionStruct Querier::mapCorrection( string seq, int len, bool drxn )
{
    CorrectionStruct c;
    if ( len < 16 ) return c;
    uint8_t query[seq.length()];
    setQuery( seq, query, seq.length(), drxn );
    QueryCorrectState q( query, len, seq.length(), min( 32, len ) );
    CharId rank, count;
    ir_->setBaseOverlap( q.q[0], q.q[1], rank, count );
    if ( !mapCorrection( q, 1, q.q[1], rank, count ) )
    {
        c.overabundant = true;
        return c;
    }
    if ( q.endCutoff < q.minOver )
    {
        c.error = true;
        return c;
    }
    if ( q.altCutoff >= q.endCutoff )
    {
        ( q.endCutoff < 40 ? c.error : c.fork ) = true;
        if ( q.endCutoff < 64 ) return c;
    }
    
    auto it = q.ends.begin();
    for ( int i ( 0 ); i < q.endOverlaps.size(); i++ )
    {
        if ( c.reads.size() + q.endCounts[i] > 300 )
        {
            if ( q.endOverlaps[i] > 60 ) c.clear();
            break;
        }
        int readCount = c.reads.size();
        qb_->getReads( c.reads, q.endRanks[i], q.endCounts[i], q.endOverlaps[i], len, drxn );
        if ( q.endOverlaps[i] <= len ) continue;
        if ( it == q.ends.end() || *it != i ) c.addReads( readCount, !drxn );
        else
        {
            c.exts[!drxn].push_back( CorrectionExt( c.reads.back()[!drxn], c.reads, readCount, drxn ) );
            it++;
        }
    }
    
    c.setExts( drxn );
    
    return c;
}

bool Querier::mapCorrection( QueryCorrectState &q, int it, uint8_t i, CharId rank, CharId count )
{
    CharCount ranks;
    CharCount counts;
    
    ir_->countRange( i, rank, count, ranks, counts );
    it++;
    
    bool fresh = false;
    for ( int j = 0; j < 4; j++ )
    {
        if ( !counts[j] || ( it < q.corrLen && j != q.q[it] ) ) continue;
        if ( it == q.corrLen && counts[j] > ( it + 10 ) * 10 ) return false;
        ( fresh ? q.fresh : fresh ) = true;
        if ( !mapCorrection( q, it, j, ranks[j], counts[j] ) ) return false;
    }
    
    if ( it < q.minOver ) return true;
    
    if ( counts.endCounts )
    {
        if ( q.fresh ) q.ends.push_back( q.endOverlaps.size() );
        if ( q.endCount < 2 && q.endCount + counts.endCounts > 1 ) q.endCutoff = it;
        q.fresh = false;
        q.endOverlaps.push_back( it );
        q.endRanks.push_back( ranks.endCounts );
        q.endCounts.push_back( counts.endCounts );
        q.endCount += counts.endCounts;
    }
    
    if ( it >= q.corrLen ) return true;
    
    for ( int j = 0; j < 4; j++ )
    {
        if ( !counts[j] || j == q.q[it] ) continue;
        if ( !q.altCutoff && ( !q.alts.empty() || counts[j] > 1 ) ) q.altCutoff = it;
        q.alts.push_back( j );
        q.altIts.push_back( it );
        q.altRanks.push_back( ranks[j] );
        q.altCounts.push_back( counts[j] );
    }
    
    return true;
}

string Querier::getConsensusExtend( QueryState &q, bool drxn )
{
    vector<Overlap> overlaps;
    for ( int i ( 0 ); i < q.endOverlaps.size() && overlaps.size() < 100; i++ )
    {
        if ( q.endCounts[i] > 300 ) break;
        qb_->getOverlaps( overlaps, q.endRanks[i], q.endCounts[i], q.endOverlaps[i], 1 );
    }
    Overlap::sortByExt( overlaps );
    
    vector<string> seqs;
    vector< vector<int> > lens;
    int cutoff = 5;
    
    for ( Overlap &ol : overlaps )
    {
        vector<int> hits;
        int maxReadCount = 0;
        for ( int i = 0; i < seqs.size(); i++ )
        {
            bool isCongruent = true;
            for ( int j = 0; j < ol.seq.length() && j < seqs[i].length(); j++ )
            {
                if ( ol.seq[j] != seqs[i][j] ) isCongruent = false;
            }
            if ( isCongruent )
            {
                hits.push_back( i );
                maxReadCount = max( maxReadCount, (int)lens[i].size() );
            }
        }
        
        if ( hits.size() > 1 && maxReadCount > 9 )
        {
            int erased = 0;
            for ( int i = 0; i < hits.size(); )
            {
                hits[i] -= erased;
                if ( lens[ hits[i] ].size() == 1 )
                {
                    seqs.erase( seqs.begin() + hits[i] );
                    lens.erase( lens.begin() + hits[i] );
                    hits.erase( hits.begin() + i );
                    erased++;
                }
                else i++;
            }
        }
        
        if ( hits.empty() )
        {
            seqs.push_back( ol.seq );
            lens.push_back( vector<int>( { (int)ol.seq.length() } ) );
        }
        else if ( hits.size() == 1 )
        {
            lens[ hits[0] ].push_back( ol.seq.length() );
        }
        else
        {
            for ( auto it = hits.rbegin(); it < hits.rend() - 1; it++ )
            {
                for ( int i = 0; i < seqs[*it].length() && i < seqs[ hits[0] ].length(); i++ )
                {
                    if ( seqs[ hits[0] ][i] != seqs[*it][i] )
                    {
                        seqs[ hits[0] ] = seqs[ hits[0] ].substr( 0, i );
                    }
                }
                lens[ hits[0] ].insert( lens[ hits[0] ].end(), lens[*it].begin(), lens[*it].end() );
                seqs.erase( seqs.begin() + *it );
                lens.erase( lens.begin() + *it );
            }
        }
    }
    
    string seq;
    vector<int> seqLens;
    for ( int i = 0; i < seqs.size(); i++ )
    {
        if ( lens[i].size() > 1 )
        {
            if ( seqLens.empty() )
            {
                seq = seqs[i];
                seqLens.insert( seqLens.end(), lens[i].begin(), lens[i].end() );
            }
            else
            {
                for ( int j = 0; j < seq.length() && j < seqs[i].length(); j++ )
                {
                    if ( seq[j] != seqs[i][j] ) seq = seq.substr( 0, j );
                }
            }
        }
        else cutoff++;
    }
    
    sort( seqLens.rbegin(), seqLens.rend() );
    int cut = 0;
    if ( seqLens.size() > cutoff ) cut = min( seqLens[cutoff], (int)seq.length() );
    
    seq = seq.substr( 0, cut );
    if ( !seq.empty() && !drxn ) revComp( seq );
    
    return seq;
}

vector<ReadId> Querier::getIds( CharId rank, CharId count, bool drxn )
{
    vector<ReadId> ids = qb_->getIds( rank, count );
    for ( ReadId& id : ids )
    {
        bool idRev = id & 0x1;
        bool isRev = idRev == drxn;
        id = id - idRev + isRev;
    }
    return ids;
}

vector<Overlap> Querier::getOverlaps( string &seq, uint16_t minOver, bool drxn )
{
    // Query index
    int seqLen = min( params.readLen, (int)seq.length() );
    uint8_t query[seqLen];
    setQuery( seq, query, seqLen, drxn );
    QueryState q( query, seqLen, max( minOver, minOver_ ) );
    CharId rank, count;
    ir_->setBaseOverlap( query[0], query[1], rank, count );
    mapReads( q, 1, rank, count );
    
    // Compile overlaps
    vector<Overlap> ols;
    for ( int i ( 0 ); i < q.endOverlaps.size() && ols.size() < maxSeqs_; i++ )
    {
        if ( q.endCounts[i] > 300 ) break;
        qb_->getOverlaps( ols, q.endRanks[i], q.endCounts[i], q.endOverlaps[i], drxn );
    }
    Overlap::sortByExt( ols );
    
    return ols;
}

vector<Overlap> Querier::getOverlaps( string &seq, uint16_t minOver, uint8_t &maxConsidered, bool drxn )
{
    // Query index
    int seqLen = min( params.readLen, (int)seq.length() );
    uint8_t query[seqLen];
    setQuery( seq, query, seqLen, drxn );
    QueryState q( query, seqLen, max( minOver, minOver_ ) );
    CharId rank, count;
    ir_->setBaseOverlap( query[0], query[1], rank, count );
    mapReads( q, 1, rank, count );
    
    // Compile overlaps
    vector<Overlap> ols;
    int iLast = 0;
    for ( int i ( 0 ); i < q.endOverlaps.size() && ols.size() < maxSeqs_; i++ )
    {
        if ( q.endCounts[i] > 300 ) break;
        qb_->getOverlaps( ols, q.endRanks[i], q.endCounts[i], q.endOverlaps[i], drxn );
        iLast = i;
    }
    Overlap::sortByExt( ols );
    
    maxConsidered = iLast + 1 < q.endOverlaps.size() ? q.endOverlaps[iLast] : minOver;
    
    return ols;
}

string Querier::getSequence( ReadId id )
{
    return qb_->getSequence( id );
}

bool Querier::isVector( Extension &ext )
{
    return false;
    if ( ext.readCount < params.queryMpCutoff ) return false;
    vector<Overlap> overlaps = getOverlaps( ext.fullSeq, minOver_, ext.drxn );
    if ( overlaps.size() < params.queryMpCutoff ) return false;
    for ( Overlap &ol : overlaps )
    {
        if ( params.isReadPe( ol.readId ) ) return false;
    }
    
    return true;
}

bool Querier::isExtendable( string& seq, bool drxn )
{
    CharId rank, count;
    CharCount ranks, counts;
    vector<uint8_t> q;
    int ol;
    ir_->primeOverlap( seq, q, rank, count, ol, drxn );
    while ( count )
    {
        if ( ++ol > seq.size() ) return true;
        ir_->countRange( q.back(), rank, count, ranks, counts );
        if ( ol > olLimits_[0] && counts.endCounts ) return true;
        
        q.push_back( drxn ? charToInt[ seq.end()[-ol] ] : charToIntComp[ seq[ol-1] ] );
        rank = ranks[q.back()];
        count = counts[q.back()];
    }
    return false;
}

ReadId Querier::isExtendable( string &seq, uint16_t minLen, bool drxn )
{
    int seqLen = min( params.readLen, (int)seq.length() );
    uint8_t query[seqLen];
    setQuery( seq, query, seqLen, drxn );
    QueryState q( query, seqLen, max( minLen, minOver_ ) );
    CharId rank, count;
    ir_->setBaseOverlap( query[0], query[1], rank, count );
    mapReads( q, 1, rank, count );
    
    ReadId extCount = 0;
    for ( int i ( 0 ); i < q.endOverlaps.size(); i++ )
    {
        extCount += q.endCounts[i];
    }
    
    return extCount;
}

vector<Extension> Querier::mapExtensions( bool &noMatches, string &seq, bool drxn, uint16_t minOver )
{
    vector<Overlap> overlaps = getOverlaps( seq, minOver, drxn );
    noMatches = overlaps.empty();
    vector<Extension> exts = compileExtensions( overlaps, drxn, true );
    mapExtensionsCull( seq, exts, 0 );
    mapExtensionsTrim( exts, overlaps );
    
    return exts;
}

vector<Extension> Querier::mapExtensions( string &seq, bool drxn, uint16_t minOver )
{
    vector<Overlap> overlaps = getOverlaps( seq, minOver, drxn );
    vector<Extension> exts = compileExtensions( overlaps, drxn, true );
    mapExtensionsCull( seq, exts, 0 );
    mapExtensionsTrim( exts, overlaps );
    
    return exts;
}

vector<Extension> Querier::mapExtensions( string &seq, bool drxn, unordered_set<SeqNum> &seeds, uint16_t minOver )
{
    vector<Overlap> overlaps = getOverlaps( seq, minOver, drxn );
    vector<Extension> exts = compileExtensions( overlaps, drxn, false );
    mapExtensionsCull( seq, exts, overlaps, seeds );
    mapExtensionsTrim( exts, overlaps );
    
    return exts;
}

void Querier::mapExtensionsCull( string &seq, vector<Extension> &exts, int base )
{
    // Remove extensions with a repeated frame
    for ( auto it = exts.begin(); it != exts.end(); )
    {
        if ( (*it).isRepeat( seq ) || ( !(*it).anyPe && isVector( *it ) ) )
        {
            it = exts.erase( it );
            continue;
        }
        it++;
    }
    
    // Validate extensions
    int maxReads = 0;
    for ( Extension &ext : exts )
    {
        maxReads = max( maxReads, ext.readCount );
    }
    
    // Give up due to too many garbage extensions with none good
    if ( exts.size() >= 15 && maxReads < 5 )
    {
        exts.clear();
    }
    
    // Remove weak extensions
    bool eraseDefault = exts.size() >= 15;
    for ( auto it1 = exts.begin(); it1 != exts.end(); )
    {
        bool doErase = eraseDefault && !(*it1).isValid();
        for ( auto it2 = exts.begin(); !doErase && it2 != exts.end(); it2++ )
        {
            doErase = it1 != it2 && (*it2).isSuperior( *it1, base, expectedPer_, constCutoff_ );
        }
        if ( doErase ) it1 = exts.erase( it1 );
        else it1++;
    }
    
    // Give up
    if ( exts.size() >= 15 )
    {
        exts.clear();
    }
    
    for ( Extension &ext : exts ) ext.cullFwd( constCutoff_ );
}

void Querier::mapExtensionsCull( string &seq, vector<Extension> &exts, vector<Overlap> &overlaps, unordered_set<SeqNum> &seeds )
{
    // Remove extensions with a repeated frame
    for ( auto it = exts.begin(); it != exts.end(); )
    {
        if ( (*it).isRepeat( seq ) || ( !(*it).anyPe && isVector( *it ) ) )
        {
            it = exts.erase( it );
            continue;
        }
        it++;
    }
    vector<Overlap> seedOverlaps;
    for ( Overlap &ol : overlaps )
    {
        if ( seeds.find( ol.readId ) != seeds.end() )
        {
            seedOverlaps.push_back( ol );
        }
    }
    
    for ( auto it = exts.begin(); it != exts.end(); )
    {
        if ( !(*it).trimToSeeds( seedOverlaps, seeds ) )
        {
            exts.erase( it );
            continue;
        }
        it++;
    }
}

void Querier::mapExtensionsTrim( vector<Extension> &exts, vector<Overlap> &overlaps )
{
    for ( Extension &ext : exts ) ext.test();
    uint16_t extCutoff = ( !overlaps.empty() ? overlaps[0].extLen * 0.9 : params.readLen );
    for ( Extension &ext : exts )
    {
        ext.trimToExtLen( extCutoff );
        ext.offsetFwdExts( 0 );
    }
    for ( Extension &ext : exts ) ext.test();
}

void Querier::mapReads( QueryState &q, uint8_t it, CharId rank, CharId count )
{
    CharCount ranks;
    CharCount counts;
    uint8_t i = q.q[it];
    uint8_t j = q.q[++it];
    
    ir_->countRange( i, rank, count, ranks, counts );
    if ( j < 4 && counts[j] && it < q.seqLen )
    {
        mapReads( q, it, ranks[j], counts[j] );
    }
    
    if ( it >= q.minOver && counts.endCounts )
    {
        q.endOverlaps.push_back( it );
        q.endRanks.push_back( ranks.endCounts );
        q.endCounts.push_back( counts.endCounts );
    }
}

vector<Overlap> Querier::mapJoin( string seq1, string seq2, uint8_t overLen )
{
    vector<Overlap> overlaps;
    uint8_t currExt = 0;
    unordered_set<ReadId> usedIds;
    
    while ( overLen + currExt < params.readLen )
    {
        uint8_t lastOverlap = max( uint8_t( params.readLen * 0.8 ), uint8_t( overLen + currExt + 1 ) );
        vector<Overlap> currOverlaps = getOverlaps( seq1, lastOverlap, lastOverlap, 1 );
        for ( Overlap &over : currOverlaps )
        {
            over.overLen -= currExt;
            over.extLen += currExt;
            if ( seq2.find( over.seq ) == 0 && usedIds.find( over.readId ) == usedIds.end() )
            {
                overlaps.push_back( over );
                usedIds.insert( over.readId );
            }
        }
        uint8_t lastExt = params.readLen - lastOverlap;
        currExt += lastExt;
        if ( lastExt >= seq2.length() || currExt + overLen + 1 == params.readLen ) break;
        seq1 += seq2.substr( 0, lastExt );
        seq2 = seq2.substr( lastExt );
    }
    
    
    Overlap::sortByExt( overlaps );
    
    return overlaps;
}

QueryJunction Querier::mapJunction( string &seq, bool drxn )
{
    return QueryJunction( ir_, qb_, seq, countLimits_, olLimits_, constCutoff_, drxn );
}

MappedSeqs Querier::mapSeed( string &seq, int errorRate, bool bestMatch )
{
    MappedSeqs ms;
    QuerySeedState qs( seq, ms, errorRate );
    uint8_t query[2][qs.len];
    qs.query[0] = query[0];
    qs.query[1] = query[1];
    setQuery( seq, query[0], qs.len, 0 );
    setQuery( seq, query[1], qs.len, 1 );
    
    for ( int i ( 0 ); i < qs.leftCount; i++ )
    {
        qs.setup( i, 0 );
        ir_->setBaseAll( query[0][qs.i-1], query[0][qs.i], qs.rank, qs.count );
        mapSeed( qs, 0, 0 );
    }
    
    for ( int i ( 0 ); i < qs.rightCount; i++ )
    {
        qs.setup( i, 1 );
        ir_->setBaseAll( query[1][qs.i-1], query[1][qs.i], qs.rank, qs.count );
        mapSeed( qs, 0, 1 );
    }
    
    ms.updateTethers( seq );
    if ( bestMatch )
    {
        ms.setBest( seq );
    }
    ms.sort();
    
    return ms;
}

bool Querier::mapSeed( QuerySeedState &qs, int errors, bool drxn )
{
    int i = qs.i;
    CharCount ranks, counts;
    ir_->countRange( qs.j, qs.rank, qs.count, ranks, counts );
    
    bool didAppend = false;
    for ( int k ( 0 ); k < 4; k++ )
    {
        int thisErrors = errors;
        if ( counts[k] && qs.advance( i, k, thisErrors, drxn ) )
        {
            qs.rank = ranks[k];
            qs.count = counts[k];
            didAppend = mapSeed( qs, thisErrors, drxn ) || didAppend;
        }
    }
    
    if ( i > qs.iAdd && counts.endCounts && ( didAppend || qs.doAdd( errors, drxn ) ) )
    {
        qs.i = i + 1;
        for ( ReadStruct &read : qb_->getReads( ranks.endCounts, counts.endCounts, drxn ) )
        {
            qs.add( read, drxn );
        }
        didAppend = true;
    }
    
    return didAppend;
}

void Querier::mapSequence( string &seq, vector<ReadId> &ids, vector<int32_t>* coords )
{
    // Query index
    int seqLen = seq.length();
    uint8_t query[seqLen];
    setQuery( seq, query, seqLen, 0 );
    CharId rank, count;
    
    for ( int i = 0; i + 50 <= seqLen; i++ )
    {
        while ( query[i] > 3 || query[i+1] > 3 && i < seqLen ) i++;
        if ( i + 50 >= seqLen ) break;
        QueryState q( &query[i], min( params.readLen, seqLen - i ), minOver_ );
        ir_->setBaseMap( query[i], query[i+1], rank, count );
        mapReads( q, 1, rank, count );
        for ( int j = 0; j < q.endOverlaps.size(); j++ )
        {
            for ( ReadId &id : qb_->getIds( q.endRanks[j], q.endCounts[j] ) )
            {
                ids.push_back( ( id & 0x1 ? id - 1: id + 1 ) );
                coords[0].push_back( i );
                coords[1].push_back( i + q.endOverlaps[j] );
            }
        }
    }
}

void Querier::setQuery( string &seq, uint8_t* query, int length, bool drxn )
{
    if ( drxn )
    {
        for ( int i ( 0 ); i < length; i++ )
        {
            query[i] = charToInt[ seq.end()[-i-1] ];
        }
    }
    else
    {
        for ( int i ( 0 ); i < length; i++ )
        {
            query[i] = charToIntComp[ seq[i] ];
        }
    }
}
