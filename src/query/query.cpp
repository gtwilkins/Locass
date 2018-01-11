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
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>

extern struct Parameters params;

Querier::Querier( IndexReader* ir, QueryBinaries* qb )
: reader_( ir ), bin_( qb )
{
    minOver_ = params.readLen * 0.45;
    maxSeqs_ = max( 300, int( params.cover * 20 ) );
    constCutoff_ = ( 5 + ( 7 * ( params.readLen / params.cover ) ) ) * ( params.readLen / 100 );
}

Querier::~Querier()
{
    delete reader_;
    delete bin_;
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
            
            bool anyInvalid = false;
            for ( auto it1 = fwdExts.begin(); it1 != fwdExts.end(); )
            {
                readCount += (*it1).readCount;
                bool rmv = false;
                for ( auto it2 = fwdExts.begin(); doTrim && it2 != fwdExts.end() && !rmv; it2++ )
                {
                    rmv = it1 != it2 && (*it2).isSuperior( *it1, constCutoff_ );
                }
                if ( rmv )
                {
                    it1 = fwdExts.erase( it1 );
                    continue;
                }
                anyInvalid = anyInvalid || !(*it1).isValid() || !(*it1).anyPe;
                it1++;
            }
            
            if ( fwdExts.size() == 1 )
            {
                fwdExts[0].addOverlap( ol );
                tmpExts.push_back( fwdExts[0] );
            }
            else
            {
                tmpExts.push_back( Extension( ol, readCount, drxn, fwdExts, !anyInvalid ) );
            }
            
            exts = tmpExts;
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
    }
    return exts;
}

int Querier::countKmer( string seq )
{
    int seqLen = min( params.readLen, (int)seq.length() );
    uint8_t query[seqLen];
    setQuery( seq, query, seqLen, 0 );
    QueryState q( query, seqLen, seqLen );
    CharId rank, count;
    reader_->setBaseAll( query[0], query[1], rank, count );
    int it = 1;
    while ( it < seqLen - 1 && count )
    {
        CharCount ranks;
        CharCount counts;
        uint8_t i = q.q[it];
        uint8_t j = q.q[++it];

        reader_->countRange( i, rank, count, ranks, counts );
        rank = ranks[j];
        count = counts[j];
    }
    
    return count;
}

void Querier::estimateCoverage( ofstream &fh, int kmerLen, int sampleSize )
{
    srand( time(NULL) );
    int maxKmer = 200;
    int counts[maxKmer]{0};
    int i = 0;
    while ( i < sampleSize )
    {
        ReadId id = ( ( rand() & 65535 ) << 16 | ( rand() & 65535 ) ) % params.seqCount;
        string seq = getSequence( id );
        if ( seq.size() < kmerLen ) continue;
        int count = countKmer( seq.substr( 0, kmerLen ) );
        if ( count > 200 ) continue;
        counts[count]++;
        i++;
    }
    
    for ( int j = 0; j < maxKmer; j++ )
    {
        fh << to_string( j ) << "," << to_string( counts[j] ) << endl;
    }
}

vector<Overlap> Querier::getOverlaps( string &seq, uint16_t minOver, bool drxn )
{
    // Query index
    int seqLen = min( params.readLen, (int)seq.length() );
    uint8_t query[seqLen];
    setQuery( seq, query, seqLen, drxn );
    QueryState q( query, seqLen, max( minOver, minOver_ ) );
    CharId rank, count;
    reader_->setBaseOverlap( query[0], query[1], rank, count );
    mapReads( q, 1, rank, count );
    
    // Compile overlaps
    vector<Overlap> overlaps;
    for ( int i ( 0 ); i < q.endOverlaps.size() && overlaps.size() < maxSeqs_; i++ )
    {
        if ( q.endCounts[i] > 300 ) break;
        bin_->getOverlaps( overlaps, q.endRanks[i], q.endCounts[i], q.endOverlaps[i], drxn );
    }
    sort( overlaps.begin(), overlaps.end(), []( const Overlap &a, const Overlap &b ) { 
        return a.extLen > b.extLen; 
    });
    
    return overlaps;
}

vector<Overlap> Querier::getOverlaps( string &seq, uint16_t minOver, uint8_t &maxConsidered, bool drxn )
{
    // Query index
    int seqLen = min( params.readLen, (int)seq.length() );
    uint8_t query[seqLen];
    setQuery( seq, query, seqLen, drxn );
    QueryState q( query, seqLen, max( minOver, minOver_ ) );
    CharId rank, count;
    reader_->setBaseOverlap( query[0], query[1], rank, count );
    mapReads( q, 1, rank, count );
    
    // Compile overlaps
    vector<Overlap> overlaps;
    int iLast = 0;
    for ( int i ( 0 ); i < q.endOverlaps.size() && overlaps.size() < maxSeqs_; i++ )
    {
        if ( q.endCounts[i] > 300 ) break;
        bin_->getOverlaps( overlaps, q.endRanks[i], q.endCounts[i], q.endOverlaps[i], drxn );
        iLast = i;
    }
    sort( overlaps.begin(), overlaps.end(), []( const Overlap &a, const Overlap &b ) { 
        return a.extLen > b.extLen; 
    });
    
    maxConsidered = iLast + 1 < q.endOverlaps.size() ? q.endOverlaps[iLast] : minOver;
    
    return overlaps;
}

string Querier::getSequence( ReadId id )
{
    return bin_->getSequence( id );
}

bool Querier::isVector( Extension &ext )
{
    if ( ext.readCount < params.queryMpCutoff )
    {
        vector<Overlap> overlaps = getOverlaps( ext.fullSeq, minOver_, ext.drxn );
        vector<Extension> fwdExts = compileExtensions( overlaps, ext.drxn, true );
        bool anyPe = false;
        for ( const Extension &fwdExt : fwdExts )
        {
            anyPe = anyPe || fwdExt.anyPe;
            ext.readCount += fwdExt.readCount;
        }
        if ( anyPe )
        {
            ext.anyPe = true;
            return true;
        }
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
    reader_->setBaseOverlap( query[0], query[1], rank, count );
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
    mapExtensionsCull( seq, exts );
    mapExtensionsTrim( exts, overlaps );
    
    return exts;
}

vector<Extension> Querier::mapExtensions( string &seq, bool drxn, uint16_t minOver )
{
    vector<Overlap> overlaps = getOverlaps( seq, minOver, drxn );
    vector<Extension> exts = compileExtensions( overlaps, drxn, true );
    mapExtensionsCull( seq, exts );
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

void Querier::mapExtensionsCull( string &seq, vector<Extension> &exts )
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
            doErase = it1 != it2 && (*it2).isSuperior( *it1, constCutoff_ );
        }
        if ( doErase )
        {
            it1 = exts.erase( it1 );
            continue;
        }
        it1++;
    }
    
    // Give up
    if ( exts.size() >= 15 )
    {
        exts.clear();
    }
    
    for ( Extension &ext : exts )
    {
        mapExtensionsCull( ext.fullSeq, ext.fwdExts );
    }
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
    uint16_t extCutoff = ( !overlaps.empty() ? overlaps[0].extLen * 0.9 : params.readLen );
    for ( Extension &ext : exts )
    {
        ext.trimToExtLen( extCutoff );
        ext.offsetFwdExts( 0 );
    }
}

void Querier::mapReads( QueryState &q, uint8_t it, CharId rank, CharId count )
{
    CharCount ranks;
    CharCount counts;
    uint8_t i = q.q[it];
    uint8_t j = q.q[++it];
    
    reader_->countRange( i, rank, count, ranks, counts );
    if ( counts[j] && it < q.seqLen )
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
    
    
    sort( overlaps.begin(), overlaps.end(), []( const Overlap &a, const Overlap &b ) { 
        return a.overLen > b.overLen; 
    });
    
    return overlaps;
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
        reader_->setBaseAll( query[0][qs.i-1], query[0][qs.i], qs.rank, qs.count );
        mapSeed( qs, 0, 0 );
    }
    
    for ( int i ( 0 ); i < qs.rightCount; i++ )
    {
        qs.setup( i, 1 );
        reader_->setBaseAll( query[1][qs.i-1], query[1][qs.i], qs.rank, qs.count );
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
    reader_->countRange( qs.j, qs.rank, qs.count, ranks, counts );
    
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
        for ( ReadStruct &read : bin_->getReads( ranks.endCounts, counts.endCounts, drxn ) )
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
    
    for ( int i = 0; i + params.readLen <= seqLen; i++ )
    {
        QueryState q( &query[i], params.readLen, minOver_ );
        reader_->setBaseMap( query[i], query[i+1], rank, count );
        mapReads( q, 1, rank, count );
        for ( int j = 0; j < q.endOverlaps.size(); j++ )
        {
            for ( ReadId &id : bin_->getIds( q.endRanks[j], q.endCounts[j] ) )
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
