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
#include "parameters.h"
#include <algorithm>
#include <cassert>

extern struct Parameters params;

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

bool CorrectionRead::congruent( string &ext, bool drxn )
{
    if ( !exts[drxn] ) return false;
    int i = drxn ? seq.length() - exts[1] : max( 0, exts[0] - (int)ext.size() );
    int j = drxn ? 0 : max( 0, (int)ext.size() - exts[0] );
    int len = min( seq.length() - i, ext.length() - j );
    while ( seq[i++] == ext[j++] && --len );
    
    return !len;
}

CorrectionExt::CorrectionExt( string seq, CorrectionRead &read, bool drxn )
: seq( seq )
{
    lens[0] = lens[1] = 0;
    ids.push_back( read.id );  
}

CorrectionExt::CorrectionExt( string seq, vector<CorrectionRead> &reads, int i, bool drxn )
: seq( seq )
{
    for ( ; i < reads.size(); i++ ) ids.push_back( reads[i].id );
    lens[0] = ids.size() > 1 ? seq.size() : 0;
    lens[1] = ids.size() > 2 ? seq.size() : 0;
}

bool CorrectionExt::addRead( CorrectionRead &read, bool drxn )
{
    if ( !read.congruent( seq, drxn ) ) return false;
    ids.push_back( read.id );
    for ( int i = 0; i < 2; i++ )
    {
        if ( read.exts[drxn] <= lens[i] ) continue;
        for ( int j = 2; --j > i; ) lens[j] = lens[j-1];
        lens[i] = read.exts[drxn];
        break;
    }
    
    return true;
}

string CorrectionExt::getCongruent( vector<CorrectionExt> &exts, bool drxn )
{
    CorrectionExt* best = NULL;
    for ( CorrectionExt &ext : exts ) if ( !best || ext.lens[0] > best->lens[0] ) best = &ext;
    if ( !best ) return "";
    int congLen = best->lens[0];
    for ( CorrectionExt &ext : exts )
    {
        for ( int i = 0; i < congLen && i < ext.lens[0]; i++ )
        {
            if ( drxn ? best->seq[i] != ext.seq[i] : best->seq.end()[-i-1] != ext.seq.end()[-i-1] ) congLen = i;
        }
    }
    
    return drxn ? best->seq.substr( 0, congLen ) : best->seq.substr( best->seq.size()-congLen );
}

void CorrectionStruct::addReads( int i, bool drxn )
{
    for ( CorrectionExt &ext : exts[drxn] )
    {
        if ( !ext.addRead( reads[i], drxn ) ) continue;
        if ( reads.size() > i+1 ) ext.lens[1] = max( reads[i].exts[drxn], ext.lens[1] );
        for ( int j = i+1; j < reads.size(); j++ ) ext.ids.push_back( reads[j].id );
    }
}

void CorrectionStruct::clear()
{
    reads.clear();
    exts[0].clear();
    exts[1].clear();
    overabundant = true;
}

void CorrectionStruct::setExts( bool drxn )
{
    sort( reads.begin(), reads.end(), [&drxn]( CorrectionRead const &a, CorrectionRead const &b ){ return a.exts[drxn] > b.exts[drxn]; } );
    unordered_map<ReadId, int> dupes;
    int cutoff = reads.empty() ? 0 : reads[0].exts[drxn]+40;
    for ( CorrectionRead &read : reads )
    {
        auto r = dupes.insert( make_pair( read.id, read.exts[drxn] ) );
        if ( !r.second ) cutoff = min( cutoff, max( r.first->second, read.exts[drxn] ) );
    }
    if ( cutoff < 40 ) reads.clear();
    for ( auto it = reads.begin(); it != reads.end(); it++ )
    {
        if ( it->exts[drxn] >= cutoff ) continue;
        reads.erase( reads.begin(), it );
        break;
    }
    
    for ( CorrectionRead &read : reads )
    {
        bool novel = true;
        for ( CorrectionExt &ext : exts[drxn] ) if ( ext.addRead( read, drxn ) ) novel = false;
        if ( novel ) exts[drxn].push_back( CorrectionExt( read[drxn], read, drxn ) );
    }
}

int CorrectionStruct::validate( int len, int kmerLen, bool drxn )
{
    if ( len >= params.readLen-10 ) return len;
    float expFwd = ( len - kmerLen + 2 ) * params.cover * ( (float)2 / (float)params.readLen );
    float expBck = ( params.readLen - len - 10 ) * params.cover * ( (float)0.25 / (float)params.readLen );
    
    int bck = 0, fwd = 0;
    for ( CorrectionRead &read : reads ) ( read.ol < len ? fwd : bck )++;
    if ( fwd < expFwd || bck > expBck ) return len;
    
    sort( reads.begin(), reads.end(), []( CorrectionRead const &a, CorrectionRead const &b ){ return a.ol > b.ol; } );
    int cutoff = 1 + sqrt( reads.size() );
    if ( cutoff < reads.size() ) return min( len, reads[cutoff].ol );
    return 0;
}

