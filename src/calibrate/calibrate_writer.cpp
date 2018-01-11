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

#include "calibrate_writer.h"
#include <cassert>
#include <algorithm>

extern Parameters params;

void CalibrateWriter::coverage()
{
    vector<int> totalCounts;
    
    int minOverlap = 1 + params.readLen / 2;
    srand( time(NULL) );
    
    vector<string> seqs;
    int attempted = 0;
    while ( seqs.size() < 1000 && attempted < 10000 )
    {
        ReadId id = ( ( rand() & 65535 ) << 16 | ( rand() & 65535 ) ) % params.seqCount;
        string seq = bwt_.getSequence( id );
        if ( seq.length() < params.readLen ) continue;
        int counts[2]{0};
        for ( bool drxn : { 0, 1 } )
        {
            vector<Extension> exts = bwt_.mapExtensions( seq, drxn, minOverlap );
            if ( exts.size() == 1 )
            {
                counts[drxn] = exts[0].readCount;
            }
        }
        int diff = abs( counts[0] - counts[1] );
        if ( counts[0] > 1 && counts[1] > 1 && diff < ( 2 + min( counts[0], counts[1] ) / 2 ) )
        {
            totalCounts.push_back( counts[0] );
            totalCounts.push_back( counts[1] );
            seqs.push_back( seq );
        }
        ++attempted;
    }
    
    vector<int> sortCounts = { totalCounts.begin(), totalCounts.end() };
    sort( sortCounts.begin(), sortCounts.end() );
    
    int midPoint = sortCounts.size() / 2;
    int i = midPoint;
    int j = i + midPoint / 2;
    int totalSum = 0;
    int totalNum = 0;
    int cutoffs[2] = { sortCounts[i], sortCounts[j] };
    
    while ( i < j )
    {
        totalSum += sortCounts[i];
        totalNum++;
        i++;
    }
    
    params.cover = ( (float)totalSum / (float)totalNum ) * 2.4;
    
    for ( int k ( 0 ); k < seqs.size(); k++ )
    {
        if ( cutoffs[0] <= totalCounts[k*2] && totalCounts[k*2+1] <= cutoffs[1] )
        {
            seqs_.push_back( seqs[k] );
        }
    }
}

void CalibrateWriter::pairing()
{
    params.locusLimits[0] = -30000;
    params.locusLimits[1] = 30000;
    int minOverlap = 1 + params.readLen / 2;
    NodeListList lociNodes;
    NodeList dummy;
    SeedLibraryCounts seedList( params );
    int i = 0, id = 0;
    while ( lociNodes.size() < 200 && i < seqs_.size() )
    {
        vector<Extension> exts = bwt_.mapExtensions( seqs_[i], 1, minOverlap );
        if ( exts.size() == 1 )
        {
            int32_t limits[2]{0};
            NodeList nodes = { new Node( seqs_[i], exts[0], seqs_[i].length(), 1 ) };
            nodes[0]->drxn_ = 2;
            ExtVars ev( nodes, dummy, limits, bwt_, false, false );
            ++i;
            if ( nodes[0]->calibrateSeed( ev ) )
            {
                SeedLibraryCount libs( params );
                int thisPairCount = Node::calibrateCount( nodes, libs );
                if ( thisPairCount >= 50 )
                {
                    lociNodes.push_back( nodes );
                    libs.id = id;
                    seedList.libs.push_back( libs );
                    ++id;
                    continue;
                }
            }
            for ( Node* node : nodes )
            {
                delete node;
            }
        }
    }
    
    vector<int> goodIds = seedList.set( params );
    
    LocusLibraryCounts locusLists;
    for ( int i : goodIds )
    {
        NodeList subGraphs[3];
        subGraphs[2].push_back( lociNodes[i][0] );
        for ( Node* fwd : lociNodes[i][0]->getDrxnNodes( 0 ) )
        {
            subGraphs[0].push_back( fwd); 
        }
        for ( Node* fwd : lociNodes[i][0]->getDrxnNodes( 1 ) )
        {
            subGraphs[1].push_back( fwd); 
        }
        
        lociNodes[i].clear();
        locusLists.libs.push_back( LocusLibraryCount( params ) );
        Locus* locus = new Locus( bwt_, subGraphs );
        locus->calibrate( locusLists.libs.back() );
        delete locus;
    }
    
    locusLists.set( params );
    
    for ( NodeList &nodes : lociNodes )
    {
        for ( Node* node : nodes )
        {
            delete node;
        }
        nodes.clear();
    }
}

void CalibrateWriter::write( Filenames* fns )
{
    cout << "Estimated effective read depth: " << to_string( params.cover ) << endl << endl;
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        cout << "Paired library " << to_string( i + 1 ) << endl;
        if ( !params.libs[i].orientation || !params.libs[i].size )
        {
            cout << "\tCould not find evidence of pairing." << endl;
            continue;
        }
        else if ( params.libs[i].orientation == 1 )
        {
            cout << "\tOrientation: Forward-Forward" << endl;
        }
        else if ( params.libs[i].orientation == 2 )
        {
            cout << "\tOrientation: Forward-Reverse" << endl;
        }
        else
        {
            cout << "\tOrientation: Reverse-Forward" << endl;
        }
        cout << "\tEstimated insert length: " << params.libs[i].size << endl;
    }
    
    params.isCalibrated = true;
    uint32_t coverage = params.cover * 100000;
    FILE* fp = fns->getBinary( true, true );
    fseek( fp, 11, SEEK_SET );
    fwrite( &params.isCalibrated, 1, 1, fp );
    fwrite( &coverage, 4, 1, fp );
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        uint16_t libMed = params.libs[i].size;
        uint16_t libMin = libMed ? params.libs[i].minDist : 0;
        uint16_t libMax = libMed ? params.libs[i].maxDist : 0;
        fseek( fp, 25 + ( i * 12 ), SEEK_SET );
        fwrite( &libMed, 2, 1, fp );
        fwrite( &libMin, 2, 1, fp );
        fwrite( &libMax, 2, 1, fp );
        fwrite( &params.libs[i].orientation, 1, 1, fp );
        fwrite( &params.libs[i].isPe, 2, 1, fp );
    }
    
    fclose( fp );
}

