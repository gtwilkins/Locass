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

CalibrateWriter::~CalibrateWriter()
{
    for ( NodeList &nodes : nodes_ )
    {
        for ( Node* node : nodes )
        {
            delete node;
        }
        nodes.clear();
    }
}

void CalibrateWriter::coverage()
{
    params.locusLimits[0] = -30000;
    params.locusLimits[1] = 30000;
    srand( time(NULL) );
    
    int attempted = 0;
    int64_t totalLen = 0;
    NodeList dummy;
    while ( nodes_.size() < 100 && attempted < 10000 )
    {
        ReadId id = ( ( rand() & 65535 ) << 16 | ( rand() & 65535 ) ) % params.seqCount;
        string seq = bwt_.getSequence( id );
        if ( seq.length() < params.readLen ) continue;
        vector<Extension> exts = bwt_.mapExtensions( seq, 0 );
        if ( exts.size() != 1 || exts[0].readCount <= 1 ) continue;
        int32_t limits[2]{0};
        seq = seq.substr( 0, exts[0].maxOverLen );
        NodeList nodes = { new Node( seq, exts[0], 0, 0 ) };
        nodes[0]->drxn_ = 2;
        ExtVars ev( nodes, dummy, limits, bwt_, false, false );
        nodes[0]->extendCount_ = 50;
        nodes[0]->extendNode( ev, 1 );
        nodes[0]->extendCount_ = 50;
        if ( nodes[0]->ends_[1] > params.readLen ) nodes[0]->extendNode( ev, 0 );
        if ( nodes[0]->seq_.length() < params.readLen * 2 || nodes[0]->reads_.size() < 15 )
        {
            for ( Node* node : nodes ) delete node;
        }
        else
        {
//            //
//            LocusLibraryCounts locusLists;
//            NodeList subGraphs[3];
//            subGraphs[2].push_back( nodes[0] );
//            for ( Node* fwd : nodes[0]->getDrxnNodes( 0 ) )
//            {
//                subGraphs[0].push_back( fwd); 
//            }
//            for ( Node* fwd : nodes[0]->getDrxnNodes( 1 ) )
//            {
//                subGraphs[1].push_back( fwd); 
//            }
//            locusLists.libs.push_back( LocusLibraryCount( params ) );
//            Locus* locus = new Locus( bwt_, subGraphs );
//            locus->calibrate( locusLists.libs.back() );
//            assert( false );
//            //
            
            nodes_.push_back( nodes );
            totalLen += nodes[0]->seq_.length();
            seqs_.push_back( seq );
        }
        ++attempted;
    }
    
    sort( nodes_.begin(), nodes_.end(), []( NodeList const &a, NodeList const &b ){
        return a[0]->coverage_ > b[0]->coverage_;
    } );
    
    double coverTotal = 0, coverCount = 0;
    int64_t cumulative = 0, cutoffs[2] = { int64_t(totalLen * 0.25), int64_t(totalLen * 0.75) };
    bool doCount = false;
    for ( NodeList &nodes : nodes_ )
    {
        cumulative += nodes[0]->seq_.length();
        if ( !doCount && cumulative > cutoffs[0] ) doCount = true;
        coverTotal += nodes[0]->coverage_ * nodes[0]->seq_.length();
        coverCount += nodes[0]->seq_.length();
        if ( doCount && cumulative > cutoffs[1] ) break;
    }
    
    sort( nodes_.begin(), nodes_.end(), []( NodeList const &a, NodeList const &b ){
        return a[0]->seq_.length() > b[0]->seq_.length();
    } );
    
    params.cover = coverTotal / coverCount;
}

void CalibrateWriter::pairing()
{
    params.locusLimits[0] = -30000;
    params.locusLimits[1] = 30000;
    NodeList dummy;
    SeedLibraryCounts seedList( params );
    for ( int i = 0; i < nodes_.size(); i++ )
    {
        int32_t limits[2]{0};
        ExtVars ev( nodes_[i], dummy, limits, bwt_, false, false );
        nodes_[i][0]->calibrateSeed( ev );
        SeedLibraryCount libs( params );
        Node::calibrateCount( nodes_[i], libs );
        libs.id = i;
        seedList.libs.push_back( libs );
    }
    
    LocusLibraryCounts locusLists;
    for ( int i : seedList.set( params ) )
    {
        NodeList subGraphs[3];
        subGraphs[2].push_back( nodes_[i][0] );
        for ( Node* fwd : nodes_[i][0]->getDrxnNodes( 0 ) )
        {
            subGraphs[0].push_back( fwd); 
        }
        for ( Node* fwd : nodes_[i][0]->getDrxnNodes( 1 ) )
        {
            subGraphs[1].push_back( fwd); 
        }
        
        nodes_[i].clear();
        locusLists.libs.push_back( LocusLibraryCount( params ) );
        Locus* locus = new Locus( bwt_, subGraphs );
        locus->calibrate( locusLists.libs.back() );
        delete locus;
    }
    
    locusLists.set( params );
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

