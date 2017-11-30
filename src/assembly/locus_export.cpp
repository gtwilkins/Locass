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

#include "locus.h"
#include <algorithm>

void Locus::blankMismatchedEnds( vector< pair<Node*, int> > &nodePath )
{
    int i = 0;
    for ( auto it = nodePath.begin(); it != nodePath.end()-1; it++ )
    {
        if ( it->second > 0 )
        {
            int bgnCoord = it->first->seq_.length() - it->second;
            string seqEnd = it->first->seq_.substr( bgnCoord );
            string seqBgn = ( it + 1 )->first->seq_.substr( 0, it->second );
            assert( seqEnd.length() == seqBgn.length() );
            if ( seqEnd != seqBgn )
            {
                for ( int i ( bgnCoord ); i < it->first->seq_.length(); i++ )
                {
                    it->first->seq_[i] = 'N';
                }
                
                for ( int i ( 0 ); i < it->second; i++ )
                {
                    ( it + 1 )->first->seq_[i] = 'N';
                }
            }
        }
        i++;
    }
}

void Locus::exportLocus( ofstream &align, ofstream &dump )
{
    int j = 0;
    int32_t ends[2] = { 0, 0 };
    
    for ( int i : { 2, 0, 1 } )
    {
        for ( Node* node : nodes_[i] )
        {
            ends[0] = min( ends[0], node->ends_[0] );
            ends[1] = max( ends[1], node->ends_[1] );
            node->id_ = to_string( j );
            j++;
        }
    }
    
    NodeList nodes;
    for ( int i : { 2, 0, 1 } )
    {
        for ( Node* node : nodes_[i] )
        {
            nodes.push_back( node );
            node->exportNodeDump( dump );
        }
    }
    
    sort( nodes.begin(), nodes.end(), []( Node* a, Node* b ){
        return a->ends_[0] < b->ends_[0];
    });
    
    for ( Node* node : nodes )
    {
        node->exportNodeAlign( align, ends[0] );
    }
}

void Locus::getContig( string &header, string &seq, uint &locusId )
{
    header = header_ + " | Locus_" + to_string( locusId );
    for ( pair<Node*, int> &nodeOverlap : getExportPath() )
    {
        seq += nodeOverlap.first->seq_.substr( 0, nodeOverlap.first->seq_.length() - max( 0, nodeOverlap.second ) );
        if ( nodeOverlap.second < 0 )
        {
            seq += string( abs( nodeOverlap.second ), 'N' );
        }
    }
}

void Locus::getExtends( string &header, string &origin, string &lSeq, string &rSeq )
{
    header = header_;
    
    string* pSeq = &lSeq;
    
    for ( pair<Node*, int> &nodeOverlap : getExportPath() )
    {
        int32_t nodeEnd = nodeOverlap.first->ends_[1] - max( 0, nodeOverlap.second );
        
        if ( nodeOverlap.first->ends_[0] < 0 )
        {
            lSeq += nodeOverlap.first->seq_.substr( 0, min( nodeEnd, 0 ) - nodeOverlap.first->ends_[0] );
            pSeq = &lSeq;
        }
        
        if ( nodeOverlap.first->ends_[0] < params.readLen && nodeEnd > 0 )
        {
            int32_t nodeBgn = max( 0, nodeOverlap.first->ends_[0] ) - nodeOverlap.first->ends_[0];
            int32_t seqLen = min( nodeEnd, params.readLen ) - max( 0, nodeOverlap.first->ends_[0] );
            origin += nodeOverlap.first->seq_.substr( nodeBgn, seqLen );
            pSeq = &origin;
        }
        
        if ( nodeEnd > params.readLen )
        {
            int32_t nodeBgn = max( params.readLen, nodeOverlap.first->ends_[0] ) - nodeOverlap.first->ends_[0];
            int32_t seqLen = nodeEnd - max( params.readLen, nodeOverlap.first->ends_[0] );
            rSeq += nodeOverlap.first->seq_.substr( nodeBgn, seqLen );
            pSeq = &rSeq;
        }
        
        if ( nodeOverlap.second < 0 )
        {
            *pSeq += string( abs( nodeOverlap.second ), 'N' );
        }
    }
}

vector< pair<Node*, int> > Locus::getExportPath()
{
    vector< pair<Node*, int> > nodePath;
    
    // Compile left path
    for ( auto it = paths_[0][0].path.rbegin(); it != paths_[0][0].path.rend() - 1; it++ )
    {
        nodePath.push_back( make_pair( *it, (*it)->getOverlap( *(it+1), true ) ) );
    }
    
    // Bridge between left and right paths
    getExportPathBridge( nodePath, paths_[0][0].path[0], paths_[1][0].path[0] );
    
    // Compile right path
    for ( auto it = paths_[1][0].path.begin(); it != paths_[1][0].path.end() - 1; it++ )
    {
        nodePath.push_back( make_pair( *it, (*it)->getOverlap( *(it+1), true ) ) );
    }
    nodePath.push_back( make_pair( paths_[1][0].path.back(), 0 ) );
    
    blankMismatchedEnds( nodePath );
    
    return nodePath;
}

void Locus::getExportPathBridge( vector< pair<Node*, int> > &nodePath, Node* lft, Node* rght )
{
    if ( lft == rght )
    {
        return;
    }
    
    for ( Edge &e : lft->edges_[1] )
    {
        if ( e.node == rght )
        {
            nodePath.push_back( make_pair( lft, e.overlap ) );
            return;
        }
    }
    
    NodeSet originSet = lft->getDrxnNodes( true, true );
    NodeSet bridgeSet;
    rght->getDrxnNodesInSet( bridgeSet, originSet, false );
    
    assert( !bridgeSet.empty() );
    
    while ( lft && lft != rght )
    {
        Node* nxt = NULL;
        int nxtHits = -1;
        int nxtOverlap = 0;
        for ( Edge &e : lft->edges_[1] )
        {
            int edgeHits = e.node->getPairHitsTotal();
            if ( bridgeSet.find( e.node ) != bridgeSet.end() && edgeHits > nxtHits )
            {
                nxtHits = edgeHits;
                nxt = e.node;
                nxtOverlap = e.overlap;
            }
        }
        nodePath.push_back( make_pair( lft, nxtOverlap ) );
        lft = nxt;
    }
}

