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

#include "node.h"
#include <algorithm>
#include "shared_functions.h"

void Node::addSeed( NodeRoll& nodes, ReadStruct& read )
{
    // Redundant read, no need to append
    bool added = false;
    for ( Node* node : nodes.nodes ) if ( node->addSeed( read ) ) added = true;
    if ( added ) return;
    
    // Check posssible overlaps and splits
    NodeIntList ols;
    NodeIntIntList splits;
    int nodeCount = nodes.size();
    for ( int i = 0; i < nodeCount; i++ ) nodes[i]->addSeed( read, nodes, ols, splits );
    
    // Remove unfavourable overlaps and splits
    NodeSet fwdSet;
    for ( int i = 0; i < ols.size(); i++ ) ols[i].first->getDrxnNodes( 0 );
    for ( int i = 0; i < splits.size(); i++ ) get<0>( splits[i] )->getDrxnNodes( 0 );
    for ( int i = 0; i < ols.size(); i++ ) if ( fwdSet.find( ols[i].first ) != fwdSet.end() ) ols.erase( ols.begin() + i-- );
    for ( int i = 0; i < splits.size(); i++ ) if ( fwdSet.find( get<0>( splits[i] ) ) != fwdSet.end() ) splits.erase( splits.begin() + i-- );
    
    int32_t maxOl = 0;
    for ( int i = 0; i < ols.size(); i++ ) maxOl = max( maxOl, ols[i].second );
    for ( int i = 0; i < splits.size(); i++ ) maxOl = max( maxOl, get<1>( splits[i] ) );
    for ( int i = 0; i < ols.size(); i++ ) if ( ols[i].second < min( maxOl, ols[i].first->getBestOverlap( 1 ) ) ) ols.erase( ols.begin() + i-- );
    
    // Do splits
    for ( int i = 0; i < splits.size(); i++ )
    {
        if ( get<1>( splits[i] ) < maxOl ) continue;
        get<0>( splits[i] )->splitNode( get<2>( splits[i] ), nodes, 1 );
        ols.push_back( make_pair( get<0>( splits[i] ), get<1>( splits[i] ) ) );
    }
    
    if ( ols.size() == 1 && ols[0].first->edges_[1].empty() )
    {
        ols[0].first->addSeed( read, ols[0].second );
        return;
    }
    
    Node* node = new Node( read );
    nodes += node;
    
    for ( int i = 0; i < ols.size(); i++ ) node->addEdge( ols[i].first, ols[i].second, 0 );
}

bool Node::addSeed( ReadStruct &read )
{
    size_t it = seq_.find( read.seq );
    if ( it == string::npos ) return false;
    
    int coords[2]{ int( ends_[0]+it ), int( ends_[0]+it+read.seq.size() ) };
    add( read.id, coords[0], coords[1], isRedundant( coords[0], coords[1] ) );
    readTest();
    ends_.init( read.tether[0], 0 );
    ends_.init( read.tether[1], 1 );
    
    return true;
}

void Node::addSeed( ReadStruct &read, int ol )
{
    int coords[2]{ ends_[1] - ol, int( ends_[1] - ol + read.seq.size() ) };
    string seq = read.seq.substr( ol );
    appendSeq( seq, 1 );
    add( read.id, coords[0], coords[1], false );
    readTest();
    ends_.init( read.tether[0], 0 );
    ends_.init( read.tether[1], 1 );
    setOrigin();
}

void Node::addSeed( ReadStruct& read, NodeRoll& nodes, NodeIntList& ols, NodeIntIntList& splits )
{
    if ( int ol = mapSeqOverlap( seq_, read.seq, params.readLen / 2 ) )
    {
        ols.push_back( make_pair( this, ol ) );
        return;
    }
    
    int32_t coords[2], splitCoords[2];
    if ( !mapSeqEnd( read.seq, seq_, params.readLen / 2, coords, 0 ) ) return;
    coords[0] += ends_[0];
    coords[1] += ends_[0];
    if ( !getSplitCoords( splitCoords, coords[1], 0 ) ) return;
    coords[1] = splitCoords[0];
    
    if ( coords[1] - coords[0] > splitCoords[0] - splitCoords[1] )
    {
        for ( Node* node : splitNode( splitCoords[1], nodes, 1 ) ) node->drxn_ = 2;
        ols.push_back( make_pair( this, coords[1] - coords[0] ) );
        assert( false );
        return;
    }
    
    splits.push_back( make_tuple( this, coords[1] - coords[0], splitCoords[1] ) );
}

void Node::seedNode( Querier& bwt, NodeRoll& nodes, vector<Node*>& added, vector<Node*>& seeded, string& seq, ReadId id, int32_t coord, bool drxn )
{
    for ( Node* node : nodes.nodes ) if ( node->add( id, seq ) ) added.push_back( node );
    if ( !added.empty() ) return;
    
    Node* seed = new Node( seq, id, coord, drxn, true );
    nodes += seed;
    
    for ( Node* node : nodes.nodes ) for ( auto& read : node->reads_ ) if ( read.second.redundant ) assert( read.second[1] - read.second[0] < params.readLen );
    bool redundant = false, good = false;
    for ( int d : { 0, 1 } )
    {
        QueryJunction qj = bwt.mapJunction( seq, d );
        seed->stop( qj.failure_, d );
        if ( qj.failure_ ) continue;
        seed->addAlts( qj.alts_, nodes, d, drxn );
        for ( QueryNode* qn : qj.alts_) if ( qn->seq.find( seq ) != string::npos ) redundant = true;
        if ( redundant )
        {
            assert( !d || ( seed->edges_[0].empty() && seed->reads_.size() == 1 ) );
            nodes.erase( seed );
            for ( Node* node : nodes.nodes ) if ( node->add( id, seq ) ) seeded.push_back( node );
            for ( Node* node : seeded ) node->extendNode( bwt, nodes, d );
            for ( Node* node : seeded ) node->extendNode( bwt, nodes, !d );
            assert( !seeded.empty() );
            return;
        }
        seed->addExtensions( qj.nodes_, nodes, d, drxn );
        seed->extendNode( bwt, nodes, d );
        good = true;
    }
    
    if ( !good ) nodes.erase( seed );
    else seeded.push_back( seed );
    for ( Node* node : nodes.nodes ) for ( auto& read : node->reads_ ) if ( read.second.redundant ) assert( read.second[1] - read.second[0] < params.readLen );
    for ( Node* node : seeded ) node->remap( bwt );
    for ( Node* node : nodes.nodes ) for ( auto& read : node->reads_ ) if ( read.second.redundant ) assert( read.second[1] - read.second[0] < params.readLen );
}

//void Node::seedNode( Querier& bwt, NodeRoll& nodes, ReadId id, int32_t coord, int drxn )
//{
//    string seq = bwt.getSequence( id );
//    
//    vector<Node*> seeded;
//    for ( Node* node : nodes.nodes ) if ( node->add( id, seq ) ) seeded.push_back( node );
//    assert( seeded.empty() );
//    
//    Node* seed = new Node( seq, id, coord, drxn, true );
//    nodes += seed;
//    for ( int d : { 0, 1 } ) seed->extendNode( bwt, nodes, d );
//    seed->extendFork( bwt, nodes, 300, 8, 0 );
//}

vector<Node*> Node::seedNode( Querier& bwt, NodeRoll& nodes, string s, int lOl, int rOl, int drxn, int32_t coord )
{
    MapNode mn;
    mn.seq = s;
    bwt.mapSequence( mn.seq, mn.ids, mn.coords );
//    mn.recoil();
    mn.setRedundant();
    
    int len = 0, mark = 0;
    bool good = false;
    struct Exts
    {
        Exts( Node* node, Coords* coords, int i, int j )
        : node( node ), created( !coords ), ended( !coords )
        {
            mapped[0][0] = i;
            mapped[0][1] = j;
            mapped[1][0] = coords ? coords->coords[0] : node->ends_[0];
            mapped[1][1] = coords ? coords->coords[1] : node->ends_[1];
        };
        Node* node;
        vector<Node*> edges[2];
        int32_t mapped[2][2];
        bool created, ended;
    };
    vector<Exts> exts;
    Coords* coords;
    for ( int i = 0; i <= mn.ids.size(); i++ ) if ( i == mn.ids.size() || len < mn.coords[1][i] )
    {
        vector<Exts> exted;
        bool bad = i == mn.ids.size();
        if ( !bad ) len = mn.coords[1][i];
        if ( !bad ) for ( Node* node : nodes.nodes ) if ( coords = node->getRead( mn.ids[i] ) )
        {
            bool added = false;
            for ( Exts& ext : exts ) if ( !ext.ended && node->isClone( ext.node ) && ( added = true ) )
            {
                ext.mapped[0][1] = mn.coords[1][i];
                ext.mapped[1][1] = coords->coords[1];
            }
            if ( !added ) exted.push_back( Exts( node, coords, mn.coords[0][i], mn.coords[1][i] ) );
            bad = true;
        }
        for ( Exts& ext : exts ) if ( ext.mapped[0][1] != len ) ext.ended = true;
        if ( good && bad )
        {
            Node* seed = new Node( mn.seq, coord, drxn );
            for ( int j = mark; j < i; j++ ) seed->add( mn.ids[j], mn.coords[0][j], mn.coords[1][j], mn.redundant[j] );
            seed->recoil();
            nodes += seed;
            exts.push_back( Exts( seed, NULL, mn.coords[0][mark], mn.coords[0][mark]+seed->size() ) );
        }
        if ( good == bad ) mark = i;
        good = !bad;
        exts.insert( exts.end(), exted.begin(), exted.end() );
    }
    
    lOl -= mn.coords[0][0];
    rOl -= mn.seq.size() - len;
    
    for ( Exts& ext : exts ) if ( !ext.created ) for ( int d : { 0, 1 } ) if ( ext.mapped[1][!d] != ext.node->ends_[!d] ) ext.node->splitNode( ext.mapped[1][!d], nodes, d )[0];;
    for ( int i = 0; i < exts.size(); i++ ) for ( int j = i+1; j < exts.size(); j++ ) for ( int d : { 0, 1 } ) assert( exts[i].mapped[0][d] < exts[j].mapped[0][d] );
//    {
//        for ( int j = i+1; j < exts.size(); j++ ) if ( exts[i].node->isClone( exts[j].node ) ) exts.erase( exts.begin() + j-- );
//        for ( int j = i+1; j < exts.size(); j++ ) for ( int d : { 0, 1 } ) assert( exts[i].mapped[0][d] < exts[j].mapped[0][d] );
//        for ( int d : { 0, 1 } ) if ( exts[i].mapped[1][!d] != exts[i].node->ends_[!d] ) exts[i].node = exts[i].node->splitNode( exts[i].mapped[1][!d], nodes, d )[0];
//        bool clone = !exts[i].node->bad_ && ( lOl > 0 || i ) && ( rOl > 0 || i+1 != exts.size() );
//        if ( clone ) exts[i].node = new Node( exts[i].node, nodes, drxn, false );
//    }
    for ( int i = 1; i < exts.size(); i++ ) exts[i-1].node->addEdge( exts[i].node, exts[i-1].mapped[0][1] - exts[i].mapped[0][0], 1, false );
    
    if ( !exts.empty() ) for ( int d : { 0, 1 } ) if ( ( d ? rOl : lOl ) > 0 )
    {
        int32_t coords[2];
        Node* base = ( d ? exts.back() : exts[0] ).node,* edge;
        Nodes block;
        for ( Node* clone : base->clones() ) for ( Edge& e : clone->edges_[d] ) block += e.node;
        for ( int i = 0; i < nodes.size(); i++ ) if ( ( edge = nodes[i] ) && !base->isClone( edge ) && mapSeqEnd( base->seq_, edge->seq_, d ? rOl : lOl, coords, d ) )
        {
            if ( !block.add( edge ) ) continue;
            int32_t cut = coords[!d] + nodes[i]->ends_[0];
            if ( !edge->getNextReadCoord( cut, !d, d ) ) continue;
            if ( d ? edge->ends_[0] < cut : cut < edge->ends_[1] ) edge = edge->splitNode( cut, nodes, d )[0];
            base->addEdge( edge, coords[1] - coords[0], d, false );
            block += edge;
        }
    }
    
    for ( int i = 0; i < exts.size(); i++ ) ( drxn ? exts[i] : exts.end()[-i-1] ).node->offset( NULL, drxn );
    if ( drxn == 2 ) for ( int i = 0; i < exts.size(); i++ ) exts[i].node->setOrigin();
    
    vector<Node*> path;
    for ( Exts& ext : exts )
    {
        for ( int d : { 0 , 1 } ) ext.node->extendNode( bwt, nodes, d );
        if ( !ext.node->bad_ ) ext.node->setVerified();
        path.push_back( ext.node );
    }
    
    return path;
    
//    bool seeded = false;
//    Nodes seeds;
//    for ( int i = 0; i < mn.ids.size(); i++ ) if ( !mn.redundant[i] ) for ( Node* node : nodes.nodes ) if ( node->getRead( mn.ids[i] ) ) seeds += node;
//    if ( !seeds.empty() )
//    {
//        for ( Node* node : seeds.nodes ) for ( int d : { 0, 1 } ) if ( node->edges_[d].empty() )
//        {
//            node->stop( 0, d );
//            node->extendNode( bwt, nodes, d );
//            if ( !node->edges_[d].empty() ) seeded = true;
//        }
//        if ( seeded ) for ( Node* node : seeds.nodes ) if ( !node->bad_ ) node->setVerified();
//        if ( seeded ) return;
//    }
//    
//    Node* seed = new Node( mn.seq );
//    nodes += seed;
//    seed->drxn_ = drxn;
//    seed->ends_[0] = 0;
//    seed->ends_[1] = seed->seq_.size();
//    if ( drxn == 2 ) seed->setOrigin();
//    for ( int i = 0; i < mn.ids.size(); i++ ) seed->add( mn.ids[i], mn.coords[0][i], mn.coords[1][i], mn.redundant[i] );
//    for ( int d : { 0, 1 } ) seed->extendNode( bwt, nodes, d );
//    if ( !seed->bad_ ) seed->setVerified();
}

void Node::trimSeed( Querier &bwt, NodeRoll &nodes )
{
    NodeList tests[2];
    int goods[2]{0};
    Nodes good[2];
    for ( Node* node : nodes.nodes )
    {
        int readCount = node->countReads( true );
        for ( int d : { 0, 1 } )
        {
            if ( readCount > 2 ) good[d].fill( node, !d, true );
            else if ( node->edges_[d].empty() )
            {
                int minOl = max( 1 + params.readLen / 2, 10 + (int)node->seq_.size() - node->getBestOverlap( !d ) );
                if ( bwt.isExtendable( node->seq_, minOl, d ) ) good[d].fill( node, !d, true );
            }
            else for ( Edge& e : node->edges_[!d] )
            {
                if ( readCount + e.node->countReads( true ) ) good[d].fill( e.node, !d, true );
                else for ( Edge& fe : e.node->edges_[!d] ) good[d].fill( fe.node, !d, true );
            }
        }
    }
    
    assert( !good[0].empty() && !good[1].empty() );
    
    for ( int i = nodes.size(); i-- > 0; ) if ( !good[0].find( nodes[i] ) || !good[1].find( nodes[i] ) ) nodes.erase( nodes[i] );
    
//    for ( Node* node : nodes.nodes )
//    {
//        int readCount = node->countReads( true );
//        for ( int i = 0; i < 2; i++ )
//        {
//            if ( !node->edges_[i].empty() ) continue;
//            if ( readCount > 2 ) goods[i]++;
//            else tests[i].push_back( node );
//        }
//        if ( node->edges_[0].empty() || node->edges_[1].empty() || readCount > 1 ) continue;
//        
//        // NYI: remove a weak and redundant node
//        NodeSet selfSet = { node }, lSet, rSet;
//        for ( Edge &e : node->edges_[0] ) e.node->getDrxnNodesNotInSet( lSet, selfSet, 1 );
//        for ( Edge &e : node->edges_[1] ) if ( lSet.find( e.node ) != lSet.end() ) assert( false );
//    }
//    
//    for ( int i = 0; i < 2; i++ )
//    {
//        for ( Node* node : tests[i] ) if ( !nodes.find( node ) && node->edges_[i].empty() )
//        {
//            bwt.isExtendable()
//            bool noMatches = false;
//            int minOl = max( 1 + params.readLen / 2, 10 + (int)node->seq_.size() - node->getBestOverlap( !i ) );
//            bwt.mapExtensions( noMatches, node->seq_, i, minOl );
//            if ( noMatches ) nodes.erase( node );
//            else assert( false );
//        }
//    }
    
    Node::merge( nodes );
    for ( Node* node : nodes.nodes ) node->setOrigin();
}
