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

bool Node::isSeed( int32_t seedLen )
{
    return validLimits_[0] < ends_[1] && 0 < ends_[1] && ends_[0] < validLimits_[3] && ends_[0] < seedLen;
}

void Node::seedAdd( ReadStruct &read )
{
    int32_t diff = read.coords[1] - ends_[1];
    if ( diff >= 0 )
    {
        validLimits_[2] = validLimits_[3] = read.tether[1];
        seq_ += read.seq.substr( read.seq.length() - diff );
        ends_[1] += diff;
    }
    else
    {
        assert( false );
    }
    reads_.insert( make_pair( read.readId, Coords( read.coords[0], read.coords[1], diff < 0 ) ) );
}

bool Node::seedCongruent( ReadStruct &read, int32_t &coord )
{
    bool congruent = read.tether[0] == read.coords[0] && validLimits_[0] <= read.tether[0] && read.tether[0] < validLimits_[3];
    bool maybe = read.tether[0] <= validLimits_[0];
    coord = min( read.tether[1], validLimits_[3] );
    
    if ( !congruent && maybe )
    {
        congruent = true;
        int32_t diff = read.coords[0] - ends_[0];
        int32_t len = max( read.tether[0], validLimits_[0] ) - read.coords[0];
        for ( int i( 0 ); i < len; i++ )
        {
            congruent = congruent && seq_[i+diff] == read.seq[i];
        }
    }
    
    while( congruent && coord != ends_[1] && seq_[coord - ends_[0]] == read.seq[coord - read.coords[0]] )
    {
        coord++;
    }
    
    if ( congruent && coord != ends_[1] )
    {
        congruent = false;
        for ( auto &read : reads_ )
        {
            congruent = congruent || read.second[1] <= coord;
        }
    }
    
    return congruent;
}

void Node::seedGetExtend( NodeList* extendNodes, NodeSet &seedSet, NodeSet &delSet, int32_t* limits )
{
    for ( int drxn : { 0, 1 } )
    {
        NodeSet currSet = seedSet;
        NodeIntMap edgeMap;
        
        while( !currSet.empty() )
        {
            NodeSet nxtSet;
            NodeSet currFwd;
            for ( Node* curr : currSet )
            {
                curr->getDrxnNodes( currFwd, drxn );
            }
            
            for ( Node* curr : currSet )
            {
                if ( currFwd.find( curr ) != currFwd.end() )
                {
                    nxtSet.insert( curr );
                }
                else
                {
                    int edgeCount = 8 + curr->edges_[!drxn].size();

                    for ( Node* bck : curr->getNextNodes( !drxn ) )
                    {
                        auto it = edgeMap.find( bck );
                        edgeCount = it != edgeMap.end() ? min( edgeCount, it->second ) : edgeCount;
                    }
                    
                    edgeCount = seedSet.find( curr ) != seedSet.end() ? curr->edges_[!drxn].size() : edgeCount;
                    edgeCount += curr->edges_[drxn].size() - curr->edges_[!drxn].size();
                    
                    if ( edgeCount < 8 && ( drxn ? curr->ends_[1] < limits[1] : limits[0] < curr->ends_[0] ) )
                    {
                        if ( curr->isContinue( drxn ) && !curr->clones_ )
                        {
                            extendNodes[drxn].push_back( curr );
                        }
                        curr->getNextNodes( nxtSet, drxn );
                    }

                    edgeMap[curr] = edgeCount;
                }
            }
            currSet = nxtSet;
        }
    }
}

bool Node::seedJoin( Node* node, int32_t coord, bool drxn )
{
//    for ( bool d : { 0, 1 } )
//    {
//        bool anyOrigin = false;
//        NodeSet bckSet = ( d == drxn ? this : node )->getDrxnNodes( !d, false, true );
//        for ( Node* bck : bckSet )
//        {
//            anyOrigin = anyOrigin || bck->drxn_ == 2;
//        }
//        if ( !anyOrigin )
//        {
//            int x = 0;
//            for ( Node* bck : bckSet )
//            {
//                bck->drxn_ = !d;
//            }
//            assert( false );
//            return;
//        }
//    }
    
    NodeSet bckSet = node->getDrxnNodes( !drxn );
    for ( Node* bck : getDrxnNodes( !drxn ) )
    {
        assert( bckSet.find( bck ) == bckSet.end() );
    }
    
    if ( drxn ? ends_[1] <= coord + 100 : coord - 100 <= ends_[0] )
    {
        Node* fwdNode = ( drxn ? node : this );
        
        NodeSet fwdSet = fwdNode->getDrxnNodes( 1, false, true );
        NodeSet revSet;

        for ( Node* fwd : fwdSet )
        {
            fwd->drxn_ = 1;
        }

        for ( Node* connected : fwdNode->getConnectedNodes( false ) )
        {
            if ( connected->drxn_ == 2 )
            {
                revSet.insert( connected );
                connected->getDrxnNodes( revSet, 0 );
            }
        }

        for ( Node* fwd : fwdNode->getDrxnNodes( 1, false, true ) )
        {
            fwdSet.insert( fwd );
            fwd->getDrxnNodesNotInSet( fwdSet, revSet, 0 );
        }

        for ( Node* fwd : fwdSet )
        {
            fwd->drxn_ = 1;
        }
        
        return true;
    }
    
    return false;
}

void Node::seedJoinLoci( Node** nodes )
{
    NodeSet fwdSets[2];
    NodeSet revSets[2];
    for ( bool drxn : { 0, 1 } )
    {
        fwdSets[drxn] = nodes[drxn]->getDrxnNodes( drxn, false, true );
    }
    
    for ( bool drxn : { 0, 1 } )
    {
        for ( Node* node : fwdSets[drxn] )
        {
            for ( Node* prv : node->getNextNodes( !drxn ) )
            {
                if ( fwdSets[drxn].find( prv ) == fwdSets[drxn].end()
                        && revSets[drxn].find( prv ) == revSets[drxn].end() )
                {
                    revSets[drxn].insert( prv );
                    prv->getDrxnNodesNotInSet( revSets[drxn], fwdSets[!drxn], !drxn );
                }
            }
        }
    }
    
    assert( false );
}

void Node::seedSetDrxnNodes( Node* fork, NodeList &nodes, bool drxn )
{
    for ( Node* fwd : fork->getDrxnNodes( drxn ) )
    {
        fwd->drxn_ = drxn;
        fwd->validLimits_[0] = fwd->validLimits_[1];
        fwd->validLimits_[3] = fwd->validLimits_[2];
        if ( find ( nodes.begin(), nodes.end(), fwd ) == nodes.end() )
        {
            nodes.push_back( fwd );
        }
    }
    vector<Edge> edges;
    for ( Edge &e : fork->edges_[drxn] )
    {
        edges.push_back( e );
        e.node->removeEdge( fork, !drxn );
    }
    fork->edges_[drxn].clear();
    
    for ( Edge &e : edges )
    {
        addEdge( e.node, e.overlap, drxn );
    }
}

Node* Node::seedSetOrigin( NodeList &forkList )
{
    Node* node = new Node();
    node->drxn_ = 2;
    node->ends_[0] = forkList[0]->ends_[0];
    node->ends_[1] = forkList[0]->ends_[1];
    node->seq_ = forkList[0]->seq_;
    node->validLimits_[0] = forkList[0]->validLimits_[1];
    node->validLimits_[1] = forkList[0]->validLimits_[1];
    node->validLimits_[2] = forkList[0]->validLimits_[2];
    node->validLimits_[3] = forkList[0]->validLimits_[2];
    
    int readCount[forkList.size()];
    readCount[0] = forkList[0]->reads_.size();
    
    for ( auto &read : forkList[0]->reads_ )
    {
        node->reads_.insert( read );
    }
    
    int32_t offset = 0;
    for ( int i( 1 ); i < forkList.size(); i++ )
    {
        int32_t overlap = forkList[i-1]->getOverlap( forkList[i], 1 );
        offset += forkList[i]->ends_[0] - forkList[i-1]->ends_[1] + overlap;
        node->seq_ += forkList[i]->seq_.substr( overlap );
        for ( auto read : forkList[i]->reads_ )
        {
            read.second.offset( -offset );
            node->reads_.insert( read );
        }
        node->ends_[1] = forkList[i]->ends_[1] - offset;
        node->validLimits_[0] = min( node->validLimits_[0], forkList[i]->validLimits_[1] - offset );
        node->validLimits_[1] = min( node->validLimits_[1], forkList[i]->validLimits_[1] - offset );
        node->validLimits_[2] = max( node->validLimits_[2], forkList[i]->validLimits_[2] - offset );
        node->validLimits_[3] = max( node->validLimits_[3], forkList[i]->validLimits_[2] - offset );
        readCount[i] = forkList[i]->reads_.size();
    }
    
    int seqLen = node->seq_.length();
    assert( node->ends_[1] - node->ends_[0] == node->seq_.length() );
    node->setCoverage();
    int x = 0;
    
    return node;
}

void Node::seedSplit( NodeList &nodes, int32_t coord )
{
    Node* node = new Node();
    node->ends_[1] = node->ends_[0] = ends_[1];
    ends_[1] = ends_[0];
    for ( auto it = reads_.begin(); it != reads_.end(); )
    {
        if ( it->second[1] <= coord )
        {
            ends_[1] = max( ends_[1], it->second[1] );
            it++;
        }
        else
        {
            node->ends_[0] = min( node->ends_[0], it->second[0] );
            node->reads_.insert( *it );
            it = reads_.erase( it );
        }
    }
    
    node->validLimits_[2] = node->validLimits_[3] = max( validLimits_[3], node->ends_[0] );
    node->validLimits_[0] = node->validLimits_[1] = max( validLimits_[0], node->ends_[0] );
    validLimits_[2] = validLimits_[3] = min( validLimits_[3], ends_[1] );
    
    node->seq_ = seq_.substr( seq_.length() - ( node->ends_[1] - node->ends_[0] ) );
    seq_.erase( seq_.begin() + ( ends_[1] - ends_[0] ), seq_.end() );
    for ( Edge &e : edges_[1] )
    {
        node->addEdge( e.node, e.overlap, 1, false );
        e.node->removeEdge( this, 0 );
    }
    edges_[1].clear();
    node->addEdge( this, 0 );
    nodes.push_back( node );
}

void Node::seedValidate( NodeSet &seedSet, NodeSet &delSet, int32_t* validLimits, int32_t* ends, bool doDel )
{
    NodeList notValid[2] = { { seedSet.begin(), seedSet.end() }, { seedSet.begin(), seedSet.end() } };
    
    for ( bool drxn : { 0, 1 } )
    {
        NodeSet currSet = { notValid[drxn].begin(), notValid[drxn].end() };
        notValid[drxn].clear();
        while ( !currSet.empty() )
        {
            NodeSet nxtSet, currFwd;
            for ( Node* curr : currSet )
            {
                curr->getDrxnNodes( currFwd, drxn );
            }

            for ( Node* curr : currSet )
            {
                if ( curr->isDeadEnd( drxn ) && doDel )
                {
                    curr->dismantleNode( delSet, drxn );
                    seedSet.erase( curr );
                }
                else if ( currFwd.find( curr ) != currFwd.end() )
                {
                    nxtSet.insert( curr );
                }
                else if ( curr->validate( drxn ) )
                {
                    curr->getNextNodes( nxtSet, drxn );
                }
                else
                {
                    notValid[drxn].push_back( curr );
                }

                validLimits[drxn] = ( drxn ? max( validLimits[1], curr->validLimits_[2] )
                                      : min( validLimits[0], curr->validLimits_[1] ) );
            }

            currSet = nxtSet;
        }
    }
    
    if ( ends[1] - ends[0] > params.maxPeMean && doDel )
    {
        int cutoff = (params.maxPeMean * params.cover * 2 ) / params.readLen;
        for ( bool drxn : { 0, 1 } )
        {
            for ( Node* node : notValid[drxn] )
            {
                int misses = node->getEndMarks( drxn ) * 5;
                int hits = 2;
                int32_t endCoord = node->ends_[drxn];
                for ( Node* fwd : node->getDrxnNodes( drxn ) )
                {
                    endCoord = drxn ? max( endCoord, fwd->ends_[1] )
                                    : min( endCoord, fwd->ends_[0] );
                    misses += fwd->reads_.size();
                    hits += fwd->getPairHitsTotal();
                }
                
                if ( misses / hits > cutoff && hits < 4 )
                {
                    node->dismantleNode( delSet, drxn );
                }
            }
        }
    }
}

bool Node::seedValidate( bool drxn )
{
    if ( validated_ )
    {
        return true;
    }
    
    NodeList tNodes = this->getTargetNodes( drxn, true );
    NodeOffsetMap fwdMap = getDrxnNodesOffset( drxn, 0, true );
    NodeOffsetMap revMap = getDrxnNodesOffset( !drxn, 0, true );
    
    while ( seedValidate( tNodes, fwdMap, revMap, drxn ) );
    for ( Node* fwd : getDrxnNodes( drxn ) )
    {
        fwd->seedValidate( tNodes, fwdMap, revMap, drxn );
    }
    
    validated_ = validated_ || ( !isContinue( 0 ) && !isContinue( 1 ) && validLimits_[1] == ends_[0] && validLimits_[2] == ends_[1] );
    
    return ( drxn ? validLimits_[2] == ends_[1] : validLimits_[1] == ends_[0] );
}

bool Node::seedValidate( NodeList &tNodes, NodeOffsetMap &fwdMap, NodeOffsetMap &revMap, bool drxn )
{
    bool didPair = false;
    unordered_set<SeqNum> usedIds;
    NodeSet tSet = getDrxnNodes( !drxn ), hitSet;
    int x = 0;

    for ( ReadMark &mark : marks_[drxn] )
    {
        if ( drxn ? validLimits_[1] < mark.mark
                  : mark.mark < validLimits_[2] )
        {
            for ( Node* t : tNodes )
            {
                auto it = t->reads_.find( mark.readId );
                if ( it != t->reads_.end()
                        && ( drxn ? it->second[1] <= t->validLimits_[2] : validLimits_[1] <= it->second[0] ) )
                {
                    didPair = true;
                    hitSet.insert( t );
                    usedIds.insert( mark.readId );
                    pushValidLimits( mark.mark, drxn );
                    t->pushValidLimits( it->second[!drxn], !drxn );
                    NodeSet midSet = t->getDrxnNodesInSet( tSet, drxn );
                    auto r = pairs_.insert( make_pair( t, 1 ) );
                    if ( !r.second ) r.first->second++;
                    if ( this != t )
                    {
                        r = t->pairs_.insert( make_pair( this, 1 ) );
                        if ( !r.second ) r.first->second++;
                        pushValidLimits( ends_[!drxn], !drxn );
                        t->pushValidLimits( t->ends_[drxn], drxn );
                        for ( Node* node : midSet )
                        {
                            node->pushValidLimits( node->ends_[0], 0 );
                            node->pushValidLimits( node->ends_[1], 1 );
                        }
                    }
                }
            }
        }
    }

    for ( Node* node : hitSet )
    {
        node->removeMarks( usedIds, false, true, !drxn );
    }
    removeMarks( usedIds, false, false, drxn );
    
    return didPair;
}
