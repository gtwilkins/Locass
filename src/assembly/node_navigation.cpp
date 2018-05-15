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

#include <algorithm>

#include "node.h"

NodeSet Node::getBetweenNodes( Node* node, bool drxn )
{
    NodeSet fwdSet = getDrxnNodes( drxn );
    NodeSet betweenSet = node->getDrxnNodesInSet( fwdSet, !drxn );
    return betweenSet;
}

NodeSet Node::getConnectedNodes( bool sameGraph )
{
    NodeSet nodes;
    getConnectedNodes( nodes, sameGraph );
    return nodes;
}

void Node::getConnectedNodes( NodeSet &nodes, bool sameGraph )
{
    nodes.insert( this );
    for ( bool drxn : { 0, 1 } )
    {
        for ( Node* nxt : getNextNodes( drxn ) )
        {
            if ( nodes.find( nxt ) == nodes.end() && ( !sameGraph || nxt->drxn_ == drxn_ ) )
            {
                nxt->getConnectedNodes( nodes, sameGraph );
            }
        }
    }
}

NodeSet Node::getDrxnNodes( bool drxn, bool sameGraph, bool inclSelf )
{
    NodeSet nodes;
    if ( inclSelf )
    {
        nodes.insert( this );
    }
    this->getDrxnNodes( nodes, drxn, sameGraph );
    return nodes;
}

void Node::getDrxnNodes( NodeSet &nodes, bool drxn, bool sameGraph )
{
    for ( Edge &edge : edges_[drxn] )
    {
        if ( nodes.find( edge.node ) == nodes.end() && ( !sameGraph || edge.node->drxn_ == drxn_ ) )
        {
            nodes.insert( edge.node );
            edge.node->getDrxnNodes( nodes, drxn, sameGraph );
        }
    }
}

void Node::getDrxnNodes( NodeSet &nodes, bool drxn, int32_t limit )
{
    if ( ( drxn ? ends_[1] < limit : ends_[0] > limit ) )
    {
        for ( Edge &edge : edges_[drxn] )
        {
            if ( nodes.find( edge.node ) == nodes.end() )
            {
                nodes.insert( edge.node );
                edge.node->getDrxnNodes( nodes, drxn, limit );
            }
        }
    }
}

NodeSet Node::getDrxnNodesInSet( NodeSet &inSet, bool drxn, bool inclSelf )
{
    NodeSet nodes;
    if ( inclSelf )
    {
        nodes.insert( this );
    }
    getDrxnNodesInSet( nodes, inSet, drxn );
    return nodes;
}

void Node::getDrxnNodesInSet( NodeSet &nodes, NodeSet &inSet, bool drxn )
{
    for ( Edge &edge : edges_[drxn] )
    {
        if ( nodes.find( edge.node ) == nodes.end() && inSet.find( edge.node ) != inSet.end() )
        {
            nodes.insert( edge.node );
            edge.node->getDrxnNodesInSet( nodes, inSet, drxn );
        }
    }
}

NodeSet Node::getDrxnNodesNotInSet( NodeSet &notSet, bool drxn, bool inclSelf )
{
    NodeSet nodes;
    if ( inclSelf )
    {
        nodes.insert( this );
    }
    getDrxnNodesNotInSet( nodes, notSet, drxn );
    return nodes;
}

void Node::getDrxnNodesNotInSet( NodeSet &nodes, NodeSet &notSet, bool drxn )
{
    for ( Edge &edge : edges_[drxn] )
    {
        if ( nodes.find( edge.node ) == nodes.end() && notSet.find( edge.node ) == notSet.end() )
        {
            nodes.insert( edge.node );
            edge.node->getDrxnNodesNotInSet( nodes, notSet, drxn );
        }
    }
}

NodeOffsetMap Node::getDrxnNodesOffset( bool drxn, int32_t limitDist, bool inclSelf )
{
    NodeOffsetMap nodes;
    if ( inclSelf )
    {
        nodes.insert( make_pair( this, make_pair( (int32_t)0, (int32_t)0 ) ) );
    }
    
    limitDist = ( limitDist > 0 
            ? ( drxn ? ends_[1] + limitDist : ends_[0] - limitDist ) 
            : ( drxn ? numeric_limits<int32_t>::max() : numeric_limits<int32_t>::min() ) );
    
    NodeSet fwdSet;
    getDrxnNodes( fwdSet, drxn, limitDist );
    NodeSet currSet = { this };
    bool didAdvance = false;
    while ( !currSet.empty() )
    {
        didAdvance = false;
        NodeSet nxtSet;
        
        for ( Node* curr : currSet )
        {
            for ( Edge &e : curr->edges_[!drxn] )
            {
                if ( fwdSet.find( e.node ) != fwdSet.end() && nodes.find( e.node ) == nodes.end() )
                {
                    nxtSet.insert( curr );
                    continue;
                }
            }
            didAdvance = true;
            int32_t currOff[2]{0};
            auto it = nodes.find( curr );
            if ( it != nodes.end() )
            {
                currOff[0] = it->second.first;
                currOff[1] = it->second.second;
            }
            for ( Edge &e : curr->edges_[drxn] )
            {
                if ( fwdSet.find( e.node ) == fwdSet.end() ) continue;
                int32_t offset = drxn ? e.node->ends_[0] + e.overlap - curr->ends_[1]
                                      : e.node->ends_[1] - e.overlap - curr->ends_[0];
                int32_t eOff[2] = { currOff[0] + offset, currOff[1] + offset };
                auto r = nodes.insert( make_pair( e.node, make_pair( eOff[0], eOff[1] ) ) );
                if ( !r.second )
                {
                    r.first->second.first = min( r.first->second.first, eOff[0] );
                    r.first->second.second = max( r.first->second.second, eOff[1] );
                }
                nxtSet.insert( e.node );
            }
        }
        
        currSet = nxtSet;
        assert( didAdvance );
    }
    
//    getDrxnNodesOffset( nodes, drxn, limitDist );
    
    return nodes;
}

void Node::getDrxnNodesOffset( NodeOffsetMap &nodes, bool drxn, int32_t &limit )
{
    if ( drxn ? ends_[1] < limit : limit < ends_[0] )
    {
        int32_t thisOffsets[2]{0};
        auto it = nodes.find( this );
        if ( it != nodes.end() )
        {
            thisOffsets[0] = it->second.first;
            thisOffsets[1] = it->second.second;
        }
        for ( Edge &e : edges_[drxn] )
        {
            int32_t offset = drxn 
                    ? e.node->ends_[0] + e.overlap - ends_[1]
                    : e.node->ends_[1] - e.overlap - ends_[0];
            int32_t offsets[2] = { offset + thisOffsets[0], offset + thisOffsets[1] };
            auto result = nodes.insert( make_pair( e.node, make_pair( offsets[0], offsets[1] ) ) );
            if ( !result.second )
            {
                if ( offsets[0] < result.first->second.first )
                {
                    result.first->second.first = offsets[0];
                    result.second = true;
                }
                if ( result.first->second.second < offsets[1] )
                {
                    result.first->second.second = offsets[1];
                    result.second = true;
                }
            }
            if ( result.second )
            {
                e.node->getDrxnNodesOffset( nodes, drxn, limit );
            }
        }
    }
}

NodeSet Node::getEndNodes( bool drxn, bool inclStopped )
{
    NodeSet ends, fwdSet = { this };
    getDrxnNodes( fwdSet, drxn );
    for ( Node* fwd : fwdSet )
    {
        if ( fwd->isContinue( drxn ) || ( inclStopped && fwd->edges_.empty() ) )
        {
            ends.insert( fwd );
        }
    }
    return ends;
}

NodeSetList Node::getNodeSetsExclusive( Node* a, Node* b, bool drxn )
{
    NodeSetList nodeSets( 2 );
    NodeSet sharedSet, fwdSet = a->getDrxnNodes( drxn, false, true );
    for ( Node* fwd : b->getDrxnNodes( drxn, false, true ) )
    {
        if ( fwdSet.find( fwd ) != fwdSet.end() )
        {
            sharedSet.insert( fwd );
        }
        else
        {
            nodeSets[1].insert( fwd );
        }
    }
    
    for ( Node* fwd : fwdSet )
    {
        if ( sharedSet.find( fwd ) == sharedSet.end() )
        {
            nodeSets[0].insert( fwd );
        }
    }
    
    return nodeSets;
}

NodeSet Node::getInvalidNodes( bool drxn )
{
    NodeSet nodes;
    if ( !isValidated() )
    {
        nodes.insert( this );
        getInvalidNodes( nodes, drxn );
    }
    return nodes;
}

void Node::getInvalidNodes( NodeSet &nodes, bool drxn )
{
    for ( Node* nxt : getNextNodes( drxn ) )
    {
        if ( !nxt->isValidated() )
        {
            nodes.insert( nxt );
            nxt->getInvalidNodes( nodes, drxn );
        }
    }
}

void Node::getNextNodes(NodeSet &nextNodes, bool drxn)
{
    for ( Edge edge: edges_[drxn] )
    {
        nextNodes.insert(edge.node);
    }
}

NodeSet Node::getNextNodes( bool drxn )
{
    NodeSet nextNodes;
    getNextNodes( nextNodes, drxn );
    return nextNodes;
}

void Node::getNextNodesInSet( NodeSet &nextSet, NodeSet &inSet, bool drxn )
{
    for ( Edge edge: edges_[drxn] )
    {
        if ( inSet.find( edge.node ) != inSet.end() )
        {
            nextSet.insert(edge.node);
        }
    }
}

NodeSet Node::getNextNodesInSet( NodeSet &inSet, bool drxn )
{
    NodeSet nextSet;
    getNextNodesInSet( nextSet, inSet, drxn );
    return nextSet;
}

NodeSet Node::getNotForwardSet( NodeSet &tmpCurrSet, bool drxn )
{
    NodeSet fwdSet, currSet;
    for ( Node* curr : tmpCurrSet )
    {
        curr->getDrxnNodes( fwdSet, drxn );
    }
    
    for ( Node* curr : tmpCurrSet )
    {
        if ( fwdSet.find( curr ) == fwdSet.end() )
        {
            currSet.insert( curr );
        }
    }
    
    return currSet;
}

NodeIntList Node::getOffsetEdges( bool drxn )
{
    NodeIntList offsets;
    for ( Edge &e : edges_[drxn] )
    {
        // Offset is of the node from its edge node
        int32_t offset = drxn ? ends_[1] - e.overlap - e.node->ends_[0]
                              : ends_[0] + e.overlap - e.node->ends_[1];
        if ( abs( offset ) > 200 )
        {
            offsets.push_back( make_pair( e.node, offset ) );
        }
    }
    return offsets;
}

void Node::getRedundantNodes( NodeSet &nodes, int32_t* coords, bool drxn )
{
    for ( Edge &e : edges_[drxn] )
    {
        int32_t offset = drxn ? ends_[1] - e.overlap - e.node->ends_[0]
                              : ends_[0] + e.overlap - e.node->ends_[1];
        int32_t offCoords[2] = { coords[0] + offset, coords[1] + offset };
        if ( e.node->ends_[0] <= coords[0] && coords[1] <= e.node->ends_[1] )
        {
            nodes.insert( e.node );
            e.node->getRedundantNodes( nodes, offCoords, drxn );
        }
    }
}

bool Node::isOffset( bool drxn )
{
    for ( int i( 0 ); i + 1 < edges_[drxn].size(); i++ )
    {
        int32_t off = drxn ? edges_[drxn][i].node->ends_[0] + edges_[drxn][i].overlap 
                           : edges_[drxn][i].node->ends_[1] - edges_[drxn][i].overlap;
        for ( int j( i + 1 ); j < edges_[drxn].size(); j++ )
        {
            int32_t off2 = drxn ? edges_[drxn][j].node->ends_[0] + edges_[drxn][j].overlap 
                                : edges_[drxn][j].node->ends_[1] - edges_[drxn][j].overlap;
            if ( abs( off - off2 ) > 100 )
            {
                return true;
            }
        }
    }
    return false;
}

bool Node::offset( int32_t off )
{
    if ( off != 0 )
    {
        assert( drxn_ != 2 );
        for ( auto &np : pairs_ )
        {
            np.first->resetFurthest( this );
        }
        
        ends_[0] += off;
        ends_[1] += off;
        for ( int i( 0 ); i < 4; i++ )
        {
            if ( validLimits_[i] != numeric_limits<int32_t>::min() && validLimits_[i] != numeric_limits<int32_t>::max() )
            {
                validLimits_[i] += off;
            }
        }
        for ( auto &read : reads_ )
        {
            read.second.offset( off );
        }
        for ( ReadMark &mark : marks_[0] )
        {
            mark.offset( off );
        }
        for ( ReadMark &mark : marks_[1] )
        {
            mark.offset( off );
        }
        return true;
    }
    return false;
}

void Node::offsetEdge( Edge &e, bool drxn )
{
    int32_t off = ends_[drxn] - e.node->ends_[!drxn] + ( drxn ? -e.overlap : e.overlap );
    e.node->offset( off );
}

void Node::offsetForward( bool drxn, bool sameGraph, bool notSelf )
{
    if ( !notSelf )
    {
        offsetNode( drxn );
    }
    
    NodeSet alreadySet = { this };
    NodeSet fwdSet = getDrxnNodes( drxn, sameGraph, true );
    NodeSet currSet = getNextNodes( drxn );
    
    Node::offsetForward( currSet, alreadySet, fwdSet, drxn );
}

void Node::offsetForward( NodeSet &tmpCurrSet, bool drxn, bool notSelf )
{
    NodeSet currSet, fwdSet, alreadySet;
    for ( Node* curr : tmpCurrSet )
    {
        curr->getDrxnNodes( fwdSet, drxn );
    }
    
    for ( Node* curr : tmpCurrSet )
    {
        if ( fwdSet.find( curr ) == fwdSet.end() )
        {
            fwdSet.insert( curr );
            curr->getNextNodes( currSet, drxn );
            if ( !notSelf )
            {
                curr->offsetNode( drxn );
            }
            alreadySet.insert( curr );
        }
    }
    
    Node::offsetForward( currSet, alreadySet, fwdSet, drxn );
}

void Node::offsetForward( NodeSet &currSet, NodeSet &alreadySet, NodeSet &fwdSet, bool drxn )
{
    while ( !currSet.empty() )
    {
        NodeSet nxtSet;
        bool didAdvance = false;
        for ( Node* curr : currSet )
        {
            bool valid = true;
            if ( alreadySet.find( curr ) == alreadySet.end() )
            {
                int32_t off = 0;
                bool first = true;
                for ( Edge &e : curr->edges_[!drxn] )
                {
                    if ( fwdSet.find( e.node ) == fwdSet.end() ) continue;
                    valid = valid && alreadySet.find( e.node ) != alreadySet.end();
                    if ( !valid ) continue;
                    int32_t edgeOffset = e.node->ends_[drxn] - curr->ends_[!drxn] + ( drxn ? -e.overlap : e.overlap );
                    off = first ? edgeOffset : ( drxn ? max( off, edgeOffset ) : min( off, edgeOffset ) );
                    first = false;
                }

                if ( valid )
                {
                    curr->offset( off );
                    alreadySet.insert( curr );
                    didAdvance = true;
                }
            }
            
            if ( valid )
            {
                for ( Node* nxt : curr->getNextNodes( drxn ) )
                {
                    if ( currSet.find( nxt ) == currSet.end() && fwdSet.find( nxt ) != fwdSet.end() )
                    {
                        nxtSet.insert( nxt );
                    }
                }
            }
            else
            {
                nxtSet.insert( curr );
            }
        }
        
        if ( !didAdvance ) break;
        assert( didAdvance );
        currSet = nxtSet;
    }
    
}

void Node::offsetIsland( NodeSet &propagated, bool drxn )
{
    if ( propagated.find( this ) == propagated.end() )
    {
        NodeSet currSet = getNextNodes( drxn );
        NodeSet alreadySet = { this };
        NodeSet fwdSet = getDrxnNodes( drxn, true, true );
        
        while ( !currSet.empty() )
        {
            Node::offsetForward( currSet, alreadySet, fwdSet, drxn );
            propagated.insert( alreadySet.begin(), alreadySet.end() );
            
            for ( Node* node : alreadySet )
            {
                for ( Node* nxt : node->getNextNodes( !drxn ) )
                {
                    if ( propagated.find( nxt ) == propagated.end() && nxt->drxn_ > 2 )
                    {
                        currSet.insert( nxt );
                        fwdSet.insert( node );
                        for ( Node* fwd : node->getDrxnNodesNotInSet( propagated, !drxn ) )
                        {
                            if ( fwd->drxn_ > 2 )
                            {
                                fwdSet.insert( fwd );
                            }
                        }
                    }
                }
            }
            
            drxn = !drxn;
        }
    }
}

void Node::setOffsetMap( NodeIntMap &offsets, NodeSet useSet, int32_t limit, bool drxn )
{
    offsets[this] = 0;
    NodeSet currSet = { this };
    while ( !currSet.empty() )
    {
        NodeSet nxtSet;
        for ( Node* curr : currSet )
        {
            int32_t currOffset = offsets[curr];
            for ( Edge &e : curr->edges_[drxn] )
            {
                if ( useSet.find( e.node ) == useSet.end() ) continue;
                if ( offsets.find( e.node ) != offsets.end() ) continue;
                int32_t offset = drxn ? e.node->ends_[0] - curr->ends_[1] + e.overlap
                                      : e.node->ends_[1] - curr->ends_[0] - e.overlap;
                offsets[e.node] = offset + currOffset;
                
                int32_t thisEnd = e.node->ends_[drxn] - offset - currOffset;
                if ( drxn ? thisEnd < limit : limit < thisEnd ) nxtSet.insert( e.node );
            }
        }
        currSet = nxtSet;
    }
}
