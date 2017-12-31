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
#include <limits>

void Node::resolveBypass( ExtVars &ev, bool drxn )
{
    NodeSet fwdSet;
    for ( auto it = ev.bypass.begin(); it != ev.bypass.end(); )
    {
        if ( find( ev.nodes.begin(), ev.nodes.end(), *it ) == ev.nodes.end() )
        {
            ev.bypass.erase( it );
            continue;
        }
        (*it)->getDrxnNodes( fwdSet, drxn );
        it++;
    }
    
    for ( Node* node : ev.bypass )
    {
        if ( fwdSet.find( node ) == fwdSet.end() )
        {
            bool isFirst = true;
            while ( node->resolveBypass( ev, isFirst, drxn ) )
            {
                node->offsetForward( drxn, true, true );
                isFirst = false;
            }
            ev.bypass.erase( node );
            break;
        }
    }
}

bool Node::resolveBypass( ExtVars &ev, bool doExtend, bool drxn )
{
    if ( doExtend )
    {
        int best = 0;
        extendForward( ev, 0, best, 8, drxn );
    }
    
    ev.ante = getDrxnNodes( !drxn, true );
    
    for ( Edge &e : edges_[!drxn] )
    {
        NodeSet fwdSet = e.node->getDrxnNodesInSet( ev.ante, drxn );
        if ( !fwdSet.empty() )
        {
            NodeSet offsetSets[2] = { fwdSet, NodeSet() };
            int hitsCount = 0;
            int32_t diffOffset = drxn ? ends_[0] + e.overlap - e.node->ends_[1]
                                      : ends_[1] - e.overlap - e.node->ends_[0];
            if ( reviewOffset( diffOffset, offsetSets, hitsCount, drxn ) )
            {
                e.node->removeEdge( this, drxn );
                e.node->edgeCount_[drxn]--;
                removeEdge( e.node, !drxn );
            }
            else
            {
                for ( auto it = edges_[!drxn].begin(); it != edges_[!drxn].end(); )
                {
                    if ( fwdSet.find( it->node ) != fwdSet.end() )
                    {
                        it->node->removeEdge( this, drxn );
                        it->node->stop_[drxn] = it->node->edges_[drxn].empty();
                        it->node->edgeCount_[drxn]--;
                        it = edges_[!drxn].erase( it );
                        continue;
                    }
                    it++;
                }
            }
            return true;
        }
    }
    
    return false;
}

bool Node::resolveOffset( ExtVars &ev, bool doExtend, bool drxn )
{
    NodeIntList offsets = getOffsetEdges( !drxn );
    bool didErase = false;
    
    if ( doExtend && edges_[!drxn].size() > 1 )
    {
        int dummy = 0;
        extendForward( ev, 0, dummy, 8, drxn );
    }
    
    if ( !ev.bypass.empty() )
    {
        return false;
    }
    
    while ( !offsets.empty() )
    {
        sort( offsets.begin(), offsets.end(), []( pair<Node*, int32_t> &a, pair<Node*, int32_t> &b ){ 
            return abs( a.second < b.second ); 
        } );
        
        NodeSet offsetSets[2];
        NodeSet diffEdges;
        int32_t diffOffset = offsets[0].second;
        
        for ( auto it = offsets.begin(); it != offsets.end(); )
        {
            diffEdges.insert( it->first );
            if ( abs( diffOffset - it->second ) <= 200 )
            {
                offsetSets[1].insert( it->first );
                it->first->getDrxnNodes( offsetSets[1], !drxn );
                it = offsets.erase( it );
                continue;
            }
            it++;
        }
        
        for ( Node* prv : getNextNodes( !drxn ) )
        {
            if ( diffEdges.find( prv ) == diffEdges.end() )
            {
                offsetSets[0].insert( prv );
                prv->getDrxnNodes( offsetSets[0], !drxn );
            }
        }
        
        int hitsCount = 0;
        bool samePref = reviewOffset( diffOffset, offsetSets, hitsCount, drxn );
        samePref = hitsCount > 0 ? samePref : reviewOffset( offsetSets, drxn );
        
        for ( auto it = edges_[!drxn].begin(); it != edges_[!drxn].end(); )
        {
            if ( diffEdges.find( it->node ) != diffEdges.end() ? samePref : !samePref )
            {
                it->node->removeEdge( this, drxn );
                it->node->stop_[drxn] = it->node->edges_[drxn].empty();
                it = edges_[!drxn].erase( it );
                didErase = true;
            }
            else it++;
        }
        
        offsetForward( drxn );
        offsets = getOffsetEdges( !drxn );
    }
    
    return true;
}

void Node::resolveOffsets( ExtVars &ev, bool drxn )
{
    while ( !ev.offset.empty() )
    {
        Node* node = *ev.offset.begin();
        if ( node->resolveOffset( ev, true, drxn ) )
        {
            node->offsetForward( drxn, true, true );
            ev.offset.erase( node );
        }
        else
        {
            Node::resolveBypass( ev, drxn );
        }
    }
}

void Node::resolveRebranch( ExtVars &ev, bool drxn )
{
    while ( !ev.rebranch.empty() )
    {
        NodeSet rebranch = ev.rebranch;
        ev.rebranch.clear();

        for ( Node* node : rebranch )
        {
            if ( find( ev.nodes.begin(), ev.nodes.end(), node ) != ev.nodes.end() )
            {
                node->rebranchNode( ev, drxn );
            }
        }
    }
}

bool Node::reviewOffset( int32_t diffOffset, NodeSet* offsetSets, int& hitsCount, bool drxn )
{
    NodeSet tCloseSet, tFarSet;
    NodeSet fwdSet = getDrxnNodes( drxn, false, true );
    for ( Node* prv : getNextNodes( !drxn ) )
    {
        tCloseSet.insert( prv );
        tFarSet.insert( prv );
        prv->getDrxnNodes( tCloseSet, !drxn, params.getFurthestPeDist( prv->ends_[!drxn], !drxn ) );
        prv->getDrxnNodes( tFarSet, !drxn, params.getFurthestMpDist( prv->ends_[!drxn], !drxn ) );
    }
    
    int onlyScores[2] = { 0, 0 };
    int32_t totalOffset = 0;
    for ( Node* t : tCloseSet )
    {
        for ( ReadMark &mark : t->getMarksBase( !drxn ) )
        {
            for ( Node* fwd : fwdSet )
            {
                auto it = fwd->reads_.find( mark.id );
                if ( it != fwd->reads_.end() )
                {
                    onlyScores[0] += offsetSets[0].find( t ) != offsetSets[0].end() 
                                  && offsetSets[1].find( t ) == offsetSets[1].end();
                    onlyScores[1] += offsetSets[0].find( t ) == offsetSets[0].end() 
                                  && offsetSets[1].find( t ) != offsetSets[1].end();
                    totalOffset += it->second[drxn] - mark.estimate;
                    hitsCount++;
                    break;
                }
            }
        }
    }
    
    NodeList tNodes( tFarSet.begin(), tFarSet.end() );
    for ( int i : { 0, 1 } )
    {
        for ( Node* node : offsetSets[i] )
        {
            for ( ReadMark &mark : node->getMarksBase( drxn ) )
            {
                for ( Node* t : tNodes )
                {
                    if ( offsetSets[i].find( t ) == offsetSets[i].end() )
                    {
                        auto it = t->reads_.find( mark.id );
                        onlyScores[i] += ( it != t->reads_.end() && mark.isValid( it->second ) );
                    }
                }
            }
        }
    }
    
    if ( hitsCount )
    {
        totalOffset /= hitsCount;
        return onlyScores[0] > onlyScores[1] || ( onlyScores[0] == onlyScores[1] && abs( totalOffset ) < abs( totalOffset - diffOffset ) );
    }
    
    return true;
}

bool Node::reviewOffset( NodeSet* offsetSets, bool drxn )
{
    NodeList tNodes, qNodes[2];
    for ( bool i : { 0, 1 } )
    {
        for ( Node* node : offsetSets[i] )
        {
            if ( offsetSets[!i].find( node ) == offsetSets[!i].end() )
            {
                qNodes[i].push_back( node );
                
            }
            else if ( !i )
            {
                tNodes.push_back( node );
            }
        }
    }
    
    int hits[2] = { 0, 0 }, reliable[2] = { 0, 0 };
    
    for ( bool i : { 0, 1 } )
    {
        for ( Node* q : qNodes[i] )
        {
            for ( ReadMark &mark : q->getMarksBase( drxn ) )
            {
                for ( Node* t : tNodes )
                {
                    if ( t->reads_.find( mark.id ) != t->reads_.end() )
                    {
                        hits[i]++;
                        reliable[i] += t->reliable_;
                        break;
                    }
                }
            }
        }
    }
    
    return !( reliable[1] > reliable[0] || ( reliable[0] == reliable[1] && hits[1] > hits[0] ) );
}
