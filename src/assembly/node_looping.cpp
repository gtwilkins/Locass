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

void Node::loop( LoopVars &lv, ExtVars &ev, bool drxn )
{
    int best = 0;
        
    while ( !lv.open.empty() )
    {
        Node* nxt = loopWalkGetNext( lv, drxn );
        nxt->loopWalk( lv, ev, best, drxn );
    }

    Node::loopReview( lv, ev, drxn );
}

void Node::loopDelete( LoopVars &lv, ExtVars &ev, bool drxn )
{
    NodeSet delSet;
    this->dismantleNode( delSet, drxn );
    for ( Node* del : delSet )
    {
        lv.open.erase( del );
        lv.branches.erase( del );
        lv.clones.erase( del );
        lv.cloned.erase( del );
        lv.scores.erase( del );
        ev.nodes.erase( remove( ev.nodes.begin(), ev.nodes.end(), del ), ev.nodes.end() );
        ev.cloneSet.erase( del );
        ev.bypass.erase( del );
    }
    ev.del.insert( delSet.begin(), delSet.end() );
}

void Node::loopReview( LoopVars &lv, ExtVars &ev, bool drxn )
{
    NodeSet fwdSet = { lv.startList.begin(), lv.startList.end() };
    NodeSet branchFwd, anteEnds;
    for ( Node* node : lv.startList )
    {
        node->getDrxnNodes( fwdSet, drxn );
    }
    
    for ( Node* branch : lv.branches )
    {
        if ( fwdSet.find( branch ) != fwdSet.end() )
        {
            branchFwd.insert( branch );
            branch->getDrxnNodes( branchFwd, drxn );
        }
    }
    
    // Rescore
    NodeSet currSet = { lv.startList.begin(), lv.startList.end() };
    while ( !currSet.empty() )
    {
        NodeSet nxtSet, currFwd;
        bool didAdvance = false;
        for ( Node* curr : currSet )
        {
            curr->getDrxnNodes( currFwd, drxn );
        }
        
        for ( Node* curr : currSet )
        {
            if ( currFwd.find( curr ) != currFwd.end() )
            {
                nxtSet.insert( curr );
                continue;
            }
            
            didAdvance = true;
            
            if ( branchFwd.find( curr ) == branchFwd.end() )
            {
                int prevScore = 0;
                bool didSet = false;

                for ( Node* prv : curr->getNextNodes( !drxn ) )
                {
                    auto it = lv.scores.find( prv );
                    if ( it != lv.scores.end() && ( !didSet || it->second.cumulScore > prevScore ) )
                    {
                        prevScore = it->second.cumulScore;
                    }
                }

                auto it = lv.scores.find( curr );
                if ( it != lv.scores.end() )
                {
                    int dummy = 0;
                    it->second.setCumul( curr, prevScore, dummy );
                    curr->getNextNodes( nxtSet, drxn );
                    if ( curr->edges_[drxn].empty() )
                    {
                        anteEnds.insert( curr );
                    }
                }
            }
        }
        
        assert( didAdvance || nxtSet.empty() );
        currSet = nxtSet;
    }
    
    while ( !anteEnds.empty() )
    {
        NodeSet nxtSet;
        for ( Node* node : anteEnds )
        {
            node->stop_[drxn] = 1;
            if ( lv.scores[node].cumulScore <= 0 && node->edges_[drxn].empty() )
            {
                for ( Node* prv : node->getNextNodes( !drxn ) )
                {
                    if ( fwdSet.find( prv ) != fwdSet.end() )
                    {
                        nxtSet.insert( prv );
                    }
                }
                node->loopDelete( lv, ev, drxn );
            }
        }
        anteEnds = nxtSet;
    }
}

void Node::loopWalk( LoopVars &lv, ExtVars &ev, int &best, bool drxn )
{
    lv.open.erase( this );
    NodeSet bckSet = getDrxnNodesInSet( lv.clones, !drxn );
    int thisCutoff = max( lv.cutoff + best + (int)bckSet.size(), -1 );
    
    if ( loopWalkIsValidClone( lv, ev, thisCutoff, drxn ) )
    {
        // Set next
        loopWalkEdges( lv, ev, drxn );
        
        // Extend and score branches
        auto itScore = lv.scores.find( this );
        int nxtBest = 0;
        for ( Node* nxt : getNextNodes( drxn ) )
        {
            if ( lv.branches.find( nxt ) != lv.branches.end() 
                    && !nxt->extendForward( ev, min( lv.target - 2, itScore->second.cumulScore ), nxtBest, lv.target, drxn ) )
            {
                nxt->loopDelete( lv, ev, drxn );
            }
        }
        
        // Update best score and reset cumulative score if target reached
        best = max( nxtBest + itScore->second.cumulScore, best );
        lv.exitFound = lv.exitFound || nxtBest > 1;
        
        // Propagate remaining open clones
        for ( Node* nxt : getNextNodes( drxn ) )
        {
            if ( lv.open.find( nxt ) != lv.open.end() && loopWalkDoWalk( lv, drxn ) )
            {
                nxt->loopWalk( lv, ev, best, drxn );
                return;
            }
        }
    }
}

bool Node::loopWalkDoWalk( LoopVars &lv, bool drxn )
{
    if ( edges_[!drxn].size() < getCloneEdgeSets( !drxn ).size() )
    {
        NodeSet fwdSet;
        for ( Node* clone : getCloneSet() )
        {
            clone->getDrxnNodes( fwdSet, drxn );
        }
        
        for ( Node* opn : lv.open )
        {
            bool anyInFwd = false;
            NodeSet opnFwdSet;
            for ( Node* clone : opn->getCloneSet() )
            {
                anyInFwd = anyInFwd || fwdSet.find( clone ) != fwdSet.end();
                clone->getDrxnNodes( opnFwdSet, drxn );
            }
            
            if ( !anyInFwd )
            {
                for ( Node* clone : getCloneSet() )
                {
                    if ( opnFwdSet.find( clone ) != opnFwdSet.end() )
                    {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

void Node::loopWalkEdges( LoopVars &lv, ExtVars &ev, bool drxn )
{
    vector<Edge> edges;
    NodeSet usedEdges = getCloneSet( true );

    loopWalkEdgesCatalog( lv, usedEdges, edges, drxn );
    loopWalkEdgesSet( lv, ev, usedEdges, edges, drxn );
}

void Node::loopWalkEdgesCatalog( LoopVars &lv, NodeSet &usedEdges, vector<Edge> &edges, bool drxn )
{
    // Block pre-existing forward edge nodes if they exist for some reason
    for ( Edge &e : edges_[drxn] )
    {
        if ( e.node->clones_ )
        {
            usedEdges.insert( e.node->clones_->begin(), e.node->clones_->end() );
        }
    }
    
    // Catalog all unique forward edge nodes for each clone;
    for ( Node* clone : getCloneSet( false ) )
    {
        int setCount = 0;
        for ( NodeSet &cloneFwd : lv.cloneFwdSets )
        {
            setCount += cloneFwd.find( clone ) != cloneFwd.end();
        }
        
        for ( Edge &e : clone->edges_[drxn] )
        {
            int edgeSetCount = 0;
            for ( NodeSet &cloneFwd : lv.cloneFwdSets )
            {
                edgeSetCount += cloneFwd.find( e.node ) != cloneFwd.end();
            }
            
            if ( usedEdges.find( e.node ) == usedEdges.end() && edgeSetCount == setCount )
            {
                edges.push_back( e );
                usedEdges.insert( e.node );
                
                if ( e.node->clones_ )
                {
                    usedEdges.insert( e.node->clones_->begin(), e.node->clones_->end() );
                }
            }
        }
    }
}

void Node::loopWalkEdgesSet( LoopVars &lv, ExtVars &ev, NodeSet &usedEdges, vector<Edge> &edges, bool drxn )
{
    ev.ante = getDrxnNodes( !drxn, true );
    
    for ( Edge &e : edges )
    {
        // Merge into already cloned convergent node
        if ( lv.cloned.find( e.node ) != lv.cloned.end() )
        {
            for ( Node* clone : e.node->getCloneSet( true ) )
            {
                if ( lv.clones.find( clone ) != lv.clones.end() )
                {
                    if ( ev.ante.find( clone ) == ev.ante.end() )
                    {
                        addEdge( clone, e.overlap, drxn );
                    }
                    break;
                }
            }
        }
        
        // Clone ante node that has not yet been cloned
        else if ( ev.ante.find( e.node ) != ev.ante.end() )
        {
            assert( ev.cloneSet.find( e.node ) == ev.cloneSet.end() );
            Node* clone = new Node( e.node, ev, drxn );
            addEdge( clone, e.overlap, drxn );
            ev.nodes.push_back( clone );
            lv.clones.insert( clone );
            lv.cloned.insert( clone );
            lv.cloned.insert( clone->clones_->begin(), clone->clones_->end() );
            lv.open.insert( clone );
        }
        
        // Contest for divergent branch
        else
        {
            addEdge( e.node, e.overlap, drxn );
            bool isFirst = true;
            while ( e.node->resolveBypass( ev, isFirst, drxn ) )
            {
                e.node->offsetForward( drxn, true, true );
                isFirst = false;
            }
            ev.ante = getDrxnNodes( !drxn, true );
            NodeSet nxtSet = getNextNodes( drxn );
            if ( nxtSet.find( e.node ) != nxtSet.end() )
            {
                lv.branches.insert( e.node );
            }
        }
    }
    
    // Rebranch if necessary
    if ( edges.size() < edgeCount_[drxn] )
    {
        NodeSet acceptable;
        NodeSet nxtSet = getNextNodes( drxn );
        usedEdges.insert( this );
        vector<Extension> exts = ev.bwt.mapExtensions( seq_, drxn );
        for ( Extension &ext : exts )
        {
            addExtension( ext, ev, acceptable, usedEdges, true, drxn );
        }
        
        for ( Node* nxt : getNextNodes( drxn ) )
        {
            if ( nxtSet.find( nxt ) == nxtSet.end() )
            {
                lv.branches.insert( nxt );
            }
        }
    }
}

Node* Node::loopWalkGetNext( LoopVars &lv, bool drxn )
{
    Node* curr = *lv.open.begin();
    NodeSet tested = { curr };
    bool updated = true;
    while ( curr->edges_[!drxn].size() < curr->getCloneEdgeSets( !drxn ).size() && updated )
    {
        updated = false;
        for ( Node* opn : lv.open )
        {
            if ( !updated && tested.find( opn ) == tested.end() )
            {
                for ( Node* opnClone : *opn->clones_ )
                {
                    NodeSet fwdSet;
                    opnClone->getDrxnNodesNotInSet( fwdSet, lv.cloned, drxn );
                    for ( Node* currClone : *curr->clones_ )
                    {
                        if ( !updated && fwdSet.find( currClone ) != fwdSet.end() )
                        {
                            curr = opn;
                            updated = true;
                        }
                    }
                }
            }
        }
    }
    return curr;
}

bool Node::loopWalkIsValidClone( LoopVars &lv, ExtVars &ev, int &cutoff, bool drxn )
{
    // Set this clone's score
    CloneScore score = getCloneScore( drxn );
    score.setScores();
    int prevScore = 0, prevMisses = 0;
    bool anyBackClones = false;
    for ( Node* node : getNextNodes( !drxn ) )
    {
        auto hit = lv.scores.find( node );
        if ( hit != lv.scores.end() )
        {
            if ( !anyBackClones 
                    || hit->second.cumulScore > prevScore 
                    || ( hit->second.cumulScore == prevScore && hit->second.cumulMisses < prevMisses ) )
            {
                anyBackClones = true;
                prevScore = hit->second.cumulScore;
                prevMisses = hit->second.cumulMisses;
            }
        }
    }
    
    // Set cumulative score
    prevScore = lv.exitFound ? prevScore - 2 : prevScore;
    score.setCumul( this, prevScore, prevMisses );
    
    lv.scores[this] = score;
    
    // Determine if well supported or insufficiently unsupported
    return score.cumulScore - score.cumulMisses >= cutoff;
}

void Node::resolveLoops( ExtVars &ev, bool drxn )
{
    NodeSet cloneFwd;
    for ( auto it = ev.cloneSet.begin(); it != ev.cloneSet.end(); )
    {
        // Erase not clone
        if ( !(*it)->clones_ )
        {
            it = ev.cloneSet.erase( it );
            continue;
        }
        
        // Hard cap on cloning
        if ( (*it)->clones_->size() >= 3 )
        {
            (*it)->dismantleNode( ev.del, drxn );
            it = ev.cloneSet.erase( it );
            continue;
        }
        
        // Compile forward nodes for determining set of clones to loop
        for ( Node* clone : (*it)->getCloneSet() )
        {
            clone->getDrxnNodes( cloneFwd, drxn );
        }
        it++;
    }
    
    for ( auto it = ev.cloneSet.begin(); it != ev.cloneSet.end(); )
    {
        if ( ev.del.find( *it ) == ev.del.end() ) it++;
        else it = ev.cloneSet.erase( it );
    }
    
    if ( ev.cloneSet.empty() ) return;
    
    LoopVars lv;
    float farness = 1;
    ev.ante.clear();
    
    for ( Node* node : ev.cloneSet )
    {
        float thisFarness = 1;
        for ( Node* clone : node->getCloneSet() )
        {
            thisFarness = min( thisFarness, (float)abs( node->ends_[drxn] - clone->ends_[drxn] ) / (float)250 );
            if ( cloneFwd.find( clone ) == cloneFwd.end() && lv.addLoopStart( node, drxn ) )
            {
                ev.ante.insert( node );
                node->getDrxnNodes( ev.ante, !drxn, true );
                farness = min( farness, thisFarness );
                break;
            }
        }
    }
    
    lv.cutoff = -3 - ( farness * 9 );
    lv.target = 13 + lv.cutoff;
    lv.exitFound = false;
    
    for ( Node* node : lv.startList )
    {
        ev.cloneSet.erase( node );
    }
    
    if ( lv.startList.empty() && !ev.cloneSet.empty() )
    {
        for ( Node* node : ev.cloneSet )
        {
            node->dismantleNode( ev.del, drxn );
        }
        ev.cloneSet.clear();
        return;
    }
    assert( !lv.startList.empty() && lv.startList.size() <= 2 );
    
    Node::loop( lv, ev, drxn );
}

