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

bool Node::anyValid( bool drxn )
{
    return ( drxn ? validLimits_[2] > ends_[0] : validLimits_[1] < ends_[1] );
}

void Node::forceValidate( int32_t* limits, NodeSetList &delSets, NodeSetList &goodSets )
{
    for ( bool drxn : { 0, 1 } )
    {
        int hits = 0;
        for ( Node* fwd : getDrxnNodes( drxn ) )
        {
            hits += fwd->getPairHitsTotal();
        }

        if ( hits < 2 )
        {
            bool doTrim = true;
            for ( Node* node : getNextNodes( drxn ) )
            {
                doTrim = doTrim && node->drxn_ != 2;
                if ( node->drxn_ != 2 )
                {
                    delSets[drxn].insert( node );
                }
            }
            if ( doTrim )
            {
                trimEnd( drxn );
            }
        }
        else
        {
            getDrxnNodes( goodSets[drxn], drxn );
        }
    }
    setValid( limits );
}

bool Node::isValidated()
{
    return validated_;
}

bool Node::isValidated( bool drxn )
{
    return !isContinue( drxn ) && ( drxn ? ends_[1] <= validLimits_[2] : validLimits_[1] <= ends_[0] );
}

int32_t Node::getValidLimit( bool drxn )
{
    return validLimits_[1 + drxn];
}

bool Node::isValidCoords( Coords &coords, int drxn )
{
    if ( validated_ )
    {
        return true;
    }
    else if ( drxn == 1 )
    {
        return coords[1] <= validLimits_[2];
    }
    else if ( drxn == 0 )
    {
        return validLimits_[1] <= coords[0];
    }
    else
    {
        return validLimits_[1] <= coords[0] && coords[1] <= validLimits_[2];
    }
}

bool Node::isValidHit( Coords* coords, bool drxn )
{
    if ( drxn ? drxn_ == 0 || (*coords)[1] <= validLimits_[2]
              : drxn_ == 1 || validLimits_[1] <= (*coords)[0] )
    {
        return true;
    }
    
    return false;
}

void Node::propagateValidation( int32_t* limits, bool drxn )
{
    NodeSet testedSet, backedSet;
    propagateValidation( limits, testedSet, backedSet, drxn );
}

void Node::propagateValidation( NodeSet &qSet, int32_t* limits, bool drxn )
{
    NodeSet testedSet, backedSet;
    for ( Node* q : qSet )
    {
        for ( Node* qClone : q->getCloneSet( true ) )
        {
            if ( testedSet.find( qClone ) == testedSet.end() && backedSet.find( qClone ) == backedSet.end() )
            {
                qClone->propagateValidation( limits, testedSet, backedSet, drxn );
            }
        }
    }
}

void Node::propagateValidation( int32_t* limits, NodeSet &testedSet, NodeSet &backedSet, bool drxn )
{
    bool doValidate = drxn_ == 2;
    
    // Navigate backwards until a node with only valid back nodes is found
    if ( drxn ? drxn_ == 1 : drxn_ == 0 )
    {
        for ( Node* bck : getNextNodes( !drxn ) )
        {
            doValidate = doValidate || bck->isValidated() || ( bck->isValidated( drxn ) && testedSet.find( bck ) != testedSet.end() );
//            doValidate = doValidate && ( bck->isValidated() || ( bck->isValidated( drxn ) && testedSet.find( bck ) != testedSet.end() ) );
            if ( !bck->isValidated() 
                    && ( drxn ? bck->drxn_ != 0 : bck->drxn_ != 1 )
                    && backedSet.find( bck ) == backedSet.end() 
                    && testedSet.find( bck ) == testedSet.end() )
            {
                backedSet.insert( bck );
                bck->propagateValidation( limits, testedSet, backedSet, drxn );
            }
        }
    }
    
    // Attempt to validate if all back nodes are valid
    if ( doValidate && testedSet.find( this ) == testedSet.end() )
    {
        validate( limits );
        testedSet.insert( this );
    }
    
    // Propagate forward if valid
    if ( isValidated( drxn ) )
    {
        for ( Node* nxt : getNextNodes( drxn ) )
        {
            if ( nxt->drxn_ <= 2 && testedSet.find( nxt ) == testedSet.end() )
            {
                nxt->propagateValidation( limits, testedSet, backedSet, drxn );
            }
        }
    }
}

void Node::pullValidLimits()
{
    validLimits_[0] = validLimits_[1];
    validLimits_[3] = validLimits_[2];
}

void Node::pushValidLimits( int32_t mark, bool drxn )
{
    if ( drxn )
    {
        if ( mark > validLimits_[3] )
        {
            validLimits_[2] = validLimits_[3];
            validLimits_[3] = mark;
        }
        else if ( mark > validLimits_[2] )
        {
            validLimits_[2] = mark;
        }
    }
    else
    {
        if ( mark < validLimits_[0] )
        {
            validLimits_[1] = validLimits_[0];
            validLimits_[0] = mark;
        }
        else if ( mark < validLimits_[1] )
        {
            validLimits_[1] = mark;
        }
    }
}

void Node::pushValidLimits( Node* t, int32_t qMark, int32_t tMark, bool drxn )
{
    pushValidLimits( qMark, drxn );
    t->pushValidLimits( tMark, !drxn );
    if ( this != t )
    {
        pushValidLimits( ends_[!drxn], !drxn );
        t->pushValidLimits( t->ends_[drxn], drxn );
        for ( Node* node : getBetweenNodes( t, !drxn ) )
        {
            node->pushValidLimits( node->ends_[0], 0 );
            node->pushValidLimits( node->ends_[1], 1 );
        }
    }
}

void Node::pushValidLimts( NodeSet &bckSet, int hits, bool drxn )
{
    NodeSet newBckSet;
    NodeSet fwdSet = getDrxnNodesInSet( bckSet, drxn );
    for ( Node* node : bckSet )
    {
        if ( fwdSet.find( node ) != fwdSet.end() )
        {
            for ( int i ( 0 ); i < min( hits, 2 ); i++ )
            {
                node->pushValidLimits( node->ends_[0], 0 );
                node->pushValidLimits( node->ends_[1], 1 );
            }
        }
        
        if ( node->ends_[0] < node->validLimits_[1] || node->validLimits_[2] < node->ends_[1] )
        {
            newBckSet.insert( node );
        }
    }
    
    bckSet = newBckSet;
}

void Node::setValid( bool drxn )
{
    if ( drxn )
    {
        validLimits_[0] = min( validLimits_[0], ends_[1] );
        validLimits_[1] = min( validLimits_[1], ends_[1] );
        validLimits_[2] = validLimits_[3] = ends_[1];
    }
    else
    {
        validLimits_[2] = max( validLimits_[2], ends_[0] );
        validLimits_[3] = max( validLimits_[3], ends_[0] );
        validLimits_[0] = validLimits_[1] = ends_[0];
    }
}

void Node::setValid()
{
    validLimits_[0] = validLimits_[1] = ends_[0];
    validLimits_[2] = validLimits_[3] = ends_[1];
    for ( int drxn : { 0, 1 } )
    {
        if ( drxn == drxn_ || drxn_ == 2 )
        {
            NodeList tNodes = getTargetNodes( drxn, true );
            setPairs( tNodes, drxn );
            
            for ( Node* fwd : getDrxnNodes( drxn ) )
            {
                NodeList tSelf = { this };
                fwd->setPairs( tSelf, drxn );
            }
        }
    }
    
    validated_ = true;
    if ( paired_ && drxn_ <= 2 )
    {
        delete paired_;
        paired_ = NULL;
    }
}

void Node::setValid( int32_t* limits )
{
    setValid();
//    this->deleteTest( drxn_ );
    
    limits[0] = min( limits[0], ends_[0] );
    limits[1] = max( limits[1], ends_[1] );
}

void Node::resetValid()
{
    validLimits_[0] = validLimits_[1] = ends_[0];
    validLimits_[2] = validLimits_[3] = ends_[1];
}

void Node::reviewEnd( ScoreMap &scores, NodeSet &validSet, NodeSet &goodSet, NodeSet &delSet, bool drxn )
{
    NodeSet endSet = { this };
    getDrxnNodesNotInSet( endSet, validSet, !drxn );
    
    Score score;
    for ( Node* node : endSet )
    {
        score += scores[node];
    }
    
    NodeSet* whichSet = ( score[1] * 3 + 10 >= score[0] ) ? &goodSet : &delSet;
    (*whichSet).insert( endSet.begin(), endSet.end() );
}

void Node::reviewFork( ScoreMap &scores, NodeSet &goodSet, NodeSet &badSet, bool drxn )
{
    float baseMisses = float( getEndMarks( drxn ) ) / min( (float)1, getMultiplicity() );
    NodeSet fwdSet = getDrxnNodes( drxn ), thisEnds;
    for ( Node* fwd : fwdSet )
    {
        if ( fwd->isContinue( drxn ) )
        {
            thisEnds.insert( fwd );
        }
    }
    
    float cutoff = max( (float)2, min( (float)10, (float)50 / max( (float)thisEnds.size(), (float)1 ) ) );
    for ( Node* fwd : fwdSet )
    {
        if ( fwd->isContinue( drxn ) )
        {
            Score score = scores[fwd];
            NodeSet endSet = fwd->getDrxnNodesInSet( fwdSet, !drxn );
            for ( Node* node : endSet )
            {
                score += scores[node];
            }
            endSet.insert( this );
            endSet.insert( fwd );
            ( score[1] * 3 + cutoff >= score[0] + baseMisses ? goodSet : badSet ).insert( endSet.begin(), endSet.end() );
        }
    }
}

bool Node::testInvalid( int32_t* limits, bool drxn )
{
    NodeSet testedSet, backedSet;
    propagateValidation( limits, testedSet, backedSet, drxn );
}

void Node::trimIsland( IslandVars &iv, NodeSet &bgnSet )
{
    NodeSet fwdSet( bgnSet.begin(), bgnSet.end() ), goodSet, islandSet;
    for ( Node* node : bgnSet )
    {
        if ( islandSet.find( node ) == islandSet.end() ) node->getConnectedNodes( islandSet, true );
        node->getDrxnNodes( fwdSet, iv.drxn );
        node->getDrxnNodes( goodSet, !iv.drxn );
    }
    
    NodeSet hitBckSet, reliBckSet;
    NodeIntMap hitMap, reliMap;
    for ( Node* fwd : fwdSet )
    {
        int hits = 0, reli = 0;
        for ( auto &np : fwd->pairs_ )
        {
            if ( np.first->drxn_ > 2 ) continue;
            hits += np.second;
            if ( !np.first->isReliable() 
                    && ( np.first->ends_[1] - np.first->ends_[0] < params.readLen * 1.5 
                            || np.first->coverage_ > params.cover * 1.2 ) ) continue;
            reli += np.second;
        }
        if ( hits )
        {
            goodSet.insert( fwd );
            fwd->getDrxnNodesInSet( goodSet, fwdSet, !iv.drxn );
            hitMap[fwd] = hits;
            fwd->getDrxnNodesInSet( hitBckSet, fwdSet, !iv.drxn );
        }
        if ( reli )
        {
            reliMap[fwd] = reli;
            fwd->getDrxnNodesInSet( reliBckSet, fwdSet, !iv.drxn );
        }
    }
    
    trimIsland( iv, fwdSet, hitBckSet, hitMap, goodSet );
    trimIsland( iv, fwdSet, reliBckSet, reliMap, goodSet );
    
    for ( Node* isl : islandSet )
    {
        if ( goodSet.find( isl ) != goodSet.end() ) continue;
        isl->dismantleNode();
        iv.ev.del.insert( isl );
    }
}

void Node::trimIsland( IslandVars &iv, NodeSet &fwdSet, NodeSet &hitBckSet, NodeIntMap &hitMap, NodeSet &goodSet )
{
    int maxHits = 0;
    Node* maxNode = NULL;
    for ( auto &hit : hitMap )
    {
        if ( hitBckSet.find( hit.first ) != hitBckSet.end() ) continue;
        int hits = hit.second;
        NodeSet hitSet = hit.first->getDrxnNodesInSet( fwdSet, !iv.drxn, false );
        for ( Node* node : hitSet )
        {
            auto it = hitMap.find( node );
            if ( it != hitMap.end() ) hits += it->second;
        }
        if ( hits > maxHits )
        {
            maxNode = hit.first;
            maxHits = hits;
        }
    }
    
    if ( maxNode )
    {
        maxNode->getDrxnNodes( goodSet, iv.drxn );
    }
}

bool Node::validate( int32_t* limits )
{
    if ( validated_ )
    {
        return validated_;
    }
    
    NodeSet pairedSet;
    if ( paired_ )
    {
        pairedSet = *paired_;
    }
    
    for ( int drxn : { 0, 1 } )
    {
        if ( drxn_ == 2 || ( drxn ? drxn_ % 3 == 1 : drxn_ % 3 == 0 ) )
        {
            NodeList tNodes = getTargetNodes( drxn ), tSelf = { this };
            for ( Node* fwd : getDrxnNodes( drxn, false, true ) )
            {
                fwd->setPairs( tNodes, drxn );
            }
            int32_t validLimit = ends_[!drxn];
            NodeSet fwdSet = getDrxnNodes( drxn, false, true );
            for ( Node* fwd : fwdSet )
            {
                if ( fwd->validated_ )
                {
                    validLimits_[1+drxn] = validLimits_[drxn * 3] = ends_[drxn];
                    break;
                }
            }
            while( validLimit != validLimits_[1+drxn] )
            {
                validLimit = validLimits_[1+drxn];
                for ( Node* fwd : fwdSet )
                {
                    fwd->setPairs( tSelf, drxn );
                }
            }
        }
    }
    
    
    if ( validLimits_[1] <= ends_[0] && ends_[1] <= validLimits_[2] && !isContinue( 0 ) && !isContinue( 1 ) )
    {
        setValid( limits );
    }
    
    return validated_;
}

bool Node::validate( bool drxn )
{
    if ( validated_ )
    {
        return validated_;
    }
    
    if ( validLimits_[1] <= validLimits_[2] )
    {
        NodeOffsetMap fwdMap = getDrxnNodesOffset( drxn, 0, false );
        NodeOffsetMap revMap = getDrxnNodesOffset( !drxn, 0, false );
        setPairs( fwdMap, revMap, drxn );
    }
    
    if ( validLimits_[1] == ends_[0] && validLimits_[2] == ends_[1] 
            && !isContinue( 0 ) && !isContinue( 1 ) )
    {
        validated_ = true;
    }
    
    return validLimits_[1+drxn] == ends_[drxn];
}

void Node::validateMerge( Node* mergeNode, bool drxn )
{
    NodeSet qSet = { this };
    NodeSet blockSet = mergeNode->getDrxnNodesNotInSet( qSet, !drxn );
    getDrxnNodesNotInSet( qSet, blockSet, !drxn );
    
    for ( Node* q : qSet )
    {
        q->validLimits_[0] = q->validLimits_[1] = q->ends_[0];
        q->validLimits_[2] = q->validLimits_[3] = q->ends_[1];
    }
    
    for ( Node* q : qSet )
    {
        q->setValid();
    }
}

