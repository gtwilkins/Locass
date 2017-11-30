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

vector< pair<int32_t, uint16_t> > Node::getCoverageMarks( bool drxn )
{
    vector< pair<int32_t, uint16_t> > coverMarks;
    int marks = max( 0, ( 2 * ( ends_[1] - ends_[0] ) / params.readLen ) - 1 );
    int spacing = marks > 0 ? max( marks, ( ends_[1] - ends_[0] - params.readLen ) ) / marks : 0;
    int remainder = ends_[1]- ends_[0] - params.readLen - ( spacing * marks );
    
    int32_t firstMark = drxn ? ends_[0] + ( params.readLen / 2 ) : ends_[1] - ( params.readLen / 2 );
    for ( int i ( 0 ); i < marks; i++ )
    {
        coverMarks.push_back( make_pair( firstMark + ( spacing * ( drxn ? i : -i ) ), 0 ) );
        firstMark += i >= marks - remainder ? 1 : 0;
    }
    coverMarks.push_back( make_pair( drxn ? ends_[1] - ( params.readLen / 2 ) : ends_[0] + ( params.readLen / 2 ), 0 ) );
    
    for ( auto &read : reads_ )
    {
        for ( pair<int32_t, uint16_t> &cover : coverMarks )
        {
            if ( read.second[0] <= cover.first && cover.first <= read.second[1] )
            {
                cover.second++;
            }
        }
    }
    
    if ( !coverMarks.empty() )
    {
        coverMarks[0].second *= 2;
        coverMarks.back().second *= coverMarks.size() > 1 ? 2 : 1;
    }
    
    for ( auto it = coverMarks.begin(); it != coverMarks.end(); )
    {
        if ( drxn ? it->first < params.readLen : it->first > 0 )
        {
            it = coverMarks.erase( it );
            continue;
        }
        it++;
    }
    
    return coverMarks;
}

int32_t Node::getLengthForCoverage()
{
    int length = ends_[1] - ends_[0];
    length -= ( ( !edges_[0].empty() ? getBestOverlap( 0 ) : min( params.readLen, length ) ) 
            + ( !edges_[1].empty() ? getBestOverlap( 1 ) : min( params.readLen, length ) ) )
            / 2;
    return max( 1, length );
}

float Node::getMultiplicity()
{
    float len = getLengthForCoverage() + params.readSpacing * 3;
    int readCount = reads_.size() + 3;
    return ( ( float( readCount * params.readLen ) / len ) / ( params.cover * 0.7 ) );
}

bool Node::isMultiple()
{
    float len = max( params.readLen, ends_[1] - ends_[0] );
    // Cutoff is 1.3x coverage, increasing towards 1.5x as node length increases
    float cutoff = params.cover * ( 1.3 + min( float(0.2), 
            len / float( 20 * params.readLen ) ) );
    return coverage_ > cutoff;
}

bool Node::isReliable( bool doSet )
{
    if ( reliable_ || unreliable_ || !doSet )
    {
        return reliable_ && !unreliable_;
    }
    setReliability();
    return reliable_;
}

void Node::resetReliability( NodeSet &nodes, bool drxn )
{
    NodeSet propagated;
    for ( Node* node : nodes )
    {
        if ( propagated.find( node ) == propagated.end() )
        {
            node->resetReliabilityPropagate( propagated, drxn );
        }
    }
}

void Node::resetReliability()
{
    this->reliable_ = false;
    this->unreliable_ = false;
}

void Node::resetReliabilityPropagate( NodeSet &propagated, bool drxn )
{
    resetReliability();
    propagated.insert( this );
    for ( Node* fwd : getDrxnNodes( drxn ) )
    {
        if ( propagated.find( this ) == propagated.end() )
        {
            fwd->resetReliability();
        }
    }
}

void Node::setCoverage()
{
    float len = getLengthForCoverage();
    float readCount = reads_.size();
    coverage_ = readCount < 3 ? params.cover : ( ( readCount - 2 ) * float( params.readLen) ) / len;
}

void Node::setCoverage( ExtVars &ev, bool subGraph, bool drxn, bool isIsland )
{
    if ( isContinue( drxn ) || ( ends_[1] - ends_[0] < params.readLen * 1.25 ) )
    {
        setCoverage();
        return;
    }
    
    vector< pair<int32_t, uint16_t> > coverMarks = getCoverageMarks( drxn );
    
    if ( coverMarks.size() < 2 || params.cover * 2 < coverMarks[0].second )
    {
        setCoverage();
        return;
    }
    
    uint16_t baseCover = max( uint16_t(params.cover * 0.9), coverMarks[0].second );
    
    int i = 1;
    while ( i < coverMarks.size() && coverMarks[i].second < baseCover * 1.2 )
    {
        baseCover = coverMarks[i].second < baseCover 
                ? max( uint16_t(params.cover * 0.9), coverMarks[i].second ) 
                : baseCover;
        i++;
    }
    
    if ( i != coverMarks.size() )
    {
        int32_t splitMark = findNextRead( coverMarks[i-1].first, drxn );
        if ( splitMark != ends_[drxn] )
        {
            if ( abs( splitMark - ends_[drxn] ) < params.readLen )
            {
                assert( false );
            }
            splitNode( splitMark, ( isIsland ? ev.island : ev.nodes ), subGraph, drxn );
            if ( edges_[drxn].size() == 1 )
            {
                edges_[drxn][0].node->setCoverage( ev, subGraph, drxn, isIsland );
            }
        }
    }
    
    setCoverage();
}

void Node::setReliability()
{
    if ( !unreliable_ && reads_.size() >= 6 && !misassembled_[0] && !misassembled_[1] )
    {
        float len = max( params.readLen, ends_[1] - ends_[0] );
        // Cutoff is 0.9x coverage, increasing to 1.3x
        float cutoff = params.cover * ( 0.9 + min( float(0.4), 
                len / float( 20 * params.readLen ) ) );
        reliable_ = coverage_ < cutoff;
    }
}

void Node::setReliable( bool doForce )
{
    reliable_ = doForce || !unreliable_;
}

void Node::setUnreliable()
{
    unreliable_ = true;
    reliable_ = false;
}
