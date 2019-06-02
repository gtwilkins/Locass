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

float Node::getAdjustedReadCount( int32_t* limits )
{
    float coeff = coverage_ / params.cover;
    float reads = reads_.size();
    for ( int i : { 0, 1 } )
    {
        for ( ReadMark const &mark : marks_[i] )
        {
            if ( mark.estimate < limits[0] || limits[1] < mark.estimate ) --reads;
        }
    }
    if ( coeff <= 1 ) return reads;
    return reads * sqrt(coeff) / coeff;
}

float Node::getAlleleCoverage( Node* forks[2], NodeList paths[2], bool drxn )
{
    
    int32_t bestOls[2]{ forks[0]->getBestOverlap( 1 ), forks[1]->getBestOverlap( 0 ) };
    int32_t lens[2]{ 1 };
    int reads[2]{0};
    for ( int i = 0; i < 2; i++ )
    {
        if ( paths[i].empty() ) continue;
        lens[i] = -( bestOls[0] + bestOls[1] ) / 2;
        lens[i] += bestOls[0] + bestOls[1] - forks[0]->getOverlap( paths[i][0], 1 ) - forks[1]->getOverlap( paths[i].back(), 0 );
        lens[i] += paths[i][0]->seq_.length();
        reads[i] = paths[i][0]->reads_.size();
        for ( int j = 1; j < paths[i].size(); j++ )
        {
            lens[i] += paths[i][j]->seq_.length() - paths[i][j]->getOverlap( paths[i][j-1], 0 );
            reads[i] += paths[i][j]->reads_.size();
        }
    }
    float cover = float( ( reads[0] + reads[1] ) * params.readLen ) 
                / float( float( max( lens[0], 2 ) + max( lens[1], 2 ) ) / (float)2 );
    
    return cover;
}

float Node::getCoverageCoeff()
{
    if ( coverage_ < params.cover ) return (float)1;
    return sqrt( params.cover / coverage_ );
}

float Node::getCoverageDrxn( int32_t dist, bool drxn, bool inclSelf )
{
    int32_t len = max( 1, dist - params.readLen );
    int32_t limits[2];
    if ( !inclSelf )
    {
        limits[!drxn] = drxn ? ends_[1] - getBestOverlap( 1 ) + params.readLen
                             : ends_[0] + getBestOverlap( 0 ) - params.readLen;
    }
    else
    {
        limits[!drxn] = drxn ? ends_[0] + params.readLen : ends_[1] - params.readLen;
    }
    limits[drxn] = drxn ? limits[0] + len : limits[1] - len;
    int covers[len]{0};
    
    NodeIntMap offsetMap = { make_pair( this, 0 ) };
    NodeSet currSet = { this };
    while ( !currSet.empty() )
    {
        NodeSet nxtSet;
        for ( Node* curr : currSet )
        {
            int currOffset = offsetMap[curr];
            for ( Edge &e : curr->edges_[drxn] )
            {
                int offset = drxn ? curr->ends_[1] - e.ol - e.node->ends_[0]
                                  : e.node->ends_[1] - e.ol - curr->ends_[0];
                offset += currOffset;
                if ( offsetMap.find( e.node ) == offsetMap.end()
                        && limits[0] < e.node->ends_[1] + offset 
                        && e.node->ends_[0] + offset < limits[1] )
                {
                    offsetMap[e.node] = offset;
                    nxtSet.insert( e.node );
                }
            }
        }
        currSet = nxtSet;
    }
    if ( !inclSelf ) offsetMap.erase( this );
    
    for ( auto &node : offsetMap )
    {
        for ( auto read : node.first->reads_ )
        {
            int i = max( limits[0], read.second[0] + node.second ) - limits[0];
            int j = min( limits[1], read.second[1] + node.second ) - limits[0];
            while ( i < j ) covers[i++]++;
        }
    }
    
    int coverSum = 0;
    for ( int i = 0; i < len; i++ ) coverSum += min( (int)params.cover, covers[i] );
    float cover = (float)coverSum / float(len);
    
    return cover;
}

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

bool Node::getReliability()
{
    if ( !unreliable_ && reads_.size() >= 6 && !misassembled_[0] && !misassembled_[1] )
    {
        float len = max( params.readLen, ends_[1] - ends_[0] );
        // Cutoff is 0.9x coverage, increasing to 1.3x
        float cutoff = params.cover * ( 1.4 + min( float(0.4), 
                len / float( 20 * params.readLen ) ) );
        return coverage_ < cutoff;
    }
    return false;
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
    float readCount = countReads();
    coverage_ = readCount < 3 ? params.cover : ( ( readCount - 2 ) * float( params.readLen) ) / len;
}

void Node::setReliability()
{
    if ( !reliable_ ) reliable_ = getReliability();
}

void Node::setReliable( bool doForce )
{
    reliable_ = doForce || !unreliable_;
}

bool Node::setUnreliable()
{
    if ( !unreliable_ || reliable_ )
    {
        unreliable_ = true;
        reliable_ = false;
        return true;
    }
    return false;
}
