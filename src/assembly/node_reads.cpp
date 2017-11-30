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

void Node::addMatchedReads( vector<Overlap> &reads )
{
    assert( false );
}

void Node::addMark( SeqNum readId, Coords &coords )
{
    Lib* lib = params.getLib( readId );
    int drxn;
    int32_t dist;
    if ( lib && (*lib).getPair( readId, dist, drxn ) )
    {
        marks_[!drxn].push_back( ReadMark( readId, coords, lib, drxn ) );
    }
}

void Node::addRead( Overlap &read, int32_t anchor, bool olDrxn )
{
    Coords coords( anchor - ( olDrxn ? read.overLen : read.extLen ), 
                   anchor + ( olDrxn ? read.extLen : read.overLen ), 
                   read.redundant );
    auto r = reads_.insert( make_pair( read.readId, coords ) );
    if ( !r.second )
    {
        removeMark( read.readId );
        r.first->second = coords;
    }
    addMark( read.readId, coords );
}

void Node::addRead( SeqNum readId, int32_t bgn, int32_t nd, bool isRedundant )
{
    Coords coords( bgn, nd, isRedundant );
    auto r = reads_.insert( make_pair( readId, coords ) );
    if ( !r.second )
    {
        removeMark( readId );
        r.first->second = coords;
    }
    addMark( readId, coords );
}

bool Node::anyReadBeyondCoord( int32_t coord, bool coordDrxn, bool drxn )
{
    for ( auto &read : reads_ )
    {
        if ( drxn ? coord < read.second[coordDrxn]
                  : read.second[coordDrxn] < coord )
        {
            return true;
        }
    }
    return false;
}

bool Node::anyReadInNode( unordered_set<SeqNum> &readIds )
{
    for ( const SeqNum &readId : readIds )
    {
        if ( reads_.find( readId ) != reads_.end() )
        {
            return true;
        }
    }
    return false;
}

int32_t Node::findNextRead( int32_t mark, bool drxn )
{
    int32_t nxt = ends_[drxn];
    if ( drxn )
    {
        for ( auto &read : reads_ )
        {
            if ( mark < read.second[0]  )
            {
                nxt = min( read.second[0], nxt );
            }
        }
    }
    else
    {
        for ( auto &read : reads_ )
        {
            if ( read.second[1] < mark )
            {
                nxt = max( read.second[1], nxt );
            }
        }
    }
    return nxt;
}

bool Node::findOverlap( Node* &hitNode, int32_t* coords, string &seq, NodeList &nodes, int minOl, bool drxn )
{
    if ( seq.length() < minOl ) return false;
    string q = drxn ? seq.substr( seq.length() - minOl ) : seq.substr( 0, minOl );
    
    bool found = false;
    for ( Node* node : nodes )
    {
        if ( minOl >= seq.length() ) break;
        size_t it = node->seq_.find( q );
        while ( it != node->seq_.npos )
        {
            int overLen = minOl;
            int iRead = drxn ? seq.length() - minOl : minOl - 1;
            int iNode = drxn ? it : node->seq_.length() - it + minOl - 1;
            if ( drxn )
            {
                while ( iRead-- > 0 && iNode-- > 0 && seq[iRead] == node->seq_[iNode] ) overLen++;
            }
            else
            {
                while ( ++iRead < seq.length() && ++iNode < node->seq_.length() && seq[iRead] == node->seq_[iNode] ) overLen++;
            }
            hitNode = node;
            found = true;
            coords[0] = drxn ? node->ends_[0] + it + minOl - overLen : node->ends_[0] + it;
            coords[1] = coords[0] + overLen;
            minOl = overLen + 1;
            if ( minOl >= seq.length() ) break;
            q = drxn ? seq.substr( seq.length() - minOl ) : seq.substr( 0, minOl );
            it = node->seq_.find( q );
        }
    }
    
    return found;
}

bool Node::findRead( SeqNum &readId, Coords *&coords, bool inclRedundant )
{
    auto hit = reads_.find( readId );
    if ( hit != reads_.end() && ( inclRedundant || !hit->second.redundant ) )
    {
        coords = &hit->second;
        return true;
    }
    return false;
}

int Node::getEndMarks( bool drxn )
{
    int endMisses = 0;
    for ( ReadMark &mark : marks_[drxn] )
    {
        if ( drxn ? mark.mark > validLimits_[3] : mark.mark < validLimits_[0] )
        {
            endMisses++;
        }
    }
    return endMisses;
}

vector<ReadMark> Node::getMarksBase( int drxn )
{
    vector<ReadMark> marks;
    for ( auto &read : reads_ )
    {
        SeqNum pairId = read.first;
        Lib* lib = params.getLib( pairId );
        int pairDrxn;
        int32_t dist;
        if ( lib && (*lib).getPair( pairId, dist, pairDrxn ) && pairDrxn != drxn )
        {
            marks.push_back( ReadMark( pairId, read.second, lib, pairDrxn ) );
        }
    }
    sortMarks( marks, drxn );
    return marks;
}

bool Node::offsetNode( bool drxn )
{
    int32_t off = 0;
    bool first = true;
    for ( Edge &e : edges_[!drxn] )
    {
        int32_t edgeOffset = e.node->ends_[drxn] - ends_[!drxn] + ( drxn ? -e.overlap : e.overlap );
        off = first ? edgeOffset : ( drxn ? max( off, edgeOffset ) : min( off, edgeOffset ) );
        first = false;
    }
    
    return offset( off );
}

void Node::reAddMark( SeqNum readId, Coords &coords )
{
    Lib* lib = params.getLib( readId );
    int drxn;
    int32_t dist;
    if ( lib && (*lib).getPair( readId, dist, drxn ) )
    {
        if ( find_if( marks_[drxn].begin(), marks_[drxn].end(), [&readId]( const ReadMark &a ){ 
            return a.readId == readId;
        } ) == marks_[drxn].end() )
        {
            marks_[!drxn].push_back( ReadMark( readId, coords, lib, drxn ) );
        }
    }
}

void Node::reAddMarks( vector<SeqNum> &readIds )
{
    for ( SeqNum &readId : readIds )
    {
        auto it = reads_.find( readId );
        if ( it != reads_.end() )
        {
            for ( Node* clone : this->getCloneSet( true ) )
            {
                clone->reAddMark( readId, it->second );
            }
        }
    }
}

void Node::removeMark( SeqNum &readId )
{
    SeqNum pairId = params.getPairId( readId );
    for ( bool drxn : { 0, 1 } )
    {
        for ( auto it = marks_[drxn].begin(); it != marks_[drxn].end(); )
        {
            if ( it->readId == pairId )
            {
                it = marks_[drxn].erase( it );
                return;
            }
            it++;
        }
    }
}

void Node::removeMarks( unordered_set<SeqNum> &readIds, bool pushLimits, bool isPair, bool drxn )
{
    if ( !readIds.empty() )
    {
        vector<ReadMark> newMarks;
        for ( ReadMark &mark : marks_[drxn] )
        {
            SeqNum readId = isPair ? params.getPairId( mark.readId ) : mark.readId;
            if ( readIds.find( readId ) == readIds.end() )
            {
                newMarks.push_back( mark );
            }
            else if ( pushLimits && !validated_ )
            {
                pushValidLimits( mark.mark, drxn );
            }
        }
        if ( newMarks.size() < marks_[drxn].size() )
        {
            marks_[drxn] = newMarks;
        }
    }
}

void Node::resetMarks()
{
    marks_[0] = getMarksBase( 0 );
    marks_[1] = getMarksBase( 1 );
}

void Node::sortMarks( vector<ReadMark> &marks, bool drxn )
{
    if ( drxn )
    {
        sort( marks.begin(), marks.end(), []( const ReadMark &a, const ReadMark &b ){
            return a.mark < b.mark;
        });
    }
    else
    {
        sort( marks.begin(), marks.end(), []( const ReadMark &a, const ReadMark &b ){
            return a.mark > b.mark;
        });
    }
}

int32_t Node::splitReads( Node* node, int32_t splitBegin, bool drxn )
{
    // Split reads
    ends_[drxn] = ends_[!drxn];
    ReadCoords selfReads;
    for ( auto &read : reads_ )
    {
        if ( drxn ? read.second[0] >= splitBegin : read.second[1] <= splitBegin )
        {
            node->reads_.insert( read );
        }
        else
        {
            selfReads.insert( read );
            if ( drxn ? read.second[1] > ends_[1] : read.second[0] < ends_[0] )
            {
                ends_[drxn] = read.second[drxn];
            }
        }
    }
    reads_ = selfReads;
    
    validLimits_[0] = max( validLimits_[1], ends_[0] );
    validLimits_[1] = max( validLimits_[1], ends_[0] );
    validLimits_[2] = min( validLimits_[2], ends_[1] );
    validLimits_[3] = min( validLimits_[2], ends_[1] );
    
    resetMarks();
    node->resetMarks();
    
    return drxn ? ends_[1] - splitBegin : splitBegin - ends_[0];
}

void Node::trimReads( int32_t endCoord, bool drxn )
{
    ends_[drxn] = endCoord;
    validLimits_[0] = max( validLimits_[0], ends_[0] );
    validLimits_[1] = max( validLimits_[1], ends_[0] );
    validLimits_[2] = min( validLimits_[2], ends_[1] );
    validLimits_[3] = min( validLimits_[3], ends_[1] );
    
    ReadCoords reads;
    vector<ReadMark> marks[2];
    unordered_set<SeqNum> removeIds;
    
    for ( auto &read : reads_ )
    {
        if ( read.second[0] >= ends_[0] && read.second[1] <= ends_[1] )
        {
            reads.insert( read );
        }
        else
        {
            removeIds.insert( read.first );
            removeIds.insert( params.getPairId( read.first ) );
        }
    }
    
    for ( bool markDrxn : { 0, 1 } )
    {
        for ( ReadMark &mark : marks_[markDrxn] )
        {
            if ( removeIds.find( mark.readId ) == removeIds.end() )
            {
                marks[drxn].push_back( mark );
            }
        }
    }
    reads_ = reads;
    marks_[0] = marks[0];
    marks_[1] = marks[1];
}
