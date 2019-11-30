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

void Node::add( ReadId id, int i, int j, bool redundant, bool ignore, int unaligned )
{
    Lib* lib = params.getLib( id );
    Coords coords( i, j, redundant );
    if ( ignore ) coords.ignore = true;
    if ( unaligned ) coords.coords[2] = unaligned;
    
    auto r = reads_.insert( make_pair( id, coords ) );
    if ( !lib ) return;
    bool drxn = lib->getPair( id );
    if ( !r.second )
    {
        rmvMark( id, drxn );
        r.first->second = coords;
    }
    if ( !r.first->second.ignore ) ( lib->isPe ? pe_ : mp_ )[drxn].push_back( NodeMark( id, i, j, lib, drxn ) );
    culled_ = false;
}

bool Node::add( ReadId id, string& seq )
{
    assert( !seq.empty() );
    size_t it = seq_.find( seq );
    if ( it == seq_.npos ) return false;
    int32_t coords[2]{ int32_t( it + ends_[0] ), int32_t( it + ends_[0] + seq.size() ) };
    add( id, coords[0], coords[1], isRedundant( coords[0], coords[1] ) );
    return true;
}

//void Node::addMark( ReadId id, int i, int j )
//{
//    Lib* lib = params.getLib( id );
//    if ( !lib ) return;
//    bool drxn = lib->getPair( id );
//    ( lib->isPe ? pe_ : mp_ )[drxn].push_back( NodeMark( id, i, j, lib, drxn ) );
//}

void Node::addMark( SeqNum readId, Coords &coords )
{
    Lib* lib = params.getLib( readId );
    int drxn;
    if ( lib && (*lib).getPair( readId, drxn ) )
    {
        marks_[!drxn].push_back( ReadMark( readId, coords, lib, drxn ) );
    }
}

void Node::addMarks( vector<ReadId> &ids )
{
    Coords* hit = NULL;
    for ( ReadId &id : ids ) if ( findRead( id, hit ) ) add( id, hit->coords[0], hit->coords[1], hit->redundant );
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

bool Node::addRead( SeqNum readId, int32_t bgn, int32_t nd, bool isRedundant, bool notClones )
{
    if ( clones_ && notClones ) return false;
    Coords coords( bgn, nd, isRedundant );
    auto r = reads_.insert( make_pair( readId, coords ) );
    if ( !r.second )
    {
        removeMark( readId );
        r.first->second = coords;
    }
    addMark( readId, coords );
    if ( clones_ )
    {
        for ( Node* clone : *clones_ )
        {
            int32_t offset = clone->ends_[0] - ends_[0];
            Coords cloneCoords( bgn + offset, nd + offset, isRedundant );
            r = clone->reads_.insert( make_pair( readId, cloneCoords ) );
            if ( !r.second )
            {
                clone->removeMark( readId );
                r.first->second = cloneCoords;
            }
            clone->addMark( readId, cloneCoords );
        }
    }
    
    return true;
}

void Node::addRead( NodeMapRead &mapRead, bool drxn )
{
    Coords coords( mapRead.coords[drxn][0]
                 , mapRead.coords[drxn][1]
                 , mapRead.seq.length() > mapRead.coords[drxn][1] - mapRead.coords[drxn][0] );
    auto r = reads_.insert( make_pair( mapRead.id, coords ) );
    if ( !r.second )
    {
        removeMark( mapRead.id );
        r.first->second = coords;
    }
    addMark( mapRead.id, coords );
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

int Node::countReads( bool notRedundant )
{
    if ( !notRedundant ) return reads_.size();
    
    int readCount = 0;
    for ( auto &read : reads_ ) if ( !read.second.redundant ) readCount++;
    assert( readCount );
    return readCount;
}

void Node::cullMarks()
{
    if ( culled_ || isContinue( 0 ) || isContinue( 1 ) ) return;
    for ( int d = 0; d < 2; d++ )
    {
        cullMarks( pe_[d], d );
        cullMarks( mp_[d], d );
    }
    culled_ = true;
    if ( cloned_ ) for ( Node* clone : cloned_->nodes ) clone->culled_ = true;
}

void Node::cullMarks( vector<NodeMark>& marks, bool drxn )
{
    int32_t limits[2]{ ends_[0] + getBestOverlap( 0 ), ends_[1] - getBestOverlap( 1 ) };
    setVerifyLimits( limits );
    
    for ( int i = 0; i < marks.size(); i++ )
    {
        if ( limits[0] < marks[i].coords[1] && marks[i].coords[0] < limits[1] ) continue;
        if ( cloned_ ) for ( Node* clone : cloned_->nodes ) clone->rmvMark( marks[i].id, drxn );
        marks.erase( marks.begin() + i-- );
    }
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

bool Node::findRead( ReadId id, Coords*& hit, bool inclRedundant )
{
    auto it = reads_.find( id );
    if ( it == reads_.end() || ( !inclRedundant && it->second.redundant ) ) return false;
    hit = &it->second;
    return true;
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

vector<NodeMark> Node::getMarks( bool drxn, bool pe )
{
    vector<NodeMark> marks;
    for ( auto &read : reads_ )
    {
        ReadId id = read.first;
        Lib* lib = params.getLib( id );
        if ( !lib || lib->isPe != pe ) continue;
        bool d = lib->getPair( id );
        if ( d == drxn ) marks.push_back( NodeMark( id, read.second[0], read.second[1], lib, d ) );
    }
    return marks;
}

vector<ReadMark> Node::getMarksBase( int drxn )
{
    vector<ReadMark> marks;
    for ( auto &read : reads_ )
    {
        SeqNum pairId = read.first;
        Lib* lib = params.getLib( pairId );
        int pairDrxn;
        if ( lib && (*lib).getPair( pairId, pairDrxn ) && pairDrxn != drxn )
        {
            marks.push_back( ReadMark( pairId, read.second, lib, pairDrxn ) );
        }
    }
    sortMarks( marks, drxn );
    return marks;
}

void Node::getMarksCount( int counts[2] )
{
    for ( auto &read : reads_ )
    {
        SeqNum pairId = read.first;
        Lib* lib = params.getLib( pairId );
        int drxn;
        if ( lib && (*lib).getPair( pairId, drxn ) )
        {
            counts[drxn]++;
        }
    }
}

Coords* Node::getRead( ReadId id )
{
    auto it = reads_.find( id );
    return it != reads_.end() ? &it->second : NULL;
}

bool Node::getSplitCoords( int32_t coords[2], int split, bool drxn )
{
    coords[0] = ends_[0];
    coords[1] = ends_[1];
    
    for ( auto &read : reads_ )
    {
        if ( read.second.redundant ) continue;
        
        // Near side
        if ( drxn ? read.second[0] < split : split < read.second[1] )
        {
            if ( drxn ? coords[0] < read.second[1] : read.second[0] < coords[1] ) coords[!drxn] = read.second[drxn];
        }
        // Far side
        else if ( drxn ? read.second[0] < coords[1] : coords[0] < read.second[1] ) coords[drxn] = read.second[!drxn];
    }
    
    return coords[0] != ends_[0] && coords[1] != ends_[1];
}

bool Node::isRedundant( int i, int j )
{
    int len = j - i;
    if ( len == params.readLen ) return false;
    for ( auto &r : reads_ ) if ( len < r.second[1]-r.second[0] && r.second[0] <= i && j <= r.second[1] ) return true;
    
    return false;
}

bool Node::isRedundant( Coords* coords )
{
    if ( coords->redundant
            || (*coords)[1] <= ends_[0] + getBestOverlap( 0 ) 
            || ends_[1] - getBestOverlap( 1 ) <= (*coords)[0] ) return true;
    int len = (*coords)[1] - (*coords)[0];
    for ( auto &read : reads_ )
    {
        if ( len >= read.second[1] - read.second[0] ) continue;
        if ( read.second[0] <= (*coords)[0] && (*coords)[1] <= read.second[1] ) return true;
    }
    
    return false;
}

bool Node::offsetNode( bool drxn )
{
    int32_t off = 0;
    bool first = true;
    for ( Edge &e : edges_[!drxn] )
    {
        int32_t edgeOffset = e.node->ends_[drxn] - ends_[!drxn] + ( drxn ? -e.ol : e.ol );
        off = first ? edgeOffset : ( drxn ? max( off, edgeOffset ) : min( off, edgeOffset ) );
        first = false;
    }
    
    return offset( off );
}

void Node::reAddMark( SeqNum readId, Coords &coords )
{
    Lib* lib = params.getLib( readId );
    int drxn;
    if ( lib && (*lib).getPair( readId, drxn ) )
    {
        if ( find_if( marks_[drxn].begin(), marks_[drxn].end(), [&readId]( const ReadMark &a ){ 
            return a.id == readId;
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

void Node::requery( Querier& bwt )
{
    vector<ReadId> ids;
    vector<int32_t> coords[2];
    bwt.mapSequence( seq_, ids, coords );
    for ( int i = 0; i < ids.size(); i++ )
    {
        if ( reads_.find( ids[i] ) != reads_.end() ) continue;
        coords[0][i] += ends_[0];
        coords[1][i] += ends_[0];
        add( ids[i], coords[0][i], coords[1][i], isRedundant( coords[0][i], coords[1][i] ) );
    }
}

void Node::remark()
{
    pe_[0].clear();
    pe_[1].clear();
    mp_[0].clear();
    mp_[1].clear();
    for ( auto &read : reads_ ) if ( !read.second.ignore && !read.second.unpaired )
    {
        ReadId id = read.first;
        Lib* lib = params.getLib( id );
        if ( !lib ) continue;
        bool drxn = lib->getPair( id );
        ( lib->isPe ? pe_ : mp_ )[drxn].push_back( NodeMark( id, read.second[0], read.second[1], lib, drxn ) );
    }
    culled_ = false;
    resort();
}

void Node::removeMark( SeqNum &readId )
{
    SeqNum pairId = params.getPairId( readId );
    for ( bool drxn : { 0, 1 } )
    {
        for ( auto it = marks_[drxn].begin(); it != marks_[drxn].end(); )
        {
            if ( it->id == pairId )
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
            SeqNum readId = isPair ? params.getPairId( mark.id ) : mark.id;
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

void Node::resetUnmarked( bool drxn )
{
    if ( getPairHitsTotal() < marks_[0].size() + marks_[1].size() )
    {
        NodeList nodes;
        for ( auto &np : pairs_ ) nodes.push_back( np.first );
        unordered_set<ReadId> ids;
        for ( ReadMark const &mark : marks_[drxn] ) ids.insert( mark.id );
        for ( ReadMark const &mark : getMarksBase( drxn ) )
        {
            bool doAdd = true;
            if ( ids.find( mark.id ) != ids.end() ) continue;
            for ( Node* t : nodes ) if ( t->reads_.find( mark.id ) != t->reads_.end() ) doAdd = false;
            if ( doAdd ) marks_[drxn].push_back( mark );
        }
        sortMarks( marks_[drxn], drxn );
    }
}

bool Node::rmvMark( ReadId id, bool drxn )
{
    vector<NodeMark> &marks = ( params.isReadPe( id ) ? pe_[drxn] : mp_[drxn] );
    for ( auto it = marks.begin(); it != marks.end(); it++ ) if ( it->id == id )
    {
        marks.erase( it );
        return true;
    }
    return false;
}

void Node::resort()
{
    int32_t limits[2]{ ends_[0] + getBestOverlap( 0, true ), ends_[1] - getBestOverlap( 1, true ) };
    for ( vector<NodeMark>* marks : { pe_, mp_ } )
    {
        sort( marks[0].begin(), marks[0].end(), []( NodeMark &a, NodeMark &b ){ return a.coords[1] > b.coords[1]; } );
        sort( marks[1].begin(), marks[1].end(), []( NodeMark &a, NodeMark &b ){ return a.coords[0] < b.coords[0]; } );
        for ( int d : { 0, 1 } )
        {
            while ( !marks[d].empty() && ( d ? marks[1][0].coords[1] <= limits[0]
                                              : limits[1] <= marks[0][0].coords[0] ) ) marks[d].erase( marks[d].begin() );
            while ( !marks[d].empty() && ( d ? limits[1] <= marks[1].back().coords[0]
                                             : marks[0].back().coords[1] <= limits[0] ) ) marks[d].pop_back();
        }
    }
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

int32_t Node::split( Node* node, int32_t cut, bool drxn )
{
    int32_t ext = cut;
    for ( auto it = reads_.begin(); it != reads_.end(); )
    {
        if ( drxn ? cut <= it->second[0] : it->second[1] <= cut )
        {
            node->reads_.insert( *it );
            it = reads_.erase( it );
        }
        else
        {
            ext = drxn ? max( ext, it->second[1] ) : min( ext, it->second[0] );
            it++;
        }
    }
    return ext;
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
            if ( removeIds.find( mark.id ) == removeIds.end() )
            {
                marks[markDrxn].push_back( mark );
            }
        }
    }
    reads_ = reads;
    marks_[0] = marks[0];
    marks_[1] = marks[1];
}
