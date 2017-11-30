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

Node::Node()
: edges_( 2 ), drxn_( 2 ), clones_( NULL ), extendCount_( 0 ), coverage_( 0 ), validated_( false ), reliable_( false ), unreliable_( false )
{
    paired_ = new NodeSet;
    stop_[0] = stop_[1] = 0;
    ends_[0] = ends_[1] = 0;
    farPairNodes_[0] = farPairNodes_[1] = NULL;
    edgeCount_[0] = edgeCount_[1] = 0;
    validLimits_[0] = validLimits_[1] = numeric_limits<int32_t>::max();
    validLimits_[2] = validLimits_[3] = numeric_limits<int32_t>::min();
    assembled_[0] = assembled_[1] = misassembled_[0] = misassembled_[1] = false;
}

Node::Node( string seq )
: Node()
{
    seq_ = seq;
}

// Seeded node
Node::Node( vector<Overlap> &reads )
: Node( reads[0].seq )
{
    ends_[0] = -reads[0].extLen;
    ends_[1] = reads[0].overLen;
    int maxOver( reads[0].overLen );
    for ( Overlap &read : reads )
    {
        addRead( read, 0, false );
        if ( read.overLen > maxOver )
        {
            seq_ += read.seq.substr( read.seq.length() -  ( read.overLen - maxOver ) );
            maxOver = read.overLen;
        }
    }
    ends_[1] = maxOver;
    validLimits_[0] = validLimits_[1] = ends_[0];
    validLimits_[2] = validLimits_[3] = ends_[1];
}

// Branch node
Node::Node( string seq, Extension &ext, int32_t prevEnd, bool drxn )
: Node( seq )
{
    drxn_ = drxn;
    ends_[!drxn] = drxn ? prevEnd - ext.maxOverLen : prevEnd + ext.maxOverLen;
    ends_[drxn] = drxn ? ends_[0] + seq.length() : ends_[1] - seq.length();
    appendNode( ext, drxn );
    setValid( !drxn );
    setCoverage();
}

// Split node
Node::Node( string seq, int32_t beginCoord, int* stop, bool drxn )
: Node( seq )
{
    drxn_ = drxn;
    stop_[0] = stop[0];
    stop_[1] = stop[1];
    ends_[!drxn] = beginCoord;
    ends_[drxn] = drxn ? ends_[0] + seq.length() : ends_[1] - seq.length();
    setValid( !drxn );
}

// Clone node
Node::Node( Node* toClone, ExtVars &ev, bool drxn )
: Node()
{
    toClone->extendComplete( ev, drxn );
    seq_ = toClone->seq_;
    reads_ = toClone->reads_;
    drxn_ = toClone->drxn_;
    ends_[0] = toClone->ends_[0];
    ends_[1] = toClone->ends_[1];
    edgeCount_[0] = toClone->edgeCount_[0];
    edgeCount_[1] = toClone->edgeCount_[1];
    toClone->addClone( this );
    resetMarks();
    coverage_ = toClone->coverage_;
}

// Seed island node
Node::Node( string seq, ReadMark &mark, Extension &ext, bool extDrxn, bool islandDrxn )
: Node( seq )
{
    drxn_ = 3 + islandDrxn;
    ends_[0] = ends_[1] = mark.estimate + ( islandDrxn ? ext.maxOverLen - params.readLen : params.readLen - ext.maxOverLen );
    ends_[!islandDrxn] = islandDrxn ? ends_[1] - seq_.length() : ends_[0] + seq_.length();
    appendNode( ext, extDrxn );
    setCoverage();
}

// Island branch
Node::Node( string seq, Extension &ext, int32_t prevEnd, bool extDrxn, bool islandDrxn )
: Node( seq )
{
    drxn_ = 3 + islandDrxn;
    ends_[!extDrxn] = extDrxn ? prevEnd - ext.maxOverLen : prevEnd + ext.maxOverLen;
    ends_[extDrxn] = extDrxn ? ends_[0] + seq.length() : ends_[1] - seq.length();
    appendNode( ext, extDrxn );
    setCoverage();
}

// Blank node
Node::Node( int32_t anchor, int32_t overlap, int subGraph, bool drxn )
: Node()
{
    drxn_ = subGraph;
    ends_[!drxn] = anchor + ( drxn ? -overlap : overlap );
    ends_[drxn] = ends_[!drxn] + ( drxn ? params.readLen : -params.readLen );
    seq_ = string( params.readLen, 'N' );
}

Node::Node( string seq, vector<Overlap> &overlaps, int drxn )
: Node( seq )
{
    drxn_ = drxn;
    ends_[1] = seq.length();
    for ( const Overlap &over : overlaps )
    {
        reads_.insert( make_pair( over.readId, Coords( overlaps[0].overLen - over.overLen, overlaps[0].overLen + over.extLen, false ) ) );
    }
}

Node::Node( MapNode* mapNode, int drxn )
: Node( mapNode->seq )
{
    drxn_ = drxn;
    ends_[0] = mapNode->estimate;
    ends_[1] = ends_[0] + seq_.length();
    for ( int i = 0; i < mapNode->ids.size(); i++ )
    {
        reads_.insert( make_pair( mapNode->ids[i], 
                Coords( ends_[0] + mapNode->coords[0][i]
                      , ends_[0] + mapNode->coords[1][i], false ) ) );
    }
    resetMarks();
}

Node::Node( ReadStruct &read )
: Node( read.seq )
{
    ends_[0] = read.coords[0];
    ends_[1] = read.coords[1];
    validLimits_[0] = validLimits_[1] = read.tether[0];
    validLimits_[2] = validLimits_[3] = read.tether[1];
    reads_.insert( make_pair( read.readId, Coords( ends_[0], ends_[1], false ) ) );
}

//void Node::nodeTest()
//{
//    NodeIntList hitScores, missScores;
//    for ( Node* fwd : getDrxnNodes( 1 ) )
//    {
//        int hits = 0;
//        int misses = 0;
//        for ( ReadMark &mark : fwd->getMarksBase( 1 ) )
//        {
//            auto it = reads_.find( mark.readId );
//            if ( it != reads_.end() )
//            {
//                if ( mark.coords[0] < it->second[0] && it->second[0] < mark.coords[1] )
//                {
//                    hits++;
//                }
//                else
//                {
//                    misses++;
//                }
//            }
//        }
//        hitScores.push_back( make_pair( fwd, hits ) );
//        missScores.push_back( make_pair( fwd, misses ) );
//    }
//    assert( false );
//}
//
//void Node::readTest()
//{
//    bool limitsGood[2] = { false, false };
//    for ( auto &read : reads_ )
//    {
//        limitsGood[0] = limitsGood[0] || read.second[0] == ends_[0];
//        limitsGood[1] = limitsGood[1] || read.second[1] == ends_[1];
//        if ( limitsGood[0] && limitsGood[1] ) break;
//    }
//    
//    assert( limitsGood[0] && limitsGood[1] );
//}
//
//void Node::offsetTest( bool drxn )
//{
//    bool didSet = false;
//    vector<int32_t> offsets;
//    for ( Edge &e : edges_[!drxn] )
//    {
//        int32_t edgeEnd = drxn ? e.node->ends_[1] - e.overlap : e.node->ends_[0] + e.overlap;
//        int32_t thisOffset = drxn ? ends_[0] - edgeEnd : edgeEnd - ends_[1];
//        offsets.push_back( thisOffset );
//        didSet = didSet || thisOffset >= 0;
//    }
//    
//    assert( didSet || edges_[!drxn].empty() );
//}

//void Node::farTest( bool drxn )
//{
//    if ( farPairNodes_[0] )
//    {
//        for ( auto &np : pairs_ )
//        {
//            if ( np.first->drxn_ <= 2 )
//            {
//                assert( drxn ? farPairCoords_[0] < np.first->ends_[1]
//                             : np.first->ends_[0] < farPairCoords_[0] );
//            }
//        }
//    }
//}

void Node::addEdge( Node* node, bool drxn, bool isLeap )
{
    int overlap = abs( ends_[drxn] - node->ends_[!drxn] );
    edges_[drxn].push_back( Edge( node, overlap, isLeap ) );
    node->edges_[!drxn].push_back( Edge( this, overlap, isLeap ) );
    edgeCount_[drxn] = max( edgeCount_[drxn], (int)edges_[drxn].size() );
    node->edgeCount_[!drxn] = max( node->edgeCount_[!drxn], (int)node->edges_[!drxn].size() );
}

void Node::addEdge( Node* node, int overlap, bool drxn, bool doOffset, bool isLeap )
{
    NodeSet bckSet = getDrxnNodes( !drxn );
    assert( bckSet.find( node ) == bckSet.end() );
    assert( overlap < 200 );
    assert( node != this );
    assert( overlap <= ends_[1] - ends_[0] );
    bool didAdd = false;
    for ( Edge &e : edges_[drxn] )
    {
        if ( node == e.node )
        {
            e.overlap = max( e.overlap, overlap );
            didAdd = true;
        }
    }
    if ( !didAdd )
    {
        edges_[drxn].push_back( Edge( node, overlap, isLeap ) );
        node->edges_[!drxn].push_back( Edge( this, overlap, isLeap ) );
    }
    if ( doOffset )
    {
        if ( drxn ? drxn_ == 0 : drxn_ == 1 )
        {
            propagateOffset( !drxn );
        }
        else
        {
            node->propagateOffset( drxn );
        }
    }
    edgeCount_[drxn] = max( edgeCount_[drxn], (int)edges_[drxn].size() );
    node->edgeCount_[!drxn] = max( node->edgeCount_[!drxn], (int)node->edges_[!drxn].size() );
}

void Node::blankEnd( int32_t len, bool drxn )
{
    int32_t unblankLen = ends_[1] - ends_[0] - len;
    int32_t i = drxn ? unblankLen : 0;
    int32_t blankEnd = drxn ? unblankLen + len : len;
    while ( i < blankEnd )
    {
        seq_[i] = 'N';
        i++;
    }
    
    for ( Edge &e : edges_[!drxn] )
    {
        if ( unblankLen < e.overlap )
        {
            e.node->blankEnd( e.overlap - unblankLen, drxn );
        }
    }
}

void Node::clearEdges( bool drxn )
{
    for ( Edge edge: edges_[drxn] )
    {
        edge.node->removeEdge( this, !drxn );
    }
    edges_[drxn].clear();
}

bool Node::deleteTest( bool drxn )
{
    NodeList tNodes = getTargetNodes( drxn, true );
    bool anyBad = false;
    for ( Node* fwd : getDrxnNodes( drxn, false, true ) )
    {
        for ( Node* t : tNodes )
        {
            int hits = 0;
            if ( t != fwd )
            {
                for ( ReadMark &mark : fwd->marks_[drxn] )
                {
                    auto it = t->reads_.find( mark.readId );
                    if ( it != t->reads_.end()
                            && !fwd->clones_
                            && mark.isValid( it->second )
                            && ( drxn ? it->second[1] <= t->validLimits_[2] 
                                      : t->validLimits_[1] <= it->second[0] ) )
                    {
                        bool isPaired = fwd->paired_ && fwd->paired_->find( t ) != fwd->paired_->end();
                        hits++;
                        anyBad = true;
                    }
                }
            }
        }
    }
    
    return anyBad;
}

bool Node::deleteTest( NodeList &tNodes, bool drxn )
{
    bool anyBad = false;
    for ( Node* t : tNodes )
    {
        int hits = 0;
        for ( ReadMark &mark : marks_[drxn] )
        {
            auto it = t->reads_.find( mark.readId );
            if ( it != t->reads_.end()
                    && !clones_
                    && mark.isValid( it->second )
                    && ( drxn ? it->second[1] <= t->validLimits_[2] 
                              : t->validLimits_[1] <= it->second[0] ) )
            {
                hits++;
                anyBad = true;
            }
        }
    }
    
    return anyBad;
}

int32_t Node::getBiggestOffset( bool drxn )
{
    int32_t biggest = 0;
    for ( Edge &e : edges_[drxn] )
    {
        biggest = max( biggest, abs( ( drxn ? ends_[1] - e.node->ends_[0] : e.node->ends_[1] - ends_[0] ) - e.overlap ) );
    }
    return biggest;
}

int Node::getBestOverlap( bool drxn )
{
    int overlap = 0;
    for ( Edge &edge : edges_[drxn] )
    {
        if ( edge.overlap > overlap )
        {
            overlap = edge.overlap;
        }
    }
    return overlap;
}

int Node::getEdgeViableCount( bool drxn )
{
    int edgeCount = 0;
    for ( Edge &e : edges_[drxn] )
    {
        edgeCount += ( e.node->stop_[drxn] == 0 );
    }
    return edgeCount;
}

int32_t Node::getFurthest( int32_t q, int32_t t, bool drxn )
{
    if ( drxn )
    {
        return max( q, t );
    }
    else
    {
        return min( q, t );
    }
}

string Node::getHeader( string header )
{
//    int digits = max( 3, (int)locusId.length() );
//    string header = ">Locus" + string( digits - locusId.length(), '0' ) + locusId;
    int digits = max( 3, (int)id_.length() );
    header += "|Node_" + string( digits - id_.length(), '0' ) + id_;
    for ( int drxn : { 0 , 1 } )
    {
        header += "|";
        bool first = true;
        for ( Edge edge : edges_[drxn] )
        {
            if ( !first )
            {
                header += "_";
            }
            header += edge.node->id_ + "x" + to_string( edge.overlap );
            first = false;
        }
    }
    return header;
}

int Node::getOverlap( Node* node, bool drxn )
{
    int overlap = 0;
    for ( Edge &edge : edges_[drxn] )
    {
        if ( edge.node == node )
        {
            overlap = edge.overlap;
            break;
        }
    }
    return overlap;
}

string Node::getSeqEnd( int len, bool drxn )
{
    assert ( len > 0 && len <= seq_.length() );
    return ( drxn ? seq_.substr( seq_.length() - len ) : seq_.substr( 0, len ) );
}

NodeList Node::getTargetNodes( bool drxn, bool inclSelf, bool inclClones )
{
    NodeSet tSet;
    if ( inclSelf )
    {
        tSet.insert( this );
    }
    getDrxnNodes( tSet, !drxn, params.getFurthestMpDist( ends_[!drxn], !drxn ) );
    if ( inclClones && clones_ )
    {
        for ( Node* clone : *clones_ )
        {
            NodeSet cloneTSet;
            clone->getDrxnNodes( cloneTSet, !drxn, params.getFurthestMpDist( clone->ends_[!drxn], !drxn ) );
            tSet.insert( cloneTSet.begin(), cloneTSet.end() );
            if ( inclSelf )
            {
                tSet.insert( clone );
            }
        }
    }
    NodeList tNodes( tSet.begin(), tSet.end() );
    return tNodes;
}

bool Node::isBeyond( int32_t bgn, int32_t nd, bool drxn )
{
    if ( drxn )
    {
        return ends_[0] > bgn && ends_[1] > nd;
    }
    else
    {
        return ends_[0] < bgn && ends_[1] < nd;
    }
}

bool Node::isContinue(bool drxn)
{
    return edges_[drxn].empty() && stop_[drxn] == 0;
}

bool Node::isDeadEnd( bool drxn )
{
    if ( edges_[drxn].empty() && stop_[drxn] && !isContinue( !drxn ) )
    {
        for ( auto &read : reads_ )
        {
            if ( read.second[drxn] != ends_[drxn] )
            {
                return false;
            }
        }
        return true;
    }
    return false;
}

bool Node::isEdge( Node* node, bool drxn )
{
    for ( Edge &e : edges_[drxn] )
    {
        if ( e.node == node )
        {
            return true;
        }
    }
    return false;
}

bool Node::isEnded(bool drxn)
{
    return ( !edges_[drxn].empty() || stop_[drxn] > 0 );
}

bool Node::isFurther( int32_t coord, bool endDrxn, bool drxn )
{
    if ( drxn )
    {
        return ends_[endDrxn] > coord;
    }
    else
    {
        return ends_[endDrxn] < coord;
    }
}

bool Node::isNearer( int32_t coord, bool endDrxn, bool drxn )
{
    if ( drxn )
    {
        return ends_[endDrxn] < coord;
    }
    else
    {
        return ends_[endDrxn] > coord;
    }
}

bool Node::pause( bool drxn )
{
    if ( isContinue( drxn ) )
    {
        stop_[drxn] = -1;
        return true;
    }
    return false;
}


void Node::propagateOffset( bool drxn )
{
    NodeSet propagated;
    propagateOffset( propagated, drxn );
}

void Node::propagateOffset( NodeSet &propagated, bool drxn )
{
    propagated.insert( this );
    if ( offsetNode( drxn ) && drxn_ != 2 )
    {
        for ( Edge &edge : edges_[drxn] )
        {
            if ( propagated.find( edge.node ) == propagated.end() )
            {
                edge.node->propagateOffset( propagated, drxn );
            }
        }
    }
}

void Node::removeEdge( Node* node, bool drxn )
{
    for ( auto it = edges_[drxn].begin(); it != edges_[drxn].end(); )
    {
        if ( it->node == node )
        {
            it = edges_[drxn].erase( it );
        }
        else
        {
            it++;
        }
    }
}

bool Node::removeEdges( NodeSet &removeSet, bool drxn )
{
    bool rVal = false;
    for ( auto it = edges_[drxn].begin(); it != edges_[drxn].end(); )
    {
        if ( removeSet.find( (*it).node ) != removeSet.end() )
        {
            (*it).node->removeEdge( this, !drxn );
            it = edges_[drxn].erase( it );
            rVal = true;
            continue;
        }
        it++;
    }
    if ( rVal && edges_[drxn].empty() )
    {
        stop_[drxn] = 1;
    }
    
    return rVal;
}

void Node::stop( int stopCode, bool drxn )
{
    stop_[drxn] = stopCode;
}

void Node::trimEnd( bool drxn )
{
    int32_t mark = ends_[!drxn] + ( drxn ? params.readLen : -params.readLen );
    if ( anyValid( drxn ) )
    {
        mark = getValidLimit( drxn );
    }
    mark = drxn ? min( mark, ends_[1] ) : max( mark, ends_[0] );
    trimSeq( mark, drxn );
}

void Node::trimSeq( int32_t coord, bool drxn, bool stop )
{
    coord = drxn ? min( coord, ends_[1] ) : max( coord, ends_[0] );
    seq_ = drxn ? seq_.substr( 0, max( coord, ends_[0] ) - ends_[0] ) : seq_.substr( min( coord, ends_[1] ) - ends_[0] );
    trimReads( coord, drxn );
    if ( stop )
    {
        stop_[drxn] = 1;
    }
}

bool Node::unpause(bool drxn)
{
    if ( stop_[drxn] == -1 )
    {
        stop_[drxn] = 0;
    }
    return isContinue( drxn );
}

void Node::dismantleNode()
{
    this->reads_.clear();
    this->stop_[0] = 3;
    this->stop_[1] = 3;
    for ( int drxn : { 0, 1 } )
    {
        for ( Edge edge : edges_[drxn] )
        {
            edge.node->removeEdge( this, !drxn );
            if ( edge.node->edges_[!drxn].empty() )
            {
                edge.node->stop_[!drxn] = 3;
            }
        }
        edges_[drxn].clear();
    }
    for ( const pair<Node*, int> &pn : pairs_ )
    {
        if ( pn.first != this )
        {
            pn.first->pairs_.erase( this );
            pn.first->resetFurthest( this );
        }
    }
    if ( clones_ )
    {
        for ( Node* clone : *clones_ )
        {
            clone->removeClone( this );
        }
        delete clones_;
        clones_ = NULL;
    }
    if ( paired_ )
    {
        delete paired_;
        paired_ = NULL;
    }
    pairs_.clear();
}

void Node::dismantleNode( NodeSet &delSet, bool drxn )
{
    this->reads_.clear();
    this->stop_[0] = 3;
    this->stop_[1] = 3;
    for ( const pair<Node*, int> &pn : pairs_ )
    {
        if ( pn.first != this )
        {
            pn.first->pairs_.erase( this );
            pn.first->resetFurthest( this );
        }
    }
    pairs_.clear();
    assert( drxn_ <= 2 );
    for ( Edge edge : edges_[!drxn] )
    {
        assert( edge.node->drxn_ <= 2 );
        edge.node->removeEdge( this, drxn );
        if ( edge.node->edges_[drxn].empty() )
        {
            edge.node->stop_[drxn] = 3;
        }
    }
    for ( Edge edge : edges_[drxn] )
    {
        assert( edge.node->drxn_ <= 2 );
        edge.node->removeEdge( this, !drxn );
        if ( edge.node->edges_[!drxn].empty() )
        {
            edge.node->dismantleNode( delSet, drxn );
        }
        else
        {
            edge.node->propagateOffset( drxn );
        }
    }
    
    if ( clones_ )
    {
        for ( Node* clone : *clones_ )
        {
            clone->removeClone( this );
        }
        delete clones_;
        clones_ = NULL;
    }
    if ( paired_ )
    {
        delete paired_;
        paired_ = NULL;
    }
    
    this->edges_[0].clear();
    this->edges_[1].clear();
    delSet.insert( this );
}
