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
: allele_( NULL ), drxn_( 2 ), hits_( this ), clones_( NULL ), extendCount_( 0 ), coverage_( 0 ), validated_( false ), reliable_( false ), unreliable_( false ), verified_( false ), dontExtend_( false ), bad_( false ), mapped_( false ), branch_( false ), culled_( false )
{
    cloned_ = NULL;
    paired_ = new NodeSet;
    pairedNodes_ = NULL;
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

Node::Node( string seq, int32_t coord, int subGraph )
: Node( seq )
{
    drxn_ = subGraph;
    ends_[0] = coord;
    ends_[1] = coord + seq.size();
    if ( subGraph < 2 ) ends_.init( ends_[subGraph] );
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
    readTest();
}

Node::Node( string seq, QueryNode* ext, int32_t prevEnd, bool drxn, int graph )
: Node( seq )
{
    drxn_ = graph;
    ends_[0] = prevEnd - ( drxn ? seq.size() : 0 );
    ends_[1] = prevEnd + ( drxn ? 0 : seq.size() );
    appendNode( ext, drxn );
    setValid( !drxn );
    setCoverage();
    ends_.init( ends_[!drxn] );
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

Node::Node( Node* toClone )
: Node()
{
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

Node::Node( Node* toClone, NodeRoll& nodes, int graph, bool bad )
: Node( toClone->seq_ )
{
    cloned_ = new NodeRoll( toClone );
    if ( toClone->cloned_ )
    {
        for ( Node* clone : toClone->cloned_->nodes )
        {
            *cloned_ += clone;
            *clone->cloned_ += this;
        }
        *toClone->cloned_ += this;
    }
    else toClone->cloned_ = new NodeRoll( this );
    drxn_ = graph;
    bad_ = bad;
    mapped_ = toClone->mapped_;
    ends_[0] = toClone->ends_[0];
    ends_[1] = toClone->ends_[1];
    ends_.origin[0] = toClone->ends_.origin[0];
    ends_.origin[1] = toClone->ends_.origin[1];
    reads_ = toClone->reads_;
    coverage_ = toClone->coverage_;
    nodes += this;
    remark();
}

Node::Node( string& seq, NodeMark& mark, bool graph )
: Node( seq )
{
    drxn_ = graph;
    ends_[0] = mark.tar[0];
    ends_[1] = ends_[0] + seq.size();
    reads_.insert( make_pair( mark.id, Coords( ends_[0], ends_[1], false ) ) );
    bad_ = true;
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
    ends_.init( read.tether[0], 0 );
    ends_.init( read.tether[1], 1 );
    add( read.readId, read.coords[0], read.coords[1], false );
    readTest();
    setOrigin();
}

Node::Node( NodeMapRead &mapRead, bool drxn )
: Node( mapRead.seq )
{
    drxn_ = drxn;
    ends_[!drxn] = mapRead.coords[!drxn][!drxn];
    ends_[drxn] = drxn ? ends_[0] + seq_.length() : ends_[1] - seq_.length();
    reads_.insert( make_pair( mapRead.id, Coords( ends_[0], ends_[1], false ) ) );
    addEdge( mapRead.nodes[0], mapRead.coords[0][1] - mapRead.coords[0][0], 0 );
    addEdge( mapRead.nodes[1], mapRead.coords[1][1] - mapRead.coords[1][0], 1 );
    mapRead.nodes[0]->removeEdge( mapRead.nodes[1], 1 );
    mapRead.nodes[1]->removeEdge( mapRead.nodes[0], 0 );
}

Node::Node( MapNode* mn, int i, int j, int drxn )
: Node()
{
    mn->setRedundant();
    drxn_ = drxn;
    ends_[0] = ends_[1] = mn->coords[0][i];
    while ( i <= j )
    {
        reads_.insert( make_pair( mn->ids[i], Coords( mn->coords[0][i], mn->coords[1][i], mn->redundant[i] ) ) );
        ends_[0] = min( ends_[0], mn->coords[0][i] );
        ends_[1] = max( ends_[1], mn->coords[1][i] );
        i++;
    }
    seq_ = mn->seq.substr( ends_[0], ends_[1] - ends_[0] );
    resetMarks();
    setCoverage();
}

Node::Node( string seq, ReadId id, int32_t estimate, int drxn )
: Node( seq )
{
    drxn_ = drxn;
    ends_[0] = ends_[1] = estimate;
    if ( drxn == 1 || drxn == 4 ) ends_[1] = ends_[0] + seq.length();
    else ends_[0] = ends_[1] -seq.length();
    reads_.insert( make_pair( id, Coords( ends_[0], ends_[1], false ) ) );
    resetMarks();
    setCoverage();
}

Node::Node( string seq, ReadId id, int32_t estimate, int drxn, bool bad )
: Node( seq )
{
    drxn_ = drxn;
    ends_[0] = estimate;
    ends_[1] = estimate + seq.size();
    ends_.init( ends_[!drxn] );
    add( id, ends_[0], ends_[1], false );
    bad_ = bad;
}

Node::~Node()
{
    for ( int d : { 0, 1 } ) for ( Edge &e : edges_[d] ) e.node->removeEdge( this, !d );
    if ( paired_ ) delete paired_;
    if ( clones_ ) delete clones_;
    if ( pairedNodes_ ) delete pairedNodes_;
}

void Node::edgeTest( bool drxn )
{
    if ( edges_[drxn].size() < 2 ) return;
    vector<string> seqs;
    for ( Edge& e : edges_[drxn] )
    {
        if ( e.ol < 0 || e.ol > e.node->size() ) continue;
        seqs.push_back( drxn ? e.node->seq_.substr( e.ol ) : string( e.node->seq_.rbegin() + e.ol, e.node->seq_.rend() ) );
    }
    for ( int i = 0; i+1 < seqs.size(); i++ )
    {
        for ( int j = i+1; j < seqs.size(); j++ )
        {
            bool match = !edges_[drxn][i].node->isClone( edges_[drxn][j].node );
            for ( int k = 0; k < min( seqs[i].size(), seqs[j].size() ); k++ ) if ( seqs[i][k] != seqs[j][k] ) match = false;
            assert( !match );
        }
    }
}

void Node::nodeTest()
{
    for ( Edge &e : edges_[0] ) assert( e.node->isEdge( this, 1 ) );
    for ( Edge &e : edges_[1] ) assert( e.node->isEdge( this, 0 ) );
}

void Node::readTest()
{
    int32_t limits[2]{ ends_[1], ends_[0] };
    for ( auto &read : reads_ )
    {
        limits[0] = min( limits[0], read.second[0] );
        limits[1] = max( limits[1], read.second[1] );
    }
    for ( NodeMark &mark : pe_[0] ) assert( mark.coords[0] >= ends_[0] && mark.coords[1] <= ends_[1] );
    for ( NodeMark &mark : pe_[1] ) assert( mark.coords[0] >= ends_[0] && mark.coords[1] <= ends_[1] );
    for ( NodeMark &mark : mp_[0] ) assert( mark.coords[0] >= ends_[0] && mark.coords[1] <= ends_[1] );
    for ( NodeMark &mark : mp_[1] ) assert( mark.coords[0] >= ends_[0] && mark.coords[1] <= ends_[1] );
    assert( limits[0] == ends_[0] && limits[1] == ends_[1] );
}
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

void Node::pairTest()
{
    for ( auto np : pairs_ )
    {
        auto it = np.first->pairs_.find( this );
        assert( it != np.first->pairs_.end() && np.second == it->second );
    }
}

void Node::addCloned( Node* clone )
{
    assert( !clone->cloned_ );
    assert( size() == clone->size() );
    clone->cloned_ = new NodeRoll( this );
    if ( !cloned_ ) cloned_ = new NodeRoll;
    for ( Node* alt : cloned_->nodes )
    {
        *alt->cloned_ += clone;
        *clone->cloned_ += alt;
    }
    *cloned_ += clone;
    for ( Node* alt : cloned_->nodes ) assert( alt->cloned_->size() == cloned_->size() );
    for ( auto& read : reads_ ) assert( clone->reads_.find( read.first ) != clone->reads_.end() );
    assert( reads_.size() == clone->reads_.size() );
}

void Node::addEdge( Node* node, bool drxn )
{
    int overlap = abs( ends_[drxn] - node->ends_[!drxn] );
    edges_[drxn].push_back( Edge( node, overlap, false ) );
    node->edges_[!drxn].push_back( Edge( this, overlap, false ) );
    edgeCount_[drxn] = max( edgeCount_[drxn], (int)edges_[drxn].size() );
    node->edgeCount_[!drxn] = max( node->edgeCount_[!drxn], (int)node->edges_[!drxn].size() );
}

void Node::addEdge( Node* node, int overlap, bool drxn, bool doOffset, bool isLeap )
{
    assert( overlap < params.readLen );
    assert( node && node != this );
    assert( overlap <= ends_[1] - ends_[0] );
    assert( overlap <= node->ends_[1] - node->ends_[0] );
    
    addEdge( Edge( node, overlap, isLeap ), drxn, true );
    edgeTest( drxn );
    node->edgeTest( !drxn );
    
    if ( doOffset && node->drxn_ != 2 )
    {
        if ( drxn ? drxn_ == 0 : drxn_ == 1 ) propagateOffset( !drxn );
        else node->propagateOffset( drxn );
    }
    edgeCount_[drxn] = max( edgeCount_[drxn], (int)edges_[drxn].size() );
    node->edgeCount_[!drxn] = max( node->edgeCount_[!drxn], (int)node->edges_[!drxn].size() );
    setCoverage();
    node->setCoverage();
    if ( bad_ && !node->bad_ && ( node->drxn_ == 2 || drxn == !node->drxn_ ) ) setNotBad( !drxn );
    if ( !bad_ && node->bad_ && ( drxn_ == 2 || drxn == drxn_ ) ) node->setNotBad( drxn );
}

void Node::addEdge( Edge edge, bool drxn, bool reciprocate )
{
    assert( edge.node );
    stop( 0, drxn );
    for ( Edge& e : edges_[drxn] )
    {
        if ( e.node != edge.node ) continue;
        assert( e.ol == edge.ol );
        e.ol = edge.ol;
        return;
    }
    edges_[drxn].push_back( edge );
    if ( reciprocate ) edge.node->addEdge( Edge( this, edge.ol, edge.leap ), !drxn, false );
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
        if ( unblankLen < e.ol )
        {
            e.node->blankEnd( e.ol - unblankLen, drxn );
        }
    }
}

void Node::clearEdges( bool drxn )
{
    for ( Edge &e: edges_[drxn] ) e.node->removeEdge( this, !drxn );
    edges_[drxn].clear();
}

vector<Node*> Node::clones( bool inclSelf )
{
    vector<Node*> nodes;
    if ( inclSelf ) nodes.push_back( this );
    if ( cloned_ ) for ( Node* node : cloned_->nodes ) nodes.push_back( node );
    return nodes;
}

//int Node::countEdges( bool drxn, bool inclClones )
//{
//    if ( !inclClones || !cloned_ ) return edges_[drxn].size();
//    Nodesx edges;
//    for ( Edge& e : edges_[drxn] ) edges += e.node;
//    for ( Node* clone : cloned_->nodes ) for ( Edge& e : clone->edges_[drxn] ) edges += e.node;
//    return edges.size();
//}

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
                    auto it = t->reads_.find( mark.id );
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
            auto it = t->reads_.find( mark.id );
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

NodeRoll Node::dismantle()
{
    NodeRoll dependent = getDependent();
    for ( int d : { 0, 1 } )
    {
        for ( Edge &e : edges_[d] )
        {
            e.node->removeEdge( this, !d );
            e.node->setBad( d );
        }
        edges_[d].clear();
    }
    for ( const pair<Node*, int> &pn : pairs_ ) if ( pn.first != this ) pn.first->pairs_.erase( this );
    pairs_.clear();
    hits_.clean();
    if ( clones_ )
    {
        for ( Node* clone : *clones_ ) clone->removeClone( this );
        delete clones_;
        clones_ = NULL;
    }
    if ( cloned_ )
    {
        for ( Node* clone : cloned_->nodes )
        {
            *clone->cloned_ -= this;
            if ( !clone->cloned_->empty() ) continue;
            delete clone->cloned_;
            clone->cloned_ = NULL;
        }
        delete cloned_;
        cloned_ = NULL;
    }
    if ( paired_ )
    {
        delete paired_;
        paired_ = NULL;
    }
    if ( pairedNodes_ )
    {
        delete pairedNodes_;
        pairedNodes_ = NULL;
    }
    
    return dependent;
}

vector<Edge> Node::edges( bool drxn )
{
    return edges_[drxn];
}

//vector<Edge> Node::getAltEdges( bool drxn )
//{
//    vector<Edge> edges;
//    for ( Node* alt : getAltForks( drxn ) )
//    {
//        for ( Edge& e : alt->edges_[drxn] )
//        {
//            bool added = false;
//            for ( int i = 0; !added && i < edges.size(); i++ ) added = e.node == edges[i].node;
//            if ( !added ) edges.push_back( e );
//        }
//    }
//    return edges;
//}

//vector<Node*> Node::getAltForks( bool drxn )
//{
//    vector<Node*> nodes{ this };
//    if ( !cloned_ ) return nodes;
//    
//    if ( edges_[drxn].empty() )
//    {
//        for ( Node* clone : cloned_->nodes ) nodes.push_back( clone );
//    }
//    else if ( edges_[!drxn].empty() )
//    {
//        for ( Node* clone : cloned_->nodes ) if ( !clone->edges_[drxn].empty() ) nodes.push_back( clone );
//    }
//    else 
//    {
//        for ( Node* clone : cloned_->nodes ) if ( clone->edges_[!drxn].empty() ) nodes.push_back( clone );
//    }
//    
//    return nodes;
//}

int32_t Node::getBiggestOffset( bool drxn )
{
    int32_t biggest = 0;
    for ( Edge &e : edges_[drxn] )
    {
        biggest = max( biggest, abs( ( drxn ? ends_[1] - e.node->ends_[0] : e.node->ends_[1] - ends_[0] ) - e.ol ) );
    }
    return biggest;
}

int Node::getBestOverlap( bool drxn, bool inclClones )
{
    int ol = edges_[drxn].empty() ? 0 : edges_[drxn][0].ol;
    for ( int i = 1; i < edges_[drxn].size(); i++ ) ol = max( ol, edges_[drxn][i].ol );
    if ( inclClones && cloned_ ) for ( Node* clone : cloned_->nodes ) ol = max( ol, clone->getBestOverlap( drxn, false ) );
    return ol;
}

NodeRoll Node::getDependent()
{
    NodeRoll nodes;
    for ( int d : { 0, 1 } )
    {
        if ( drxn_ != 2 && d != drxn_ ) continue;
        for ( Edge &e : edges_[d] ) if ( e.node->drxn_ == d && e.node->edges_[!d].size() == 1 ) nodes.add( e.node );
    }
    return nodes;
}

Edge Node::getEdge( Node* node, bool drxn )
{
    for ( Edge e : edges_[drxn] ) if ( e.node == node ) return e;
    return Edge( NULL, 0, false );
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
            header += edge.node->id_ + "x" + to_string( edge.ol );
            first = false;
        }
    }
    return header;
}

int32_t Node::getLen( vector<Node*>& path )
{
    int32_t len = path[0]->size();
    for ( int i = 1; i < path.size(); i++ )
    {
        Edge e = path[i]->getEdge( path[i-1], 0 );
        assert( e.node && e.ol > 0 );
        len += path[i]->size() - e.ol;
    }
    return len;
}

int32_t Node::getLength()
{
    int32_t len = ends_[1] - ends_[0];
    int seqLen = seq_.length();
    assert( len == seq_.length() );
    return len;
}

int Node::getOverlap( Node* node, bool drxn )
{
    int overlap = 0;
    for ( Edge &edge : edges_[drxn] )
    {
        if ( edge.node == node )
        {
            overlap = edge.ol;
            break;
        }
    }
    return overlap;
}

string Node::getSeq( vector<Node*>& path )
{
    string seq = path[0]->seq_;
    for ( int i = 1; i < path.size(); i++ )
    {
        Edge e = path[i]->getEdge( path[i-1], 0 );
        assert( e.node && e.ol > 0 );
        seq += path[i]->seq_.substr( e.ol );
    }
    return seq;
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

bool Node::inheritEdges( bool drxn )
{
    if ( drxn_ >= 2 ) return false;
    NodeSet nxtSet = getNextNodes( drxn );
    NodeSet prvSet;
    for ( Node* nxt : nxtSet ) nxt->getNextNodes( prvSet, !drxn );
    if ( prvSet.size() <= 1 ) return false;
    bool didEdge = false;
    for ( Node* prv : prvSet )
    {
        if ( prv == this ) continue;
        NodeSet prvNxtSet = prv->getNextNodes( drxn );
        for ( Node* nxt : nxtSet )
        {
            if ( prvNxtSet.find( nxt ) == nxtSet.end() )
            {
                int ol = ( drxn ? mapSeqOverlap( prv->seq_, nxt->seq_, 23 ) 
                                : mapSeqOverlap( nxt->seq_, prv->seq_, 23 ) );
                if ( ol && ol < params.readLen )
                {
                    Node* bck = ( drxn_ == drxn ? prv : nxt );
                    Node* fwd = ( drxn_ == drxn ? nxt : prv );
                    if ( bck == fwd ) continue;
                    NodeSet bckFwdSet = bck->getDrxnNodes( drxn_ );
                    NodeSet fwdFwdSet = fwd->getDrxnNodes( drxn_ );
                    if ( bckFwdSet.find( fwd ) != bckFwdSet.end() ) continue;
                    if ( fwdFwdSet.find( bck ) != fwdFwdSet.end() ) continue;
                    bck->addEdge( fwd, ol, drxn_ );
                    NodeSet qSet = bck->getDrxnNodesNotInSet( bckFwdSet, drxn_ );
                    for ( Node* q : qSet ) if ( q->validated_ ) q->updatePairs();
                    didEdge = true;
                }
            }
        }
    }
    
    return didEdge;
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

bool Node::isClone( Node* node )
{
    return node == this || ( cloned_ && cloned_->find( node ) );
}

bool Node::isClone( vector<Node*>& path, int i, bool drxn )
{
    if ( !isClone( drxn ? path[i] : path.end()[-i-1] ) ) return false;
    if ( ++i < path.size() ) for ( Edge& e : edges_[drxn] ) if ( e.node->isClone( path, i, drxn ) ) return true;
    return i >= path.size();
}

bool Node::isContinue( bool drxn )
{
    return edges_[drxn].empty() && stop_[drxn] == 0 && !cloned_;
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

bool Node::isEdge( Node* node, bool drxn, bool inclClone )
{
    for ( Edge &e : edges_[drxn] )
    {
        if ( e.node == node ) return true;
        if ( e.node->cloned_ && inclClone ) for ( Node* clone : e.node->cloned_->nodes ) if ( clone == node ) return true;
    }
    return false;
}

bool Node::isEnded( bool drxn )
{
    return ( edges_[drxn].empty() && stop_[drxn] );
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
            if ( propagated.find( edge.node ) != propagated.end() || edge.node->drxn_ == 2 ) continue;
            edge.node->propagateOffset( propagated, drxn );
        }
    }
}

bool Node::removeEdge( Node* node, bool drxn, bool reciprocate )
{
    bool removed = false;
    for ( auto it = edges_[drxn].begin(); it != edges_[drxn].end(); )
    {
        if ( it->node == node )
        {
            it = edges_[drxn].erase( it );
            removed = true;
        }
        else it++;
    }
    if ( edges_[drxn].empty() ) stop( BLUNT_END, drxn );
    if ( removed && reciprocate ) node->removeEdge( this, !drxn );
    return removed;
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

void Node::setOrigin()
{
    ends_.origin[0] = ends_[0];
    ends_.origin[1] = ends_[1];
}

int Node::size()
{
    return ends_[1] - ends_[0];
}

void Node::sortEdges( bool drxn )
{
    sort( edges_[drxn].begin(), edges_[drxn].end(), []( const Edge &a, const Edge &b ){
        return a.ol > b.ol;
    }  );
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

void Node::trimEnd( int32_t coord, NodeList &nodes, bool drxn )
{
    NodeSet delSet;
    for ( Node* nxt : getNextNodes( drxn ) ) nxt->dismantleNode( delSet, drxn );
    for ( int i = 0; i < nodes.size(); )
    {
        if ( delSet.find( nodes[i] ) != delSet.end() ) nodes.erase( nodes.begin() + i );
        else i++;
    }
    for ( Node* del : delSet ) delete del;
    trimSeq( coord, drxn );
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

bool Node::unpause( bool drxn )
{
    if ( unpauseable( drxn ) ) stop_[drxn] = 0;
    return isContinue( drxn );
}

bool Node::unpauseable( bool drxn )
{
    return stop_[drxn] == PAUSED || stop_[drxn] == BACK_END;
}

void Node::deleteNodes( NodeSet &delSet, NodeList &nodes )
{
    for ( int i = 0; i < nodes.size(); i++ )
    {
        if ( delSet.find( nodes[i] ) == delSet.end() ) continue;
        nodes[i]->dismantleNode();
        delete nodes[i];
        nodes.erase( nodes.begin() + i-- );
    }
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
    assert( drxn_ != 2 );
    for ( Edge edge : edges_[!drxn] )
    {
        edge.node->removeEdge( this, drxn );
        if ( edge.node->edges_[drxn].empty() )
        {
            edge.node->stop_[drxn] = 3;
        }
    }
    for ( Edge edge : edges_[drxn] )
    {
        assert( edge.node->drxn_ != 2 );
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
