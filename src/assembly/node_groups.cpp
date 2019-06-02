/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "node_groups.h"
#include "node.h"
#include <algorithm>
#include <iostream>
#include <fstream>

Nodes::Nodes( Node* node )
{
    nodes.insert( node );
}

Nodes::Nodes( vector<Node*> base )
: nodes( base.begin(), base.end() )
{ }

Nodes::Nodes( Node* node, bool drxn, bool inclNode, bool inclClone )
{
    fill( node, drxn, inclNode, inclClone );
}

Nodes::Nodes( Node* node, int32_t limit, bool drxn, bool inclNode )
{
    assert( false );
}

Nodes::Nodes( Node* node, int32_t limit, bool drxn, bool inclNode, bool inclClone )
{
    fill( node, limit, drxn, inclNode, inclClone );
}

Nodes::Nodes( vector<Node*> &base, bool drxn, bool inclBase )
{
    if ( inclBase ) nodes.insert( base.begin(), base.end() );
    for ( Node* node : base ) node->getDrxnNodes( nodes, drxn );
}

void Nodes::operator += ( Node* node )
{
    nodes.insert( node );
}

void Nodes::operator += ( Nodes &rhs )
{
    for ( Node* n : rhs.nodes ) nodes.insert( n );
}

void Nodes::operator -= ( Node* node )
{
    nodes.erase( node );
}

void Nodes::operator -= ( Nodes &rhs )
{
    for ( Node* n : rhs.nodes ) nodes.erase( n );
}

bool Nodes::add( Node* node )
{
    if ( !node ) return false;
    return nodes.insert( node ).second;
}

Nodes Nodes::between( Node* a, Node* b, bool inclNodes )
{
    Nodes fwd( a, 1, inclNodes, false );
    return Nodes::inSet( b, fwd, 0, inclNodes );
}

void Nodes::cancel( Nodes &a, Nodes &b )
{
    for ( auto it = a.nodes.begin(); it != a.nodes.end(); )
    {
        if ( b.erase( *it ) ) it = a.nodes.erase( it );
        else it++;
    }
}

void Nodes::clear()
{
    nodes.clear();
}

Nodes Nodes::connected( Node* node )
{
    Nodes nodes;
    nodes.fill( node );
    return nodes;
}

void Nodes::dumpBad( bool bad )
{
    for ( auto it = nodes.begin(); it != nodes.end(); )
    {
        if ( (*it)->bad_ == bad ) it = nodes.erase( it );
        else it++;
    }
}

void Nodes::dumpGraph( int graph, bool dump )
{
    for ( auto it = nodes.begin(); it != nodes.end(); )
    {
        if ( (*it)->drxn_ == graph ? dump : !dump ) it = nodes.erase( it );
        else it++;
    }
}

bool Nodes::empty()
{
    return nodes.empty();
}

bool Nodes::erase( Node* node )
{
    auto it = nodes.find( node );
    if ( it == nodes.end() ) return false;
    nodes.erase( it );
    return true;
}

Nodes Nodes::exclusive( Nodes& ignore )
{
    Nodes ex;
    for ( Node* node : nodes ) if ( !ignore.find( node ) ) ex.add( node );
    return ex;
}

void Nodes::fill( Node* node )
{
    nodes.insert( node );
    for ( Edge &e : node->edges_[0] ) if ( !find( e.node ) ) fill( e.node );
    for ( Edge &e : node->edges_[1] ) if ( !find( e.node ) ) fill( e.node );
}

void Nodes::fill( Node* node, bool drxn, bool inclNode, bool inclClone )
{
    if ( inclNode && !add( node ) ) return;
    for ( Edge &e : node->edges_[drxn] ) if ( !find( e.node ) ) fill( e.node, drxn, true, inclClone );
    if ( !inclClone || !node->cloned_ || node->edges_[!drxn].empty() ) return;
    for ( Node* clone : node->cloned_->nodes )
    {
        if ( clone->edges_[drxn].empty() ) continue;
        if ( node->edges_[drxn].empty() || clone->edges_[!drxn].empty() ) fill( clone, drxn, false, true );
    }
}

void Nodes::fill( Node* node, int32_t limit, bool drxn, bool inclNode, bool inclClone )
{
    if ( inclNode && !add( node ) ) return;
    if ( limit <= 0 ) return;
    for ( Edge &e : node->edges_[drxn] ) fill( e.node, limit - e.node->size() + e.ol, drxn, true, inclClone );
    if ( !inclClone || !node->cloned_ || node->edges_[!drxn].empty() ) return;
    for ( Node* clone : node->cloned_->nodes )
    {
        if ( clone->edges_[drxn].empty() ) continue;
        if ( node->edges_[drxn].empty() || clone->edges_[!drxn].empty() ) fill( clone, limit, drxn, false, true );
    }
}

void Nodes::fillBad( Node* node, bool bad, bool drxn )
{
    if ( node->bad_ == bad ) nodes.insert( node );
    for ( Edge &e : node->edges_[drxn] ) if ( e.node->bad_ == bad && !find( e.node ) ) fillBad( e.node, bad, drxn );
}

void Nodes::fillBranch( Node* node )
{
    if ( node->branch_ ) nodes.insert( node );
    for ( Edge &e : node->edges_[0] ) if ( e.node->branch_ && !find( e.node ) ) fillBranch( e.node );
    for ( Edge &e : node->edges_[1] ) if ( e.node->branch_ && !find( e.node ) ) fillBranch( e.node );
}

void Nodes::fillIn( Node* node, Nodes &include, bool drxn, bool inclNode )
{
    if ( inclNode && ( !include.find( node ) || !add( node ) ) ) return;
    for ( Edge &e : node->edges_[drxn] ) if ( !find( e.node ) && include.find( e.node ) ) fillIn( e.node, include, drxn, true );
//    if ( node->cloned_ && inclNode ) for ( Node* clone : node->cloned_->nodes ) fillIn( clone, include, drxn, false );
}

void Nodes::fillIn( Node* node, Nodes &include, int32_t limit, bool drxn, bool inclNode )
{
    if ( inclNode ) nodes.insert( node );
    if ( limit <= 0 ) return;
    for ( Edge &e : node->edges_[drxn] ) if ( !find( e.node ) && include.find( e.node ) ) fillIn( e.node, include, limit - e.node->size() + e.ol, drxn, true );
//    if ( node->cloned_ && inclNode ) for ( Node* clone : node->cloned_->nodes ) fillIn( clone, include, limit, drxn, false );
}

void Nodes::fillNext( Node* node, bool drxn )
{
    for ( Edge &e : node->edges_[drxn] ) nodes.insert( e.node );
}

void Nodes::fillNot( Node* node, Nodes &ignore, bool drxn, bool inclNode )
{
    if ( inclNode && ( ignore.find( node ) || !add( node ) ) ) return;
    for ( Edge &e : node->edges_[drxn] ) fillNot( e.node, ignore, drxn, true );
}

void Nodes::fillNot( Node* node, Nodes &ignore, int32_t limit, bool drxn, bool inclNode )
{
    if ( inclNode ) nodes.insert( node );
    if ( limit <= 0 ) return;
    for ( Edge &e : node->edges_[drxn] ) if ( !find( e.node ) && !ignore.find( e.node ) ) fillNot( e.node, ignore, limit - e.node->size() + e.ol, drxn, true );
//    if ( node->cloned_ && inclNode ) for ( Node* clone : node->cloned_->nodes ) fillNot( clone, ignore, limit, drxn, false );
}

bool Nodes::find( Node* node )
{
    return nodes.find( node ) != nodes.end();
}

void Nodes::forkSets( Node* fork, Node* branch, Nodes fwds[2], Nodes bcks[2], int32_t limit, bool drxn )
{
    Nodes fwd( fork, limit, drxn, false ), bck( branch, limit, !drxn, false );
    for ( Edge& e : fork->edges_[drxn] ) fwds[ branch != e.node ].fillIn( e.node, fwd, drxn, true );
    for ( Edge& e : branch->edges_[!drxn] ) bcks[ fork != e.node ].fillIn( e.node, bck, !drxn, true );
    Nodes::cancel( fwds[0], fwds[1] );
    Nodes::cancel( bcks[0], bcks[1] );
}

Nodes Nodes::inSet( Node* node, Nodes &include, bool drxn, bool inclNode )
{
    Nodes nodes;
    nodes.fillIn( node, include, drxn, inclNode );
    return nodes;
}

Nodes Nodes::inSet( Node* node, Nodes &include, int32_t limit, bool drxn, bool inclNode )
{
    Nodes nodes;
    nodes.fillIn( node, include, limit, drxn, inclNode );
    return nodes;
}

Nodes Nodes::isBad( Node* node, bool bad, bool drxn )
{
    Nodes nodes;
    nodes.fillBad( node, bad, drxn );
    return nodes;
}

Nodes Nodes::isBranch( Node* node )
{
    Nodes nodes;
    nodes.fillBranch( node );
    return nodes;
}

Nodes Nodes::notSet( Node* node, Nodes &ignore, bool drxn, bool inclNode )
{
    Nodes nodes;
    nodes.fillNot( node, ignore, drxn, inclNode );
    return nodes;
}

Nodes Nodes::notSet( Node* node, Nodes &ignore, int32_t limit, bool drxn, bool inclNode )
{
    Nodes nodes;
    nodes.fillNot( node, ignore, limit, drxn, inclNode );
    return nodes;
}

int Nodes::size()
{
    return nodes.size();
}

NodeRoll::NodeRoll( Node* node )
: nodes( 1, node )
{ }

NodeRoll::NodeRoll( Nodes base )
: nodes( base.nodes.begin(), base.nodes.end() )
{
    
}

Node* NodeRoll::operator[]( int i )
{
    assert( i < nodes.size() );
    return nodes[i];
}

void NodeRoll::operator += ( NodeRoll &rhs )
{
    for ( Node* node : rhs.nodes ) nodes.push_back( node );
}

void NodeRoll::operator += ( Nodes rhs )
{
    nodes.insert( nodes.end(), rhs.nodes.begin(), rhs.nodes.end() );
}

void NodeRoll::operator -= ( Nodes& rhs )
{
    for ( int i = 0; i < nodes.size(); i++ ) if ( rhs.find( nodes[i] ) ) nodes.erase( nodes.begin() + i-- );
}

void NodeRoll::operator += ( Node* node )
{
    nodes.push_back( node );
}

void NodeRoll::operator -= ( Node* node )
{
    nodes.erase( std::remove( nodes.begin(), nodes.end(), node ), nodes.end() );
}

bool NodeRoll::add( Node* node )
{
    if ( find( node ) ) return false;
    nodes.push_back( node );
    return true;
}

void NodeRoll::add( NodeRoll &rhs )
{
    for ( Node* node : rhs.nodes ) add( node );
}

void NodeRoll::clear()
{
    nodes.clear();
}

NodeRoll NodeRoll::clones( Node* node )
{
    NodeRoll nodes( node );
    if ( node->cloned_ ) for ( Node* clone : node->cloned_->nodes ) nodes += clone;
    return nodes;
}

bool NodeRoll::find( Node* node )
{
    return std::find( nodes.begin(), nodes.end(), node ) != nodes.end();
}

void NodeRoll::erase( Node* node )
{
    assert( find( node ) );
    node->dismantle();
    delete node;
    remove( node );
}

void NodeRoll::erase( Node* node, int &i )
{
    for ( int j = 0; j < size(); j++ )
    {
        if ( nodes[j] != node ) continue;
        NodeRoll dependent = nodes[j]->dismantle();
        nodes.erase( nodes.begin() + j );
        delete node;
        if ( j <= i ) i--;
        for ( Node* dep : dependent.nodes ) erase( dep, i );
        return;
    }
}

bool NodeRoll::empty()
{
    return nodes.empty();
}

NodeRoll NodeRoll::getGraph( int drxn )
{
    NodeRoll graph;
    for ( Node* node : nodes ) if ( node->drxn_ == drxn ) graph.add( node );
    return graph;
}

NodeRoll NodeRoll::next( Node* node, bool drxn )
{
    NodeRoll nodes;
    for ( Edge& e : node->edges_[drxn] ) nodes += e.node;
    return nodes;
}

void NodeRoll::remove( Node* node )
{
    nodes.erase( std::remove( nodes.begin(), nodes.end(), node ), nodes.end() );
}

void NodeRoll::print( string fn, int32_t i, int32_t j )
{
    Nodes used[2];
    for ( Node* node : nodes ) node->id_ = to_string( 0 );
    for ( Node* node : nodes ) if ( ( i < node->ends_[0] && node->ends_[0] < j ) || ( i < node->ends_[1] && node->ends_[1] < j ) ) used[0] += node;
    vector<Node*> tar( used[0].nodes.begin(), used[0].nodes.end() );
    sort( tar.begin(), tar.end(), []( Node* a, Node* b ){ return a->ends_[0] < b->ends_[0]; } );
    assert( !tar.empty() );
    int32_t limit = tar[0]->ends_[0];
    for ( int k = 0; k < tar.size(); k++ ) tar[k]->id_ = to_string( k+1 );
    
    ofstream ofs( fn );
    for ( Node* base : tar ) if ( !base->bad_ && !used[1].find( base ) ) print( ofs, base, used, limit );
    ofs.close();
}

void NodeRoll::print( ofstream& ofs, Node* node, Nodes used[2], int32_t limit )
{
    if ( !used[1].add( node ) ) return;
    string header = ">NODE " + string( 5 - node->id_.size(), '0' ) + node->id_ + " |";
    if ( node->edges_[0].empty() ) header += " END";
    for ( Edge& e : node->edges_[0] )
    {
        if ( used[0].find( e.node ) && !used[1].find( e.node ) ) print( ofs, e.node, used, limit );
        header += " " + e.node->id_;
    }
    header += " |";
    if ( node->edges_[1].empty() ) header += " END";
    for ( Edge& e : node->edges_[1] ) header += " " + e.node->id_;
    header += " | READS " + to_string( node->countReads( true ) );
    if ( node->cloned_ )
    {
        header += " | CLONED";
        for ( Node* clone : node->cloned_->nodes ) header += " " + clone->id_;
    }
    
    ofs << header << "\n" << string( node->ends_[0] - limit, '-' ) << node->seq_ << "\n";
    for ( Edge& e : node->edges_[1] )
    {
        if ( used[0].find( e.node ) && !used[1].find( e.node ) ) print( ofs, e.node, used, limit );
    }
}

int NodeRoll::size()
{
    return nodes.size();
}

void NodeRoll::test( bool loop )
{
    if ( !loop ) return;
    Nodes drxn[2];
    for ( Node* n : nodes ) if ( n->drxn_ == 2 ) for ( int d : { 0, 1 } ) drxn[d].fill( n, d, true, false );
    for ( int d : { 0, 1 } ) for ( Node* n : drxn[d].nodes ) if ( n->drxn_ < 2 ) assert( n->drxn_ == d );
    for ( Node* n : nodes )
    {
        assert( n->pairs_.empty() );
        n->hits_.test();
        n->pairTest();
        bool bad = n->drxn_ < 2;
        for ( Edge &e : n->edges_[0] )
        {
            string lSeq = e.node->seq_.substr( e.node->seq_.size() - e.ol );
            string rSeq = n->seq_.substr( 0, e.ol );
            assert( lSeq == rSeq || e.isLeap );
            assert( e.node->isEdge( n, 1 ) );
            if ( n->drxn_ == 1 && ( e.node->drxn_ == 1 || e.node->drxn_ == 2 ) && !e.node->bad_ ) bad = false;
            assert( find( e.node ) );
        }
        for ( Edge &e : n->edges_[1] )
        {
            string lSeq = n->seq_.substr( n->seq_.size() - e.ol );
            string rSeq = e.node->seq_.substr( 0, e.ol );
            assert( lSeq == rSeq || e.isLeap );
            assert( e.node->isEdge( n, 0 ) );
            if ( n->drxn_ == 0 && ( e.node->drxn_ == 0 || e.node->drxn_ == 2 ) && !e.node->bad_ ) bad = false;
            assert( find( e.node ) );
        }
        if ( n->drxn_ < 2 && !n->bad_ ) assert( drxn[n->drxn_].find( n ) );
        assert( !bad || n->bad_ );
        assert( bad == n->bad_ );
        if ( loop && !n->bad_ ) assert( !Nodes( n, true, false, false ).find( n ) );
    }
}

bool NodeOffset::add( int32_t dist, bool minOnly )
{
    if ( dist == dists[0] || dist == dists[1] || dist == dists[2] ) return false;
    else if ( abs( dist ) < abs( dists[0] ) )
    {
        dists[0] = dist;
        return true;
    }
    else if ( minOnly ) return false;
    else if ( abs( dists[2] ) < abs( dist ) )
    {
        dists[2] = dist;
        return true;
    }
    else if ( min( abs( dists[0] - dists[1] ), abs( dists[1] - dists[2] ) ) < min( abs( dists[0] - dist ), abs( dist - dists[2] ) ) )
    {
        dists[1] = dist;
        return true;
    }
    return false;
}

int32_t NodeOffset::diff( int32_t est )
{
    int32_t best = abs( est - abs( dists[0] ) );
    for ( int i = 1; i < 3; i++ ) if ( dists[i-1] != dists[i] ) best = min( best, abs( est - abs( dists[i] ) ) );
    return best;
}

int32_t NodeOffset::diff( NodeOffset& rhs, int32_t est, bool drxn )
{
    int32_t best = abs( est - ( drxn ? rhs[0] - dists[0] : dists[0] - rhs[0] ) );
    for ( int i = 1; i < 3; i++ )
    {
        if ( dists[i-1] == dists[i] ) continue;
        for ( int j = 1; j < 3; j++ )
        {
            if ( rhs[j-1] == rhs[j] ) continue;
            best = min( best, abs( est - ( drxn ? dists[i] - rhs[j] : rhs[j] - dists[i] ) ) );
        }
    }
    return best;
}

NodeOffsets::NodeOffsets( Node* node, bool drxn, bool inclNode, bool inclClone )
{
    fill( node, 0, 0, drxn, drxn, inclNode, inclClone );
}

NodeOffsets::NodeOffsets( Node* node, int32_t limit, bool orient, bool drxn, bool inclNode )
{
    fill( node, 0, limit, orient, drxn, inclNode, true );
}

bool NodeOffsets::add( Node* node, int32_t dist )
{
    auto it = map.find( node );
    if ( it != map.end() ) return it->second.add( dist );
    map.insert( make_pair( node, NodeOffset( dist ) ) );
    return true;
}

bool NodeOffsets::add( Node* fork, Edge& e, int32_t base, bool orient, bool drxn )
{
    int32_t dist = ( ( orient == drxn ? e.node : fork )->size() - e.ol ) * ( drxn ? 1 : -1 );
    return add( e.node, dist + base );
}

void NodeOffsets::clear()
{
    map.clear();
}

bool NodeOffsets::doLoop( Node* clone, int32_t dist, bool drxn )
{
    NodeOffset* off = get( clone );
    return !off || ( drxn ? (*off)[0] < dist : dist < (*off)[0] );
}

bool NodeOffsets::empty()
{
    return map.empty();
}

void NodeOffsets::erase( Node* node )
{
    map.erase( node );
}

int32_t NodeOffsets::extend( Node* fork, Edge& e, int32_t dist, bool orient, bool drxn )
{
    return dist + ( ( drxn ? 1 : -1 ) * ( ( orient == drxn ? e.node : fork )->size() - e.ol ) );
}

bool NodeOffsets::find( Node* node )
{
    return map.find( node ) != map.end();
}

void NodeOffsets::fill( Node* node, int32_t dist, int32_t limit, bool orient, bool drxn, bool inclNode, bool inclClone, bool looped )
{
    if ( inclNode && looped && ( looped = ( !doLoop( node, dist, drxn ) ) ) ) return;
    if ( inclNode && !add( node, dist ) ) return;
    if ( limit && ( drxn ? abs( limit ) < dist : dist < -abs( limit ) ) ) return;
    for ( Edge &e : node->edges_[drxn] ) fill( e.node, extend( node, e, dist, orient, drxn ), limit, orient, drxn, true, inclClone, looped );
    if ( !inclClone || looped || !node->cloned_ || node->edges_[!drxn].empty() ) return;
    for ( Node* clone : node->cloned_->nodes )
    {
        if ( node->edges_[drxn].empty() || clone->edges_[!drxn].empty() ) fill( clone, dist, limit, orient, drxn, inclNode, true, true );
    }
}

void NodeOffsets::fillIn( Node* node, Nodes& include, int32_t dist, bool orient, bool drxn, bool inclNode, bool looped )
{
    if ( inclNode && looped && ( looped = ( !doLoop( node, dist, drxn ) ) ) ) return;
    if ( inclNode && ( !include.find( node ) || !add( node, dist ) ) ) return;
    for ( Edge &e : node->edges_[drxn] ) fillIn( e.node, include, extend( node, e, dist, orient, drxn ), orient, drxn, true, looped );
    if ( !node->cloned_ || node->edges_[!drxn].empty() ) return;
    for ( Node* clone : node->cloned_->nodes )
    {
        if ( node->edges_[drxn].empty() || clone->edges_[!drxn].empty() ) fillIn( clone, include, dist, orient, drxn, inclNode, true );
    }
}

void NodeOffsets::fillNot( Node* node, Nodes& ignore, int32_t dist, int32_t limit, bool orient, bool drxn, bool inclNode, bool looped )
{
    if ( inclNode && looped && ( looped = ( !doLoop( node, dist, drxn ) ) ) ) return;
    if ( inclNode && ( ignore.find( node ) || !add( node, dist ) ) ) return;
    if ( limit && ( drxn ? abs( limit ) < dist : dist < -abs( limit ) ) ) return;
    for ( Edge &e : node->edges_[drxn] ) fillNot( e.node, ignore, extend( node, e, dist, orient, drxn ), limit, orient, drxn, true, looped );
    if ( looped || !node->cloned_ || node->edges_[!drxn].empty() ) return;
    for ( Node* clone : node->cloned_->nodes )
    {
        if ( node->edges_[drxn].empty() || clone->edges_[!drxn].empty() ) fillNot( clone, ignore, dist, orient, drxn, inclNode, true );
    }
}

//void NodeOffsets::fillNot( Node* node, Nodes& ignore, int32_t dist, bool orient, bool drxn, bool inclNode, bool looped )
//{
//    if ( inclNode && looped && ( looped = ( !doLoop( node, dist, drxn ) ) ) ) return;
//    if ( inclNode && ( ignore.find( node ) || !add( node, dist ) ) ) return;
//    for ( Edge &e : node->edges_[drxn] ) fillNot( e.node, ignore, extend( node, e, dist, orient, drxn ), orient, drxn, true, looped );
//    if ( !node->cloned_ || node->edges_[!drxn].empty() ) return;
//    for ( Node* clone : node->cloned_->nodes )
//    {
//        if ( node->edges_[drxn].empty() || clone->edges_[!drxn].empty() ) fillNot( clone, ignore, dist, orient, drxn, inclNode, true );
//    }
//}

NodeOffset* NodeOffsets::get( Node* node )
{
    auto it = map.find( node );
    return it != map.end() ? &it->second : NULL;
}
