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

Nodes::Nodes( Node* node, bool drxn, bool inclNode )
{
    fill( node, drxn, inclNode);
}

Nodes::Nodes( Node* node, int32_t limit, bool drxn, bool inclNode )
{
    fill( node, limit, drxn, inclNode );
}

//Nodes::Nodes( vector<Node*> &base, bool drxn, bool inclBase )
//{
//    if ( inclBase ) nodes.insert( base.begin(), base.end() );
//    for ( Node* node : base ) node->getDrxnNodes( nodes, drxn );
//}

void Nodes::operator += ( Node* node )
{
    nodes.insert( node );
}

void Nodes::operator += ( Nodes &rhs )
{
    for ( Node* n : rhs.nodes ) nodes.insert( n );
}

void Nodes::operator += ( vector<Node*> &adds )
{
    nodes.insert( adds.begin(), adds.end() );
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

//Nodesx Nodesx::between( Node* a, Node* b, bool inclNodes )
//{
//    Nodesx fwd( a, 1, inclNodes, false );
//    return Nodesx::inSet( b, fwd, 0, inclNodes );
//}

//void Nodesx::cancel( Nodesx &a, Nodesx &b )
//{
//    for ( auto it = a.nodes.begin(); it != a.nodes.end(); )
//    {
//        if ( b.erase( *it ) ) it = a.nodes.erase( it );
//        else it++;
//    }
//}

void Nodes::clear()
{
    nodes.clear();
}

//Nodesx Nodesx::connected( Node* node )
//{
//    Nodesx nodes;
//    nodes.fill( node );
//    return nodes;
//}

//void Nodesx::dumpBad( bool bad )
//{
//    for ( auto it = nodes.begin(); it != nodes.end(); )
//    {
//        if ( (*it)->bad_ == bad ) it = nodes.erase( it );
//        else it++;
//    }
//}

//void Nodesx::dumpGraph( int graph, bool dump )
//{
//    for ( auto it = nodes.begin(); it != nodes.end(); )
//    {
//        if ( (*it)->drxn_ == graph ? dump : !dump ) it = nodes.erase( it );
//        else it++;
//    }
//}

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

//Nodesx Nodesx::exclusive( Nodesx& ignore )
//{
//    Nodesx ex;
//    for ( Node* node : nodes ) if ( !ignore.find( node ) ) ex.add( node );
//    return ex;
//}

void Nodes::fill( Node* node )
{
    nodes.insert( node );
    for ( Edge &e : node->edges_[0] ) if ( !find( e.node ) ) fill( e.node );
    for ( Edge &e : node->edges_[1] ) if ( !find( e.node ) ) fill( e.node );
}

void Nodes::fill( Node* node, bool drxn, bool inclNode )
{
    if ( inclNode && !add( node ) ) return;
    for ( Edge &e : node->edges_[drxn] ) if ( !find( e.node ) ) fill( e.node, drxn, true );
//    if ( !inclClone || !node->cloned_ || node->edges_[!drxn].empty() ) return;
//    for ( Node* clone : node->cloned_->nodes )
//    {
//        if ( clone->edges_[drxn].empty() ) continue;
//        if ( node->edges_[drxn].empty() || clone->edges_[!drxn].empty() ) fill( clone, drxn, false, true );
//    }
}

void Nodes::fill( Node* node, int32_t limit, bool drxn, bool inclNode )
{
    if ( inclNode && !add( node ) ) return;
    if ( limit <= 0 ) return;
    for ( Edge &e : node->edges_[drxn] ) fill( e.node, limit - e.node->size() + e.ol, drxn, true );
}

void Nodes::fillAll( Node* node )
{
    if ( add( node ) ) for ( int d : { 0, 1 } ) for ( Edge& e : node->edges_[d] ) fillAll( e.node );
}

//void Nodes::fillBad( Node* node, bool bad, bool drxn )
//{
//    if ( node->bad_ == bad ) nodes.insert( node );
//    for ( Edge &e : node->edges_[drxn] ) if ( e.node->bad_ == bad && !find( e.node ) ) fillBad( e.node, bad, drxn );
//}

//void Nodes::fillBranch( Node* node )
//{
//    if ( node->branch_ ) nodes.insert( node );
//    for ( Edge &e : node->edges_[0] ) if ( e.node->branch_ && !find( e.node ) ) fillBranch( e.node );
//    for ( Edge &e : node->edges_[1] ) if ( e.node->branch_ && !find( e.node ) ) fillBranch( e.node );
//}

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
}

//void Nodes::fillNext( Node* node, bool drxn )
//{
//    for ( Edge &e : node->edges_[drxn] ) nodes.insert( e.node );
//}

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
}

bool Nodes::find( Node* node )
{
    return nodes.find( node ) != nodes.end();
}

Nodes Nodes::looped( Node* node )
{
    Nodes fwd( node, 1, false ), loop;
    if ( fwd.find( node ) ) loop.fillIn( node, fwd, 0, true );
    return loop;
}

vector<Node*> Nodes::meet( Node* node )
{
    vector<Node*> met;
    Nodes pathed;
    if ( find( node ) ) met.push_back( node );
    else meet( node, pathed, met );
    return met;
}

void Nodes::meet( Node* node, Nodes& pathed, vector<Node*>& met )
{
    for ( int d: { 0, 1 } ) for ( Edge& e : node->edges_[d] ) if ( pathed.add( e.node ) )
    {
        if ( find( e.node ) ) met.push_back( e.node );
        else meet( e.node, pathed, met );
    }
}

//void Nodesx::forkSets( Node* fork, Node* branch, Nodesx fwds[2], Nodesx bcks[2], int32_t limit, bool drxn )
//{
//    Nodesx fwd( fork, limit, drxn, false ), bck( branch, limit, !drxn, false );
//    for ( Edge& e : fork->edges_[drxn] ) fwds[ branch != e.node ].fillIn( e.node, fwd, drxn, true );
//    for ( Edge& e : branch->edges_[!drxn] ) bcks[ fork != e.node ].fillIn( e.node, bck, !drxn, true );
//    Nodesx::cancel( fwds[0], fwds[1] );
//    Nodesx::cancel( bcks[0], bcks[1] );
//}

//Nodesx Nodesx::inSet( Node* node, Nodesx &include, bool drxn, bool inclNode )
//{
//    Nodesx nodes;
//    nodes.fillIn( node, include, drxn, inclNode );
//    return nodes;
//}

//Nodesx Nodesx::inSet( Node* node, Nodesx &include, int32_t limit, bool drxn, bool inclNode )
//{
//    Nodesx nodes;
//    nodes.fillIn( node, include, limit, drxn, inclNode );
//    return nodes;
//}

//Nodesx Nodesx::isBad( Node* node, bool bad, bool drxn )
//{
//    Nodesx nodes;
//    nodes.fillBad( node, bad, drxn );
//    return nodes;
//}

//Nodesx Nodesx::isBranch( Node* node )
//{
//    Nodesx nodes;
//    nodes.fillBranch( node );
//    return nodes;
//}

//Nodesx Nodesx::notSet( Node* node, Nodesx &ignore, bool drxn, bool inclNode )
//{
//    Nodesx nodes;
//    nodes.fillNot( node, ignore, drxn, inclNode );
//    return nodes;
//}
//
//Nodesx Nodesx::notSet( Node* node, Nodesx &ignore, int32_t limit, bool drxn, bool inclNode )
//{
//    Nodesx nodes;
//    nodes.fillNot( node, ignore, limit, drxn, inclNode );
//    return nodes;
//}

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

bool NodeRoll::cull( string s )
{
    for ( Node* node : nodes ) if ( node->seq_ == s )
    {
        erase( node );
        return true;
    }
    return false;
}

bool NodeRoll::find( Node* node )
{
    return std::find( nodes.begin(), nodes.end(), node ) != nodes.end();
}

void NodeRoll::flip( Nodes& flips )
{
    NodeRoll flipping;
    for ( int i = 0; i < size(); i++ ) if ( flips.find( nodes[i] ) )
    {
        flipping += nodes[i];
        nodes.erase( nodes.begin() + i-- );
    }
    vector<Node*> origins, fos;
    for ( Node* f : flipping.getGraph( 2 ).nodes )
    {
        cout << ">Origin "<< f->coverage_ << endl << f->seq_ << endl;
        fos.push_back( f );
        cout << ">rOrigin "<< f->coverage_ << endl << revCompNew( f->seq_ ) << endl;
    }
    for ( Node* f : flipping.getGraph( 2 ).nodes ) for ( auto& read : f->reads_ ) if ( !read.second.redundant )
    {
        ReadId id = params.getRevId( read.first );
        bool found = false;
        for ( Node* node : nodes ) if ( node->drxn_ < 2 && node->getRead( id ) )
        {
            node->setOrigin();
            origins.push_back( node );
        }
    }
    for ( int i = 0; i < flipping.size(); i++ )
    {
        flipping[i]->dismantle();
        delete flipping[i];
    }
    return;
    
    unordered_map<Node*, Node*> fresh;
    Nodes forks;
    
    for ( int i = 0; i < flipping.size(); i++ )
    {
        Node* q = flipping[i];
        bool done = false;
        string rc = revCompNew( q->seq_ );
        int32_t coords[2];
        for ( int d : { 0, 1 } ) if ( !done ) for ( Node* t : nodes ) if ( done = mapSeqEnd( rc, t->seq_, params.readLen, coords, d ) )
        {
            int32_t tCut[2]{ t->ends_[0] + coords[0], t->ends_[0] + coords[1] };
            int32_t qCut = d ? q->ends_[0] + coords[1] - coords[0] : q->ends_[1] - coords[1] + coords[0];
            for ( int j : { 0, 1 } ) if ( tCut[j] != t->ends_[j] )
            {
                assert( t->getNextReadCoord( tCut[j], j, !j ) );
                t = t->splitNode( tCut[j], *this, !j )[0];
            }
            if ( qCut != q->ends_[d] )
            {
                assert( q->getNextReadCoord( qCut, d, !d ) );
                for ( Node* node : q->splitNode( qCut, flipping, !d ) ) forks += node;
                i--;
            }
            else forks += q;
            
            break;
        }
        if ( !done ) for ( Node* t : nodes ) if ( rc.find( t->seq_ ) != string::npos ) assert( false );
        if ( !done )
        {
            Node* freshed = new Node( rc, q, 1 );
            if ( q->drxn_ == 2 ) freshed->setOrigin();
            fresh[q] = freshed;
        }
    }
    
    for ( const pair<Node*, Node*>& dupe : fresh )
    {
        for ( int d : { 0, 1 } ) for ( Edge& e : dupe.first->edges_[d] )
        {
            auto it = fresh.find( e.node );
            if ( forks.find( e.node ) )
            {
                assert( it == fresh.end() );
                string rc = revCompNew( e.node->seq_ );
                vector<Node*> hits;
                for ( Node* t : nodes ) if ( mapSeqOverlap( rc, t->seq_, min( (int)rc.size(), params.readLen ), !d ) ) hits.push_back( t );
                for ( Node* t : hits ) dupe.second->addEdge( t, e.ol, !d, false, e.leap );
                assert( !hits.empty() );
            }
            else
            {
                assert( it != fresh.end() );
                dupe.second->addEdge( it->second, e.ol, !d, false, e.leap );
            }
        }
    }
    
    for ( int i = 0; i < flipping.size(); i++ )
    {
        flipping[i]->dismantle();
        delete flipping[i];
    }
    for ( int d : { 0, 1 } ) for ( Node* node : nodes ) node->edgeTest( d );
    Node::updateStates( *this );
    Node::reverify( *this );
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

bool NodeRoll::getLine( string& s, vector<Node*>& line )
{
    for ( Edge& e : line.back()->edges_[1] ) 
    {
        if ( mapSeqOverlap( s, e.node->seq_, e.ol+1, 1 ) )
        {
            line.push_back( e.node );
            string seq = Node::getSeq( line );
            assert( seq.find( s ) != string::npos );
            if ( seq.find( s ) != string::npos ) return true;
            line.pop_back();
        }
        if ( s.find( e.node->seq_ ) != string::npos )
        {
            line.push_back( e.node );
            if ( getLine( s, line ) ) return true;
            line.pop_back();
        }
    }
    return false;
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

void NodeRoll::print( bool inclBad )
{
    
    struct NodePrint
    {
        string head, seq;
        int32_t coord;
    };
    
    vector<NodePrint> prints;
    int id = 0;
    Nodes used;
    for ( Node* node : nodes ) if ( used.add( node ) )
    {
        bool ignore = node->bad_;
        if ( ignore && inclBad ) for ( int d : { 0, 1 } ) for ( Edge& e : node->edges_[d] ) if ( !e.node->bad_ ) ignore = false;
        if ( ignore ) continue;
        
        vector<Node*> path = node->getPath();
        NodePrint np;
        int readCount = 0;
        for ( Node* n : path ) used += n ;
        for ( Node* n : path ) readCount += n->countReads( true );
        np.head = path.size() == 1 ? ">Node:" : ">Nodes:";
        for ( Node* n : path ) np.head += " " + to_string( ++id );
        np.head += "   ||   Reads: " + to_string( readCount );
        np.seq = Node::getSeq( path );
        np.coord = path[0]->ends_[0];
        prints.push_back( np );
    }
    
    sort( prints.begin(), prints.end(), []( NodePrint& a, NodePrint& b ){ return a.coord < b.coord; } );
    for ( NodePrint& np : prints ) cout << np.head << endl << string( np.coord - prints[0].coord, '-') + np.seq << endl; 
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

bool NodeRoll::purify( string s )
{
    int culls = 0, ol = 0;
    int32_t coords[2];
    for ( int i = 0; i < size(); i++ ) if ( ol = mapSeqOverlap( s, nodes[i]->seq_, params.readLen*2/3, 0 ) )
    {
        vector<Node*> line{ nodes[i] };
        if ( !getLine( s, line ) ) continue;
        while ( line.size() > 1 && mapSeqOverlap( s, line[1]->seq_, ol, 0 ) ) line.erase( line.begin() );
        for ( int j = 0; j < line.size(); j++ )
        {
            vector<Node*> bad[2];
            if ( j ) for ( Edge& e : line[j]->edges_[0] ) if ( e.node != line[j-1] ) bad[0].push_back( e.node );
            if ( j+1 < line.size() ) for ( Edge& e : line[j]->edges_[1] ) if ( e.node != line[j+1] ) bad[1].push_back( e.node );
            for ( int d : { 0, 1 } ) for ( Node* node : bad[d] ) line[j]->removeEdge( node, d, true );
            culls += bad[0].size() + bad[1].size();
        }
        for ( Node* f : Nodes( line.back(), 1, true ).nodes ) f->setVerified();
        return true;
//        Node* base = nodes[i];
//        Node* good[2]{ NULL, NULL };
//        vector<Node*> bad[2];
//        for ( int d : { 0, 1 } ) for ( Edge& e : nodes[i]->edges_[d] )
//        {
//            if ( s.find( e.node->seq_ ) != string::npos || mapSeqOverlap( s, e.node->seq_, params.readLen, d ) )
//            {
//                assert( !good[d] );
//                good[d] = e.node;
//            }
//            else bad[d].push_back( e.node );
//        }
//        assert( good[0] && good[1] );
//        for ( int d : { 0, 1 } ) if ( good[d] ) for ( Node* node : bad[d] ) nodes[i]->removeEdge( node, d, true );
//        culls += bad[0].size() + bad[1].size();
    }
    return false;
}

void NodeRoll::purify( vector<Node*>& base, vector<Node*>& tar, bool drxn )
{
    Nodes fwd, bck, good[2];
    for ( Node* node : base ) fwd.fill( node, drxn, true );
    for ( Node* node : base ) bck.fill( node, !drxn, true );
    for ( Node* node : tar ) good[0].fillIn( node, fwd, !drxn, true );
    for ( Node* node : tar ) good[0].fill( node, drxn, true );
    for ( Node* f : fwd.nodes )
    {
        if ( !good[0].find( f ) && !good[1].find( f ) && !bck.find( f ) ) erase( f );
        else for ( Edge e : f->edges( !drxn ) )
        {
            if ( !good[0].find( e.node ) && !good[1].find( e.node ) && !bck.find( e.node ) ) erase( e.node );
        }
    }
}

int NodeRoll::size()
{
    return nodes.size();
}

bool NodeRoll::solidify( string s )
{
    for ( Node* node : nodes ) if ( mapSeqOverlap( s, node->seq_, params.readLen, 0 ) )
    {
        vector<Node*> line{ node };
        assert( getLine( s, line ) );
        solidify( line );
        return true;
//        Node* clone = new Node( node, *this, node->drxn_, node->bad_ );
//        assert( false );
    }
    assert( false );
    return false;
}

bool NodeRoll::solidify( vector<Node*>& line )
{
    while ( line.size() > 1 && line[1]->edges_[0].size() == 1 ) line.erase( line.begin() );
    while ( line.size() > 1 && line.end()[-2]->edges_[1].size() == 1 ) line.pop_back();
    assert( line.size() > 2 );
    for ( int i = 1; i < line.size(); i++ ) assert( line[i]->getEdge( line[i-1], 0 ).node );
    vector<Node*> clones{ new Node( line[1], *this, line[1]->drxn_, line[1]->bad_ ) };
    for ( int i = 2; i+1 < line.size(); i++ )
    {
        Edge edge = line[i]->getEdge( line[i-1], 0 );
        edge.node = clones.back();
        clones.push_back( new Node( line[i], *this, line[i]->drxn_, line[i]->bad_ ) );
        clones.back()->addEdge( edge, 0, true );
    }
    clones[0]->addEdge( line[1]->getEdge( line[0], 0 ), 0, true );
    clones.back()->addEdge( line.end()[-2]->getEdge( line.back(), 1 ), 1, true );
    line[0]->removeEdge( line[1], 1, true );
    line.back()->removeEdge( line.end()[-2], 0, true );
    for ( Node* clone : clones ) clone->setVerified();
    test( true );
}

bool NodeRoll::solidify( string& s, vector<Node*>& line )
{
    for ( Edge& e : line.back()->edges_[1] ) if ( s.find( e.node->seq_ ) != string::npos )
    {
        line.push_back( e.node );
        string seq = Node::getSeq( line );
        if ( seq == s )
        {
            while ( line.size() > 1 && line[1]->edges_[0].size() == 1 ) line.erase( line.begin() );
            while ( line.size() > 1 && line.end()[-2]->edges_[1].size() == 1 ) line.pop_back();
            assert( line.size() == 3 );
            Edge edges[2]{ line[0]->getEdge( line[1], 1 ), line.back()->getEdge( line[1], 0 ) };
            assert( edges[0].node && edges[1].node );
            Node* clone = new Node( line[1], *this, line[1]->drxn_, line[1]->bad_ );
            line[0]->addEdge( clone, edges[0].ol, 1, false, edges[0].leap );
            line[2]->addEdge( clone, edges[1].ol, 0, false, edges[1].leap );
            line[0]->removeEdge( line[1], 1, true );
            line[2]->removeEdge( line[1], 0, true );
            clone->setVerified();
            return true;
        }
        else if ( solidify( s, line ) ) return true;
        line.pop_back();
    }
    return false;
}

void NodeRoll::summarise()
{
    CharId edges = 0, reads = 0, verifieds = 0, bads = 0;
    for ( Node* node : nodes )
    {
        for ( int d : { 0, 1 } ) edges += node->edges_[d].size();
        reads += node->reads_.size();
        if ( node->verified_ ) verifieds++;
        if ( node->bad_ ) bads++;
    }
    cout << "Node count:" << nodes.size() << endl;
    cout << "Edge count: " << edges << endl;
    cout << "Read count: " << reads << endl;
    cout << "Verified count: " << verifieds << endl;
    cout << "Bad count: " << bads << endl;
}

void NodeRoll::test( bool loop )
{
    if ( !loop ) return;
    Nodes drxn[2];
    for ( Node* n : nodes ) if ( n->drxn_ == 2 ) for ( int d : { 0, 1 } ) drxn[d].fill( n, d, true );
    for ( Node* n : nodes )
    {
        assert( n->pairs_.empty() );
        n->hits_.test();
        n->pairTest();
//        bool bad = n->drxn_ < 2;
        for ( Edge &e : n->edges_[0] ) if ( !e.leap )
        {
            string lSeq = e.node->seq_.substr( e.node->seq_.size() - e.ol );
            string rSeq = n->seq_.substr( 0, e.ol );
            assert( lSeq == rSeq || e.leap );
            assert( e.node->isEdge( n, 1 ) );
//            if ( n->drxn_ == 1 && ( e.node->drxn_ == 1 || e.node->drxn_ == 2 ) && !e.node->bad_ ) bad = false;
            assert( find( e.node ) );
        }
        for ( Edge &e : n->edges_[1] ) if ( !e.leap )
        {
            string lSeq = n->seq_.substr( n->seq_.size() - e.ol );
            string rSeq = e.node->seq_.substr( 0, e.ol );
            assert( lSeq == rSeq || e.leap );
            assert( e.node->isEdge( n, 0 ) );
//            if ( n->drxn_ == 0 && ( e.node->drxn_ == 0 || e.node->drxn_ == 2 ) && !e.node->bad_ ) bad = false;
            assert( find( e.node ) );
        }
//        if ( n->drxn_ < 2 && !n->bad_ ) assert( drxn[n->drxn_].find( n ) );
//        assert( !bad || n->bad_ );
//        assert( bad == n->bad_ );
    }
}

void NodeRoll::updateBad()
{
    Nodes good[2];
    for ( Node* node : nodes ) if ( node->drxn_ == 2 ) for ( int d : { 0, 1 } ) good[d].fill( node, d, true );
    for ( Node* node : nodes ) if ( !node->bad_ )
    {
        bool bad[2]{ !good[0].find( node ), !good[1].find( node ) };
        if ( bad[0] && bad[1] )
        {
            node->bad_ = true;
            node->verified_ = false;
            node->branch_ = false;
        }
        else for ( int d : { 0, 1 } ) if ( bad[d] && !bad[!d] && node->drxn_ == d )
        {
            node->drxn_ = !d;
            if ( !node->verified_ )
            {
                node->clearPaired( true );
                node->ends_.reset( node->drxn_ );
            }
        }
    }
//    for ( Node* node : nodes ) if ( !node->bad_ && !good[0].find( node ) && !good[1].find( node ) )
//    {
//        node->bad_ = true;
//        node->verified_ = false;
//        node->branch_ = false;
//    }
}

void NodeRoll::verify( int drxn )
{
    Nodes tar;
    for ( Node* node : nodes ) if ( node->drxn_ == 2 )
    {
        if ( drxn < 2 ) tar.fill( node, drxn, true );
        else assert( false );
    }
    for ( Node* node : tar.nodes ) node->setVerified();
}

NodeOffset::NodeOffset( int32_t dist )
: dists{ make_pair( dist, dist ) }
{}

bool NodeOffset::add( int32_t dist )
{
    for ( pair<int32_t, int32_t>& d : dists ) if ( d.first <= dist && dist <= d.second ) return false;
    
    if ( ( dists.size() > 1 ) && ( dists.back().second - dists[0].first > params.maxMpMean*1.5 ) ) return false;
    if ( dist < dists.back().second + 100 ) for ( int i = 0; i < dists.size(); i++ )
    {
        if ( dist < dists[i].first )
        {
            if ( dists[i].first - 100 < dist && dists[i].second - dists[i].first > 500 ) return false;
            if ( dists[i].first - 100 < dist ) dists[i].first = dist;
            else dists.insert( dists.begin() + i, make_pair( dist, dist ) );
            return true;
        }
        
        if ( dist < dists[i].second + 100 )
        {
            dists[i].second = dist;
            if ( i+1 < dists.size() && dists[i+1].first - 100 < dist )
            {
                dists[i].second = dists[i+1].second;
                dists.erase( dists.begin() + i+1 );
            }
            if ( dists[i].second - dists[i].first > 500 ) return false;
            return true;
        }
    }
    else if ( ( dists.size() > 1 ) && ( dists[0].first + params.maxMpMean*1.5 < dist ) ) return false;
    
    dists.push_back( make_pair( dist, dist ) );
    return true;
}

//bool NodeOffset::add( int32_t dist, bool minOnly )
//{
//    if ( dist == dists[0] || dist == dists[1] || dist == dists[2] ) return false;
//    else if ( abs( dist ) < abs( dists[0] ) )
//    {
//        dists[0] = dist;
//        return true;
//    }
//    else if ( minOnly ) return false;
//    else if ( abs( dists[2] ) < abs( dist ) )
//    {
//        dists[2] = dist;
//        return true;
//    }
//    else if ( min( abs( dists[0] - dists[1] ), abs( dists[1] - dists[2] ) ) < min( abs( dists[0] - dist ), abs( dist - dists[2] ) ) )
//    {
//        dists[1] = dist;
//        return true;
//    }
//    return false;
//}

int32_t NodeOffset::diff( int32_t est, bool drxn )
{
    int32_t best = drxn ? dists[0].first - est : est - dists[0].second;
    for ( pair<int32_t, int32_t>& d : dists )
    {
        if ( d.first <= est && est <= d.second ) return 0;
        int32_t diffs[2]{ drxn ? d.first - est : est - d.first, drxn ? d.second - est : est - d.second };
        for ( int i : { 0, 1 } ) if ( abs( diffs[i] ) < abs( best ) ) best = diffs[i];
    }
    return best;
}

//int32_t NodeOffset::min( bool drxn )
//{
//    return drxn ? dists[0].first : dists.back().second;
//}

//NodeOffsets::NodeOffsets( Node* node, bool drxn, bool inclNode )
//{
//    fill( node, 0, 0, drxn, drxn, inclNode );
//}

NodeOffsets::NodeOffsets( Node* node, int32_t limit, bool orient, bool drxn, bool inclNode )
{
    fill( node, 0, limit, orient, drxn, inclNode );
}

bool NodeOffsets::add( Node* node, int32_t dist )
{
    auto it = map.find( node );
    if ( it != map.end() ) return it->second.add( dist );
    map.insert( make_pair( node, NodeOffset( dist ) ) );
    return true;
}

//bool NodeOffsets::add( Node* fork, Edge& e, int32_t base, bool orient, bool drxn )
//{
//    int32_t dist = ( ( orient == drxn ? e.node : fork )->size() - e.ol ) * ( drxn ? 1 : -1 );
//    return add( e.node, dist + base );
//}

void NodeOffsets::clear()
{
    map.clear();
}

bool NodeOffsets::discontinue( int32_t dist, int32_t limit, bool drxn )
{
    if ( limit && ( drxn ? abs( limit ) < dist : dist < -abs( limit ) ) ) return true;
    return false;
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

void NodeOffsets::fill( Node* node, int32_t dist, int32_t limit, bool orient, bool drxn, bool inclNode )
{
    if ( inclNode && !add( node, dist ) ) return;
    if ( discontinue( dist, limit, drxn ) ) return;
    for ( Edge &e : node->edges_[drxn] ) fill( e.node, extend( node, e, dist, orient, drxn ), limit, orient, drxn, true );
}

void NodeOffsets::fillIn( Node* node, Nodes& include, int32_t dist, int32_t limit, bool orient, bool drxn, bool inclNode )
{
    if ( inclNode && ( !include.find( node ) || !add( node, dist ) ) ) return;
    if ( discontinue( dist, limit, drxn ) ) return;
    for ( Edge &e : node->edges_[drxn] ) fillIn( e.node, include, extend( node, e, dist, orient, drxn ), limit, orient, drxn, true );
}

void NodeOffsets::fillNot( Node* node, Nodes& ignore, int32_t dist, int32_t limit, bool orient, bool drxn, bool inclNode )
{
    if ( inclNode && ( ignore.find( node ) || !add( node, dist ) ) ) return;
    if ( discontinue( dist, limit, drxn ) ) return;
    for ( Edge &e : node->edges_[drxn] ) fillNot( e.node, ignore, extend( node, e, dist, orient, drxn ), limit, orient, drxn, true );
}

NodeOffset* NodeOffsets::get( Node* node )
{
    auto it = map.find( node );
    return it != map.end() ? &it->second : NULL;
}

NodeDists::NodeDists( Node* node, bool drxn, bool inclNode )
{
    fill( node, 0, 0, drxn, drxn, inclNode );
}

NodeDists::NodeDists( Node* node, int32_t limit, bool orient, bool drxn, bool inclNode )
{
    fill( node, 0, limit, orient, drxn, inclNode );
}

bool NodeDists::add( Node* node, int32_t dist, bool drxn )
{
    auto it = map.insert( make_pair( node, dist ) );
    if ( !it.second )
    {
        if ( drxn ? it.first->second <= dist : dist <= it.first->second ) return false;
        it.first->second = dist;
    }
    return true;
}

bool NodeDists::discontinue( int32_t dist, int32_t limit, bool drxn )
{
    if ( limit && ( drxn ? abs( limit ) < dist : dist < -abs( limit ) ) ) return true;
    return false;
}

bool NodeDists::empty()
{
    return map.empty();
}

void NodeDists::fill( Node* node, int32_t dist, int32_t limit, bool orient, bool drxn, bool inclNode )
{
    if ( inclNode && !add( node, dist, drxn ) ) return;
    if ( discontinue( dist, limit, drxn ) ) return;
    for ( Edge& e : node->edges_[drxn] ) fill( e.node, NodeOffsets::extend( node, e, dist, orient, drxn ), limit, orient, drxn, true );
}

void NodeDists::fillNot( Node* node, Nodes& ignore, int32_t dist, int32_t limit, bool orient, bool drxn, bool inclNode )
{
    if ( inclNode && ( ignore.find( node ) || !add( node, dist, drxn ) ) ) return;
    if ( discontinue( dist, limit, drxn ) ) return;
    for ( Edge& e : node->edges_[drxn] ) fillNot( e.node, ignore, NodeOffsets::extend( node, e, dist, orient, drxn ), limit, orient, drxn, true );
}

void NodeDists::fillReads( Node* node, int32_t dist, int readCount, int readLimit, bool orient, bool drxn )
{
    readCount += node->countReads( true );
    if ( !add( node, dist, drxn ) || readCount > readLimit ) return;
    for ( Edge& e : node->edges_[drxn] ) fillReads( e.node, NodeOffsets::extend( node, e, dist, orient, drxn ), readCount, readLimit, orient, drxn );
}

bool NodeDists::find( Node* node )
{
    return map.find( node ) != map.end();
}

int32_t* NodeDists::get( Node* node )
{
    auto it = map.find( node );
    return it != map.end() ? &it->second : NULL;
}

void NodeCounts::countOverlaps( Node* node, int ol, int minOl, bool drxn, bool inclNode )
{
    if ( ol <= 0 || ol < minOl || ( inclNode && !map.insert( make_pair( node, ol ) ).second ) ) return;
    ol -= node->size();
    for ( Edge& e : node->edges_[drxn] ) countOverlaps( e.node, ol + e.ol, minOl, drxn, true );
}

void NodeCounts::countReads( Node* node, bool drxn )
{
    int cur = 0;
    for ( Edge& e : node->edges_[!drxn] )
    {
        auto it = map.find( e.node );
        if ( it == map.end() ) return;
        cur = max( cur, it->second );
    }
    
    cur += node->countReads( true );
    if ( map.insert( make_pair( node, cur ) ).second ) for ( Edge& e : node->edges_[drxn] ) countReads( e.node, drxn );
}

int32_t* NodeCounts::get( Node* node )
{
    auto it = map.find( node );
    return it != map.end() ? &it->second : NULL;
}

