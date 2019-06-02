/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "node_pairs.h"
#include "node.h"

int32_t NodePair::estimate()
{
    assert( count );
    return sum / count;
}

void NodePairs::add( Node* node, int32_t est, int32_t len, bool isPe, bool drxn )
{
    auto it = pairs[drxn].find( node );
    if ( it != pairs[drxn].end() )
    {
        it->second.maxLen = max( it->second.maxLen, len );
        it->second.sum += est;
        it->second.count++;
    }
    else pairs[drxn].insert( make_pair( node, NodePair( est, len ) ) );
}

void NodePairs::clean()
{
    for ( int d : { 0, 1 } )
    {
        for ( auto& np : pairs[d] ) np.first->hits_.pairs[!d].erase( self );
        pairs[d].clear();
    }
    count = 0;
}

bool NodePairs::empty()
{
    return !count && pairs[0].empty() && pairs[1].empty();
}

void NodePairs::erase( Node* node, bool drxn )
{
    pairs[drxn].erase( node );
}

int NodePairs::get()
{
    int hits = count;
    for ( int d : { 0, 1 } ) for ( auto& np : pairs[d] ) hits += np.second.count;
    return hits;
}

int NodePairs::get( bool drxn )
{
    int hits = 0;
    for ( auto& np : pairs[drxn] ) hits += np.second.count;
    return hits;
}

int NodePairs::get( Node* node, bool drxn )
{
    auto it = pairs[drxn].find( node );
    return it == pairs[drxn].end() ? 0 : it->second.count;
}

int NodePairs::get( Node* node, NodeOffsets& selfOffs, NodeOffsets& tarOffs, bool drxn )
{
    auto it = pairs[drxn].find( node );
    if ( it == pairs[drxn].end() ) return 0;
    NodeOffset* off[2]{ selfOffs.get( self ), tarOffs.get( node ) };
    assert( off[0] && off[1] );
    if ( !off[0] || !off[1] ) return 0;
    int32_t est = it->second.sum / max( 1, it->second.count );
    int32_t cutoff = min( 200 + ( it->second.maxLen / 3 ), 1000 );
    int32_t dists[2]{ abs( (*off[0])[0] - (*off[1])[1] ), abs( (*off[0])[1] - (*off[1])[0] ) };
    for ( int i = 0; i < 2; i++ ) if ( abs( est - dists[i] ) < cutoff ) return it->second.count;
    return 0;
}

int NodePairs::get( Nodes& nodes, bool drxn )
{
    int score = nodes.find( self ) ? count : 0;
    for ( Node* node : nodes.nodes ) score += get( node, drxn );
    return score;
}

int NodePairs::get( Nodes& nodes, NodeOffsets& selfOffs, NodeOffsets& tarOffs, bool drxn )
{
    int score = 0;
    for ( Node* node : nodes.nodes ) score += get( node, selfOffs, tarOffs, drxn );
    return score;
}

NodePair* NodePairs::getPair( Node* node, bool drxn )
{
    auto it = pairs[drxn].find( node );
    return it == pairs[drxn].end() ? NULL : &it->second;
}

void NodePairs::reset()
{
    pairs[0].clear();
    pairs[1].clear();
    count = 0;
}

void NodePairs::setRedundant( Nodes& q, Nodes& t, Node* fwd, bool drxn )
{
    NodePair* np = getPair( fwd, drxn );
    assert( np );
    int32_t est = np->estimate(), cutoff = 200 + np->maxLen / 10;
    vector<Node*> tars = fwd->clones();
    for ( Node* b : self->clones() )
    {
        for ( Node* f : tars )
        {
            if ( !( np = b->hits_.getPair( f, drxn ) ) || abs( est - np->estimate() ) > cutoff ) continue;
            q += b;
            t += f;
        }
    }
}

void NodePairs::test()
{
    for ( int d : { 0, 1 } )
    {
        for ( auto& np : pairs[d] )
        {
            auto it = np.first->hits_.pairs[!d].find( self );
            assert( it != np.first->hits_.pairs[!d].end() && it->second.count == np.second.count );
        }
    }
}