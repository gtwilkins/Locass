/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "node_pairs.h"
#include "node.h"

NodePair::NodePair( ReadId id, int32_t est, bool mp )
: libs( params.libs.size(), 0 ), count( 1 )
{
    pairs[mp].push_back( make_pair( id, est ) );
    libs[ getLib( id ) ]++;
}

void NodePair::add( ReadId id, int32_t est, bool mp )
{
    int i = 0;
    while ( i < pairs[mp].size() && pairs[mp][i].second <= est ) if ( pairs[mp][i++].first == id ) return;
    pairs[mp].insert( pairs[mp].begin() + i, make_pair( id, est ) );
    libs[ getLib( id ) ]++;
    count++;
}

int32_t NodePair::estimate()
{
    assert( count );
    bool mp = pairs[0].empty();
    int i = ( pairs[mp].size()-1 ) / 2, j = pairs[mp].size() / 2;
    return ( pairs[mp][i].second + pairs[mp][j].second ) / 2;
}

int NodePair::getLib( ReadId id )
{
    for ( int i = 0; i < params.libs.size(); i++ ) if ( id < params.libs[i].endCount ) return i;
    assert( false );
    return 0;
}

int32_t NodePair::margin()
{
    if ( !pairs[0].empty() ) return 100 + ( 100 * max( 0, 10 - (int)pairs[0].size() ) ) + ( params.maxPeMean / 3 );
    int best = 0;
    for ( int i = 1; i < libs.size(); i++ ) if ( libs[i] > libs[best] ) best = i;
    return 200 + params.libs[best].size / 3;
}

void NodePairs::add( Node* node, ReadId id, int32_t est, int32_t len, bool mp, bool drxn )
{
    auto it = pairs[drxn].find( node );
    if ( it != pairs[drxn].end() ) it->second.add( id, est, mp );
    else pairs[drxn].insert( make_pair( node, NodePair( id, est, mp ) ) );
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

//int NodePairs::get( Node* node, NodeOffsets& selfOffs, NodeOffsets& tarOffs, bool drxn )
//{
//    auto it = pairs[drxn].find( node );
//    if ( it == pairs[drxn].end() ) return 0;
//    NodeOffset* off[2]{ selfOffs.get( self ), tarOffs.get( node ) };
//    assert( off[0] && off[1] );
//    if ( !off[0] || !off[1] ) return 0;
//    int32_t est = it->second.sum / max( 1, it->second.count );
//    int32_t cutoff = min( 200 + ( it->second.maxLen / 3 ), 1000 );
//    int32_t dists[2]{ abs( (*off[0])[0] - (*off[1])[1] ), abs( (*off[0])[1] - (*off[1])[0] ) };
//    for ( int i = 0; i < 2; i++ ) if ( abs( est - dists[i] ) < cutoff ) return it->second.count;
//    return 0;
//}

int NodePairs::get( Nodes& nodes, bool drxn )
{
    int score = nodes.find( self ) ? count : 0;
    for ( Node* node : nodes.nodes ) score += get( node, drxn );
    return score;
}

//int NodePairs::get( Nodes& nodes, NodeOffsets& selfOffs, NodeOffsets& tarOffs, bool drxn )
//{
//    int score = 0;
//    for ( Node* node : nodes.nodes ) score += get( node, selfOffs, tarOffs, drxn );
//    return score;
//}

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
    unordered_set<ReadId> ids;
    for ( int mp : { 0, 1 } ) for ( pair<ReadId, int32_t> rp : np->pairs[mp] ) ids.insert( rp.first );
    vector<Node*> tars = fwd->clones();
    for ( Node* b : self->clones() ) for ( Node* f : tars ) if ( np = b->hits_.getPair( f, drxn ) )
    {
        for ( int mp : { 0, 1 } ) for ( pair<ReadId, int32_t> rp : np->pairs[mp] ) if ( ids.find( rp.first ) != ids.end() )
        {
            q += b;
            t += f;
            break;
        }
    }
//    int32_t est = np->estimate(), cutoff = 200 + np->maxLen / 10;
//    vector<Node*> tars = fwd->clones();
//    for ( Node* b : self->clones() )
//    {
//        for ( Node* f : tars )
//        {
//            if ( !( np = b->hits_.getPair( f, drxn ) ) || abs( est - np->estimate() ) > cutoff ) continue;
//            q += b;
//            t += f;
//        }
//    }
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