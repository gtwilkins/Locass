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

#include "node_path.h"
#include "correct_read.h"
#include "seed_fork.h"
#include "leap.h"
#include <algorithm>

PathPairing::PathPairing( NodePath* l, NodePath* r, int32_t diff )
: diffs{ diff }, missed{ 0 }, score( 0 )
{
    node[0] = l;
    node[1] = r;
    for ( int d : { 0, 1 } ) node[d]->pairs_[!d].push_back( this );
}

PathPairing::~PathPairing()
{
    for ( int d : { 0, 1 } ) node[d]->pairs_[!d].erase( remove( node[d]->pairs_[!d].begin(), node[d]->pairs_[!d].end(), this ), node[d]->pairs_[!d].end() );
}

PathEdge::PathEdge( NodePath* fork, NodePath* branch, Edge& e, int32_t diff, bool drxn )
: diff( diff ), ol( e.ol ), score( -1 ), multi( 1 ), leap( e.leap )
{
    edge[!drxn] = fork;
    edge[drxn] = branch;
    for ( int d : { 0, 1 } ) edge[d]->edges_[!d].push_back( this );
}

PathEdge::PathEdge( NodePath* l, NodePath* r, PathEdge* pe )
: diff( pe->diff ), ol( pe->ol ), score( pe->score ), multi( 1 ), leap( pe->leap )
{
    edge[0] = l;
    edge[1] = r;
    for ( int d : { 0, 1 } ) edge[d]->edges_[!d].push_back( this );
    edge[0]->path_.back()->addEdge( edge[1]->path_[0], ol, 1, false, leap );
}

PathEdge::~PathEdge()
{
    for ( int d : { 0, 1 } ) edge[d]->edges_[!d].erase( remove( edge[d]->edges_[!d].begin(), edge[d]->edges_[!d].end(), this ), edge[d]->edges_[!d].end() );
    for ( int d : { 0, 1 } ) edge[d]->breaks_[!d].erase( remove( edge[d]->breaks_[!d].begin(), edge[d]->breaks_[!d].end(), this ), edge[d]->breaks_[!d].end() );
}

void PathEdge::claim( NodePath* np, bool drxn )
{
    edge[0]->path_.back()->removeEdge( edge[1]->path_[0], 1, true );
    edge[drxn]->edges_[!drxn].erase( remove( edge[drxn]->edges_[!drxn].begin(), edge[drxn]->edges_[!drxn].end(), this ), edge[drxn]->edges_[!drxn].end() );
    edge[drxn] = np;
    np->edges_[!drxn].push_back( this );
    edge[0]->path_.back()->addEdge( edge[1]->path_[0], ol, 1, false, leap );
}

void PathEdge::downgrade()
{
    for ( int d : { 0, 1 } )
    {
        edge[d]->edges_[!d].erase( remove( edge[d]->edges_[!d].begin(), edge[d]->edges_[!d].end(), this ), edge[d]->edges_[!d].end() );
        edge[d]->breaks_[!d].push_back( this );
    }
}

bool PathEdge::sever( NodePath* l, NodePath* r )
{
    for ( PathEdge* pe : l->edges_[1] ) if ( pe->edge[1] == r )
    {
        l->path_.back()->removeEdge( r->path_[0], 1, true );
        delete pe;
        return true;
    }
    assert( false );
    return false;
}

bool PathPair::add( Node* base, NodeRoll* tar, unordered_map<Node*, NodePath*>& mapped, ReadId id, bool d )
{
    Coords* coords;
    int32_t* off;
    if ( tar ) for ( Node* t : tar->nodes ) if ( t != base && ( coords = t->getRead( id ) ) )
    {
        auto it = mapped.find( t );
        if ( it == mapped.end() ) return false;
        assert( off = it->second->get( t ) );
        int32_t mark = *off + (*coords)[d] - t->ends_[1];
        bool marked = false;
        for ( pair<NodePath*, int32_t>& m : marks_[d] ) if ( m.first == it->second && m.second == mark ) marked = true;
        if ( !marked ) marks_[d].push_back( make_pair( it->second, mark ) );
    }
    return true;
}

bool PathPair::confirm( bool cull )
{
    int32_t cutoff = dist_ * 0.2 + 300;
    for ( int i = 0; i < marks_[0].size(); i++ ) for ( int j = 0; j < marks_[1].size(); j++ ) if ( marks_[0][i].first == marks_[1][j].first )
    {
        if ( marks_[1][j].second - marks_[0][i].second - dist_ < cutoff ) return false;
    }
    best_ = dist_ * 1.2 + 300;
    vector<int32_t> good[2]{ vector<int32_t>( marks_[0].size(), false ), vector<int32_t>( marks_[1].size(), false ) };
    vector<PathPairs*> pairs, unpairs;
    for ( int i = 0; i < good[0].size(); i++ ) for ( int j = 0; j < good[1].size(); j++ )
    {
        for ( PathPairs* pp : marks_[0][i].first->paired_[1] ) if ( pp->node_[1] == marks_[1][j].first )
        {
            int32_t dist = marks_[1][j].second - pp->diff_ - marks_[0][i].second - dist_;
            if ( dist <= cutoff )
            {
                good[0][i] = good[1][j] = true;
                if ( abs( dist ) < abs( best_ ) ) best_ = dist;
                pairs.push_back( pp );
            }
            else if ( cull ) unpairs.push_back( pp );
            break;
        }
    }
    
    for ( PathPairs* pp : unpairs ) if ( find( pairs.begin(), pairs.end(), pp ) == pairs.end() ) pp->shared_.erase( remove( pp->shared_.begin(), pp->shared_.end(), this ), pp->shared_.end() );
    
    if ( pairs.empty() ) return false;
    
    for ( int d : { 0, 1 } ) for ( int i = good[d].size(); i-- > 0; ) if ( !good[d][i] ) marks_[d].erase( marks_[d].begin() + i );
    
    assert( marks_[0].empty() == marks_[1].empty() );
    
    if ( marks_[0].size() == 1 && marks_[1].size() == 1 )
    {
        assert( pairs.size() == 1 );
        if ( cull ) pairs[0]->shared_.erase( remove( pairs[0]->shared_.begin(), pairs[0]->shared_.end(), this ), pairs[0]->shared_.end() );
        pairs[0]->uniques_.push_back( PathUnique( id_, marks_[0][0].second, marks_[1][0].second, dist_, best_ ) );
        return false;
    }
    
    if ( !cull ) for ( PathPairs* pp : pairs ) pp->shared_.push_back( this );
    
    return true;
}

//vector<PathPairs*> PathPair::get()
//{
//    vector<PathPairs*> pairs
//    for ( int i = 0; i < marks_[0].size(); i++ );
//    return pairs;
//}

PathPairs::PathPairs( NodePath* l, NodePath* r, int32_t diff, bool forked )
: diff_( diff ), score_( 0 ), forked_( forked )
{
    node_[0] = l;
    node_[1] = r;
    for ( int d : { 0, 1 } ) node_[d]->paired_[!d].push_back( this );
}

PathPairs::PathPairs( PathPairs* pp, NodePath* np, bool drxn )
: diff_( pp->diff_ ), score_( 0 ), forked_( pp->forked_ )
{
    node_[!drxn] = pp->node_[!drxn];
    node_[drxn] = np;
    for ( int d : { 0, 1 } ) node_[d]->paired_[!d].push_back( this );
    for ( PathPair* p : pp->shared_ ) for ( int i = 0; i < p->marks_[drxn].size(); i++ ) if ( p->marks_[drxn][i].first == pp->node_[drxn] )
    {
        p->marks_[drxn].push_back( make_pair( np, p->marks_[drxn][i].second ) );
        shared_.push_back( p );
    }
    for ( PathUnique& pu : pp->uniques_ )
    {
        PathPair* p = new PathPair( pu.id_, pu.dist_ );
        for ( int d : { 0, 1 } ) p->marks_[d].push_back( make_pair( pp->node_[d], d ? pu.r_ : pu.l_ ) );
        p->marks_[drxn].push_back( make_pair( np, drxn ? pu.r_ : pu.l_ ) );
        pp->shared_.push_back( p );
        shared_.push_back( p );
    }
    pp->uniques_.clear();
}

PathPairs::~PathPairs()
{
    for ( int d : { 0, 1 } ) node_[d]->paired_[!d].erase( remove( node_[d]->paired_[!d].begin(), node_[d]->paired_[!d].end(), this ), node_[d]->paired_[!d].end() );
}

bool PathPairs::add( NodePath* l, NodePath* r, int32_t diff, bool forked )
{
    for ( PathPairs* pp : l->paired_[1] ) if ( pp->node_[1] == r )
    {
        if ( forked ) pp->forked_ = true;
        if ( diff <= pp->diff_ ) return false;
        pp->diff_ = diff;
        return true;
    }
    new PathPairs( l, r, diff, forked );
    return true;
}

void PathPairs::discard()
{
    for ( PathPair* p : shared_ )
    {
        vector<bool> good( p->marks_[1].size(), false );
        for ( int i = 0; i < p->marks_[0].size(); i++ )
        {
            bool bad = true;
            for ( PathPairs* pp : p->marks_[0][i].first->paired_[1] ) if ( pp != this )
            {
                for ( int j = 0; j < p->marks_[1].size(); j++ ) if ( pp->node_[1] == p->marks_[1][j].first )
                {
                    good[j] = true;
                    bad = false;
                }
            }
            if ( bad ) p->marks_[0].erase( p->marks_[0].begin() + i-- );
        }
        for ( int i = p->marks_[1].size(); i-- > 0; ) if ( !good[i] ) p->marks_[1].erase( p->marks_[1].begin() + i );
        assert( p->marks_[0].empty() == p->marks_[1].empty() );
        if ( p->marks_[0].empty() ) delete p;
    }
    delete this;
}

NodePath::NodePath( NodePath* np, vector<NodePath*>& paths, NodeRoll& nodes, int multi )
: offs_( np->offs_ ), verified_( np->verified_ ), id_( np->id_ ), multi_( multi ), loop_( np->loop_ )
{
    for ( int i = 0; i < np->path_.size(); i++ )
    {
        Node* node = np->path_[i];
        path_.push_back( new Node( node, nodes, node->drxn_, node->bad_ ) );
        if ( i ) for ( Edge& e : node->edges_[0] ) path_.back()->addEdge( e.node == np->path_[i-1] ? path_[i-1] : e.node, e.ol, 0, false, e.leap );
        if ( path_.size() < np->path_.size() ) for ( Edge& e : node->edges_[1] ) if ( e.node != np->path_[i+1] ) path_.back()->addEdge( e, 1, true );
        if ( i ) assert( path_[i]->getEdge( path_[i-1], 0 ).node );
    }
    NodeRoll edges[2]{ NodeRoll::next( np->path_[0], 0 ), NodeRoll::next( np->path_.back(), 1 ) };
    for ( Edge& e : np->path_[0]->edges_[0] ) if ( !edges[0].find( e.node ) ) assert( false );
    for ( Edge& e : np->path_[0]->edges_[0] ) if ( !edges[0].find( e.node ) ) path_[0]->addEdge( e, 0, true );
    for ( Edge& e : np->path_.back()->edges_[1] ) if ( !edges[1].find( e.node ) ) assert( false );
    for ( Edge& e : np->path_.back()->edges_[1] ) if ( !edges[1].find( e.node ) ) path_.back()->addEdge( e, 1, true );
    
    paths.push_back( this );
    for ( int d : { 0, 1 } ) coords_[d] = np->coords_[d];
    for ( int d : { 0, 1 } ) ends_[d] = np->ends_[d];
    for ( int d : { 0, 1 } ) ended_[d] = np->ended_[d];
    for ( int d : { 0, 1 } ) for ( PathPairs* pp : np->paired_[d] ) new  PathPairs( pp, this, !d );
    np->multi_ = max( 1, np->multi_ - multi );
}

NodePath::NodePath( Node* seed, vector<NodePath*>& paths, int32_t coord )
: offs_{ make_pair( seed, coord ) }, path_{ seed }, verified_( seed->verified_ ), multi_( 0 ), loop_( 0 )
{
    paths.push_back( this );
    coords_[0] = coords_[1] = ends_[1] = coord;
    ends_[0] = ends_[1] - seed->size();
    ended_[0] = ended_[1] = false;
}

NodePath::~NodePath()
{
    for ( int d : { 0, 1 } )
    {
        while ( !edges_[d].empty() ) delete edges_[d].back();
        while ( !breaks_[d].empty() ) delete breaks_[d].back();
        while ( !pairs_[d].empty() ) delete pairs_[d].back();
        while ( !paired_[d].empty() ) delete paired_[d].back();
    }
}

void NodePath::addEdge( Edge& e, vector<NodePath*>& paths, bool branch, bool drxn )
{
    ended_[drxn] = true;
    for ( PathEdge* pe : edges_[drxn] ) if ( e.node == ( drxn ? pe->edge[1]->path_[0] : pe->edge[0]->path_.back() ) ) return;
    
    bool added = false;
    for ( NodePath* np : paths ) if ( np->findFork( e.node, !drxn ) && ( added = true ) )
    {
        new PathEdge( this, np, e, ( drxn ? np->coords_[0] - getCoord( e, 1 ) : getCoord( e, 0 ) - np->coords_[1] ), drxn );
    }
    if ( !added )
    {
        NodePath* edge = new NodePath( e.node, paths, getCoord( e, drxn ) );
        new PathEdge( this, edge, e, 0, drxn );
        edge->extend( paths, branch, drxn );
    }
}

bool NodePath::addPair( NodePath* np, int32_t diff )
{
    for ( PathPairing* pp : pairs_[1] ) if ( pp->node[1] == np )
    {
        int i = 0;
        while ( i < pp->diffs.size() && pp->diffs[i] < diff ) i++;
        if ( pp->diffs[i] == diff ) return false;
        pp->diffs.insert( pp->diffs.begin() + i, diff );
        pp->missed.insert( pp->missed.begin() + i, 0 );
        return true;
    }
    new PathPairing( this, np, diff );
    return true;
}

void NodePath::cleanPairs( vector<NodePath*>& paths )
{
    for ( NodePath* np : paths ) for ( int i = 0; i < np->paired_[1].size(); i++ )
    {
        if ( !np->paired_[1][i]->score_ && np->paired_[1][i]->uniques_.empty() && np->paired_[1][i]->shared_.empty() ) delete np->paired_[1][i--];
    }
}

void NodePath::create( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& seeds, vector<NodePath*>& paths )
{
    for ( NodePath* np : paths ) delete np;
    seeds.clear();
    paths.clear();
    
    nodes.summarise();
    for ( Node* node : nodes.getGraph( 2 ).nodes )
    {
        bool used = false;
        for ( NodePath* np : paths ) if ( used = np->findNode( node ) )
        {
            if ( find( seeds.begin(), seeds.end(), np ) == seeds.end() ) seeds.push_back( np );
            break;
        }
        if ( !used )
        {
            seeds.push_back( new NodePath( node, paths, node->ends_[1] ) );
            seeds.back()->extend( paths, false, 1 );
            for ( int d : { 0, 1 } ) seeds.back()->extend( paths, true, d );
        }
    }
    
    for ( int d : { 0, 1 } ) for ( NodePath* np : paths ) if ( np->edges_[d].empty() ) np->fill( d ? np->path_.back() : np->path_[0], np->coords_[d], 500, d );
    for ( NodePath* np : paths )
    {
        vector<NodePath*> path;
        np->setPairs( np, path, 0, params.maxMpMean * 1.2 + 300, false, false );
    }
    for ( NodePath* np : paths )
    {
        vector<NodePath*> path;
        np->setPairs( np, path, 0, false );
    }
    for ( NodePath* np : paths ) np->setScore();
    for ( NodePath* np : paths ) np->setScores();
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) np->reduce( d );
    for ( NodePath* np : paths ) np->setMulti();
    setReads( bwt, nodes, paths );
}

bool NodePath::doesReach( unordered_set<NodePath*>& include, int32_t coord, bool drxn )
{
    if ( include.find( this ) == include.end() ) return false;
    if ( drxn ? coord < coords_[1] : coords_[0] - path_[0]->size() < coord ) return true;
    for ( PathEdge* pe : edges_[drxn] ) if ( pe->edge[drxn]->doesReach( include, coord, drxn ) ) return true;
    return false;
}

void NodePath::extend( vector<NodePath*>& paths, bool branch, bool drxn )
{
    if ( ended_[drxn] ) return;
    vector<Edge> edges = ( drxn ? path_.back() : path_[0] )->edges_[drxn];
    bool verified = ( drxn ? path_.back() : path_[0] )->verified_;
    while ( edges.size() == 1 && edges_[drxn].empty() && edges[0].node->edges_[!drxn].size() == 1 )
    {
//        if ( edges[0].node->edges_[!drxn].size() > 1 )
//        {
//            if ( branch && verified ) addEdge( edges[0], paths, branch, drxn );
//            return;
//        }
        if ( !( verified = edges[0].node->verified_ ) ) verified_ = false;
        coords_[drxn] = getCoord( edges[0], drxn );
        ends_[drxn] = drxn ? coords_[1] : coords_[0] - edges[0].node->size();
        path_.insert( drxn ? path_.end() : path_.begin(), edges[0].node );
        assert( offs_.insert( make_pair( edges[0].node, coords_[drxn] ) ).second );
        edges = edges[0].node->edges_[drxn];
    }
    
    if ( edges.empty() ) ended_[drxn] = true;
    if ( !branch || !verified || edges.empty() ) return;
    
    for ( Edge& e : edges ) addEdge( e, paths, branch, drxn );
    for ( PathEdge* pe : edges_[drxn] ) pe->edge[drxn]->extend( paths, true, !drxn );
}

void NodePath::fill( Node* node, int32_t dist, int32_t limit, bool drxn )
{
    for ( Edge& e : node->edges_[drxn] )
    {
        int32_t coord = getCoord( node, e, dist, drxn );
        auto it = offs_.insert( make_pair( e.node, coord ) );
        if ( !it.second )
        {
            if ( drxn ? it.first->second <= coord : coord <= it.first->second ) continue;
            it.first->second = coord;
        }
        if ( limit > 0 ) fill( e.node, coord, limit - e.node->size() + e.ol, drxn );
    }
}


bool NodePath::findFork( Node* node, bool drxn )
{
    Node* fork = drxn ? path_.back() : path_[0];
    if ( node == fork ) return true;
    if ( node->edges_[!drxn].empty() && fork->isClone( node ) ) return true;
    return false;
}

bool NodePath::findNode( Node* node )
{
    if ( find( path_.begin(), path_.end(), node ) != path_.end() ) return true;
    assert( !path_.empty() );
    if ( node->cloned_ ) for ( int d : { 0, 1 } ) if ( node->edges_[d].empty() )
    {
        for ( Node* p : path_ ) if ( p->isClone( node ) && p->edges_[!d].empty() ) return true;
        for ( Node* p : path_ ) if ( p->isClone( node ) && p->edges_[!d].empty() ) assert( false );
        if ( findFork( node, !d ) ) assert( false );
        if ( findFork( node, !d ) ) return true;
    }
    return false;
}

int32_t* NodePath::get( Node* node )
{
    auto it = offs_.find( node );
    return it != offs_.end() ? &it->second : NULL;
}

//int32_t NodesPath::getCoord( bool drxn )
//{
//    return drxn ? coords_[1] : coords_[0] - path_[0]->size() ;
//}

int32_t NodePath::getLen( NodePath* rhs, int32_t diff, bool drxn )
{
    return ( drxn ? rhs : this )->ends_[1] - ( drxn ? this : rhs )->ends_[0] - diff;
}

int32_t NodePath::getCoord( Edge& e, bool drxn )
{
    return drxn ? coords_[1] + e.node->size() - e.ol : coords_[0] - path_[0]->size() + e.ol;
}

int32_t NodePath::getCoord( Node* node, Edge& e, int32_t dist, bool drxn )
{
    return dist + ( drxn ? e.node->size() - e.ol : -node->size() + e.ol );
}

Node* NodePath::getEnd( Querier& bwt, NodeRoll& nodes, unordered_set<NodePath*>& used, bool drxn )
{
    assert( edges_[drxn].empty() );
    for ( auto& no : offs_ ) if ( no.first->isContinue( drxn ) ) return drxn ? path_.back() : path_[0];
    unordered_map<Node*, int32_t> offs;
    getOffsets( offs, used, ends_[drxn], 0, max( params.maxPeMax, params.maxMpMean ), drxn );
    return Leap::leapEnd( bwt, nodes, drxn ? path_.back() : path_[0], offs, drxn );
}

int NodePath::getMulti( bool drxn )
{
    int multi = 0, verified = verified_;
    for ( PathEdge* pe : edges_[drxn] )
    {
        if ( pe->score > 3 || pe->edge[drxn]->verified_ ) verified = 1;
        if ( pe->edge[drxn]->edges_[!drxn].size() == 1 ) multi += pe->edge[drxn]->getMulti( drxn );
    }
    return max( verified, multi );
}

int32_t NodePath::getOffset( NodePath* np, bool drxn )
{
    for ( PathEdge* pe : edges_[drxn] ) if ( pe->edge[drxn] == np ) return pe->diff;
    return 0;
}

void NodePath::getOffsets( unordered_map<Node*, int32_t>& offs, unordered_set<NodePath*>& include, int32_t base, int32_t diff, int32_t limit, bool drxn )
{
    if ( include.find( this ) == include.end() ) return;
    
    bool added = false;
    for ( auto& no : offs_ )
    {
        int32_t off = ( drxn ? base - no.second : no.second - no.first->size() - base ) - diff;
        auto it = offs.insert( make_pair( no.first, off ) );
        if ( !it.second && off < it.first->second && ( it.second = true ) ) it.first->second = off;
        if ( it.second && off < limit ) added = true;
    }
    
    if ( added ) for ( PathEdge* pe : edges_[!drxn] ) pe->edge[!drxn]->getOffsets( offs, include, base, diff+pe->diff, limit, drxn );
}

int32_t NodePath::getOverlap( NodePath* np, bool drxn )
{
    for ( PathEdge* pe : edges_[drxn] ) if ( pe->edge[drxn] == np ) return pe->ol;
    return 0;
}

bool NodePath::isBranchable( bool drxn )
{
    if ( edges_[drxn].size() < 2 ) return false;
    for ( PathEdge* pe : edges_[drxn] ) if ( pe->score < 2 || ( pe->score < 4 && !pe->edge[drxn]->verified_ ) ) return false;
    return true;
}

void NodePath::print( vector<NodePath*>& paths )
{
    int32_t minCoord = 0;
    for ( int i = 0; i < paths.size(); i++ ) paths[i]->id_ = i+1;
    for ( NodePath* np : paths ) minCoord = min( minCoord, np->ends_[0] );
//    for ( NodePath* np : paths )
//    {
//        string seq = Node::getSeq( np->path_ );
//        float cover = Node::getCoverage( np->path_ );
//        cout << ">" << np->id_ << ") Coords: (" << np->ends_[0] << ")-(" << np->ends_[1] << "), Coverage: " << cover;
//        cout << ", Left edges:" << ( np->edges_[0].empty() ? " NONE" : "" );
//        for ( PathEdge* pe : np->edges_[0] ) cout << " (" << pe->edge[0]->id_ << ")";
//        cout << ", Right edges:" << ( np->edges_[1].empty() ? " NONE" : "" );
//        for ( PathEdge* pe : np->edges_[1] ) cout << " (" << pe->edge[1]->id_ << ")";
//        if ( !np->verified_ ) cout << " UNVERIFIED";
//        cout << endl << string( np->ends_[0] - minCoord, '-' )<< seq << endl;
//    }
}

void NodePath::reduce( bool drxn )
{
    if ( edges_[drxn].size() < 2 ) return;
    int limits[2]{ 10, 0 };
    unordered_set<NodePath*> t;
    setBranch( t, params.maxPeMean, !drxn );
    for ( PathEdge* pe : edges_[drxn] ) if ( pe->score < 0 )
    {
        pe->score = 0;
        unordered_set<NodePath*> q;
        pe->edge[drxn]->setBranch( q, params.maxPeMean, drxn );
        for ( NodePath* np : q ) for ( PathPairing* pp : np->pairs_[!drxn] ) if ( t.find( pp->node[!drxn] ) != t.end() ) pe->score += pp->score;
    }
    for ( PathEdge* pe : edges_[drxn] ) limits[0] = min( limits[0], pe->score );
    for ( PathEdge* pe : edges_[drxn] ) limits[1] = max( limits[1], pe->score );
    if ( limits[1] > 9 ) for ( int i = 0; i < edges_[drxn].size(); i++ ) if ( !edges_[drxn][i]->score ) edges_[drxn][i--]->downgrade();
}

void NodePath::reset( vector<NodePath*>& paths )
{
    vector<PathPairs*> repairs;
    for ( NodePath* np : paths )
    {
        np->loop_ = 0;
        unordered_map<NodePath*, int32_t> diffs;
        vector<NodePath*> pathed;
        np->resetDiffs( np, pathed, diffs, 0, params.maxMpMean*1.2 + 300 );
        for ( int i = 0; i < np->paired_[1].size(); i++ )
        {
            auto it = diffs.find( np->paired_[1][i]->node_[1] );
            if ( it != diffs.end() )
            {
                if ( abs( np->paired_[1][i]->diff_ - it->second ) > 100 )
                {
                    repairs.push_back( np->paired_[1][i] );
                }
                np->paired_[1][i]->diff_ = it->second;
            }
            else np->paired_[1][i--]->discard();
        }
    }
    
    unordered_set<ReadId> repaired;
    for ( PathPairs* pp : repairs )
    {
        vector<PathPair*> shared;
        for ( PathPair* p : pp->shared_ ) if ( repaired.insert( p->id_ ).second ) shared.push_back( p );
        for ( PathPair* p : shared ) if ( !p->confirm( true ) ) delete p;
    }
    
    cleanPairs( paths );
}

void NodePath::resetDiffs( NodePath* np, vector<NodePath*>& path, unordered_map<NodePath*, int32_t>& diffs, int32_t diff, int32_t limit )
{
    if ( !path.empty() && np == this && ( !loop_ || -diff < loop_ ) ) loop_ = -diff;
    if ( find( path.begin(), path.end(), np ) != path.end() ) return;
    
    auto ins = diffs.insert( make_pair( np, diff ) );
    if ( !ins.second )
    {
        if ( diff <= ins.first->second ) return;
        ins.first->second = diff;
    }
    
    if ( !path.empty() ) limit -= np->size();
    path.push_back( np );
    if ( limit > 0 ) for ( PathEdge* pe : np->edges_[1] ) resetDiffs( pe->edge[1], path, diffs, diff + pe->diff, limit + pe->ol );
    path.pop_back();
}

void NodePath::setBranch( unordered_set<NodePath*>& branch, int32_t dist, bool drxn )
{
    if ( !branch.insert( this ).second ) return;
    dist -= ends_[1] - ends_[0];
    if ( dist > 0 ) for ( PathEdge* pe : edges_[drxn] ) pe->edge[drxn]->setBranch( branch, dist + pe->ol, drxn );
}

void NodePath::setBranch( unordered_set<NodePath*>& used, Nodes& extable, Nodes& leapable, bool drxn )
{
    Node* branch = drxn ? path_[0] : path_.back();
    int32_t cutoff = max( 1, params.maxPeMean - branch->size() );
    for ( auto& no : NodeDists( branch, cutoff, drxn, drxn, true ).map ) if ( no.second <= cutoff && no.first->isContinue( drxn ) )
    {
        extable += no.first;
        return;
    }
    int32_t coord = drxn ? ends_[0] + params.maxPeMean : ends_[1] - params.maxPeMean;
    for ( PathEdge* pe : edges_[!drxn] ) if ( pe->edge[!drxn]->doesReach( used, coord, drxn ) ) leapable += branch;
}

void NodePath::setEnds( unordered_set<NodePath*>& pathed, vector<NodePath*>& ends, bool drxn )
{
    if ( !pathed.insert( this ).second ) return;
    
    int limits[2]{ 10, 0 }, i = 0;
    vector<PathEdge*> edges[3];
    if ( verified_ ) for ( PathEdge* pe : edges_[drxn] )
    {
        limits[0] = min( limits[0], pe->score );
        limits[1] = max( limits[1], pe->score );
        if ( pe->edge[drxn]->verified_ && ( i = 2 ) ) edges[2].push_back( pe );
        else edges[ pe->score > 2 ].push_back( pe );
    }
    if ( !i && !edges[1].empty() && ( edges[0].empty() || ( limits[1]+5 >= limits[0]*10 ) ) ) i = 1;
    
    if ( i ) for ( PathEdge* pe : edges[i] ) pe->edge[drxn]->setEnds( pathed, ends, drxn );
    else ends.push_back( this );
}

//bool NodesPath::setEnds( bool drxn )
//{
//    if ( !edges_[drxn].empty() ) return false;
//    fill( drxn ? path_.back() : path_[0], coords_[drxn], drxn ? ends_[1]+500 : ends_[0]-500+params.readLen, drxn );
//    return true;
//}

//void NodesPath::setMates( bool drxn )
//{
//    Lib* lib;
//    for ( Node* node : path_ ) for ( auto& read : node->reads_ ) if ( ( lib = params.getLib( read.first ) ) && !lib->isPe )
//    {
//        ReadId id = read.first;
//        int d;
//        if ( lib->getPair( id, d ) && d == drxn ) mp_[d].insert( make_pair( d ? id : read.first, lib->size ) );
//    }
//}

void NodePath::setMulti()
{
    multi_ = max( getMulti( 0 ), getMulti( 1 ) );
    
    if ( verified_ ) return;
    for ( int d : { 0, 1 } ) for ( PathEdge* pe : edges_[d] ) if ( pe->score > 3 ) return;
    
    multi_ = 0;
}

void NodePath::setPairs( NodePath* np, vector<NodePath*>& path, int32_t diff, int32_t limit, bool lFork, bool rFork )
{
    if ( !path.empty() && np == this && ( !loop_ || -diff < loop_ ) ) loop_ = -diff;
    if ( find( path.begin(), path.end(), np ) != path.end() ) return;
    if ( np != this && !PathPairs::add( this, np, diff, lFork && rFork ) ) return;
    
    if ( !lFork ) for ( PathEdge* pe : np->edges_[1] ) if ( pe->edge[1]->edges_[0].size() > 1 ) lFork = true;
    if ( lFork && np->edges_[1].size() > 1 ) rFork = true;
    
    path.push_back( this );
    if ( np != this ) limit -= np->size();
    if ( limit > 0 ) for ( PathEdge* pe : np->edges_[1] ) setPairs( pe->edge[1], path, diff + pe->diff, limit + pe->ol, lFork, rFork );
    path.pop_back();
}

void NodePath::setPairs( NodePath* np, vector<NodePath*>& path, int32_t diff, bool forked )
{
    if ( ends_[0] - diff - np->ends_[1] > params.maxMpMean ) return;
    if ( np == this && !path.empty() ) addPair( this, diff );
    if ( count( path.begin(), path.end(), this ) ) return;
    path.push_back( this );
    
    bool failed = false;
    if ( forked ) failed = !np->addPair( this, diff );
    
    if ( edges_[1].size() > 1 ) forked = true;
    if ( edges_[1].size() == 1 && edges_[1][0]->edge[1]->edges_[0].size() > 1 ) forked = true;
    if ( !failed ) for ( PathEdge* pe : edges_[1] ) pe->edge[1]->setPairs( np, path, diff + pe->diff, forked );
    path.pop_back();
}

void NodePath::setReachable( unordered_map<NodePath*, vector<int32_t> >& reached, int32_t diff, int32_t limit, bool drxn )
{
    auto it = reached.insert( make_pair( this, vector<int32_t>{ diff } ) );
    if ( !it.second )
    {
        for ( int32_t d : it.first->second ) if ( d == diff ) return;
        it.first->second.push_back( diff );
        sort( it.first->second.begin(), it.first->second.end() );
    }
    
    if ( limit > 0 ) for ( PathEdge* pe : edges_[drxn] )
    {
        int32_t dist = limit - ( drxn ? pe->edge[1]->coords_[1] - pe->diff - coords_[1] : coords_[0] - pe->diff - pe->edge[0]->coords_[0] );
        assert( dist < limit );
        pe->edge[drxn]->setReachable( reached, diff + pe->diff, dist, drxn );
    }
}

void NodePath::setReachable( vector< pair<NodePath*, int32_t> >& reached, unordered_set<NodePath*>& block, int32_t diff, int32_t limit, bool drxn, bool ignore )
{
    if ( !ignore )
    {
        if ( block.find( this ) != block.end() ) return;
        int i = 0;
        while ( i < reached.size() && reached[i].first != this ) i++;
        if ( i == reached.size() ) reached.push_back( make_pair( this, diff ) );
        else if ( diff >= reached[i].second ) return;
        else reached[i].second = diff;
//        auto it = reached.insert( make_pair( this, diff ) );
//        if ( !it.second && diff >= it.first->second ) return;
//        if ( !it.second ) it.first->second = diff;
    }
    
    if ( limit > 0 ) for ( PathEdge* pe : edges_[drxn] )
    {
        int32_t dist = limit - ( drxn ? pe->edge[1]->coords_[1] - pe->diff - coords_[1] : coords_[0] - pe->diff - pe->edge[0]->coords_[0] );
        assert( dist < limit );
        pe->edge[drxn]->setReachable( reached, block, diff + pe->diff, dist, drxn, false );
    }
}

void NodePath::setReachable( unordered_set<NodePath*>& reached, bool drxn )
{
    if ( !reached.insert( this ).second ) return;
    for ( PathEdge* pe : edges_[drxn] ) pe->edge[drxn]->setReachable( reached, drxn );
}

void NodePath::setReachable( unordered_set<NodePath*>& reached, unordered_set<NodePath*>& blocked, int32_t limit, bool drxn, bool ignore )
{
    if ( !ignore && ( blocked.find( this ) != blocked.end() || !reached.insert( this ).second ) ) return;
    
    if ( limit > 0 ) for ( PathEdge* pe : edges_[drxn] )
    {
        int32_t dist = limit - ( drxn ? pe->edge[1]->coords_[1] - pe->diff - coords_[1] : coords_[0] - pe->diff - pe->edge[0]->coords_[0] );
        assert( dist < limit );
        pe->edge[drxn]->setReachable( reached, blocked, dist, drxn, false );
    }
}

void NodePath::setReads( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& paths )
{
    unordered_map<Node*, NodePath*> mapped;
    for ( NodePath* np : paths ) for ( Node* node : np->path_ ) mapped.insert( make_pair( node, np ) );
    for ( NodePath* np : paths ) for ( Node* node : np->path_ ) node->remap( bwt );
    
    vector<PathPair*> reads;
    Lib* lib;
    Coords* coords[2];
    unordered_set<ReadId> used;
    for ( NodePath* np : paths )
    {
        bool forked = false;;
        int32_t minDist = params.maxMpMean * 2 + 200;
        for ( PathPairs* pp : np->paired_[1] ) if ( pp->forked_ && ( forked = true ) ) minDist = min( minDist, pp->node_[1]->ends_[0] - pp->diff_ - np->ends_[1] );
        if ( !forked ) continue;
        
        for ( pair<Node*, int32_t> a : np->offs_ ) for ( auto& read : a.first->reads_ ) if ( lib = params.getLib( read.first ) )
        {
            ReadId id[2]{ read.first, read.first };
            if ( lib->size * 1.2 + 300 < np->ends_[1] - a.second + minDist ) continue;
            if ( lib->getPair( id[1] ) != 1 ) continue;
            if ( used.find( id[0] ) != used.end() ) continue;
            coords[0] = &read.second;
            
            bool found = false, failed = false;
            for ( PathPairs* pp : np->paired_[1] ) if ( pp->forked_ )
            {
                for ( pair<Node*, int32_t> b : pp->node_[1]->offs_ ) if ( ( coords[1] = b.first->getRead( id[1] ) ) && ( found = true ) )
                {
                    assert( used.insert( id[0] ).second );
                    PathPair* p = new PathPair( id[0], lib->size );
                    Node* node[2]{ a.first, b.first };
                    for ( int d : { 0, 1 } ) p->marks_[d].push_back( make_pair( pp->node_[d], ( d ? b : a ).second + (*coords[d])[d] - node[d]->ends_[1] ) );
                    for ( int d : { 0, 1 } ) if ( !p->add( node[d], coords[d]->ignore ? &nodes : node[d]->cloned_, mapped, id[d], d ) ) failed = true;
                    if ( failed || !p->confirm( false ) ) delete p;
                    else reads.push_back( p );
                    break;
                }
                if ( found ) break;
            }
        }
    }
    cleanPairs( paths );
}

void NodePath::setScore()
{
    if ( paired_[1].empty() ) return;
    
    for ( const pair<Node*, int32_t>& no : offs_ ) for ( auto& np : no.first->hits_.pairs[1] )
    {
        int32_t* off, cutoff = np.second.estimate() + np.second.margin();
        for ( PathPairs* pp : paired_[1] ) if ( off = pp->node_[1]->get( np.first ) )
        {
            int32_t dist = *off - no.second - pp->diff_;
            if ( dist < cutoff ) pp->score_ += np.second.count;
        }
    }
}

void NodePath::setScores()
{
    if ( pairs_[1].empty() ) return;
    
    for ( const pair<Node*, int32_t>& no : offs_ )
    {
        for ( auto& np : no.first->hits_.pairs[1] )
        {
            int32_t* off, est = np.second.estimate();
            for ( PathPairing* pp : pairs_[1] ) if ( off = pp->node[1]->get( np.first ) )
            {
                int32_t dist = *off - no.second;
                int32_t cutoff = np.second.margin();
                int32_t best = cutoff;
                vector<int32_t> diffs;
                for ( int32_t diff : pp->diffs ) diffs.push_back( abs( dist - diff - est ) );
                for ( int32_t diff : diffs ) best = min( best, diff );
                if ( best == cutoff ) continue;
                cutoff += best;
                pp->score += np.second.count;
                for ( int i = 0; i < diffs.size(); i++ ) if ( diffs[i] > cutoff ) pp->missed[i] += np.second.count;
            }
        }
    }
    for ( int i = 0; i < pairs_[1].size(); i++ ) if ( !pairs_[1][i]->score ) delete pairs_[1][i--];
}

void NodePath::setTarget( unordered_set<NodePath*>& reached, int32_t dist, int32_t limit, bool drxn, bool ignore )
{
    if ( !ignore && !reached.insert( this ).second ) return;
    if ( !ignore ) dist += ends_[1] - ends_[0];
    if ( dist < limit ) for ( PathEdge* pe : edges_[drxn] ) pe->edge[drxn]->setTarget( reached, dist - pe->ol, limit, drxn );
}

int32_t NodePath::size()
{
    return ends_[1] - ends_[0];
}

PathEdgeScore::PathEdgeScore( NodePath* np, int32_t diff, int32_t limit, bool drxn )
: node( np ), diff( diff ), hits( 0 ), miss( 0 ), unique( 0 )
{
    np->setReachable( fwd, diff, limit, drxn );
}

//bool PathEdgeScore::add( PathPairing* pp, vector< pair<NodesPath*, int32_t> >& tar, vector<int32_t>& diffs, bool missed, bool drxn )
//{
////    auto it = tar.find( pp->node[!drxn] );
////    if ( it == tar.end() ) return false;
////    ( missed ? miss : score ) += pp->hits( diffs, it->second );
//    return true;
//}

void PathEdgeScore::add( NodePath* q, vector<int32_t>& qDiffs, unordered_set<NodePath*>& alts, vector< pair<NodePath*, int32_t> >& tar, bool drxn )
{
    bool used = false;
    for ( auto& nd : tar ) if ( nd.first == q ) used = true;
    
    int scores[2]{0};
    
    for ( PathPairing* pp : q->pairs_[!drxn] )
    {
        vector<int32_t> tDiffs;
        for ( auto& nd : tar ) if ( nd.first == pp->node[!drxn] ) tDiffs.push_back( nd.second );
        
        // This pair is not part of the target
        if ( tDiffs.empty() )
        {
            // Score if this is a valid alt
            if ( alts.find( pp->node[!drxn] ) != alts.end() ) scores[0] += pp->score;
            continue;
        }
        
        int maxScore = 0, minScore = pp->score, sum = 0;
        for ( int32_t td : tDiffs )
        {
            bool scored = false;
            int best = 0;
            for ( int32_t qd : qDiffs ) for ( int i = 0; i < pp->diffs.size(); i++ ) if ( td + qd == pp->diffs[i] )
            {
                best = max( best, pp->score - pp->missed[i] );
                scored = true;
            }
            assert( scored );
            maxScore = max( maxScore, best );
            minScore = min( minScore, best );
            sum += best;
        }
        
        assert( sum <= pp->score );
        assert( minScore <= pp->score );
        scores[1] += used ? min( sum, minScore ) : min( sum, pp->score );
    }
    
    miss += scores[0];
    hits += scores[1];
    if ( !scores[0] || scores[1] > scores[0]*5 ) unique += scores[1] - scores[0];
}

vector<int32_t>* PathEdgeScore::get( NodePath* np )
{
    auto it = fwd.find( np );
    return it != fwd.end() ? &it->second : NULL;
}

void PathEdgeScore::getKeep( unordered_set<NodePath*>& keep, int32_t limit, bool drxn )
{
    for ( auto& nd : fwd ) if ( drxn ? nd.first->ends_[0] - nd.second.back() < limit 
                                     : limit < nd.first->ends_[1] + nd.second.back() ) keep.insert( nd.first );
}

void PathEdgeScore::setKeep( unordered_set<NodePath*>& keep )
{
    for ( auto it = fwd.begin(); it != fwd.end(); )
    {
        if ( it->first != node && keep.find( it->first ) == keep.end() ) it = fwd.erase( it );
        else it++;
    }
}

PathMapping::PathMapping( vector<NodePath*>& path )
: path_( path ), diff_( 0 )
{
    for ( NodePath* np : path_ ) diffs_[1].push_back( make_pair( np, 0 ) );
    for ( NodePath* np : path_ ) diffs_[0].insert( diffs_[0].begin(), make_pair( np, 0 ) );
    for ( int d : { 0, 1 } ) edge( d );
    while ( advance() || crossroad() );
}

void PathMapping::add( NodePath* np, int32_t diff, bool drxn )
{
    diff_ += diff;
    if ( diff ) for ( auto& nd : diffs_[drxn] ) nd.second += diff;
    if ( diff ) for ( auto& nd : alts_[!drxn] ) nd.second += diff;
    for ( PathEdge* pe : np->edges_[!drxn] ) if ( pe->edge[!drxn] != getFork( drxn ) ) alts_[!drxn].push_back( make_pair( pe->edge[!drxn], pe->diff ) );
    path_.insert( drxn ? path_.end() : path_.begin(), np );
    diffs_[drxn].push_back( make_pair( np, 0 ) );
    diffs_[!drxn].insert( diffs_[!drxn].begin(), make_pair( np, diff_ ) );
    exts_[drxn].clear();
    edge( drxn );
}

bool PathMapping::advance()
{
    int scores[2][2];
    for ( int d : { 0, 1 } )
    {
        while ( exts_[d].size() == 1 ) add( exts_[d][0].node, exts_[d][0].diff, d );
        sort( exts_[d].begin(), exts_[d].end(), []( PathEdgeScore& a, PathEdgeScore& b ){ return a.hits+a.unique > b.hits+b.unique; } );
        scores[d][0] = exts_[d].empty() ? 0 : exts_[d][0].hits + exts_[d][0].unique;
        scores[d][1] = exts_[d].size() < 2 ? 0 : exts_[d][1].hits + exts_[d][1].unique;
    }
    
    if ( !scores[0][0] && !scores[1][0] ) return false;
    int d = scores[1][0] - scores[1][1] > scores[0][0] - scores[0][1];
    assert( !exts_[d].empty() );
    for ( int i = 1; i < exts_[d].size(); i++ ) alts_[d].push_back( make_pair( exts_[d][i].node, diff_ + exts_[d][i].diff ) );
    add( exts_[d][0].node, exts_[d][0].diff, d );
    return true;
}

bool PathMapping::crossroad()
{
    bool pathable[2]{ false, false };
    for ( int d : { 0, 1 } ) for ( PathEdgeScore& pes : exts_[d] ) if ( pes.node->verified_ ) pathable[d] = true;
    if ( !pathable[0] || pathable[1] ) return false;
    assert( false );
    return true;
}

void PathMapping::edge( bool drxn )
{
    if ( !getFork( drxn )->verified_ ) return;
    
    // instantiate forward edges with nodes up 1000 distance away
    for ( PathEdge* pe : getFork( drxn )->edges_[drxn] ) exts_[drxn].push_back( PathEdgeScore( pe->edge[drxn], pe->diff, 1000, drxn ) );
    
    unordered_set<NodePath*> base, keep, alts;
    
    // Ensure that foward nodes are not unevenly included due minor differences in distance
    for ( PathEdgeScore& pes : exts_[drxn] ) pes.getKeep( keep, getFork( drxn )->ends_[drxn] + ( drxn ? 400 : -400 ), drxn );
    for ( NodePath* np : path_ ) keep.erase( np );
    for ( PathEdgeScore& pes : exts_[drxn] ) pes.setKeep( keep );
    
    // Block for the alts
    for ( NodePath* np : path_ ) base.insert( np );
    for ( int d : { 0, 1 } ) for ( PathEdge* pe : getFork( d )->edges_[d] ) pe->edge[d]->setReachable( base, d );
    
    for ( auto& nd : alts_[!drxn] ) nd.first->setReachable( alts, base, 500, !drxn );
    
    vector<int32_t>* diffs;
    for ( int i = path_.size(); i-- > 0; )
    {
        NodePath* np = diffs_[drxn][i].first;
        int32_t diff = diffs_[drxn][i].second;
        int32_t dist = np->getLen( diffs_[drxn].back().first, diff, drxn );
        for ( PathPairing* pp : np->pairs_[drxn] ) for ( PathEdgeScore& pes : exts_[drxn] ) if ( diffs = pes.get( pp->node[drxn] ) )
        {
            int best = 0;
            for ( int32_t d : *diffs ) for ( int j = 0; j < pp->diffs.size(); j++ ) if ( d + diff == pp->diffs[j] ) best = max( best, pp->score - pp->missed[j] );
            pes.hits += best;
        }
        if ( i > 1 && exts_[drxn].size() > 1 && np->edges_[!drxn].size() > 1 )
        {
            int x = 0;
        }
    }
    int x = 0;
//    for ( PathEdgeScore& pes : edges_[drxn] )
//    {
//        vector< pair<NodesPath*, int32_t> > rev = alts;
//        pes.node->setReachable( rev, base, pes.diff, !drxn, false );
//        
//        for ( auto& nd : pes.fwd ) for ( PathPairing* pp : nd.first->pairs_[!drxn] )
//        {
//            if ( pes.add( pp, diffs_[drxn], nd.second, false, drxn ) ) continue;
//            if ( pes.add( pp, rev, nd.second, true, drxn ) ) continue;
//        }
//    }
}

NodePath* PathMapping::getFork( bool drxn )
{
    assert( !path_.empty() );
    return drxn ? path_.back() : path_[0];
}

vector<NodePath*> PathMapping::map( NodePath* seed )
{
    vector<NodePath*> path = { seed };
    PathMapping pm( path );
    return path;
}

//void PathScores::add( Node* node )
//{
//    float coef = 1;
//    added += node;
//    if ( node->drxn_ < 2 )
//    {
//        float bad[2]{ 0, 0 };
//        for ( auto& np : node->hits_.pairs[!node->drxn_] ) bad[np.first->bad_] += np.second.count;
//        if ( bad[1] >= bad[0] ) return;
//        coef = ( bad[0] - bad[1] ) / ( bad[0] + bad[1] );
//    }
//    for ( int d : { 0, 1 } )
//    {
//        for ( auto& np : node->hits_.pairs[d] ) add( np.first, np.second.count * coef, d );
//        add( node, node->hits_.count * coef, d );
//    }
//}
//
//void PathScores::add( Node* node, float score, bool drxn )
//{
//    if ( score <= 0 ) return;
//    auto it = scores[drxn].insert( make_pair( node, score ) );
//    if ( !it.second ) it.first->second += score;
//}
//
//float PathScores::get( Node* node, int drxn )
//{
//    float score = 0;
//    for ( int d = drxn == 1; d < min( drxn+1, 2 ); d++ )
//    {
//        auto it = scores[drxn].find( node );
//        if ( it != scores[drxn].end() ) score += it->second;
//    }
//    return score;
//}
//
//float PathScores::get( Nodesx& nodes, int drxn )
//{
//    float score = 0;
//    for ( Node* node : nodes.nodes ) score += get( node, drxn );
//    return score;
//}
//
//PathScores PathScores::unadded()
//{
//    PathScores r;
//    for ( int d : { 0, 1 } ) for ( auto& np : scores[d] ) if ( !added.find( np.first ) ) r.scores[d].insert( np );
//    return r;
//}
//
//NodeBranch::NodeBranch()
//: fork( NULL ), branch( NULL ), score( 0 )
//{
//    
//}
//
//NodeBranch::NodeBranch( Node* fork, Node* branch )
//: fork( fork ), branch( branch ), score( 0 )
//{
//    
//}
//
//NodeBranch::NodeBranch( Node* fork, Node* branch, Nodesx& fwdIn, NodeScores& scores, bool drxn )
//: fork( fork ), branch( branch ), fwd( Nodesx::inSet( branch, fwdIn, drxn, true ) ), score( 0 )
//{
//    for ( Node* f : fwd.nodes ) score += scores.get( f );
//}
//
//NodeBranch::NodeBranch( Node* fork, Node* branch, Nodesx& fwdIn, PathScores& scores, bool drxn )
//: fork( fork ), branch( branch ), fwd( Nodesx::inSet( branch, fwdIn, drxn, true ) ), score( 0 )
//{
//    score = scores.get( fwd, drxn );
//}
//
//vector<NodeBranch> NodeBranch::get( Node* fork, Nodesx fwd, NodeScores &scores, bool drxn )
//{
//    vector<NodeBranch> branches;
//    if ( !fork ) return branches;
//    for ( Edge &e : fork->edges_[drxn] ) if ( fwd.find( e.node ) ) branches.push_back( NodeBranch( fork, e.node, fwd, scores, drxn ) );
//    sort( branches );
//    return branches;
//}
//
//vector<NodeBranch> NodeBranch::get( Node* fork, Nodesx& fwd, PathScores& scores, bool drxn )
//{
//    vector<NodeBranch> branches;
//    for ( Edge &e : fork->edges_[drxn] ) if ( fwd.find( e.node ) ) branches.push_back( NodeBranch( fork, e.node, fwd, scores, drxn ) );
//    sort( branches );
//    return branches;
//}
//
//void NodeBranch::sort( vector<NodeBranch> &branches )
//{
//    std::sort( branches.begin(), branches.end(), []( NodeBranch &a, NodeBranch &b ){ return a.score > b.score; } );
//}
//
//PathChunk::PathChunk( NodeList &pilot, Node* alt, int &i )
//{
//    hits[0] = hits[1] = 0;
//    edges[0][0] = edges[0][1] = edges[1][0] = edges[1][1] = false;
//    path[0].insert( path[0].end(), pilot.begin() + i, pilot.end() );
//    pilot.erase( pilot.begin() + i, pilot.end() );
//    i = 0;
//    if ( alt ) path[1].push_back( alt );
//    else i = 1;
//}
//
//void PathChunk::add( vector<PathChunk> &chunks, Node* alt, int &i )
//{
//    for ( int j = 0; j < path[0].size(); j++ )
//    {
//        if ( path[0][j] != alt ) continue;
//        i = j;
//        chunks.push_back( PathChunk( path[0], NULL, i ) );
//        return;
//    }
//    path[1].push_back( alt );
//}
//
//NodeBranch PathChunk::alt( NodeScores& scores, Nodesx& added, Nodesx& btw, int& cumul, bool drxn )
//{
//    NodeBranch best;
//    best.score = 9;
//    for ( int i : { 0, 1 } )
//    {
//        int miss = 0;
//        for ( Node* node : path[i] )
//        {
//            Nodesx fwd = Nodesx::notSet( node, added, params.maxPeMean, drxn, false );
//            fwd -= btw;
//            for ( Edge& e : node->edges_[drxn] )
//            {
//                if ( !fwd.find( e.node ) ) continue;
//                NodeBranch alt( node, e.node, fwd, scores, drxn );
//                alt.score -= cumul + min( miss, hits[!i] );
//                if ( alt.score > best.score )
//                {
//                    best = alt;
//                }
//            }
//            miss += scores.get( node );
//        }
//    }
//    
//    cumul += min( hits[0], hits[1] );
//    return best;
//}
//
//Nodesx PathChunk::between( vector<PathChunk> &chunks )
//{
//    Nodesx nodes[3];
//    
//    for ( PathChunk& chunk : chunks )
//    {
//        for ( int i = 0; i < 2 && !chunk.path[i].empty(); i++ )
//        {
//            nodes[0].fill( chunk.path[i][0], 1, true, false );
//            nodes[1].fill( chunk.path[i].back(), 0, true, false );
//        }
//    }
//    for ( Node* node : nodes[0].nodes ) if ( nodes[1].find( node ) ) nodes[2].add( node );
//    
//    return nodes[2];
//}
//
////void PathChunk::branch( NodeScores& scores, Nodes& added, Nodes& btw, bool drxn )
////{
////    branches[drxn].clear();
////    for ( int i : { 0, 1 } )
////    {
////        for ( Node* node : path[i] ) NodeBranch::add( branches[drxn], node, added, btw, scores, params.maxPeMean, drxn );
////    }
////    NodeBranch::sort( branches[drxn] );
////}
//
////vector<PathChunk> PathChunk::fill( NodeList &pilot, NodeScores &scores, Nodes &added )
////{
////    Nodes allNodes = Nodes::between( pilot[0], pilot.back(), true ), altNodes;
////    for ( Node* node : allNodes.nodes ) if ( !added.find( node ) ) altNodes.add( node );
////    NodeScores altScores = scores.split( altNodes );
////    
////    Node* cur = pilot[0];
////    int i = 0;
////    vector<PathChunk> chunks( 1, PathChunk( pilot, NULL, i ) );
////    while ( cur )
////    {
////        vector<NodeBranch> branches = NodeBranch::get( cur, allNodes, altScores, 1 );
////        if ( branches.empty() ) break;
////        cur = branches[0].branch;
////        if ( !added.find( cur ) ) scores.add2( cur, true );
////        added.add( cur );
////        for ( auto &np : cur->pairs_ ) if ( altNodes.find( np.first ) ) altScores.add2( np );
////        if ( i )
////        {
////            if ( cur == chunks.back().path[0][i] ) i++;
////            else chunks.push_back( PathChunk( chunks.back().path[0], cur, i ) );
////        }
////        else chunks.back().add( chunks, cur, i );
////    }
////    for ( int i = 0; i < chunks.size(); i++ )
////    {
////        chunks[i].edges[0][0] = chunks[i].edges[0][1] = i > 0;
////        chunks[i].edges[1][0] = chunks[i].edges[1][1] = i < chunks.size()-1;
////        if ( chunks[i].path[1].empty() ) chunks[i].edges[0][1] = chunks[i].edges[1][1] = false;
////        for ( int j : { 0, 1 } ) for ( Node* node : chunks[i].path[j] ) chunks[i].hits[j] += scores.get( node );
////    }
////    
////    return chunks;
////}
//
//bool PathChunk::setFork( Node* branches[2], bool drxn )
//{
//    for ( int i = 0; i < 2 && !path[i].empty(); i++ )
//    {
//        Nodesx nots;
//        for ( int j = 0; j < path[i].size() && !branches[i]; j++ )
//        {
//            Node* b = drxn ? path[i].end()[-j-1] : path[i][j];
//            for ( Node* f : Nodesx::notSet( b, nots, drxn, true ).nodes )
//            {
//                if ( f->isContinue( drxn ) ) branches[i] = b;
//                nots.add( f );
//            }
//        }
//    }
//    return branches[0] || branches[1];
//}
//
//bool PathChunk::take( Node* fork, PathChunk& tar, Nodesx &ignore, bool drxn )
//{
//    edges[!drxn][0] = true;
//    for ( int i : { 0, 1 } )
//    {
//        auto it = find( tar.path[i].begin(), tar.path[i].end(), fork );
//        if ( it == tar.path[i].end() ) continue;
//        int j = tar.path[!i].empty();
//        if ( !j ) path[1].insert( drxn ? path[1].begin() : path[1].end(), tar.path[!i].begin(), tar.path[!i].end() );
//        path[j].insert( drxn ? path[j].begin() : path[j].end(), drxn ? it+1 : tar.path[i].begin(), drxn ? tar.path[i].end() : it );
//        tar.path[i].erase( drxn ? it+1 : tar.path[i].begin(), drxn ? tar.path[i].end() : it );
//        edges[!drxn][1] = true;
//        tar.edges[drxn][i] = true;
//        if ( !j ) tar.edges[drxn][!i] = true;
//        return !j;
//    }
//    
//    int j = tar.hits[1] > tar.hits[0];
//    path[1].insert( drxn ? path[1].begin() : path[1].end(), tar.path[j].begin(), tar.path[j].end() );
//    
//    return true;
//}
//
//NodePath::NodePath( Node* seed, PathScores& scores )
//{
//    path_[0].push_back( seed );
//    scores.add( seed );
//    
//    for ( int again = 1; again-- > 0; )
//    {
//        while ( NodePath::path( path_[0], scores, 0 ) );
//        while ( NodePath::path( path_[0], scores, 1 ) );
//        again = unpause( path_[0], scores );
//    }
//}
//
//NodePath::NodePath( vector<Node*> a, vector<Node*> b )
//{
//    path_[0] = a;
//    path_[1] = b;
//}
//
//void NodePath::branch( vector<NodePath>& path, PathScores& scores, bool drxn )
//{
//    if ( ( drxn ? path.back().path_[0].back() : path[0].path_[0][0] )->cloned_ ) return;
//    Nodesx btw = Nodesx::between( path[0].path_[0][0], path.back().path_[0].back(), true );
//    int cumul = 0;
//    PathScores altScores = scores.unadded();
//    NodeBranch best;
//    for ( int i = 0; i < path.size() && cumul < 20; i++ ) ( drxn ? path.end()[-i-1] : path[i] ).branch( best, scores, altScores, btw, cumul, drxn );
//    if ( !best.branch ) return;
//    vector<Node*> ext;
//    while ( path.size() > 1 && !( drxn ? path.back() : path[0] ).branch( best, path, scores, ext, drxn ) );
//}
//
//void NodePath::branch( NodeBranch& best, PathScores& scores, PathScores& altScores, Nodesx& btw, int& cumul, bool drxn )
//{
//    float miss[2]{ (float)cumul, (float)cumul };
//    for ( int i = 0; i < 2; i++ )
//    {
//        for ( int j = 0; j < path_[i].size(); j++ )
//        {
//            Node* node = drxn ? path_[i].end()[-j-1] : path_[i][j];
//            for ( Edge& e : node->edges_[drxn] )
//            {
//                if ( scores.added.find( e.node ) ) continue;
//                NodeBranch nb( node, e.node );
//                nb.fwd = Nodesx::notSet( e.node, scores.added, drxn, true );
//                nb.fwd -= btw;
//                nb.score = scores.get( nb.fwd, drxn ) - ( path_[!i].empty() ? 0 : miss[i] );
//                if ( nb.score >= 10 && nb.score > best.score ) best = nb;
//            }
//            
//            miss[i] += node->hits_.get();
//        }
//    }
//    cumul = min( miss[0], miss[1] );
//}
//
//bool NodePath::branch( NodeBranch& best, vector<NodePath>& path, PathScores& scores, vector<Node*>& ext, bool drxn )
//{
//    for ( int i = 0; i < 2; i++ )
//    {
//        auto it = find( path_[i].begin(), path_[i].end(), best.fork );
//        if ( it == path_[i].end() ) continue;
//        scores.add( best.branch );
//        best.branch->offset( best.fork, drxn );
//        if ( path_[!i].empty() )
//        {
//            ext.insert( drxn ? ext.end() : ext.begin(), drxn ? it+1 : path_[i].begin(), drxn ? path_[i].end() : it );
//            path_[i].erase( drxn ? it+1 : path_[i].begin(), drxn ? path_[i].end() : it );
//            NodePath np( ext, vector<Node*>( 1, best.branch ) );
//            while( NodePath::path( np.path_[1], scores, drxn ) );
//            path.insert( drxn ? path.end() : path.begin(), np );
//        }
//        else
//        {
//            path_[i].erase( drxn ? it+1 : path_[i].begin(), drxn ? path_[i].end() : it );
//            path_[i].insert( drxn ? path_[i].end() : path_[i].begin(), best.branch );
//            path_[!i].insert( drxn ? path_[!i].end() : path_[!i].begin(), ext.begin(), ext.end() );
//            if ( !i ) swap( path_[0], path_[1] );
//            while( NodePath::path( path_[1], scores, drxn ) );
//            assert( false );
//        }
//        return true;
//    }
//    
//    ext.insert( drxn ? ext.begin() : ext.end(), path_[0].begin(), path_[0].end() );
//    path.erase( drxn ? path.end()-1 : path.begin() );
//    return false;
//}
//
//vector< vector<NodePath> > NodePath::create( NodeRoll& nodes )
//{
//    vector< vector<NodePath> > paths;
//    vector< pair<Node*, int> > seeds;
//    for ( Node* node : nodes.getGraph( 2 ).nodes )
//    {
//        seeds.push_back( make_pair( node, 0 ) );
//        Nodesx tars[2]{ Nodesx( node, params.maxPeMean, 0, true, false ), Nodesx( node, params.maxPeMean, 1, true, false ) };
//        for ( Node* t : tars[0].nodes ) seeds.back().second += t->hits_.get( tars[1], 1 );
//    }
//    sort( seeds.begin(), seeds.end(), []( pair<Node*, int> &a, pair<Node*, int> &b ){ return a.second > b.second; } );
//    
//    Nodesx used;
//    for ( int i = 0; i < seeds.size(); i++ )
//    {
//        if ( used.find( seeds[i].first ) ) continue;
//        PathScores scores;
//        vector<NodePath> path( 1, NodePath( seeds[i].first, scores ) );
//        while( path[0].fill( path, scores, 0 ) );
//        while( path.back().fill( path, scores, 1 ) );
//        NodePath::branch( path, scores, 0 );
//        NodePath::branch( path, scores, 1 );
//        paths.push_back( path );
//        break;
//    }
//    return paths;
//}
//
//bool NodePath::fill( vector<NodePath>& path, PathScores& scores, bool drxn )
//{
//    assert( !path_[0].empty() && path_[1].empty() );
//    Node* cur = NULL;
//    if ( !drxn ) for ( int i = path_[0].size(); !cur && --i >= 0; ) if ( path_[0][i]->drxn_ != 1 ) cur = path_[0][i];
//    if ( drxn ) for ( int i = 0; !cur && i < path_[0].size(); i++ ) if ( path_[0][i]->drxn_ != 0 ) cur = path_[0][i];
//    
//    PathScores altScores = scores.unadded();
//    Nodesx btw = Nodesx::between( path_[0][0], path_[0].back(), true );
//    vector<Node*> alt, parts[3];
//    Node* forks[2]{ NULL, NULL };
//    while ( cur && ( !forks[0] || !forks[1] ) )
//    {
//        vector<NodeBranch> branches = NodeBranch::get( cur, btw, altScores, drxn );
//        if ( branches.empty() || ( alt.empty() && !branches[0].score ) ) return false;
//        if ( scores.added.find( branches[0].branch ) )
//        {
//            if ( !alt.empty() ) forks[drxn] = branches[0].branch;
//        }
//        else
//        {
//            if ( alt.empty() ) forks[!drxn] = cur;
//            branches[0].branch->offset( cur, drxn );
//            alt.insert( drxn ? alt.end() : alt.begin(), branches[0].branch );
//        }
//        cur = branches[0].branch;
//    }
//    assert( forks[0] && forks[1] );
//    for ( Node* node : alt ) scores.add( node );
//    for ( Node* node : path_[0] )
//    {
//        if ( node == forks[1] ) forks[1] = NULL;
//        parts[ ( forks[0] ? 0 : ( forks[1] ? 1 : 2 ) ) ].push_back( node );
//        if ( node == forks[0] ) forks[0] = NULL;
//    }
//    path_[0] = parts[int(!drxn)*2];
//    path.insert( drxn ? path.end() : path.begin(), NodePath( parts[1], alt ) );
//    path.insert( drxn ? path.end() : path.begin(), NodePath( parts[int(drxn)*2], vector<Node*>() ) );
//    return true;
//}
//
//bool NodePath::path( NodeList& path, PathScores& scores, bool drxn )
//{
//    assert( !path.empty() );
//    Node* cur = drxn ? path.back() : path[0];
//    Nodesx fwd = Nodesx::notSet( cur, scores.added, params.maxPeMean, drxn, false );
//    vector<NodeBranch> branches = NodeBranch::get( cur, fwd, scores, drxn );
//    if ( branches.empty() || !branches[0].score ) return false;
//    branches[0].branch->offset( drxn ? path.back() : path[0] , drxn );
//    if ( drxn ) assert( branches[0].branch->ends_[0] < path.back()->ends_[1] && path.back()->ends_[1] < branches[0].branch->ends_[1] );
//    if ( !drxn ) assert( path[0]->ends_[0] < branches[0].branch->ends_[1] && branches[0].branch->ends_[0] < path[0]->ends_[0] );
//    path.insert( drxn ? path.end() : path.begin(), branches[0].branch );
//    scores.add( branches[0].branch );
//    return branches[0].branch->verified_;
//}
//
//bool NodePath::unpause( NodeList& path, PathScores& scores )
//{
//    if ( !path[0]->verified_ && !path.back()->verified_ ) return false;
//    assert( false ); 
//    return true;
//}
