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

//PathPairing::PathPairing( NodePath* l, NodePath* r, int32_t diff )
//: diffs{ diff }, missed{ 0 }, score( 0 )
//{
//    node[0] = l;
//    node[1] = r;
//    for ( int d : { 0, 1 } ) node[d]->pairs_[!d].push_back( this );
//}
//
//PathPairing::~PathPairing()
//{
//    for ( int d : { 0, 1 } ) node[d]->pairs_[!d].erase( remove( node[d]->pairs_[!d].begin(), node[d]->pairs_[!d].end(), this ), node[d]->pairs_[!d].end() );
//}

PathEdge::PathEdge( NodePath* fork, NodePath* branch, Edge& e, int32_t diff, bool drxn )
: diff( diff ), ol( e.ol ), score( -1 ), multi( 1 ), leap( e.leap ), bad( false )
{
    edge[!drxn] = fork;
    edge[drxn] = branch;
    node[!drxn] = fork->getTerminus( drxn );
    node[drxn] = branch->getTerminus( !drxn );
    for ( int d : { 0, 1 } ) edge[d]->edges_[!d].push_back( this );
}

PathEdge::PathEdge( NodePath* l, NodePath* r, PathEdge* pe )
: diff( pe->diff ), ol( pe->ol ), score( pe->score ), multi( 1 ), leap( pe->leap ), bad( false )
{
    edge[0] = l;
    edge[1] = r;
    for ( int d : { 0, 1 } ) edge[d]->edges_[!d].push_back( this );
    edge[0]->path_.back()->addEdge( edge[1]->path_[0], ol, 1, false, leap );
    for ( int d : { 0, 1 } ) node[d] = edge[d]->getTerminus( !d );
}

PathEdge::PathEdge( NodePath* np, Node* ne, PathEdge* pe, bool drxn )
: diff( pe->diff ), ol( pe->ol ), score( pe->score ), leap( pe->leap ), bad( pe->bad )
{
    edge[drxn] = np;
    edge[!drxn] = pe->edge[!drxn];
    node[!drxn] = pe->node[!drxn];
    node[drxn] = ne;
    node[!drxn]->addEdge( node[drxn], ol, drxn, false, leap );
}

PathEdge::~PathEdge()
{
    for ( int d : { 0, 1 } ) edge[d]->edges_[!d].erase( remove( edge[d]->edges_[!d].begin(), edge[d]->edges_[!d].end(), this ), edge[d]->edges_[!d].end() );
    for ( int d : { 0, 1 } ) edge[d]->alts_[!d].erase( remove( edge[d]->alts_[!d].begin(), edge[d]->alts_[!d].end(), this ), edge[d]->alts_[!d].end() );
}

void PathEdge::claim( NodePath* np, bool drxn )
{
    edge[0]->path_.back()->removeEdge( edge[1]->path_[0], 1, true );
    edge[drxn]->edges_[!drxn].erase( remove( edge[drxn]->edges_[!drxn].begin(), edge[drxn]->edges_[!drxn].end(), this ), edge[drxn]->edges_[!drxn].end() );
    edge[drxn] = np;
    np->edges_[!drxn].push_back( this );
    edge[0]->path_.back()->addEdge( edge[1]->path_[0], ol, 1, false, leap );
}

void PathEdge::downgrade( int drxn )
{
    assert( drxn < 2 );
    for ( int d : { 0, 1 } ) if ( d != drxn )
    {
        edge[d]->edges_[!d].erase( remove( edge[d]->edges_[!d].begin(), edge[d]->edges_[!d].end(), this ), edge[d]->edges_[!d].end() );
        edge[d]->alts_[!d].push_back( this );
    }
    bad = true;
}

bool PathEdge::isAlt()
{
    for ( int d : { 0, 1 } ) for ( PathEdge* pe : edge[d]->edges_[!d] ) if ( pe == this ) return false;
    assert( false );
    return true;
}

bool PathEdge::sever( NodePath* l, NodePath* r )
{
    for ( PathEdge* pe : l->edges_[1] ) if ( pe->edge[1] == r )
    {
        l->path_.back()->removeEdge( r->path_[0], 1, true );
        delete pe;
        return true;
    }
    for ( PathEdge* pe : l->alts_[1] ) if ( pe->edge[1] == r )
    {
        l->path_.back()->removeEdge( r->path_[0], 1, true );
        delete pe;
        return true;
    }
    assert( false );
    return false;
}

bool PathEdge::sever()
{
    return node[0]->removeEdge( node[1], 1, true );
}

void PathEdge::sort( vector<PathEdge*>& edges )
{
    std::sort( edges.begin(), edges.end(), []( PathEdge* a, PathEdge* b ){ return a->score > b->score; } );
    assert( edges.size() < 2 || edges[1]->score <= edges[0]->score );
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
: diff_( pp->diff_ ), score_( pp->score_ ), forked_( pp->forked_ ), ids_( pp->ids_.begin(), pp->ids_.end() )
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
    ids_.insert( pp->ids_.begin(), pp->ids_.end() );
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

bool PathPairs::steal( NodePath* np, bool drxn )
{
    assert( shared_.empty() && uniques_.empty() );
    
    if ( node_[drxn] == np )
    {
        np->score_ += ids_.size();
        return false;
    }
    
    for ( PathPairs* alt : np->paired_[drxn] ) if ( alt->node_[drxn] == node_[drxn] )
    {
        alt->ids_.insert( ids_.begin(), ids_.end() );
        alt->score_ += score_;
        return false;
    }
    
    node_[!drxn] = np;
    np->paired_[drxn].push_back( this );
    return true;
}

void PathPairs::test()
{
    for ( int d : { 0, 1 } ) assert( find( node_[d]->paired_[!d].begin(), node_[d]->paired_[!d].end(), this ) != node_[d]->paired_[!d].end() );
}

Branch::Branch( PathEdge* pe, bool drxn )
: edge_( pe ), match_( NULL ), base_( 0 ), hits_( 0 ), miss_( 0 ), used_( false )
{
    pe->edge[drxn]->setBranch( branch_, params.maxPeMean*2, drxn );
}

Branch::~Branch()
{
    for ( pair<Branch*, int>& match : matched_ ) match.first->cull( this );
}

vector<Branch*> Branch::create( NodePath* fork, bool drxn )
{
    vector<Branch*> branches;
    for ( PathEdge* pe : fork->edges_[drxn] ) branches.push_back( new Branch( pe, drxn ) );
    
    for ( int i = 0; i+1 < branches.size(); i++ ) for ( auto it = branches[i]->branch_.begin(); it != branches[i]->branch_.end(); )
    {
        bool bad = true;
        for ( int j = i+1; bad && j < branches.size(); j++ ) if ( branches[j]->branch_.find( *it ) == branches[j]->branch_.end() ) bad = false;
        if ( bad ) for ( int j = i+1;  j < branches.size(); j++ ) branches[j]->branch_.erase( *it );
        if ( bad ) it = branches[i]->branch_.erase( it );
        else it++;
    }
    
    return branches;
}

bool Branch::cull()
{
    bool culled = false;
    for ( int i = 1; i < matched_.size(); i++ )
    {
        int hits[3]{ matched_[i].second, matched_[0].second, matched_[i].first->matched_[0].second };
        if ( hits[0]-1 >= hits[1] / 5 ||  hits[0]-1 >= hits[2] / 5 || hits[0] * 10 > hits[1] + hits[2] ) continue;
        assert( culled = matched_[i].first->cull( this ) );
        matched_.erase( matched_.begin() + i-- );
    }
    return culled;
}

bool Branch::cull( Branch* b )
{
    bool culled = false;
    for ( int i = 0; i < matched_.size(); i++ ) if ( matched_[i].first == b && ( culled = true ) ) matched_.erase( matched_.begin() + i-- );
    return culled;
}

int Branch::getMax( bool inclBase )
{
    return ( inclBase ? base_ : 0 ) + ( matched_.empty() ? 0 : matched_[0].second );
}

bool Branch::match( vector<NodePath*>& path, vector<Branch*> branches[2] )
{
    vector< unordered_set<ReadId> > used[2]{ vector< unordered_set<ReadId> >( branches[0].size() ), vector< unordered_set<ReadId> >( branches[1].size() ) };
    vector< vector< unordered_set<ReadId> > > hits( branches[0].size(), vector< unordered_set<ReadId> >( branches[1].size() ) );
    for ( int d : { 0, 1 } ) setBase( branches[d], path, d );
    for ( int d : { 0, 1 } ) for ( Branch* b : branches[d] ) b->matched_.clear();
    
    int matched = 0;
    for ( int i = 0; i < branches[0].size(); i++ ) for ( NodePath* np : branches[0][i]->branch_ ) for ( PathPairs* pp : np->paired_[1] )
    {
        for ( int j = 0; j < branches[1].size(); j++ ) if ( branches[1][j]->branch_.find( pp->node_[1] ) != branches[1][j]->branch_.end() )
        {
            used[0][i].insert( pp->ids_.begin(), pp->ids_.end() );
            used[1][j].insert( pp->ids_.begin(), pp->ids_.end() );
            hits[i][j].insert( pp->ids_.begin(), pp->ids_.end() );
            matched += pp->ids_.size();
        }
    }
    for ( int d : { 0, 1 } ) for ( int i = 0; i < branches[d].size(); i++ ) branches[d][i]->hits_ = used[d][i].size();
    if ( !matched ) return false;
    
    for ( int d : { 0, 1 } ) if ( used[d].size() > 1 ) for ( auto it = used[d][0].begin(); it != used[d][0].end(); )
    {
        bool bad = true;
        for ( int i = 1; bad && i < used[d].size(); i++ ) if ( used[d][i].find( *it ) == used[d][i].end() ) bad = false;
        if ( bad ) for ( int i = 0; i < used[0].size(); i++ ) for ( int j = 0; j < used[1].size(); j++ ) hits[i][j].erase( *it );
        if ( bad ) for ( int i = 1; bad && i < used[d].size(); i++ ) used[d][i].erase( *it );
        if ( bad ) it = used[d][0].erase( it );
        else it++;
    }
    
    matched = 0;
    for ( int i = 0; i < branches[0].size(); i++ ) for ( int j = 0; j < branches[1].size(); j++ ) if ( !hits[i][j].empty() )
    {
        branches[0][i]->matched_.push_back( make_pair( branches[1][j], hits[i][j].size() ) );
        branches[1][j]->matched_.push_back( make_pair( branches[0][i], hits[i][j].size() ) );
        branches[0][i]->sort();
        branches[1][j]->sort();
        matched += hits[i][j].size();
    }
    if ( !matched ) return false;
    for ( int d : { 0, 1 } ) sort( branches[d] );
    return true;
}

bool Branch::match( vector<NodePath*>& path, vector<Branch*>& l, vector<Branch*>& r )
{
    vector< unordered_set<ReadId> > used[2]{ vector< unordered_set<ReadId> >( l.size() ), vector< unordered_set<ReadId> >( r.size() ) };
    vector< vector< unordered_set<ReadId> > > hits( l.size(), vector< unordered_set<ReadId> >( r.size() ) );
    setBase( l, path, 0 );
    setBase( r, path, 1 );
    for ( Branch* b : l ) b->hits_ = b->miss_ = 0;
    for ( Branch* b : r ) b->hits_ = b->miss_ = 0;
    for ( Branch* b : l ) b->match_ = NULL;
    for ( Branch* b : r ) b->match_ = NULL;
    
    
    int matched = 0;
    for ( int i = 0; i < l.size(); i++ ) for ( NodePath* np : l[i]->branch_ ) for ( PathPairs* pp : np->paired_[1] )
    {
        for ( int j = 0; j < r.size(); j++ ) if ( r[j]->branch_.find( pp->node_[1] ) != r[j]->branch_.end() )
        {
            used[0][i].insert( pp->ids_.begin(), pp->ids_.end() );
            used[1][j].insert( pp->ids_.begin(), pp->ids_.end() );
            hits[i][j].insert( pp->ids_.begin(), pp->ids_.end() );
            matched += pp->ids_.size();
        }
    }
    if ( !matched ) return false;
    for ( int d : { 0, 1 } ) for ( auto it = used[d][0].begin(); it != used[d][0].end(); )
    {
        bool bad = true;
        for ( int i = 1; bad && i < used[d].size(); i++ ) if ( used[d][i].find( *it ) == used[d][i].end() ) bad = false;
        if ( bad ) for ( int i = 0; i < l.size(); i++ ) for ( int j = 0; j < r.size(); j++ ) hits[i][j].erase( *it );
        if ( bad ) for ( int i = 1; bad && i < used[d].size(); i++ ) used[d][i].erase( *it );
        if ( bad ) it = used[d][0].erase( it );
        else it++;
    }
    for ( int i = 0; i < l.size(); i++ ) for ( int j = 0; j < r.size(); j++ ) if ( !hits[i][j].empty() )
    {
        bool best[2]{ hits[i][j].size() > l[i]->hits_, hits[i][j].size() > r[j]->hits_ };
        if ( best[0] ) l[i]->match_ = r[j];
        if ( best[1] ) r[j]->match_ = l[i];
        l[i]->miss_ += best[0] ? l[i]->hits_ : hits[i][j].size();
        r[j]->miss_ += best[1] ? r[j]->hits_ : hits[i][j].size();
        if ( best[0] ) l[i]->hits_ = hits[i][j].size();
        if ( best[1] ) r[j]->hits_ = hits[i][j].size();
    }
    sort( l );
    sort( r );
    return true;
}

bool Branch::matched( vector<PathEdge*> taken[2][2], bool drxn )
{
    if ( matched_.size() != 1 ) return false;
    for ( int d : { 0, 1 } ) for ( int i : { 0, 1 } ) taken[d][i].clear();
    taken[drxn][1].push_back( edge_ );
    taken[!drxn][ matched_[0].first->matched_.size() == 1 ].push_back( matched_[0].first->edge_ );
    return true;
}

void Branch::setBase( vector<Branch*>& branches, vector<NodePath*>& path, bool drxn )
{
    for ( Branch* b : branches ) b->base_ = 0;
    for ( NodePath* np : path ) for ( PathPairs* pp : np->paired_[drxn] ) for ( Branch* b : branches )
    {
        if ( b->branch_.find( pp->node_[drxn] ) != b->branch_.end() ) b->base_ += pp->ids_.size();
    }
    std::sort( branches.begin(), branches.end(), []( Branch* a, Branch* b ){ return a->base_ > b->base_; } );
}

void Branch::setUsed( vector<Branch*>& branches, vector< vector<NodePath*> >& loci, bool drxn )
{
    for ( Branch* b : branches ) for ( vector<NodePath*>& path : loci )
    {
        if ( b->used_ = ( find( path.begin(), path.end(), b->edge_->edge[drxn] ) != path.end() ) ) break;
    }
}

void Branch::sort( vector<Branch*>& branches )
{
    std::sort( branches.begin(), branches.end(), []( Branch* a, Branch* b )
    { return ( a->matched_.empty() ? 0 : a->matched_[0].second ) > ( b->matched_.empty() ? 0 : b->matched_[0].second ); } );
}

void Branch::sort()
{
    std::sort( matched_.begin(), matched_.end(), []( pair<Branch*, int> a, pair<Branch*, int> b ){ return a.second > b.second; } );
}

NodePath::NodePath( NodePath* np, vector<NodePath*>& paths, NodeRoll& nodes, int multi )
: offs_( np->offs_ ), verified_( np->verified_ ), id_( np->id_ ), multi_( multi ), loop_( np->loop_ ), score_( 0 ), bridge_( false ), blank_( false )
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
    for ( int d : { 0, 1 } ) for ( PathPairs* pp : np->paired_[d] ) new PathPairs( pp, this, !d );
    np->multi_ = max( 1, np->multi_ - multi );
    cover_ = Node::getCoverage( path_ );
}

NodePath::NodePath( Node* seed, vector<NodePath*>& paths, int32_t coord )
: offs_{ make_pair( seed, coord ) }, path_{ seed }, verified_( seed->verified_ ), cover_( 0 ), multi_( 0 ), loop_( 0 ), score_( 0 ), blank_( false )
{
    paths.push_back( this );
    coords_[0] = coords_[1] = ends_[1] = coord;
    ends_[0] = ends_[1] - seed->size();
    ended_[0] = ended_[1] = false;
}

NodePath::NodePath( NodePath* np, vector<NodePath*>& paths, NodeRoll& nodes )
: verified_( np->verified_ ), id_( np->id_ ), cover_( np->cover_ ), loop_( np->loop_ ), score_( np->score_ ), bridge_( false ), blank_( false )
{
    Nodes block[2];
    for ( int d : { 0, 1 } ) for ( PathEdge* pe : np->edges_[d] ) block[d] += pe->edge[d]->getTerminus( !d );
    for ( int d : { 0, 1 } ) ends_[d] = np->ends_[d];
    for ( int d : { 0, 1 } ) coords_[d] = np->coords_[d];
    
    for ( Node* node : np->path_ ) path_.push_back( new Node( node, nodes, node->drxn_, node->bad_ ) );
    
    for ( int i = 0; i < path_.size(); i++ )
    {
        int32_t* off = np->get( np->path_[i] );
        assert( off );
        offs_.insert( make_pair( path_[i], *off ) );
        for ( int d : { 0, 1 } ) for ( Edge& e : np->path_[i]->edges_[d] )
        {
            if ( !d && i && e.node == np->path_[i-1] ) continue;
            if ( ( d ? i+1 == path_.size() : !i ) && block[d].find( e.node ) ) continue;
            if ( d && i+1 < path_.size() && e.node == np->path_[i+1] ) path_[i]->addEdge( path_[i+1], e.ol, 1, false, e.leap );
            else path_[i]->addEdge( e.node, e.ol, d, false, e.leap );
        }
    }
    
    for ( int d : { 0, 1 } ) for ( PathEdge* pe : np->alts_[d] )
    {
        PathEdge* alt = NULL;
        for ( int i = 0; !alt && i < path_.size(); i++ ) if ( np->path_[i] == pe->node[!d] ) alt = new PathEdge( this, path_[i], pe, !d );
        assert( alt );
        alts_[d].push_back( alt );
        ( find( pe->edge[d]->edges_[!d].begin(), pe->edge[d]->edges_[!d].end(), pe ) != pe->edge[d]->edges_[!d].end()
                ? pe->edge[d]->edges_[!d] : pe->edge[d]->alts_[!d] ).push_back( alt );
    }
    
    paths.push_back( this );
}

//NodePath::NodePath( NodePath* np, PathEdge* l, PathEdge* r, NodeRoll& nodes )
//: verified_( np->verified_ ), id_( np->id_ ), loop_( np->loop_ ), score_( np->score_ ), bridge_( false ), blank_( false )
//{
//    Nodes block[2];
//    for ( int d : { 0, 1 } ) for ( PathEdge* pe : np->edges_[d] ) if ( pe != ( d ? r : l ) ) block[d] += pe->edge[d]->getTerminus( !d );
//    for ( int d : { 0, 1 } ) ends_[d] = np->ends_[d];
//    for ( int d : { 0, 1 } ) coords_[d] = np->coords_[d];
//    
//    for ( Node* node : np->path_ ) path_.push_back( new Node( node, nodes, node->drxn_, node->bad_ ) );
//    
//    for ( int i = 0; i < path_.size(); i++ )
//    {
//        int32_t* off = np->get( np->path_[i] );
//        assert( off );
//        offs_.insert( make_pair( path_[i], *off ) );
//        for ( int d : { 0, 1 } ) for ( Edge& e : np->path_[i]->edges_[d] )
//        {
//            if ( !d && i && e.node == np->path_[i-1] ) continue;
//            if ( ( d ? i+1 == path_.size() : !i ) && block[d].find( e.node ) ) continue;
//            if ( d && i+1 < path_.size() && e.node == np->path_[i+1] ) path_[i]->addEdge( path_[i+1], e.ol, 1, false, e.leap );
//            else path_[i]->addEdge( e.node, e.ol, d, false, e.leap );
//        }
//    }
//    
//    for ( int d : { 0, 1 } ) np->edges_[d].erase( remove( np->edges_[d].begin(), np->edges_[d].end(), ( d ? r : l ) ), np->edges_[d].end() );
//    for ( int d : { 0, 1 } ) np->getTerminus( d )->removeEdge( ( d ? r : l )->edge[d]->getTerminus( !d ), d, true );
//    l->edge[1] = r->edge[0] = this;
//    edges_[0].push_back( l );
//    edges_[1].push_back( r );
//    unordered_set<NodePath*> branch[2][2];
//    for ( int d : { 0, 1 } ) np->setBranch( branch[d][0], np->size()+params.maxMpMean, d );
//    for ( int d : { 0, 1 } ) setBranch( branch[d][1], np->size()+params.maxMpMean, d );
//    for ( int d : { 0, 1 } ) for ( int i = 0; i < np->paired_[d].size(); i++ )
//    {
//        bool found[2];
//        for ( int j : { 0, 1 } ) found[j] = branch[d][j].find( np->paired_[d][i]->node_[d] ) != branch[d][j].end();
//        if ( found[0] && !found[1] ) continue;
//        if ( !found[0] && found[1] )
//        {
//            paired_[d].push_back( np->paired_[d][i] );
//            np->paired_[d][i]->node_[!d] = this;
//            np->paired_[d].erase( np->paired_[d].begin() + i-- );
//        }
//        else new PathPairs( np->paired_[d][i], this, !d );
//    }
//    
//    for ( int d : { 0, 1 } ) for ( PathEdge* pe : np->alts_[d] )
//    {
//        PathEdge* alt = new PathEdge( this, pe, !d );
//        alts_[d].push_back( alt );
//        ( find( pe->edge[d]->edges_[!d].begin(), pe->edge[d]->edges_[!d].end(), pe ) != pe->edge[d]->edges_[!d].end()
//                ? pe->edge[d]->edges_[!d] : pe->edge[d]->alts_[!d] ).push_back( alt );
//    }
//    for ( int i = 0; i < path_.size(); i++ ) if ( np->path_[i]->verified_ )
//    {
//        np->path_[i]->reverify();
//        path_[i]->setVerified();
//    }
//    
//    cover_ = Node::getCoverage( path_ );
//}


NodePath::~NodePath()
{
    for ( int d : { 0, 1 } )
    {
        while ( !edges_[d].empty() ) delete edges_[d].back();
        while ( !alts_[d].empty() ) delete alts_[d].back();
//        while ( !pairs_[d].empty() ) delete pairs_[d].back();
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

//bool NodePath::addPair( NodePath* np, int32_t diff )
//{
//    for ( PathPairing* pp : pairs_[1] ) if ( pp->node[1] == np )
//    {
//        int i = 0;
//        while ( i < pp->diffs.size() && pp->diffs[i] < diff ) i++;
//        if ( pp->diffs[i] == diff ) return false;
//        pp->diffs.insert( pp->diffs.begin() + i, diff );
//        pp->missed.insert( pp->missed.begin() + i, 0 );
//        return true;
//    }
//    new PathPairing( this, np, diff );
//    return true;
//}

void NodePath::clonePairs( NodePath* base )
{
    unordered_set<NodePath*> branch[2][2];
    for ( int d : { 0, 1 } ) base->setBranch( branch[d][0], base->size()+params.maxMpMean, d );
    for ( int d : { 0, 1 } ) setBranch( branch[d][1], base->size()+params.maxMpMean, d );
    for ( int d : { 0, 1 } ) for ( int i = 0; i < base->paired_[d].size(); i++ )
    {
        bool found[2];
        for ( int j : { 0, 1 } ) found[j] = branch[d][j].find( base->paired_[d][i]->node_[d] ) != branch[d][j].end();
        if ( found[0] && !found[1] ) continue;
        if ( !found[0] && found[1] )
        {
            paired_[d].push_back( base->paired_[d][i] );
            base->paired_[d][i]->node_[!d] = this;
            base->paired_[d].erase( base->paired_[d].begin() + i-- );
        }
        else new PathPairs( base->paired_[d][i], this, !d );
    }
}

void NodePath::cleanPairs( vector<NodePath*>& paths )
{
    for ( NodePath* np : paths ) for ( int i = 0; i < np->paired_[1].size(); i++ )
    {
        if ( !np->paired_[1][i]->score_ 
                && np->paired_[1][i]->uniques_.empty() 
                && np->paired_[1][i]->shared_.empty()
                && np->paired_[1][i]->ids_.empty() 
                && np->paired_[1][i]->mates_.empty()) delete np->paired_[1][i--];
    }
}

void NodePath::create( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& seeds, vector<NodePath*>& paths )
{
    for ( NodePath* np : paths ) delete np;
    seeds.clear();
    paths.clear();
    
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
    
    setFills( paths );
    for ( NodePath* np : paths ) np->cover_ = Node::getCoverage( np->path_ );
    for ( NodePath* np : paths )
    {
        vector<NodePath*> path;
        np->setPairs( np, path, 0, params.maxMpMean * 1.2 + 300, false, false );
    }
    
    for ( NodePath* np : paths ) np->seq_ = Node::getSeq( np->path_ );
    for ( NodePath* np : paths ) np->setScore();
    for ( NodePath* np : paths ) np->setMates( bwt, nodes, paths );
//    for ( NodeMark& nm : seeds[5]->path_[0]->pe_[0] ) cout << ">" << nm.id << endl << bwt.getSequence( nm.id ) << endl;
//    setReads( bwt, nodes, paths );
    cleanPairs( paths );
    setEdges( paths );
}

void NodePath::cull( vector<NodePath*>& seed, vector<NodePath*>& paths, SeedExtend exts[2], NodeRoll& nodes )
{
    Lib* lib;
    int d;
    Coords* coords;
    for ( NodePath* a : paths ) for ( Node* node : a->path_ ) for ( auto& read: node->reads_ ) if ( ( lib = params.getLib( read.first ) ) && !lib->isPe )
    {
        ReadId id = read.first;
        if ( lib->getPair( id, d ) && d ) for ( NodePath* b : paths ) if ( a != b ) for ( Node* t : b->path_ ) if ( coords = t->getRead( id ) )
        {
            auto ins = a->tmp_.insert( make_pair( b, unordered_set<ReadId>{ id } ) );
            if ( !ins.second ) ins.first->second.insert( id );
            auto ins2 = b->tmp_.insert( make_pair( a, unordered_set<ReadId>{ id } ) );
            if ( !ins2.second ) ins.first->second.insert( id );
        }
    }
    for ( int i = 0; i < paths.size(); i++ ) if ( paths[i]->id_ = i+1 );
    for ( int d : { 0, 1 } ) exts[d].reset();
    
    for ( int d : { 0, 1 } ) if ( d )
    {
        unordered_set<NodePath*> pathed, base;
        for ( NodePath* np : seed ) if ( np->id_ > 1 ) for ( PathEdge* pe : np->edges_[!d] ) pe->edge[!d]->setReachable( pathed, !d );
        for ( NodePath* np : seed ) if ( np->id_ == 1 ) for ( int i : { 0, 1 } ) np->setReachable( base, i );
        base.insert( pathed.begin(), pathed.end() );
        for ( NodePath* np : seed ) if ( np->id_ > 1 ) np->cull( paths, pathed, base, exts[d], nodes, d );
        
        pathed.clear();
        for ( NodePath* np : seed ) for ( PathEdge* pe : np->edges_[!d] ) pe->edge[!d]->setReachable( pathed, !d );
        for ( NodePath* np : seed ) if ( np->id_ > 1 ) np->cull( paths, pathed, nodes, d );
        
        unordered_set<NodePath*> bad, good;
        for ( NodePath* np : seed )
        {
            unordered_set<NodePath*> bck;
            np->setReachable( bck, !d );
            pathed.insert( bck.begin(), bck.end() );
        }
        for ( NodePath* np : paths ) for ( auto& ids : np->tmp_ ) if ( pathed.find( ids.first ) != pathed.end() ) np->setReachable( good, d );
        for ( NodePath* np : pathed ) for ( PathEdge* pe : np->edges_[!d] ) if ( pathed.find( pe->edge[!d] ) == pathed.end() )
        {
            if ( good.find( pe->edge[!d] ) == good.end() )  np->setReachable( bad, !d );
        }
        for ( NodePath* np : pathed ) bad.erase( np );
        for ( NodePath* np : good ) bad.erase( np );
        for ( NodePath* np : bad )
        {
            paths.erase( remove( paths.begin(), paths.end(), np ), paths.end() );
            for ( Node* node : np->path_ ) nodes.erase( node );
            delete np;
        }
    }
}

void NodePath::cull( vector<NodePath*>& seed, vector<NodePath*>& paths, NodeRoll& nodes )
{
    unordered_set<NodePath*> sets[2], trash;
    Coords* coords;
    Lib* lib;
    for ( int d : { 0, 1 } ) for ( NodePath* np : seed ) np->setReachable( sets[d], d );
    for ( NodePath* np : paths ) if ( sets[0].find( np ) == sets[0].end() && sets[1].find( np ) == sets[1].end() )
    {
        for ( int d : { 0, 1 } ) if ( np->edges_[!d].empty() )
        {
            bool good = false;
            for ( Node* node : np->path_ ) if ( !good ) for ( auto& read : node->reads_ ) if ( ( lib = params.getLib( read.first ) ) && !lib->isPe )
            {
                bool bad = false;
                ReadId id = params.getPairId( read.first );
                for ( NodePath* alt : paths ) if ( !bad && np != alt ) for ( Node* t : alt->path_ ) if ( bad = ( coords = t->getRead( read.first ) ) ) break;
                if ( !bad ) for ( NodePath* alt : paths ) if ( !good && np != alt ) for ( Node* t : alt->path_ ) if ( good = ( coords = t->getRead( id ) ) ) break;
            }
            if ( !good ) trash.insert( np );
        }
    }
    for ( NodePath* np : trash )
    {
        paths.erase( remove( paths.begin(), paths.end(), np ), paths.end() );
        paths.erase( remove( paths.begin(), paths.end(), np ), paths.end() );
        delete np;
    }
}

void NodePath::cull( vector<NodePath*>& paths, unordered_set<NodePath*>& pathed, unordered_set<NodePath*>& base, SeedExtend& ext, NodeRoll& nodes, bool drxn )
{
    if ( !pathed.insert( this ).second ) return;
    struct CullEdge
    {
        PathEdge* edge;
        unordered_set<NodePath*> reached;
        unordered_set<ReadId> hits;
    };
    
    score_ = 0;
    
    if ( edges_[drxn].size() > 1 )
    {
        vector<CullEdge> scores;
        unordered_set<NodePath*> tar( base.begin(), base.end() );
        setReachable( tar, !drxn );
        for ( PathEdge* pe : edges_[drxn] )
        {
            CullEdge edge;
            edge.edge = pe;
            pe->score = pe->shared = 0;
            pe->edge[drxn]->setReachable( edge.reached, drxn );
            for ( NodePath* np : edge.reached ) for ( const pair<NodePath*, unordered_set<ReadId> >& ids : np->tmp_ )
            {
                if ( tar.find( ids.first ) != tar.end() ) edge.hits.insert( ids.second.begin(), ids.second.end() );
            }
            scores.push_back( edge );
        }
        
        sort( scores.begin(), scores.end(), []( CullEdge& a, CullEdge& b ){ return a.hits.size() > b.hits.size(); } );
        
        if ( ( score_ = scores[0].hits.size() ) > 2 ) for ( int i = 1; i < scores.size(); i++ )
        {
            for ( ReadId id : scores[i].hits ) ( scores[0].hits.find( id ) == scores[0].hits.end() ? scores[i].edge->score : scores[i].edge->shared )++;
            int prime = scores[0].hits.size() - scores[i].edge->shared;
            if ( i == 1 || prime < scores[0].edge->score ) scores[0].edge->score = prime;
            
            cout << scores[0].edge->edge[drxn]->id_ << " [" << prime << "]   vs   " << scores[i].edge->edge[drxn]->id_ << " [" << scores[i].edge->score << "]" << endl;
            if ( scores[i].edge->score || scores[i].edge->shared > 2 || prime < 3 )
            {
                if ( scores[i].edge->shared || prime < 50 || scores[0].edge->score > 1 ) continue;
            }
            unordered_set<NodePath*> safe;
            for ( int j = 0; j < scores.size(); j++ ) if ( i != j ) safe.insert( scores[j].reached.begin(), scores[j].reached.end() );
            for ( NodePath* np : scores[i].reached ) if ( safe.find( np ) == safe.end() )
            {
                
                for ( Node* node : np->path_ ) nodes.erase( node );
                paths.erase( remove( paths.begin(), paths.end(), np ), paths.end() );
                pathed.erase( np );
                delete np;
            }
            scores.erase( scores.begin() + i-- );
        }
        PathEdge::sort( edges_[drxn] );
        
        if ( scores[0].hits.size() < 2 ) ext.addFork( drxn ? path_.back() : path_[0] );
        else for ( CullEdge& ce : scores ) if ( ce.hits.size() > 1 ) ce.edge->edge[drxn]->cull( paths, pathed, base, ext, nodes, drxn );
    }
    else if ( !edges_[drxn].empty() ) edges_[drxn][0]->edge[drxn]->cull( paths, pathed, base, ext, nodes, drxn );
}

void NodePath::cull( vector<NodePath*>& paths, unordered_set<NodePath*>& pathed, NodeRoll& nodes, bool drxn )
{
    if ( !pathed.insert( this ).second ) return;
    
    if ( edges_[drxn].size() > 1 && edges_[drxn][0]->score && score_ > 3 )
    {
        unordered_set<NodePath*> base;
        edges_[drxn][0]->edge[drxn]->setReachable( base, drxn );
        for ( int i = 1; i < edges_[drxn].size(); i++ ) if ( edges_[drxn][i]->score < edges_[drxn][0]->score )
        {
            int best = 0, alts = 0;
            unordered_set<NodePath*> q;
            edges_[drxn][i]->edge[drxn]->setReachable( base, drxn );
            for ( NodePath* np : base ) if ( np->edges_[drxn].size() > 1 && q.find( np ) == q.end() ) best = max( best, np->edges_[drxn][1]->score );
            if ( edges_[drxn][i]->score < min( best, 2 ) || ( alts && !edges_[drxn][i]->score ) )
            {
                PathEdge::sever( edges_[drxn][i]->edge[0], edges_[drxn][i]->edge[1] );
                i--;
            }
            else if ( best < min( edges_[drxn][i]->score, 2 ) )
            {
                for ( NodePath* np : base ) if ( np->edges_[drxn].size() > 1 && q.find( np ) == q.end() ) for ( int j = 1; j < np->edges_[drxn].size(); j++ )
                {
                    if ( np->edges_[drxn][0]->score > np->edges_[drxn][j]->score )
                    {
                        PathEdge::sever( edges_[drxn][i]->edge[0], edges_[drxn][i]->edge[1] );
                        j--;
                    }
                }
            }
            else
            {
                int x = 0;
            }
        }
    }
    for ( PathEdge* pe : edges_[drxn] ) pe->edge[drxn]->cull( paths, pathed, nodes, drxn );
}

void NodePath::destroy( NodeRoll& nodes )
{
    for ( Node* node : path_ ) for ( int d : { 0, 1 } ) node->clearEdges( d );
    for ( Node* node : path_ ) nodes.erase( node );
    delete this;
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
    if ( !branch || ( !params.rna && !verified ) || edges.empty() ) return;
    
    for ( Edge& e : edges ) addEdge( e, paths, branch, drxn );
    for ( PathEdge* pe : edges_[drxn] ) pe->edge[drxn]->extend( paths, true, !drxn );
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

//void NodePath::forceCut( vector<NodePath*>& paths, NodeRoll& nodes, int i, int j, bool drxn )
//{
//    NodePath* edges[2]{ NULL, NULL };
//    for ( NodePath* np : paths )
//    {
//        if ( np->id_ == i ) edges[0] = np;
//        if ( np->id_ == j ) edges[1] = np;
//    }
//    
//    PathEdge* edge = NULL;
//    assert( edges[0] && edges[1] );
//    for ( PathEdge* pe : edges[0]->edges_[drxn] ) if ( pe->edge[drxn] == edges[1] ) edge = pe;
//    assert( edge );
//    
//    Node* node = drxn ? edges[0]->path_.back() : edges[0]->path_[0];
//    Node* clone = new Node( node, nodes, node->drxn_, true );
//    clone->addEdge( drxn ? edges[1]->path_[0] : edges[1]->path_.back(), edge->ol, drxn, false, edge->leap );
//    PathEdge::sever( edge->edge[0], edge->edge[1] );
//}

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

int32_t NodePath::getLen( vector<NodePath*>& path )
{
    int32_t len = path[0]->size();
    for ( int i = 1; i < path.size(); i++ )
    {
        int ol = params.readLen;
        for ( PathEdge* pe : path[i]->edges_[0] ) if ( pe->edge[0] == path[i-1] ) ol = min( ol, pe->ol );
        len += path[i]->size() - ol;
    }
    return len;
}

int32_t NodePath::getCoord( Edge& e, bool drxn )
{
    return drxn ? coords_[1] + e.node->size() - e.ol : coords_[0] - path_[0]->size() + e.ol;
}

int32_t NodePath::getCoord( Node* node, Edge& e, int32_t dist, bool drxn )
{
    return dist + ( drxn ? e.node->size() - e.ol : -node->size() + e.ol );
}

float NodePath::getCoverage( vector<NodePath*>& path, int minLen, int drxn )
{
    vector<Node*> pathed;
    for ( NodePath* np : path ) pathed.insert( drxn ? pathed.end() : pathed.begin(), np->path_.begin(), np->path_.end() );
    return Node::getCoverage( pathed, minLen );
}

Node* NodePath::getEnd( Querier& bwt, NodeRoll& nodes, unordered_set<NodePath*>& used, bool drxn )
{
    for ( auto& no : offs_ ) if ( no.first->isContinue( drxn ) ) return drxn ? path_.back() : path_[0];
    vector< pair<Node*, int32_t> > exts = ( drxn ? path_.back() : path_[0] )->getExtendable( 300, drxn );
    if ( !exts.empty() )
    {
        return drxn ? path_.back() : path_[0];
    }
    unordered_map<Node*, int32_t> offs;
    getOffsets( offs, used, ends_[drxn], 0, max( params.maxPeMax, params.maxMpMean ), drxn );
    return NULL;
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

int32_t NodePath::getOffset( Node* node, bool drxn )
{
    int32_t coord = 0;
    for ( int i = 0; i < path_.size(); i++ )
    {
        if ( i ) coord += path_[i-1]->size() - path_[i]->getOverlap( path_[i-1], 0 );
        if ( path_[i] == node ) return drxn ? size() - path_[i]->size() - coord : coord;
    }
    return 0;
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

string NodePath::getSeq( vector<NodePath*>& path )
{
    vector<Node*> pathed;
    for ( NodePath* np : path ) pathed.insert( pathed.end(), np->path_.begin(), np->path_.end() );
    return Node::getSeq( pathed );
}

Node* NodePath::getTerminus( bool drxn )
{
    assert( !path_.empty() );
    return drxn ? path_.back() : path_[0];
}

bool NodePath::isBranchable( bool drxn )
{
    if ( edges_[drxn].size() < 2 ) return false;
    for ( PathEdge* pe : edges_[drxn] ) if ( pe->score < 2 || ( pe->score < 4 && !pe->edge[drxn]->verified_ ) ) return false;
    return true;
}

void NodePath::naturalise( vector<NodePath*>& paths )
{
    for ( NodePath* np : paths )
    {
        int32_t diff = np->path_.back()->ends_[1] - np->ends_[1];
        if ( !diff ) continue;
        for ( auto it = np->offs_.begin(); it != np->offs_.end(); it++ ) it->second += diff;
        for ( int d : { 0, 1 } )
        {
            np->ends_[d] += diff;
            np->coords_[d] += diff;
            for ( PathEdge* pe : np->edges_[d] ) pe->diff += d ? -diff : diff;
        }
    }
}

void NodePath::print( vector<NodePath*>& paths, vector< vector<NodePath*> >& loci )
{
    int32_t minCoord = 0;
    for ( int i = 0; i < paths.size(); i++ ) paths[i]->id_ = i+1;
    for ( NodePath* np : paths ) minCoord = min( minCoord, np->ends_[0] );
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) while ( !np->paired_[d].empty() ) delete np->paired_[d].back();
    for ( NodePath* np : paths ) for ( Node* node : np->path_ ) node->setVerified();
    for ( NodePath* np : paths )
    {
        vector<NodePath*> path;
        np->setPairs( np, path, 0, params.maxMpMean * 1.2 + 300, false, false );
    }
    
    for ( NodePath* np : paths ) np->setScore();
    cleanPairs( paths );
    setEdges( paths );
    
    unordered_set<NodePath*> used;
    int count = 0;
    for ( vector<NodePath*>& path : loci )
    {
        int32_t coords[3]{ path[0]->ends_[0], path.back()->ends_[1], NodePath::getLen( path ) };
        cout << ">Allele " << ++count << ",  Coords: " << coords[0] << " to " << coords[1] << ",  Length: " << coords[2] << endl << NodePath::getSeq( path ) << endl;
        for ( NodePath* np : path ) np->print( minCoord );
        used.insert( path.begin(), path.end() );
    }
    
    for ( NodePath* np : paths ) if ( used.find( np ) == used.end() ) np->print( minCoord );
}

void NodePath::print( vector<NodePath*>& paths )
{
    ofstream ofs( "/home/glen/Genomes/HpDump/PathAlign.fa" );
    
    int32_t minCoord = 0;
    for ( int i = 0; i < paths.size(); i++ ) paths[i]->id_ = i+1;
    for ( NodePath* np : paths ) minCoord = min( minCoord, np->ends_[0] );
    
    for ( NodePath* np : paths ) np->print( minCoord, ofs );
    for ( NodePath* np : paths ) np->printUnpaired( paths );
}

void NodePath::print( int32_t minCoord, ofstream& ofs )
{
    bool seeded = false;
    for ( Node* node : path_ ) if ( node->drxn_ == 2 ) seeded = true;
    string seq = Node::getSeq( path_ );
//    float cover = Node::getCoverage( path_ );
    float cover = Node::getCoverageMedian( path_ );
    int reads = 0;
    for ( Node* node : path_ ) reads += node->countReads( true );
    ofs << ">" << ( seeded ? "SEED   " : "" ) << id_ << ") Coords: (" << ends_[0] << ")-(" << ends_[1] << "), Cover: " << (int)cover << "[" << reads << "]";
    ofs << ", L:" << ( edges_[0].empty() ? " NONE" : "" );
    for ( PathEdge* pe : edges_[0] ) ofs << " (" << pe->edge[0]->id_ << ")";
    ofs << ", R:" << ( edges_[1].empty() ? " NONE" : "" );
    for ( PathEdge* pe : edges_[1] ) ofs << " (" << pe->edge[1]->id_ << ")";
    
    for ( int d : { 0, 1 } ) if ( !paired_[d].empty() )
    {
        ofs << ( d ? ",   Right pairs:" : ",   Left pairs:" );
        for ( PathPairs* pp : paired_[d] )
        {
            int count = pp->ids_.size();
            for ( ReadId id : pp->mates_ ) if ( pp->ids_.find( id ) == pp->ids_.end() ) count++;;
            ofs << "  " << pp->node_[d]->id_ << "(" << pp->ids_.size() << ")";
        }
    }
    ofs << endl << string( ends_[0] - minCoord, '-' )<< seq << endl;
    
    
}

void NodePath::print( int32_t minCoord )
{
    string seq = Node::getSeq( path_ );
    float cover = Node::getCoverage( path_ );
    cout << ">" << id_ << ") Coords: (" << ends_[0] << ")-(" << ends_[1] << "), Cover: " << cover;
    cout << ", L:" << ( edges_[0].empty() ? " NONE" : "" );
    for ( PathEdge* pe : edges_[0] ) cout << " (" << pe->edge[0]->id_ << ")";
    cout << ", R:" << ( edges_[1].empty() ? " NONE" : "" );
    for ( PathEdge* pe : edges_[1] ) cout << " (" << pe->edge[1]->id_ << ")";
//    if ( !verified_ ) cout << " UNVERIFIED";

    for ( int d : { 0, 1 } ) if ( !paired_[d].empty() )
    {
        cout << ( d ? ",   Right pairs:" : ",   Left pairs:" );
        for ( PathPairs* pp : paired_[d] ) cout << "  " << pp->node_[d]->id_ << "(" << pp->ids_.size() << ")";
    }
    cout << endl << string( ends_[0] - minCoord, '-' )<< seq << endl;
}

void NodePath::printUnpaired( vector<NodePath*>& paths )
{
    struct Unpaired
    {
        ReadId id;
        int32_t est;
    };
    unordered_map<NodePath*, vector<Unpaired> > pairs[2];
    Lib* lib;
    Coords* coords;
    Unpaired up;
    for ( Node* node : path_ ) for ( auto& read : node->reads_ ) if ( ( lib = params.getLib( read.first ) ) && !lib->isPe )
    {
        int d;
        up.id = read.first;
        lib->getPair( up.id, d );
        up.est = lib->size - abs( node->ends_[d] - read.second[!d] ) - getOffset( node, d );
        for ( NodePath* np : paths ) if ( np != this ) for ( Node* t : np->path_ ) if ( coords = t->getRead( up.id ) )
        {
            Unpaired match = up;
            match.est -= abs( t->ends_[!d] - coords->coords[d] ) + np->getOffset( t, !d );
            auto ins = pairs[d].insert( make_pair( np, vector<Unpaired>{ match } ) );
            if ( !ins.second ) ins.first->second.push_back( match );
        }
    }
//    vector< pair<int, vector<int> > > pairs[2][2];
//    Coords* coords;
//    for ( Node* node : path_ ) for ( int d : { 0, 1 } )
//    {
//        vector< vector<NodeMark>* > marks{ &node->mp_[d] };
//        if ( edges_[d].empty() ) marks.push_back( &node->pe_[d] );
//        for ( vector<NodeMark>* m : marks ) for ( NodeMark& nm : *m )
//        {
//            ReadId id[2]{ nm.id, params.getRevId( nm.id ) };
//            for ( int i : { 0, 1 } ) for ( NodePath* np : paths ) for ( Node* t : np->path_ ) if ( coords = t->getRead( id[i] ) )
//            {
//                bool added = false;
//                for ( pair<int, vector<int> >& p : pairs[d][i] ) if ( added = p.first == np->id_ )
//                {
//                    p.second.push_back( nm.dist );
//                    break;
//                }
//                if ( !added ) pairs[d][i].push_back( make_pair( np->id_, vector<int>{ nm.dist } ) );
//            }
//        }
//    }
    
    for ( int d : { 0, 1 } ) if ( !pairs[d].empty() )
    {
        cout << "Node " << id_ << ( d ? "   RIGHT   " : "   LEFT   " ) << endl;
        for ( pair<NodePath*, vector<Unpaired> > p : pairs[d] )
        {
            cout << "    " << p.first->id_ << ":";
            sort( p.second.begin(), p.second.end(), []( Unpaired& a, Unpaired& b ){ return a.est < b.est; } );
            unordered_set<ReadId> used;
            for ( Unpaired& u : p.second ) if ( used.insert( u.id ).second ) cout << "  " << u.est << " [" << u.id << "]";
            cout << endl;
        }
//        for ( pair<int, vector<int> >& p : pairs[d][i] )
//        {
//            cout << "   " << p.first << ":";
//            for ( int j : p.second ) cout << "  " << j;
//            cout << endl;
//        }
    }
}

bool NodePath::isContinue( int32_t dist, bool drxn )
{
    dist -= size();
    if ( dist < 0 ) return false;
    for ( PathEdge* pe : edges_[drxn] ) if ( pe->edge[drxn]->isContinue( dist + pe->ol, drxn ) ) return true;
    if ( edges_[drxn].empty() ) for ( auto& no : offs_ ) if ( ( drxn ? no.second - ends_[1] : ends_[0] - ( no.second - no.first->size() ) ) < dist )
    {
        if ( no.first->isContinue( drxn ) ) return true;
    }
    return false;
}

bool NodePath::isEnding( PathEdge* edge, float cover, bool drxn )
{
    if ( edges_[drxn].empty() ) return true;
    if ( size() > params.readLen * 1.5 ) return cover_ * 5 < cover;
    return getTerminus( !drxn )->getCoverage( drxn ) * 5 < cover;
}

bool NodePath::isEnding( int dist, bool drxn )
{
    dist -= max( 1, size()-params.readLen );
    if ( dist < 0 ) return false;
    for ( PathEdge* pe : edges_[drxn] ) if ( !pe->edge[drxn]->isEnding( dist, drxn ) ) return false;
    return blank_ = true;
}

bool NodePath::isForked( PathEdge* base, bool drxn )
{
    for ( PathEdge* pe : edges_[drxn] ) if ( pe != base ) return true;
    return false;
}

bool NodePath::isPathed( Node* node )
{
    return find( path_.begin(), path_.end(), node ) != path_.end();
}

bool NodePath::merge( vector<NodePath*>& paths )
{
    bool merged = false;
    for ( NodePath* np : paths ) while ( np->merge() ) merged = true;
    int deleted = 0, culled = 0;
    for ( int i = 0; i < paths.size(); i++ ) if ( paths[i]->path_.empty() )
    {
        delete paths[i];
        paths.erase( paths.begin() + i-- );
        deleted++;
    }
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) PathEdge::sort( np->edges_[d] );
    return merged || deleted;
}

bool NodePath::merge()
{
    if ( edges_[1].size() != 1 || edges_[1][0]->bad || edges_[1][0]->edge[1]->edges_[0].size() != 1 ) return false;
    PathEdge* edge = edges_[1][0];
    NodePath* np = edge->edge[1];
    
    path_.insert( path_.end(), np->path_.begin(), np->path_.end() );
    
    for ( auto& no : np->offs_ ) offs_.insert( make_pair( no.first, no.second - edge->diff ) );
    ends_[1] = np->ends_[1] - edge->diff;
    coords_[1] = np->coords_[1] - edge->diff;
    score_ += np->score_;
    
    for ( int d : { 0, 1 } ) while ( !np->paired_[d].empty() )
    {
        if ( np->paired_[d].back()->steal( this, d ) ) np->paired_[d].pop_back();
        else delete np->paired_[d].back();
    }
    
    delete edge;
    for ( PathEdge* pe : np->edges_[1] )
    {
        pe->edge[0] = this;
        edges_[1].push_back( pe );
    }
    for ( int d : { 0, 1 } ) for ( PathEdge* pe : np->alts_[d] )
    {
        pe->edge[!d] = this;
        alts_[d].push_back( pe );
    }
    np->edges_[1].clear();
    for ( int d : { 0, 1 } ) np->alts_[d].clear();
    np->path_.clear();
    cover_ = Node::getCoverage( path_ );
    
    return true;
}

void NodePath::reduce( bool drxn )
{
    if ( edges_[drxn].size() < 2 || edges_[drxn][0]->score < 10 ) return;
    
    static int culls = 0;
    int base = edges_[drxn][0]->score;
    float ploidy = cover_ * 2 / params.cover;
    
    for ( int i = 1; i < edges_[drxn].size(); i++ ) if ( edges_[drxn][i]->score < min( 3, base / 10 )  )
    {
        PathEdge* pe = edges_[drxn][i];
        bool blunt = pe->node[drxn]->isBlunt( 0, 5, drxn );
        cout << ++culls << ( blunt ? "   BLUNT" : "    CONT" ) << "    " << ploidy << endl;
        int x = 0;
    }
    
    
    int limits[2]{ 10, 0 };
//    unordered_set<NodePath*> t;
//    setBranch( t, params.maxPeMean, !drxn );
//    for ( PathEdge* pe : edges_[drxn] ) if ( pe->score < 0 )
//    {
//        pe->score = 0;
//        unordered_set<NodePath*> q;
//        pe->edge[drxn]->setBranch( q, params.maxPeMean, drxn );
//        for ( NodePath* np : q ) for ( PathPairs* pp : np->paired_[!drxn] ) if ( t.find( pp->node_[!drxn] ) != t.end() ) pe->score += pp->score_;
//    }
    for ( PathEdge* pe : edges_[drxn] ) limits[0] = min( limits[0], pe->score );
    for ( PathEdge* pe : edges_[drxn] ) limits[1] = max( limits[1], pe->score );
    if ( limits[1] > 9 ) for ( int i = 0; i < edges_[drxn].size(); i++ ) if ( !edges_[drxn][i]->score ) edges_[drxn][i--]->downgrade( drxn );
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

vector< vector<NodePath*> > NodePath::resolve( Querier& bwt, vector<NodePath*>& paths, NodeRoll& nodes )
{
    setBlunts( bwt, paths, nodes );
    setCrosses( paths, nodes );
    setUnforked( paths, nodes );
//    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) np->unfork( d );
    merge( paths );
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) PathEdge::sort( np->edges_[d] );
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) np->setBridges( d );
    merge( paths );
    setCrosses( paths, nodes );
    setBlanks( paths );
    vector< vector<NodePath*> > loci = setPaths( paths );
//    setCulled( loci, paths, nodes );
    
    int count = 0;
//    for ( vector<NodePath*>& path : loci )
//    {
//        int32_t coords[3]{ path[0]->ends_[0], path.back()->ends_[1], NodePath::getLen( path ) };
//        cout << ">Allele " << ++count << ",  Coords: " << coords[0] << " to " << coords[1] << ",  Length: " << coords[2] << endl << NodePath::getSeq( path ) << endl;
//    }
    
    unordered_set<NodePath*> used;
    vector<NodePath*> blocked;
    for ( vector<NodePath*>& path : loci ) used.insert( path.begin(), path.end() );
    for ( NodePath* np : paths ) if ( used.find( np ) == used.end() && np->size() > 200 ) blocked.push_back( np );
    bool verified = false;
    for ( vector<NodePath*>& path : loci ) for ( NodePath* np : path ) if ( !np->verified_ ) for ( Node* node : np->path_ ) if ( node->setVerified() ) verified = true ;
    if ( verified ) Node::verify( nodes );
    
    static bool cheat = true;
//    if ( cheat )
//    {
//        vector<NodePath*> path = { loci[13][1] };
//        vector<PathEdge*> taken[2][2];
//        taken[0][1].push_back( path[0]->edges_[0][1] );
//        taken[1][1].push_back( path[0]->edges_[1][1] );
//        setCloned( path, taken, paths, nodes );
//    }
//    if ( cheat )
//    {
//        vector<NodePath*> path = { loci[3][0] };
//        vector<PathEdge*> taken[2][2];
//        taken[1][1].push_back( path[0]->edges_[1][1] );
//        setCloned( path, taken, paths, nodes );
//    }
    cheat = false;
    return loci;
}

vector< vector<NodePath*> > NodePath::resolve( Querier& bwt, vector<NodePath*>& paths, vector<NodePath*>& seeds, NodeRoll& nodes )
{
//    vector<int> libs( params.libs.size(), 0 );
//    for ( Node* node : nodes.nodes )
//    {
//        for ( int i = 0; i < params.libs.size(); i++ ) for ( int d : { 0, 1 } ) for ( auto& np : node->hits_.pairs[d] ) libs[i] += np.second.libs[i];
//    }
    
    vector< vector<NodePath*> > loci;
    vector<int> lens;
    for ( NodePath* seed : seeds )
    {
        vector<NodePath*> locus{ seed };
        vector< vector<NodePath*> > haps;
        if ( abs( seed->cover_ - params.cover ) < abs( seed->cover_ - ( params.cover / 2 ) ) )
        {
            
        }
        else haps.push_back( vector<NodePath*>{ seed } );
        
        for ( int d : { 0, 1 } ) for ( NodePath* ext = ( d ? locus.back() : locus[0] ); ext && !ext->edges_[d].empty(); )
        {
            PathEdge::sort( ext->edges_[d] );
            if ( !ext->edges_[d][0]->score ) break;
            ext = ext->edges_[d][0]->edge[d];
            locus.insert( d ? locus.end() : locus.begin(), ext );
        }
        lens.push_back( NodePath::getLen( locus ) );
        loci.push_back( locus );
    }
    return loci;
}

void NodePath::setBlanks( vector<NodePath*>& paths )
{
    sort( paths.begin(), paths.end(), []( NodePath* a, NodePath* b ) { return a->score_ > b->score_; } );
    
    unordered_map<NodePath*, float> pathed;
    for ( NodePath* np : paths ) if ( !np->blank_ && pathed.find( np ) == pathed.end() )
    {
        np->setBlank( pathed, np->size() > params.readLen * 1.5 ? np->cover_ : Node::getCoverage( np->path_, params.readLen * 1.5 ), false );
    }
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) if ( !np->blank_ && np->isEnding( params.readLen / 2, d ) ) break;
}

void NodePath::setBlank( unordered_map<NodePath*, float>& pathed, float cover, bool blank )
{
    auto it = pathed.insert( make_pair( this, cover ) );
    if ( !it.second )
    {
        if ( !blank ) return;
        if ( it.first->second >= cover ) return;
        it.first->second = cover;
    }
    
    if ( blank && size() > params.readLen + 5 && ( cover_ * 5 > cover ) ) return;
    
    if ( blank_ ) assert( blank );
    blank_ = blank;
    
    for ( int d : { 0, 1 } ) for ( PathEdge* pe : alts_[d] ) pe->edge[d]->setBlank( pathed, cover, true );
    if ( blank ) for ( int d : { 0, 1 } ) for ( PathEdge* pe : edges_[d] ) if ( !pe->bad ) pe->edge[d]->setBlank( pathed, cover, true );
    if ( !blank ) for ( int d : { 0, 1 } ) for ( PathEdge* pe : edges_[d] ) if ( !pe->bad )
    {
        vector<NodePath*> path{ pe->edge[0], pe->edge[1] };
        float ext = path[d]->size() > params.readLen * 1.5 ? path[d]->cover_ : NodePath::getCoverage( path, params.readLen * 1.5 );
        blank = path[d]->cover_ * 5 < cover;
        if ( !blank || pe->edge[d]->isEnding( pe, cover, d ) ) pe->edge[d]->setBlank( pathed, blank ? cover : ext, blank );
    }
    
}

bool NodePath::setBlunts( Querier& bwt, vector<NodePath*>& paths, NodeRoll& nodes )
{
    int blunted = 0, culled = 0;
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) PathEdge::sort( np->edges_[d] );
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) for ( int i = 1; i < np->edges_[d].size(); i++ ) if ( np->edges_[d][i]->score < 3 )
    {
        PathEdge* base = np->edges_[d][0],* pe = np->edges_[d][i];
        if ( base->score / 5 <= pe->score ) continue;
        int cover = pe->node[!d]->getCoverage( d );
        if ( cover < 200 ) continue;
        int ol = params.readLen - max( 6, params.readLen * 300 / ( cover+500 ) );
        int diff = params.readLen - pe->ol;
        int exts = 0;
        bool bad = pe->ol < ol ;
        if ( !bad )
        {
            ol = min( ol + ( params.readLen - pe->ol ), params.readLen-5 );
            string seq = pe->node[d]->getSeqEnd( min( params.readLen, pe->node[d]->size() ),!d );
            exts = bwt.getExtendable( seq, ol, d );
            bad = exts < 2 + sqrt( cover / 200 );
        }
        if ( !bad ) continue;
        
        assert( pe->sever() );
        delete pe;
        i--;
        blunted++;
    }
    if ( !blunted ) return false;
    Node::updateStates( nodes );
    for ( int i = 0; i < paths.size(); i++ )
    {
        bool bad = true;
        for ( Node* node : paths[i]->path_ ) if ( !( bad = node->bad_ ) ) break;
        if ( bad ) for ( int d : { 0, 1 } ) for ( PathEdge* pe : paths[i]->edges_[d] ) if ( !pe->edge[d]->getTerminus( !d )->bad_ ) bad = false;
        if ( !bad ) continue;
        delete paths[i];
        paths.erase( paths.begin() + i-- );
        culled++;
    }
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

bool NodePath::setBridges( bool drxn )
{
    if ( edges_[drxn].size() < 2 ) return false;
    int score = edges_[drxn][0]->score;
    vector<NodePath*> path;
    for ( int i = 1; i < edges_[drxn].size(); i++ ) if ( !edges_[drxn][i]->bad && edges_[drxn][i]->edge[drxn]->setBridge( score, edges_[drxn][i]->score, drxn ) )
    {
        edges_[drxn][i--]->downgrade( drxn );
    }
}

bool NodePath::setBridge( int base, int score, bool drxn )
{
    if ( edges_[drxn].empty() ) return false;
    score = max( score, max( edges_[0][0]->score, edges_[1][0]->score ) );
    if ( base / 3 <= score ) return false;
    PathEdge* pe = edges_[drxn][0];
    bool bridged = false;
    if ( pe->bad ) bridged = bridge_ = true;
    else if ( pe->edge[drxn]->edges_[!drxn][0] != pe )
    {
        int alt = pe->edge[drxn]->edges_[!drxn][0]->score;
        if ( alt / 3 <= score ) return false;
        if ( ( alt + base ) / 8 <= score ) return false;
        pe->downgrade( !drxn );
        bridged = bridge_ = true;
    }
    else if ( pe->edge[drxn]->setBridge( base, score, drxn ) ) bridged = bridge_ = true;
    
    return bridged;
}

void NodePath::setCloned( vector<NodePath*>& path, vector<PathEdge*> taken[2][2], vector<NodePath*>& paths, NodeRoll& nodes )
{
    vector<NodePath*> cloned;
    for ( NodePath* np : path ) cloned.push_back( new NodePath( np, paths, nodes ) );
    
    for ( int i = 1; i < path.size(); i++ ) for ( PathEdge* pe : path[i-1]->edges_[1] ) if ( pe->edge[1] == path[i] ) new PathEdge( cloned[i-1], cloned[i], pe);
    
    NodePath* ends[2] = { cloned[0], cloned.back() };
    for ( int d : { 0, 1 } ) for ( PathEdge* pe : taken[d][1] ) pe->claim( ends[d], !d );
    for ( int d : { 0, 1 } ) for ( PathEdge* pe : taken[d][0] )
    {
        pe = new PathEdge( ends[d], ends[d]->getTerminus( d ), pe, !d );
        ends[d]->edges_[d].push_back( pe );
        pe->edge[d]->edges_[!d].push_back( pe );
    }
    for ( NodePath* np : cloned ) for ( Node* node : np->path_ ) node->setVerified();
    
    for ( int i = 0; i < path.size(); i++ ) cloned[i]->clonePairs( path[i] );
}

bool NodePath::setCrosses( vector<NodePath*>& paths, NodeRoll& nodes )
{
    sort( paths.begin(), paths.end(), []( NodePath* a, NodePath* b ) { return a->score_ > b->score_; } );
    
    bool crossed = false;
    for ( int i = 0; i < paths.size(); i++ ) if ( paths[i]->edges_[0].size() > 1 && paths[i]->edges_[1].size() > 1 )
    {
        NodePath* np = paths[i];
        vector<NodePath*> path{ paths[i] };
        vector<Branch*> branches[2]{ Branch::create( paths[i], 0 ), Branch::create( paths[i], 1 ) };
        if ( Branch::match( path, branches ) && paths[i]->setCross( branches, paths, nodes ) ) crossed = true;
    }
    if ( crossed && merge( paths ) ) setCrosses( paths, nodes );
    return crossed;
}

bool NodePath::setCross( vector<Branch*> branches[2], vector<NodePath*>& paths, NodeRoll& nodes )
{
    if ( branches[0].size() < 2 || branches[1].size() < 2 || branches[0][0]->getMax() < 2 ) return false;
    
    vector<NodePath*> path = { this };
    vector<PathEdge*> taken[2][2];
    vector<Branch*> test[2]{ branches[0], branches[1] };
    bool unended[2]{ false, false };
    bool culled = false;
    for ( int d : { 0, 1 } ) for ( Branch* b : branches[d] ) b->cull();
    for ( int d : { 0, 1 } ) for ( Branch* b : branches[d] ) if ( b->matched( taken, d ) )
    {
        bool unique = b->matched_[0].first->matched_.size() == 1;
        unordered_set<Branch*> used = { b->matched_[0].first }, cull[2];
        for ( pair<Branch*, int> bp : b->matched_[0].first->matched_ ) used.insert( bp.first );
        if ( !unique && b->matched_[0].first->matched_.back().second < 5 ) continue;
        
        bool good = false, bad = false;
        for ( int i : { 0, 1 } ) for ( Branch* alt : branches[i] ) if ( used.find( alt ) == used.end() )
        {
            if ( alt->matched_.size() == 1 && alt->matched_[0].first->matched_.size() == 1 ) good = true;
            else if ( alt->matched_.empty() ) cull[i].insert( alt );
            else if ( min( b->matched_[0].second, alt->matched_[0].first->matched_.back().second ) < 5 ) bad = true;
            else good = true;
        }
        
        if ( bad ) continue;
        for ( int i : { 0, 1 } ) for ( Branch* b : cull[i] ) if ( culled = true ) b->edge_->downgrade( i );
        if ( !unique ) continue;
        if ( good ) setCloned( path, taken, paths, nodes );
        return true;
    }
    
    return culled;
//    for ( int d : { 0, 1 } ) for ( Branch* b : branches[d] )
//    {
//        if ( b->matched_.size() == 1 && b->matched_[0].first->matched_.size() == 1 ) matched[d].push_back( b );
//        else if ( b->getMax( true ) < branches[d][0]->getMax( true ) / 10 ) trash[d].push_back( b );
//        else unmatched[d].push_back( b );
//    }
//    
//    for ( int d : { 0, 1 } ) for ( Branch* b : trash[d] )
//    {
//        branches[d].erase( remove( branches[d].begin(), branches[d].end(), b ), branches[d].end() );
//        b->edge_->downgrade( d );
//        delete b;
//    }
//    for ( int d : { 0, 1 } ) if ( branches[d].size() != 2 ) return true;
//    
//    
//    vector<PathEdge*> taken[2][2];
//    vector<NodePath*> path = { this };
//    if ( !matched[0].empty() && !unmatchable[0] && !unmatchable[1] )
//    {
//        assert( matched[0][0]->matched( taken, 0 ) );
//        setCloned( path, taken, paths, nodes );
//        return true;
//    }
//    
//    for ( int d : { 0, 1 } ) for ( Branch* b : unmatched[d] ) if ( b->matched( taken, d ) )
//    {
//        setCloned( path, taken, paths, nodes );
//        return true;
//    }
//    
//    return false;
}

//NodePath* NodePath::setCross( vector<Branch*> branches[2], NodeRoll& nodes )
//{
//    int hits = branches[0][0]->hits_, base = min( branches[0][0]->base_, branches[1][0]->base_ );
//    vector<PathEdge*> culls[2];
//    for ( int d : { 0, 1 } ) if ( ( hits = branches[!d][1]->hits_ ) > 2 ) for ( Branch* b : branches[d] )
//    {
//        if ( ( b->match_ ? b->match_->hits_ : hits ) > b->hits_ * 5 ) culls[d].push_back( b->edge_ );
//    }
//    if ( culls[0].empty() && culls[1].empty() && ( hits = branches[0][0]->hits_ ) ) for ( int d : { 0, 1 } ) for ( Branch* b : branches[d] )
//    {
//        if ( ( b->match_ ? b->match_->hits_ : hits ) > ( b->hits_ + b->base_ ) * 5 ) culls[d].push_back( b->edge_ );
//    }
//    
//    if ( ( !culls[0].empty() || !culls[1].empty() ) && ( base + hits < 10 ) ) for ( int d : { 0, 1 } ) culls[d].clear();
//    
//    for ( int d : { 0, 1 } ) for ( PathEdge* pe : culls[d] ) pe->downgrade( d );
//
//    if ( !culls[0].empty() || !culls[1].empty() )
//    {
//        return NULL;
//    }
//
//    vector<PathEdge*> edges[2];
//    for ( int d : { 0, 1 } ) for ( Branch* b : branches[d] ) if ( b->match_ == ( d ? branches[0][0] : branches[0][0]->match_ ) ) edges[d].push_back( b->edge_ );
//
//    bool dupe = edges[0].size() == 1 && edges[1].size() == 1 && ( branches[0][0]->hits_ > max( 2, branches[0][0]->miss_ ) );
//    for ( int d : { 0, 1 } ) for ( Branch* b : branches[d] ) if ( b->hits_ <= b->miss_ * 5 );
//    if ( dupe )
//    {
//        return new NodePath( this, edges[0][0], edges[1][0], nodes );
//    }
//    return NULL;
//}

bool NodePath::setCulled( vector< vector<NodePath*> >& loci, vector<NodePath*>& paths, NodeRoll& nodes )
{
    sort( paths.begin(), paths.end(), []( NodePath* a, NodePath* b ){ return a->cover_ < b->cover_; } );
    
    int unedged = 0, kept = 0, extras = 0;
    unordered_set<NodePath*> goods[2];
    
    float covers = 0;
    int32_t lens = 0;
    
    vector<float> covered;
    for ( vector<NodePath*>& path : loci ) goods[1].insert( path.begin(), path.end() );
    for ( vector<NodePath*>& path : loci )
    {
        for ( PathEdge* pe : path[0]->edges_[0] ) pe->edge[0]->setGood( goods[1], params.readLen+pe->ol, 0 );
        for ( PathEdge* pe : path.back()->edges_[1] ) pe->edge[1]->setGood( goods[1], params.readLen+pe->ol, 1 );
    }
    for ( vector<NodePath*>& path : loci )
    {
        float cover = NodePath::getCoverage( path, params.readLen*2 );
        int32_t len = max( 1, NodePath::getLen( path ) - params.readLen );
        covers += cover * len;
        lens += len;
        covered.push_back( cover );
        
        vector<PathEdge*> culled;
//        for ( NodePath* np : path ) for ( int d : { 0, 1 } ) for ( PathEdge* pe : np->alts_[d] ) if ( !pe->edge[d]->isContinue( params.readLen + pe->ol, d ) )
//        {
//            bool forked = pe->edge[d]->isForked( pe, !d );
//            vector<NodePath*> pathed;
//            if ( pe->isAlt() || pe->edge[d]->setCulled( pathed, goods, pe->edge[d]->size() - pe->ol, cover, forked, d ) ) culled.push_back( pe );
//            else kept++;
//        }
//        for ( PathEdge* pe : culled )
//        {
//            pe->sever();
//        }
        unedged += culled.size();
    }
    
    covers /= lens;
    
    if ( covers > 100 ) for ( NodePath* np : paths ) if ( goods[0].find( np ) == goods[0].end() &&  goods[1].find( np ) == goods[1].end() )
    {
        vector<NodePath*> pathed;
        for ( int d : { 0, 1 } ) if ( !np->isContinue( params.readLen*1.5, d ) && np->setCulled( pathed, 0, covers / 100, d ) )
        {
            extras++;
            goods[0].insert( np );
            break;
        }
    }
    
    if ( !unedged && goods[0].empty() ) return false;
    
    for ( int i = 0; i < paths.size(); i++ ) if ( goods[0].find( paths[i] ) != goods[0].end() )
    {
        int x = 0;
//        paths[i]->destroy( nodes );
//        paths.erase( paths.begin() + i-- );
    }
    Node::updateStates( nodes );
    return true;
}

bool NodePath::setCulled( vector<NodePath*>& pathed, unordered_set<NodePath*> goods[2], int32_t dist, float base, bool forked, bool drxn )
{
    if ( goods[1].find( this ) != goods[1].end() ) return true;
    if ( find( pathed.begin(), pathed.end(), this ) != pathed.end() ) return true;
    
    bool far = dist > params.readLen / 2, bad = false, good = false, advanced = false;
    
    pathed.push_back( this );
    for ( PathEdge* pe : edges_[drxn] ) if ( !pe->bad && goods[1].find( pe->edge[drxn] ) == goods[1].end() )
    {
        bool reforked = pe->edge[drxn]->isForked( pe, !drxn );
        if ( far && reforked ) continue;
        advanced = true;
        if ( !pe->edge[drxn]->setCulled( pathed, goods, dist + pe->edge[drxn]->size() - pe->ol, base, reforked, drxn ) ) good = true;
    }
    float cover = NodePath::getCoverage( pathed );
    if ( advanced && !good ) bad = true;
    else if ( far || !good ) bad = ( NodePath::getCoverage( pathed ) * 12 < base );
    pathed.pop_back();
    
    if ( bad && !forked )
    {
        goods[0].insert( this );
        goods[0].insert( pathed.begin(), pathed.end() );
    }
    
    return bad;
}

bool NodePath::setCulled( vector<NodePath*>& pathed, int32_t len, float cutoff, bool drxn )
{
    len += size();
    bool tested = len > params.readLen * 1.5, bad = false, good = false;
    
    pathed.push_back( this );
    if ( tested ) bad = NodePath::getCoverage( pathed, 1, drxn ) < cutoff;
    else for ( PathEdge* pe : edges_[drxn] ) if ( !pe->edge[drxn]->setCulled( pathed, len - pe->ol, cutoff, drxn ) ) good = true;
    pathed.pop_back();
    return tested ? bad : !good;
}

void NodePath::setEnds( unordered_set<NodePath*>& pathed, vector<NodePath*>& ends, bool drxn )
{
    if ( !pathed.insert( this ).second ) return;
    
    int limits[2]{ 10, 0 }, i = 0;
    vector<PathEdge*> edges[3];
    if ( verified_ || edges_[!drxn].empty() ) for ( PathEdge* pe : edges_[drxn] )
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

void NodePath::setEdges( vector<NodePath*>& paths )
{
    for ( NodePath* np : paths ) for ( PathEdge* pe : np->edges_[1] ) pe->score = 0;
    for ( NodePath* np : paths ) if ( !np->paired_[1].empty() )
    {
        int32_t limit = 0;
        for ( PathPairs* pp : np->paired_[1] ) limit = max( limit, pp->node_[1]->ends_[1] - pp->diff_  - np->ends_[1] );
        limit = limit * 1.1 + 100;
        
        unordered_set<PathEdge*> used;
        vector< unordered_set<PathEdge*> > edged( np->paired_[1].size() );
        for ( PathEdge* pe : np->edges_[1] )
        {
            vector<PathEdge*> path = { pe };
            np->setEdges( pe->edge[1], limit + pe->ol, path, edged, used );
        }
        for ( int i = 0; i < np->paired_[1].size(); i++ )
        {
            for ( PathEdge* pe : edged[i] ) pe->score += np->paired_[1][i]->ids_.size();
        }
    }
}

void NodePath::setEdges( NodePath* np, int32_t limit, vector<PathEdge*>& path, vector< unordered_set<PathEdge*> >& edged, unordered_set<PathEdge*>& used )
{
    for ( int i = 0; i < paired_[1].size(); i++ ) if ( paired_[1][i]->node_[1] == np ) edged[i].insert( path.begin(), path.end() );
    
    limit -= np->size();
    if ( limit > 0 ) for ( PathEdge* pe : np->edges_[1] ) if ( find( path.begin(), path.end(), pe ) == path.end() )
    {
        if ( !used.insert( pe ).second )
        {
            for ( unordered_set<PathEdge*>& e : edged ) if ( e.find( pe ) != e.end() ) e.insert( path.begin(), path.end() );
        }
        else
        {
            path.push_back( pe );
            setEdges( pe->edge[1], limit + pe->ol, path, edged, used );
            path.pop_back();
        }
    }
//    if ( limit > 0 ) for ( PathEdge* pe : np->edges_[1] ) if ( find( path.begin(), path.end(), pe ) == path.end() )
//    {
//        path.push_back( pe );
//        setEdges( pe->edge[1], limit + pe->ol, path, edged );
//        path.pop_back();
//    }
}

void NodePath::setExclusive( vector< unordered_set<NodePath*> >& sets )
{
    for ( int i = 0; i+1 < sets.size(); i++ ) for ( auto it = sets[i].begin(); it != sets[i].end(); )
    {
        bool bad = true;
        for ( int j = i+1; bad && j < sets.size(); j++ ) if ( sets[j].find( *it ) == sets[j].end() ) bad = false;
        if ( bad ) for ( int j = i+1;  j < sets.size(); j++ ) sets[j].erase( *it );
        if ( bad ) it = sets[i].erase( it );
        else it++;
    }
}

void NodePath::setExtends( vector< vector<NodePath*> >& loci, SeedExtend seed[2] )
{
    unordered_set<NodePath*> pathed;
    for ( vector<NodePath*>& path : loci ) pathed.insert( path.begin(), path.end() );
    for ( vector<NodePath*>& path : loci )
    {
        NodePath* ends[2]{ path[0], path.back() };
        Node* node[2]{ NULL, NULL };
        float coverage = NodePath::getCoverage( path, params.readLen*2 );
        for ( int d : { 0, 1 } ) for ( auto& no : ends[d]->offs_ ) if ( no.first->isContinue( d ) && ( node[d] = ends[d]->getTerminus( d ) ) ) break;
        for ( int d : { 0, 1 } ) if ( !node[d] && !ends[d]->getTerminus( d )->getExtendable( 300, d ).empty() ) node[d] = ends[d]->getTerminus( d );
        
        for ( int d : { 0, 1 } ) if ( node[d] ) seed[d].addFork( node[d] );
        for ( int d : { 0, 1 } ) for ( NodePath* np : path ) for ( PathEdge* pe : np->alts_[d] )
        {
            if ( pathed.find( pe->edge[d] ) == pathed.end() ) seed[d].addAlt( pe->node[d] );
        }
    }
}

void NodePath::setFills( vector<NodePath*>& paths )
{
    Nodes block;
    for ( NodePath* np : paths ) for ( Node* node : np->path_ ) block += node;
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) if ( np->edges_[d].empty() ) np->setFill( np->getTerminus( d ), block, np->coords_[d], 300, d );
}

void NodePath::setFill( Node* node, Nodes& block, int32_t dist, int32_t limit, bool drxn )
{
    for ( Edge& e : node->edges_[drxn] ) if ( !block.find( node ) )
    {
        int32_t coord = getCoord( node, e, dist, drxn );
        auto it = offs_.insert( make_pair( e.node, coord ) );
        if ( !it.second )
        {
            if ( drxn ? it.first->second <= coord : coord <= it.first->second ) continue;
            it.first->second = coord;
        }
        if ( limit > 0 ) setFill( e.node, block, coord, limit - e.node->size() + e.ol, drxn );
    }
}

void NodePath::setGood( unordered_set<NodePath*>& goods, int32_t len, bool drxn )
{
    if ( !goods.insert( this ).second || len < 0 ) return;
    len -= size();
    for ( PathEdge* pe : edges_[drxn] ) if ( !pe->bad ) pe->edge[drxn]->setGood( goods, ( pe->score ? params.readLen : len ) + pe->ol, drxn );
}

void NodePath::setMates( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& paths )
{
    Lib* lib;
    int32_t coords[2];
    int added = 0;
    for ( int d : { 0, 1 } ) for ( Node* node : path_ ) for ( NodeMark& nm : node->mp_[d] )
    {
        string seq = bwt.getSequence( nm.id );
        int ol = min( (int)seq.size(), 40 );
        vector<PathPairs*> bests;
        for ( PathPairs* pp : paired_[d] ) if ( mapSeqEnd( seq, pp->node_[d]->seq_, ol, coords, d ) )
        {
            if ( coords[1] - coords[0] > ol ) bests.clear();
            ol = coords[1] - coords[0];
            bests.push_back( pp );
        }
        if ( !bests.empty() && ol < seq.size() ) for ( Node* node : nodes.nodes ) if ( mapSeqEnd( seq, node->seq_, ol+1, coords, d ) )
        {
            bests.clear();
            break;
        }
        for ( PathPairs* pp : bests ) pp->mates_.insert( nm.id );
        if ( !bests.empty() ) added++;
    }
}

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

vector< vector<NodePath*> > NodePath::setPaths( vector<NodePath*>& paths )
{
    vector< vector<NodePath*> > loci;
    sort( paths.begin(), paths.end(), []( NodePath* a, NodePath* b ) { return a->score_ > b->score_; } );
    
    unordered_set<NodePath*> used, ends[2];
    unordered_map<NodePath*, int> dupes;
    for ( NodePath* np : paths ) if ( np->score_ > 1 && !np->blank_ && used.find( np ) == used.end() )
    {
        vector<NodePath*> path = np->setPath( loci );
        loci.push_back( path );
        used.insert( path.begin(), path.end() );
        ends[0].insert( path[0] );
        ends[1].insert( path.back() );
    }
    for ( vector<NodePath*>& path : loci ) for ( NodePath* np : path )
    {
        auto it = dupes.insert( make_pair( np, 1 ) );
        if ( !it.second ) it.first->second++;
    }
    for ( int i = loci.size(); i-- > 0; )
    {
        vector<int> duped;
        for ( NodePath* np : loci[i] ) duped.push_back( dupes[np] );
    }
    return loci;
}

vector<NodePath*> NodePath::setPath( vector< vector<NodePath*> >& loci )
{
    vector<NodePath*> path{ this };
    bool ended[2]{ edges_[0].empty(), edges_[1].empty() };
    NodePath* ends[2]{ this, this };
    while ( !ended[0] || !ended[1] )
    {
        for ( int d : { 0, 1 } ) while ( !ended[d] && !( ended[d] = ends[d]->edges_[d].empty() ) && ends[d]->edges_[d].size() == 1 )
        {
            PathEdge* edge = ends[d]->edges_[d][0];
            if ( !( ended[d] = edge->bad ) ) path.insert( d ? path.end() : path.begin(), ends[d] = edge->edge[d] );
        }
        
        int score[2]{ ended[0] ? 0 : ends[0]->edges_[0][0]->score, ended[1] ? 0 : ends[1]->edges_[1][0]->score };
        int d = ended[0] || ( !ended[1] && score[1] > score[0] );
        if ( ended[d] ) break;
        
        PathEdge* edge = NULL;
        vector<Branch*> branches[2];
        branches[d] = Branch::create( ends[d], d );
        for ( int i = 0; i < path.size() && !edge && branches[d].size() > 1; i++ )
        {
            NodePath* fork = d ? path.end()[-i-1] : path[i],* branch = i+1 == path.size() ? NULL : ( d ? path.end()[-i-2] : path[i+1] );
            if ( fork->edges_[!d].size() < 2 || !fork->edges_[!d][1]->score ) continue;
            
            vector<NodePath*> cross( d ? path.end()-i-1 : path.begin(), d ? path.end() : path.begin()+i+1 );
            
            for ( Branch* b : branches[!d] ) delete b;
            branches[!d] = Branch::create( fork, !d );
            if ( !Branch::match( cross, branches[0], branches[1] ) ) continue;
            for ( int j : { 0, 1 } ) Branch::setUsed( branches[j], loci, j );
            
            Branch* base = NULL;
            if ( branch ) for ( Branch* b : branches[!d] ) if ( base = ( b->edge_->edge[!d] == branch ? b : NULL ) ) break;
            assert( (bool)branch == (bool)base );
            if ( !base ) for ( Branch* b : branches[!d] ) if ( base = ( !b->used_ ? b : NULL ) ) break;
            
            if ( !branch )
            {
                Branch* crossed[2]{ NULL, NULL };
                for ( int d : { 0, 1 } ) for ( Branch* b : branches[d] ) if ( b->match_ && !b->used_ && ( crossed[d] = b ) ) break;
                if ( crossed[0] && crossed[1] && crossed[0] != branches[0][0] && crossed[1] != branches[1][0] )
                {
                    assert( false );
                    crossed[0] = crossed[1] = NULL;
                }
                for ( int d : { 0, 1 } ) if ( crossed[!d] && !crossed[d] ) for ( Branch* b : branches[d] ) if ( b == crossed[!d]->match_ ) crossed[d] = b;
                assert( (bool)crossed[0] == (bool)crossed[1] );
                if ( crossed[!d] ) path.insert( d ? path.begin() : path.end(), ends[!d] = crossed[!d]->edge_->edge[!d] );
                if ( crossed[d] ) edge = crossed[d]->edge_;
                assert( crossed[0] && crossed[1] );
            }
            else for ( int i = 0; i < branches[d].size(); i++ ) if ( base->match_ != branches[d][i] && branches[d][i]->match_ && branches[d][i]->match_ != base )
            {
                delete branches[d][i];
                branches[d].erase( branches[d].begin() + i-- );
            }
        }
        
        if ( !edge && !ended[d] && ( ended[d] = branches[d].empty() ) ) continue;
        if ( !edge && branches[d].size() == 1 ) edge = branches[d][0]->edge_;
        
        if ( !edge )
        {
            Branch::setBase( branches[d], path, d );
            Branch::setUsed( branches[d], loci, d );
            for ( Branch* b : branches[d] ) if ( !edge && b->base_ && !b->edge_->bad && !b->used_ ) edge = b->edge_;
            if ( !edge ) for ( Branch* b : branches[d] ) if ( !edge && b->base_ && !b->edge_->bad ) edge = b->edge_;
            
            if ( ended[d] = ( !edge || edge->edge[d]->blank_ || edge->edge[d]->isEnding( params.readLen / 2, d ) ) ) edge = NULL;
        }
        
        for ( int i : { 0, 1 } ) for ( Branch* b : branches[i] ) delete b;
        if ( !ended[d] ) path.insert( d ? path.end() : path.begin(), ends[d] = edge->edge[d] );
    }
    
    return path;
}

//void NodePath::setPairs( NodePath* np, vector<NodePath*>& path, int32_t diff, bool forked )
//{
//    if ( ends_[0] - diff - np->ends_[1] > params.maxMpMean ) return;
//    if ( np == this && !path.empty() ) addPair( this, diff );
//    if ( count( path.begin(), path.end(), this ) ) return;
//    path.push_back( this );
//    
//    bool failed = false;
//    if ( forked ) failed = !np->addPair( this, diff );
//    
//    if ( edges_[1].size() > 1 ) forked = true;
//    if ( edges_[1].size() == 1 && edges_[1][0]->edge[1]->edges_[0].size() > 1 ) forked = true;
//    if ( !failed ) for ( PathEdge* pe : edges_[1] ) pe->edge[1]->setPairs( np, path, diff + pe->diff, forked );
//    path.pop_back();
//}

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

void NodePath::setReachable( unordered_map<NodePath*, int32_t>& reached, int32_t dist, int32_t limit, bool drxn )
{
    auto it = reached.insert( make_pair( this, dist ) );
    if ( !it.second )
    {
        if ( it.first->second <= dist ) return;
        it.first->second = dist;
    }
    
    if ( dist < limit ) for ( PathEdge* pe : edges_[drxn] ) pe->edge[drxn]->setReachable( reached, dist - pe->ol + pe->edge[drxn]->size(), limit, drxn );
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
    reached.insert( this );
    for ( PathEdge* pe : edges_[drxn] ) if ( reached.insert( pe->edge[drxn] ).second ) pe->edge[drxn]->setReachable( reached, drxn );
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
    for ( NodePath* np : paths ) for ( Node* node : np->path_ ) if ( node->verified_ ) node->remap( bwt );
    
    vector<PathPair*> reads;
    Lib* lib;
    Coords* coords[2];
    unordered_set<ReadId> used;
    for ( NodePath* np : paths )
    {
        bool forked = false;
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

bool NodePath::setUnforked( vector<NodePath*>& paths, NodeRoll& nodes )
{
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } )
    {
        np->setUnforked( d );
    }
    return false;
}

bool NodePath::setUnforked( bool drxn )
{
    if ( edges_[drxn].size() < 2 ) return false;
    int base = edges_[drxn][0]->score;
    for ( int i = 1; i < edges_[drxn].size(); i++ ) if ( edges_[drxn][i]->score*5 < edges_[drxn][0]->score-5 )
    {
        vector<NodePath*> path = { this, edges_[drxn][i]->edge[drxn] };
        for ( NodePath* np = path.back(); np->edges_[0].size() == 1 && np->edges_[1].size() == 1; )
        {
            np = np->edges_[drxn][0]->edge[drxn];
            path.push_back( np );
        }
        int edge = edges_[drxn][i]->score;
        if ( path.back()->edges_[!drxn].size() < 2 ) continue;
        int x = 0;
    }
    return false;
}

void NodePath::setScore()
{
    for ( const pair<Node*, int32_t>& no : offs_ )
    {
        bool pathed[2]{ isPathed( no.first ), false };
        for ( auto& np : no.first->hits_.pairs[1] )
        {
            int32_t* off, cutoff = np.second.estimate() + np.second.margin();
            for ( PathPairs* pp : paired_[1] ) if ( off = pp->node_[1]->get( np.first ) )
            {
                int32_t dist = *off - no.second - pp->diff_;
                if ( dist < cutoff )
                {
                    pathed[1] = pp->node_[1]->isPathed( np.first );
                    pp->score_ += np.second.count;
                    for ( int d : { 0, 1 } ) if ( pathed[d] ) pp->node_[d]->score_ += np.second.count;
                    for ( int i : { 0, 1 } ) for ( pair<ReadId, int32_t>& id : np.second.pairs[i] ) pp->ids_.insert( id.first );
                }
            }
            if ( pathed[0] && isPathed( np.first ) ) score_ += np.second.count;
        }
    }
    for ( Node* node : path_ ) score_ += node->hits_.count;
}

void NodePath::setTarget( unordered_set<NodePath*>& reached, int32_t dist, int32_t limit, bool drxn, bool ignore )
{
    if ( !ignore && !reached.insert( this ).second ) return;
    if ( !ignore ) dist += ends_[1] - ends_[0];
    if ( dist < limit ) for ( PathEdge* pe : edges_[drxn] ) pe->edge[drxn]->setTarget( reached, dist - pe->ol, limit, drxn );
}

//bool NodePath::setUnbad( vector<NodePath*>& paths, NodeRoll& nodes )
//{
//    int unbads = 0;
//    if ( !params.rna ) return false;
//    for ( NodePath* np : paths ) if ( np->path_[0]->bad_ ) for ( int d : { 0, 1 } ) if ( !np->edges_[d].empty() )
//    {
//        PathEdge::sort( np->edges_[d] );
//    }
//    if ( unbads ) Node::updateStates( nodes );
//    return unbads;
//}

int32_t NodePath::size()
{
    return ends_[1] - ends_[0];
}

//void NodePath::sweep( NodeRoll& nodes, bool drxn )
//{
//    vector<PathEdge*> trash;
//    for ( int d : { 0, 1 } ) trash.insert( trash.end(), breaks_[d].begin(), breaks_[d].end() );
//    
//    Node* node[2]{ drxn ? path_[0] : path_.back(), NULL };
//    for ( PathEdge* pe : edges_[!drxn] ) if ( ( node[1] = drxn ? pe->edge[!drxn]->path_.back() : pe->edge[!drxn]->path_[0] ) && node[1]->bad_ )
//    {
//        Node* clone = new Node( node[0], nodes, node[0]->drxn_, true );
//        clone->addEdge( node[1], pe->ol, !drxn, false, pe->leap );
//        trash.push_back( pe );
//    }
//    for ( PathEdge* pe : trash ) PathEdge::sever( pe->edge[0], pe->edge[1] );
//    assert( edges_[!drxn].size() > 0 );
//    if ( edges_[0].size() == 1 && edges_[1].size() == 1 ) edges_[drxn][0]->edge[drxn]->sweep( nodes, drxn );
//}

bool NodePath::unfork( bool drxn )
{
    bool unforked = false;
    for ( int i = 1; i < edges_[drxn].size(); i++ ) if ( edges_[drxn][i]->score + max( 0, edges_[drxn][i]->score-1 ) * 5 < ( edges_[drxn][0]->score-2 ) / 2 )
    {
        edges_[drxn][i--]->downgrade( drxn );
    }
    
    return unforked;
}

HapBranch::HapBranch( NodePath* branch, int32_t dist, bool drxn )
: branch_( branch )
{
    branch->setReachable( dists_, dist, params.maxMpMean*1.3+200, drxn );
}

int32_t* HapBranch::get( NodePath* np )
{
    auto it = dists_.find( np );
    return it != dists_.end() ? &it->second : NULL;
}

HapFork::HapFork( NodePath* fork, Haplo* base, vector<HapFork*>& forks, bool drxn )
: fork_( fork ), base_{ base }
{
    for ( PathEdge* pe : fork->edges_[drxn] ) branches_.push_back( new HapBranch( pe->edge[drxn], pe->edge[drxn]->size()-pe->ol, drxn ) );
    for ( HapBranch* hp : branches_ ) for ( pair<NodePath*, int32_t> npd : hp->dists_) reached_.insert( npd.first );
    
    for ( HapFork* hf : forks ) if ( hf->fork_ )
    {
        if ( reached_.find( hf->fork_ ) != reached_.end() ) hf->bck_.push_back( this );
        if ( hf->reached_.find( fork_ ) != hf->reached_.end() ) bck_.push_back( hf );
    }
    forks.push_back( this );
    
    //NYI
    assert( reached_.find( fork ) == reached_.end() );
}

HapFork::~HapFork()
{
    for ( HapBranch* hp : branches_ ) delete hp;
    for ( pair<HapFork*, vector<HapScore*> >& c : cross_ )
    {
        for ( int i = 0; i < c.first->cross_.size(); i++ ) if ( c.first->cross_[i].first == this ) c.first->cross_.erase( c.first->cross_.begin() +i-- );
        for ( HapScore* hs : c.second ) delete hs;
    }
}

vector<HapScore*> HapFork::cross( HapFork* alt, bool drxn )
{
    for ( pair<HapFork*, vector<HapScore*> >& c : cross_ ) if ( c.first == alt ) return c.second;
    vector<HapScore*> scores;
    for ( HapBranch* f : branches_ ) for ( HapBranch* b : branches_ ) scores.push_back( new HapScore( drxn ? b : f, drxn ? f : b ) );
    cross_.push_back( make_pair( alt, scores ) );
    alt->cross_.push_back( make_pair( this, scores ) );
    return scores;
}

bool HapFork::cull( vector<HapFork*>& forks, bool drxn )
{
    for ( HapFork* hf : forks ) hf->bck_.erase( remove( hf->bck_.begin(), hf->bck_.end(), this ), hf->bck_.end() );
    forks.erase( remove( forks.begin(), forks.end(), this ), forks.end() );
    delete this;
    return true;
}

bool HapFork::resolve( vector<HapFork*> forks[2], bool drxn )
{
    assert( base_.size() <= 2 );
    if ( resolveBranch( forks[drxn], drxn ) ) return true;
        
    for ( int i = 0; i+1 < base_.size(); i++ ) for ( int j = i+1; j < base_.size(); j++ ) if ( base_[i]->fork_[!drxn] == base_[j]->fork_[!drxn] )
    {
        HapFork* alts[2]{ drxn ?  base_[i]->fork_[0] : this, drxn ? this : base_[i]->fork_[1] };
        if ( alts[!drxn]->resolveBranch( forks[!drxn], !drxn ) ) return false;
        vector<HapScore*> scores = cross( base_[i]->fork_[!drxn], drxn );
        for ( HapScore* hs : scores )
        {
            hs->pref_ = hs->ids_;
            for ( HapScore* alt : scores ) if ( alt != hs ) for ( ReadId id : alt->ids_ ) hs->pref_.erase( id );
            assert( hs->pref_.empty() );
        }
        
        assert( alts[0]->base_.size() == 2 && alts[1]->base_.size() == 2 );
        for ( int d : { 0, 1 } ) alts[d]->fork_ = NULL;
        for ( int i = 0; i < base_.size(); i++ )
        {
            for ( int d : { 0, 1 } ) base_[i]->unfork( d );
            for ( int d : { 0, 1 } ) base_[i]->setFork( alts[d]->branches_[i]->branch_, forks[d], d );
        }
        for ( int d : { 0, 1 } ) alts[d]->cull( forks[d], d );
        return true;
    }
//    for ( Haplo* h : base_ ) h->setFork( forks, drxn );
    int x = 0;
    return false;
}

bool HapFork::resolveBranch( vector<HapFork*>& forks, bool drxn )
{
    vector< pair<HapBranch*, unordered_set<ReadId> > > ids;
    for ( HapBranch* hb : branches_ ) ids.push_back( make_pair( hb, unordered_set<ReadId>{} ) );
    for ( Haplo* h : base_ ) h->setBranches( drxn );
    for ( Haplo* h : base_ ) for ( HapScore* hs : h->branch_[drxn] ) for ( pair<HapBranch*, unordered_set<ReadId> >& branch : ids )
    {
        if ( hs->branch_[drxn] != branch.first ) continue;
        branch.second.insert( hs->ids_.begin(), hs->ids_.end() );
        break;
    }
    sort( ids.begin(), ids.end(), []( pair<HapBranch*, unordered_set<ReadId> >& a, pair<HapBranch*, unordered_set<ReadId> >& b ){ return a.second.size() > b.second.size(); } );
    
    bool good = false, bad = false;
    for ( int i = 0; i < ids.size(); i++ ) if ( ids[i].second.empty() )
    {
        if ( !i && ( bad = true ) ) break;
    }
    if ( bad ) return cull( forks, drxn );
    
    for ( Haplo* h : base_ ) for ( HapScore* hs : h->branch_[drxn] ) hs->unique_ = hs->miss_ = 0;
    for ( Haplo* h : base_ )
    {
        unordered_set<ReadId> used;
        for ( Haplo* alt : base_ ) if ( h != alt ) for ( HapScore* hs : alt->branch_[drxn] ) used.insert( hs->ids_.begin(), hs->ids_.end() );
        for ( HapScore* hs : h->branch_[drxn] ) hs->unique_ = 0;
        for ( HapScore* hs : h->branch_[drxn] ) for ( ReadId id : hs->pref_ ) if ( used.find( id ) == used.end() ) hs->unique_++;
        for ( HapScore* hs : h->branch_[drxn] ) if ( hs->unique_ ) for ( Haplo* alt : base_ ) if ( h != alt )
        {
            for ( HapScore* hsa : alt->branch_[drxn] ) if ( hsa->branch_[drxn] == hs->branch_[drxn] ) hsa->miss_ += hs->unique_;
        }
    }
    
    
    assert( branches_.size() == 2 );
    
    for ( Haplo* h : base_ ) sort( h->branch_[drxn].begin(), h->branch_[drxn].end(), []( HapScore* a, HapScore* b )
    { return a->unique_ == b->unique_ ? a->miss_ < b->miss_ : a->unique_ > b->unique_; } );
    
    for ( Haplo* h : base_ ) if ( h->branch_[drxn][0]->unique_ && !h->branch_[drxn][1]->unique_ ) good = true;
    if ( good )
    {
        fork_ = NULL;
        for ( Haplo* h : base_ ) h->unfork( drxn );
        for ( Haplo* h : base_ ) h->setFork( h->branch_[drxn][0]->branch_[drxn]->branch_, forks, drxn );
        return cull( forks, drxn );
    }
    for ( Haplo* h : base_ ) for ( HapScore* hs : h->branch_[drxn] ) assert( !hs->unique_ );
    int x = 0;
    return false;
}

HapScore::HapScore( Haplo* h, HapBranch* branch, bool drxn )
{
    branch_[drxn] = branch;
    int32_t* dist;
    for ( NodePath* np : h->path_ ) for ( PathPairs* pp : np->paired_[drxn] ) if ( dist = branch->get( pp->node_[drxn] ) )
    {
        ids_.insert( pp->ids_.begin(), pp->ids_.end() );
        mates_.insert( pp->mates_.begin(), pp->mates_.end() );
    }
}

HapScore::HapScore( HapBranch* l, HapBranch* r )
{
    branch_[0] = l;
    branch_[1] = r;
    int32_t* dist;
    for ( pair<NodePath*, int32_t> d : l->dists_ ) for ( PathPairs* pp : d.first->paired_[1] ) if ( dist = r->get( pp->node_[1] ) )
    {
        ids_.insert( pp->ids_.begin(), pp->ids_.end() );
        mates_.insert( pp->mates_.begin(), pp->mates_.end() );
    }
}

Haplo::Haplo( NodePath* seed )
: path_{ seed }
{ }

void Haplo::add( NodePath* np, bool drxn )
{
    assert( find( path_.begin(), path_.end(), np ) == path_.end() );
    path_.insert( drxn ? path_.end() : path_.begin(), np );
    for ( int d : { 0, 1 } ) assert( branch_[d].empty() );
}

void Haplo::setBranches( bool drxn )
{
    if ( !branch_[drxn].empty() || !fork_[drxn] ) return;
    for ( HapBranch* hb : fork_[drxn]->branches_ ) branch_[drxn].push_back( new HapScore( this, hb, drxn ) );
    for ( HapScore* hs : branch_[drxn] ) hs->pref_ = hs->ids_;
    for ( HapScore* hs : branch_[drxn] ) for ( HapScore* alt : branch_[drxn] ) if ( hs != alt ) for ( ReadId id : alt->ids_ ) hs->pref_.erase( id );
    sort( branch_[drxn].begin(), branch_[drxn].end(), []( HapScore* a, HapScore* b ){ return a->ids_.size() > b->ids_.size(); } );
}

void Haplo::setFork( NodePath* fork, vector<HapFork*>& forks, bool drxn )
{
    fork_[drxn] = NULL;
    if ( !fork ) fork = drxn ? path_.back() : path_[0];
    else add( fork, drxn );
    while ( !fork->edges_[drxn].empty() )
    {
        PathEdge::sort( fork->edges_[drxn] );
        if ( fork->edges_[drxn].size() > 1 && ( !fork->edges_[drxn][0]->score || fork->edges_[drxn][1]->score ) ) break;
        add( fork = fork->edges_[drxn][0]->edge[drxn], drxn );
    }
    if ( fork->edges_[drxn].empty() ) return;
    for ( HapFork* hp : forks ) if ( hp->fork_ == fork && ( fork_[drxn] = hp ) )
    {
        hp->base_.push_back( this );
        return;
    }
    fork_[drxn] = new HapFork( fork, this, forks, drxn );
}

void Haplo::unfork( bool drxn )
{
    for ( HapScore* hs : branch_[drxn] ) delete hs;
    branch_[drxn].clear();
    
}

vector<Haplo*> Haplo::create( Querier& bwt, vector<NodePath*>& paths, vector<NodePath*>& seeds, NodeRoll& nodes )
{
    vector<Haplo*> haps;
    for ( NodePath* seed : seeds )
    {
        haps.push_back( new Haplo( seed ) );
        if ( abs( seed->cover_ - params.cover ) < abs( seed->cover_ - ( params.cover / 2 ) ) ) haps.push_back( new Haplo( seed ) );
    }
    vector<HapFork*> forks[2];
    for ( int d : { 0, 1 } ) for ( Haplo* h : haps ) h->setFork( NULL, forks[d], d );
    for ( int advanced = 1; advanced-- > 0; )
    {
        for ( int d : { 0, 1 } ) for ( int i = 0; i < forks[d].size(); i++ ) if ( forks[d][i]->bck_.empty() && forks[d][i]->resolve( forks, d ) )
        {
            advanced = true;
            i--;
        }
        assert( advanced );
    }
    
    return haps;
}

//PathEdgeScore::PathEdgeScore( NodePath* np, int32_t diff, int32_t limit, bool drxn )
//: node( np ), diff( diff ), hits( 0 ), miss( 0 ), unique( 0 )
//{
//    np->setReachable( fwd, diff, limit, drxn );
//}
//
////bool PathEdgeScore::add( PathPairing* pp, vector< pair<NodesPath*, int32_t> >& tar, vector<int32_t>& diffs, bool missed, bool drxn )
////{
//////    auto it = tar.find( pp->node[!drxn] );
//////    if ( it == tar.end() ) return false;
//////    ( missed ? miss : score ) += pp->hits( diffs, it->second );
////    return true;
////}
//
//void PathEdgeScore::add( NodePath* q, vector<int32_t>& qDiffs, unordered_set<NodePath*>& alts, vector< pair<NodePath*, int32_t> >& tar, bool drxn )
//{
//    bool used = false;
//    for ( auto& nd : tar ) if ( nd.first == q ) used = true;
//    
//    int scores[2]{0};
//    
//    for ( PathPairing* pp : q->pairs_[!drxn] )
//    {
//        vector<int32_t> tDiffs;
//        for ( auto& nd : tar ) if ( nd.first == pp->node[!drxn] ) tDiffs.push_back( nd.second );
//        
//        // This pair is not part of the target
//        if ( tDiffs.empty() )
//        {
//            // Score if this is a valid alt
//            if ( alts.find( pp->node[!drxn] ) != alts.end() ) scores[0] += pp->score;
//            continue;
//        }
//        
//        int maxScore = 0, minScore = pp->score, sum = 0;
//        for ( int32_t td : tDiffs )
//        {
//            bool scored = false;
//            int best = 0;
//            for ( int32_t qd : qDiffs ) for ( int i = 0; i < pp->diffs.size(); i++ ) if ( td + qd == pp->diffs[i] )
//            {
//                best = max( best, pp->score - pp->missed[i] );
//                scored = true;
//            }
//            assert( scored );
//            maxScore = max( maxScore, best );
//            minScore = min( minScore, best );
//            sum += best;
//        }
//        
//        assert( sum <= pp->score );
//        assert( minScore <= pp->score );
//        scores[1] += used ? min( sum, minScore ) : min( sum, pp->score );
//    }
//    
//    miss += scores[0];
//    hits += scores[1];
//    if ( !scores[0] || scores[1] > scores[0]*5 ) unique += scores[1] - scores[0];
//}
//
//vector<int32_t>* PathEdgeScore::get( NodePath* np )
//{
//    auto it = fwd.find( np );
//    return it != fwd.end() ? &it->second : NULL;
//}
//
//void PathEdgeScore::getKeep( unordered_set<NodePath*>& keep, int32_t limit, bool drxn )
//{
//    for ( auto& nd : fwd ) if ( drxn ? nd.first->ends_[0] - nd.second.back() < limit 
//                                     : limit < nd.first->ends_[1] + nd.second.back() ) keep.insert( nd.first );
//}
//
//void PathEdgeScore::setKeep( unordered_set<NodePath*>& keep )
//{
//    for ( auto it = fwd.begin(); it != fwd.end(); )
//    {
//        if ( it->first != node && keep.find( it->first ) == keep.end() ) it = fwd.erase( it );
//        else it++;
//    }
//}

//PathMapping::PathMapping( vector<NodePath*>& path )
//: path_( path ), diff_( 0 )
//{
//    for ( NodePath* np : path_ ) diffs_[1].push_back( make_pair( np, 0 ) );
//    for ( NodePath* np : path_ ) diffs_[0].insert( diffs_[0].begin(), make_pair( np, 0 ) );
//    for ( int d : { 0, 1 } ) edge( d );
//    while ( advance() || crossroad() );
//}
//
//void PathMapping::add( NodePath* np, int32_t diff, bool drxn )
//{
//    diff_ += diff;
//    if ( diff ) for ( auto& nd : diffs_[drxn] ) nd.second += diff;
//    if ( diff ) for ( auto& nd : alts_[!drxn] ) nd.second += diff;
//    for ( PathEdge* pe : np->edges_[!drxn] ) if ( pe->edge[!drxn] != getFork( drxn ) ) alts_[!drxn].push_back( make_pair( pe->edge[!drxn], pe->diff ) );
//    path_.insert( drxn ? path_.end() : path_.begin(), np );
//    diffs_[drxn].push_back( make_pair( np, 0 ) );
//    diffs_[!drxn].insert( diffs_[!drxn].begin(), make_pair( np, diff_ ) );
//    exts_[drxn].clear();
//    edge( drxn );
//}
//
//bool PathMapping::advance()
//{
//    int scores[2][2];
//    for ( int d : { 0, 1 } )
//    {
//        while ( exts_[d].size() == 1 ) add( exts_[d][0].node, exts_[d][0].diff, d );
//        sort( exts_[d].begin(), exts_[d].end(), []( PathEdgeScore& a, PathEdgeScore& b ){ return a.hits+a.unique > b.hits+b.unique; } );
//        scores[d][0] = exts_[d].empty() ? 0 : exts_[d][0].hits + exts_[d][0].unique;
//        scores[d][1] = exts_[d].size() < 2 ? 0 : exts_[d][1].hits + exts_[d][1].unique;
//    }
//    
//    if ( !scores[0][0] && !scores[1][0] ) return false;
//    int d = scores[1][0] - scores[1][1] > scores[0][0] - scores[0][1];
//    assert( !exts_[d].empty() );
//    for ( int i = 1; i < exts_[d].size(); i++ ) alts_[d].push_back( make_pair( exts_[d][i].node, diff_ + exts_[d][i].diff ) );
//    add( exts_[d][0].node, exts_[d][0].diff, d );
//    return true;
//}
//
//bool PathMapping::crossroad()
//{
//    bool pathable[2]{ false, false };
//    for ( int d : { 0, 1 } ) for ( PathEdgeScore& pes : exts_[d] ) if ( pes.node->verified_ ) pathable[d] = true;
//    if ( !pathable[0] || pathable[1] ) return false;
//    assert( false );
//    return true;
//}
//
//void PathMapping::edge( bool drxn )
//{
//    if ( !getFork( drxn )->verified_ ) return;
//    
//    // instantiate forward edges with nodes up 1000 distance away
//    for ( PathEdge* pe : getFork( drxn )->edges_[drxn] ) exts_[drxn].push_back( PathEdgeScore( pe->edge[drxn], pe->diff, 1000, drxn ) );
//    
//    unordered_set<NodePath*> base, keep, alts;
//    
//    // Ensure that foward nodes are not unevenly included due minor differences in distance
//    for ( PathEdgeScore& pes : exts_[drxn] ) pes.getKeep( keep, getFork( drxn )->ends_[drxn] + ( drxn ? 400 : -400 ), drxn );
//    for ( NodePath* np : path_ ) keep.erase( np );
//    for ( PathEdgeScore& pes : exts_[drxn] ) pes.setKeep( keep );
//    
//    // Block for the alts
//    for ( NodePath* np : path_ ) base.insert( np );
//    for ( int d : { 0, 1 } ) for ( PathEdge* pe : getFork( d )->edges_[d] ) pe->edge[d]->setReachable( base, d );
//    
//    for ( auto& nd : alts_[!drxn] ) nd.first->setReachable( alts, base, 500, !drxn );
//    
//    vector<int32_t>* diffs;
//    for ( int i = path_.size(); i-- > 0; )
//    {
//        NodePath* np = diffs_[drxn][i].first;
//        int32_t diff = diffs_[drxn][i].second;
//        int32_t dist = np->getLen( diffs_[drxn].back().first, diff, drxn );
//        for ( PathPairing* pp : np->pairs_[drxn] ) for ( PathEdgeScore& pes : exts_[drxn] ) if ( diffs = pes.get( pp->node[drxn] ) )
//        {
//            int best = 0;
//            for ( int32_t d : *diffs ) for ( int j = 0; j < pp->diffs.size(); j++ ) if ( d + diff == pp->diffs[j] ) best = max( best, pp->score - pp->missed[j] );
//            pes.hits += best;
//        }
//        if ( i > 1 && exts_[drxn].size() > 1 && np->edges_[!drxn].size() > 1 )
//        {
//            int x = 0;
//        }
//    }
//    int x = 0;
////    for ( PathEdgeScore& pes : edges_[drxn] )
////    {
////        vector< pair<NodesPath*, int32_t> > rev = alts;
////        pes.node->setReachable( rev, base, pes.diff, !drxn, false );
////        
////        for ( auto& nd : pes.fwd ) for ( PathPairing* pp : nd.first->pairs_[!drxn] )
////        {
////            if ( pes.add( pp, diffs_[drxn], nd.second, false, drxn ) ) continue;
////            if ( pes.add( pp, rev, nd.second, true, drxn ) ) continue;
////        }
////    }
//}
//
//NodePath* PathMapping::getFork( bool drxn )
//{
//    assert( !path_.empty() );
//    return drxn ? path_.back() : path_[0];
//}
//
//vector<NodePath*> PathMapping::map( NodePath* seed )
//{
//    vector<NodePath*> path = { seed };
//    PathMapping pm( path );
//    return path;
//}

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
