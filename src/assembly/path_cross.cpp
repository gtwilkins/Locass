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

#include "path_cross.h"
#include "locus_port.h"
#include <limits>
#include <algorithm>

CrossPair::CrossPair( CrossBranch* l, CrossBranch* r )
: score_( 0 ), miss_( 0 )
{
    branch[0] = l;
    branch[1] = r;
    l->pairs_.push_back( this );
    r->pairs_.push_back( this );
}

CrossPair::~CrossPair()
{
    for ( int d : { 0, 1 } ) branch[d]->pairs_.erase( remove( branch[d]->pairs_.begin(), branch[d]->pairs_.end(), this ), branch[d]->pairs_.end() );
}

CrossBranch::CrossBranch( CrossBranch* fork )
: fork_( fork ), loop_( 0 ), cover( 0 )
{
    if ( fork ) path_ = fork->path_;
    if ( fork ) diffs_ = fork->diffs_;
}

CrossBranch::CrossBranch( NodePath* fork, NodePath* block, int32_t off, bool drxn )
: fork_( NULL ), loop_( 0 ), cover( 0 )
{
    for ( int d : { 0, 1 } ) ends_[d] = fork->ends_[d]+off;
    
    int32_t limit = params.maxMpMean * 1.3 + 300;
    if ( block )
    {
        path_ = { fork };
        unordered_set<NodePath*> used;
        for ( PathEdge* pe : fork->edges_[drxn] ) if ( pe->edge[drxn] != block )
        {
            fill( pe->edge[drxn], used, off + ( drxn ? -pe->diff : pe->diff ), limit, drxn );
        }
    }
    else for ( NodePath* np = fork; np; )
    {
        path_.insert( drxn ? path_.end() : path_.begin(), np );
        assert( diffs_.insert( make_pair( np, off ) ).second );
        limit -= np->size();
        ends_[drxn] = np->ends_[drxn] + off;
        if ( np->edges_[drxn].size() == 1 )
        {
            limit += np->edges_[drxn][0]->ol;
            off += np->edges_[drxn][0]->diff * ( drxn ? -1 : 1 );
            np = np->edges_[drxn][0]->edge[drxn];
        }
        else
        {
            unordered_set<NodePath*> used;
            for ( PathEdge* pe : np->edges_[drxn] ) fill( pe->edge[drxn], used, off + ( drxn ? -pe->diff : pe->diff ), limit+pe->ol, drxn );
            np = NULL;
        }
    }
}

void CrossBranch::fill( NodePath* np, unordered_set<NodePath*>& used, int32_t diff, int32_t limit, bool drxn )
{
    if ( np == ( drxn ? path_[0] : path_.back() ) )
    {
        int32_t looped = abs( np->ends_[!drxn] + diff - ends_[!drxn] );
        if ( !loop_ || looped < loop_ ) loop_ = looped;
        return;
    }
    if ( !used.insert( np ).second ) return;
    
    auto it = diffs_.insert( make_pair( np, diff ) );
    bool redundant = false;
    if ( !it.second )
    {
        if ( drxn ? diff >= it.first->second : diff <= it.first->second ) redundant = true;
        else it.first->second = diff;
    }
    
    limit -= np->size();
    if ( limit > 0 && !redundant ) for ( PathEdge* pe : np->edges_[drxn] ) fill( pe->edge[drxn], used, diff + ( drxn ? -pe->diff : pe->diff ), limit+pe->ol, drxn );
    used.erase( np );
}

int32_t* CrossBranch::get( NodePath* np )
{
    auto it = diffs_.find( np );
    return it != diffs_.end() ? &it->second : NULL;
}

vector<PathEdge*> CrossBranch::getEdges( bool drxn )
{
    int limits[2]{ 0, 0 };
    vector<PathEdge*> edges = ( drxn ? path_.back() : path_[0] )->edges_[drxn];
    bool culled = false;
    for ( int i = 0; i < edges.size(); i++ )
    {
        limits[1] = max( limits[1], edges[i]->score );
        if ( edges[i]->score > 2 || ( edges[i]->score && edges[i]->edge[drxn]->verified_ ) ) continue;
        limits[0] = max( limits[0], edges[i]->score + 1 );
        edges.erase( edges.begin() + i-- );
        culled = true;
    }
    if ( culled && ( limits[1] < limits[0] * 10 ) ) edges.clear();
    return edges;
}

PathCross::PathCross( vector<NodePath*>& cross, int32_t off, bool crossed )
: cross_( cross ), crossed_( crossed )
{
    ends_[0] = cross[0]->ends_[0];
    ends_[1] = cross.back()->ends_[1] + off;
    
//    vector<Node*> path;
//    for ( NodePath* np : cross ) for ( Node* node : np->path_ ) path.push_back( node );
//    cover_ = Node::getCoverage( path );
    
    if ( crossed )
    {
        for ( PathEdge* pe : cross[0]->edges_[0] ) extend( new CrossBranch( NULL ), pe->edge[0], pe->diff, true, 0 );
        for ( PathEdge* pe : cross.back()->edges_[1] ) extend( new CrossBranch( NULL ), pe->edge[1], -pe->diff, true, 1 );
    }
    else
    {
        assert( cross.size() == 2 );
        branch_[0].push_back( new CrossBranch( cross[0], NULL, 0, 0 ) );
        branch_[0].push_back( new CrossBranch( cross[1], cross[0], 0, 0 ) );
        branch_[1].push_back( new CrossBranch( cross[1], NULL, off, 1 ) );
        branch_[1].push_back( new CrossBranch( cross[0], cross[1], off, 1 ) );
    }
    
    for ( int d : { 0, 1 } )
    {
        int32_t* diff = NULL;
        vector<NodePath*> shares;
        for ( auto& nd : branch_[d][0]->diffs_ )
        {
            bool shared = true;
            for ( int i = 1; shared && i < branch_[d].size(); i++ )
            {
                if ( !( diff = branch_[d][i]->get( nd.first ) ) || abs( nd.second - *diff ) > 500 ) shared = false;
                if ( shared ) assert( diff );
            }
            if ( shared ) shares.push_back( nd.first );
        }
        for ( NodePath* np : shares ) for ( CrossBranch* cb : branch_[d] ) cb->diffs_.erase( np );
    }
    
    unordered_set<ReadId> used;
    for ( CrossBranch* l : branch_[0] ) for ( CrossBranch* r : branch_[1] ) new CrossPair( l, r );
    for ( CrossBranch* l : branch_[0] )
    {
        for ( const pair<NodePath*, int32_t>& nd : l->diffs_ ) for ( PathPairs* pp : nd.first->paired_[1] )
        {
            for ( PathPair* p : pp->shared_ ) if ( used.insert( p->id_ ).second ) match( p->marks_, p->id_, p->dist_ );
            for ( PathUnique& pu : pp->uniques_ ) if ( used.insert( pu.id_ ).second )
            {
                vector< pair<NodePath*, int32_t> > tars[2]{ { make_pair( pp->node_[0], pu.l_ ) }, { make_pair( pp->node_[1], pu.r_ ) } };
                match( tars, pu.id_, pu.dist_ );
            }
        }
    }
    for ( CrossBranch* l : branch_[0] ) for ( CrossPair* cp : l->pairs_ )
    {
        cp->score_ = cp->unique_.size();
        for ( int d : { 0, 1 } ) for ( ReadId id : cp->branch[d]->unique_ ) if ( cp->unique_.find( id ) == cp->unique_.end() && cp->hits_.find( id ) == cp->hits_.end() ) cp->miss_++;
    }
    for ( int d : { 0, 1 } ) for ( CrossBranch* fork : forks_[d] )
    {
        vector<CrossBranch*> pile[2];
        for ( CrossBranch* branch : branch_[d] )
        {
            bool forked = false;
            for ( CrossBranch* cb = branch; cb->fork_ && !( forked = cb->fork_ == fork ); ) cb = cb->fork_;
            pile[forked].push_back( branch );
        }
        unordered_set<ReadId> used;
        for ( CrossBranch* cb : pile[0] ) used.insert( cb->hits_.begin(), cb->hits_.end() );
        for ( CrossBranch* cb : pile[1] ) for ( ReadId id : cb->hits_ )
        {
            fork->hits_.insert( id );
            if ( used.find( id ) == used.end() ) fork->unique_.insert( id );
        }
    }
}

PathCross::~PathCross()
{
    for ( int d : { 0, 1 } ) for ( CrossBranch* cb : branch_[d] ) delete cb;
}

bool PathCross::claim( vector<NodePath*>& paths, NodeRoll& nodes )
{
    if ( !crossed_ )
    {
        assert( cross_.size() == 2 && branch_[0].size() == 2 && branch_[1].size() == 2 && branch_[0][0]->pairs_.size() == 2 && branch_[1][0]->pairs_.size() == 2 );
        bool good = false;
        for ( PathEdge* pe : cross_[0]->edges_[1] ) if ( pe->edge[1] == cross_[1] ) good = true;
        if ( !good || cross_[0]->edges_[1].size() < 2 || cross_[1]->edges_[0].size() < 2 ) return false;
        int score[2]{0};
        for ( int d : { 0, 1 } ) if ( ( score[d] = branch_[d][0]->pairs_[1]->score_ ) < branch_[0][0]->pairs_[0]->score_ * 4 + 4 ) return false;
        return PathEdge::sever( cross_[0], cross_[1] );
    }
    
    vector<CrossPair*> paired;
    
    bool overpaired = false, matched = false;
    unordered_map<CrossBranch*, int> matches[2];
    for ( CrossBranch* cb : branch_[0] ) for ( CrossPair* cp : cb->pairs_ ) if ( cp->score_ >= max( 1, cp->miss_ * 4 ) )
    {
        paired.push_back( cp );
        for ( int d : { 0, 1 } )
        {
            auto ins = matches[d].insert( make_pair( cp->branch[d], 1 ) );
            if ( ins.second ) continue;
            ins.first->second++;
            overpaired = true;
            assert( false );
        }
    }
    
    for ( int i = 0; i < paired.size(); i++ ) if ( paired[i]->score_ < 2 && ( matches[0].size() < branch_[0].size() || matches[1].size() < branch_[1].size() ) )
    {
        for ( int d : { 0, 1 } )
        {
            auto it = matches[d].find( paired[i]->branch[d] );
            assert( it != matches[d].end() );
            if ( it->second < 2 ) matches[d].erase( it );
            else it->second--;
        }
        paired.erase( paired.begin() + i-- );
    }
    
    for ( CrossPair* cp : paired ) if ( matches[0][ cp->branch[0] ] < 2 && matches[1][ cp->branch[1] ] < 2 && separate( cp->branch, paths, nodes ) ) matched = true;
    
    if ( matched || ( forks_[0].empty() && forks_[1].empty()) ) return matched;
    
    int best = 4;
    CrossBranch* branches[2]{ NULL, NULL };
    for ( int d : { 0, 1 } ) for ( CrossBranch* fork : forks_[d] ) if ( fork->unique_.size() > 3 )
    {
        for ( CrossBranch* branch : branch_[!d] )
        {
            int score[3]{0};
            for ( ReadId id : branch->unique_ )
            {
                if ( fork->unique_.find( id ) != fork->unique_.end() ) score[1]++;
                else if ( fork->hits_.find( id ) == fork->hits_.end() ) score[0]++;
            }
            for ( ReadId id : fork->unique_ ) if ( branch->hits_.find( id ) == branch->hits_.end() ) score[2]++;
            if ( score[0] || score[1] < best ) continue;
            branches[d] = fork;
            branches[!d] = branch;
            best = score[1];
        }
    }
    
    if ( branches[0] && branches[1] ) separate( branches, paths, nodes );
    
    return matched;
}

vector<PathCross*> PathCross::create( vector<NodePath*>& paths )
{
    vector<PathCross*> crosses;
    for ( NodePath* np : paths ) if ( np->edges_[0].size() > 1 )
    {
        vector<NodePath*> cross{ np };
        int32_t off = 0;
        while ( !cross.empty() && np->edges_[1].size() == 1 )
        {
            off -= np->edges_[1][0]->diff;
            np = np->edges_[1][0]->edge[1];
            if ( np->edges_[0].size() != 1 ) cross.clear();
            else cross.push_back( np );
        }
        if ( !cross.empty() && cross[0]->isBranchable( 0 ) && cross.back()->isBranchable( 1 ) ) crosses.push_back( new PathCross( cross, off, true ) );
    }
    for ( NodePath* np : paths ) if ( np->edges_[1].size() > 1 )
    {
        for ( PathEdge* pe : np->edges_[1] ) if ( pe->edge[1]->edges_[0].size() > 1 )
        {
            vector<NodePath*> cross{ np, pe->edge[1] };
            crosses.push_back( new PathCross( cross, -pe->diff, false ) );
        }
    }
    
    return crosses;
}

void PathCross::extend( CrossBranch* cb, NodePath* branch, int32_t diff, bool split, bool drxn )
{
    if ( cb->path_.empty() || ( cb->fork_ && cb->fork_->path_.size() == cb->path_.size() ) ) cb->ends_[!drxn] = branch->ends_[!drxn] + diff;
    cb->path_.insert( drxn ? cb->path_.end() : cb->path_.begin(), branch );
    assert( cb->diffs_.insert( make_pair( branch, diff ) ).second );
    cb->ends_[drxn] = branch->ends_[drxn] + diff;
    
    vector<PathEdge*> edges = cb->getEdges( drxn );
    
    if ( abs( cb->ends_[drxn] - ends_[!drxn] ) > params.maxMpMean * 1.5 ) split = false;
    if ( branch->edges_[!drxn].size() > 1 ) split = false;
    if ( edges.size() > 1 ) for ( PathEdge* pe : edges ) if ( pe->edge[drxn]->edges_[!drxn].size() > 1 ) split = false;
    
    if ( edges.size() == 1 ) extend( cb, edges[0]->edge[drxn], diff + ( drxn ? -edges[0]->diff : edges[0]->diff ), split, drxn );
    else if ( edges.size() > 1 && split )
    {
        for ( PathEdge* pe : edges ) extend( new CrossBranch( cb ), pe->edge[drxn], diff + ( drxn ? -pe->diff : pe->diff ), true, drxn );
        forks_[drxn].push_back( cb );
    }
    else
    {
        branch_[drxn].push_back( cb );
        int32_t limit = params.maxMpMean * 1.3 + 300 - abs( cb->ends_[drxn] - ends_[!drxn] );
        unordered_set<NodePath*> used( cb->path_.begin(), cb->path_.end() );
        for ( PathEdge* pe : edges ) cb->fill( pe->edge[drxn], used, diff + ( drxn ? -pe->diff : pe->diff ), limit, drxn );
    }
}

void PathCross::match( vector< pair<NodePath*, int32_t> > tars[2], ReadId id, int32_t dist )
{
    vector< vector<int32_t> > hits[2];
    vector<bool> used[2]{ vector<bool>( tars[0].size(), false ), vector<bool>( tars[1].size(), false ) };
    for ( int d : { 1, 0 } )
    {
        bool good = false;
        int32_t* off;
        for ( CrossBranch* cb : branch_[d] )
        {
            vector<int32_t> hit;
            for ( int i = 0; i < tars[d].size(); i++ ) if ( ( off = cb->get( tars[d][i].first ) ) && ( used[d][i] = true ) ) hit.push_back( tars[d][i].second + (*off) );
            if ( !hit.empty() ) good = true;
            hits[d].push_back( hit );
        }
        if ( !good ) return;
    }
    
    int32_t best = dist * 2.6 + 600, minDist = dist * 0.7 - 200, cutoff = dist * 1.3 + 300;
    
    for ( int i = 0; i < tars[0].size(); i++ ) for ( int j = 0; j < tars[1].size(); j++ )
    {
        if ( !used[0][i] || !used[1][j] ) for ( PathPairs* pp : tars[0][i].first->paired_[1] )
        {
            if ( pp->node_[1] == tars[1][j].first && tars[1][j].second - pp->diff_ - tars[0][i].second < cutoff ) return;
        }
        else
        {
            unordered_set<NodePath*> pathed( cross_.begin(), cross_.end() );
            if ( reach( tars[0][i].first, pathed, tars[1], 0, tars[0][i].second + cutoff ) ) return;
        }
    }
    
    vector< vector<int32_t> > dists( branch_[0].size(), vector<int32_t>( branch_[1].size(), best ) );
    for ( int i = 0; i < branch_[0].size(); i++ ) for ( int32_t l : hits[0][i] ) for ( int j = 0; j < branch_[1].size(); j++ ) for ( int32_t r : hits[1][j] )
    {
        int32_t est = r - l;
        if ( est < minDist ) return;
        best = min( best, est );
        dists[i][j] = min( dists[i][j], est );
    }
    if ( cutoff < best ) return;
    best -= dist;
    
    unordered_set<CrossBranch*> matched[2];
    for ( int i = 0; i < dists.size(); i++ ) for ( int j = 0; j < dists[i].size(); j++ )
    {
        if ( cutoff < dists[i][j] && dists[i][j] < cutoff + abs( best ) ) return;
        if ( dists[i][j] < cutoff ) for ( int d : { 0, 1 } ) matched[d].insert( branch_[d][ d ? j : i ] );
    }
    for ( int d : { 0, 1 } ) if ( matched[d].size() == branch_[d].size() ) return;
    for ( int d : { 0, 1 } ) if ( matched[d].size() == 1 ) for ( CrossBranch* cb : matched[d] ) cb->unique_.insert( id );
    for ( int d : { 0, 1 } ) for ( CrossBranch* cb : branch_[d] ) ( matched[d].find( cb ) != matched[d].end() ? cb->hits_ : cb->miss_ ).insert( id );
    
    for ( int i = 0; i < dists.size(); i++ ) for ( int j = 0; j < dists[i].size(); j++ ) if ( dists[i][j] < cutoff )
    {
        if ( matched[0].size() == 1 && matched[1].size() == 1 ) branch_[0][i]->pairs_[j]->unique_.insert( id );
        branch_[0][i]->pairs_[j]->hits_.insert( id );
    }
}

bool PathCross::reach( NodePath* np, unordered_set<NodePath*>& pathed, vector< pair<NodePath*, int32_t> >& tars, int32_t diff, int32_t cutoff )
{
    for ( pair<NodePath*, int32_t>& t : tars ) if ( t.first == np && t.second - diff < cutoff ) return true;
    if ( cutoff < np->ends_[1] - cutoff || !pathed.insert( np ).second ) return false;
    for ( PathEdge* pe : np->edges_[1] ) if ( reach( pe->edge[1], pathed, tars, diff + pe->diff, cutoff) ) return true;
    pathed.erase( np );
    return false;
}

bool PathCross::resolve( vector<NodePath*>& paths, NodeRoll& nodes )
{
    bool resolved = false;
    
    for ( int again = 1; again-- > 0; )
    {
        vector<PathCross*> crosses = create( paths );
//        LocusExport( nodes, "/home/glen/LocassDump/dump2" );
        if ( crosses.empty() ) return resolved;
        for ( PathCross* pc : crosses ) if ( pc->claim( paths, nodes ) ) again = 1;
        for ( PathCross* pc : crosses ) delete pc;
        
        if ( again ) NodePath::reset( paths );
        if ( again ) resolved = true;
    }
    
    return resolved;
    AllelePaths::create( paths, nodes );
    int x = 0;
}

bool PathCross::separate( CrossBranch* paired[2], vector<NodePath*>& paths, NodeRoll& nodes )
{
    bool claimable[2][2]{ { false, false }, { false, false } };
    int len[2]{ paired[0]->fork_ ? (int)paired[0]->fork_->path_.size() : 0, paired[1]->fork_ ? (int)paired[1]->fork_->path_.size() : 0 };
    NodePath* fork[2]{ len[0] ? paired[0]->fork_->path_[0] : cross_[0], len[1] ? paired[1]->fork_->path_.back() : cross_.back() };
    NodePath* branch[2]{ len[0] ? paired[0]->path_.rbegin()[ len[0] ] : paired[0]->path_.back(), len[1] ? paired[1]->path_[ len[1] ] : paired[1]->path_[0] };
    for ( int d : { 0, 1 } ) for ( PathEdge* pe : fork[d]->edges_[d] ) claimable[d][ pe->edge[d] == branch[d] ] = true;
    for ( int d : { 0, 1 } ) for ( int i : { 0, 1 } ) if ( !claimable[d][i] ) return false;
    
    for ( NodePath* np : paths )
    {
        if ( np->path_[0]->cloned_ ) assert( !np->path_[0]->edges_[0].empty() );
        if ( np->path_.back()->cloned_ ) assert( !np->path_.back()->edges_[1].empty() );
    }
    
    NodePath* ends[2][2]{ { cross_[0], NULL }, { cross_.back(), NULL } };
    for ( int i = 0; i < cross_.size(); i++ )
    {
        NodePath* np = new NodePath( cross_[i], paths, nodes, 1 );
        if ( i ) for ( PathEdge* pe : cross_[i-1]->edges_[1] ) if ( pe->edge[1] == cross_[i] ) new PathEdge( ends[1][1], np, pe );
        for ( int d : { 0, 1 } ) if ( d || !i ) ends[d][1] = np;
    }
    
    for ( int d : { 0, 1 } ) for ( int i = 0; i < len[d]; i++ )
    {
        NodePath* np[2]{ d ? paired[1]->path_[i] : paired[0]->path_.end()[-i-1], NULL };
        np[1] = new NodePath( np[0], paths, nodes, 1 );
        for ( PathEdge* pe : np[0]->edges_[!d] ) if ( pe->edge[!d] == ends[d][0] ) new PathEdge( d ? ends[1][1] : np[1], d ? np[1] : ends[0][1], pe );
        for ( int j : { 0, 1 } ) ends[d][j] = np[j];
    }
    
    for ( int d : { 0, 1 } ) for ( PathEdge* pe : ends[d][0]->edges_[d] ) if ( pe->edge[d] == branch[d] )
    {
        if ( find( branch_[d].begin(), branch_[d].end(), paired[d] ) != branch_[d].end() ) pe->claim( ends[d][1], !d );
        else new PathEdge( d ? ends[1][1] : branch[0], d ? branch[1] : ends[0][1], pe );
        break;
    }
    
    for ( NodePath* np : paths )
    {
        if ( np->path_[0]->cloned_ ) assert( !np->path_[0]->edges_[0].empty() );
        if ( np->path_.back()->cloned_ ) assert( !np->path_.back()->edges_[1].empty() );
    }
    
    return true;
}
