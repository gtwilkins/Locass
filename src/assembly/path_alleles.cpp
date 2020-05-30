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

#include "path_alleles.h"
#include <algorithm>
#include <set>

AlleleBlock::AlleleBlock( vector<NodePath*>& path )
: homo_( true )
{
    for( NodePath* np : path ) path_[0].insert( path_[0].end(), np->path_.begin(), np->path_.end() );
    assert( !path_[0].empty() );
}

AlleleBlock::AlleleBlock( vector<Node*> paths[2], bool drxn )
: homo_( true )
{
    int len = 1, limit = min( paths[0].size(), paths[1].size() );
    while ( len < limit && ( drxn ? paths[0].end()[-len-1]->isClone( paths[1].end()[-len-1] ) : paths[0][len]->isClone( paths[1][len] ) ) ) len++;
    for ( int d : { 0, 1 } ) path_[d].insert( path_[d].end(), drxn ? paths[d].end()-len : paths[d].begin(), drxn ? paths[d].end() : paths[d].begin()+len );
    for ( int d : { 0, 1 } ) paths[d].erase( drxn ? paths[d].end()-len : paths[d].begin(), drxn ? paths[d].end() : paths[d].begin()+len );
}

AlleleBlock::AlleleBlock( vector<Node*> paths[2], vector< pair<int, int> >& matches, int& i, int& j, int& k )
{
    for ( ; k < matches.size() && i == matches[k].first && j == matches[k].second; k++ )
    {
        path_[0].push_back( paths[0][i++] );
        path_[1].push_back( paths[1][j++] );
    }
    if ( homo_ = !path_[0].empty() ) return;
    
    while ( i < ( k < matches.size() ? matches[k].first : paths[0].size() ) ) path_[0].push_back( paths[0][i++] );
    while ( j < ( k < matches.size() ? matches[k].second : paths[1].size() ) ) path_[1].push_back( paths[1][j++] );
}

int AlleleBlock::getScore( int i, int& count, bool drxn )
{
    int score = 0, f = min( count, (int)path_[ homo_ ? 0 : i ].size() );
    count -= f;
    
    if ( !homo_ || !f ) return 0;
    if ( f == path_[0].size() ) return Node::getLen( path_[0] );
    
    vector<Node*> counted( drxn ? path_[0].begin() : path_[0].end()-f, drxn ? path_[0].begin()+f : path_[0].end() );
    score = Node::getLen( counted );
    assert( false );
    return score;
}

int AlleleBlock::getScore( int i, int& uncount, int& count, bool drxn )
{
    int b = min( uncount, (int)path_[ homo_ ? 0 : i ].size() ), f = min( count, (int)path_[ homo_ ? 0 : i ].size() );
    uncount -= b;
    count -= f;
    
    if ( !homo_ || b == path_[0].size() || !f ) return 0;
    if ( homo_ && !b && f == path_[0].size() ) return Node::getLen( path_[0] );
    
    vector<Node*> counted( drxn ? path_[0].begin()+b : path_[0].end()-f, drxn ? path_[0].begin()+f : path_[0].end()-b );
    return Node::getLen( counted );
}

AlleleNode::AlleleNode( vector<NodePath*>& path )
: homo_( true )
{
    len[0] = len[1] = 0;
    for ( int i : { 0, 1 } ) for ( int j : { 0, 1 } ) ol[i][j] = 0;
    bested_[0] = bested_[1] = false;
    path_[0] = path;
    blocks_.push_back( new AlleleBlock( path_[0] ) );
    for ( int i = 0; i < path_[0].size(); i++ ) len[0] += path_[0][i]->size() - ( i ? path_[0][i]->getOverlap( path_[0][i-1], 0 ) : 0 );
}

AlleleNode::AlleleNode( vector<NodePath*> path[3], AlleleNode* base, bool drxn )
: homo_( false )
{
    vector<Node*> paths[2];
    for ( int p : { 0, 1 } )
    {
        bested_[0] = bested_[1] = false;
        path_[p] = path[p];
        len[p] = 0;
        for ( int i = 0; i < path_[p].size(); i++ ) len[p] += path_[p][i]->size() - ( i ? path_[p][i]->getOverlap( path_[p][i-1], 0 ) : 0 );
        ol[!drxn][p] = ( drxn ? base->path_[0].back() : base->path_[0][0] )->getOverlap( drxn ? path_[p][0] : path_[p].back(), drxn );
        for ( NodePath* np : path_[p] ) paths[p].insert( paths[p].end(), np->path_.begin(), np->path_.end() );
    }
    
    AlleleBlock* ends[2]{ NULL, NULL };
    if ( drxn ? paths[0][0]->isClone( paths[1][0] ) : paths[0].back()->isClone( paths[1].back() ) ) ends[!drxn] = new AlleleBlock( paths, !drxn );
    if ( !path[2].empty() && !paths[0].empty() && !paths[1].empty() && 
            ( drxn ? paths[0].back()->isClone( paths[1].back() ) : paths[0][0]->isClone( paths[1][0] ) ) ) ends[drxn] = new AlleleBlock( paths, drxn );
    
    if ( ends[0] ) blocks_.push_back( ends[0] );
    vector< pair<int, int> > cur, matches;
    int32_t best = 0;
    for ( Node* node : paths[0] ) best += node->size();
    if ( !paths[0].empty() && !paths[1].empty() ) match( paths, cur, matches, best, 0, 0 );
    int i = 0, j = 0, k = 0;
    while ( i < paths[0].size() || j < paths[1].size() ) blocks_.push_back( new AlleleBlock( paths, matches, i, j, k ) );
    if ( ends[1] ) blocks_.push_back( ends[1] );
}

AlleleNode::~AlleleNode()
{
    for ( AlleleBlock* ab : blocks_ ) delete ab;
}

//AlleleNode::AlleleNode( NodesPath* seed, vector<AlleleNode*>& path )
//{
//    path_[0].push_back( seed );
//    for ( int d : { 0, 1 } ) for ( NodesPath* fork = seed; fork && fork->edges_[d].size() == 1 && fork->edges_[d][0]->edge[d]->edges_[!d].size() == 1; )
//    {
//        fork = fork->edges_[d][0]->edge[d];
//        path_[0].insert( d ? path_[0].end() : path_[0].begin(), fork );
//    }
//    
//    blocks_.push_back( new AlleleBlock( path_[0] ) );
//    path.push_back( this );
//    for ( int d : { 0, 1 } ) extend( path, d );
//}
//
//AlleleNode::AlleleNode( NodesPath* fork, vector<AlleleNode*>& path, bool drxn )
//{
//    for ( int d : { 0, 1 } ) for ( NodesPath* np = fork->edges_[drxn][d]->edge[drxn]; np; )
//    {
//        path_[d].insert( drxn ? path_[d].end() : path_[d].begin(), np );
//        np = np->edges_[drxn].size() == 1 ? np->edges_[drxn][0]->edge[drxn] : NULL;
//    }
//    
//    path.insert( drxn ? path.end() : path.begin(), this );
//    
//    int i = 0;
//    while ( drxn ? path_[0].end()[-i-1] == path_[1].end()[-i-1] : path_[0][i] == path_[1][i] ) i++;
//    if ( i )
//    {
//        vector<NodesPath*> merge( drxn ? path_[0].end()-i : path_[0].begin(), drxn ? path_[0].end() : path_[0].begin()+i );
//        for ( int d : { 0, 1 } ) path_[d].erase( drxn ? path_[d].end()-i : path_[d].begin(), drxn ? path_[d].end() : path_[d].begin()+i );
//        new AlleleNode( merge, path, drxn );
//    }
//    
//    fill( i ? 2 : (int)drxn );
//}
//
//AlleleNode::AlleleNode( vector<NodesPath*>& merge, vector<AlleleNode*>& path, bool drxn )
//{
//    path_[0] = merge;
//    blocks_.push_back( new AlleleBlock( path_[0] ) );
//    path.insert( drxn ? path.end() : path.begin(), this );
//    extend( path, drxn );
//}
//
//void AlleleNode::extend( vector<AlleleNode*>& path, bool drxn )
//{
//    assert( !path_[0].empty() && path_[1].empty() );
//    NodesPath* fork = drxn ? path_[0].back() : path_[0][0];
//    if ( fork->edges_[drxn].size() != 2 || fork->multi_ != 2 ) return;
//    new AlleleNode( fork, path, drxn );
//}
//
//void AlleleNode::fill()
//{
//    vector<Node*> paths[2];
//    for ( int d : { 0, 1 } ) for ( NodesPath* np : path_[d] ) paths[d].insert( paths[d].end(), np->path_.begin(), np->path_.end() );
//    
//    AlleleBlock* ends[2]{ NULL, NULL };
//    if ( !paths[0].empty() && !paths[1].empty() && paths[0][0]->isClone( paths[1][0] ) ) ends[0] = new AlleleBlock( paths, 0 );
//    if ( !paths[0].empty() && !paths[1].empty() && paths[0].back()->isClone( paths[1].back() ) ) ends[1] = new AlleleBlock( paths, 1 );
//    
//    if ( ends[0] ) blocks_.push_back( ends[0] );
//    
//    vector< pair<int, int> > cur, matches;
//    int32_t best = 0;
//    for ( Node* node : paths[0] ) best += node->size();
//    if ( !paths[0].empty() && !paths[1].empty() ) match( paths, cur, matches, best, 0, 0 );
//    int i = 0, j = 0, k = 0;
//    while ( i < paths[0].size() || j < paths[1].size() ) blocks_.push_back( new AlleleBlock( paths, matches, i, j, k ) );
//    
//    if ( ends[1] ) blocks_.push_back( ends[1] );
//}

void AlleleNode::match( vector<Node*> paths[2], vector< pair<int, int> >& cur, vector< pair<int, int> >& matches, int32_t& best, int i, int j )
{
    bool matched = cur.empty();
    int limit = paths[1].size();
    for ( ; i < paths[0].size(); i++ ) if ( paths[0][i]->cloned_ )
    {
        for ( int k = limit; k-- > j; )
        {
            if ( paths[0][i]->isClone( paths[1][k] ) )
            {
                cur.push_back( make_pair( i, k ) );
                match( paths, cur, matches, best, i+1, k+1 );
                matched = true;
                cur.pop_back();
                limit = k;
            }
        }
    }
    if ( matched ) return;
    
    int32_t score = 0;
    auto it = cur.begin();
    vector<Node*> unmatched;
    for ( int k = 0; k < paths[0].size(); k++ )
    {
        if ( matched = ( it != cur.end() && k == it->second ) ) it++;
        else unmatched.push_back( paths[0][k] );
        if ( !unmatched.empty() && ( matched || k+1 == paths[0].size() ) )
        {
            score += Node::getLen( unmatched );
            unmatched.clear();
        }
    }
    
    if ( score < best )
    {
        matches = cur;
        best = score;
    }
}

int AlleleNode::getScore()
{
    int score = 0;
    for ( AlleleBlock* ab : blocks_ ) if ( ab->homo_ ) score += Node::getLen( ab->path_[0] );
    return score;
}

int AlleleNode::getScore( int i, int& j, bool drxn )
{
    if ( homo_ ) i = 0;
    if ( j <= 0 ) return 0;
    
    int scored = 0, score = 0;
    for ( int k = 0; k < path_[i].size() && j && j--; k++ ) scored += ( drxn ? path_[i][k] : path_[i].end()[-k-1] )->path_.size();
    
    for ( int k = 0; k < blocks_.size() && scored >= 0; k++ ) score += ( drxn ? blocks_[k] : blocks_.end()[-k-1] )->getScore( i, scored, drxn );
    
    return score;
}

int AlleleNode::getScore( int i, int& b, int& f, bool drxn )
{
    if ( homo_ ) i = 0;
    if ( b <= 0 && f <= 0 ) return 0;
    
    int unscored = 0, scored = 0, score = 0;
    for ( int k = 0; k < path_[i].size() && ( ( b && b-- ) || ( f && f-- ) ); k++ )
    {
        ( b ? unscored : scored ) += ( drxn ? path_[i][k] : path_[i].end()[-k-1] )->path_.size();
    }
    
    for ( int k = 0; k < blocks_.size() && scored > 0; k++ ) score += ( drxn ? blocks_[k] : blocks_.end()[-k-1] )->getScore( i, unscored, scored, drxn );
    
    return score;
}

//AlleleMatch::AlleleMatch( AlleleBranch* a, AlleleBranch* b )
//: hits( 0 ), miss( 0 )
//{
//    branch[0] = a;
//    branch[1] = b;
//    a->pairs_.push_back( this );
//    b->pairs_.push_back( this );
//    
//    vector< pair<int, int> > cur;
//    match( cur, 0, 0 );
//}
//
//bool AlleleMatch::advance( bool drxn )
//{
//    int len[2]{ branch[0]->dists_.back().second, branch[1]->dists_.back().second };
//    int limit[2]{ matched.empty() ? 0 : branch[0]->dists_[ matched.back().first ].second
//                , matched.empty() ? 0 : branch[1]->dists_[ matched.back().second ].second };
//    for ( int d : { 0, 1 } ) if ( len[d] > len[!d] ) assert( len[d] - limit[d] > len[!d] - limit[!d] );
//    
//    vector<AlleleBranch*> advanced[2];
//    for ( int d : { 0, 1 } )
//    {
//        if ( len[d] > len[!d] ) advanced[d].push_back( branch[d] );
//        else for ( PathEdge* pe : ( drxn ? branch[d]->path_.back() : branch[d]->path_[0] )->edges_[drxn] ) advanced[d].push_back( new AlleleBranch( branch[d], pe, drxn ) );
//    }
//    
//    return false;
//}
//
//void AlleleMatch::match( vector< pair<int, int> >& cur, int i, int j )
//{
//    bool record = !cur.empty();
//    int limit = branch[1]->dists_.size();
//    for ( ; i < branch[0]->dists_.size(); i++ ) if ( branch[0]->dists_[i].first->cloned_ )
//    {
//        for ( int k = limit; k-- > j; )
//        {
//            if ( branch[0]->dists_[i].first->isClone( branch[1]->dists_[k].first ) )
//            {
//                cur.push_back( make_pair( i, k ) );
//                match( cur, i+1, k+1 );
//                record = false;
//                cur.pop_back();
//                limit = k;
//            }
//        }
//    }
//    if ( !record ) return;
//    
//    int score[2]{0}, limits[2]{ -params.readLen, -params.readLen };
//    i = j = 0;
//    for ( pair<int, int>& paired : cur )
//    {
//        int32_t dist = min( branch[0]->dists_[ paired.first ].second - limits[0], branch[1]->dists_[ paired.second ].second - limits[1] );
//        int32_t base[2]{ branch[0]->dists_[i].second - branch[0]->dists_[i].first->size(), branch[1]->dists_[j].second - branch[1]->dists_[j].first->size() };
//        int32_t misses[2]{0};
//        while ( ++i < paired.first ) misses[0] = branch[0]->dists_[i].second - base[0];
//        while ( ++j < paired.second ) misses[1] = branch[1]->dists_[j].second - base[1];
//        score[0] += misses[0] + misses[1];
//        score[1] += min( branch[0]->dists_[ paired.first ].first->size(), dist );
//        limits[0] = branch[0]->dists_[i].second;
//        limits[1] = branch[1]->dists_[j].second;
//    }
//    
//    if ( score[1] > hits || ( score[1] == hits && score[0] < miss ) )
//    {
//        matched = cur;
//        hits = score[1];
//        miss = score[0];
//    }
//}

//AlleleBranch::AlleleBranch( vector<NodesPath*>& path, int32_t ol, bool drxn )
//: path_( path )
//{
//    int32_t dist = -ol;
//    vector<Node*> nodes;
//    for ( NodesPath* np : path ) nodes.insert( nodes.end(), np->path_.begin(), np->path_.end() );
//    if ( !drxn ) reverse( nodes.begin(), nodes.end() );
//    
//    for ( int i = 0; i < nodes.size(); i++ )
//    {
//        dist += nodes[i]->size() - ( i ? nodes[i-1]->getOverlap( nodes[i], drxn ) : 0 );
//        dists_.push_back( make_pair( nodes[i], dist ) );
//    }
//}
//
//AlleleBranch::AlleleBranch( AlleleBranch* fork, PathEdge* edge, bool drxn )
//: path_( fork->path_ ), dists_( fork->dists_ )
//{
//    int32_t dist = dists_.back().second - edge->ol;
//    
//    vector<Node*> nodes;
//    for ( PathEdge* pe = edge; pe; )
//    {
//        path_.insert( drxn ? path_.end() : path_.begin(), pe->edge[drxn] );
//        nodes.insert( nodes.end(), pe->edge[drxn]->path_.begin(), pe->edge[drxn]->path_.end() );
//        pe = pe->edge[drxn]->edges_[drxn].size() == 1 ? pe->edge[drxn]->edges_[drxn][0] : NULL;
//    }
//    if ( !drxn ) reverse( nodes.begin(), nodes.end() );
//    
//    for ( int i = 0; i < nodes.size(); i++ )
//    {
//        dist += nodes[i]->size() - ( i ? nodes[i-1]->getOverlap( nodes[i], drxn ) : 0 );
//        dists_.push_back( make_pair( nodes[i], dist ) );
//    }
//}
//
//vector<AlleleBranch*> AlleleBranch::split( bool drxn )
//{
//    vector<AlleleBranch*> branches;
//    
//    for ( PathEdge* pe : ( drxn ? path_.back() : path_[0] )->edges_[drxn] ) branches.push_back( new AlleleBranch( this, pe, drxn ) );
//    
//    return branches;
//}

AlleleGraph::AlleleGraph( NodePath* branch, unordered_set<NodePath*> pathed[2], vector<AlleleGraph*> graphs[2], int i, bool drxn )
: paired( NULL ), len( 0 ), multi( branch->multi_ )
{
    graphs[i].push_back( this );
    
    vector<PathEdge*> pes;
    while ( branch )
    {
        path.push_back( branch );
        len += branch->size();
        for ( int j = 0; j < branch->path_.size(); j++ ) nodes.push_back( drxn ? branch->path_[j] : branch->path_.end()[-j-1] );
        pes.clear();
        bool again = pathed[!i].find( branch ) == pathed[!i].end();
        if ( !again ) for ( AlleleGraph* ag : graphs[!i] ) if ( ag->path[0] == branch )
        {
            paired = ag;
            ag->paired = this;
        }
        for ( PathEdge* pe : branch->edges_[drxn] ) if ( pathed[i].find( pe->edge[drxn] ) != pathed[i].end() )
        {
            pes.push_back( pe );
            if ( pathed[!i].find( pe->edge[drxn] ) != pathed[!i].end() ) again = false;
            for ( PathEdge* re : pe->edge[drxn]->edges_[!drxn] ) if ( re->edge[!drxn] != branch && pathed[i].find( re->edge[!drxn] ) != pathed[i].end() ) again = false;
        }
        branch = again && pes.size() == 1 ? pes[0]->edge[drxn] : NULL;
        if ( branch ) len -= pes[0]->ol;
        if ( branch ) multi = min( multi, branch->multi_ );
    }
    
    for ( PathEdge* pe : pes )
    {
        AlleleGraph* edge = NULL;
        for ( AlleleGraph* ag : graphs[i] ) if ( ag->path[0] == pe->edge[drxn] && ( edge = ag ) ) break;
        if ( !edge ) edge = new AlleleGraph( pe->edge[drxn], pathed, graphs, i, drxn );
        edges[1].push_back( make_pair( edge, pe->ol ) );
        edge->edges[0].push_back( make_pair( this, pe->ol ) );
    }
}

void AlleleGraph::set()
{
    for ( pair<AlleleGraph*, int32_t>& edge : edges[1] ) setGaps( edge.first, -edge.second, false );
    sort( branches.begin(), branches.end(), []( pair<AlleleGraph*, int32_t>& a, pair<AlleleGraph*, int32_t>& b ){ return a.second < b.second; } );
    for ( int i = 0; i+1 < branches.size(); i++ ) for ( int j = i+1; j < branches.size(); j++ )
    {
        if ( branches[i].first == branches[j].first ) branches.erase( branches.begin() + j-- );
    }
}

void AlleleGraph::setGaps( AlleleGraph* ag, int32_t dist, bool branched )
{
    auto it = gaps.insert( make_pair( ag, dist ) );
    if ( !it.second )
    {
        if ( it.first->second <= dist ) return;
        it.first->second = dist;
    }
    
    if ( !branched && ( ag->paired || !ag->scores.empty() ) && ( branched = true ) ) branches.push_back( make_pair( ag, dist ) );
    
    for ( pair<AlleleGraph*, int>& edge : ag->edges[1] ) setGaps( edge.first, dist + ag->len - edge.second, branched );
}

bool AlleleTar::accept( int32_t est, int32_t dist )
{
    return est <= dist * 1.2 + 300 && dist * 0.7 - 200 <= est;
}

void AlleleTar::add( NodePath* np, ReadId id, bool hit )
{
    auto ins = ( hit ? hits_ : miss_ ).insert( make_pair( np, unordered_set<ReadId>{ id } ) );
    if ( !ins.second ) assert( ins.first->second.insert( id ).second );
}

bool AlleleTar::confirmHit( PathPair* p, bool drxn )
{
    pair<int32_t, int>* t;
    for ( pair<NodePath*, int32_t> mark : p->marks_[!drxn] ) if ( tar_.find( mark.first ) == tar_.end() ) return false;
    for ( pair<NodePath*, int32_t> mark : p->marks_[drxn] )
    {
        auto it = que_.find( mark.first );
        if ( it == que_.end() ) return false;
        bool good = false;
        for ( pair<NodePath*, int32_t> tar : p->marks_[!drxn] ) for ( PathPairs* pp : mark.first->paired_[!drxn] )
        {
            if ( pp->node_[!drxn] == tar.first && ( t = get( tar.first ) ) )
            {
                int32_t est = ( ( mark.second + it->second ) - ( tar.second + t->first ) ) * ( drxn ? 1 : -1 );
                if ( ignoreHit( tar.first, mark.first, est, p->dist_, drxn ) ) return false;
                if ( accept( est, p->dist_ ) ) good = true;
            }
        }
        if ( !good ) return false;
    }
    return true;
}

bool AlleleTar::confirmMiss( PathPair* p, bool drxn )
{
    for ( pair<NodePath*, int32_t> mark : p->marks_[!drxn] ) if ( tar_.find( mark.first ) != tar_.end() ) return false;
    for ( pair<NodePath*, int32_t> mark : p->marks_[drxn] ) if ( mark.first->multi_ > 1 || que_.find( mark.first ) != que_.end() ) return false;
    for ( pair<NodePath*, int32_t> q : p->marks_[drxn] ) for ( pair<NodePath*, int32_t> t : p->marks_[drxn] )
    {
        if ( ignoreMiss( t.first, q.first, drxn ) ) assert( false );
        if ( ignoreMiss( t.first, q.first, drxn ) ) return false;
    }
    assert( false );
    return true;
}

bool AlleleTar::extend( bool drxn )
{
    vector< pair<NodePath*, bool> > branches;
    for ( PathEdge* pe : path_.back()->edges_[drxn] ) if ( que_.find( pe->edge[drxn] ) != que_.end() ) branches.push_back( make_pair( pe->edge[drxn], true ) );
    
    if ( branches.empty() ) return false;
    if ( branches.size() == 1 )
    {
        path_.push_back( branches[0].first );
        return true;
    }
    
    vector< unordered_set<ReadId> > ids( branches.size() );
    vector<int> miss( branches.size(), 0 );
    for ( int i = 0; i < branches.size(); i++ )
    {
        unordered_set<NodePath*> pathed;
        miss[i] = setScore( branches[i].first, pathed, ids[i], drxn );
    }
    
    for ( int i = 0; i < branches.size(); i++ ) if ( branches[i].second && !ids[i].empty() )
    {
        for ( int j = 0; j < branches.size(); j++ ) if ( i != j && branches[j].second && ids[i].size() > ids[j].size() )
        {
            int alts[2]{ miss[j], 0 }, base[2]{ miss[i], 0 };
            for ( ReadId id : ids[j] ) if ( ids[i].find( id ) == ids[i].end() ) alts[1]++;
            base[1] = (int)ids[i].size() + alts[1] - (int)ids[j].size();
            if ( ( base[0] + alts[1] ) * 4 < base[1] + alts[0] - 1 || alts[1] * 4 < base[1] - 1 )
            {
                branches[j].second = false;
            }
        }
    }
    
    NodePath* ext = NULL;
    int count = 0;
    for ( int i = 0; i < branches.size(); i++ ) if ( branches[i].second && ( ext = branches[i].first )) count++;
    for ( int i = 0; i < branches.size(); i++ ) if ( !branches[i].second )
    {
        unordered_set<NodePath*> keep;
        for ( int j = 0; j < branches.size(); j++ ) if ( i != j && branches[j].second ) secure( branches[j].first, branches[i].first, keep, drxn );
        assert( !keep.empty() );
        for ( auto q = que_.begin(); q != que_.end(); )
        {
            if ( keep.find( q->first ) != keep.end() ) q++;
            else
            {
                hits_.erase( q->first );
                q = que_.erase( q );
            }
        }
    }
    if ( ext && count == 1 )
    {
        path_.push_back( ext );
        return true;
    }
    return false;
}

void AlleleTar::fill( NodePath* np, unordered_set<NodePath*>& ignores, bool complete, bool drxn )
{
    auto it = tar_.find( np );
    if ( it == tar_.end() || ( !complete && it->second.second >= np->multi_ ) || !ignores.insert( np ).second ) return;
    for ( PathEdge* pe : np->edges_[drxn] ) fill( pe->edge[drxn], ignores, complete, drxn );
}

pair<int32_t, int>* AlleleTar::get( NodePath* np )
{
    auto it = tar_.find( np );
    return it != tar_.end() ? &it->second : NULL;
}

bool AlleleTar::ignoreHit( NodePath* t, NodePath* q, int32_t est, int32_t dist, bool drxn )
{
    auto it = ignores_.find( t );
    if ( it == ignores_.end() ) return false;
    auto it2 = it->second.find( q );
    if ( it2 == it->second.end() ) return false;
    auto it3 = que_.find( q );
    assert( it3 != que_.end() );
    int32_t ext = drxn ? it2->second - it3->second : it3->second - it2->second;
    int32_t cutoff = dist * 1.2 + 300 + abs( dist - est );
    return est + ext < dist * 1.2 + 300 + abs( dist - est );
}

bool AlleleTar::ignoreMiss( NodePath* t, NodePath* q, bool drxn )
{
    unordered_set<NodePath*> pathed;
    return reach( t, q, pathed, drxn );
}

void AlleleTar::insert( NodePath* np, int32_t off, int multi )
{
    ends_[0] = min( ends_[0], np->ends_[0] + off );
    ends_[1] = max( ends_[1], np->ends_[1] + off );
    auto it = tar_.insert( make_pair( np, make_pair( off, multi ) ) );
    if ( it.second ) return;
    it.first->second.second += multi;
    auto ins = loop_.insert( make_pair( np, 1 ) );
    if ( !ins.second ) ins.first->second++;
}

void AlleleTar::match( bool drxn )
{
    prepare( drxn );
    pair<int32_t, int>* t;
    unordered_set<ReadId> used;
    for ( const pair<NodePath*, int32_t>& q : que_ ) for ( PathPairs* pp : q.first->paired_[!drxn] )
    {
        t = get( pp->node_[!drxn] );
        for ( PathUnique& pu : pp->uniques_ )
        {
            if ( t )
            {
                int32_t est = ( pu.r_ + ( drxn ? q.second : t->first ) ) - ( pu.l_ + ( drxn ? t->first : q.second ) );
                if ( accept( est, pu.dist_ ) && !ignoreHit( pp->node_[!drxn], q.first, est, pu.dist_, drxn ) ) add( q.first, pu.id_, true );
            }
            else if ( q.first->multi_ == 1 && accept( pu.best_ + pu.dist_, pu.dist_ ) && !ignoreMiss( pp->node_[!drxn], q.first, drxn ) )
            {
                int32_t est = ( drxn ? pu.r_ : pu.l_ ) + q.second + ( ( drxn ? -1 : 1 ) * ( pu.dist_ * 1.2 + 300 ) );
                if ( drxn ? ends_[0] < est : est < ends_[1] ) add( q.first, pu.id_, false );
            }
        }
        for ( PathPair* p : pp->shared_ ) if ( used.insert( p->id_ ).second && ( t ? confirmHit( p, drxn ) : confirmMiss( p, drxn ) ) )
        {
            for ( pair<NodePath*, int32_t> mark : p->marks_[drxn] ) add( mark.first, p->id_, t );
        }
    }
    
    for ( PathEdge* pe : path_.back()->edges_[drxn] ) if ( que_.find( pe->edge[drxn] ) != que_.end() ) setBranch( pe->edge[drxn], drxn );
}

void AlleleTar::prepare( bool drxn )
{
    ignores_.clear();
    for ( const pair<NodePath*, pair<int32_t, int> >& t : tar_ )
    {
        for ( PathEdge* pe : t.first->edges_[drxn] ) if ( tar_.find( pe->edge[drxn] ) == tar_.end() )
        {
            unordered_set<NodePath*> ignores;
            fill( t.first, ignores, que_.find( pe->edge[drxn] ) == que_.end(), !drxn );
            if ( !ignores.empty() )
            {
                vector<NodePath*> path;
                reach( pe->edge[drxn], path, ignores, t.second.first + ( drxn ? -pe->diff : pe->diff ), params.maxMpMean*1.2 + 300 - pe->ol, drxn );
            }
        }
    }
    
}

bool AlleleTar::reach( NodePath* np, NodePath* q, unordered_set<NodePath*>& pathed, bool drxn )
{
    if ( np == q ) return true;
    if ( que_.find( np ) == que_.end() || !pathed.insert( np ).second ) return false;
    for ( PathEdge* pe : np->edges_[drxn] ) if ( reach( pe->edge[drxn], q, pathed, drxn ) ) return true;
    pathed.erase( np );
    return false;
}

void AlleleTar::reach( NodePath* np, vector<NodePath*>& path, unordered_set<NodePath*>& ignores, int32_t diff, int32_t limit, bool drxn )
{
    if ( find( path.begin(), path.end(), np ) != path.end() ) return;
    auto it = tar_.find( np );
    if ( it != tar_.end() && it->second.second >= np->multi_ ) return;
    
    if ( que_.find( np ) != que_.end() ) for ( NodePath* t : ignores )
    {
        auto found = ignores_.find( t );
        if ( found != ignores_.end() )
        {
            auto ins = found->second.insert( make_pair( np, diff ) );
            if ( !ins.second )
            {
                if ( drxn ? ins.first->second <= diff : diff <= ins.first->second ) return;
                ins.first->second = diff;
            }
        }
        else ignores_.insert( make_pair( t, unordered_map<NodePath*, int32_t>{ make_pair( np, diff ) } ) );
    }
    
    path.push_back( np );
    limit -= np->size();
    if ( limit > 0 ) for ( PathEdge* pe : np->edges_[drxn] ) reach( pe->edge[drxn], path, ignores, diff + ( drxn ? -pe->diff : pe->diff ), limit + pe->ol, drxn );
    path.pop_back();
}

void AlleleTar::secure( NodePath* np, NodePath* block, unordered_set<NodePath*>& keep, bool drxn )
{
    if ( np == block || que_.find( np ) == que_.end() || !keep.insert( np ).second ) return;
    for ( PathEdge* pe : np->edges_[drxn] ) secure( pe->edge[drxn], block, keep, drxn );
}

int AlleleTar::setBranch( NodePath* np, bool drxn )
{
    auto it = branches_.find( np );
    if ( it != branches_.end() ) return it->second.second;
    
    branches_.insert( make_pair( np, make_pair( (NodePath*)NULL, 0 ) ) );
    
    vector<NodePath*> branched;
    int fork = setForks( np, branched, drxn );
    
    unordered_set<ReadId> ids[2];
    for ( const pair<NodePath*, unordered_set<ReadId> >& miss : miss_ )
    {
        int i = find( branched.begin(), branched.end(), miss.first ) != branched.end();
        ids[i].insert( miss.second.begin(), miss.second.end() );
    }
    
    it = branches_.find( np );
    assert( it != branches_.end() && !branched.empty() );
    
    if ( forks_.find( branched.back() ) != forks_.end() ) it->second.first = branched.back();
    for ( ReadId id : ids[1] ) if ( ids[0].find( id ) == ids[0].end() ) it->second.second++;
    return fork + it->second.second;
}

int AlleleTar::setForks( NodePath* np, vector<NodePath*>& branched, bool drxn )
{
    branched.push_back( np );
    vector< pair<NodePath*, int> > branches;
    for ( PathEdge* pe : np->edges_[drxn] ) if ( que_.find( pe->edge[drxn] ) != que_.end() ) branches.push_back( make_pair( pe->edge[drxn], 0 ) );
    
    if ( branches.size() == 1 ) return setForks( branches[0].first, branched, drxn );
    else if ( branches.size() > 1 )
    {
        auto ins = forks_.insert( make_pair( np, branches ) );
        if ( !ins.second ) return ins.first->second.empty() ? 0 : ins.first->second[0].second;
        for ( pair<NodePath*, int>& branch : branches ) branch.second = setBranch( branch.first, drxn );
        
        int limit = branches[0].second;
        for ( int i = 1; i < branches.size(); i++ ) limit = min( limit, branches[i].second );
        for ( int i = 0; i < branches.size(); i++ ) if ( !limit || limit < branches[0].second ) branches.erase( branches.begin() + i-- );
        auto it = forks_.find( np );
        assert( it != forks_.end() );
        if ( branches.empty() ) forks_.erase( it );
        else it->second = branches;
        return branches.empty() ? 0 : branches[0].second;
    }
    return 0;
}

int AlleleTar::setScore( NodePath* np, unordered_set<NodePath*>& pathed, unordered_set<ReadId>& hits, bool drxn )
{
    if ( que_.find( np ) == que_.end() || !pathed.insert( np ).second ) return 0;
    
    auto it = hits_.find( np );
    if ( it != hits_.end() ) for ( ReadId id : it->second ) hits.insert( id );
    
    for ( PathEdge* pe : np->edges_[drxn] ) setScore( pe->edge[drxn], pathed, hits, drxn );
    
    auto branch = branches_.find( np );
    if ( branch == branches_.end() ) return 0;
    if ( !branch->second.first ) return branch->second.second;
    auto fork = forks_.find( branch->second.first );
    assert( fork != forks_.end() && !fork->second.empty() );
    return branch->second.second + fork->second[0].second;
}

void AlleleTar::survey( NodePath* np, NodePath* alt, unordered_map<NodePath*, int>& block, int32_t diff, int32_t limit, bool drxn )
{
    auto it = block.find( np );
    if ( it != block.end() && it->second >= np->multi_ ) return;
    
    auto ins = que_.insert( make_pair( np, drxn ? -diff : diff ) );
    if ( !ins.second )
    {
        if ( drxn ? ins.first->second <= -diff : diff <= ins.first->second ) return;
        ins.first->second = drxn ? -diff : diff;
    }
    
    limit -= np->size();
    if ( np == alt ) merged_ = true;
    else if ( limit < 0 || np->edges_[drxn].empty() ) fanned_ = true;
    else for ( PathEdge* pe : np->edges_[drxn] ) survey( pe->edge[drxn], alt, block, diff + pe->diff, limit + pe->ol, drxn );
}

AlleleAlign::AlleleAlign( AlleleGraph* branch[2], bool best )
: score( 0 )
{
    for ( int i : { 0, 1 } ) path[i].push_back( branch[i] );
    if ( best )
    {
        for ( int i : { 0, 1 } ) while ( path[i].back()->edges[1].size() == 1 )
        {
            if ( find( path[i].begin(), path[i].end(), path[i].back()->edges[1][0].first ) != path[i].end() ) assert( false );
            path[i].push_back( path[i].back()->edges[1][0].first );
        }
    }
    else for ( int i : { 0, 1 } ) align[i].push_back( make_pair( branch[i], -1 ) );
}

bool AlleleAlign::add( AlleleGraph* a, AlleleGraph* b, int i, int j, bool drxn )
{
    if ( score < b->scores[j].second ) return false;
    
    vector<AlleleGraph*> added[2];
    for ( int k : { 0, 1 } ) if ( path[k].back() != ( k ? b : a ) )
    {
        for ( pair<AlleleGraph*, int32_t>& edge : path[k].back()->edges[1] ) if ( edge.first == ( k ? b : a ) ) added[k] = vector<AlleleGraph*>{ edge.first };
        if ( !added[k].empty() ) continue;
        unordered_set<AlleleGraph*> joins, pathed;
        for ( pair<AlleleGraph*, int32_t>& edge : path[k].back()->edges[1] ) join( edge.first, ( k ? b : a ), joins, pathed );
        
        for ( AlleleGraph* ag = path[k].back(); ag && ag != ( k ? b : a ); )
        {
            vector<AlleleGraph*> edges;
            for ( pair<AlleleGraph*, int32_t>& edge : ag->edges[1] ) if ( joins.find( edge.first ) != joins.end() ) edges.push_back( edge.first );
            ag = edges.size() == 1 ? edges[0] : NULL;
            if ( find( added[k].begin(), added[k].end(), ag ) != added[k].end() ) assert( false );
            if ( ag ) added[k].push_back( ag );
        }
        if ( !added[k].empty() && added[k].back() == ( k ? b : a ) ) continue;
        
        vector<AlleleGraph*> backs = { k ? b : a };
        for ( AlleleGraph* ag = backs[0]; ag; )
        {
            vector<AlleleGraph*> edges;
            for ( pair<AlleleGraph*, int32_t>& edge : ag->edges[0] ) if ( joins.find( edge.first ) != joins.end() ) edges.push_back( edge.first );
            ag = edges.size() == 1 ? edges[0] : NULL;
            if ( find( added[k].begin(), added[k].end(), ag ) != added[k].end() ) assert( false );
            if ( find( added[k].begin(), added[k].end(), ag ) != added[k].end() ) break;
            if ( ag ) backs.push_back( ag );
        }
        added[k].insert( added[k].end(), backs.rbegin(), backs.rend() );
    }
    
    for ( int k : { 0, 1 } ) for ( AlleleGraph* ag : added[k] )
    {
        int count = multi( ag, k, true );
        if ( ag->paired ) for ( AlleleGraph* alt : added[!k] ) if ( alt == ag->paired ) count++;
        if ( count >= ag->multi ) return false;
    }
    
    for ( int k : { 0, 1 } ) path[k].insert( path[k].end(), added[k].begin(), added[k].end() );
    
    b->scores[j].second = score;
    gaps[0].push_back( gap( align[0].back().first, align[0].back().second, a, i, drxn ) );
    gaps[1].push_back( gap( align[1].back().first, align[1].back().second, b, j, drxn ) );
    align[0].push_back( make_pair( a, i ) );
    align[1].push_back( make_pair( b, j ) );
    lens.push_back( a->nodes[ a->scores[i].first ]->size() );
    score += lens.back() + min( 0, min( gaps[0].back(), gaps[1].back() ) );
    return true;
}

bool AlleleAlign::advance( int i, AlleleGraph* b, int j )
{
    if ( align[1][i].first == b ) return align[1][i].second < j;
    return align[1][i].first->gaps.find( b ) != align[1][i].first->gaps.end();
}

AlleleAlign AlleleAlign::branch( AlleleGraph* b, int j )
{
    int i = lens.size() - bool( b );
    if ( b ) while ( i && !advance( i, b, j ) ) i--;
    
    AlleleAlign aa;
    for ( int k : { 0, 1 } ) aa.align[k].insert( aa.align[k].end(), align[k].begin(), align[k].begin() + i + 1 );
    for ( int k : { 0, 1 } ) aa.gaps[k].insert( aa.gaps[k].end(), gaps[k].begin(), gaps[k].begin() + i );
    aa.lens.insert( aa.lens.end(), lens.begin(), lens.begin() + i );
    aa.rescore();
    
    return aa;
}

void AlleleAlign::challenge( AlleleAlign& challenger )
{
    if ( !challenger.score || challenger.score < score ) return;
    for ( int i : { 0, 1 } ) align[i] = challenger.align[i];
    for ( int i : { 0, 1 } ) path[i] = challenger.path[i];
    for ( int i : { 0, 1 } ) gaps[i] = challenger.gaps[i];
    lens = challenger.lens;
    score = challenger.score;
    
}

bool AlleleAlign::join( AlleleGraph* ag, AlleleGraph* tar, unordered_set<AlleleGraph*>& joins, unordered_set<AlleleGraph*>& pathed )
{
    bool joined = ( ag == tar ) || ( joins.find( ag ) != joins.end() );
    if ( !joined && !pathed.insert( ag ).second ) return false;
    if ( !joined ) for ( pair<AlleleGraph*, int32_t>& edge : ag->edges[1] ) if ( join( edge.first, tar, joins, pathed ) ) joined = true;
    if ( joined ) joins.insert( ag );
    return joined;
}

int32_t AlleleAlign::gap( AlleleGraph* base, int i, AlleleGraph* ext, int j, bool drxn )
{
    i = i >= 0 ? base->scores[i].first : -1;
    j = ext->scores[j].first;
    int32_t len = ( i >= 0 && ext == base ) ? -base->nodes[i]->getOverlap( base->nodes[i+1], drxn ) : 0;
    
    if ( base != ext )
    {
        if ( i < 0 ) len += base->nodes[++i]->size();
        while ( ++i < base->nodes.size() ) len += base->nodes[i]->size() - base->nodes[i-1]->getOverlap( base->nodes[i], drxn );
        auto it = base->gaps.find( ext );
        assert( it != base->gaps.end() );
        len += it->second;
    }
    
    for ( int k = ( base == ext ? i+1 : 0 ); k < j; ) len += ext->nodes[k]->size() - ext->nodes[k]->getOverlap( ext->nodes[++k], drxn );
    
    return len;
}

int AlleleAlign::multi( AlleleGraph* ag, int i, bool paired )
{
    int count = 0;
    for ( AlleleGraph* pathed : path[i] ) if ( pathed == ag ) count++;
    if ( paired && ag->paired ) for ( AlleleGraph* pathed : path[!i] ) if ( pathed == ag->paired ) count++;
    return count;
}

bool AlleleAlign::output( vector<NodePath*> ext[3], bool drxn )
{
    for ( int i : { 0, 1 } ) for ( int j = 1; j < path[i].size(); j++ )
    {
        bool joined = false;
        for ( pair<AlleleGraph*, int32_t>& edge : path[i][j-1]->edges[1] ) if ( edge.first == path[i][j] ) joined = true;
        if ( joined ) continue;
        
        bool erased = false;
        for ( int k = 1; k < align[i].size(); k++ ) if ( find( path[i].begin(), path[i].begin()+j, align[i][k].first ) == path[i].begin()+j )
        {
            int cut[2];
            for ( int d : { 0, 1 } ) cut[d] = find( path[d].begin(), path[d].end(), align[d][k].first ) - path[d].begin();
            if ( cut[!i] < path[!i].size() && align[!i][k].first == align[!i][k-1].first ) cut[!i]++;
            for ( int d : { 0, 1 } ) path[d].erase( path[d].begin() + cut[d], path[d].end() );
            for ( int d : { 0, 1 } ) align[d].erase( align[d].begin()+k, align[d].end() );
            erased = true;
        }
        assert( erased );
        path[i].erase( path[i].begin() + j, path[i].end() );
    }
    
    int merge = 0;
    for ( ; merge < min( path[0].size(), path[1].size() ) && path[0].end()[-merge-1]->paired == path[1].end()[-merge-1]; merge++ );
    for ( int i : { 0, 1 } ) for ( int j = 0; j + merge < path[i].size(); j++ ) ext[i].insert( ext[i].end(), path[i][j]->path.begin(), path[i][j]->path.end() );
    for ( int i = path[0].size()-merge; i < path[0].size(); i++ ) ext[2].insert( ext[2].end(), path[0][i]->path.begin(), path[0][i]->path.end() );
    if ( !drxn ) for ( int i = 0; i < 3; i++ ) reverse( ext[i].begin(), ext[i].end() );
    return true;
//    return !ext[2].empty();
}

void AlleleAlign::rescore()
{
    score = 0;
    for ( int i = 0; i < lens.size(); i++ ) score += lens[i] + min( 0, min( gaps[0][i], gaps[1][i] ) );
    for ( int i : { 0, 1 } ) path[i] = vector<AlleleGraph*>{ align[i][0].first };
    for ( int i : { 0, 1 } ) for ( int j = 1; j < align[i].size(); j++ ) if ( align[i][j].first != path[i].back() ) path[i].push_back( align[i][j].first );
}

bool AlleleAlign::terminate( AlleleGraph* a, AlleleGraph* b, int i, int j )
{
    for ( int k = 1; k-1 < lens.size(); k++ )
    {
        if ( align[0][k].first == a && align[0][k].second == i && align[1][k].first == b && align[1][k].second == j ) return true;
    }
    return false;
}

void AlleleAlign::trim( int i )
{
    if ( i == lens.size() ) return;
    for ( int j : { 0, 1 } ) align[j].erase( align[j].begin() + i + 1, align[j].end() );
    for ( int j : { 0, 1 } ) gaps[j].erase( gaps[j].begin() + i, gaps[j].end() );
    lens.erase( lens.begin() + i, lens.end() );
    rescore();
}

AllelePaths::AllelePaths( NodePath* seed, unordered_set<NodePath*>& used )
: id_( seed->id_ ), score_( 0 ), contested_( false )
{
    ambiguous_[0] = ambiguous_[1] = true;
    
    vector<NodePath*> path = { seed };
    for ( int d : { 0, 1 } ) for ( NodePath* fork = seed; fork && fork->edges_[d].size() == 1 && fork->edges_[d][0]->edge[d]->edges_[!d].size() == 1; )
    {
        fork = fork->edges_[d][0]->edge[d];
        path.insert( d ? path.end() : path.begin(), fork );
    }
    path_.push_back( new AlleleNode( path ) );
    for ( int d : { 0, 1 } ) while ( extend( d ) );
    
    for ( AlleleNode* an : path_ ) if ( an->homo_ ) for ( NodePath* np : an->path_[0] ) used.insert( np );
    for ( AlleleNode* an : path_ ) for ( AlleleBlock* ab : an->blocks_ ) if ( ab->homo_ ) score_ += Node::getLen( ab->path_[0] );
    int x = 0;
}

AllelePaths::~AllelePaths()
{
    for ( AlleleNode* an : path_ ) delete an;
}

void AllelePaths::align( vector<NodePath*> path[2], int32_t len[2], int32_t ol[2], bool drxn )
{
    unordered_set<NodePath*> pathed[2], base;
    for ( int i : { 0, 1 } ) for ( NodePath* np : path[i] ) pathed[i].insert( np );
    for ( AlleleNode* an : path_ ) for ( int i : { 0, 1 } ) for ( NodePath* np : an->path_[i] ) base.insert( np );
    
    int ext = len[1] > len[0];
    setPathable( path[!ext].back(), pathed[!ext], pathed[ext], base, len[ext]-len[!ext] + 500, drxn );
    if ( pathed[!ext].find( path[ext].back() ) == pathed[!ext].end() ) setPathable( path[ext].back(), pathed[ext], pathed[!ext], base, 500, drxn );
    
    vector<AlleleGraph*> graphs[2];
    AlleleGraph* branch[2]{ new AlleleGraph( path[0][0], pathed, graphs, 0, drxn ), new AlleleGraph( path[1][0], pathed, graphs, 1, drxn ) };
    setAligns( pathed, graphs, drxn );
    for ( int i : { 0, 1 } ) for ( AlleleGraph* ag : graphs[i] ) ag->set();
    
    AlleleAlign aligned( branch, false ), best( branch, true );
    plot( branch[0], 0, aligned, best, drxn );
    
    vector<NodePath*> exts[3];
    if ( best.output( exts, drxn ) )
    {
        path_.insert( drxn ? path_.end() : path_.begin(), new AlleleNode( exts, drxn ? path_.back() : path_[0], drxn ) );
        if ( !exts[2].empty() ) path_.insert( drxn ? path_.end() : path_.begin(), new AlleleNode( exts[2] ) );
    }
    
    for ( int i : { 0, 1 } ) for ( AlleleGraph* ag : graphs[i] ) delete ag;
}

void AllelePaths::align( bool drxn )
{
    if ( ( drxn ? path_.back() : path_[0] )->homo_ || path_.size() < 2 ) return;
    
    AlleleTar tars[2];
    int32_t len[2]{0}, diff[2]{0};
    setQuery( tars, len, diff, drxn );
    setTarget( tars[0], false, !drxn );
    setTarget( tars[1], false, !drxn );
    unordered_map<NodePath*, int> blocked = getMultis();
    for ( int p : { 0, 1 } ) for ( PathEdge* pe : tars[p].path_.back()->edges_[drxn] )
    {
        tars[p].survey( pe->edge[drxn], tars[!p].path_.back(), blocked, diff[p] + pe->diff, max( len[0], len[1] ) - len[p] + 1000 + pe->ol, drxn );
    }
    for ( int p : { 0, 1 } ) if ( tars[p].merged_ )
    {
        tars[p].match( drxn );
        while ( tars[p].path_.back() != tars[!p].path_.back() && tars[p].extend( drxn ) );
        if ( setMerged( tars, drxn ) ) return;
    }
    assert( false );
}

void AllelePaths::challenge( AllelePaths* ap )
{
    vector<NodePath*> branches[2][2][2];
    int len[2][2];
    setBranches( branches[0], len[0] );
    ap->setBranches( branches[1], len[1] );
    
    bool matched = false;
    for ( int d : { 0, 1 } ) if ( len[!d][1] && len[d][0] ) for ( int i : { 0, 1 } ) for ( int j : { 0, 1 } )
    {
        int ii = find( branches[!d][1][i].begin(), branches[!d][1][i].end(), branches[d][0][j][0] ) - branches[!d][1][i].begin() + 1;
        int jj = find( branches[d][0][j].begin(), branches[d][0][j].end(), branches[!d][1][i][0] ) - branches[d][0][j].begin() + 1;
        assert( ii >= 0 && jj >= 0 && ii <= branches[!d][1][i].size()+1 && jj <= branches[d][0][j].size()+1 );
        if ( ( ii >= branches[!d][1][i].size() ) && ( jj >= branches[d][0][j].size() ) ) continue;
        assert( !matched );
        matched = true;
        
        AllelePaths* l = d ? this : ap,* r = d ? ap : this;
        int score[2]{0}, ib = jj >= branches[d][0][j].size() ? ii : 0, jb = ii >= branches[!d][1][i].size() ? jj : 0;
        bool diverged[2]{ false, false }, engulfed[2]{ false, false };
        vector< pair<AlleleNode*, int> > matches[2];
        for ( int k = 0; k < len[d][0]; k++ ) matches[1].push_back( make_pair( r->path_[k], r->path_[k]->homo_ ? 0 : j ) );
        for ( int k = 0; k < len[!d][1]; k++ ) matches[0].push_back( make_pair( l->path_.end()[-k-1], l->path_.end()[-k-1]->homo_ ? 0 : i ) );
        
        if ( jb )
        {
            ii = 0;
            for ( ; jb && ii < branches[!d][1][i].size() && branches[d][0][j][jb-1] == branches[!d][1][i][ii]; ii++ ) jb--;
            if ( ii < branches[!d][1][i].size() ) diverged[0] = true;
        }
        else for ( int k = len[d][0]; !diverged[0] && k < r->path_.size(); k++ )
        {
            int p = !r->path_[k]->path_[1].empty() && r->path_[k]->path_[1][0] == branches[!d][1][i][ii], pp = 0;
            for ( ; pp < r->path_[k]->path_[p].size() && ii < branches[!d][1][i].size() && r->path_[k]->path_[p][pp] == branches[!d][1][i][ii]; pp++ ) ii++;
            if ( k+1 == r->path_.size() && pp ) engulfed[1] = true;
            if ( pp == r->path_[k]->path_[p].size() || ii == branches[!d][1][i].size() ) score[1] += r->path_[k]->getScore();
            else if ( diverged[0] = true ) score[1] += r->path_[k]->getScore( p, pp, 1 );
            matches[1].push_back( make_pair( r->path_[k], p ) );
        }
        
        if ( ib )
        {
            jj = 0;
            for ( ; ib && jj < branches[d][0][j].size() && branches[!d][1][i][ib-1] == branches[d][0][j][jj]; jj++ ) ib--;
            if ( jj < branches[d][0][j].size() ) diverged[1] = true;
        }
        else for ( int k = l->path_.size() - len[!d][1]; !diverged[1] && k-- > 0; )
        {
            int p = !l->path_[k]->path_[1].empty() && l->path_[k]->path_[1].back() == branches[d][0][j][jj], pp = 0;
            for ( ; pp < l->path_[k]->path_[p].size() && jj < branches[d][0][j].size() && l->path_[k]->path_[p].end()[-1-pp] == branches[d][0][j][jj]; pp++ ) jj++;
            if ( !k && pp ) engulfed[0] = true;
            if ( pp == l->path_[k]->path_[p].size() || jj == branches[d][0][j].size() ) score[0] += l->path_[k]->getScore();
            else if ( diverged[1] = true ) score[0] += l->path_[k]->getScore( p, pp, 0 );
            matches[0].push_back( make_pair( l->path_[k], p ) );
        }
        
        for ( int k = l->path_.size() - len[!d][1]; k < l->path_.size(); k++ ) score[0] += l->path_[k]->getScore( i, ib, ii, 1 );
        for ( int k = len[d][0]; k-- > 0; ) score[1] += r->path_[k]->getScore( j, jb, jj, 0 );
        for ( int k : { 0, 1 } ) if ( engulfed[k] && score[k] < score[!k] - min( 500, max( 200, score[k] / 3 ) ) )
        {
            for ( pair<AlleleNode*, int>& match : matches[k] ) match.first->bested_[ match.second ] = true;
            int x = 0;
        }
        int x = 0;
    }
    
    int x = 0;
}

void AllelePaths::create( vector<NodePath*>& paths, NodeRoll& nodes )
{
    vector<AllelePaths*> alleles;
    unordered_set<NodePath*> used;
    for ( NodePath* np : paths ) if ( np->multi_ == 2 && used.insert( np ).second ) alleles.push_back( new AllelePaths( np, used ) );
    sort( alleles.begin(), alleles.end(), []( AllelePaths* a, AllelePaths* b ){ return a->score_ > b->score_; } );
    for ( int again = 1; again-- > 0; ) for ( int i = 0; !again && i < alleles.size(); i++ )
    {
        for ( int d : { 0, 1 } ) alleles[i]->align( d );
        unordered_map<NodePath*, int> multis = alleles[i]->getMultis();
        for ( int j = i+1; j < alleles.size(); j++ ) if ( alleles[j]->discard( alleles[i], multis ) ) alleles.erase( alleles.begin() + j-- );
        alleles[i]->setClones( paths, nodes );
        for ( int d : { 0, 1 } ) if ( again = alleles[i]->merge( alleles, d ) ) break;
        int x = 0;
    }
    assert( false );
//    reduce( alleles );
    assert( false );
}

bool AllelePaths::discard( AllelePaths* ap, unordered_map<NodePath*, int>& multis )
{
    bool unabsorbed = false, absorbed = false;
    for ( AlleleNode* an : path_ ) for ( int p : { 0, 1 } ) for ( NodePath* np : an->path_[p] )
    {
        auto it = multis.find( np );
        if ( it == multis.end() && ( unabsorbed = true ) ) continue;
        
        if ( np->multi_ < it->second + ( an->homo_ ? 2 : 1 ) ) absorbed = true;
    }
    if ( absorbed && unabsorbed && score_ > ap->score_ - 1000 )
    {
        absorbed = false;
        ap->contested_ = contested_ = true;
        assert( false );
    }
    if ( absorbed || !unabsorbed )
    {
        delete this;
        return true;
    }
    return false;
}

bool AllelePaths::dump()
{
    for ( AlleleNode* an : path_ )
    {
        bool dumped = an->bested_[0] || an->bested_[1], kept = !an->homo_ && an != path_[0] && an != path_.back() && ( !an->bested_[0] || an->bested_[1] );
        if ( !dumped || kept ) return false;
    }
    delete this;
    return true;
}

bool AllelePaths::extend( bool drxn )
{
    assert( !path_.empty() );
    if ( !( drxn ? path_.back() : path_[0] )->path_[1].empty() ) return false;
    NodePath* fork = drxn ? path_.back()->path_[0].back() : path_[0]->path_[0][0]; 
    if ( fork->edges_[drxn].size() != 2 || fork->multi_ != 2 ) return false;
    
    vector<NodePath*> paths[3];
    int32_t len[2]{0}, ol[2]{ fork->edges_[drxn][0]->ol, fork->edges_[drxn][1]->ol };
    for ( int i : { 0, 1 } ) for ( PathEdge* pe = fork->edges_[drxn][i]; pe; )
    {
        len[i] += pe->edge[drxn]->size() - pe->ol;
        paths[i].push_back( pe->edge[drxn] );
        pe = pe->edge[drxn]->edges_[drxn].size() == 1 ? pe->edge[drxn]->edges_[drxn][0] : NULL;
    }
    
    while ( !paths[0].empty() && !paths[1].empty() && paths[0].back() == paths[1].back() )
    {
        paths[2].insert( paths[2].begin(), paths[0].back() );
        for ( int i : { 0, 1 } ) paths[i].pop_back();
    }
    
    if ( !drxn ) for ( int i = 0; i < 3; i++ ) reverse( paths[i].begin(), paths[i].end() );
    path_.insert( drxn ? path_.end() : path_.begin(), new AlleleNode( paths, drxn ? path_.back() : path_[0], drxn ) );
    if ( paths[2].empty() ) return false;
    path_.insert( drxn ? path_.end() : path_.begin(), new AlleleNode( paths[2] ) );
    return true;
}

unordered_map<NodePath*, int> AllelePaths::getMultis()
{
    unordered_map<NodePath*, int> multis;
    for ( AlleleNode* an : path_ ) for ( int p : { 0, 1 } ) for ( NodePath* np : an->path_[p] )
    {
        auto ins = multis.insert( make_pair( np, an->homo_ ? 2 : 1 ) );
        if ( !ins.second ) ins.first->second += an->homo_ ? 2 : 1;
    }
    return multis;
}

bool AllelePaths::merge( vector<AllelePaths*>& alleles, bool drxn )
{
    if ( !( drxn ? path_.back() : path_[0] )->homo_ ) return false;
    
    NodePath* fork = drxn ? path_.back()->path_[0].back() : path_[0]->path_[0][0];
    vector<AllelePaths*> candidates;
    for ( AllelePaths* ap : alleles ) if ( ap != this && fork == ( drxn ? ap->path_[0]->path_[0][0] : ap->path_.back()->path_[0].back() ) ) candidates.push_back( ap );
    if ( candidates.empty() ) return false;
    
    AlleleTar tar;
    setTarget( tar, true, drxn );
    unordered_map<NodePath*, int> blocked = getMultis();
    for ( PathEdge* pe : fork->edges_[drxn] ) tar.survey( pe->edge[drxn], NULL, blocked, pe->diff, params.maxMpMean + fork->size(), drxn );
    assert( false );
    
    return true;
}

void AllelePaths::plot( AlleleGraph* a, int i, AlleleAlign& aligned, AlleleAlign& best, bool drxn )
{
    int trim = aligned.lens.size();
    if ( i < a->marks.size() ) for ( pair<AlleleGraph*, int>& mark : a->marks[i] )
    {
        aligned.trim( trim );
        if ( aligned.advance( aligned.lens.size(), mark.first, mark.second ) )
        {
            if ( aligned.add( a, mark.first, i, mark.second, drxn ) ) plot( a, i+1, aligned, best, drxn );
        }
        else if ( !aligned.terminate( a, mark.first, i, mark.second ) )
        {
            AlleleAlign realigned = aligned.branch( mark.first, mark.second );
            if ( realigned.add( a, mark.first, i, mark.second, drxn ) ) plot( a, i+1, realigned, best, drxn );
        }
    }
    else if ( !a->branches.empty() ) for ( pair<AlleleGraph*, int32_t>& branch : a->branches )
    {
        if ( branch.first->multi > aligned.multi( branch.first, 0, false ) )
        {
            aligned.trim( trim );
            plot( branch.first, 0, aligned, best, drxn );
        }
    }
    else best.challenge( aligned );
}

void AllelePaths::reduce( vector<AllelePaths*>& alleles )
{
    for ( int i = 0; i < alleles.size(); i++ ) for ( int j = i+1; j < alleles.size(); j++ ) alleles[i]->challenge( alleles[j] );
    for ( int i = 0; i < alleles.size(); i++ ) if ( alleles[i]->dump() ) alleles.erase( alleles.begin() + i-- );
}

void AllelePaths::setAligns( unordered_set<NodePath*> pathed[2], vector<AlleleGraph*> graphs[2], bool drxn  )
{
    Nodes base, shared;
    unordered_map<Node*, pair<AlleleGraph*, int> > marks;
    
    for ( AlleleGraph* ag : graphs[0] ) for ( Node* node : ag->nodes ) if ( node->cloned_ ) base += node;
    
    for ( AlleleGraph* ag : graphs[1] ) for ( int i = 0; i < ag->nodes.size(); i++ ) if ( ag->paired || ag->nodes[i]->cloned_ )
    {
        bool marked = false;
        if ( ag->nodes[i]->cloned_ ) for ( Node* clone : ag->nodes[i]->cloned_->nodes ) if ( base.find( clone ) && ( marked = true ) ) shared += clone;
        if ( marked ) assert( marks.insert( make_pair( ag->nodes[i], make_pair( ag, ag->scores.size() ) ) ).second );
        if ( marked || ag->paired ) ag->scores.push_back( make_pair( i, 0 ) );
    }
    
    for ( AlleleGraph* ag : graphs[0] ) for ( int i = 0; i < ag->nodes.size(); i++ ) if ( ag->paired || ( ag->nodes[i]->cloned_ && shared.find( ag->nodes[i] ) ) )
    {
        vector< pair<AlleleGraph*, int> > marked;
        if ( ag->paired ) marked.push_back( make_pair( ag->paired, i ) );
        if ( ag->nodes[i]->cloned_ ) for ( Node* clone : ag->nodes[i]->cloned_->nodes )
        {
            auto it = marks.find( clone );
            if ( it != marks.end() ) marked.push_back( it->second );
        }
        assert( !marked.empty() );
        ag->scores.push_back( make_pair( i, 0 ) );
        ag->marks.push_back( marked );
    }
}

void AllelePaths::setBranches( vector<NodePath*> branches[2][2], int len[2] )
{
    len[0] = path_[0]->homo_ ? ( path_.size() < 2 ? 0 : 2 ) : 1;
    len[1] = path_.back()->homo_ ? ( path_.size() < 2 ? 0 : 2 ) : 1;
    
    if ( len[0] ) for ( int i : { 0, 1 } ) branches[0][i] = { path_[ len[0]-1 ]->path_[i].rbegin(), path_[ len[0]-1 ]->path_[i].rend() };
    if ( len[0] == 2 ) for ( int i : { 0, 1 } ) branches[0][i].insert( branches[0][i].end(), path_[0]->path_[0].rbegin(), path_[0]->path_[0].rend() );
    
    if ( len[1] ) for ( int i : { 0, 1 } ) branches[1][i] = path_.end()[ -len[1] ]->path_[i];
    if ( len[1] == 2 ) for ( int i : { 0, 1 } ) branches[1][i].insert( branches[1][i].end(), path_.back()->path_[0].begin(), path_.back()->path_[0].end() );
    
}

void AllelePaths::setClones( vector<NodePath*>& paths, NodeRoll& nodes )
{
    if ( contested_ ) return;
    
    vector<NodePath*> clones;
    NodePath* ends[2]{ path_[0]->path_[0][0], path_.back()->path_[0].back() };
    NodePath* np[2][2]{ { NULL, NULL }, { NULL, NULL } },* cloned[2][2]{ { NULL, NULL }, { NULL, NULL } };
    for ( AlleleNode* an : path_ ) for ( int p : { 0, 1 } ) for ( int i = 0; i < an->path_[p].size(); i++ )
    {
        np[p][1] = an->path_[p][i];
        cloned[p][1] = NULL;
        if ( np[p][1] != ends[0] && np[p][1] != ends[1] && np[p][1]->multi_ > ( 1 + an->homo_ ) )
        {
            cloned[p][1] = new NodePath( np[p][1], paths, nodes, 1+an->homo_ );
            an->path_[p][i] = cloned[p][1];
            clones.push_back( cloned[p][1] );
        }
        
        if ( np[p][0] ) for ( int d : { 0, 1 } ) if ( ( !i && an->homo_ ) || d == p )
        {
            NodePath* edge[2]{ cloned[d][0] ? : np[d][0], cloned[p][1] ? : np[p][1] };
            bool retain[2]{ cloned[d][0] || ends[0] == np[d][0], cloned[p][1] || ends[1] == np[p][1] };
            if ( edge[0] == cloned[d][0] || edge[1] == cloned[p][1] ) for ( int j = 0; j < np[d][0]->edges_[1].size(); j++ )
            {
                if ( np[d][0]->edges_[1][j]->edge[1] != np[p][1] ) continue;

                if ( retain[0] && retain[1] ) new PathEdge( edge[0], edge[1], np[d][0]->edges_[1][j] );
                else if ( retain[0] ) np[d][0]->edges_[1][j--]->claim( edge[0], 0 );
                else if ( retain[1] ) np[d][0]->edges_[1][j]->claim( edge[1], 1 );
            }
        }
        
        for ( int d : { 0, 1 } ) if ( an->homo_ || d == p )
        {
            np[d][0] = np[d][1];
            cloned[d][0] = cloned[d][1];
        }
    }
    
    NodePath::reset( paths );
}

void AllelePaths::setDiffs()
{
    for ( int i : { 0, 1 } ) diffs_[i].clear();
    NodePath* fork[2]{ NULL, NULL };
    unordered_set<NodePath*> pathed, path;
    for ( AlleleNode* an : path_ )
    {
        for ( int i : { 0, 1 } ) for ( NodePath* np : an->path_[i] ) pathed.insert( np );
        if ( an->homo_ )
        {
            if ( !fork[0] ) fork[0] = an->path_[0][0];
            fork[1] = an->path_[0].back();
        }
    }
    for ( int i : { 0, 1 } ) setDiffs( fork[i], 0, path, pathed, i, i );
    for ( int i : { 0, 1 } ) for ( PathEdge* pe : fork[i]->edges_[!i] ) setDiffs( pe->edge[!i], pe->diff, path, pathed, i, !i );
}

void AllelePaths::setDiffs( NodePath* np, int32_t diff, unordered_set<NodePath*>& path, unordered_set<NodePath*>& pathed, int i, bool drxn )
{
    if ( pathed.find( np ) == pathed.end() || !path.insert( np ).second ) return;
    auto it = diffs_[i].insert( make_pair( np, diff ) );
    
    bool complete = false;
    if ( !it.second && !( complete = ( diff <= it.first->second ) ) ) it.first->second = diff;
    if ( !complete ) for ( PathEdge* pe : np->edges_[drxn] ) setDiffs( pe->edge[drxn], diff, path, pathed, i, drxn );
    path.erase( np );
}

bool AllelePaths::setMerged( AlleleTar tar[2], bool drxn )
{
    if ( tar[0].path_.back() == tar[1].path_.back() )
    {
        delete ( drxn ? path_.back() : path_[0] );
        path_.erase( drxn ? path_.end()-1 : path_.begin() );
        int len = 1, maxLen = min( tar[0].path_.size(), tar[1].path_.size() );
        while ( len < maxLen && tar[0].path_.rbegin()[len] == tar[1].path_.rbegin()[len] ) len++;
        
        vector<NodePath*> ext[3];
        for ( int p : { 0, 1 } ) ext[p].insert( ext[p].end(), tar[p].path_.begin(), tar[p].path_.end()-len );
        ext[2].insert( ext[2].end(), tar[0].path_.end()-len, tar[0].path_.end() );
        if ( !drxn ) for ( int i = 0; i < 3; i++ ) reverse( ext[i].begin(), ext[i].end() );
        AlleleNode* bridge = new AlleleNode( ext, drxn ? path_.back() : path_[0], drxn );
        AlleleNode* merge = new AlleleNode( ext[2] );
        path_.insert( drxn ? path_.end() : path_.begin(), bridge );
        path_.insert( drxn ? path_.end() : path_.begin(), merge );
        return true;
    }
    assert( false );
    return false;
}

void AllelePaths::setMultis( unordered_map<NodePath*, int>& multis )
{
    for ( AlleleNode* an : path_ ) for ( int p : { 0, 1 } ) for ( NodePath* np : an->path_[p] )
    {
        auto ins = multis.insert( make_pair( np, an->homo_ ? 2 : 1 ) );
        if ( !ins.second ) ins.first->second += an->homo_ ? 2 : 1;
    }
}

void AllelePaths::setPathable( NodePath* np, unordered_set<NodePath*>& pathed, unordered_set<NodePath*>& merge, unordered_set<NodePath*>& base, int32_t limit, bool drxn )
{
    if ( np->multi_ < 3 && base.find( np ) != base.end() ) return;
    
    bool merged = merge.find( np ) != merge.end();
    for ( PathEdge* pe : np->edges_[drxn] ) if ( ( !merged || merge.find( pe->edge[drxn] ) != merge.end() ) && pathed.insert( pe->edge[drxn] ).second )
    {
        if ( limit > 0 ) setPathable( pe->edge[drxn], pathed, merge, base, limit + pe->ol - pe->edge[drxn]->size(), drxn );
    }
}

void AllelePaths::setQuery( AlleleTar tar[2], int32_t len[2], int32_t diff[2], bool drxn )
{
    NodePath* fork = drxn ? path_.end()[-2]->path_[0].back() : path_[1]->path_[0][0];
    AlleleNode* an = drxn ? path_.back() : path_[0];
    for ( int d : { 0, 1 } ) tar[0].ends_[d] = tar[1].ends_[d] = fork->ends_[d];
    for ( int p : { 0, 1 } ) for ( int i = 0; i < an->path_[p].size(); i++ )
    {
        NodePath* np = drxn ? an->path_[p][i] : an->path_[p].end()[-i-1];
        len[p] += np->size() - ( tar[p].path_.empty() ? fork : tar[p].path_.back() )->getOverlap( np, drxn );
        diff[p] += ( tar[p].path_.empty() ? fork : tar[p].path_.back() )->getOffset( np, drxn );
        tar[p].path_.push_back( np );
        tar[p].insert( np, drxn ? -diff[p] : diff[p], 1 );
    }
}

void AllelePaths::setTarget( AlleleTar& tar, bool homo, bool drxn )
{
    int32_t diff[2]{0};
    NodePath* np[2][2]{ { NULL, NULL }, { NULL, NULL } };
    for ( int i = !homo; i < path_.size(); i++ )
    {
        AlleleNode* an = drxn ? path_[i] : path_.end()[-i-1];
        
        for ( int p : { 0, 1 } ) for ( int j = 0; j < an->path_[p].size(); j++ )
        {
            np[p][1] = drxn ? an->path_[p][j] : an->path_[p].end()[-j-1];
            if ( np[p][0] )
            {
                if ( an->homo_ ) diff[0] = diff[1] = max( diff[0] + np[0][0]->getOffset( np[0][1], drxn ), diff[1] + np[0][0]->getOffset( np[0][1], drxn ) );
                else diff[p] += np[p][0]->getOffset( np[p][1], drxn );
            }
            tar.insert( np[p][1], drxn ? -diff[p] : diff[p], an->homo_ ? 2 : 1 );
            for ( int t : { 0, 1 } ) if ( p == t || an->homo_ ) np[t][0] = np[p][1];
        }
    }
}

void AllelePaths::setTarget( NodePath* np, vector<NodePath*>& pathed, unordered_map<NodePath*, int32_t>& t, int32_t dist, int32_t limit, bool drxn )
{
    if ( find( pathed.begin(), pathed.end(), np ) != pathed.end() ) return;
    
    auto it = t.insert( make_pair( np, dist ) );
    if ( !it.second )
    {
        if ( it.first->second >= dist ) return;
        it.first->second = dist;
    }
    
    limit -= np->size();
    if ( limit < 0 ) return;
    
    pathed.push_back( np );
    for ( PathEdge* pe : np->edges_[drxn] ) setTarget( pe->edge[drxn], pathed, t, dist, limit + pe->ol, drxn );
    pathed.pop_back();
}