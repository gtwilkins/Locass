/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "node_path.h"
#include "correct_read.h"
#include "seed_fork.h"
#include <algorithm>

void PathScores::add( Node* node )
{
    float coef = 1;
    added += node;
    if ( node->drxn_ < 2 )
    {
        float bad[2]{ 0, 0 };
        for ( auto& np : node->hits_.pairs[!node->drxn_] ) bad[np.first->bad_] += np.second.count;
        if ( bad[1] >= bad[0] ) return;
        coef = ( bad[0] - bad[1] ) / ( bad[0] + bad[1] );
    }
    for ( int d : { 0, 1 } )
    {
        for ( auto& np : node->hits_.pairs[d] ) add( np.first, np.second.count * coef, d );
        add( node, node->hits_.count * coef, d );
    }
}

void PathScores::add( Node* node, float score, bool drxn )
{
    if ( score <= 0 ) return;
    auto it = scores[drxn].insert( make_pair( node, score ) );
    if ( !it.second ) it.first->second += score;
}

float PathScores::get( Node* node, int drxn )
{
    float score = 0;
    for ( int d = drxn == 1; d < min( drxn+1, 2 ); d++ )
    {
        auto it = scores[drxn].find( node );
        if ( it != scores[drxn].end() ) score += it->second;
    }
    return score;
}

float PathScores::get( Nodes& nodes, int drxn )
{
    float score = 0;
    for ( Node* node : nodes.nodes ) score += get( node, drxn );
    return score;
}

PathScores PathScores::unadded()
{
    PathScores r;
    for ( int d : { 0, 1 } ) for ( auto& np : scores[d] ) if ( !added.find( np.first ) ) r.scores[d].insert( np );
    return r;
}

NodeBranch::NodeBranch()
: fork( NULL ), branch( NULL ), score( 0 )
{
    
}

NodeBranch::NodeBranch( Node* fork, Node* branch )
: fork( fork ), branch( branch ), score( 0 )
{
    
}

NodeBranch::NodeBranch( Node* fork, Node* branch, Nodes& fwdIn, NodeScores& scores, bool drxn )
: fork( fork ), branch( branch ), fwd( Nodes::inSet( branch, fwdIn, drxn, true ) ), score( 0 )
{
    for ( Node* f : fwd.nodes ) score += scores.get( f );
}

NodeBranch::NodeBranch( Node* fork, Node* branch, Nodes& fwdIn, PathScores& scores, bool drxn )
: fork( fork ), branch( branch ), fwd( Nodes::inSet( branch, fwdIn, drxn, true ) ), score( 0 )
{
    score = scores.get( fwd, drxn );
}

vector<NodeBranch> NodeBranch::get( Node* fork, Nodes fwd, NodeScores &scores, bool drxn )
{
    vector<NodeBranch> branches;
    if ( !fork ) return branches;
    for ( Edge &e : fork->edges_[drxn] ) if ( fwd.find( e.node ) ) branches.push_back( NodeBranch( fork, e.node, fwd, scores, drxn ) );
    sort( branches );
    return branches;
}

vector<NodeBranch> NodeBranch::get( Node* fork, Nodes& fwd, PathScores& scores, bool drxn )
{
    vector<NodeBranch> branches;
    for ( Edge &e : fork->edges_[drxn] ) if ( fwd.find( e.node ) ) branches.push_back( NodeBranch( fork, e.node, fwd, scores, drxn ) );
    sort( branches );
    return branches;
}

void NodeBranch::sort( vector<NodeBranch> &branches )
{
    std::sort( branches.begin(), branches.end(), []( NodeBranch &a, NodeBranch &b ){ return a.score > b.score; } );
}

PathChunk::PathChunk( NodeList &pilot, Node* alt, int &i )
{
    hits[0] = hits[1] = 0;
    edges[0][0] = edges[0][1] = edges[1][0] = edges[1][1] = false;
    path[0].insert( path[0].end(), pilot.begin() + i, pilot.end() );
    pilot.erase( pilot.begin() + i, pilot.end() );
    i = 0;
    if ( alt ) path[1].push_back( alt );
    else i = 1;
}

void PathChunk::add( vector<PathChunk> &chunks, Node* alt, int &i )
{
    for ( int j = 0; j < path[0].size(); j++ )
    {
        if ( path[0][j] != alt ) continue;
        i = j;
        chunks.push_back( PathChunk( path[0], NULL, i ) );
        return;
    }
    path[1].push_back( alt );
}

NodeBranch PathChunk::alt( NodeScores& scores, Nodes& added, Nodes& btw, int& cumul, bool drxn )
{
    NodeBranch best;
    best.score = 9;
    for ( int i : { 0, 1 } )
    {
        int miss = 0;
        for ( Node* node : path[i] )
        {
            Nodes fwd = Nodes::notSet( node, added, params.maxPeMean, drxn, false );
            fwd -= btw;
            for ( Edge& e : node->edges_[drxn] )
            {
                if ( !fwd.find( e.node ) ) continue;
                NodeBranch alt( node, e.node, fwd, scores, drxn );
                alt.score -= cumul + min( miss, hits[!i] );
                if ( alt.score > best.score )
                {
                    best = alt;
                }
            }
            miss += scores.get( node );
        }
    }
    
    cumul += min( hits[0], hits[1] );
    return best;
}

Nodes PathChunk::between( vector<PathChunk> &chunks )
{
    Nodes nodes[3];
    
    for ( PathChunk& chunk : chunks )
    {
        for ( int i = 0; i < 2 && !chunk.path[i].empty(); i++ )
        {
            nodes[0].fill( chunk.path[i][0], 1, true, false );
            nodes[1].fill( chunk.path[i].back(), 0, true, false );
        }
    }
    for ( Node* node : nodes[0].nodes ) if ( nodes[1].find( node ) ) nodes[2].add( node );
    
    return nodes[2];
}

//void PathChunk::branch( NodeScores& scores, Nodes& added, Nodes& btw, bool drxn )
//{
//    branches[drxn].clear();
//    for ( int i : { 0, 1 } )
//    {
//        for ( Node* node : path[i] ) NodeBranch::add( branches[drxn], node, added, btw, scores, params.maxPeMean, drxn );
//    }
//    NodeBranch::sort( branches[drxn] );
//}

//vector<PathChunk> PathChunk::fill( NodeList &pilot, NodeScores &scores, Nodes &added )
//{
//    Nodes allNodes = Nodes::between( pilot[0], pilot.back(), true ), altNodes;
//    for ( Node* node : allNodes.nodes ) if ( !added.find( node ) ) altNodes.add( node );
//    NodeScores altScores = scores.split( altNodes );
//    
//    Node* cur = pilot[0];
//    int i = 0;
//    vector<PathChunk> chunks( 1, PathChunk( pilot, NULL, i ) );
//    while ( cur )
//    {
//        vector<NodeBranch> branches = NodeBranch::get( cur, allNodes, altScores, 1 );
//        if ( branches.empty() ) break;
//        cur = branches[0].branch;
//        if ( !added.find( cur ) ) scores.add2( cur, true );
//        added.add( cur );
//        for ( auto &np : cur->pairs_ ) if ( altNodes.find( np.first ) ) altScores.add2( np );
//        if ( i )
//        {
//            if ( cur == chunks.back().path[0][i] ) i++;
//            else chunks.push_back( PathChunk( chunks.back().path[0], cur, i ) );
//        }
//        else chunks.back().add( chunks, cur, i );
//    }
//    for ( int i = 0; i < chunks.size(); i++ )
//    {
//        chunks[i].edges[0][0] = chunks[i].edges[0][1] = i > 0;
//        chunks[i].edges[1][0] = chunks[i].edges[1][1] = i < chunks.size()-1;
//        if ( chunks[i].path[1].empty() ) chunks[i].edges[0][1] = chunks[i].edges[1][1] = false;
//        for ( int j : { 0, 1 } ) for ( Node* node : chunks[i].path[j] ) chunks[i].hits[j] += scores.get( node );
//    }
//    
//    return chunks;
//}

bool PathChunk::setFork( Node* branches[2], bool drxn )
{
    for ( int i = 0; i < 2 && !path[i].empty(); i++ )
    {
        Nodes nots;
        for ( int j = 0; j < path[i].size() && !branches[i]; j++ )
        {
            Node* b = drxn ? path[i].end()[-j-1] : path[i][j];
            for ( Node* f : Nodes::notSet( b, nots, drxn, true ).nodes )
            {
                if ( f->isContinue( drxn ) ) branches[i] = b;
                nots.add( f );
            }
        }
    }
    return branches[0] || branches[1];
}

bool PathChunk::take( Node* fork, PathChunk& tar, Nodes &ignore, bool drxn )
{
    edges[!drxn][0] = true;
    for ( int i : { 0, 1 } )
    {
        auto it = find( tar.path[i].begin(), tar.path[i].end(), fork );
        if ( it == tar.path[i].end() ) continue;
        int j = tar.path[!i].empty();
        if ( !j ) path[1].insert( drxn ? path[1].begin() : path[1].end(), tar.path[!i].begin(), tar.path[!i].end() );
        path[j].insert( drxn ? path[j].begin() : path[j].end(), drxn ? it+1 : tar.path[i].begin(), drxn ? tar.path[i].end() : it );
        tar.path[i].erase( drxn ? it+1 : tar.path[i].begin(), drxn ? tar.path[i].end() : it );
        edges[!drxn][1] = true;
        tar.edges[drxn][i] = true;
        if ( !j ) tar.edges[drxn][!i] = true;
        return !j;
    }
    
    int j = tar.hits[1] > tar.hits[0];
    path[1].insert( drxn ? path[1].begin() : path[1].end(), tar.path[j].begin(), tar.path[j].end() );
    
    return true;
}

NodePath::NodePath( Node* seed, PathScores& scores )
{
    path_[0].push_back( seed );
    scores.add( seed );
    
    for ( int again = 1; again-- > 0; )
    {
        while ( NodePath::path( path_[0], scores, 0 ) );
        while ( NodePath::path( path_[0], scores, 1 ) );
        again = unpause( path_[0], scores );
    }
}

NodePath::NodePath( vector<Node*> a, vector<Node*> b )
{
    path_[0] = a;
    path_[1] = b;
}

void NodePath::branch( vector<NodePath>& path, PathScores& scores, bool drxn )
{
    if ( ( drxn ? path.back().path_[0].back() : path[0].path_[0][0] )->cloned_ ) return;
    Nodes btw = Nodes::between( path[0].path_[0][0], path.back().path_[0].back(), true );
    int cumul = 0;
    PathScores altScores = scores.unadded();
    NodeBranch best;
    for ( int i = 0; i < path.size() && cumul < 20; i++ ) ( drxn ? path.end()[-i-1] : path[i] ).branch( best, scores, altScores, btw, cumul, drxn );
    if ( !best.branch ) return;
    vector<Node*> ext;
    while ( path.size() > 1 && !( drxn ? path.back() : path[0] ).branch( best, path, scores, ext, drxn ) );
}

void NodePath::branch( NodeBranch& best, PathScores& scores, PathScores& altScores, Nodes& btw, int& cumul, bool drxn )
{
    float miss[2]{ (float)cumul, (float)cumul };
    for ( int i = 0; i < 2; i++ )
    {
        for ( int j = 0; j < path_[i].size(); j++ )
        {
            Node* node = drxn ? path_[i].end()[-j-1] : path_[i][j];
            for ( Edge& e : node->edges_[drxn] )
            {
                if ( scores.added.find( e.node ) ) continue;
                NodeBranch nb( node, e.node );
                nb.fwd = Nodes::notSet( e.node, scores.added, drxn, true );
                nb.fwd -= btw;
                nb.score = scores.get( nb.fwd, drxn ) - ( path_[!i].empty() ? 0 : miss[i] );
                if ( nb.score >= 10 && nb.score > best.score ) best = nb;
            }
            
            miss[i] += node->hits_.get();
        }
    }
    cumul = min( miss[0], miss[1] );
}

bool NodePath::branch( NodeBranch& best, vector<NodePath>& path, PathScores& scores, vector<Node*>& ext, bool drxn )
{
    for ( int i = 0; i < 2; i++ )
    {
        auto it = find( path_[i].begin(), path_[i].end(), best.fork );
        if ( it == path_[i].end() ) continue;
        scores.add( best.branch );
        best.branch->offset( best.fork, drxn );
        if ( path_[!i].empty() )
        {
            ext.insert( drxn ? ext.end() : ext.begin(), drxn ? it+1 : path_[i].begin(), drxn ? path_[i].end() : it );
            path_[i].erase( drxn ? it+1 : path_[i].begin(), drxn ? path_[i].end() : it );
            NodePath np( ext, vector<Node*>( 1, best.branch ) );
            while( NodePath::path( np.path_[1], scores, drxn ) );
            path.insert( drxn ? path.end() : path.begin(), np );
        }
        else
        {
            path_[i].erase( drxn ? it+1 : path_[i].begin(), drxn ? path_[i].end() : it );
            path_[i].insert( drxn ? path_[i].end() : path_[i].begin(), best.branch );
            path_[!i].insert( drxn ? path_[!i].end() : path_[!i].begin(), ext.begin(), ext.end() );
            if ( !i ) swap( path_[0], path_[1] );
            while( NodePath::path( path_[1], scores, drxn ) );
            assert( false );
        }
        return true;
    }
    
    ext.insert( drxn ? ext.begin() : ext.end(), path_[0].begin(), path_[0].end() );
    path.erase( drxn ? path.end()-1 : path.begin() );
    return false;
}

vector< vector<NodePath> > NodePath::create( NodeRoll& nodes )
{
    vector< vector<NodePath> > paths;
    vector< pair<Node*, int> > seeds;
    for ( Node* node : nodes.getGraph( 2 ).nodes )
    {
        seeds.push_back( make_pair( node, 0 ) );
        Nodes tars[2]{ Nodes( node, params.maxPeMean, 0, true, false ), Nodes( node, params.maxPeMean, 1, true, false ) };
        for ( Node* t : tars[0].nodes ) seeds.back().second += t->hits_.get( tars[1], 1 );
    }
    sort( seeds.begin(), seeds.end(), []( pair<Node*, int> &a, pair<Node*, int> &b ){ return a.second > b.second; } );
    
    Nodes used;
    for ( int i = 0; i < seeds.size(); i++ )
    {
        if ( used.find( seeds[i].first ) ) continue;
        PathScores scores;
        vector<NodePath> path( 1, NodePath( seeds[i].first, scores ) );
        while( path[0].fill( path, scores, 0 ) );
        while( path.back().fill( path, scores, 1 ) );
        NodePath::branch( path, scores, 0 );
        NodePath::branch( path, scores, 1 );
        paths.push_back( path );
        break;
    }
    return paths;
}

bool NodePath::fill( vector<NodePath>& path, PathScores& scores, bool drxn )
{
    assert( !path_[0].empty() && path_[1].empty() );
    Node* cur = NULL;
    if ( !drxn ) for ( int i = path_[0].size(); !cur && --i >= 0; ) if ( path_[0][i]->drxn_ != 1 ) cur = path_[0][i];
    if ( drxn ) for ( int i = 0; !cur && i < path_[0].size(); i++ ) if ( path_[0][i]->drxn_ != 0 ) cur = path_[0][i];
    
    PathScores altScores = scores.unadded();
    Nodes btw = Nodes::between( path_[0][0], path_[0].back(), true );
    vector<Node*> alt, parts[3];
    Node* forks[2]{ NULL, NULL };
    while ( cur && ( !forks[0] || !forks[1] ) )
    {
        vector<NodeBranch> branches = NodeBranch::get( cur, btw, altScores, drxn );
        if ( branches.empty() || ( alt.empty() && !branches[0].score ) ) return false;
        if ( scores.added.find( branches[0].branch ) )
        {
            if ( !alt.empty() ) forks[drxn] = branches[0].branch;
        }
        else
        {
            if ( alt.empty() ) forks[!drxn] = cur;
            branches[0].branch->offset( cur, drxn );
            alt.insert( drxn ? alt.end() : alt.begin(), branches[0].branch );
        }
        cur = branches[0].branch;
    }
    assert( forks[0] && forks[1] );
    for ( Node* node : alt ) scores.add( node );
    for ( Node* node : path_[0] )
    {
        if ( node == forks[1] ) forks[1] = NULL;
        parts[ ( forks[0] ? 0 : ( forks[1] ? 1 : 2 ) ) ].push_back( node );
        if ( node == forks[0] ) forks[0] = NULL;
    }
    path_[0] = parts[int(!drxn)*2];
    path.insert( drxn ? path.end() : path.begin(), NodePath( parts[1], alt ) );
    path.insert( drxn ? path.end() : path.begin(), NodePath( parts[int(drxn)*2], vector<Node*>() ) );
    return true;
}

bool NodePath::path( NodeList& path, PathScores& scores, bool drxn )
{
    assert( !path.empty() );
    Node* cur = drxn ? path.back() : path[0];
    Nodes fwd = Nodes::notSet( cur, scores.added, params.maxPeMean, drxn, false );
    vector<NodeBranch> branches = NodeBranch::get( cur, fwd, scores, drxn );
    if ( branches.empty() || !branches[0].score ) return false;
    branches[0].branch->offset( drxn ? path.back() : path[0] , drxn );
    if ( drxn ) assert( branches[0].branch->ends_[0] < path.back()->ends_[1] && path.back()->ends_[1] < branches[0].branch->ends_[1] );
    if ( !drxn ) assert( path[0]->ends_[0] < branches[0].branch->ends_[1] && branches[0].branch->ends_[0] < path[0]->ends_[0] );
    path.insert( drxn ? path.end() : path.begin(), branches[0].branch );
    scores.add( branches[0].branch );
    return branches[0].branch->verified_;
}

bool NodePath::unpause( NodeList& path, PathScores& scores )
{
    if ( !path[0]->verified_ && !path.back()->verified_ ) return false;
    assert( false ); 
    return true;
}
