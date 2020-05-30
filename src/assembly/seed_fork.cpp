/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "seed_fork.h"
#include "node_claim.h"
#include <algorithm>

ExtendBranch::ExtendBranch( Node* node, Nodes& fwd, bool drxn )
: node( node ), score( 0 )
{
    branch.fillIn( node, fwd, drxn, true );
}

ExtendBranch::ExtendBranch( Node* node, bool drxn )
: node( node ), score( 0 )
{
    branch.fill( node, params.maxPeMean - node->size(), drxn, true );
}

ExtendScores::ExtendScores( Node* node, bool drxn )
{
    for ( Node* b : Nodes( node, params.maxPeMean+100 - node->size(), !drxn, true ).nodes ) if ( !b->bad_ ) add( b, drxn );
}

ExtendScores::ExtendScores( Nodes& loop, bool drxn )
{
    Nodes bck = loop;
    for ( Node* node : loop.nodes ) bck.fill( node, max( params.readLen, params.maxPeMean - params.readLen ), !drxn, false );
    for ( Node* node : bck.nodes ) if ( !node->bad_ ) add( node, drxn );
}

void ExtendScores::add( Node* node, bool drxn )
{
    if ( !added.add( node ) ) return;
    for ( auto &np : node->hits_.pairs[drxn] )
    {
        auto it = scores.insert( make_pair( np.first, np.second.count ) );
        if ( !it.second ) it.first->second += np.second.count;
    }
}

vector<ExtendBranch> ExtendScores::branch( Node* node, bool drxn )
{
    vector<ExtendBranch> branches;
    Nodes fwd( node, max( params.readLen, params.maxPeMean - params.readLen ), drxn, false );
    for ( Edge& e : node->edges_[drxn] ) branches.push_back( ExtendBranch( e.node, fwd, drxn ) );
    for ( ExtendBranch& eb : branches ) for ( Node* f : eb.branch.nodes ) eb.score += get( f );
    sort( branches.begin(), branches.end(), []( ExtendBranch& a, ExtendBranch& b ){ return a.score > b.score; } );
    return branches;
}

vector<ExtendBranch> ExtendScores::branch( Nodes& loop, bool drxn )
{
    vector<ExtendBranch> branches;
    Nodes branchable;
    for ( Node* node : loop.nodes ) for ( Edge& e : node->edges_[drxn] ) if ( !loop.find( e.node ) ) branchable += e.node;
    for ( Node* node : branchable.nodes ) branches.push_back( ExtendBranch( node, drxn ) );
    for ( ExtendBranch& eb : branches ) for ( Node* f : eb.branch.nodes ) eb.score += get( f );
    sort( branches.begin(), branches.end(), []( ExtendBranch& a, ExtendBranch& b ){ return a.score > b.score; } );
    return branches;
}

int ExtendScores::get( Node* node )
{
    auto it = scores.find( node );
    return it != scores.end() ? it->second : 0;
}

ExtendAlt::ExtendAlt( Node* branch )
: branch( branch ), lastExt( 0 ), gen( 0 ), stopped( false )
{
    assert( !branch->bad_ );
}

void ExtendAlt::add( vector<ExtendAlt>& alts, Node* node )
{
    if ( node->bad_ ) return;
    for ( ExtendAlt& alt : alts ) if ( alt.branch == node ) return;
    alts.push_back( ExtendAlt( node ) );
}

void ExtendAlt::confirm( vector<ExtendAlt>& alts, NodeRoll& nodes, Nodes& bck, bool drxn )
{
    Nodes fwd[2];
    for ( int i = 0; i < alts.size(); i++ )
    {
        bool bad = true;
        if ( !nodes.find( alts[i].branch ) );
        else if ( bck.find( alts[i].branch ) ) fwd[0].fillIn( alts[i].branch, bck, drxn, true );
        else if ( alts[i].branch->isBlunt( 0, 2, drxn ) ) nodes.erase( alts[i].branch );
        else if ( !alts[i].branch->bad_ )
        {
            if ( bck.empty() ) bad = false;
            for ( Edge& e : alts[i].branch->edges_[!drxn] ) if ( bck.find( e.node ) ) bad = false;
            if ( bad ) fwd[1].fillNot( alts[i].branch, bck, !drxn, false );
        }
        if ( bad ) alts.erase( alts.begin() + i-- );
    }
    
    for ( Node* f : fwd[0].nodes ) for ( Edge& e : f->edges_[drxn] ) if ( !bck.find( e.node ) ) ExtendAlt::add( alts, e.node );
    for ( Node* f : fwd[1].nodes ) for ( Edge& e : f->edges_[!drxn] ) if ( bck.find( e.node ) ) ExtendAlt::add( alts, f );
}

bool ExtendAlt::get( Querier &bwt, NodeRoll& nodes, NodeRoll& ext, bool drxn )
{
    if ( branch->isContinue( drxn ) ) ext += branch;
    if ( stopped || branch->isContinue( drxn ) ) return true;
    
    int32_t cutoff =  params.maxPeMean - branch->size();
    vector< pair<Node*, int> > ends = get( bwt, nodes, cutoff, drxn );
    
    // This is a dead end
    if ( stopped = ends.empty() ) return true;
    
    // This is a splayed branch
    if ( stopped = ( ends.size() > 5 && ends[5].second < ( cutoff / 2 ) ) ) return true;
    
    // Prune if this is an overly substantial branch
    if ( ends[0].second > 2000 ) stopped = true;
    if ( !branch->verified_ && ends[0].second > cutoff && !branch->verify() ) stopped = true;
    if ( stopped ) return true;
    
    for ( int i = 0; i < ends.size(); i++ ) if ( ends[i].second <= 2000 ) ext += ends[i].first;
    
    return true;
}

vector< pair<Node*, int> > ExtendAlt::get( Querier &bwt, NodeRoll &nodes, int32_t cutoff, bool drxn )
{
    vector< pair<Node*, int> > ends;
    
    for ( const pair<Node*, int32_t>& nd : NodeDists( branch, drxn, false ).map )
    {
        if ( nd.first->isContinue( drxn ) ) ends.push_back( make_pair( nd.first, nd.second ) );
    }
    sort( ends.begin(), ends.end(), []( pair<Node*, int>& a, pair<Node*, int>& b ){ return a.second < b.second; } );
    if ( ends.empty() ) return ends;
    
    if ( ends.size() > 5 ) for ( int i = 0; i < ends.size(); i++ ) if ( !ends[i].first->extendable( bwt, drxn ) ) ends.erase( ends.begin() + i-- );
    if ( ends.size() > 2 ) for ( int i = 0; i < ends.size(); i++ ) if ( ends[i].second > cutoff ) ends.erase( ends.begin() + i, ends.end() );
    
    return ends;
}

bool ExtendEdge::get( NodeRoll& nodes, NodeRoll& ext, bool drxn )
{
    for ( int d : { 0, 1 } ) if ( !nodes.find( edge[d] ) ) return false;
//    for ( int d : { 0, 1 } ) if ( !nodes.find( edge[d] ) )
//    {
//        edge[d] = NULL;
//        Coords* coords;
//        for ( Node* node : nodes.nodes ) if ( ( coords = node->getRead( id[d] ) ) && coords->coords[!d] == node->ends_[!d] )
//        {
//            for ( Edge& e : node->edges_[!d] ) if ( ( coords = e.node->getRead( id[!d] ) ) && coords->coords[d] == node->ends_[d] )
//            {
//                edge[d] = node;
//                edge[!d] = e.node;
//                break;
//            }
//            if ( edge[d] ) break;
//        }
//    }
//    if ( !edge[drxn] ) return false;
    
    vector<Node*> adds;
    if ( !get( edge[drxn], adds, 0, drxn ) )
    {
        assert( nodes.find( edge[!drxn] ) && nodes.find( edge[drxn] ) );
        edge[!drxn]->removeEdge( edge[drxn], drxn, true );
        return false;
    }
    for ( Node* node : adds ) ext.add( node );
    
    return !ext.empty();
}

bool ExtendEdge::get( Node* node, vector<Node*>& adds, int readCount, bool drxn )
{
    readCount += node->countReads( true );
    if ( readCount >= 15 ) return true;
    if ( node->isContinue( drxn ) )
    {
        adds.push_back( node );
        return true;
    }
    else if ( node->edges_[drxn].empty() && readCount < 5 ) return false;
    bool good = false;
    for ( Edge& e : node->edges_[drxn] ) if ( get( e.node, adds, readCount, drxn ) ) good = true;
    return good;
}

ExtendFork::ExtendFork( Node* a, Node* b )
: lastExt( 0 ), impatience( 0 ), gen( 0 )
{
    branch[0] = a;
    branch[1] = b;
    loop[0] = loop[1] = NULL;
}

bool ExtendFork::advance( vector<ExtendAlt>& alts, bool drxn )
{
    bool advanced = false;
    Nodes tested, pathed, bck;
    vector< pair<Node*, int> > missed;
    for ( int i : { 0, 1 } ) if ( branch[i] ) branch[i]->verify( tested, drxn );
    for ( int i : { 0, 1 } ) if ( branch[i] )
    {
        ExtendScores scores( branch[i], drxn );
        for ( int again = 1; branch[i] && again-- > 0; )
        {
            vector<ExtendBranch> branches = scores.branch( branch[i], drxn );
            int cutoff = params.branchMinHits / ( advanced ? 1 : 6 );
            if ( again = ( !branches.empty() && ( branches.size() == 1 || branches[0].score >= cutoff ) ) )
            {
                branch[i] = pathed.add( branches[0].node ) ? branches[0].node : NULL;
                for ( int j = 1; j < branches.size(); j++ ) missed.push_back( make_pair( branches[j].node, branches[j].score ) );
                if ( branch[i] ) advanced = true;
            }
        }
    }
    
    if ( !advanced ) return false;
    if ( isLoop( drxn ) ) return true;
    if ( missed.empty() ) return branch[0] || branch[1];
    
    sort( missed.begin(), missed.end(), []( const pair<Node*, int>& a, const pair<Node*, int>& b ){ return a.second > b.second; } );
    assert( branch[0] || branch[1] );
    for ( int i : { 0, 1 } ) if ( !branch[i] && !loop[i] && missed[0].second >= 2 ) branch[i] = missed[0].first;
    for ( pair<Node*, int>& alt : missed ) ExtendAlt::add( alts, alt.first );
    
    return true;
}

bool ExtendFork::get( NodeRoll &nodes, NodeRoll& ext, vector<ExtendAlt>& alts, Nodes& bck, bool priming, bool drxn, bool first )
{
    if ( getLoop( ext, alts, bck, drxn ) ) return true;
    
    NodeRoll extable[2];
    int counts[2]{0};
    for ( int i : { 0, 1 } ) if ( branch[i] && !nodes.find( branch[i] ) ) branch[i] = NULL;
    for ( int i : { 0, 1 } ) if ( branch[i] && ( branch[i]->bad_ || bck.find( branch[i] ) ) ) branch[i] = NULL;
    for ( int i : { 0, 1 } ) if ( branch[i] )
    {
        Nodes fwd( branch[i], params.maxPeMean - params.readLen, drxn, true );
        for ( Node* f : fwd.nodes )
        {
            if ( f->isContinue( drxn ) && !ext.find( f ) ) extable[i] += f;
            else for ( Edge &e : f->edges_[drxn] ) if ( !fwd.find( e.node ) ) counts[i]++;
        }
        if ( extable[i].empty() && !counts[i] ) branch[i] = NULL;
    }
    if ( !branch[0] && !branch[1] ) return false;
    extable[0].add( extable[1] );
    counts[0] += counts[1];
    
    if ( first && ( ( lastExt + counts[0] > 5 ) || ( counts[0] > extable[0].size() ) ) )
    {
        if ( priming ) return true;
        first = advance( alts, drxn );
        return get( nodes, ext, alts, bck, false, drxn, first );
    }
    
    lastExt = max( 0, extable[0].size() + counts[0] - 3 );
    if ( impatience > 14 && priming ) return false;
    if ( impatience + lastExt / 2 > 35 ) return false;
    impatience = max( 0, impatience + extable[0].size() - 3 );
    ext += extable[0];
    for ( Node* node : extable[0].nodes ) bck.fill( node, !drxn, true );
    for ( int i : { 0, 1 } ) if ( branch[i] ) bck.fill( branch[i], !drxn, true );
    
    return true;
}

bool ExtendFork::getLoop( NodeRoll &ext, vector<ExtendAlt>& alts, Nodes& bck, bool drxn )
{
    if ( branch[0] || branch[1] || ( !loop[0] && !loop[1] ) ) return false;
    
    for ( int i : { 0, 1 } ) if ( loop[i] )
    {
        Nodes looped = Nodes::looped( loop[i] );
        if ( looped.empty() )
        {
            assert( false );
            branch[i] = loop[i];
            loop[i] = NULL;
            continue;
        }
        
        ExtendScores scores( looped, drxn );
        vector<ExtendBranch> branches = scores.branch( looped, drxn );
        int cutoff = params.branchMinHits / 2;
        for ( int j = 0; j < branches.size(); j++ )
        {
            if ( j == 0 && branches[0].score >= cutoff ) branch[i] = branches[0].node;
            else if ( j == 1 && branches[1].score >= cutoff && !branch[!i] && !loop[!i] ) branch[!i] = branches[1].node;
            else ExtendAlt::add( alts, branches[i].node );
        }
        if ( branches.empty() || branch[i] ) loop[i] = NULL;
    }
    
    return !branch[0] && !branch[1];
}

bool ExtendFork::isLoop( bool drxn )
{
    for ( int i : { 0, 1 } ) if ( branch[i] && Nodes( branch[i], drxn, false ).find( branch[i] ) )
    {
        loop[i] = branch[i];
        branch[i] = NULL;
    }
    
    return !branch[0] && !branch[1] && ( loop[0] || loop[1] );
}

void SeedExtend::addAlt( Node* node )
{
    ExtendAlt::add( alts_, node );
}

void SeedExtend::addEdge( Node* fork, Node* branch, bool drxn )
{
    ExtendEdge edge;
    edge.edge[!drxn] = fork;
    edge.edge[drxn] = branch;
    edge.id[!drxn] = fork->getTerminalRead( drxn );
    edge.id[drxn] = branch->getTerminalRead( !drxn );
    edges_.push_back( edge );
}

void SeedExtend::addFork( Node* node )
{
    if ( !node ) return;
    for ( ExtendFork& ef : forks_ ) if ( ef.branch[0] == node || ef.branch[1] == node ) return;
    forks_.push_back( ExtendFork( node, NULL ) );
    gen_ = 0;
}

void SeedExtend::cull( Querier& bwt, NodeRoll& nodes, bool drxn )
{
    Nodes forks, good, tested;
    for ( ExtendFork& ef : forks_ ) for ( int i : { 0, 1 } ) if ( ef.branch[i] && !nodes.find( ef.branch[i] ) ) ef.branch[i] = NULL;
    for ( ExtendFork& ef : forks_ ) for ( int i : { 0, 1 } ) if ( ef.branch[i] ) good.fill( ef.branch[i], !drxn, true );
    for ( ExtendFork& ef : forks_ ) for ( int i : { 0, 1 } ) if ( ef.branch[i] ) ef.branch[i]->pruneFwd( bwt, nodes, &good, &tested, drxn );
    
    good.clear();
    for ( ExtendFork& ef : forks_ ) for ( int i : { 0, 1 } ) if ( ef.branch[i] ) forks += ef.branch[i];
    for ( ExtendFork& ef : forks_ ) for ( int i : { 0, 1 } ) if ( ef.branch[i] ) good.fill( ef.branch[i], drxn, true );
    int x = 0;
    for ( Node* node : good.nodes ) if ( !forks.find( node ) ) node->pruneFwd( bwt, nodes, &good, NULL, !drxn );
}

bool SeedExtend::empty()
{
    return forks_.empty() && alts_.empty();
}

bool SeedExtend::extend( Querier& bwt, NodeRoll& nodes, bool priming, bool drxn )
{
    NodeRoll ext;
    Nodes bck;
    cull( bwt, nodes, drxn );
    for ( int i = 0; i < forks_.size(); i++ ) if ( !forks_[i].get( nodes, ext, alts_, bck, priming, drxn ) ) forks_.erase( forks_.begin() + i-- );
    int extCount = ext.size(), altCount = 0;
    
    ExtendAlt::confirm( alts_, nodes, bck, drxn );
    for ( int i = 0; i < alts_.size(); i++ ) if ( !alts_[i].get( bwt, nodes, ext, drxn ) ) alts_.erase( alts_.begin() + i-- );
//    for ( int i = 0; i < edges.size(); i++ ) if ( !edges[i].get( nodes, ext, drxn ) ) edges.erase( edges.begin() + i-- );
    
    if ( !forks_.empty() ) gen_ = 0;
    else if ( ++gen_ >= 5 ) return false;
    
    for ( ExtendAlt& af : alts_ ) if ( !af.stopped ) altCount++;
    cout << ( drxn ? "Right" : "Left" ) << "; Branch extends: " << extCount << ", Alt extends: " << ext.size() - extCount << ", Alt count: " << altCount << ", Alt stopped: " << alts_.size() - altCount << endl;
    
    for ( Node* node : ext.nodes ) node->extendNode( bwt, nodes, drxn );
    
    return !ext.empty();
}

void SeedExtend::reset()
{
    forks_.clear();
    alts_.clear();
    gen_ = 0;
}

