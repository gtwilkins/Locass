/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "seed_fork.h"
#include "node_claim.h"
#include <algorithm>

AltFork::AltFork( Node* branch )
: branch( branch ), lastExt( 0 ), gen( 0 ), stopped( false )
{
    assert( !branch->bad_ );
}

void AltFork::add( vector<AltFork>& alts, Node* node )
{
    if ( node->bad_ ) return;
    for ( AltFork& alt : alts ) if ( alt.branch == node ) return;
    alts.push_back( AltFork( node ) );
}

void AltFork::confirm( vector<AltFork>& alts, NodeRoll& nodes, Nodes& bck, bool drxn )
{
    Nodes fwd[2];
    for ( int i = 0; i < alts.size(); i++ )
    {
        bool bad = true;
        if ( bck.find( alts[i].branch ) ) fwd[0].fillIn( alts[i].branch, bck, drxn, true );
        else if ( alts[i].branch->isBlunt( 0, 2 , drxn ) ) nodes.erase( alts[i].branch );
        else if ( !alts[i].branch->bad_ )
        {
            if ( bck.empty() ) bad = false;
            for ( Edge& e : alts[i].branch->edges_[!drxn] ) if ( bck.find( e.node ) ) bad = false;
            if ( bad ) fwd[1].fillNot( alts[i].branch, bck, !drxn, false );
        }
        if ( bad ) alts.erase( alts.begin() + i-- );
    }
    
    for ( Node* f : fwd[0].nodes ) for ( Edge& e : f->edges_[drxn] ) if ( !bck.find( e.node ) ) AltFork::add( alts, e.node );
    for ( Node* f : fwd[1].nodes ) for ( Edge& e : f->edges_[!drxn] ) if ( bck.find( e.node ) ) AltFork::add( alts, f );
    
//    for ( AltFork& af : alts )
//    {
//        if ( af.stopped || af.branch->isContinue( drxn ) ) continue;
//        for ( Edge& e : af.branch->edges_[!drxn] ) if ( e.node->isForkSplayed( params.readLen, 12, drxn ) )
//        {
//            af.stopped = true;
//        }
//    }
}

bool AltFork::get( Querier &bwt, NodeRoll &nodes, NodeRoll& ext, bool drxn )
{
    if ( branch->isContinue( drxn ) ) ext += branch;
    if ( stopped || branch->isContinue( drxn ) ) return true;
    
    int32_t cutoff =  params.maxPeMean - branch->size();
    vector< pair<Node*, int> > ends = get( bwt, nodes, cutoff, drxn );
    
    // This is a dead end
    if ( stopped = ends.empty() )
    {
        int x = 0;
        return true;
    }
    
    // This is a splayed branch
    if ( stopped = ( ends.size() > 5 && ends[5].second < ( cutoff / 2 ) ) )
    {
        int x = 0;
        return true;
    }
    
    // Prune if this is an overly substantial branch
    if ( ends[0].second > 2000 ) stopped = true;
    if ( !branch->verified_ && ends[0].second > cutoff && !branch->verify() )
    {
        stopped = true;
    }
    if ( stopped ) return true;
    
    for ( int i = 0; i < ends.size(); i++ ) if ( ends[i].second <= 2000 ) ext += ends[i].first;
    
    return true;
}

vector< pair<Node*, int> > AltFork::get( Querier &bwt, NodeRoll &nodes, int32_t cutoff, bool drxn )
{
    vector< pair<Node*, int> > ends;
    for ( const pair<Node*, NodeOffset>& no : NodeOffsets( branch, drxn, false, false ).map )
    {
        if ( no.first->isContinue( drxn ) ) ends.push_back( make_pair( no.first, abs( no.second.dists[0] ) ) );
    }
    sort( ends.begin(), ends.end(), []( pair<Node*, int>& a, pair<Node*, int>& b ){ return a.second < b.second; } );
    if ( ends.empty() ) return ends;
    
    if ( ends.size() > 5 ) for ( int i = 0; i < ends.size(); i++ ) if ( !ends[i].first->extendable( bwt, drxn ) ) ends.erase( ends.begin() + i-- );
    if ( ends.size() > 2 ) for ( int i = 0; i < ends.size(); i++ ) if ( ends[i].second > cutoff ) ends.erase( ends.begin() + i, ends.end() );
    
    return ends;
}

//bool AltFork::get( Querier &bwt, NodeRoll &nodes, NodeRoll& ext, bool drxn )
//{
//    if ( stopped ) return true;
//    
//    // Ensure that this branch is verified or attempt to prune
//    if ( !branch->verified_ )
//    {
//        Nodes tested;
//        branch->verify( tested, drxn );
//        if ( !branch->verified_ )
//        {
//            branch->sortEdges( !drxn );
//            for ( int i = branch->edges_[!drxn].size(); --i >= 0 && branch->edges_[!drxn].size() > 1; )
//            {
//                if ( !branch->edges_[!drxn][i].node->bad_ ) branch->pruneEdges( bwt, nodes, branch->edges_[!drxn][i], !drxn );
//            }
//            return false;
//        }
//    }
//    
//    // Ensure that this branch is not a dead end
//    int count = 0;
//    Nodes fwd( branch, drxn, true );
//    for ( Node* f : fwd.nodes ) if ( f->isContinue( drxn ) ) count++;
//    if ( !count ) return false;
//    
//    if ( ( min( lastExt, count ) / 2 ) * gen++ < 15 )
//    {
//        for ( Node* f : fwd.nodes ) if ( f->isContinue( drxn ) ) ext.add( f );
//        lastExt = count;
//    }
//    else
//    {
//        Nodes tested;
//        branch->verify( tested, drxn );
//        stopped = true;
//    }
//    return true;
//}

//bool AltFork::prune( Querier &bwt, NodeRoll &nodes, bool drxn )
//{
//    assert( false );
////    branch->extendFork( bwt, nodes, max( params.readLen, params.maxPeMean - params.readLen ), 6, !drxn );
////    branch->verifyFork( params.maxPeMean, false, drxn );
////    ClaimFork( branch, nodes, false, !drxn );
////    return branch->bad_;
//}

SeedFork::SeedFork( Node* a, Node* b )
: lastExt( 0 ), impatience( 0 ), gen( 0 )
{
    branch[0] = a;
    branch[1] = b;
}

bool SeedFork::advance( bool drxn )
{
    Nodes tested, bck;
    for ( int i = !branch[0]; i < 2 && branch[i]; i++ )
    {
        branch[i]->verify( tested, drxn );
        bck.fill( branch[i], params.maxPeMean - branch[i]->size() + params.readLen*2, !drxn, true, false );
    }
    NodeScores scores;
    for ( Node* b : bck.nodes ) scores.add( b, drxn, false );
    
    bool advanced = false;
    for ( int again = 1; again-- > 0; )
    {
        for ( int i = !branch[0]; i < 2 && branch[i]; i++ )
        {
            // Advance fork until branched
            while ( !branch[!i] && advance( branch[i], branch[!i], scores, advanced, drxn ) );
            
            // Advance branch until fork
            while ( branch[!i] && !Nodes( branch[!i], drxn, false, false ).find( branch[i] ) && advance( branch[i], branch[!i], scores, advanced, drxn ) ) again = 1;
        }
        if ( again ) advanced = true;
    }
    return advanced;
}

bool SeedFork::advance( Node*& node, Node*& alt, NodeScores &scores, bool advanced, bool drxn )
{
    if ( !node->verified_ ) return false;
    vector<NodeBranch> branches = NodeBranch::get( node, Nodes( node, params.avgPeMean, drxn, false, false ), scores, drxn );
    int cutoff = params.branchMinHits / ( advanced ? 1 : 6 );
    if ( branches.empty() || ( branches.size() > 1 && branches[0].score < cutoff ) ) return false;
    if ( node == alt ) alt = NULL;
    cutoff = max( 2, int( params.branchMinHits / 6 ) );
    for ( int i = 1; i < branches.size(); i++ )
    {
        if ( !alt && branches[i].score > cutoff )
        {
            alt = branches[i].branch;
            scores.add( alt, drxn, false );
        }
        else AltFork::add( alts, branches[i].branch );
    }
    
    scores.erase( node );
    node = branches[0].branch;
    assert( !node->bad_ );
    if ( node != alt ) scores.add( node, drxn, false );
    impatience = 0;
    return true;
}

//bool SeedFork::cull( Querier &bwt, NodeRoll &nodes, Nodes& bck, bool priming, bool drxn )
//{
//    Nodes fwd[2];
//    for ( int i = 0; i < alts.size(); i++ )
//    {
//        bool bad = true;
//        if ( bck.find( alts[i].branch ) ) fwd[0].fillIn( alts[i].branch, bck, drxn, true );
//        else if ( alts[i].branch->isBlunt( 0, 2 , drxn ) ) nodes.erase( alts[i].branch );
//        else if ( !alts[i].branch->bad_ )
//        {
//            for ( Edge& e : alts[i].branch->edges_[!drxn] ) if ( bck.find( e.node ) ) bad = false;
//            if ( bad ) fwd[1].fillNot( alts[i].branch, bck, !drxn, false );
//        }
//        alts.erase( alts.begin() + i-- );
//    }
//    
//    for ( Node* f : fwd[0].nodes ) for ( Edge& e : f->edges_[drxn] ) if ( !bck.find( e.node ) ) AltFork::add( alts, e.node );
//    for ( Node* f : fwd[1].nodes ) for ( Edge& e : f->edges_[!drxn] ) if ( bck.find( e.node ) ) AltFork::add( alts, f );
//    
//    for ( AltFork& af : alts )
//    {
//        if ( af.stopped || af.branch->isContinue( drxn ) ) continue;
//        for ( Edge& e : af.branch->edges_[!drxn] ) if ( e.node->isForkSplayed( params.readLen, 12, drxn ) ) af.stopped = true;
//    }
//}

bool SeedFork::extend( vector<SeedFork>& forks, Querier& bwt, NodeRoll& nodes, bool priming, bool drxn )
{
    NodeRoll ext;
    Nodes bck;
    nodes.test();
    for ( int i = 0; i < forks.size(); i++ ) if ( !forks[i].getExt( ext, bck, priming, drxn ) ) forks.erase( forks.begin() + i-- );
    for ( int i = 0; i < forks.size(); i++ ) if ( !forks[i].getAlt( bwt, nodes, ext, bck, drxn ) ) forks.erase( forks.begin() + i-- );
    
    for ( Node* node : ext.nodes ) node->extendNode( bwt, nodes, drxn );
    return !ext.empty();
}

//bool SeedFork::getAlts( Querier& bwt, NodeRoll& nodes, NodeRoll& ext, Nodes& bck, bool drxn )
//{
//    if ( !bck.empty() ) AltFork::confirm( alts, nodes, bck, drxn );
//}

bool SeedFork::getAlt( Querier& bwt, NodeRoll& nodes, NodeRoll& ext, Nodes& bck, bool drxn )
{
    int extCount = ext.size(), altCount = 0;
    AltFork::confirm( alts, nodes, bck, drxn );
    for ( int i = 0; i < alts.size(); i++ ) if ( !alts[i].get( bwt, nodes, ext, drxn ) ) alts.erase( alts.begin() + i-- );
    
    if ( branch[0] || branch[1] ) gen = 0;
    else if ( ++gen >= 5 ) return false;
    
    for ( AltFork& af : alts ) if ( !af.stopped ) altCount++;
    cout << ( drxn ? "Right" : "Left" ) << "; Branch extends: " << extCount << ", Alt extends: " << ext.size() - extCount << ", Alt count: " << altCount << ", Alt stopped: " << alts.size() - altCount << endl;
    
    for ( AltFork& af : alts ) if ( !af.stopped ) return true;
    return branch[0] || branch[1];
}

//Nodes SeedFork::getBack( vector<SeedFork> &forks, bool drxn )
//{
//    Nodes fwd, bck;
//    for ( SeedFork& sf : forks ) for ( int d : { 0, 1 } ) if ( sf.branch[d] ) fwd.fill( sf.branch[d], drxn, true, false );
//    for ( Node* f : fwd.nodes ) if ( f->edges_[drxn].empty() ) bck.fill( f, !drxn, true, true );
//    return bck;
//}

bool SeedFork::getExt( NodeRoll& ext, Nodes& bck, bool priming, bool drxn, bool first )
{
    NodeRoll extable[2];
    int counts[2]{0};
    for ( int i = !branch[0]; i < 2 && branch[i]; i++ ) if ( bck.find( branch[i] ) ) branch[i] = NULL;
    for ( int i = !branch[0]; i < 2 && branch[i]; i++ )
    {
        Nodes fwd( branch[i], params.maxPeMean - params.readLen, drxn, true, false );
        for ( Node* f : fwd.nodes )
        {
            if ( f->isContinue( drxn ) && !ext.find( f ) ) extable[i] += f;
            else for ( Edge &e : f->edges_[drxn] ) if ( !fwd.find( e.node ) ) counts[i]++;
        }
        if ( extable[i].empty() && !counts[i] ) branch[i] = NULL;
    }
    extable[0].add( extable[1] );
    counts[0] += counts[1] + extable[0].size();
    
    if ( first && ( ( lastExt + counts[0] > 5 ) || ( counts[0] > extable[0].size() ) ) )
    {
        if ( priming ) return true;
        first = advance( drxn );
        return getExt( ext, bck, false, drxn, first );
    }
    
    lastExt = max( 0, extable[0].size() + counts[0] - 3 );
    if ( impatience > 14 && priming ) return false;
    if ( !branch[0] && !branch[1] && alts.empty() ) return false;
    if ( impatience + lastExt / 2 > 35 ) assert( false );
    impatience = max( 0, impatience + extable[0].size() - 3 );
    ext += extable[0];
    for ( Node* node : extable[0].nodes ) bck.fill( node, !drxn, true, true );
    for ( int i = 0; i < 2; i++ ) if ( branch[i] ) bck.fill( branch[i], !drxn, true, true );
    
    return true;
}

//bool SeedFork::getSide( NodeRoll& ext, Nodes& bck, bool drxn )
//{
//    Nodes tested;
//    for ( int i = 0; i < alts.size(); i++ )
//    {
//        if ( !bck.find( alts[i].branch ) ) continue;
//        alts[i].branch->verify( tested, drxn );
//        for ( Node* f : Nodes::inSet( alts[i].branch, bck, drxn, true ).nodes ) for ( Edge& e : f->edges_[drxn] ) if ( !bck.find( e.node ) ) sides.add( e.node );
//        alts.erase( alts.begin() + i-- );
//    }
//    for ( Node* node : sides.nodes ) assert( !node->bad_ );
//    for ( int i = 0; i < sides.size(); i++ )
//    {
//        int count = 0;
//        
//        // Advance this side if it has merged back into the fork
//        if ( bck.find( sides[i] ) )
//        {
//            sides[i]->verify( tested, drxn );
//            for ( Node* f : Nodes::inSet( sides[i], bck, drxn, true ).nodes ) for ( Edge& e : f->edges_[drxn] ) if ( !bck.find( e.node ) ) sides.add( e.node );
//        }
//        // Count how many open branches from this side
//        else for ( Node* f : Nodes( sides[i], drxn, true ).nodes ) if ( f->isContinue( drxn ) ) count++;
//        
//        // Extend this side far enough that it may be scrutinised
//        if ( count )
//        {
//            int dist = ( params.maxPeMean * 3 + params.readLen - sides[i]->size() ) / count;
//            count = 0;
//            for ( Node* f : Nodes( sides[i], dist, drxn, true ).nodes )
//            {
//                if ( !f->isContinue( drxn ) ) continue;
//                ext.add( f );
//                count++;
//            }
//            
//            // Ready to be scrutinised
//            if ( !count ) AltFork::add( alts, sides[i] );
//        }
//        
//        // Remove blunt side or one ready to be scrutinised
//        if ( !count ) sides -= sides[i--];
//    }
//    for ( Node* node : sides.nodes ) assert( !node->bad_ );
//    
//    if ( !branch[0] && !branch[1] && sides.empty() && alts.empty() )
//    {
//        return false;
//    }
//    return true;
//}

bool SeedFork::restart( Querier& bwt, NodeRoll& nodes, vector<Node*> path[2], bool drxn )
{
    bool ended[2]{ !path[0].empty(), !path[1].empty() };
    for ( int i = 0; i < 2; i++ )
    {
        int cumul = 0;
        for ( int j = 0; j < path[i].size() && ended[i] && cumul < 5; j++ )
        {
            Node* node = path[i][( drxn ? path[i].size()-j-1 : j )];
            assert( !node->cloned_ || !node->edges_[drxn].empty() );
            for ( Node* f : Nodes( node, drxn, true, false ).nodes )
            {
                if ( f->isContinue( drxn ) ) ended[i] = false;
                assert( !f->cloned_ );
            }
            if ( !branch[i] || !ended[i] ) branch[i] = node;
            cumul += node->hits_.get();
        }
        assert( path[i].empty() || branch[i] );
    }
    for ( int i = !ended[0]; i < 2 && ended[i]; i++ )
    {
        if ( !ended[!i] && branch[!i] && ( drxn ? branch[!i]->ends_[1] - params.maxPeMean < branch[i]->ends_[1]
                                                : branch[i]->ends_[0] < branch[!i]->ends_[0] + params.maxPeMean ) ) return true;
        if ( Node::leap( bwt, nodes, path[i], drxn ) ) return true;
        assert( false );
    }
    if ( !ended[0] && !ended[1] ) return true;
    return false;
}

//void SeedFork::test( Querier &bwt, NodeRoll &nodes, bool drxn )
//{
//    for ( int i = 0; i < alts.size(); i++ ) if ( !alts[i].test( ) ) alts.erase( alts.begin() + i-- );
//}
