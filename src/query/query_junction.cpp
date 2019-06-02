/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "query_junction.h"
#include "parameters.h"
#include <cassert>

QueryJunction::QueryJunction( IndexReader* ir, QueryBinaries* qb, string &seq, ReadId* counts, int* ols, int cutoff, bool drxn )
: seq_( seq ), failure_( 0 ), drxn_( drxn ), branch_( false )
{
    QState qs;
    ir->primeOverlap( seq_, qs.q, qs.rank, qs.count, qs.ol, drxn );
    query( ir, qs, counts, ols );
    if ( failure_ = qs.failure() ) return;
//    nodes_ = QueryNode::graph( qb, qs, "", drxn );
    vector<QueryNode*> nodes;
    for ( QueryNode* node : QueryNode::graph( qb, qs, "", drxn ) ) if ( node->confirm() ) nodes.push_back( node );
    int maxExt = 0, maxReads = 0;
    for ( QueryNode* node : nodes ) maxExt = max( maxExt, node->setExt() );
    for ( QueryNode* node : nodes ) maxReads = max( maxReads, node->setReads( maxExt - cutoff ) );
    for ( QueryNode* node : nodes ) if ( !add( node, maxExt > cutoff && maxReads >= 10, true ) ) delete node;
    assert( !nodes_.empty() );
//    for ( QueryNode* node : nodes ) node->setSeq( drxn );
//    if ( maxExt > cutoff ) QueryNode::cull( nodes_, maxExt - cutoff );
//    for ( QueryNode* node : nodes ) node->setSeq( drxn );
//    setAlts();
}

QueryJunction::~QueryJunction()
{
//    for ( QueryNode* node : nodes_ ) for ( int i = 0; i < node->edges[0].size(); i++ ) if ( !node->edges[0][i] ) node->edges[0].erase( node->edges[0].begin() + i-- );
//    for ( QueryNode* node : nodes_ ) if ( node->edges[0].empty() ) node->discard();
    for ( QueryNode* node : nodes_ ) delete node;
    for ( QueryNode* node : alts_ ) delete node;
}

bool QueryJunction::add( QueryNode* node, bool cull, bool first )
{
    if ( cull && node->readCount < 2 ) return false;
    while ( node->merge() );
    assert( node->setSeq( drxn_ ) );
    if ( first ) for ( QueryNode* e : node->edges[0] ) if ( e->len ) first = false;
    if ( first ) ( node->len ? nodes_ : alts_ ).push_back( node );
    for ( QueryNode* e : node->edges[1] ) if ( !e->fixed ) add( e, false, !node->len && e->len );
    return true;
}

uint8_t QueryJunction::getChar( int ol )
{
    if ( ol >= seq_.size() ) return 4;
    return ( drxn_ ? charToInt[ seq_.end()[-ol-1] ] : charToIntComp[ seq_[ol] ] );
}

//void QueryJunction::setAlts()
//{
//    for ( int i = 0; i < nodes_.size(); i++ )
//    {
//        if ( nodes_[i]->len == nodes_[i]->seq.size() ) continue;
//        alts_.push_back( nodes_[i] );
//        nodes_.erase( nodes_.begin() + i-- );
//    }
//    
//    if ( nodes_.empty() ) for ( QueryNode* alt : alts_ ) alt->setNonAlt( nodes_ );
//    
//    if ( nodes_.empty() && !failure_ ) failure_ = 3;
//    
//    for ( QueryNode* node : alts_ )
//    {
//        assert( !seq_.empty() );
//        assert( !node->len );
//    }
//}

void QueryJunction::query( IndexReader* ir, QState &qs, ReadId* counts, int* ols )
{
    while ( qs.count && qs.edges.empty() && !failure_ )
    {
        ir->countRange( qs.q.back(), qs.rank, qs.count, qs.ranks, qs.counts );
        
        if ( branch_ ) qs.record();
        else if ( qs.ol >= ols[0] ) branch_ = qs.count <= counts[qs.ol];
        else if ( qs.ol >= ols[1] ) failure_ = 2;
        
        uint8_t i = qs.perfect ? getChar( qs.ol ) : 4;
        if ( !branch_ && i < 4 && !qs.counts[i] ) failure_ = failure_ ? : 1;
        if ( failure_ ) return;
        
        CharId branchMin = branch_ ? qs.perfect && i > 3 : 0;
        if ( branch_ && !branchMin )
        {
            int branchCount[2]{ qs.counts[0] > 0, qs.counts[0] > 1 }, iMax = 0;
            for ( int j = 1; j < 4; j++ )
            {
                if ( !qs.counts[j] ) continue;
                if ( qs.counts[j] > qs.counts[iMax] ) iMax = j; 
                if ( qs.counts[j] > 1 ) branchCount[1]++;
                branchCount[0]++;
            }
            
            bool strong = qs.counts[iMax] >= 15 && qs.ol < ols[1];
            if ( i < 4 && i != iMax ) branchMin = 1;
            else branchMin = branchCount[strong] == 1 ? 0 : strong + 1;
            if ( !branchMin ) i = iMax;
        }
        
        if ( branchMin ) qs.branch( i, branchMin );
        else if ( !qs.advance( i ) ) assert( false );
//        else if ( !qs.advance( i ) ) failure_ = 1;
    }
    
    for ( QState &qsEdge : qs.edges ) if ( qs.gen < 3 ) query( ir, qsEdge, counts, ols );
}