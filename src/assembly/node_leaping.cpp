/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "node.h"
#include <algorithm>

bool Node::leap( Querier& bwt, NodeRoll& nodes, NodeList& path, bool drxn )
{
    int32_t ends[2]{ path.back()->ends_[0], path[0]->ends_[1] };
    for ( Node* node : path )
    {
        node->branch_ = true;
        ends[0] = min( ends[0], node->ends_.limits[0][2] );
        ends[1] = max( ends[1], node->ends_.limits[1][2] );
    }
    Node* head = ( drxn ? path.back() : path[0] );
    if ( ends[1] - ends[0] - params.maxPeMean < 0 )
    {
        for ( Node* b : Nodes( ( drxn ? path[0] : path.back() ), params.maxPeMean - ( ends[1] - ends[0] ), !drxn, true, false ).nodes )
        {
            b->branch_ = !b->bad_;
        }
    }
    return ( drxn ? path.back() : path[0] )->leap( bwt, nodes, drxn );
}

bool Node::leap( Querier& bwt, NodeRoll& nodes, bool drxn )
{
    assert( abs( ends_[drxn] - ends_.limits[drxn][0] ) < params.readLen );
    
    Nodes island, ignore;
    for ( Node* n : leapTarget( true, drxn ).nodes ) ignore.fillBad( n, true, !drxn );
    for ( Node* b : Nodes::isBranch( this ).nodes ) for ( Node* n : nodes.nodes ) if ( b->verify( n, b->pe_[drxn], 0, drxn ) && !ignore.find( n )  ) island.add( n );
    
    nodes.test();
    for ( int lept = 1; lept-- > 0; )
    {
        for ( Node* b : Nodes::isBranch( this ).nodes )
        {
            if ( lept ) break;
            for ( Node* n : nodes.nodes ) if ( !b->isPaired( n ) ) b->verify( n, b->pe_[drxn], 0, drxn );
            while ( !b->pe_[drxn].empty() )
            {
                NodeMark mark = b->pe_[drxn].back();
                b->pe_[drxn].pop_back();
                string seq = bwt.getSequence( mark.id );
                if ( drxn ) cout << seq << endl;
                if ( b->add( nodes, seq, mark, drxn ) ) continue;
                Node* node = new Node( seq, mark, drxn );
                nodes += node;
                if ( lept = node->leapSeed( bwt, nodes ) )
                {
                    leapExtend( bwt, node, nodes, drxn );
                    Nodes seeds = Nodes::isBad( node, true, drxn );
                    seeds.fillBad( node, true, !drxn );
                    island += seeds;
                    b->add( node, mark, drxn );
                    break;
                }
                nodes.erase( node );
            }
        }
    }
    
    nodes.test();
    NodeScores scores[2]; // [0] = base, [1] = islands
    Nodes tar = leapTarget( true, drxn );
    for ( Node* f : island.nodes )
    {
        for ( auto &np : f->hits_.pairs[!drxn] )
        {
            if ( !tar.find( np.first ) ) continue;
            scores[0].append( np.first, np.second.count );
            scores[1].append( f, np.second.count );
        }
    }
    NodeIntList hits[2]{ leapScore( scores[0], tar, !drxn ), leapScore( scores[1], island, drxn ) };
    for ( int d : { 0, 1 } )
    {
        Nodes tried;
        for ( int i = 1; i < hits[d].size(); i++ ) if ( hits[d][i].second < hits[d][0].second-3 ) hits[d].erase( hits[d].begin()+i, hits[d].end() );
        for ( int i = 0; i < hits[d].size(); i++ )
        {
            Nodes fwd;
            hits[d][i].first->leapTarget( fwd, 10, false, d ? drxn : !drxn );
            for ( Node* f : fwd.nodes ) if ( tried.add( f ) ) hits[d].insert( hits[d].begin()+i+1, make_pair( f, hits[d][i].second ) );
        }
    }
    
    Node* best[2]{ NULL, NULL };
    int maxOl = 0, maxScore = 0;
    for ( pair<Node*, int>& b : hits[0] )
    {
        for ( pair<Node*, int>& f : hits[1] )
        {
            int ol = mapSeqOverlap( ( drxn ? b.first : f.first )->seq_, ( drxn ? f.first : b.first )->seq_, 25 );
            int score = b.second + f.second + bool(ol) + ( ol > maxOl ? 2 : 0 ) + ( ol > maxOl*2 );
            if ( score <= maxScore ) continue;
            best[0] = b.first;
            best[1] = f.first;
            maxOl = ol;
            maxScore = b.second + f.second + bool(ol);
        }
    }
    assert( maxOl );
    Nodes branches = Nodes::isBranch( this ), tested;
    for ( Node* b : branches.nodes )
    {
        b->clearPaired( true );
        b->branch_ = false;
    }
    best[0]->addEdge( best[1], maxOl, drxn );
    for ( Node* b : branches.nodes ) if ( b->verified_ ) b->setVerified();
    for ( Node* b : branches.nodes )
    {
        bool doVerify = false;
        if ( !b->verified_ ) for ( Edge& e : b->edges_[!drxn] ) if ( e.node->verified_ ) doVerify = true;
        if ( doVerify ) b->verify( tested, drxn );
    }
    nodes.test();
    return true;
}

void Node::leapExtend( Querier& bwt, Node* seed, NodeRoll& nodes, bool drxn )
{
    for ( int d : { drxn, !drxn } )
    {
        int lastExt = 1;
        for ( int again = 1; again-- > 0; )
        {
            NodeRoll ext;
            Nodes tar = leapTarget( true, drxn ), extable, fwd = Nodes::isBad( seed, true, d );
            NodeScores scores;
            for ( Node* f : fwd.nodes ) for ( Node* t : tar.nodes ) if ( !t->isPaired( f ) ) t->verify( f, t->pe_[drxn], 0, drxn );
            for ( Node* f : fwd.nodes ) scores.append( f, f->hits_.get( !d ) );
            NodeIntList hits = leapScore( scores, fwd, !d );
            int dist = d == drxn ? params.readLen * 4 / lastExt : params.readLen * 8 / lastExt;
            for ( int i = 0; i < hits.size() && hits[0].second - hits[i].second < 3; i++ )
            {
                extable.fillIn( hits[i].first, fwd, dist / ( 1 + hits[0].second - hits[i].second ), d, true );
            }
            lastExt = 0;
            for ( Node* f : extable.nodes )
            {
                for ( Edge& e : f->edges_[d] ) if ( !extable.find( e.node ) ) lastExt++;
                if ( f->isContinue( d ) ) ext += f;
            }
            again = !ext.empty();
            lastExt += ext.size();
            for ( Node* node : ext.nodes ) node->extendNode( bwt, nodes, d );
        }
    }
}

NodeIntList Node::leapScore( NodeScores& scores, Nodes& allowed, bool drxn )
{
    NodeIntList hits;
    for ( auto & ns : scores.scores )
    {
        int score = ns.second;
        for ( Node* b : Nodes::inSet( ns.first, allowed, drxn, false ).nodes ) score += scores.get( b );
        hits.push_back( make_pair( ns.first, score ) );
    }
    sort( hits.begin(), hits.end(), []( pair<Node*, int>& a, pair<Node*, int>& b ) { return a.second > b.second; } );
    return hits;
}

bool Node::leapSeed( Querier& bwt, NodeRoll& nodes )
{
    bool good = false;
    for ( int d : { 0, 1 } )
    {
        QueryJunction qj = bwt.mapJunction( seq_, d );
        stop( qj.failure_, d );
        if ( qj.failure_ ) continue;
        good = true;
        addAlts( qj.alts_, nodes, d, drxn_ );
        addExtensions( qj.nodes_, nodes, d, drxn_ );
        extendNode( bwt, nodes, d );
    }
    
    if ( good ) requery( bwt );
    
    return good;
}

Nodes Node::leapTarget( bool resetGood, bool drxn )
{
    Nodes tar( this, params.maxPeMean, !drxn, true );
    tar.dumpBad( true );
    for ( Edge &e : edges_[drxn] ) e.node->leapTarget( tar, 10, resetGood, drxn );
    return tar;
}

void Node::leapTarget( Nodes& tar, int extReadMax, bool resetGood, bool drxn )
{
    if ( tar.find( this ) ) return;
    bool good = resetGood && branch_;
    if ( !good && resetGood ) for ( auto &np : hits_.pairs[!drxn] ) if ( np.first->verified_ ) good = true;
    extReadMax = branch_ ? 10 : extReadMax - countReads( true );
    if ( !good && extReadMax < 0 ) return;
    tar.add( this );
    for ( Edge &e : edges_[drxn] ) e.node->leapTarget( tar, extReadMax, resetGood, drxn );
}
