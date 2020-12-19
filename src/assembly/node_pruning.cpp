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

#include "node.h"
#include "node_claim.h"
#include "locus_fill.h"
#include "prune_bubble.h"
#include <algorithm>

bool Node::isBlunt( int readCount, int readLimit, bool drxn )
{
    readCount += countReads( true );
    if ( readCount > readLimit || verified_ || cloned_ ) return false;
    for ( Edge& e : edges_[drxn] ) if ( !e.node->isBlunt( readCount, readLimit, drxn ) ) return false;
    return !edges_[drxn].empty() || stop_[drxn] == DEAD_END || stop_[drxn] == SHORT_END || stop_[drxn] == BLUNT_END;
}

bool Node::isBranchComplete( int readCount, int readMin, bool drxn )
{
    readCount += countReads( true );
    if ( readCount >= readMin ) return true;
    for ( Edge& e : edges_[drxn] ) if ( !e.node->isBranchComplete( readCount, readMin, drxn ) ) return false;
    return !isContinue( drxn );
}

bool Node::isForkComplete( int32_t dist, int readMin, bool drxn )
{
    bool good = true;
    if ( dist <= 0 ) return true;
    for ( Node* f : Nodes( this, dist, drxn, true ).nodes ) if ( f->isContinue( drxn ) ) good = false;
    if ( !good && readMin ) for ( Edge& e : edges_[drxn] ) if ( !e.node->isBranchComplete( 0, readMin, drxn ) ) return false;
    return good || readMin;
}

//bool Node::isForkSplayed( int32_t dist, int maxBranches, bool drxn )
//{
//    Nodesx fwd( this, dist, drxn, false, false );
//    int branches = 0;
//    for ( Node* f : fwd.nodes )
//    {
//        bool ended = !f->stop_[drxn];
//        for ( Edge& e : f->edges_[drxn] ) if ( fwd.find( e.node ) ) ended = false;
//        if ( ended ) branches++;
//    }
//    return branches <= maxBranches;
//}

bool Node::isMergeable( Node* alt, vector<Node*> merges[2], bool drxn )
{
    if ( !cloned_ ) return !alt->cloned_;
    Nodes pairs;
    for ( Node* clone : cloned_->nodes )
    {
        if ( clone->edges_[drxn].size() > 1 ) return false;
        Node* edge = clone->edges_[drxn].empty() ? NULL : clone->edges_[drxn][0].node;
        merges[0].push_back( clone );
        merges[1].push_back( edge );
        if ( !pairs.add( edge ) );
        if ( !alt->isClone( edge ) || edge->edges_[!drxn].size() != 1 ) return false;
    }
    return !alt->cloned_ || alt->cloned_->size() == pairs.size();
}

bool Node::isSubstantial( int readCount, int readLimit, bool drxn )
{
    readCount += countReads( true );
    if ( readCount >= readLimit ) return true;
    for ( Edge& e : edges_[drxn] ) if ( e.node->isSubstantial( readCount, readLimit, drxn ) ) return true;
    return false;
}

bool Node::merge( NodeRoll &nodes )
{
    updateStates( nodes );
    int merged = 0;
    Nodes tested;
    nodes.test( true );
    for ( Node* node : nodes.nodes ) node->readTest();
    for ( Node* node : nodes.getGraph( 2 ).nodes )
    {
        if ( !nodes.find( node ) || !tested.add( node ) ) continue;
        for ( int d : { 0, 1 } ) while ( node->merge( nodes, d ) ) merged++;
        node->merge( nodes, tested, merged );
    }
    nodes.test( true );
    if ( merged ) Node::verify( nodes );
    cout << "Merged: " << merged << endl;
    return merged;
}

void Node::merge( NodeRoll &nodes, Nodes& tested, int& merged )
{
    for ( int d : { 0, 1 } )
    {
        vector<Edge> nxt;
        for ( Edge& e : edges_[d] ) if ( tested.add( e.node ) ) nxt.push_back( e );
        for ( Edge& e : nxt ) while ( e.node->merge( nodes, d ) ) merged++;
        for ( Edge& e : nxt ) e.node->merge( nodes, tested, merged );
    }
}

bool Node::merge( NodeRoll& nodes, bool drxn, bool verify )
{
    assert( !verify );
    if ( edges_[drxn].size() != 1 || edges_[drxn][0].node->edges_[!drxn].size() != 1 ) return false;
    if ( edges_[drxn][0].leap || ( edges_[drxn][0].leap = edges_[drxn][0].ol <= 0 ) ) return false;
    
    Node* alt = edges_[drxn][0].node;
    int ol = edges_[drxn][0].ol;
    vector<Node*> merges[2]{ { this }, { alt } };
    if ( alt->drxn_ == 2 && drxn_ < 2 ) return false;
    if ( !isMergeable( alt, merges, drxn ) ) return false;
    for ( int i = 0; i < merges[0].size(); i++ ) for ( int d = 0; d < 2 && merges[d][i]; d++ ) merges[d][i]->clearPaired( false );
    
    int32_t ext = alt->size() - ol;
    if ( alt->drxn_ == 2 ) alt->drxn_ = drxn;
    alt->offset( drxn ? ends_[1] - alt->ends_[0] - ol : ends_[0] - alt->ends_[1] + ol );
    seq_ = drxn ? seq_ + alt->getSeqEnd( ext, 1 ) : alt->getSeqEnd( ext, 0 ) + seq_;
    reads_.insert( alt->reads_.begin(), alt->reads_.end() );
    
    for ( int i = 1; i < merges[0].size(); i++ )
    {
        int32_t off = merges[0][i]->ends_[0] - ends_[0];
        merges[0][i]->seq_ = seq_;
        merges[0][i]->reads_ = reads_;
        for ( auto& read : merges[0][i]->reads_ ) read.second.offset( off );
    }
    
    for ( int i = 0; i < merges[0].size(); i++ )
    {
        int stop = merges[ bool( merges[1][i] ) ][i]->stop_[drxn];
        merges[0][i]->ends_[drxn] += drxn ? ext : -ext;
        merges[0][i]->bad_ = merges[0][i]->bad_ && merges[1][i]->bad_;
        merges[0][i]->mapped_ = merges[0][i]->mapped_ && merges[1][i]->mapped_;
        if ( merges[1][i] ) for ( Edge& e : merges[1][i]->edges( drxn ) ) merges[0][i]->addEdge( e, drxn, true );
        if ( merges[1][i] ) nodes.erase( merges[1][i] );
        merges[0][i]->stop_[drxn] = stop;
        merges[0][i]->setCoverage();
        merges[0][i]->remark();
        merges[0][i]->verified_ = false;
        merges[0][i]->ends_.reset( merges[0][i]->drxn_ );
        merges[0][i]->readTest();
    }
    
    return true;
}

bool Node::updateBad()
{
    Nodes tested[2], goods[2];
    for ( int d : { 0, 1 } ) if ( !updateBad( tested[d], d ) ) return bad_ = false;
    for ( int d : { 0, 1 } ) for ( Node* node : tested[d].nodes ) if ( !node->bad_ )
    {
        bool good = false;
        for ( Edge& e : node->edges_[!d] ) if ( !tested[d].find( e.node ) ) good = true;
        if ( good ) goods[d].fillIn( node, tested[d], d, true );
    }
    for ( int d : { 0, 1 } ) for ( Node* node : tested[d].nodes ) if ( !goods[d].find( node ) ) node->bad_ = true;
    return bad_ = true;
}

bool Node::updateBad( Nodes& tested, bool drxn )
{
    if ( drxn_ == 2 ) return false;
    if ( !tested.add( this ) ) for ( Edge& e : edges_[drxn] ) if ( !e.node->updateBad( tested, drxn ) ) return bad_ = false;
    return true;
}

void Node::updateStates( NodeRoll& nodes )
{
    Nodes drxn[2], tested[2];
    for ( Node* node : nodes.nodes ) if ( node->drxn_ == 2 ) for ( int d : { 0, 1 } ) drxn[d].fill( node, d, true );
    for ( Node* node : nodes.nodes ) if ( !node->bad_ && !drxn[0].find( node ) && !drxn[1].find( node ) )
    {
        node->bad_ = true;
        node->setUnverified();
    }
    for ( Node* node : nodes.nodes ) if ( node->drxn_ == 2 ) for ( int d : { 0, 1 } ) node->updateState( tested[d], drxn[!d], d );
    for ( int d : { 0, 1 } ) for ( Node* node : drxn[d].nodes ) assert( !node->bad_ );
}

void Node::updateState( Nodes& tested, Nodes& oppose, bool drxn )
{
    if ( !tested.add( this ) ) return;
    
    if ( !bad_ && drxn_ == !drxn && !oppose.find( this ) )
    {
        drxn_ = drxn;
        setUnverified();
    }
    
    if ( bad_ )
    {
        drxn_ = drxn;
        bad_ = false;
        setUnverified();
    }
    
    for ( Edge& e : edges_[drxn] ) e.node->updateState( tested, oppose, drxn );
}

void Node::pruneBad( NodeRoll& nodes )
{
    if ( nodes.size() < 1000 || params.rna ) return;
    
    updateStates( nodes );
    Nodes keep[3], tested[2];
    for ( Node* node : nodes.getGraph( 2 ).nodes ) node->pruneBad( keep, tested, 0, 0, 1 );
    for ( Node* node : keep[2].nodes ) for ( int d : { 0, 1 } )
    {
        node->pruneBad( keep, d );
        for ( Node* clone : node->clones( false ) ) clone->pruneBad( keep, d );
    }
    
    int initial = nodes.size(), cleared[2]{0}, culled = 0;
    for ( int d : { 0, 1 } ) for ( Node* node : keep[d].nodes ) if ( !keep[!d].find( node ) )
    {
        assert( !keep[2].find( node) );
        for ( Edge& e : node->edges_[!d] ) assert( !keep[2].empty() );
        assert( node->bad_ );
        node->clearEdges( !d, false );
        cleared[!d]++;
    }
    for ( int d : { 0, 1 } ) keep[2] += keep[d];
    
    for ( int i = 0; i < nodes.nodes.size(); i++ ) if ( !keep[2].find( nodes[i] ) )
    {
        assert( !keep[2].find( nodes[i] ) );
        nodes.erase( nodes[i], i );
        culled++;
    }
    
    cout << "Initial node count: " << initial << ", culled: " << culled << ". kept: " << keep[2].size() << ", cleared LEFT:" << cleared[0] << ", cleared RIGHT: " << cleared[1] << endl;
}

void Node::pruneBad( Nodes keep[3], Nodes tested[2], int32_t rev, int32_t dblrev, bool drxn )
{
    if ( !bad_ && ( rev || dblrev ) ) return;
    if ( !keep[2].add( this ) && bad_ ) return;
    
    if ( dblrev )
    {
        if ( !bad_ ) return;
        if ( dblrev < 300 ) for ( Edge& e : edges_[drxn] ) e.node->pruneBad( keep, tested, 0, e.node->size()-e.ol+dblrev, drxn );
    }
    else if ( rev )
    {
        if ( !bad_ ) return;
        for ( Edge& e : edges_[drxn] ) e.node->pruneBad( keep, tested, e.node->size()-e.ol+rev, 0, drxn );
        if ( rev < 200 ) for ( Edge& e : edges_[!drxn] ) e.node->pruneBad( keep, tested, rev, rev, !drxn );
    }
    else
    {
        if ( !tested[drxn].add( this ) ) return;
        for ( int d : { 0, 1 } ) if ( d == drxn || drxn_ > 1 ) for ( Edge& e : edges_[d] ) e.node->pruneBad( keep, tested, 0, 0, d );
        for ( Edge& e : edges_[!drxn] ) if ( e.node->bad_ ) e.node->pruneBad( keep, tested, max( 1, e.node->size()-e.ol ), 0, !drxn );
    }
}

void Node::pruneBad( Nodes keep[3], bool drxn )
{
    for ( Edge& e : edges_[drxn] ) if ( !keep[2].find( e.node ) && keep[!drxn].add( e.node ) ) e.node->pruneBad( keep, !drxn );
}

bool Node::prepFork( Querier &bwt, NodeRoll& nodes, int32_t dist, bool drxn )
{
    // Asure that this is a fork, and that the opposite direction is completed
    if ( edges_[drxn].size() < 2 || !isForkComplete( dist, 20, !drxn ) || isBlunt( 0, 20, !drxn ) ) return false;
    
    // Complete branches of this fork
    if ( !extendFork( bwt, nodes, dist, 8, drxn ) ) return false;
    
    // Asure that at least two branches are substantial;
    Nodes branches;
    for ( Edge& e : edges_[drxn] ) if ( !e.node->isBlunt( 0, 20, drxn ) ) branches += e.node;
    if ( branches.size() < 2 ) return false;
    
    // If any branches fork back, complete those reverse forks
    for ( Node* node : branches.nodes ) if ( !node->extendFork( bwt, nodes, dist, 6, !drxn ) && !node->isForkComplete( dist, 20, !drxn ) ) return false;
    
    // Complete all the pairings
    for ( Node* node : branches.nodes ) node->verifyFork( params.maxPeMean, false, !drxn );
    
    // Valid if on the branches has a fork in the opposite direction
    for ( Node* node : branches.nodes ) for ( Edge& e : node->edges_[!drxn] ) if ( !isClone( e.node ) && !e.node->isBlunt( 0, 20, drxn ) ) return true;
    
    // This fork must be accompanied by a fork in the opposite direction
    vector<Node*> fwd;
    for ( Edge& e : edges_[!drxn] ) if ( !e.node->isBlunt( 0, 20, drxn ) ) fwd.push_back( e.node );
    if ( fwd.size() > 1 ) return true;
    
    // Extend any other forks in this direction before the first fork in the opposite direction
    for ( Node* prv = this; fwd.size() == 1; )
    {
        if ( !fwd[0]->extendFork( bwt, nodes, dist, 8, drxn ) ) return false;
        
        for ( Edge& e : fwd[0]->edges_[drxn] ) if ( e.node != prv && !e.node->isBlunt( 0, 3, drxn ) )
        {
            if ( !e.node->extendFork( bwt, nodes, dist, 6, !drxn ) && !e.node->isForkComplete( dist, 20, !drxn ) ) return false;
        }
        prv = fwd[0];
        fwd.clear();
        if ( !prv->isForkComplete( dist, 20, !drxn ) || prv->isBlunt( 0, 20, !drxn ) ) return false;
        for ( Edge& e : prv->edges_[!drxn] ) if ( !e.node->isBlunt( 0, 3, !drxn ) ) fwd.push_back( e.node );
    }
    
    if ( fwd.size() < 2 ) return false;
    
    // Do all the necessary verifying
    fwd = vector<Node*>{ this };
    for ( Node* prv = this; fwd.size() == 1 && ( prv = fwd[0] ); )
    {
        fwd.clear();
        for ( Edge& e : prv->edges_[drxn] ) e.node->verifyFork( params.maxPeMean, false, !drxn );
        for ( Edge& e : prv->edges_[!drxn] ) if ( !e.node->isBlunt( 0, 3, !drxn ) ) fwd.push_back( e.node );
    }
    
    return true;
}

void Node::prune( Querier& bwt, NodeRoll& nodes )
{
    verify( nodes );
    pruneBlunt( bwt, nodes );
    merge( nodes );
    prunePaths( bwt, nodes );
    recoordinate( nodes );
    pruneIslands( bwt, nodes );
//    Bubble::prune( nodes );
    remap( bwt, nodes );
    reverify( nodes );
//    LocusFill::fill( bwt, nodes );
//    prunePaths( bwt, nodes );
    pruneBlunt( bwt, nodes );
    pruneIslands( bwt, nodes );
    nodes.test( true );
//    nodes.print( "/home/glen/locas_test.fasta", -3000, 0 );
}

void Node::pruneBlunt( Querier &bwt, NodeRoll& nodes )
{
    updateStates( nodes );
    int nodeCount[2]{ nodes.size(), 0 }, cloneCulled = 0;
    
//    for ( int d : { 0, 1 } ) for ( int again = 1; again-- > 0; )
//    {
//        for ( int i = 0; i < nodes.size(); i++ ) if ( nodes[i]->cloned_ && nodes[i]->edges_[d].empty() )
//        {
//            nodes.erase( nodes[i], i );
//            again = 1;
//        }
//    }
    
    nodeCount[1] = nodes.size();
    for ( bool d : { 0, 1 } )
    {
        Nodes strong;
        for ( Node* node : nodes.nodes ) if ( node->edges_[d].empty() ) node->setStrong( strong, 0, 10, !d );
        cout << strong.size() << endl;
        for ( int i = 0; i < nodes.size(); i++ ) if ( !strong.find( nodes[i] ) && nodes[i]->isBlunt( 0, 2, d ) )
        {
            bool verified = false, strength = true;
            for ( Edge& e : nodes[i]->edges_[!d] )
            {
                if ( e.node->verified_ ) verified = true;
                if ( !( strength = ( strong.find( e.node ) || e.node->isBlunt( 0, 2, !d ) ) ) ) break;
            }
            if ( verified && strength ) for ( Node* f : Nodes( nodes[i], d, true ).nodes ) nodes.erase( f, i );
        }
    }
    
    int loops[2]{0};
    for ( Node* node : nodes.nodes )
    {
        if ( !node->cloned_ ) continue;
        assert( !node->edges_[0].empty() || !node->edges_[1].empty() );
        for ( int d : { 0, 1 } ) if ( node->edges_[d].empty() ) loops[d]++;
    }
    nodes.test( true );
    cout << "Start with: " << nodeCount[0] << ", added: " << nodeCount[1] - nodeCount[0] << ", erased: " << nodeCount[1] - nodes.size() << " blunt nodes, culled: " << cloneCulled << " clones." << endl;
    cout << "Left loops: " << loops[0] << ", right loops: " << loops[1] << endl;
    
    NodeCounts counts[2];
    Nodes safe[2], bad;
    for ( Node* node : nodes.nodes ) for ( bool d : { 0, 1 } ) if ( node->edges_[d].empty() )
    {
        counts[d].countReads( node, !d );
        int minOl = max( params.readLen / 2 + 1, node->size() - node->getBestOverlap( !d, true ) + 10 );
        if ( !node->stop_[d] && bwt.isExtendable( node->seq_, minOl, d ) ) safe[d].fill( node, !d, true );
    }
    
    for ( bool d : { 0, 1 } ) for ( auto& nc : counts[d].map ) if ( nc.second < 10 && !safe[d].find( nc.first ) )
    {
        int32_t* p, cutoff = ( nc.second-1 ) * 10 + max( 0, nc.second-2 ) * 20 + 5;
        for ( Edge& e : nc.first->edges_[!d] ) if ( ( p = counts[d].get( e.node ) ) && *p > cutoff )
        {
            if ( e.node->getCoverage( d ) > cutoff ) bad += nc.first;
        }
    }
    for ( int i = 0; i < nodes.size(); i++ ) if ( bad.find( nodes[i] ) ) nodes.erase( nodes[i], i );
    pruneBad( nodes );
}

bool Node::pruneBranch( Querier& bwt, NodeRoll& nodes, int readsLeft, bool drxn )
{
    int32_t dist = max( params.readLen, params.maxPeMean - params.readLen );
    if ( verified_ ) return false;
    
    Nodes tried[2];
    if ( prepFork( bwt, nodes, dist, !drxn ) )
    {
        if ( ClaimNode::claim( this, nodes, tried, !drxn ) ) return true;
//        if ( ClaimNode::claim( this, nodes, !drxn ) ) return true;
    }
    
    readsLeft -= countReads( true );
    if ( readsLeft >= 0 ) for ( Edge& e : edges( drxn ) ) if ( e.node->pruneBranch( bwt, nodes, readsLeft, drxn ) ) return true;
    return false;
}

bool Node::pruneFwd( Querier& bwt, NodeRoll& nodes, Nodes* good, Nodes* tested, bool drxn )
{
    if ( tested && !tested->add( this ) ) return false;
    
    bool pruned = false;
    if ( edges_[drxn].size() > 1 )
    {
        int ol = params.readLen * 0.8;
        pair<string, ReadId> best = bwt.getCoverage( seq_, ol, drxn );
        float cover = min( coverage_, (float)best.second * params.readLen / ( params.readLen - ol ) );

        int testDist = min( params.readLen - 10 , 1 + ( 5000 / min( 1000, max( 100, (int)cover ) ) ) );
        vector< pair<Node*, float> > edges;
        
        bool kept = false;
        if ( !best.first.empty() && cover > 70 ) for ( Edge& e : edges_[drxn] )
        {
            int diff = 0;
            int len = min( (int)best.first.size(), e.node->size() - e.ol ), dist = e.node->size() - e.ol;
            assert( e.leap || e.ol > 0 );
            for ( int i = 0; !diff && i < len; i++ ) if ( drxn ? e.node->seq_[i+e.ol] != best.first[i] : e.node->seq_.end()[-i-1-e.ol] != best.first.end()[-i-1] ) diff = i+1;
            if ( !diff ) kept = true;
            
            if ( !diff ) edges.push_back( make_pair( e.node, cover ) );
            else if ( dist > testDist ) edges.push_back( make_pair( e.node, float( e.node->countReads( true )-1 ) * params.readLen / (float)dist ) );
            else
            {
                int lens[2]{ min( size(), params.readLen-diff ), diff }, cut = drxn ? e.node->ends_[0] + e.ol + diff : e.node->ends_[1] - e.ol - diff;
                string seq = drxn ? getSeqEnd( lens[0], 1 ) + e.node->seq_.substr( e.ol, lens[1] ) 
                                  : e.node->seq_.substr( e.node->size()-lens[1], lens[1] ) + getSeqEnd( lens[0], 0 );
                pair<string, ReadId> alt = bwt.getCoverage( seq, ol, drxn );
                for ( auto& read : e.node->reads_ ) if ( !read.second.redundant && ( drxn ? read.second[1] <= cut : cut <= read.second[0] ) ) alt.second++;
                edges.push_back( make_pair( e.node, float( alt.second - 1 ) * params.readLen / float( dist + params.readLen - ol ) ) );
            }
        }
        if ( !kept ) return false;
        
        for ( pair<Node*, float> edge : edges ) if ( edge.second < cover / 20 && ( !good || !good->find( edge.first ) ) )
        {
            if ( drxn_ < 2 && drxn != drxn_ && !edge.first->bad_ );
            else if ( edge.first->edges_[!drxn].size() == 1 ) nodes.erase( edge.first );
            else ( removeEdge( edge.first, drxn, true ) );
        }
    }
    
    if ( tested ) for ( Edge& e : edges_[drxn] ) if ( e.node->pruneFwd( bwt, nodes, good, tested, drxn) ) pruned = true;
    return pruned;
}

void Node::pruneIslands( Querier &bwt, NodeRoll& nodes )
{
    Nodes tested;
    for ( Node* node : nodes.nodes ) if ( node->drxn_ == 2 ) tested.fill( node );
    
    int culled = 0;
    for ( int i = 0; i < nodes.size(); i++ ) if ( nodes[i]->bad_ && !tested.find( nodes[i] ) )
    {
        Nodes island;
        island.fill( nodes[i] );
        bool good = false;
        for ( Node* node : island.nodes ) if ( node->edges_[0].empty() && node->isSubstantial( 0, 20, 1 ) ) good = true;
        if ( !good ) culled += island.size();
        if ( !good ) for ( Node* node : island.nodes ) nodes.erase( node, i );
        else tested += island;
    }
    
    cout << "Pruned " << culled << " island nodes, " << nodes.size() << " total nodes remain." << endl;
}

void Node::pruneLoops( NodeRoll& nodes )
{
    Nodes tested[2];
    for ( Node* node : nodes.nodes ) for ( int d : { 0, 1 } ) if ( tested[d].add( node ) )
    {
        vector<Node*> path;
        while ( node && node->edges_[d].size() == 1 )
        {
            node = node->edges_[d][0].node;
            tested[d] += node;
            for ( int i = 0; i < path.size(); i++ ) if ( path[i] == node )
            {
                Node* edge = path[i];
                node = path.back();
                int minOl = node->getOverlap( edge, d );
                for ( int j = i; j+1 < path.size(); j++ )
                {
                    int ol = path[j]->getOverlap( path[j+1], d );
                    if ( ol >= minOl ) break;
                    node = path[j];
                    edge = path[j+1];
                    minOl = ol;
                }
                assert( false );
                node = NULL;
                break;
            }
            if ( node ) path.push_back( node );
        }
    }
}

void Node::prunePaths( Querier& bwt, NodeRoll& nodes )
{
    int counted = 0;
    pruneBad( nodes );
    for ( int again = 1; counted++ < 3 && again-- > 0; )
    {
        pruneBlunt( bwt, nodes );
        pruneLoops( nodes );
        Nodes tested[2], tried[2];
        for ( int d : { 0, 1 } ) for ( Node* node : nodes.getGraph( 2 ).nodes ) if ( nodes.find( node ) && node->prunePaths( bwt, nodes, tested[d], tried, 0, d ) ) again = 1;
        reverify( nodes );
        if ( !again && Bubble::prune( nodes ) ) again = 1;
    }
    Node::merge( nodes );
}

bool Node::prunePaths( Querier& bwt, NodeRoll& nodes, Nodes& tested, Nodes tried[2], int branched, bool drxn )
{
    if ( bad_ || !tested.add( this ) ) return false;
    
    bool claimed = false;
    
    int32_t dist = max( params.readLen, params.maxPeMean - params.readLen );
    if ( prepFork( bwt, nodes, dist, !drxn ) )
    {
        nodes.test( true );
        claimed = ClaimNode::claim( this, nodes, tried, !drxn );
        nodes.test( true );
    }
    
    if ( branched > 2 )
    {
        return false;
    }
    
    int branches[2]{0};
    for ( int d : { 0, 1 } ) for ( Edge& e : edges( d ) ) if ( e.node->verified_ ) branches[d]++;
    branched = max( branched + branches[drxn] - max( 1, branches[!drxn] ), 0 );
    
    if ( !bad_ && verified_ ) for ( Edge& e : edges( drxn ) ) if ( e.node->prunePaths( bwt, nodes, tested, tried, branched, drxn ) ) claimed = true;
    
    return claimed;
}

void Node::setStrong( Nodes& strong, int readCount, int readLimit, bool drxn )
{
    readCount += countReads( true );
    if ( readCount >= readLimit ) strong.fill( this, drxn, true );
    else for ( Edge& e : edges_[drxn] ) if ( !strong.find( e.node ) ) e.node->setStrong( strong, readCount, readLimit, drxn );
}

