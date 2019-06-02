/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "node.h"
#include "node_claim.h"
#include <algorithm>

void Node::absorbClone( Node* clone )
{
    for ( int d : { 0, 1 } )
    {
        for ( Edge& e : clone->edges_[d] ) if ( !isEdge( e.node, d, true ) ) addEdge( e.node, e.ol, d, false, e.isLeap );
        clone->clearEdges( d );
        if ( !bad_ && ( drxn_ == 2 || drxn_ == d ) ) for ( Edge& e : edges_[d] ) e.node->setNotBad( d );
    }
    
}

bool Node::absorbClone( Node* clone, bool drxn )
{
    assert( edges_[drxn].empty() || clone->edges_[drxn].empty() );
    Nodes fwd[2]{ Nodes( this, drxn, false, false ), Nodes( clone, drxn, false, false ) };
    if ( !clone->bad_ || fwd[0].find( clone ) || fwd[1].find( this ) ) return false;
    Nodes bck[2]{ Nodes( this, !drxn, false, false ), Nodes( clone, !drxn, false, false ) };
    for ( int i = 0; i < 2 && !bad_ && !clone->bad_; i++ )
    {
        for ( Node* f : fwd[i].nodes )
        {
            if ( !f->cloned_ || !f->edges_[drxn].empty() ) continue;
            for ( Node* node : f->cloned_->nodes )
            {
                if ( !bck[i].find( node ) && bck[!i].find( node ) ) assert( false );
            }
        }
    }
    absorbClone( clone );
    return true;
}

bool Node::claimBad( Node* node, bool drxn )
{
    bool removed = false, test = false;
    for ( int i = 0; i < node->edges_[!drxn].size(); i++ )
    {
        if ( node->edges_[!drxn][i].node == this ) test = true;
        if ( node->edges_[!drxn][i].node->bad_ ) continue;
        node->edges_[!drxn][i].node->removeEdge( node, drxn );
        node->edges_[!drxn].erase( node->edges_[!drxn].begin() + i-- );
        removed = true;
    }
    assert( removed && test );
    if ( removed ) node->setBad( drxn );
    return removed;
}

void Node::claimBads( NodeRoll &nodes )
{
    for ( Node* node : nodes.nodes )
    {
        if ( node->bad_ ) continue;
        for ( int d : { 0, 1 } )
        {
            if ( node->drxn_ < 2 && d != node->drxn_ ) continue;
            for ( Edge &e : node->edges_[d] ) if ( e.node->bad_ )e.node->setNotBad( d );
        }
    }
}

bool Node::claimGood( bool drxn )
{
    if ( bad_ ) return false;
    if ( drxn_ < 2 && drxn_ != drxn ) return false;
    bool claimed = false;
    for ( Edge &e : edges_[drxn] )
    {
        if ( e.node->drxn_ >= 2 || e.node->drxn_ == drxn || !e.node->bad_ ) continue;
        int32_t offset = drxn ? ends_[1] - e.ol - e.node->ends_[0] : ends_[0] - e.node->ends_[1] + e.ol;
        e.node->offset( offset );
        e.node->drxn_ = drxn;
        e.node->bad_ = false;
        e.node->claimGood( drxn );
        claimed = true;
    }
    return claimed;
}

void Node::discardBlunt( NodeRoll& nodes, bool drxn )
{
    for ( Edge& e : edges_[drxn] ) e.node->discardBlunt( nodes, drxn );
    nodes.erase( this );
}

Edge Node::getFoldTarget( bool drxn )
{
    Edge f( NULL, 0, false );
    for ( Edge& e : edges_[drxn] )
    {
        if ( e.node->isBlunt( 0, 2, drxn ) ) continue;
        if ( f.node || !e.node->verified_ ) return Edge( NULL, 0, false );
        f = e;
    }
    return f;
}

void Node::getFoldable( NodeLine& fold, vector<NodeLine>& folds, bool drxn )
{
    for ( int i = 1; i < edges_[drxn].size(); i++ )
    {
        NodeLine alt = fold;
        alt.add( edges_[drxn][i], drxn );
        edges_[drxn][0].node->getFoldable( alt, folds, drxn );
        folds.push_back( alt );
    }
    if ( !edges_[drxn].empty() )
    {
        fold.add( edges_[drxn][0], drxn );
        edges_[drxn][0].node->getFoldable( fold, folds, drxn );
    }
}

bool Node::isBadRelevant()
{
    if ( !bad_ ) return false;
    assert( drxn_ < 2 );
    for ( Edge &e : edges_[drxn_] )
    {
        if ( e.node->bad_ ) continue;
        if ( e.node->verified_ ) return true;
        for ( Edge &re : e.node->edges_[!drxn_] ) if ( re.node->verified_ ) return true;
    }
    return false;
}

bool Node::isBlunt( int readCount, int readLimit, bool drxn )
{
    readCount += countReads( true );
    if ( readCount > readLimit || verified_ || cloned_ ) return false;
    for ( Edge& e : edges_[drxn] ) if ( !e.node->isBlunt( readCount, readLimit, drxn ) ) return false;
    return !edges_[drxn].empty() || stop_[drxn] == DEAD_END || stop_[drxn] == SHORT_END || stop_[drxn] == BLUNT_END;
}

//bool Node::isExtendLen( int32_t dist, bool drxn )
//{
//    dist -= size();
//    if ( dist <= 0 ) return true;
//    NodeOffsets offs( this, dist, drxn );
//    for ( auto &no : offs.map )
//    {
//        if ( !no.first->edges_[drxn].empty() || dist < no.second[0] ) continue;
//        if ( no.first->stop_[drxn] == DEAD_END || no.first->stop_[drxn] == SHORT_END ) continue;
//        if ( no.first->isContinue( drxn ) || no.first->unpauseable( drxn ) ) return false;
//    }
//    return true;
//}

bool Node::isBranchComplete( int readCount, int readMin, bool drxn )
{
    readCount += countReads( true );
    if ( readCount >= readMin ) return true;
    for ( Edge& e : getAltEdges( drxn ) ) if ( !e.node->isBranchComplete( readCount, readMin, drxn ) ) return false;
    return !isContinue( drxn );
}

bool Node::isFoldable( Node* from, bool drxn )
{
    for ( Edge& b : edges_[!drxn] )
    {
        if ( b.node == from ) continue;
        for ( Edge& f : b.node->edges_[drxn] )
        {
            if ( f.node != this && !f.node->isBlunt( 0, 2, drxn ) ) return false;
        }
    }
    return true;
}

bool Node::isForkComplete( int32_t dist, int readMin, bool drxn )
{
    bool good = true;
    if ( dist <= 0 ) return true;
    for ( Node* f : Nodes( this, dist, drxn, true, true ).nodes ) if ( f->isContinue( drxn ) ) good = false;
    if ( !good && readMin ) for ( Edge& e : getAltEdges( drxn ) ) if ( !e.node->isBranchComplete( 0, readMin, drxn ) ) return false;
    return good || readMin;
}

//bool Node::isForkPrepped( int32_t dist, bool drxn )
//{
//    if ( !isForkComplete( dist, drxn ) ) return false;
//    for ( Edge& e : edges_[drxn] ) if ( !e.node->isForkComplete( dist, !drxn ) ) return false;
//    return true;
//}

bool Node::isForkSplayed( int32_t dist, int maxBranches, bool drxn )
{
    Nodes fwd( this, dist, drxn, false, false );
    int branches = 0;
    for ( Node* f : fwd.nodes )
    {
        bool ended = !f->stop_[drxn];
        for ( Edge& e : f->edges_[drxn] ) if ( fwd.find( e.node ) ) ended = false;
        if ( ended ) branches++;
    }
    return branches <= maxBranches;
}

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

bool Node::merge( NodeRoll &nodes )
{
    int merged = 0;
    Nodes tested;
    for ( Node* node : nodes.getGraph( 2 ).nodes )
    {
        if ( !nodes.find( node ) || !tested.add( node ) ) continue;
        for ( int d : { 0, 1 } ) while ( node->merge( nodes, d ) ) merged++;
        node->merge( nodes, tested, merged );
    }
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
    if ( edges_[drxn][0].isLeap || ( edges_[drxn][0].isLeap = edges_[drxn][0].ol <= 0 ) ) return false;
    
    Node* alt = edges_[drxn][0].node;
    int ol = edges_[drxn][0].ol;
    vector<Node*> merges[2]{ { this }, { alt } };
    if ( !isMergeable( alt, merges, drxn ) ) return false;
    for ( int i = 0; i < merges[0].size(); i++ ) for ( int d = 0; d < 2 && merges[d][i]; d++ ) merges[d][i]->clearPaired( false );
    
    int32_t ext = alt->size() - ol;
    alt->offset( drxn ? ends_[1] - alt->ends_[0] - ol : ends_[0] - alt->ends_[1] + ol );
    seq_ = drxn ? seq_ + alt->getSeqEnd( ext, 1 ) : alt->getSeqEnd( ext, 0 ) + seq_;
    reads_.insert( alt->reads_.begin(), alt->reads_.end() );
    
    for ( int i = 1; i < merges[0].size(); i++ )
    {
        int32_t off = ends_[0] - merges[0][i]->ends_[0];
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

bool Node::prepBranch( int32_t dist, int maxHits, bool drxn )
{
    if ( !isForkComplete( dist - size() + params.readLen, 20, drxn ) || !isForkComplete( dist, 0, !drxn ) || isBlunt( 0, 20, drxn )  ) return false;
    verifyFork( params.maxPeMean, false, !drxn );
    if ( maxHits < 0 ) return true;
    int hits = 0;
    Nodes fwd( this, dist, drxn, true, true ), bck( this, !drxn, false, true );
    for ( Node* f : fwd.nodes ) for ( auto &np : f->hits_.pairs[!drxn] ) if ( bck.find( np.first ) ) hits += np.second.count;
    if ( hits <= maxHits )
    {
        int x = 0;
    }
    return hits <= maxHits && !bck.find( this );
}

bool Node::prepFork( Querier &bwt, NodeRoll& nodes, int32_t dist, bool drxn )
{
    // Asure that this is a fork, and that the opposite direction is completed
    if ( getAltEdges( drxn ).size() < 2 || !isForkComplete( dist, 20, !drxn ) || isBlunt( 0, 20, !drxn ) ) return false;
    
    // Complete branches of this fork
    if ( !extendFork( bwt, nodes, dist, 8, drxn ) ) return false;
    
    // Asure that at least two branches are substantial;
    Nodes branches;
    for ( Edge& e : getAltEdges( drxn ) ) if ( !e.node->isBlunt( 0, 20, drxn ) ) branches += e.node;
    if ( branches.size() < 2 ) return false;
    
    // If any branches fork back, complete those reverse forks
    for ( Node* node : branches.nodes ) if ( !node->extendFork( bwt, nodes, dist, 6, !drxn ) && !node->isForkComplete( dist, 20, !drxn ) ) return false;
    
    // Complete all the pairings
    for ( Node* node : branches.nodes ) node->verifyFork( params.maxPeMean, false, !drxn );
    
    
    for ( Node* node : branches.nodes ) for ( Edge& e : node->edges_[!drxn] ) if ( !isClone( e.node ) && !e.node->isBlunt( 0, 20, drxn ) ) return true;
    
    vector<Node*> fwd;
    for ( Edge& e : getAltEdges( !drxn ) ) if ( !e.node->isBlunt( 0, 20, drxn ) ) fwd.push_back( e.node );
    if ( fwd.size() > 1 ) return true;
    
    for ( Node* prv = this; fwd.size() == 1; )
    {
        if ( !fwd[0]->extendFork( bwt, nodes, dist, 8, drxn ) ) return false;
        
        for ( Edge& e : fwd[0]->getAltEdges( drxn ) ) if ( e.node != prv && !e.node->isBlunt( 0, 3, drxn ) )
        {
            if ( !e.node->extendFork( bwt, nodes, dist, 6, !drxn ) && !e.node->isForkComplete( dist, 20, !drxn ) ) return false;
        }
        prv = fwd[0];
        fwd.clear();
        if ( !prv->isForkComplete( dist, 20, !drxn ) || prv->isBlunt( 0, 20, !drxn ) ) return false;
        for ( Edge& e : prv->getAltEdges( !drxn ) ) if ( !e.node->isBlunt( 0, 3, !drxn ) ) fwd.push_back( e.node );
    }
    
    if ( fwd.size() < 2 ) return false;
    
    fwd = vector<Node*>{ this };
    for ( Node* prv = this; fwd.size() == 1 && ( prv = fwd[0] ); )
    {
        fwd.clear();
        for ( Edge& e : prv->getAltEdges( drxn ) ) e.node->verifyFork( params.maxPeMean, false, !drxn );
        for ( Edge& e : prv->getAltEdges( !drxn ) ) if ( !e.node->isBlunt( 0, 3, !drxn ) ) fwd.push_back( e.node );
    }
    
    return true;
}

void Node::prune( Querier& bwt, NodeRoll& nodes )
{
    Node::verify( nodes );
    Node::pruneBlunt( bwt, nodes );
    Node::merge( nodes );
    Node::prunePaths( bwt, nodes );
    recoordinate( nodes );
    remap( bwt, nodes );
    reverify( nodes );
    Node::prunePaths( bwt, nodes );
    Node::pruneBlunt( bwt, nodes );
    nodes.test( true );
//    nodes.print( "/home/glen/locas_test.fasta", -3000, 0 );
}

//void Node::pruneBad( Querier &bwt, NodeRoll &nodes )
//{
//    for ( int d : { 0, 1 } )
//    {
//        Nodes tested;
//        for ( Node* node : nodes.getGraph( 2 ).nodes ) if ( nodes.find( node ) ) node->pruneBad( bwt, nodes, tested, d );
//    }
//}

//void Node::pruneBad( Querier &bwt, NodeRoll& nodes, Nodes& tested, bool drxn )
//{
//    for ( Edge& e : edges_[!drxn] ) if ( drxn_ < 2 && !e.node->bad_ && !tested.find( e.node ) ) return;
//    if ( !verified_ || !tested.add( this ) ) return;
//    assert( !bad_ );
//    bool bad = false;
//    for ( Edge& e : edges_[!drxn] ) if ( e.node->bad_ && !e.node->isBlunt( 0, 5, !drxn ) ) bad = true;
//    if ( bad ) ClaimBad( this, nodes, drxn );
//    for ( Edge& e : edges_[drxn] ) e.node->pruneBad( bwt, nodes, tested, drxn );
//}

void Node::pruneBlunt( Querier &bwt, NodeRoll& nodes )
{
    int nodeCount[2]{ nodes.size(), 0 }, cloneCulled = 0;
    
    for ( int i = 0; i < nodes.size(); i++ )
    {
        if ( !nodes[i]->cloned_ ) continue;
        if ( !nodes[i]->pruneClone( 0 ) && !nodes[i]->pruneClone( 1 ) ) continue;
        vector<Node*> clones = nodes[i]->cloned_->nodes;
        nodes.erase( nodes[i], i );
        for ( Node* clone : clones ) clone->setState();
        clones[0]->reverify();
        cloneCulled++;
    }
    
    for ( int i = 0; i < nodes.size(); i++ )
    {
        if ( !nodes[i]->verified_ ) continue;
        for ( int d : { 0, 1 } ) for ( Node* node : NodeRoll::next( nodes[i], d ).nodes ) if ( node->isContinue( d ) ) node->extendNode( bwt, nodes, d );
    }
    nodeCount[1] = nodes.size();
    for ( bool d : { 0, 1 } )
    {
        Nodes strong;
        for ( Node* node : nodes.nodes ) if ( node->edges_[d].empty() ) node->setStrong( strong, 0, 10, !d );
        cout << strong.size() << endl;
        for ( int i = 0; i < nodes.size(); i++ )
        {
            if ( strong.find( nodes[i] ) || !nodes[i]->isBlunt( 0, 2, d ) ) continue;
            bool verified = false, strength = true;
            for ( Edge& e : nodes[i]->edges_[!d] )
            {
                if ( e.node->verified_ ) verified = true;
                if ( !( strength = strong.find( e.node ) || e.node->isBlunt( 0, 2, !d ) ) ) break;
            }
            if ( verified && strength ) for ( Node* f : Nodes( nodes[i], d, true, false ).nodes ) nodes.erase( f, i );
        }
    }
    
    int loops[2]{0};
    for ( Node* node : nodes.nodes )
    {
        if ( !node->cloned_ ) continue;
        assert( !node->edges_[0].empty() || !node->edges_[1].empty() );
        for ( int d : { 0, 1 } ) if ( node->edges_[d].empty() )
        {
            loops[d]++;;
        }
    }
    cout << "Start with: " << nodeCount[0] << ", added: " << nodeCount[1] - nodeCount[0] << ", erased: " << nodeCount[1] - nodes.size() << " blunt nodes, culled: " << cloneCulled << " clones." << endl;
    cout << "Left loops: " << loops[0] << ", right loops: " << loops[1] << endl;
}

void Node::pruneBlunt( NodeRoll &nodes )
{
    for ( int d : { 0, 1 } )
    {
        Nodes tested;
        for ( Node* node : nodes.getGraph( 2 ).nodes ) if ( nodes.find( node ) ) node->pruneBlunt( nodes, tested, d );
    }
}

void Node::pruneBlunt( NodeRoll& nodes, Nodes& tested, bool drxn )
{
    if ( !verified_ || !tested.add( this ) ) return;
    pruneBlunt( nodes, 0 );
    pruneBlunt( nodes, 1 );
    for ( Edge& e : edges_[drxn] ) e.node->pruneBlunt( nodes, tested, drxn );
}

void Node::pruneBlunt( NodeRoll& nodes, bool drxn )
{
    vector<Edge> edges[3];
    
    // Ensure that there are any substantial paths from this node
    for ( Edge& e : edges_[drxn] ) if ( e.node->verified_ ) edges[0].push_back( e );
    if ( edges[0].empty() ) return;
    
    // Catalog unsubstantial edges as foldable or not
    for ( Edge& e : edges_[drxn] )
    {
        if ( e.node->verified_ || !e.node->isBlunt( 0, 2, drxn ) ) continue;
        edges[ 1 + ( ( edges[0].size() != 1 ) || !e.node->isFoldable( this, drxn ) ) ].push_back( e );
    }
    
    // Discard unfoldable edges
    for ( Edge& e : edges[2] ) e.node->discardBlunt( nodes, drxn );
    if ( edges[1].empty() ) return;
    assert( edges[0].size() == 1 );
    
    // Attempt to fold any remaining unsubstantial edges
    NodeLine tar( this, drxn );
    Nodes mapped( this );
    for ( Edge e = edges[0][0]; tar.add( e, drxn ); ) e = e.node->getFoldTarget( drxn );
    for ( Edge& e : edges[1] ) tar.map( e, mapped, 0, 0, drxn );
    tar.fold( nodes, mapped, drxn );
    nodes.test();
    
//    Edge f = getFoldTarget( drxn );
//    if ( !f.node ) return;
//    vector<Edge> foldable;
//    for ( int i = 0; i < edges_[drxn].size(); i++ )
//    {
//        if ( !edges_[drxn][i].node->isBlunt( 0, 2, drxn ) ) continue;
//        if ( edges_[drxn][i].node->isFoldable( this, drxn ) ) foldable.push_back( edges_[drxn][i] );
//        else edges_[drxn][i--].node->discardBlunt( nodes, drxn );
//    }
//    if ( foldable.empty() ) return;
//    NodeLine tar( this, drxn );
//    Nodes mapped( this );
//    while ( tar.add( f, drxn ) ) f = f.node->getFoldTarget( drxn );
//    for ( Edge& e : foldable ) tar.map( e, mapped, 0, 0, drxn );
//    tar.fold( nodes, mapped, drxn );
//    nodes.test();
}

void Node::pruneBubble( Querier& bwt, NodeRoll& nodes )
{
    for ( int d : { 0, 1 } )
    {
        Nodes tested;
        for ( Node* node : nodes.getGraph( 2 ).nodes ) node->pruneBubble( nodes, tested, d );
    }
}

void Node::pruneBubble( NodeRoll& nodes, Nodes& tested, bool drxn )
{
    if ( !tested.add( this ) || !verified_ ) return;
    if ( drxn_ < 2 && edges_[!drxn].size() > 1 )
    {
        ClaimFork( this, nodes, !drxn );
    }
    for ( Edge& e : edges_[drxn] ) e.node->pruneBubble( nodes, tested, drxn );
}

bool Node::pruneClone( bool drxn )
{
    if ( !cloned_ || !edges_[drxn].empty() ) return false;
    if ( edges_[!drxn].empty() ) return true;
    int cloned = 0;
    for ( Node* clone : cloned_->nodes ) if ( !clone->edges_[0].empty() || !clone->edges_[1].empty() ) cloned++;
    if ( !cloned ) return false;
    for ( Edge& e : edges_[!drxn] )
    {
        for ( Edge& re : e.node->edges_[drxn] ) assert( !cloned_->find( re.node ) );
        for ( Edge& re : e.node->edges_[drxn] ) if ( cloned_->find( re.node ) ) return true;
    }
    for ( Node* clone : cloned_->nodes ) if ( !clone->edges_[0].empty() && !clone->edges_[1].empty() ) return false;
    for ( Node* clone : cloned_->nodes ) if ( Nodes( clone, drxn, false, false ).find( this ) ) return false;
    for ( Node* clone : cloned_->nodes )
    {
        if ( clone->edges_[0].empty() && clone->edges_[1].empty() ) continue;
        for ( Edge& e : edges_[!drxn] ) clone->addEdge( e, !drxn, true );
        clone->setCoverage();
    }
    for ( Edge& e : edges( !drxn ) ) e.node->removeEdge( this, drxn, true );
    return true;
}

void Node::pruneClones( Querier& bwt, NodeRoll& nodes )
{
    for ( int again = 1; again-- > 0; )
    {
        for ( int i = 0; i < nodes.size(); i++ )
        {
            if ( nodes[i]->cloned_ && nodes[i]->edges_[0].empty() && nodes[i]->pruneClones( nodes, i, 0 ) ) again = 1;
            if ( nodes[i]->cloned_ && nodes[i]->edges_[1].empty() && nodes[i]->pruneClones( nodes, i, 1 ) ) again = 1;
        }
    }
}

bool Node::pruneClones( NodeRoll& nodes, int& i, bool drxn )
{
    assert( !cloned_->empty() );
    for ( Node* clone : cloned_->nodes )
    {
        if ( absorbClone( clone, drxn ) ) nodes.erase( clone, i );
        else if ( clone->absorbClone( this, drxn ) ) nodes.erase( this, i );
        else continue;
        return true;
    }
    if ( bad_ ) return false;
    
    for ( Node* clone : cloned_->nodes ) if ( !clone->bad_ ) return false;
    Nodes bads, forks;
    for ( Node* clone : cloned_->nodes ) bads.fillBad( clone, true, drxn );
    for ( Node* bad : bads.nodes ) for ( Edge& e : bad->edges_[drxn] ) if ( !e.node->bad_ ) forks += e.node;
    assert( !forks.empty() );
    for ( Node* fork : forks.nodes )
    {
        Node* clone = new Node( fork, nodes, drxn_, true );
        for ( auto it = fork->edges_[!drxn].begin(); it != fork->edges_[!drxn].end(); it++ )
        {
            if ( !bads.find( it->node ) ) continue;
            clone->addEdge( it->node, it->ol, !drxn, false, it->isLeap );
            it->node->removeEdge( fork, drxn );
            it = fork->edges_[!drxn].erase( it )-1;
        }
    }
    while ( cloned_ )
    {
        assert( !cloned_->empty() );
        absorbClone( cloned_->nodes.back() );
        nodes.erase( cloned_->nodes.back(), i );
    }
    
    return true;
}

//void Node::pruneEdges( Querier &bwt, NodeRoll &nodes )
//{
//    NodeRoll seed = nodes.getGraph( 2 );
//    Nodes tested( seed.nodes );
//    for ( Node* node : seed.nodes )
//    {
//        for ( Edge& e : node->edges_[0] ) e.node->pruneEdges( bwt, nodes, tested, 0 );
//        for ( Edge& e : node->edges_[1] ) e.node->pruneEdges( bwt, nodes, tested, 1 );
//    }
//}

//void Node::pruneEdges( Querier& bwt, NodeRoll& nodes, Nodes& tested, bool drxn )
//{
//    for ( Edge& e : edges_[!drxn] ) if ( drxn_ < 2 && !e.node->bad_ && !tested.find( e.node ) ) return;
//    if ( tested.find( this ) ) return;
//    tested.add( this );
//    
//    sortEdges( 0 );
//    sortEdges( 1 );
//    for ( int i = edges_[drxn].size(); verified_ && --i >= 0 && edges_[drxn].size() > 1; )
//    {
//        pruneEdges( bwt, nodes, edges_[drxn][i], drxn );
//    }
//    for ( int i = edges_[!drxn].size(); --i >= 0 && edges_[!drxn].size() > 1; )
//    {
//        if ( edges_[!drxn][i].node->bad_ ) pruneEdges( bwt, nodes, edges_[!drxn][i], !drxn );
//    }
//    
//    if ( verified_ ) for ( int i = 0; i < edges_[drxn].size(); i++ ) edges_[drxn][i].node->pruneEdges( bwt, nodes, tested, drxn );
//}

//void Node::pruneEdges( Querier& bwt, NodeRoll& nodes, Edge& e, bool drxn )
//{
//    if ( e.node->edges_[!drxn].size() < 2 ) return;
//    
//    extendFork( bwt, nodes, params.maxPeMean-params.readLen, 10, drxn );
//    e.node->extendFork( bwt, nodes, params.maxPeMean-params.readLen, 10, !drxn );
//    verifyFork( params.maxPeMean, false, drxn );
//    e.node->verifyFork( params.maxPeMean, false, !drxn );
//    Nodes fwds[2], bcks[2]; // fwds[0] = branch, fwds[1] = rest, bcks[0] = this, bcks[1] = rest
//    NodeOffsets offs[2];
//    NodeOffsets::forkSets( this, e, offs, fwds, bcks, params.maxPeMax, drxn );
//    
//    int pen = getBestOverlap( drxn ) + e.node->getBestOverlap( !drxn ) - ( e.ol * 2 );
//    int hits[2][2]{ {0,0}, {0,0} };
//    int cutoff = 1 + ( params.branchMinHits * 30 / ( ( 30 + pen ) * 4 ) );
//    
//    for ( int i : {0,1} ) for ( int j : {0,1} ) for ( Node* b : bcks[i].nodes ) hits[i][j] += b->hits_.get( fwds[j], offs[0], offs[1], drxn );
//    
//    if ( hits[1][0] > cutoff && !hits[0][0] && hits[0][1] > cutoff && !hits[1][1] )
//    {
//        Node* cut = e.node;
//        cut->removeEdge( this, !drxn );
//        removeEdge( cut, drxn );
//        cut->setBad( drxn );
//        setBad( !drxn );
//    }
//}

bool Node::pruneFork( Querier& bwt, NodeRoll& nodes, bool drxn )
{
    int32_t dist = max( params.readLen, params.maxPeMean - params.readLen );
    if ( !prepFork( bwt, nodes, dist, drxn ) ) return false;
    ClaimFork cf( this, nodes, !drxn );
    return cf.claimed_;
}

void Node::pruneLoops( Querier& bwt, NodeRoll& nodes )
{
    for ( int again = 1; again-- > 0; )
    {
        for ( Node* node : nodes.nodes )
        {
            if ( node->cloned_ ) for ( int d = 0; !again && d < 2; d++ ) again = node->pruneLoops( bwt, nodes, d );
            if ( again ) break;
        }
    }
}

bool Node::pruneLoops( Querier& bwt, NodeRoll& nodes, bool drxn )
{
    if ( !edges_[drxn].empty() || bad_ ) return false;
    
    for ( Node* clone : cloned_->nodes ) if ( pruneLoops( clone, nodes, drxn ) ) return true;
    
    verifyFork( params.maxPeMax, true, drxn );
//    NodeOffsets offs[2];
//    NodeOffsets::forkSets( this, offs, params.maxPeMax, drxn );
//    int score = 0;
//    Nodes fwd;
//    for ( auto& f : offs[1].map ) fwd += f.first;
//    for ( auto& b : offs[0].map ) score += b.first->hits_.get( fwd, offs[0], offs[1], drxn );
    
    assert( false );
    ClaimFork( this, nodes, drxn );
    return false;
}

bool Node::pruneLoops( Node* clone, NodeRoll& nodes, bool drxn )
{
    Nodes bck[2]{ Nodes::isBad( this, false, !drxn ), Nodes::isBad( clone, false, !drxn ) };
    Nodes fwd[3]{ Nodes::inSet( clone, bck[0], drxn, true ), Nodes(), Nodes() };
    NodeRoll forks[2];
    for ( Node* b : Nodes::notSet( this, bck[1], !drxn, true ).nodes )
    {
        for ( Edge& e : b->edges_[!drxn] )
        {
            if ( !bck[0].find( e.node ) || fwd[0].find( e.node ) ) continue;
            forks[0].add( e.node );
            fwd[1].fillIn( b, bck[0], drxn, true );
            fwd[2].fill( e.node, drxn, false, false );
        }
    }
    forks[0] -= fwd[2];
    if ( forks[0].empty() ) return false;
    fwd[0] -= fwd[1];
    for ( Node* f : fwd[0].nodes ) for ( Edge& e : f->edges_[drxn] ) if ( fwd[1].find( e.node ) ) forks[1].add( e.node );
    NodeOffsets offs( clone, drxn, true );
    NodeOffset* off = offs.get( this );
    int best[2]{ 200, 200 };
    for ( int i = 0; i < forks[0].size(); i++ )
    {
        NodeOffsets forkOffs( forks[0][i], drxn, true );
        NodeOffset* forkOff[2]{ forkOffs.get( this ), forkOffs.get( clone ) };
        int diff = abs( (*forkOff[0])[!drxn] ) - abs( (*forkOff[1])[!drxn] );
        int same = off ? abs( (*off)[!drxn] ) - diff : 200;
        diff = abs( diff );
        if ( !i || diff < best[0] ) best[0] = diff;
        if ( !i || same < best[1] ) best[1] = same;
    }
    if ( best[1] < 30 && best[0] >= 200 ) return false;
//    if ( best[0] >= 30 || best[1] <= 200 ) return false;
    for ( Node* fork : forks[1].nodes )
    {
        fork->clearPaired( true );
        Node* forkClone = new Node( fork, nodes, fork->drxn_, false );
        vector<Edge> edges = fork->edges_[!drxn];
        for ( Edge& e : edges )
        {
            if ( !fwd[0].find( e.node ) ) continue;
            e.node->addEdge( forkClone, e.ol, drxn, false, e.isLeap );
            e.node->removeEdge( fork, drxn );
            fork->removeEdge( e.node, !drxn );
        }
    }
    if ( best[0] < 30 && best[1] >= 200 )
    {
        absorbClone( clone );
        nodes.erase( clone );
        nodes.test( true );
    }
    else if ( forks[1].empty() ) return false;
    for ( Node* fork : forks[1].nodes ) if ( fork->verified_ ) fork->setVerified();
    return true;
}

void Node::pruneMismatched( NodeRoll& nodes )
{
    NodeRoll seed = nodes.getGraph( 2 );
    Nodes tested( seed.nodes );
    for ( Node* node : seed.nodes )
    {
        for ( int d : { 0, 1 } ) for ( Edge& e : node->edges_[d] ) e.node->pruneMismatched( nodes, tested, d );
    }
}

void Node::pruneMismatched( NodeRoll& nodes, Nodes& tested, bool drxn )
{
    if ( tested.find( this ) ) return;
    for ( Edge& e : edges_[!drxn] ) if ( !e.node->bad_ && !tested.find( e.node ) ) return;
    tested += this;
    vector<int32_t> offsets;
    for ( Edge& e : edges_[!drxn] ) if ( !e.node->bad_ ) offsets.push_back( getOffset( e.node, !drxn ) );
    sort( offsets.begin(), offsets.end(), [&]( int32_t a, int32_t b ){ return drxn ? a < b : a > b; } );
    assert( abs( offsets[0] - offsets.back() ) < 200 );
    if ( offsets[0] ) offset( offsets[0] );
    for ( Edge& e : edges_[drxn] ) e.node->pruneMismatched( nodes, tested, drxn );
}

void Node::prunePaths( Querier& bwt, NodeRoll& nodes )
{
    for ( int again = 1; again-- > 0; )
    {
        pruneBlunt( bwt, nodes );
        Nodes tested[2];
        for ( int d : { 0, 1 } ) for ( Node* node : nodes.getGraph( 2 ).nodes ) if ( nodes.find( node ) && node->prunePaths( bwt, nodes, tested[d], d ) ) again = 1;
        reverify( nodes );
        
    }
    Node::merge( nodes );
}

bool Node::prunePaths( Querier& bwt, NodeRoll& nodes, Nodes& tested, bool drxn )
{
    if ( bad_ || !tested.add( this ) ) return false;
    
    bool claimed = false;
    
    int32_t dist = max( params.readLen, params.maxPeMean - params.readLen );
    if ( prepFork( bwt, nodes, dist, !drxn ) )
    {
        bool tryAgain = !( claimed = ClaimNode::claim( this, nodes, !drxn ) );
        if ( ends_[1] == -3845 ) tryAgain = false;
        if ( tryAgain )
        {
            ClaimFork cf( this, nodes, !drxn );
            if ( cf.claimed_ )
            {
                assert( false );
                int x = 0;
            }
//            if ( cf.claimed_ ) claimed = true;
//            vector<Node*> clones;
//            if ( cloned_ ) for ( Node* clone : cloned_->nodes ) if ( cf.cloned_.find( clone ) ) clones.push_back( clone );
//            for ( Node* clone : clones ) if ( clone->prunePaths( bwt, nodes, tested, drxn ) ) claimed = true;
        }
    }
    
    vector<Edge> altEdges = getAltEdges( !drxn );
    
    if ( !bad_ && verified_ ) for ( Edge& e : edges( drxn ) ) if ( e.node->prunePaths( bwt, nodes, tested, drxn ) ) claimed = true;
    
    if ( !bad_ && !verified_ && !claimed && prepBranch( dist, 2, drxn ) )
    {
        cout << "Branch: " << ends_[!drxn] << endl;
        int x = 0;
    }
    
    return claimed;
}

void Node::prunePrep( Querier& bwt, NodeRoll& nodes )
{
    NodeRoll seed = nodes.getGraph( 2 );
    Nodes tested( seed.nodes );
    for ( Node* node : seed.nodes )
    {
        for ( Edge& e : node->edges_[0] ) e.node->prunePrep( bwt, nodes, tested, 0 );
        for ( Edge& e : node->edges_[1] ) e.node->prunePrep( bwt, nodes, tested, 1 );
    }
    tested = Nodes( seed.nodes );
    for ( Node* node : seed.nodes )
    {
        for ( Edge& e : node->edges_[0] ) e.node->prunePrep( nodes, tested, 0 );
        for ( Edge& e : node->edges_[1] ) e.node->prunePrep( nodes, tested, 1 );
    }
}

void Node::prunePrep( Querier &bwt, NodeRoll &nodes, Nodes& tested, bool drxn )
{
    if ( !tested.add( this ) ) return;
    for ( Edge& e : edges_[!drxn] )
    {
        if ( !e.node->bad_ ) continue;
        e.node->extendFork( bwt, nodes, params.maxPeMean, 10, !drxn );
        e.node->extendFork( bwt, nodes, params.maxPeMean, 10, drxn );
        if ( e.node->bad_ ) e.node->verifyFork( params.maxPeMean, false, !drxn );
        else e.node->verifyFill();
    }
    if ( !verified_ ) return;
    for ( int i = 0; i < edges_[drxn].size(); i++ )
    {
        if ( !edges_[drxn][i].node->verified_ ) edges_[drxn][i].node->extendFork( bwt, nodes, params.readLen, 10, drxn );
    }
    for ( int i = 0; i < edges_[drxn].size(); i++ ) edges_[drxn][i].node->prunePrep( bwt, nodes, tested, drxn );
}

void Node::prunePrep( NodeRoll &nodes, Nodes& tested, bool drxn )
{
    if ( !verified_ || !tested.add( this ) ) return;
    for ( Edge& e : edges_[!drxn] ) if ( !e.node->verified_ && !e.node->bad_ ) e.node->verifyFill();
    for ( Edge& e : edges_[drxn] ) e.node->prunePrep( nodes, tested, drxn );
//    for ( int d : { 0, 1 } ) for ( Edge& e : edges_[d] ) assert( !e.node->isContinue( d ) );
}

void Node::setStrong( Nodes& strong, int readCount, int readLimit, bool drxn )
{
    readCount += countReads( true );
    if ( edges_[!drxn].empty() && cloned_ )
    {
        int maxFwd = 0;
        for ( Edge& e : getAltEdges( !drxn ) ) maxFwd = max( maxFwd, e.node->setStrong( readLimit - readCount, !drxn ) );
        readCount += maxFwd;
    }
    if ( readCount >= readLimit ) strong.fill( this, drxn, true, true );
    else for ( Edge& e : edges_[drxn] ) if ( !strong.find( e.node ) ) e.node->setStrong( strong, readCount, readLimit, drxn );
}

int Node::setStrong( int readLimit, bool drxn )
{
    int readCount = countReads( true ), maxFwd = 0;
    if ( ( readLimit -= readCount ) <= 0 ) return readCount;
    for ( Edge& e : getAltEdges( drxn ) ) if ( ( maxFwd = max( maxFwd, setStrong( readLimit, drxn ) ) ) >= readLimit ) break;
    return readCount + maxFwd;
}

