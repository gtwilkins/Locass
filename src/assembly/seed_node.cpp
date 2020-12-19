/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "seed_node.h"

SeedNode::SeedNode( Node* node, Coords* coords, int i, int j )
: node_( node ), created_( !coords ), ended_( !coords )
{
    mapped_[0][0] = i;
    mapped_[0][1] = j;
    mapped_[1][0] = coords ? coords->coords[0] : node->ends_[0];
    mapped_[1][1] = coords ? coords->coords[1] : node->ends_[1];
}

bool SeedNode::setEdges( string& s, vector<Node*>& ends, SeedNode& seed, NodeRoll& nodes, int ol, bool drxn )
{
    assert( !ends.empty() );
    bool edged = false;
    int32_t coords[2];
    Nodes block;
    for ( Node* node : ends ) for ( Node* clone : node->clones() ) block += clone;
    for ( Node* node : ends ) for ( Edge& e : node->edges_[drxn] ) block += e.node;
    
    vector<Edge> edges;
    int diff = drxn ? s.size() - seed.mapped_[0][1] : seed.mapped_[0][0];
    for ( int i = 0; i < nodes.size(); i++ ) if ( !block.find( nodes[i] ) && mapSeqEnd( s, nodes[i]->seq_, ol, coords, drxn ) )
    {
        int32_t cut = coords[!drxn] + nodes[i]->ends_[0];
        if ( !nodes[i]->getNextReadCoord( cut, !drxn, drxn ) ) continue;
        coords[!drxn] = cut - nodes[i]->ends_[0];
        vector<Node*> edged = ( nodes[i]->ends_[!drxn] != cut ? nodes[i]->splitNode( cut, nodes, drxn ) : vector<Node*>{ nodes[i] } );
        for ( Node* node : edged ) if ( block.add( node ) ) edges.push_back( Edge( node, coords[1] - coords[0] - diff, false ) );
    }
    for ( Node* node : ends ) for ( Edge& e : edges ) node->addEdge( e.node, e.ol, drxn, false, e.ol < 1 );
    cout << ( edges.empty() ? "Unedged seed" : to_string( edges.size() ) + " edges seed" ) << endl;
    return !edges.empty();
}

//bool SeedNode::setEdges( vector<Node*>& ends, NodeRoll& nodes, int ol, bool drxn )
//{
//    assert( !ends.empty() );
//    bool edged = false;
//    int32_t coords[2];
//    Nodes block;
//    for ( Node* node : ends ) for ( Node* clone : node->clones() ) block += clone;
//    for ( Node* node : ends ) for ( Edge& e : node->edges_[drxn] ) block += e.node;
//    
//    vector<Edge> edges;
//    for ( int i = 0; i < nodes.size(); i++ ) if ( !block.find( nodes[i] ) && mapSeqEnd( ends[0]->seq_, nodes[i]->seq_, ol, coords, drxn ) )
//    {
//        int32_t cut = coords[!drxn] + nodes[i]->ends_[0];
//        if ( !nodes[i]->getNextReadCoord( cut, !drxn, drxn ) ) continue;
//        coords[!drxn] = cut - nodes[i]->ends_[0];
//        vector<Node*> edged = ( nodes[i]->ends_[!drxn] != cut ? nodes[i]->splitNode( cut, nodes, drxn ) : vector<Node*>{ nodes[i] } );
//        for ( Node* node : edged ) if ( block.add( node ) ) edges.push_back( Edge( node, coords[1] - coords[0], false ) );
//    }
//    for ( Node* node : ends ) for ( Edge& e : edges ) node->addEdge( e.node, e.ol, drxn, false );
//    cout << ( edges.empty() ? "Unedged seed" : to_string( edges.size() ) + " edges seed" ) << endl;
//    return !edges.empty();
//}

vector<SeedNode> SeedNode::setNodes( Querier& bwt, NodeRoll& nodes, string& s, int drxn )
{
    vector<SeedNode> seeds;
    vector<ReadId> ids;
    vector<int32_t> mapped[2];
    vector<bool> redundant;
    bwt.mapSequence( s, ids, mapped );
    for ( int i = 0; i < ids.size(); )
    {
        for ( int j = i; j < ids.size() && mapped[1][j] <= mapped[1][i]; j++ ) assert( mapped[0][i] <= mapped[0][j] );
        for ( int j = i; j < ids.size() && mapped[1][j] <= mapped[1][i]; j++ ) redundant.push_back( mapped[0][i] < mapped[0][j] || mapped[1][j] < mapped[1][i] );
        i = redundant.size();
    }
    
    int len = 0, mark = 0;
    bool good = false;
    Coords* coords;
    
    for ( int i = 0; i <= ids.size(); i++ ) if ( i == ids.size() || len < mapped[1][i] )
    {
        vector<SeedNode> seeded;
        bool bad = i == ids.size();
        if ( !bad ) len = mapped[1][i];
        
        // This read is already used in a node
        if ( !bad ) for ( Node* node : nodes.nodes ) if ( ( coords = node->getRead( ids[i] ) ) && !coords->redundant )
        {
            bool added = false;
            for ( SeedNode& sn : seeded ) if ( node->isClone( sn.node_ ) && ( added = true ) ) break;
            for ( SeedNode& sn : seeds ) if ( !added && !sn.ended_ && node->isClone( sn.node_ ) && ( added = true ) )
            {
                if ( node != sn.node_ ) continue;
                sn.mapped_[0][1] = mapped[1][i];
                sn.mapped_[1][1] = coords->coords[1];
            }
            if ( !added ) seeded.push_back( SeedNode( node, coords, mapped[0][i], mapped[1][i] ) );
            bad = true;
        }
        for ( SeedNode& sn : seeds ) if ( sn.mapped_[0][1] != len ) sn.ended_ = true;
        
        // This run of reads requires a new node
        if ( good && bad )
        {
            Node* seed = new Node( s, 0, drxn );
            for ( int j = mark; j < i; j++ ) seed->add( ids[j], mapped[0][j], mapped[1][j], redundant[j] );
            seed->recoil();
            nodes += seed;
            seeds.push_back( SeedNode( seed, NULL, mapped[0][mark], mapped[0][mark]+seed->size() ) );
            if ( drxn == 2 ) seed->setOrigin();
        }
        
        if ( good == bad ) mark = i;
        good = !bad;
        seeds.insert( seeds.end(), seeded.begin(), seeded.end() );
    }
    
    // Split already used nodes where necessary
    for ( SeedNode& sn : seeds ) assert( sn.node_->ends_[0] <= sn.mapped_[1][0] && sn.mapped_[1][1] <= sn.node_->ends_[1] );
    for ( int d : { 0, 1 } ) for ( SeedNode& sn : seeds ) if ( !sn.created_ ) sn.node_ = sn.node_->splitNode( sn.mapped_[1][!d], nodes, d )[0];
//    for ( SeedNode& sn : seeds ) if ( !sn.created_ ) for ( int d : { 0, 1 } )
//    {
//        if ( sn.mapped_[1][!d] != sn.node_->ends_[!d] ) sn.node_ = sn.node_->splitNode( sn.mapped_[1][!d], nodes, d )[0];
//    }
    
    return seeds;
}

void SeedNode::setPath( vector<SeedNode>& seeds, vector<Node*> ends[2], NodeRoll& nodes, int drxn )
{
    assert( !seeds.empty() );
    if ( seeds.size() == 1 )
    {
        for ( int d : { 0, 1 } ) ends[d] = vector<Node*>{ seeds[0].node_ };
        return;
    }
    
    for ( int i = 0; i+1 < seeds.size(); i++ ) assert( seeds[i].mapped_[0][0] < seeds[i+1].mapped_[0][0] && seeds[i].mapped_[0][1] < seeds[i+1].mapped_[0][1] );
    
    Nodes good;
    vector<int> ols, breaks;
    vector< vector<Node*> > edges[2];
    vector<bool> cloneable;
    int merges = 0;
    for ( int i = 0; i < seeds.size(); )
    {
        breaks.push_back( i );
        if ( i ) ols.push_back( seeds[i-1].mapped_[0][1] - seeds[i].mapped_[0][0] );
        
        if ( seeds[i].created_ )
        {
            for ( int d : { 0, 1 } ) edges[d].push_back( vector<Node*>{ seeds[i].node_ } );
            cloneable.push_back( false );
            i++;
            continue;
        }
        
        Nodes reached;
        for ( Node* node : seeds[i].node_->clones() ) reached += node;
        
        int j = i;
        for ( int good = 1; good-- > 0 && ++j < seeds.size(); )
        {
            for ( Node* node : seeds[j].node_->clones() ) for ( Edge& e : node->edges_[0] ) if ( reached.find( e.node ) && ( good = 1 ) ) reached += node;
        }
        for ( int k = j-1; --k > i; ) for ( Node* node : seeds[k].node_->clones() ) if ( reached.erase( node ) )
        {
            for ( Edge& e : node->edges_[1] ) if ( reached.find( e.node ) ) reached += node;
        }
        vector<Node*> edged[2];
        for ( Node* node : seeds[i].node_->clones() ) if ( reached.find( node ) ) edged[0].push_back( node );
        for ( Node* node : seeds[j-1].node_->clones() ) if ( reached.find( node ) ) edged[1].push_back( node );
        bool clone = true;
        for ( int d : { 0, 1 } ) edges[d].push_back( edged[d] );
        for ( int d : { 0, 1 } ) for ( Node* node : edged[d] ) if ( node->bad_ || node->edges_[d].empty() ) clone = false;
        cloneable.push_back( clone );
        if ( clone ) merges++;
        assert( reached.size() >= j-i );
        assert( j > i );
        i = j;
        good += reached;
    }
    breaks.push_back( seeds.size() );
    
    if ( drxn == 2 && !cloneable[0] && !cloneable.back() ) merges++;
    for ( int i = 0; merges > 1 && i < cloneable.size(); i++ ) if ( drxn ? cloneable[i] : cloneable.end()[-i-1] )
    {
        int j = drxn ? breaks[i] : breaks.end()[-i-2], k = drxn ? breaks[i+1] : breaks.end()[-i-1];
        for ( int ii = j; ii < k; ii++ )
        {
            seeds[ii].node_ = new Node( seeds[ii].node_, nodes, drxn, false );
            if ( drxn == 2 ) seeds[ii].node_->setOrigin();
        }
        for ( int ii = j+1; ii < k; ii++ ) seeds[ii-1].node_->addEdge( seeds[ii].node_, seeds[ii-1].mapped_[0][1] - seeds[ii].mapped_[0][0], 1, false );
        ( drxn ? edges[0][i] : edges[0].end()[-i-1] ) = vector<Node*>{ seeds[j].node_ };
        ( drxn ? edges[1][i] : edges[1].end()[-i-1] ) = vector<Node*>{ seeds[k-1].node_ };
        merges--;
    }
    
    assert( merges < 2 );
    for ( int i = 0; i < ols.size(); i++ ) for ( Node* l : edges[1][i] ) for ( Node* r : edges[0][i+1] ) l->addEdge( r, ols[i], 1, false );
    ends[0] = edges[0][0];
    ends[1] = edges[1].back();
}

bool SeedNode::seed( Querier& bwt, NodeRoll& nodes, string s, int lOl, int rOl, int drxn )
{
    vector<Node*> ends[2];
    vector<SeedNode> seeds = setNodes( bwt, nodes, s, drxn );
    setPath( seeds, ends, nodes, drxn );
//    lOl -= seeds[0].mapped_[0][0];
//    rOl -= s.size() - seeds.back().mapped_[0][1];
    if ( lOl > 0 ) setEdges( s, ends[0], seeds[0], nodes, lOl, 0 );
    if ( rOl > 0 ) setEdges( s, ends[1], seeds.back(), nodes, rOl, 1 );
//    if ( lOl > 0 ) setEdges( ends[0], nodes, lOl, 0 );
//    if ( rOl > 0 ) setEdges( ends[1], nodes, rOl, 1 );
    if ( drxn == 2 );
    else if ( drxn ) for ( int i = 0; i < seeds.size(); i++ ) seeds[i].node_->offsetNode( 1 );
    else for ( int i = seeds.size(); --i > 0; ) seeds[i].node_->offsetNode( 0 );
    
    for ( SeedNode& sn : seeds ) sn.node_->setCoverage();
    for ( SeedNode& sn : seeds ) sn.node_->bad_ = false;
    for ( SeedNode& sn : seeds ) sn.node_->setVerified();
    for ( SeedNode& sn : seeds ) sn.node_->verify( drxn );
//    if ( !rOl ) seeds.back().node_->extendFork( bwt, nodes, 1000, 10, 1 );
//    if ( !lOl ) seeds.back().node_->extendFork( bwt, nodes, 1000, 10, 0 );
    nodes.verify( 1 );
    
    
    return true;
}
