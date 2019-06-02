/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "node_claim.h"
#include "node.h"
#include <algorithm>
#include <unordered_set>

ClaimJoin::ClaimJoin( ClaimNode* fork, ClaimNode* branch, int32_t dist, bool drxn )
: dist( dist )
{
    node[!drxn] = fork;
    node[drxn] = branch;
    for ( int d : { 0, 1 } ) node[d]->joins_[!d].push_back( this );
}

ClaimJoin::~ClaimJoin()
{
    for ( int d : { 0, 1 } ) node[d]->joins_[!d].erase( remove( node[d]->joins_[!d].begin(), node[d]->joins_[!d].end(), this ), node[d]->joins_[!d].end() );
}

int ClaimJoin::diff( bool drxn )
{
    return drxn ? dist : -dist;
}

void ClaimJoin::fill( Node* node, NodeOffsets& offs, bool drxn )
{
    if ( offs.find( node ) && bridge.add( node ) )
    {
        for ( Edge& e : node->getAltEdges( drxn ) ) fill( e.node, offs, drxn );
    }
}

ClaimRepair::ClaimRepair( Node* l, Node* r )
{
    fork[0] = l;
    fork[1] = r;
    for ( int d : { 0, 1 } ) fill( fork[d], -fork[d]->size(), d );
}

void ClaimRepair::dupe()
{
    ClaimRepair cr( fork[0], fork[1] );
    unordered_map<Node*, int32_t> diffs[2]{ get( cr, 0 ), get( cr, 1 ) };
    
    Nodes repairs;
    for ( int d : { 0, 1 } ) for ( auto& no : offs[d] ) if ( diffs[d].find( no.first ) == diffs[d].end() ) for ( auto& np : no.first->hits_.pairs[!d] )
    {
        if ( !no.first->cloned_ && !np.first->cloned_ ) continue;
        auto it = diffs[!d].find( np.first );
        if ( it != diffs[!d].end() && ( no.second + it->second ) < np.second.maxLen + 500 )
        {
            if ( no.first->cloned_ ) repairs += no.first;
            if ( np.first->cloned_ ) repairs += np.first;
            break;
        }
    }
    
    for ( Node* node : repairs.nodes ) node->reverify();
}

void ClaimRepair::fill( Node* node, int32_t off, bool drxn )
{
    auto it = offs[drxn].insert( make_pair( node, off ) );
    if ( !it.second )
    {
        if ( it.first->second <= off ) return;
        it.first->second = off;
    }
    
    if ( node->cloned_ ) for ( Node* clone : node->cloned_->nodes ) if ( !clone->edges_[0].empty() || !clone->edges_[1].empty() )
    {
        if ( node->edges_[drxn].empty() || clone->edges_[!drxn].empty() ) fill( clone, off, drxn );
    }
    
    for ( Edge& e : node->edges_[drxn] ) fill( e.node, off + node->size() - e.ol, drxn );
    
}

unordered_map<Node*, int32_t> ClaimRepair::get( ClaimRepair& cr, bool drxn )
{
    unordered_map<Node*, int32_t> diffs;
    for ( auto& no : offs[drxn] )
    {
        auto it = cr.offs[drxn].find( no.first );
        if ( it == cr.offs[drxn].end() || abs( no.second - it->second ) > 200 ) diffs.insert( no );
    }
    return diffs;
}

void ClaimRepair::trim()
{
    ClaimRepair cr( fork[0], fork[1] );
    unordered_map<Node*, int32_t> diffs[2]{ get( cr, 0 ), get( cr, 1 ) };
    
    Nodes repairs;
    for ( int d : { 0, 1 } ) for ( auto& no : diffs[d] ) if ( no.first->cloned_ ) for ( auto& np : no.first->hits_.pairs[!d] )
    {
        auto it = diffs[!d].find( np.first );
        if ( it != diffs[!d].end() && ( no.second + it->second ) < np.second.maxLen + 500 )
        {
            repairs += no.first;
            break;
        }
    }
    for ( Node* node : repairs.nodes ) node->reverify();
}

ClaimRedundant::ClaimRedundant( vector<ClaimNode*> alts[2], Node* paired[2], int hits )
: score( hits )
{
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : alts[d] )
    {
        cn->redundant_[!d].push_back( this );
        nodes[d].push_back( cn );
    }
    used.push_back( make_pair( paired[0], paired[1] ) );
}

ClaimRedundant::~ClaimRedundant()
{
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : nodes[d] )
    {
        cn->redundant_[!d].erase( remove( cn->redundant_[!d].begin(), cn->redundant_[!d].end(), this ), cn->redundant_[!d].end() );
    }
}

bool ClaimRedundant::add( vector<ClaimNode*> alts[2], Node* paired[2], int hits )
{
    for ( int d : { 0, 1 } )
    {
        if ( alts[d].size() != nodes[d].size() ) return false;
        for ( ClaimNode* cn : nodes[d] ) if ( find( alts[d].begin(), alts[d].end(), cn ) == alts[d].end() ) return false;
    }
    used.push_back( make_pair( paired[0], paired[1] ) );
    score += hits;
    
    return true;
}

bool ClaimRedundant::disregard( Node* paired[2] )
{
    for ( pair<Node*, Node*>& up : used )
    {
        if ( up.first->isClone( paired[0] ) && up.second->isClone( paired[1] ) ) return true;
    }
    return false;
}


ClaimPairing::ClaimPairing( ClaimNode* l, ClaimNode* r, int32_t dist )
: diffs{ dist }, missed{ 0 }, paths( 1 ), score( 0 )
{
    node[0] = l;
    node[1] = r;
    for ( int d : { 0, 1 } ) node[d]->pairs_[!d].push_back( this );
}

ClaimPairing::~ClaimPairing()
{
    for ( int d : { 0, 1 } ) node[d]->pairs_[!d].erase( remove( node[d]->pairs_[!d].begin(), node[d]->pairs_[!d].end(), this ), node[d]->pairs_[!d].end() );
}

int ClaimPairing::hit( vector<int32_t>& hitDiffs )
{
    int miss = -score;
    for ( int32_t d : hitDiffs ) for ( int i = 0; i < diffs.size(); i++ ) if ( diffs[i] == d ) miss = max( miss, -missed[i] );
    assert( !miss );
    return miss + score;
}

ClaimScore::ClaimScore( Node* l, Node* r, ClaimPairing* cp, int32_t dist, int32_t est )
: pairs{ cp }, diffs( 1 )
{
    node[0] = l;
    node[1] = r;
    for ( int32_t diff : cp->diffs ) assert( ( dist - diff ) > 0 );
    for ( int32_t diff : cp->diffs ) diffs[0].push_back( abs( ( dist - diff ) - est ) );
}

void ClaimScore::addRedundant( vector<ClaimNode*> alts[2], int hits )
{
    assert( alts[0].size() == 1 || alts[1].size() == 1 );
    assert( alts[0].size() > 1 || alts[1].size() > 1 );
    for ( ClaimRedundant* cr : pairs[0]->node[0]->redundant_[1] ) if ( cr->add( alts, node, hits ) ) return;
    new ClaimRedundant( alts, node, hits );
}

bool ClaimScore::cull( int32_t cutoff )
{
    bool incomplete = false;
    for ( int i = 0; i < pairs.size(); i++ )
    {
        bool good = false, bad = false;
        for ( int32_t diff : diffs[i] ) ( diff < cutoff ? good : bad ) = true;
        assert( good );
        if ( good && bad ) incomplete = true;
        if ( !good ) diffs.erase( diffs.begin() + i );
        if ( !good ) pairs.erase( pairs.begin() + i-- );
    }
    assert( !incomplete );
    return !incomplete;
}

bool ClaimScore::redundant( vector<ClaimNode*>& claims, int32_t est, int32_t best, int32_t cutoff, int hits )
{
    for ( ClaimRedundant* cr : pairs[0]->node[0]->redundant_[1] ) if ( cr->disregard( node ) ) return true;
    if ( !node[0]->cloned_ && !node[1]->cloned_ ) return false;
    
    Nodes clones[2];
    node[0]->hits_.setRedundant( clones[0], clones[1], node[1], 1 );
    if ( clones[0].size() < 2 && clones[1].size() < 2 ) return false;
    if ( clones[0].size() > 1 && clones[1].size() > 1 ) return true;
    
    int32_t* off;
    vector<ClaimNode*> cns[2];
    vector<Node*> nodes[2];
    vector<int32_t> offs[2];
    for ( int d : { 0, 1 } ) for ( Node* clone : clones[d].nodes ) for ( ClaimNode* cn : claims ) if ( ( off = cn->get( clone ) ) && cn->isValid( clone, d ) )
    {
        cns[d].push_back( cn );
        nodes[d].push_back( clone );
        offs[d].push_back( *off );
    }
    
    vector<ClaimScore> scores;
    for ( int i = 0; i < cns[0].size(); i++ ) for ( int j = 0; j < cns[1].size(); j++ ) if ( nodes[0][i] != node[0] || nodes[1][j] != node[1] )
    {
        if ( cns[0][i]->isShared( nodes[0][i], cns[1][j], 1 ) ) return true;
        if ( cns[1][j]->isShared( nodes[1][j], cns[0][i], 0 ) ) return true;
        
        ClaimPairing* cp = cns[0][i]->getPairing( cns[1][j], 1 );
        if ( !cp ) continue;
        
        assert( !cp->node[0]->isInvalid( nodes[0][i], offs[1][j], cp->node[1], 1 ) );
        assert( !cp->node[1]->isInvalid( nodes[1][j], offs[0][i], cp->node[0], 0 ) );
        if ( cp->node[0]->isInvalid( nodes[0][i], offs[1][j], cp->node[1], 1 ) ) return true;
        if ( cp->node[1]->isInvalid( nodes[1][j], offs[0][i], cp->node[0], 0 ) ) return true;
        
        for ( ClaimScore& cs : scores ) if ( cs.node[0] == nodes[0][i] && cs.node[1] == nodes[1][j] )
        {
            cs.pairs.push_back( cp );
            cs.diffs.push_back( vector<int32_t>() );
            for ( int32_t diff : cp->diffs ) cs.diffs.back().push_back( abs( ( offs[1][j] - offs[0][i] - diff ) - est ) );
            cp = NULL;
        }
        
        if ( cp ) scores.push_back( ClaimScore( nodes[0][i], nodes[1][j], cp, offs[1][j] - offs[0][i], est ) );
    }
    
    for ( ClaimScore& cs : scores ) for ( int i = 0; i < cs.diffs.size(); i++ ) for ( int32_t diff : cs.diffs[i] ) best = min( best, diff );
    
    cutoff += best;
    if ( !cull( cutoff ) ) return true;
    for ( ClaimScore& cs : scores ) if ( !cs.cull( cutoff ) ) return true;
    
    vector<ClaimNode*> alts[2];
    for ( int d : { 0, 1 } ) setClaims( alts[d], d );
    for ( ClaimScore& cs : scores ) for ( int d : { 0, 1 } ) cs.setClaims( alts[d], d );
    if ( alts[0].size() != 1 && alts[1].size() != 1 ) return true;
    
    bool distal[2]{ false, false };
    if ( setDistal( distal ) ) return true;
    for ( ClaimScore& cs : scores ) if ( cs.setDistal( distal ) ) return true;
    
    Nodes usedNodes;
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : alts[d] ) cn->setBase( usedNodes, true );
//    for ( int d : { 0, 1 } ) for ( Node* clone : clones[d].nodes ) if ( !usedNodes.find( clone ) ) assert( alts[0].size() == 1 && alts[1].size() == 1 );
    for ( int d : { 0, 1 } ) for ( Node* clone : clones[d].nodes ) if ( !usedNodes.find( clone ) ) return true;
    
    if ( alts[0].size() == 1 && alts[1].size() == 1 ) return false;
    
    addRedundant( alts, hits );
    
    return true;
}

void ClaimScore::score( vector<ClaimNode*>& claims, int32_t est, int32_t cutoff, int hits )
{
    if ( !setAlts( est ) ) return;
    
    int32_t best = diffs[0][0];
    for ( int i = 0; i < diffs.size(); i++ ) for ( int32_t diff : diffs[i] ) best = min( best, diff );
    
    if ( redundant( claims, est, best, cutoff, hits ) )
    {
        return;
    }
    
    if ( best > cutoff )
    {
        return;
    }
    cutoff += best;
    if ( !cull( cutoff ) && pairs.size() > 1 ) return;
    
    bool distal[2]{ false, false };
    if ( setDistal( distal ) ) return;
    
    if ( pairs.empty() ) assert( false );
    if ( pairs.empty() ) return;
    if ( pairs.size() > 1 )
    {
        vector<ClaimNode*> alts[2];
        for ( int d : { 0, 1 } ) setClaims( alts[d], d );
        addRedundant( alts, hits );
        return;
    }
    
    pairs[0]->score += hits;
    assert( pairs[0]->diffs.size() == diffs[0].size() );
    for ( int i = 0; i < diffs[0].size(); i++ ) if ( cutoff < diffs[0][i] ) pairs[0]->missed[i] += hits;
    for ( int i = 0; i < diffs[0].size(); i++ ) if ( cutoff < diffs[0][i] ) assert( false );
}

bool ClaimScore::setAlts( int32_t est )
{
    vector<ClaimNode*> cns[2]{ { pairs[0]->node[0] }, { pairs[0]->node[1] } };
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : cns[d][0]->alts_ ) if ( cn->get( node[d] ) ) cns[d].push_back( cn );
    for ( int i = 0; i < cns[0].size(); i++ ) for ( int j = !i; j < cns[1].size(); j++ )
    {
        ClaimPairing* paired = cns[0][i]->getPairing( cns[1][j], 1 );
        if ( !paired ) continue;
        int32_t* offs[2]{ cns[0][i]->get( node[0] ), cns[1][j]->get( node[1] ) };
        assert( offs[0] && offs[1] );
        int32_t dist = *offs[1] - *offs[0];
        pairs.push_back( paired );
        diffs.push_back( vector<int32_t>() );
        for ( int32_t diff : paired->diffs ) diffs.back().push_back( abs( ( dist - diff ) - est ) );
    }
    return true;
}

bool ClaimScore::setDistal( bool distal[2] )
{
    for ( ClaimPairing* cp : pairs ) for ( int d : { 0, 1 } ) if ( cp->node[d]->distal_.find( node[d] ) ) distal[d] = true;
    return distal[0] && distal[1];
}

void ClaimScore::setClaims( vector<ClaimNode*>& claims, bool drxn )
{
    for ( ClaimPairing* cp : pairs ) if ( find( claims.begin(), claims.end(), cp->node[drxn] ) == claims.end() ) claims.push_back( cp->node[drxn] );
}

ClaimDupe::ClaimDupe( ClaimNode* seed, vector<ClaimDupe*>& paths, bool drxn )
: path{ seed }, duped( false ), ignored( false )
{
    paths.push_back( this );
    extend( paths, drxn );
}

ClaimDupe::ClaimDupe( ClaimDupe* seed, ClaimNode* branch, vector<ClaimDupe*>& paths, bool drxn )
: path( seed->path ), duped( false ), ignored( false )
{
    paths.push_back( this );
    path.insert( drxn ? path.end() : path.begin(), branch );
    extend( paths, drxn );
}

bool ClaimDupe::create( ClaimNode* fork, vector<ClaimDupe*> paths[2], bool drxn )
{
    new ClaimDupe( fork, paths[!drxn], !drxn );
    for ( ClaimEdge& ce : fork->edges_[drxn] ) new ClaimDupe( ce.node, paths[drxn], drxn );
    
    unordered_set<ClaimNode*> shared;
    for ( int d : { 0, 1 } ) for ( int i = 0; i+1 < paths[d].size(); i++ ) for ( int j = i+1; j < paths[d].size(); j++ )
    {
        for ( ClaimNode* cn : paths[d][i]->pathed ) if ( paths[d][j]->pathed.find( cn ) != paths[d][j]->pathed.end() ) shared.insert( cn );
    }
    
    return paths[0].size() > 1 && paths[1].size() > 1;
}

bool ClaimDupe::dump( ClaimNode* fork, ClaimDupe* l, ClaimDupe* r, vector<ClaimDupe*> paths[2], int& pairing )
{
    ClaimDupe* cds[2]{ l, r };
    unordered_set<ClaimNode*> shared, tar[2], alt[2];
    int misses = 0, hits[2]{0};
    
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : paths[d] ) if ( cd != cds[d] ) for ( ClaimNode* cn : cd->pathed )
    {
        ( cds[d]->pathed.find( cn ) != cds[d]->pathed.end() ? shared : alt[d] ).insert( cn );
    }
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : cds[d]->path ) if ( shared.find( cn ) == shared.end() ) tar[d].insert( cn );
    
    for ( ClaimNode* cn : tar[0] ) for ( ClaimPairing* cp : cn->pairs_[1] ) if ( tar[1].find( cp->node[1] ) != tar[1].end() )
    {
        vector<int32_t> diffs;
        if ( !cn->reach( cp->node[1], fork, 0, diffs, false, 1 ) ) continue;
        assert( cp->hit( diffs ) );
        pairing += cp->hit( diffs );
    }
    if ( pairing ) return false;
    
    unordered_set<ClaimRedundant*> used[2];
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : tar[d] )
    {
        for ( ClaimPairing* cp : cn->pairs_[!d] ) if ( alt[!d].find( cp->node[!d] ) != alt[!d].end() )
        {
            vector<int32_t> diffs;
            if ( !cn->reach( cp->node[!d], fork, 0, diffs, false, !d ) ) continue;
            hits[d] += cp->hit( diffs );
            assert( cp->hit( diffs ) );
        }
        for ( ClaimRedundant* cr : cn->redundant_[!d] ) if ( used[d].insert( cr ).second )
        {
            bool good = true;
            for ( ClaimNode* cnp : cr->nodes[d] ) if ( tar[d].find( cnp ) == tar[d].end() ) good = false;
            for ( ClaimNode* cnp : cr->nodes[!d] ) if ( alt[!d].find( cnp ) == alt[!d].end() ) good = false;
            if ( !good ) continue;
            hits[d] += cr->score;
        }
    }
    
    if ( !hits[0] || !hits[1] ) return false;
    assert( !misses );
    
    misses /= 2;
    bool failed = min( hits[0], hits[1] ) - misses < 2;
    if ( !( hits[0] + hits[1] - misses > 7 ) ) failed = true;
    
    cout << "        " << ( failed ? "NOT DUPED: " : "DUPED: " ) << hits[0] << "-" << hits[1] << endl;
    
    if ( min( hits[0], hits[1] ) - misses < 2 ) return false;
    
    return hits[0] + hits[1] - misses > 7;
}

bool ClaimDupe::dupe( ClaimNode* fork, NodeRoll& nodes, bool drxn )
{
    vector<ClaimDupe*> paths[2];
    
    bool duped = false;
    
    if ( ClaimDupe::create( fork, paths, drxn ) )
    {
        cout << "    Duplicate attempt, coords: " << fork->getFork( 1 )->ends_[1] << ", paths: " << paths[0].size() << "-" << paths[1].size() << endl;
        duped = setDumps( fork, paths ) && setIgnores( paths ) && setDupes( fork->getFork( drxn ), paths, nodes, drxn );
        if ( !duped ) setSplits( paths );
        cout << "    " << ( duped ? "DUPLICATED!" : "FAILED!" ) << endl << endl;
    }
    else assert( false );
    
    if ( duped ) for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : paths[d] ) for ( ClaimNode* cn : cd->path ) if ( cn->setState( d ) ) break;
    
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : paths[d] ) delete cd;
    
    
    return duped;
}

void ClaimDupe::dupeNodes( unordered_set<ClaimNode*>& dupes, NodeRoll& cloned, Nodes& base, bool drxn )
{
    for ( int i = 0; i < path.size(); i++ )
    {
        ClaimNode* cn = drxn ? path[i] : path.end()[-i-1];
        if ( dupes.find( cn ) == dupes.end() ) return;
        if ( cn == cn->clone_ ) cn->clone_ = new ClaimNode( cn, cloned, base );
    }
}

bool ClaimDupe::dupePath( Node* fork, NodeRoll& cloned, Nodes& base, bool drxn )
{
    vector<ClaimDupe*> claim[2];
    unordered_set<ClaimNode*> alts;
    if ( !isDupe( claim, alts, drxn ) ) return false;
    
    ClaimRepair cr( fork, fork );
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : claim[d] ) cd->dupeNodes( alts, cloned, base, d );
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : claim[d] ) cd->dupeEdges( alts, d );
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : claim[d] ) cd->duped = true;
    cr.dupe();
    
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : claim[d] ) for ( ClaimNode* cn : cd->path ) if ( cn->clone_ != cn )
    {
        delete cn->clone_;
        cn->clone_ = cn;
    }
    
    return true;
}

void ClaimDupe::dupeEdges( unordered_set<ClaimNode*>& trims, bool drxn )
{
    ClaimNode* cn = NULL;
    for ( int i = 0; ( !cn || cn != cn->clone_ ) && i < path.size(); i++ )
    {
        cn = drxn ? path[i] : path.end()[-i-1];
        for ( ClaimEdge& ce : cn->edges_[!drxn] ) if ( trims.find( ce.node->clone_ ) == trims.end() )
        {
            cn->clone_->getFork( !drxn )->addEdge( ce.node->clone_->getFork( drxn ), ce.ol, !drxn, false, ce.isLeap );
            cn->clone_->addEdge( ce.node->clone_, ce.ol, ce.diff, ce.isLeap, !drxn );
        }
    }
    assert( cn == cn->clone_ );
    for ( int i = 0; i < cn->edges_[!drxn].size(); i++ ) if ( trims.find( cn->edges_[!drxn][i].node ) != trims.end() )
    {
        assert( cn->removeEdge( cn->edges_[!drxn][i--].node, !drxn ) );
    }
    cn = drxn ? path[0] : path.back();
    if ( cn == cn->clone_ || !cn->clone_->edges_[!drxn].empty() ) return;
//    for ( ClaimEdge& ce : cn->edges_[!drxn] ) if ( trims.find( ce.node->clone_ ) == trims.end() )
//    {
//        
//    }
    assert( false );
}

void ClaimDupe::extend( vector<ClaimDupe*>& paths, bool drxn )
{
    for ( ClaimNode* cn = drxn ? path.back() : path[0]; !cn->edges_[drxn].empty(); )
    {
        for ( int i = 1; i < cn->edges_[drxn].size(); i++ ) new ClaimDupe( this, cn->edges_[drxn][i].node, paths, drxn );
        cn = cn->edges_[drxn][0].node;
        path.insert( drxn ? path.end() : path.begin(), cn );
    }
    
    for ( ClaimNode* cn : path )
    {
        pathed.insert( cn );
        for ( ClaimJoin* cj : cn->joins_[drxn] ) fill( cj->node[drxn], drxn );
    }
}

void ClaimDupe::fill( ClaimNode* cn, bool drxn )
{
    if ( !pathed.insert( cn ).second ) return;
    for ( ClaimEdge& ce : cn->edges_[drxn] ) fill( ce.node, drxn );
    for ( ClaimJoin* cj : cn->joins_[drxn] ) fill( cj->node[drxn], drxn );
}

bool ClaimDupe::isDupe( vector<ClaimDupe*>& paths )
{
    assert( !kept.empty() );
    if ( dumped.empty() || duped ) return false;
    
    for ( ClaimDupe* cd : paths ) if ( cd != this && cd->dumped.size() >= kept.size() )
    {
        bool doDupe = false;
        for ( ClaimDupe* cdp : dumped ) if ( !cdp->duped ) doDupe = true;
        for ( ClaimDupe* cdp : kept ) if ( find( cd->dumped.begin(), cd->dumped.end(), cdp ) == cd->dumped.end() ) doDupe = false;
        if ( doDupe ) return true;
    }
    
    return false;
}

bool ClaimDupe::isDupe( vector<ClaimDupe*> claim[2], unordered_set<ClaimNode*>& alts, bool drxn )
{
    if ( kept.empty() || dumped.empty() || duped ) return false;
    claim[drxn] = kept;
    claim[!drxn] = kept[0]->kept;
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : claim[d] )
    {
        if ( cd->kept.size() != claim[!d].size() ) return false;
        for ( int i = 0; i < cd->kept.size(); i++ ) if ( cd->kept[i] != claim[!d][i] ) return false;
    }
    
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : claim[d][0]->dumped ) if ( !cd->duped ) for ( ClaimNode* cn : cd->path ) alts.insert( cn );
    
    return true;
}

bool ClaimDupe::setDumps( ClaimNode* fork, vector<ClaimDupe*> paths[2] )
{
    bool anyDumped = false;
    for ( ClaimDupe* l : paths[0] ) for ( ClaimDupe* r : paths[1] )
    {
        int pairing = 0;
        bool dumped = dump( fork, l, r, paths, pairing );
        ( dumped ? l->dumped : l->kept ).push_back( r );
        ( dumped ? r->dumped : r->kept ).push_back( l );
        if ( !dumped ) l->pairings.push_back( pairing );
        if ( !dumped ) r->pairings.push_back( pairing );
        if ( dumped ) anyDumped = true;
    }
    return anyDumped;
}

bool ClaimDupe::setDupes( Node* fork, vector<ClaimDupe*> paths[2], NodeRoll& nodes, bool drxn )
{
    Nodes base;
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : paths[d] ) if ( !cd->kept.empty() ) for ( ClaimNode* cn : cd->path ) cn->setBase( base, false );
    
    NodeRoll cloned;
    bool duped = false;
    for ( ClaimDupe* cd : paths[!drxn] ) if ( cd->dupePath( fork, cloned, base, drxn ) ) duped = true;
    
    for ( Node* clone : cloned.nodes )
    {
        assert( nodes.add( clone ) );
        clone->reverify();
    }
    
    if ( duped ) assert( !cloned.empty() );
    
    return duped;
}

bool ClaimDupe::setIgnores( vector<ClaimDupe*> paths[2] )
{
    unordered_set<ClaimDupe*> ignore[2], include[2];
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : paths[d] ) ( cd->dumped.empty() ? ignore : include )[d].insert( cd );
    
    if ( include[0].size() < 2 || include[1].size() < 2 ) return false;
    if ( ignore[0].empty() && ignore[1].empty() ) return true;
    
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : ignore[d] ) for ( int pairing : cd->pairings ) if ( pairing > 1 ) assert( false );
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : ignore[d] ) for ( int pairing : cd->pairings ) if ( pairing > 1 ) return false;
    
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : ignore[d] )
    {
        int maxLen = 0, minLen = cd->path.size();
        for ( ClaimDupe* alt : include[d] ) 
        {
            int i = 0, limit = min( cd->path.size(), alt->path.size() );
            while ( i < limit && ( d ? cd->path[i] == alt->path[i] : cd->path.end()[-i-1] == alt->path.end()[-i-1] ) ) i++;
            assert( i < limit );
            minLen = min( minLen, i );
            maxLen = max( maxLen, i );
        }
        if ( maxLen != minLen )
        {
            if ( cd->path[maxLen]->getFork( !d )->isBlunt( 0, 5, d ) )
            {
                continue;
            }
            return false;
        }
    }
    
    for ( int d : { 0, 1 } ) if ( !ignore[d].empty() ) for ( ClaimDupe* cd : include[!d] )
    {
        bool unclaimed = true;
        for ( int i = 0; i < cd->kept.size(); i++ ) if ( include[d].find( cd->kept[i] ) != include[d].end() && cd->pairings[i] > 4 ) unclaimed = false;
        if ( unclaimed )
        {
            return false;
        }
    }
    
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : ignore[d] ) cd->kept.clear();
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : include[d] ) for ( int i = 0; i < cd->kept.size(); i++ )
    {
        if ( ignore[!d].find( cd->kept[i] ) == ignore[!d].end() ) continue;
        cd->kept.erase( cd->kept.begin() + i );
        cd->pairings.erase( cd->pairings.begin() + i-- );
    }
    
    return true;
}

void ClaimDupe::setSplits( vector<ClaimDupe*> paths[2] )
{
    for ( int d : { 0, 1 } ) for ( ClaimDupe* cd : paths[d] )
    {
        int pairings = 0;
        for ( int pairing : cd->pairings ) if ( pairing > 5 ) pairings++;
        if ( pairings > 1 )
        {
            cout << "        SPLITABLE: ";
            for ( int i = 0; i < cd->pairings.size(); i++ ) cout << ( i ? "-" : "" ) << cd->pairings[i];
            cout << endl;
            ( d ? cd->path.back() : cd->path[0] )->split_ = true;
        }
    }
}

ClaimTrim::ClaimTrim( ClaimNode* fork, ClaimNode* branch, ClaimNode* bridge, bool drxn )
{
    forks[0] = drxn ? fork : branch;
    forks[1] = drxn ? branch : fork;
    for ( int d : { 0, 1 } )
    {
        ols[d] = 0;
        fill( forks[d], groups[d][0], d );
        for ( ClaimEdge& re : forks[!d]->edges_[d] ) if ( re.node != forks[d] && re.node != bridge )
        {
            fill( re.node, groups[d][1], d );
            ols[d] = max( ols[d], re.ol );
        }
        
        for ( ClaimNode* cn : groups[d][0] ) if ( groups[d][1].find( cn ) != groups[d][1].end() ) shared[d].insert( cn );
    }
}

void ClaimTrim::fill( ClaimNode* cn, unordered_set<ClaimNode*>& group, bool drxn )
{
    if ( !group.insert( cn ).second ) return;
    for ( ClaimEdge& ce : cn->edges_[drxn] ) fill( ce.node, group, drxn );
    for ( ClaimJoin* cj : cn->joins_[drxn] ) fill( cj->node[drxn], group, drxn );
}

void ClaimTrim::score( int l, int r, int& hits, int& unique )
{
//    unordered_set<ClaimNode*> pathable[2];
//    for ( int i : { 0, 1 } ) for ( ClaimNode* cn : groups[0][i] ) pathable[ i == l ].insert( cn );
//    for ( int i : { 0, 1 } ) for ( ClaimNode* cn : groups[1][i] ) pathable[ i == r ].insert( cn );
    bool path[2]{ !l, !r };
    for ( ClaimNode* cn : groups[0][l] ) for ( ClaimPairing* cp : cn->pairs_[1] ) if ( groups[1][r].find( cp->node[1] ) != groups[1][r].end() )
    {
        if ( shared[0].find( cn ) != shared[0].end() && shared[1].find( cp->node[1] ) != shared[1].end() ) continue;
        vector<int32_t> diffs;
        if ( !cn->reach( cp->node[1], forks, 0, diffs, path, false, 1 ) ) continue;
//        if ( !cn->reach( cp->node[1], 0, pathable, diffs, false, false, 1 ) ) continue;
        assert( !diffs.empty() );
        int pairScore = cp->hit( diffs );
        hits += pairScore;
        unique += pairScore;
    }
    unordered_set<ClaimRedundant*> shares;
    for ( ClaimNode* cn : groups[0][l] ) for ( ClaimRedundant* cr : cn->redundant_[1] )
    {
        bool good = true;
        for ( ClaimNode* cnp : cr->nodes[0] ) if ( groups[0][l].find( cnp ) == groups[0][l].end() ) good = false;
        for ( ClaimNode* cnp : cr->nodes[1] ) if ( groups[1][r].find( cnp ) == groups[1][r].end() ) good = false;
        if ( good ) for ( ClaimNode* cnl : cr->nodes[0] ) for ( ClaimNode* cnr : cr->nodes[1] )
        {
            vector<int32_t> diffs;
            if ( !cnl->reach( cnr, forks, 0, diffs, path, false, 1 ) ) good = false;
        }
        if ( good ) shares.insert( cr );
    }
    for ( ClaimRedundant* cr : shares )
    {
        hits += cr->score;
        unique += cr->score;
    }
}

bool ClaimTrim::trim( ClaimNode* fork, bool drxn )
{
    for ( ClaimEdge& ce : fork->edges_[drxn] )
    {
        ClaimEdge* edge[2]{ &ce, &ce };
        ClaimNode* branch = ce.node,* bridge = ce.node->isBridge() ? ce.node : NULL;
        if ( bridge )
        {
            edge[1] = &branch->edges_[drxn][0];
            branch = edge[1]->node;
        }
        if ( branch->edges_[!drxn].size() < 2 ) continue;
        
        cout << "    Trim attempt, coords: " << fork->getFork( drxn )->ends_[1] << ", edges: " << fork->edges_[drxn].size() << endl; 
        
        ClaimTrim ct( fork, branch, bridge, drxn );
        
        if ( ct.trim( min( edge[0]->ol, edge[1]->ol ) ) )
        {
            if ( bridge )
            {
                ( edge[1]->ol < edge[0]->ol ? fork : branch ) = bridge;
            }
            
            ClaimRepair cr( ( drxn ? branch : fork )->getFork( 1 ), ( drxn ? fork : branch )->getFork( 0 ) );
            assert( fork->removeEdge( branch, drxn ) );
            cr.trim();
        
            fork->setState( !drxn );
            branch->setState( drxn );
            if ( fork->edges_[drxn].size() > 1 ) trim( fork, drxn );
            cout << "    TRIMMED!" << endl << endl;
            return true;
        }
        else cout << "    FAILED!" << endl << endl;
    }
//    for ( ClaimEdge& ce : fork->edges_[drxn] ) if ( ce.node->edges_[!drxn].size() > 1 )
//    {
//        ClaimTrim ct( fork, ce.node, drxn );
//        if ( !ct.trim( ce ) ) continue;
//        
//        ClaimNode* branch = ce.node;
//        ClaimRepair cr( ( drxn ? branch : fork )->getFork( 1 ), ( drxn ? fork : branch )->getFork( 0 ) );
//        assert( fork->removeEdge( branch, drxn ) );
//        cr.trim();
//        
//        fork->setState( !drxn );
//        branch->setState( drxn );
//        if ( fork->edges_[drxn].size() > 1 ) trim( fork, drxn );
//        cout << "    TRIMMED!" << endl << endl;
//        return true;
//    }
    
    return false;
}

bool ClaimTrim::trim( int ol )
{
    for ( int d : { 0, 1 } ) if ( groups[d][0].empty() || groups[d][1].empty() ) return false;
    
    int hits[3]{0}, unique[3]{0}, diff = min( 100, min( 60, ols[0] - ol ) + min( 60, ols[1] - ol ) ) / 10;
    score( 0, 1, hits[0], unique[0] );
    score( 1, 0, hits[1], unique[1] );
    score( 0, 0, hits[2], unique[2] );
    
    // Filter out contended paths
    if ( !unique[0] || !unique[1] || unique[2] ) return false;
    
    bool failed = min( hits[0] + unique[0], hits[1] + unique[1] ) - hits[2] < 4;
    if ( !( hits[0] + hits[1] + unique[0] + unique[1] + diff - hits[2] > 8 ) ) failed = true;
    
    cout << "        " << ( failed ? "NOT TRIMMED: " : "TRIMMED: " ) << unique[0] << "(" << hits[0] << ")-" << unique[1] << "(" << hits[1] << ")" << endl;
    
    if ( min( hits[0] + unique[0], hits[1] + unique[1] ) - hits[2] < 4 ) return false;
    
    return hits[0] + hits[1] + unique[0] + unique[1] + diff - hits[2] > 8;
}

ClaimNode::ClaimNode( ClaimNode* cn, NodeRoll& cloned, Nodes& base )
: clone_( this ), split_( false )
{
    for ( int d : { 0, 1 } )
    {
        ends_[d] = cn->ends_[d];
        edged_[d] = cn->edged_[d];
    }
    
    for ( Node* node : cn->path_ )
    {
        vector<Edge> edges[2];
//        for ( int d : { 0 , 1 } ) for ( Edge& e : getEdges( node, d ) ) if ( !base.find( e.node ) ) edges[d].push_back( e );
        for ( int d : { 0 , 1 } ) for ( Edge& e : node->getAltEdges( d ) ) if ( !base.find( e.node ) ) edges[d].push_back( e );
        path_.push_back( new Node( node, cloned, node->drxn_, false ) );
        for ( int d : { 0 , 1 } ) for ( Edge& e : edges[d] ) path_.back()->addEdge( e, d, true );
    }
    for ( int i = 0; i+1 < path_.size(); i++ ) path_[i]->addEdge( cn->path_[i]->getEdge( path_[i+1], 1, true, true ), 1, true );
}

ClaimNode::ClaimNode( Node* node, vector<ClaimNode*>& claims, int32_t coord, bool orient, bool drxn )
: clone_( this ), path_{ node }, split_( false )
{
    claims.push_back( this );
    for ( int d : { 0, 1 } )
    {
        ends_[d] = coord;
        edged_[d] = false;
        for ( Node* fork : getForks( node, d ) ) offs_.insert( make_pair( fork, coord ) );
    }
    extend( claims, orient, drxn );
}

ClaimNode::~ClaimNode()
{
    reset();
    for ( int d : { 0, 1 } ) for ( ClaimEdge& ce : edges_[d] ) ce.node->removeEdge( this, !d, false );
}

void ClaimNode::addDistal( Node* node, bool drxn )
{
    if ( !get( node ) || !distal_.add( node ) ) return;
    for ( Edge& e : node->edges_[drxn] ) addDistal( e.node, drxn );
}

void ClaimNode::addEdge( ClaimNode* cn, int ol, int32_t diff, bool isLeap, bool drxn, bool reciprocate )
{
    for ( ClaimEdge& ce : edges_[drxn] ) if ( ce.node == cn ) return;
    edges_[drxn].push_back( ClaimEdge( cn, ol, diff, isLeap ) );
    if ( reciprocate ) cn->addEdge( this, ol, -diff, isLeap, !drxn, false );
}

bool ClaimNode::addEdge( Edge& e, vector<ClaimNode*>& claims, bool orient, bool drxn )
{
    int32_t coord = getCoord( NULL, e, ends_[drxn], orient, drxn );
    bool edged = false;
    for ( ClaimNode* cn : claims ) if ( cn != this && cn->findFork( e.node, !drxn ) )
    {
        addEdge( cn, e.ol, cn->ends_[!drxn] - coord, e.isLeap, drxn );
        edged = true;
    }
    if ( !edged ) addEdge( new ClaimNode( e.node, claims, coord, orient, drxn ), e.ol, 0, e.isLeap, drxn );
    
    return edged_[drxn] = true;
}

void ClaimNode::addJoin( ClaimNode* cn, Node* hang, NodeOffsets& offs, int32_t dist, bool drxn )
{
    if ( dist < -500 ) return;
    if ( abs( dist ) > 500 )
    {
        assert( false );
    }
    ClaimJoin* joined = NULL;
    for ( ClaimJoin* cj : joins_[drxn] ) if ( cj->node[drxn] == cn ) joined = cj;
    if ( joined ) joined->dist = min( joined->dist, dist );
    else joined = new ClaimJoin( this, cn, dist, drxn );
    joined->fill( hang, offs, !drxn );
}

//void ClaimNode::addJoin( ClaimNode* cn, int32_t diff, bool drxn, bool reciprocate )
//{
//    bool added = false;
//    for ( ClaimJoin& ce : joins_[drxn] ) if ( added = ( ce.node == cn ) )
//    {
//        if ( drxn ? ce.diff <= diff : diff <= ce.diff ) return;
//        ce.diff = diff;
//        break;
//    }
//    if ( !added ) joins_[drxn].push_back( ClaimEdge( cn, 0, diff, true ) );
//    
//    if ( reciprocate ) cn->addJoin( this, -diff, !drxn, false );
//}

bool ClaimNode::branch( vector<ClaimNode*>& claims, vector<ClaimNode*>& forks, bool orient, bool test )
{
    if ( test )
    {
        if ( forks.size() > 1 || forks[0]->isBranched( !orient ) ) return true;
        for ( ClaimEdge& ce : forks[0]->edges_[orient] ) if ( ce.node->edges_[!orient].size() > 1 ) return true;
        if ( !forks[0]->split( claims, orient, !orient ) ) return false;
        return true;
    }
    assert( !forks.empty() );
    bool success = false;;
    for ( ClaimNode* cn : forks )
    {
        bool branched = false;
        for ( ClaimEdge& ce : cn->edges_[orient] ) ce.node->setBranched( claims, cn->edges_[orient].size() > 1, orient, orient );
        for ( ClaimEdge& ce : cn->edges_[orient] ) if ( ce.node->edges_[!orient].size() > 1 ) branched = true;
        if ( cn->setBranched( claims, branched, orient, !orient ) || branched ) success = true;
    }
    return success;
    
//    if ( forks.size() > 1 || forks[0]->isBranched( !orient ) ) return true;
//    for ( ClaimEdge& ce : forks[0]->edges_[orient] ) if ( ce.node->edges_[!orient].size() > 1 ) return true;
//    if ( !forks[0]->split( claims, orient, !orient ) ) return false;
//    return true;
}

bool ClaimNode::claim( Node* fork, NodeRoll& nodes, bool drxn )
{
    cout << "Claim coords: " << fork->ends_[0] << " " << fork->ends_[1] << endl << endl;
    
    bool claimed = false;
    
    for ( int again = 1; again-- > 0; )
    {
        if ( fork->ends_[1] == -3844 )
        {
            int x = 0;
        }
        if ( fork->ends_[1] == -3322 )
        {
            int x = 0;
        }
        if ( fork->ends_[1] == -3954 )
        {
            int x = 0;
        }
        if ( fork->ends_[1] == -4794 )
        {
            int x = 0;
        }
        bool success = false;
        vector<ClaimNode*> claims, forks;
        if ( create( fork, claims, forks, drxn ) ) for ( int retry = 1; retry-- > 0; )
        {
//            branch( claims, forks, drxn, fork->ends_[1] == -3322 );
            branch( claims, forks, drxn, false );
            complete( claims, drxn );
            score( claims, drxn );
            for ( ClaimNode* cn : forks )
            {
                if ( success = cn->trim( 2 ) ) break;
                if ( success = cn->dupe( nodes, 2 ) ) break;
            }
            
            if ( success )
            {
                for ( ClaimNode* cn : forks ) if ( cn->edges_[drxn].size() > 1 ) again = 1;
                claimed = true;
            }
            else retry = resplit( forks, claims, drxn );
        }
        
        for ( ClaimNode* cn : claims ) delete cn;
    }
    
    cout << ( claimed ? "CLAIMED!" : "NOT CLAIMED!" ) << endl << endl;
    
    return claimed;
}

bool ClaimNode::complete( vector<ClaimNode*>& claims, bool orient )
{
    Nodes base;
    for ( ClaimNode* cn : claims ) cn->setBase( base, true );
    
    for ( int d : { 0, 1 } )
    {
        int32_t limit = 0;
        vector<ClaimNode*> ends;
        for ( ClaimNode* cn : claims ) if ( cn->edges_[d].empty() ) ends.push_back( cn );
        for ( ClaimNode* cn : ends ) limit = max( limit, abs( cn->ends_[d] ) );
        
        setEnds( ends, base, orient, d );
        setJoins( claims, base, limit, orient, d );
        for ( ClaimNode* cn : ends ) cn->setDistal( d );
    }
    for ( ClaimNode* cn : claims ) cn->testLoop( NULL );
}

bool ClaimNode::create( Node* fork, vector<ClaimNode*>& claims, vector<ClaimNode*>& forks, bool orient )
{
    bool looped = fork->cloned_ && fork->edges_[!orient].empty();
    if ( fork->cloned_ ) for ( Node* node : fork->cloned_->nodes ) if ( node->edges_[!orient].empty() && !node->edges_[orient].empty() ) looped = true;
    
    if ( looped ) for ( Node* node : fork->clones() )
    {
        if ( node->edges_[!orient].empty() ) continue;
        bool found = false;
        for ( ClaimNode* cn : claims ) if ( cn->findFork( node, orient ) ) found = true;
        assert( !found );
        if ( !found ) forks.push_back( new ClaimNode( node, claims, 0, orient, !orient ) );
    }
    else ( forks.push_back( new ClaimNode( fork, claims, 0, orient, !orient ) ) );
    
    for ( ClaimNode* cn : forks ) cn->edge( claims, orient, orient );
    for ( ClaimNode* cn : forks ) for ( ClaimEdge& ce : cn->edges_[orient] ) ce.node->edge( claims, orient, !orient );
    for ( int i = 0; i < forks.size(); i++ ) if ( forks[i]->edges_[orient].size() < 2 ) forks.erase( forks.begin() + i-- );
    
    for ( int again = 1; again-- > 0; )
    {
        Nodes branches[2];
        vector<ClaimNode*> unended[2];
        for ( int d : { 0, 1 } ) for ( ClaimNode* cn : claims ) if ( !cn->edged_[d] ) unended[d].push_back( cn );
        for ( int d : { 0, 1 } ) for ( ClaimNode* cn : unended[d] ) if ( !cn->edges_[d].empty() ) for ( Edge& e : cn->getEdges( NULL, d ) ) branches[d] += e.node;
        for ( int d : { 0, 1 } ) for ( ClaimNode* cn : unended[d] )
        {
            bool doEdge = false;
            for ( Node* node : cn->getForks( NULL, d ) ) if ( branches[!d].find( node ) ) again = doEdge = true;
            for ( Edge& e: cn->getEdges( NULL, d ) ) if ( branches[!d].find( e.node ) ) again = doEdge = true;
            if ( doEdge )
            {
                cn->edge( claims, orient, d );
            }
        }
    }
//    for ( int again = 1; again-- > 0; ) for ( int d : { 0, 1 } )
//    {
//        Nodes branches[2];
//        vector<ClaimNode*> ends[2];
//        unordered_set<ClaimNode*> added;
//        bool extended[2]{ false, false };
//        for ( ClaimNode* cn : claims ) if ( !cn->edged_[d] ) ends[ cn->edges_[d].empty() ].push_back( cn );
//        for ( int i : { 0, 1 } ) for ( ClaimNode* cn : ends[i] ) for ( Edge& e: cn->getEdges( NULL, orient ) ) branches[i] += e.node;
//        for ( int i : { 0, 1 } ) for ( ClaimNode* cn : ends[i] ) for ( Edge& e: cn->getEdges( NULL, orient ) ) if ( branches[!i].find( e.node ) )
//        {
//            cn->edge( claims, orient, d );
//            for ( ClaimEdge& ce : cn->edges_[d] ) added.insert( ce.node );
//            extended[i] = true;
//            again = 1;
//            break;
//        }
//        int x = 0;
//        for ( ClaimNode* cn : added ) cn->edge( claims, orient, !d );
//    }
    
    return !forks.empty();
}

bool ClaimNode::dupe( NodeRoll& nodes, int drxn )
{
    for ( int d : { 0, 1 } ) if ( drxn == 2 || d == drxn )
    {
        for ( ClaimEdge& ce : edges_[d] ) if ( ce.node->dupe( nodes, d ) ) return true;
    }
    
    for ( int d : { 0, 1 } ) if ( (bool)d == (bool)drxn && edges_[d].size() > 1 )
    {
        if ( isBranched( !d ) && ClaimDupe::dupe( this, nodes, d ) ) return true;
    }
    
    return false;
}

void ClaimNode::edge( vector<ClaimNode*>& claims, bool orient, bool drxn )
{
    if ( edged_[drxn] ) return;
    for ( Node* fork : getForks( NULL, drxn ) ) if ( !fork->isForkComplete( params.shortLen(), 20, drxn ) ) return;
    for ( Edge& e : getEdges( NULL, drxn ) ) addEdge( e, claims, orient, drxn );
    edged_[drxn] = true;
}

void ClaimNode::extend( vector<ClaimNode*>& claims, bool orient, bool drxn )
{
    vector<Edge> edges = getEdges( NULL, drxn );
    while ( edges.size() == 1 && edges_[drxn].empty() )
    {
        vector<Edge> rEdges = getEdges( edges[0].node, !drxn );
        if ( rEdges.size() > 1 && addEdge( edges[0], claims, orient, drxn ) ) return;
        ends_[drxn] = getCoord( NULL, edges[0], ends_[drxn], orient, drxn );
        path_.insert( drxn ? path_.end() : path_.begin(), edges[0].node );
        for ( Node* fork : getForks( NULL, drxn ) ) offs_.insert( make_pair( fork, ends_[drxn] ) );
        edges = getEdges( NULL, drxn );
        if ( edges.empty() ) edged_[drxn] = true;
    }
}

void ClaimNode::fill( Node* node, Nodes& block, int32_t dist, int32_t limit, bool orient, bool drxn )
{
    if ( block.find( node ) || abs( dist ) > abs( limit ) ) return;
    for ( Node* fork : getForks( node ) )
    {
        assert( !block.find( fork ) );
        auto it = offs_.insert( make_pair( fork, dist ) );
        if ( it.second ) continue;
        if ( it.first->second <= dist ) assert( fork == node );
        if ( it.first->second > dist ) it.first->second = dist;
        else return;
    }
    for ( Edge& e : getEdges( node, drxn ) ) fill( e.node, block, getCoord( node, e, dist, orient, drxn ), limit, orient, drxn );
}

void ClaimNode::fill( Node* node, Nodes& include, int32_t dist, bool orient, bool drxn )
{
    for ( Edge& e : getEdges( node, drxn ) ) if ( include.find( e.node ) )
    {
        int32_t coord = getCoord( node, e, dist, orient, drxn );
        auto it = offs_.insert( make_pair( e.node, coord ) );
        if ( !it.second && it.first->second <= coord ) continue;
        for ( Node* fork : getForks( e.node ) ) offs_[fork] = coord;
        fill( e.node, include, dist, orient, drxn );
    }
}

bool ClaimNode::findFork( Node* q, bool drxn )
{
    for ( Node* fork : getForks( NULL, drxn ) ) if ( q == fork ) return true;
    return false;
}

int32_t* ClaimNode::get( Node* node )
{
    auto it = offs_.find( node );
    return it == offs_.end() ? NULL : &it->second;
}

int32_t ClaimNode::getCoord( Node* fork, Edge& e, int32_t dist, bool orient, bool drxn )
{
    if ( !fork ) fork = getFork( drxn );
    return dist + ( ( ( drxn ? e.node : fork )->size() - e.ol ) * ( drxn ? 1 : -1 ) );
}

vector<Edge> ClaimNode::getEdges( Node* fork, bool drxn, bool blunt )
{
    vector<Edge> edges;
    bool isEnd = !fork;
    for ( Node* node : getForks( fork, drxn ) ) for ( Edge& e : node->edges_[drxn] )
    {
        for ( Node* branch : getForks( e.node, drxn ) )
        {
            bool added = !blunt && branch->isBlunt( 0, 3, drxn );
            for ( Edge& edge : edges ) if ( edge.node == branch ) added = true;
            if ( isEnd ) for ( ClaimEdge& ce : edges_[drxn] ) if ( ce.node->findFork( e.node, !drxn ) ) added = true;
            if ( !added ) edges.push_back( Edge( branch, e.ol, e.isLeap ) );
        }
    }
    return edges;
}

Nodes ClaimNode::getExts()
{
    Nodes exts;
    for ( const pair<Node*, int32_t>& no : offs_ ) if ( no.second < ends_[0] || ends_[1] < no.second ) exts += no.first;
    return exts;
}

Node* ClaimNode::getFork( bool drxn )
{
    return drxn ? path_.back() : path_[0];
}

vector<Node*> ClaimNode::getForks( Node* fork, bool drxn )
{
    if ( !fork ) fork = getFork( drxn );
    vector<Node*> forks;
    if ( !fork->cloned_ || !fork->edges_[drxn].empty() ) forks.push_back( fork );
    if ( fork->cloned_ ) for ( Node* node : fork->cloned_->nodes )
    {
        if ( node->edges_[drxn].empty() ) continue;
        if ( node->edges_[0].empty() && node->edges_[1].empty() ) continue;
        if ( fork->edges_[0].empty() || fork->edges_[1].empty() || node->edges_[!drxn].empty() ) forks.push_back( node );
    }
    assert( !forks.empty() );
    return forks;
}

vector<Node*> ClaimNode::getForks( Node* fork )
{
    vector<Node*> forks{ fork };
    if ( fork->cloned_ ) for ( Node* clone : fork->cloned_->nodes )
    {
        if ( clone->edges_[0].empty() && clone->edges_[1].empty() ) continue;
        if ( clone->edges_[0].empty() || clone->edges_[1].empty() || fork->edges_[0].empty() || fork->edges_[1].empty() ) forks.push_back( clone );
    }
    return forks;
}

vector< pair<Node*, int32_t> > ClaimNode::getHangs( bool orient, bool drxn )
{
    vector< pair<Node*, int32_t> > hangs;
    
    if ( !edged_[drxn] ) for ( Edge& e : getEdges( NULL, drxn, true ) )
    {
        bool used = false;
        for ( ClaimEdge& ce : edges_[drxn] ) if ( ce.node->findFork( e.node, !drxn ) ) used = true;
        if ( !used ) hangs.push_back( make_pair( e.node, getCoord( NULL, e, ends_[drxn], orient, drxn ) ) );
    }
    
    return hangs;
}

ClaimPairing* ClaimNode::getPairing( ClaimNode* cn, bool drxn )
{
    for ( ClaimPairing* cp : pairs_[drxn] ) if ( cp->node[drxn] == cn ) return cp;
    return NULL;
}

bool ClaimNode::isBranched( bool drxn )
{
    if ( !edged_[drxn] ) return false;
    if ( edges_[drxn].size() > 1 ) return true;
    for ( ClaimEdge& ce : edges_[drxn] ) if ( ce.node->isBranched( drxn ) ) return true;
    return false;
}

bool ClaimNode::isBridge()
{
    if ( path_.size() != 1 || !edged_[0] || !edged_[1] || edges_[0].size() != 1 || edges_[1].size() != 1 ) return false;
    return path_[0]->countReads( true ) < 2;
}

bool ClaimNode::isExtend( int32_t coord )
{
    return coord < ends_[0] || ends_[1] < coord;
}

bool ClaimNode::isForked( bool drxn )
{
    return edges_[drxn].size() > 1 || ( isBridge() && edges_[!drxn][0].node->edges_[drxn].size() > 1 );
}

bool ClaimNode::isInvalid( Node* node, int32_t off, ClaimNode* cn, bool drxn )
{
    if ( isShared( node, cn, drxn ) ) return true;
    if ( joins_[drxn].empty() || !edges_[drxn].empty() ) return false;
    if ( ends_[0] <= off && off <= ends_[1] ) return false;
    for ( ClaimJoin* cj : joins_[drxn] ) if ( cj->node[drxn]->reach( cn, drxn ) ) if ( cj->bridge.find( node ) ) return false;
    return true;
}

bool ClaimNode::isShared( Node* node, ClaimNode* cn, bool drxn )
{
    for ( ClaimShared& cs : shared_ ) if ( cs.shared.find( node ) && cs.node->reach( cn, drxn ) ) return true;
    return false;
}

bool ClaimNode::isValid( Node* node, bool drxn )
{
    if ( !node->cloned_ || ( !node->edges_[0].empty() && node->edges_[1].empty() ) ) return true;
    
    for ( Node* clone : node->cloned_->nodes ) if ( get( clone ) )
    {
        if ( !clone->edges_[0].empty() && !clone->edges_[1].empty() ) return false;
        if ( node->edges_[drxn].empty() && !clone->edges_[drxn].empty() ) return false;
    }
    
    return true;
}

bool ClaimNode::reach( ClaimNode* tar, bool drxn )
{
    if ( tar == this ) return true;
    for ( ClaimEdge& ce : edges_[drxn] ) if ( ce.node->reach( tar, drxn ) ) return true;
    for ( ClaimJoin* cj : joins_[drxn] ) if ( cj->node[drxn]->reach( tar, drxn ) ) return true;
    return false;
}

bool ClaimNode::reach( ClaimNode* tar, ClaimNode* forks[2], int32_t diff, vector<int32_t>& diffs, bool path[2], bool forked, bool drxn )
{
    if ( this == tar )
    {
        if ( !forked ) return false;
        if ( find( diffs.begin(), diffs.end(), diff ) == diffs.end() ) diffs.push_back( diff );
        return true;
    }
    
    for ( ClaimEdge& ce : edges_[drxn] )
    {
        bool pathed = forked;
        if ( forks[!drxn] == this ) pathed = path[!drxn] && ( ( forks[drxn] == ce.node ) == path[drxn] );
        if ( forks[drxn] == ce.node ) pathed = path[drxn] && ( ( forks[!drxn] == this ) == path[!drxn] );
        if ( !ce.node->reach( tar, forks, diff + ce.diff, diffs, path, pathed, drxn ) ) return false;
    }
    for ( ClaimJoin* cj : joins_[drxn] )
    {
        bool pathed = forked;
        assert( forks[!drxn] != this );
        assert( forks[drxn] != cj->node[drxn] );
        if ( forks[!drxn] == this ) pathed = path[!drxn] && !path[drxn] && ( forks[drxn] != cj->node[drxn] );
        if ( forks[drxn] == cj->node[drxn] ) pathed = path[drxn] && !path[!drxn] && ( forks[!drxn] != this );
        if ( !cj->node[drxn]->reach( tar, forks, diff + cj->diff( drxn ), diffs, path, pathed, drxn ) ) return false;
    }
    
    return true;
}

bool ClaimNode::reach( ClaimNode* tar, ClaimNode* fork, int32_t diff, vector<int32_t>& diffs, bool forked, bool drxn )
{
    if ( this == fork ) forked = true;
    if ( this == tar )
    {
        if ( !forked ) return false;
        if ( find( diffs.begin(), diffs.end(), diff ) == diffs.end() ) diffs.push_back( diff );
        return true;
    }
    
    for ( ClaimEdge& ce : edges_[drxn] ) if ( !ce.node->reach( tar, fork, diff + ce.diff, diffs, forked, drxn ) ) return false;
    for ( ClaimJoin* cj : joins_[drxn] ) if ( !cj->node[drxn]->reach( tar, fork, diff + cj->diff( drxn ), diffs, forked, drxn ) ) return false;
    
    return true;
}

bool ClaimNode::reach( ClaimNode* tar, int32_t diff, vector<int32_t>& diffs, bool drxn )
{
    if ( tar == this )
    {
        if ( find( diffs.begin(), diffs.end(), diff ) == diffs.end() ) diffs.push_back( diff );
        return true;
    }
    
    bool reached = false;
    for ( ClaimEdge& ce : edges_[drxn] ) if ( ce.node->reach( tar, diff + ce.diff, diffs, drxn ) ) reached = true;
    for ( ClaimJoin* cj : joins_[drxn] ) if ( cj->node[drxn]->reach( tar, diff + cj->diff( drxn ), diffs, drxn ) ) reached = true;
    return reached;
}

bool ClaimNode::removeEdge( ClaimNode* cn, bool drxn, bool reciprocate )
{
    bool removed = false;
    for ( int i = 0; i < edges_[drxn].size(); i++ ) if ( edges_[drxn][i].node == cn )
    {
        edges_[drxn].erase( edges_[drxn].begin() + i-- );
        removed = true;
    }
    if ( removed && reciprocate )
    {
        bool trimmed = false;
        Nodes forks( getForks( NULL, drxn ) ), branches( cn->getForks( NULL, !drxn ) );
        if ( clone_ != this ) forks -= clone_->getFork( drxn );
        for ( ClaimEdge& ce : edges_[drxn] ) branches -= ce.node->getFork( !drxn );
        for ( Node* node : forks.nodes ) for ( Node* branch : branches.nodes ) if ( node->removeEdge( branch, drxn, true ) ) trimmed = true;
        assert( trimmed );
        cn->removeEdge( this, !drxn, false );
    }
    return removed;
}

bool ClaimNode::removeJoin( ClaimNode* cn, bool drxn, bool reciprocate )
{
    bool removed = false;
    for ( int i = 0; i < joins_[drxn].size(); i++ ) if ( joins_[drxn][i]->node[drxn] == cn )
    {
        joins_[drxn].erase( joins_[drxn].begin() + i-- );
        removed = true;
    }
    if ( removed && reciprocate ) cn->removeJoin( this, !drxn, false );
    return removed;
}

void ClaimNode::reset()
{
    for ( int d : { 0, 1 } )
    {
        while ( !pairs_[d].empty() ) delete pairs_[d].back();
        while ( !redundant_[d].empty() ) delete redundant_[d].back();
        while ( !joins_[d].empty() ) delete joins_[d].back();
    }
    for ( ClaimNode* cn : alts_ ) cn->alts_.erase( remove( cn->alts_.begin(), cn->alts_.end(), this ), cn->alts_.end() );
    alts_.clear();
    distal_.clear();
}

bool ClaimNode::resplit( vector<ClaimNode*>& forks, vector<ClaimNode*>& claims, bool orient )
{
    bool splitable = false, didSplit = false;
    for ( ClaimNode* cn : claims ) if ( cn->split_ ) splitable = true;
    if ( !splitable )
    {
        for ( ClaimNode* cn : forks ) if ( !cn->isBranched( !orient ) && cn->split( claims, orient, !orient ) )
        {
            didSplit = true;
        }
    }
    
    if ( !splitable && !didSplit ) return false;
    
    for ( ClaimNode* cn : claims )
    {
        cn->reset();
        for ( auto it = cn->offs_.begin(); it != cn->offs_.end(); )
        {
            if ( cn->ends_[0] <= it->second && it->second <= cn->ends_[1] ) it++;
            else it = cn->offs_.erase( it );
        }
    }
    
    if ( !didSplit ) for ( int i = 0; i < claims.size(); i++ ) if ( claims[i]->split_ )
    {
        assert( claims[i]->edges_[0].empty() || claims[i]->edges_[1].empty() );
        for ( int d : { 0, 1 } ) if ( claims[i]->edges_[d].empty() )
        {
            assert( claims[i]->split_ );
            claims[i]->split_ = false;
            claims[i]->edge( claims, orient, d );
            if ( !claims[i]->edges_[d].empty() ) didSplit = true;
        }
        assert( !claims[i]->split_ );
    }
    
    return didSplit;
}

void ClaimNode::score( vector<ClaimNode*>& claims, bool orient )
{
    Nodes base;
    for ( ClaimNode* cn : claims ) for ( pair<Node*, int32_t> no : cn->offs_ ) base += no.first;
    
    for ( ClaimNode* cn : claims ) cn->setPairs( cn, 0, false, false );
    
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : claims ) if ( cn->edges_[d].size() > 1 || !cn->joins_[d].empty() )
    {
        for ( Node* node : cn->getForks( NULL, d ) ) node->verifyFork( params.maxPeMean, false, d );
    }
    
    for ( ClaimNode* cn : claims ) cn->setScores( claims, base );
//    for ( ClaimNode* cn : claims ) cn->setScores( claims, base, orient );
}

void ClaimNode::setBase( Nodes& base, bool distal )
{
    if ( distal ) for ( const pair<Node*, int32_t>& no : offs_ ) base += no.first;
    else for ( Node* node : path_ )
    {
        base += node;
        if ( node->cloned_ ) for ( Node* clone : node->cloned_->nodes ) if ( offs_.find( clone ) != offs_.end() ) base += clone;
    }
}

void ClaimNode::setBlocked( Node* node, bool drxn )
{
    if ( offs_.find( node ) == offs_.end() ) return;
    for ( Node* fork : getForks( node ) ) offs_.erase( fork );
    for ( Edge& e : getEdges( node, drxn ) ) setBlocked( e.node, drxn );
}

//void ClaimNode::setBlocks( Nodes& block, bool drxn )
//{
//    vector<ClaimNode*> forks;
//    setForks( forks, drxn );
//    assert( !forks.empty() );
//    int32_t limit = ends_[drxn];
//    for ( ClaimNode* cn : forks ) limit = drxn ? max( limit, cn->ends_[1]+500 ) : min( limit, cn->ends_[0]-500 );
//    for ( const pair<Node*, int32_t>& no : offs_ ) if ( drxn ? limit < no.second : no.second < limit ) block += no.first;
//}

bool ClaimNode::setBranched( vector<ClaimNode*>& claims, bool branched, bool orient, bool drxn )
{
    // branched = path has already branched, forked = path branches forward
    if ( !branched && !edged_[drxn] ) edge( claims, orient, drxn );
    if ( edges_[drxn].size() > 1 ) setForked( claims, orient, drxn );
    
    bool forked = edges_[drxn].size() > 1;
    for ( ClaimEdge& ce : edges_[drxn] ) if ( ce.node->setBranched( claims, branched || forked, orient, drxn ) ) forked = true;
    if ( forked ) setForked( claims, orient, !drxn );
    return forked;
}

void ClaimNode::setDistal( bool drxn )
{
    for ( const pair<Node*, int32_t>& no : offs_ ) if ( drxn ? ends_[1] < no.second : no.second < ends_[0] )
    {
        for ( Edge& e : no.first->edges_[!drxn] ) if ( !get( e.node ) && !e.node->isBlunt( 0, 3, !drxn ) )
        {
            addDistal( no.first, drxn );
        }
    }
}

void ClaimNode::setEnds( vector<ClaimNode*>& ends, Nodes& base, bool orient, bool drxn )
{
    if ( ends.size() < 2 ) return;
    
    int32_t maxExt = 0,* off;
    vector<int32_t> limits;
    for ( ClaimNode* cn : ends )
    {
        vector<ClaimNode*> forks;
        cn->setForks( forks, drxn );
        assert( !forks.empty() );
        int32_t limit = cn->ends_[drxn];
        for ( ClaimNode* fork : forks ) limit = drxn ? max( limit, fork->ends_[1]+500 ) : min( limit, fork->ends_[0]-500 );
        maxExt = max( maxExt, abs( limit - cn->ends_[drxn] ) );
        limits.push_back( limit );
    }
    if ( !maxExt ) return;
    
    Nodes exts;
    for ( int i = 0; i < ends.size(); i++ ) for ( pair<Node*, int32_t>& no : ends[i]->getHangs( orient, drxn ) )
    {
        ends[i]->fill( no.first, base, no.second, abs( ends[i]->ends_[drxn] - limits[i] )+500, orient, drxn );
        for ( pair<Node*, int32_t> no : ends[i]->offs_ ) if ( ends[i]->isExtend( no.second ) ) exts += no.first;
    }
    
    Nodes shared;
    for ( Node* node : exts.nodes ) for ( ClaimNode* cn : ends ) if ( off = cn->get( node ) ) cn->fill( node, exts, *off, orient, drxn );
    for ( Node* node : exts.nodes )
    {
        int found = 0;
        for ( ClaimNode* cn : ends ) if ( off = cn->get( node ) ) found++;
        if ( found == ends.size() ) for ( ClaimNode* cn : ends ) cn->offs_.erase( node );
        if ( found < ends.size() && found > 1 ) shared += node;
    }
    
    vector<ClaimNode*> alts;
    for ( int i = 0; i < ends.size(); i++ ) if ( ends[i]->setExts( shared, limits[i], drxn ) ) alts.push_back( ends[i] );
    for ( ClaimNode* cn : alts ) for ( ClaimNode* alt : alts ) if ( cn != alt ) cn->alts_.push_back( alt );
//    for ( int i = 0; i < ends.size(); i++ ) if ( !ends[i]->setExts( shared, limits[i], drxn ) ) ends.erase( ends.begin() + i-- );
//    for ( ClaimNode* cn : ends ) for ( ClaimNode* alt : ends ) if ( cn != alt ) cn->alts_.push_back( alt );
    
//    // Extend further
//    for ( ClaimNode* cn : ends ) for ( pair<Node*, int32_t>& no : cn->getHangs( orient, drxn ) ) cn->fill( no.first, base, no.second, limit+500, orient, drxn );
//    for ( ClaimNode* cn : ends ) for ( pair<Node*, int32_t> no : cn->offs_ ) exts += no.first;
//    for ( Node* node : exts.nodes ) for ( ClaimNode* cn : ends )
//    {
//        auto it = cn->offs_.find( node );
//        if ( it != cn->offs_.end() ) cn->fill( node, base, exts, it->second, orient, drxn );
//    }
//
//    // Catalogue all nodes remaining nodes to block
//    for ( ClaimNode* cn : ends ) cn->setBlocks( block, drxn );
//    for ( ClaimNode* cn : ends ) for ( pair<Node*, int32_t> no : cn->offs_ ) if ( !base.find( no.first ) )
//    {
//        for ( ClaimNode* alt : ends ) if ( cn != alt && alt->offs_.find( no.first ) != alt->offs_.end() ) block += no.first;
//    }
//
//    // Remove blocked
//    for ( Node* node : block.nodes ) for ( ClaimNode* cn : ends ) if ( cn->offs_.find( node ) != cn->offs_.end() ) cn->setBlocked( node, drxn );
}

bool ClaimNode::setExts( Nodes& shared, int32_t limit, bool drxn )
{
    Nodes keep;
    for ( const pair<Node*, int32_t>& no : offs_ ) if ( shared.find( no.first ) ) setKeep( no.first, keep, !drxn );
    for ( auto it = offs_.begin(); it != offs_.end(); )
    {
        if ( keep.find( it->first ) || ( drxn ? it->second <= limit : limit <= it->second ) ) it++;
        else it = offs_.erase( it );
    }
    return !keep.empty();
}

void ClaimNode::setForked( vector<ClaimNode*>& claims, bool orient, bool drxn )
{
    edge( claims, orient, drxn );
    for ( ClaimEdge& ce : edges_[drxn] )
    {
        ce.node->edge( claims, orient, !drxn );
        if ( ce.node->isBridge() ) ce.node->edges_[drxn][0].node->edge( claims, orient, !drxn );
    }
}

void ClaimNode::setForks( vector<ClaimNode*>& forks, bool drxn )
{
    if ( find( forks.begin(), forks.end(), this ) != forks.end() ) return;
    if ( edges_[drxn].size() > 1 ) forks.push_back( this );
    else for ( ClaimEdge& ce : edges_[!drxn] ) ce.node->setForks( forks, drxn );
}

void ClaimNode::setJoins( vector<ClaimNode*>& claims, Nodes& base, int32_t limit, bool orient, bool drxn )
{
    vector< pair<ClaimNode*, Nodes> > ends;
    for ( ClaimNode* cn : claims ) if ( cn->edges_[drxn].empty() )
    {
        Nodes ext = cn->getExts();
        if ( !ext.empty() ) ends.push_back( make_pair( cn, ext ) );
    }
    
    for ( ClaimNode* cn : claims ) if ( !cn->edged_[drxn] && !cn->edges_[0].empty() )
    {
        NodeOffsets offs;
        NodeOffset* off;
        for ( pair<Node*, int32_t>& hang : cn->getHangs( orient, drxn ) ) offs.fillNot( hang.first, base, hang.second, limit+500, orient, drxn, true );
        
        for ( ClaimNode* alt : claims ) if ( cn != alt && !alt->edged_[!drxn] )
        {
            for ( pair<Node*, int32_t>& hang : alt->getHangs( orient, !drxn ) )
            {
//                if ( off = offs.get( hang.first ) ) cn->addJoin( alt, hang.first, offs, hang.second - (*off)[0], drxn );
                if ( off = offs.get( hang.first ) ) cn->addJoin( alt, hang.first, offs, drxn ? hang.second - (*off)[0] : (*off)[0] - hang.second, drxn );
            }
        }
        
        if ( !offs.empty() && !cn->edges_[drxn].empty() ) for ( pair<ClaimNode*, Nodes>& ext : ends )
        {
            ClaimShared cs( cn );
            for ( Node* node : ext.second.nodes ) if ( offs.find( node ) ) cs.shared.fillIn( node, ext.second, drxn, true );
            if ( !cs.shared.empty() ) ext.first->shared_.push_back( cs );
        }
    }
}

void ClaimNode::setKeep( Node* node, Nodes& keep, bool drxn )
{
    int32_t* off[2]{ get( node ), NULL };
    if ( !off[0] || !keep.add( node ) ) return;
    for ( Edge& e : node->edges_[drxn] ) setKeep( e.node, keep, drxn );
    for ( Node* clone : node->clones( false ) ) if ( ( off[1] = get( clone ) ) && off[0] == off[1] ) setKeep( clone, keep, drxn );
}

void ClaimNode::setPairs( ClaimNode* cn, int32_t diff, bool lFork, bool rFork )
{
    if ( !lFork ) for ( ClaimEdge& ce : cn->edges_[1] ) if ( ce.node->isForked( 0 ) ) lFork = true;
    if ( !joins_[1].empty() ) lFork = true;
    
    bool added = false;
    for ( ClaimPairing* cp : pairs_[1] ) if ( !added && ( added = ( cp->node[1] == cn ) ) )
    {
        for ( int i = 0; i <= cp->diffs.size(); i++ )
        {
            if ( i == cp->diffs.size() || diff < cp->diffs[i] )
            {
                cp->diffs.insert( cp->diffs.begin() + i, diff );
                cp->missed.insert( cp->missed.begin() + i, 0 );
            }
            if ( cp->diffs[i] == diff ) break;
        }
    }
    if ( !added && lFork && rFork ) new ClaimPairing( this, cn, diff );
    
    if ( cn->edges_[1].size() > 1 ) rFork = true;
    for ( ClaimEdge& ce : cn->edges_[1] ) for ( ClaimEdge& re : ce.node->edges_[0] ) if ( re.node->isForked( 1 ) ) rFork = true;
    if ( !joins_[1].empty() ) rFork = true;
    
    for ( ClaimEdge& ce : cn->edges_[1] ) setPairs( ce.node, diff+ce.diff, lFork, lFork && rFork );
    for ( ClaimJoin* cj : cn->joins_[1] ) setPairs( cj->node[1], diff+cj->diff( 1 ), lFork, lFork && rFork );
}

bool ClaimNode::setRedundant( vector<ClaimNode*>& claims, ClaimPairing* cp, Node* paired[2], int32_t cutoff, int32_t est, int hits, bool orient )
{
    if ( !paired[0]->cloned_ && !paired[1]->cloned_ ) return false;
    for ( ClaimRedundant* cr : redundant_[orient] ) if ( cr->disregard( paired ) ) return true;
    
    Nodes clones[2];
    paired[0]->hits_.setRedundant( clones[0], clones[1], paired[1], 1 );
    if ( clones[0].size() < 2 && clones[1].size() < 2 ) return false;
    if ( clones[0].size() > 1 && clones[1].size() > 1 ) return true;
    
    int32_t* off;
    vector< pair<ClaimNode*, int32_t> > offs[2];
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : claims ) if ( cn != cp->node[d] )
    {
        for ( Node* clone : clones[d].nodes ) if ( off = cn->get( clone ) ) break;
        if ( off ) offs[d].push_back( make_pair( cn, *off ) );
    }
    
    assert( offs[0].empty() || offs[1].empty() );
    if ( offs[0].empty() && offs[1].empty() ) return true;
    
    bool drxn = !offs[1].empty();
    assert( off = cp->node[!drxn]->get( paired[!drxn] ) );
    est = *off + ( drxn ? est : -est );
    vector<ClaimNode*> alts[2]{ { cp->node[0] }, { cp->node[1] } };
    
    for ( pair<ClaimNode*, int32_t>& no : offs[drxn] )
    {
        ClaimPairing* matched = NULL;
        for ( ClaimPairing* acp : cp->node[!drxn]->pairs_[drxn] ) if ( acp->node[drxn] == no.first ) matched = acp;
        if ( !matched ) return true;
        bool success = false, fail = false;
        for ( int32_t diff : matched->diffs ) if ( abs( no.second - diff - est ) >= cutoff ) assert( false );
        for ( int32_t diff : matched->diffs ) ( abs( no.second - diff - est ) < cutoff ? success : fail ) = true;
        if ( success && fail ) assert( false );
        if ( success && fail ) return true;
        if ( !fail ) alts[drxn].push_back( no.first );
    }
    
    Nodes base;
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : alts[d] ) cn->setBase( base, true );
    for ( int d : { 0, 1 } ) for ( Node* clone : clones[d].nodes ) if ( !base.find( clone ) ) return true;
    
    for ( ClaimRedundant* cr : redundant_[orient] ) if ( cr->add( alts, paired, hits ) ) return true;
    
    new ClaimRedundant( alts, paired, hits );
    
    return true;
}

void ClaimNode::setScores( vector<ClaimNode*>& claims, Nodes& base )
{
    if ( pairs_[1].empty() ) return;
    
    for ( pair<Node*, int32_t> no : offs_ ) if ( isValid( no.first, 0 ) )
    {
        for ( auto& np : no.first->hits_.pairs[1] ) if ( base.find( np.first ) && !get( np.first ) )
        {
            int32_t* off, est = np.second.estimate();
            for ( ClaimPairing* cp : pairs_[1] ) if ( off = cp->node[1]->get( np.first ) )
            {
                if ( !cp->node[1]->isValid( np.first, 1 ) ) break;
                if ( isInvalid( no.first, no.second, cp->node[1], 1 ) ) break;
                if ( cp->node[1]->isInvalid( np.first, *off, this, 0 ) ) break;
                if ( path_[0]->ends_[0] == -4858 && cp->node[1]->path_.back()->ends_[1] == -4466 )
                {
                    int x = 0;
                }
                ClaimScore cs( no.first, np.first, cp, *off - no.second, est );
                cs.score( claims, est, 200 + np.second.maxLen / 3, np.second.count );
                break;
            }
        }
    }
    
    for ( int i = 0; i < pairs_[1].size(); i++ ) if ( !pairs_[1][i]->score ) delete pairs_[1][i--];
}

void ClaimNode::setScores( vector<ClaimNode*>& claims, Nodes& base, bool orient )
{
    if ( pairs_[orient].empty() ) return;
    
    for ( pair<Node*, int32_t> no : offs_ ) if ( isValid( no.first, !orient ) )
    {
        bool bad = false;
        for ( ClaimNode* alt : alts_ ) if ( alt->get( no.first ) ) bad = true;
        if ( bad ) continue;
        for ( auto& np : no.first->hits_.pairs[orient] ) if ( base.find( np.first ) && !get( np.first ) )
        {
            int32_t* off, est = np.second.estimate();
            for ( ClaimPairing* cp : pairs_[orient] ) if ( off = cp->node[orient]->get( np.first ) )
            {
                if ( !cp->node[orient]->isValid( np.first, orient ) ) break;
                if ( isInvalid( no.first, no.second, cp->node[orient], orient ) ) break;
                if ( cp->node[orient]->isInvalid( np.first, *off, this, !orient ) ) break;
                if ( isShared( no.first, cp->node[orient], orient ) ) break;
                if ( cp->node[orient]->isShared( np.first, this, !orient ) ) break;
                if ( distal_.find( no.first ) && cp->node[orient]->distal_.find( np.first ) ) break;
                for ( ClaimNode* alt : cp->node[orient]->alts_ ) if ( alt->get( np.first ) ) bad = true;
                if ( bad ) break;
                
                int32_t dist = orient ? *off - no.second : no.second - *off;
                int32_t cutoff = 200 + np.second.maxLen / 3;
                int32_t best = -1;
                vector<int32_t> diffs;
                for ( int32_t diff : cp->diffs ) diffs.push_back( abs( abs( dist - diff ) - est ) );
                for ( int32_t diff : diffs ) if ( best == -1 || diff < best ) best = diff;
                bool good = dist > 0;
                for ( int32_t diff : cp->diffs ) if ( dist - diff > 0 ) good = true;
                assert( good && !diffs.empty() );
                if ( cutoff < best )
                {
                    break;
                }
                cutoff += best;
                Node* paired[2]{ orient ? no.first : np.first, orient ? np.first : no.first };
                if ( setRedundant( claims, cp, paired, cutoff, est, np.second.count, orient ) ) break;
                cp->score += np.second.count;
                for ( int i = 0; i < diffs.size(); i++ ) if ( cutoff < diffs[i] ) cp->missed[i] += np.second.count;
                for ( int i = 0; i < diffs.size(); i++ ) if ( cutoff < diffs[i] ) assert( false );
                break;
            }
        }
    }
    
    for ( int i = 0; i < pairs_[orient].size(); i++ ) if ( !pairs_[orient][i]->score ) delete pairs_[orient][i--];
}

bool ClaimNode::setState( bool drxn )
{
    for ( int i = 0; i < path_.size(); i++ )
    {
        if ( ( drxn ? path_[i] : path_.end()[-i-1] )->setState() ) return true;
    }
    return false;
}

bool ClaimNode::split( vector<ClaimNode*>& claims, bool orient, bool drxn )
{
    if ( !edged_[!drxn] )
    {
        for ( Node* fork : getForks( NULL, !drxn ) ) if ( !fork->isForkComplete( params.shortLen(), 20, !drxn ) )
        {
            return false;
        }
        edge( claims, orient, !drxn );
        for ( ClaimEdge& ce : edges_[!drxn] )
        {
            for ( Node* fork : ce.node->getForks( NULL, drxn ) ) if ( !fork->isForkComplete( params.shortLen(), 20, drxn ) )
            {
                return false;
            }
            ce.node->edge( claims, orient, drxn );
        }
    }
    
    if ( !edged_[drxn] )
    {
        assert( edges_[drxn].empty() );
        edge( claims, orient, drxn );
        for ( ClaimEdge& ce : edges_[drxn] ) ce.node->edge( claims, orient, !drxn );
        return edges_[drxn].size() > 1;
    }
    
    if ( edges_[drxn].size() == 1 ) return edges_[drxn][0].node->split( claims, orient, drxn );
    
    return false;
}

void ClaimNode::testLoop( ClaimNode* cn )
{
    assert( cn != this );
    if ( !cn ) cn = this;
    for ( ClaimEdge& ce : edges_[1] ) ce.node->testLoop( cn );
    for ( ClaimJoin* cj : joins_[1] ) cj->node[1]->testLoop( cn );
}

bool ClaimNode::trim( int drxn )
{
    for ( int d : { 0, 1 } ) if ( drxn == 2 || d == drxn )
    {
        for ( ClaimEdge& ce : edges_[d] ) if ( ce.node->trim( d ) ) return true;
    }
    
    for ( int d : { 0, 1 } ) if ( edges_[d].size() > 1 && ClaimTrim::trim( this, d ) ) return true;
//    for ( int d : { 0, 1 } ) if ( edges_[d].size() > 1 )
//    {
//        bool good = false;
//        for ( ClaimEdge& ce : edges_[d] ) if ( ce.node->edges_[!d].size() > 1 ) good = true;
//        if ( good && ClaimTrim::trim( this, d ) ) return true;
//    }
    
    return false;
}

Claim::Claim( Node* fork, int32_t dist, bool orient, bool drxn )
: path{ fork }, coord( dist )
{
    for ( Edge& e : fork->edges_[drxn] ) if ( !e.node->isBlunt( 0, 3, drxn ) ) edges.push_back( e );
    offs.add( fork, dist );
    advance( fork, orient, drxn );
}

Claim::Claim( Claim* fork, Edge& e, bool orient, bool drxn )
: path( fork->path ), edges{ e }, offs( fork->offs ), coord( fork->coord )
{
    advance( terminus( drxn ), orient, drxn );
}

void Claim::advance( Node* fork, bool orient, bool drxn )
{
    while ( fork && edges.size() == 1 )
    {
        path.insert( drxn ? path.end() : path.begin(), edges[0].node );
        offs.add( fork, edges[0], coord, orient, drxn );
        fork = edges[0].node;
        coord = ( *offs.get( fork ) )[0];
        edges.clear();
        for ( Edge& e : fork->edges_[drxn] ) if ( !e.node->isBlunt( 0, 3, drxn ) ) edges.push_back( e );
        if ( !fork->cloned_ ) continue;
        for ( Node* clone : fork->cloned_->nodes )
        {
            if ( !fork->edges_[drxn].empty() && !clone->edges_[!drxn].empty() ) continue;
            offs.add( clone, coord );
            for ( Edge& e : clone->edges_[drxn] )
            {
                bool good = !e.node->isBlunt( 0, 3, drxn );
                for ( Edge& f : edges ) if ( e.node == f.node ) good = false;
                if ( good ) edges.push_back( e );
            }
        }
    }
}

void Claim::block( vector<Claim*>& claims, bool drxn )
{
    Nodes base, blocked;
    for ( Claim* c : claims ) for ( Node* node : c->path ) base += node;
    for ( auto& no : claims[0]->offs.map )
    {
        bool bad = !base.find( no.first );
        for ( int i = 1; bad && i < claims.size(); i++ )
        {
            NodeOffset* off = claims[i]->offs.get( no.first );
            if ( !off || abs( no.second[0] - (*off)[0] ) > 500 ) bad = false;
        }
        if ( bad ) blocked += no.first;
    }
    for ( Claim* c : claims )
    {
        Nodes keep;
        c->block( c->terminus( !drxn ), blocked, keep, drxn );
        c->cull( keep );
    }
}

void Claim::block( Node* node, Nodes& blocked, Nodes& keep, bool drxn )
{
    if ( !offs.find( node ) || blocked.find( node ) || !keep.add( node ) ) return;
    for ( Edge& e : node->edges_[drxn] ) block( e.node, blocked, keep, drxn );
    if ( !node->cloned_ || node->edges_[!drxn].empty() ) return;
    for ( Node* clone : node->cloned_->nodes )
    {
        if ( node->edges_[drxn].empty() || clone->edges_[!drxn].empty() ) block( clone, blocked, keep, drxn );
    }
}

void Claim::cull( Nodes& keep )
{
    for ( auto it = offs.map.begin(); it != offs.map.end(); )
    {
        if ( keep.find( it->first ) ) it++;
        else it = offs.map.erase( it );
    }
}

//int32_t Claim::estimate( Claim* f, int32_t& best, int32_t est, bool drxn )
//{
//    if ( !valid || !f->valid ) return -1;
//    int32_t diff = off.diff( f->off, est, drxn );
//    best = min( best, diff );
//    return diff;
//}

bool Claim::identical( Claim* rhs, bool drxn )
{
    if ( path.size() < rhs->path.size() ) return false;
    for ( int i = 0; i < rhs->path.size(); i++ ) if ( place( i, drxn ) != rhs->place( i, drxn ) ) return false;
    return true;
}

Node* Claim::place( int i, bool drxn )
{
    assert( i < path.size() );
    if ( i >= path.size() ) return NULL;
    return drxn ? path[i] : path.end()[-i-1];
}

bool Claim::redundant( Claim* rhs, bool drxn, bool reciprocate )
{
    Node* fork = rhs->terminus( !drxn );
    if ( fork == terminus( !drxn ) ) return false;
    NodeOffset* no[2]{ offs.get( fork ), rhs->offs.get( fork ) };
    if ( no[0] && no[1] && abs( (*no[0])[0] - (*no[1])[0] ) < 100 ) return true;
    if ( reciprocate && rhs->redundant( this, drxn, false ) ) return true;
    return false;
}

void Claim::retract()
{
    Nodes keep( path );
    for ( auto it = offs.map.begin(); it != offs.map.end(); )
    {
        if ( keep.find( it->first ) )
        {
            it->second[2] = it->second[1] = it->second[0];
            it++;
        }
        else it = offs.map.erase( it );
    }
    assert( offs.map.size() == path.size() );
}

vector<Claim*> Claim::split( bool init, bool orient, bool drxn )
{
    vector<Claim*> branches;
    int edgeCount = 0;
    for ( Edge& e : edges ) if ( e.node->isBranchComplete( 0, 15, drxn ) ) edgeCount++;
    if ( edgeCount < 2 ) return branches;
    if ( !init ) for ( Edge& e : edges ) if ( !offs.find( e.node ) ) return branches;
    
    retract();
    while ( edges.size() > 1 )
    {
        branches.push_back( new Claim( this, edges.back(), orient, drxn ) );
        edges.pop_back();
    }
    advance( terminus( drxn ), orient, drxn );
    for ( Claim* c : branches ) assert( !c->identical( this, drxn ) );
    return branches;
}

Nodes Claim::target( vector<Claim*>& branches )
{
    Nodes tar;
    for ( Claim* cb : branches ) for ( auto& no : cb->offs.map ) tar += no.first;
    for ( auto& no : branches[0]->offs.map )
    {
        bool bad = true;
        for ( int j = 1; bad && j < branches.size(); j++ ) if ( !branches[j]->offs.get( no.first ) ) bad = false;
        if ( bad ) tar -= no.first;
    }
    return tar;
}

Node* Claim::terminus( bool drxn )
{
    assert( !path.empty() );
    if ( path.empty() ) return NULL;
    return drxn ? path.back() : path[0];
}

ClaimPair::ClaimPair( Claim* fwd, Edge e, int altCount )
: path( fwd ), edge( e ), alts( altCount, 0 ), uniques( 0 ), drop( false )
{}

ClaimBlocks::ClaimBlocks( vector<Claim*> test[2] )
{
    Nodes used[3], inter[3];
    for ( int d : { 0, 1 } ) for ( Claim* c : test[d] ) for ( auto& no : c->offs.map ) used[d] += no.first;
    for ( int d : { 0, 1 } ) used[2] += used[d];
    
    NodeRoll branches[2], forks[2];
    for ( int d : { 0, 1 } ) for ( Node* node : used[d].nodes ) for ( Edge& e : node->edges_[!d] ) if ( !used[2].find( e.node ) ) branches[d].add( e.node );
    for ( Node* node : branches[0].nodes ) inter[0].fillNot( node, used[2], 1, true );
    for ( Node* node : branches[1].nodes ) if ( inter[0].find( node ) ) inter[1].fillNot( node, used[2], 0, true );
    for ( Node* node : inter[1].nodes ) if ( inter[0].find( node ) ) inter[2] += node;
    for ( int d : { 0, 1 } ) for ( int i = 0; i < branches[d].size(); i++ ) if ( !inter[2].find( branches[d][i] ) ) branches[d].remove( branches[d][i--] );
    for ( int d : { 0, 1 } ) for ( Node* node : branches[d].nodes ) for ( Edge& e : node->edges_[d] ) if ( used[2].find( e.node  ) ) forks[d].add( e.node );
    if ( inter[2].empty() ) return;
    
    bool pairs[forks[0].size()][forks[1].size()], paired[forks[0].size()]{false};
    for ( int i = 0; i < forks[0].size(); i++ )
    {
        Nodes btw = Nodes::inSet( forks[0][i], inter[2], 1, false );
        for ( int j = 0; j < forks[1].size(); j++ ) pairs[i][j] = false;
        for ( int j = 0; j < forks[1].size(); j++ ) for ( Edge& e : forks[1][j]->edges_[0] ) if ( btw.find( e.node ) ) pairs[i][j] = true;
    }
    for ( int i = 0; i < forks[0].size(); i++ )
    {
        if ( paired[i] ) continue;
        vector<int> same{ i };
        for ( int j = i+1; j < forks[0].size(); j++ )
        {
            bool diff = false;
            for ( int k = 0; k < forks[1].size(); k++ ) if ( pairs[i][k] != pairs[j][k] ) diff = true;
            if ( !diff ) same.push_back( j );
        }
        Nodes block[2];
        for ( int j : same ) block[0].fillIn( forks[0][j], used[0], 0, true );
        for ( int j = 0; j < forks[1].size(); j++ ) if ( pairs[i][j] ) block[1].fillIn( forks[1][j], used[1], 1, true );
        for ( int j : same ) paired[j] = true;
        block[0] += block[1];
        blocks.push_back( block[0] );
    }
    for ( int i = 0; i < blocks.size(); i++ )
    {
        for ( int j = 0; j < blocks.size(); j++ )
        {
            if ( i == j || blocks[i].size() < blocks[j].size() ) continue;
            bool bad = true;
            for ( Node* node : blocks[j].nodes ) if ( !( bad = blocks[i].find( node ) ) ) break;
            if ( bad && j < i ) i--;
            if ( bad ) blocks.erase( blocks.begin() + j-- );
        }
    }
}

bool ClaimBlocks::ambiguous( vector<Claim*> test[2], Node* bck, Node* fwd, bool drxn )
{
    Nodes clones[2];
    bck->hits_.setRedundant( clones[0], clones[1], fwd, drxn );
    if ( clones[0].size() < 2 && clones[1].size() < 2 ) return false;
    for ( int d : { 0, 1 } )
    {
        if ( clones[d].size() < 2 ) continue;
        for ( Claim* c : test[drxn ? d : !d] )
        {
            bool found[2]{ false, false };
            for ( Node* node : clones[d].nodes )
            {
                if ( !node->edges_[0].empty() && !node->edges_[1].empty() ) found[ c->offs.find( node ) ] = true;
            }
            if ( found[0] && found[1] ) return true;
        }
    }
    
    for ( Nodes& block : blocks )
    {
        bool bad[2]{ false, false };
        for ( int d : { 0, 1 } ) for ( Node* node : clones[d].nodes ) if ( block.find( node ) ) bad[d] = true;
        if ( bad[0] && bad[1] ) return true;
    }
    
    used.push_back( make_pair( clones[0], clones[1] ) );
    return false;
}

bool ClaimBlocks::disregard( vector<Claim*> test[2], Node* bck, Node* fwd, bool drxn )
{
    for ( pair<Nodes, Nodes>& block : used ) if ( block.first.find( bck ) && block.second.find( fwd ) ) return true;
    
    if ( ( bck->cloned_ || fwd->cloned_ ) && ambiguous( test, bck, fwd, drxn ) ) return true;
    
    for ( Nodes& block : blocks )
    {
        bool bad[2]{ block.find( bck ), block.find( fwd ) };
        if ( bck->cloned_ ) for ( Node* node : bck->cloned_->nodes ) if ( block.find( bck ) ) bad[0] = true;
        if ( fwd->cloned_ ) for ( Node* node : fwd->cloned_->nodes ) if ( block.find( fwd ) ) bad[1] = true;
        if ( bad[0] && bad[1] ) return true;
    }
    return false;
}

ClaimScores::ClaimScores( Claim* path, vector<Claim*> test[2], bool drxn )
: path( path ), claimed( false ), dropped( false )
{
    for ( Claim* c : test[drxn] )
    {
        Edge edge( NULL, 0, false );
        Nodes branches;
        for ( Node* branch : c->terminus( !drxn )->getAltForks( !drxn ) ) branches += branch;
        
        for ( Node* fork : path->terminus( drxn )->getAltForks( drxn ) )
        {
            for ( Edge& e : fork->edges_[drxn] ) if ( branches.find( e.node ) && edge.node != c->terminus( !drxn ) ) edge = e;
        }
        
        if ( edge.node ) pairs.push_back( ClaimPair( c, edge, test[!drxn].size()-1 ) );
    }
    
    for ( ClaimPair& cp : pairs ) cp.pairs = vector<int>( pairs.size(), 0 );
}

void ClaimScores::create( vector<Claim*> test[2], vector<ClaimScores*>& scores, bool orient, bool drxn )
{
    for ( ClaimScores* cs : scores ) delete cs;
    scores.clear();
    for ( Claim* c : test[!drxn] ) scores.push_back( new ClaimScores( c, test, drxn ) );
    for ( ClaimScores* a : scores ) for ( ClaimScores* b : scores ) if ( b != a ) a->alts.push_back( b );
    
    ClaimScores::score( test, scores, orient, drxn );
    for ( ClaimScores* cs : scores ) cs->trim( drxn );
}

bool ClaimScores::disregard( Node* bck, Node* fwd, Nodes& used, bool drxn )
{
    if ( !bck->cloned_ ) return false;
    for ( Node* b : bck->cloned_->nodes )
    {
        if ( !used.find( b ) ) continue;
        for ( Node* f : fwd->clones() ) if ( b->hits_.get( f, drxn ) ) return true;
    }
    return false;
}

int ClaimScores::get( ClaimScores* alt, Claim* branch )
{
    for ( ClaimPair& cp : pairs ) if ( cp.path == branch ) for ( int i = 0; i < alts.size(); i++ ) if ( alts[i] == alt ) return cp.alts[i];
    return 0;
}

bool ClaimScores::paired( Claim* c, bool exists )
{
    for ( ClaimPair& cp : pairs ) if ( cp.path == c ) return exists ? true : cp.work;
    return false;
}

bool ClaimScores::preferred( int hits[2][3], int ols[2], int altHits[2][3], int altOl, bool dual )
{
    // Filter out contended paths
    if ( hits[0][0] || !hits[1][0] || !altHits[1][0] ) return false;
    if ( min( hits[1][0] + hits[1][1], altHits[1][0] + altHits[1][2] ) - min( hits[0][1], hits[0][2] ) < 4 )return false;
//    if ( min( hits[0], hits[1] ) * 5 >= min( hits[2], altHits[3] )-1 ) return false;
//    if ( hits[0] * 5 >= min( hits[1], altHits[1] )-1 ) return false;
    
    if ( altHits[0][0] ) return min( hits[1][0], altHits[1][0] ) - max( hits[0][1], hits[0][2] ) > 8;
    
    // Assign bonus score for preferential overlap if different forks
    int diff = dual ? min( 100, min( 60, ols[1] - ols[0] ) + min( 60, altOl - ols[0] ) ) / 10 : 0;
    int hit = hits[1][1] + altHits[1][2] + hits[1][0] + altHits[1][0];
//    int hit = hits[2] + altHits[3] + diff;
//    if ( hits[0] || hits[1] || altHits[0] || altHits[1] ) hit = min( hits[2], altHits[3] );
//    int hit = hits[1] + altHits[1] + diff;
//    if ( max( hits[0], altHits[0] ) * 5 >= min( hits[1], altHits[1] )-1 ) hit = min( hits[1], altHits[1] );
    
    return hit + diff - min( hits[0][1], hits[0][2] ) > 8;
}

bool ClaimScores::remove( Claim* branch )
{
    for ( int i = 0; i < pairs.size(); i++ )
    {
        if ( pairs[i].path != branch ) continue;
        pairs.erase( pairs.begin() + i );
        dropped = true;
        for ( ClaimPair& cp : pairs ) cp.pairs.erase( cp.pairs.begin() + i );
        return true;
    }
    return false;
}

void ClaimScores::score( vector<Claim*> test[2], vector<ClaimScores*>& scores, bool orient, bool drxn )
{
    ClaimBlocks cb( test );
    Nodes tars[2]{ Claim::target( test[0] ), Claim::target( test[1] ) };
    int d = tars[0].nodes.size() < tars[1].nodes.size();
    
    for ( Node* b : tars[!d].nodes )
    {
        for ( auto& np : b->hits_.pairs[d] )
        {
            if ( !tars[d].find( np.first ) ) continue;
            if ( cb.disregard( test, b, np.first, d ) ) continue;
            
            // Prepare
            Node* bck = d == drxn ? b : np.first,* fwd = d == drxn ? np.first : b;
            int32_t est = np.second.estimate() + ( d == orient ? 0 : b->size() - np.first->size() ); 
            int32_t cutoff = 200 + np.second.maxLen / 3;
            int32_t best = cutoff * 2;
            
            // Set distance estimates for each branch pair
            for ( ClaimScores* cs : scores ) cs->score( bck, fwd, est, best, drxn );
            
            // Abort if estimates look dubious
            assert( best >= 0 );
            if ( best > cutoff ) continue;
            
            // Convert estimates to boolean if better than cutoff
            bool unique = score( scores, cutoff+best, drxn );
            
            // Record where one back branch matches, but another doesn't
            for ( ClaimScores* cs : scores ) cs->score( np.second.count, unique );
        }
    }
}

bool ClaimScores::score( vector<ClaimScores*>& scores, int32_t cutoff, bool drxn )
{
    vector<Claim*> claims[2];
    for ( ClaimScores* cs : scores )
    {
        for ( ClaimPair& cp : cs->pairs )
        {
            if ( !( cp.work = ( cp.work >= 0 && cp.work < cutoff ) ) ) continue;
            claims[!drxn].push_back( cs->path );
            claims[drxn].push_back( cp.path );
        }
    }
    
    for ( int d : { 0, 1 } ) if ( claims[d].size() == 1 && claims[!d].size() != 1 ) return false;
    for ( int d : { 0, 1 } ) for ( int i = 0; i+1 < claims[d].size(); i++ ) for ( int j = i+1; j < claims[d].size(); j++ )
    {
        if ( claims[d][i] == claims[d][j] ) return false;
    }
    for ( int d : { 0, 1 } ) for ( int i = 0; i+1 < claims[d].size(); i++ ) for ( int j = i+1; j < claims[d].size(); j++ )
    {
        if ( !claims[d][i]->redundant( claims[d][j], d ) ) return false;
    }
    return true;
}

void ClaimScores::score( Node* bck, Node* fwd, int32_t est, int32_t& best, bool drxn )
{
    for ( ClaimPair& cp : pairs ) cp.work = -1;
    NodeOffset* offs[2];
    for ( Node* b : bck->clones() )
    {
        if ( !( offs[0] = path->offs.get( b ) ) ) continue;
        for ( ClaimPair& cp : pairs )
        {
            for ( Node* f : fwd->clones() )
            {
                if ( !( offs[1] = cp.path->offs.get( f ) ) ) continue;
                int32_t diff = offs[0]->diff( *offs[1], est, drxn );
                cp.work = cp.work < 0 || diff < cp.work ? diff : cp.work;
                best = min( best, diff );
            }
        }
    }
}

void ClaimScores::score( int hit, bool unique )
{
    for ( ClaimPair& cp : pairs )
    {
        if ( !cp.work ) continue;
        for ( int i = 0; i < alts.size(); i++ ) if ( !alts[i]->paired( cp.path ) ) cp.alts[i] += hit;
        for ( int i = 0; i < pairs.size(); i++ ) if ( !pairs[i].work ) cp.pairs[i] += hit;
        if ( unique ) cp.uniques += hit;
    }
}

bool ClaimScores::set( ClaimScores* alt, Claim* weak, Claim* strong, int hits[2][3], int ols[2], bool& paired )
{
    int i = 0, j = 0, k = 0;
    while ( i < pairs.size() && pairs[i].path != weak ) i++;
    while ( j < pairs.size() && pairs[j].path != strong ) j++;
    while ( k < alts.size() && alts[k] != alt ) k++;
    assert( k < alts.size() );
    if ( j == pairs.size() || !pairs[j].uniques ) return false;
    paired = i < pairs.size();
    hits[0][0] = paired ? pairs[i].uniques : 0;
    hits[0][1] = paired ? pairs[i].pairs[j] : 0;
    hits[0][2] = paired ? pairs[i].alts[k] : 0;
    hits[1][0] = pairs[j].uniques;
    hits[1][1] = paired ? pairs[j].pairs[i] : 0;
    hits[1][2] = pairs[j].alts[k];
    ols[0] = paired ? pairs[i].edge.ol : 0;
    ols[1] = pairs[j].edge.ol;
    return true;
}

void ClaimScores::trim( bool drxn )
{
    int hits[2][3], altHits[2][3];
    for ( int i = 0; i < pairs.size(); i++ )
    {
        // Only path pairs that share no unique read pairs are able to be trimmed
        if ( hits[0][0] = pairs[i].uniques ) continue;
        
        for ( int j = 0; !pairs[i].drop && j < pairs.size(); j++ )
        {
            // Stronger branch must have at least some unique read pairs
            if ( i == j || !( hits[1][0] = pairs[j].uniques ) ) continue;
            // Stronger branch must be favoured over the weaker branch by this path
            if ( ( hits[0][1] = pairs[i].pairs[j] ) >= ( hits[1][1] = pairs[j].pairs[i] ) ) continue;
            
            for ( int k = 0; !pairs[i].drop && k < alts.size(); k++ )
            {
//                // Ensure that either this path is redundant with the alt, or pair [i] is redundant with pair [j]
//                if ( pairs[i].alts[j] && pairs[i].alts[k] ) continue;
                // Stronger branch must be favoured over the weaker branch by this path
                if ( ( hits[0][2] = pairs[i].alts[k] ) >= ( hits[1][2] = pairs[j].alts[k] ) ) continue;
                
                bool dual = path->terminus( drxn ) != alts[k]->path->terminus( drxn ), paired;
                int ols[2][2]{ { pairs[i].edge.ol, pairs[j].edge.ol }, { 0, 0 } };
                
                if ( !alts[k]->set( this, pairs[j].path, pairs[i].path, altHits, ols[1], paired ) ) continue;
                
                // Ensure that the target is weaker than the preferred pair and that the alt is strongly paired to the target
                if ( !preferred( hits, ols[0], altHits, ols[1][1], dual ) ) continue;
                
                cout << hits[0][0] << "-" << max( hits[1][1], hits[0][2] ) << " & " 
                        << altHits[0][0] << "-" << max( altHits[1][1], altHits[1][2] ) << "    ";
                cout << path->terminus( drxn )->ends_[drxn] << endl;
                
                if ( !paired ) alts[k]->dropped = true;
                else if ( preferred( altHits, ols[1], hits, ols[0][1], dual ) )
                {
                    for ( ClaimPair& cp : alts[k]->pairs ) if ( cp.path == pairs[j].path ) alts[k]->dropped = cp.drop = true;
                }
                
                pairs[i].drop = dropped = true;
            }
        }
    }
}

ClaimMap::ClaimMap( ClaimScores* cs, Nodes& bases, NodeRoll& cloned, bool drxn )
: bases( bases ), cloned( cloned ), fork( cs->path->terminus( drxn ) ), claimed( false )
{
    if ( cs->claimed || !cs->dropped ) return;
    scores[0].push_back( cs );
    claims[!drxn].insert( cs->path );
    for ( ClaimPair& cp : cs->pairs ) claims[drxn].insert( cp.path );
    
    // Search alts for identical pairing, or else nodes that will need to be cloned
    for ( ClaimScores* alt : cs->alts ) scores[alt->claimed || !match( alt, drxn )].push_back( alt );
    
    assert( confirm( drxn ) );
    for ( int d : { 0, 1 } ) for ( Claim* c : claims[d] ) claim( c, d );
    
    Nodes branches;
    for ( Claim* c : claims[drxn] ) branches += get( c->terminus( !drxn ) );
    for ( Edge& e : ( fork = get( fork ) )->edges( drxn ) )
    {
        bool discard = bases.find( e.node ) && !branches.find( e.node );
        if ( e.node->cloned_ ) for ( Node* node : safe[drxn].nodes ) if ( e.node->cloned_->find( node ) ) discard = false;
        if ( discard ) assert( claimed = fork->removeEdge( e.node, drxn, true ) );
    }
    
    for ( ClaimPair& cp : cs->pairs ) disconnect( cp, drxn );
    
    for ( ClaimPair& cp : cs->pairs ) connect( cp, drxn );
    
    if ( claimed ) for ( ClaimScores* cs : scores[0] ) cs->claimed = true;
}

void ClaimMap::claim( Claim* c, bool drxn )
{
    Node* node[3] = { NULL, NULL, NULL };
    for ( int i = 0; i < c->path.size(); i++ )
    {
        node[0] = node[1] = c->place( i, drxn );
        node[1] = get( node[0] );
        if ( node[2] ) node[1]->addEdge( node[0]->getEdge( node[2], !drxn, true, true ), !drxn, true );
        if ( node[0] == node[1] ) break;
        for ( int d : { 0, 1 } ) for ( Edge& e : node[0]->edges_[!d] ) if ( !bases.find( e.node ) ) node[1]->addEdge( e, !d, true );
        
        node[2] = node[1];
    }
    assert( node[0] == node[1] );
    
    if ( !node[2] ) return;
    
    // If any of the path was cloned, remove obsolete edge to the last cloned, or else invalid edge to different pair
    for ( Edge e : node[0]->edges( !drxn ) )
    {
        // Only mistmatched edges
        if ( !alts[0].find( e.node ) ) continue;
        // Don't remove fork in shared path nodes
        if ( alts[1].find( node[0] ) && alts[1].find( e.node ) ) continue;
        
        assert( claimed = node[0]->removeEdge( e.node, !drxn, true ) );
    }
}

bool ClaimMap::confirm( bool drxn )
{
    Nodes branches[2], breaks[2];
    bool broken = false, looped = false;
    for ( int d : { 0, 1 } )
    {
        for ( Claim* c : claims[d] )
        {
            bool branched = false;
            for ( int i = 0; !branched && i < c->path.size(); i++ )
            {
                Node* node = c->place( i, d );
                if ( branched = !alts[0].find( node ) ) branches[d] += node;
                if ( node->cloned_ ) assert( !node->edges_[d].empty() );
                if ( !node->cloned_ || !node->edges_[!d].empty() ) continue;
                if ( !branched ) broken = true;
                else if ( !i ) breaks[d] += node;
            }
        }
    }
    
    if ( !broken && breaks[0].empty() && breaks[1].empty() ) return true;
    Nodes loop;
    for ( Node* node : branches[1].nodes ) loop.fill( node, 1, false, false );
    for ( Node* node : branches[0].nodes ) if ( loop.find( node ) ) looped = true;
    if ( looped && broken ) assert( false );
    if ( looped && broken ) return false;
    if ( looped ) for ( int d : { 0, 1 } ) for ( Node* node : breaks[d].nodes ) safe[d].add( node );
    return true;
}

void ClaimMap::connect( ClaimPair& cp, bool drxn )
{
    if ( safe[!drxn].find( fork ) || safe[drxn].find( cp.path->terminus( !drxn ) ) ) return;
    cp.edge.node = get( cp.path->terminus( !drxn ) );
    fork->addEdge( cp.edge, drxn, true );
    for ( ClaimScores* alt : scores[1] )
    {
        if ( alt->claimed ) continue;
        for ( ClaimPair& acp : alt->pairs )
        {
            if ( acp.path != cp.path ) continue;
            acp.edge.node = cp.edge.node;
            alt->path->terminus( drxn )->addEdge( acp.edge, drxn, true );
        }
    }
}

void ClaimMap::disconnect( ClaimPair& cp, bool drxn )
{
    Node* branch = cp.path->terminus( !drxn );
    if ( branch != get( branch ) ) return;
    Nodes forks( fork );
    for ( ClaimScores* alt : scores[1] ) if ( alt->paired( cp.path, true ) ) forks += alt->path->terminus( drxn );
    for ( Edge& e : branch->edges( !drxn ) )
    {
        bool discard = bases.find( e.node ) && !forks.find( e.node );
        if ( e.node->cloned_ ) for ( Node* node : safe[!drxn].nodes ) if ( e.node->cloned_->find( node ) ) discard = false;
        if ( discard ) assert( claimed = branch->removeEdge( e.node, !drxn, true ) );
    }
}

Node* ClaimMap::get( Node* node )
{
    if ( !alts[0].find( node ) ) return node;
    auto it = used.find( node );
    if ( it != used.end() ) return it->second;
    Node* clone = new Node( node, cloned, node->drxn_, false );
    used.insert( make_pair( node, clone ) );
    bases += clone;
    claimed = true;
    return clone;
}

bool ClaimMap::match( ClaimScores* cs, bool drxn )
{
    vector<Claim*> matched[2];
    for ( ClaimPair& cp : cs->pairs ) matched[ claims[drxn].find( cp.path ) != claims[drxn].end() ].push_back( cp.path );
    
    for ( Claim* c : matched[0] ) unmatch( c, drxn );
    
    if ( cs->path->terminus( drxn ) != fork ) return false;
    
    if ( matched[0].empty() && matched[1].size() == claims[drxn].size() )
    {
        claims[!drxn].insert( cs->path );
        return true;
    }
    
    for ( Claim* c : matched[1] ) for ( Node* node : c->path ) alts[1] += node;
    unmatch( cs->path, !drxn );
    return false;
}

void ClaimMap::unmatch( Claim* path, bool drxn )
{
    int i = 0;
    for ( Claim* c : claims[drxn] )
    {
        for ( int j = 0; j < min( path->path.size(), c->path.size() ); j++ )
        {
            if ( path->place( j, drxn ) != c->place( j, drxn ) ) break;
            i = max( i, j+1 );
        }
    }
    if ( i )
    {
        int x = 0;
    }
    assert( i < path->path.size() );
    for ( int j = 0; j < i; j++ ) alts[0] += path->place( j, drxn );
}

ClaimFork::ClaimFork( Node* fork, NodeRoll& nodes, bool drxn )
: fork_( fork ), orient_( drxn ), drxn_( drxn ), claimed_( false )
{
    if ( fork->cloned_ && fork->edges_[!drxn].empty() )
    {
        int x = 0;
    }
    
    bool looped = fork->cloned_ && fork->edges_[!drxn].empty();
    if ( fork->cloned_ ) for ( Node* clone : fork->cloned_->nodes ) if ( clone->edges_[!drxn].empty() ) looped = true;
    for ( Edge& e : fork->edges_[drxn] ) add( fork, e, 0, drxn );
    if ( looped ) for ( Node* clone : fork->cloned_->nodes ) for ( Edge& e : clone->edges_[drxn] ) add( clone, e, 0, drxn );
    
//    for ( Node* f : fork->getAltForks( drxn ) ) for ( Edge& e : f->edges_[!drxn] ) add( f, e, 0, false, !drxn );
//    for ( Node* f : fork->getAltForks( drxn ) ) forks_ += f;
//    for ( Node* f : forks_.nodes ) for ( Edge& e : f->edges_[!drxn] ) add( f, e, 0, false, !drxn );
    for ( int d : { 0, 1 } ) if ( !complete( d ) ) return;
    claim( nodes );
}

ClaimFork::ClaimFork( Node* fork, Node* branch, bool drxn )
: fork_( fork ), orient_( drxn ), drxn_( drxn ), claimed_( false )
{
    Claim* c = new Claim( fork, 0, orient_, !drxn );
    paths_[!drxn].push_back( c );
    for ( Edge& e : fork->getAltEdges( drxn ) )
    {
        paths_[drxn].push_back( new Claim( e.node, NodeOffsets::extend( fork, e, 0, drxn, drxn ), orient_, drxn ) );
    }
    for ( int d : { 0, 1 } ) for ( Claim* c : paths_[d] ) added_.push_back( c );
    
    c->offs.fill( c->terminus( !drxn ), c->coord, 500, drxn, !drxn, false, true );
    complete( drxn );
    ClaimScores::create( paths_, scores_, drxn, drxn );
    int x = 0;
}

ClaimFork::~ClaimFork()
{
    for ( Claim* c : added_ ) delete c;
    for ( ClaimScores* cs : scores_ ) delete cs;
}

void ClaimFork::add( Node* fork, Edge edge, int32_t dist, bool drxn )
{
    if ( edge.node->isBlunt( 0, 3, drxn ) ) return;
    for ( Claim* c : paths_[drxn] ) if ( edge.node == c->terminus( !drxn ) ) return;
    if ( edge.node->cloned_ )
    {
        for ( Node* clone : edge.node->cloned_->nodes )
        {
            if ( clone->edges_[drxn].empty() ) continue;
            if ( edge.node->edges_[drxn].empty() || clone->edges_[!drxn].empty() ) add( fork, Edge( clone, edge.ol, edge.isLeap ), dist, drxn );
        }
    }
    
    dist = NodeOffsets::extend( fork, edge, dist, orient_, drxn );
    
    if ( !edge.node->cloned_ || !edge.node->edges_[drxn].empty() )
    {
        paths_[drxn].push_back( new Claim( edge.node, dist, orient_, drxn ) );
        added_.push_back( paths_[drxn].back() );
    }
    
    for ( Edge& e : edge.node->edges_[!drxn] ) add( edge.node, e, dist, !drxn );
}

bool ClaimFork::claim( NodeRoll& nodes )
{
    if ( paths_[0].size() < 2 || paths_[1].size() < 2 ) return false;
    
    
    ClaimScores::create( paths_, scores_, orient_, drxn_ );
    if ( !confirm() ) return claim( nodes );
    
    if ( fork_->ends_[1] == -3822 )
    {
        int x = 0;
    }
    
    bool claimed = false, claimable = false;
    for ( ClaimScores* cs : scores_ ) if ( cs->dropped && fork_->isClone( cs->path->terminus( drxn_ ) ) ) claimable = true;
    if ( claimable )
    {
        Nodes bases;
        for ( int d : { 0, 1 } ) for ( Claim* c : paths_[d] ) for ( Node* node : c->path ) for ( Node* clone : node->getAltForks( !d ) ) bases += clone;
        for ( ClaimScores* cs : scores_ ) if ( ClaimMap( cs, bases, cloned_, drxn_ ).claimed ) claimed = true;
    }
    if ( !claimed )
    {
        vector<Nodes> fwds;
        for ( Claim* c : paths_[!drxn_] ) fwds.push_back( Nodes( c->terminus( !drxn_ ), !drxn_, false, false ) );
        cout << "FAILED ";
        cout << fork_->ends_[drxn_] << " ";
        cout << "paths: " << paths_[!drxn_].size() << "-" << paths_[drxn_].size() << " scores: ";
        int i = 0;
        for ( ClaimScores* cs : scores_ )
        {
            cout << "(" << ++i << ") ";
            for ( ClaimPair& cp : cs->pairs )
            {
                for ( int j = 0; j < cp.alts.size(); j++ ) cout << ( j ? "-" : "" ) << cp.alts[j];
                cout << " ";
            }
        }
        cout << endl;
        int x = 0;
        for ( ClaimScores* cs : scores_ ) if ( cs->dropped && fork_->isClone( cs->path->terminus( drxn_ ) ) ) return false;
        
        return !claimable && retry() && claim( nodes );
    }
    
    for ( Node* clone : cloned_.nodes ) if ( nodes.add( clone ) ) clone->reverify();
    
    for ( int d : { 0, 1 } ) for ( Claim* c : paths_[d] ) for ( int i = 0; i < c->path.size(); i++ ) if ( c->place( i, !d )->setState() ) break;
    
    for ( int d : { 0, 1 } )
    {
        for ( Claim* c : paths_[d] )
        {
            Node* branch = c->terminus( !d );
            if ( !branch->bad_ || branch->drxn_ != d ) continue;
            for ( Edge& e : branch->edges_[!d] ) assert( e.node->drxn_ == !d || e.node->bad_ );
        }
    }
    
    if ( claimed ) claimed_ = true;
    return claimed_;
}

bool ClaimFork::complete( bool drxn )
{
    if ( paths_[drxn].size() == 1 )
    {
        for ( Claim* c : paths_[drxn][0]->split( true, orient_, drxn ) )
        {
            paths_[drxn].push_back( c );
            added_.push_back( c );
        }
    }
    
    if ( paths_[drxn].size() < 2 ) return false;
    
    Node* same = paths_[drxn][0]->terminus( drxn );
    for ( int i = 1; same && i < paths_[drxn].size(); i++ ) if ( paths_[drxn][i]->terminus( drxn ) != same ) same = NULL;
    if ( same ) return true;
    
    int32_t dist = 0;
    for ( Claim* c : paths_[drxn] ) dist = max( dist, abs( c->coord ) );
    for ( Claim* c : paths_[drxn] ) c->offs.fill( c->terminus( drxn ), c->coord, dist + 500, orient_, drxn, false, true );
    dist = max( dist, params.maxPeMean + params.readLen );
    Claim::block( paths_[drxn], drxn );
    for ( Claim* c : paths_[drxn] )
    {
        Nodes keep, alts[2];
        for ( Node* node : c->path ) keep += node;
        for ( Claim* alt : paths_[drxn] ) for ( auto& no : alt->offs.map ) alts[ c != alt ] += no.first;
        for ( auto& no : c->offs.map ) if ( abs( no.second[0] ) <= dist || alts[1].find( no.first ) ) keep.fillIn( no.first, alts[0], !drxn, true );
        c->cull( keep );
    }
    
    return true;
}

bool ClaimFork::confirm()
{
    unordered_set<Claim*> drop[2], unique[2];
    for ( ClaimScores* cs: scores_ )
    {
        if ( cs->dropped ) drop[!drxn_].insert( cs->path );
        for ( ClaimPair& cp : cs->pairs )
        {
            if ( cp.drop ) drop[drxn_].insert( cp.path );
            else if ( cp.uniques < 5 ) continue;
            unique[!drxn_].insert( cs->path );
            unique[drxn_].insert( cp.path );
        }
    }
    
    for ( int d : { 0, 1 } )
    {
        for ( Claim* c : paths_[d] )
        {
            if ( drop[d].find( c ) != drop[d].end() || unique[d].find( c ) != unique[d].end() ) continue;
            Node* fork = c->terminus( !d );
            vector<Claim*> contested;
            for ( Claim* c : drop[d] ) if ( fork == c->terminus( !d ) ) contested.push_back( c );
            if ( contested.empty() ) continue;
            if ( d ) assert( false );
            contest( c, contested, d );
            return false;
        }
    }
    
    for ( ClaimScores* cs: scores_ )
    {
        for ( int i = 0; i < cs->pairs.size(); i++ ) if ( cs->pairs[i].drop ) cs->pairs.erase( cs->pairs.begin() + i-- );
    }
    
    return true;
}

void ClaimFork::contest( Claim* branch, vector<Claim*>& contested, bool drxn )
{
    int j = branch->path.size();
    for ( int i = 0; !contested.empty() && i < branch->path.size(); i++ )
    {
        Node* node = branch->place( i, drxn );
        for ( int k = 0; k < contested.size(); k++ )
        {
            if ( contested[k]->place( i, drxn ) == node ) continue;
            contested.erase( contested.begin() + k-- );
            j = min( j, i );
        }
        
        if ( !contested.empty() ) continue;
        
        if ( i != j )
        {
            ClaimFork cf( node, branch->place( i-1, drxn ), !drxn );
            if ( cf.claimed_ ) j = i;
        }
        
        if ( i == j ) paths_[drxn].erase( remove( paths_[drxn].begin(), paths_[drxn].end(), branch ), paths_[drxn].end() );
        else for ( int d : { 0, 1 } ) paths_[d].clear();
    }
    assert( contested.empty() );
}

bool ClaimFork::retry()
{
    bool split[2]{ false, false };
    unordered_set<Claim*> splitable[2];
    
    for ( ClaimScores* cs : scores_ )
    {
        if ( !fork_->isClone( cs->path->terminus( drxn_ ) ) ) continue;
        
        int uniques = 0;
        for ( ClaimPair& cp : cs->pairs )
        {
            if ( cp.uniques < 5 ) continue;
            
            if ( uniques++ ) splitable[!drxn_].insert( cs->path );
            
            for ( ClaimScores* alt : cs->alts ) for ( ClaimPair& acp : alt->pairs ) if ( acp.path == cp.path && acp.uniques > 4 ) splitable[drxn_].insert( cp.path );
        }
    }
    
    for ( bool d : { 0, 1 } )
    {
        for ( Claim* c : splitable[d] )
        {
            for ( Claim* s : c->split( false, orient_, d ) )
            {
                paths_[d].push_back( s );
                added_.push_back( s );
                split[d] = true;
            }
        }
        
        if ( !split[d] ) continue;
        for ( Claim* c : paths_[d] ) c->retract();
        complete( d );
    }
    
    if ( split[0] || split[1] )
    {
        cout << "Split!" << endl;
    }
    return split[0] || split[1];
}
