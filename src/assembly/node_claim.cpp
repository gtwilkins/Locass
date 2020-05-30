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

void ClaimJoin::fill( Node* node, NodeDists& dists, bool drxn )
{
    if ( dists.find( node ) && bridge.add( node ) ) for ( Edge& e : node->edges_[drxn] ) fill( e.node, dists, drxn );
}

ClaimRepair::ClaimRepair( Node* l, Node* r )
{
    forks[0].push_back( l );
    forks[1].push_back( r );
    for ( int d : { 0, 1 } ) for ( Node* node : forks[d] ) fill( node, -node->size(), d );
}

ClaimRepair::ClaimRepair( vector<Node*>& l, vector<Node*>& r )
{
    forks[0] = l;
    forks[1] = r;
    for ( int d : { 0, 1 } ) for ( Node* node : forks[d] ) fill( node, -node->size(), d );
}

void ClaimRepair::dupe( vector< pair<ClaimNode*, ClaimNode*> >& cloned, NodeRoll& nodes )
{
    Node* cloneForks[2]{ cloned[0].second->path_[0], cloned.back().second->path_.back() };
    ClaimRepair cr[2]{ ClaimRepair( forks[0], forks[1] ), ClaimRepair( cloneForks[0], cloneForks[1] ) };
    unordered_map<Node*, int32_t> diffs[2][2]{ { get( cr[0], 0 ), get( cr[1], 0 ) }, { get( cr[0], 1 ), get( cr[1], 1 ) } };
    
    for ( int d : { 0, 1 } ) cloneForks[d]->setState();
    nodes.updateBad();
    
    bool verified[2]{ false, false };
    for ( int d : { 0, 1 } ) for ( Edge& e : cloneForks[d]->edges_[d] ) if ( e.node->verified_ ) verified[d] = true;
    for ( pair<ClaimNode*, ClaimNode*>& clones : cloned ) if ( verified[0] && verified[1] ) for ( Node* node : clones.second->path_ ) node->verified_ = true;
    
    Nodes repairs;
    for ( pair<ClaimNode*, ClaimNode*>& clones : cloned )
    {
        for ( Node* node : clones.first->path_ ) for ( int d : { 0, 1 } ) for ( int i : { 0, 1 } ) diffs[d][i].erase( node );
        for ( Node* node : clones.second->path_ ) for ( int d : { 0, 1 } ) for ( int i : { 0, 1 } ) diffs[d][i].erase( node );
        for ( Node* node : clones.first->path_ ) node->clearPaired( true );
        for ( Node* node : clones.second->path_ ) node->clearPaired( true );
        for ( Node* node : clones.second->path_ ) repairs += node;
    }
    for ( Node* node : repairs.nodes ) node->reverify();
    
    for ( int d : { 0, 1 } ) for ( int i : { 0, 1 } ) for ( auto& no : diffs[d][i] ) sever( no.first, diffs[!d][!i], !d );
}

void ClaimRepair::fill( Node* node, int32_t off, bool drxn )
{
    auto it = offs[drxn].insert( make_pair( node, off ) );
    if ( !it.second )
    {
        if ( it.first->second <= off ) return;
        it.first->second = off;
    }
    
    if ( off < params.maxMpMean ) for ( Edge& e : node->edges_[drxn] ) fill( e.node, off + node->size() - e.ol, drxn );
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

void ClaimRepair::sever( Node* node, unordered_map<Node*, int32_t>& tar, bool drxn )
{
    Nodes used;
    vector< vector< pair<Node*, int32_t> > > repairs;
    for ( auto& np : node->hits_.pairs[drxn] ) if ( np.first->cloned_ && tar.find( np.first ) != tar.end() && !used.find( np.first ) )
    {
        vector< pair<Node*, int32_t> > ests = { make_pair( np.first, np.second.estimate() ) };
        for ( Node* clone : np.first->cloned_->nodes )
        {
            used += clone;
            auto it = node->hits_.pairs[drxn].find( clone );
            if ( it != node->hits_.pairs[drxn].end() ) ests.push_back( make_pair( clone, it->second.estimate() ) );
        }
        if ( ests.size() > 1 ) repairs.push_back( ests );
    }
    
    if ( repairs.empty() ) return;
    
    int32_t limit = 0;
    for ( vector< pair<Node*, int32_t> >& ests : repairs ) for ( pair<Node*, int32_t>& est : ests ) limit = max( limit, abs( est.second ) );
    
    NodeDists dists( node, limit+500, drxn, drxn, true );
    int32_t* dist = NULL;
    for ( vector< pair<Node*, int32_t> >& ests : repairs )
    {
        for ( pair<Node*, int32_t>& est : ests ) est.second = ( dist = dists.get( est.first ) ) ? abs( est.second - *dist ) : params.maxMpMean+1;
        sort( ests.begin(), ests.end(), []( pair<Node*, int32_t>& a, pair<Node*, int32_t>& b ){ return a.second < b.second; } );
        for ( pair<Node*, int32_t>& est : ests ) if ( est.second > min( ests[0].second*2 + 200, params.maxMpMean ) )
        {
            node->hits_.erase( est.first, drxn );
            est.first->hits_.erase( node, !drxn );
        }
    }
}

void ClaimRepair::trim()
{
    ClaimRepair cr( forks[0], forks[1] );
    unordered_map<Node*, int32_t> diffs[2]{ get( cr, 0 ), get( cr, 1 ) };
    for ( int d : { 0, 1 } ) for ( auto& no : diffs[d] ) sever( no.first, diffs[!d], !d );
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

int ClaimRedundant::get( unordered_set<ClaimNode*>& lTar, unordered_set<ClaimNode*>& rTar, ClaimNode* l, ClaimNode* r, bool lEdge, bool rEdge )
{
    for ( ClaimNode* cn : nodes[0] ) if ( lTar.find( cn ) == lTar.end() ) return 0;
    for ( ClaimNode* cn : nodes[1] ) if ( rTar.find( cn ) == rTar.end() ) return 0;
    
    vector<bool> tars( nodes[1].size(), false );
    unordered_set<ClaimNode*> pathed;
    for ( ClaimNode* cn : nodes[0] )
    {
        bool reached = false;
        if ( !reach( cn, l, r, tars, pathed, reached, lEdge, rEdge ) || !reached ) return 0;
    }
    for ( bool reached : tars ) if ( !reached ) return 0;
    return score;
}

bool ClaimRedundant::reach( ClaimNode* cn, ClaimNode* l, ClaimNode* r, vector<bool>& tars, unordered_set<ClaimNode*>& pathed, bool& reached, bool lEdge, bool rEdge )
{
    for ( int i = 0; i < nodes[1].size(); i++ ) if ( cn == nodes[1][i] )
    {
        if ( l || r ) return false;
        tars[i] = reached = true;
        return true;
    }
    
    if ( !pathed.insert( cn ).second ) return true;
    if ( !lEdge && !rEdge && cn == l ) l = NULL;
    
    for ( ClaimEdge& ce : cn->edges_[1] )
    {
        if ( ce.node == r && !lEdge && !rEdge && l ) return false;
        
        bool edged = lEdge == ( cn == l ) && rEdge == ( ce.node == r );
        if ( !lEdge && !rEdge ) edged = ce.node == r;
        
        if ( !reach( ce.node, edged ? NULL : l, edged ? NULL : r, tars, pathed, reached, lEdge, rEdge ) ) return false;
    }
    
    for ( ClaimJoin* cj : cn->joins_[1] ) if ( !reach( cj->node[1], l, !l && cj->node[1] == r ? NULL : r, tars, pathed, reached, lEdge, rEdge ) ) return false;
    pathed.erase( cn );
    return true;
}

ClaimPairing::ClaimPairing( ClaimNode* l, ClaimNode* r, int32_t dist )
: diffs{ dist }, missed{ 0 }, hits( 0 )
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
    int miss = -hits;
    for ( int32_t d : hitDiffs ) for ( int i = 0; i < diffs.size(); i++ ) if ( diffs[i] == d ) miss = max( miss, -missed[i] );
    assert( !miss );
    return miss + hits;
}

int ClaimPairing::get( ClaimNode* l, ClaimNode* r, bool lEdge, bool rEdge )
{
    unordered_set<ClaimNode*> pathed;
    int score = 0;
    return reach( node[0], l, r, pathed, 0, score, lEdge, rEdge ) ? score : 0;
}

int ClaimPairing::get( vector<ClaimNode*>& path, int mode )
{
    unordered_set<ClaimNode*> pathed;
    int score = 0;
    return reach( node[0], path, 0, mode, pathed, 0, score, false, false ) ? score : 0;
}

bool ClaimPairing::reach( ClaimNode* cn, vector<ClaimNode*>& path, int i, int mode, unordered_set<ClaimNode*>& pathed, int32_t diff, int& score, bool success, bool fail )
{
    if ( cn == path[0] ) ( mode == 0 ? fail : success ) = true;
    if ( cn == path.back() ) ( mode == 1 ? fail : success ) = true;
    if ( i < path.size() ) i = ( cn == path[i] ? i+1 : 0 );
    
    if ( cn == node[1] )
    {
        if ( mode == 2 ? ( i != path.size() ) : ( !success || fail ) ) return false;
        for ( int i = 0; i < diffs.size(); i++ ) if ( diff == diffs[i] ) score = max( score, hits - missed[i] );
        return true;
    }
    
    if ( !pathed.insert( cn ).second ) return true;
    for ( ClaimEdge& ce : cn->edges_[1] ) if ( !reach( ce.node, path, i, mode, pathed, diff+ce.diff, score, success, fail ) ) return false;
    for ( ClaimJoin* cj : cn->joins_[1] ) if ( !reach( cj->node[1], path, i, mode, pathed, diff+cj->dist, score, success, fail ) ) return false;
    pathed.erase( cn );
    
    return true;
}

bool ClaimPairing::reach( ClaimNode* cn, ClaimNode* l, ClaimNode* r, unordered_set<ClaimNode*>& pathed, int32_t diff, int& score, bool lEdge, bool rEdge )
{
    if ( cn == node[1] )
    {
        if ( l || r ) return false;
        for ( int i = 0; i < diffs.size(); i++ ) if ( diff == diffs[i] ) score = max( score, hits - missed[i] );
        return true;
    }
    
    if ( !pathed.insert( cn ).second ) return true;
    if ( !lEdge && !rEdge && cn == l ) l = NULL;
    
    for ( ClaimEdge& ce : cn->edges_[1] )
    {
        if ( ce.node == r && !lEdge && !rEdge && l ) return false;
        
        bool edged = lEdge == ( cn == l ) && rEdge == ( ce.node == r );
        if ( !lEdge && !rEdge ) edged = ce.node == r;
        
        if ( !reach( ce.node, edged ? NULL : l, edged ? NULL : r, pathed, diff + ce.diff, score, lEdge, rEdge ) ) return false;
    }
    
    for ( ClaimJoin* cj : cn->joins_[1] ) if ( !reach( cj->node[1], l, !l && cj->node[1] == r ? NULL : r, pathed, diff + cj->dist, score, lEdge, rEdge ) ) return false;
    pathed.erase( cn );
    
    return true;
}

ClaimScore::ClaimScore( Node* l, Node* r, ClaimPairing* cp, int32_t dist, int32_t est )
: pairs{ cp }, diffs( 1 )
{
    node[0] = l;
    node[1] = r;
    for ( int32_t diff : cp->diffs ) diffs[0].push_back( abs( ( dist - diff ) - est ) );
}

void ClaimScore::addRedundant( vector<ClaimNode*> alts[2], int hits )
{
    if ( alts[0].size() > 1 && alts[1].size() > 1 ) return;
    assert( alts[0].size() > 1 || alts[1].size() > 1 );
    for ( ClaimRedundant* cr : pairs[0]->node[0]->redundant_[1] ) if ( cr->add( alts, node, hits ) ) return;
    new ClaimRedundant( alts, node, hits );
}

void ClaimScore::cull( int32_t cutoff )
{
    for ( int i = 0; i < pairs.size(); i++ )
    {
        bool good = false, bad = false;
        for ( int32_t diff : diffs[i] ) ( diff < cutoff ? good : bad ) = true;
        if ( !good ) diffs.erase( diffs.begin() + i );
        if ( !good ) pairs.erase( pairs.begin() + i-- );
    }
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
        
        if ( cp->node[0]->isInvalid( nodes[0][i], offs[0][i], cp->node[1], 1 ) ) return true;
        if ( cp->node[1]->isInvalid( nodes[1][j], offs[1][j], cp->node[0], 0 ) ) return true;
        
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
    cull( cutoff );
    for ( ClaimScore& cs : scores ) cs.cull( cutoff );
    
    vector<ClaimNode*> alts[2];
    for ( int d : { 0, 1 } ) setClaims( alts[d], d );
    for ( ClaimScore& cs : scores ) for ( int d : { 0, 1 } ) cs.setClaims( alts[d], d );
    if ( alts[0].size() != 1 && alts[1].size() != 1 ) return true;
    
    bool distal[2]{ false, false };
    if ( setDistal( distal ) ) return true;
    for ( ClaimScore& cs : scores ) if ( cs.setDistal( distal ) ) return true;
    
    Nodes usedNodes;
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : alts[d] ) cn->setBase( usedNodes, true );
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
    cull( cutoff );
    
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
    
    pairs[0]->hits += hits;
    assert( pairs[0]->diffs.size() == diffs[0].size() );
    for ( int i = 0; i < diffs[0].size(); i++ ) if ( cutoff < diffs[0][i] ) pairs[0]->missed[i] += hits;
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

ClaimBranch::ClaimBranch( ClaimEdge& ce, bool shared, bool drxn )
: edge( ce ), path{ ce.node }, claimed( false ), shared( shared )
{
    if ( !drxn ) while ( path[0]->edges_[0].size() == 1 )
    {
        path.insert( path.begin(), path[0]->edges_[0][0].node );
    }
    if ( drxn ) while ( path.back()->edges_[1].size() == 1 )
    {
        path.push_back( path.back()->edges_[1][0].node );
    }
    fill( ce.node, drxn );
}

bool ClaimBranch::claim( unordered_set<ClaimBranch*> claims[2] )
{
    if ( dumped[1].empty() ) return false;
    assert( !dumped[0].empty() );
    
    // Catalog the paths to keep or discard
    unordered_set<ClaimBranch*> discard[2];
    for ( int i : { 0, 1 } ) for ( pair<ClaimBranch*, int>& cbp : dumped[i] ) ( i ? discard : claims )[1].insert( cbp.first );
    for ( int i : { 0, 1 } ) for ( pair<ClaimBranch*, int>& cbp : dumped[0][0].first->dumped[i] ) ( i ? discard : claims )[0].insert( cbp.first );
    
    // Ensure all path pairs can be cleanly claimed
    for ( int d : { 0, 1 } ) for ( ClaimBranch* cb : claims[d] )
    {
        if ( cb->dumped[0].size() != claims[!d].size() ) return false;
        for ( pair<ClaimBranch*, int>& cbp : cb->dumped[0] ) if ( claims[!d].find( cbp.first ) == claims[!d].end() ) return false;
    }
    
    return true;
}

void ClaimBranch::fill( ClaimNode* cn, bool drxn )
{
    if ( !pathed.insert( cn ).second ) return;
    for ( ClaimEdge& ce : cn->edges_[drxn] ) fill( ce.node, drxn );
    for ( ClaimJoin* cj : cn->joins_[drxn] ) fill( cj->node[drxn], drxn );
}

bool ClaimBranch::isShared( unordered_set<ClaimBranch*>& alts )
{
    if ( shared ) return true;
    for ( pair<ClaimBranch*, int>& cbp : dumped[0] ) if ( !cbp.first->claimed && alts.find( cbp.first ) == alts.end() ) return true;
    return false;
}

ClaimDupe::ClaimDupe( ClaimNode* seed )
: path{ seed }, gen( 0 )
{
    while ( path[0]->edges_[0].size() == 1 ) path.insert( path.begin(), path[0]->edges_[0][0].node );
    while ( path.back()->edges_[1].size() == 1 ) path.push_back( path.back()->edges_[1][0].node );
    forks[0] = path[0];
    forks[1] = path.back();
    edge[0] = edge[1] = NULL;
    for ( int d : { 0, 1 } ) for ( ClaimEdge& ce : forks[d]->edges_[d] ) branches[d].push_back( new ClaimBranch( ce, false, d ) );
    for ( int d : { 0, 1 } ) for ( ClaimBranch* cb : branches[d] ) edged[d].insert( cb->edge.node );
}

ClaimDupe::ClaimDupe( ClaimBranch* seed, ClaimDupe* cd, bool drxn )
: path( seed->path ), gen( cd->gen+1 )
{
    forks[drxn] = drxn ? path.back() : path[0];
    forks[!drxn] = cd->forks[!drxn];
    while ( cd->edge[drxn] ) cd = cd->edge[drxn];
    edge[drxn] = NULL;
    edge[!drxn] = cd;
    cd->edge[drxn] = this;
    for ( ClaimEdge& ce : ( drxn ? path.back() : path[0] )->edges_[drxn] ) branches[drxn].push_back( new ClaimBranch( ce, false, drxn ) );
    for ( ClaimBranch* cb : branches[drxn] ) edged[drxn].insert( cb->edge.node );
//    for ( int d : { 0, 1 } ) edged[d].insert( cd->edged[d].begin(), cd->edged[d].end() );
    edged[!drxn].insert( drxn ? cd->path.back() : cd->path[0] );
}

ClaimDupe::~ClaimDupe()
{
    for ( int d : { 0, 1 } ) for ( ClaimBranch* cb : branches[d] ) delete cb;
    for ( int d : { 0, 1 } ) if ( edge[d] ) edge[d]->edge[!d] = NULL;
}

bool ClaimDupe::advance( ClaimBranch* cb, vector<ClaimNode*>& claims, NodeRoll& nodes, Nodes tried[2], bool& split, bool drxn )
{
    int branched = 0;
    for ( pair<ClaimBranch*, int>& cbp : cb->dumped[0] ) if ( cbp.second > 4 ) branched++;
    if ( branched < 2 ) return false;
    
    ClaimNode* cn = drxn ? cb->path.back() : cb->path[0];
    if ( !cn->edged_[drxn] || cn->edges_[drxn].empty() )
    {
        cn->split_ = split = true;
        return false;
    }
    
    ClaimDupe cd( cb, this, drxn );
    for ( pair<ClaimBranch*, int>& cbp : cb->dumped[0] )
    {
        cd.branches[!drxn].push_back( new ClaimBranch( cbp.first->edge, cbp.first->shared || cbp.first->dumped[0].size() > 1, !drxn ) );
    }
    cout << string( ( gen+2 )*4, ' ' ) << "ADVANCING " << ( drxn ? "RIGHT" : "LEFT" ) << endl;
    if ( cd.resolve( claims, nodes, tried ) ) return true;
    return false;
}

bool ClaimDupe::claim( unordered_set<ClaimBranch*> claimed[2], vector<ClaimNode*>& claims, NodeRoll& nodes )
{
    assert( !edge[0] || !edge[1] );
    vector<ClaimDupe*> dupes = { this };
    while ( dupes[0]->edge[0] ) dupes.insert( dupes.begin(), dupes[0]->edge[0] );
    while ( dupes.back()->edge[1] ) dupes.insert( dupes.end(), dupes.back()->edge[1] );
    
    unordered_set<ClaimNode*> edges[2];
    for ( int d : { 0, 1 } ) for ( ClaimEdge& ce : ( d ? dupes.back()->path.back() : dupes[0]->path[0] )->edges_[d] ) edges[d].insert( ce.node );
    for ( int d : { 0, 1 } ) for ( ClaimBranch* cb : claimed[d] ) if ( edges[d].find( cb->edge.node ) == edges[d].end() ) return false;
    
    bool unduped[2]{ edges[0].size() == claimed[0].size(), edges[1].size() == claimed[1].size() };
    for ( int d : { 0, 1 } ) for ( ClaimBranch* cb : claimed[d] ) if ( cb->shared ) unduped[d] = false;
    
    if ( dupes.size() == 1 && unduped[0] && unduped[1] ) return true;
    if ( dupes.size() == 1 ) assert( !unduped[0] && !unduped[1] );
    if ( unduped[0] && unduped[1] ) assert( dupes.size() != 2 );
    
    vector< pair<ClaimNode*, ClaimNode*> > cloned;
    for ( int i = unduped[0]; i + unduped[1] < dupes.size(); i++ ) dupes[i]->setPath( cloned, claims, nodes );
    
    ClaimRepair cr( dupes[ unduped[0] ]->path[0]->getFork( 0 ), dupes.end()[ -1-unduped[1] ]->path.back()->getFork( 1 ) );
    for ( int d : { 0, 1 } ) if ( unduped[d] ) ( d ? dupes.back() : dupes[0] )->setDettached( cloned, d );
    for ( int d : { 0, 1 } ) if ( !unduped[d] ) for ( ClaimBranch* cb : claimed[d] )
    {
        ( d ? cloned.back() : cloned[0] ).second->addEdge( cb->edge.node, cb->edge.ol, cb->edge.diff, cb->edge.isLeap, d, true, true );
        if ( !cb->isShared( claimed[!d] ) ) ( d ? cloned.back() : cloned[0] ).first->removeEdge( cb->edge.node, d );
    }
    cr.dupe( cloned, nodes );
    
    return true;
}

bool ClaimDupe::dupe( ClaimNode* fork, vector<ClaimNode*>& claims, NodeRoll& nodes, Nodes tried[2] )
{
    ClaimDupe cd( fork );
    bool duped = cd.resolve( claims, nodes, tried );
    cout << endl;
    return duped;
}

bool ClaimDupe::resolve( vector<ClaimNode*>& claims, NodeRoll& nodes, Nodes tried[2] )
{
    cout << string( ( gen+1 )*4, ' ' ) << "Duplicate attempt, coords: " << forks[0]->getFork( 0 )->ends_[0] << " to " << forks[1]->getFork( 1 )->ends_[1] << ", paths: " << branches[0].size() << "-" << branches[1].size() << endl;
    
    bool duped = false, split = false, announced = false;
    
    // Pair branches and dump false pairings
    for ( ClaimBranch* l : branches[0] ) for ( ClaimBranch* r : branches[1] )
    {
        int unique = setUniques( l, r );
        bool dumped = !unique && setDump( l, r );
        l->dumped[dumped].push_back( make_pair( r, unique ) );
        r->dumped[dumped].push_back( make_pair( l, unique ) );
    }
    
    // Claim any clean pairings
    for ( ClaimBranch* l : branches[0] )
    {
        unordered_set<ClaimBranch*> claimed[2];
        if ( l->claim( claimed ) && claim( claimed, claims, nodes ) ) duped = true;
    }
    if ( duped && ( announced = true ) ) cout << string( ( gen+1 )*4, ' ' ) << "DUPED CLEANLY" << endl;
    
    // Extend any branch that pairs strongly to more than one branch
    if ( !duped ) for ( int d : { 0, 1 } ) for ( ClaimBranch* cb : branches[d] ) if ( !duped && ( duped = advance( cb, claims, nodes, tried, split, d ) ) ) break;
    if ( duped && !announced && ( announced = true ) ) cout << string( ( gen+1 )*4, ' ' ) << "DUPED AHEAD" << endl;
    
    // Claim any branch than has only one pairing, but that pair has multiple in return
    if ( !duped && !split )
    {
        vector<ClaimBranch*> claimable[2];
        for ( int d : { 0, 1 } ) for ( ClaimBranch* cb : branches[d] ) if ( cb->dumped[0].size() == 1 ) claimable[d].push_back( cb );
        for ( int d : { 0, 1 } ) for ( ClaimBranch* cb : claimable[d] )
        {
            unordered_set<ClaimBranch*> claimed[2];
            claimed[d].insert( cb );
            claimed[!d].insert( cb->dumped[0][0].first );
//            if ( claim( claimed, claims, nodes ) ) cb->claimed = duped = true;
        }
    }
    if ( duped && !announced && ( announced = true ) ) cout << string( ( gen+1 )*4, ' ' ) << "DUPED SOLO" << endl;
    if ( !announced && split ) cout << string( ( gen+1 )*4, ' ' ) << "SET SPLIT" << endl;
    else if ( !announced && !duped ) cout << string( ( gen+1 )*4, ' ' ) << "NOT DUPED!" << endl;
    
    if ( !duped && !split && !gen ) for ( ClaimNode* cn : path ) for ( Node* node : cn->path_ ) tried[0] += node;
    
    return duped;
}

void ClaimDupe::setDettached( vector< pair<ClaimNode*, ClaimNode*> >& cloned, bool drxn )
{
    ClaimEdge ce = drxn ? path[0]->getEdge( cloned.back().first, 0 ) : path.back()->getEdge( cloned[0].first, 1 );
    assert( ce.node );
    ( drxn ? path[0] : path.back() )->addEdge( ( drxn ? cloned.back() : cloned[0] ).second, ce.ol, ce.diff, ce.isLeap, !drxn, true, true );
    ( drxn ? path[0] : path.back() )->removeEdge( ce.node, !drxn );
}

bool ClaimDupe::setDump( ClaimBranch* l, ClaimBranch* r )
{
    unordered_set<ClaimNode*> alts[2];
    for ( ClaimBranch* cb : branches[0] ) if ( cb != l ) for ( ClaimNode* cn : cb->pathed ) if ( l->pathed.find( cn ) == l->pathed.end() ) alts[0].insert( cn );
    for ( ClaimBranch* cb : branches[1] ) if ( cb != r ) for ( ClaimNode* cn : cb->pathed ) if ( r->pathed.find( cn ) == r->pathed.end() ) alts[1].insert( cn );
    for ( ClaimBranch* ll : branches[0] ) if ( ll != l ) for ( ClaimBranch* rr : branches[1] ) if ( rr != r )
    {
        ClaimBranch* cb[2][2]{ { l, ll }, { r, rr } };
        int hits[2]{ setScore( cb, alts, 0, 1 ), setScore( cb, alts, 1, 0 ) };
        if ( !hits[0] || !hits[1] ) continue;
        bool success = min( hits[0], hits[1] ) > 1 && hits[0] + hits[1] > 6;
        cout << string( ( gen+2 )*4, ' ' ) << ( success ? "DUPED: " : "NOT DUPED: " ) << hits[0] << "-" << hits[1] << endl;
        if ( min( hits[0], hits[1] ) > 1 && hits[0] + hits[1] > 6 ) return true;
    }
    return false;
}

void ClaimDupe::setPath( vector< pair<ClaimNode*, ClaimNode*> >& cloned, vector<ClaimNode*>& claims, NodeRoll& nodes )
{
    for ( int i = 0; i < path.size(); i++ )
    {
        ClaimNode* clone = new ClaimNode( path[i], claims, nodes );
        for ( int d : { 0, 1 } ) for ( ClaimEdge& ce : path[i]->edges_[d] )
        {
            if ( !d && !cloned.empty() && ce.node == cloned.back().first )
                clone->addEdge( cloned.back().second, ce.ol, ce.diff, ce.isLeap, 0, true, true );
            if ( !d && ( i ? ce.node != path[i-1] : edged[0].find( ce.node ) == edged[0].end() ) )
                clone->addEdge( ce.node, ce.ol, ce.diff, ce.isLeap, 0, true, true );
            if ( d && ( i+1 < path.size() ? ce.node != path[i+1] : edged[1].find( ce.node ) == edged[1].end() ) )
                clone->addEdge( ce.node, ce.ol, ce.diff, ce.isLeap, 1, true, true );
        }
        cloned.push_back( make_pair( path[i], clone ) );
    }
}

int ClaimDupe::setScore( ClaimBranch* cb[2][2], unordered_set<ClaimNode*> alts[2], int l, int r )
{
    int hits = 0;
    unordered_set<ClaimNode*> tar[2];
    unordered_set<ClaimRedundant*> used;
    for ( ClaimNode* cn : cb[0][l]->pathed ) if ( cb[0][!l]->pathed.find( cn ) == cb[0][!l]->pathed.end() ) tar[0].insert( cn );
    for ( ClaimNode* cn : cb[1][r]->pathed ) if ( cb[1][!r]->pathed.find( cn ) == cb[1][!r]->pathed.end() ) tar[1].insert( cn );
    
    for ( ClaimNode* cn : tar[0] )
    {
        for ( ClaimPairing* cp : cn->pairs_[1] ) if ( tar[1].find( cp->node[1] ) != tar[1].end() ) 
            hits += cp->get( cb[0][l]->path.back(), cb[1][r]->path[0], false, false );
        for ( ClaimRedundant* cr : cn->redundant_[1] ) if ( used.insert( cr ).second )
        {
            bool good = !r;
            if ( r ) for ( ClaimNode* cnp : cr->nodes[1] ) if ( tar[1].find( cnp ) != tar[1].end() ) good = true;
            if ( good ) hits += cr->get( ( l ? alts : tar )[0], ( r ? alts : tar )[1], cb[0][l]->path.back(), cb[1][r]->path[0], false, false );
        }
    }
    return hits;
}

int ClaimDupe::setUniques( ClaimBranch* l, ClaimBranch* r )
{
    int unique = 0;
    unordered_set<ClaimNode*> tar[2]{ l->pathed, r->pathed };
    for ( ClaimBranch* cd : branches[0] ) if ( cd != l ) for ( ClaimNode* cn : cd->pathed ) tar[0].erase( cn );
    for ( ClaimBranch* cd : branches[1] ) if ( cd != r ) for ( ClaimNode* cn : cd->pathed ) tar[1].erase( cn );
    unordered_set<ClaimRedundant*> redundants;
    for ( ClaimNode* cn : tar[0] )
    {
        for ( ClaimPairing* cp : cn->pairs_[1] ) if ( tar[1].find( cp->node[1] ) != tar[1].end() )
            unique += cp->get( l->path.back(), r->path[0], false, false );
        for ( ClaimRedundant* cr : cn->redundant_[1] ) if ( redundants.insert( cr ).second )
            unique += cr->get( tar[0], tar[1], l->path.back(), r->path[0], false, false );
    }
    return unique;
}

ClaimTrim::ClaimTrim( ClaimNode* l, ClaimNode* r, ClaimNode* bridge )
{
    forks[0] = l;
    forks[1] = r;
    path = bridge ? vector<ClaimNode*>{ l, bridge, r } : vector<ClaimNode*>{ l, r };
    for ( int d : { 0, 1 } )
    {
        ols[d] = 0;
        fill( forks[d], groups[d][0], d );
        for ( ClaimEdge& re : forks[!d]->edges_[d] ) if ( re.node != forks[d] && re.node != bridge )
        {
            fill( re.node, groups[d][1], d );
            ols[d] = max( ols[d], re.ol );
        }
        
        vector<ClaimNode*> shared;
        for ( ClaimNode* cn : groups[d][0] ) if ( groups[d][1].find( cn ) != groups[d][1].end() ) shared.push_back( cn );
        for ( ClaimNode* cn : shared ) for ( int i : { 0, 1 } ) groups[d][i].erase( cn );
    }
}

void ClaimTrim::fill( ClaimNode* cn, unordered_set<ClaimNode*>& group, bool drxn )
{
    if ( !group.insert( cn ).second ) return;
    for ( ClaimEdge& ce : cn->edges_[drxn] ) fill( ce.node, group, drxn );
    for ( ClaimJoin* cj : cn->joins_[drxn] ) fill( cj->node[drxn], group, drxn );
}

int ClaimTrim::score( int l, int r )
{
    int hits = 0;
    unordered_set<ClaimRedundant*> used;
    for ( ClaimNode* cn : groups[0][l] )
    {
        for ( ClaimPairing* cp : cn->pairs_[1] ) if ( groups[1][r].find( cp->node[1] ) != groups[1][r].end() )
            hits += cp->get( path, !l && !r ? 2 : (int)r );
        for ( ClaimRedundant* cr : cn->redundant_[1] ) if ( used.insert( cr ).second )
            hits += cr->get( groups[0][l], groups[1][r], forks[0], forks[1], !l, !r );
    }
    return hits;
}

void ClaimTrim::score( int l, int r, int& hits, int& unique )
{
//    unordered_set<ClaimNode*> pathable[2];
//    for ( int i : { 0, 1 } ) for ( ClaimNode* cn : groups[0][i] ) pathable[ i == l ].insert( cn );
//    for ( int i : { 0, 1 } ) for ( ClaimNode* cn : groups[1][i] ) pathable[ i == r ].insert( cn );
    bool bools[2]{ !l, !r };
    for ( ClaimNode* cn : groups[0][l] ) for ( ClaimPairing* cp : cn->pairs_[1] ) if ( groups[1][r].find( cp->node[1] ) != groups[1][r].end() )
    {
        int test = cp->get( path, !l && !r ? 2 : (int)r );
        vector<int32_t> diffs;
        unordered_set<ClaimNode*> pathed;
        if ( !cn->reachViaEdgeOnly( cp->node[1], forks, pathed, 0, diffs, bools, false, 1 ) ) continue;
//        if ( !cn->reach( cp->node[1], 0, pathable, diffs, false, false, 1 ) ) continue;
        assert( !diffs.empty() );
        int pairScore = cp->hit( diffs );
        assert( pairScore == test );
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
            unordered_set<ClaimNode*> pathed;
            if ( !cnl->reachViaEdgeOnly( cnr, forks, pathed, 0, diffs, bools, false, 1 ) ) good = false;
        }
        if ( good ) shares.insert( cr );
    }
    for ( ClaimRedundant* cr : shares )
    {
        hits += cr->score;
        unique += cr->score;
    }
}

bool ClaimTrim::trim( ClaimNode* fork, NodeRoll& nodes, bool drxn )
{
    for ( ClaimEdge& ce : fork->edges_[drxn] )
    {
        if ( ce.node->edges_[!drxn].size() > 1 )
        {
            ClaimTrim ct( drxn ? fork : ce.node, drxn ? ce.node : fork, NULL );
            
            cout << "    Trim attempt, coords: " << fork->getFork( drxn )->ends_[1] << ", edges: " << fork->edges_[drxn].size() << endl; 
            
            if ( ct.trim( ce.ol ) )
            {
                cout << "    TRIMMED!" << endl << endl;
                ClaimRepair cr( ( drxn ? ce.node : fork )->getFork( 1 ), ( drxn ? fork : ce.node )->getFork( 0 ) );
                assert( fork->removeEdge( ce.node, drxn ) );
                cr.trim();
                nodes.updateBad();
                if ( fork->edges_[drxn].size() > 1 ) trim( fork, nodes, drxn );
                return true;
            }
            
            cout << "    FAILED!" << endl << endl;
        }
        
        if ( ce.node->isBridge() )
        {
            cout << "    Trim attempt, coords: " << fork->getFork( drxn )->ends_[1] << ", bridge with: " << ce.node->edges_[0].size() << "-" << ce.node->edges_[1].size() <<  " edges" << endl; 
            
            bool failed = false;
            for ( ClaimEdge& l : ce.node->edges_[0] ) for ( ClaimEdge& r : ce.node->edges_[1] )
            {
                ClaimTrim ct( l.node, r.node, ce.node );
                if ( !ct.trim( min( l.ol, r.ol ) ) ) failed = true;
            }
            
            if ( !failed )
            {
                cout << "    TRIMMED!" << endl << endl;
                vector<Node*> forks[2];
                for ( int d : { 0, 1 } ) for ( ClaimEdge& re : ce.node->edges_[d] ) forks[!d].push_back( re.node->getFork( !d ) );
                ClaimRepair cr( forks[0], forks[1] );
                ce.node->dismantle();
                cr.trim();
                nodes.updateBad();
                return true;
            }
            
            cout << "    FAILED!" << endl << endl;
        }
    }
    
    return false;
}

bool ClaimTrim::trim( int ol )
{
    for ( int d : { 0, 1 } ) if ( groups[d][0].empty() || groups[d][1].empty() ) return false;
    
    int hits[3]{ score( 0, 1 ), score( 1, 0 ), score( 0, 0 ) }, diff = min( 100, min( 60, ols[0] - ol ) + min( 60, ols[1] - ol ) ) / 10;
    
    // Filter out contended paths
    if ( !hits[0] || !hits[1] || min( hits[0], hits[1] ) / 15 < hits[2] ) return false;
    
    bool failed = min( hits[0], hits[1] ) - hits[2] < 2;
    if ( !( hits[0] + hits[1] + diff - hits[2] > 8 ) ) failed = true;
    
    cout << "        " << ( failed ? "NOT TRIMMED: " : "TRIMMED: " ) << hits[0] << "-" << hits[1] << endl;
    
    if ( min( hits[0], hits[1] ) - hits[2] < 2 ) return false;
    
    return hits[0] + hits[1] + diff - hits[2] > 8;
}

ClaimNode::ClaimNode( ClaimNode* cn, vector<ClaimNode*>& claims, NodeRoll& nodes )
: clone_( this ), split_( false )
{
    for ( int d : { 0, 1 } )
    {
        ends_[d] = cn->ends_[d];
        edged_[d] = cn->edged_[d];
    }
    
    Nodes base( cn->path_ );
    for ( int d : { 0, 1 } ) for ( ClaimEdge& ce : cn->edges_[d] ) base += ce.node->getFork( !d );
    
    for ( Node* node : cn->path_ )
    {
        path_.push_back( new Node( node, nodes, node->drxn_, node->bad_ ) );
        for ( int d : { 0 , 1 } ) for ( Edge& e : node->edges_[d] ) if ( !base.find( e.node ) ) path_.back()->addEdge( e, d, true );
    }
    
    for ( int i = 1; i < path_.size(); i++ ) for ( Edge& e : cn->path_[i]->edges_[0] ) if ( e.node == cn->path_[i-1] )
    {
        path_[i]->addEdge( path_[i-1], e.ol, 0, false, e.leap );
    }
    claims.push_back( this );
}

ClaimNode::ClaimNode( ClaimNode* cn, NodeRoll& cloned, Nodes& base )
: clone_( this ), split_( false )
{
    for ( Node* node : cn->path_ )
    {
        path_.push_back( new Node( node, cloned, node->drxn_, node->bad_ ) );
        for ( int d : { 0 , 1 } ) for ( Edge& e : node->edges_[d] ) if ( !base.find( e.node ) ) path_.back()->addEdge( e, d, true );
    }
    
    for ( int i = 1; i < path_.size(); i++ ) path_[i]->addEdge( cn->path_[i]->getEdge( path_[i-1], 0 ), 0, true );
}

ClaimNode::ClaimNode( Node* node, vector<ClaimNode*>& claims, int32_t coord, bool orient, bool drxn )
: clone_( this ), path_{ node }, offs_{ make_pair( node, coord ) }, split_( false )
{
    claims.push_back( this );
    ends_[0] = ends_[1] = coord;
    edged_[0] = edged_[1] = false;
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

void ClaimNode::addEdge( ClaimNode* cn, int ol, int32_t diff, bool isLeap, bool drxn, bool reciprocate, bool nodeEdge )
{
    for ( ClaimEdge& ce : edges_[drxn] ) if ( ce.node == cn ) return;
    edges_[drxn].push_back( ClaimEdge( cn, ol, diff, isLeap ) );
    if ( nodeEdge ) getFork( drxn )->addEdge( cn->getFork( !drxn ), ol, drxn, false, isLeap );
    if ( reciprocate ) cn->addEdge( this, ol, diff, isLeap, !drxn, false );
}

bool ClaimNode::addEdge( Edge& e, vector<ClaimNode*>& claims, bool orient, bool drxn )
{
    int32_t coord = getCoord( NULL, e, ends_[drxn], orient, drxn );
    bool edged = false;
    for ( ClaimNode* cn : claims ) if ( cn != this && cn->findFork( e.node, !drxn ) )
    {
        addEdge( cn, e.ol, drxn ? cn->ends_[0] - coord : coord - cn->ends_[1], e.leap, drxn );
        edged = true;
    }
    if ( !edged ) addEdge( new ClaimNode( e.node, claims, coord, orient, drxn ), e.ol, 0, e.leap, drxn );
    
    return edged_[drxn] = true;
}

void ClaimNode::addJoin( ClaimNode* cn, Node* hang, NodeDists& dists, int32_t dist, bool drxn )
{
    if ( dist < -500 ) return;
    ClaimJoin* joined = NULL;
    for ( ClaimJoin* cj : joins_[drxn] ) if ( cj->node[drxn] == cn ) joined = cj;
    if ( joined ) joined->dist = min( joined->dist, dist );
    else joined = new ClaimJoin( this, cn, dist, drxn );
    joined->fill( hang, dists, !drxn );
}

bool ClaimNode::branch( vector<ClaimNode*>& claims, vector<ClaimNode*>& forks, bool orient )
{
    assert( !forks.empty() );
    bool success = false;
    unordered_set<ClaimNode*> tested[2];
    for ( ClaimNode* cn : forks )
    {
        bool branched = false;
        for ( ClaimEdge& ce : cn->edges_[orient] ) ce.node->setBranched( claims, tested[orient], cn->edges_[orient].size() > 1, orient, orient );
        for ( ClaimEdge& ce : cn->edges_[orient] ) if ( ce.node->edges_[!orient].size() > 1 ) branched = true;
        if ( cn->setBranched( claims, tested[!orient], branched, orient, !orient ) || branched ) success = true;
    }
    return success;
    
//    if ( forks.size() > 1 || forks[0]->isBranched( !orient ) ) return true;
//    for ( ClaimEdge& ce : forks[0]->edges_[orient] ) if ( ce.node->edges_[!orient].size() > 1 ) return true;
//    if ( !forks[0]->split( claims, orient, !orient ) ) return false;
//    return true;
}

bool ClaimNode::branch( vector<ClaimNode*>& claims, unordered_set<ClaimNode*>& tested, bool branched, bool orient, bool drxn )
{
    split_ = false;
    
    if ( !tested.insert( this ).second ) return branched;
    
    if ( !edged_[drxn] && ( !branched || !edges_[drxn].empty() ) ) edge( claims, orient, drxn );
    
    if ( edges_[drxn].empty() ) return branched;
    
    if ( edges_[drxn].size() > 1 )
    {
        for ( ClaimNode* cn = this; cn->edges_[!drxn].size() == 1; )
        {
            ClaimNode* nxt = cn->edges_[!drxn][0].node;
            cn->edge( claims, orient, !drxn );
            for ( ClaimEdge& ce : cn->edges_[!drxn] ) if ( ce.node != nxt ) ce.node->setBranched( claims, false, orient, !drxn );
            cn = nxt;
        }
        
        for ( ClaimEdge& ce : edges_[drxn] ) ce.node->setBranched( claims, false, orient, drxn );
        for ( ClaimEdge& ce : edges_[drxn] ) ce.node->branch( claims, tested, true, orient, drxn );
        
        return true;
    }
    
    return edges_[drxn][0].node->branch( claims, tested, branched, orient, drxn );
}

bool ClaimNode::claim( Node* node, NodeRoll& nodes, Nodes tried[2], bool drxn )
{
    if ( tried[0].find( node ) && tried[1].find( node ) ) return false;
    
    cout << "Claim coords: " << node->ends_[0] << " " << node->ends_[1] << endl << endl;
    
    bool claimed = false;
    
    for ( int again = 1; again-- > 0; )
    {
        vector<ClaimNode*> claims;
        ClaimNode* fork = new ClaimNode( node, claims, 0, drxn, !drxn );
        for ( int retry = 1; retry-- > 0; )
        {
            if ( !fork->create( claims, drxn ) ) break;;
            complete( claims, drxn );
            score( claims, drxn );
            unordered_set<ClaimNode*> trimmed[2], duped[2];
            if ( deloop( claims ) || fork->trim( trimmed, nodes, tried, 2 ) || fork->dupe( duped, claims, nodes, tried, 2 ) )
            {
                if ( fork->edges_[drxn].size() > 1 ) again = 1;
                tried[1].clear();
                claimed = true;
            }
            else retry = resplit( claims, drxn );
        }
        
        for ( ClaimNode* cn : claims ) delete cn;
    }
    
    cout << ( claimed ? "CLAIMED!" : "NOT CLAIMED!" ) << endl << endl;
    
    return claimed;
}

bool ClaimNode::deloop( vector<ClaimNode*>& claims )
{
    for ( ClaimNode* cn : claims ) for ( int d : { 0, 1 } ) if ( cn->edges_[d].size() == 1 )
    {
        vector<ClaimEdge> path;
        while ( cn && cn->edges_[d].size() == 1 )
        {
            ClaimEdge ce = cn->edges_[d][0];
            for ( int i = 0; i < path.size(); i++ ) if ( path[i].node == ce.node )
            {
                for ( int j = i+1; j < path.size(); j++ ) if ( path[j].ol < ce.ol )
                {
                    ce = path[j];
                    cn = path[j-1].node;
                }
                assert( cn->getFork( d )->removeEdge( ce.node->getFork( !d ), d, true ) );
                return true;
            }
            path.push_back( cn->edges_[d][0] );
            cn = path.back().node;
        }
    }
    return false;
}

void ClaimNode::complete( vector<ClaimNode*>& claims, bool orient )
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
}

bool ClaimNode::create( vector<ClaimNode*>& claims, bool orient )
{
    for ( int again = 1; again-- > 0; )
    {
        unordered_set<ClaimNode*> tested[2];
        for ( int d : { 0, 1 } ) branch( claims, tested[d], false, orient, d );
    
        Nodes branches[2];
        vector<ClaimNode*> unended[2];
        for ( int d : { 0, 1 } ) for ( ClaimNode* cn : claims ) if ( !cn->edged_[d] ) unended[d].push_back( cn );
        for ( int d : { 0, 1 } ) for ( ClaimNode* cn : unended[d] ) if ( !cn->edges_[d].empty() ) for ( Edge& e : cn->getEdges( NULL, d ) ) branches[d] += e.node;
        for ( int d : { 0, 1 } ) for ( ClaimNode* cn : unended[d] )
        {
            bool doEdge = branches[!d].find( cn->getFork( d ) );
            for ( Edge& e: cn->getEdges( NULL, d ) ) if ( branches[!d].find( e.node ) ) doEdge = true;
            if ( doEdge && cn->edge( claims, orient, d ) ) again = 1;
        }
    }
    
    return true;
}

void ClaimNode::dismantle()
{
    for ( int d : { 0, 1 } ) while ( !edges_[d].empty() ) removeEdge( edges_[d].back().node, d );
    for ( Node* node : path_ ) node->dismantle();
}

bool ClaimNode::dupe( unordered_set<ClaimNode*> duped[2], vector<ClaimNode*>& claims, NodeRoll& nodes, Nodes tried[2], int drxn )
{
    for ( int d : { 0, 1 } ) if ( drxn == 2 || d == drxn )
    {
        if ( !duped[d].insert( this ).second ) return false;
        for ( ClaimEdge& ce : edges_[d] ) if ( ce.node->dupe( duped, claims, nodes, tried, d ) ) return true;
    }
    
    if ( tried[0].find( path_[0] ) && tried[0].find( path_.back() ) ) return false;
    
    for ( int d : { 0, 1 } ) if ( (bool)d == (bool)drxn && edges_[d].size() > 1 )
    {
        if ( isBranched( !d ) && ClaimDupe::dupe( this, claims, nodes, tried ) ) return true;
    }
    
    return false;
}

bool ClaimNode::edge( vector<ClaimNode*>& claims, bool orient, bool drxn )
{
    if ( edged_[drxn] ) return false;
    if ( !getFork( drxn )->isForkComplete( params.shortLen(), 20, drxn ) ) return false;
    for ( Edge& e : getFork( drxn )->edges_[drxn] ) addEdge( e, claims, orient, drxn );
    return edged_[drxn] = true;
}

void ClaimNode::extend( vector<ClaimNode*>& claims, bool orient, bool drxn )
{
    for ( vector<Edge> edges = getEdges( NULL, drxn ); edges.size() == 1 && edges_[drxn].empty(); )
    {
        if ( getEdges( edges[0].node, !drxn ).size() > 1 && addEdge( edges[0], claims, orient, drxn ) ) return;
        ends_[drxn] = getCoord( NULL, edges[0], ends_[drxn], orient, drxn );
        path_.insert( drxn ? path_.end() : path_.begin(), edges[0].node );
        offs_.insert( make_pair( edges[0].node, ends_[drxn] ) );
        edges = getEdges( NULL, drxn );
        if ( edges.empty() ) edged_[drxn] = true;
    }
}

void ClaimNode::fill( Node* node, Nodes& block, int32_t dist, int32_t limit, bool orient, bool drxn )
{
    if ( block.find( node ) || abs( dist ) > abs( limit ) ) return;
    
    auto it = offs_.insert( make_pair( node, dist ) );
    if ( !it.second && ( drxn ? it.first->second <= dist : dist <= it.first->second ) ) return;
    if ( !it.second ) it.first->second = dist;
    
    for ( Edge& e : getEdges( node, drxn ) ) fill( e.node, block, getCoord( node, e, dist, orient, drxn ), limit, orient, drxn );
}

void ClaimNode::fill( Node* node, Nodes& include, int32_t dist, bool orient, bool drxn )
{
    for ( Edge& e : getEdges( node, drxn ) ) if ( include.find( e.node ) )
    {
        int32_t coord = getCoord( node, e, dist, orient, drxn );
        auto it = offs_.insert( make_pair( e.node, coord ) );
        if ( !it.second && ( drxn ? it.first->second <= coord : coord <= it.first->second ) ) continue;
        if ( !it.second ) it.first->second = coord;
        fill( e.node, include, dist, orient, drxn );
    }
}

bool ClaimNode::findFork( Node* q, bool drxn )
{
    return q == getFork( drxn );
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

ClaimEdge ClaimNode::getEdge( ClaimNode* cn, bool drxn )
{
    ClaimEdge edge( NULL, 0, 0, false );
    for ( ClaimEdge& ce : edges_[drxn] ) if ( ce.node == cn ) edge = ce;
    return edge;
}

vector<Edge> ClaimNode::getEdges( Node* fork, bool drxn, bool blunt )
{
    vector<Edge> edges;
    bool isEnd = !fork;
    if ( !fork ) fork = getFork( drxn );
    for ( Edge& e : fork->edges_[drxn] ) if ( blunt || !e.node->isBlunt( 0, 3 , drxn ) )
    {
        bool added = false;
        if ( isEnd ) for ( ClaimEdge& ce : edges_[drxn] ) if ( e.node == ce.node->getFork( !drxn ) ) added = true;
        if ( !added ) edges.push_back( Edge( e.node, e.ol, e.leap ) );
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
    if ( !edged_[0] || !edged_[1] || edges_[0].empty() || edges_[1].empty() ) return false;
    for ( int d : { 0, 1 } ) for ( ClaimEdge& ce : edges_[d] ) if ( ce.node->edges_[!d].size() == 1 ) return false;
    
    int readCount = 0;
    for ( Node* node : path_ ) readCount += node->countReads( true );
    if ( readCount < 2 ) return true;
    
    for ( Node* node : path_ ) if ( !node->cloned_ ) return false;
    
    for ( int d : { 0, 1 } ) for ( ClaimEdge& ce : edges_[d][0].node->edges_[!d] ) if ( ce.node != this )
    {
        if ( ( d ? ce.node->path_.back() : ce.node->path_[0] )->isClone( path_, 0, !d ) ) return true;
    }
    return false;
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
    vector<ClaimNode*> path;
    for ( ClaimJoin* cj : joins_[drxn] ) if ( cj->node[drxn]->reach( cn, path, drxn ) ) if ( cj->bridge.find( node ) ) return false;
    return true;
}

bool ClaimNode::isShared( Node* node, ClaimNode* cn, bool drxn )
{
    vector<ClaimNode*> path;
    for ( ClaimShared& cs : shared_ ) if ( cs.shared.find( node ) && cs.node->reach( cn, path, drxn ) ) return true;
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

bool ClaimNode::reach( ClaimNode* tar, vector<ClaimNode*>& path, bool drxn )
{
    if ( tar == this ) return true;
    if ( find( path.begin(), path.end(), this ) != path.end() ) return false;
    bool reached = false;
    path.push_back( this );
    for ( ClaimEdge& ce : edges_[drxn] ) if ( !reached && ce.node->reach( tar, path, drxn ) ) reached = true;
    for ( ClaimJoin* cj : joins_[drxn] ) if ( !reached && cj->node[drxn]->reach( tar, path, drxn ) ) reached = true;
    path.pop_back();
    return reached;
}

bool ClaimNode::reachViaEdgeOnly( ClaimNode* tar, ClaimNode* forks[2], unordered_set<ClaimNode*>& pathed, int32_t diff, vector<int32_t>& diffs, bool path[2], bool forked, bool drxn )
{
    if ( this == tar )
    {
        if ( !forked ) return false;
        if ( find( diffs.begin(), diffs.end(), diff ) == diffs.end() ) diffs.push_back( diff );
        return true;
    }
    
    if ( !pathed.insert( this ).second ) return true;
    
    for ( ClaimEdge& ce : edges_[drxn] )
    {
        bool forkPath = forked;
        if ( forks[!drxn] == this ) forkPath = path[!drxn] && ( ( forks[drxn] == ce.node ) == path[drxn] );
        if ( forks[drxn] == ce.node ) forkPath = path[drxn] && ( ( forks[!drxn] == this ) == path[!drxn] );
        if ( !ce.node->reachViaEdgeOnly( tar, forks, pathed, diff + ce.diff, diffs, path, forkPath, drxn ) ) return false;
    }
    for ( ClaimJoin* cj : joins_[drxn] )
    {
        bool forkPath = forked;
        if ( forks[!drxn] == this ) forkPath = path[!drxn] && !path[drxn] && ( forks[drxn] != cj->node[drxn] );
        if ( forks[drxn] == cj->node[drxn] ) forkPath = path[drxn] && !path[!drxn] && ( forks[!drxn] != this );
        if ( !cj->node[drxn]->reachViaEdgeOnly( tar, forks, pathed, diff + cj->diff( drxn ), diffs, path, forkPath, drxn ) ) return false;
    }
    
    pathed.erase( this );
    return true;
}

bool ClaimNode::reachViaForkOnly( ClaimNode* tar, ClaimNode* fork, unordered_set<ClaimNode*>& pathed, int32_t diff, vector<int32_t>& diffs, bool forked, bool drxn )
{
    if ( this == fork ) forked = true;
    if ( this == tar )
    {
        if ( !forked ) return false;
        if ( find( diffs.begin(), diffs.end(), diff ) == diffs.end() ) diffs.push_back( diff );
        return true;
    }
    
    if ( !pathed.insert( this ).second ) return true;
    
    for ( ClaimEdge& ce : edges_[drxn] ) if ( !ce.node->reachViaForkOnly( tar, fork, pathed, diff + ce.diff, diffs, forked, drxn ) ) return false;
//    for ( ClaimJoin* cj : joins_[drxn] ) if ( !cj->node[drxn]->reach( tar, fork, diff + cj->diff( drxn ), diffs, forked, drxn ) ) return false;
    for ( ClaimJoin* cj : joins_[drxn] ) if ( !cj->node[drxn]->reachViaForkOnly( tar, fork, pathed, diff + cj->dist, diffs, forked, drxn ) ) return false;
    
    pathed.erase( this );
    return true;
}

//bool ClaimNode::reach( ClaimNode* tar, int32_t diff, vector<int32_t>& diffs, bool drxn )
//{
//    if ( tar == this )
//    {
//        if ( find( diffs.begin(), diffs.end(), diff ) == diffs.end() ) diffs.push_back( diff );
//        return true;
//    }
//    
//    bool reached = false;
//    for ( ClaimEdge& ce : edges_[drxn] ) if ( ce.node->reach( tar, diff + ce.diff, diffs, drxn ) ) reached = true;
//    for ( ClaimJoin* cj : joins_[drxn] ) if ( cj->node[drxn]->reach( tar, diff + cj->diff( drxn ), diffs, drxn ) ) reached = true;
//    return reached;
//}

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
        assert( getFork( drxn )->removeEdge( cn->getFork( !drxn ), drxn, true ) );
        cn->removeEdge( this, !drxn, false );
    }
    return removed;
}

//bool ClaimNode::removeJoin( ClaimNode* cn, bool drxn, bool reciprocate )
//{
//    bool removed = false;
//    for ( int i = 0; i < joins_[drxn].size(); i++ ) if ( joins_[drxn][i]->node[drxn] == cn )
//    {
//        joins_[drxn].erase( joins_[drxn].begin() + i-- );
//        removed = true;
//    }
//    if ( removed && reciprocate ) cn->removeJoin( this, !drxn, false );
//    return removed;
//}

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

bool ClaimNode::resplit( vector<ClaimNode*>& claims, bool orient )
{
    bool didSplit = false;
    
    unordered_set<ClaimNode*> splitable[2], tested[2];
    for ( ClaimNode* cn : claims ) if ( cn->split_ ) for ( int d : { 0, 1 } ) if ( cn->edges_[d].empty() ) splitable[d].insert( cn );
    
    if ( splitable[0].empty() && splitable[1].empty() ) return false;
    
    for ( ClaimNode* cn : claims )
    {
        cn->reset();
        for ( auto it = cn->offs_.begin(); it != cn->offs_.end(); )
        {
            if ( cn->ends_[0] <= it->second && it->second <= cn->ends_[1] ) it++;
            else it = cn->offs_.erase( it );
        }
    }
    
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : splitable[d] ) if ( cn->branch( claims, tested[d], false, orient, d ) ) didSplit = true;
    
    return didSplit;
}

void ClaimNode::score( vector<ClaimNode*>& claims, bool orient )
{
    Nodes base;
    for ( ClaimNode* cn : claims ) for ( pair<Node*, int32_t> no : cn->offs_ ) base += no.first;
    
    for ( ClaimNode* cn : claims )
    {
        unordered_set<ClaimNode*> pathed;
        cn->setPairs( cn, pathed, 0, false, false );
        assert( pathed.empty() );
    }
    
    for ( int d : { 0, 1 } ) for ( ClaimNode* cn : claims ) if ( cn->edges_[d].size() > 1 || !cn->joins_[d].empty() )
    {
        cn->getFork( d )->verifyFork( params.maxPeMean, false, d );
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

//void ClaimNode::setBlocked( Node* node, bool drxn )
//{
//    if ( offs_.find( node ) == offs_.end() ) return;
//    for ( Node* fork : getForks( node ) ) offs_.erase( fork );
//    for ( Edge& e : getEdges( node, drxn ) ) setBlocked( e.node, drxn );
//}

void ClaimNode::setBranched( vector<ClaimNode*>& claims, bool reversed, bool orient, bool drxn )
{
    if ( !edged_[!drxn] ) edge( claims, orient, !drxn );
    
    ClaimNode* cn = isBridge() ? edges_[drxn][0].node : this;
    if ( !cn->edged_[!drxn] ) cn->edge( claims, orient, !drxn );
    
    if ( !reversed ) for ( ClaimEdge& ce : cn->edges_[!drxn] ) ce.node->setBranched( claims, true, orient, !drxn );
}

bool ClaimNode::setBranched( vector<ClaimNode*>& claims, unordered_set<ClaimNode*>& tested, bool branched, bool orient, bool drxn )
{
    if ( !tested.insert( this ).second ) return false;
    
    // branched = path has already branched, forked = path branches forward
    if ( !branched && !edged_[drxn] ) edge( claims, orient, drxn );
    if ( edges_[drxn].size() > 1 ) setForked( claims, orient, drxn );
    
    bool forked = edges_[drxn].size() > 1;
    for ( ClaimEdge& ce : edges_[drxn] ) if ( ce.node->setBranched( claims, tested, branched || forked, orient, drxn ) ) forked = true;
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
    
    for ( ClaimNode* cn : claims ) if ( !cn->edged_[drxn] )
    {
        NodeDists dists;
        int32_t* dist;
        for ( pair<Node*, int32_t>& hang : cn->getHangs( orient, drxn ) ) dists.fillNot( hang.first, base, hang.second, limit+500, orient, drxn, true );
        
        for ( ClaimNode* alt : claims ) if ( cn != alt && !alt->edged_[!drxn] )
        {
            for ( pair<Node*, int32_t>& hang : alt->getHangs( orient, !drxn ) )
            {
                if ( cn->findFork( hang.first, drxn ) )
                {
                    cn->addJoin( alt, hang.first, dists, drxn ? hang.second - cn->ends_[1] : cn->ends_[0] - hang.second, drxn );
                }
                else if ( dist = dists.get( hang.first ) ) cn->addJoin( alt, hang.first, dists, drxn ? hang.second - *dist : *dist - hang.second, drxn );
            }
        }
        
        if ( !dists.empty() && !cn->edges_[drxn].empty() ) for ( pair<ClaimNode*, Nodes>& ext : ends )
        {
            ClaimShared cs( cn );
            for ( Node* node : ext.second.nodes ) if ( dists.find( node ) ) cs.shared.fillIn( node, ext.second, drxn, true );
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

void ClaimNode::setPairs( ClaimNode* cn, unordered_set<ClaimNode*>& pathed, int32_t diff, bool lFork, bool rFork )
{
    if ( cn->ends_[0] - cn->path_[0]->size() - diff - ends_[1] > params.maxMpMean ) return;
    if ( !pathed.insert( cn ).second ) return;
    
    if ( !lFork ) for ( ClaimEdge& ce : cn->edges_[1] ) if ( ce.node->isForked( 0 ) ) lFork = true;
    if ( !joins_[1].empty() ) lFork = true;
    
    bool added = false, failed = false;
    for ( ClaimPairing* cp : pairs_[1] ) if ( !added && ( added = ( cp->node[1] == cn ) ) )
    {
        for ( int i = 0; i <= cp->diffs.size(); i++ )
        {
            if ( failed = ( i < cp->diffs.size() && diff == cp->diffs[i] ) ) break;
            if ( i < cp->diffs.size() && cp->diffs[i] < diff ) continue;
//            if ( !i && ( diff < cp->diffs.back() - int32_t( params.maxMpMean / max( 1, (int)cp->diffs.size() - 1 ) ) ) ) return;
            cp->diffs.insert( cp->diffs.begin() + i, diff );
            cp->missed.insert( cp->missed.begin() + i, 0 );
            break;
        }
        break;
    }
    if ( !added && lFork && rFork ) new ClaimPairing( this, cn, diff );
    
    if ( cn->edges_[1].size() > 1 ) rFork = true;
    for ( ClaimEdge& ce : cn->edges_[1] ) for ( ClaimEdge& re : ce.node->edges_[0] ) if ( re.node->isForked( 1 ) ) rFork = true;
    if ( !joins_[1].empty() ) rFork = true;
    
    if ( !failed ) for ( ClaimEdge& ce : cn->edges_[1] ) setPairs( ce.node, pathed, diff+ce.diff, lFork, lFork && rFork );
    if ( !failed ) for ( ClaimJoin* cj : cn->joins_[1] ) setPairs( cj->node[1], pathed, diff+cj->diff( 1 ), lFork, lFork && rFork );
    pathed.erase( cn );
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
                ClaimScore cs( no.first, np.first, cp, *off - no.second, est );
                cs.score( claims, est, np.second.margin(), np.second.count );
                break;
            }
        }
    }
    
    for ( int i = 0; i < pairs_[1].size(); i++ ) if ( !pairs_[1][i]->hits ) delete pairs_[1][i--];
}

//bool ClaimNode::setState( bool drxn )
//{
//    for ( int i = 0; i < path_.size(); i++ )
//    {
//        Node* node = ( drxn ? path_[i] : path_.end()[-i-1] );
//        if ( node->setState() ) return true;
//    }
//    return false;
//}

bool ClaimNode::split( vector<ClaimNode*>& claims, bool orient, bool drxn )
{
    if ( !edged_[!drxn] )
    {
        if ( !getFork( !drxn )->isForkComplete( params.shortLen(), 20, !drxn ) ) return false;
        edge( claims, orient, !drxn );
        for ( ClaimEdge& ce : edges_[!drxn] )
        {
            if ( !ce.node->getFork( drxn )->isForkComplete( params.shortLen(), 20, drxn ) ) return false;
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

bool ClaimNode::trim( unordered_set<ClaimNode*> trimmed[2], NodeRoll& nodes, Nodes tried[2], int drxn )
{
    for ( int d : { 0, 1 } ) if ( drxn == 2 || d == drxn )
    {
        if ( !trimmed[d].insert( this ).second ) return false;
        for ( ClaimEdge& ce : edges_[d] ) if ( ce.node->trim( trimmed, nodes, tried, d ) ) return true;
    }
    
    if ( tried[1].find( path_[0] ) && tried[1].find( path_.back() ) ) return false;
    
    for ( int d : { 0, 1 } ) if ( edges_[d].size() > 1 && ClaimTrim::trim( this, nodes, d ) ) return true;
    
    if ( drxn == 2 ) for ( int d : { 0, 1 } ) for ( ClaimNode* cn : trimmed[d] )
    {
        bool complete = cn->edged_[0] && cn->edged_[1];
        for ( int i : { 0, 1 } ) for ( ClaimEdge& ce : edges_[d] ) if ( !ce.node->edged_[!i] ) complete = false;
        if ( complete ) for ( Node* node : cn->path_ ) tried[1] += node;
    }
    
    return false;
}
