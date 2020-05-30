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
#include <algorithm>

void Node::add( Node* tar, NodeMark& mark, bool drxn )
{
    Coords* hit;
    tar->findRead( mark.id, hit );
    add( tar, mark, hit, drxn );
}

void Node::add( Node* tar, NodeMark& mark, Coords* hit, bool drxn )
{
    if ( hit->ignore ) return;
    
    ReadId id[2]{ params.getPairId( mark.id ), mark.id };
    if ( tar != this || ( drxn ? (*hit)[0] < mark[0] : mark[1] < (*hit)[1] ) )
    {
        int32_t dists[2][2]{ { mark[!drxn] - ends_[0], ends_[1] - mark[!drxn] }, { (*hit)[drxn] - tar->ends_[0], tar->ends_[1] - (*hit)[drxn] } };
        int32_t len = max( mark.dist, max( dists[0][drxn], dists[1][!drxn] ) );
        int32_t ests[2]{ dists[1][drxn] - dists[0][drxn] + len, dists[0][!drxn] - dists[1][!drxn] + len };
        bool mp = params.isReadMp( mark.id );
        hits_.add( tar, id[0], ests[0], mark.dist, mp, drxn );
        tar->hits_.add( this, id[1], ests[1], mark.dist, mp, !drxn );
    }
    else hits_.count++;
    tar->rmvMark( params.getPairId( mark.id ), !drxn );
}

bool Node::add( NodeRoll& nodes, string& seq, NodeMark& mark, bool drxn )
{
    for ( Node* node : nodes.nodes ) if ( node->reads_.find( mark.id ) != node->reads_.end() ) return true;
    for ( Node* node : nodes.nodes )
    {
        size_t it = node->seq_.find( seq );
        if ( it == node->seq_.npos ) continue;
        int32_t coords[2]{ int32_t( it + node->ends_[0] ), int32_t( it + node->ends_[0] + seq.size() ) };
        node->add( mark.id, coords[0], coords[1], node->isRedundant( coords[0], coords[1] ) );
        add( node, mark, drxn );
        assert( false );
        return true;
    }
    return false;
}

bool Node::bridgeVerified( Nodes& verified, bool drxn )
{
    if ( verified.find( this ) ) return true;
    
    bool bridged = false;
    for ( Edge& e : edges_[drxn] ) if ( e.node->bridgeVerified( verified, drxn ) ) bridged = true;
    if ( bridged ) verified += this;
    return bridged;
}

bool Node::canSetVerified()
{
    if ( isContinue( 0 ) || isContinue( 1 ) ) return false;
    return !bad_ && ends_.verified();
}

bool Node::canVerify()
{
    if ( drxn_ >= 2 || verified_ ) return true;
    if ( bad_ ) return false;
    for ( Edge& e : edges_[!drxn_] ) if ( !e.node->bad_ && e.node->verified_ ) return true;
    return false;
}

void Node::confirmVerified()
{
    if ( drxn_ >= 2 || !verified_ ) return;
    if ( !bad_ ) for ( Edge& e : edges_[!drxn_] ) if ( e.node->verified_ ) return;
    verified_ = false;
    reverify();
    for ( Edge& e : edges_[drxn_] ) e.node->confirmVerified();
}

void Node::getVerified( Nodes& verified, bool bridge, bool drxn )
{
    if ( !verified_ || !verified.add( this ) ) return;
    
    for ( int d : { 0, 1 } ) if ( d == drxn || drxn_ >= 2 )
    {
        for ( Edge& e : edges_[d] ) e.node->getVerified( verified, bridge, d );
        for ( Edge& e : edges_[!d] ) if ( !verified.find( this ) ) e.node->bridgeVerified( verified, !d );
    }
}

void Node::reverify( NodeRoll& nodes )
{
    updateStates( nodes );
    Nodes verified;
    for ( Node* node : nodes.nodes ) if ( node->drxn_ >= 2 ) node->getVerified( verified, true, 1 );
//    Nodes unverified[2];
    for ( Node* node : nodes.nodes )
    {
        if ( node->pairedNodes_ ) delete node->pairedNodes_;
        node->pairedNodes_ = NULL;
        node->branch_ = false;
        node->hits_.reset();
        node->remark();
        node->verified_ = verified.find( node );
    }
//    for ( int d : { 0, 1 } ) for ( Node* node : unverified[d].nodes ) assert( !node->verified_ );
//    for ( int d : { 0, 1 } ) for ( Node* node : unverified[d].nodes ) if ( node->drxn_ < 2 ) node->verified_ = false;
    for ( Node* node : nodes.nodes ) if ( !node->verified_ ) node->ends_.reset( node->drxn_ );
    for ( Node* node : nodes.nodes ) if ( node->verified_ ) node->setVerified();
    verify( nodes );
    NodeRoll ends[2];
    int count = 0;
    for ( Node* node : nodes.nodes ) if ( node->verified_ )
    {
        for ( int d : { 0, 1 } )
        {
            bool ended = !node->edges_[d].empty();
            for ( Edge& e : node->edges_[d] ) if ( e.node->verified_ ) ended = false;
            if ( ended ) ends[d] += node;
        }
        count++;
    }
    cout << "Reverified nodes: " << count << ", started with: " << verified.size() << endl;
    cout << "Left verified ends: " << ends[0].size() << ", right verified ends: " << ends[1].size() << endl;
}

bool Node::reverify()
{
    vector<Node*> test = clones();
    for ( Node* node : test )
    {
        node->branch_ = false;
        node->clearPaired( true );
        if ( node->bad_ ) node->verified_ = false;
        if ( !node->verified_ ) node->ends_.reset( node->drxn_ );
    }
    for ( Node* node : test )
    {
        if ( node->verified_ ) node->setVerified();
        else if ( !node->bad_ ) node->verify();
    }
    return verified_;
}

void Node::setAllPairs()
{
    if ( !verified_ ) branch_ = true;
    
    for ( bool d : { 0, 1 } )
    {
        if ( !cloned_ )
        {
            for ( Node* node : Nodes( this, params.maxPeMax, d, true ).nodes ) verify( node, pe_[d], 3, d );
            for ( Node* node : Nodes( this, params.maxMpMax, d, true ).nodes ) verify( node, mp_[d], 3, d );
            continue;
        }
        
        NodeRoll clones[2]{ NodeRoll::clones( this ), NodeRoll() };
        vector<NodeOffsets> offs;
        for ( Node* clone : clones[0].nodes ) offs.push_back( NodeOffsets( clone, params.maxMpMax, d, d, true ) );
        for ( int i = 0; i < 2; i++ )
        {
            int32_t limit = i ? params.maxMpMax : params.maxPeMax;
            Nodes fwd( this, limit, d, true );
            for ( Node* clone : cloned_->nodes ) fwd.fill( clone, limit, d, true );
            for ( Node* f : fwd.nodes )
            {
                clones[1] = NodeRoll::clones( f );
                verifyClone( clones, offs, i ? mp_[d] : pe_[d], 3, d );
            }
        }
    }
}

void Node::setUnverified()
{
    clearPaired( true );
    verified_ = false;
    ends_.reset( drxn_ );
}


bool Node::setVerified()
{
    ends_.init( ends_[0], 0 );
    ends_.init( ends_[1], 1 );
    if ( !bad_ ) verified_ = true;
    setAllPairs();
    if ( bad_ ||  isContinue( 0 ) || isContinue( 1 ) ) return verified_ = false;
    
    return true;
}

void Node::setVerifyLimits( int32_t limits[2] )
{
    int32_t ols[2]{ getBestOverlap( 0 ), getBestOverlap( 1 ) };
    if ( cloned_ )
    {
        for ( Node* clone : cloned_->nodes )
        {
            ols[0] = max( ols[0], clone->getBestOverlap( 0 ) );
            ols[1] = max( ols[1], clone->getBestOverlap( 1 ) );
        }
    }
    limits[0] = ends_[0] + ( ols[0] ? : params.readLen );
    limits[1] = ends_[1] - ( ols[1] ? : params.readLen );
}

int Node::verify( Node* tar, vector<NodeMark> &marks, int fwdHits, bool drxn )
{
    int added = 0;
    Coords* hit;
    if ( !culled_ ) cullMarks();
    int32_t limits[2];
    tar->setVerifyLimits( limits );
    
    for ( int i = 0; i < marks.size(); i++ )
    {
        if ( !bad_ && !branch_ && !tar->bad_ && !ends_.isValidMark( marks[i].coords[drxn], fwdHits, drxn ) ) return added;
        if ( !tar->findRead( marks[i].id, hit ) || hit->ignore ) continue;
        if ( cloned_ || tar->cloned_ ) return verifyClone( tar, marks, fwdHits, drxn );
        if ( limits[0] < hit->coords[1] && hit->coords[0] < limits[1] )
        {
            if ( drxn != drxn_ && !bad_ ) ends_.push( marks[i].coords[0], !drxn );
            if ( !bad_ && !tar->bad_ && !branch_ && !tar->branch_ ) tar->ends_.push( (*hit)[drxn], drxn );
            add( tar, marks[i], hit, drxn );
            if ( tar != this ) fwdHits++;
            added++;
        }
        marks.erase( marks.begin() + i-- );
    }
    setPaired( tar );
    return added;
}

void Node::verify( NodeRoll& nodes )
{
    updateStates( nodes );
    NodeRoll seed = nodes.getGraph( 2 );
    Nodes tested;
    for ( Node* node : seed.nodes ) if ( tested.add( node ) ) node->verify();
    for ( Node* node : seed.nodes )
    {
        if ( node->ends_.verified( 0 ) ) for ( Edge &e : node->edges_[0] ) e.node->verify( tested, 0 );
        if ( node->ends_.verified( 1 ) ) for ( Edge &e : node->edges_[1] ) e.node->verify( tested, 1 );
    }
}

void Node::verify( bool drxn )
{
    if ( !verify() ) return;
    Nodes tested( this );
    for ( Edge& e : edges_[drxn] ) e.node->verify( tested, drxn );
}

void Node::verify( Nodes& tested, bool drxn )
{
    if ( verified_ && tested.find( this ) ) return;
    if ( bad_ )
    {
        bad_ = false;
        drxn_ = drxn;
        setUnverified();
    }
    if ( !verify() && drxn_ < 2 ) return;
    tested += this;
    
    for ( Edge &e : edges_[drxn] ) e.node->verify( tested, drxn );
    for ( Edge &e : edges_[!drxn] ) if ( e.node->bad_ && tested.add( e.node ) ) e.node->verifyFork( params.maxPeMean*1.3+200, false, drxn );
}

bool Node::verify()
{
    if ( verified_ ) return true;
    bool proceed = drxn_ >= 2, d = drxn_;
    if ( !proceed ) for ( Edge &e : edges_[!d] ) if ( e.node->verified_ ) proceed = true;
    assert( !bad_ );
    if ( bad_ || !proceed ) return false;
    if ( branch_ && ( !hits_.empty() ) ) return reverify();
    branch_ = false;
    
    NodeRoll bck( Nodes( this, params.maxPeMax, !d, false ) ), fwd( this );
    bck += this;
    fwd += Nodes( this, params.maxPeMax, d, false );
    int hits[2]{0};
    for ( Node* b : bck.nodes ) if ( b->ends_.targetable( d ) )
    {
        b->resort();
        for ( Node* f : fwd.nodes )
        {
            if ( !b->isPaired( f ) ) b->verify( f, b->pe_[d], b == this ? hits[d] : 0, d );
            if ( b->bad_ ) continue;
            int pairHit = b->hits_.get( f, d );
            if ( b != this ) hits[!d] += pairHit;
            if ( f != this ) hits[d] += pairHit;
        }
    }
    
    if ( hits[d] > 2 && hits[!d] > 2 ) return setVerified();
    if ( hits[0] > 2 ) ends_.init( ends_[0], 0 );
    if ( hits[1] > 2 ) ends_.init( ends_[1], 1 );
    return false;
}

int Node::verifyClone( Node* tar, vector<NodeMark>& marks, int fwdHits, bool drxn )
{
    int added = 0;
    NodeRoll clones[2]{ NodeRoll( this ), NodeRoll( tar ) };
    if ( cloned_ ) for ( Node* clone : cloned_->nodes ) clones[0] += clone;
    if ( tar->cloned_ ) for ( Node* clone : tar->cloned_->nodes ) clones[1] += clone;
    vector<NodeOffsets> offs;
    for ( Node* clone : clones[0].nodes ) offs.push_back( NodeOffsets( clone, params.maxPeMean, drxn, drxn, true ) );
    
    return verifyClone( clones, offs, marks, fwdHits, drxn );
}

int Node::verifyClone( NodeRoll clones[2], vector<NodeOffsets>& offs, vector<NodeMark>& marks, int fwdHits, bool drxn )
{
    struct Match
    {
        Match( Node* a, Node* b, int32_t diff ): diff( diff ){ node[0] = a; node[1] = b; pref[0] = pref[1] = true; };
        Node* node[2];
        int32_t diff;
        bool pref[2];
    };
    
    int added = 0;
    Coords* hit;
    assert( !clones[0].empty() && !clones[1].empty() );
    Node* tar = clones[1].nodes[0];
    if ( !culled_ ) cullMarks();
    int32_t limits[2]{ tar->ends_[0] + tar->getBestOverlap( 0, true ), tar->ends_[1] - tar->getBestOverlap( 1, true ) };
    NodeOffset* off = NULL;
    
    for ( int i = 0; i < marks.size(); i++ ) if ( tar->findRead( marks[i].id, hit, true ) && !hit->ignore )
    {
        // Establish estimated distance between nodes and cutoff deviation from that estimate
        int32_t dists[2]{ abs( marks[i][!drxn] - ends_[drxn] ), abs( (*hit)[drxn] - tar->ends_[drxn] ) };
        int32_t est = ( dists[1] - dists[0] + marks[i].dist ) * ( drxn ? 1 : -1 );
        int32_t cutoffs[2]{ -marks[i].dist / ( params.isReadPe( marks[i].id ) ? 1 : 2 ), 200 + marks[i].dist / 5 };
        int32_t best = marks[i].dist*2;
        int32_t cutoff = 100 + marks[i].dist / 5;
        int32_t self = -1;
        
        // Record all pairs of clones that match the estimate within its cutoff
        vector<Match> matches;
        for ( int j = 0; j < clones[0].size(); j++ ) for ( int k = 0; k < clones[1].size(); k++ ) if ( off = offs[j].get( clones[1][k] ) )
        {
            Match match( clones[0][j], clones[1][k], off->diff( est, drxn ) );
            best = min( best, abs( match.diff ) );
            if ( match.node[0] == match.node[1] && abs( marks[i][drxn] - ends_[drxn] ) < abs( hit->coords[drxn] - tar->ends_[drxn] ) ) continue;
            if ( match.node[0] == match.node[1] && ( self < 0 || abs( match.diff ) < self ) ) self = abs( match.diff );
            if ( cutoffs[0] <= match.diff && match.diff <= cutoffs[1] ) matches.push_back( match );
        }
        cutoff += best;
        
        if ( hit->coords[1] <= limits[0] || limits[1] <= hit->coords[0] ) matches.clear();
        sort( matches.begin(), matches.end(), []( Match& a, Match& b ){ return abs( a.diff ) < abs( b.diff ); } );
        for ( int j = 0; j < matches.size(); j++ ) for ( int k = j+1; k < matches.size(); k++ ) for ( int d = 0; d < 2; d++ )
        {
            if ( matches[j].node[d] == matches[k].node[d] && abs( matches[j].diff ) + 20 < abs( matches[k].diff ) ) matches[k].pref[d] = false;
        }
        
        // Remove obviously inferior pairings
        if ( self >= 0 ) for ( int j = 0; j < matches.size(); j++ ) if ( matches[j].node[0] != matches[j].node[1] ) assert( abs( matches[j].diff ) <= self );
        if ( self >= 0 ) for ( int j = 0; j < matches.size(); j++ ) if ( matches[j].node[0] != matches[j].node[1] && self < abs( matches[j].diff ) ) matches.erase( matches.begin() + j-- );
        for ( int j = 1; j < matches.size(); j++ ) if ( cutoff < abs( matches[j].diff ) ) matches.erase( matches.begin() + j-- );
        for ( int j = 1; j < matches.size(); j++ ) if ( !matches[j].pref[0] && !matches[j].pref[1] ) matches.erase( matches.begin() + j-- );
        
        Nodes used[2];
        for ( Match& match : matches )
        {
            // Set hit
            NodeMark mark = marks[i];
            mark.offset( match.node[0]->ends_[0] - ends_[0] );
            match.node[0]->add( match.node[1], mark, drxn );
            
            // Push back limits
            int32_t diffs[2]{ match.node[0]->ends_[0] - ends_[0], match.node[1]->ends_[0] - tar->ends_[0] };
            if ( !match.node[0]->bad_ && !match.node[1]->bad_ && used[0].add( match.node[0] ) ) match.node[0]->ends_.push( mark[!drxn] + diffs[0], !drxn );
            if ( !match.node[0]->bad_ && !match.node[1]->bad_ && used[1].add( match.node[1] ) ) match.node[1]->ends_.push( (*hit)[drxn] + diffs[1], drxn );
            if ( match.node[0] == this ) added++;
        }
        ReadId ids[2]{ marks[i].id, params.getPairId( marks[i--].id ) };
        for ( int d : { 0, 1 } ) for ( int j = 0; j < clones[d].size(); j++ ) clones[d][j]->rmvMark( ids[d], d ? !drxn : drxn );
    }
    
    return added;
}

void Node::verifyFill()
{
    if ( verified_ || bad_ || drxn_ >= 2 ) return;
    setVerified();
    for ( Edge& e : edges_[!drxn_] ) e.node->verifyFill();
}

void Node::verifyFork( int32_t dist, bool force, bool drxn )
{
    NodeOffsets offs[2]{ NodeOffsets( this, dist - size(), drxn, !drxn, true ), NodeOffsets( this, dist, drxn, drxn, false ) };
    for ( const pair<Node*, NodeOffset>& b : offs[0].map )
    {
        if ( !b.first->verified_ && !b.first->bad_ ) b.first->branch_ = true;
        for ( const pair<Node*, NodeOffset>& f : offs[1].map )
        {
            if ( force || !b.first->isPaired( f.first ) ) b.first->verify( f.first, b.first->pe_[drxn], 0, drxn );
        }
    }
}
