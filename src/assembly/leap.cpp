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

#include "leap.h"
#include "shared_functions.h"
#include <algorithm>

LeapHit::LeapHit( Node* node, ReadId id, int32_t coord, int mp )
: node ( node )
{
    hit[mp].push_back( make_pair( id, coord ) );
}

bool LeapHit::add( Node* q, vector<LeapHit> hits[2], NodeRoll& nodes, ReadId ids[2], int32_t coords[2], int mp )
{
    Node* tar = NULL;
    for ( Node* node : nodes.nodes ) if ( tar = ( node->reads_.find( ids[1] ) != node->reads_.end() ? node : NULL ) ) break;
    if ( !tar ) return false;
    if ( tar->bad_ )
    {
        add( hits[0], q, ids[0], coords[0], mp );
        add( hits[1], tar, ids[1], coords[1], mp );
    }
    return true;
}

void LeapHit::add( vector<LeapHit>& hits, Node* node, ReadId id, int32_t coord, int mp )
{
    for ( LeapHit& lh : hits ) if ( lh.node == node )
    {
        for ( pair<ReadId, int32_t>& read : lh.hit[mp] ) if ( read.first == id ) return;
        lh.hit[mp].push_back( make_pair( id, coord ) );
        return;
    }
    hits.push_back( LeapHit( node, id, coord, mp ) );
}

void LeapHit::affirm( vector<LeapHit>& hits, NodeRoll& nodes )
{
    vector< pair<ReadId, int32_t> > dropped[2];
    for ( LeapHit& lh : hits ) for ( int mp : { 0, 1 } ) for ( int i = 0; i < lh.hit[mp].size(); i++ )
    {
        if ( lh.node->reads_.find( lh.hit[mp][i].first ) == lh.node->reads_.end() )
        {
            dropped[mp].push_back( lh.hit[mp][i] );
            lh.hit[mp].erase( lh.hit[mp].begin() + i-- );
        }
    }
    for ( int i = 0; i < hits.size(); i++ ) if ( hits[i].hit[0].empty() && hits[i].hit[1].empty() ) hits.erase( hits.begin() + i-- );
    for ( int mp : { 0, 1 } ) for ( pair<ReadId, int32_t>& read : dropped[mp] ) for ( Node* node : nodes.nodes )
    {
        auto it = node->reads_.find( read.first );
        if ( it != node->reads_.end() ) LeapHit::add( hits, node, read.first, read.second, mp );
    }
}

LeapMark::LeapMark( Node* node, string seq, ReadId qId, ReadId tId, int32_t qCoord, int32_t tCoord, bool mp )
: node( node ), seq( seq ), mp( mp ), redundant( false )
{
    ids[0] = qId;
    ids[1] = tId;
    coords[0] = qCoord;
    coords[1] = tCoord;
    ols[0] = ols[1] = ols[2] = joined = 0;
}

bool LeapMark::add( NodeRoll& nodes )
{
    bool added = false;
    for ( Node* node : nodes.nodes ) if ( node->add( ids[1], seq ) ) added = true;
    return added;
}

void LeapMark::join( Node* node, bool drxn )
{
    int32_t coords[2];
    if ( !mapSeqEnd( seq, node->seq_, 40, coords, drxn ) ) return;
    joined = max( joined, coords[1] - coords[0] );
    ols[drxn] = max( ols[drxn], coords[1] - coords[0] );
    assert( joined < seq.size() );
}

void LeapMark::overlap( LeapMark* lm )
{
    if ( seq == lm->seq ) ols[2] = lm->ols[2] = seq.size();
    else if ( seq.size() < lm->seq.size() && lm->seq.find( seq ) != string::npos )
    {
        if ( !lm->mp ) redundant = true;
        lm->ols[2] = max( lm->ols[2], (int)seq.size() );
    }
    else if ( lm->seq.size() < seq.size() && seq.find( lm->seq ) != string::npos )
    {
        if ( !mp ) lm->redundant = true;
        ols[2] = max( ols[2], (int)lm->seq.size() );
    }
    else
    {
        int ol[2]{ mapSeqOverlap( lm->seq, seq, 40 ), mapSeqOverlap( seq, lm->seq, 40 ) };
        for ( int d : { 0, 1 } ) if ( ol[d] )
        {
            ols[d] = max( ols[d], ol[d] );
            lm->ols[!d] = max( lm->ols[!d], ol[d] );
        }
    }
}

bool LeapMark::seed( Querier& bwt )
{
    if ( mp || redundant ) return false;
    if ( !ols[0] && !ols[1] && !ols[2] ) return false;
    if ( ols[0] && ols[1] ) return true;
    bool good = true;
    for ( int d : { 0, 1 } ) if ( !ols[d] && !bwt.isExtendable( seq, d ) ) good = false;
    return good;
}

LeapEnd::LeapEnd( LeapHit& lh, vector<LeapHit>& hits, bool orient, bool drxn )
: fork( lh.node ), dists( lh.node, params.maxMpMean, orient, !drxn, true ), score( 0 )
{
    for ( int mp : { 0, 1 } ) for ( pair<ReadId, int32_t>& read : lh.hit[mp] )
    {
        auto it = fork->reads_.find( read.first );
        assert( it != fork->reads_.end() );
        add( read.first, read.second - abs( fork->ends_[orient] - it->second[orient] ), mp );
        base.insert( read.first );
    }
//    x;
//    estCount[0] = estCount[1] = estSum[0] = estSum[1] = 0;
//    int32_t* dist;
//    for ( LeapHit& alt : hits ) if ( dist = dists.get( alt.node ) ) for ( int d : { 0, 1 } )
//    {
//        score += alt.hit[d].size();
//        estCount[d] += alt.hit[d].size();
//        for ( pair<ReadId, int32_t>& read : alt.hit[d] ) estSum[d] += read.second - *dist;
//    }
}

bool LeapEnd::add( ReadId id, int32_t coord, bool mp )
{
    if ( !reads.insert( id ).second ) return false;
    estimates[mp].push_back( make_pair( id, coord ) );
    sort( estimates[mp].begin(), estimates[mp].end(), []( pair<ReadId, int32_t>& a, pair<ReadId, int32_t>& b ){ return a.second < b.second; } );
    score++;
    return true;
}

void LeapEnd::add( LeapEnd* le, int32_t off )
{
    for ( int mp : { 0, 1 } ) for ( auto read : le->estimates[mp] ) if ( le->base.find( read.first ) != le->base.end() )
    {
        add( read.first, read.second + off, mp );
    }
}

int32_t LeapEnd::estimate()
{
    for ( int d : { 0, 1 } ) if ( !estimates[d].empty() )
    {
        int i = estimates[d].size() / 2, j = ( estimates[d].size()+1 ) % 2;
        return ( estimates[d][i-j].second + estimates[d][i].second ) / 2;
    }
    return 0;
}

int LeapEnd::getBranchReads( Node* node, int readCount, int readLimit, bool drxn )
{
    readCount += node->countReads( true );
    if ( readCount > readLimit ) return 0;
    int maxReads = readCount;
    for ( Edge& e : node->edges_[drxn] )
    {
        int branchReads = getBranchReads( e.node, readCount, readLimit, drxn );
        if ( !branchReads ) return 0;
        maxReads = max( maxReads, branchReads );
    }
    return maxReads;
}

NodeDists LeapEnd::getForks( bool drxn )
{
    NodeDists forks;
    int branchReads = getBranchReads( fork, 0, 5, drxn );
    if ( branchReads ) setForks( fork, forks, 5 - branchReads + fork->countReads( true ), !drxn );
    return forks;
}

bool LeapEnd::obsolete( LeapEnd* le )
{
    for ( auto& read : reads ) if ( le->reads.find( read ) == le->reads.end() ) return false;
    return true;
}

bool LeapEnd::retract( LeapEnd* le, bool drxn )
{
    NodeDists forks = getForks( drxn );
    return forks.get( le->fork );
}

void LeapEnd::setForks( Node* node, NodeDists& forks, int readsLeft, bool drxn )
{
    readsLeft -= node->countReads( true );
    int32_t* dist = dists.get( node );
    assert( dist );
    if ( !dist ) return;
    forks.add( node, *dist, drxn );
    if ( readsLeft >= 0 ) for ( Edge& e: node->edges_[drxn] ) setForks( e.node, forks, readsLeft, drxn );
}

LeapEnds::LeapEnds( vector<LeapHit>& hits, NodeRoll& nodes, bool consolidate, bool orient, bool drxn )
{
    LeapHit::affirm( hits, nodes );
    for ( LeapHit& lh : hits ) ends.push_back( new LeapEnd( lh, hits, orient, drxn ) );
    
    unordered_map<LeapEnd*, unordered_set<LeapEnd*> > forks[2];
    for ( LeapEnd* le : ends ) for ( int d : { 0, 1 } ) forks[d].insert( make_pair( le, unordered_set<LeapEnd*>() ) );
    for ( int i = 0; i < ends.size(); i++ ) for ( int j = i+1; j < ends.size(); j++ )
    {
        int32_t* dist[2]{ ends[i]->dists.get( ends[j]->fork ), ends[j]->dists.get( ends[i]->fork ) };
        for ( int d : { 0, 1 } ) if ( dist[d] )
        {
            int32_t ests[2]{ ends[d?j:i]->estimate(), ends[d?i:j]->estimate() };
            int32_t est = drxn == orient ? ests[1] - ests[0] : ests[0] - ests[1];
            if ( abs( *dist[d] ) > est + 300 ) dist[d] = NULL;
        }
        if ( dist[0] && dist[1] )
        {
            assert( false );
            int32_t ests[2]{ ends[i]->estimate(), ends[j]->estimate() };
            int32_t diffs[2]{ abs( *dist[0] - ( ests[0] - ests[1] ) ), abs( *dist[1] - ( ests[1] - ests[0] ) ) };
            assert( drxn && orient );
            for ( int d : { 0, 1 } ) if ( dist[d] && diffs[d] < 300 && diffs[d] < diffs[!d] - 300 ) dist[!d] = NULL;
            if ( dist[0] && dist[1] ) assert( false );
            if ( dist[0] && dist[1] ) dist[0] = dist[1] = NULL;
        }
        for ( int d : { 0, 1 } ) if ( dist[d] )
        {
            int b = d ? i : j, f = d ? j : i;
            auto iti = forks[d].find( ends[i] );
            auto itj = forks[!d].find( ends[j] );
            assert( iti != forks[d].end() && itj != forks[!d].end() );
            iti->second.insert( ends[j] );
            itj->second.insert( ends[i] );
            int32_t off = drxn == orient ? -abs( *dist[d] ) : abs( *dist[d] );
            ends[f]->add( ends[b], off );
            if ( ends[f]->retract( ends[b], drxn ) ) ends[b]->add( ends[f], -off );
        }
    }
    sort( ends.begin(), ends.end(), []( LeapEnd* a, LeapEnd* b ){ return a->score > b->score; } );
    
    unordered_set<LeapEnd*> culled;
    for ( auto& fork : forks[0] ) for ( LeapEnd* le : fork.second ) if ( fork.first->obsolete( le ) && culled.insert( fork.first ).second ) break;
    cull( culled, forks );
    if ( !consolidate ) return;
    for ( auto& fork : forks[1] ) for ( LeapEnd* le : fork.second ) if ( fork.first->obsolete( le ) && culled.insert( fork.first ).second ) break;
    cull( culled, forks );
    for ( LeapEnd* le : ends ) if ( le->score < 2 || le->estimates[0].empty() ) culled.insert( le );
    cull( culled, forks );
}

LeapEnds::~LeapEnds()
{
    for ( LeapEnd* le : ends ) delete le;
}

void LeapEnds::cull( unordered_set<LeapEnd*>& culled, unordered_map<LeapEnd*, unordered_set<LeapEnd*> > forks[2] )
{
    for ( int i : { 0, 1 } ) for ( auto it = forks[i].begin(); it != forks[i].end(); )
    {
        if ( culled.find( it->first ) != culled.end() ) it = forks[i].erase( it );
        else
        {
            for ( auto it2 = it->second.begin(); it2 != it->second.end(); )
            {
                if ( culled.find( *it2 ) != culled.end() ) it2 = it->second.erase( it2 );
                else it2++;
            }
            it++;
        }
    }
    for ( int i = 0; i < ends.size(); i++ ) if ( culled.find( ends[i] ) != culled.end() )
    {
        delete ends[i];
        ends.erase( ends.begin() + i-- );
    }
    culled.clear();
}

vector< pair<Node*, int> > LeapEnds::get( bool orient, bool drxn )
{
    unordered_map<Node*, int> scores;
    for ( LeapEnd* le : ends )
    {
        NodeDists dists[2];
        dists[0].fillReads( le->fork, 0, 0, 5, orient, !drxn );
        for ( auto& no : dists[0].map ) dists[1].fill( no.first, no.second, max( 200, abs( le->estimate() ) ), orient, drxn, true );
        for ( auto& no : dists[1].map )
        {
            auto it = scores.insert( make_pair( no.first, le->score ) );
            if ( !it.second ) it.first->second = max( it.first->second, le->score );
        }
    }
    vector< pair<Node*, int> > branches;
    for ( const pair<Node*, int>& score : scores ) branches.push_back( score );
    return branches;
}

Leap::Leap( Querier& bwt, NodeRoll& nodes, unordered_map<Node*, int32_t>& offs, int32_t coord, bool drxn )
: coord_( coord ), seeded_( false )
{
    Lib* lib;
    for ( const pair<Node*, int32_t>& no : offs )
    {
        if ( no.first->cloned_ ) for ( Node* clone : no.first->cloned_->nodes ) clone->setAllPairs();
        if ( no.second < params.maxPeMax ) for ( auto& read : no.first->reads_ ) if ( ( lib = params.getLib( read.first ) ) && lib->isPe )
        {
            ReadId ids[2]{ read.first, read.first };
            if ( lib->getPair( ids[1] ) != drxn ) continue;
            int32_t coords[2]{ no.second + abs( no.first->ends_[drxn] - read.second[drxn] )
                             , no.second - lib->size + abs( no.first->ends_[drxn] - read.second[!drxn] ) };
            bool ignore = no.first->reads_.find( ids[1] ) != no.first->reads_.end();
            for ( Node* clone : no.first->clones() ) for ( auto& np : clone->hits_.pairs[drxn] ) if ( np.first->reads_.find( ids[1] ) != no.first->reads_.end() ) ignore = true;
            if ( ignore || LeapHit::add( no.first, hits_, nodes, ids, coords, 0 ) ) continue;
            marks_.push_back( new LeapMark( no.first, bwt.getSequence( ids[1] ), ids[0], ids[1], coords[0], coords[1], false ) );
        }
        if ( no.second < params.maxMpMean ) for ( NodeMark& nm : no.first->mp_[drxn] )
        {
            ReadId ids[2]{ params.getPairId( nm.id ), nm.id };
            int32_t coords[2]{ no.second + abs( no.first->ends_[drxn] - nm.coords[drxn] )
                             , no.second - nm.dist + abs( no.first->ends_[drxn] - nm.coords[!drxn] ) };
            if ( coords[1] < 0 || LeapHit::add( no.first, hits_, nodes, ids, coords, 1 ) ) continue;
            marks_.push_back( new LeapMark( no.first, bwt.getSequence( ids[1] ), ids[0], ids[1], coords[0], coords[1], true ) );
        }
    }
    
    add( nodes );
    for ( int i = 0; i+1 < marks_.size(); i++ ) for ( int j = i+1; j < marks_.size(); j++ ) marks_[i]->overlap( marks_[j] );
    for ( auto& no : offs ) if ( no.second < 100 ) for ( LeapMark* lm : marks_ ) lm->join( no.first, !drxn );
}

Leap::~Leap()
{
    for ( LeapMark* lm : marks_ ) delete lm;
}

void Leap::add( NodeRoll& nodes )
{
    for ( int i = 0; i < marks_.size(); i++ ) if ( marks_[i]->add( nodes ) )
    {
        LeapHit::add( marks_[i]->node, hits_, nodes, marks_[i]->ids, marks_[i]->coords, marks_[i]->mp );
        delete marks_[i];
        marks_.erase( marks_.begin() + i-- );
    }
}

Node* Leap::bridge( NodeRoll& nodes, bool drxn )
{
    LeapEnds ends[2]{ LeapEnds( hits_[0], nodes, true, drxn, drxn ), LeapEnds( hits_[1], nodes, true, drxn, !drxn ) };
    vector< pair<Node*, int> > branches[2]{ ends[0].get( drxn, drxn ), ends[1].get( drxn, !drxn ) };
    Node* best[2]{ NULL, NULL };
    int maxScore = 0, maxOl = 0;
    for ( pair<Node*, int>& b : branches[0] ) if ( !b.first->bad_ ) for ( pair<Node*, int>& f : branches[1] )
    {
        int ol = mapSeqOverlap( ( drxn ? b.first : f.first )->seq_, ( drxn ? f.first : b.first )->seq_, 25 );
        if ( b.second + f.second + ( ol / 10 ) + ( ol > maxOl ) <= maxScore + ( maxOl / 10 )  ) continue;
        best[0] = b.first;
        best[1] = f.first;
        maxOl = ol;
        maxScore = b.second + f.second;
    }
    assert( false );
    if ( !maxOl ) return NULL;
    best[0]->addEdge( best[1], maxOl, drxn );
    best[1]->setState();
    for ( Node* node : nodes.nodes ) if ( node->ends_.ends[0] == -20324 && node->ends_[1] == -19611 ) node->setVerified();
    Node::verify( nodes );
    return best[1];
}

void Leap::confirmBranch( Node* node, unordered_map<Node*, int32_t>& offs, Nodes& branch, bool drxn )
{
    if ( offs.find( node ) == offs.end() || !branch.add( node ) ) return;
    for ( Edge& e : node->edges_[drxn] ) confirmBranch( e.node, offs, branch, drxn );
}

void Leap::extend( Querier& bwt, NodeRoll& nodes, bool drxn )
{
    for ( LeapEnd* le : LeapEnds( hits_[1], nodes, false, drxn, !drxn ).ends ) if ( le->score > 1 )
    {
        le->fork->extendFork( bwt, nodes, params.maxPeMean*2, 16, !drxn );
        le->fork->extendFork( bwt, nodes, params.maxPeMean*3, 8, drxn );
    }
}

void Leap::extend( Querier& bwt, NodeRoll& nodes, Node* node, bool drxn )
{
    int counts[2]{ nodes.size(), nodes.size() };
    node->extendFork( bwt, nodes, params.maxPeMean, 8, !drxn );
    node->extendFork( bwt, nodes, params.maxPeMean, 8, drxn );
    counts[1] = nodes.size();
}

void Leap::getBranch( Node* node, unordered_map<Node*, int32_t>& offs, int32_t dist, int32_t limit, bool neg, bool drxn )
{
    auto it = offs.insert( make_pair( node, dist ) );
    if ( !it.second )
    {
        if ( dist < 0 && it.first->second > 0 ) assert( false );
        if ( dist < 0 && it.first->second > 0 ) return;
        bool update = abs( dist ) < abs( it.first->second );
        if ( !update && dist ) return;
        if ( update ) it.first->second = dist;
    }
    if ( abs( dist ) < limit ) for ( Edge& e : node->edges_[drxn] )
    {
        int32_t coord = dist + ( neg ? e.ol - e.node->size() : node->size() - e.ol );
        if ( neg || coord < limit ) getBranch( e.node, offs, coord, limit, neg, drxn );
    }
}

Node* Leap::join( Querier& bwt, NodeRoll& nodes, bool drxn )
{
    assert( joins_.empty() );
    LeapEnds ends[2]{ LeapEnds( hits_[0], nodes, true, drxn, drxn ), LeapEnds( hits_[1], nodes, true, drxn, !drxn ) };
    for ( LeapEnd* le : ends[1].ends ) if ( le->score > 1 && !le->fork->isBlunt( 0, 3, drxn ) )
    {
        Nodes island( le->fork, params.maxPeMean + 200, !drxn, true );
        for ( LeapEnd* re : ends[0].ends ) if ( island.find( re->fork ) )
        {
            return le->fork;
        }
    }
    
    return NULL;
}

bool Leap::leapBranch( Querier& bwt, NodeRoll& nodes, Node* node, bool drxn )
{
    while ( !node->bad_ && node->pruneBranch( bwt, nodes, 80, drxn ) );
    if ( node->bad_ ) return false;
    unordered_map<Node*, int32_t> offs;
    getBranch( node, offs, 0, params.maxPeMean, false, drxn );
    for ( Edge& e : node->edges_[!drxn] ) if ( !e.node->bad_ ) getBranch( e.node, offs, e.ol - e.node->size(), params.maxPeMean, true, !drxn );
    Nodes branch;
    for ( auto& no : offs ) if ( no.second >= 0 )
    {
        no.first->setAllPairs();
        for ( auto& np : no.first->hits_.pairs[!drxn] ) confirmBranch( np.first, offs, branch, drxn );
    }
    for ( auto it = offs.begin(); it != offs.end(); )
    {
        if ( it->second >= 0 || branch.find( it->first ) ) it++;
        else it = offs.erase( it );
    }
    
    Leap l( bwt, nodes, offs, node->ends_[!drxn], !drxn );
    bool lept = l.leap( bwt, nodes, !drxn );
    while ( !node->bad_ && node->pruneBranch( bwt, nodes, 80, drxn ) );
    
    return lept;
}

Node* Leap::leapEnd( Querier& bwt, NodeRoll& nodes, Node* node, unordered_map<Node*, int32_t>& offs, bool drxn )
{
    Leap l( bwt, nodes, offs, node->ends_[drxn], drxn );
    return l.leap( bwt, nodes, drxn );
}

Node* Leap::leap( Querier& bwt, NodeRoll& nodes, bool drxn )
{
    Node* bridged = NULL;
    extend( bwt, nodes, drxn );
    if ( bridged = bridge( nodes, drxn ) ) return bridged;
    if ( !seed( bwt, nodes, drxn ) ) return bridged;
    if ( bridged = join( bwt, nodes, drxn ) ) return bridged;
    extend( bwt, nodes, drxn );
    remap( bwt, nodes, drxn );
    if ( bridged = bridge( nodes, drxn ) ) return bridged;
    assert( false );
    return bridge( nodes, drxn );
}

void Leap::remap( Querier& bwt, NodeRoll& nodes, bool drxn )
{
//    LeapEnds ends( hits_[1], nodes, drxn, !drxn );
//    Nodes tar;
//    for ( LeapEnd* le : ends.ends ) for ( auto& no : le->offs.map ) tar += no.first;
//    bool remapped = false;
////    for ( Node* node : tar.nodes ) if ( node->remap( bwt ) ) remapped = true;
//    for ( int mp : { 0, 1 } ) for ( LeapMark& lm : marks_[mp] )
//    {
//        if ( LeapHit::add( lm.node, hits_, nodes, lm.ids, lm.coords, mp ) )
//        {
//            int x = 0;
//        }
//    }
}

bool Leap::seed( Querier& bwt, NodeRoll& nodes, bool drxn )
{
    sort( marks_.begin(), marks_.end(), []( LeapMark* a, LeapMark* b ){ return a->ols[0] + a->ols[1] > b->ols[0] + b->ols[1]; } );
    Nodes seeds;
    for ( int i = 0; i < marks_.size(); i++ )
    {
        bool used = LeapHit::add( marks_[i]->node, hits_, nodes, marks_[i]->ids, marks_[i]->coords, marks_[i]->mp );
        if ( !used && marks_[i]->seed( bwt ) )
        {
            int32_t coord = coord_ + ( drxn ? -marks_[i]->coords[1] - marks_[i]->seq.size() : marks_[i]->coords[1] );
            vector<Node*> added, seeded;
            Node::seedNode( bwt, nodes, added, seeded, marks_[i]->seq, marks_[i]->ids[1], coord, drxn );
            for ( Node* node : seeded ) extend( bwt, nodes, node, drxn );
            for ( Node* node : added ) extend( bwt, nodes, node, drxn );
            for ( Node* node : seeded ) seeds += node;
            if ( !seeded.empty() ) LeapHit::add( marks_[i]->node, hits_, nodes, marks_[i]->ids, marks_[i]->coords, 0 );
            if ( !seeded.empty() || !added.empty() ) used = true;
        }
        if ( used )
        {
            delete marks_[i];
            marks_.erase( marks_.begin() + i-- );
        }
    }
    add( nodes );
    
    return !seeds.empty();
}

//LeapBranch::LeapBranch( Node* node, bool drxn )
//{
//    getBranch( offs_[0], node, 0, params.maxPeMean, !drxn, drxn );
//    getBranch( offs_[1], node, 0, params.maxPeMean, !drxn, drxn );
//}
//
//bool LeapBranch::claim( Querier& bwt, NodeRoll& nodes, Node* node, bool drxn )
//{
//    if ( node->verified_ ) return false;
//    
//    Nodes branch;
//    getBranch( branch, 15, drxn );
//    for ( Node* node : branch.nodes ) if ( node != this )
//    {
//        if ( node->pruneFork( bwt, nodes, !drxn ) ) assert( false );
////        if ( node->pruneFork( bwt, nodes, !drxn ) ) return true;
//    }
//    
//    return false;
//}
//
//void LeapBranch::getBranch( Nodes& branch, Node* node, int readsLeft, bool drxn )
//{
//    branch += node;
//    readsLeft -= node->countReads( true );
//    if ( readsLeft >= 0 ) for ( Edge& e : node->edges_[drxn] ) getBranch( branch, e.node, readsLeft, drxn );
//}
//
//void LeapBranch::getBranch( unordered_map<Node*, int32_t>& offs, Node* node, int32_t dist, int32_t limit, bool orient, bool drxn )
//{
//    auto it = offs.insert( make_pair( node, dist ) );
//    if ( !it.second )
//    {
//        if ( it.first->second <= dist ) return;
//        it.first->second = dist;
//    }
//    if ( dist < limit ) for ( Edge& e : node->edges_[drxn] )
//    {
//        getBranch( offs, e.node, dist + ( drxn == orient ? e.node : node )->size() - e.ol, limit, orient, drxn );
//    }
//}
//
//bool LeapBranch::leap( Querier& bwt, NodeRoll& nodes, Node* node, bool drxn )
//{
//    if ( claim( bwt, nodes, node, drxn ) ) return true;
//    
//    LeapBranch();
//    
//    return false;
//}
