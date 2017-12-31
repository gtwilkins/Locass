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

#include "path_review.h"
#include "path_reassembly.h"
#include "path_merge.h"
#include <algorithm>
#include <cassert>

extern Parameters params;

bool AltPath::doesConflict( AltPath &bck, bool drxn )
{
    return ( drxn ? limits[0] < bck.limits[1] : bck.limits[0] < limits[1] );
}

bool ConPath::doesConflict( AltPath &bck, bool drxn )
{
    int32_t limit = drxn ? max( paths[0][0]->ends_[0], paths[1][0]->ends_[0] ) 
                         : min( paths[0].back()->ends_[1], paths[1].back()->ends_[1] );
    return ( drxn ? limit < bck.limits[1] : bck.limits[1] < limit );
}

PathReview::PathReview( PathVars &pv, NodeList &path, int32_t* reliable, int32_t* forkLimits, bool anteFinished )
: pv_( pv ), fork_( NULL ), truncate_( NULL ), path_( path ), pathSet_( path.begin(), path.end() )
, drxn_( pv.drxn ), anteFinished_( anteFinished )
{
    path[0]->getDrxnNodes( pathSet_, !drxn_, params.getFurthestMpDist( path[0]->ends_[!drxn_], !drxn_ ) );
    NodeSet bckSet = { path_[0] };
    path_.back()->getDrxnNodes( bckSet, !drxn_, true );
    for ( Node* bck : bckSet )
    {
        if ( pathSet_.find( bck ) == pathSet_.end() ) convSet_.insert( bck );
    }
    
    reliLimits_[0] = reliable[0] + ( params.maxPeMean / 2 );
    reliLimits_[1] = reliable[1] - ( params.maxPeMean / 2 );
    forkLimits_[0] = forkLimits[0];
    forkLimits_[1] = forkLimits[1];
}

bool PathReview::review( Path &path, NodeList &sideNodes, NodeSet &delSet )
{
    bool assembled = resolveForks( delSet );
    assembled = resolveEnd( delSet ) && assembled;
    assembled = resolveMisassembly( delSet ) && assembled;
    if ( resolveUnspanned( path, delSet ) ) assembled = false;
    
    // Set new fork
    if ( fork_ && delSet.find( path.fork ) == delSet.end() ) path.fork = fork_;
    
    // Update status of spans and remove invalid ones
    if ( assembled && path.fork ) reviewSpans( path );
    
    // Assign divergent branches to be extended
    for ( Node* node : sideSet_ )
    {
        if ( delSet.find( node ) == delSet.end() ) sideNodes.push_back( node );
    }
    
    for ( int i = 0; i < path.spans.size(); )
    {
        if ( delSet.find( path.spans[i].node ) != delSet.end() ) path.spans.erase( path.spans.begin() + i );
        else i++;
    }
    
    for ( Node* node : path.path ) if ( delSet.find( node ) != delSet.end() ) assembled = false;
    
    
    return assembled;
}

bool PathReview::resolveAlleles( ConPath &con, NodeSet &delSet, bool &diverged )
{
    if ( find( con.paths[0].begin(), con.paths[0].end(), truncate_ ) != con.paths[0].end()
            || find( con.paths[1].begin(), con.paths[1].end(), truncate_ ) != con.paths[1].end() )
    {
        return false;
    }
    
    if ( con.diverged || diverged )
    {
        diverged = true;
        return true;
    }
    
    for ( Node* node : con.paths[0] )
    {
        if ( find( con.paths[1].begin(), con.paths[1].end(), node ) != con.paths[1].end() )
        {
            return false;
        }
    }
    
    float cover = Node::getAlleleCoverage( con.forks, con.paths, drxn_ );
    
    // Attempt to fold a diminutive path
    if ( cover < params.cover * 1.2 && Node::foldAlleles( pv_.nds[pv_.drxn], con.forks, con.paths, delSet, drxn_ ) ) return false;
    
    // Detect misassembly
    NodeSet forkSet( con.paths[0].begin(), con.paths[0].end() );
    forkSet.insert( con.paths[1].begin(), con.paths[1].end() );
    bool misassembled = false;
    if ( cover < params.cover * 1.5 && con.forks[!drxn_]->isMisassembled( pv_, con.forks[drxn_], forkSet, drxn_ ) )
    {
        Reassemble re( con.forks[!drxn_], pv_, false );
        if ( re.reassemble( pv_, delSet, true ) ) return false;
        misassembled = true;
    }
    
    // Characterize alleles
    NodeSet tSets[2];
    tSets[drxn_] = con.forks[drxn_]->getDrxnNodesInSet( pathSet_, drxn_, true );
    tSets[!drxn_] = con.forks[!drxn_]->getDrxnNodesInSet( pathSet_, !drxn_, true );
    int hits[2][2], reli[2][2], marks[2][2];
    int32_t farDist[2] = { con.forks[!drxn_]->ends_[drxn_], con.forks[!drxn_]->ends_[drxn_] };
    for ( int i : { 0, 1 } )
    {
        hits[i][0] = 1;
        hits[i][1] = 1;
        reli[i][0] = 1;
        reli[i][1] = 1;
        marks[i][0] = 1;
        marks[i][1] = 1;
        for ( int j = 0; j < con.paths[i].size(); j++ )
        {
            if ( con.paths[i][j]->farPairNodes_[0] )
            {
                farDist[i] = drxn_ ? min( farDist[i], con.paths[i][j]->farPairCoords_[0] ) 
                                   : max( farDist[i], con.paths[i][j]->farPairCoords_[0] );
            }
            con.paths[i][j]->getMarksCount( marks[i] );
            for ( auto &np : con.paths[i][j]->pairs_ )
            {
                for ( int k : { 0, 1 } )
                {
                    if ( tSets[k].find( np.first ) != tSets[k].end() )
                    {
                        hits[i][k] += np.second;
                        if ( np.first->isReliable( false ) ) reli[i][k] += np.second;
                    }
                }
            }
        }
        farDist[i] = min( abs( farDist[i] - con.forks[!drxn_]->ends_[drxn_] ), params.maxPeMean );
    }
    
    int pref;
    float fwdRatios[2], bckRatios[2];
    fwdRatios[0] = min( (float)1, float(hits[0][drxn_]) / float( max( 1, marks[0][drxn_] ) ) );
    fwdRatios[1] = min( (float)1, float(hits[1][drxn_]) / float( max( 1, marks[1][drxn_] ) ) );
    if ( reli[0][!drxn_] == reli[1][!drxn_] )
    {
        pref = hits[1][!drxn_] > hits[0][!drxn_];
        bckRatios[0] = min( (float)1, float(reli[0][!drxn_]) / float( max( 1, marks[0][!drxn_] ) ) );
        bckRatios[1] = min( (float)1, float(reli[1][!drxn_]) / float( max( 1, marks[1][!drxn_] ) ) );
    }
    else
    {
        pref = reli[1][!drxn_] > reli[0][!drxn_];
        bckRatios[0] = min( (float)1, float(hits[0][!drxn_]) / float( max( 1, marks[0][!drxn_] ) ) );
        bckRatios[1] = min( (float)1, float(hits[1][!drxn_]) / float( max( 1, marks[1][!drxn_] ) ) );
    }
    
    if ( cover > params.cover * 1.5 )
    {
        int scores[2] = { ( reli[0][!drxn_] - 1 ) * farDist[0], ( reli[1][!drxn_] - 1 ) * farDist[1] };
        if ( !scores[0] && ! scores[1] )
        {
            scores[0] = ( hits[0][!drxn_] - 1 ) * farDist[0];
            scores[1] = ( hits[1][!drxn_] - 1 ) * farDist[0];
        }
        if ( !scores[0] && ! scores[1] )
        {
            scores[0] = hits[0][drxn_];
            scores[1] = hits[1][drxn_];
        }
        if ( scores[0] > scores[1] * 1.5 )
        {
            con.paths[1][0]->dismantleNode( delSet, drxn_ );
        }
        else if ( scores[1] > scores[0] * 1.5 )
        {
            con.paths[0][0]->dismantleNode( delSet, drxn_ );
        }
        else
        {
            bool allUnreliable = true;
            for ( int i : { 0, 1 } )
            {
                for ( Node* node : con.paths[i] )
                {
                    allUnreliable = !node->setUnreliable() && allUnreliable;
                }
            }
            return allUnreliable;
        }
        return false;
    }
    
    float coverMod = max( float(1), cover / float( params.cover * .8 ) );
    float badness = ( ( bckRatios[pref] / ( bckRatios[0] + bckRatios[1] ) ) - 0.5 ) * coverMod;
    if ( badness > 0.25 )
    {
        badness += ( ( ( fwdRatios[!pref] / ( fwdRatios[0] + fwdRatios[1] ) ) - 0.5 ) * coverMod );
        if ( badness > 0.5 )
        {
            pv_.misassEst[drxn_] = con.paths[!pref][0]->ends_[!drxn_];
            pv_.misassEst[!drxn_] = pv_.misassEst[drxn_] + ( drxn_ ? params.readLen - params.maxPeMean 
                                                                   : params.maxPeMean - params.readLen );
            Reassemble re( con.paths[!pref][0], pv_, false );
            if ( re.reassemble( pv_, delSet ) ) return false;
        }
    }
    
    return true;
}

AltPath PathReview::resolveBranch( Node* node, NodeSet &tSet, NodeSet &delSet )
{
    AltPath ds( node );
    ds.far = node->ends_[drxn_];
    NodeFloatMap selfScores;
    Node* curr = node;
    while ( curr )
    {
        Node* best = NULL;
        NodeList delBranches, sideBranches;
        int bestScore = 0;
        float bestSelfScore = 0;
        
        for ( Node* nxt : curr->getNextNodes( drxn_ ) )
        {
            if ( usedSet_.find( nxt ) != usedSet_.end() ) continue;
            if ( pathSet_.find( nxt ) != pathSet_.end() ) continue;
            float score = 0;
            float selfScore = 0;
            float fwdReads = 0;
            bool doesContinue = false;
            for ( Node* fwd : nxt->getDrxnNodes( drxn_, true, true ) )
            {
                doesContinue = doesContinue || fwd->isContinue( drxn_ );
                fwdReads += fwd->getAdjustedReadCount( reliLimits_ );
                auto it = scores_.find( fwd );
                if ( it != scores_.end() ) score += it->second;
                auto it2 = selfScores.find( fwd );
                if ( it2 != selfScores.end() ) selfScore += it2->second;
            }
            if ( fwdReads >= 6 || !doesContinue ) delBranches.push_back( nxt );
            else
            {
                sideBranches.push_back( nxt );
                ds.doesContinue = ds.doesContinue || doesContinue;
            }
            if ( score > bestScore || ( score == bestScore && selfScore > bestSelfScore ) )
            {
                best = nxt;
                bestScore = score;
                bestSelfScore = selfScore;
            }
        }
        
        for ( Node* del : delBranches )
        {
            bool doDel = drxn_ ? del->ends_[1] < path_.back()->ends_[1] - 100
                               : path_.back()->ends_[0] + 100 < del->ends_[0];
            if ( del != best && doDel ) del->dismantleNode( delSet, drxn_ );
        }
        for ( Node* side : sideBranches )
        {
            if ( side != best ) sideSet_.insert( side );
        }
        
        if ( !best ) break;
        curr = best;
        ds.path.push_back( curr );
        ds.doesBridge = ds.doesBridge || bridgeSet_.find( curr ) != bridgeSet_.end();
        auto it = scores_.find( curr );
        auto it2 = selfScores.find( curr );
        if ( it != scores_.end() ) ds.hits += it->second;
        if ( it2 != selfScores.end() ) ds.hits += it2->second;
        ds.reads += curr->getAdjustedReadCount( reliLimits_ );
        if ( curr->farPairNodes_[0] )
            ds.far = drxn_ ? min( ds.far, curr->farPairCoords_[0] ) 
                           : max( ds.far, curr->farPairCoords_[0] );
        
        
        if ( tSet.find( curr ) != tSet.end() )
        {
            float coeff = curr->getCoverageCoeff();
            for ( auto &np : curr->pairs_ )
            {
                float adjScore = (float)np.second * coeff;
                auto r = selfScores.insert( make_pair( np.first, adjScore ) );
                if ( !r.second ) r.first->second += adjScore;
            }
        }
    }
    
    if ( !ds.path.empty() )
    {
        ds.far = max( 1, abs( ds.far - ds.path[0]->ends_[!drxn_] ) );
        int32_t dist = ( abs( ds.path[0]->ends_[!drxn_] - ds.path.back()->ends_[drxn_] ) * 0.9 ) - params.readLen;
        ds.limits[!drxn_] = ds.path[0]->ends_[!drxn_];
        ds.limits[drxn_] = drxn_ ? ds.limits[0] + dist : ds.limits[1] - dist;
        for ( Node* fwd : ds.path.back()->getDrxnNodes( drxn_, true, true ) )
            ds.doesContinue = ds.doesContinue || fwd->isContinue( drxn_ );
    }
    
    return ds;
}

void PathReview::resolveConverge( Node* forks[2], NodeSet &convSet, NodeSet &delSet )
{
    ConPath con;
    con.forks[!drxn_] = forks[0];
    con.forks[drxn_] = forks[1];
    NodeSet usedSet;
    NodeSet divSet;
    NodeSet tSet = forks[0]->getDrxnNodesInSet( pathSet_, !drxn_, true );
    tSet.insert( forks[1] );
    forks[1]->getDrxnNodesInSet( tSet, pathSet_, drxn_ );
    for ( int i : { 0, 1 } )
    {
        NodeSet pathSet, notSet, onlySet;
        Node* curr = forks[0];
        con.score = 0;
        while ( curr )
        {
            Node* branch = NULL;
            int bestPathScore;
            int bestTargScore;
            NodeSet currFwd = curr->getDrxnNodesInSet( convSet, drxn_ );
            for ( auto it = onlySet.begin(); it != onlySet.end(); )
            {
                if ( currFwd.find( *it ) == currFwd.end() 
                        || usedSet.find( *it ) != usedSet.end() ) it = onlySet.erase( it );
                else it++;
            }
            for ( Node* nxt : curr->getNextNodes( drxn_ ) )
            {
                if ( nxt == forks[1] ) continue;
                if ( convSet.find( nxt ) == convSet.end() )
                {
                    divSet.insert( nxt );
                    nxt->getDrxnNodes( divSet, drxn_ );
                }
                else
                {
                    int pathScore = 0;
                    int targScore = 0;
                    NodeSet nxtFwdSet = nxt->getDrxnNodesInSet( convSet, drxn_, true );
                    for ( Node* fwd : nxtFwdSet )
                    {
                        for ( auto &np : fwd->pairs_ )
                        {
                            if ( pathSet.find( np.first ) != pathSet.end() ) pathScore += np.second;
                            if ( notSet.find( np.first ) != notSet.end() ) pathScore -= np.second;
                            if ( usedSet.find( fwd ) != usedSet.end() ) continue;
                            if ( tSet.find( np.first ) != tSet.end() ) targScore += np.second;
                        }
                    }
                    for ( Node* fwd : onlySet )
                    {
                        if ( nxtFwdSet.find( fwd ) == nxtFwdSet.end() )
                        {
                            targScore -= fwd->getPairHitsTotal();
                        }
                    }
                    if ( !branch || ( pathScore > bestPathScore ) || ( pathScore == bestPathScore && targScore > bestTargScore ) )
                    {
                        branch = nxt;
                        bestPathScore = pathScore;
                        bestTargScore = targScore;
                    }
                }
            }
            if ( branch )
            {
                con.paths[i].push_back( branch );
                con.score += branch->getPairHitsTotal();
                pathSet.insert( branch );
                NodeSet nxtSet = branch->getDrxnNodesInSet( convSet, drxn_, true );
                NodeSet notNxtSet;
                for ( Node* nxt : curr->getNextNodes( drxn_ ) )
                {
                    if ( nxt != branch ) nxt->getDrxnNodesInSet( notNxtSet, convSet, drxn_ );
                }
                for ( Node* fwd : currFwd )
                {
                    if ( nxtSet.find( fwd ) == nxtSet.end() ) notSet.insert( fwd );
                }
                for ( Node* fwd : nxtSet )
                {
                    if ( fwd != branch && notNxtSet.find( fwd ) == notNxtSet.end() ) onlySet.insert( fwd );
                }
            }
            curr = branch;
        }
        
        usedSet.insert( con.paths[i].begin(), con.paths[i].end() );
        tSet = forks[0]->getDrxnNodesInSet( pathSet_, !drxn_, true );
        tSet.insert( forks[1] );
        forks[1]->getDrxnNodes( tSet, drxn_ );
    }
    
    pathSet_.insert( con.paths[0].begin(), con.paths[0].end() );
    pathSet_.insert( con.paths[1].begin(), con.paths[1].end() );
    
    
    vector<AltPath> conDivs;
    for ( int i : { 0, 1 } )
    {
        for ( int j = 0; j < con.paths[i].size(); j++ )
        {
            bool isFork = false;
            for ( Node* nxt : con.paths[i][j]->getNextNodes( drxn_ ) )
            {
                if ( pathSet_.find( nxt ) == pathSet_.end() ) isFork = true;
            }
            if ( !isFork ) continue;
            AltPath ap = resolveDiverge( con.paths[i][j], delSet );
            if ( ap.path.empty() ) continue;
            ap.score = ap.hits;
            Node* divNode = NULL;
            for ( Node* node : ap.path )
            {
                if ( convSet.find( node ) == convSet.end() )
                {
                    divNode = node;
                    break;
                }
            }
            if ( !divNode ) continue;

            for ( int k = j+1; k < con.paths[i].size(); k++ )
                ap.score -= con.paths[i][k]->getPairHitsTotal();
            
            if ( !ap.doesContinue ) ap.path[0]->dismantleNode( delSet, drxn_ );
            else if ( ap.reads <= 5 )
            {
                sideSet_.insert( divNode );
                con.diverged = true;
            }
            else if ( ap.doesContinue && ap.reads > 5 && ap.score >= 0 )
            {
                if ( !conDivs.empty() && ap.score < conDivs[0].score ) continue;
                if ( !conDivs.empty() && ap.score > conDivs[0].score ) conDivs.clear();
                conDivs.push_back( ap );
                con.diverged = true;
            }
        }
    }
    
    for ( AltPath &conDiv : conDivs ) divs_.push_back( conDiv );
    
    if ( !con.diverged )
    {
        for ( Node* node : convSet )
        {
            if ( pathSet_.find( node ) == pathSet_.end() ) node->dismantleNode( delSet, drxn_ );
        }
    
        for ( int i : { 0, 1 } )
        {
            for ( Node* node : convSet )
            {
                if ( pathSet_.find( node ) == pathSet_.end() ) node->dismantleNode( delSet, drxn_ );
            }
            for ( int j = 0; j < con.paths[i].size()-1; j++ )
            {
                NodeSet goodSet = { con.paths[i][j+1] };
                for ( int k = 0; k < con.paths[!i].size()-1; k++ )
                {
                    if ( con.paths[!i][k] == con.paths[i][j] ) goodSet.insert( con.paths[!i][k+1] );
                }
                for ( Node* nxt : con.paths[i][j]->getNextNodes( drxn_ ) )
                {
                    if ( sideSet_.find( nxt ) != sideSet_.end() ) continue;
                    if ( goodSet.find( nxt ) != goodSet.end() ) continue;
                    con.paths[i][j]->removeEdge( nxt, drxn_ );
                    nxt->removeEdge( con.paths[i][j], !drxn_ );
                    if ( nxt->edges_[!drxn_].empty() ) nxt->dismantleNode( delSet, drxn_ );
                }
            }
        }
    }
    
    if ( !drxn_ )
    {
        reverse( con.paths[0].begin(), con.paths[0].end() );
        reverse( con.paths[1].begin(), con.paths[1].end() );
    }
    
    if ( !con.diverged )
    {
        while ( !con.paths[0].empty() && !con.paths[1].empty()
                && con.paths[0][0] == con.paths[1][0] )
        {
            con.forks[0] = con.paths[0][0];
            con.paths[0].erase( con.paths[0].begin() );
            con.paths[1].erase( con.paths[1].begin() );
        }
        
        while ( !con.paths[0].empty() && !con.paths[1].empty()
                && con.paths[0].back() == con.paths[1].back() )
        {
            con.forks[1] = con.paths[0].back();
            con.paths[0].pop_back();
            con.paths[1].pop_back();
        }
        
        if ( con.paths[0].empty() || con.paths[1].empty() )
        {
            return;
        }
        
        if ( !divs_.empty() && con.doesConflict( divs_.back(), drxn_ ) && con.score >= 5 
                && ( !divs_.back().doesBridge && divs_.back().score < con.score ) )
        {
            divs_.back().path[0]->dismantleNode( delSet, drxn_ );
            divs_.pop_back();
        }
    }
    
    cons_.push_back( con );
}

AltPath PathReview::resolveDiverge( Node* fork, NodeSet &delSet )
{
    NodeSet tSet;
    NodeSet fwdSet = fork->getDrxnNodesNotInSet( pathSet_, drxn_, false );
    NodeSet qSet = fork->getDrxnNodesNotInSet( fwdSet, drxn_, false );
    NodeList qNodes( qSet.begin(), qSet.end() );
    for ( Node* fwd : fwdSet )
    {
        float fwdScore = 0;
        for ( ReadMark const &mark : fwd->marks_[!drxn_] )
        {
            for ( Node* q : qNodes )
            {
                if ( q->reads_.find( mark.id ) != q->reads_.end() )
                {
                    fwdScore++;
                    auto r = scores_.insert( make_pair( q, 1 ) );
                    if ( !r.second ) r.first->second++;
                    bridgeSet_.insert( fwd );
                    bridgeSet_.insert( q );
                    break;
                }
            }
        }
        for ( auto &np : fwd->pairs_ )
        {
            if ( fwdSet.find( np.first ) == fwdSet.end() )
                fwdScore += np.second * np.first->getCoverageCoeff();
        }
        
        auto r = scores_.insert( make_pair( fwd, fwdScore ) );
        if ( !r.second ) r.first->second += fwdScore;
        if ( r.first->second > 0.5 )
        {
            tSet.insert( fwd );
            fwd->getDrxnNodesInSet( tSet, fwdSet, !drxn_ );
        }
    }
    
    return resolveBranch( fork, tSet, delSet );
}

bool PathReview::resolveEnd( NodeSet &delSet )
{
    bool doesContinue = false;
    for ( Node* fwd : path_.back()->getDrxnNodes( drxn_, true, true ) )
    {
        doesContinue = doesContinue || fwd->isContinue( drxn_ );
    }
    if ( !doesContinue )
    {
        for ( AltPath &ap : divs_ )
        {
            if ( drxn_ ? path_.back()->ends_[1] < ap.path.back()->ends_[1]
                       : ap.path.back()->ends_[0] < path_.back()->ends_[0] )
            {
                NodeSet divSet = ap.path[0]->getDrxnNodes( drxn_, true, true );
                NodeSet mainSet = ap.fork->getDrxnNodesNotInSet( divSet, drxn_ );
                bool doMerge = false;
                for ( Node* div : divSet )
                {
                    for ( Node* node : mainSet )
                    {
                        for ( ReadMark &mark : div->marks_[drxn_] )
                        {
                            doMerge = doMerge || node->reads_.find( mark.id ) != node->reads_.end();
                            if ( doMerge ) break;
                        }
                    }
                }
                auto it = find( path_.begin(), path_.end(), ap.fork );
                if ( doMerge && it != path_.end() )
                {
                    NodeList path( it+1, path_.end() );
                    PathMerge merge( ap.fork, path, divSet, drxn_ );
                    if ( merge.merge( pv_, delSet, drxn_ ) )
                    {
                        truncate_ = ap.fork;
                        return false;
                    }
                    else if ( abs( ap.path.back()->ends_[drxn_] - path_.back()->ends_[drxn_] ) > params.maxPeMean )
                    {
                        ap.path[0]->dismantleNode( delSet, drxn_ );
                        return false;
                    }
                }
            }
        }
    }
    
    return true;
}

bool PathReview::resolveForks( NodeSet &delSet )
{
    float expFar = ( params.avgPeMean * .9 ) - params.readLen;
    
    Node* forks[2] = { NULL, NULL };
    NodeSet convSet, bckSet;
    for ( int i = 0; i < path_.size(); i++ )
    {
        // Check if this completes a convergence
        if ( forks[0] )
        {
            for ( Node* prv : path_[i]->getNextNodes( !drxn_ ) )
            {
                if ( convSet_.find( prv ) == convSet_.end() ) continue;
                forks[1] = path_[i];
                path_[i]->getDrxnNodesInSet( bckSet, convSet_, !drxn_ );
                for ( Node* con : convSet )
                {
                    if ( bckSet.find( con ) != bckSet.end() ) continue;
                    if ( pathSet_.find( con ) != pathSet_.end() ) continue;
                    forks[1] = NULL;
                }
                break;
            }
            
            if ( forks[1] )
            {
                resolveConverge( forks, convSet, delSet );
                forks[0] = NULL;
                forks[1] = NULL;
                convSet.clear();
                bckSet.clear();
            }
            else convSet.insert( path_[i] );
        }
        
        // Check if this is a forward fork
        for ( Node* nxt : path_[i]->getNextNodes( drxn_ ) )
        {
            if ( convSet_.find( nxt ) == convSet_.end() ) continue;
            if ( !forks[0] ) forks[0] = path_[i];
            path_[i]->getDrxnNodesInSet( convSet, convSet_, drxn_ );
            break;
        }
        
        if ( forks[0] || i == path_.size()-1 || path_[i]->edges_[drxn_].size() <= 1 ) continue;
        
        // Compete divergent against others
        AltPath ds = resolveDiverge( path_[i], delSet );
        if ( ds.path.empty() ) continue;

        float mod = max( (float)0.3, min( ds.far / expFar, (float)1.1 ) );
        if ( ds.far > params.maxPeMean ) 
            mod += min( (float)1, float( ds.far  - params.maxPeMean ) / float( params.maxMpMean ) );
        if ( ds.doesBridge ) mod = max( 1.5, mod + 0.5 );
        ds.score = ( ds.hits * 2 * mod ) - ( ds.reads / 1 + mod );

        bool doErase = false;
        if ( !divs_.empty() && ds.doesConflict( divs_.back(), drxn_ ) )
        {
            if ( ds.score > divs_.back().score )
            {
                divs_.back().path[0]->dismantleNode( delSet, drxn_ );
                divs_.pop_back();
            } else doErase = true;
        }

        if ( doErase )
        {
            ds.path[0]->dismantleNode( delSet, drxn_ );
        }

        divs_.push_back( ds );
        usedSet_.insert( ds.path.begin(), ds.path.end() );
    }
    
    for ( int i = 0; i < divs_.size(); )
    {
        int32_t dist = drxn_ ? path_.back()->ends_[1] - divs_[i].path.back()->ends_[1]
                             : divs_[i].path.back()->ends_[0] - path_.back()->ends_[0];
        dist -= params.maxPeMean;
        if ( dist < 0 ) divs_[i].score += divs_[i].hits;
        
        bool doExtend = divs_[i].score > -30 && divs_[i].doesContinue;
        bool doErase = !doExtend && dist > 0 && !i && divs_[i].doesContinue;
        if ( !divs_[i].doesContinue && ( divs_[i].score < 3 || divs_[i].hits < 3 ) ) doErase = true;
        
        if ( doExtend )
        {
            doErase = true;
            for ( Node* node : divs_[i].path )
            {
                if ( convSet_.find( node ) == convSet.end() )
                {
                    sideSet_.insert( node );
                    doErase = false;
                    break;
                }
            }
            
        }
        if ( !drxn_ ) for ( Node* node : pv_.nds[1] ) assert( !node->edges_[0].empty() );
        if ( doErase || delSet.find( divs_[i].path[0] ) != delSet.end() )
        {
            divs_[i].path[0]->dismantleNode( delSet, drxn_ );
            divs_.erase( divs_.begin() + i );
        } 
        else if ( !divs_[i].doesContinue && dist > 0 && delSet.find( divs_[i].path.back() ) == delSet.end() )
        {
            PathMerge merge( divs_[i].fork, divs_[i].path, pathSet_, drxn_ );
            if ( merge.merge( pv_, delSet, drxn_ ) )
            {
                truncate_ = divs_[i].fork;
                return false;
            }
            else divs_[i].path[0]->dismantleNode( delSet, drxn_ );
            divs_.erase( divs_.begin() + i );
        } else i++;
    }
    
    return true;
}

bool PathReview::resolveMisassembly( NodeSet &delSet )
{
    bool diverged = false;
    bool valid = true;
    Node* conBgn = !cons_.empty() ? cons_[0].forks[!drxn_] : NULL;
    Node* conEnd = NULL;
    int conCurr = 0;
    int32_t limit = path_.back()->ends_[drxn_];
    if ( !divs_.empty() ) limit = divs_[0].fork->ends_[drxn_];
    limit = limit + ( drxn_ ? - params.maxPeMean / 2 : params.maxPeMean / 2 );
    for ( int i = 0; i < path_.size(); i++ )
    {
        if ( delSet.find( path_[i] ) != delSet.end() ) break;
        if ( drxn_ ? limit < path_[i]->ends_[1] : path_[i]->ends_[0] < limit ) break;
        if ( truncate_ == path_[i] ) return false;
        if ( path_[i] == pv_.unspanned ) break;
        valid = valid && path_[i]->isValidated();
        
        
        // Conclude convergence
        if ( conEnd == path_[i] ) conEnd = NULL;
        
        // Scrutinise this homozygous node
        if ( !conEnd && path_[i]->drxn_ != 2 && !diverged && valid )
        {
            if ( path_[i] == pv_.weak || path_[i]->isMisassembled( pv_, drxn_ ) )
            {
                if ( abs( path_[i]->ends_[drxn_] - path_.back()->ends_[drxn_] ) < params.maxPeMean ) break;
                valid = valid && path_[i] != pv_.weak;
                Reassemble re( path_[i], pv_, path_[i] == pv_.weak );
                if ( re.reassemble( pv_, delSet ) ) return false;
            }
        }
        
        if ( !conEnd && !diverged && valid ) fork_ = path_[i];
        
        if ( path_[i] == conBgn )
        {
            if ( !resolveAlleles( cons_[conCurr], delSet, diverged ) ) return false;
            conEnd = cons_[conCurr].forks[drxn_];
            if ( ++conCurr < cons_.size() ) conBgn = cons_[conCurr].forks[!drxn_];
        }
        
        diverged = diverged || ( !divs_.empty() && path_[i] ==  divs_[0].fork );
    }
    
    return true;
}

bool PathReview::resolveUnspanned( Path &path, NodeSet &delSet )
{
    if ( !pv_.unspanned || pv_.unspanned->isEnded( drxn_ ) ) return false;
    
    for ( Node* nxt : pv_.unspanned->getNextNodes( drxn_ ) ) nxt->dismantleNode( delSet, drxn_ );
    pv_.unspanned->stop( 5, drxn_ );
    
    if ( !path.spans.empty() )
    {
        NodeSet fwdSet = path.spans.back().node->getDrxnNodes( drxn_, true, true );
        NodeSet bckSet;
        for ( Node* fwd : fwdSet )
        {
            if ( !fwd->farPairNodes_[0] ) continue;
            if ( fwdSet.find( fwd->farPairNodes_[0] ) != fwdSet.end() ) continue;
            bckSet.insert( fwd );
            fwd->getDrxnNodesInSet( bckSet, fwdSet, !drxn_ );
        }
        for ( Node* fwd : fwdSet )
        {
            if ( bckSet.find( fwd ) != bckSet.end() ) continue;
            fwd->dismantleNode( delSet, drxn_ );
        }
    }
    
    return true;
}

void PathReview::reviewSpans( Path &path )
{
    NodeSet forkFwd = path.fork->getDrxnNodes( drxn_ );
    
    for ( int i = 0; i < path.spans.size(); )
    {
        if ( forkFwd.find( path.spans[i].node ) != forkFwd.end() )
        {
            path.spans.erase( path.spans.begin() + i );
            continue;
        }
        else if ( !path.spans[i].complete || !path.spans[i].spanned )
        {
            NodeSet reliFwd;
            for ( Node* fwd : path.spans[i].node->getDrxnNodesNotInSet( forkFwd, drxn_ ) )
            {
                if ( fwd->isReliable() ) reliFwd.insert( fwd );
            }
            if ( !reliFwd.empty() ) path.spans[i].complete = true;
            
            NodeSet farFwd;
            for ( Node* fwd : reliFwd )
            {
                if ( !fwd->farPairNodes_[1] || !path.spans[i].complete ) continue;
                farFwd.insert( fwd->farPairNodes_[1] );
                fwd->farPairNodes_[1]->getDrxnNodesNotInSet( farFwd, forkFwd, drxn_ );
            }
            
            path.spans[i].spanned = farFwd.find( path.spans[i].node ) != farFwd.end();
        }
        i++;
    }
}

