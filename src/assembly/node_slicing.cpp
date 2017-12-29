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
#include "shared_functions.h"
#include <algorithm>

//void Node::slice( PathVars &pv, NodeSet &delSet, bool drxn )
//{
//    NodeSet bckSet = { this };
//    getDrxnNodes( bckSet, !drxn, params.getFurthestPeMean( pv.misassMark[!drxn], !drxn ) );
//    vector<int32_t> marks;
//    int32_t limits[2];
//    limits[0] = drxn ? pv.misassMark[0] - params.maxPeMean : pv.misassMark[1];
//    limits[1] = drxn ? pv.misassMark[0] : pv.misassMark[1] + params.maxPeMean;
//    
//    for ( Node* node : bckSet )
//    {
//        node->resetUnmarked( !drxn );
//        for ( ReadMark const &mark : node->marks_[!drxn] )
//            if ( limits[0] <= mark.mark && mark.mark <= limits[1] )
//                marks.push_back( mark.mark );
//    }
//    sort( marks.begin(), marks.end() );
//    
//    int32_t peWindow = max( params.readLen, params.avgPeMean / 2 );
//    int32_t peWindowCutoff = max( float( 3 ), float( peWindow * params.peCover ) / float( params.readLen * 6 ) );
//    bool misassembled = false;
//    for ( int i = 0; i + peWindowCutoff < marks.size(); i++ )
//        if ( abs( marks[i] - marks[i + peWindowCutoff] ) < peWindow )
//            misassembled = true;
//    
//    if ( !misassembled )
//    {
//        vector<PathSeq> pss = getCombinedSeqs( params.maxPeMean + params.readLen, !drxn );
//        PathSeq::setWeakspot( pss, pv.misassEst );
//        PathSeq::map( pv, pss, this, drxn );
//        unordered_set<ReadId> usedIds;
//        for ( PathSeq &ps : pss )  ps.tryComplete( pv, usedIds );
//        if ( !usedIds.empty() ) return;
//    }
//    
//    if ( slice( pv, misassembled, drxn ) ) return;
//    Node* thisNode = this;
//    if ( abs( pv.misassMark[!drxn] - ends_[!drxn] ) > params.readLen * 2 )
//    {
//        int32_t splitCoord = drxn ? max( ends_[0], pv.misassMark[0] - int( 1.5 * params.readLen) )
//                                  : min( ends_[1], pv.misassMark[1] + int( 1.5 * params.readLen ) );
//        if ( splitCoord != ends_[!drxn] )
//        {
//            splitCoord = findNextRead( splitCoord, drxn );
//            if ( splitCoord != ends_[drxn] )
//            {
//                thisNode = splitNode( splitCoord, pv.nodes, drxn, drxn );
//            }
//        }
//    }
//    thisNode->dismantleNode( delSet, drxn );
//}

bool Node::slice( PathVars &pv, bool misassembled, bool drxn )
{
    NodeSet tSet = getDrxnNodes( drxn, true, true );
    int32_t limits[2] = { pv.misassMark[0], pv.misassMark[1] };
    limits[drxn] += ( drxn ? params.maxPeMean : -params.maxPeMean );
    getDrxnNodes( tSet, drxn, limits[drxn] );
    vector<string> seqs;
    vector<ReadId> ids;
    vector<int32_t> markEsts;
    for ( Node* t : tSet )
    {
        if ( t->coverage_ < params.cover * 1.5 )
        {
            for ( ReadMark const &mark : t->marks_[drxn] )
            {
                if ( pv.usedIds.find( mark.id ) != pv.usedIds.end() ) continue;
                if ( limits[0] <= mark.mark && mark.mark <= limits[1] )
                {
                    seqs.push_back( pv.bwt.getSequence( mark.id ) );
                    ids.push_back( mark.id );
                    markEsts.push_back( mark.estimate );
                }
            }
        }
    }
    
    unordered_set<int> edges[2];
    for ( int i = 0; i < seqs.size(); i++ )
    {
        for ( int j = 0; j < seqs.size(); j++ )
        {
            if ( i != j && mapSeqOverlap( seqs[i], seqs[j], params.readLen * 0.4 ) )
            {
                edges[0].insert( i );
                edges[1].insert( j );
            }
        }
    }
    
    NodeList islands;
    int32_t dummy[2]{0}, dummyLimits[2] = { params.locusLimits[0], params.locusLimits[1] };
    ExtVars ev( pv.nds[pv.drxn], pv.nds[pv.drxn+3], dummy, pv.bwt );
    IslandVars iv( ev, drxn );
    bool didBridge = false;
    
    // Attempt to merge forward into the target nodes, if that succeeds, attempt to merge backwards - this is a bridge
    for ( int i = 0; i < seqs.size(); i++ )
    {
        bool doAdd = true;
        for ( Node* node : islands )
        {
            doAdd = doAdd && node->reads_.find( ids[i] ) == node->reads_.end();
        }
        for ( Node* node : pv.nds[pv.drxn] )
        {
            doAdd = doAdd && node->reads_.find( ids[i] ) == node->reads_.end();
        }
        if ( !doAdd ) continue;
        
        Node* bgn = new Node( seqs[i], ids[i], markEsts[i], drxn+3 );
        islands.push_back( bgn );
        int32_t extLimits[2] = { bgn->ends_[0] - params.maxPeMean, bgn->ends_[1] + params.maxPeMean };
        if ( bgn->sliceExtend( iv, extLimits[pv.drxn], pv.drxn ) && bgn->sliceExtend( iv, extLimits[pv.drxn], !pv.drxn ) )
        {
            didBridge = true;
            NodeSet thisSet = iv.merged[!pv.drxn];
            for ( Node* merge : iv.merged[!pv.drxn] ) merge->getDrxnNodes( thisSet, pv.drxn, true );
            for ( Node* merge : thisSet )
            {
                merge->drxn_ = pv.drxn;
                pv.nds[pv.drxn].push_back( merge );
                pv.newSet.insert( merge );
                islands.erase( remove( islands.begin(), islands.end(), merge ), islands.end() );
                for ( ReadId id : ids )
                {
                    if ( merge->reads_.find( id ) != merge->reads_.end() ) pv.usedIds.insert( id );
                }
                for ( Node* bck : merge->getNextNodes( !pv.drxn ) )
                {
                    if ( thisSet.find( bck ) == thisSet.end() && bck->drxn_ > 2 )
                    {
                        merge->removeEdge( bck, !pv.drxn );
                        bck->removeEdge( merge, pv.drxn );
                    }
                }
            }
            for ( Node* merge : iv.merged[!pv.drxn] )
            {
                merge->offsetForward( pv.drxn );
                merge->propagateValidation( dummyLimits, pv.drxn );
            }
            iv.merged[!pv.drxn].clear();
        }
    }
    
    bool doSplit = false, didSplit = false;
    int maxHitCount = 0;
    NodeSet usedSet;
    for ( Node* isl : islands )
    {
        if ( isl->edges_[drxn].empty() && usedSet.find( isl ) == usedSet.end() )
        {
            // Tally back pairs
            int thisMaxHitCount = 0, bckCount = 0;
            for ( Node* conn : isl->getConnectedNodes( true ) )
                for ( ReadMark const &mark : conn->marks_[drxn] )
                    for ( Node* node : pv.nds[pv.drxn] )
                        bckCount += node->reads_.find( mark.id ) != node->reads_.end();
                
            // Tally forward pairs
            NodeSet bckSet = isl->getDrxnNodes( !drxn, true, true );
            for ( Node* bck : bckSet )
            {
                if ( !bck->edges_[!drxn].empty() ) continue;
                int thisReadCount = 0, thisHitCount = 0;
                for ( Node* fwd : bck->getDrxnNodesInSet( bckSet, drxn, true ) )
                {
                    usedSet.insert( fwd );
                    thisReadCount += fwd->reads_.size();
                    for ( ReadMark const &mark : fwd->marks_[!drxn] )
                    {
                        for ( Node* t : tSet )
                        {
                            auto it = t->reads_.find( mark.id );
                            if ( it == t->reads_.end() ) continue;
                            ReadId pairId = params.getPairId( mark.id );
                            bool doCount = true;
                            for ( Node* node : pv.nds[pv.drxn] )
                            {
                                doCount = doCount && node->reads_.find( pairId ) == node->reads_.end();
                            }
                            thisHitCount += doCount;
                        }
                    }
                }
                thisMaxHitCount = max( thisMaxHitCount, thisHitCount - ( thisReadCount <= 3 ) );
            }
            
            // If multiple forward pairs and no back, try split
            if ( thisMaxHitCount > 1 && !bckCount )
            {
                maxHitCount = max( maxHitCount, thisMaxHitCount );
                bool thisSplit = false;
                for ( Node* conn : isl->getConnectedNodes( true ) )
                {
                    usedSet.insert( conn );
                    thisSplit = ( conn->edges_[drxn].empty() && conn->splitExtend( pv.bwt, this, pv.nds[pv.drxn], drxn ) ) || thisSplit;
                }
                if ( thisSplit ) didSplit = true;
                else doSplit = true;
            }
        }
    }
    
    for ( Node* node : islands )
    {
        node->dismantleNode();
        delete node;
    }
    
    if ( didSplit || didBridge ) return true;
    if ( !didSplit && !doSplit && !misassembled )
    {
        vector<ReadMark> newMarks;
        for ( ReadMark const &mark : marks_[drxn] )
        {
            if ( find( ids.begin(), ids.end(), mark.id ) == ids.end() )
            {
                newMarks.push_back( mark );
            }
        }
        marks_[drxn] = newMarks;
        return true;
    }
    
    return false;
}

bool Node::sliceOrBridge( PathVars &pv, Node* target, int32_t coords[2], NodeSet &delSet )
{
    drxn_ = pv.drxn + 3;
    NodeList islands = { this };
    int32_t dummy[2]{0};
    ExtVars ev( pv.nds[pv.drxn], islands, dummy, pv.bwt );
    IslandVars iv( ev, pv.drxn );
    NodeSet tSet = target->getDrxnNodes( pv.drxn, false, true );
    NodeSet hitSet, notAgain;
    
    bool didBridge = false;
    Node* curr = this;
    while ( curr )
    {
        int32_t limit = curr->ends_[!pv.drxn] + ( pv.drxn ? -params.maxPeMax : params.maxPeMax );
        if ( !curr->isContinue( !pv.drxn ) ) notAgain.insert( curr );
        curr->sliceExtend( iv, limit, !pv.drxn );
        
        // A merge has occurred
        if ( !iv.merged[!pv.drxn].empty() )
        {
            NodeSet mergeSet;
            for ( Node* merge : iv.merged[!pv.drxn] )
            {
                if ( merge->getBiggestOffset( !pv.drxn ) < params.readLen * 2 ) mergeSet.insert( merge );
                else
                {
                    assert( false );
                    merge->clearEdges( !pv.drxn );
                    merge->stop( 1, pv.drxn );
                }
            }
            iv.merged[!pv.drxn] = mergeSet;
            if ( !mergeSet.empty() )
            {
                didBridge = true;
                break;
            }
        }
        
        for ( Node* isl : islands )
        {
            if ( hitSet.find( isl ) != hitSet.end() ) continue;
            for ( ReadMark &mark : isl->marks_[pv.drxn] )
            {
                for ( Node* node : pv.nds[pv.drxn] )
                {
                    if ( tSet.find( node ) != tSet.end() ) continue;
                    if ( node->reads_.find( mark.id ) != node->reads_.end() )
                    {
                        hitSet.insert( isl );
                    }
                }
            }
        }
        Node* nxtCurr = NULL;
        NodeSet hitFwdSet;
        for ( Node* hit : hitSet ) hit->getDrxnNodes( hitFwdSet, pv.drxn, true );
        for ( Node* hit : hitSet )
        {
            if ( hitFwdSet.find( hit ) != hitFwdSet.end() 
                    || notAgain.find( hit ) != notAgain.end() ) continue;
            bool doesContinue = false;
            for ( Node* bck : hit->getDrxnNodes( !pv.drxn, true, true ) )
            {
                doesContinue = doesContinue || bck->isContinue( !pv.drxn );
            }
            if ( !doesContinue ) continue;
            if ( !nxtCurr || ( pv.drxn ? hit->ends_[0] < nxtCurr->ends_[0] 
                                       : nxtCurr->ends_[1] < hit->ends_[1] ) )
            {
                nxtCurr = hit;
            }
        }
        curr = nxtCurr;
    }
    
    // No need to slice and no bridge
    if ( !hitSet.empty() && !didBridge )
    {
        for ( Node* isl : islands )
        {
            isl->dismantleNode();
            delete isl;
        }
        
        return false;
    }
    
    if ( didBridge )
    {
        NodeSet mergeSet = iv.merged[!pv.drxn];
        for ( Node* merge : iv.merged[!pv.drxn] ) merge->getDrxnNodes( mergeSet, pv.drxn, true );

        for ( Node* merge : mergeSet )
        {
            islands.erase( remove( islands.begin(), islands.end(), merge ), islands.end() );
            merge->drxn_ = pv.drxn;
            pv.nds[pv.drxn].push_back( merge );
        }

        for ( Node* isl : islands )
        {
            isl->dismantleNode();
            delete isl;
        }
        
        for ( Node* merge : mergeSet ) merge->setValid();
        
        return true;
    }
    
    // Do slice
    if ( hitSet.empty() )
    {
        bool doSetUnreliable = false;
        int32_t splitCoord = coords[!pv.drxn];
        if ( target->getNextReadCoord( splitCoord, !pv.drxn, pv.drxn ) )
        {
            target = target->splitNode( splitCoord, pv.nds[pv.drxn], pv.drxn, pv.drxn );
            doSetUnreliable = target->ends_[1] - target->ends_[0] < params.readLen * 2;
        }
        else
        {
            NodeList hitNodes;
            vector<int32_t> hitCoords[2];
            target->overlapExtend( pv.nds[pv.drxn], coords, hitNodes, hitCoords, pv.drxn, pv.drxn );
            for ( int i = 0; i < hitNodes.size(); i++ )
            {
                addEdge( hitNodes[i], hitCoords[1][i] - hitCoords[0][i], pv.drxn );
            }
        }
        
        tSet.clear();
        target->getDrxnNodes( tSet, pv.drxn );
        
        NodeSet badSet = getNextNodes( pv.drxn );
        for ( Node* t : tSet )
        {
            for ( Node* isl : islands )
            {
                for ( ReadMark &mark : isl->marks_[!pv.drxn] )
                {
                    if ( t->reads_.find( mark.id ) == t->reads_.end() ) continue;
                    badSet.insert( t );
                    break;
                }
            }
        }
        
        clearEdges( pv.drxn );
        for ( Node* bad : badSet )
        {
            if ( bad->drxn_ != pv.drxn ) continue;
            bad->dismantleNode( delSet, pv.drxn );
        }
        
        for ( Node* isl : islands )
        {
            isl->dismantleNode();
            delete isl;
        }
        
        if ( doSetUnreliable ) target->setUnreliable();
    }
    
    return true;
}

bool Node::sliceExtend( IslandVars &iv, int32_t limit, bool drxn )
{
    int mergeCount = iv.merged[drxn].size();
    NodeSet extSet = { this };
    while ( !extSet.empty() )
    {
        for ( Node* node : extSet )
        {
            node->extendCount_ = 20;
            node->extendIsland( iv, drxn );
        }
        extSet.clear();

        NodeSet currSet = { this };
        NodeSet usedSet;
        NodeIntMap edgeCounts;
        edgeCounts[this] = 0;
        while ( !currSet.empty() )
        {
            NodeSet nxtSet;
            for ( Node* curr : currSet )
            {
                int currCount = edgeCounts[curr];
                if ( usedSet.find( curr ) == usedSet.end() && curr->coverage_ < params.cover * 1.3
                        && ( drxn ? curr->ends_[1] < limit : limit < curr->ends_[0] ) )
                {
                    for ( Node* nxt : curr->getNextNodes( drxn ) )
                    {
                        int nxtCount = currCount + ( curr->edges_[drxn].size() > nxt->edges_[!drxn].size() );
                        if ( nxtCount <= 2 && nxt->coverage_ < params.cover * 1.3 
                                && ( drxn ? nxt->ends_[1] < limit : limit < nxt->ends_[0] ) )
                        {
                            edgeCounts[nxt] = nxtCount;
                            if ( nxt->isContinue( drxn ) ) extSet.insert( nxt );
                            else nxtSet.insert( nxt );
                        }
                    }
                }
            }
            currSet = nxtSet;
        }
    }
    
    return iv.merged[drxn].size() > mergeCount;
}

bool Node::splitExtend( Querier &bwt, Node* target, NodeList &nodes, bool drxn )
{
    vector<Extension> exts = bwt.mapExtensions( seq_, drxn );
    bool didSplit = false;
    for ( Extension &ext : exts )
    {
        MergeHit merge;
        target->checkExtensionMerge( ext, merge );
        if ( merge.node && (*merge.coords)[!drxn] != merge.node->ends_[!drxn] )
        {
            didSplit = merge.node->splitNode( (*merge.coords)[!drxn], nodes, drxn, drxn )|| didSplit;
        }
    }
    return didSplit;
}

