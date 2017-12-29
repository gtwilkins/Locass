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
#include "node_structs.h"
#include <algorithm>

void Node::addExtension( Extension &ext, IslandVars &iv, bool doesBranch, bool drxn )
{
    MergeHit merge;
    checkExtension( ext, iv, merge );
    
    if ( merge.node && merge.node != this && iv.ev.ante.find( merge.node ) == iv.ev.ante.end() )
    {
        addExtensionMerge( merge, ext, iv, doesBranch, drxn );
    }
    
    if ( !merge.node )
    {
        if ( doesBranch )
        {
            Node* node = new Node( getSeqEnd( ext.maxOverLen, drxn ), ext, ends_[drxn], drxn, iv.drxn );
            addEdge( node, ext.maxOverLen, drxn );
            iv.ev.island.push_back( node );
            if ( !ext.fwdExts.empty() )
            {
                iv.ev.ante.insert( this );
                node->addExtensions( ext.fwdExts, iv, drxn );
                iv.ev.ante.erase( this );
            }
            node->setCoverage();
        }
        else 
        {
            appendNode( ext, drxn );
            if ( !ext.fwdExts.empty() )
            {
                addExtensions( ext.fwdExts, iv, drxn );
            }
        }
        setCoverage();
    }
}

void Node::addExtensionMerge( MergeHit &merge, Extension &ext, IslandVars &iv, bool doesBranch, bool drxn )
{
    Node* thisNode = this;
    
    if ( !ext.overlaps.empty() )
    {
        if ( doesBranch )
        {
            thisNode = new Node( getSeqEnd( ext.maxOverLen, drxn ), ext, ends_[drxn], drxn, iv.drxn );
            addEdge( thisNode, ext.maxOverLen, drxn );
            iv.ev.island.push_back( thisNode );
            thisNode->setCoverage();
        }
        else
        {
            appendNode( ext, drxn );
            setCoverage();
        }
    }
    
    if ( merge.node->drxn_ != !iv.drxn )
    {
        bool isIslandMerge = ( merge.node->drxn_ > 2 );
        Node* mergeNode = NULL;
        
        if ( drxn == iv.drxn && isIslandMerge )
        {
            mergeNode = merge.node->splitNode( (*merge.coords)[!drxn], iv.ev.island, iv.drxn, drxn );
            thisNode->addEdge( mergeNode, merge.overlap, drxn );
        }
        else if ( drxn != iv.drxn )
        {
            mergeNode = merge.node->mergeNode( ( isIslandMerge ? iv.ev.island : iv.ev.nodes ), merge.coords, iv.drxn, drxn );
            mergeNode->addEdge( thisNode, merge.overlap, !drxn );
        }
        
        if ( mergeNode && mergeNode != merge.node && iv.merged[drxn].find( merge.node ) != iv.merged[drxn].end() )
        {
            iv.merged[drxn].erase( merge.node );
            iv.merged[drxn].insert( mergeNode );
        }

        if ( !isIslandMerge )
        {
            iv.merged[drxn].insert( thisNode );
        }
    }
}

void Node::addExtensions( vector<Extension> &exts, IslandVars &iv, bool drxn )
{
    int32_t endCoord = ends_[drxn];
    
    bool doesBranch = exts.size() > 1 || !edges_[drxn].empty() || clones_;
    for ( Extension &ext : exts )
    {
        addExtension( ext, iv, doesBranch, drxn );
    }
    
    stop_[drxn] = ( edges_[drxn].empty() && ends_[drxn] == endCoord ) ? 1 : stop_[drxn];
}

bool Node::canValidate( IslandVars &iv )
{
    if ( validLimits_[1] <= validLimits_[2] )
    {
        return true;
    }
    
    int hits = 0;
    NodeSet tSet = getDrxnNodes( !iv.drxn );
    for ( Node* q : getDrxnNodes( iv.drxn ) )
    {
        for ( auto &pn : q->pairs_ )
        {
            if ( tSet.find( pn.first ) != tSet.end() )
            {
                hits += pn.second;
            }
        }
    }
    for ( int i( 0 ); i < min( hits, 2 ); i++ )
    {
        pushValidLimits( ends_[0], 0 );
        pushValidLimits( ends_[1], 1 );
    }
    
    return validLimits_[1] <= validLimits_[2];
}

void Node::checkExtension( Extension &ext, IslandVars &iv, MergeHit &merge )
{
    for ( Node* &node : iv.origin )
    {
        if ( node->checkExtensionMerge( ext, merge ) ) return;
    }
    for ( Node* &node : iv.ev.nodes )
    {
        if ( node->checkExtensionMerge( ext, merge ) ) return;
    }
    for ( Node* &node : iv.ev.island )
    {
        if ( node->checkExtensionMerge( ext, merge ) ) return;
    }
    for ( Node* const &node : iv.ante )
    {
        if ( node->checkExtensionMerge( ext, merge ) ) return;
    }
}

void Node::extendIsland( IslandVars &iv, bool drxn )
{
    iv.ev.ante = getDrxnNodes( !drxn );
    while ( extendCount_ && isContinue( drxn ) )
    {
        vector<Extension> exts = iv.ev.bwt.mapExtensions( seq_, drxn );
        addExtensions( exts, iv, drxn );
        extendCount_--;
    }
    setCoverage();
//    setCoverage( iv.ev, iv.drxn, drxn, true );
}

void Node::islandDelete( IslandVars &iv, Node* node )
{
    node->dismantleNode();
    assert ( find( iv.ev.island.begin(), iv.ev.island.end(), node ) != iv.ev.island.end() );
    iv.ev.island.erase( remove( iv.ev.island.begin(), iv.ev.island.end(), node ), iv.ev.island.end() );
    iv.merged[0].erase( node );
    iv.merged[1].erase( node );
    delete node;
}

bool Node::islandPairToBase( IslandVars &iv, NodeList &tNodes )
{
    for ( Node* t : tNodes )
    {
        if ( !paired_ || paired_->find( t ) == paired_->end() )
        {
            for ( auto it = marks_[iv.drxn].begin(); it != marks_[iv.drxn].end(); )
            {
                auto hit = t->reads_.find( (*it).id );
                if ( hit != t->reads_.end() )
                {
                    auto selfResult = pairs_.insert( make_pair( t, 1 ) );
                    if ( !selfResult.second )
                    {
                        selfResult.first->second++;
                    }
                    auto pairResult = t->pairs_.insert( make_pair( this, 1 ) );
                    if ( !pairResult.second )
                    {
                        pairResult.first->second++;
                    }
                    
                    SeqNum pairId = params.getPairId( (*it).id );
                    islandUpdateValid( pairId );
                    if ( params.isReadPe( pairId ) )
                    {
                        iv.peReads.insert( pairId );
                        if ( t->reliable_ )
                        {
                            iv.peReadsReliable.insert( pairId );
                        }
                    }
                    else
                    {
                        iv.mpReads.insert( pairId );
                    }
                    
                    it = marks_[iv.drxn].erase( it );
                    continue;
                }
                it++;
            }
            
            if ( paired_ )
            {
                paired_->insert( t );
            }
        }
    }
    
    for ( auto &pn : pairs_ )
    {
        if ( pn.first->drxn_ <= 2 )
        {
            return true;
        }
    }
    
    return false;
}

void Node::islandUpdateValid( SeqNum &readId )
{
    auto it = reads_.find( readId );
    assert( it != reads_.end() );
    if ( it != reads_.end() && ( it->second[0] < validLimits_[1] || validLimits_[2] < it->second[1] ) )
    {
        // Update own valid limits
        if ( it->second[0] < validLimits_[0] )
        {
            validLimits_[1] = min( validLimits_[0], it->second[1] );
            validLimits_[0] = it->second[0];
        }
        else if ( it->second[0] != validLimits_[0] )
        {
            validLimits_[1] = min( validLimits_[1], it->second[0] );
        }
        
        if ( validLimits_[3] < it->second[1] )
        {
            validLimits_[2] = max( validLimits_[3], it->second[0] );
            validLimits_[3] = it->second[1];
        }
        else if ( validLimits_[3] != it->second[1] )
        {
            validLimits_[2] = max( validLimits_[2], it->second[1] );
        }
        
        // If there are any hits in either direction, update all intervening valid limits
        for ( bool drxn : { 0, 1 } )
        {
            NodeSet fwdSet = getDrxnNodes( drxn, true, false );
            NodeSet bckSet;
            bool anyFwd = false;
            
            for ( Node* fwd : fwdSet )
            {
                if ( drxn ? fwd->ends_[0] < fwd->validLimits_[3]
                          : fwd->validLimits_[0] < fwd->ends_[1] )
                {
                    fwd->getDrxnNodesInSet( bckSet, fwdSet, !drxn );
                    fwd->validLimits_[(!drxn) * 2] = fwd->validLimits_[(!drxn) * 2 + 1] = fwd->ends_[!drxn];
                    anyFwd = true;
                }
            }
            
            if ( anyFwd )
            {
                validLimits_[drxn * 2] = validLimits_[drxn * 2 + 1] = ends_[drxn];
            }
            
            for ( Node* bck : bckSet )
            {
                bck->validLimits_[0] = bck->validLimits_[1] = bck->ends_[0];
                bck->validLimits_[2] = bck->validLimits_[3] = bck->ends_[1];
            }
        }
    }
}

bool Node::islandReview( IslandVars &iv, NodeSet &islandSet )
{
    int readCount = 0;
    NodeIntMap pairMap;
    NodeSet fwdSet;
    
    for ( Node* node : islandSet )
    {
        int thisHits = 0;
        for ( unordered_set<SeqNum> *readIds : { &iv.peReads, &iv.mpReads } )
        {
            for ( const SeqNum &readId : *readIds )
            {
                if ( node->reads_.find( readId ) != node->reads_.end() )
                {
                    thisHits++;
                }
            }
        }
        
        if ( thisHits )
        {
            node->getDrxnNodes( fwdSet, iv.drxn );
            pairMap[node] = thisHits;
            readCount += thisHits;
        }
    }
    
    Node* best = NULL;
    int bestHits = 0;
    
    for ( auto &np : pairMap )
    {
        if ( fwdSet.find( np.first ) == fwdSet.end() )
        {
            int thisHits = np.second;
            for ( Node* fwd : np.first->getDrxnNodes( iv.drxn ) )
            {
                auto it = pairMap.find( fwd );
                thisHits += it != pairMap.end() ? it->second : 0;
            }
            if ( thisHits > bestHits )
            {
                best = np.first;
                bestHits = thisHits;
            }
        }
    }
    
    bool doDelete = true;
    
    if ( best && readCount > 1 )
    {
        NodeSet propagated;
        best->offsetIsland( propagated, iv.drxn );
        best->bridgeIslandOffset( iv, islandSet, iv.drxn );
        for ( Node* node : islandSet )
        {
            doDelete = doDelete && ( iv.drxn ? node->ends_[1] < iv.limit : iv.limit < node->ends_[0] );
        }
    }
    
    if ( doDelete )
    {
        iv.ev.del.insert( islandSet.begin(), islandSet.end() );
    }
}

void Node::islandReview( IslandVars &iv, NodeSet &islandSet, NodeSetList &peIslands, NodeSetList &mpIslands )
{
    int nodeCount = 0, readCount = 0;
    bool isPe = false;
    
    NodeSet mainlandSet;
    
    for ( Node* node : islandSet )
    {
        bool thisBridge = false;
        for ( const SeqNum &readId : iv.peReads )
        {
            if ( node->reads_.find( readId ) != node->reads_.end() )
            {
                isPe = true;
                readCount++;
            }
        }
        for ( const SeqNum &readId : iv.mpReads )
        {
            if ( node->reads_.find( readId ) != node->reads_.end() )
            {
                readCount++;
            }
        }
        
        for ( auto &np : node->pairs_ )
        {
            if ( np.first->drxn_ <= 2 )
            {
                mainlandSet.insert( np.first );
                thisBridge = true;
            }
        }
        
        nodeCount += thisBridge;
    }
    
    if ( isPe > 0 && readCount > 1 )
    {
        peIslands.push_back( islandSet );
    }
    else if ( readCount > 1 && nodeCount > 1 && mainlandSet.size() > 1 )
    {
        mpIslands.push_back( islandSet );
    }
}

bool Node::islandReviewUnbridged( IslandVars &iv, NodeSet &islandSet )
{
    int peCount = 0, peReliCount = 0, mpCount = 0, nodeCount = 0, readCount = 0;
    for ( Node* node : islandSet )
    {
        int thisCount = 0;
        for ( const SeqNum &readId : iv.peReads )
        {
            if ( node->reads_.find( readId ) != node->reads_.end() )
            {
                peCount++;
                peReliCount += find( iv.peReadsReliable.begin(), iv.peReadsReliable.end(), readId ) != iv.peReadsReliable.end();
                thisCount++;
            }
        }
        for ( const SeqNum &readId : iv.mpReads )
        {
            if ( node->reads_.find( readId ) != node->reads_.end() )
            {
                mpCount++;
                thisCount++;
            }
        }
        
        if ( thisCount > 0 )
        {
            nodeCount++;
            readCount += thisCount;
        }
    }
    
    if ( readCount >= 8 && ( ( peCount > 0 && mpCount > 0 ) || ( peReliCount >= 4 ) ) )
    {
        return true;
    }
    
    return false;
}

void Node::islandSetExtend( IslandVars &iv, NodeSet &islandSet, NodeSet* &extSets )
{
    NodeSet fwdSet;
    NodeList peNodes, mpNodes;
    NodeIntMap pairMap = Node::islandSetExtendGetBridges( iv, islandSet, peNodes, mpNodes, fwdSet );
    
    if ( pairMap.empty() )
    {
        iv.ev.del.insert( islandSet.begin(), islandSet.end() );
        return;
    }
    
    NodeSet tested;
    Node::islandSetValidate( iv, peNodes, tested );
    Node::islandSetValidate( iv, mpNodes, tested );
    
    int totalHits = 0;
    NodeList anchorList = Node::islandSetExtendGetAnchors( iv, pairMap, fwdSet, totalHits );
    
    if ( !anchorList.empty() )
    {
        
        
        int32_t cutoffs[2];
        cutoffs[0] = iv.drxn ? anchorList[0]->ends_[0] - ( params.maxPeMean / 2 )
                             : anchorList[0]->ends_[1] - params.maxPeMax;
        cutoffs[1] = iv.drxn ? anchorList[0]->ends_[0] + params.maxPeMax
                             : anchorList[0]->ends_[1] + ( params.maxPeMean / 2 );
        
        NodeSet usedBckSet, usedFwdSet;
        int bestHits = 0;
        for ( int i( 0 ); i < anchorList.size(); i++ )
        {
            if ( usedBckSet.find( anchorList[i] ) == usedBckSet.end()
                    && !anchorList[i]->islandSetExtendMulti( iv, extSets, usedBckSet, cutoffs, !i, bestHits, !iv.drxn ) 
                    && i == 0 && totalHits < 30 )
            {
                anchorList[i]->islandSetExtendMulti( iv, extSets, usedFwdSet, cutoffs, false, bestHits, iv.drxn );
            }
        }
        
        if ( totalHits < 30 )
        {
            for ( NodeList* nodes : { &peNodes, &mpNodes } )
            {
                for ( Node* node : *nodes )
                {
                    // Set extend for under-explored, paired nodes
                    for ( bool drxn : { 0, 1 } )
                    {
                        if ( node->isContinue( drxn ) && node->ends_[1] - node->ends_[0] < params.maxPeMean )
                        {
                            extSets[drxn].insert( node );
                        }
                        for ( Node* nxt : node->getNextNodes( drxn ) )
                        {
                            if ( nxt->isContinue( drxn ) && nxt->ends_[1] - nxt->ends_[0] < params.maxPeMean )
                            {
                                extSets[drxn].insert( nxt );
                            }
                        }
                    }
                }
            }
        }
    }
    else if ( iv.round < 5 )
    {
        for ( Node* node : peNodes )
        {
            if ( iv.ev.del.find( node ) == iv.ev.del.end() )
            {
                int32_t cutoffs[2] = { node->ends_[1] - params.maxPeMean, node->ends_[0] + params.maxPeMean };
                node->islandSetExtendSingle( iv, extSets, cutoffs, !iv.drxn );
            }
        }
        
        for ( Node* node : mpNodes )
        {
            if ( iv.ev.del.find( node ) == iv.ev.del.end() && 
                    ( node->isReliable( true ) || node->isValidated( !iv.drxn ) ) )
            {
                int32_t cutoffs[2] = { node->ends_[1] - params.maxPeMean, node->ends_[0] + params.maxPeMean };
                int bestHits = 0;
                NodeSet usedSet;
                if ( node->isValidated( !iv.drxn ) )
                {
                    node->islandSetExtendMulti( iv, extSets, usedSet, cutoffs, false, bestHits, !iv.drxn );
                }
                else
                {
                    node->islandSetExtendSingle( iv, extSets, cutoffs, 2 );
                }
            }
        }
    }
}

NodeList Node::islandSetExtendGetAnchors( IslandVars &iv, NodeIntMap &pairMap, NodeSet &fwdSet, int &totalHits )
{
    NodeList anchorList;
    NodeSet usedSet;
    
    bool first = true;
    int usedHits = 0;
    
    while ( anchorList.size() <= 2 )
    {
        Node* bestNode = NULL;
        int bestHits = 1;
        
        for ( auto &np : pairMap )
        {
            if ( first )
            {
                totalHits += np.second;
            }
            
            if ( fwdSet.find( np.first ) == fwdSet.end() && usedSet.find( np.first ) == usedSet.end() )
            {
                int hits = np.second;
                for ( Node* fwd : np.first->getDrxnNodes( iv.drxn ) )
                {
                    if ( usedSet.find( fwd ) == usedSet.end() )
                    {
                        auto it = pairMap.find( fwd );
                        if ( it != pairMap.end() )
                        {
                            hits += it->second;
                        }
                    }
                }
                
                if ( hits > bestHits )
                {
                    bestNode = np.first;
                    bestHits = hits;
                }
            }
        }
        
        if ( !bestNode ) break;
        anchorList.push_back( bestNode );
        usedSet.insert( bestNode );
        bestNode->getDrxnNodes( usedSet, iv.drxn );
        usedHits += bestHits;
        first = false;
    }
    
    if ( totalHits - usedHits > (int)params.cover && !anchorList.empty() )
    {
        anchorList[0]->getConnectedNodes( iv.ev.del, true );
        anchorList.clear();
    }
    
    // Backtrack two nodes if best anchor is in a fragmented group of nodes
    if ( !anchorList.empty() && pairMap[anchorList[0]] < 2 )
    {
        for ( Node* nxt : anchorList[0]->getNextNodes( iv.drxn ) )
        {
            if ( pairMap.find( nxt ) == pairMap.end() && nxt->reads_.size() <= 3 )
            {
                for ( Node* nxtNxt : nxt->getNextNodes( iv.drxn ) )
                {
                    bool doAdd = false;
                    for ( Node* bck : nxtNxt->getNextNodes( !iv.drxn ) )
                    {
                        doAdd = doAdd || usedSet.find( bck ) == usedSet.end();
                    }
                    if ( doAdd )
                    {
                        anchorList.push_back( nxtNxt );
                    }
                }
            }
        }
    }
    
    return anchorList;
}

NodeIntMap Node::islandSetExtendGetBridges( IslandVars &iv, NodeSet &islandSet, NodeList &peNodes, NodeList &mpNodes, NodeSet &fwdSet )
{
    NodeIntMap pairMap;
    
    for ( Node* node : islandSet )
    {
        bool isPaired = false;
        for ( auto &np : node->pairs_ )
        {
            if ( np.first->drxn_ <= 2 )
            {
                isPaired = true;
                auto r = pairMap.insert( make_pair( node, np.second ) );
                if ( !r.second )
                {
                    r.first->second += np.second;
                }
                node->getDrxnNodes( fwdSet, iv.drxn );
            }
        }
        
        if ( isPaired )
        {
            if ( node->anyReadInNode( iv.peReads ) )
            {
                peNodes.push_back( node );
            }
            else if ( node->anyReadInNode( iv.mpReads ) )
            {
                mpNodes.push_back( node );
            }
        }
    }
    
    return pairMap;
}

bool Node::islandSetExtendMulti( IslandVars &iv, NodeSet* &extSets, NodeSet &usedSet, int32_t* cutoffs, bool isBestAnchor, int &bestHits, int drxn )
{
    NodeSet fwdSet = getDrxnNodes( drxn, false, true );
    getNextNodes( fwdSet, !drxn );
    int32_t limits[2] = { ends_[0] - params.maxPeMax, ends_[1] + params.maxPeMax };
    ScoreMap scores = getScoreMap( fwdSet, limits, drxn );
    NodeList forks;
    
    for ( int i( 0 ); i < 1 + isBestAnchor; i++ )
    {
        NodeSet currSet = getNextNodes( !drxn );
        if ( currSet.empty() )
        {
            currSet.insert( this );
        }
        Node* bestBranch = NULL;
        Node* currBranch = NULL;
        Score currScore;
        Score cumulThisScore;
        Score cumulRealScore;
        
        while ( !currSet.empty() )
        {
            for ( Node* curr : currSet )
            {
                for ( Node* nxt : curr->getNextNodes( drxn ) )
                {
                    Score thisScore;
                    Score realScore;

                    for ( Node* fwd : nxt->getDrxnNodes( drxn, false, true ) )
                    {
                        if ( usedSet.find( fwd ) == usedSet.end() )
                        {
                            thisScore += scores[fwd];
                        }
                        realScore += scores[fwd];
                    }

                    if ( thisScore[1] > currScore[1] || ( thisScore[1] == currScore[1] && thisScore[0] < currScore[0] ) )
                    {
                        bestBranch = nxt;
                        currScore = thisScore;
                    }
                }
            }
            
            currSet.clear();
            
            if ( bestBranch )
            {
                currSet.insert( bestBranch );

                if ( usedSet.find( bestBranch ) == usedSet.end() )
                {
                    cumulThisScore += scores[bestBranch];
                }
                cumulRealScore += scores[bestBranch];
                
                currBranch = bestBranch;
                bestBranch = NULL;
                currScore.clear();
            }
        }
        
        if ( !currBranch ) break;
        if ( cumulThisScore[1] < 2 )
        {
            if ( isBestAnchor && cumulRealScore[1] < bestHits - 2 ) break;
            if ( !isBestAnchor && cumulThisScore[0] > 3 ) break;
        }
        
        bestHits = max( (int)cumulThisScore[1], bestHits ); 
        
        if ( currScore[1] < 2 && currScore[0] - ( currScore[1] * 3 ) > 12 ) break;
        
        if ( usedSet.find( currBranch ) == usedSet.end() )
        {
            forks.push_back( currBranch );
        }
        
        usedSet.insert( currBranch );
        currBranch->getDrxnNodes( usedSet, drxn );
        currBranch->getDrxnNodes( usedSet, !drxn );
    }
    
    bool didExtend = false;
    
    for ( Node* fork : forks )
    {
        didExtend = fork->islandSetExtendSingle( iv, extSets, cutoffs, drxn ) || didExtend;
    }
    
    return didExtend;
}

bool Node::islandSetExtendSingle( IslandVars &iv, NodeSet* &extSets, int32_t* cutoffs, int drxns )
{
    bool didExtend = false;
    
    for ( bool drxn : { 0, 1 } )
    {
        if ( drxn == drxns || drxn == 2 )
        {
            NodeIntMap edgeMap;
            NodeSet currSet = { this };
            edgeMap[this] = 0;
            while ( !currSet.empty() )
            {
                NodeSet nxtSet;
                for ( Node* curr : currSet )
                {
                    if ( edgeMap[curr] < 5 && ( drxn ? curr->ends_[1] < cutoffs[1] : cutoffs[0] < curr->ends_[0] ) )
                    {
                        if ( curr->isContinue( drxn ) )
                        {
                            extSets[drxn].insert( curr );
                            didExtend = true;
                        }
                        else
                        {
                            for ( Node* nxt : curr->getNextNodes( drxn ) )
                            {
                                if ( edgeMap.find( nxt ) == edgeMap.end() )
                                {
                                    edgeMap[nxt] = edgeMap[curr] + min( 4, (int)curr->edges_[drxn].size() - (int)nxt->edges_[!drxn].size() );
                                    nxtSet.insert( nxt );
                                }
                            }
                        }
                    }
                }
                currSet = nxtSet;
            }
        }
    }
    
    return didExtend;
}

void Node::islandSetValidate( IslandVars &iv, NodeList &testNodes, NodeSet &tested )
{
    NodeSet vSets[2];
    for ( Node* node : testNodes )
    {
        node->validate( iv.drxn );
        tested.insert( node );
    }
    
    for ( bool drxn : { 0, 1 } )
    {
        for ( Node* node : testNodes )
        {
            if ( node->isValidated( drxn ) )
            {
                vSets[drxn].insert( node );
            }
        }
        
        NodeSet nxtSet;
        for ( Node* v : vSets[drxn] )
        {
            for ( Node* nxt : v->getNextNodes( drxn ) )
            {
                if ( tested.find( nxt ) == tested.end() && nxt->validate( drxn ) )
                {
                    nxtSet.insert( nxt );
                }
                tested.insert( nxt );
            }
        }
        vSets[drxn] = nxtSet;
    }
}

void Node::mergeIsland( ExtVars &ev, bool drxn, bool trimBack )
{
    // Set coordinates compatible with mainland
    offsetForward( drxn, true );
    
    if ( ends_[0] < -40000 )
    {
        int x = 0;
    }
    
    if ( trimBack )
    {
        // Trim island nodes that will not become part of mainland if required
        NodeSet fwdSet = getDrxnNodes( drxn, true, true );
        for ( Node* fwd : fwdSet )
        {
            for ( Node* prv : fwd->getNextNodes( !drxn ) )
            {
                if ( fwdSet.find( prv ) == fwdSet.end() && prv->drxn_ > 2 )
                {
                    fwd->removeEdge( prv, !drxn );
                    prv->removeEdge( fwd, drxn );
                }
            }
        }
    }
    
    int32_t limits[2] = { 100000, -100000 };
    for ( Node* fwd : getDrxnNodes( drxn, true, true ) )
    {
        if ( fwd->drxn_ > 2 )
        {
            fwd->drxn_ = drxn;
            ev.island.erase( remove( ev.island.begin(), ev.island.end(), fwd ), ev.island.end() );
            ev.nodes.push_back( fwd );
            fwd->clearPairs();
            fwd->resetMarks();
            fwd->reliable_ = false;
            fwd->unreliable_ = false;
            fwd->validated_ = false;
            fwd->validLimits_[0] = fwd->validLimits_[1] = fwd->validLimits_[2] = fwd->validLimits_[3] = fwd->ends_[!drxn];
            fwd->stop_[0] = 0;
            fwd->stop_[1] = 0;
            limits[0] = min( limits[0], fwd->ends_[0] );
            limits[1] = max( limits[1], fwd->ends_[1] );
        }
    }
}

void Node::pushValidLimits( IslandVars &iv, Node* markNode, Node* hitNode, int32_t &markCoord, Coords* coords )
{
    markNode->pushValidLimits( markCoord, iv.drxn );
    if ( hitNode->drxn_ <= 2 )
    {
        validLimits_[0] = min( validLimits_[0], markCoord );
        validLimits_[1] = min( validLimits_[1], markCoord );
        validLimits_[2] = max( validLimits_[2], markCoord );
        validLimits_[3] = max( validLimits_[3], markCoord );
    }
    else
    {
        hitNode->pushValidLimits( (*coords)[iv.drxn], !iv.drxn );
    }
    
    if ( markNode != this )
    {
        pushValidLimits( ends_[iv.drxn], iv.drxn );
        markNode->pushValidLimits( markNode->ends_[!iv.drxn], !iv.drxn );
    }
    if ( hitNode != this && hitNode->drxn_ > 2 )
    {
        pushValidLimits( ends_[!iv.drxn], !iv.drxn );
        hitNode->pushValidLimits( hitNode->ends_[iv.drxn], iv.drxn );
    }
}

bool Node::overlapExtend( NodeList &nodes, int32_t* coords, NodeList &hitNodes, vector<int32_t>* hitCoords, bool subGraph, bool drxn )
{
    NodeList currNodes = { this };
    vector<int32_t> currCoords[2];
    currCoords[0].push_back( coords[0] );
    currCoords[1].push_back( coords[1] );
    NodeSet usedSet;
    
    while ( !currNodes.empty() )
    {
        NodeList nxtNodes;
        vector<int32_t> nxtCoords[2];
        NodeSet currFwd;
        for ( Node* const &curr : currNodes )
        {
            curr->getDrxnNodes( currFwd, drxn );
        }
        for ( int i = 0; i < currNodes.size(); i++ )
        {
            if ( currCoords[0][i] >= currCoords[1][i] ) continue;
            if ( usedSet.find( currNodes[i] ) != usedSet.end() ) continue;
            if ( currFwd.find( currNodes[i] ) == currFwd.end() )
            {
                int32_t nxtCoord[2] = { currCoords[0][i], currCoords[1][i] };
                if ( currNodes[i]->getNextReadCoord( nxtCoord[!drxn], !drxn, drxn ) )
                {
                    if ( nxtCoord[!drxn] != currNodes[i]->ends_[!drxn] )
                    {
                        if ( drxn == subGraph )
                        {
                            currNodes[i] = currNodes[i]->splitNode( nxtCoord[!drxn], nodes, drxn, drxn );
                        }
                        else
                        {
                            int32_t splitCoord = currNodes[i]->ends_[!drxn];
                            for ( auto &read : currNodes[i]->reads_ )
                            {
                                if ( drxn ? read.second[0] < nxtCoord[0] : nxtCoord[1] < read.second[1] )
                                {
                                    splitCoord = drxn ? max( splitCoord, read.second[1] ) : min( splitCoord, read.second[0] );
                                }
                            }
                            currNodes[i]->splitNode( splitCoord, nodes, !drxn, !drxn );
                        }
                    }
                    hitNodes.push_back( currNodes[i] );
                    usedSet.insert( currNodes[i] );
                    currNodes[i]->getDrxnNodes( usedSet, drxn );
                    
                    assert( currNodes[i]->ends_[0] <= nxtCoord[0] );
                    assert( currNodes[i]->ends_[1] >= nxtCoord[1] );
                    hitCoords[0].push_back( nxtCoord[0] );
                    hitCoords[1].push_back( nxtCoord[1] );
                }
                else
                {
                    for ( Edge &e : currNodes[i]->edges_[drxn] )
                    {
                        if ( currNodes[i]->drxn_ == 2 ) continue;
                        int32_t diff = drxn ? e.node->ends_[0] - currNodes[i]->ends_[1] + e.overlap
                                            : e.node->ends_[1] - currNodes[i]->ends_[0] - e.overlap;
                        nxtNodes.push_back( e.node );
                        nxtCoords[0].push_back( currCoords[0][i] + diff );
                        nxtCoords[1].push_back( currCoords[1][i] + diff );
                    }
                }
            }
            else
            {
                nxtNodes.push_back( currNodes[i] );
                nxtCoords[0].push_back( currCoords[0][i] );
                nxtCoords[1].push_back( currCoords[1][i] );
            }
        }
        
        currNodes = nxtNodes;
        currCoords[0] = nxtCoords[0];
        currCoords[1] = nxtCoords[1];
    }
    
    return !hitNodes.empty();
}

void Node::reviewMerged( ExtVars &ev, NodeSet &mergeSet, bool drxn )
{
    ev.rebranch.clear();
    ev.cloneSet.clear();
    
    Node::offsetForward( mergeSet, drxn, false );
    NodeSet currSet = Node::getNotForwardSet( mergeSet, drxn );
    
    bool didAdvance = false;
    while ( !currSet.empty() )
    {
        didAdvance = false;
        NodeSet currFwd, nxtSet;
        for ( Node* curr : currSet )
        {
            curr->getDrxnNodes( currFwd, drxn );
        }

        for ( Node* curr : currSet )
        {
            if ( currFwd.find( curr ) == currFwd.end() )
            {
                while ( curr->resolveBypass( ev, false, drxn ) )
                {
                    Node::offsetForward( currSet, drxn, false );
                }
                
                if ( curr->resolveOffset( ev, false, drxn ) )
                {
                    Node::offsetForward( currSet, drxn, false );
                }
                
                didAdvance = true;
                curr->getNextNodes( nxtSet, drxn );
            }
            else
            {
                nxtSet.insert( curr );
            }
        }

        currSet = nxtSet;
        assert( didAdvance );
    }
    
    NodeSet testedSet, backedSet;
    for ( Node* node : mergeSet )
    {
        if ( testedSet.find( node ) == testedSet.end() && backedSet.find( node ) == backedSet.end() )
        {
            node->propagateValidation( ev.limits, testedSet, backedSet, drxn );
        }
    }
    
    ev.offset.clear();
    ev.bypass.clear();
}

bool Node::seedIslandsCheckRead( IslandVars &iv, ReadMark &mark )
{
    for ( Node* node : iv.ev.island )
    {
        if ( node->reads_.find( mark.id ) != node->reads_.end() )
        {
            return false;
        }
    }
    return true;
}

//bool Node::seedIslandBridge( IslandVars &iv, vector<IslandRead*> &path, Node* hitNode, int32_t* coords, bool drxn )
//{
//    NodeList hitNodes;
//    vector<int32_t> hitCoords[2];
//    Node::overlapExtend( hitNode, coords, hitNodes, hitCoords, drxn );
//    if ( !hitNodes.empty() )
//    {
//        string seq = path[0]->seq;
//        vector<ReadId> seqIds;
//        vector<int32_t> seqCoords[2];
//        for ( int i = 1; i < path.size(); i++ )
//        {
//            seq += path[i]->seq.substr( path[i]->overlaps[0] );
//        }
//        
//        for ( int i = 0; i < hitNodes.size(); i++ )
//        {
//            int32_t coord = hitCoords[drxn][i];
//            hitNodes[i]->getNextReadCoord( coord, drxn, drxn );
//            string extra = drxn ? hitNodes[i]->seq_.substr( hitCoords[1][i] - hitNodes[i]->ends_[0], coord - hitCoords[1][i] )
//                                : hitNodes[i]->seq_.substr( coord - hitNodes[i]->ends_[0], hitCoords[0][i] - coord ) ;
//            extra = extra.empty() ? extra : ( drxn ? extra.substr( 0, extra.length() - 1 ) : extra.substr( 1 ) );
//            hitCoords[drxn][i] = hitCoords[drxn][i] + ( drxn ? extra.size() : -extra.size() );
//            if ( hitNodes.size() == 1 )
//            {
//                seq = drxn ? seq + extra : extra + seq;
//                iv.ev.bwt.mapSequence( seq, seqIds, seqCoords );
//                Node* node = new Node( seq, seqIds, seqCoords, hitCoords[drxn][i], iv.drxn + 3, drxn );
//                hitNodes[i]->addEdge( node, abs( node->ends_[drxn] - hitCoords[!drxn][i] ), !drxn, false, false );
//                iv.ev.island.push_back( node );
//                node->setCoverage();
//            }
//            else
//            {
//                string branch = drxn ? path.back()->seq + extra : extra + path[0]->seq;
//                branch = branch.substr( 1, branch.length() - 2 );
//                vector<ReadId> branchIds;
//                vector<int32_t> branchCoords[2];
//                iv.ev.bwt.mapSequence( branch, branchIds, branchCoords );
//                assert( false );
//            }
//        }
//        for ( Node* node : iv.ev.nodes )
//        {
//            for ( ReadId &id : seqIds )
//            {
//                assert( node->reads_.find( id ) == node->reads_.end() );
//            }
//        }
//        return true;
//    }
//    return false;
//}

void Node::seedIslandsClump( IslandVars &iv, vector<ReadMark> &marks, unordered_set<SeqNum> &seeds, bool drxn )
{
    for ( ReadMark &mark : marks )
    {
        if ( !Node::seedIslandsCheckRead( iv, mark ) ) continue;
        
        string seq = iv.ev.bwt.getSequence( mark.id );
        vector<Extension> exts = iv.ev.bwt.mapExtensions( seq, !drxn, seeds );
        
        for ( Extension &ext : exts )
        {
            MergeHit merge;
            checkExtension( ext, iv, merge );
            
            if ( !ext.overlaps.empty() )
            {
                string extSeq = ( !drxn ? seq.substr( seq.length() - ext.maxOverLen ) : seq.substr( 0, ext.maxOverLen ) );
                Node* node = new Node( extSeq, mark, ext, !drxn, drxn );
                iv.ev.island.push_back( node );
                
                if ( merge.node && merge.node->drxn_ == !iv.drxn )
                {
                    Node::islandDelete( iv, node );
                    iv.ev.del.erase( node );
                }
                else if ( merge.node && merge.node->drxn_ <= 2 )
                {
                    Node* mergeNode = merge.node->mergeNode( iv.ev.nodes, merge.coords, iv.drxn, !drxn );
                    mergeNode->addEdge( node, merge.overlap, drxn );
                    iv.merged[!drxn].insert( node );
                }
                else if ( merge.node )
                {
                    Node* mergeNode = merge.node->mergeNode( iv.ev.island, merge.coords, iv.drxn, !drxn );
                    mergeNode->addEdge( node, merge.overlap, drxn );
                }
                else if ( !node->seedIslandsConfirm( iv, seeds, !drxn ) )
                {
                    Node::islandDelete( iv, node );
                    iv.ev.del.erase( node );
                }
            }
        }
    }
}

bool Node::seedIslandsConfirm( IslandVars &iv, unordered_set<SeqNum> &seeds, bool drxn )
{
    bool didExtend = true;
    while ( didExtend )
    {
        didExtend = false;
        NodeSet nodes = { this };
        getDrxnNodes( nodes, drxn );
        for ( Node* node : nodes )
        {
            for ( const SeqNum &readId : seeds )
            {
                if ( node->reads_.find( readId ) != node->reads_.end() )
                {
                    return true;
                }
            }
            if ( node->isContinue( drxn ) )
            {
                iv.ev.ante = node->getDrxnNodes( !drxn, true, true );
                vector<Extension> exts = iv.ev.bwt.mapExtensions( node->seq_, drxn, seeds, ( params.readLen * 2 ) / 3 );
                didExtend = didExtend || !exts.empty();
                node->addExtensions( exts, iv, drxn );
            }
        }
    }
    
    return false;
}

//void Node::seedIslandsClumps( IslandVars &iv, vector<ReadMark> &clumps )
//{
//    vector<IslandRead> reads;
//    for ( ReadMark &mark : clumps )
//    {
//        reads.push_back( IslandRead( iv.ev.bwt.getSequence( mark.id ), mark ) );
//        for ( int i = 0; i < reads.size() - 1; i++ )
//        {
//            if ( reads[i].seq == reads.back().seq )
//            {
//                reads.pop_back();
//                break;
//            }
//        }
//    }
//    
//    for ( int i = 0; i < reads.size(); i++ )
//    {
//        for ( int j = 0; j < reads.size(); j++ )
//        {
//            if ( i != j ) reads[i].setOverlap( reads[j] );
//        }
//    }
//    
//    unordered_set<IslandRead*> forks[2], branches[2], ends;
//    for ( IslandRead &read : reads )
//    {
//        if ( !read.edges[0] && !read.edges[1] ) continue;
//        if ( !read.edges[0] ) ends.insert( &read );
//        for ( int drxn : { 0, 1 } )
//        {
//            if ( read.edges[drxn] && read.edges[drxn]->edges[!drxn] != &read )
//            {
//                branches[drxn].insert( &read );
//                branches[drxn].insert( read.edges[drxn]->edges[!drxn] );
//                forks[!drxn].insert( read.edges[drxn] );
//            }
//        }
//    }
//    
//    vector< vector<IslandRead*> > paths;
//    unordered_set<IslandRead*> pathed;
//    for ( const unordered_set<IslandRead*> &readList : { ends, forks[0], branches[0] } )
//    {
//        for ( IslandRead* read : readList )
//        {
//            if ( pathed.find( read ) != pathed.end() ) continue;
//            vector<IslandRead*> path = { read };
//            while ( path.back()->edges[1] 
//                    && find( path.begin(), path.end(), path.back()->edges[1] ) == path.end()
//                    && forks[1].find( path.back() ) == forks[1].end()
//                    && branches[1].find( path.back() ) == branches[1].end() )
//            {
//                path.push_back( path.back()->edges[1] );
//            }
//            paths.push_back( path );
//            pathed.insert( path.begin(), path.end() );
//        }
//    }
//    
//    vector<MapNode*> mns, mes;
//    for ( vector<IslandRead*> &path : paths )
//    {
//        mns.push_back( new MapNode() );
//        mns.back()->seq = path[0]->seq;
//        int32_t sumCoord = path[0]->coord;
//        for ( int i = 1; i < path.size(); i++ )
//        {
//            mns.back()->seq += path[i]->seq.substr( path[i]->overlaps[0] );
//            sumCoord += path[i]->coord;
//        }
//        iv.ev.bwt.mapSequence( mns.back()->seq, mns.back()->ids, mns.back()->coords );
//        mns.back()->estimate = sumCoord / (int)path.size();
//        Node::seedIslandsClumpsCheck( iv, mns, mns.back() );
//    }
//    
//    for ( int i = 0; i < paths.size(); i++ )
//    {
//        if ( paths[i][0]->edges[0] )
//        {
//            mes.push_back( new MapNode() );
//            mns[i]->addEdge( mes.back(), paths[i][0]->overlaps[0], 0 );
//            int32_t sumCoord = mns[i]->estimate, coordCount = 1;
//            for ( int j = 0; j < paths.size(); j++ )
//            {
//                if ( paths[j].back() == paths[i][0]->edges[0] )
//                {
//                    if ( paths[j].back()->edges[1] == paths[i][0] ) paths[j].back()->edges[1] = NULL;
//                    mns[j]->addEdge( mes.back(), paths[i][0]->overlaps[0], 1 );
//                    sumCoord += mns[j]->estimate;
//                    coordCount++;
//                }
//            }
//            mes.back()->estimate = sumCoord / coordCount;
//            assert( !mes.back()->edges[0].empty() && !mes.back()->edges[1].empty() );
//        }
//        if ( paths[i].back()->edges[1] )
//        {
//            mes.push_back( new MapNode() );
//            mns[i]->addEdge( mes.back(), paths[i].back()->overlaps[1], 1 );
//            int32_t sumCoord = mns[i]->estimate, coordCount = 1;
//            for ( int j = 0; j < paths.size(); j++ )
//            {
//                if ( paths[j][0] == paths[i].back()->edges[1] )
//                {
//                    if ( paths[j][0]->edges[0] == paths[i].back() ) paths[j][0]->edges[0] = NULL;
//                    mns[j]->addEdge( mes.back(), paths[i].back()->overlaps[1], 0 );
//                    sumCoord += mns[j]->estimate;
//                    coordCount++;
//                }
//            }
//            mes.back()->estimate = sumCoord / coordCount;
//            assert( !mes.back()->edges[0].empty() && !mes.back()->edges[1].empty() );
//        }
//    }
//    
//    for ( int i = 0; i < mes.size(); )
//    {
//        if ( !mes[i]->edges[0].empty() && !mes[i]->edges[1].empty() )
//        {
//            mes[i]->setEdgeSeq();
//            iv.ev.bwt.mapSequence( mes[i]->seq, mes[i]->ids, mes[i]->coords );
//            if ( !mes[i]->recoil() || !Node::seedIslandsClumpsCheck( iv, mes, mes[i] ) )
//            {
//                delete mes[i];
//                mes.erase( mes.begin() + i );
//                continue;
//            }
//        }
//        i++;
//    }
//    
//    MapNode::collapse( mns, mes );
//    
//    for ( int i = 0; i < mns.size(); i++ )
//    {
//        if ( mns[i]->edges[!iv.drxn].empty() && mns[i]->bridges[!iv.drxn].empty() )
//        {
//            Node* hitNode = NULL;
//            NodeList hitNodes;
//            int32_t coords[2];
//            vector<int32_t> hitCoords[2];
//            string secondSeq;
//            Node::findOverlap( hitNode, coords, mns[i]->seq, iv.ev.nodes, 16, !iv.drxn );
//            if ( !hitNode && mns[i]->getSecondSeq( secondSeq, !iv.drxn ) )
//            {
//                Node::findOverlap( hitNode, coords, secondSeq, iv.ev.nodes, 32, !iv.drxn );
//                if ( hitNode ) mns[i]->setSecondSeq( !iv.drxn );
//            }
//            if ( hitNode && hitNode->overlapExtend( iv.ev.nodes, coords, hitNodes, hitCoords, iv.drxn, !iv.drxn ) )
//            {
//                string seq = iv.drxn ? mns[i]->seq.substr( 0, *min_element( mns[i]->coords[1].begin(), mns[i]->coords[1].end() ) - 1 )
//                                     : mns[i]->seq.substr( *max_element( mns[i]->coords[0].begin(), mns[i]->coords[0].end() ) + 1 );
//                for ( int j = 0; j < hitNodes.size(); j++ )
//                {
//                    int32_t coord = hitCoords[!iv.drxn][j];
//                    hitNodes[j]->getNextReadCoord( coord, !iv.drxn, !iv.drxn );
//                    int32_t extra = abs( coord - hitCoords[!iv.drxn][j] ) - 1;
//                    assert( extra >= 0 );
//                    MapNode* mn = new MapNode();
//                    mn->seq = iv.drxn ? hitNodes[j]->seq_.substr( coord - hitNodes[j]->ends_[0], extra ) + seq 
//                                      : seq + hitNodes[j]->seq_.substr( hitCoords[1][j] - hitNodes[j]->ends_[0], extra );
//                    hitCoords[!iv.drxn][j] = hitCoords[!iv.drxn][j] + ( iv.drxn ? -extra : extra );
//                    iv.ev.bwt.mapSequence( mn->seq, mn->ids, mn->coords );
//                    mn->addEdge( mns[i], seq.length(), iv.drxn );
//                    mn->addEdge( hitNodes[j], hitCoords[1][j]- hitCoords[0][j], !iv.drxn );
//                    if ( mn->recoil() )
//                    {
//                        mes.push_back( mn );
//                        if ( hitNodes.size() == 1 )
//                        {
//                            MapNode::fold( mn, mes, iv.drxn );
//                        }
//                    }
//                    else
//                    {
//                        delete mn;
//                    }
//                }
//                MapNode::collapse( mns, mes );
//            }
//        }
//        
//        mns[i]->checkLoop();
//        mns[i]->node = new Node( mns[i], 3 + iv.drxn );
//        iv.ev.island.push_back( mns[i]->node );
//        for ( int j = 0; j < mns[i]->bridges[!iv.drxn].size(); j++ )
//        {
//            mns[i]->bridges[!iv.drxn][j]->addEdge( mns[i]->node, mns[i]->bridgeOverlaps[!iv.drxn][j], iv.drxn );
//            iv.merged[!iv.drxn].insert( mns[i]->node );
//        }
//        
//        int32_t nodeLimits[2] = { mns[i]->node->ends_[1], mns[i]->node->ends_[0] };
//        for ( auto &read : mns[i]->node->reads_ )
//        {
//            nodeLimits[0] = min( nodeLimits[0], read.second[0] );
//            nodeLimits[1] = max( nodeLimits[1], read.second[1] );
//        }
//        assert( nodeLimits[0] == mns[i]->node->ends_[0] );
//        assert( nodeLimits[1] == mns[i]->node->ends_[1] );
//        assert( nodeLimits[1] - nodeLimits[0] == mns[i]->node->seq_.length() );
//    }
//    
//    for ( int i = 0; i < mns.size(); i++ )
//    {
//        for ( int j = 0; j < mns[i]->edges[iv.drxn].size(); j++ )
//        {
//            NodeSet fwdSet = mns[i]->edges[iv.drxn][j]->node->getDrxnNodes( iv.drxn );
//            if ( fwdSet.find( mns[i]->node ) == fwdSet.end() 
//                    && mns[i]->node != mns[i]->edges[iv.drxn][j]->node )
//            {
//                mns[i]->node->addEdge( mns[i]->edges[iv.drxn][j]->node, mns[i]->edgeOverlaps[iv.drxn][j], iv.drxn );
//            }
//        }
//        mns[i]->node->setCoverage();
//    }
//}

bool Node::seedIslandsClumpsCheck( IslandVars &iv, vector<MapNode*> &mns, MapNode* mn )
{
    int iBegin = -1, iEnd = mn->ids.size();
    Node* edgeNode = NULL;
    Coords* edgeCoords = NULL;
    for ( int i = 0; i < mn->ids.size(); i++ )
    {
        for ( NodeList const &nodes : { iv.ev.nodes, iv.ev.island } )
        {
            for ( Node* const &node : nodes )
            {
                auto it = node->reads_.find( mn->ids[i] );
                if ( it != node->reads_.end() )
                {
                    if ( iBegin <= 0 )
                    {
                        edgeNode = node;
                        edgeCoords = &it->second;
                    }
                    if ( iBegin == -1 ) iBegin = i;
                    iEnd = i;
                    break;
                }
            }
            if ( iEnd == i ) break;
        }
        if ( edgeNode && iEnd < i ) break;
    }
    if ( edgeNode )
    {
        if ( iBegin )
        {
            if ( iv.drxn ) edgeNode = edgeNode->splitNode( (*edgeCoords)[0], ( edgeNode->drxn_ <= 2 ? iv.ev.nodes : iv.ev.island ), 1, 1 );
            else edgeNode->mergeNode( ( edgeNode->drxn_ <= 2 ? iv.ev.nodes : iv.ev.island ), edgeCoords, 0, 1 );
            int overlap = mn->coords[1][iBegin-1] - mn->coords[0][iBegin];
            if ( iEnd < mn->ids.size() - 1 )
            {
                MapNode* newMn = new MapNode();
                newMn->seq = mn->seq;
                newMn->ids.insert( newMn->ids.end(), mn->ids.begin() + iEnd, mn->ids.end() );
                newMn->coords[0].insert( newMn->coords[0].end(), mn->coords[0].begin() + iEnd, mn->coords[0].end() );
                newMn->coords[1].insert( newMn->coords[1].end(), mn->coords[1].begin() + iEnd, mn->coords[1].end() );
                newMn->recoil();
                newMn->bridges[1].insert( newMn->bridges[1].end(), mn->bridges[1].begin(), mn->bridges[1].end() );
                newMn->bridgeOverlaps[1].insert( newMn->bridgeOverlaps[1].end(), mn->bridgeOverlaps[1].begin(), mn->bridgeOverlaps[1].end() );
                for ( int i = 0; i < mn->edges[1].size(); i++ )
                {
                    newMn->addEdge( mn->edges[1][i], mn->edgeOverlaps[1][i], 1 );
                }
                mns.push_back( newMn );
                Node::seedIslandsClumpsCheck( iv, mns, newMn );
            }
            mn->ids.erase( mn->ids.begin() + iBegin, mn->ids.end() );
            mn->coords[0].erase( mn->coords[0].begin() + iBegin, mn->coords[0].end() );
            mn->coords[1].erase( mn->coords[1].begin() + iBegin, mn->coords[1].end() );
            mn->recoil();
            mn->removeEdges( 1 );
            mn->addEdge( edgeNode, overlap, 1 );
        }
        else if ( iEnd < mn->ids.size() - 1 )
        {
            if ( iv.drxn ) edgeNode->mergeNode( ( edgeNode->drxn_ <= 2 ? iv.ev.nodes : iv.ev.island ), edgeCoords, 1, 0 );
            else edgeNode = edgeNode->splitNode( (*edgeCoords)[1], ( edgeNode->drxn_ <= 2 ? iv.ev.nodes : iv.ev.island ), 0, 0 );
            int overlap = mn->coords[1][iEnd] - mn->coords[0][iEnd+1];
            mn->ids.erase( mn->ids.begin(), mn->ids.begin()+iEnd+1 );
            mn->coords[0].erase( mn->coords[0].begin(), mn->coords[0].begin()+iEnd+1 );
            mn->coords[1].erase( mn->coords[1].begin(), mn->coords[1].begin()+iEnd+1 );
            mn->recoil();
            mn->removeEdges( 0 );
            mn->addEdge( edgeNode, overlap, 0 );
        }
        else
        {
            mn->removeEdges( 0 );
            mn->removeEdges( 1 );
            return false;
        }
    }
    return true;
}

bool Node::seedIslandsSingle( IslandVars &iv, ReadMark &mark, unordered_set<SeqNum> &seeds, bool drxn )
{
    if ( !Node::seedIslandsCheckRead( iv, mark ) ) return false;
    
    string seq = iv.ev.bwt.getSequence( mark.id );
    vector<Extension> exts = iv.ev.bwt.mapExtensions( seq, drxn, ( params.readLen * 2 ) / 3 );

    for ( Extension &ext : exts )
    {
        MergeHit merge;
        checkExtension( ext, iv, merge );

        if ( !ext.overlaps.empty() && ( !merge.node || merge.node->drxn_ > 2 ) )
        {
            string extSeq = ( drxn ? seq.substr( seq.length() - ext.maxOverLen ) : seq.substr( 0, ext.maxOverLen ) );
            Node* node = new Node( extSeq, mark, ext, drxn, drxn );
            iv.ev.island.push_back( node );

            // Merge with overlapping island node if applicable
            if ( merge.node )
            {
                Node* mergeNode = merge.node->splitNode( (*merge.coords)[!drxn], iv.ev.island, iv.drxn, drxn );
                node->addEdge( mergeNode, merge.overlap, drxn );
                return true;
            }

            // Ensure that node is extendable in both directions
            else
            {
                if ( node->seedIslandsConfirm( iv, seeds, !drxn ) )
                {
                    node->extendCount_ = 10;
                    node->extendIsland( iv, 0 );
                    node->extendCount_ = 10;
                    node->extendIsland( iv, 1 );
                    return true;
                }
                else
                {
                    Node::islandDelete( iv, node );
                    iv.ev.del.erase( node );
                    return false;
                }
            }
        }
    }
}

bool Node::setBlank( IslandVars &iv, NodeSet &foldable, bool drxn )
{
    sort( edges_[!drxn].begin(), edges_[!drxn].end(), []( Edge &a, Edge &b ){ 
        return a.overlap > b.overlap;
    } );
    
    seq_ = ( drxn ? edges_[!drxn][0].node->seq_.substr( edges_[!drxn][0].node->seq_.length() - edges_[!drxn][0].overlap ) 
                  : edges_[!drxn][0].node->seq_.substr( 0, edges_[!drxn][0].overlap ) );
    ends_[drxn] = ends_[!drxn] + ( drxn ? seq_.length() : -seq_.length() );
    
    NodeSet folded;
    NodeIntMap offMap;
    
    vector<SeqNum> readIds;
    vector<string> seqs;
    vector<int32_t> coords;
    
    for ( Edge &re : edges_[!drxn] )
    {
        NodeSet currSet;
        int32_t reOff = edges_[!drxn][0].overlap - re.overlap;
        for ( Edge &e : re.node->edges_[drxn] )
        {
            if ( folded.find( e.node ) == folded.end()
                    && foldable.find( e.node ) != foldable.end()
                    && e.overlap <= re.overlap )
            {
                offMap[e.node] = reOff + e.overlap;
                currSet.insert( e.node );
            }
        }
        
        while ( !currSet.empty() )
        {
            NodeSet nxtSet;
            for ( Node* curr : currSet )
            {
                int32_t currOff = offMap[curr];
                folded.insert( curr );
                assert( curr != this );
                for ( auto &read : curr->reads_ )
                {
                    int32_t readBgnOff = currOff - abs( read.second[!drxn] - curr->ends_[!drxn] );
                    int32_t readEndOff = abs( read.second[drxn] - curr->ends_[!drxn] );
                    if ( readBgnOff > 0 && readBgnOff < curr->seq_.length() )
                    {
                        string readSeq = drxn ? curr->seq_.substr( currOff, readEndOff - currOff )
                                              : curr->seq_.substr( curr->seq_.length() - readEndOff, readEndOff - currOff );
                        int32_t coord = drxn ? ends_[0] + ( edges_[!drxn][0].overlap - readBgnOff )
                                             : ends_[1] - ( edges_[!drxn][0].overlap - readBgnOff );
                        
                        if ( !drxn )
                        {
                            reverse( readSeq.begin(), readSeq.end() );
                        }
                        
                        readIds.push_back( read.first );
                        seqs.push_back( readSeq );
                        coords.push_back( coord );
                    }
                }
                for ( Edge &e : curr->edges_[drxn] )
                {
                    int32_t eOff = currOff + e.overlap - ( e.node->ends_[1] - e.node->ends_[0] );
                    if ( folded.find( e.node ) == folded.end()
                            && foldable.find( e.node ) != foldable.end()
                            && eOff >= 0 )
                    {
                        offMap[e.node] = eOff;
                        nxtSet.insert( e.node );
                    }
                }
            }
            currSet = nxtSet;
        }
    }
    
    string bestString;
    int bestScore = 0;
    
    for ( auto it = seqs.begin(); it != seqs.end() && it+1 != seqs.end(); it++ )
    {
        vector<string*> pSeqs;
        for ( auto it2 = it + 1; it2 != seqs.end(); it2++ )
        {
            pSeqs.push_back( &(*it2) );
        }
        
        int hits = 0, swaps = 1, len = 0;
        for ( int i( 0 ); i < it->length() && !pSeqs.empty(); i++ )
        {
            for ( auto it2 = pSeqs.begin(); it2 != pSeqs.end(); )
            {
                if ( (**it2)[i] != (*it)[i] )
                {
                    it2 = pSeqs.erase( it2 );
                    continue;
                }
                it2++;
            }
            if ( !pSeqs.empty() )
            {
                hits += pSeqs.size();
                swaps += i > 0 && (*it)[i-1] != (*it)[i];
                len = i;
            }
        }
        
        int score = swaps * hits;
        if ( score > bestScore )
        {
            bestScore = score;
            bestString = it->substr( 0, len );
        }
    }
    
    if ( seqs.size() == 1 )
    {
        bestString = seqs[0];
    }
    
    for ( int i( 0 ); i < seqs.size(); i++ )
    {
        int len = 0;
        while ( seqs[i][len] == bestString[len] 
                && len < seqs[i].length() 
                && len < bestString.length() )
        {
            len++;
        }
        int32_t coordsBgn = drxn ? coords[i] : ends_[0] - len;
        int32_t coordsEnd = drxn ? ends_[1] + len : coords[i];
        reads_.insert( make_pair( readIds[i], Coords( coordsBgn, coordsEnd, len != seqs[i].length() ) ) );
    }
    
    if ( reads_.empty() )
    {
        folded.clear();
        return false;
    }
    
    if ( !drxn )
    {
        reverse( bestString.begin(), bestString.end() );
    }
    seq_ = drxn ? seq_ + bestString : bestString + seq_;
    ends_[drxn] += ( drxn ? bestString.length() : -bestString.length() );
    
    NodeSet tSet;
    
    assert( folded.find( this ) == folded.end() );
    
    for ( Node* node : folded )
    {
        for ( auto &np : node->pairs_ )
        {
            np.first->pairs_.erase( node );
            np.first->resetFurthest( node );
            tSet.insert( np.first );
        }
        node->pairs_.clear();
        node->dismantleNode();
        iv.ev.del.insert( node );
    }
    
    for ( bool d : { 0, 1 } )
    {
        vector<ReadMark> marks = getMarksBase( d );
        for ( Node* t : tSet )
        {
            for ( auto it = marks.begin(); it != marks.end(); )
            {
                if ( t->reads_.find( it->id ) != t->reads_.end() )
                {
                    if ( d == iv.drxn )
                    {
                        pushValidLimits( it->mark, d );
                    }
                    
                    auto r = t->pairs_.insert( make_pair( this, 1 ) );
                    if ( !r.second )
                    {
                        r.first->second++;
                    }
                    r = pairs_.insert( make_pair( t, 1 ) );
                    if ( !r.second )
                    {
                        r.first->second++;
                    }
                    SeqNum pairId = params.getPairId( it->id );
                    islandUpdateValid( pairId );
                    
                    it = marks.erase( it );
                    continue;
                }
                it++;
            }
        }
        
        marks_[d] = marks;
    }
    
    validLimits_[0] = validLimits_[1] = ends_[0];
    validLimits_[2] = validLimits_[3] = ends_[1];
    
    return true;
}
