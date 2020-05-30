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

#include "locus.h"
#include "locus_pathing_structs.h"
#include <algorithm>

void Locus::resetIslandVars( bool drxn )
{
    ivs_[drxn]->merged[0].clear();
    ivs_[drxn]->merged[1].clear();
    for ( unordered_set<SeqNum>* readIds : { &ivs_[drxn]->peReads, &ivs_[drxn]->mpReads, &ivs_[drxn]->peReadsReliable } )
    {
        for ( auto it = readIds->begin(); it != readIds->end(); )
        {
            bool doErase = true;;
            for ( Node* node : nodes_[3 + drxn] )
            {
                if ( node->reads_.find( *it ) != node->reads_.end() )
                {
                    doErase = false;
                    break;
                }
            }
            if ( doErase )
            {
                it = readIds->erase( it );
                continue;
            }
            it++;
        }
    }
}

void Locus::setIslandVars( Node* pathEnd, bool drxn )
{
    if ( !ivs_[drxn] )
    {
        ExtVars ev( nodes_[drxn], nodes_[drxn + 3], validLimits_, bwt_ );
        ivs_[drxn] = new IslandVars( ev, drxn );
    }
    ivs_[drxn]->pathEnd = pathEnd;
    ivs_[drxn]->origin.clear();
    ivs_[drxn]->ante.clear();
    ivs_[drxn]->limit = pathEnd->ends_[drxn];
    for ( Node* origin : originEnds_[drxn] )
    {
        ivs_[drxn]->origin.push_back( origin );
        origin->getDrxnNodes( ivs_[drxn]->ante, !drxn );
    }
}

void Locus::leap()
{
    clock_t startTime = clock();
    for ( bool drxn : { 0, 1 } )
    {
        for ( Path &path : paths_[drxn] )
        {
            if ( path.path.empty() || finished_[drxn] || stopCodes_[drxn] ) continue;
            if ( !leap( path, drxn ) ) continue;
            path.completed = false;
            completed_[drxn] = false;
            desperation_[drxn] = false;
            resetIslandVars( drxn );
        }
    }
    leapTime_ += (double)( clock() - startTime ) / (double) CLOCKS_PER_SEC;
}

bool Locus::leap( Path &path, bool drxn )
{
    
    int32_t limits[3];
    NodeList reliableNodes;
    leapGetNodes( path.path.back(), reliableNodes, limits, drxn );
    setIslandVars( path.path.back(), drxn );
    ivs_[drxn]->limit = limits[2];
    
    vector<ReadMark> peClumps, peMarks, mpMarks;
    unordered_set<SeqNum> readIds;
    leapGetReads( reliableNodes, peClumps, peMarks, mpMarks, readIds, limits, drxn );
    
    if ( peClumps.empty() || !leapSetIslandClumps( *ivs_[drxn], peClumps, readIds ) )
    {
        leapSetIslandSingles( *ivs_[drxn], peMarks, mpMarks, readIds, limits );
    }
    
    if ( !leapReview( *ivs_[drxn] ) || params.haploid )
    {
        return leapBridge( *ivs_[drxn] );
    }
    return true;
}

bool Locus::leapBridge( IslandVars &iv )
{
    NodeSet usedSet;
    NodeSetList peIslands, mpIslands;

    for ( Node* node : iv.ev.island )
    {
        if ( usedSet.find( node ) == usedSet.end() )
        {
            NodeSet islandSet = node->getConnectedNodes( true );
            usedSet.insert( islandSet.begin(), islandSet.end() );
            Node::islandReview( iv, islandSet, peIslands, mpIslands );
        }
    }
    
    if ( !peIslands.empty() && Node::bridgeIsland( iv, peIslands ) )
    {
        leapCleanup( iv );
        leapSetExtendMerge( iv );
        return leapReview( iv );
    }
    leapCleanup( iv );
    
    return false;
}

void Locus::leapCleanup( IslandVars &iv )
{
    NodeSet delSet;
    for ( Node* del : iv.ev.del )
    {
        if ( find( iv.ev.island.begin(), iv.ev.island.end(), del ) != iv.ev.island.end() )
        {
            Node::islandDelete( iv, del );
        }
        else if ( find( iv.ev.nodes.begin(), iv.ev.nodes.end(), del ) != iv.ev.nodes.end() )
        {
            delSet.insert( del );
        }
    }
    deleteNodes( delSet, iv.drxn );
    iv.ev.del.clear();
}

void Locus::leapCleanup( IslandVars &iv, NodeSet &goodSet )
{
    for ( Node* del : iv.ev.del )
    {
        if ( goodSet.find( del ) == goodSet.end() )
        {
            Node::islandDelete( iv, del );
        }
    }
    iv.ev.del.clear();
}

bool Locus::leapExtend( IslandVars &iv )
{
    NodeSet extSets[2];
    iv.round = 0;
    
    while ( leapSetExtend( iv, extSets ) )
    {
        for ( bool drxn : { iv.drxn, !iv.drxn } )
        {
            for ( Node* node : extSets[drxn] )
            {
                node->extendCount_ = 20;
                node->extendIsland( iv, drxn );
            }
        }
    }
    
    leapCleanup( iv );
}

void Locus::leapGetNodes( Node* pathEnd, NodeList &reliableNodes, int32_t* limits, bool drxn )
{
    limits[drxn] = pathEnd->ends_[drxn];
    
    for ( Node* fwd : pathEnd->getDrxnNodes( drxn ) )
    {
        if ( fwd->isReliable( false ) )
        {
            limits[drxn] = drxn ? max( limits[1], fwd->ends_[1] ) : min( limits[0], fwd->ends_[0] );
            reliableNodes.push_back( fwd );
        }
    }
    limits[2] = limits[drxn];
    limits[drxn] += drxn ? 5000 : -5000;
    limits[!drxn] = drxn ? max( params.readLen, limits[2] - 2500 ) : min( 0, limits[2] + 2500 );
    
    NodeSet nodes = { pathEnd };
    pathEnd->getDrxnNodes( nodes, !drxn, params.getFurthestMpMean( pathEnd->ends_[drxn], !drxn ) );
    
    for ( Node* node : nodes )
    {
        if ( node->isReliable( false ) )
        {
            reliableNodes.push_back( node );
        }
    }
    
    sort( reliableNodes.begin(), reliableNodes.end(), [&drxn]( Node* &a, Node* &b ){
        return drxn ? a->ends_[1] > b->ends_[1] : a->ends_[0] < b->ends_[0];
    });
}

void Locus::leapGetPeClumps( vector<ReadMark> &peClumps, vector<ReadMark> &peMarks )
{
    int32_t cutoff = params.readLen * 1.5;
    int minCount = max( 3, int( params.peCover / 6 ) );
    vector< pair<int32_t, int32_t> > ranges;
    if ( peMarks.size() > minCount )
    {
        for ( int i( 0 ); i < peMarks.size() - minCount; i++ )
        {
            if ( peMarks[i + minCount].estimate - peMarks[i].estimate <= cutoff )
            {
                if ( !ranges.empty() && peMarks[i].estimate <= ranges.back().second )
                {
                    ranges.back().second = peMarks[i + minCount].estimate + cutoff;
                }
                else
                {
                    ranges.push_back( make_pair( peMarks[i].estimate - cutoff, peMarks[i + minCount].estimate + cutoff ) );
                }
            }
        }
    }
    
    for ( ReadMark &mark : peMarks )
    {
        for ( pair<int32_t, int32_t> &rng : ranges )
        {
            if ( rng.first <= mark.estimate && mark.estimate <= rng.second )
            {
                peClumps.push_back( mark );
                break;
            }
        }
    }
}

void Locus::leapGetReads( NodeList &reliableNodes, vector<ReadMark> &peClumps, vector<ReadMark> &peMarks, vector<ReadMark> &mpMarks, unordered_set<SeqNum> &readIds, int32_t* limits, bool drxn )
{
    NodeList allNodes = getAllNodes();
    allNodes.insert( allNodes.end(), nodes_[3 + drxn].begin(), nodes_[3 + drxn].end() );
    
    for ( Node* node : reliableNodes )
    {
        node->getUnpairedMarks( allNodes, peMarks, mpMarks, leptReads_[drxn], readIds, limits, true, drxn );
    }
    
    sort( peMarks.begin(), peMarks.end(), []( ReadMark &a, ReadMark &b ){
        return a.estimate < b.estimate;
    });
    
    leapGetPeClumps( peClumps, peMarks );
    
    if ( !peClumps.empty() )
    {
        vector<ReadMark> peRevMarks;
        
        for ( Node* node : reliableNodes )
        {
            node->getUnpairedMarks( allNodes, peRevMarks, mpMarks, leptReads_[drxn], readIds, limits, false, !drxn );
        }
        
        sort( peRevMarks.begin(), peRevMarks.end(), []( ReadMark &a, ReadMark &b ){
            return a.estimate < b.estimate;
        });
        
        leapGetPeClumps( peClumps, peRevMarks );
    }
}

bool Locus::leapReview( IslandVars &iv )
{
    for ( Node* merged : iv.merged[iv.drxn] )
    {
        if ( merged->drxn_ > 2 )
        {
            merged->getConnectedNodes( iv.ev.del, true );
        }
        else
        {
            iv.ev.ante.clear();
            merged->rebranchNode( iv.ev, iv.drxn );
        }
    }
    
    leapCleanup( iv );
    
    NodeSet usedSet;

    for ( Node* node : iv.ev.island )
    {
        if ( usedSet.find( node ) == usedSet.end() )
        {
            NodeSet islandSet = node->getConnectedNodes( true );
            usedSet.insert( islandSet.begin(), islandSet.end() );
            
            Node::islandReview( iv, islandSet );
        }
    }
    
    leapCleanup( iv );
    
    leptReads_[iv.drxn].insert( iv.peReads.begin(), iv.peReads.end() );
    leptReads_[iv.drxn].insert( iv.mpReads.begin(), iv.mpReads.end() );
    
    if ( !iv.merged[!iv.drxn].empty() )
    {
        clock_t startTime = clock();
        Node::reviewMerged( iv.ev, iv.merged[!iv.drxn], iv.drxn );
        debriefExtend( iv.ev, iv.drxn, true );
        revTime_ += (double)( clock() - startTime ) / (double) CLOCKS_PER_SEC;
        
        return true;
    }
    
    return false;
}

bool Locus::leapSetExtend( IslandVars &iv, NodeSet* extSets )
{
    extSets[0].clear();
    extSets[1].clear();
    leapCleanup( iv );
    
    // Set merged island nodes and trim hanging nodes
    leapSetExtendMerge( iv );
    
    // Set pairs between island and mainland
    leapSetExtendBridge( iv );
    
    // Set island nodes to extend for this round
    leapSetExtendEnds( iv, extSets );
    
    if ( extSets[0].empty() && extSets[1].empty() )
    {
        return false;
    }
    return true;
}

void Locus::leapSetExtendBridge( IslandVars &iv )
{
    NodeSet tSet = { iv.pathEnd };
    iv.pathEnd->getDrxnNodes( tSet, !iv.drxn, params.getFurthestMpMean( iv.pathEnd->ends_[iv.drxn], !iv.drxn ) );
    iv.pathEnd->getDrxnNodes( tSet, iv.drxn, true );
    NodeList tNodes( tSet.begin(), tSet.end() );
    for ( Node* node : iv.ev.island )
    {
        node->islandPairToBase( iv, tNodes );
    }
}

void Locus::leapSetExtendEnds( IslandVars &iv, NodeSet* &extSets )
{
    NodeSet usedSet;
    
    for ( Node* node : iv.ev.island )
    {
        if ( usedSet.find( node ) == usedSet.end() )
        {
            NodeSet islandSet = node->getConnectedNodes( true );
            usedSet.insert( islandSet.begin(), islandSet.end() );
            Node::islandSetExtend( iv, islandSet, extSets );
        }
    }
    
    leapCleanup( iv );
    
    for ( bool drxn : { 0, 1 } )
    {
        for ( auto it = extSets[drxn].begin(); it != extSets[drxn].end(); )
        {
            if ( find( iv.ev.island.begin(), iv.ev.island.end(), *it ) == iv.ev.island.end() )
            {
                Node* node = *it;
                it = extSets[drxn].erase( it );
                continue;
            }
            it++;
        }
    }
    
    iv.round++;
}

void Locus::leapSetExtendMerge( IslandVars &iv )
{
    NodeSet mergedSet, reachableSet, hangingSet, branchSet, usedSet;
    
    for ( auto it = iv.merged[!iv.drxn].begin(); it != iv.merged[!iv.drxn].end(); )
    {
        if ( (*it)->drxn_ > 2 )
        {
            (*it)->getConnectedNodes( reachableSet, true );
            it = iv.merged[!iv.drxn].erase( it );
            continue;
        }
        it++;
    }
    
    if ( reachableSet.empty() ) return;
    
    for ( Node* node : reachableSet )
    {
        for ( Node* nxt : node->getNextNodes( !iv.drxn ) )
        {
            if ( nxt->drxn_ <= 2 )
            {
                iv.merged[!iv.drxn].insert( node );
                branchSet.insert( node );
                mergedSet.insert( node );
                node->getDrxnNodes( mergedSet, iv.drxn );
            }
        }
    }
    
    for ( Node* node : reachableSet )
    {
        if ( mergedSet.find( node ) == mergedSet.end() && node->removeEdges( mergedSet, iv.drxn ) )
        {
            hangingSet.insert( node );
        }
    }
    
    for ( Node* node : hangingSet )
    {
        bool isMerged = usedSet.find( node ) != usedSet.end();
        
        if ( !isMerged )
        {
            NodeSet islandSet = node->getConnectedNodes( true );
            isMerged = node->islandReview( iv, islandSet );
            ( isMerged ? usedSet : iv.ev.del ).insert( islandSet.begin(), islandSet.end() );
        }
        
        if ( isMerged )
        {
            iv.merged[iv.drxn].insert( node );
        }
    }
    
    for ( Node* node : branchSet )
    {
        node->mergeIsland( iv.ev, iv.drxn );
    }
}

bool Locus::leapSetIslandClumps( IslandVars &iv, vector<ReadMark> &peClumps, unordered_set<SeqNum> &readIds )
{
//    Node::seedIslandsClumps( iv, peClumps );
//    for ( ReadMark &mark : peClumps )
//    {
//        for ( Node* node : iv.ev.island )
//        {
//            if ( node->reads_.find( mark.readId ) != node->reads_.end() )
//            {
//                leptReads_[iv.drxn].insert( mark.readId );
//            }
//        }
//    }
    for ( ReadMark &mark : peClumps )
    {
        Node::seedIslandsSingle( iv, mark, readIds, iv.drxn );
        leptReads_[iv.drxn].insert( mark.id );
    }
    leapCleanup( iv );
    leapExtend( iv );
    
    return !iv.merged[!iv.drxn].empty();
}

void Locus::leapSetIslandSingles( IslandVars &iv, vector<ReadMark> &peMarks, vector<ReadMark> &mpMarks, unordered_set<SeqNum> &readIds, int32_t* limits )
{
    for ( ReadMark &mark : peMarks )
    {
        if ( leptReads_[iv.drxn].find( mark.id ) == leptReads_[iv.drxn].end() 
                && ( iv.drxn ? mark.estimate >= limits[2] 
                             : mark.estimate <= limits[2] ) )
        {
            Node::seedIslandsSingle( iv, mark, readIds, iv.drxn );
        }
        leptReads_[iv.drxn].insert( mark.id );
    }
    
    for ( ReadMark &mark : mpMarks )
    {
        Node::seedIslandsSingle( iv, mark, readIds, iv.drxn );
        leptReads_[iv.drxn].insert( mark.id );
    }
    
    leapCleanup( iv );
    leapExtend( iv );
}

