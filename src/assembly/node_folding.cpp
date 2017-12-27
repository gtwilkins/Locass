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

void Node::bluntEnd( bool drxn )
{
    if ( clones_ ) return;
    
    int32_t minEnd = ends_[!drxn] + ( drxn ? getBestOverlap( 0 ) : -getBestOverlap( 1 ) );
    if ( edges_[!drxn].empty() )
    {
        minEnd = minEnd + ( drxn ? params.readLen : -params.readLen );
    }
    if ( drxn ? ends_[1] <= minEnd : minEnd <= ends_[0] ) return;
    
    int32_t endCoords[3] = { minEnd, minEnd, minEnd };
    
    for ( auto &read : reads_ )
    {
        for ( int i( 0 ); i < 3; i++ )
        {
            if ( drxn ? endCoords[i] <= read.second.coords[1] : read.second.coords[0] <= endCoords[i] )
            {
                for ( int j ( 3 ); i < --j; )
                {
                    endCoords[j] = endCoords[j-1];
                }
                endCoords[i] = read.second.coords[drxn];
                break;
            }
        }
    }
    
    if ( endCoords[2] != ends_[drxn] )
    {
        truncateNode( endCoords[2], drxn );
    }
}

bool Node::foldAlleles( NodeList &nodes, Node* forks[2], NodeList paths[2], NodeSet &delSet, bool drxn )
{
    int reads[2]{0};
    int minOl[2] = { params.readLen, params.readLen };
    int pref;
    bool anyClones = false;
    for ( int i : { 0, 1 } )
    {
        for ( int j = 0; j < paths[i].size(); j++ )
        {
            reads[i] += paths[i][j]->reads_.size();
            minOl[i] = min( minOl[i], paths[i][j]->getBestOverlap( 0 ) );
            minOl[i] = min( minOl[i], paths[i][j]->getBestOverlap( 1 ) );
            anyClones = anyClones || paths[i][j]->clones_;
        }
        if ( reads[i] > reads[!i] ) pref = i;
    }
    
    if ( anyClones || ( minOl[0] < 0 && minOl[1] < 0 ) ) return false;
    
    int bestOls[2] = { max( forks[0]->getOverlap( paths[0][0], 1 )
                          , forks[0]->getOverlap( paths[1][0], 1 ) )
                     , max( forks[1]->getOverlap( paths[0].back(), 0 )
                          , forks[1]->getOverlap( paths[1].back(), 0 ) ) };
    if ( reads[!pref] <= 2 || reads[pref] > reads[!pref] * 5 )
    {
        if ( minOl[pref] < 0 && minOl[!pref] >= 0 )
        {
            pref = !pref;
        }
        Node* node = new Node();
        nodes.push_back( node );
        node->drxn_ = drxn;
        node->seq_ = forks[0]->seq_.substr( forks[0]->seq_.length() - bestOls[0] );
        
        // Get length of both paths and set sequence
        int lens[2]{0};
        for ( int i : { 0, 1 } )
        {
            int len = bestOls[0] - forks[0]->getOverlap( paths[i][0], 1 );
            for ( int j = 0; j < paths[i].size(); j++ )
            {
                Node* prv = j ? paths[i][j-1] : forks[0];
                if ( i == pref ) node->seq_ += paths[i][j]->seq_.substr( prv->getOverlap( paths[i][j], 1 ) );
                len += paths[i][j]->seq_.length();
                if ( j ) len -= prv->getOverlap( paths[i][j], 1 );
            }
            len += bestOls[1] - forks[1]->getOverlap( paths[i].back(), 0 );
            lens[i] = len;
        }
        int lastOl = paths[pref].back()->getOverlap( forks[1], 1 );
        node->seq_ += forks[1]->seq_.substr( lastOl, bestOls[1] - lastOl );
        node->ends_[1] = node->seq_.length();
        
        // Add reads to node
        // Reads are added as offset from their closest path end
        for ( int i : { 0, 1 } )
        {
            int32_t offset = bestOls[0] - forks[0]->getOverlap( paths[i][0], 1 );
            for ( int j = 0; j < paths[i].size(); j++ )
            {
                int32_t fromEnd = lens[i] - offset - paths[i][j]->seq_.length();
                assert( fromEnd >= 0 );
                for ( auto &read : paths[i][j]->reads_ )
                {
                    int32_t fromEnds[2] = { offset + read.second[0] - paths[i][j]->ends_[0]
                                         , paths[i][j]->ends_[1] - read.second[1] + fromEnd };
                    int32_t readCoords[2] = { fromEnds[0], fromEnds[0] + read.second[1] - read.second[0] };
                    if ( fromEnds[1] < fromEnds[0] )
                    {
                        readCoords[1] = node->seq_.length() - fromEnds[1];
                        readCoords[0] = readCoords[1] - ( read.second[1] - read.second[0] );
                    }
                    node->addRead( read.first, readCoords[0], readCoords[1], i != pref );
                    assert( fromEnds[0] >= 0 && fromEnds[1] >= 0 );
                }
                if ( j+1 < paths[i].size() ) offset += paths[i][j]->seq_.length() - paths[i][j]->getOverlap( paths[i][j+1], 1 );
            }
        }
//        int32_t coord = abs( paths[pref][0]->getOverlap( forks[0], 0 ) - bestOls[0] );
//        for ( int i = 0; i < paths[pref].size(); i++ )
//        {
//            if ( !i )
//            {
//                int ol = forks[0]->getOverlap( paths[pref][0], 1 );
//                if ( ol < 0 )
//                {
//                    node->seq_ += string( -ol, 'N' );
//                    ol = 0;
//                }
//                node->seq_ += paths[pref][0]->seq_.substr( ol );
//            }
//            else
//            {
//                coord -= paths[pref][i]->getOverlap( paths[pref][i-1], 0 );
//                node->seq_ += paths[pref][i]->seq_.substr( paths[pref][i-1]->getOverlap( paths[pref][i], 1 ) );
//            }
//            for ( auto &read : paths[pref][i]->reads_ )
//            {
//                int32_t coords[2] = { coord + read.second[0] - paths[pref][i]->ends_[0]
//                                    , coord + read.second[1] - paths[pref][i]->ends_[0] };
//                node->addRead( read.first, coords[0], coords[1], false );
//            }
//            coord += paths[pref][i]->seq_.length();
//        }
        
//        coord = abs( paths[!pref][0]->getOverlap( forks[0], 0 ) - bestOls[0] );
//        for ( int i = 0; i < paths[!pref].size(); i++ )
//        {
//            if ( i ) coord -= paths[!pref][i]->getOverlap( paths[!pref][i-1], 0 );
//            for ( auto &read : paths[!pref][i]->reads_ )
//            {
//                int32_t coords[2] = { coord + read.second[0] - paths[!pref][i]->ends_[0]
//                                    , coord + read.second[1] - paths[!pref][i]->ends_[0] };
//                int diff = coords[1] - node->ends_[1];
//                if ( diff > 0 )
//                {
//                    coords[0] -= diff;
//                    coords[1] -= diff;
//                }
//                node->addRead( read.first, coords[0], coords[1], false );
//            }
//            coord += paths[!pref][i]->seq_.length();
//        }
        paths[0][0]->clearEdges( 0 );
        paths[0].back()->clearEdges( 1 );
        paths[1][0]->clearEdges( 0 );
        paths[1].back()->clearEdges( 1 );
        forks[!drxn]->addEdge( node, bestOls[!drxn], drxn );
        node->addEdge( forks[drxn], bestOls[drxn], drxn );
        for ( int i : { 0, 1 } )
        {
            for ( int j = 0; j < paths[i].size(); j++ )
            {
                paths[i][j]->dismantleNode();
                delSet.insert( paths[i][j] );
            }
        }
        node->setValid();
        return true;
    }
    return false;
}

void Node::foldBranch( NodeList &nodes, Node* merges[2], int32_t coords[2], int ol, bool drxn )
{
    merges[drxn] = merges[drxn]->splitNode( coords[drxn], nodes, drxn, drxn );
    merges[!drxn]->mergeNode( nodes, coords[!drxn], drxn );
    merges[!drxn]->addEdge( merges[drxn], ol, drxn );
}

Node* Node::foldEdge( ExtVars &ev, Node* targetNode, bool drxn )
{
    assert( targetNode->ends_[0] <= targetNode->validLimits_[0] && targetNode->validLimits_[3] <= targetNode->ends_[1] );
    validLimits_[3*(!drxn)] = drxn ? max( min( validLimits_[0], ends_[1] - params.readLen ), ends_[0] )
                                   : min( max( validLimits_[3], ends_[0] + params.readLen ), ends_[1] );
    vector<MapStruct> qStructs = getMapStructQuery( validLimits_[3*(!drxn)], drxn );
    vector<MapStruct> tStructs = getMapStructTarget( targetNode, validLimits_[3*(!drxn)], !drxn );
    MapResult result;
    
    for ( MapStruct &q : qStructs )
    {
        q.isEnd = true;
        for ( MapStruct &t : tStructs )
        {
            t.isEnd = true;
            mapFold( result, q, t, params.readLen, 12, false, drxn );
        }
    }
    
    if ( result.l && result.r )
    {
        return ( drxn ? result.r : result.l )->foldEdgeHit( ev, result, drxn );
    }
    else if ( drxn ? targetNode->validLimits_[3] - validLimits_[0] < params.readLen * 3
                   : validLimits_[3] - targetNode->validLimits_[0] < params.readLen * 3 )
    {
        return foldEdgeMiss( ev, targetNode, validLimits_[(!drxn) * 3], drxn );
    }
    return NULL;
}

Node* Node::foldEdgeHit( ExtVars &ev, MapResult &result, bool drxn )
{
    // Truncate island end if necessary
    assert( foldPrep( ev, result.coords[drxn], result.len, true, !drxn ) );
    
    // Prep and edge target
    if ( ( drxn ? result.l : result.r )->foldPrep( ev, result.coords[!drxn], result.len, true, drxn ) )
    {
        ( drxn ? result.l : result.r )->addEdge( this, result.len, drxn, false, result.len <= 0 );
    }
    else
    {
        for ( Edge &e : ( drxn ? result.l : result.r )->edges_[drxn] )
        {
            int len = result.len - abs( result.coords[!drxn] - ( drxn ? result.l : result.r )->ends_[drxn] ) + e.overlap;
            e.node->addEdge( this, len, drxn, false, len <= 0 || e.isLeap );
        }
    }
    
    return ( edges_[!drxn].empty() ? NULL : this );
}

Node* Node::foldEdgeMiss( ExtVars &ev, Node* targetNode, int32_t limit, bool drxn )
{
    // Prep if necessary
    int dummy = 0;
    assert( foldPrep( ev, ends_[!drxn], dummy, true, !drxn ) );
    if ( drxn ? ends_[0] <= targetNode->ends_[1] : targetNode->ends_[0] <= ends_[1] )
    {
        assert( foldPrep( ev, limit, dummy, true, !drxn ) );
        assert( targetNode->foldPrep( ev, targetNode->validLimits_[drxn*3], dummy, true, drxn ) );
    }
    
    // Blunt if necessary
//    if ( drxn ? ends_[0] <= targetNode->ends_[1] : targetNode->ends_[0] <= ends_[1] )
//    {
//        bluntEnd( !drxn );
//        if ( targetNode->edges_[drxn].empty() )
//        {
//            targetNode->bluntEnd( drxn );
//        }
//    }
    
    // Add leap edge
    if ( drxn ? targetNode->ends_[1] < ends_[0] : ends_[1] < targetNode->ends_[0] )
    {
        int overlap = drxn ? targetNode->ends_[1] - ends_[0] : ends_[1] - targetNode->ends_[0];
        targetNode->addEdge( this, overlap, drxn, false, true );
    }
    else if ( abs( targetNode->ends_[drxn] - ends_[!drxn] ) < params.readLen * 2 && !targetNode->clones_ )
    {
        int32_t qLen = seq_.length() - getBestOverlap( drxn );
        int32_t tLen = targetNode->seq_.length() - targetNode->getBestOverlap( !drxn );
        int overlap = min( min( qLen, tLen ), max( 30, abs( targetNode->ends_[drxn] - ends_[!drxn] ) ) );
        if ( overlap < 30 )
        {
            blankEnd( min( 30, qLen ), !drxn );
            targetNode->blankEnd( min( 30, tLen ), drxn );
            overlap = max( 0, 30 - min( qLen, tLen ) );
        }

        targetNode->addEdge( this, overlap, drxn, false, true );
    }
    
    return ( edges_[!drxn].empty() ? NULL : this );
}

NodeListList Node::foldEdgeGetTargets( Node* targetNode, int32_t* limits, bool drxn )
{
    NodeSet tSet = { targetNode };
    NodeSet tFwdSet;
    
    for ( Node* node : getDrxnNodes( drxn, true, true ) )
    {
        for ( auto &np : node->pairs_ )
        {
            if ( np.first->drxn_ <= 2 )
            {
                tFwdSet.insert( np.first );
                np.first->getDrxnNodes( tFwdSet, drxn );
            }
        }
    }
    
    targetNode->getDrxnNodesInSet( tSet, tFwdSet, !drxn );
    
    NodeSet usedSet;
    NodeListList targetPaths;
    NodeList targetPath;
    bool allowFork = drxn ? ends_[0] < targetNode->ends_[1] : targetNode->ends_[0] < ends_[1];
    targetNode->getFoldPath( targetPaths, targetPath, usedSet, tSet, limits, allowFork, !drxn );
    
    return targetPaths;
}

Node* Node::foldEnd( ExtVars &ev, Node* altNode, bool drxn )
{
    NodeSet altBckSet = altNode->getDrxnNodes( !drxn, false, true );
    NodeSet tSet = altNode->getDrxnNodes( drxn, false, true );
    NodeIntMap limitMap;
    
    int32_t target = drxn ? min( max( validLimits_[3], ends_[0] + params.readLen ), ends_[1] )
                          : max( min( validLimits_[0], ends_[1] - params.readLen ), ends_[0] );
    int32_t limits[2] = { target - params.readLen, target + params.readLen };
    
    // Plot best fold path
    NodeList foldPath = foldEndGetFoldPath( altBckSet, drxn );
    
    // Get all paired target nodes and best target fork
    Node* tFork = foldEndGetPairs( limitMap, foldPath, tSet, drxn );
    
    if ( foldPath.empty() || !tFork ) return NULL;
    
    // Plot best target path toward fold end
    NodeList altPath;
    NodeListList targetPaths = foldEndGetAltPath( foldPath.back(), tFork, altPath, limits, drxn );
    
    vector<MapStruct> qStructs = getMapStructQuery( target, !drxn );
    vector<MapStruct> tStructs = getMapStructTarget( targetPaths, target, limits, drxn );
    MapResult result;
    for ( MapStruct &q : qStructs )
    {
        q.isEnd = true;
        for ( MapStruct &t : tStructs )
        {
            t.isEnd = false;
            mapFold( result, q, t, params.readLen, 12, q.nodes[0]->drxn_ == t.nodes[0]->drxn_, !drxn );
        }
    }
    
    if ( result.l && result.r )
    {
        return ( drxn ? result.l : result.r )->foldEndHit( ev, ( drxn ? result.r : result.l ), result, drxn );
    }
    else if ( !altPath.empty() )
    {
        return foldEndMiss( ev, altPath, target, limitMap[tFork], drxn );
    }
    
    return NULL;
}

Node* Node::foldEndHit( ExtVars &ev, Node* hitNode, MapResult &result, bool drxn )
{
    if ( drxn_ <= 2 && !edges_[drxn].empty() )
    {
        for ( Node* nxt : getNextNodes( drxn ) )
        {
            nxt->dismantleNode( ev.del, drxn );
        }
    }
    
    assert( foldPrep( ev, result.coords[!drxn], result.len, true, drxn ) );
    if ( hitNode->foldPrep( ev, result.coords[drxn], result.len, false, !drxn ) )
    {
        if ( result.len < ends_[1] - ends_[0] )
        {
            addEdge( hitNode, result.len, drxn, false, result.len <= 0 );
        }
        else if ( result.len == ends_[1] - ends_[0] )
        {
            interEdge( ev, hitNode, drxn );
        }
    }
    else
    {
        for ( Edge &e : hitNode->edges_[drxn] )
        {
            int overlap = result.len + e.overlap - abs( hitNode->ends_[drxn] - result.coords[drxn] );
            e.node->addEdge( this, overlap, !drxn, false, e.isLeap || result.len <= 0 );
        }
    }
    
    return ( edges_[!drxn].empty() ? NULL : this );
}

Node* Node::foldEndMiss( ExtVars &ev, NodeList &altPath, int32_t qLimit, int32_t tLimit, bool drxn )
{
    int dummy = 0;
    assert( foldPrep( ev, qLimit, dummy, true, drxn ) );
    
    for ( int i( 0 ); i < altPath.size(); i++ )
    {
        if ( drxn ? qLimit < altPath[i]->ends_[0] : altPath[i]->ends_[1] < qLimit )
        {
            if ( i > 0 )
            {
                for ( Node* nxt : altPath[i]->getNextNodes( drxn ) )
                {
                    int overlap = drxn ? ends_[1] - nxt->ends_[0] : nxt->ends_[1] - ends_[0];
                    if ( overlap < 0 )
                    {
                        nxt->addEdge( this, overlap, !drxn, false, true );
                    }
                }
            }
            else
            {
                int overlap = drxn ? ends_[1] - altPath[i]->ends_[0] : altPath[i]->ends_[1] - ends_[0];
                altPath[i]->addEdge( this, overlap, !drxn, false, true );
            }
            
            return this;
        }
    }
    
    if ( drxn ? qLimit < tLimit : tLimit < qLimit )
    {
        assert( altPath.back()->foldPrep( ev, tLimit, dummy, false, !drxn ) );
        altPath.back()->addEdge( this, -abs( qLimit - tLimit ), !drxn, false, true );
        return this;
    }
    return NULL;
}

NodeListList Node::foldEndGetAltPath( Node* foldEnd, Node* tFork, NodeList &altPath, int32_t* limits, bool drxn )
{
    NodeList targetPath;
    NodeList currList;
    NodeSet altSet;
    NodeSet forkBckSet = tFork->getDrxnNodes( !drxn, false, true );
    NodeSet foldBckSet = foldEnd->getDrxnNodes( !drxn, false, true );
    
    for ( Node* prv : foldEnd->getNextNodes( !drxn ) )
    {
        if ( forkBckSet.find( prv ) != forkBckSet.end() )
        {
            for ( Node* nxt : prv->getNextNodes( drxn ) )
            {
                if ( forkBckSet.find( nxt ) != forkBckSet.end() )
                {
                    currList.push_back( nxt );
                    altSet.insert( nxt );
                    nxt->getDrxnNodesInSet( altSet, forkBckSet, drxn );
                }
            }
        }
    }
    
    NodeIntMap scores;
    for ( Node* alt : altSet )
    {
        scores[alt] = alt->getPairHitsTotal();
    }
    
    bool doAdd = false;
    bool dontAdd = false;
    
    while ( !currList.empty() )
    {
        int iBest = -1;
        int bestScore = -1;
        
        for ( int i ( 0 ); i < currList.size(); i++ )
        {
            if ( altSet.find( currList[i] ) != altSet.end() )
            {
                int currScore = scores[currList[i]];
                for ( Node* fwd : currList[i]->getDrxnNodesInSet( altSet, !drxn ) )
                {
                    currScore += scores[fwd];
                }
                if ( currScore > bestScore )
                {
                    iBest = i;
                    bestScore = currScore;
                }
            }
        }
        
        if ( iBest < 0 ) break;
        Node* curr = currList[iBest];
        doAdd = ( doAdd || !dontAdd ) && drxn ? limits[0] < curr->ends_[1] - getBestOverlap( 1 ) && curr->ends_[0] < limits[1]
                                              : limits[0] < curr->ends_[1] && curr->ends_[0] + getBestOverlap( 0 ) < limits[1];
        dontAdd = dontAdd || doAdd;
        currList.clear();
        
        altPath.push_back( curr );
        if ( doAdd )
        {
            targetPath.push_back( curr );
        }
        
        for ( Node* nxt : curr->getNextNodes( drxn ) )
        {
            currList.push_back( nxt );
        }
    }
    
    NodeListList targetPaths;
    if ( !targetPath.empty() )
    {
        reverse( targetPath.begin(), targetPath.end() );
        targetPaths.push_back( targetPath );
    }
    
    NodeSet usedSet = { altPath.begin(), altPath.end() };
    for ( Node* bck : forkBckSet )
    {
        if ( bck->edges_[!drxn].empty() 
                && foldBckSet.find( bck ) == foldBckSet.end()
                && usedSet.find( bck ) == usedSet.end()
                && bck->ends_[0] <= limits[1] && limits[0] <= bck->ends_[1] )
        {
            NodeList thisPath;
            bck->getFoldPath( targetPaths, thisPath, usedSet, forkBckSet, limits, false, drxn );
        }
    }
    
    bool isComplex = false;
    for ( Node* node : altPath )
    {
        for ( Node* nxt : node->getNextNodes( drxn ) )
        {
            if ( altSet.find( nxt ) != altSet.end() && usedSet.find( nxt ) == usedSet.end() )
            {
                NodeList thisPath;
                nxt->getFoldPath( targetPaths, thisPath, usedSet, altSet, limits, false, drxn );
                isComplex = true;
            }
        }
    }
    
    if ( isComplex )
    {
        altPath.clear();
    }
    
    return targetPaths;
}

NodeList Node::foldEndGetFoldPath( NodeSet &altBckSet, bool drxn )
{
    NodeList foldPath;
    NodeSet foldSet;
    NodeSet bckSet = getDrxnNodesNotInSet( altBckSet, !drxn, true );
    for ( Node* bck : bckSet )
    {
        for ( Node* prv : bck->getNextNodes( !drxn ) )
        {
            if ( altBckSet.find( prv ) != altBckSet.end() )
            {
                foldSet.insert( bck );
                bck->getDrxnNodesInSet( foldSet, bckSet, drxn );
            }
        }
    }
    
    NodeIntMap scores;
    for ( Node* fold : foldSet )
    {
        scores[fold] = fold->getPairHitsTotal();
    }
    
    int x = 0;
    
    Node* curr = this;
    while ( curr )
    {
        foldPath.push_back( curr );
        curr = NULL;
        int currHits = -1;
        for ( Node* prv : foldPath.back()->getNextNodes( !drxn ) )
        {
            if ( foldSet.find( prv ) != foldSet.end() )
            {
                int prvHits = scores[prv];
                for ( Node* bck : prv->getDrxnNodesInSet( foldSet, !drxn ) )
                {
                    prvHits += scores[bck];
                }
                if ( prvHits > currHits )
                {
                    currHits = prvHits;
                    curr = prv;
                }
            }
        }
    }
    
    return foldPath;
}

Node* Node::foldEndGetPairs( NodeIntMap &limitMap, NodeList &foldPath, NodeSet &tSet, bool drxn )
{
    NodeSet tHits, tFwdSet;
    NodeIntMap pairMap;
    
    // Find all pairs between the fold path and the target set
    for ( Node* fold : foldPath )
    {
        NodeList tNodes;
        for ( auto &np : fold->pairs_ )
        {
            if ( tSet.find( np.first ) != tSet.end() )
            {
                tNodes.push_back( np.first );
                tHits.insert( np.first );
                np.first->getDrxnNodes( tFwdSet, drxn );
                auto r = pairMap.insert( np );
                if ( !r.second )
                {
                    r.second += np.second;
                }
            }
        }
        
        if ( !tNodes.empty() )
        {
            for ( ReadMark &mark : fold->getMarksBase( !drxn ) )
            {
                for ( Node* t : tNodes )
                {
                    auto it = t->reads_.find( mark.id );
                    if ( it != t->reads_.end() )
                    {
                        auto r = limitMap.insert( make_pair( t, it->second[!drxn] ) );
                        if ( !r.second )
                        {
                            r.first->second = drxn ? min( r.first->second, it->second[0] ) 
                                                   : max( r.first->second, it->second[1] );
                        }
                        auto it2 = fold->reads_.find( params.getPairId( mark.id ) );
                        if ( it2 != fold->reads_.end() )
                        {
                            r = limitMap.insert( make_pair( fold, it2->second[drxn] ) );
                            if ( r.second )
                            {
                                r.first->second = drxn ? max( r.first->second, it2->second[1] )
                                                       : min( r.first->second, it2->second[0] );
                            }
                        }
                    }
                }
            }
        }
    }
    
    if ( limitMap.find( this ) == limitMap.end() )
    {
        limitMap[this] = drxn ? min( ends_[1], max( min( ends_[0] + params.readLen, ends_[1] ), validLimits_[3] ) ) 
                              : max( ends_[0], min( max( ends_[1] - params.readLen, ends_[0] ), validLimits_[0] ) );
    }
    
    Node* tBestFork = NULL;
    int tBestScore = -1;
    
    // Set best fork among hit target nodes
    for ( Node* t : tHits )
    {
        if ( tFwdSet.find( t ) == tFwdSet.end() )
        {
            int thisScore = pairMap[t];
            for ( Node* fwd : t->getDrxnNodes( drxn ) )
            {
                auto it = pairMap.find( fwd );
                thisScore += it != pairMap.end() ? it->second : 0;
            }
            
            if ( thisScore > tBestScore )
            {
                tBestFork = t;
                tBestScore = thisScore;
            }
        }
    }
    
    return tBestFork;
}

//void Node::foldNodesIn( NodeSet &folded, NodeSet &foldable, NodeSet &acceptors, bool drxn )
//{
//    Node* pNodes[2] = { this, NULL };
//    int32_t accOffset;
//    
//    // See if there is a single forward acceptor
//    for ( Edge &e : edges_[drxn] )
//    {
//        if ( acceptors.find( e.node ) != acceptors.end() )
//        {
//            if ( pNodes[1] )
//            {
//                pNodes[1] = NULL;
//                break;
//            }
//            pNodes[1] = e.node;
//            accOffset = drxn ? e.node->ends_[0] + e.overlap - ends_[1]
//                             : e.node->ends_[1] - e.overlap - ends_[0];
//        }
//    }
//    
//    NodeSet tSet;
//    NodeIntList foldOffsets;
//    NodeSet usedFes;
//    
//    // Catalog all nodes to fold and their offsets
//    for ( Edge &e : edges_[drxn] )
//    {
//        int readCount = e.node->reads_.size();
//        if ( foldable.find( e.node ) != foldable.end() && readCount <= 3 )
//        {
//            int32_t eOffset = drxn ? e.node->ends_[0] + e.overlap - ends_[1]
//                                   : e.node->ends_[1] - e.overlap - ends_[0];
//            foldOffsets.push_back( make_pair( e.node, eOffset ) );
//            for ( Edge &fe : e.node->edges_[drxn] )
//            {
//                if ( foldable.find( fe.node ) != foldable.end() 
//                        && usedFes.find( fe.node ) == usedFes.end()
//                        && readCount + fe.node->reads_.size() <= 3 )
//                {
//                    int32_t feOffset = drxn ? fe.node->ends_[0] + fe.overlap - e.node->ends_[1] + eOffset
//                                            : fe.node->ends_[1] - fe.overlap - e.node->ends_[0] + eOffset;
//                    foldOffsets.push_back( make_pair( e.node, feOffset ) );
//                    usedFes.insert( fe.node );
//                }
//            }
//        }
//    }
//    
//    // Transfer reads from foldable nodes
//    vector<SeqNum> readIds[2];
//    for ( pair<Node*, int32_t> &fold : foldOffsets )
//    {
//        for ( auto &read : fold.first->reads_ )
//        {
//            int32_t readCoords[2] = { read.second[0] - fold.second, read.second[1] - fold.second };
//            Node* pNode = this;
//            vector<SeqNum>* pIds = &readIds[0];
//            if ( pNodes[1] && ( drxn ? pNodes[1]->ends_[0] < readCoords[0] : readCoords[1] < pNodes[1]->ends_[1] ) )
//            {
//                pNode = pNodes[1];
//                pIds = &readIds[1];
//                readCoords[0] += accOffset;
//                readCoords[1] += accOffset;
//            }
//            
//            readCoords[0] = max( readCoords[0], pNode->ends_[0] );
//            readCoords[1] = min( readCoords[1], pNode->ends_[1] );
//
//            if ( readCoords[0] < readCoords[1] )
//            {
//                pNode->addRead( read.first, readCoords[0], readCoords[1], true, bool(drxn_ % 3) );
//                pIds->push_back( params.getPairId( read.first ) );
//            }
//        }
//    }
//    
//    // Transfer pair node hits
//    for ( Node* t : tSet )
//    {
//        for ( int i( 0 ); i < 1 + bool( pNodes[1] ); )
//        {
//            for ( SeqNum &readId : readIds[i] )
//            {
//                if ( t->reads_.find( readId ) != t->reads_.end() )
//                {
//                    auto r = pNodes[i]->pairs_.insert( make_pair( t, 1 ) );
//                    if ( !r.second )
//                    {
//                        r.first->second++;
//                    }
//
//                    r = t->pairs_.insert( make_pair( this, 1 ) );
//                    if ( !r.second )
//                    {
//                        r.first->second++;
//                    }
//                    
//                    SeqNum pairId = params.getPairId( readId );
//                    pNodes[i]->islandUpdateValid( pairId );
//                }
//            }
//        }
//    }
//}

void Node::getFoldPath( NodeListList &targetPaths, NodeList &thisPath, NodeSet &usedSet, NodeSet &allowSet, int32_t* limits, bool allowFork, bool drxn )
{
    usedSet.insert( this );
    
    if ( drxn ? limits[0] < ends_[1] - getBestOverlap( 1 ) && ends_[0] < limits[1]
              : limits[0] < ends_[1] && ends_[0] + getBestOverlap( 0 ) < limits[1] )
    {
        thisPath.push_back( this );
    }
    
    bool didGoNext = false;
    
    if ( drxn ? ends_[1] - 12 < limits[1] : limits[0] < ends_[0] + 12 )
    {
        for ( Node* nxt : getNextNodes( drxn ) )
        {
            if ( usedSet.find( nxt ) == usedSet.end() 
                    && allowSet.find( nxt ) != allowSet.end() )
            {
                if ( allowFork || edges_[drxn].size() <= 1 )
                {
                    if ( !didGoNext )
                    {
                        nxt->getFoldPath( targetPaths, thisPath, usedSet, allowSet, limits, allowFork, drxn );
                        didGoNext = true;
                    }
                    else
                    {
                        NodeList newPath;
                        nxt->getFoldPath( targetPaths, newPath, usedSet, allowSet, limits, allowFork, drxn );
                    }
                }
                else
                {
                    if ( !didGoNext )
                    {
                        thisPath.push_back( nxt );
                        reverse( thisPath.begin(), thisPath.end() );
                        targetPaths.push_back( thisPath );
                        didGoNext = true;
                    }
                    else
                    {
                        targetPaths.push_back( NodeList( 1, nxt ) );
                    }
                }
            }
        }
    }
    
    if ( !didGoNext && !thisPath.empty() )
    {
        reverse( thisPath.begin(), thisPath.end() );
        targetPaths.push_back( thisPath );
    }
}

bool Node::foldPrep( ExtVars &ev, int32_t coord, int &overlap, bool isEnd, bool drxn )
{
    if ( isEnd && drxn_ > 2 )
    {
        for ( Node* nxt : getDrxnNodes( drxn ) )
        {
            nxt->removeEdge( this, !drxn );
        }
        edges_[drxn].clear();
    }
    
    if ( coord != ends_[drxn] )
    {
        int32_t splitCoord = ends_[drxn];
        for ( auto &read : reads_ )
        {
            if ( drxn ? coord < read.second[1] && read.second[0] < splitCoord
                      : read.second[0] < coord && splitCoord < read.second[1] )
            {
                splitCoord = read.second[!drxn];
            }
        }
        
        if ( splitCoord != ends_[drxn] )
        {
            splitNode( splitCoord, ( drxn_ <= 2 ? ev.nodes : ev.island ), ( drxn_ <= 2 ? drxn : !drxn ), drxn );
            overlap -= abs( coord - ends_[drxn] );
        }
        else
        {
            return false;
        }
    }
    return true;
}

MapStruct Node::getMapStructEnd( int32_t target, int32_t* limits, bool drxn )
{
    int32_t coord = drxn ? min( limits[1], ends_[1] ) : max( limits[0], ends_[0] );
    int32_t minLen = max( 0, getBestOverlap( drxn ) - max( 0, ( drxn ? ends_[1] - coord : coord - ends_[0] ) ) );
    int32_t seqLen = drxn ? coord - max( limits[0], ends_[0] ) 
                          : min( limits[1], ends_[1] ) - coord;
    seqLen = min( max( minLen, seqLen ), (int)seq_.length() - abs( coord - ends_[drxn] ) );
    string seq = drxn ? seq_.substr( coord - ends_[0] - seqLen, seqLen )
                      : seq_.substr( coord - ends_[0], seqLen );
    
    MapStruct ms( seq, minLen, coord, target, drxn );
    ms.nodes.push_back( this );
    
    return ms;
}

vector<MapStruct> Node::getMapStructQuery( int32_t target, bool drxn )
{
    int32_t limits[2] = { target - params.readLen, target + params.readLen };
    offsetForward( !drxn, true, true );
    
    vector<MapStruct> qStructs;
    MapStruct qStruct = getMapStructEnd( target, limits, drxn );
    getMapStructQuery( qStructs, qStruct, limits, drxn );
    return qStructs;
}

void Node::getMapStructQuery( vector<MapStruct> &qStructs, MapStruct &qStruct, int32_t* limits, bool drxn )
{
    if ( !edges_[!drxn].empty() && ( drxn ? limits[0] < ends_[0] : ends_[1] < limits[1] ) )
    {
        for ( Edge &e : edges_[!drxn] )
        {
            MapStruct q = qStruct;
            int seqLen = e.node->seq_.length() - e.overlap;
            seqLen = min( params.readLen * 2 - q.len, seqLen );
            q.len += seqLen;
            int32_t seqBgn = drxn ? e.node->seq_.length() - seqLen - e.overlap
                                  : e.overlap;
            string seq = e.node->seq_.substr( seqBgn, seqLen );
            q.seq = drxn ? seq + q.seq : q.seq + seq;
            q.nodes.push_back( e.node );
            q.minLen = drxn ? q.minLen + seqLen : q.minLen;
            
            e.node->getMapStructQuery( qStructs, q, limits, drxn );
        }
    }
    else
    {
        qStructs.push_back( qStruct );
    }
}

vector<MapStruct> Node::getMapStructTarget( Node* targetNode, int32_t target, bool drxn )
{
    vector<MapStruct> tStructs;
    int32_t targetLimit = drxn ? max( target, targetNode->validLimits_[0] )
                               : min( target, targetNode->validLimits_[3] ) ;
    int32_t limits[2] = { targetLimit - 200, targetLimit + 200 };
    
    NodeListList targetPaths = foldEdgeGetTargets( targetNode, limits, !drxn );
    getMapStructTargetCheckPaths( targetPaths, limits, drxn );
    for ( NodeList &targetPath : targetPaths )
    {
        MapStruct ms = getMapStructTarget( targetPath, targetLimit, limits, drxn );
        tStructs.push_back( ms );
    }
    
    return tStructs;
}

vector<MapStruct> Node::getMapStructTarget( NodeListList &targetPaths, int32_t target, int32_t* limits, bool drxn )
{
    vector<MapStruct> tStructs;
    getMapStructTargetCheckPaths( targetPaths, limits, drxn );
    for ( NodeList &targetPath : targetPaths )
    {
        MapStruct ms = getMapStructTarget( targetPath, target, limits, drxn );
        tStructs.push_back( ms );
    }
    return tStructs;
}

MapStruct Node::getMapStructTarget( NodeList &targetPath, int32_t target, int32_t* limits, bool drxn )
{
    MapStruct ms = targetPath[0]->getMapStructEnd( target, limits, drxn );
    
    for ( int i ( 1 ); i < targetPath.size(); i++ )
    {
        if ( ms.len < params.readLen * 2 )
        {
            int overlap = targetPath[i-1]->getOverlap( targetPath[i], !drxn );
            int seqLen = targetPath[i]->seq_.length() - overlap;
            seqLen = min( params.readLen * 2 - ms.len, seqLen );
            ms.len += seqLen;
            int32_t seqBgn = drxn ? targetPath[i]->seq_.length() - seqLen - overlap
                                  : overlap;
            string seq = targetPath[i]->seq_.substr( seqBgn, seqLen );
            ms.seq = drxn ? seq + ms.seq : ms.seq + seq;
            ms.nodes.push_back( targetPath[i] );
            ms.minLen = drxn ? ms.minLen + seqLen : ms.minLen;
        }
    }
    
    return ms;
}

void Node::getMapStructTargetCheckPaths( NodeListList &targetPaths, int32_t* limits, bool drxn )
{
    for ( int i ( 0 ); i < targetPaths.size(); )
    {
        if ( !targetPaths.empty() && 
                ( drxn ? targetPaths[i][0]->ends_[1] <= limits[0] 
                       : limits[1] <= targetPaths[i][0]->ends_[0] ) )
        {
            targetPaths[i].clear();
        }
        int32_t len;
        for ( int j ( 0 ); j < targetPaths[i].size(); )
        {
            // Check if path is not yet within limits
            if ( drxn ? limits[1] <= targetPaths[i][j]->ends_[0]
                      : targetPaths[i][j]->ends_[1] <= limits[0] )
            {
                targetPaths[i].erase( targetPaths[i].begin(), targetPaths[i].begin() + j + 1 );
                continue;
            }
            
            // Check if path has gone beyond limits
            if ( drxn ? targetPaths[i][j]->ends_[0] <= limits[0]
                      : limits[1] <= targetPaths[i][j]->ends_[1] )
            {
                targetPaths[i].erase( targetPaths[i].begin() + j + 1, targetPaths[i].end() );
            }
            
            // Check path edge is a leap
            if ( j > 0 )
            {
                bool doSplit = false;
                int overlap;
                for ( Edge &e : targetPaths[i][j]->edges_[drxn] )
                {
                    if ( e.node == targetPaths[i][j-1] )
                    {
                        doSplit = doSplit || e.isLeap;
                        overlap = e.overlap;
                    }
                }
                
                len += targetPaths[i][j]->ends_[1] - targetPaths[i][j]->ends_[0] - overlap;
                doSplit = doSplit || len > params.readLen * 2;
                        
                if ( doSplit )
                {
                    targetPaths.push_back( NodeList( targetPaths[i].begin() + j, targetPaths[i].end() ) );
                    targetPaths[i].erase( targetPaths[i].begin() + j, targetPaths[i].end() );
                    break;
                }
            }
            else
            {
                len = drxn ? min( limits[1], targetPaths[i][j]->ends_[1] ) - targetPaths[i][j]->ends_[0]
                           : targetPaths[i][j]->ends_[1] - max( limits[0], targetPaths[i][j]->ends_[0] );
            }
            
            j++;
        }
        if ( targetPaths[i].empty() )
        {
            targetPaths.erase( targetPaths.begin() + i );
            continue;
        }
        i++;
    }
}

bool Node::getNextReadCoord( int32_t &coord, bool coordDrxn, bool readDrxn )
{
    int32_t nxtCoord = ends_[readDrxn];
    for ( auto &read : reads_ )
    {
        if ( readDrxn ? read.second[coordDrxn] < nxtCoord && coord <= read.second[coordDrxn]
                      : nxtCoord < read.second[coordDrxn] && read.second[coordDrxn] <= coord )
        {
            nxtCoord = read.second[coordDrxn];
        }
    }
    coord = nxtCoord;
    return nxtCoord != ends_[readDrxn];
}

void Node::interEdge( ExtVars &ev, Node* node, bool drxn )
{
    node->clearPairs();
    
    int32_t offset = node->ends_[!drxn] - ends_[!drxn];
    for ( auto &read : reads_ )
    {
        node->addRead( read.first, read.second[0] + offset, read.second[1] + offset, true );
    }
    
    NodeSet nodeFwdSet = node->getDrxnNodes( !drxn );
    
    for ( Edge &e : edges_[!drxn] )
    {
        if ( nodeFwdSet.find( e.node ) == nodeFwdSet.end() )
        {
            e.node->addEdge( node, e.overlap, drxn );
        }
    }
    
    assert( false );
    
    dismantleNode();
    ev.del.insert( this );
    
    if ( node->drxn_ <= 2 && node->validated_ )
    {
        node->setValid();
    }
}

//void Node::joinEnds( IslandVars &iv, Node** nodes, int overlap )
//{
//    vector<Overlap> overlaps = iv.ev.bwt.mapJoin( nodes[0]->seq_, nodes[1]->seq_.substr( overlap ), overlap );
//    
//    if ( !overlaps.empty() )
//    {
//        uint16_t coords[2]{0};
////        NodeSet fwdSet = nodes[!iv.drxn]->getDrxnNodes( iv.drxn );
////        NodeSet bckSet = nodes[!iv.drxn]->getDrxnNodes( !iv.drxn, true, true );
////        for ( int i = 0; i < overlaps.size(); i++ )
////        {
////            coords[0] = max( coords[0], overlaps[i].overLen );
////            coords[1] = max( coords[1], overlaps[i].extLen );
////            for ( const NodeList &nodeList : { iv.ev.nodes, iv.ev.island } )
////            {
////                for ( Node* node : nodeList )
////                {
////                    auto it = node->reads_.find( overlaps[i].readId );
////                    if ( it != node->reads_.end() )
////                    {
////                        int j = i;
////                        while ( ++j < overlaps.size() && node->reads_.find( overlaps[j].readId ) != node->reads_.end() );
////                        if ( !i )
////                        {
////                            assert( false );
////                        }
////                        else
////                        {
////                            assert( false );
////                        }
////                    }
////                }
////            }
////        }
//        joinEnds( iv, nodes, overlaps, overlap );
//        
//        assert( iv.ev.del.find( nodes[0] ) == iv.ev.del.end() && iv.ev.del.find( nodes[1] ) == iv.ev.del.end() );
//        string seq = nodes[0]->seq_.substr( nodes[0]->seq_.length() - coords[0] ) + nodes[1]->seq_.substr( 0, coords[1] );
//        int edges[2] = { coords[0], coords[1] + overlap };
//        coords[0] = nodes[0]->ends_[1] - coords[0];
//        coords[1] = nodes[0]->ends_[1] + coords[1];
//        Node* node = new Node( seq, overlaps, iv.drxn + 3 );
//        iv.ev.island.push_back( node );
//        iv.merged[!iv.drxn].insert( node );
//        nodes[0]->addEdge( node, edges[0], 1, iv.drxn );
//        nodes[1]->addEdge( node, edges[1], 0, !iv.drxn );
//        node->setCoverage();
//        node->resetMarks();
//    }
//    else
//    {
//        nodes[!iv.drxn]->addEdge( nodes[iv.drxn], overlap, iv.drxn );
//        iv.merged[!iv.drxn].insert( nodes[iv.drxn] );
//    }
//    
//}

//bool Node::joinEnds( IslandVars &iv, Node** nodes, vector<Overlap> &overlaps, int overlap )
//{
//    int32_t limits[2]{0};
//    for ( int i = 0; i < overlaps.size(); i++ )
//    {
//        limits[0] = max( limits[0], overlaps[i].overLen );
//        limits[1] = max( limits[1], overlaps[i].extLen );
//        
//        for ( const NodeList &nodeList : { iv.ev.nodes, iv.ev.island } )
//        {
//            for ( Node* node : nodeList )
//            {
//                auto it = node->reads_.find( overlaps[i].readId );
//                if ( it != node->reads_.end() )
//                {
//                    int j = i;
//                    while ( ++j < overlaps.size() && node->reads_.find( overlaps[j].readId ) != node->reads_.end() );
//                    auto it2 = node->reads_.find( overlaps[j-1].readId );
//                    if ( !i )
//                    {
//                        NodeSet nxtSet = nodes[0]->getNextNodes( 1 );
//                        if ( nxtSet.find( node ) != nxtSet.end() )
//                        {
//                            
//                        }
//                        assert( false );
//                    }
//                    else
//                    {
//                        assert( false );
//                    }
//                }
//            }
//        }
//    }
//}

void Node::mapFold( MapResult &result, MapStruct &q, MapStruct &t, int32_t distWiggle, int minLen, bool fixedDist, bool drxn )
{
    MapStruct* l = drxn ? &t : &q;
    MapStruct* r = drxn ? &q : &t;
    int iBest = -1, jBest = -1, bestLen = 0;
    bool isSet = result.l && result.r;
    
    mapSeqs( *l, *r, result.score, iBest, jBest, bestLen, minLen, distWiggle, fixedDist, isSet );
    
    if ( bestLen > 0 )
    {
        mapSetResult( result, *l, *r, iBest, jBest, bestLen );
    }
}

void Node::mapSeqs( MapStruct &l, MapStruct &r, int &bestScore, int &iBest, int &jBest, int &bestLen, int minLen, int32_t distWiggle, bool fixedDist, bool &isSet )
{
    int32_t offset = r.coord - r.len - l.coord;
    
    int m[l.len][r.len];
    
    for ( int i ( 0 ); i < l.len; i++ )
    {
        m[i][0] = l.seq[i] == r.seq[0];
    }
    
    for ( int j ( 1 ); j < r.len; j++ )
    {
        m[0][j] = l.seq[0] == r.seq[j];
    }
    
    for ( int i ( 1 ); i < l.len; i++ )
    {
        for ( int j ( 1 ); j < r.len; j++ )
        {
            m[i][j] = l.seq[i] == r.seq[j] ? m[i-1][j-1] + 1 : 0;
            
            // Better than min length; goes beyond min i; starts before min j
            if ( m[i][j] >= minLen && i + 1 >= l.minLen && j + 1 - m[i][j] <= r.minLen )
            {
                int len = m[i][j];
                int dist = abs( i - j - offset ) / ( 1 + ( !fixedDist ) * 2 );
                int lCut = l.isEnd ? max( l.cutLen - i - 1, 0 ) : 0;
                int rCut = r.isEnd ? max( ( j - len + 1 ) - r.cutLen, 0 ) : 0;
                int lGap = l.isEnd ? max( ( i + 1 - len - l.cutLen ) / 7, 0 ) : 0;
                int rGap = r.isEnd ? max( ( r.cutLen - j - 1 ) / 7, 0 ) : 0;
                int score = len - dist - lCut - rCut - lGap - rGap;
                int adjScore = len - ( ( lCut + rCut ) / 2 ) - max( abs( i - offset - j ) - distWiggle, 0 );
                if ( adjScore >= 0 && ( score > bestScore || !isSet ) )
                {
                    iBest = i;
                    jBest = j;
                    bestScore = score;
                    bestLen = len;
                    isSet = true;
                }
            }
        }
    }
}

void Node::mapSetResult( MapResult &result, MapStruct &l, MapStruct &r, int iBest, int jBest, int len )
{
    result.l = l.nodes[0];
    result.r = r.nodes[0];
    result.coords[0] = l.coord + iBest + 1;
    result.coords[1] = r.coord - r.len + ( jBest - len + 1 );
    result.len = len;
    
    for ( int i( 1 ); i < l.nodes.size(); i++ )
    {
        int32_t coords[2] = { result.coords[0] - result.len, result.coords[0] };
        int32_t offset = l.nodes[i-1]->ends_[1] - l.nodes[i-1]->getOverlap( l.nodes[i], 1 ) - l.nodes[i]->ends_[0];
        bool mustTrim = !l.nodes[i]->anyReadBeyondCoord( coords[1] + 1 + offset, 1, 0 );
        
        if ( coords[1] <= l.nodes[i-1]->ends_[1] || coords[0] < l.nodes[i]->ends_[0] || mustTrim )
        {
            break;
        }
        else
        {
            result.coords[0] += offset;
            result.l = l.nodes[i];
        }
    }
    
    if ( result.l->ends_[1] < result.coords[0] )
    {
        result.len -= result.coords[0] - result.l->ends_[1];
        result.coords[0] = result.l->ends_[1];
    }
    
    for ( int i ( 1 ); i < r.nodes.size(); i++ )
    {
        int32_t coords[2] = { result.coords[1], result.coords[1] + result.len };
        int32_t offset = r.nodes[i]->ends_[1] - r.nodes[i]->getOverlap( r.nodes[i-1], 1 ) - r.nodes[i-1]->ends_[0];
        bool mustTrim = !r.nodes[i]->anyReadBeyondCoord( coords[0] - 1 + offset, 0, 1 );
        
        if ( r.nodes[i-1]->ends_[0] <= coords[0] || r.nodes[i]->ends_[1] < coords[1] || mustTrim )
        {
            break;
        }
        else
        {
            result.coords[1] += offset;
            result.r = r.nodes[i];
        }
    }
    
    if ( result.coords[1] < result.r->ends_[0] )
    {
        result.len -= result.r->ends_[0] - result.coords[1];
        result.coords[1] = result.r->ends_[0];
    }
}

void Node::truncateNode( int32_t coord, bool drxn )
{
    while ( !edges_[drxn].empty() )
    {
        edges_[drxn][0].node->removeEdge( this, !drxn );
        edges_[drxn].erase( edges_[drxn].begin() );
    }
    
    int trimCount = 0;
    for ( auto &read : reads_ )
    {
        trimCount += ( drxn ? coord < read.second[1] : read.second[0] < coord );
    }
    
    int32_t seqLen = drxn ? coord - ends_[0] : ends_[1] - coord;
    assert( seqLen >= getBestOverlap( !drxn ) );
    seq_ = drxn ? seq_.substr( 0, seqLen ) : seq_.substr( seq_.length() - seqLen );
    ends_[drxn] = coord;
    int32_t maxEnds[2] = { ends_[!drxn], ends_[!drxn] };
    for ( auto it = reads_.begin(); it != reads_.end(); )
    {
        if ( it->second[!drxn] != ends_[!drxn] )
        {
            if ( it->second[1] < ends_[0] || ends_[1] < it->second[0] 
                    || ( trimCount > 2 && ( it->second[0] < ends_[0] || ends_[1] < it->second[1] ) ) )
            {
                it = reads_.erase( it );
                continue;
            }
        }
        
        if ( it->second[0] < ends_[0] )
        {
            it->second[0] = ends_[0];
            it->second.redundant = true;
        }
        
        if ( ends_[1] < it->second[1] )
        {
            it->second[1] = ends_[1];
            it->second.redundant = true;
        }
        
        maxEnds[0] = min( maxEnds[0], it->second[0] );
        maxEnds[1] = max( maxEnds[1], it->second[1] );
        
        it++;
    }
    
    if ( drxn )
    {
        ends_[1] = min( ends_[1], maxEnds[1] );
        validLimits_[2] = min( ends_[1], validLimits_[2] );
        validLimits_[3] = min( ends_[1], validLimits_[3] );
        seq_ = seq_.substr( 0, ends_[1] - ends_[0] );
    }
    else
    {
        ends_[0] = max( ends_[0], maxEnds[0] );
        validLimits_[0] = max( ends_[0], validLimits_[0] );
        validLimits_[1] = max( ends_[0], validLimits_[1] );
        seq_ = seq_.substr( seq_.length() - ( ends_[1] - ends_[0] ) );
    }
    
    this->getLength();
}
