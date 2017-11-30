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

#include "seed.h"
#include <algorithm>

Seed::Seed( string &header, string &seq, Querier &bwt, int errorCount )
: header_( header ), seq_( seq ), bwt_( bwt )
{
    validLimits_[0] = seq.length();
    validLimits_[1] = 0;
    tether_[0] = seq.length();
    tether_[1] = 0;
    ends_[0] = seq.length();
    ends_[1] = 0;
    assert( errorCount <= 20 );
    MappedSeqs ms = bwt_.mapSeed( seq, errorCount, true );
    
    for ( ReadStruct &read : ms.reads )
    {
        vector<int32_t> hits, coords;
        NodeSet fwdSet;
        tether_[0] = min( tether_[0], read.tether[0] );
        tether_[1] = max( tether_[1], read.tether[1] );
        
        for ( int i( 0 ); i < nodes_.size(); i++ )
        {
            int32_t coord;
            if ( nodes_[i]->seedCongruent( read, coord ) )
            {
                hits.push_back( i );
                coords.push_back( coord );
                nodes_[i]->getDrxnNodes( fwdSet, 0 );
            }
        }
        
        for ( int i ( 0 ); i < hits.size() && hits.size() > 1; )
        {
            if ( fwdSet.find( nodes_[hits[i]] ) != fwdSet.end() )
            {
                hits.erase( hits.begin() + i );
                coords.erase( coords.begin() + i );
            }
            else
            {
                i++;
            }
        }
        
        if ( hits.size() > 1 || hits.size() == 1 && ( nodes_[ hits[0] ]->ends_[1] != coords[0] || !nodes_[ hits[0] ]->edges_[1].empty() ) )
        {
            Node* node = new Node( read );
            nodes_.push_back( node );
            for ( int i ( 0 ); i < hits.size(); i++ )
            {
                if ( coords[i] != nodes_[ hits[i] ]->ends_[1] )
                {
                    nodes_[ hits[i] ]->seedSplit( nodes_, coords[i] );
                }
                node->addEdge( nodes_[ hits[i] ], 0 );
            }
        }
        else if ( hits.size() == 1 )
        {
            nodes_[ hits[0] ]->seedAdd( read );
        }
        else
        {
            Node* node = new Node( read );
            nodes_.push_back( node );
        }
    }
    
    for ( Node* node : nodes_ )
    {
        ends_[0] = min( ends_[0], node->ends_[0] );
        ends_[1] = max( ends_[1], node->ends_[1] );
        validLimits_[0] = min( validLimits_[0], node->validLimits_[1] );
        validLimits_[1] = max( validLimits_[1], node->validLimits_[2] );
    }
}

void Seed::assemble()
{
    for ( Node* node : nodes_ )
    {
        node->resetMarks();
        node->setCoverage();
    }
    
    NodeList extendNodes[2] = { nodes_, nodes_ }, dummy;
    int32_t limits[2] = { (int32_t)seq_.length() - params.maxPeMean - params.readLen, params.maxPeMean + params.readLen };
    
    while ( !extendNodes[0].empty() || !extendNodes[1].empty() )
    {
        extendNodes[0].clear();
        extendNodes[1].clear();
        NodeSet seedSet, delSet;
        for ( Node* node : nodes_ )
        {
            if ( node->isSeed( seq_.length() ) )
            {
                seedSet.insert( node );
            }
        }
        
        Node::seedGetExtend( extendNodes, seedSet, delSet, limits );
        
        for ( int drxn : { 0, 1 } )
        {
            for ( Node* node : extendNodes[drxn] )
            {
                if ( !node->clones_ )
                {
                    ExtVars ev( nodes_, dummy, validLimits_, bwt_, false, false );
                    ev.ante = node->getDrxnNodes( !drxn );
                    node->extendCount_ = 1 + ( ( params.maxPeMean * 3 ) / params.readLen );
                    node->extendNode( ev, drxn );
                    ends_[drxn] = drxn ? max( ends_[1], node->ends_[1] ) : min( ends_[0], node->ends_[0] );
                }
            }
        }
        
        Node::seedValidate( seedSet, delSet, validLimits_, ends_ );
        for ( Node* del : delSet )
        {
            nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
            delete del;
        }
    }
}

void Seed::checkDivergent( NodeList &path )
{
    NodeSet drxnSets[2] = {  path.back()->getDrxnNodes( 0, false, true ), path[0]->getDrxnNodes( 1, false, true ) };
    NodeSet pathSet = path[0]->getDrxnNodesInSet( drxnSets[0], 1, true );
    NodeList divNodes[2];
    
    for ( bool drxn : { 0, 1 } )
    {
        for ( Node* node : pathSet )
        {
            for ( Node* nxt : node->getNextNodes( drxn ) )
            {
                if ( pathSet.find( nxt ) == pathSet.end() )
                {
                    divNodes[drxn].push_back( nxt );
                }
            }
        }
        
        checkDivergentBack( divNodes[drxn], drxnSets[drxn], drxn );
    }
    
    NodeIntMap divScores;
    for ( bool drxn : { 0, 1 } )
    {
        for ( Node* div : divNodes[drxn] )
        {
            int score = div->getPairHitsTotal();
            for ( Node* fwd : div->getDrxnNodes( drxn ) )
            {
                score += fwd->getPairHitsTotal();
            }
            assert( score < 2 );
            divScores[div] = score;
        }
    }
}

void Seed::checkDivergentBack( NodeList &div, NodeSet &pathSet, bool drxn )
{
    NodeSet delSet;
    for ( int i ( 0 ); i < div.size(); )
    {
        NodeSet fwdSet = div[i]->getDrxnNodes( drxn, false, true );
        NodeList branches;
        for ( Node* fwd : fwdSet )
        {
            for ( Node* prv : fwd->getNextNodes( !drxn ) )
            {
                if ( find( branches.begin(), branches.end(), prv ) == branches.end()
                        && pathSet.find( fwd ) == pathSet.end() )
                {
                    branches.push_back( prv );
                }
            }
        }
        
        for ( Node* node : branches )
        {
            NodeSet tSet = node->getDrxnNodesInSet( pathSet, drxn );
            NodeSet qSet = node->getDrxnNodesNotInSet( pathSet, !drxn, true );
            int score[2] = { 0, 0 };
            for ( Node* t : tSet )
            {
                for ( auto &np : t->pairs_ )
                {
                    if ( qSet.find( np.first ) != qSet.end() )
                    {
                        score[0] += np.second;
                    }
                    else if ( tSet.find( np.first ) == tSet.end() && pathSet.find( np.first ) != pathSet.end() )
                    {
                        score[1] += np.second;
                    }
                }
            }
            
            if ( score[0] > score[1] && score[0] > 2 )
            {
                for ( Node* nxt : node->getNextNodes( drxn ) )
                {
                    for ( Node* prv : nxt->getNextNodes( !drxn ) )
                    {
                        if ( pathSet.find( nxt ) != pathSet.end() && pathSet.find( prv ) != pathSet.end() )
                        {
                            nxt->removeEdge( prv, !drxn );
                            prv->removeEdge( nxt, drxn );
                            if ( prv->edges_[drxn].empty() )
                            {
                                prv->dismantleNode( delSet, !drxn );
                            }
                        }
                    }
                }
            }
            else
            {
                for ( Node* nxt : node->getNextNodes( drxn ) )
                {
                    if ( pathSet.find( nxt ) != pathSet.end() )
                    {
                        node->removeEdge( nxt, drxn );
                        nxt->removeEdge( node, !drxn );
                    }
                }
                
                if ( node->edges_[drxn].empty() )
                {
                    node->dismantleNode( delSet, !drxn );
                }
            }
        }
        
        i++;
    }
    
    for ( Node* del : delSet )
    {
        delete del;
        div.erase( remove( div.begin(), div.end(), del ), div.end() );
        nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
    }
}

vector<Locus*> Seed::getLoci()
{
    resolveBackForks();
    vector<Locus*> loci;
    
    NodeIntMap scores;
    NodeList bestPath = getLociGetPath( loci.empty() );
//    checkDivergent( bestPath );
    
    Node* forks[2];
    while ( getLociSetDivergent( bestPath, forks ) )
    {
        getLociResolveDivergent( scores, bestPath, forks );
    }
    
    while ( getLociSetConvergent( bestPath, forks ) )
    {
        if ( abs( forks[0]->ends_[0] + forks[0]->getBestOverlap( 1 ) ) < abs( forks[1]->ends_[1] - forks[1]->getBestOverlap( 0 ) ) )
        {
            bestPath.erase( bestPath.begin(), find( bestPath.begin(), bestPath.end(), forks[0] ) );
        }
        else
        {
            bestPath.erase( find( bestPath.begin(), bestPath.end(), forks[1] ) + 1, bestPath.end() );
        }
    }
    
    assert( !getLociSetDivergent( bestPath, forks ) );
    
    if ( !bestPath.empty() )
    {
        forks[0] = forks[0] ? forks[0] : bestPath[0];
        forks[1] = forks[1] ? forks[1] : bestPath.back();
        NodeList nodes[3];
        NodeList forkList( find( bestPath.begin(), bestPath.end(), forks[0] ), find( bestPath.begin(), bestPath.end(), forks[1] ) + 1 );
        
        Node* origin = NULL;
        int32_t bestTether = 0;
        for ( Node* node : forkList )
        {
            int32_t nodeTether = min( tether_[0] - node->ends_[0], node->ends_[1] - tether_[1] );
            if ( !origin || nodeTether > bestTether )
            {
                origin = node;
                bestTether = nodeTether;
            }
        }
        
        NodeSet locusSet = { origin };
        nodes[2].push_back( origin );
        for ( int drxn : { 0, 1 } )
        {
            for ( Node* node : origin->getDrxnNodes( drxn ) )
            {
                node->drxn_ = drxn;
                nodes[drxn].push_back( node );
                locusSet.insert( node );
            }
            origin->offsetForward( drxn, false, true );
        }
        
        Locus* locus = new Locus( bwt_, nodes );
        locus->header_ = header_;
        loci.push_back( locus );
        
        for ( Node* node : locusSet )
        {
            nodes_.erase( remove( nodes_.begin(), nodes_.end(), node ), nodes_.end() );
        }
    }
    
    for ( Node* node : nodes_ )
    {
        node->dismantleNode();
        delete node;
    }
    
    return loci;
}

NodeList Seed::getLociGetPath( bool doForce )
{
    NodeList startNodes = getLociGetPathGetStarts( doForce );
    NodeList bestPath;
    int bestPathScore = -1;
    
    for ( Node* node : startNodes )
    {
        NodeList path = { node };
        int score = getLociGetPathCurrScore( node, path, 1 );
        getLociGetPath( node, path, score, 1 );
        
        if ( score > bestPathScore )
        {
            bestPathScore = score;
            bestPath = path;
        }
    }
    
    if ( !bestPath.empty() && !bestPath[0]->edges_[0].empty() )
    {
        reverse( bestPath.begin(), bestPath.end() );
        getLociGetPath( bestPath.back(), bestPath, bestPathScore, 0 );
        reverse( bestPath.begin(), bestPath.end() );
    }
    
    assert( !bestPath.empty() );
    
    return bestPath;
}

void Seed::getLociGetPath( Node* curr, NodeList &path, int &score, bool drxn )
{
    while ( !curr->edges_[drxn].empty() )
    {
        Node* bestNxt = NULL;
        int bestNxtScore = -1;
        for ( Node* nxt : curr->getNextNodes( drxn ) )
        {
            int nxtScore = 0;
            for ( Node* fwd : nxt->getDrxnNodes( drxn, false, true ) )
            {
                for ( Node* node : path )
                {
                    nxtScore += fwd->getPairHits( node );
                }
            }

            if ( nxtScore > bestNxtScore )
            {
                bestNxt = nxt;
                bestNxtScore = nxtScore;
            }
        }

        assert( bestNxt );
        curr = bestNxt;
        path.push_back( curr );
        score += getLociGetPathCurrScore( curr, path, drxn );
    }
    
}

int Seed::getLociGetPathCurrScore( Node* curr, NodeList &path, bool drxn )
{
    int score = 0;
    
    if ( drxn ? tether_[0] <= curr->ends_[0] || tether_[1] <= curr->ends_[1] 
              : curr->ends_[0] <= tether_[0] || curr->ends_[1] <= tether_[1] )
    {
        for ( Node* node : path )
        {
            if ( drxn ? node->ends_[0] <= tether_[0] || node->ends_[1] <= tether_[1]
                      : tether_[0] <= node->ends_[0] || tether_[1] <= curr->ends_[1] )
            {
                score += node->getPairHits( curr );
            }
        }
    }
    
    return score;
}

NodeList Seed::getLociGetPathGetStarts( bool doForce )
{
    NodeList spanNodes, startNodes, endNodes;
    NodeSet spanFwd;
    for ( Node* node : nodes_ )
    {
        if ( node->ends_[0] <= tether_[0] || node->ends_[1] <= tether_[1] )
        {
            bool isSpan = false;
            for ( auto &np : node->pairs_ )
            {
                isSpan = isSpan || tether_[0] <= np.first->ends_[0] || tether_[1] <= np.first->ends_[1];
            }
            if ( isSpan )
            {
                spanNodes.push_back( node );
                node->getDrxnNodes( spanFwd, 1 );
            }
        }
        
        if ( node->edges_[0].empty() )
        {
            endNodes.push_back( node );
        }
    }
    
    for ( Node* node : spanNodes )
    {
        if ( spanFwd.find( node ) == spanFwd.end() )
        {
            startNodes.push_back( node );
        }
    }
    
    return ( startNodes.empty() && doForce ? endNodes : startNodes );
}

void Seed::getLociResolveDivergent( NodeIntMap &scores, NodeList &path, Node** forks )
{
    NodeSet pathSet( path.begin(), path.end() );
    int divScores[2] = { 0, 0 };
    NodeSet divSets[2];
    
    for ( int drxn : { 0, 1 } )
    {
        divSets[drxn] = forks[drxn]->getDrxnNodesNotInSet( pathSet, drxn );
        for ( Node* fwd : divSets[drxn] )
        {
            divScores[drxn] += scores[fwd];
        }
    }
    
    NodeSet delSet;
    for ( int drxn : { 0, 1 } )
    {
        if ( divScores[drxn] < divScores[!drxn] || divScores[drxn] < 2 )
        {
            for ( Node* node : divSets[drxn] )
            {
                node->dismantleNode( delSet, drxn );
            }
        }
    }
    
    for ( Node* del : delSet )
    {
        delete del;
        nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
    }
}

bool Seed::getLociSetConvergent( NodeList &path, Node** forks )
{
    forks[0] = NULL;
    forks[1] = NULL;
    NodeSet pathSet = { path.begin(), path.end() };
    
    for ( Node* node : path )
    {
        if ( node->edges_[0].size() > 1 )
        {
            forks[0] = node;
        }
        if ( node->edges_[1].size() > 1 && !forks[1] )
        {
            forks[1] = node;
        }
    }
    
    if ( forks[0] && forks[1] )
    {
        NodeSet fwdSet = forks[1]->getDrxnNodesNotInSet( pathSet, 1 );
        NodeSet bckSet = forks[0]->getDrxnNodesInSet( fwdSet, 0 );
        return !bckSet.empty();
    }
    return false;
    
}

bool Seed::getLociSetDivergent( NodeList &path, Node** forks )
{
    forks[0] = NULL;
    forks[1] = NULL;
    bool problem = false;
    NodeSet pathSet = { path.begin(), path.end() };
    NodeSet drxnSets[2] = { path.back()->getDrxnNodes( 0 ), path[0]->getDrxnNodes( 1 ) };
    
    for ( Node* node : path )
    {
        for ( int drxn : { 0, 1 } )
        {
            for ( Node* nxt : node->getNextNodes( drxn ) )
            {
                if ( pathSet.find( nxt ) == pathSet.end()
                        && drxnSets[!drxn].find( nxt ) == drxnSets[!drxn].end() )
                {
                    problem = problem || ( !drxn && forks[1] );
                    forks[drxn] = !drxn || !forks[drxn] ? node : forks[drxn];
                }
            }
        }
    }
    
    return problem;
}

bool Seed::resolveBackFork( Node** forks, NodeSet &delSet )
{
    NodeSet qSets[2], tSets[2], forkSets[2];
    for ( bool drxn : { 0, 1 } )
    {
        forkSets[!drxn] = forks[!drxn]->getDrxnNodes( drxn, false, true );
        forks[!drxn]->getDrxnNodes( forkSets[!drxn], !drxn );
        qSets[drxn] = forks[drxn]->getDrxnNodesNotInSet( forkSets[!drxn], !drxn );
        tSets[drxn] = forks[drxn]->getDrxnNodes( drxn, false, true );
    }
    
    int midScore = 0;
    for ( Node* fwd : forks[0]->getDrxnNodesNotInSet( qSets[1], 1 ) )
    {
        bool isMid = tSets[1].find( fwd ) != tSets[1].end();
        for ( auto &np: fwd->pairs_ )
        {
            if ( tSets[0].find( np.first ) != tSets[0].end() || ( isMid && tSets[1].find( np.first ) != tSets[1].end() ) )
            {
                midScore += np.second;
            }
        }
    }
    
    if ( !qSets[0].empty() || !qSets[1].empty() )
    {
        int scores[2] = { 0, 0 };
        for ( bool drxn : { 0, 1 } )
        {
            for ( Node* q : qSets[drxn] )
            {
                for ( auto &np : q->pairs_ )
                {
                    if ( tSets[drxn].find( np.first ) != tSets[drxn].end() )
                    {
                        scores[drxn] += np.second;
                    }
                }
            }
        }
        
        bool delMid = midScore < scores[0] && midScore < scores[1];
        for ( bool drxn : { 0, 1 } )
        {
            bool delQ = !delMid && scores[drxn] <= scores[!drxn];
            for ( Node* nxt : forks[drxn]->getNextNodes( !drxn ) )
            {
                bool isQ = qSets[drxn].find( nxt ) != qSets[drxn].end();
                if ( ( isQ && delQ ) || ( !isQ && delMid ) )
                {
                    forks[drxn]->removeEdge( nxt, !drxn );
                    nxt->removeEdge( forks[drxn], drxn );
                    if ( nxt->edges_[!drxn].empty() )
                    {
                        nxt->dismantleNode( delSet, drxn );
                    }
                }
            }
        }
        
        return true;
    }
    return false;
}

bool Seed::resolveBackForkBypass( Node* fork, NodeSet &delSet )
{
    NodeSet fwdSet = fork->getDrxnNodes( 1 );
    NodeSet tSet = fork->getDrxnNodes( 0, false, true );
    bool didErase = true;
    
    while ( didErase && fork->edges_.size() > 1 )
    {
        didErase = false;
        for ( Node* nxt : fork->getNextNodes( 1 ) )
        {
            NodeSet qSet = nxt->getDrxnNodesInSet( fwdSet, 0 );
            if ( !qSet.empty() )
            {
                int fwdScore = 0, qScore = 0, qReads = 0;

                for ( Node* q : qSet )
                {
                    qReads += q->reads_.size();
                    for ( auto &np : q->pairs_ )
                    {
                        if ( tSet.find( np.first ) != tSet.end() || fwdSet.find( np.first ) != fwdSet.end() )
                        {
                            qScore += np.second;
                        }
                    }
                }

                for ( Node* fwd : fwdSet )
                {
                    for ( auto &np : fwd->pairs_ )
                    {
                        if ( tSet.find( np.first ) != tSet.end() )
                        {
                            fwdScore += np.second;
                        }
                    }
                }
                
                if ( qScore == 0 && fwdScore > 0 && qReads > params.peCover )
                {
                    assert( false );
                    for ( Node* q : fork->getNextNodes( 1 ) )
                    {
                        if ( qSet.find( q ) != qSet.end() )
                        {
                            fork->removeEdge( q, 1 );
                            q->removeEdge( fork, 0 );
                            if ( q->edges_[0].empty() )
                            {
                                q->dismantleNode( delSet, 1 );
                            }
                        }
                    }
                }
                else
                {
                    fork->removeEdge( nxt, 1 );
                    nxt->removeEdge( fork, 0 );
                }
                
                didErase = true;
                break;
            }
        }
    }
    
    return fork->edges_[1].size() <= 1;
}

bool Seed::resolveBackForkDouble( Node* fork, NodeSet &delSet )
{
    bool didResolve = true;
    while ( didResolve && fork->edges_[1].size() > 1 )
    {
        didResolve = false;
        Node* nodes[4] = { fork, NULL, NULL, NULL };
        for ( int i ( 0 ); i < fork->edges_[1].size() - 1 && !didResolve; i++ )
        {
            nodes[2] = fork->edges_[1][i].node;
            NodeSet prvSet = nodes[2]->getNextNodes( 0 );
            for ( int j ( i + 1 ); j < fork->edges_[1].size() && !didResolve; j++ )
            {
                nodes[3] = fork->edges_[1][j].node;
                for ( Node* prv : nodes[3]->getNextNodes( 0 ) )
                {
                    if ( prv != nodes[0] && prvSet.find( prv ) != prvSet.end() )
                    {
                        nodes[1] = prv;
                        resolveBackForkDouble( nodes, delSet );
                        didResolve = true;
                    }
                }
            }
        }
    }
    
    return fork->edges_[1].size() <= 1;
}

void Seed::resolveBackForkDouble( Node** forks, NodeSet &delSet )
{
    NodeSetList bckSets = Node::getNodeSetsExclusive( forks[0], forks[1], 0 );
    NodeSetList fwdSets = Node::getNodeSetsExclusive( forks[2], forks[3], 1 );
    int pairScores[2][2];
    int pref[2] = { -1, -1 };
    for ( int i : { 0, 1 } )
    {
        for ( int j : { 0, 1 } )
        {
            pairScores[i][j] = 0;
            for ( Node* bck : bckSets[i] )
            {
                for ( Node* fwd : fwdSets[j] )
                {
                    auto it = bck->pairs_.find( fwd );
                    if ( it != bck->pairs_.end() )
                    {
                        pairScores[i][j] += it->second;
                    }
                }
            }
        }
    }
    
    bool didResolve = false;
    for ( int i ( 0 ); i < 2 && !didResolve; i++ )
    {
        for ( int j ( 0 ); j < 2 && !didResolve; j++ )
        {
            if ( pairScores[i][j] > pairScores[i][!j] )
            {
                if ( pairScores[!i][j] <= pairScores[!i][!j] )
                {
                    forks[i]->removeEdge( forks[2+!j], 1 );
                    forks[2+!j]->removeEdge( forks[i], 0 );
                    forks[!i]->removeEdge( forks[2+j], 1 );
                    forks[2+j]->removeEdge( forks[!i], 0 );
                    didResolve = true;
                }
            }
        }
    }
    
    if ( !didResolve )
    {
        int minPairs, minReads, iMinPairs = -1, iMinReads = -1;
        bool anyPairs = false;
        
        for ( int i ( 0 ); i < 4; i++ )
        {
            int pairCount = 0, readCount = 0;
            
            for ( Node* node : ( i < 2 ? bckSets[i] : fwdSets[i-2] ) )
            {
                pairCount += node->getPairHitsTotal();
                readCount += node->reads_.size();
            }
            
            anyPairs = anyPairs || pairCount > 0;
            
            if ( iMinPairs = -1 || pairCount < minPairs )
            {
                iMinPairs = i;
                minPairs = pairCount;
            }
            
            if ( iMinReads = -1 || readCount < minReads )
            {
                iMinReads = i;
                minReads = pairCount;
            }
        }
        
        int iMinSet = ( anyPairs ? iMinPairs : iMinReads );
        for ( Node* node : ( iMinSet < 2 ? bckSets[iMinSet] : fwdSets[iMinSet-2] ) )
        {
            node->dismantleNode();
            delSet.insert( node );
            assert( false );
        }
    }
}

void Seed::resolveBackForks()
{
    NodeList forks;
    NodeSet delSet;
    for ( Node* node : nodes_ )
    {
        if ( node->edges_[1].size() > 1 
                && !resolveBackForkBypass( node, delSet )
                && !resolveBackForkDouble( node, delSet ) )
        {
            forks.push_back( node );
        }
    }
    
    for ( Node* del : delSet )
    {
        nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
        forks.erase( remove( forks.begin(), forks.end(), del ), forks.end() );
    }
    
    for ( Node* fork : forks )
    {
        bool didResolve = true;
        while ( didResolve && fork->edges_[1].size() > 1 && delSet.find( fork ) == delSet.end() )
        {
            didResolve = false;
            
            for ( Node* branch : fork->getNextNodes( 1 ) )
            {
                Node* forkPair[2] = { fork, branch };

                while ( forkPair[1] && !didResolve )
                {
                    Node* nxt = NULL;
                    if ( forkPair[1]->edges_[0].size() > 1 )
                    {
                        didResolve = resolveBackFork( forkPair, delSet );
                    }
                    else if ( forkPair[1]->edges_[1].size() == 1 )
                    {
                        nxt = forkPair[1]->edges_[1][0].node;
                    }
                    forkPair[1] = nxt;
                }
            }
        }
    }
}

Seed::~Seed()
{
//    for ( Node* node : nodes_ )
//    {
//        delete node;
//    }
}

