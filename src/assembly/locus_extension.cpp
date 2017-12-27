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
#include <algorithm>

void Locus::extendLocus()
{
    clock_t startTime = clock();
    
    for ( bool drxn : { 0, 1 } )
    {
        endNodes_[drxn] = originEnds_[drxn];
        setExtend( drxn );
    }
    
    while ( ( canExtend( 0 ) || canExtend( 1 ) ) )
    {
        for ( int drxn : { 0, 1 } )
        {
            while ( canExtend( drxn ) )
            {
                extendNodes( drxn );
                if ( !updateExtension( drxn ) ) break;
            }
        }
        locusTest();
        if ( !plot() ) leap();
        locusTest();
    }
    
    duration_ = (double)( clock() - startTime ) / (double) CLOCKS_PER_SEC;
}

Node* Locus::advanceEnd( Node* node, ScoreMap &scores, NodeList &sideNodes, NodeSet &delSet, bool drxn, bool isEnd )
{
    node->propagateValidation( validLimits_, drxn );
    
    
    Node* prevNode = NULL;
    while ( node != prevNode )
    {
        prevNode = node;
        NodeList tNodes = node->getTargetNodes( drxn, true );
        scores = getScoreMap( node, tNodes, drxn );
        
        // Compute a score for each next node
        ScoreList nxtScores, stopScores;
        advanceEndScoreNext( node, scores, nxtScores, stopScores, drxn );
        
        // Advance to best next node if possible
        if ( !nxtScores.empty() )
        {
            sort( nxtScores.begin(), nxtScores.end(), []( const pair<Node*, Score> &a, const pair<Node*, Score> &b )
            {
                return a.second.hits == b.second.hits ? a.second.misses < b.second.misses : a.second.hits > b.second.hits;
            });
            
            if ( nxtScores[0].first->isValidated() && nxtScores[0].second[1] > 2 )
            {
                node = advanceEndSetNext( nxtScores, stopScores, sideNodes, delSet, isEnd, drxn );
            }
        }
    }
    
//    if ( node->deleteTest( drxn ) )
//    {
//        node->propagateValidation( validLimits_, drxn );
//        assert( false );
//    }
    
    return node;
}

void Locus::advanceEndScoreNext( Node* node, ScoreMap &scores, ScoreList &nxtScores, ScoreList &stopScores, bool drxn )
{
    for ( Edge &edge : node->edges_[drxn] )
    {
        Score thisScore;
        thisScore[0] += float( 10 * ( params.readLen - edge.overlap ) ) / float( params.readLen );
        NodeSet nxtFwdSet = { edge.node };
        edge.node->getDrxnNodes( nxtFwdSet, drxn );
        bool doesContinue = false;
        for ( Node* fwd : nxtFwdSet )
        {
            thisScore += scores[fwd];
            doesContinue = doesContinue ? : fwd->isContinue( drxn );
        }
        if ( doesContinue )
        {
            nxtScores.push_back( make_pair( edge.node, thisScore ) );
        }
        else
        {
            stopScores.push_back( make_pair( edge.node, thisScore ) );
        }
    }
}

Node* Locus::advanceEndSetNext( ScoreList &nxtScores, ScoreList &stopScores, NodeList &sideNodes, NodeSet &delSet, bool isEnd, bool drxn )
{
    // Parse inferior branches
    float cutoff = float( 8 + ( isEnd * 12 ) ) / float( nxtScores.size() );
    for ( int i( 1 ); i < nxtScores.size(); i++ )
    {
        if ( nxtScores[i].second[1] < 2 && nxtScores[i].second[0] > cutoff + ( nxtScores[i].second[1] * 3 ) )
        {
            delSet.insert( nxtScores[i].first );
        }
        else if ( isEnd || nxtScores[i].second[1] - nxtScores[i].second[0] >= 2 )
        {
            sideNodes.push_back( nxtScores[i].first );
        }
    }

    // Delete bad stopped branches
    for ( pair<Node*, Score> &stopScore : stopScores )
    {
        if ( stopScore.second[1] <= 2 || stopScore.second[1] < nxtScores[0].second[1] / 3 )
        {
            delSet.insert( stopScore.first );
        }
    }

    // Set new advanced node
    return nxtScores[0].first;
}

bool Locus::canExtend( bool drxn )
{
    harsh_[drxn] = forkCount_[drxn] > ( abs( ends_[drxn] ) + 1000 ) / 10;
    bool doTerminate = forkCount_[drxn] > ( abs( ends_[drxn] ) + 1000 ) / 3;
    if ( doTerminate )
    {
        stopCodes_[drxn] = 6;
    }
    for ( Node* node : toExtend_[drxn] )
    {
        if ( doTerminate )
        {
            node->stop_[drxn] = 6;
        }
    }
    return !doTerminate && !completed_[drxn];
}

bool Locus::debriefExtend( ExtVars &ev, bool drxn, bool rePlot )
{
    rePlot = rePlot || !ev.bypass.empty() || !ev.offset.empty() || !ev.cloneSet.empty();
    
    while ( !ev.rebranch.empty() || !ev.bypass.empty() || !ev.offset.empty() || !ev.cloneSet.empty() )
    {
        // Resolve rebranches
        if ( !ev.rebranch.empty() )
        {
            Node::resolveRebranch( ev, drxn );
        }
        
        // Resolve bypass nodes
        while ( !ev.bypass.empty() )
        {
            Node::resolveBypass( ev, drxn );
        }
        
        // Resolve offset nodes
        if ( !ev.offset.empty() )
        {
            Node::resolveOffsets( ev, drxn );
            locusTest();
            continue;
        }
        
        // Resolve cloned nodes
        if ( !ev.cloneSet.empty() )
        {
            Node::resolveLoops( ev, drxn );
            deleteNodes( ev.del, drxn );
            ev.del.clear();
            continue;
        }
    }
    
    if ( rePlot && doReplot( drxn ) )
    {
        plotPrep( drxn );
        for ( Path &path : paths_[drxn] )
        {
            while ( !plotPath( path, drxn ) );
        }
        
        return false;
    }
    
    return true;
}

bool Locus::doExtendNode( Node* &node, bool drxn )
{
    bool doExtend = ( drxn ? node->ends_[1] < limits_[1] : node->ends_[0] > limits_[0] );
    return doExtend;
}

bool Locus::doReplot( bool drxn )
{
    if ( !updateForks( drxn ) )
    {
        bool doReplot = false;
        for ( NodeList* nodes : { &toExtend_[drxn], &endNodes_[drxn], &sideNodes_[drxn] } )
        {
            for ( Node* node : *nodes )
            {
                doReplot = doReplot || find( nodes_[drxn].begin(), nodes_[drxn].end(), node ) == nodes_[drxn].end();
            }
        }
        
        return doReplot;
    }
    return true;
}
        
bool Locus::extendNodes( bool drxn )
{
    for ( Node* node : toExtend_[drxn] )
    {
        if ( doExtendNode( node, drxn ) )
        {
            ExtVars ev( nodes_[drxn], nodes_[drxn + 3], validLimits_, bwt_ );
            node->extendNode( ev, drxn );
            forkCount_[drxn] += node->edges_[drxn].size();
            locusTest();
            if ( !debriefExtend( ev, drxn ) ) break;
        }
    }
}

bool Locus::getForwardLoops( Node* &loop, ExtVars &ev, bool drxn )
{
    if ( !ev.del.empty() )
    {
        deleteNodes( ev.del, drxn );
        ev.del.clear();
    }
    
    Node::propagateValidation( ev.cloneSet, validLimits_, drxn );
    reviewSet( ev.cloneSet, drxn );
    
    NodeList fwdLoops;
    for ( auto it = ev.cloneSet.begin(); it != ev.cloneSet.end(); )
    {
        if ( find( nodes_[drxn].begin(), nodes_[drxn].end(), *it ) == nodes_[drxn].end() || !(*it)->clones_ )
        {
            it = ev.cloneSet.erase( it );
            continue;
        }
        fwdLoops.push_back( *it );
        it++;
    }
    
    sort( fwdLoops.begin(), fwdLoops.end(), [&drxn]( Node* &a, Node* &b ){
        return ( drxn ? a->ends_[0] < b->ends_[0] : a->ends_[1] > b->ends_[1] );
    } );
    
    loop = fwdLoops.empty() ? NULL : fwdLoops[0];
    
    return !fwdLoops.empty();
}

//void Locus::resolveLoop( ExtVars &ev, bool drxn )
//{
//    Node* nxtLoop = NULL;
//    while ( getForwardLoops( nxtLoop, ev, drxn ) )
//    {
//        LoopVars lv( nxtLoop );
//        nxtLoop->loop( lv, ev, drxn );
//    }
//}

void Locus::reviewSet( NodeSet &qSet, bool drxn )
{
    NodeSet vSet = getValidSet(), endsSet, reviewSet, goodSet, badSet, delSet, testSet;
    NodeList reviewForks;
    
    // Get all nodes forward of a query node
    for ( Node* q : qSet )
    {
        testSet.insert( q );
        NodeSet anteSet = q->getDrxnNodes( !drxn, true );
        for ( Node* qClone : q->getCloneSet( true ) )
        {
            if ( anteSet.find( qClone ) == anteSet.end() )
            {
                testSet.insert( qClone );
            }
            
            for ( Node* fwdEnd : qClone->getEndNodes( drxn, true ) )
            {
                endsSet.insert( fwdEnd );
                fwdEnd->getDrxnNodesNotInSet( endsSet, vSet, !drxn );
            }
        }
    }
    
    // Get review forks; Forks are review nodes that are preceded by a valid node
    for ( Node* node : endsSet )
    {
        for ( Node* bck : node->getNextNodes( !drxn ) )
        {
            if ( bck->isValidated() )
            {
                reviewForks.push_back( node );
                reviewSet.insert( node );
                node->getDrxnNodes( reviewSet, drxn );
                break;
            }
        }
    }
    
    // Score all review nodes
    ScoreMap scores = Node::getScoreMap( reviewSet, validLimits_, drxn );
    
    // Review remaining unreviewed ends; This should not actually occur
    for ( Node* node : reviewForks )
    {
        node->reviewFork( scores, goodSet, badSet, drxn );
    }
    
    for ( Node* t : testSet )
    {
        if ( badSet.find( t ) != badSet.end() && goodSet.find( t ) == goodSet.end() )
        {
            delSet.insert( t );
        }
    }
    
    deleteNodes( delSet, drxn );
    
    for ( Node* del : delSet )
    {
        qSet.erase( del );
    }
}

bool Locus::setBestEnds( Node* bgn, ScoreMap &scores, NodeSet &extSet, NodeSet &delSet, bool drxn, bool isEnd )
{
    float endCutoff = -14 + ( isEnd ? 0 : 6 );
//    float endCutoff = -14;
    float forkCutoff = ( isEnd ? 5 : 2 ) * endCutoff;
    
    vector< pair<Node*, float> > endScores;
    
    for ( Node* fork : setBestEndsGetForks( bgn, scores, drxn ) )
    {
        float forkScore;
        bool didAdd = setBestEndsSetFork( fork, forkScore, scores, endScores, endCutoff, delSet, drxn ) || didAdd;
        
        NodeSet forkFwdSet = fork->getDrxnNodes( drxn );
        NodeListList endLists = getEndLists( forkFwdSet, drxn );
        float cutoff = endLists.empty() ? endCutoff : min( endCutoff, forkCutoff / float( endLists.size() ) );
        if ( endLists.size() > 80 )
        {
            fork->getDrxnNodes( delSet, drxn );
            continue;
        }
        
        for ( NodeList &endList : endLists )
        {
            float endScore = 10 * ( ( (float)getWeakestEdge( fork, endList[0], forkFwdSet, drxn ) / (float)params.readLen ) - (float)0.8 );
            didAdd = setBestEndsSetEnd( endList, forkScore, endScore, scores, endScores, cutoff, delSet, drxn ) || didAdd;
        }
    }
    
    return setBestEndsSetExtend( endScores, extSet, isEnd );
}

NodeSet Locus::setBestEndsGetForks( Node* bgn, ScoreMap &scores, bool drxn )
{
    NodeSet forkSet, currSet = { bgn };
    while ( !currSet.empty() )
    {
        NodeSet nxtSet;
        for ( Node* curr : currSet )
        {
            auto it = scores.find( curr );
            if ( curr->isContinue( drxn ) || !curr->isValidated() || it->second[1] < it->second[0] )
            {
                forkSet.insert( curr );
            }
            else
            {
                curr->getNextNodes( nxtSet, drxn );
            }
        }
        currSet = nxtSet;
    }
    return forkSet;
}

bool Locus::setBestEndsSetEnd( NodeList &endList, float &forkScore, float &endScore, ScoreMap &scores, NodeFloatList &endScores, float &cutoff, NodeSet &delSet, bool drxn )
{
    for ( Node* node : endList )
    {
        auto it = scores.find( node );
        endScore += it->second[1] - it->second[0];
    }
    
    if ( endScore >= cutoff )
    {
        endScores.push_back( make_pair( endList[0], endScore + forkScore ) );
        return true;
    }
    
    delSet.insert( endList.begin(), endList.end() );
    
    return false;
}

bool Locus::setBestEndsSetExtend( NodeFloatList &endScores, NodeSet &extSet, bool isEnd )
{
    sort( endScores.begin(), endScores.end(), []( const pair<Node*, int> &a, const pair<Node*, int> &b ){
        return a.second < b.second;
    });
    
    if ( !endScores.empty() )
    {
        int extendPool = isEnd ? 120 / endScores.size() : 30 / endScores.size();
        for ( pair<Node*, float> &score : endScores )
        {
            score.second -= endScores.back().second;
            int thisExtend = extendPool * ( score.second + 1 ) / ( endScores[0].second + 1 );
            score.first->extendCount_ = max( isEnd ? 5 : 3, min( thisExtend, isEnd ? 60 : 15 ) );
            extSet.insert( score.first );
        }
        return true;
    }
    return false;
}

bool Locus::setBestEndsSetFork( Node* fork, float &forkScore, ScoreMap &scores, NodeFloatList &endScores, float &endCutoff, NodeSet &delSet, bool drxn )
{
    forkScore = scores[fork][1] * 2 - scores[fork][0] 
        + ( 10 * ( ( float( fork->getBestOverlap( !drxn ) ) / float(params.readLen) ) - float(0.8) ) );
    
    if ( fork->edges_[drxn].empty() )
    {
        if ( forkScore < endCutoff && fork->drxn_ != 2 )
        {
            if ( scores[fork][1] >= 2 && !fork->clones_ )
            {
                fork->trimEnd( drxn );
            }
            else
            {
                delSet.insert( fork );
            }
        }
        else if ( fork->isContinue( drxn ) )
        {
            endScores.push_back( make_pair( fork, forkScore ) );
            return true;
        }
    }
    return false;
}

void Locus::setExtend( bool drxn )
{
    NodeSet extSet, goodSet, delSet;
    toExtend_[drxn].clear();
    
    setExtendEnds( extSet, delSet, drxn );
    setExtendSides( extSet, delSet, drxn );
    
    for ( Node* node : extSet )
    {
        goodSet.insert( node );
        node->getDrxnNodes( goodSet, !drxn );
        toExtend_[drxn].push_back( node );
        ends_[drxn] = drxn ? max( ends_[1], node->ends_[1] ) : min( ends_[0], node->ends_[0] );
    }
    
    sort( toExtend_[drxn].begin(), toExtend_[drxn].end(), []( Node* &a, Node* &b )
    { 
        return a->extendCount_ > b->extendCount_; 
    } );
    
    deleteNodes( delSet, goodSet, drxn );
//    summarise( drxn );
}

void Locus::setExtendEnds( NodeSet &extSet, NodeSet &delSet, bool drxn )
{
    for ( auto it = endNodes_[drxn].begin(); it != endNodes_[drxn].end(); )
    {
        ScoreMap scores;
        *it = advanceEnd( *it, scores, sideNodes_[drxn], delSet, drxn, true );
        bool doErase = false;
        for ( auto it2 = endNodes_[drxn].begin(); it2 != it; it2++ )
        {
            doErase = doErase ? : *it == *it2;
        }
        doErase = doErase ? : !setBestEnds( *it, scores, extSet, delSet, drxn, true );
        if ( doErase )
        {
            it = endNodes_[drxn].erase( it );
        }
        else
        {
            it++;
        }
    }
}

void Locus::setExtendLimits()
{
    limits_[0] = ends_[0] - max( params.readLen * 5, params.libs[0].size );
    limits_[1] = ends_[1] + max( params.readLen * 5, params.libs[0].size );
}

void Locus::setExtendSides( NodeSet &extSet, NodeSet &delSet, bool drxn )
{
    for ( int i( 0 ); i < sideNodes_[drxn].size(); )
    {
        bool doErase = extSet.find( sideNodes_[drxn][i] ) != extSet.end() || sideNodes_[drxn][i]->stop_[drxn] != 0;
        
        // Erase side node if it merges to an end node's end
        for ( Node* fwd : sideNodes_[drxn][i]->getDrxnNodes( drxn ) )
        {
            doErase = doErase ? : extSet.find( fwd ) != extSet.end();
            doErase = doErase ? : find( endNodes_[drxn].begin(), endNodes_[drxn].end(), fwd ) != endNodes_[drxn].end();
        }
        
        // Process side node
        ScoreMap scores;
        if ( !doErase )
        {
            Node* newEnd = advanceEnd( sideNodes_[drxn][i], scores, sideNodes_[drxn], delSet, drxn, false );
            sideNodes_[drxn][i] = newEnd;
        }
        
        for ( int j( 0 ); j != i; j++ )
        {
            doErase = doErase ? : sideNodes_[drxn][i] == sideNodes_[drxn][j];
        }
        
        doErase = doErase ? : !setBestEnds( sideNodes_[drxn][i], scores, extSet, delSet, drxn, false );
        
        if ( doErase )
        {
            sideNodes_[drxn].erase( sideNodes_[drxn].begin() + i );
        }
        else
        {
            i++;
        }
    }
}

bool Locus::updateExtension( bool drxn )
{
    if ( endNodes_[drxn].empty() )
    {
        return false;
    }
    setExtend( drxn );
//    summarise( drxn );
    int overLimit = 0;
    for ( Node* node : toExtend_[drxn] )
    {
        overLimit += !doExtendNode( node, drxn );
    }
    return overLimit >= 2 && toExtend_[drxn].size() != overLimit;
}
