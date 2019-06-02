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
#include <limits>

void Node::addClone( Node* node, bool isFirst )
{
    if ( !clones_ ) clones_ = new NodeList;
    if ( isFirst )
    {
        node->addClone( this, false );
        for ( Node* clone : *clones_ )
        {
            node->addClone( clone, false );
            clone->addClone( node, false );
        }
    }
    if ( find( clones_->begin(), clones_->end(), node ) == clones_->end() )
    {
        clones_->push_back( node );
    }
}

bool Node::anyCloneInSet( NodeSet &nodeSet )
{
    if ( nodeSet.find( this ) != nodeSet.end() )
    {
        return true;
    }
    for ( Node* clone : getCloneSet() )
    {
        if ( nodeSet.find( clone ) != nodeSet.end() )
        {
            return true;
        }
    }
    return false;
}

//void Node::cloneNode( NodeList &nodes, bool drxn )
//{
//    Node* clone = new Node( this );
//    nodes.push_back( clone );
//    for ( Edge &e : edges_[!drxn] )
//    {
//        e.node->addEdge( clone, e.overlap, drxn, true, e.isLeap );
//    }
//}

//void Node::extendLoopNode( LoopVars &lv, ExtVars &ev, bool drxn )
//{
//    ev.ante.clear();
//    getDrxnNodes( ev.ante, !drxn, true );
//    
//    while ( isContinue( drxn ) && extendCount_ > 0 && lv.rebranch.find( this ) == lv.rebranch.end() )
//    {
//        vector<Extension> exts = ev.bwt.mapExtensions( seq_, drxn );
//        bool doesBranch = exts.size() > 1;
//        for ( auto it = exts.begin(); it != exts.end(); )
//        {
//            it->fwdExts.clear();
//            MergeHit merge;
//            
//            bool doAdd = true;
//            for ( Node* clone : lv.clones )
//            {
//                if ( clone->checkExtensionMerge( *it, merge ) )
//                {
//                    it = exts.erase( it );
//                    lv.rebranch.insert( this );
//                    doAdd = false;
//                    break;
//                }
//            }
//            for ( Node* node : ev.island )
//            {
//                if ( node->checkExtensionMerge( *it, merge ) )
//                {
//                    it = exts.erase( it );
//                    lv.rebranch.insert( this );
//                    doAdd = false;
//                    break;
//                }
//            }
//            it += doAdd;
//        }
//        
//        extendCount_ = max( 0, extendCount_ - (int)exts.size() );
//        addExtensions( exts, ev, doesBranch, drxn );
//    }
//    
//    lv.extended.insert( this );
//    setCoverage();
//}

CloneScore Node::getCloneComparison( Node* clone, bool drxn )
{
    NodeSet tSet;
    getDrxnNodes( tSet, !drxn, params.getFurthestMpDist( ends_[!drxn], !drxn ) );
    clone->getDrxnNodes( tSet, !drxn, params.getFurthestMpDist( clone->ends_[!drxn], !drxn ) );
    NodeList tNodes( tSet.begin(), tSet.end() );
    
    NodeIntMap pairs;
    vector<ReadMark> marks = getMarksBase( drxn );
    unordered_set<SeqNum> hitIds;
    vector<int32_t> offsets = { clone->ends_[drxn] - ends_[drxn] };
    CloneScore score = setCloneScore( pairs, tNodes, offsets, marks, hitIds, drxn, false );
    return score;
}

NodeSetList Node::getCloneEdgeSets( bool drxn )
{
    NodeSet cloneSet = getCloneSet();
    NodeSetList cloneEdgeSet;
    for ( Node* node : cloneSet )
    {
        for ( Edge &edge : node->edges_[drxn] )
        {
            bool found = false;
            for ( NodeSet edgeSet : cloneEdgeSet )
            {
                if ( edgeSet.find( edge.node ) != edgeSet.end() )
                {
                    found = true;
                    break;
                }
            }
            if ( !found )
            {
                cloneEdgeSet.push_back( edge.node->getCloneSet() );
            }
        }
    }
    return cloneEdgeSet;
}

//vector<int32_t> Node::getCloneOffsets()
//{
//    vector<int32_t> offsets;
//    if ( clones_ )
//    {
//        for ( Node* clone : *clones_ )
//        {
//            offsets.push_back( clone->ends_[0] - this->ends_[0] );
//        }
//    }
//    return offsets;
//}
//
//vector<int32_t> Node::getCloneOffsets( NodeList &tNodes )
//{
//    vector<int32_t> offsets;
//    if ( clones_ )
//    {
//        for ( Node* clone : *clones_ )
//        {
//            if ( find( tNodes.begin(), tNodes.end(), clone ) == tNodes.end() )
//            {
//                offsets.push_back( clone->ends_[0] - this->ends_[0] );
//            }
//        }
//    }
//    return offsets;
//}

CloneScore Node::getCloneScore( bool drxn )
{
    PairingVars pv;
    CloneTargetVars ctv = getCloneTargetVars( pv.tNodes, drxn );
    vector<ReadMark> marks = getMarksBase( drxn );
    pv.marks = &marks;
    return setCloneScore( pv, ctv, drxn );
}

NodeSet Node::getCloneSet( bool inclSelf )
{
    NodeSet cloneSet;
    if ( inclSelf )
    {
        cloneSet.insert( this );
    }
    if ( clones_ )
    {
        cloneSet.insert( clones_->begin(), clones_->end() );
    }
    return cloneSet;
}

CloneTargetVars Node::getCloneTargetVars( bool drxn )
{
    CloneTargetVars ctv;
    
    ctv.tSetThis.insert( this );
    getDrxnNodes( ctv.tSetThis, !drxn, params.getFurthestMpDist( ends_[!drxn], !drxn ) );
    
    if ( clones_ )
    {
        for ( Node* node : *clones_ )
        {
            ctv.offsets.push_back( node->ends_[1] -ends_[1] );

            NodeSet tSetClone = { node };
            node->getDrxnNodes( tSetClone, !drxn, params.getFurthestMpDist( node->ends_[!drxn], !drxn ) );
            ctv.tSetAlt.insert( tSetClone.begin(), tSetClone.end() );
        }
    }
    
    return ctv;
}

CloneTargetVars Node::getCloneTargetVars( NodeList &tNodes, bool drxn )
{
    CloneTargetVars ctv = getCloneTargetVars( drxn );
    tNodes.insert( tNodes.end(), ctv.tSetThis.begin(), ctv.tSetThis.end() );
    for ( Node* t : ctv.tSetAlt )
    {
        if ( ctv.tSetThis.find( t ) == ctv.tSetThis.end() )
        {
            tNodes.push_back( t );
        }
    }
    return ctv;
}

//CloneScore Node::getLoopBranchScore( LoopVars &lv, bool drxn )
//{
//    CloneScore branchScore = lv.scores[this];
//    branchScore.cumulScore = std::numeric_limits<int>::min();
//    for ( Node* fwd : getDrxnNodes( drxn ) )
//    {
//        auto it = lv.scores.find( fwd );
//        if ( it != lv.scores.end() )
//        {
//            branchScore += it->second;
//            branchScore.mergeFurthest( it->second, drxn );
//            if ( fwd->edges_[drxn].empty() )
//            {
//                branchScore.cumulScore = max( branchScore.cumulScore, it->second.cumulScore );
//            }
//        }
//    }
//    branchScore.furthest[0] = abs( branchScore.furthest[0] - ends_[!drxn] );
//    branchScore.furthest[1] = abs( branchScore.furthest[1] - ends_[!drxn] );
//    return branchScore;
//}

//void Node::loopDelete( LoopVars &lv, ExtVars &ev, NodeSet &delSet, NodeSet &goodSet, bool drxn )
//{
//    for ( Node* del : delSet )
//    {
//        if ( ev.del.find( del ) == ev.del.end() && goodSet.find( del ) == goodSet.end() )
//        {
//            del->loopDelete( lv, ev, drxn );
//        }
//    }
//}

//void Node::loopReview( LoopVars &lv, ExtVars &ev, bool drxn )
//{
//    NodeSet branches, cloneBranches;
//    NodeIntList endScores;
//    
//    loopReviewCatalog( lv, ev, branches, cloneBranches, drxn );
//    loopReviewRescore( lv, ev, branches, cloneBranches, endScores, drxn );
//    loopReviewBranches( lv, ev, branches, cloneBranches, endScores, drxn );
//    
//    NodeSet fwdSet = getDrxnNodes( drxn, true, true );
//    NodeSet anteSet = getDrxnNodes( !drxn, true, true );
//    NodeSet goodSet, anteEnds;
//    
//    if ( ev.del.find( this ) == ev.del.end() )
//    {
//        for ( Node* fwd : fwdSet )
//        {
//            for ( Node* clone : fwd->getCloneSet() )
//            {
//                assert( anteSet.find( clone ) != anteSet.end() );
//            }
//            if ( fwd->edges_[drxn].empty() )
//            {
//                if ( !fwd->clones_ )
//                {
//                    goodSet.insert( fwd );
//                    fwd->getDrxnNodesInSet( goodSet, fwdSet, !drxn );
//                }
//                else
//                {
//                    anteEnds.insert( fwd );
//                }
//            }
//        }
//    }
//    
//    while ( !anteEnds.empty() )
//    {
//        NodeSet nxtSet;
//        for ( Node* ante : anteEnds )
//        {
//            if ( fwdSet.find( ante ) != fwdSet.end() && goodSet.find( ante ) == goodSet.end() )
//            {
//                ante->getNextNodes( nxtSet, !drxn );
//                ante->loopDelete( lv, ev, drxn );
//            }
//        }
//        anteEnds = nxtSet;
//    }
//    
//    // Re-branch
//    for ( Node* node : lv.rebranch )
//    {
//        if ( ev.del.find( node ) == ev.del.end() )
//        {
//            node->rebranchNode( ev, drxn );
//        }
//    }
//}
//
//void Node::loopReviewBranches( LoopVars &lv, ExtVars &ev, NodeSet branches, NodeSet &cloneBranches, NodeIntList &endScores, bool drxn )
//{
//    NodeSet fwdSet = getDrxnNodes( drxn, true, true );
//    NodeSet goodSet, badSet, thisSet;
//    branches.insert( cloneBranches.begin(), cloneBranches.end() );
//    for ( Node* branch : branches )
//    {
//        thisSet.insert( branch );
//        branch->getDrxnNodes( thisSet, drxn );
//    }
//    
//    sort( endScores.begin(), endScores.end(), [&fwdSet]( pair<Node*, int32_t> &a, pair<Node*, int32_t> &b ) {
//        return a.second == b.second 
//                ? ( fwdSet.find( a.first ) == fwdSet.end() && fwdSet.find( b.first ) != fwdSet.end() )
//                : a.second > b.second;
//    });
//    
//    for ( pair<Node*, int32_t> &endScore : endScores )
//    {
//        if ( goodSet.find( endScore.first ) == goodSet.end() && badSet.find( endScore.first ) == badSet.end()  )
//        {
//            endScore.first->loopReviewBranchesVerify( lv, branches, goodSet, badSet, thisSet, drxn );
//        }
//    }
//    
//    loopReviewBranchesScavenge( lv, branches, goodSet, badSet, fwdSet, drxn );
//    loopReviewBranchesCleanup( lv, ev, branches, goodSet, badSet, drxn );
//    
//    NodeSet branchFwd;
//    for ( Node* node : branches )
//    {
//        node->getDrxnNodes( branchFwd, drxn );
//    }
//    
//    NodeSet anteSet = getDrxnNodes( !drxn, true, true );
//    for ( Node* fwd : getDrxnNodes( drxn, true, true ) )
//    {
//        for ( Node* clone : fwd->getCloneSet() )
//        {
//            assert( anteSet.find( clone ) != anteSet.end() );
//        }
//    }
//}
//
//void Node::loopReviewBranchesCleanup( LoopVars &lv, ExtVars &ev, NodeSet &branches, NodeSet &goodSet, NodeSet &badSet, bool drxn )
//{
//    for ( Node* branch : branches )
//    {
//        if ( goodSet.find( branch ) == goodSet.end() )
//        {
//            branch->loopDelete( lv, ev, drxn );
//        }
//    }
//    
//    NodeSet anteSet = getDrxnNodes( !drxn, true );
//    
//    for ( Node* fwd : getDrxnNodes( drxn, true, true ) )
//    {
//        if ( badSet.find( fwd ) != badSet.end() )
//        {
//            fwd->loopDelete( lv, ev, drxn );
//            continue;
//        }
//        for ( Node* clone : fwd->getCloneSet() )
//        {
//            if ( anteSet.find( clone ) == anteSet.end() )
//            {
//                if ( goodSet.find( fwd ) != goodSet.end() )
//                {
//                    clone->loopDelete( lv, ev, drxn );
//                }
//            }
//        }
//    }
//}
//
//void Node::loopReviewBranchesScavenge( LoopVars &lv, NodeSet &branches, NodeSet &goodSet, NodeSet &badSet, NodeSet &fwdSet, bool drxn )
//{
//    NodeSet anteSet = getDrxnNodes( !drxn, true, true );
//    
//    for ( Node* branch : branches )
//    {
//        if ( branch->clones_ 
//                && fwdSet.find( branch ) != fwdSet.end() 
//                && goodSet.find( branch ) == goodSet.end() 
//                && badSet.find( branch ) == badSet.end() )
//        {
//            int thisBest = -1, cloneBest = -1;
//            
//            for ( Node* fwd : branch->getDrxnNodes( drxn, true, true ) )
//            {
//                if ( fwd->edges_[drxn].empty() && fwd->getCloneEdgeSets( drxn ).empty() )
//                {
//                    thisBest = max( thisBest, lv.scores[fwd].cumulScore );
//                }
//            }
//            
//            NodeSet cloneSet;
//            
//            for ( Node* clone : branch->getCloneSet() )
//            {
//                if ( anteSet.find( clone ) == anteSet.end() )
//                {
//                    cloneSet.insert( clone );
//                }
//            }
//            
//            for ( Node* clone : cloneSet )
//            {
//                for ( Node* fwd : clone->getDrxnNodes( drxn, true, true ) )
//                {
//                    if ( fwd->edges_[drxn].empty() && fwd->getCloneEdgeSets( drxn ).empty() )
//                    {
//                        cloneBest = max( thisBest, lv.scores[fwd].cumulScore );
//                    }
//                }
//            }
//            
//            if ( thisBest == cloneBest )
//            {
//                for ( Node* bck : branch->getNextNodes( !drxn ) )
//                {
//                    auto it = lv.scores.find( bck );
//                    thisBest += ( it != lv.scores.end() && it->second.cumulScore - it->second.cumulMisses < 0 );
//                }
//            }
//            
//            if ( thisBest > 0 && thisBest > cloneBest )
//            {
//                goodSet.insert( branch );
//                badSet.insert( cloneSet.begin(), cloneSet.end() );
//            }
//            else
//            {
//                goodSet.insert( cloneSet.begin(), cloneSet.end() );
//                badSet.insert( branch );
//            }
//        }
//    }
//    
//    // Propagate good set then bad set backwards the forwards
//    NodeSet* sets[2] = { &goodSet, &badSet };
//    for ( bool j : { !drxn, drxn } )
//    {
//        for ( int i : { 0, 1 } )
//        {
//            NodeSet currSet;
//
//            for ( Node* fwd : fwdSet )
//            {
//                if ( sets[i]->find( fwd ) != sets[i]->end() )
//                {
//                    currSet.insert( fwd );
//                }
//            }
//
//            while ( !currSet.empty() )
//            {
//                NodeSet nxtSet;
//
//                for ( Node* curr : currSet )
//                {
//                    for ( Node* nxt : curr->getNextNodes( j ) )
//                    {
//                        if ( sets[i]->find( nxt ) == sets[i]->end() 
//                                && sets[!i]->find( nxt ) == sets[!i]->end() 
//                                && fwdSet.find( nxt ) != fwdSet.end() )
//                        {
//                            sets[i]->insert( nxt );
//                            for ( Node* clone : nxt->getCloneSet() )
//                            {
//                                if ( anteSet.find( clone ) == anteSet.end() )
//                                {
//                                    sets[!i]->insert( clone );
//                                }
//                            }
//                            nxtSet.insert( nxt );
//                        }
//                    }
//                }
//
//                currSet = nxtSet;
//            }
//        }
//    }
//}
//
//void Node::loopReviewBranchesVerify( LoopVars &lv, NodeSet &branches, NodeSet &goodSet, NodeSet &badSet, NodeSet &thisSet, bool drxn )
//{
//    NodeSet thisBranches;
//    NodeSet currSet = { this };
//    
//    while ( !currSet.empty() )
//    {
//        NodeSet nxtSet;
//        for ( Node* curr : currSet )
//        {
//            if ( branches.find( curr ) != branches.end() )
//            {
//                thisBranches.insert( curr );
//            }
//
//            int prevScore = numeric_limits<int>::min();
//            NodeSet prevSet;
//
//            for ( Node* bck : curr->getNextNodes( !drxn ) )
//            {
//                auto it = lv.scores.find( bck );
//                if ( it != lv.scores.end() )
//                {
//                    if ( it->second.cumulScore > prevScore )
//                    {
//                        prevSet.clear();
//                        prevScore = it->second.cumulScore;
//                    }
//                    if ( it->second.cumulScore == prevScore )
//                    {
//                        prevSet.insert( bck );
//                    }
//                }
//            }
//
//            for ( Node* prev : prevSet )
//            {
//                if ( goodSet.find( prev ) != goodSet.end() )
//                {
//                    thisBranches.insert( curr );
//                }
//                else if ( thisSet.find( prev ) != thisSet.end() && badSet.find( prev ) == badSet.end() )
//                {
//                    nxtSet.insert( prev );
//                }
//            }
//        }
//        currSet = nxtSet;
//    }
//    
//    NodeSet branchFwdSet = thisBranches;
//    NodeSet endPathSet = { this };
//    for ( Node* branch : thisBranches )
//    {
//        branch->getDrxnNodes( branchFwdSet, drxn );
//    }
//    getDrxnNodesInSet( endPathSet, branchFwdSet, !drxn );
//    
//    for ( Node* node : endPathSet )
//    {
//        goodSet.insert( node );
//        for ( Node* clone : node->getCloneSet() )
//        {
//            badSet.insert( clone );
//        }
//    }
//}
//
//void Node::loopReviewCatalog( LoopVars &lv, ExtVars &ev, NodeSet &branches, NodeSet &cloneBranches, bool drxn )
//{
//    NodeSet fwdSet = getDrxnNodes( drxn, true, true );
//    NodeSet anteSet = getDrxnNodes( !drxn, true, true );
//    NodeSet divSets[2];
//    
//    // Identify clone branches; Identify and classify ends; Classify divergent or antecedent
//    for ( Node* fwd : fwdSet )
//    {
//        if ( lv.scores.find( fwd ) == lv.scores.end() )
//        {
//            fwd->loopDelete( lv, ev, drxn );
//        }
//        else
//        {
//            for ( Node* clone : fwd->getCloneSet() )
//            {
//                if ( anteSet.find( clone ) == anteSet.end() )
//                {
//                    divSets[0].insert( fwd );
//                    divSets[1].insert( clone );
//                }
//            }
//        }
//    }
//    
//    for ( int i : { 0, 1 } )
//    {
//        for ( Node* node : divSets[i] )
//        {
//            bool isBranch = false;
//            for ( Node* prv : node->getNextNodes( !drxn ) )
//            {
//                isBranch = isBranch || divSets[i].find( prv ) == divSets[i].end();
//            }
//            
//            if ( isBranch )
//            {
//                ( i ? cloneBranches : branches ).insert( node );
//            }
//        }
//    }
//}
//
//void Node::loopReviewRescore( LoopVars &lv, ExtVars &ev, NodeSet &branches, NodeSet &cloneBranches, NodeIntList &endScores, bool drxn )
//{
//    CloneScoreMap scores;
//    
//    for ( Node* branch : branches )
//    {
//        branch->loopReviewRescoreBranch( scores, drxn );
//    }
//    
//    for ( Node* branch : cloneBranches )
//    {
//        branch->loopReviewRescoreBranch( scores, drxn );
//    }
//    
//    for ( auto &score : scores )
//    {
//        lv.scores[score.first] = score.second;
//    }
//    
//    
//    loopReviewRescoreEnds( lv, ev, branches, endScores, drxn );
//    loopReviewRescoreEnds( lv, ev, cloneBranches, endScores, drxn );
//}
//
//void Node::loopReviewRescoreBranch( CloneScoreMap &scores, bool drxn )
//{
//    PairingVars pv;
//    getCloneTargetVars( pv.tNodes, drxn );
//    
//    for ( Node* fwd : getDrxnNodes( drxn, true, true ) )
//    {
//        CloneScore score = fwd->getCloneScore( pv, drxn );
//        score.setScores();
//        auto r = scores.insert( make_pair( fwd, score ) );
//        if ( !r.second )
//        {
//            if ( r.first->second.selfScore - r.first->second.altScore < score.selfScore - score.altScore )
//            {
//                r.first->second = score;
//            }
//        }
//    }
//}
//
//void Node::loopReviewRescoreEnds( LoopVars &lv, ExtVars &ev, NodeSet &branches, NodeIntList &endScores, bool drxn )
//{
//    NodeSet currSet, endSet;
//    NodeSet fwdSet = branches;
//    
//    for ( Node* branch : branches )
//    {
//        branch->getDrxnNodes( fwdSet, drxn );
//        branch->getNextNodes( currSet, drxn );
//        if ( branch->edges_[drxn].empty() )
//        {
//            endSet.insert( branch );
//        }
//    }
//    
//    while ( !currSet.empty() )
//    {
//        NodeSet nxtSet, currFwd;
//        for ( Node* curr : currSet )
//        {
//            curr->getDrxnNodes( currFwd, drxn );
//        }
//        
//        for ( Node* curr : currSet )
//        {
//            int prevScore = numeric_limits<int>::min();
//            
//            for ( Node* prv : curr->getNextNodes( !drxn ) )
//            {
//                if ( fwdSet.find( prv ) != fwdSet.end() )
//                {
//                    assert( lv.scores.find( prv ) != lv.scores.end() );
//                    prevScore = max( prevScore, lv.scores[prv].cumulScore );
//                }
//            }
//            
//            if ( currFwd.find( curr ) == currFwd.end() )
//            {
//                lv.scores[curr].cumulScore += prevScore;
//                curr->getNextNodes( nxtSet, drxn );
//                
//                if ( curr->edges_[drxn].empty() 
//                        && ( ev.cloneSet.find( curr ) != ev.cloneSet.end() 
//                                || curr->getCloneEdgeSets( drxn ).empty() ) )
//                {
//                    endSet.insert( curr );
//                }
//            }
//            else
//            {
//                nxtSet.insert( curr );
//            }
//        }
//        currSet = nxtSet;
//    }
//    
//    for ( Node* node : endSet )
//    {
//        endScores.push_back( make_pair( node, lv.scores[node].cumulScore ) );
//    }
//}

//void Node::loopWalk( LoopVars &lv, ExtVars &ev, PairingVars &pv, int &cutoff, bool isAnte, bool isBranch, bool drxn )
//{
//    lv.open.erase( this );
//    
//    if ( loopWalkIsValid( drxn ) && loopWalkIsValidClone( lv.scores, pv, cutoff, isAnte, isBranch, drxn ) )
//    {
//        // Set next
//        loopWalkEdges( lv, ev, isAnte, drxn );
//        
//        // Propagate divergent clones if this is antecedent and reset score if good branch
//        for ( Node* nxt : getNextNodes( drxn ) )
//        {
//            if ( isAnte && nxt->clones_ && lv.branches.find( nxt ) != lv.branches.end() 
//                    && nxt->loopWalkDivergent( lv, ev, pv, cutoff, drxn ) > 0 )
//            {
//                auto itScore = lv.scores.find( this );
//                itScore->second.cumulScore = 0;
//                itScore->second.cumulMisses = 0;
//            }
//        }
//        
//        // Propagate remaining open clones
//        for ( Node* nxt : getNextNodes( drxn ) )
//        {
//            if ( nxt->clones_ && lv.open.find( nxt ) != lv.open.end() && loopWalkDoWalk( lv, drxn ) )
//            {
//                nxt->loopWalk( lv, ev, pv, cutoff, isAnte, isBranch, drxn );
//                return;
//            }
//        }
//    }
//}

//int Node::loopWalkDivergent( LoopVars &lv, ExtVars &ev, PairingVars &pv, int cutoff, bool drxn )
//{
//    lv.open.erase( this );
//    NodeSet anteOpenSet = lv.open;
//    lv.open.clear();
//    lv.open.insert( this );
//    
//    while ( !lv.open.empty() )
//    {
//        Node* nxt = loopWalkGetNext( lv, drxn );
//        nxt->loopWalk( lv, ev, pv, cutoff, false, nxt == this, drxn );
//    }
//    
//    lv.open = anteOpenSet;
//    
//    if ( !lv.merged.empty() && loopWalkDivergentMerge( lv, ev, drxn ) )
//    {
//        lv.branches.erase( this );
//        return 0;
//    }
//    
//    int scoreLimits[2] = { numeric_limits<int>::max(), numeric_limits<int>::min() };
//    bool anyValid = false; 
//    
//    for ( Node* fwd : getDrxnNodes( drxn, true, true ) )
//    {
//        if ( fwd->edges_[drxn].empty() && fwd->getCloneEdgeSets( drxn ).empty() )
//        {
//            auto it = lv.scores.find( fwd );
//            if ( it != lv.scores.end() )
//            {
//                anyValid = true;
//                scoreLimits[0] = min( scoreLimits[0], it->second.cumulScore );
//                scoreLimits[1] = max( scoreLimits[1], it->second.cumulScore );
//            }
//        }
//    }
//    
//    if ( ( scoreLimits[1] <= 0 && scoreLimits[0] < 0 ) || !anyValid )
//    {
//        return 0;
//    }
//    
//    return scoreLimits[0] + scoreLimits[1];
//}

//bool Node::loopWalkDivergentMerge( LoopVars &lv, ExtVars &ev, bool drxn )
//{
//    NodeSet fwdSet = getDrxnNodes( drxn, true, true );
//    NodeSet mergeSet;
//    
//    for ( auto it = lv.merged.begin(); it != lv.merged.end(); )
//    {
//        if ( fwdSet.find( *it ) != fwdSet.end() )
//        {
//            mergeSet.insert( *it );
//            (*it)->getDrxnNodesInSet( mergeSet, fwdSet, !drxn );
//            lv.open.insert( *it );
//            it = lv.merged.erase( it );
//            continue;
//        }
//        it++;
//    }
//    
//    if ( mergeSet.empty() )
//    {
//        return false;
//    }
//    
//    loopWalkRescore( lv, mergeSet, true, drxn );
//    
//    for ( Node* fwd : fwdSet )
//    {
//        if ( mergeSet.find( fwd ) == mergeSet.end() )
//        {
//            for ( Node* bck : fwd->getNextNodes( !drxn ) )
//            {
//                if ( mergeSet.find( bck ) != mergeSet.end() )
//                {
//                    lv.branches.insert( fwd );
//                    break;
//                }
//            }
//        }
//    }
//    
//    return true;
//}

//void Node::loopWalkEdges( LoopVars &lv, ExtVars &ev, bool isAnte, bool drxn )
//{
//    if ( lv.merged.find( this ) == lv.merged.end() )
//    {
//        Node* furthest = NULL;
//        vector<Edge> edges;
//        NodeSet usedEdges;
//
//        loopWalkEdgesCatalog( lv, ev, furthest, usedEdges, edges, isAnte, drxn );
//        assert( furthest != NULL );
//        furthest->loopWalkEdgesRebranch( lv, ev, usedEdges, edges, isAnte, drxn );
//        loopWalkEdgesClone( lv, ev, edges, isAnte, drxn );
//    }
//}

//void Node::loopWalkEdgesCatalog( LoopVars &lv, ExtVars &ev, Node* &furthest, NodeSet &usedEdges, vector<Edge> &edges, bool isAnte, bool drxn )
//{
//    // Block pre-existing forward edge nodes if they exist for some reason
//    for ( Edge &e : edges_[drxn] )
//    {
//        if ( e.node->clones_ )
//        {
//            usedEdges.insert( e.node->clones_->begin(), e.node->clones_->end() );
//            lv.open.insert( e.node );
//        }
//    }
//    
//    // Catalog all unique forward edge nodes for each clone; Extend any that are continuing
//    for ( Node* clone : getCloneSet( false ) )
//    {
//        furthest = !furthest || ( drxn ? furthest->ends_[1] < clone->ends_[1] : clone->ends_[0] < furthest->ends_[0] ) ? clone : furthest;
//        
//        for ( Edge &e : clone->edges_[drxn] )
//        {
//            if ( usedEdges.find( e.node ) == usedEdges.end() )
//            {
//                edges.push_back( e );
//                usedEdges.insert( e.node );
//                
//                if ( e.node->clones_ )
//                {
//                    usedEdges.insert( e.node->clones_->begin(), e.node->clones_->end() );
//                }
//                else if ( e.node->isContinue( drxn ) && lv.extended.find( e.node ) == lv.extended.end() )
//                {
//                    e.node->extendCount_ = max( isAnte ? 3 : 0, extendCount_ - 1 );
//                    e.node->extendLoopNode( lv, ev, drxn );
//                }
//            }
//        }
//    }
//}

//void Node::loopWalkEdgesClone( LoopVars &lv, ExtVars &ev, vector<Edge> &edges, bool isAnte, bool drxn )
//{
//    ev.ante.clear();
//    ev.ante.insert( this );
//    getDrxnNodes( ev.ante, !drxn, true );
//    
//    for ( Edge &e : edges )
//    {
//        if ( ( !e.node->clones_ || ( !loopWalkEdgesMerge( lv, ev, e, drxn ) ) && lv.cloned.find( e.node ) == lv.cloned.end()  ) )
//        {
//            Node* clone = new Node( e.node, ev, drxn );
//            clone->extendCount_ = max( 0, lv.extended.find( e.node ) != lv.extended.end() ? e.node->extendCount_ : extendCount_ - 1 );
//            addEdge( clone, e.overlap, drxn );
//            ev.nodes.push_back( clone );
//            lv.clones.insert( clone );
//            lv.cloned.insert( clone );
//            lv.cloned.insert( clone->clones_->begin(), clone->clones_->end() );
//            lv.open.insert( clone );
//            if ( isAnte && ev.ante.find( e.node ) == ev.ante.end() )
//            {
//                lv.branches.insert( clone );
//            }
//            if ( !isAnte && ev.ante.find( e.node ) != ev.ante.end() )
//            {
//                lv.merged.insert( clone );
//            }
//        }
//    }
//}

//bool Node::loopWalkEdgesMerge( LoopVars &lv, ExtVars &ev, Edge &edge, bool drxn )
//{
//    if ( lv.cloned.find( edge.node ) != lv.cloned.end() )
//    {
//        int32_t bestOffset;
//        Node* bestNode = NULL;
//        for ( Node* clone : edge.node->getCloneSet( true ) )
//        {
//            if ( ev.ante.find( clone ) == ev.ante.end() && lv.clones.find( clone ) != lv.clones.end() )
//            {
//                int32_t offset = abs( ( drxn ? ends_[1] - clone->ends_[0] : clone->ends_[1] - ends_[0] ) - edge.overlap );
//                if ( !bestNode || offset < bestOffset )
//                {
//                    bestOffset = offset;
//                    bestNode = clone;
//                }
//            }
//        }
//        
//        if ( bestNode  )
//        {
//            if ( bestOffset <= edge.node->getBiggestOffset( !drxn ) )
//            {
//                addEdge( bestNode, edge.overlap, drxn );
//            }
//            else
//            {
//                lv.rebranch.insert( bestNode );
//            }
//            
//            return true;
//        }
//    }
//    
//    return false;
//}

//void Node::loopWalkEdgesRebranch( LoopVars &lv, ExtVars &ev, NodeSet &usedEdges, vector<Edge> &edges, bool isAnte, bool drxn )
//{
//    if ( edges.size() < edgeCount_[drxn] )
//    {
//        ev.ante.clear();
//        getDrxnNodes( ev.ante, !drxn, true );
//        vector<Extension> exts = ev.bwt.mapExtensions( seq_, drxn );
//        bool doesBranch = true;
//        for ( auto it = exts.begin(); it != exts.end(); )
//        {
//            it->fwdExts.clear();
//            MergeHit merge;
//            
//            bool doAdd = true;
//            for ( Node* node : usedEdges )
//            {
//                if ( node->checkExtensionMerge( *it, merge ) )
//                {
//                    it = exts.erase( it );
//                    doAdd = false;
//                    break;
//                }
//            }
//            
//            for ( Node* clone : lv.clones )
//            {
//                if ( doAdd && clone->checkExtensionMerge( *it, merge ) )
//                {
//                    it = exts.erase( it );
//                    lv.rebranch.insert( this );
//                    doAdd = false;
//                    break;
//                }
//            }
//            it += doAdd;
//        }
//        
//        int i = edges_[drxn].size();
//        addExtensions( exts, ev, doesBranch, drxn );
//        
//        while ( i < edges_[drxn].size() )
//        {
//            edges.push_back( edges_[drxn][i] );
//            if ( edges_[drxn][i].node->isContinue( drxn ) )
//            {
//                edges_[drxn][i].node->extendCount_ = max( isAnte ? 3 : 0, extendCount_ - 1 );
//                edges_[drxn][i].node->extendLoopNode( lv, ev, drxn );
//            }
//            i++;
//        }
//    }
//}

//bool Node::loopWalkIsValid( bool drxn )
//{
//    int32_t limits[2] = { std::numeric_limits<int32_t>::min(), std::numeric_limits<int32_t>::min() };
//    propagateValidation( limits, drxn );
//    NodeSet invSet = getInvalidNodes( !drxn );
//    ScoreMap scores = Node::getScoreMap( invSet, limits, drxn );
//    Score score;
//    for ( pair<Node*, Score> nodeScore : scores )
//    {
//        score += nodeScore.second;
//    }
//    
//    return score[1] * 3 + 10 >= score[0];
//}

//bool Node::loopWalkRescore( LoopVars &lv, NodeSet &fwdSet, bool isAnte, bool drxn )
//{
//    CloneScoreMap scores;
//    CloneScore thisScore = lv.scores[this];
//    thisScore.setCumul( this, 0, 0, isAnte );
//    scores[this] = thisScore;
//    NodeSet currSet;
//    for ( Node* nxt : getNextNodes( drxn ) )
//    {
//        if ( fwdSet.find( nxt ) != fwdSet.end() )
//        {
//            currSet.insert( nxt );
//        }
//    }
//    
//    while ( !currSet.empty() )
//    {
//        NodeSet nxtSet;
//        for ( Node* curr : currSet )
//        {
//            bool doScore = true;
//            int prevScore = numeric_limits<int>::min();
//            int prevMisses = numeric_limits<int>::max();
//            
//            for ( Node* bck : curr->getNextNodes( !drxn ) )
//            {
//                auto it = scores.find( bck );
//                doScore = doScore && ( it != scores.end() || fwdSet.find( bck ) == fwdSet.end() );
//                if ( it != scores.end() )
//                {
//                    prevScore = max( prevScore, it->second.cumulScore );
//                    prevMisses = prevScore == it->second.cumulScore ? min( prevMisses, it->second.cumulMisses ) : prevMisses;
//                }
//            }
//            
//            if ( doScore )
//            {
//                CloneScore currScore = lv.scores[curr];
//                currScore.setCumul( curr, prevScore, prevMisses, isAnte );
//                scores[curr] = currScore;
//                
//                for ( Node* nxt : curr->getNextNodes( drxn ) )
//                {
//                    if ( fwdSet.find( nxt ) != fwdSet.end() && currSet.find( nxt ) == currSet.end() )
//                    {
//                        nxtSet.insert( nxt );
//                    }
//                }
//            }
//            else
//            {
//                nxtSet.insert( curr );
//            }
//        }
//        currSet = nxtSet;
//    }
//    
//    for ( auto &score : scores )
//    {
//        lv.scores[score.first] = score.second;
//    }
//}

void Node::removeClone( Node* node )
{
    if ( clones_ )
    {
        clones_->erase( remove( clones_->begin(), clones_->end(), node ), clones_->end() );
        if ( clones_->empty() )
        {
            delete clones_;
            clones_ = NULL;
        }
    }
}

//void Node::setLoopBranchScore( LoopVars &lv, PairingVars &pv, CloneTargetVars &ctv, bool drxn )
//{
//    NodeSet branchSet = getDrxnNodes( drxn );
//    NodeSet currSet = { this };
//    while ( !currSet.empty() )
//    {
//        NodeSet nxtSet;
//        for ( Node* curr : currSet )
//        {
//            bool doScore = true;
//            
//            // Ensure that all previous nodes in this branch have been scored
//            for ( Node* bck : curr->getNextNodes( !drxn ) )
//            {
//                if ( lv.scores.find( bck ) == lv.scores.end() && branchSet.find( bck ) != branchSet.end() )
//                {
//                    nxtSet.insert( curr );
//                    doScore = false;
//                }
//            }
//            
//            if ( doScore )
//            {
//                curr->setLoopNodeScore( lv, pv, ctv, drxn, curr != this );
//                curr->getNextNodes( nxtSet, drxn );
//            }
//        }
//        currSet = nxtSet;
//    }
//}
//
//void Node::setLoopNodeScore( LoopVars &lv, PairingVars &pv, CloneTargetVars &ctv, bool drxn, bool addCumulative )
//{
//    vector<ReadMark> marks = getMarksBase( drxn );
//    pv.marks = &marks;
//    CloneScore score = setCloneScore( pv, ctv, drxn );
//    score.setScores();
//    
//    if ( addCumulative )
//    {
//        int cumulScore = 0;
//        bool firstEdge = true;
//        for ( Node* bck : getNextNodes( !drxn ) )
//        {
//            auto hit = lv.scores.find( bck );
//            if ( hit != lv.scores.end() )
//            {
//                score.mergeFurthest( hit->second, drxn );
//                if ( ( firstEdge || hit->second.cumulScore > cumulScore ) )
//                {
//                    cumulScore = hit->second.cumulScore;
//                }
//            }
//        }
//        score.cumulScore += cumulScore;
//    }
//    
//    lv.scores[this] = score;
//}

