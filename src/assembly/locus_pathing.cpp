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
#include <math.h>
#include <algorithm>
#include <complex>
#include <cassert>
#include <iostream>

void Locus::finalise()
{
    finalise_ = true;
    forkNodes_[0] = originEnds_[0];
    forkNodes_[1] = originEnds_[1];
    
    for ( bool drxn : { 0, 1 } )
    {
        for ( Path &path : paths_[drxn] )
        {
            path.spans.clear();
            path.alleles.clear();
            while( !plotPath( path, drxn ) );
        }
    }
}

bool Locus::plot()
{
    locusTest();
    for ( bool drxn : { 0, 1 } )
    {
        if ( !completed_[drxn] )
        {
            plotPrep( drxn );
            completed_[drxn] = true;
            for ( Path &path : paths_[drxn] )
            {
                if ( !path.completed )
                {
                    while( !plotPath( path, drxn ) );
                }
            }
        }
    }
    
    for ( bool drxn : { 0, 1 } )
    {
        if ( !completed_[drxn] )
        {
            setExtend( drxn );
        }
    }
    
    setExtendLimits();
    return !completed_[0] || !completed_[1];
}

bool Locus::plotPath( Path &path, bool drxn )
{
    endNodes_[drxn].clear();
    sideNodes_[drxn].clear();
    toExtend_[drxn].clear();
    path.reset( forkNodes_[drxn], drxn );
    setForkLimits();
    vector<Span> initSpans = path.spans;
    locusTest();
    
    // Prepare branches for first loop
    BranchList branches;
    NodeList notReliableNodes;
    PathBranch best, last;
    PathVars pv( bwt_, nodes_, remappedReads_[drxn], finalise_, drxn );
    plotPathGetFirst( pv, best, path, branches, drxn );
    
    while ( best.branch )
    {
        // Parse converging paths
        if ( !plotPathConverge( pv, best, path, drxn ) ) break;
//        if ( !plotPathConverge( pv, best.branch, path, drxn ) ) break;
        
        // Update current multiplicity state and update spans if changing state
        plotPathUpdateMultiplicity( best, path, notReliableNodes, drxn );
        
        // Update current reliability state and parse dubious but not multiplicitous path nodes
        plotPathUpdateReliability( best, path, notReliableNodes, drxn );

        // Check if branch pairs beyond spans
        plotPathUpdateSpans( best, path, drxn );
        
        // Add best to path and catalog branches
        if ( !plotPathAdd( pv, best, branches, path ) ) break;

        // Prepare next iteration
        plotPathGetNext( pv, best, last, path, branches, drxn );
    }
    
    if ( !pv.rerun )
    {
        // Prune weak divergent branches and set side nodes
        NodeSet goodSet, delSet;
    //    plotPathTrimPrep( pv, path, drxn );
//        plotPathTrimDivergent( pv, path, goodSet, delSet, drxn );
//        plotPathTrimEnds( pv, best, branches, goodSet, delSet, drxn );
        deleteNodes( delSet, goodSet, drxn );
        delSet.clear();
        locusTest();
        PathReview review( pv, path.path, reliable_, forkLimits_, finished_[!drxn] );
        pv.rerun = !review.review( path, sideNodes_[drxn], delSet );
        deleteNodes( delSet, drxn );
        locusTest();
    }
    
    // Set furthest fork
    plotPathUpdateFork( pv, path, drxn );
    locusTest();
    
    if ( pv.rerun )
    {
        path.spans = initSpans;
        return false;
    }
     
    // Set new end and delta nodes or determine if locus is complete in this direction
    plotPathSetEnds( pv, best, last, path, drxn );
    
    assert( path.path.back()->drxn_ == 2 || find( nodes_[drxn].begin(), nodes_[drxn].end(), path.path.back() ) != nodes_[drxn].end() );
    
    // Update spans for next round of pathing
    reviewSpans( drxn );
    
    assert( path.path.back()->drxn_ == 2 || find( nodes_[drxn].begin(), nodes_[drxn].end(), path.path.back() ) != nodes_[drxn].end() );
    
    return true;
}

bool Locus::plotPathAdd( PathVars &pv, PathBranch &best, BranchList &branches, Path &path )
{
    if ( best.branch->inheritEdges( pv.drxn ) )
    {
        pv.rerun = true;
        return false;
    }
    
    bool isReliable = best.branch->isReliable() && !pv.weak && !pv.unspanned;
    if ( !best.valid && !pv.unspanned )
    {
        pv.unspanned = best.branch;
        isReliable = false;
    }
    
    float farMod = min( float( best.farCoords[0] ) / float( params.maxPeMean - params.readLen ), 
                        ( best.reliable == 1 ? (float)1.5 : best.reliable ) );
    farMod = farMod < 1 ? ( farMod - 0.5 ) / 0.5 : min( float( 3 ), farMod );
    float reliFarMod = min( float( best.farCoords[1] ) / float( params.maxPeMean - params.readLen ), (float)3 );
    if ( best.farNodes[0] && best.farNodes[0]->drxn_ != best.branch->drxn_ )
    {
        farMod = max( (float)1, farMod );
        reliFarMod = max( (float)1, reliFarMod );
    }
    float score = min( params.branchMinHits, best.hits ) * farMod 
                + min( params.branchMinHits, best.reliable ) * reliFarMod;
    
    bool doScrutinise = abs( best.branch->ends_[!pv.drxn] - reliable_[!pv.drxn] ) > params.maxPeMean 
            && best.branch->drxn_ != 2;
    if ( score < params.branchMinHits )
    {
        float bckCover = best.branch->getCoverageDrxn( params.maxPeMean, !pv.drxn, false );
        float fwdCover = best.branch->getCoverageDrxn( params.maxPeMean, pv.drxn, true );
        float ratio = max( 0.3, min( 1.0, ( min( fwdCover, bckCover ) / params.cover ) / 0.9 ) );
        score /= ratio;
        isReliable = false;
    }
    
    if ( score < params.branchMinHits && doScrutinise && !pv.weak )
    {
        pv.weak = best.branch;
        isReliable = false;
    }
    
    // Catalogue all divergent branches for later processing
    for ( PathBranch &branch : branches )
    {
        if ( branch.state == 2 ) path.divergent.push_back( branch );
    }

    path.path.push_back( best.branch );
    
    if ( isReliable && !pv.unspanned && !pv.weak
            && pv.drxn ? reliable_[1] < best.branch->ends_[1] 
                       : best.branch->ends_[0] < reliable_[0] )
    {
        reliable_[pv.drxn] = best.branch->ends_[pv.drxn];
    }
    
    return true;
}

bool Locus::plotPathConverge( PathVars &pv, PathBranch &convFork, Path &path, bool drxn )
{
    if ( convFork.branch->edges_[!drxn].size() >= 2 && convFork.branch->drxn_ != 2 )
    {
        if ( convFork.branch->inheritEdges( !drxn ) )
        {
            pv.rerun = true;
            return false;
        }
        
        NodeSet currSet = { convFork.branch };
        NodeSet fwdSet = convFork.fwdSet;
        while ( !currSet.empty() )
        {
            NodeSet nxtSet;
            for ( Node* curr : currSet )
            {
                for ( Node* nxt : curr->getNextNodes( !drxn ) )
                {
                    if ( pv.tSet.find( nxt ) == pv.tSet.end() )
                    {
                        if ( path.isReliable ) nxt->isReliable( true );
                        nxt->getDrxnNodes( fwdSet, drxn );
                        pv.addTarget( nxt, fwdSet );
                        nxtSet.insert( nxt );
                    }
                }
            }
            currSet = nxtSet;
        }
        
        convFork.setScores( pv );
        convFork.farCoords[0] = max( 1, abs( convFork.farCoords[0] - convFork.branch->ends_[!drxn] ) );
        convFork.farCoords[1] = max( 1, abs( convFork.farCoords[1] - convFork.branch->ends_[!drxn] ) );
    }
    
    return true;
}

//bool Locus::plotPathConverge( PathVars &pv, Node* convEnd, Path &path, bool drxn )
//{
//    if ( convEnd->edges_[!drxn].size() >= 2 && convEnd->drxn_ != 2 )
//    {
//        NodeSet bckSet = convEnd->getDrxnNodesNotInSet( pv.tSet, !drxn );
//        BranchList newDivergent;
//        for ( PathBranch &div : path.divergent )
//        {
//            if ( div.state == 1 && bckSet.find( div.branch ) != bckSet.end() )
//            {
//                div.state = 0;
//                NodeSet convSet = div.branch->getDrxnNodesInSet( bckSet, drxn );
//                if ( !plotPathConverge( pv, div, convEnd, newDivergent, convSet, path, drxn ) )
//                {
//                    pv.rerun = true;
//                    return false;
//                }
//            }
//        }
//        
//        if ( !newDivergent.empty() )
//        {
//            path.divergent.insert( path.divergent.end(), newDivergent.begin(), newDivergent.end() );
//        }
//    }
//    return true;
//}
//
//bool Locus::plotPathConverge( PathVars &pv, PathBranch &convBranch, Node* convEnd, BranchList &newDiv, NodeSet &convSet, Path &path, bool drxn )
//{
//    // Compile prime and alt path
//    PathScore scores[2];
//    scores[1].bgnFork = convBranch.fork;
//    scores[1].endFork = convEnd;
//    auto it = convBranch.fork ? find( path.path.begin(), path.path.end(), convBranch.fork ) + 1 : path.path.begin();
//    scores[0].path.insert( scores[0].path.end(), ( it == path.path.end() ? path.path.begin() : it ), path.path.end() );
//    scores[1].path.push_back( convBranch.branch );
//    
//    if ( !plotPathConvergePath( pv, convBranch, convEnd, newDiv, path, scores[1].path, convSet, drxn ) )
//    {
//        return false;
//    }
//    
//    // Score prime and alt path
//    int32_t bgnCoord = scores[0].path[0]->ends_[!drxn];
//    int32_t bestOls[2] = { max( scores[0].path[0]->getBestOverlap( !drxn ), scores[1].path[0]->getBestOverlap( !drxn ) )
//                         , max( scores[0].path.back()->getOverlap( convEnd, drxn ), scores[1].path.back()->getOverlap( convEnd, drxn ) ) };
//    plotPathConvergeScore( pv, bgnCoord, bestOls, convEnd, scores[0], drxn );
//    plotPathConvergeScore( pv, bgnCoord, bestOls, convEnd, scores[1], drxn );
//    
//    if ( plotPathConvergeRedundant( scores[1], path, drxn ) )
//    {
//        return true;
//    }
//    
//    if ( plotPathConvergeConflict( scores[0], scores[1], path, drxn ) )
//    {
//        return false;
//    }
//    
//    if ( plotPathConvergeClones( scores, drxn ) )
//    {
//        return false;
//    }
//    
//    if ( plotPathConvergeCompare( scores, path, drxn ) )
//    {
//        return false;
//    }
//    
//    // Accept convergent path
//    path.alleles.push_back( Allele( scores ) );
//    scores[1].altPath = scores[0].path;
//    path.convergents.push_back( scores[1] );
//    pv.resetPathReliable( scores[0].path );
//    
//    if ( !convBranch.state )
//    {
//        pv.addTarget( scores[1].path );
//    }
//    
//    return true;
//}
//
//bool Locus::plotPathConvergeClones( PathScore* scores, bool drxn )
//{
//    bool best[2];
//    
//    // Determine if each path is not bettered by clone path
//    for ( int i : { 0, 1 } )
//    {
//        Node* cloneBgn = NULL;
//        NodeSet cloneSet;
//        for ( auto it = scores[i].path.begin(); it != scores[i].path.end(); it++ )
//        {
//            if ( (*it)->clones_ )
//            {
//                cloneBgn = *it;
//                cloneSet.insert( it, scores[i].path.end() );
//                break;
//            }
//        }
//        best[i] = !cloneBgn || plotPathConvergeClonesAnyBetter( cloneBgn, cloneSet, drxn );
//    }
//    
//    // Delete a path if it is bettered by a clone path while the other path is not
//    for ( int i : { 0, 1 } )
//    {
//        if ( !best[i] && best[!i] )
//        {
//            NodeSet delSet = { scores[i].path[0] };
//            deleteNodes( delSet, drxn );
//            return true;
//        }
//    }
//    
//    return false;
//}
//
//bool Locus::plotPathConvergeClonesAnyBetter( Node* origBgn, NodeSet &origSet, bool drxn )
//{
//    for ( Node* cloneBgn : origBgn->getCloneSet())
//    {
//        CloneScore score = origBgn->getCloneComparison( cloneBgn, drxn );
//        
//        NodeSet thisOrigSet = { origBgn };
//        Node* currClone = cloneBgn;
//        
//        while ( currClone )
//        {
//            Node* nxtLoopClone = NULL;
//            for ( Node* currCloneNxt : currClone->getNextNodes( drxn ) )
//            {
//                for ( Node* nxtClone : currCloneNxt->getCloneSet() )
//                {
//                    if ( origSet.find( nxtClone ) != origSet.end() && thisOrigSet.find( nxtClone ) == thisOrigSet.end() )
//                    {
//                        thisOrigSet.insert( nxtClone );
//                        score += nxtClone->getCloneComparison( currCloneNxt, drxn );
//                        nxtLoopClone = currCloneNxt;
//                    }
//                }
//            }
//            currClone = nxtLoopClone;
//        }
//        
//        if ( score.altPref > score.selfPref )
//        {
//            return false;
//        }
//    }
//    return true;
//}
//
//bool Locus::plotPathConvergeCompare( PathScore* scores, Path &path, bool drxn )
//{
//    int32_t furthest = min( params.maxPeMean - params.readLen, max( scores[0].furthest, scores[1].furthest ) ) * 0.8;
//    
//    bool badFar[2] = { scores[0].furthest < furthest, scores[1].furthest < furthest };
//    bool bestFar[2] = { scores[0].furthest > scores[1].furthest, scores[0].furthest < scores[1].furthest };
//    bool onlyReliable[2] = { false, false };
//    bool moreReliable[2] = { false, false };
//    plotPathConvergeCompareReliable( scores, path, onlyReliable, moreReliable, drxn );
//    
//    // Is one well supported while the other isn't? -- Delete other
//    for ( int i : { 0, 1 } )
//    {
//        bool thisBetter = onlyReliable[i] && !onlyReliable[!i];
//        thisBetter = thisBetter || ( bestFar[i] 
//                && ( !scores[!i].furthest 
//                || ( scores[i].reliable && !scores[!i].reliable ) ) );
//        thisBetter = thisBetter || ( badFar[!i]
//                && ( scores[i].reliable > scores[!i].reliable 
//                || ( scores[i].furthest > scores[!i].furthest * 2 ) ) );
//        
//        if ( thisBetter )
//        {
//            NodeSet thisDelSet = { scores[!i].path[0] };
//            deleteNodes( thisDelSet, drxn );
//            return true;
//        }
//    }
//            
//    float combLen = float( scores[0].len + scores[1].len ) / float(2);
//    float coCover =  float( params.readLen * ( scores[0].readCount + scores[1].readCount ) ) / combLen;
//    bool notDiploid = coCover > params.cover * 1.4;
//    
//    // Is at least one not haploid?
//    if ( notDiploid )
//    {
//        for ( int i : { 0, 1 } )
//        {
//            if ( !badFar[i] 
//                    && ( ( moreReliable[i] && !moreReliable[!i] ) 
//                    || ( scores[i].reliable >= 3 && scores[!i].reliable < 1 ) ) )
//            {
//                NodeSet thisDelSet = { scores[!i].path[0] };
//                deleteNodes( thisDelSet, drxn );
//                return true;
//            }
//        }
//        
//        // Is one better, reliable and haploid? -- Set other unreliable, else set both unreliable
//        Node* spanBgn = NULL;
//        bool unreliable[2] = { 
//            !scores[0].reliable || !bestFar[0] || !scores[0].isHaploid, 
//            !scores[1].reliable || !bestFar[1] || !scores[1].isHaploid
//        };
//        for ( int i : { 0, 1 } )
//        {
//            if ( unreliable[i] || ( !unreliable[0] && !unreliable[1] ) )
//            {
//                for ( Node* node : scores[i].path )
//                {
//                    spanBgn = !spanBgn || ( node->isFurther( spanBgn->ends_[!drxn], !drxn, !drxn ) ) ? node : spanBgn;
//                    node->setUnreliable();
//                }
//            }
//        }
//        
//        // Update span and current path state
//        path.isReliable = false;
//        path.isMulti = true;
//        path.addSpan( spanBgn, path.path.back(), drxn );
//    }
//    
//    return false;
//}
//
//void Locus::plotPathConvergeCompareReliable( PathScore* scores, Path &path, bool* onlyReliable, bool* moreReliable, bool drxn )
//{
//    int32_t limits[2];
//    int32_t offset = ( params.maxPeMean - params.readLen ) * ( drxn ? -1 : 1 );
//    Node::getLimits( limits, scores[0].path, true );
//    Node::getLimits( limits, scores[1].path, false );
//    limits[0] += offset;
//    limits[1] += offset;
//    
//    NodeSet reliSet[2];
//    vector<Allele> allelePairs[2];
//    vector< pair<int32_t, int32_t> > reliableSpans[2];
//    
//    for ( int i : { 0, 1 } )
//    {
//        Node::getReliablePairNodes( scores[i].path, reliSet[i], !drxn );
//        path.getAllelesInSet( allelePairs[i], reliSet[i] );
//        for ( Path &revPath : paths_[!drxn] )
//        {
//            revPath.getAllelesInSet( allelePairs[i], reliSet[i] );
//        }
//        
//        for ( Node* node : reliSet[i] )
//        {
//            int32_t ends[2] = { 
//                max( limits[0], node->ends_[0] + node->getBestOverlap( 0 ) / 2 ), 
//                min( limits[1], node->ends_[1] - node->getBestOverlap( 1 ) / 2 )
//            };
//            
//            if ( ends[0] < ends[1] )
//            {
//                auto hit = reliableSpans[i].end();
//                for ( auto it = reliableSpans[i].begin(); it != reliableSpans[i].end(); )
//                {
//                    if ( ( ends[0] <= it->first && it->first <= ends[1] ) || ( ends[0] <= it->second && it->second <= ends[1] ) )
//                    {
//                        ends[0] = min( ends[0], it->first );
//                        ends[1] = max( ends[1], it->second );
//                        hit = hit == reliableSpans[i].end() ? it : hit;
//                        hit->first = ends[0], hit->second = ends[1];
//                        if ( hit != it )
//                        {
//                            it = reliableSpans[i].erase( it );
//                            continue;
//                        }
//                    }
//                    it++;
//                }
//                if ( hit == reliableSpans[i].end() )
//                {
//                    reliableSpans[i].push_back( make_pair( ends[0], ends[1] ) );
//                }
//            }
//        }
//    }
//    
//    for ( int i : { 0, 1 } )
//    {
//        for ( Allele &allele : allelePairs[i] )
//        {
//            onlyReliable[i] = onlyReliable[i] || !allele.anyInSet( reliSet[!i] );
//        }
//        
//        if ( reliSet[i].size() > reliSet[!i].size() && reliSet[!i].size() <= 1 )
//        {
//            int32_t len = reliableSpans[!i].empty() ? 0 : reliableSpans[!i][0].first - reliableSpans[!i][0].second;
//            vector< pair<int32_t, int32_t> > altSpans = reliableSpans[!i];
//            for ( pair<int32_t, int32_t> &span : reliableSpans[i] )
//            {
//                len += span.second - span.first;
//                for ( int j( 0 ); j < altSpans.size(); )
//                {
//                    if ( span.first <= altSpans[j].first && altSpans[j].first < span.second )
//                    {
//                        altSpans[j].first = span.second;
//                    }
//                    if ( span.first < altSpans[j].second && altSpans[j].second <= span.second )
//                    {
//                        altSpans[j].second = span.first;
//                    }
//                    if ( altSpans[j].first < span.first && span.second < altSpans[j].second )
//                    {
//                        altSpans.push_back( make_pair( altSpans[j].first, span.first ) );
//                        altSpans.push_back( make_pair( altSpans[j].second, span.second ) );
//                        altSpans.erase( altSpans.begin() + j );
//                        continue;
//                    }
//                    if ( altSpans[j].second <= altSpans[j].first )
//                    {
//                        altSpans.erase( altSpans.begin() + j );
//                        continue;
//                    }
//                    j++;
//                }
//            }
//            moreReliable[i] = altSpans.empty() && len >= params.readSpacing * 10 + 20;
//        }
//    }
//}
//
//bool Locus::plotPathConvergeConflict( PathScore &prime, PathScore &alt, Path &path, bool drxn )
//{
//    for ( auto it = path.convergents.begin(); it != path.convergents.end(); it++ )
//    {
//        if ( find( (*it).altPath.begin(), (*it).altPath.end(), prime.path[0] ) != (*it).altPath.end() 
//                || find( prime.path.begin(), prime.path.end(), (*it).altPath[0] ) != prime.path.end()
//                || find( (*it).path.begin(), (*it).path.end(), alt.path[0] ) != (*it).path.end() )
//        {
//            bool isBetter = (*it).furthest != alt.furthest 
//                    ? (*it).furthest > alt.furthest
//                    : (*it).reliable != alt.reliable
//                    ? (*it).reliable > alt.reliable
//                    : (*it).readCount > alt.readCount;
//            
//            NodeSet convSet( (*it).path.begin(), (*it).path.end() );
//            NodeSet altSet( alt.path.begin(), alt.path.end() );
//            
//            if ( isBetter && (*it).endFork != alt.endFork )
//            {
//                alt.path.back()->removeEdge( alt.endFork, drxn );
//                alt.endFork->removeEdge( alt.path.back(), !drxn );
//            }
//            else if ( isBetter && alt.bgnFork && find( (*it).path.begin()+1, (*it).path.end(), alt.path[0] ) != (*it).path.end() )
//            {
//                alt.bgnFork->removeEdge( alt.path[0], drxn );
//                alt.path[0]->removeEdge( alt.bgnFork, !drxn );
//            }
//            else if ( !isBetter && (*it).endFork != alt.endFork )
//            {
//                (*it).path.back()->removeEdge( (*it).endFork, drxn );
//                (*it).endFork->removeEdge( (*it).path.back(), !drxn );
//            }
//            else if ( !isBetter && (*it).bgnFork && find( alt.path.begin()+1, alt.path.end(), (*it).path[0] ) != alt.path.end() )
//            {
//                (*it).bgnFork->removeEdge( (*it).path[0], drxn );
//                (*it).path[0]->removeEdge( (*it).bgnFork, !drxn );
//            }
//            
//            deleteNodes( ( isBetter ? altSet : convSet ), ( isBetter ? convSet : altSet ), drxn );
//            return true;
//        }
//    }
//    return false;
//}
//
//bool Locus::plotPathConvergePath( PathVars pv, PathBranch &convBranch, Node* convEnd, BranchList &newDiv, Path &path, NodeList &convPath, NodeSet &convSet, bool drxn )
//{
//    NodeSet delSet;
//    BranchList branches = plotPathGetBranches( convPath[0], drxn );
//    pv.addTarget( convPath[0], convBranch.fwdSet );
//    PathBranch::setScores( pv, branches, path.spans );
//    PathBranch best = plotPathGetBestNext( convBranch, convEnd, branches, newDiv, convSet, path, delSet, drxn );
//    while( best.branch )
//    {
//        convPath.push_back( best.branch );
//        pv.addTarget( best.branch, best.fwdSet );
//        branches = plotPathGetBranches( best.branch, drxn );
//        PathBranch::setScores( pv, branches, path.spans );
//        best = plotPathGetBestNext( convBranch, convEnd, branches, newDiv, convSet, path, delSet, drxn );
//    }
//    
//    if ( !delSet.empty() )
//    {
//        deleteNodes( delSet, drxn );
//        return false;
//    }
//    return true;
//}
//
//bool Locus::plotPathConvergeRedundant( PathScore &alt, Path &path, bool drxn )
//{
//    NodeSet altSet( alt.path.begin(), alt.path.end() );
//    for ( auto it = path.convergents.begin(); it != path.convergents.end(); )
//    {
//        // Previous convergent set is redundant
//        if ( altSet.find( (*it).path[0] ) != altSet.end() )
//        {
//            bool doErase = true;
//            for ( Node* node : (*it).path )
//            {
//                doErase = doErase && altSet.find( node ) != altSet.end();
//            }
//
//            if ( doErase )
//            {
//                it = path.convergents.erase( it );
//                continue;
//            }
//        }
//        
//        // This convergent set is redundant
//        NodeSet convSet( (*it).path.begin(), (*it).path.end() );
//        if ( convSet.find( alt.path[0] ) != convSet.end() )
//        {
//            bool isRedundant = true;
//            for ( Node* node : alt.path )
//            {
//                isRedundant = isRedundant && convSet.find( node ) != convSet.end();
//            }
//            
//            if ( isRedundant )
//            {
//                return true;
//            }
//        }
//        
//        it++;
//    }
//    return false;
//}
//
//void Locus::plotPathConvergeScore( PathVars &pv, int32_t bgnCoord, int32_t bestOls[2], Node* postFork, PathScore &score, bool drxn )
//{
//    score.len = -( ( bestOls[0] + bestOls[1] ) / 2 );
//    score.len += bestOls[0] - score.path[0]->getBestOverlap( !drxn ) + bestOls[1] - score.path.back()->getOverlap( postFork, drxn );
//    
//    Node* prev = NULL;
//    int32_t furthest = bgnCoord;
//    for ( Node* node : score.path )
//    {
//        score.readCount += node->reads_.size() / ( node->clones_ ? 1 + node->clones_->size() : 1 );
//        score.len += node->seq_.length() - ( prev ? node->getOverlap( prev, !drxn ) : 0 );
//        auto it = pv.reliable.find( node );
//        score.reliable += it != pv.reliable.end() ? it->second : 0;
//        furthest = node->farPairNodes_[0] ? ( drxn ? min( furthest, node->farPairCoords_[0] ) 
//                                                   : max( furthest, node->farPairCoords_[0] ) ) 
//                                          : furthest;
//        prev = node;
//    }
//    
//    score.furthest = abs( bgnCoord - furthest );
//    
//    int32_t minLen = 5 + float( 5 * params.readLen ) / float( params.cover );
//    float coverage = float( score.readCount * params.readLen ) / float( max( minLen, score.len ) );
//    score.isHaploid = coverage < 0.7 * params.cover && score.len >= minLen;
//    score.notHaploid = coverage >= 0.7 * params.cover;
//}

PathBranch Locus::plotPathGetBestNext( BranchList &branches, bool isFirst )
{
    auto pBest = branches.end();
    
    for ( auto it = branches.begin(); it != branches.end(); it++ )
    {
        // Set new best
        if ( pBest == branches.end() || (*it).score[1] > (*pBest).score[1] )
        {
            pBest = it;
        }
    }
    
    PathBranch best;
    if ( pBest != branches.end() && ( isFirst || (*pBest).branch->isValidated() || (*pBest).hits > 4 ) )
    {
        best = *pBest;
        branches.erase( pBest );
        
        // Annotate converging vs diverging branches
        for ( PathBranch &branch : branches )
        {
            branch.altBranch = best.branch;
            branch.altScore = best.score;
            branch.state = 2;
            branch.spanned -= best.spanned;
            for ( Node* fwd : branch.fwdSet )
            {
                if ( best.fwdSet.find( fwd ) != best.fwdSet.end() )
                {
                    branch.state = 1;
                    break;
                }
            }
        }
    }
    return best;
}

PathBranch Locus::plotPathGetBestNext( PathBranch &convBranch, Node* convEnd, BranchList &branches, BranchList &newDiv, NodeSet &convSet, Path &path, NodeSet &delSet, bool drxn )
{
    PathBranch best;
    for ( PathBranch &branch : branches )
    {
        // Is this best convergent branch so far?
        if ( convSet.find( branch.branch ) != convSet.end() && ( !best.branch || branch.score[1] > best.score[1] ) )
        {
            best = branch;
        }
    }
    
    for ( PathBranch &branch : branches )
    {
        // This is a divergent branch
        if ( branch.branch != convEnd && convSet.find( branch.branch ) == convSet.end() )
        {
            branch.altBranch = convBranch.branch;
            branch.altScore = convBranch.score;
            branch.state = 2;
            
            bool isMainBranch = false;
            for ( Node* prv : branch.branch->getNextNodes( !drxn ) )
            {
                isMainBranch = isMainBranch || find( path.path.begin(), path.path.end(), prv ) != path.path.end();
            }
            
            if ( isMainBranch ) break;
            
            NodeSet fwdSet = convEnd->getDrxnNodes( drxn );
            for ( Node* fwd : branch.fwdSet )
            {
                if ( fwdSet.find( fwd ) != fwdSet.end() )
                {
                    branch.state = 3;
                    convBranch.state = 1;
                    break;
                }
            }
            newDiv.push_back( branch );
        }
        // This is a convergent path that is not the best
        else if ( branch.branch != convEnd && branch.branch != best.branch )
        {
            delSet.insert( branch.branch );
        }
    }
    return best;
}

BranchList Locus::plotPathGetBranches( PathVars &pv, bool drxn )
{
    BranchList branches;
    if ( forkNodes_[drxn].empty() || finalise_ )
    {
        setOriginEnds();
        forkNodes_[drxn] = originEnds_[drxn];
        reviewSpans( drxn );
    }
    NodeSet fwdSet;
    for ( Node* node : forkNodes_[drxn] )
    {
        branches.push_back( PathBranch( node, drxn ) );
        fwdSet.insert( node );
        node->getDrxnNodes( fwdSet, drxn, params.getFurthestMpDist( node->ends_[drxn], drxn ) );
        node->getDrxnNodes( pv.tSet, !drxn, params.getFurthestMpDist( node->ends_[!drxn], !drxn ) );
    }
    pv.setTarget( fwdSet, nodes_[drxn + 3] );
    return branches;
}

BranchList Locus::plotPathGetBranches( Node* target, bool drxn )
{
    BranchList branches;
    for ( Node* node : target->getNextNodes( drxn ) )
    {
        branches.push_back( PathBranch( node, target, drxn ) );
    }
    return branches;
}

void Locus::plotPathGetFirst( PathVars &pv, PathBranch &best, Path &path, BranchList &branches, bool drxn )
{
    branches = plotPathGetBranches( pv, drxn );
    PathBranch::setScores( pv, branches, path.spans );
    best = plotPathGetBestNext( branches, true );
    if ( best.branch )
    {
        path.fork = best.branch;

        // Set initial state of multiplicity and reliability
        path.isMulti = ( !path.spans.empty() && !path.spans.back().ended ) || ( best.branch && best.branch->isMultiple() );
        if ( path.isMulti && ( path.spans.empty() || path.spans.back().ended ) )
        {
            path.spans.push_back( Span( best.branch ) );
        }
        path.isReliable = path.isMulti || best.branch->isReliable( true );
    }
}

void Locus::plotPathGetNext( PathVars &pv, PathBranch &best, PathBranch &last, Path &path, BranchList &branches, bool drxn )
{
    last = best;
    pv.addTarget( last.branch, last.fwdSet );
    
    branches = plotPathGetBranches( best.branch, drxn );
    PathBranch::setScores( pv, branches, path.spans );
    best = plotPathGetBestNext( branches, false );
}

void Locus::plotPathSetEnds( PathVars &pv, PathBranch &best, PathBranch &last, Path &path, bool drxn )
{
    // Terminate path if invalid, irretrievably problematic or beyond extend limit
    if ( ( last.branch && ( drxn ? last.branch->ends_[1] >= params.locusLimits[1] 
                                 : last.branch->ends_[0] <= params.locusLimits[0] ) ) )
    {
        path.completed = true;
        finished_[drxn] = true;
    }
    else if ( pv.isInvalid )
    {
        path.completed = true;
    }
    // Set an end if last path node continues
    else if ( last.doesContinue )
    {
        completed_[drxn] = false;
        endNodes_[drxn].push_back( last.branch );
        desperation_[drxn] = 0;
    }
    // Terminate path if no more continuing and no divergent paths
    else if ( sideNodes_[drxn].empty() )
    {
        path.completed = true;
    }
    else if ( desperation_[drxn] < 3 )
    {
        completed_[drxn] = false;
        desperation_[drxn]++;
    }
    
    assert( sideNodes_[drxn].empty() || desperation_[drxn] >= 3 || !completed_[drxn] || finished_[drxn] || pv.weak );
}

//bool Locus::plotPathSetForks( PathVars &pv, Path &path, NodeSet &goodSet, NodeList &notReliableNodes, bool drxn )
//{
//    NodeSet badBranches( notReliableNodes.begin(), notReliableNodes.end() );
//    for ( PathBranch &div : path.divergent )
//    {
//        if ( div.state == 1 || goodSet.find( div.branch ) != goodSet.end() )
//        {
//            badBranches.insert( div.altBranch );
//        }
//    }
//    
//    if ( pv.misassembled && abs( path.path.back()->ends_[pv.drxn] - pv.misassembled->ends_[pv.drxn] ) > params.maxPeMean )
//    {
//        NodeSet delSet;
//        locusTest();
//        pv.misassembled->slice( pv, delSet, pv.drxn );
//        deleteNodes( delSet, drxn );
//        locusTest();
//        return false;
//    }
//    
////    if ( pv.weak )
////    {
////        for ( Node* node : path.path )
////        {
////            if ( !node->isContinue( drxn ) && node->stop_[pv.drxn] != 5
////                    && ( drxn ? pv.invalidCoord < node->ends_[1]
////                              : node->ends_[0] < pv.invalidCoord ) )
////            {
////                if ( !node->isMisassembled()
////                        && node->reassemble( pv, pv.drxn ) )
////                {
////                    locusTest();
////                    return false;
////                }
////                NodeSet delSet = pv.weak->getDrxnNodes( pv.drxn );
////                deleteNodes( delSet, pv.drxn );
////                pv.weak->stop_[pv.drxn] = 5;
////                return false;
////            }
////        }
////    }
//
//    for ( auto it = path.path.begin(); it != path.path.end(); it++ )
//    {
//        if ( ( !finalise_ && badBranches.find( *it ) != badBranches.end() ) || *it == pv.weak || *it == pv.misassembled ) 
//        {
//            path.path.erase( it, path.path.end() );
//            break;
//        }
//        else if ( (*it)->isReliable( false ) )
//        {
//            reliable_[pv.drxn] = pv.drxn ? max( reliable_[1], (*it)->ends_[1] ) 
//                                         : min( reliable_[0], (*it)->ends_[0] );
//        }
//    }
//    
//    for ( auto it = path.alleles.begin(); it != path.alleles.end(); )
//    {
//        if ( !(*it).complete && !(*it).paths[0].empty() )
//        {
//            if ( find( path.path.begin(), path.path.end(), (*it).paths[0].back() ) == path.path.end() )
//            {
//                it = path.alleles.erase( it, path.alleles.end() );
//                continue;
//            }
//            (*it).complete = true;
//        }
//        it++;
//    }
//    
//    path.fork = path.path.empty() ? path.fork : path.path.back();
//    if ( path.fork )
//    {
//        forkNodes_[drxn].clear();
//        forkNodes_[drxn].push_back( path.fork );
//        if ( path.path.empty() )
//        {
//            path.path.push_back( path.fork );
//        }
//    }
//    
//    return true;
//}

void Locus::plotPathTrimBranch( PathVars &pv, PathBranch &branch, NodeSet &goodSet, NodeSet &delSet, bool isEnd, bool drxn )
{
    NodeSetList endSets;
    NodeSet bestEndSet = { branch.branch };
    float bestScore = 2;
    
    for ( Node* fwd : branch.fwdSet )
    {
        if ( fwd->edges_[drxn].empty() )
        {
            NodeSet endSet = { fwd };
            fwd->getDrxnNodesInSet( endSet, branch.fwdSet, !drxn );
            endSets.push_back( endSet );
            float score = 0;
            for ( Node* node : endSet )
            {
                if ( isEnd )
                {
                    auto it = pv.hits.find( node );
                    score += ( it != pv.hits.end() ? it->second : 0 );
                }
                else
                {
                    auto it = pv.adjusted.find( node );
                    score += ( it != pv.adjusted.end() ? it->second : 0 );
                }
            }
            
            if ( score > bestScore )
            {
                bestEndSet = endSet;
                bestScore = score;
            }
        }
    }
    
    NodeSet bckSet;
    for ( Node* node : bestEndSet )
    {
        if ( pv.hits.find( node ) != pv.hits.end() )
        {
            bckSet.insert( node );
            node->getDrxnNodesInSet( bckSet, bestEndSet, !drxn );
        }
    }
    bestEndSet = bckSet;
    
    for ( NodeSet &endSet :endSets )
    {
        float readCount = 0;
        float score = 0, totalScore = 0;
        for ( Node* node : endSet )
        {
            float nodeScore;
            if ( isEnd )
            {
                auto it = pv.hits.find( node );
                nodeScore = ( it != pv.hits.end() ? it->second : 0 );
            }
            else
            {
                auto it = pv.adjusted.find( node );
                nodeScore = ( it != pv.adjusted.end() ? it->second : 0 );
            }
            totalScore += nodeScore;
            if ( bestEndSet.find( node ) == bestEndSet.end() )
            {
                readCount += node->reads_.size();
                score += nodeScore;
            }
        }
        readCount = min( readCount / 2, params.branchMinHits ) + bestScore - totalScore;
        score = max( (float)0, score - ( bestScore - totalScore ) );
        
        ( ( readCount - score * 3 > ( isEnd ? 22 : 12 ) ) ? delSet : goodSet ).insert( endSet.begin(), endSet.end() );
    }
}

void Locus::plotPathTrimDivergent( PathVars &pv, Path &path,NodeSet &goodSet, NodeSet &delSet, bool drxn )
{
    sort( path.divergent.begin(), path.divergent.end(), [&]( const PathBranch &a, const PathBranch &b ){
        return ( drxn ? a.branch->ends_[0] + a.overlap > b.branch->ends_[0] + b.overlap
            : a.branch->ends_[1] - a.overlap < b.branch->ends_[1] - b.overlap ); 
    } );
    
    int extended = 0;
    for ( PathBranch &div : path.divergent )
    {
        if ( div.state < 2 ) continue;
        
        float badScore = ( div.score[0] * ( ( div.altScore[1] + 10 ) / ( div.score[1] + 10 ) ) )                // base misses
                       + max( (float)0, (!div.doesContinue) * 20 - div.score[1] );                              // penalty for dead end
                       + ( 10 * ( float( params.readLen - div.overlap ) / float( params.readLen ) ) );          // penalty for weak edge
                       + ( harsh_[drxn] ? min( (float)0, div.score[1] - 10 - ( extended * 10 ) ) : 0 );         // desperate penalty
        
        // Delete bad branches
        if ( div.invalid || badScore > 16 )
        {
            delSet.insert( div.branch );
        }
        // Extend reasonable branches
        else
        {
            if ( div.fwdSet.size() > 1 )
            {
                plotPathTrimBranch( pv, div, goodSet, delSet, false, drxn );
            }
            sideNodes_[drxn].push_back( div.branch );
            goodSet.insert( div.branch );
            extended++;
        }
    }
}

void Locus::plotPathTrimEnds( PathVars &pv, PathBranch &best, BranchList &branches, NodeSet &goodSet, NodeSet &delSet, bool drxn )
{
    if ( best.branch )
    {
        branches.push_back( best );
    }
    
    for ( PathBranch &branch : branches )
    {
        if ( branch.score[1] * 3 + 12 > min( params.branchMinHits * 2, branch.score[0] ) )
        {
            goodSet.insert( branch.fwdSet.begin(), branch.fwdSet.end() );
        }
        else
        {
            if ( branch.fwdSet.size() == 1 )
            {
                delSet.insert( branch.branch );
            }
            else
            {
                goodSet.insert( branch.branch );
                plotPathTrimBranch( pv, branch, goodSet, delSet, true, drxn );
            }
        }
    }
}

void Locus::plotPathTrimPrep( PathVars &pv, Path &path, bool drxn )
{
    if ( path.path.empty() ) return;
    
    NodeSet divSet = path.path.back()->getDrxnNodes( drxn );
    for ( PathBranch &div : path.divergent )
    {
        divSet.insert( div.branch );
        div.branch->getDrxnNodes( divSet, drxn );
    }
    
    NodeSet fwdSet = path.path[0]->getDrxnNodesNotInSet( divSet, drxn, true );
    NodeSet bckSet = path.path.back()->getDrxnNodesInSet( fwdSet, !drxn, true );
    for ( Node* bck : bckSet )
    {
        for ( Node* nxt : bck->getNextNodes( drxn ) )
        {
            if ( bckSet.find( nxt ) == bckSet.end() && divSet.find( nxt ) == divSet.end() )
            {
                PathBranch div( nxt, drxn );
                div.setScores( pv );
            }
        }
    }
    
    for ( auto it = path.divergent.begin(); it != path.divergent.end(); )
    {
        if ( (*it).state == 3 )
        {
            for ( Node* fwd : (*it).fwdSet )
            {
                if ( pv.tSet.find( fwd ) != pv.tSet.end() )
                {
                    (*it).state = 1;
                    break;
                }
            }
            (*it).state = (*it).state == 1 ? 1 : 2;
        }
        if ( (*it).state == 0 )
        {
            it = path.divergent.erase( it );
            continue;
        }
        it++;
    }
}

void Locus::plotPathUpdateFork( PathVars &pv, Path &path, bool drxn )
{
    if ( path.fork )
    {
        forkNodes_[drxn].clear();
        forkNodes_[drxn].push_back( path.fork );
        if ( find( nodes_[drxn].begin(), nodes_[drxn].end(), path.fork ) == nodes_[drxn].end() )
        {
            forkNodes_[drxn].clear();
        }
    }
    if ( forkNodes_[drxn].empty() )
    {
        forkNodes_[drxn] = originEnds_[drxn];
        return;
    }
    
    if ( pv.newSet.empty() ) return;
    
    NodeSet goodSet( forkNodes_[drxn].begin(), forkNodes_[drxn].end() );
    for ( Node* fork : forkNodes_[drxn] )
    {
        if ( fork->drxn_ == 2 ) continue;
        fork->getDrxnNodes( goodSet, !drxn, true );
    }
    
    for ( Node* node : pv.newSet )
    {
        NodeSet bckSet = node->getDrxnNodes( !drxn, true, true );
        for ( auto it = goodSet.begin(); it != goodSet.end(); )
        {
            if ( bckSet.find( *it ) == bckSet.end() ) it = goodSet.erase( it );
            else it++;
        }
    }
    
    for ( auto it = forkNodes_[drxn].begin(); it != forkNodes_[drxn].end(); )
    {
        if ( goodSet.find( *it ) == goodSet.end() ) it = forkNodes_[drxn].erase( it );
        else it++;
    }
    
    if ( !forkNodes_[drxn].empty() ) return;
    if ( originEnds_[drxn].empty() ) setOriginEnds();
    goodSet.insert( originEnds_[drxn].begin(), originEnds_[drxn].end() );
    Node* forks[2] = { originEnds_[drxn][0], originEnds_[drxn][originEnds_[drxn].size() > 1] };
    Node* fork = originEnds_[drxn][0];
    bool complete = false;
    while ( !complete )
    {
        while ( forks[0] == forks[1] && !complete )
        {
            if ( goodSet.find( forks[0] ) == goodSet.end() ) complete = true;
            else if ( forks[0]->edges_[drxn].size() == 1 )
            {
                fork = forks[0];
                forks[0] = forks[0]->edges_[drxn][0].node;
                forks[1] = forks[0];
            }
            else if ( forks[0]->edges_[drxn].size() == 2 )
            {
                fork = forks[0];
                forks[0] = forks[0]->edges_[drxn][0].node;
                forks[1] = forks[1]->edges_[drxn][1].node;
            }
            else complete = true;
        }
        
        if ( complete ) break;
        
        for ( int i : { 0, 1 } )
        {
            while ( forks[i]->edges_[0].size() == 1 && forks[i]->edges_[1].size() == 1 )
            {
                forks[i] = forks[i]->edges_[drxn][0].node;
            }
        }
        complete = complete || forks[0] != forks[1];
    }
    
    forkNodes_[drxn].push_back( fork );
}

void Locus::plotPathUpdateMultiplicity( PathBranch &best, Path &path, NodeList &notReliableNodes, bool drxn )
{
    if ( !path.isMulti && best.branch->isMultiple() )
    {
        notReliableNodes.clear();
        path.spans.push_back( Span( best.branch ) );
        path.isMulti = true;
    }
    else if ( path.isMulti )
    {
        path.completeSpan( best.branch, drxn );
    }
}

void Locus::plotPathUpdateReliability( PathBranch &best, Path &path, NodeList &notReliableNodes, bool drxn )
{
    path.isReliable = path.isMulti || best.branch->isReliable( true );
    if ( path.isReliable && !notReliableNodes.empty() )
    {
        setReliable( notReliableNodes, drxn );
        notReliableNodes.clear();
    }
    else if ( !path.isReliable )
    {
        notReliableNodes.push_back( best.branch );
    }
}

void Locus::plotPathUpdateSpans( PathBranch &best, Path &path, bool drxn )
{
    for ( Span &span : path.spans )
    {
        span.spanned = span.spanned || ( span.ended && best.valid );
        if ( !best.valid  && !span.spanned
                && abs( best.branch->ends_[drxn] - span.node->ends_[drxn] ) > params.maxPeMax ) best.invalid = true;
    }
}

void Locus::plotPrep( bool drxn )
{
    NodeSet qSet( forkNodes_[drxn].begin(), forkNodes_[drxn].end() );
    if( qSet.empty() )
    {
        qSet.insert( endNodes_[drxn].begin(), endNodes_[drxn].end() );
    }
    if( qSet.empty() )
    {
        qSet.insert( sideNodes_[drxn].begin(), sideNodes_[drxn].end() );
    }
    if( qSet.empty() )
    {
        qSet.insert( originEnds_[drxn].begin(), originEnds_[drxn].end() );
    }
    Node::propagateValidation( qSet, validLimits_, drxn );
}

void Locus::reviewSpans( bool drxn )
{
    NodeSet forkFwd = { forkNodes_[drxn].begin(), forkNodes_[drxn].end() };
    
    for ( Path &path : paths_[drxn] )
    {
        for ( int i = 0; i < path.spans.size(); )
        {
            if ( forkFwd.find( path.spans[i].node ) != forkFwd.end()
                    || find( nodes_[drxn].begin(), nodes_[drxn].end(), path.spans[i].node ) == nodes_[drxn].end() )
            {
                path.spans.erase( path.spans.begin() + i );
            } else i++;
        }
    }
}

void Locus::setReliable( NodeList &path, bool drxn )
{
    float len = path[0]->seq_.length() - ( path[0]->getBestOverlap( !drxn ) + path.back()->getBestOverlap( drxn ) ) / 2;
    float readCount = path[0]->reads_.size();
    for ( int i( 1 ); i < path.size(); i++ )
    {
        len += path[i]->seq_.length() - path[i]->getOverlap( path[i - 1], !drxn );
        readCount += path[i]->reads_.size();
    }
    float coverage = float( readCount * params.readLen ) / max( len, float( 1 ) );
    if ( coverage < params.cover * 1.3 + min( float(0.4), len / float( 20 * params.readLen ) ) )
    {
        for ( Node* node : path )
        {
            node->setReliable( false );
        }
    }
}

bool Locus::updateForks( bool drxn )
{
    bool newBranch = !forkNodes_[drxn].empty();
    NodeList branches;
    NodeSet forkFwdSet;
    for ( Node* fork : forkNodes_[drxn] )
    {
        if ( find( nodes_[drxn].begin(), nodes_[drxn].end(), fork ) != nodes_[drxn].end() )
        {
            branches.push_back( fork );
            forkFwdSet.insert( fork );
            fork->getDrxnNodes( forkFwdSet, drxn );
        }
    }
    
    // Find ends that are not reachable by current forks
    for ( Node* node : nodes_[drxn] )
    {
        if ( node->edges_[drxn].empty() && forkFwdSet.find( node ) == forkFwdSet.end() )
        {
            newBranch = true;
            branches.push_back( node );
        }
    }
    
    if ( !newBranch )
    {
        return false;
    }
    
    NodeList tmp = branches;
    
    while ( branches.size() > 1 )
    {
        bool doErase = false;
        NodeSet fwdSet = branches[0]->getDrxnNodes( drxn );
        
        // Delete any branches that are forwards reachable, or delete this if any are backwards reachable
        for ( int i ( 1 ); i < branches.size(); )
        {
            if ( fwdSet.find( branches[i] ) != fwdSet.end() )
            {
                NodeSet bckSet = branches[i]->getDrxnNodesNotInSet( fwdSet, !drxn, true );
                bool iErase = false;
                for ( Node* nxt : branches[0]->getNextNodes( drxn ) )
                {
                    iErase = iErase || bckSet.find( nxt ) == bckSet.end();
                }
                if ( !iErase )
                {
                    doErase = true;
                    break;
                }
                branches.erase( branches.begin() + i );
                continue;
            }
            i++;
        }
        
        // Find nodes back from this that are not backwards reachable from this
        if ( !doErase )
        {
            NodeSet bckSet = branches[0]->getDrxnNodes( !drxn, true );
            for ( Node* bck : bckSet )
            {
                for ( Node* nxt : bck->getNextNodes( drxn ) )
                {
                    if ( bckSet.find( nxt ) == bckSet.end() && nxt != branches[0] )
                    {
                        doErase = true;
                        if ( find( branches.begin(), branches.end(), bck ) == branches.end() )
                        {
                            branches.push_back( bck );
                        }
                    }
                }
            }
        }
        
        if ( doErase )
        {
            branches.erase( branches.begin() );
        }
        else if ( branches.size() > 1 )
        {
            branches.clear();
        }
    }
    
    if ( !branches.empty() )
    {
        NodeSet fwdSet = branches[0]->getDrxnNodes( drxn, false, true );
        NodeSet bckSet = branches[0]->getDrxnNodes( !drxn, true, true );
        for ( Node* node : tmp )
        {
            assert( fwdSet.find( node ) != fwdSet.end() || bckSet.find( node ) != bckSet.end() || node->drxn_ == 2 );
        }
    }
    
    
    forkNodes_[drxn] = branches.empty() ? originEnds_[drxn] : branches;
    
    reviewSpans( drxn );
    
    return true;
}
