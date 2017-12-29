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

#include <algorithm>

#include "locus_pathing_structs.h"

Span::Span( Node* bgn )
: node( bgn ), spanned( false ), complete( false ), ended( false )
{}

PathBranch::PathBranch()
{
    branch = NULL;
    altBranch = NULL;
    fork = NULL;
    overlap = 0;
    reliable = 0;
    hits = 0;
    adjusted = 0;
    valid = true;
    invalid = false;
    doesContinue = false;
    state = 0;
    score[0] = 0, score[1] = 0, altScore[0] = 0, altScore[1] = 0;
}

PathBranch::PathBranch( Node* node, bool drxn )
: PathBranch()
{
    fwdSet.insert( node );
    node->getDrxnNodes( fwdSet, drxn );
    branch = node;
}

PathBranch::PathBranch( Node* node, Node* prev, bool drxn )
: PathBranch( node, drxn )
{
    fork = prev;
}

void PathBranch::setAlt( PathBranch &alt )
{
    altBranch = alt.branch;
    altScore[0] = alt.score[0];
    altScore[1] = alt.score[1];
}

void PathBranch::setScores( PathVars &pv, BranchList &branches, vector<Span> &spans )
{
    int32_t minStart = pv.drxn ? numeric_limits<int32_t>::max() : numeric_limits<int32_t>::min();
    int32_t maxFar[2] = { int32_t( ( params.maxPeMean - params.readLen ) * 0.8 ), params.readLen / 2 };
    int maxSpans = 0;
    float maxReliable = 0;
    
    for ( PathBranch &branch : branches )
    {
        branch.setScores( pv );
        branch.setSpans( spans, pv.drxn );
        minStart = pv.drxn ? min( minStart, branch.branch->ends_[0] ) : max( minStart, branch.branch->ends_[1] );
        maxReliable = max( maxReliable, branch.reliable );
    }
    
    for ( PathBranch &branch : branches )
    {
        branch.farCoords[0] = max( 0, ( pv.drxn ? minStart - branch.farCoords[0] : branch.farCoords[0] - minStart ) );
        branch.farCoords[1] = max( 1, ( pv.drxn ? minStart - branch.farCoords[1] : branch.farCoords[1] - minStart ) );
        if ( branch.reliable + 2 > maxReliable / 2 )
        {
            maxFar[0] = max( maxFar[0], branch.farCoords[0] );
            maxFar[1] = max( maxFar[1], branch.farCoords[1] );
            maxSpans = max( maxSpans, branch.spanned );
        }
    }
    
    maxFar[0] = min( maxFar[0], ( params.maxPeMean - params.readLen ) );
    maxFar[1] = min( maxFar[1], int32_t( ( params.maxPeMean - params.readLen ) * 0.8 ) );
    
    for ( PathBranch &branch : branches )
    {
        float farMod = min( float( branch.farCoords[0] ) / float( maxFar[0] ), 
                            float( branch.farCoords[1] + params.readLen ) / float( maxFar[1] ) );
        farMod = farMod > 1 ? sqrt( farMod ) : pow( farMod, 2 );
        float spanMod = sqrt( float( 1 + maxSpans ) / float( 1 + branch.spanned ) );
        float reliMod = 0.3 + ( sqrt( ( branch.reliable + 1.0 ) / ( maxReliable + 1.0 ) ) * 0.9 ); 
        float totalMod = max( reliMod * farMod * spanMod, float( 10 ) / float( 10 + branch.reads ) );
        if ( totalMod < 1 && branch.farNodes[0] && branch.farNodes[0]->drxn_ != branch.branch->drxn_ ) totalMod = 1;
        branch.score[1] = branch.adjusted * totalMod + branch.islands;
        branch.score[0] = max( (float)0, ( branch.reads - branch.score[1] ) ) / totalMod;
    }
}

void PathBranch::setScores( PathVars &pv )
{
    hits = reliable = adjusted = reads = 0;
    farNodes[0] = farNodes[1] = branch;
    farCoords[0] = farCoords[1] = branch->ends_[!pv.drxn];
    
    for ( Node* fwd : fwdSet )
    {
        doesContinue = doesContinue || fwd->isContinue( pv.drxn );
        
        if ( reads < params.branchMinHits * 4 )
        {
            reads += ( fwd->reads_.size() / ( fwd->coverage_ < params.cover ? (float)1 : fwd->coverage_ / params.cover ) );
        }
        
        // Set hits and furthest
        auto it = pv.hits.find( fwd );
        if ( it != pv.hits.end() ) 
        {
            hits += it->second;
            if ( !fwd->farPairNodes_[0] || fwd->drxn_ == 2 )
            {
                fwd->setFurthest( pv.tSet, pv.drxn );
            }
            
            if ( fwd->farPairNodes_[0] && ( pv.drxn ? fwd->farPairCoords_[0] < farCoords[0]
                                                    : farCoords[0] < fwd->farPairCoords_[0] ) )
            {
                farNodes[0] = fwd->farPairNodes_[0];
                farCoords[0] = fwd->farPairCoords_[0];
            }
        }
        
        // Set reliable hits and furthest reliable
        it = pv.reliable.find( fwd );
        if ( it != pv.reliable.end() ) 
        {
            reliable += it->second;
            if ( !fwd->farPairNodes_[1] )
            {
                fwd->setFurthest( pv.tSet, pv.drxn );
            }
            
            if ( fwd->farPairNodes_[1] && ( pv.drxn ? fwd->farPairCoords_[1] < farCoords[1]
                                                    : farCoords[1] < fwd->farPairCoords_[1] ) )
            {
                farNodes[1] = fwd->farPairNodes_[1];
                farCoords[1] = fwd->farPairCoords_[1];
            }
        }
        
        // Set adjusted hits
        auto it2 = pv.adjusted.find( fwd );
        if ( it2 != pv.adjusted.end() ) adjusted += it2->second;
        
        // Set island hits
        it = pv.islandHits.find( fwd );
        if ( it != pv.islandHits.end() ) islands += it->second;
    }
    
    reads = min( reads / 2, params.branchMinHits * 2 );
}

void PathBranch::setSpans( vector<Span> &spans, bool drxn )
{
    spanned = 0;
    if ( farNodes[0] )
    {
        NodeSet farFwdSet = farNodes[0]->getDrxnNodesNotInSet( fwdSet, drxn );
        for ( Span &span : spans )
        {
            if ( farFwdSet.find( span.node ) != farFwdSet.end() )
            {
                spanned++;
            }
            else
            {
                valid = valid && span.spanned;
            }
        }
    }
}

Allele::Allele( PathScore* scores )
{
    paths[0] = scores[0].path;
    paths[1] = scores[1].path;
    complete = false;
}

bool Allele::anyInSet( NodeSet &nodes )
{
    bool found[2] = { false, false };
    anyInSet( nodes, found );
    return found[0] && found[1];
}

void Allele::anyInSet( NodeSet &nodes, bool* found )
{
    for ( int i : { 0, 1 } )
    {
        for ( Node* node : paths[i] )
        {
            if ( nodes.find( node ) != nodes.end() )
            {
                found[i] = true;
                break;
            }
        }
    }
}

//void Path::addSpan( Node* bgn, Node* nd, bool drxn )
//{
//    bool doAdd = true;
//    for ( auto it = spans.begin(); it != spans.end(); )
//    {
//        if ( bgn->isFurther( (*it).begin->ends_[!drxn], !drxn, !drxn ) )
//        {
//            spans.erase( it, spans.end() );
//            break;
//        }
//        if ( !(*it).end || bgn->isFurther( (*it).end->ends_[!drxn], !drxn, !drxn ) )
//        {
//            doAdd = false;
//            if ( (*it).end && nd->isFurther( (*it).end->ends_[drxn], drxn, drxn ) )
//            {
//                (*it).end = NULL;
//            }
//            break;
//        }
//        it++;
//    }
//    if ( doAdd )
//    {
//        spans.push_back( Span( bgn ) );
//    }
//}

void Path::completeSpan( Node* node, bool drxn )
{
    bool doEnd = node->ends_[1] - node->ends_[0] > params.readLen * 2;
    NodeSet reliSet, reliFwd;
    while ( node )
    {
        if ( !node->getReliability() && ( reliSet.empty() || node->coverage_ > params.cover * 1.1 ) ) break;
        reliSet.insert( node );
        if ( node->edges_[drxn].size() == 1 ) node = node->edges_[drxn][0].node;
        else
        {
            for ( Node* nxt : node->getNextNodes( drxn ) )
            {
                if ( nxt->getReliability() ) reliSet.insert( node );
            }
            node = NULL;
        }
    }
    
    if ( reliSet.empty() ) return;
    
    for ( Node* reli : reliSet )
    {
        if ( reli->farPairNodes_[0] ) reli->farPairNodes_[0]->getDrxnNodes( reliFwd, drxn );
    }
    
    if ( reliFwd.empty() ) return;
    
    bool allEnded = true;
    for ( Span &span : spans )
    {
        if ( span.ended && span.spanned ) continue;
        if ( reliFwd.find( span.node ) != reliFwd.end() ) span.spanned = true;
        if ( span.spanned ) span.ended = true;
        if ( !span.ended ) span.ended = doEnd;
        allEnded = allEnded && span.ended;
    }
    
    if ( allEnded ) isMulti = false;
}

void Path::getAllelesInSet( vector<Allele> &rAlleles, NodeSet &nodes )
{
    for ( Allele &allele : alleles )
    {
        bool found[2] = { false, false };
        allele.anyInSet( nodes, found );
        if ( found[0] && found[1] )
        {
            rAlleles.push_back( allele );
        }
    }
}

void Path::reset( NodeList forks, bool drxn )
{
    fork = NULL;
    convergents.clear();
    convFork.clear();
    divergent.clear();
    path.clear();
    
    NodeSet bckSet;
    for ( Node* node : forks )
    {
        if ( node->drxn_ == 2 ) break;
        node->getDrxnNodes( bckSet, !drxn );
    }
    
    for ( int i = 0; i < spans.size(); )
    {
        if ( bckSet.find( spans[i].node ) == bckSet.end() )
        {
            spans.erase( spans.begin() + i );
        }
        else i++;
    }
}