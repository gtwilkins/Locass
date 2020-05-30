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

//bool Node::checkMisassembly( PathVars &pv )
//{
//    if ( drxn_ == 2 ) return false;
//    if ( !misassembled_[pv.drxn] && !isContinue( pv.drxn ) && isReliable( false )
//            && isMisassembled( pv, pv.drxn ) )
//    {
//        int32_t dist = params.maxPeMean + params.readLen;
//        if ( reassemble( pv, this, pv.misassEst, dist, false, pv.drxn ) ) return true;
//        pv.misassembled = this;
//    }
//    
//    return false;
//}

//bool Node::checkMisassembly( PathVars &pv, NodeList* paths )
//{
////    if ( isMisassembled( pv, paths, pv.drxn ) )
////    {
////        int32_t dist = params.maxPeMean + params.readLen;
////        if ( reassemble( pv, this, pv.misassEst, dist, false, pv.drxn ) ) return true;
////        
////        if ( isMisassembled( pv, paths, pv.drxn ) )
////        {
////            pv.misassembled = this;
////        }
////    }
//    
//    return false;
//}

//vector<PathSeq> Node::getCombinedSeqs( int dist, bool drxn )
//{
//    vector<PathSeq> nss( 1, PathSeq( this ) );
//    NodeSet usedNodes, branchNodes;
//    int i = 0;
//    while ( i < nss.size() )
//    {
//        while ( nss[i].dist < dist && nss[i].nodes.back()->drxn_ < 2 )
//        {
//            vector<Edge> unusedEdges;
//            nss[i].nodes.back()->sortEdges( drxn );
//            for ( Edge &e : nss[i].nodes.back()->edges_[drxn] )
//            {
//                if ( usedNodes.find( e.node ) == usedNodes.end() ) unusedEdges.push_back( e );
//            }
//            Edge* be = !unusedEdges.empty() ? &unusedEdges[0] : NULL;
//            if ( unusedEdges.size() > 1 )
//            {
//                nss.push_back( PathSeq( nss[i] ) );
//            }
//            else if ( unusedEdges.empty() && !nss[i].nodes.back()->edges_[drxn].empty() )
//            {
//                be = &nss[i].nodes.back()->edges_[drxn][0];
//                for ( int j = 0; j < nss[i].nodes.back()->edges_[drxn].size(); j++ )
//                {
//                    for ( Node* fwd : nss[i].nodes.back()->getDrxnNodes( drxn, true, false ) )
//                    {
//                        if ( branchNodes.find( fwd ) != branchNodes.end() )
//                        {
//                            be = &nss[i].nodes.back()->edges_[drxn][j];
//                            break;
//                        }
//                    }
//                }
//                assert( be );
//            }
//
//            if ( be )
//            {
//                usedNodes.insert( be->node );
//                branchNodes.erase( be->node );
//                int extLen = be->node->ends_[1] - be->node->ends_[0];
//                int seqLen = be->node->seq_.length();
//                if ( be->overlap > 0 ) extLen -= be->overlap;
//                string ext = drxn ? be->node->seq_.substr( max( 0, be->overlap ) )
//                                  : be->node->seq_.substr( 0, extLen );
//                if ( be->overlap < 0 )
//                {
//                    extLen -= be->overlap;
//                    ext = ( drxn ? string( -be->overlap, 'N' ) + ext : ext + string( -be->overlap, 'N' ) );
//                }
//                int32_t nodeEnds[2];
//                nodeEnds[0] = nss[i].ends[drxn] - ( drxn ? be->overlap : extLen );
//                nodeEnds[1] = nss[i].ends[drxn] + ( drxn ? extLen : be->overlap );
//                nss[i].seq = drxn ? nss[i].seq + ext : ext + nss[i].seq;
//                nss[i].nodes.push_back( be->node );
//                nss[i].nodeCoords.push_back( Coords( nodeEnds[0], nodeEnds[1], false ) );
//                nss[i].dist += extLen;
//                nss[i].ends[drxn] = drxn ? nss[i].ends[1] + extLen : nss[i].ends[0] - extLen;
//            }
//            else
//            {
//                nss[i].dist + dist;
//                assert( false );
//            }
//            if ( nss[i].dist > dist )
//            {
//                int trim = min( (int)nss[i].seq.length(), nss[i].dist - dist );
//                nss[i].seq = drxn ? nss[i].seq.substr( 0 , (int)nss[i].seq.length() - trim )
//                                   : nss[i].seq.substr( trim );
//                nss[i].ends[drxn] += drxn ? -trim : trim;
//                nss[i].dist -= trim;
//            }
//        }
//        i++;
//        bool nxtRead = false;
//        while ( !nxtRead && i < nss.size() )
//        {
//            for ( Node* nxt : nss[i].nodes.back()->getNextNodes( drxn ) )
//            {
//                nxtRead = nxtRead || usedNodes.find( nxt ) == usedNodes.end();
//            }
//            if ( !nxtRead )
//            {
//                nss.erase( nss.begin() + i );
//            }
//        }
//    }
//    
//    for ( PathSeq &ns : nss )
//    {
//        if ( !drxn )
//        {
//            reverse( ns.nodes.begin(), ns.nodes.end() );
//            reverse( ns.nodeCoords.begin(), ns.nodeCoords.end() );
//        }
//    }
//    
//    return nss;
//}

bool Node::getMisassembled( PathVars &pv, bool drxn )
{
    NodeSet currSet = { this };
    Node* curr = this;
    vector<int32_t> marks, ests;
    while ( curr && curr->isReliable() && abs( curr->ends_[pv.drxn] - ends_[!pv.drxn] ) < params.maxPeMean )
    {
        for ( ReadMark &mark : curr->marks_[pv.drxn] )
        {
            if ( !params.isReadPe( mark.id ) ) continue;
            if ( curr != this && ( drxn ? mark.estimate < ends_[0] 
                                        : ends_[1] < mark.estimate ) ) continue;
            ests.push_back( mark.estimate );
        }
        curr = curr->edges_[pv.drxn].size() == 1 ? curr->edges_[pv.drxn][0].node : NULL;
        currSet.insert( curr );
    }
    
    if ( ests.empty() ) return false;
    
    sort( ests.begin(), ests.end() );
    int32_t coord = ests[0];
    int maxCount = ests.size() > 1;
    int len = ests.size() > 1 ? ests[1] - ests[0] : 1;
    for ( int i = 0; i + maxCount < ests.size(); i++ )
    {
        if ( ests[i+maxCount] - ests[i] )
        {
            coord = ests[i];
            len = ests[i+maxCount] - ests[i];
        }
        int j = maxCount + 1;
        while ( ests[i+j] - ests[i] <= params.readLen  )
        {
            coord = ests[i];
            maxCount = j;
            len = ests[i+j] - ests[i];
            j++;
        }
    }
    pv.misassEst[0] = coord;
    pv.misassEst[1] = coord + len;
    
    return true;
}

bool Node::isMisassembled()
{
    return misassembled_[0] || misassembled_[1];
}

bool Node::isMisassembled( PathVars &pv, bool drxn )
{
    if ( assembled_[drxn] || !reliable_ ) return false;
    
    int32_t peWindow = max( params.readLen, params.avgPeMean / 2 );
    vector<int32_t> marks, ests;
    
    NodeSet currSet;
    Node* curr = this;
    while ( curr && curr->isReliable() )
    {
        currSet.insert( curr );
        curr = ( abs( curr->ends_[drxn] - ends_[drxn] ) < peWindow 
                && curr->edges_[drxn].size() == 1 ) ? curr->edges_[drxn][0].node : NULL;
    }
    
    int32_t limits[2] = { ends_[0] - peWindow, ends_[1] + peWindow };
    for ( Node* node : currSet )
    {
        for ( ReadMark &mark : node->marks_[drxn] )
        {
            if ( params.isReadPe( mark.id ) && limits[0] < mark.mark && mark.mark < limits[1] )
            {
                marks.push_back( mark.mark );
                ests.push_back( mark.estimate );
            }
        }
    }
    
    if ( isMisassembled( marks, ests, pv.misassMark, pv.misassEst, pv.drxn )
            && ( drxn ? pv.misassEst[0] <= ends_[1] : ends_[0] <= pv.misassEst[1] ) )
    {
        return true;
    }
    
    assembled_[drxn] = true;
    return false;
}

bool Node::isMisassembled( PathVars &pv, Node* farFork, NodeSet &forkSet, bool drxn )
{
    if ( !reliable_ || assembled_[drxn] ) return false;
    
    vector<int32_t> marks, ests;
    int32_t cutoff = ends_[drxn] + ( drxn ? params.maxPeMean : -params.maxPeMean );
    Node* curr = farFork->isReliable() ? farFork : NULL;
    while ( curr && ( drxn ? farFork->ends_[0] < cutoff : cutoff < farFork->ends_[1] ) )
    {
        forkSet.insert( curr );
        if ( curr->edges_[drxn].size() == 1 ) curr = curr->edges_[drxn][0].node;
        else
        {
            for ( Node* nxt : curr->getNextNodes( drxn ) )
            {
                if ( drxn ? nxt->ends_[0] < cutoff 
                          : cutoff < nxt->ends_[1] ) forkSet.insert( nxt );
            }
            break;
        }
    }
    
    for ( Node* node : forkSet )
    {
        if ( drxn ? cutoff < node->ends_[0] : node->ends_[1] < cutoff ) continue;
        for ( ReadMark &mark : node->marks_[drxn] )
        {
            if ( drxn ? cutoff < mark.mark : mark.mark < cutoff ) continue;
            if ( params.isReadPe( mark.id ) )
            {
                marks.push_back( mark.mark );
                ests.push_back( mark.estimate );
            }
        }
    }
    
    if ( isMisassembled( marks, ests, pv.misassMark, pv.misassEst, pv.drxn )
            && ( drxn ? pv.misassMark[0] - params.readLen < farFork->ends_[0] 
                      : farFork->ends_[1] < pv.misassMark[1] + params.readLen ) )
    {
        return true;
    }
    
    return false;
}

bool Node::isMisassembled( vector<int32_t> &marks, vector<int32_t> &ests, int32_t markLimits[2], int32_t estLimits[2], bool drxn )
{
    int32_t peWindow = max( params.readLen, params.avgPeMean / 2 );
    int32_t peWindowCutoff = max( float( 3 ), float( peWindow * params.peCover ) / float( params.readLen * 6 ) );
    
    bool misassembled[2]{false};
    for ( int i : { 0, 1 } )
    {
        vector<int32_t>* reads = ( i ? &marks : &ests );
        ( drxn ? sort( reads->begin(), reads->end() ) :  sort( reads->rbegin(), reads->rend() ) );
        int32_t** estMark = ( i ? &markLimits : &estLimits );
        for ( int j = 0; j + peWindowCutoff < reads->size() && !misassembled[i]; j++ )
        {
            int k = j + peWindowCutoff;
            while ( k < reads->size() 
                    && abs( (*reads)[k] - (*reads)[j] ) < peWindow )
            {
                (*estMark)[!drxn] = (*reads)[j];
                (*estMark)[drxn] = (*reads)[k++];
                misassembled[i] = true;
            }
        }
    }
    
    if ( misassembled[0] && !misassembled[1] )
    {
        estLimits[0] = markLimits[0] + ( drxn ? -params.avgPeMean : params.avgPeMean );
        estLimits[1] = markLimits[1] + ( drxn ? -params.avgPeMean : params.avgPeMean );
    }
    if ( !misassembled[0] && misassembled[1] )
    {
        markLimits[0] = estLimits[0] + ( drxn ? params.avgPeMean : -params.avgPeMean );
        markLimits[1] = estLimits[1] + ( drxn ? params.avgPeMean : -params.avgPeMean );
    }
    
    
    return misassembled[0] || misassembled[1];
}

bool Node::isMisassembledRev( int32_t markCoords[2], bool drxn )
{
    int32_t limits[2];
    limits[0] = drxn ? markCoords[0] - params.maxPeMean : markCoords[0] + params.readLen;
    limits[1] = drxn ? markCoords[1] - params.readLen : markCoords[1] + params.maxPeMean;
    NodeSet bckSet = { this };
    getDrxnNodes( bckSet, !drxn, limits[!drxn] );
    vector<int32_t> marks;
    
    for ( Node* node : bckSet )
    {
        node->resetUnmarked( !drxn );
        for ( ReadMark const &mark : node->marks_[!drxn] )
            if ( limits[0] <= mark.mark && mark.mark <= limits[1] )
                marks.push_back( mark.mark );
    }
    sort( marks.begin(), marks.end() );
    
    int32_t peWindow = max( params.readLen, params.avgPeMean / 2 );
    int32_t peWindowCutoff = max( float( 3 ), float( peWindow * params.peCover ) / float( params.readLen * 6 ) );
    
    for ( int i = 0; i + peWindowCutoff < marks.size(); i++ )
        if ( abs( marks[i] - marks[i + peWindowCutoff] ) < peWindow )
            return true;
    
    return false;
}

void Node::mapRead( NodeMapRead &nmr, NodeSet &tSet )
{
    for ( bool drxn : { 0, 1 } )
    {
        int seqLen = 16;
        if ( seqLen > nmr.seq.length() ) break;
        string q = drxn ? nmr.seq.substr( nmr.seq.length() - seqLen ) : nmr.seq.substr( 0, seqLen );
        for ( Node* node : tSet )
        {
            size_t it = node->seq_.find( q );
            while ( it != node->seq_.npos )
            {
                if ( drxn )
                {
                    int i = q.length();
                    nmr.coords[1][0] = node->ends_[0] + it;
                    it += q.length();
                    while ( i < nmr.seq.length() && it < node->seq_.length() && nmr.seq[i++] == node->seq_[it++] )
                    {
                        seqLen++;
                    }
                    nmr.coords[1][1] = nmr.coords[1][0] + seqLen;
                    nmr.nodes[1] = node ;
                    nmr.lens[1] = seqLen++;
                }
                else
                {
                    int i = nmr.seq.length() - q.length();
                    nmr.coords[0][1] = node->ends_[0] + it + q.length();
                    while ( i >= 0 && it >= 0 && nmr.seq[--i] == node->seq_[--it] )
                    {
                        seqLen++;
                    }
                    nmr.coords[0][0] = nmr.coords[0][1] - seqLen;
                    nmr.nodes[0] = node ;
                    nmr.lens[0] = seqLen++;
                }
                if ( seqLen > nmr.seq.length() ) break;
                q = drxn ? nmr.seq.substr( nmr.seq.length() - seqLen ) : nmr.seq.substr( 0, seqLen );
                it = node->seq_.find( q );
            }
        }
    }
}

//bool Node::reassemble( PathVars &pv, bool drxn )
//{
//    if ( misassembled_[drxn] || !pv.weak ) return false;
//    int32_t dist = abs( ends_[!drxn] - pv.weak->ends_[!drxn] ) + params.maxPeMean + params.readLen;
//    int32_t estimates[2];
//    
//    return reassemble( pv, pv.weak, estimates, dist, true, drxn );
//}
//
//bool Node::reassemble( PathVars &pv, Node* targ, int32_t estimates[2], int32_t dist, bool remap, bool drxn )
//{
//    vector<PathSeq> pss = getCombinedSeqs( dist, !drxn );
//    if ( remap )
//    {
//        estimates[0] = pss[0].ends[0] + params.readLen;
//        estimates[1] = pss[0].ends[1] - params.readLen;
//    }
//    PathSeq::setWeakspot( pss, estimates );
//    PathSeq::map( pv, pss, targ, drxn );
//
//    unordered_set<ReadId> usedIds;
//    int bestCount, bestOffset;
//    int iBest = PathSeq::getBest( pv, pss, bestCount, bestOffset, remap );
//    while ( iBest != -1 )
//    {
//        pss[iBest].setEdges();
//        if ( pss[iBest].tryBridge( pv, bestCount, bestOffset, drxn ) ) return true;
//        if ( !pss[iBest].tryMap( pv, usedIds, remap, drxn ) ) pss[iBest].exhausted = true;
//        else
//        {
//            pss = getCombinedSeqs( dist, !drxn );
//            PathSeq::map( pv, pss, targ, drxn );
//            PathSeq::setWeakspot( pss, estimates );
//        }
//        iBest = PathSeq::getBest( pv, pss, bestCount, bestOffset, remap );
//    }
//
//    return !usedIds.empty();
//}
