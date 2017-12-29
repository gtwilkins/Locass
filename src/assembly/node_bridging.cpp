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
#include <string.h>

bool Node::bridgeIsland( IslandVars &iv, NodeSetList &islandSets )
{
    NodeList islandEnds;
    NodeSet peSet;
    NodeIntMap qLimits, tLimits;
    
    Node::bridgeIslandGetPairs( iv, islandSets, peSet, qLimits, tLimits );
    
    if ( peSet.empty() )
    {
        NodeList peNodes;
        vector< vector<SeqNum> > peReadLists;
        for ( NodeSet isl : islandSets )
        {
            for ( Node* node : isl )
            {
                vector<SeqNum> thisPe;
                for ( const SeqNum &readId : iv.peReads )
                {
                    SeqNum pairId = params.getPairId( readId );
                    if ( node->reads_.find( readId ) != node->reads_.end() )
                    {
                        thisPe.push_back( readId );
                    }
                    if ( node->reads_.find( pairId ) != node->reads_.end() )
                    {
                        thisPe.push_back( pairId );
                    }
                }
                if ( !thisPe.empty() )
                {
                    peNodes.push_back( node );
                    peReadLists.push_back( thisPe );
                }
            }
        }
        return false;
    }
    
    for ( NodeSet &islandSet : islandSets )
    {
        Node::bridgeIslandGetIslandEnds( iv, peSet, islandEnds, qLimits, islandSet );
    }
    
    if ( !islandEnds.empty() )
    {
        NodeSet endSet( islandEnds.begin(), islandEnds.end() );
        if ( iv.limit > 25000 )
        {
            Node::trimIsland( iv, endSet );
        }
        NodeList mainEnds = Node::bridgeIslandGetMainEnds( iv, peSet, tLimits );
        
        if ( !mainEnds.empty() )
        {
            NodeList tNodes;
            for ( auto &t : tLimits )
            {
                if ( iv.ev.del.find( t.first ) == iv.ev.del.end() )
                {
                    tNodes.push_back( t.first );
                }
            }

            vector<bool> iAttempted( islandEnds.size(), false );
            return Node::bridgeIslandSet( iv, mainEnds, islandEnds, iAttempted, tNodes );
        }
    }
    
    return false;
}

//void Node::bridgeIslandDump( IslandVars &iv, NodeSetList &islandSets, Querier &bwt )
//{
//    NodeList islandEnds;
//    NodeSet peSet;
//    NodeIntMap qLimits, tLimits;
//    
//    Node::bridgeIslandGetPairs( iv, islandSets, peSet, qLimits, tLimits );
//    
//    if ( peSet.empty() )
//    {
//        NodeList peNodes;
//        vector< vector<SeqNum> > peReadLists;
//        for ( NodeSet isl : islandSets )
//        {
//            for ( Node* node : isl )
//            {
//                vector<SeqNum> thisPe;
//                for ( const SeqNum &readId : iv.peReads )
//                {
//                    SeqNum pairId = params.getPairId( readId );
//                    if ( node->reads_.find( readId ) != node->reads_.end() )
//                    {
//                        thisPe.push_back( readId );
//                    }
//                    if ( node->reads_.find( pairId ) != node->reads_.end() )
//                    {
//                        thisPe.push_back( pairId );
//                    }
//                }
//                if ( !thisPe.empty() )
//                {
//                    peNodes.push_back( node );
//                    peReadLists.push_back( thisPe );
//                }
//            }
//        }
//    }
//    
//    for ( NodeSet &islandSet : islandSets )
//    {
//        Node::bridgeIslandGetIslandEnds( iv, peSet, islandEnds, qLimits, islandSet );
//    }
//    
//    NodeList mainEnds = Node::bridgeIslandGetMainEnds( iv, peSet, tLimits );
//    NodeList tNodes;
//    for ( auto &t : tLimits )
//    {
//        if ( iv.ev.del.find( t.first ) == iv.ev.del.end() )
//        {
//            tNodes.push_back( t.first );
//        }
//    }
//    Node::bridgeIslandSetOffsets( iv, islandEnds, tNodes );
//    NodeSet dumpSet, endSet;
//    int32_t coords[3] = { std::numeric_limits<int32_t>::min(), std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::max() };
//    for ( bool drxn : { 0, 1 } )
//    {
//        for ( Node* node : ( drxn == iv.drxn ? islandEnds : mainEnds ) )
//        {
//            coords[drxn] = drxn ? min( coords[1], node->ends_[0] ) : max( coords[0], node->ends_[1] );
//            dumpSet.insert( node );
//            node->getDrxnNodes( dumpSet, 0 );
//            node->getDrxnNodes( dumpSet, 1 );
//            for ( Node* nxt : node->getNextNodes( drxn ) )
//            {
//                nxt->getDrxnNodes( endSet, !drxn );
//            }
//            endSet.insert( node );
//            node->getDrxnNodes( endSet, !drxn );
//            node->offsetForward( !drxn, false, true );
//            coords[2] = min( coords[2], node->ends_[0] );
//        }
//    }
//    
//    coords[0] -= 500;
//    coords[1] += 500;
//    
//    vector< pair<SeqNum, int32_t> > readIds;
//    for ( Node* node : dumpSet )
//    {
//        for ( bool drxn : {0, 1} )
//        {
//            if ( drxn ? coords[0] < node->ends_[1] : node->ends_[0] < coords[1] )
//            {
//                for ( ReadMark &mark : node->marks_[drxn] )
//                {
//                    if ( coords[0] <= mark.estimate && mark.estimate <= coords[1] )
////                    if ( coords[0] <= mark.estimate && mark.estimate <= coords[1] && params.isReadPe( mark.readId ) )
//                    {
//                        bool doAdd = true;
//                        for ( Node* a : dumpSet )
//                        {
//                            if ( a->reads_.find( mark.id ) != a->reads_.end() )
//                            {
//                                doAdd = false;
//                                break;
//                            }
//                        }
//                        if ( doAdd )
//                        {
//                            readIds.push_back( make_pair( mark.id, mark.estimate ) );
//                        }
//                    }
//                }
//            }
//        }
//    }
//    
//    ofstream dump2( "/home/glen/PythonProjects/BioJunkyard/data/Export/dump2" );
//    int x = 0;
//    for ( Node* node : endSet )
//    {
//        if ( node->drxn_ > 2 )
//        {
//            node->id_ = to_string( x );
//            x++;
//        }
//        node->exportNodeDump( dump2 );
//    }
//    dump2.close();
//    
//    sort( readIds.begin(), readIds.end(), []( pair<SeqNum, int32_t> &a, pair<SeqNum, int32_t> &b ) {
//        return a.second < b.second;
//    });
//    
//    coords[2] = readIds[0].second;
//    
//    
//    if ( !readIds.empty() )
//    {
//        ofstream fh( ( iv.drxn ? "/home/glen/PythonProjects/BioJunkyard/data/Export/align1.fa" : "/home/glen/PythonProjects/BioJunkyard/data/Export/align0.fa" ) );
//        for ( int i(0); i < mainEnds.size(); i++ )
//        {
//            fh << ">Main" << i << endl;
//            fh << mainEnds[i]->seq_ << endl;
//        }
//        
//        for ( Node* node : mainEnds )
//        {
//            node->bridgeIslandDump( bwt, fh, iv.drxn );
//        }
//        
//        for ( int i(0); i < mainEnds.size(); i++ )
//        {
//            fh << ">Island" << i << endl;
//            fh << islandEnds[i]->seq_ << endl;
//        }
//        
//        for ( Node* node : islandEnds )
//        {
//            node->bridgeIslandDump( bwt, fh, !iv.drxn );
//        }
//        
//        for ( pair<SeqNum, int32_t> readId : readIds )
//        {
//            bool doAdd = false;
//            int32_t offset = max( 0, readId.second - readIds[0].second );
//            string seq = bwt.getSequence( readId.first );
//            vector<string> seqs;
//            vector<int> offsets;
//            size_t it = seq.find( "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG" );
//            int hit = 0;
//            if ( it != seq.npos )
//            {
//                seqs.push_back( seq.substr( 0, it ) );
//                seqs.push_back( seq.substr( it + 38 ) );
//                offsets.push_back( 0 );
//                offsets.push_back( it );
//            }
//            else if ( params.isReadPe( readId.first ) )
//            {
//                int j = 0, len = 15;
//                while ( j < seq.length() )
//                {
//                    seqs.push_back( seq.substr( j, len ) );
//                    offsets.push_back( j );
//                    j += len;
//                }
//                seqs.push_back( seq.substr( seq.length() - len, len ) );
//                offsets.push_back( seq.length() - len );
//            }
//            
//            for ( Node* node : endSet )
//            {
//                for ( int k ( 0 ); k < seqs.size(); k++ )
//                {
//                    size_t it = node->seq_.find( seqs[k] );
//                    if ( it != node->seq_.npos && seqs[k].length() >= 15 )
//                    {
//                        offset = max( 0, node->ends_[0] + int32_t(it) - coords[2] - offsets[k] );
//                        doAdd = true;
//                        hit = k;
//                        break;
//                    }
//                }
//                if ( doAdd ) break;
//            }
//            
//            if ( doAdd || params.isReadPe( readId.first ) )
//            {
//                fh << ">" << readId.first << endl;
//                if ( params.isReadPe( readId.first ) )
//                {
//                    fh << string( offset, '-' ) << seq << endl;
//                }
//                else
//                {
//                    fh << string( offset, '-' ) << seqs[hit] << endl;
//                }
//            }
//        }
//        fh.close();
//    }
//    
//    int y = 0;
//}

//void Node::bridgeIslandDump( Querier &bwt, ofstream &fh, bool drxn )
//{
//    vector< pair<SeqNum, Coords> > reads;
//    for ( auto &read : reads_ )
//    {
//        reads.push_back( read );
//    }
//    sort( reads.begin(), reads.end(), [&drxn]( pair<SeqNum, Coords> &a, pair<SeqNum, Coords> &b ) {
//        return ( drxn ? a.second[1] > b.second[1] 
//                      : a.second[0] < b.second[0] );
//    });
//    
//    for ( int i ( 0 ); i < min( 3, (int)reads.size() ); i++ )
//    {
//        fh << ">" << reads[i].first << endl;
//        fh << bwt.getSequence( reads[i].first ) << endl;
//    }
//    
//    for ( Node* fwd : getDrxnNodes( drxn ) )
//    {
//        for ( auto &read : fwd->reads_ )
//        {
//            fh << ">" << read.first << endl;
//            fh << string( read.second[0] - fwd->ends_[0], '-' ) << bwt.getSequence( read.first ) << endl;
//        }
//    }
//}

void Node::bridgeIslandGetEnds( IslandVars &iv, NodeIntMap &limitMap, NodeList &endList, NodeSet &endSet, bool isIsland, bool drxn )
{
//    NodeSet foldable;
//    NodeList forkList = Node::bridgeIslandGetForks( iv, limitMap, endList, endSet, foldable, isIsland, drxn );
//    
//    NodeSet multiFork;
//    for ( Node* node : forkList )
//    {
//        int endCount = 0;
//        for ( Node* fwd : node->getDrxnNodes( drxn ) )
//        {
//            endCount += endSet.find( fwd ) != endSet.end();
//        }
//        if ( endCount > 1 )
//        {
//            multiFork.insert( node );
//        }
//    }
//    
//    for ( Node* node : endSet )
//    {
//        bool dontFold = node->reads_.size() > 2 || !node->edges_[drxn].empty();
//        NodeSet bckSet = node->getDrxnNodes( !drxn );
//        for ( Node* bck : bckSet )
//        {
//            dontFold = dontFold && multiFork.find( bck ) == multiFork.end();
//        }
//        if ( dontFold && find( endList.begin(), endList.end(), node ) == endList.end() )
//        {
//            endList.push_back( node );
//            for ( auto it = forkList.begin(); it != forkList.end(); )
//            {
//                if ( bckSet.find( *it ) != bckSet.end() )
//                {
//                    it = forkList.erase( it );
//                    continue;
//                }
//                it++;
//            }
//        }
//    }
//    
//    // Convert forks into new blank ends
//    for ( auto it = forkList.begin(); it != forkList.end(); )
//    {
//        vector<Edge> edges = (*it)->edges_[drxn];
//        sort( edges.begin(), edges.end(), []( Edge &a, Edge &b ){
//            return a.overlap > b.overlap;
//        });
//        
//        if ( edges.empty() )
//        {
//            endList.push_back( *it );
//            it = forkList.erase( it );
//            continue;
//        }
//        
//        while ( edges.size() > 1 )
//        {
//            bool didAdvance = false;
//            NodeSet prvSet;
//            for ( Node* prv : edges[0].node->getNextNodes( !drxn ) )
//            {
//                prvSet.insert( prv );
//            }
//            while ( edges.size() > 1 && !didAdvance )
//            {
//                for ( Node* prv : edges[1].node->getNextNodes( !drxn ) )
//                {
//                    if ( prvSet.find( prv ) == prvSet.end() )
//                    {
//                        edges.erase( edges.begin() );
//                        didAdvance = true;
//                        break;
//                    }
//                }
//                if ( !didAdvance )
//                {
//                    edges.erase( edges.begin() + 1 );
//                }
//            }
//        }
//        
//        Node* node = new Node( (*it)->ends_[drxn], edges[0].overlap, ( isIsland ? 3 + iv.drxn : drxn ), drxn );
//        ( isIsland ? iv.ev.island : iv.ev.nodes ).push_back( node );
//        endList.push_back( node );
//        limitMap[node] = node->ends_[drxn];
//        NodeSet edgedSet;
//        
//        for ( Edge &re : edges[0].node->edges_[!drxn] )
//        {
//            for ( Node* nxt : re.node->getNextNodes( drxn ) )
//            {
//                if ( nxt->reads_.size() <= 2 && nxt->edges_[drxn].empty() && nxt != node && !nxt->clones_ )
//                {
//                    foldable.insert( nxt );
//                }
//            }
//            re.node->addEdge( node, re.overlap, drxn, false );
//            edgedSet.insert( re.node );
//        }
//        
//        if ( !node->setBlank( iv, foldable, drxn ) )
//        {
//            node->dismantleNode();
//            iv.ev.del.insert( node );
//            endList.push_back( *it );
//            it = forkList.erase( it );
//            continue;
//        }
//        
//        limitMap[node] = node->ends_[drxn];
//        assert( !node->reads_.empty() );
//        
//        for ( auto it2 = it; it2 != forkList.end(); )
//        {
//            if ( edgedSet.find( *it2 ) != edgedSet.end() )
//            {
//                it2 = forkList.erase( it2 );
//                continue;
//            }
//            it2++;
//        }
//    }
//    
//    for ( auto it = endList.begin(); it != endList.end(); )
//    {
//        NodeSet endFwdSet = (*it)->getDrxnNodes( drxn );
//        bool doErase = false;
//        for ( auto it2 = endList.begin(); it2 != endList.end(); it2++ )
//        {
//            doErase = doErase || endFwdSet.find( *it2 ) != endFwdSet.end();
//        }
//        if ( doErase )
//        {
//            it = endList.erase( it );
//            continue;
//        }
//        it++;
//    }
}

NodeList Node::bridgeIslandGetForks( IslandVars &iv, NodeIntMap &limitMap, NodeSet &endSet, bool drxn )
{
    NodeSet forkSet, forkBckSet;
    NodeList endList;
    
    // Find potential forks and definite ends
    for ( Node* node : endSet )
    {
        int readCount = 1;
        int32_t limit = limitMap[node];
        for ( auto &read : node->reads_ )
        {
            readCount += ( drxn ? read.second[1] < limit : limit < read.second[0] );
        }
        
        if ( readCount <= 2 && node->edges_[!drxn].size() == 1 )
        {
            Node* fork = node;
            Node* prv = node->edges_[!drxn][0].node;
            fork = prv->edges_[drxn].size() > 1 ? prv : fork;
            if ( prv->edges_[!drxn].size() == 1 && prv->edges_[drxn].size() > 1 && readCount + prv->reads_.size() <= 2 )
            {
                fork = prv->edges_[drxn][0].node;
            }
            forkSet.insert( fork );
            fork->getDrxnNodes( forkBckSet, !drxn );
        }
        else
        {
            endList.push_back( node );
            forkBckSet.insert( node );
            node->getDrxnNodes( forkBckSet, !drxn );
        }
    }
    
    for ( Node* fork : forkSet )
    {
        if ( forkBckSet.find( fork ) == forkBckSet.end() )
        {
            endList.push_back( fork );
        }
    }
    
    for ( Node* node : endList )
    {
        for ( int i ( 0 ); i < node->edges_[drxn].size(); i++ )
        {
            if ( node->edges_[drxn][i].isLeap )
            {
                node->edges_[drxn][i].node->dismantleNode( iv.ev.del, drxn );
                i--;
            }
        }
    }
    
    return endList;
}

void Node::bridgeIslandGetPairs( IslandVars &iv, NodeSetList &islandSets, NodeSet &peSet, NodeIntMap &qLimits, NodeIntMap &tLimits )
{
    NodeSet tSet = iv.pathEnd->getDrxnNodes( iv.drxn, true, true );
    if ( abs( iv.limit - iv.pathEnd->ends_[!iv.drxn] ) < params.maxMpMax )
    {
        iv.pathEnd->getDrxnNodes( tSet, !iv.drxn, params.getFurthestMpDist( iv.limit, !iv.drxn ) );
    }
    NodeList tNodes( tSet.begin(), tSet.end() );
    
    for ( NodeSet islandSet : islandSets )
    {
        for ( Node* node : islandSet )
        {
            for ( ReadMark &mark : node->getMarksBase( iv.drxn ) )
            {
                for ( Node* t : tNodes )
                {
                    auto it = t->reads_.find( mark.id );
                    if ( it != t->reads_.end() )
                    {
                        SeqNum pairId = params.getPairId( mark.id );
                        auto it2 = node->reads_.find( pairId );
                        int32_t markCoord = it2->second[!iv.drxn];
                        auto r = qLimits.insert( make_pair( node, markCoord ) );
                        if ( !r.second )
                        {
                            r.first->second = ( iv.drxn ? min( r.first->second, markCoord ) 
                                                        : max( r.first->second, markCoord ) );
                        }
                        int32_t tLimit = iv.drxn ? max( it->second[1], t->validLimits_[3] )
                                                 : min( it->second[0], t->validLimits_[0] );
                        r = tLimits.insert( make_pair( t, tLimit ) );
                        if ( !r.second )
                        {
                            r.first->second = ( iv.drxn ? max( r.first->second, tLimit ) 
                                                        : min( r.first->second, tLimit ) );
                        }
                        t->pushValidLimits( it->second[iv.drxn], iv.drxn );
                        
                        if ( params.isReadPe( mark.id ) )
                        {
                            iv.peReads.insert( pairId );
                            peSet.insert( node );
                            peSet.insert( t );
                        }
                        else
                        {
                            iv.mpReads.insert( pairId );
                        }
                    }
                }
            }
        }
    }
}

void Node::bridgeIslandOffset( IslandVars &iv, NodeSet &islandSet, bool drxn )
{
    int32_t peOffSum = 0, mpOffSum = 0, peLimits[2], mpLimits[2];
    int32_t islandEndCoord = ends_[drxn], baseEndCoord = ends_[!drxn];
    Coords* baseCoords = NULL;
    int peCount = 0, mpCount = 0;
    
    for ( Node* node : islandSet )
    {
        NodeList tNodes;
        for ( auto &np : node->pairs_ )
        {
            if ( np.first->drxn_ <= 2 )
            {
                tNodes.push_back( np.first );
            }
        }
        
        for ( ReadMark &mark : node->getMarksBase( drxn ) )
        {
            for ( Node* t : tNodes )
            {
                auto hit = t->reads_.find( mark.id );
                if ( hit != t->reads_.end() )
                {
                    if ( node == this )
                    {
                        islandEndCoord = drxn ? min( islandEndCoord, mark.mark ) : max( islandEndCoord, mark.mark );
                    }
                    if ( t == this )
                    {
                        baseEndCoord = drxn ? max( baseEndCoord, hit->second[1] ) : min( baseEndCoord, hit->second[0] );
                        if ( hit->second[drxn] == baseEndCoord )
                        {
                            baseCoords = &hit->second;
                        }
                    }
                    
                    int32_t thisOff = drxn ? hit->second[0] - mark.estimate : hit->second[1] - mark.estimate;
                    if ( params.isReadPe( mark.id ) )
                    {
                        peLimits[0] = peCount > 0 ? min( peLimits[0], thisOff ) : thisOff;
                        peLimits[1] = peCount > 0 ? max( peLimits[1], thisOff ) : thisOff;
                        peOffSum += thisOff;
                        peCount++;
                    }
                    else
                    {
                        mpLimits[0] = mpCount > 0 ? min( mpLimits[0], thisOff ) : thisOff;
                        mpLimits[1] = mpCount > 0 ? max( mpLimits[1], thisOff ) : thisOff;
                        mpOffSum += thisOff;
                        mpCount++;
                    }
                    break;
                }
            }
        }
    }
    
    islandEndCoord = drxn ? min( ends_[0], islandEndCoord - params.readLen )
                          : max( ends_[1], islandEndCoord + params.readLen );
    
    if ( !peCount && !mpCount ) return;
    if ( peCount > 2 )
    {
        peOffSum -= ( peLimits[0] + peLimits[1] );
        peCount -= 2;
    }
    if ( mpCount > 2 )
    {
        mpOffSum -= ( mpLimits[0] + mpLimits[1] );
        mpCount -= 2;
    }
    
    int32_t off = peCount ? peOffSum / peCount : mpOffSum / mpCount;
    
    NodeSet propagated;
    offset( off );
    offsetIsland( propagated, drxn );
}

void Node::bridgeIslandConsolidateEnds( NodeIntMap &limitMap, NodeList &candidates, NodeSet &bckSet, bool drxn )
{
    bool didUpdate = true;
    while ( didUpdate )
    {
        didUpdate = false;
        for ( auto it = candidates.begin(); it != candidates.end(); )
        {
            if ( bckSet.find( *it ) != bckSet.end() )
            {
                it = candidates.erase( it );
                continue;
            }
            it++;
        }

        for ( int i ( 0 ); i + 1 < candidates.size(); )
        {
            bool iObsolete = false;
            NodeSet candiFwdSet = candidates[i]->getDrxnNodes( drxn );
            for ( int j ( i + 1 ); j < candidates.size(); )
            {
                bool jObsolete = false;
                NodeSet currSet = candidates[j]->getNextNodes( drxn );
                while ( !currSet.empty() )
                {
                    NodeSet nxtSet;
                    for ( Node* curr : currSet )
                    {
                        if ( candiFwdSet.find( curr ) != candiFwdSet.end() )
                        {
                            if ( find( candidates.begin(), candidates.end(), curr ) == candidates.end() )
                            {
                                candidates.push_back( curr );
                                curr->getDrxnNodes( bckSet, !drxn );
                                int32_t limit = curr->ends_[!drxn];
                                curr->getNextReadCoord( limit, drxn, drxn );
                                limitMap[curr] = limit;
                            }
                            iObsolete = true;
                            jObsolete = true;
                            didUpdate = true;
                        }
                        else
                        {
                            curr->getNextNodes( nxtSet, drxn );
                        }
                    }
                    currSet = nxtSet;
                }
                if ( jObsolete )
                {
                    candidates.erase( candidates.begin() + j );
                    continue;
                }
                j++;
            }
            if ( iObsolete )
            {
                candidates.erase( candidates.begin() + i );
                continue;
            }
            i++;
        }
    }
}

void Node::bridgeIslandGetIslandEnds( IslandVars &iv, NodeSet &peSet, NodeList &endList, NodeIntMap &limitMap, NodeSet &islandSet )
{
    NodeSet fwdSet, bckSet;
    NodeList candidates;
    
    for ( Node* node : islandSet )
    {
        if ( peSet.find( node ) != peSet.end() )
        {
            candidates.push_back( node );
            node->getDrxnNodes( fwdSet, iv.drxn );
        }
    }
    
    for ( Node* node : candidates )
    {
        if ( fwdSet.find( node ) == fwdSet.end() )
        {
            node->getDrxnNodes( bckSet, !iv.drxn );
        }
    }
    
    for ( Node* bck : bckSet )
    {
        if ( limitMap.find( bck ) != limitMap.end() )
        {
            candidates.push_back( bck );
            bck->getDrxnNodes( fwdSet, iv.drxn );
        }
    }
    
    if ( candidates.empty() ) return;
    
    Node::bridgeIslandConsolidateEnds( limitMap, candidates, fwdSet, !iv.drxn );
    NodeSet endSet( candidates.begin(), candidates.end() );
    
    assert( !endSet.empty() );
    
    NodeList thisEndList = Node::bridgeIslandGetForks( iv, limitMap, endSet, !iv.drxn );
    endList.insert( endList.end(), thisEndList.begin(), thisEndList.end() );
    
    assert( !endList.empty() );
}

NodeList Node::bridgeIslandGetMainEnds( IslandVars &iv, NodeSet &peSet, NodeIntMap &limitMap )
{
    NodeSet bckSet, fwdSet;
    NodeList candidates;
    
    for ( auto &t : limitMap )
    {
        if ( peSet.find( t.first ) != peSet.end() )
        {
            candidates.push_back( t.first );
            t.first->getDrxnNodes( bckSet, !iv.drxn );
            t.first->getDrxnNodes( fwdSet, iv.drxn );
        }
    }
    
    for ( Node* fwd : fwdSet )
    {
        if ( bckSet.find( fwd ) == bckSet.end() && !fwd->pairs_.empty()
                && find( candidates.begin(), candidates.end(), fwd ) == candidates.end() )
        {
            candidates.push_back( fwd );
            fwd->getDrxnNodes( bckSet, !iv.drxn );
            if ( limitMap.find( fwd ) == limitMap.end() )
            {
                limitMap[fwd] = fwd->validLimits_[iv.drxn*3];
            }
        }
    }
    
    assert( !candidates.empty() );
    
    Node::bridgeIslandConsolidateEnds( limitMap, candidates, bckSet, iv.drxn );
    NodeSet endSet = { candidates.begin(), candidates.end() };
    assert( !endSet.empty() );
    
    NodeList endList = Node::bridgeIslandGetForks( iv, limitMap, endSet, iv.drxn );
    
    return Node::bridgeIslandGetMainEndsSort( iv, endList, iv.drxn );
}

NodeList Node::bridgeIslandGetMainEndsSort( IslandVars &iv, NodeList &endList, bool drxn )
{
    NodeIntMap pairMap;
    NodeSet usedSet;
    
    for ( Node* node : endList )
    {
        for ( Node* bck : node->getDrxnNodesNotInSet( usedSet, !drxn, true ) )
        {
            pairMap[bck] = bck->getPairHitsTotal();
            usedSet.insert( bck );
        }
    }
    
    usedSet.clear();
    NodeList sortedEnds;
    
    while ( !endList.empty() )
    {
        int bestEnd = 0;
        int bestScore = 0;
        int32_t furthest = 0;
        NodeSet bestUsed;
        for ( int i ( 0 ); i < endList.size(); i++ )
        {
            int pairCount = 0;
            NodeSet currUsed = endList[i]->getDrxnNodesNotInSet( usedSet, !drxn, true );
            for ( Node* bck : currUsed )
            {
                auto it = pairMap.find( bck );
                pairCount += it != pairMap.end() ? it->second : 0;
            }
            int32_t farness = abs( endList[i]->ends_[drxn] );
            if ( i == bestEnd || pairCount > bestScore || ( pairCount == bestScore && farness > furthest ) )
            {
                bestEnd = i;
                bestScore = pairCount;
                furthest = farness;
                bestUsed = currUsed;
            }
        }
        sortedEnds.push_back( endList[bestEnd] );
        usedSet.insert( bestUsed.begin(), bestUsed.end() );
        endList.erase( endList.begin() + bestEnd );
    }
    
    return sortedEnds;
}

//bool Node::bridgeIslandPerfect( IslandVars &iv, NodeList &mainEnds, NodeList &islandEnds )
//{
//    int maxHit = 0;
//    int32_t maxCoords[2];
//    NodeSet nodeSets[2];
//    Node* nodes[2]{ NULL };
//    
//    for ( Node* isl : islandEnds )
//    {
//        nodeSets[iv.drxn].insert( isl );
//        isl->getDrxnNodes( nodeSets[iv.drxn], !iv.drxn );
//    }
//    for ( Node* mn : mainEnds )
//    {
//        nodeSets[!iv.drxn].insert( mn );
//        mn->getDrxnNodes( nodeSets[!iv.drxn], iv.drxn );
//    }
//    for ( int i : { 0, 1 } )
//    {
//        for ( Node* n1 : nodeSets[i] )
//        {
//            string seq = i ? n1->seq_.substr( 0, max( 15, maxHit + 1 ) )
//                           : n1->seq_.substr( n1->seq_.length() - max( 15, maxHit + 1 ) );
//            for ( Node* n2 : nodeSets[!i] )
//            {
//                size_t it = n2->seq_.find( seq );
//                if ( it != n2->seq_.npos )
//                {
//                    int j = seq.length();
//                    if ( i )
//                    {
//                        while ( it + j < n2->seq_.length() 
//                                && j < n1->seq_.length()
//                                && n2->seq_[it + j] == n1->seq_[j] )
//                        {
//                            j++;
//                        }
//                    }
//                    else
//                    {
//                        while ( j < n1->seq_.length() && it > 0 
//                                && n1->seq_[ n1->seq_.length() - j - 1 ] == n2->seq_[ it - 1 ] )
//                        {
//                            j++;
//                            it--;
//                        }
//                    }
//                    int32_t hitCoords[2] = { ( i ? n2->ends_[0] + (int32_t)it + j : n1->ends_[1] ) 
//                                           , ( i ? n1->ends_[0] : n2->ends_[0] + (int32_t)it ) };
//                    if ( j > maxHit && abs( hitCoords[1] - hitCoords[0] + j ) < j + 100 )
//                    {
//                        maxHit = j;
//                        nodes[i] = n1;
//                        nodes[!i] = n2;
//                        memcpy( &maxCoords, &hitCoords, 8 );
//                    }
//                }
//            }
//        }
//    }
//    if ( maxHit && abs( maxCoords[1] - maxCoords[0] + maxHit ) < maxHit + 100 )
//    {
//        if ( maxCoords[0] == nodes[0]->ends_[1] 
//                && maxCoords[1] == nodes[1]->ends_[0] )
//        {
//            Node::joinEnds( iv, nodes, maxHit );
//            return true;
//        }
//        else if ( maxCoords[0] == nodes[0]->ends_[1] )
//        {
//            return nodes[0]->bridgeIslandPerfect( iv, nodes[1], maxCoords[1], maxHit, 1 );
//        }
//        else if ( maxCoords[1] == nodes[1]->ends_[0] )
//        {
//            return nodes[1]->bridgeIslandPerfect( iv, nodes[0], maxCoords[0], maxHit, 0 );
//        }
//    }
//    return false;
//}

//bool Node::bridgeIslandPerfect( IslandVars &iv, Node* node, int32_t coord, int overlap, bool drxn )
//{
//    NodeList currNodes = { node };
//    vector<int32_t> currCoords = { coord }, currOverlaps = { overlap };
//    NodeSet usedSet;
//    Node* nodes[2] = { this, this };
//    
//    while ( !currNodes.empty() )
//    {
//        NodeList nxtNodes;
//        vector<int32_t> nxtCoords, nxtOverlaps;
//        NodeSet currFwd;
//        for ( Node* const &curr : currNodes )
//        {
//            curr->getDrxnNodes( currFwd, drxn );
//        }
//        for ( int i = 0; i < currNodes.size(); i++ )
//        {
//            if ( currOverlaps[i] <= 0 ) continue;
//            if ( usedSet.find( currNodes[i] ) != usedSet.end() ) continue;
//            if ( currFwd.find( currNodes[i] ) == currFwd.end() )
//            {
//                int32_t nxtCoord = currCoords[i];
//                if ( currNodes[i]->ends_[!drxn] == currCoords[i] )
//                {
//                    nodes[drxn] = currNodes[i];
//                    usedSet.insert( currNodes[i] );
//                    currNodes[i]->getDrxnNodes( usedSet, drxn );
//                    Node::joinEnds( iv, nodes, currOverlaps[i] );
//                }
//                else if ( currNodes[i]->getNextReadCoord( nxtCoord, !drxn, drxn ) )
//                {
//                    assert( false );
//                }
//                else
//                {
//                    for ( Edge &e : currNodes[i]->edges_[drxn] )
//                    {
//                        int32_t diff = abs( currNodes[i]->ends_[drxn] - currCoords[i] ) - e.overlap;
//                        nxtCoord = e.node->ends_[!drxn] + ( drxn ? -min( diff, 0 ) : min( diff, 0 ) );
//                        nxtNodes.push_back( e.node );
//                        nxtCoords.push_back( nxtCoord );
//                        nxtOverlaps.push_back( min( currOverlaps[i] - diff, currOverlaps[i] ) );
//                    }
//                }
//            }
//            else
//            {
//                nxtNodes.push_back( currNodes[i] );
//                nxtCoords.push_back( currCoords[i] );
//                nxtOverlaps.push_back( currOverlaps[i] );
//            }
//        }
//        
//        currNodes = nxtNodes;
//        currCoords = nxtCoords;
//        currOverlaps = nxtOverlaps;
//    }
//    
//    return !usedSet.empty();
//}

bool Node::bridgeIslandSet( IslandVars &iv, NodeList &mainEnds, NodeList &islandEnds, vector<bool> &iAttempted, NodeList &tNodes )
{
    if ( !Node::bridgeIslandSetOffsets( iv, islandEnds, tNodes ) )
    {
        return false;
    }
    
//    if ( Node::bridgeIslandPerfect( iv, mainEnds, islandEnds ) )
//    {
//        return true;
//    }
    
    NodeSet islandShared;
    NodeSetList mainSets;
    Node::bridgeIslandSetEndSets( iv, mainEnds, islandEnds, mainSets, islandShared );
    
    vector< vector<int> > setsHits;
    vector< vector<int32_t> > setsOverlaps;
    int iBestMain = 0;
    int hitsBestMain = std::numeric_limits<int>::min();
    int32_t overlapBestMain = std::numeric_limits<int32_t>::min();
    vector<int> hitsBestMains;
    int unattempted = 0;
    
    for ( int i ( 0 ); i < islandEnds.size(); i++ )
    {
        vector<int> setHits;
        vector<int32_t> setOverlaps;
        int bestHits = 0;
        
        for ( int j ( 0 ); j < mainEnds.size(); j++ )
        {
            int hits = 0;
            for ( Node* node : islandEnds[i]->getDrxnNodesNotInSet( islandShared, iv.drxn, true ) )
            {
                for ( auto &np : node->pairs_ )
                {
                    hits += mainSets[j].find( np.first ) != mainSets[j].end() ? np.second : 0;
                }
                
            }
            if ( j > 0 )
            {
                bestHits = max( hits, bestHits );
            }
            setHits.push_back( hits );
            
            int32_t overlap = iv.drxn ? mainEnds[j]->ends_[1] - islandEnds[i]->ends_[0] 
                                      : islandEnds[i]->ends_[1] - mainEnds[j]->ends_[0];
            setOverlaps.push_back( overlap );
        }
        setsHits.push_back( setHits );
        setsOverlaps.push_back( setOverlaps );
        
        bestHits = setHits[0] - bestHits;
        hitsBestMains.push_back( bestHits );
        
        if ( iAttempted[i] ) continue;
        unattempted++;
        if ( bestHits > hitsBestMain 
             || ( bestHits == hitsBestMain 
                  && abs( setOverlaps[0] - ( params.readLen / 2 ) ) < abs( overlapBestMain - ( params.readLen / 2 ) ) ) )
        {
            hitsBestMain = bestHits;
            overlapBestMain = setOverlaps[0];
            iBestMain = i;
        }
    }
    if ( iAttempted[iBestMain] ) return false;
    
    vector<bool>islandsEdged( islandEnds.size(), false );
    vector<bool>mainsEdged( mainEnds.size(), false );
    int islandsEdgedCount = 0, mainsEdgedCount = 0;
    
    NodeSet islandSet = islandEnds[iBestMain]->getConnectedNodes( true );
    
//    for ( int i ( iBestMain ); --i >= 0; )
//    {
//        int32_t thisOffset = abs( setsOverlaps[i][0] - ( params.readLen / 2 ) );
//        int32_t bestOffset = abs( setsOverlaps[iBestMain][0] - ( params.readLen / 2 ) );
//        bool doSwap;
//        if ( islandSet.find( islandEnds[i] ) != islandSet.end() )
//        {
//            doSwap = thisOffset < bestOffset && hitsBestMains[i] > 0;
//            for ( int j ( 1 ); j < mainEnds.size(); j++ )
//            {
//                int32_t altOffset = abs( setsOverlaps[i][0] - ( params.readLen / 2 ) ) - ( params.readLen / 2 );
//                doSwap = doSwap && !( setsHits[i][j] == setsHits[i][0] && altOffset < bestOffset );
//            }
//        }
//        else
//        {
//            doSwap = thisOffset < bestOffset && setsOverlaps[iBestMain][0] < 0;
//        }
//        
//        if ( doSwap && !iAttempted[i] )
//        {
//            iBestMain = i;
//            islandSet = islandEnds[iBestMain]->getConnectedNodes( true );
//        }
//    }
    
    if ( setsOverlaps[iBestMain][0] > params.readLen * 2.5 )
    {
        for ( int i( iBestMain ); i < islandEnds.size(); i++ )
        {
            int32_t thisOffset = abs( setsOverlaps[i][0] - ( params.readLen / 2 ) );
            int32_t bestOffset = abs( setsOverlaps[iBestMain][0] - ( params.readLen / 2 ) );
            if ( !iAttempted[i] 
                    && islandSet.find( islandEnds[i] ) != islandSet.end() 
                    && thisOffset < params.readLen * 2 
                    && thisOffset < bestOffset )
            {
                iBestMain = i;
                islandSet = islandEnds[iBestMain]->getConnectedNodes( true );
            }
        }
    }
    
    NodeSet islandMerged = islandEnds[iBestMain]->foldEdge( iv.ev, mainEnds[0], iv.drxn );
    iAttempted[iBestMain] = true;
    if ( !islandMerged.empty() )
    {
        NodeSet propagated;
        for ( Node* merge : islandMerged )
        {
            iv.merged[!iv.drxn].insert( merge );
            merge->offsetNode( iv.drxn );
            merge->offsetIsland( propagated, iv.drxn );
        }
    }
    else if ( unattempted > 1 )
    {
        return Node::bridgeIslandSet( iv, mainEnds, islandEnds, iAttempted, tNodes );
    }
    else
    {
        return false;
    }
    
    mainsEdged[0] = true;
    islandsEdged[iBestMain] = true;
    islandsEdgedCount = 1;
    mainsEdgedCount = 1;

    // Attempt to edge remaining main ends
    for ( int j( 1 ); j < mainEnds.size(); j++ )
    {
        int iBest = iBestMain;
        int iBestScore = 0;
        for ( int i( 0 ); i < islandEnds.size(); i++ )
        {
            if ( i != iBestMain )
            {
                int score = setsHits[i][j] - setsHits[i][0];
                if ( score > iBestScore || ( score == iBestScore && iBest == iBestMain ) )
                {
                    iBestScore = score;
                    iBest = i;
                }
            }
        }
        
        if ( iBest != iBestMain )
        {
            if ( islandSet.find( islandEnds[iBest] ) != islandSet.end() )
            {
                islandMerged = islandEnds[iBest]->foldEdge( iv.ev, mainEnds[j], iv.drxn );
                for ( Node* merge : islandMerged )
                {
                    iv.merged[!iv.drxn].insert( merge );
                    mainsEdgedCount++;
                    mainsEdged[j] = true;
                    islandsEdged[iBest] = true;
                }
            }
            break;
        }
//        if ( iBest == iBestMain )
//        {
//            Node* mainFolded = mainEnds[j]->foldEnd( iv.ev, mainEnds[0], iv.drxn );
//            if ( mainFolded )
//            {
//                for ( Node* nxt : mainFolded->getNextNodes( iv.drxn ) )
//                {
//                    if ( nxt->drxn_ > 2 )
//                    {
//                        iv.merged[!iv.drxn].insert( nxt );
//                    }
//                }
//            }
//        }
    }

//    for ( int i( 0 ); i < islandEnds.size(); i++ )
//    {
//        if ( !islandsEdged[i] && islandSet.find( islandEnds[i] ) != islandSet.end() )
//        {
//            Node* islandFolded = islandEnds[i]->foldEnd( iv.ev, islandEnds[iBestMain], !iv.drxn );
//            if ( islandFolded )
//            {
//                for ( Node* nxt : islandFolded->getNextNodes( !iv.drxn ) )
//                {
//                    if ( nxt->drxn_ <= 2 )
//                    {
//                        iv.merged[!iv.drxn].insert( nxt );
//                        break;
//                    }
//                }
//                break;
//            }
//        }
//    }
    
    return mainsEdgedCount > 0;
}

void Node::bridgeIslandSetEndSets( IslandVars &iv, NodeList &mainEnds, NodeList &islandEnds, NodeSetList &mainSets, NodeSet &islandShared )
{
    NodeSet mainBckSet = mainEnds[0]->getDrxnNodes( !iv.drxn );
    NodeSet mainShared;
    if ( mainEnds.size() > 2 )
    {
        mainEnds.erase( mainEnds.begin() + 2, mainEnds.end() );
    }
    if ( mainEnds.size() == 2 )
    {
        for ( Node* bck : mainEnds[1]->getDrxnNodes( !iv.drxn ) )
        {
            if ( mainBckSet.find( bck ) != mainBckSet.end() )
            {
                mainShared.insert( bck );
            }
        }
    }
    for ( Node* node : mainEnds )
    {
        mainSets.push_back( node->getDrxnNodesNotInSet( mainShared, !iv.drxn, true ) );
    }
    
    NodeSet islandFwd;
    for ( Node* node : islandEnds )
    {
        for ( Node* fwd : node->getDrxnNodes( iv.drxn, true ) )
        {
            if ( islandFwd.find( fwd ) != islandFwd.end() )
            {
                islandShared.insert( fwd );
            }
            else
            {
                islandFwd.insert( fwd );
            }
        }
    }
}

bool Node::bridgeIslandSetOffsets( IslandVars &iv, NodeList &islandEnds, NodeList &tNodes )
{
    for ( auto it = islandEnds.begin(); it != islandEnds.end(); )
    {
        (*it)->offsetForward( iv.drxn, false, true );
        int32_t offsetSum = 0;
        int32_t offsetLimits[2];
        int offsetCount = 0;
        bool pe = false, mp = false, reliable = false;
        
        for ( Node* fwd : (*it)->getDrxnNodes( iv.drxn, false, true ) )
        {
            for ( ReadMark &mark : fwd->getMarksBase( iv.drxn ) )
            {
                SeqNum pairId = params.getPairId( mark.id );
                for ( Node* t : tNodes )
                {
                    auto it = t->reads_.find( mark.id );
                    if ( it != t->reads_.end() )
                    {
                        if ( params.isReadPe( pairId ) )
                        {
                            int32_t thisOff = iv.drxn ? it->second[0] - mark.estimate : it->second[1] - mark.estimate;
                            offsetSum += thisOff;
                            offsetLimits[0] = offsetCount > 0 ? min( offsetLimits[0], thisOff ) : thisOff;
                            offsetLimits[1] = offsetCount > 0 ? max( offsetLimits[1], thisOff ) : thisOff;
                            offsetCount++;
                            pe = true;
                        }
                        else
                        {
                            mp = true;
                        }
                        reliable = reliable || t->isReliable( false );
                    }
                }
            }
        }
        
        if ( !pe || ( !mp && !reliable ) )
        {
            islandEnds.erase( it );
            continue;
        }
        
        assert( offsetCount );
        
        if ( offsetCount > 2 )
        {
            offsetSum -= ( offsetLimits[0] + offsetLimits[1] );
            offsetCount -= 2;
        }
        
        int32_t off = offsetSum / offsetCount;
        (*it)->offset( off );
        (*it)->offsetForward( iv.drxn, false, true );
        it++;
    }
    
    sort( islandEnds.begin(), islandEnds.end(), [&iv]( Node* &a, Node* &b ){
        return ( iv.drxn ? a->ends_[1] < b->ends_[1] : a->ends_[0] > b->ends_[0] );
    } );
    
    return !islandEnds.empty();
}

bool Node::mapBridge( Node* target, PathVars &pv, MapNode* mn )
{
    pv.bwt.mapSequence( mn->seq, mn->ids, mn->coords );
    mn->recoil();
    
    int32_t coords[2] = { mn->bridgeCoords[0][0], mn->bridgeCoords[1][0] };
    NodeList hitNodes[2];
    vector<int32_t> hitCoords[2][2];
    NodeSet cloneSet, reachacble;
    for ( int i : { 0, 1 } )
    {
        int32_t iCoords[2]{ coords[i], coords[i] };
        iCoords[i] = iCoords[!i] + ( i ? mn->bridgeOverlaps[i][0] : -mn->bridgeOverlaps[i][0] );
        mn->bridges[i][0]->overlapExtend( pv.nds[pv.drxn], iCoords, hitNodes[i], hitCoords[i], pv.drxn, i );
        reachacble = target->getDrxnNodes( pv.drxn, true, true );
        target->getDrxnNodes( reachacble, !pv.drxn, true );
        for ( int j = 0; j < hitNodes[i].size(); )
        {
            if ( hitNodes[i][j]->drxn_ == !pv.drxn )
            {
                assert( false );
                hitNodes[i].erase( hitNodes[i].begin() + j );
                hitCoords[i][0].erase( hitCoords[i][0].begin() + j );
                hitCoords[i][1].erase( hitCoords[i][1].begin() + j );
                if ( hitNodes[i].empty() ) return false;
            }
            else
            {
                cloneSet.insert( hitNodes[i][j] );
                hitNodes[i][j]->getDrxnNodes( cloneSet, i );
                j++;
            }
        }
    }
    
    NodeList nodes;
    NodeSet newSet;
    vector<int> overlaps;
    setBridge( pv, newSet, nodes, overlaps, mn, pv.drxn );
    if ( !nodes.empty() )
    {
        if ( !pv.drxn )
        {
            reverse( nodes.begin(), nodes.end() );
            reverse( overlaps.begin(), overlaps.end() );
        }
        
        for ( Node* &node : nodes )
        {
            if ( newSet.find( node ) == newSet.end() )
            {
                if ( cloneSet.find( node ) != cloneSet.end() )
                {
                    node = new Node( node );
                    newSet.insert( node );
                    assert( false );
                }
                else
                {
                    node->clearEdges( !pv.drxn );
                    reachacble.insert( node );
                    if ( node->drxn_ > 2 )
                    {
                        node->clearEdges( pv.drxn );
                        node->clearPairs();
                        node->drxn_ = pv.drxn ;
                        newSet.insert( node );
                        pv.nds[pv.drxn+3].erase( remove( pv.nds[pv.drxn+3].begin(), pv.nds[pv.drxn+3].end(), node ), pv.nds[pv.drxn+3].end() );
                        pv.nds[pv.drxn].push_back( node );
                    }
                }
            }
            else pv.nds[pv.drxn].push_back( node );
        }
        
        for ( int i = 0; i < hitNodes[!pv.drxn].size(); i++ )
        {
            hitNodes[!pv.drxn][i]->addEdge( nodes[0], hitCoords[!pv.drxn][1][i] - hitCoords[!pv.drxn][0][i], pv.drxn );
            for ( int j = 0; j < hitNodes[pv.drxn].size(); j++ )
            {
                hitNodes[!pv.drxn][i]->removeEdge( hitNodes[pv.drxn][j], pv.drxn );
                hitNodes[pv.drxn][j]->removeEdge( hitNodes[!pv.drxn][i], !pv.drxn );
            }
        }
        for ( int i = 0; i < overlaps.size(); i++ )
        {
            nodes[i]->addEdge( nodes[i+1], overlaps[i], pv.drxn );
        }
        for ( int i = 0; i < hitNodes[pv.drxn].size(); i++ )
        {
            NodeSet fwdSet = nodes.back()->getDrxnNodes( pv.drxn );
            if ( fwdSet.find( hitNodes[pv.drxn][i] ) == fwdSet.end() )
            {
                nodes.back()->addEdge( hitNodes[pv.drxn][i], hitCoords[pv.drxn][1][i] - hitCoords[pv.drxn][0][i], pv.drxn );
            }
        }
        for ( Node* node : nodes )
        {
            node->setValid();
        }
        pv.newSet.insert( nodes.begin(), nodes.end() );
        return true;
    }
    return false;
}

//bool Node::mapBridge( PathVars &pv, Node* target, MapNode* mn, int32_t* coords, bool drxn )
//{
//    pv.bwt.mapSequence( mn->seq, mn->ids, mn->coords );
//    mn->recoil();
//    
//    NodeList hitNodes[2];
//    vector<int32_t> hitCoords[2][2];
//    NodeSet cloneSet;
//    for ( int i : { 0, 1 } )
//    {
//        int32_t iCoords[2]{ coords[i], coords[i] };
//        iCoords[i] = iCoords[!i] + ( i ? mn->bridgeOverlaps[i][0] : -mn->bridgeOverlaps[i][0] );
//        mn->bridges[i][0]->overlapExtend( pv.nodes, iCoords, hitNodes[i], hitCoords[i], pv.drxn, i );
//        cloneSet = target->getDrxnNodes( drxn, true, true );
//        target->getDrxnNodes( cloneSet, !drxn, true );
//        for ( int j = 0; j < hitNodes[i].size(); )
//        {
//            if ( cloneSet.find( hitNodes[i][j] ) == cloneSet.end() )
//            {
//                hitNodes[i].erase( hitNodes[i].begin() + j );
//                hitCoords[i][0].erase( hitCoords[i][0].begin() + j );
//                hitCoords[i][1].erase( hitCoords[i][1].begin() + j );
//                if ( hitNodes[i].empty() ) return false;
//            }
//            else j++;
//        }
//    }
//    
//    NodeSet fwdSet( hitNodes[drxn].begin(), hitNodes[drxn].end() );
//    for ( Node* node : hitNodes[drxn] )
//    {
//        node->getDrxnNodes( fwdSet, drxn );
//    }
//    for ( Node* node : hitNodes[!drxn] )
//    {
//        assert( fwdSet.find( node ) == fwdSet.end() );
//    }
//    
//    NodeList nodes;
//    NodeSet newSet;
//    vector<int> overlaps;
//    setBridge( pv, newSet, nodes, overlaps, mn, pv.drxn );
//    if ( !nodes.empty() )
//    {
//        if ( !pv.drxn )
//        {
//            reverse( nodes.begin(), nodes.end() );
//            reverse( overlaps.begin(), overlaps.end() );
//        }
//        
//        for ( Node* &node : nodes )
//        {
//            if ( newSet.find( node ) == newSet.end() )
//            {
//                if ( cloneSet.find( node ) != cloneSet.end() )
//                {
//                    node = new Node( node );
//                    newSet.insert( node );
//                }
//                else
//                {
//                    node->clearEdges( !pv.drxn );
//                    cloneSet.insert( node );
//                    if ( node->drxn_ > 2 )
//                    {
//                        node->clearEdges( pv.drxn );
//                        node->clearPairs();
//                        node->drxn_ = pv.drxn ;
//                        newSet.insert( node );
//                        pv.islands.erase( remove( pv.islands.begin(), pv.islands.end(), node ), pv.islands.end() );
//                        pv.nodes.push_back( node );
//                    }
//                }
//            }
//            else pv.nodes.push_back( node );
//        }
//        
//        for ( int i = 0; i < hitNodes[!pv.drxn].size(); i++ )
//        {
//            hitNodes[!pv.drxn][i]->addEdge( nodes[0], hitCoords[!pv.drxn][1][i] - hitCoords[!pv.drxn][0][i], pv.drxn );
//            for ( int j = 0; j < hitNodes[pv.drxn].size(); j++ )
//            {
//                hitNodes[!pv.drxn][i]->removeEdge( hitNodes[pv.drxn][j], pv.drxn );
//                hitNodes[pv.drxn][j]->removeEdge( hitNodes[!pv.drxn][i], !pv.drxn );
//            }
//        }
//        for ( int i = 0; i < overlaps.size(); i++ )
//        {
//            nodes[i]->addEdge( nodes[i+1], overlaps[i], pv.drxn );
//        }
//        for ( int i = 0; i < hitNodes[pv.drxn].size(); i++ )
//        {
//            NodeSet fwdSet = nodes.back()->getDrxnNodes( pv.drxn );
//            if ( fwdSet.find( hitNodes[pv.drxn][i] ) == fwdSet.end() )
//            {
//                nodes.back()->addEdge( hitNodes[pv.drxn][i], hitCoords[pv.drxn][1][i] - hitCoords[pv.drxn][0][i], pv.drxn );
//            }
//        }
//        for ( Node* node : nodes )
//        {
//            node->setValid();
//        }
//        return true;
//    }
//    return false;
//}

void Node::setBridge( PathVars &pv, NodeSet &newSet, NodeList &nodes, vector<int> &overlaps, MapNode* mn, bool drxn )
{
    int i = 0;
    int j = 0;
    int32_t currLimits[2] = { 0, 0 }, prevLimits[2];
    while ( j < mn->ids.size() )
    {
        bool found = false;
        for ( int k = pv.drxn; k < 5; k += 3 )
        {
            for ( Node* n : pv.nds[k] )
            {
                auto it = n->reads_.find( mn->ids[j] );
                if ( it != n->reads_.end() && !it->second.redundant )
                {
                    if ( i != j )
                    {
                        if ( !nodes.empty() ) overlaps.push_back( prevLimits[1] - currLimits[0] );
                        nodes.push_back( new Node( mn, i, j - 1, drxn ) );
                        newSet.insert( nodes.back() );
                        prevLimits[0] = currLimits[0];
                        prevLimits[1] = currLimits[1];
                        currLimits[0] = mn->coords[0][j];
                        currLimits[1] = mn->coords[1][j];
                    }
                    int32_t splitCoords[2] = { it->second[0], it->second[1] };
                    while ( it != n->reads_.end() && ++j < mn->ids.size() && ( splitCoords[1] <= it->second[1] || it->second.redundant ) )
                    {
                        splitCoords[1] = it->second[1];
                        currLimits[0] = min( currLimits[0], mn->coords[0][j-1] );
                        currLimits[1] = max( currLimits[1], mn->coords[1][j-1] );
                        it = n->reads_.find( mn->ids[j] );
                    }
                    if ( !nodes.empty() ) overlaps.push_back( prevLimits[1] - currLimits[0] );
                    nodes.push_back( n->splitNodeDual( splitCoords, pv.nds[k], drxn ) );
                    prevLimits[0] = currLimits[0];
                    prevLimits[1] = currLimits[1];
                    if ( j < mn->ids.size() )
                    {
                        currLimits[0] = mn->coords[0][j];
                        currLimits[1] = mn->coords[1][j];
                    }
                    i = j;
                    found = true;
                    break;
                }
            }
            if ( found ) break;
        }
        if ( found ) continue;
        currLimits[0] = min( currLimits[0], mn->coords[0][j] );
        currLimits[1] = max( currLimits[1], mn->coords[1][j] );
        j++;
    }
    if ( i < j )
    {
        if ( !nodes.empty() ) overlaps.push_back( prevLimits[1] - currLimits[0] );
        nodes.push_back( new Node( mn, i, j - 1, drxn ) );
        newSet.insert( nodes.back() );
    }
}

