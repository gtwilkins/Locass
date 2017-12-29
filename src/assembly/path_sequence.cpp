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

#include "path_sequence.h"
#include <cassert>
#include <algorithm>

ReadEndMap::ReadEndMap( string seq, ReadMark &mark, Node* hit, int32_t* hitCoords, int32_t offset, bool mapDrxn )
: seq ( seq ), node( hit ), edge( NULL ), id( mark.id ), off( offset )
{
    coords[mapDrxn] = hitCoords[0];
    coords[mapDrxn+1] = hitCoords[1];
    coords[(!mapDrxn)*2] = mapDrxn ? coords[2] - seq.length(): coords[0] + seq.length();
    ol = mapDrxn ? coords[2] - coords[1] : coords[1] - coords[0];
    edgeOl = 0;
    drxn = mark.estimate < mark.mark;
    doMap = edged = false;
}

string ReadEndMap::getSeq( unordered_set<ReadId> &usedIds, bool drxn )
{
    vector<ReadEndMap*> path = { this };
    while ( path.back()->edge ) path.push_back( path.back()->edge );
    if ( !drxn ) reverse( path.begin(), path.end() );
    string s = path[0]->seq;
    usedIds.insert( id );
    for ( int i = 1; i < path.size(); i++ )
    {
        usedIds.insert( path[i]->id );
        s += path[i]->seq.substr( path[i - drxn]->edgeOl );
    }
    
    return s;
}

void ReadEndMap::move( bool mapDrxn )
{
    if ( !edge ) return;
    node = edge->node;
    coords[1] = edge->coords[1];
    int32_t offset = mapDrxn ? edge->coords[0] + edgeOl - coords[2]
                             : edge->coords[2] - edgeOl - coords[0];
    if ( !offset ) return;
    for ( int i : { 0, 2 } ) coords[i] += offset;
    off += drxn ? -offset : offset;
}

void ReadEndMap::setEdge( vector<ReadEndMap*> &reads, bool drxn )
{
    for ( ReadEndMap* read : reads )
    {
        if ( read == this || read->seq == seq || edgeOl == seq.length() ) continue;
        int readOl = mapSeqOverlap( ( drxn ? seq : read->seq ), ( drxn ? read->seq : seq ), max( ol + 1, edgeOl ) );
        if ( readOl )
        {
            if ( readOl == ol && edge->id > read->id ) continue;
            edge = read;
            edgeOl = readOl;
        }
    }
}

SeqPath::SeqPath(NodeList& path, Node* anchor, int32_t offset)
: nodes( path )
{
    seq = nodes[0]->seq_;
    nodeCoords[0].push_back( 0 );
    nodeCoords[1].push_back( nodes[0]->seq_.length() );
    ends[0] = nodes[0]->ends_[0] - ( nodes[0] == anchor ? offset : 0 );
    for ( int i = 1; i < nodes.size(); i++ )
    {
        int ol = nodes[i-1]->getOverlap( nodes[i], 1 );
        nodeCoords[0].push_back( seq.length() - ol );
        if ( ol < 0 )
        {
            seq += string( -ol, 'N' );
            ol = 0;
        }
        seq += nodes[i]->seq_.substr( ol );
        nodeCoords[1].push_back( seq.length() );
        if ( nodes[i] == anchor )
        {
            ends[0] = nodes[i]->ends_[1] -seq.length() - offset;
        }
    }
    ends[1] = ends[0] + seq.length();
    assert( nodeCoords[1].back() == seq.length() );
}

SeqPathReassemble::SeqPathReassemble( NodeList &path, Node* anchor, int32_t offset )
: SeqPath( path, anchor, offset )
{}

SeqPathReassemble::~SeqPathReassemble()
{
    for ( int i : { 0, 1 } )
    {
        for ( int j = 0; j < reads[i].size(); j++ )
        {
            delete reads[i][j];
        }
    }
}

bool SeqPathReassemble::doMap( PathVars &pv )
{
    NodeList edgeNodes;
    for ( int i = 0; i < nodes.size() - 1; i++ )
    {
        edgeNodes.push_back( new Node() );
        edgeNodes.back()->ends_[0] = ends[1];
        edgeNodes.back()->ends_[1] = ends[0];
        edgeNodes.back()->drxn_ = nodes[i]->drxn_ > 0;
    }
    
    NodeSet delSet, added;
    for ( int i : { 0, 1 } )
    {
        for ( ReadEndMap* read : reads[i] )
        {
            bool redundant = read->ol < read->seq.length();
            if ( !read->doMap ) continue;
            for ( Node* node : pv.nds[pv.drxn+3] )
            {
                if ( node->reads_.find( read->id ) == node->reads_.end() ) continue;
                node->clearEdges( 0 );
                node->clearEdges( 1 );
                node->dismantleNode();
                delSet.insert( node );
            }
            
            for ( int j = 0; j < nodes.size(); j++ )
            {
                int32_t readCoords[2] = { read->coords[0] + ends[0], read->coords[2] + ends[0] };
                int offset = nodes[j]->ends_[0]- ends[0] - nodeCoords[0][j];
                if ( nodeCoords[0][j] <= read->coords[0] && read->coords[2] <= nodeCoords[1][j] )
                {
                    nodes[j]->addRead( read->id, readCoords[0] + offset, readCoords[1] + offset, redundant );
                    added.insert( nodes[j] );
                    break;
                }
                if ( read->coords[2] <= nodeCoords[1][j] )
                {
                    if ( nodeCoords[0][j] <= read->coords[0] )
                    {
                        nodes[j]->addRead( read->id, max( nodes[j]->ends_[0], readCoords[0] )
                                                   , min( nodes[j]->ends_[1], readCoords[1] ), redundant );
                        added.insert( nodes[j] );
                    }
                    else if ( j > 0 )
                    {
                        edgeNodes[j-1]->addRead( read->id, readCoords[0], readCoords[1], redundant );
                        edgeNodes[j-1]->ends_[0] = min( edgeNodes[j-1]->ends_[0], readCoords[0] );
                        edgeNodes[j-1]->ends_[1] = max( edgeNodes[j-1]->ends_[1], readCoords[1] );
                    }
                    break;
                }
            }
        }
    }
    
    for ( Node* del : delSet )
    {
        delete del;
        pv.nds[pv.drxn].erase( remove( pv.nds[pv.drxn].begin(), pv.nds[pv.drxn].end(), del ), pv.nds[pv.drxn].end() );
        pv.nds[pv.drxn+3].erase( remove( pv.nds[pv.drxn+3].begin(), pv.nds[pv.drxn+3].end(), del ), pv.nds[pv.drxn+3].end() );
    }
    
    bool didMap = false;
    for ( int i = 0; i < edgeNodes.size(); i++ )
    {
        if ( edgeNodes[i]->reads_.empty() ) delete edgeNodes[i];
        else
        {
            pv.nds[pv.drxn].push_back( edgeNodes[i] );
            nodes[i]->removeEdge( nodes[i+1], 1 );
            nodes[i+1]->removeEdge( nodes[i], 0 );
            int overlaps[2] = { ends[0] + nodeCoords[1][i] - edgeNodes[i]->ends_[0]
                              , edgeNodes[i]->ends_[1] - nodeCoords[0][i+1] - ends[0] };
            edgeNodes[i]->seq_ = seq.substr( edgeNodes[i]->ends_[0] - ends[0], edgeNodes[i]->ends_[1] - edgeNodes[i]->ends_[0] );
            nodes[i + !pv.drxn]->addEdge( edgeNodes[i], overlaps[!pv.drxn], pv.drxn );
            edgeNodes[i]->addEdge( nodes[i + pv.drxn], overlaps[pv.drxn], pv.drxn );
            edgeNodes[i]->setCoverage();
            edgeNodes[i]->setValid();
            didMap = true;
        }
    }
    
    for ( Node* node : added )
    {
        node->setCoverage();
        node->updatePairs();
        didMap = true;
    }
    
    return didMap;
}

int32_t SeqPathReassemble::getCoord( ReadEndMap* read )
{
    int i = 0;
    while ( nodes[i] != read->node ) i++;
    assert( i < nodes.size() );
    
    return nodes[i]->ends_[0] + read->coords[1] - nodeCoords[0][i];
}

string SeqPathReassemble::getExtraSeq( ReadEndMap* read, int &ol, bool drxn )
{
//    int32_t endCoord;
//    for ( int i = 0; i < nodes.size(); i++ )
//    {
//        if ( nodes[i] == read->node )
//        {
//            endCoord 
//        }
//    }
    ol = max( ol, read->ol );
    if ( drxn ) ol = min( ol, (int)seq.length() - read->coords[1] );
    else ol = min( ol, read->coords[1] );
    int i = drxn ? read->coords[2] : read->coords[1] - ol;
    int j = drxn ? read->coords[1] + ol : read->coords[0];
    assert( i <= j );
    string rSeq = seq.substr( i, j - i );
    if ( drxn )
    {
        size_t it = rSeq.find( 'N' );
        if ( it != rSeq.npos )
        {
            ol -= rSeq.length() - it;
            rSeq = rSeq.substr( 0, it );
        }
    }
    else
    {
        size_t it = rSeq.find_last_of( 'N' );
        if ( it != rSeq.npos )
        {
            ol -= 1 + it;
            rSeq = rSeq.substr( it + 1 );
        }
    }
    
    assert( rSeq.length() < 85 );
    return rSeq;
}

void SeqPathReassemble::map( string &q, ReadMark &mark, int minOls[2] )
{
    int32_t coords[2];
    int32_t estimates[2]{ mark.estimate, mark.estimate };
    bool markDrxn = mark.estimate < mark.mark;
    estimates[markDrxn] = markDrxn ? estimates[0] + q.length() : estimates[1] - q.length();
    
    for ( bool i : { 0, 1 } )
    {
        if ( usedIds[i].find( mark.id ) != usedIds[i].end() ) continue;
        if ( mapSeqEnd( q, seq, minOls[i], coords, i ) )
        {
            int32_t mapped[2];
            mapped[i] = coords[i] + ends[0];
            mapped[!i] = i ? mapped[1] - q.length() : mapped[0] + q.length();
            if ( mapped[0] < ends[0] || ends[1] < mapped[1] ) continue;
            usedIds[i].insert( mark.id );
            
            int j = i ? 0 : nodes.size()-1;
            while ( i ? ( j+1 < nodes.size() && ( nodeCoords[1][j] < coords[1] || nodeCoords[0][j+1] <= coords[0] ) )
                      : ( j > 0 && ( coords[0] < nodeCoords[0][j] || coords[1] <= nodeCoords[1][j-1] ) ) ) j += ( i ? 1 : -1 );
            Node* node = nodes[j];
            assert( j >= 0 && j < nodes.size() );
            int32_t offset = markDrxn ? estimates[0] - mapped[0] : mapped[1] - estimates[1];
            reads[i].push_back( new ReadEndMap( q, mark, node, coords, offset, i ) );
            assert( reads[i].back()->ol > 0 );
            complete = false;
        }
    }
}

void SeqPathReassemble::remap( string &q, ReadMark &mark )
{
    int32_t estimates[2]{ mark.estimate, mark.estimate };
    bool markDrxn = mark.estimate < mark.mark;
    estimates[markDrxn] = markDrxn ? estimates[0] + q.length() : estimates[1] - q.length();
    
    for ( bool i : { 0, 1 } )
    {
        if ( usedIds[i].find( mark.id ) != usedIds[i].end() ) continue;
        {
            int minOl = params.readLen / 3;
            int32_t coords[2]{0};
            int bestOl = 0;
            Node* node = NULL;
            for ( ReadEndMap* const &read : reads[i] )
            {
                int ol = mapSeqOverlap( ( i ? q : read->seq ), ( i ? read->seq : q ), minOl );
                if ( ol )
                {
                    coords[0] = i ? read->coords[1] : read->coords[2] - ol;
                    coords[1] = i ? read->coords[0] + ol : read->coords[1];
                    minOl = ol + 1;
                    bestOl = ol;
                    node = read->node;
                }
            }
            
            if ( bestOl )
            {
                int32_t mapped[2] = { coords[0] + ends[0], coords[1] + ends[0] };
                mapped[!i] = i ? mapped[1] - q.length() : mapped[0] + q.length();
                int32_t offset = markDrxn ? estimates[1] - mapped[1] : mapped[0] - estimates[0];
                reads[i].push_back( new ReadEndMap( q, mark, node, coords, offset, i ) );
                usedIds[i].insert( mark.id );
                complete = false;
            }
        }
    }
}

void SeqPathReassemble::setBridges( vector<SeqPathReassemble*> seqs, MapNode* mn, unordered_set<ReadId> &usedIds, float &bestScore, bool drxn, int fromDrxn )
{
    for ( vector<ReadEndMap*> &path : paths[drxn] )
    {
        int iMax = path.size()-2;
        for ( SeqPathReassemble* s : seqs )
        {
            for ( ReadEndMap* read : s->reads[!drxn] )
            {
                if ( read->edge ) continue;
                int i = 0;
                while ( i < iMax && path[i]->id != read->id ) i++;
                if ( i >= iMax ) continue;
                ReadEndMap* hitRead = path[i];
                int ols[3];
                ols[drxn] = path.back()->ol;
                ols[!drxn] = read->ol;
                ols[2] = min( ols[0], ols[1] );
                int len = abs( path[i]->coords[!drxn*2] - path.back()->coords[drxn*2] ) - path.back()->ol - read->ol;
                int dist = drxn ? ( path.back()->coords[1] + ends[0] ) - ( read->coords[1] + s->ends[0] )
                                : ( read->coords[1] + s->ends[0] ) - ( path.back()->coords[1] + ends[0] );
                int totalOl = ols[0] + ols[1];
                int drxns[2] = { 0 };
                int offsets[2] = { 0 };
                while ( i < path.size() - 1 )
                {
                    totalOl += path[i]->edgeOl;
                    ols[2] = min( ols[2], path[i]->edgeOl );
                    drxns[path[i]->drxn]++;
                    offsets[path[i]->drxn] += path[i]->off;
                    i++;
                }
                drxns[path[i]->drxn]++;
                offsets[path[i]->drxn] += path[i]->off;
                offsets[0] /= max( drxns[0], 1 );
                offsets[1] /= max( drxns[1], 1 );
                
                // If offsets are positive, then current is too long and a short cut is implicated
                // If offsets are negative, then an expanded bridge is implicated
                
                float score = totalOl * ( float( ols[0] + ols[1] + ols[2] ) / float( 3 * params.readLen ) );
                
                // Penalize disagreement between length expansion and offset
                score -= abs( offsets[drxn] ) * drxns[drxn];
                score -= abs( len - dist + offsets[!drxn] ) * drxns[!drxn];
                
                // Bonus or penalty for expansion depending on support
                float bonus = totalOl * min( (float) 1, float( abs( len - dist ) ) / (float)params.readLen );
                float olMod = min( (float)2, ( (float)totalOl / float(params.readLen) ) - 3 );
                for ( int i : { 0, 1, 2 } ) olMod += min( (float)1, ( float(ols[i]*4) / (float)params.readLen ) - 2 );
                float ratioMod = 1 - abs( float( len - dist + offsets[!drxn] ) 
                               / float( 50 + abs( float(len-dist) / (float)10 ) ) );
                float mod = min( olMod, ratioMod );
                mod = min( mod, sqrt( mod ) );
                if ( abs( bonus * mod ) > abs( score ) )
                {
                    int x = 0;
                }
                score += bonus * mod;
                
                if ( ( fromDrxn == 2 || drxns[fromDrxn] > 2 ) && score > bestScore )
                {
                    ReadEndMap* olReads[2] = { drxn ? read : path.back(), drxn ? path.back() : read };
                    SeqPathReassemble* hitSeqs[2] = { drxn ? s : this, drxn ? this : s };
                    int seqOls[2]{ params.readLen - 1, params.readLen - 1 };
                    usedIds.clear();
                    mn->seq = hitSeqs[0]->getExtraSeq( olReads[0], seqOls[0], 0 );
                    mn->seq += hitRead->getSeq( usedIds, drxn );
                    mn->seq += hitSeqs[1]->getExtraSeq( olReads[1], seqOls[1], 1 );
                    for ( int i : { 0, 1 } )
                    {
                        mn->bridges[i].clear();
                        mn->bridgeOverlaps[i].clear();
                        mn->bridgeCoords[i].clear();
                        mn->bridges[i].push_back( olReads[i]->node );
                        mn->bridgeOverlaps[i].push_back( seqOls[i] );
                        mn->bridgeCoords[i].push_back( hitSeqs[i]->getCoord( olReads[i] ) );
                    }
                    bestScore = score;
                }
            }
        }
    }
}

void SeqPathReassemble::setEdges()
{
    for ( bool drxn : { 0, 1 } )
    {
        // Set redundant
        unordered_set<ReadEndMap*> edged;
        for ( ReadEndMap* r1 : reads[drxn] )
        {
            for ( ReadEndMap* r2 : reads[drxn] )
            {
                if ( r1 == r2 ) continue;
                if ( r1->seq == r2->seq && r1->id < r2->id )
                {
                    if ( r2->edge && r2->edge->id > r1->id ) continue;
                    r2->edge = r1;
                    r2->edgeOl = r2->seq.length();
                }
            }
        }
        
        // Set edges
        for ( ReadEndMap* const &read : reads[drxn] )
        {
            read->setEdge( reads[drxn], drxn );
            if ( read->edge ) edged.insert( read->edge );
        }
        
        // Trim possible loops
        for ( ReadEndMap* read : reads[drxn] )
        {
            if ( edged.find( read ) == edged.end() || !read->edge ) continue;
            bool loop = false;
            ReadEndMap* weak = read;
            ReadEndMap* curr = read->edge;
            unordered_set<ReadEndMap*> pathedReads;
            while ( curr->edge && curr != read && pathedReads.find( curr ) == pathedReads.end() )
            {
                if ( curr->edgeOl < weak->edgeOl ) weak = curr;
                if ( curr->edge == read ) loop = true;
                pathedReads.insert( curr );
                curr = curr->edge;
            }
            assert( pathedReads.find( curr ) == pathedReads.end() );
            if ( loop )
            {
                weak->edge = NULL;
                edged.erase( weak );
            }
        }
        
        // Set path ends
        for ( ReadEndMap* r1 : reads[drxn] )
        {
            r1->edged = edged.find( r1 ) != edged.end();
            if ( edged.find( r1 ) == edged.end() && r1->edge )
            {
                vector<ReadEndMap*> path = { r1 };
                while ( path.back()->edge ) path.push_back( path.back()->edge );
                for ( int i = path.size()-1; --i >= 0; ) path[i]->move( drxn );
                if ( path.size() > 2 ) paths[drxn].push_back( path );
            }
        }
    }
    
    sortReads();
}

void SeqPathReassemble::setHalves( MapNode* mn, unordered_set<ReadId> &usedIds, int &bestScore, bool drxn )
{
    for ( vector<ReadEndMap*> &path : paths[drxn] )
    {
        int thisScore = path.back()->ol;
        for ( int i = 0; i < path.size()-1; i++ )
        {
            thisScore += path[i]->edgeOl;
        }
        
        if ( thisScore > bestScore )
        {
            int seqOl = params.readLen - 1;
            usedIds.clear();
            mn->seq = drxn ? "" : getExtraSeq( path.back(), seqOl, 0 );
            mn->seq += path[0]->getSeq( usedIds, drxn );
            mn->seq += drxn ? getExtraSeq( path.back(), seqOl, 1 ) : "";
            mn->bridges[drxn].clear();
            mn->bridgeOverlaps[drxn].clear();
            mn->bridgeCoords[drxn].clear();
            mn->bridges[drxn].push_back( path.back()->node );
            mn->bridgeOverlaps[drxn].push_back( seqOl );
            mn->bridgeCoords[drxn].push_back( getCoord( path.back() ) );
            bestScore = thisScore;
        }
    }
}

void SeqPathReassemble::sortReads()
{
    // Sorted left to right then by extension length
    sort( reads[0].begin(), reads[0].end(), []( ReadEndMap* const &a, ReadEndMap* const &b ){
        return a->coords[1] == b->coords[1] ? a->coords[2] > b->coords[2]
                                          : a->coords[1] < b->coords[1];
    } );
    sort( reads[1].begin(), reads[1].end(), []( ReadEndMap* const &a, ReadEndMap* const &b ){
        return a->coords[1] == b->coords[1] ? a->coords[0] < b->coords[0]
                                          : a->coords[1] < b->coords[1];
    } );
}

SeqPathMerge::SeqPathMerge( NodeList &path, Node* anchor, int32_t offset )
: SeqPath( path, anchor, offset ), hitSeq( NULL )
{
    
}

bool SeqPathMerge::doMerge( PathVars &pv, NodeSet &delSet, bool drxn )
{
    if ( ol >= params.readLen )
    {
        int diff = ol + 1 - params.readLen;
        selfCoord += drxn ? -diff : diff;
        ol -= diff;
    }
    
    Node* mergeNodes[2] = { NULL, NULL };
    int32_t splitCoords[2];
    splitCoords[drxn] = hitCoord;
    splitCoords[!drxn] = selfCoord;
    mergeNodes[drxn] = hitSeq->getSplitCoord( splitCoords[drxn], drxn );
    mergeNodes[!drxn] = getSplitCoord( splitCoords[!drxn], !drxn );
    if ( mergeNodes[0] && mergeNodes[1] )
    {
        if ( splitCoords[!drxn] != mergeNodes[!drxn]->ends_[drxn] )
        {
            int32_t splitCoord = mergeNodes[!drxn]->ends_[drxn];
            for ( auto &read : mergeNodes[!drxn]->reads_ )
            {
                if ( drxn ? splitCoords[!drxn] <= read.second[1] && read.second[0] < splitCoord
                          : read.second[0] <= splitCoords[!drxn] && splitCoord < read.second[1] )
                {
                    splitCoord = read.second[!drxn];
                }
            }
            
            if ( splitCoord == mergeNodes[!drxn]->ends_[!drxn] )
            {
                int diff = abs( mergeNodes[!drxn]->ends_[!drxn] - splitCoords[!drxn] );
                bool didFind = false;
                for ( Edge &e : mergeNodes[!drxn]->edges_[!drxn] )
                {
                    if ( find( nodes.begin(), nodes.end(), e.node ) == nodes.end() ) continue;
                    if ( e.overlap > diff );
                    mergeNodes[!drxn] = e.node;
                    splitCoords[!drxn] = e.node->ends_[drxn];
                    ol -= ( diff - e.overlap );
                    didFind = true;
                }
                if ( !didFind ) return false;
            }
            else if ( splitCoord != mergeNodes[!drxn]->ends_[drxn] )
            {
                mergeNodes[!drxn]->splitNode( splitCoord, pv.nds[pv.drxn], drxn, drxn );
                for ( Node* nxt : mergeNodes[!drxn]->getNextNodes( drxn ) )
                {
                    nxt->dismantleNode( delSet, drxn );
                }
                ol -= abs( splitCoords[!drxn] - mergeNodes[!drxn]->ends_[drxn] );
            }
            else
            {
                assert( false );
                return false;
            }
        }
        
        
        for ( Node* nxt : mergeNodes[!drxn]->getNextNodes( drxn ) )
        {
            nxt->dismantleNode( delSet, drxn );
        }
        
        int32_t splitCoord = splitCoords[drxn];
        if ( !mergeNodes[drxn]->getNextReadCoord( splitCoord, !drxn, drxn ) )
        {
            int endDist = abs( mergeNodes[drxn]->ends_[drxn] - splitCoords[drxn] );
            bool didMerge = false;
            for ( Edge &e : mergeNodes[drxn]->edges_[drxn] )
            {
                int diff = endDist - e.overlap;
                if ( diff < 0 ) continue;
                mergeNodes[!drxn]->addEdge( e.node, ol - diff, drxn );
                didMerge = true;
            }
            if ( !didMerge ) return false;
        }
        else
        {
            ol -= abs( splitCoord - splitCoords[drxn] );
            mergeNodes[drxn] = mergeNodes[drxn]->splitNode( splitCoord, pv.nds[pv.drxn], drxn, drxn );
            mergeNodes[!drxn]->addEdge( mergeNodes[drxn], ol, drxn );
        }
        
        bool doValidate = drxn;
        for ( Node* node : nodes )
        {
            if ( doValidate || node == mergeNodes[!drxn] ) node->setValid();
            if ( node == mergeNodes[!drxn] ) doValidate = !doValidate;
        }
        return true;
    }
    
    assert( false );
    return false;
}

Node* SeqPathMerge::getSplitCoord( int32_t &coord, bool drxn )
{
    Node* node = NULL;
    for ( int i = 0; i < nodes.size(); i++ )
    {
        if ( drxn ? nodeCoords[0][i] <= coord && ( i + 1 == nodes.size() || coord < nodeCoords[0][i+1] )
                  : coord <= nodeCoords[1][i] )
        {
            node = nodes[i];
            coord = node->ends_[0] + coord - nodeCoords[0][i];
//            int32_t splitCoord = coord;
//            if ( !node->getNextReadCoord( splitCoord, !drxn, drxn ) )
//            {
//                int endDist = node->ends_[]
//                assert( false );
//            }
//            diff = abs( splitCoord - coord );
//            coord = splitCoord;
            return node;
        }
    }
    
    return node;
}

void SeqPathMerge::merge( vector<SeqPathMerge*> seqs, bool endOnly, bool drxn )
{
    int len = min( params.readLen, (int)seq.length() - params.readLen );
    if ( endOnly ) len = drxn ? min( len, nodeCoords[1][1] - nodeCoords[1][0] )
                              : min( len, nodeCoords[0][1] - nodeCoords[0][0] );
    int attempts = 1 + len / 8;
    for ( int i = 0; i < attempts; i++ )
    {
        int32_t endCoord = drxn ? seq.length() - ( i * 8 ) : i * 8;
        int32_t thisCoords[2];
        string q = drxn ? seq.substr( 0, endCoord ) : seq.substr( endCoord );
        for ( SeqPathMerge* s : seqs )
        {
            if ( mapSeqEnd( q, s->seq, 16, thisCoords, drxn ) )
            {
                while ( drxn ? ( endCoord + 1 < seq.length() && thisCoords[1] < s->seq.length() 
                                    && seq[endCoord+1] == s->seq[thisCoords[1]] )
                             : ( endCoord > 0 && thisCoords[0] > 0 
                                    && seq[endCoord-1] == s->seq[thisCoords[0]] ) )
                {
                    endCoord += drxn ? 1 : -1;
                    thisCoords[drxn] += drxn ? 1 : -1;
                }
                int thisOl = thisCoords[1] - thisCoords[0];
                int offset = endCoord + ends[0] - ( thisCoords[drxn] + s->ends[0] );
                int thisScore = thisOl - abs( offset );
                if ( thisScore > 0 && ( !hitSeq || thisScore > hitScore ) )
                {
                    selfCoord = endCoord;
                    nodeCoord = endCoord + ends[0];
                    hitCoord = thisCoords[!drxn];
                    ol = thisOl;
                    hitScore = thisScore;
                    hitSeq = s;
                }
            }
        }
        if ( hitSeq ) break;
    }
}

