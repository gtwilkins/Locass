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

#include "node_structs.h"
#include "node.h"
#include <algorithm>
#include "shared_functions.h"

extern struct Parameters params;

using namespace std;

Coords::Coords( int32_t start, int32_t end, bool isRedundant )
: redundant( isRedundant )
{
    coords[0] = start;
    coords[1] = end;
};

void Coords::offset( int32_t offset )
{
    coords[0] += offset;
    coords[1] += offset;
}

int32_t& Coords::operator []( bool i )
{
    return coords[i];
}

Score::Score()
{
    hits = 0;
    misses = 0;
}

void Score::clear()
{
    hits = 0;
    misses = 0;
}

void Score::operator +=( const Score& rhs )
{
    hits += rhs.hits;
    misses += rhs.misses; 
}

void Score::operator -=( const Score& rhs )
{
    hits = max( hits - rhs.hits, (float)0 );
    misses = max( hits - rhs.misses, (float)0 ); 
}

bool Score::operator > ( const Score &rhs )
{
    return hits == rhs.hits ? misses < rhs.misses : hits > rhs.hits;
}

float& Score::operator [](int i)
{
    if( i ) return hits;
    else return misses;
}

float Score::ratio( float hits, float misses )
{
    return hits + hits ? ( hits + hits ) / ( hits + misses + hits + misses ) : 0;
}

void PathVars::addTarget( Node* t, NodeSet &fwdSet )
{
    tSet.insert( t );
    float cloneCount = t->clones_ ? 1 + t->clones_->size() : 1;
    float coeff = max( ( t->coverage_ / params.cover ) / float( cloneCount + 0.3 ), float( cloneCount ) );
    for ( auto &np : t->pairs_ )
    {
        if ( fwdSet.find( np.first ) != fwdSet.end() )
        {
            furthest = drxn ? max( furthest, t->ends_[0] + params.readLen )
                            : min( furthest, t->ends_[1] - params.readLen );
            
            // Set hits
            auto r = hits.insert( np );
            if ( !r.second ) r.first->second += np.second;
            
            // Set reliable hits
            if ( t->isReliable() )
            {
                r = reliable.insert( np );
                if ( !r.second ) r.first->second += np.second;
            }
            
            // Set adjustedHits
            float adjHits = np.second / ( t->drxn_ == np.first->drxn_ ? coeff : 1 );
            auto r2 = adjusted.insert( make_pair( np.first, adjHits ) );
            if ( !r2.second ) r2.first->second += adjHits;
            
            // Update furthest if necessary
            if ( !np.first->farPairNodes_[0]
                    || ( np.first->farPairNodes_[0] && np.first->farPairNodes_[0] != t
                            && ( drxn ? t->ends_[0] < np.first->farPairCoords_[0]
                                      : np.first->farPairCoords_[0] < t->ends_[1] ) ) )
            {
                np.first->setFurthest( tSet, drxn );
            }
        }
    }
}

void PathVars::addTarget( NodeList &path )
{
    if ( !path.empty() )
    {
        NodeSet fwdSet = path.back()->getDrxnNodes( drxn );
        for ( Node* node : path )
        {
            addTarget( node, fwdSet );
        }
    }
}

void PathVars::resetPathReliable( NodeList &path )
{
    NodeSet resetSet;
    for ( Node* node : path )
    {
        if ( !node->isReliable() )
        {
            NodeSet fwdSet = node->getDrxnNodes( drxn );
            for ( auto &np : node->pairs_ )
            {
                if ( fwdSet.find( np.first ) != fwdSet.end() )
                {
                    resetSet.insert( np.first );
                    if ( np.first->farPairNodes_[1] && np.first->farPairNodes_[1] == node  )
                    {
                        np.first->farPairNodes_[1] = NULL;
                    }
                }
            }
        }
    }
    
    for ( Node* node : resetSet )
    {
        int reliableCount = 0;
        for ( auto &np : node->pairs_ )
        {
            if ( np.first->isReliable() && tSet.find( np.first ) != tSet.end() )
            {
                reliableCount += np.second;
            }
        }
        
        reliable[node] = reliableCount;
    }
}

void PathVars::setTarget( NodeSet &fwdSet, NodeList &islands )
{
    for ( Node* t : tSet )
    {
        addTarget( t, fwdSet );
    }
    
    for ( Node* node : islands )
    {
        for ( auto &np : node->pairs_ )
        {
            if ( fwdSet.find( np.first ) != fwdSet.end() )
            {
                auto r = islandHits.insert( np );
                if ( !r.second ) r.first->second += np.second;
            }
        }
    }
}

void CloneScore::addHit( Coords* coords, bool drxn )
{
    if ( drxn )
    {
        if ( (*coords)[0] < furthest[0] )
        {
            furthest[1] = furthest[0];
            furthest[0] = (*coords)[0];
        }
        else
        {
            furthest[1] = min( furthest[1], (*coords)[0] );
        }
    }
    else
    {
        if ( (*coords)[1] > furthest[0] )
        {
            furthest[1] = furthest[0];
            furthest[0] = (*coords)[1];
        }
        else
        {
            furthest[1] = max( furthest[1], (*coords)[1] );
        }
    }
}

void CloneScore::mergeFurthest( const CloneScore &rhs, bool drxn )
{
    vector<int32_t> combined = { furthest[0], furthest[1], rhs.furthest[0], rhs.furthest[1] };
    sort( combined.begin(), combined.end(), [&drxn]( int32_t &a, int32_t &b ){
        return ( drxn ? a < b : a > b );
    });
    furthest[0] = combined[0];
    furthest[1] = combined[1];
}

void CloneScore::operator += ( const CloneScore &rhs )
{
    selfHits += rhs.selfHits;
    selfPref += rhs.selfPref;
    selfOnly += rhs.selfOnly;
    selfScore += rhs.selfScore;
    altPref += rhs.altPref;
    altOnly += rhs.altOnly;
    altScore += rhs.altScore;
}

bool CloneScore::operator > ( const CloneScore &rhs )
{
    if ( selfOnly > 0 && rhs.selfOnly > 0 )
    {
        if ( furthest[1] > rhs.furthest[0] )
        {
            return true;
        }
        if ( rhs.furthest[1] > furthest[0] )
        {
            return false;
        }
    }
    if ( selfOnly > 0 && rhs.selfOnly < 1 && furthest[0] > rhs.furthest[0] )
    {
        return true;
    }
    if ( rhs.selfOnly > 0 && selfOnly < 1 && rhs.furthest[0] > furthest[0] )
    {
        return false;
    }
    if ( cumulScore > rhs.cumulScore )
    {
        return true;
    }
    return false;
}

void CloneScore::setScores()
{
    selfScore = selfPref + selfOnly * 4;
    altScore = ( altPref + altOnly * 4 );
    cumulScore = selfScore - altScore;
}

void CloneScore::setCumul( Node* node, int prevScore, int prevMisses )
{
    float multiplicity = min( (float)4, node->getMultiplicity() );
    int len = max( 1, node->ends_[1] - node->ends_[0] - params.readLen );
    int tmpAltScore = ( multiplicity > 1 ? max( int( altScore > 0 ), int( altScore / ( multiplicity + 0.5 ) ) )
                                         : altScore );
    
    cumulMisses = selfScore > tmpAltScore ? 0 : prevMisses;
    
    cumulScore = selfScore - tmpAltScore + prevScore;
    
    if ( multiplicity < 1 )
    {
        cumulMisses += 1 + ( ( len * ( (float)1.5 - multiplicity ) ) / 5 );
    }
    else if ( !selfScore )
    {
        cumulMisses += 1 + ( len / 15 );
    }
}

ReadMark::ReadMark( SeqNum &readId, Coords &inCoords, Lib* lib, bool drxn )
: id( readId ), mark( inCoords[!drxn] )
{
    estimate = drxn ? inCoords[0] + lib->size : inCoords[1] - lib->size;
    coords[0] = mark + ( drxn ? (*lib).minDist : -(*lib).maxDist );
    coords[1] = mark + ( drxn ? (*lib).maxDist : -(*lib).minDist );
}

void ReadMark::getMissScore( float &missScore, int32_t* limits )
{
    int maxEnds[2] = { max( limits[0], coords[0] ), min( limits[1], coords[1] ) };
    missScore += max( (float)0, ( (float)( maxEnds[1] - maxEnds[0] ) / (float)( coords[1] - coords[0] ) ) );
}

bool ReadMark::isValid( Coords &hitCoords )
{
    return coords[0] <= hitCoords[0] && hitCoords[1] <= coords[1];
}

bool ReadMark::isValid( Coords &hitCoords, pair<int32_t, int32_t> &markOffset, pair<int32_t, int32_t> &hitOffset, bool drxn )
{
    int32_t limits[2] = { coords[0] - markOffset.first, coords[1] - markOffset.second };
    return limits[0] <= hitCoords[0] - hitOffset.second 
            && hitCoords[1] - hitOffset.first <= limits[1];
}

void ReadMark::offset( int32_t offset )
{
    coords[0] += offset;
    coords[1] += offset;
    estimate += offset;
    mark += offset;
}

bool LoopVars::addLoopStart( Node* node, bool drxn )
{
    
    for ( Node* clone : node->getCloneSet( true ) )
    {
        cloned.insert( clone );
        if ( startList.empty() )
        {
            cloneFwdSets.push_back( clone->getDrxnNodes( drxn, true, true ) );
        }
        else
        {
            NodeSet fwdSet = clone->getDrxnNodes( drxn, true, true );
            bool didAdd = false;
            for ( NodeSet &cloneFwd : cloneFwdSets )
            {
                bool converges = false;
                bool divergent = false;
                if ( cloneFwd.find( clone ) == cloneFwd.end() )
                {
                    for ( Node* fwd : cloneFwd )
                    {
                        ( fwdSet.find( fwd ) != fwdSet.end() ? converges : divergent ) = true;
                    }
                }
                
                if ( converges && divergent )
                {
                    cloneFwd.insert( fwdSet.begin(), fwdSet.end() );
                    didAdd = true;
                    break;
                }
            }
            if ( !didAdd ) return false;
        }
    }
    
    startList.push_back( node );
    clones.insert( node );
    open.insert( node );
    
    return true;
}

float Edge::getOverlapWeakness()
{
    return min( float(1), float( max( 0, params.readLen - overlap ) ) / float( max( 1, params.readLen ) ) ); 
}

MapStruct::MapStruct( string &inSeq, int32_t inMinLen, int32_t inCoord, int32_t target, bool drxn )
: seq( inSeq ), len( inSeq.length() ), coord( inCoord )
{
    // drxn is which side of the mapping this seq is
    minLen = drxn ? len - inMinLen : inMinLen;
    cutLen = drxn ? len - coord + target : target - coord;
}

IslandRead::IslandRead( string inSeq, ReadMark &mark )
: seq( inSeq ), id( mark.id )
{
    edges[0] = edges[1] = NULL;
    coord = mark.estimate < mark.mark ? mark.estimate : mark.estimate - seq.length();
}

void IslandRead::setOverlap( IslandRead &read )
{
    int overLen = ( ( edges[1] && read.edges[0] ) ? min( overlaps[1], read.overlaps[0] ) : 30 );
    size_t it1 = seq.length() - overLen;
    string q = seq.substr( it1 );
    size_t it2 = read.seq.find( q );
    while ( it2 != read.seq.npos && it2 < it1 )
    {
        while ( it2 > 0 && seq[it1-1] == read.seq[it2-1] )
        {
            it1--;
            it2--;
            overLen++;
        }
        if ( !it2 )
        {
            if ( !edges[1] || overLen > overlaps[1] )
            {
                edges[1] = &read;
                overlaps[1] = overLen;
            }
            if ( !read.edges[0] || overLen > read.overlaps[0] )
            {
                read.edges[0] = this;
                read.overlaps[0] = overLen;
            }
        }
        it1 = seq.length() - ++overLen;
        q = seq.substr( it1 );
        it2 = read.seq.find( q );
    }
}

void MapNode::addEdge( MapNode* mapNode, int overlap, bool drxn )
{
    edges[drxn].push_back( mapNode );
    edgeOverlaps[drxn].push_back( overlap );
    mapNode->edges[!drxn].push_back( this );
    mapNode->edgeOverlaps[!drxn].push_back( overlap );
}

void MapNode::addEdge( Node* edgeNode, int overlap, bool drxn )
{
    bridges[drxn].push_back( edgeNode );
    bridgeOverlaps[drxn].push_back( overlap );
}

void MapNode::checkLoop()
{
    vector<ReadId>::iterator found = find( ids.begin() + 1, ids.end(), ids[0] );
    while( found != ids.end() )
    {
        int i = found - ids.begin();
        ids.erase( ids.begin(), found );
        coords[0].erase( coords[0].begin(), coords[0].begin()+i );
        coords[1].erase( coords[1].begin(), coords[1].begin()+i );
        found = find( ids.begin() + 1, ids.end(), ids[0] );
        int trim = coords[0][0];
        for ( int j = 0; j < ids.size(); j++ )
        {
            coords[0][j] -= trim;
            coords[1][j] -= trim;
        }
        seq = seq.substr( trim );
    }
    
    found = find( ids.begin(), ids.end()-1, ids.back() );
    while ( found != ids.end()-1 )
    {
        int i = found - ids.begin();
        ids.erase( found+1, ids.end() );
        coords[0].erase( coords[0].begin()+i+1, coords[0].end() );
        coords[1].erase( coords[1].begin()+i+1, coords[1].end() );
        found = find( ids.begin(), ids.end()-1, ids.back() );
        int limit = 0;
        for ( int j = 0; j < ids.size(); j++ )
        {
            limit = max( limit, coords[1][j] );
        }
        seq = seq.substr( 0, limit );
    }
    
    for ( int i = 0; i < ids.size(); i++ )
    {
        found = find( ids.begin()+i+1, ids.end(), ids[i] );
        assert( found == ids.end() );
    }
}

void MapNode::collapse( vector<MapNode*> &mapNodes, vector<MapNode*> &mapEdges )
{
    while ( !mapEdges.empty() )
    {
        MapNode* forks[2] = { mapEdges[0]->edges[0].empty() ? NULL : mapEdges[0]->edges[0][0]
                            , mapEdges[0]->edges[1].empty() ? NULL : mapEdges[0]->edges[1][0] };
        assert( forks[0] || forks[1] );
        bool drxn = forks[1] && ( forks[1]->edges[0].size() > 1 || !forks[0] );
        if ( forks[drxn]->edges[!drxn].size() == 1 
                && forks[drxn]->bridges[!drxn].empty() 
                && mapEdges[0]->edges[!drxn].empty() )
        {
            MapNode::fold( mapEdges[0], mapEdges, drxn );
            continue;
        }
        vector<MapNode*> branches;
        for ( MapNode* edge : forks[drxn]->edges[!drxn] )
        {
            if ( find( mapEdges.begin(), mapEdges.end(), edge ) != mapEdges.end() ) branches.push_back( edge );
        }
        bool doTransfer = branches.size() == forks[drxn]->edges[!drxn].size() && forks[drxn]->bridgeOverlaps[!drxn].empty();
        
        vector<int> hits;
        int i = 0;
        while ( i < branches[0]->ids.size() )
        {
            for ( int j = 1; j < branches.size(); j++ )
            {
                if ( branches[j]->ids.size() > i &&
                        ( drxn ? branches[j]->ids.end()[-i-1] == branches[0]->ids.end()[-i-1]
                               : branches[j]->ids[i] == branches[0]->ids[i] ) )
                {
                    hits.push_back( j );
                }
            }
            doTransfer = doTransfer && ( hits.size() == branches.size() - 1 );
            if ( !doTransfer ) break;
            hits.clear();
            i++;
        }

        if ( i )
        {
            vector<ReadId> transIds = { ( drxn ? branches[0]->ids.end()-i : branches[0]->ids.begin() )
                                      , ( drxn ? branches[0]->ids.end() : branches[0]->ids.begin()+i ) };
            vector<int32_t> transCoords[2];
            transCoords[0] = { ( drxn ? branches[0]->coords[0].end()-i : branches[0]->coords[0].begin() )
                             , ( drxn ? branches[0]->coords[0].end() : branches[0]->coords[0].begin()+i ) };
            transCoords[1] = { ( drxn ? branches[0]->coords[1].end()-i : branches[0]->coords[1].begin() )
                             , ( drxn ? branches[0]->coords[1].end() : branches[0]->coords[1].begin()+i ) };
            int32_t transLimits[2] = { *min_element( transCoords[0].begin(), transCoords[0].end() )
                                     , *max_element( transCoords[1].begin(), transCoords[1].end() ) };
            int32_t extra = transLimits[1] - transLimits[0] - branches[0]->edgeOverlaps[drxn][0];
            int32_t offset = ( drxn ? -transLimits[0] : forks[0]->seq.length() - branches[0]->edgeOverlaps[0][0] );
            forks[drxn]->seq = drxn ? branches[0]->seq.substr( transLimits[0], extra ) + forks[1]->seq
                                    : forks[0]->seq + branches[0]->seq.substr( transLimits[1] - extra, extra );
            for ( int j = 0; j < transIds.size(); j++ )
            {
                transCoords[0][j] += offset;
                transCoords[1][j] += offset;
            }

            if ( drxn )
            {
                for ( int j = 0; j < forks[1]->ids.size(); j++ )
                {
                    forks[1]->coords[0][j] += extra;
                    forks[1]->coords[1][j] += extra;
                }
                forks[1]->ids.insert( forks[1]->ids.begin(), transIds.begin(), transIds.end() );
                forks[1]->coords[0].insert( forks[1]->coords[0].begin(), transCoords[0].begin(), transCoords[0].end() );
                forks[1]->coords[1].insert( forks[1]->coords[1].begin(), transCoords[1].begin(), transCoords[1].end() );
            }
            else
            {
                forks[0]->ids.insert( forks[0]->ids.end(), transIds.begin(), transIds.end() );
                forks[0]->coords[0].insert( forks[0]->coords[0].end(), transCoords[0].begin(), transCoords[0].end() );
                forks[0]->coords[1].insert( forks[0]->coords[1].end(), transCoords[1].begin(), transCoords[1].end() );
            }

            for ( MapNode* b : branches )
            {
                b->ids.erase( ( drxn ? b->ids.end()-i : b->ids.begin() ), ( drxn ? b->ids.end() : b->ids.begin()+i ) );
                b->coords[0].erase( ( drxn ? b->coords[0].end()-i : b->coords[0].begin() )
                                  , ( drxn ? b->coords[0].end() : b->coords[0].begin()+i ) );
                b->coords[1].erase( ( drxn ? b->coords[1].end()-i : b->coords[1].begin() )
                                  , ( drxn ? b->coords[1].end() : b->coords[1].begin()+i ) );
                if ( !b->ids.empty() )
                {
                    int32_t limit = drxn ? *max_element( b->coords[1].begin(), b->coords[1].end() ) 
                                         : *min_element( b->coords[0].begin(), b->coords[0].end() );
                    b->edgeOverlaps[drxn][0] = transLimits[1] - transLimits[0] - ( drxn ? b->seq.length() - limit : limit );
                    b->seq = drxn ? b->seq.substr( 0, limit ) : b->seq.substr( limit );
                    for ( int j = 0; j < forks[drxn]->edges[!drxn].size(); j++ )
                    {
                        if ( forks[drxn]->edges[!drxn][j] == b )
                        {
                            forks[drxn]->edgeOverlaps[!drxn][j] = b->edgeOverlaps[drxn][0];
                        }
                    }

                    if ( !drxn )
                    {
                        for ( int j = 0; j < b->ids.size(); j++ )
                        {
                            b->coords[0][j] -= limit;
                            b->coords[1][j] -= limit;
                        }
                    }
                }
            }
        }

        if ( hits.empty() )
        {
            if ( branches[0]->edges[!drxn].empty() || !branches[0]->bridges[!drxn].empty() )
            {
                if ( branches[0]->ids.empty() )
                {
                    MapNode::fold( branches[0], mapEdges, drxn );
                }
                else
                {
                    mapNodes.push_back( branches[0] );
                    mapEdges.erase( remove( mapEdges.begin(), mapEdges.end(), branches[0] ), mapEdges.end() );
                }
            }
            else
            {
                MapNode::fold( branches[0], mapEdges, !drxn );
            }
            branches.erase( branches.begin() );
        }
        else
        {
            mapNodes.push_back( new MapNode() );
            mapNodes.back()->addEdge( forks[drxn], branches[0]->edgeOverlaps[drxn][0], drxn );
            mapNodes.back()->estimate = forks[drxn]->estimate;
            hits.insert( hits.begin(), 0 );
            for ( int hit : hits )
            {
                branches[hit]->edges[drxn][0] = mapNodes.back();
                branches[hit]->edgeOverlaps[drxn][0] = 0;
                mapNodes.back()->edges[!drxn].push_back( branches[hit] );
                mapNodes.back()->edgeOverlaps[!drxn].push_back( 0 );
                for ( int j = forks[drxn]->edges[!drxn].size(); --j >= 0 ; )
                {
                    if ( forks[drxn]->edges[!drxn][j] == branches[hit] )
                    {
                        forks[drxn]->edges[!drxn].erase( forks[drxn]->edges[!drxn].begin() + j );
                        forks[drxn]->edgeOverlaps[!drxn].erase( forks[drxn]->edgeOverlaps[!drxn].begin() + j );
                    }
                }
            }
        }
    }
}

void MapNode::fold( MapNode* mapNode, vector<MapNode*> &mapNodes, bool drxn )
{
    MapNode* mn[2] = { ( drxn ? mapNode : mapNode->edges[0][0] ), ( drxn ? mapNode->edges[1][0] : mapNode ) };
    for ( int i = 0; i < mn[!drxn]->bridges[!drxn].size(); i++ )
    {
        mn[drxn]->bridges[!drxn].push_back( mn[!drxn]->bridges[!drxn][i] );
        mn[drxn]->bridgeOverlaps[!drxn].push_back( mn[!drxn]->bridgeOverlaps[!drxn][i] );
    }
    
    int ol = mapNode->edgeOverlaps[drxn][0];
    if ( !mapNode->ids.empty() )
    {
        int32_t extra = mn[0]->seq.length() - mn[0]->edgeOverlaps[1][0];
        mn[drxn]->seq = mn[0]->seq.substr( 0, extra ) + mn[1]->seq;
        for ( int i = 0; i < mn[1]->ids.size(); i++ )
        {
            mn[0]->ids.push_back( mn[1]->ids[i] );
            mn[0]->coords[0].push_back( mn[1]->coords[0][i] + extra );
            mn[0]->coords[1].push_back( mn[1]->coords[1][i] + extra );
        }
        mn[1]->ids = mn[0]->ids;
        mn[1]->coords[0] = mn[0]->coords[0];
        mn[1]->coords[1] = mn[0]->coords[1];
        ol = mapNode->edges[!drxn].empty() ? ol : mapNode->edgeOverlaps[!drxn][0];
    }
    
    if ( !mapNode->edges[0].empty() && !mapNode->edges[1].empty() )
    {
        for ( int i : { 0, 1 } )
        {
            for ( int j = 0; j < mapNode->edges[i][0]->edges[!i].size(); j++ )
            {
                if ( mapNode->edges[i][0]->edges[!i][j] == mapNode )
                {
                    mapNode->edges[i][0]->edges[!i][j] = mapNode->edges[!i][0];
                    mapNode->edges[i][0]->edgeOverlaps[!i][j] = ol;
                }
            }
        }
    }
    else
    {
        mn[drxn]->removeEdge( mn[!drxn], !drxn );
    }
    
    mapNodes.erase( remove( mapNodes.begin(), mapNodes.end(), mapNode ), mapNodes.end() );
    delete mapNode;
}

bool MapNode::getSecondSeq( string &secondSeq, bool drxn )
{
    for ( int i = 0; i < ids.size(); i++ )
    {
        if ( drxn && coords[1].end()[-i-1] != seq.length() )
        {
            secondSeq = seq.substr( 0, coords[1].end()[-i-1] );
            return true;
        }
        if ( !drxn && coords[0][i] != 0 )
        {
            secondSeq = seq.substr( coords[0][i] );
            return true;
        }
    }
    return false;
}

bool MapNode::recoil()
{
    if ( ids.empty() )
    {
        for ( int drxn : { 0, 1 } )
        {
            if ( !edges[drxn].empty() )
            {
                for ( int i = 0; i < edges[drxn][0]->edges[!drxn].size(); i++ )
                {
                    if ( edges[drxn][0]->edges[!drxn][i] == this )
                    {
                        if ( !edges[!drxn].empty() )
                        {
                            edges[drxn][0]->edges[!drxn][i] = edges[!drxn][0];
                        }
                        else if ( !bridges[!drxn].empty() )
                        {
                            edges[drxn][0]->edges[!drxn].erase( edges[drxn][0]->edges[!drxn].begin() + i );
                            edges[drxn][0]->edgeOverlaps[!drxn].erase( edges[drxn][0]->edgeOverlaps[!drxn].begin() + i );
                            if ( !bridges[!drxn].empty() )
                            {
                                bridgeOverlaps[!drxn][0] -= seq.length() - edgeOverlaps[drxn][0];
                                edges[drxn][0]->bridges[!drxn].push_back( bridges[!drxn][0] );
                                edges[drxn][0]->bridgeOverlaps[!drxn].push_back( bridgeOverlaps[!drxn][0] );
                            }
                        }
                        break;
                    }
                }
            }
        }
        
        return false;
    }
    
    int32_t limits[2] = { *min_element( coords[0].begin(), coords[0].end() )
                        , *max_element( coords[1].begin(), coords[1].end() ) };
    int32_t trims[2] = { limits[0], (int32_t)seq.length() - limits[1] };
    seq = seq.substr( limits[0], limits[1] - limits[0] );
    for ( int drxn : { 0, 1 } )
    {
        for ( int i = 0; i < edges[drxn].size(); i++ )
        {
            edgeOverlaps[drxn][i] -= trims[drxn];
            for ( int j = 0; j < edges[drxn][i]->edges[!drxn].size(); j++ )
            {
                if ( edges[drxn][i]->edges[!drxn][j] == this )
                {
                    edges[drxn][i]->edgeOverlaps[!drxn][j] = edgeOverlaps[drxn][i];
                }
            }
        }
        for ( int i = 0; i < coords[drxn].size(); i++ )
        {
            coords[drxn][i] -= trims[0];
        }
        for ( int i = 0; i < bridges[drxn].size(); i++ )
        {
            bridgeOverlaps[drxn][i] -= trims[drxn];
        }
    }
    return true;
}

void MapNode::removeEdge( MapNode* mn, bool drxn )
{
    for ( int i = 0; i < edges[drxn].size(); )
    {
        if ( edges[drxn][i] == mn )
        {
            edges[drxn].erase( edges[drxn].begin() + i );
            edgeOverlaps[drxn].erase( edgeOverlaps[drxn].begin() + i );
            continue;
        }
        i++;
    }
}

void MapNode::removeEdges( bool drxn )
{
    for ( MapNode* e : edges[drxn] )
    {
        e->removeEdge( this, !drxn );
    }
    edges[drxn].clear();
    edgeOverlaps[drxn].clear();
}

void MapNode::setEdgeSeq()
{
    int32_t edgeCoords[2] = { *max_element( edges[0][0]->coords[0].begin(), edges[0][0]->coords[0].end() )
                            , *min_element( edges[1][0]->coords[1].begin(), edges[1][0]->coords[1].end() ) };
    
    seq = edges[0][0]->seq.substr( ++edgeCoords[0] ) 
        + edges[1][0]->seq.substr( edgeOverlaps[1][0], --edgeCoords[1] - edgeOverlaps[1][0] );
    
    edgeOverlaps[0][0] = edges[0][0]->seq.length() - edgeCoords[0];
    edgeOverlaps[1][0] = edgeCoords[1];
}

void MapNode::setSecondSeq( bool drxn )
{
    for ( int i = 0; i < ids.size(); i++ )
    {
        if ( drxn && coords[1].end()[-i-1] != seq.length() )
        {
            seq = seq.substr( 0, coords[1].end()[-i-1] );
            ids.erase( ids.end()-i, ids.end() );
            coords[0].erase( coords[0].end()-i, coords[0].end() );
            coords[1].erase( coords[1].end()-i, coords[1].end() );
            return;
        }
        if ( !drxn && coords[0][i] != 0 )
        {
            seq = seq.substr( coords[0][i] );
            ids.erase( ids.begin(), ids.begin() + i );
            coords[0].erase( coords[0].begin(), coords[0].begin() + i );
            coords[1].erase( coords[1].begin(), coords[1].begin() + i );
            int32_t trim = coords[0][0];
            for ( int j = 0; j < ids.size(); j++ )
            {
                coords[0][j] -= trim;
                coords[1][j] -= trim;
            }
            return;
        }
    }
}

void MapNode::splitEdge( Node* node, int i, bool drxn )
{
    MapNode* mn = new MapNode();
}

bool NodeMapRead::checkDoubleHit()
{
    if ( nodes[0] && nodes[1] )
    {
        bool drxn = lens[1] > lens[0];
        bool doUnset = true;
        if ( coords[0][1] == nodes[0]->ends_[1] && nodes[1]->ends_[0] == coords[1][0] )
        {
            for ( Edge &e : nodes[0]->edges_[1] )
            {
                if ( e.node == nodes[1] 
                        && ( coords[1][1] - coords[1][0] > e.overlap || e.isLeap ) )
                {
                    doUnset = false;
                }
            }
        }
        if ( doUnset )
        {
            nodes[!drxn] = NULL;
            lens[!drxn] = 0;
            return false;
        }
        return true;
    }
    return false;
}

bool NodeMapRead::checkHit( ReadMark &mark, bool markDrxn, bool drxn )
{
    if ( nodes[drxn] )
    {
        int estimate = drxn != markDrxn ? ( drxn ? mark.estimate - params.readLen 
                                                 : mark.estimate + params.readLen ) 
                                        : mark.estimate;
        if ( nodes[drxn]->reads_.find( id ) != nodes[drxn]->reads_.end()
                || abs( estimate - coords[drxn][markDrxn] ) > 150 )
        {
            nodes[drxn] = NULL;
            lens[drxn] = 0;
        }
        else
        {
            return true;
        }
    }
    return false;
}

void NodeMapReadHits::add( NodeSet hitNodes, int32_t* hitCoords, bool doAdd, bool drxn )
{
    int iHome = -1;
    for ( int i = 0; i < coords[0].size(); )
    {
        if ( hitCoords[0] <= coords[1][i] && coords[0][i] <= hitCoords[1]  )
        {
            iHome = iHome == -1 ? i : iHome;
            coords[0][iHome] = min( coords[0][iHome], min( coords[0][i], hitCoords[0] ) );
            coords[1][iHome] = max( coords[1][iHome], max( coords[1][i], hitCoords[1] ) );
            nodes[iHome].insert( hitNodes.begin(), hitNodes.end() );
            hits[drxn][iHome]++;
            doAdd = false;
            if ( i != iHome )
            {
                nodes[iHome].insert( nodes[i].begin(), nodes[i].end() );
                nodes.erase( nodes.begin() + i );
                hits[0][iHome] += hits[0][i];
                hits[1][iHome] += hits[1][i];
                hits[0].erase( hits[0].begin() + i );
                hits[1].erase( hits[1].begin() + i );
                coords[0].erase( coords[0].begin() + i );
                coords[1].erase( coords[1].begin() + i );
                continue;
            }
        }
        i++;
    }
    
    if ( doAdd )
    {
        nodes.push_back( hitNodes );
        hits[drxn].push_back( 1 );
        hits[!drxn].push_back( 0 );
        coords[0].push_back( hitCoords[0] );
        coords[1].push_back( hitCoords[1] );
    }
}

MappedReadEnd::MappedReadEnd( string inSeq, ReadId inId, int32_t* inCoords, int32_t inOffset, bool endDrxn, bool fromDrxn )
: seq( inSeq ), id( inId ), offset( inOffset ), drxn( fromDrxn ), doMap( false )
{
    coords[0 + endDrxn] = inCoords[0];
    coords[1 + endDrxn] = inCoords[1];
    coords[(!endDrxn)*2] = endDrxn ?  coords[2] - seq.length() : coords[0]  + seq.length();
    counts[0] = counts[1] = 0;
    counts[fromDrxn] = 1;
}

bool MappedReadEnd::checkMap( int score, bool overDrxn )
{
    if ( doMap )
    {
        int overlap = abs( coords[overDrxn*2] - coords[1] );
        return abs( offset ) - ( overlap / 2 ) - score < 100;
    }
    return false;
}

string MappedReadEnd::getExtSeq( bool drxn )
{
    int len = min( (int)seq.length(), abs( coords[1] - coords[drxn*2] ) );
    return ( drxn ? seq.substr( seq.length() - len ) : seq.substr( 0, len ) );
}

int MappedReadEnd::getMinOverlap( int minlen, bool overDrxn )
{
    int i = 0;
    while ( i < seq.length() &&
            ( overDrxn ? seq[i] == seq[i+1] : seq[seq.length()-i-1] == seq[seq.length()-i-2] ) ) i++;
    return minlen + i;
}

int MappedReadEnd::getOverlap( bool overDrxn )
{
    return abs( coords[overDrxn*2] - coords[1] );
}

void MappedReadEnd::set( MappedReadEnd &read )
{
    seq = read.seq ;
    id = read.id;
    coords[0] = read.coords[0];
    coords[1] = read.coords[1];
    coords[2] = read.coords[2];
    counts[0] = counts[1] = 0;
    drxn = read.drxn;
    counts[drxn] = 1;
}

bool MappedReadEnd::withinLimits( int32_t* limits )
{
    return ( limits[0] < coords[2] && coords[0] < limits[1] );
}

void PathRead::offset( vector<PathRead*> reads, int offset )
{
    for ( PathRead* read : reads )
    {
        read->coord += offset;
    }
}

void PathRead::setOverlaps( vector<PathRead*> &reads )
{
    for ( int i = 0; i < reads.size(); i++ )
    {
        for ( int j = i+1; j < reads.size(); )
        {
            if ( reads[i]->seq != reads[j]->seq )
            {
                j++;
                continue;
            }
            reads[i]->ids.insert( reads[i]->ids.end(), reads[j]->ids.begin(), reads[j]->ids.end() );
            reads[i]->coord = ( reads[i]->coord + reads[j]->coord ) / 2;
            delete reads[j];
            reads.erase( reads.begin() + j );
        }
    }
    for ( PathRead* r1 : reads )
    {
        for ( PathRead* r2 : reads )
        {
            if ( abs( r1->coord - r2->coord ) > 500 ) continue;
            if ( r1 == r2 ) continue;
            int ols[2]{0};
            ols[0] = mapSeqOverlap( r2->seq, r1->seq, max( r1->ols[0]+1, 25 ) );
            ols[1] = mapSeqOverlap( r1->seq, r2->seq, max( r1->ols[1]+1, 25 ) );
            if ( ols[0] )
            {
                r1->edges[0] = r2;
                r1->ols[0] = ols[0];
                if ( ols[0] > r2->ols[1] )
                {
                    r2->edges[1] = r1;
                    r2->ols[1] = ols[0];
                }
            }
            if ( ols[1] )
            {
                r1->edges[1] = r2;
                r1->ols[1] = ols[1];
                if ( ols[1] > r2->ols[0] )
                {
                    r2->edges[0] = r1;
                    r2->ols[0] = ols[1];
                }
            }
        }
    }
}