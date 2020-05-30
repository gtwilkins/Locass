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

#include "node.h"

bool Node::calibrateSeed( ExtVars &ev )
{
    int32_t limits[2] = { ends_[0], ends_[1] };
    for ( bool drxn : { 0, 1 } )
    {
        while ( isContinue( drxn ) && ends_[1] - ends_[0] < 20000 )
        {
            extendCount_ = 99;
            ev.ante = getDrxnNodes( !drxn );
            extendNode( ev, drxn );
            limits[0] = min( limits[0], ends_[0] );
            limits[1] = max( limits[1], ends_[1] );
        }
    }
    
    NodeSet delSet;
    for ( bool drxn : { 0, 1 } )
    {
        NodeSet currSet = getNextNodes( drxn );
        NodeSet safeSet;
        if ( currSet.size() == 1 ) safeSet = currSet;
        
        while ( !currSet.empty() && currSet.size() < 15 )
        {
            NodeSet nxtSet;
            for ( Node* curr : currSet )
            {
                while ( curr->isContinue( drxn ) && abs( curr->ends_[drxn] - limits[!drxn] ) < 20000 )
                {
                    curr->extendCount_ = 99;
                    ev.ante = curr->getDrxnNodes( !drxn );
                    curr->extendNode( ev, drxn );
                    limits[0] = min( limits[0], curr->ends_[0] );
                    limits[1] = max( limits[1], curr->ends_[1] );
                }
                
                for ( Node* nxt : curr->getNextNodes( drxn ) )
                {
                    if ( nxt->edges_[!drxn].size() > 1 || ( curr->edges_[drxn].size() == 1 && safeSet.find( curr ) != safeSet.end() ) )
                    {
                        safeSet.insert( nxt );
                    }
                    if ( safeSet.find( nxt ) != safeSet.end() || safeSet.find( curr ) != safeSet.end() )
                    {
                        nxtSet.insert( nxt );
                    }
                }
            }
            
            currSet = nxtSet;
        }
        
        for ( Node* fwd : getDrxnNodes( drxn ) )
        {
            if ( fwd->coverage_ > params.cover * 2 )
            {
                fwd->dismantleNode( delSet, drxn );
            }
        }
    }
    
    for ( auto it = ev.nodes.begin(); it != ev.nodes.end(); )
    {
        if ( delSet.find( *it ) != delSet.end() )
        {
            delete *it;
            it = ev.nodes.erase( it );
            continue;
        }
        it++;
    }
    
    return coverage_ < params.cover * 2;
}

int Node::calibrateCount( NodeList &nodes, SeedLibraryCount &libs )
{
    int totalCount = 0;
    
    for ( Node* n1 : nodes )
    {
        for ( auto &read : n1->reads_ )
        {
            for ( int i ( 0 ); i < params.libs.size(); i++ )
            {
                if ( read.first < params.libs[i].endCount )
                {
                    libs.libCounts[i]++;
                    break;
                }
            }
            if ( !( read.first & 0x1 ) )
            {
                bool isRight = read.first & 0x3;
                ReadId pairId[2] = { read.first + ( isRight ? -2 : 2 )
                                   , read.first + ( isRight ? -1 : 3 ) };
                for ( Node* n2 : nodes )
                {
                    auto it = n2->reads_.find( pairId[0] );
                    if ( it != n2->reads_.end() )
                    {
                        libs.add( read.first, max( abs( read.second[0] - it->second[1] ), abs( read.second[1] - it->second[0] ) ), true, true );
                        totalCount++;
                    }
                    it = n2->reads_.find( pairId[1] );
                    if ( it != n2->reads_.end() )
                    {
                        if ( read.second[0] < it->second[0] && read.second[1] < it->second[1] )
                        {
                            libs.add( read.first, abs( read.second[0] - it->second[1] ), true, false );
                        }
                        else if( it->second[0] < read.second[0] && it->second[1] < read.second[1] )
                        {
                            libs.add( read.first, abs( it->second[0] - read.second[1] ), false, true );
                        }
                        totalCount++;
                    }

                }
            }
        }
    }
    
    return totalCount;
}

void Node::calibrate( NodeList &nodes, LocusLibraryCount &lib )
{
    int32_t limits[ params.libs.size() ][2];
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        limits[i][0] = std::numeric_limits<int32_t>::max();
        limits[i][1] = std::numeric_limits<int32_t>::min();
    }
    
    for ( Node* n1 : nodes )
    {
        if ( n1->coverage_ > params.cover * 2 ) continue;
        for ( auto &read : n1->reads_ )
        {
            if ( read.first & 0x1 ) continue;
            for ( int i ( 0 ); i < params.libs.size(); i++ )
            {
                if ( params.libs[i].endCount <= read.first ) continue;
                bool isRight = read.first & 0x3;
                ReadId pairId[2] = { read.first + ( isRight ? -2 : 2 )
                                   , read.first + ( isRight ? -1 : 3 ) };
                for ( Node* n2 : nodes )
                {
                    auto it = n2->reads_.find( pairId[0] );
                    if ( it != n2->reads_.end() )
                    {
                        lib.ffPairs[i].push_back( max( abs( read.second[0] - it->second[1] ), abs( read.second[1] - it->second[0] ) ) );
                        limits[i][0] = min( limits[i][0], min( it->second[0], read.second[0] ) );
                        limits[i][1] = max( limits[i][1], max( it->second[1], read.second[1] ) );
                        break;
                    }
                    it = n2->reads_.find( pairId[1] );
                    if ( it != n2->reads_.end() )
                    {
                        if ( read.second[0] < it->second[0] && read.second[1] < it->second[1] )
                        {
                            lib.frPairs[i].push_back( abs( read.second[0] - it->second[1] ) );
                        }
                        else if( it->second[0] < read.second[0] && it->second[1] < read.second[1] )
                        {
                            lib.rfPairs[i].push_back( abs( it->second[0] - read.second[1] ) );
                        }
                        limits[i][0] = min( limits[i][0], min( it->second[0], read.second[0] ) );
                        limits[i][1] = max( limits[i][1], max( it->second[1], read.second[1] ) );
                        break;
                    }
                }
                break;
            }
        }
    }
    
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        lib.lens[i] = max( 0, limits[i][1] ) - min( 0, limits[i][0] );
    }
}

