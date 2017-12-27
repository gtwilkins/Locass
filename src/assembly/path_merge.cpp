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

#include "path_merge.h"
#include <algorithm>
#include <cassert>

PathMerge::PathMerge( Node* fork, NodeList &altPath, NodeSet &pathSet, bool drxn )
{
    branch = altPath[0];
    int32_t limit = altPath.back()->ends_[drxn];
    NodeListList qPaths;
    qPaths.push_back( NodeList( altPath.begin(), altPath.end() ) );
    for ( Node* nxt : altPath.back()->getNextNodes( drxn ) )
    {
        limit = drxn ? max( nxt->ends_[1], limit ) : min( nxt->ends_[0], limit );
        NodeList qPath = { altPath.back() };
        qPath.push_back( nxt );
        qPaths.push_back( qPath );
    }
    limit = drxn ? limit + params.readLen * 1.5 : limit - params.readLen * 1.5;
    
    for ( NodeList &qPath : qPaths )
    {
        Node* qAnchor = qPath[0];
        if ( !drxn ) reverse( qPath.begin(), qPath.end() );
        qSeqs_.push_back( new SeqPathMerge( qPath, qAnchor, 0 ) );
    }
    
    NodeListList tPaths;
    for ( Node* nxt : fork->getNextNodes( drxn ) )
    {
        if ( pathSet.find( nxt ) == pathSet.end() ) continue;
        NodeList path = { nxt };
        while ( drxn ? path.back()->ends_[1] < limit : limit < path.back()->ends_[0] )
        {
            NodeList branches;
            for ( Node* pathNxt : path.back()->getNextNodes( drxn ) )
            {
                if ( pathSet.find( pathNxt ) != pathSet.end() ) branches.push_back( pathNxt );
            }
            if ( branches.empty() ) break;
            path.push_back( branches[0] );
            if ( branches.size() > 1 )
            {
                for ( int i = 1; i < branches.size(); i++ )
                {
                    NodeList tPath = { path.end()[-2] };
                    tPath.push_back( branches[i] );
                    tPaths.push_back( tPath );
                }
                break;
            }
        }
        tPaths.push_back( path );
    }
    
    for ( NodeList &tPath : tPaths )
    {
        Node* tAnchor = tPath[0];
        if ( !drxn ) reverse( tPath.begin(), tPath.end() );
        tSeqs_.push_back( new SeqPathMerge( tPath, tAnchor, 0 ) );
    }
}

PathMerge::~PathMerge()
{
    for ( SeqPathMerge* s : qSeqs_ )
    {
        delete s;
    }
    for ( SeqPathMerge* s : tSeqs_ )
    {
        delete s;
    }
}

bool PathMerge::merge( PathVars &pv, NodeSet &delSet, bool drxn )
{
    SeqPathMerge* bestHit = NULL;
    
    bool endOnly = false;
    for ( SeqPathMerge* q : qSeqs_ )
    {
        q->merge( tSeqs_, endOnly, drxn );
        if ( q->hitSeq && ( !bestHit || ( drxn ? q->nodeCoord > bestHit->nodeCoord
                                               : q->nodeCoord < bestHit->nodeCoord ) ) ) 
        {
            bestHit = q;
        }
        endOnly = false;
    }
    
    bool didMerge = false;
    if ( bestHit ) didMerge = bestHit->doMerge( pv, delSet, drxn );
    
    return didMerge ;
}
