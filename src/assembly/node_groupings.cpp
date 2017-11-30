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
#include <limits>
#include <algorithm>

void Node::getLimits( int32_t* limits, NodeList &nodes, bool doReset )
{
    if ( doReset )
    {
        limits[0] = std::numeric_limits<int32_t>::max();
        limits[1] = std::numeric_limits<int32_t>::min();
    }
    for ( Node* node : nodes )
    {
        limits[0] = min( limits[0], node->ends_[0] );
        limits[1] = max( limits[1], node->ends_[1] );
    }
}

NodeList Node::getNodeListIntersection( NodeList &a, NodeList &b )
{
    NodeList c( a.begin(), a.end() );
    for ( Node* node : b )
    {
        if ( find( c.begin(), c.end(), node ) == c.end() )
        {
            c.push_back( node );
        }
    }
    return c;
}

void Node::sortNodeListByFurthest( NodeList &nodes, bool endDrxn, bool drxn )
{
    sort( nodes.begin(), nodes.end(), [&]( Node* &a, Node* &b ){
        return ( drxn ? a->ends_[endDrxn] > b->ends_[endDrxn] : a->ends_[endDrxn] < b->ends_[endDrxn] );
    } );
}

