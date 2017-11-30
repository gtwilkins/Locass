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

void Node::setFurthest( NodeSet &tSet, bool drxn )
{
    if ( !pairs_.empty() )
    {
        NodeList tTmp;
        for ( auto &np : pairs_ )
        {
            if ( tSet.find( np.first ) != tSet.end() )
            {
                tTmp.push_back( np.first );
            }
        }
        
        for ( ReadMark &mark : getMarksBase( drxn ) )
        {
            for ( Node* t : tTmp )
            {
                auto it = t->reads_.find( mark.readId );
                if ( it != t->reads_.end() )
                {
                    if ( !farPairNodes_[0] || ( drxn ? it->second[0] < farPairCoords_[0]
                                                     : farPairCoords_[0] < it->second[1] ) )
                    {
                        farPairNodes_[0] = t;
                        farPairCoords_[0] = it->second[!drxn];
                    }
                    if ( t->isReliable() && ( !farPairNodes_[1] || ( drxn ? it->second[0] < farPairCoords_[1]
                                                                          : farPairCoords_[1] < it->second[1] ) ) )
                    {
                        farPairNodes_[1] = t;
                        farPairCoords_[1] = it->second[!drxn];
                    }
                }
            }
        }
    }
}

void Node::resetFurthest( Node* node )
{
    if ( farPairNodes_[0] == node || farPairNodes_[1] == node )
    {
        farPairNodes_[0] = farPairNodes_[1] = NULL;
    }
}

void Node::resetFurthest( NodeSet &resetSet )
{
    if ( farPairNodes_[0] && resetSet.find( farPairNodes_[0] ) != resetSet.end() )
    {
        farPairNodes_[0] = farPairNodes_[1] = NULL;
    }
    else if ( farPairNodes_[1] && resetSet.find( farPairNodes_[1] ) != resetSet.end() )
    {
        farPairNodes_[1] = NULL;
    }
}
