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
#include <algorithm>

void Locus::calibrate( LocusLibraryCount &lib )
{
    calibrate_ = true;
    for ( bool drxn : { 0, 1 } )
    {
        endNodes_[drxn] = originEnds_[drxn];
        setExtend( drxn );
    }
    
    while ( ( canExtend( 0 ) || canExtend( 1 ) ) )
    {
        for ( int drxn : { 0, 1 } )
        {
            while ( canExtend( drxn ) )
            {
                extendNodes( drxn );
                if ( !updateExtension( drxn ) ) break;
            }
        }
        plot();
    }
    
    finalise();
    for ( int i : { 0, 1 } )
    {
        assert( paths_[i][0].path[0]->drxn_ == 2 );
    }
    
    NodeList nodes = getAllNodes();
    Node::calibrate( nodes, lib );
    for ( int i : { 0, 1 } )
    {
        assert( paths_[i][0].path[0]->drxn_ == 2 );
    }
    
    vector< pair<float, int> > coverages;
    
    int totalLen = 0;
    if ( !paths_[0][0].path.empty() && !paths_[1][0].path.empty() )
    {
        NodeSet alleleSet;
        NodeSet pathSet;
        for ( int i : { 0, 1 } )
        {
            assert( paths_[i][0].path[0]->drxn_ == 2 );
            pathSet.insert( paths_[i][0].path.begin(), paths_[i][0].path.end() );
            pathSet.insert( paths_[i][0].alleleSet.begin(), paths_[i][0].alleleSet.end() );
            alleleSet.insert( paths_[i][0].alleleSet.begin(), paths_[i][0].alleleSet.end() );
        }
        
        for ( Node* node : pathSet )
        {
            int nodeLen = max( 1, node->ends_[1] - node->ends_[0] - params.readLen );
            if ( alleleSet.find( node ) != alleleSet.end() )
            {
                nodeLen = max( 1, nodeLen / 2 );
                coverages.push_back( make_pair( node->coverage_ * 2, nodeLen ) );
                totalLen += nodeLen;
            }
            else
            {
                coverages.push_back( make_pair( node->coverage_, nodeLen ) );
                totalLen += nodeLen;
            }
        }
    }
    
    sort( coverages.begin(), coverages.end(), []( pair<float, int> &a, pair<float, int> &b ){
        return a.first < b.first;
    });
    
    lib.len = 0;
    int cutoffLen = totalLen * 0.67;
    for ( int i ( 0 ); i < coverages.size(); i++ )
    {
        if ( cutoffLen < coverages[i].second )
        {
            lib.len = totalLen;
            lib.coverage = coverages[i].first;
            break;
        }
        cutoffLen -= coverages[i].second;
    }
    
    calibrate_ = false;
}