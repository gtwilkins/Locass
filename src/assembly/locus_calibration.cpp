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
    
    NodeList nodes = getAllNodes();
    Node::calibrate( nodes, lib );
    
    vector< pair<float, int> > coverages;
    
    int totalLen = 0;
    if ( !paths_[0][0].path.empty() && !paths_[1][0].path.empty() )
    {
        NodeSet fwdSet = paths_[0][0].path.back()->getDrxnNodes( 1 );
        NodeSet bckSet = paths_[1][0].path.back()->getDrxnNodes( 0 );
        NodeSet alleleSet;
        for ( bool drxn : { 0, 1 } )
        {
            for ( Allele &allele : paths_[drxn][0].alleles )
            {
                alleleSet.insert( allele.paths[0].begin(), allele.paths[0].end() );
                alleleSet.insert( allele.paths[1].begin(), allele.paths[1].end() );
            }
        }
        
        for ( Node* node : fwdSet )
        {
            int nodeLen = max( 1, node->ends_[1] - node->ends_[0] - params.readLen );
            if ( alleleSet.find( node ) != alleleSet.end() )
            {
                coverages.push_back( make_pair( node->coverage_ * 2, nodeLen ) );
                totalLen += nodeLen;
            }
            else if ( bckSet.find( node ) != bckSet.end() )
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
}