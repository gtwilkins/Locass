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

#include "locus_path.h"
#include "leap.h"

//vector<LocusPath*> LocusPath::create( Querier& bwt, NodeRoll& nodes )
//{
//    Node::verify( nodes );
//    vector<LocusPath*> loci;
//    vector<NodesPath*> seeds, paths;
//    NodesPath::create( nodes, seeds, paths );
//    vector< vector<NodesPath*> > alleles;
//    for ( NodesPath* seed : seeds ) alleles.push_back( PathMapping::map( seed ) );
//    
//    unordered_set<NodesPath*> used, ends[2], branches[2];
//    for ( vector<NodesPath*>& allele : alleles )
//    {
//        for ( NodesPath* np : allele ) used.insert( np );
//        ends[0].insert( allele[0] );
//        ends[1].insert( allele.back() );
//    }
//    for ( int d : { 0, 1 } ) for ( NodesPath* np : used ) for ( PathEdge* pe : np->edges_[d] ) if ( used.find( pe->edge[d] ) == used.end() )
//    {
//        branches[d].insert( pe->edge[d] );
//    }
//    vector<Node*> forks[2];
//    Nodesx extable[2], leapable[2];
//    for ( int d : { 0, 1 } ) for ( NodesPath* np : branches[d] ) np->setBranch( used, extable[d], leapable[d], d );
//    for ( int d : { 0, 1 } ) for ( NodesPath* np : ends[d] ) forks[d].push_back( np->getEnd( bwt, nodes, used, d ) );
//    for ( int d : { 0, 1 } ) for ( Node* node : leapable[d].nodes ) Leap::leapBranch( bwt, nodes, node, d );
//    assert( false );
//    return loci;
//}
