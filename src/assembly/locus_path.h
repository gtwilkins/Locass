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


#ifndef LOCUS_PATH_H
#define LOCUS_PATH_H

#include "node.h"
#include "node_path.h"

class LocusPath;

struct LocusEdge
{
    LocusPath* edge;
    PathEdge* pe;
};

class LocusPath
{
    LocusPath( PathEdge* edge, unordered_set<NodePath*> forks[2], unordered_set<PathEdge*>& edged );
    bool extend( unordered_set<NodePath*> forks[2], unordered_set<PathEdge*>& edged, bool drxn );
    static bool isBest( PathEdge* edge );
    static bool isBridge( NodePath* np );
    bool isBridge( unordered_set<NodePath*> forks[2], unordered_set<NodePath*>& bridges );
    static bool isFork( NodePath* np, unordered_set<NodePath*>& bridges, bool drxn );
    bool isRedundant( LocusPath* t );
    static bool setBranches( vector<LocusPath*>& paths, unordered_set<NodePath*> forks[2], unordered_set<PathEdge*>& edged );
    static bool setBridges( vector<LocusPath*>& paths, unordered_set<NodePath*> forks[2], unordered_set<NodePath*>& bridges, unordered_set<PathEdge*>& edged );
    static bool setRedundant( vector<LocusPath*>& paths );
    vector<NodePath*> path_;
    bool blocked_[2], fork_[2], branch_[2];
    vector<LocusEdge*> edges_[2], bridges_[2], alts_[2];
    int score_, len_;
    vector<Node*> nodes_;
    bool bridge_;
    
public:
    static vector<LocusPath*> create( vector<NodePath*>& paths );
};


class Diploid
{
    void easyExtend( bool drxn );
    vector<NodePath*> getEdges( bool i, bool drxn );
    vector<NodePath*> alleles_[2];
    bool homo_;
public:
    Diploid( NodePath* a, NodePath* b );
    static void create( NodePath* a, NodePath* b, vector<NodePath*> alleles[2] );
};

#endif /* LOCUS_PATH_H */

