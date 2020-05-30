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

#ifndef PATH_CROSS_H
#define PATH_CROSS_H

#include "node_path.h"
#include "path_alleles.h"

struct PathCross;
struct CrossBranch;

struct CrossPair
{
    CrossPair( CrossBranch* l, CrossBranch* r );
    ~CrossPair();
    CrossBranch* branch[2];
    unordered_set<ReadId> unique_, hits_;
    int score_, miss_;
};

struct CrossBranch
{
    CrossBranch( CrossBranch* fork );
    CrossBranch( NodePath* fork, NodePath* block, int32_t off, bool drxn );
    void fill( NodePath* np, unordered_set<NodePath*>& used, int32_t diff, int32_t limit, bool drxn );
    int32_t* get( NodePath* np );
    vector<PathEdge*> getEdges( bool drxn );
    CrossBranch* fork_;
    vector<NodePath*> path_;
    unordered_map<NodePath*, int32_t> diffs_;
    vector<CrossPair*> pairs_;
    unordered_set<ReadId> hits_, miss_, unique_;
    int32_t ends_[2], loop_;
    float cover;
};

struct PathCross
{
    PathCross( vector<NodePath*>& cross, int32_t off, bool crossed );
    ~PathCross();
    bool claim( vector<NodePath*>& paths, NodeRoll& nodes );
    static vector<PathCross*> create( vector<NodePath*>& paths );
    void extend( CrossBranch* cb, NodePath* branch, int32_t diff, bool split, bool drxn );
    void match( vector< pair<NodePath*, int32_t> > tars[2], ReadId id, int32_t dist );
    bool reach( NodePath* np, unordered_map<NodePath*, int32_t>& pathed, vector< pair<NodePath*, int32_t> >& tars, int32_t diff, int32_t cutoff );
    static bool resolve( vector<NodePath*>& paths, NodeRoll& nodes );
    bool separate( CrossBranch* paired[2], vector<NodePath*>& paths, NodeRoll& nodes );
    vector<NodePath*> cross_;
    vector<CrossBranch*> branch_[2], forks_[2];
    int32_t ends_[2];
    bool crossed_;
};


#endif /* PATH_CROSS_H */

