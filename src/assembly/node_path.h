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

#ifndef NODE_PATH_H
#define NODE_PATH_H

#include "node.h"
#include "node_structs.h"
#include "query.h"

class NodePath;

struct PathPairing
{
    PathPairing( NodePath* l, NodePath* r, int32_t diff );
    ~PathPairing();
    int hits( vector<int32_t>& hitDiffs, int32_t diff );
    NodePath* node[2];
    vector<int32_t> diffs;
    vector<int> missed;
    int score;
};

struct PathEdge
{
    PathEdge( NodePath* fork, NodePath* branch, Edge& e, int32_t diff, bool drxn );
    PathEdge( NodePath* l, NodePath* r, PathEdge* pe );
    ~PathEdge();
    void claim( NodePath* np, bool drxn );
    void downgrade();
    static bool sever( NodePath* l, NodePath* r );
    NodePath* edge[2];
    int32_t diff;
    int ol, score, multi;
    bool leap;
};

struct PathUnique
{
    PathUnique( ReadId id, int32_t l, int32_t r, int32_t dist, int32_t best ): id_( id ), l_( l ), r_( r ), dist_( dist ), best_( best ){};
    ReadId id_;
    int32_t l_, r_, dist_, best_;
};

struct PathPair
{
    PathPair( ReadId id, int32_t dist ) : id_( id ), dist_( dist ){};
    bool add( Node* base, NodeRoll* tar, unordered_map<Node*, NodePath*>& mapped, ReadId id, bool d );
    bool confirm( bool cull );
//    vector<PathPairs*> get();
    ReadId id_;
    vector< pair<NodePath*, int32_t> > marks_[2];
    int32_t dist_, best_;
};

struct PathPairs
{
    PathPairs( NodePath* l, NodePath* r, int32_t diff, bool forked );
    PathPairs( PathPairs* pp, NodePath* np, bool drxn );
    ~PathPairs();
    static bool add( NodePath* l, NodePath* r, int32_t diff, bool forked );
    void discard();
    NodePath* node_[2];
    vector<PathUnique> uniques_;
    vector<PathPair*> shared_;
    int32_t diff_;
    int score_;
    bool forked_;
};

class NodePath
{
public:
    NodePath( NodePath* np, vector<NodePath*>& paths, NodeRoll& nodes, int multi );
    ~NodePath();
    static void create( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& seeds, vector<NodePath*>& paths );
    int32_t* get( Node* node );
//    int32_t getCoord( bool drxn );
    int32_t getLen( NodePath* rhs, int32_t diff, bool drxn );
    int32_t getOffset( NodePath* np, bool drxn );
    int32_t getOverlap( NodePath* np, bool drxn );
    Node* getEnd( Querier& bwt, NodeRoll& nodes, unordered_set<NodePath*>& used, bool drxn );
    bool isBranchable( bool drxn );
    static void print( vector<NodePath*>& paths );
    static void reset( vector<NodePath*>& paths );
    void setBranch( unordered_set<NodePath*>& branch, int32_t dist, bool drxn );
    void setBranch( unordered_set<NodePath*>& used, Nodes& extable, Nodes& leapable, bool drxn );
    void setEnds( unordered_set<NodePath*>& pathed, vector<NodePath*>& ends, bool drxn );
//    void setMates( bool drxn );
    void setReachable( unordered_map<NodePath*, vector<int32_t> >& reached, int32_t diff, int32_t limit, bool drxn );
    void setReachable( vector< pair<NodePath*, int32_t> >& reached, unordered_set<NodePath*>& block, int32_t diff, int32_t limit, bool drxn, bool ignore=false );
    void setReachable( unordered_set<NodePath*>& reached, bool drxn );
    void setReachable( unordered_set<NodePath*>& reached, unordered_set<NodePath*>& blocked, int32_t limit, bool drxn, bool ignore=false );
    static void setReads( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& paths );
    void setTarget( unordered_set<NodePath*>& reached, int32_t dist, int32_t limit, bool drxn, bool ignore=false );
    int32_t size();
    
    vector<PathEdge*> edges_[2], breaks_[2];
    vector<PathPairing*> pairs_[2];
    vector<PathPairs*> paired_[2];
    unordered_map<Node*, int32_t> offs_;
    unordered_map<ReadId, int32_t> mp_[2];
//    unordered_map<ReadId, CrossRead> pe_[2], mp_[2];
    vector<Node*> path_;
    int32_t coords_[2], ends_[2], loop_;
    int id_, multi_;
    bool verified_, ended_[2];
    
private:
    NodePath( Node* seed, vector<NodePath*>& paths, int32_t coord );
    
    void addEdge( Edge& e, vector<NodePath*>& paths, bool branch, bool drxn );
    bool addPair( NodePath* np, int32_t diff );
    static void cleanPairs( vector<NodePath*>& paths );
    bool doesReach( unordered_set<NodePath*>& include, int32_t coord, bool drxn );
    void extend( vector<NodePath*>& paths, bool branch, bool drxn );
    void fill( Node* node, int32_t dist, int32_t limit, bool drxn );
    bool findFork( Node* node, bool drxn );
    bool findNode( Node* node );
    int32_t getCoord( Edge& e, bool drxn );
    int32_t getCoord( Node* node, Edge& e, int32_t dist, bool drxn );
    int getMulti( bool drxn );
    void getOffsets( unordered_map<Node*, int32_t>& offs, unordered_set<NodePath*>& include, int32_t base, int32_t diff, int32_t limit, bool drxn );
    void reduce( bool drxn );
    void resetDiffs( NodePath* np, vector<NodePath*>& path, unordered_map<NodePath*, int32_t>& diffs, int32_t diff, int32_t limit );
    void setMulti();
    void setPairs( NodePath* np, vector<NodePath*>& path, int32_t diff, int32_t limit, bool lFork, bool rFork );
    void setPairs( NodePath* np, vector<NodePath*>& path, int32_t diff, bool forked );
//    void setReads( bool drxn );
    void setScore();
    void setScores();
    
};

struct PathEdgeScore
{
    PathEdgeScore( NodePath* np, int32_t diff, int32_t limit, bool drxn );
//    bool add( PathPairing* pp, vector< pair<NodesPath*, int32_t> >& tar, vector<int32_t>& diffs, bool missed, bool drxn );
    void add( NodePath* q, vector<int32_t>& qDiffs, unordered_set<NodePath*>& base, vector< pair<NodePath*, int32_t> >& tar, bool drxn );
    vector<int32_t>* get( NodePath* np );
    void getKeep( unordered_set<NodePath*>& keep, int32_t limit, bool drxn );
    void setKeep( unordered_set<NodePath*>& keep );
    NodePath* node;
    unordered_map<NodePath*, vector<int32_t> > fwd;
    int32_t diff;
    int hits, miss, unique;
};

class PathMapping
{
public:
    static vector<NodePath*> map( NodePath* seed );
private:
    PathMapping( vector<NodePath*>& path );
    void add( NodePath* np, int32_t diff, bool drxn );
    bool advance();
    bool crossroad();
    void edge( bool drxn );
    NodePath* getFork( bool drxn );
    vector<PathEdgeScore> exts_[2];
    vector<NodePath*>& path_;
    vector< pair<NodePath*, int32_t> > diffs_[2], alts_[2];
    int32_t diff_;
};


#endif /* NODE_PATH_H */

