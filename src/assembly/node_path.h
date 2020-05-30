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
#include "seed_fork.h"
#include "query.h"

class NodePath;

struct PathEdge
{
    PathEdge( NodePath* fork, NodePath* branch, Edge& e, int32_t diff, bool drxn );
    PathEdge( NodePath* l, NodePath* r, PathEdge* pe );
    PathEdge( NodePath* np, Node* ne, PathEdge* pe, bool drxn );
    ~PathEdge();
    void claim( NodePath* np, bool drxn );
    void downgrade( int drxn );
    bool isAlt();
    static bool sever( NodePath* l, NodePath* r );
    bool sever();
    static void sort( vector<PathEdge*>& edges );
    NodePath* edge[2];
    Node* node[2];
    int32_t diff;
    int ol, score, multi, shared;
    bool leap, bad;
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
    bool steal( NodePath* np, bool drxn );
    void test();
    NodePath* node_[2];
    vector<PathUnique> uniques_;
    vector<PathPair*> shared_;
    unordered_set<ReadId> ids_, mates_;
    int32_t diff_;
    int score_;
    bool forked_;
};

struct Branch
{
    Branch( PathEdge* pe, bool drxn );
    ~Branch();
    static vector<Branch*> create( NodePath* fork, bool drxn );
    bool cull();
    bool cull( Branch* b );
    int getMax( bool inclBase=false );
    static bool match( vector<NodePath*>& path, vector<Branch*> branches[2] );
    static bool match( vector<NodePath*>& path, vector<Branch*>& l, vector<Branch*>& r );
    bool matched( vector<PathEdge*> taken[2][2], bool drxn );
    static void setBase( vector<Branch*>& branches, vector<NodePath*>& path, bool drxn );
    static void setUsed( vector<Branch*>& branches, vector< vector<NodePath*> >& loci, bool drxn  );
    static void sort( vector<Branch*>& branches );
    void sort();
    PathEdge* edge_;
    vector< pair<Branch*, int> > matched_;
    Branch* match_;
    unordered_set<NodePath*> branch_;
    int base_, hits_, miss_;
    bool used_;
};

class NodePath
{
public:
    NodePath( NodePath* np, vector<NodePath*>& paths, NodeRoll& nodes, int multi );
    ~NodePath();
    static void create( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& seeds, vector<NodePath*>& paths );
    static void cull( vector<NodePath*>& seed, vector<NodePath*>& paths, SeedExtend exts[2], NodeRoll& nodes );
    static void cull( vector<NodePath*>& seed, vector<NodePath*>& paths, NodeRoll& nodes );
    void cull( vector<NodePath*>& paths, unordered_set<NodePath*>& pathed, unordered_set<NodePath*>& base, SeedExtend& ext, NodeRoll& nodes, bool drxn );
    void cull( vector<NodePath*>& paths, unordered_set<NodePath*>& pathed, NodeRoll& nodes, bool drxn );
    static void finalise( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& seeds, vector<NodePath*>& paths );
    static void forceCut( vector<NodePath*>& paths, NodeRoll& nodes, int i, int j, bool drxn );
    int32_t* get( Node* node );
//    int32_t getCoord( bool drxn );
    static int32_t getLen( vector<NodePath*>& path );
    int32_t getLen( NodePath* rhs, int32_t diff, bool drxn );
    int32_t getOffset( Node* node, bool drxn );
    int32_t getOffset( NodePath* np, bool drxn );
    int32_t getOverlap( NodePath* np, bool drxn );
    Node* getEnd( Querier& bwt, NodeRoll& nodes, unordered_set<NodePath*>& used, bool drxn );
    Node* getTerminus( bool drxn );
    bool isBranchable( bool drxn );
    static void naturalise( vector<NodePath*>& paths );
    static void print( vector<NodePath*>& paths, vector< vector<NodePath*> >& loci );
    static void print( vector<NodePath*>& paths );
    void print( int32_t minCoord, ofstream& ofs );
    void print( int32_t minCoord );
    void printUnpaired( vector<NodePath*>& paths );
    static void reset( vector<NodePath*>& paths );
    static vector< vector<NodePath*> > resolve( Querier& bwt, vector<NodePath*>& paths, NodeRoll& nodes );
    static vector< vector<NodePath*> > resolve( Querier& bwt, vector<NodePath*>& paths, vector<NodePath*>& seeds, NodeRoll& nodes );
    void setBranch( unordered_set<NodePath*>& branch, int32_t dist, bool drxn );
    void setBranch( unordered_set<NodePath*>& used, Nodes& extable, Nodes& leapable, bool drxn );
    bool setBridges( bool drxn );
    bool setBridge( int base, int score, bool drxn );
    void setEnds( unordered_set<NodePath*>& pathed, vector<NodePath*>& ends, bool drxn );
    static void setExtends( vector< vector<NodePath*> >& loci, SeedExtend seed[2] );
    void setMates( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& paths );
    void setReachable( unordered_map<NodePath*, vector<int32_t> >& reached, int32_t diff, int32_t limit, bool drxn );
    void setReachable( unordered_map<NodePath*, int32_t>& reached, int32_t dist, int32_t limit, bool drxn );
    void setReachable( vector< pair<NodePath*, int32_t> >& reached, unordered_set<NodePath*>& block, int32_t diff, int32_t limit, bool drxn, bool ignore=false );
    void setReachable( unordered_set<NodePath*>& reached, bool drxn );
    void setReachable( unordered_set<NodePath*>& reached, unordered_set<NodePath*>& blocked, int32_t limit, bool drxn, bool ignore=false );
    static void setReads( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& paths );
    void setTarget( unordered_set<NodePath*>& reached, int32_t dist, int32_t limit, bool drxn, bool ignore=false );
//    static bool setUnbad( vector<NodePath*>& paths, NodeRoll& nodes);
    int32_t size();
    void sweep( NodeRoll& nodes, bool drxn );
    
    vector<PathEdge*> edges_[2], alts_[2];
    vector<PathPairs*> paired_[2];
    unordered_map<Node*, int32_t> offs_;
    unordered_map<NodePath*, unordered_set<ReadId> > tmp_;
    vector<Node*> path_;
    string seq_;
    int32_t coords_[2], ends_[2], loop_;
    int id_, multi_, score_;
    float cover_;
    bool verified_, ended_[2], bridge_, blank_;
    
private:
    NodePath( Node* seed, vector<NodePath*>& paths, int32_t coord );
    NodePath( NodePath* np, vector<NodePath*>& paths, NodeRoll& nodes );
//    NodePath( NodePath* np, PathEdge* l, PathEdge* r, NodeRoll& nodes );
    
    void addEdge( Edge& e, vector<NodePath*>& paths, bool branch, bool drxn );
    bool addPair( NodePath* np, int32_t diff );
    void clonePairs( NodePath* base );
    static void cleanPairs( vector<NodePath*>& paths );
    void destroy( NodeRoll& nodes );
    bool doesReach( unordered_set<NodePath*>& include, int32_t coord, bool drxn );
    void extend( vector<NodePath*>& paths, bool branch, bool drxn );
    bool findFork( Node* node, bool drxn );
    bool findNode( Node* node );
    int32_t getCoord( Edge& e, bool drxn );
    int32_t getCoord( Node* node, Edge& e, int32_t dist, bool drxn );
    static float getCoverage( vector<NodePath*>& path, int minLen=0, int drxn=1 );
    int getMulti( bool drxn );
    void getOffsets( unordered_map<Node*, int32_t>& offs, unordered_set<NodePath*>& include, int32_t base, int32_t diff, int32_t limit, bool drxn );
    static string getSeq( vector<NodePath*>& path );
    bool isContinue( int32_t dist, bool drxn );
    bool isEnding( PathEdge* edge, float cover, bool drxn );
    bool isEnding( int dist, bool drxn );
    bool isForked( PathEdge* base, bool drxn );
    bool isPathed( Node* node );
    static bool merge( vector<NodePath*>& paths );
    bool merge();
    void reduce( bool drxn );
    void resetDiffs( NodePath* np, vector<NodePath*>& path, unordered_map<NodePath*, int32_t>& diffs, int32_t diff, int32_t limit );
    static void setBlanks( vector<NodePath*>& paths );
    void setBlank( unordered_map<NodePath*, float>& pathed, float cover, bool blank );
//    void setBlanks( unordered_set<NodePath*>& blanks, unordered_set<NodePath*>& good, float cover, bool blank );
    static bool setBlunts( Querier& bwt, vector<NodePath*>& paths, NodeRoll& nodes );
    bool setBlunt( Querier& bwt, bool drxn );
    static void setCloned( vector<NodePath*>& path, vector<PathEdge*> taken[2][2], vector<NodePath*>& paths, NodeRoll& nodes );
    static bool setCrosses( vector<NodePath*>& paths, NodeRoll& nodes );
    bool setCross( vector<Branch*> branches[2], vector<NodePath*>& paths, NodeRoll& nodes );
    NodePath* setCross( vector<Branch*> branches[2], NodeRoll& nodes );
    static bool setCulled( vector< vector<NodePath*> >& loci, vector<NodePath*>& paths, NodeRoll& nodes );
    bool setCulled( vector<NodePath*>& pathed, unordered_set<NodePath*> goods[2], int32_t dist, float base, bool forked, bool drxn );
    bool setCulled( vector<NodePath*>& pathed, int32_t len, float cutoff, bool drxn );
    static void setEdges( vector<NodePath*>& paths );
    void setEdges( NodePath* np, int32_t limit, vector<PathEdge*>& path, vector< unordered_set<PathEdge*> >& edged, unordered_set<PathEdge*>& used );
    static void setExclusive( vector< unordered_set<NodePath*> >& sets );
    static void setFills( vector<NodePath*>& paths );
    void setFill( Node* node, Nodes& block, int32_t dist, int32_t limit, bool drxn );
    void setGood( unordered_set<NodePath*>& goods, int32_t len, bool drxn );
    void setMulti();
    void setPairs( NodePath* np, vector<NodePath*>& path, int32_t diff, int32_t limit, bool lFork, bool rFork );
    static vector< vector<NodePath*> > setPaths( vector<NodePath*>& paths );
    vector<NodePath*> setPath( vector< vector<NodePath*> >& loci );
//    void setPairs( NodePath* np, vector<NodePath*>& path, int32_t diff, bool forked );
//    void setReads( bool drxn );
    static bool setUnforked( vector<NodePath*>& paths, NodeRoll& nodes );
    bool setUnforked( bool drxn );
    void setScore();
    void setScores();
    bool unfork( bool drxn );
    
};

class Haplo;
struct HapScore;

struct HapBranch
{
    HapBranch( NodePath* branch, int32_t dist, bool drxn );
    int32_t* get( NodePath* np );
    NodePath* branch_;
    unordered_map<NodePath*, int32_t> dists_;
};

struct HapFork
{
    HapFork( NodePath* fork, Haplo* base, vector<HapFork*>& forks, bool drxn );
    ~HapFork();
    vector<HapScore*> cross( HapFork* alt, bool drxn );
    bool cull( vector<HapFork*>& forks, bool drxn );
    bool resolve( vector<HapFork*> forks[2], bool drxn );
    bool resolveBranch( vector<HapFork*>& forks, bool drxn );
    NodePath* fork_;
    vector< pair<HapFork*, vector<HapScore*> > > cross_;
    vector<HapFork*> bck_;
    vector<HapBranch*> branches_;
    vector<Haplo*> base_;
    unordered_set<NodePath*> reached_;
};

struct HapScore
{
    HapScore( Haplo* h, HapBranch* branch, bool drxn );
    HapScore( HapBranch* l, HapBranch* r );
//    HapScore( HapBranch* l, bool drxn );
    unordered_set<ReadId> ids_, mates_, pref_;
    HapBranch* branch_[2];
    int unique_, miss_;
};

class Haplo
{
    Haplo( NodePath* seed );
    void add( NodePath* np, bool drxn );
public:
    void setBranches( bool drxn );
    void setFork( NodePath* fork, vector<HapFork*>& forks, bool drxn );
    void unfork( bool drxn );
    static vector<Haplo*> create( Querier& bwt, vector<NodePath*>& paths, vector<NodePath*>& seeds, NodeRoll& nodes );
    vector<NodePath*> path_;
    HapFork* fork_[2];
    vector<HapScore*> branch_[2];
};

//struct PathEdgeScore
//{
//    PathEdgeScore( NodePath* np, int32_t diff, int32_t limit, bool drxn );
////    bool add( PathPairing* pp, vector< pair<NodesPath*, int32_t> >& tar, vector<int32_t>& diffs, bool missed, bool drxn );
//    void add( NodePath* q, vector<int32_t>& qDiffs, unordered_set<NodePath*>& base, vector< pair<NodePath*, int32_t> >& tar, bool drxn );
//    vector<int32_t>* get( NodePath* np );
//    void getKeep( unordered_set<NodePath*>& keep, int32_t limit, bool drxn );
//    void setKeep( unordered_set<NodePath*>& keep );
//    NodePath* node;
//    unordered_map<NodePath*, vector<int32_t> > fwd;
//    int32_t diff;
//    int hits, miss, unique;
//};

//class PathMapping
//{
//public:
//    static vector<NodePath*> map( NodePath* seed );
//private:
//    PathMapping( vector<NodePath*>& path );
//    void add( NodePath* np, int32_t diff, bool drxn );
//    bool advance();
//    bool crossroad();
//    void edge( bool drxn );
//    NodePath* getFork( bool drxn );
//    vector<PathEdgeScore> exts_[2];
//    vector<NodePath*>& path_;
//    vector< pair<NodePath*, int32_t> > diffs_[2], alts_[2];
//    int32_t diff_;
//};


#endif /* NODE_PATH_H */

