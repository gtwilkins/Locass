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

#ifndef PATH_ALLELES_H
#define PATH_ALLELES_H

#include "node_path.h"

struct AlleleBlock
{
    AlleleBlock( vector<NodePath*>& path );
    AlleleBlock( vector<Node*> paths[2], bool drxn );
    AlleleBlock( vector<Node*> paths[2], vector< pair<int, int> >& matches, int& i, int& j, int& k );
    int getScore( int i, int& count, bool drxn );
    int getScore( int i, int& uncount, int& count, bool drxn );
    vector<Node*> path_[2];
    bool homo_;
};

struct AlleleNode
{
    AlleleNode( vector<NodePath*>& path );
    AlleleNode( vector<NodePath*> path[3], AlleleNode* base, bool drxn );
    ~AlleleNode();
//    AlleleNode( NodesPath* seed, vector<AlleleNode*>& path );
//    AlleleNode( NodesPath* fork, vector<AlleleNode*>& path, bool drxn );
//    AlleleNode( vector<NodesPath*>& merge, vector<AlleleNode*>& path, bool drxn );
//    void extend( vector<AlleleNode*>& path, bool drxn );
//    void fill();
    void match( vector<Node*> paths[2], vector< pair<int, int> >& cur, vector< pair<int, int> >& matches, int32_t& best, int i, int j );
    int getScore();
    int getScore( int i, int& j, bool drxn );
    int getScore( int i, int& b, int& f, bool drxn );
    vector<NodePath*> path_[2];
    vector<AlleleBlock*> blocks_;
    int32_t len[2], ol[2][2];
    bool homo_, bested_[2];
};

//struct AlleleBranch;
//
//struct AlleleMatch
//{
//    AlleleMatch( AlleleBranch* a, AlleleBranch* b );
//    bool advance( bool drxn );
//    void match( vector< pair<int, int> >& cur, int i, int j );
//    AlleleBranch* branch[2];
//    vector< pair<int, int> > matched;
//    int hits, miss;
//};
//
//struct AlleleBranch
//{
//    AlleleBranch( vector<NodesPath*>& path, int32_t ol, bool drxn );
//    AlleleBranch( AlleleBranch* fork, PathEdge* edge, bool drxn );
//    vector<AlleleBranch*> split( bool drxn );
//    vector<NodesPath*> path_;
//    vector< pair<Node*, int32_t> > dists_;
//    vector<AlleleMatch*> pairs_;
//};

struct AlleleGraph
{
    AlleleGraph( NodePath* branch, unordered_set<NodePath*> pathed[2], vector<AlleleGraph*> graphs[2], int i, bool drxn );
    void set();
    void setGaps( AlleleGraph* ag, int32_t dist, bool branched );
    vector<NodePath*> path;
    vector<Node*> nodes;
    vector< pair<int, int> > scores;
    vector< vector< pair<AlleleGraph*, int> > > marks;
    vector< pair<AlleleGraph*, int32_t> > edges[2], branches;
    unordered_map<AlleleGraph*, int32_t> gaps;
    AlleleGraph* paired;
    int32_t len, multi;
};

struct AlleleTar
{
    AlleleTar(): merged_( false ), fanned_( false ){};
    bool accept( int32_t est, int32_t dist );
    void add( NodePath* np, ReadId id, bool hit );
    bool confirmHit( PathPair* p, bool drxn );
    bool confirmMiss( PathPair* p, bool drxn );
    bool discard( NodePath* np, bool drxn );
    bool extend( bool drxn );
    void fill( NodePath* np, unordered_set<NodePath*>& ignores, bool complete, bool drxn );
    pair<int32_t, int>* get( NodePath* np );
    bool ignoreHit( NodePath* t, NodePath* q, int32_t est, int32_t dist, bool drxn );
    bool ignoreMiss( NodePath* t, NodePath* q, bool drxn );
//    bool ignore( NodesPath* np, NodesPath* t, int32_t limit, bool good, bool bad, bool drxn );
    void insert( NodePath* np, int32_t off, int multi );
    void match( bool drxn );
    void prepare( bool drxn );
    bool reach( NodePath* np, NodePath* q, unordered_set<NodePath*>& pathed, bool drxn );
    void reach( NodePath* np, vector<NodePath*>& path, unordered_set<NodePath*>& ignores, int32_t diff, int32_t limit, bool drxn );
    void secure( NodePath* np, NodePath* block, unordered_set<NodePath*>& keep, bool drxn );
    int setBranch( NodePath* np, bool drxn );
    int setForks( NodePath* np, vector<NodePath*>& branched, bool drxn );
    int setScore( NodePath* np, unordered_set<NodePath*>& pathed, unordered_set<ReadId>& hits, bool drxn );
    void survey( NodePath* np, NodePath* alt, unordered_map<NodePath*, int>& block, int32_t diff, int32_t limit, bool drxn );
    unordered_map< NodePath*, pair<int32_t, int> > tar_;
    unordered_map< NodePath*, int32_t> que_;
    unordered_map< NodePath*, int> loop_;
    unordered_map< NodePath*, unordered_set<ReadId> > hits_, miss_;
    unordered_map< NodePath*, unordered_map<NodePath*, int32_t> > ignores_;
    unordered_map< NodePath*, vector< pair< NodePath*, int> > > forks_;
    unordered_map< NodePath*, pair<NodePath*, int> > branches_;
    vector<NodePath*> path_;
    int32_t ends_[2];
    bool merged_, fanned_;
};

struct AlleleAlign
{
    AlleleAlign(): score( 0 ){};
    AlleleAlign( AlleleGraph* branch[2], bool best );
    bool add( AlleleGraph* a, AlleleGraph* b, int i, int j, bool drxn );
    bool advance( int i, AlleleGraph* b, int j );
    void append( AlleleGraph* ag, int i );
    AlleleAlign branch( AlleleGraph* b, int j );
    void challenge( AlleleAlign& challenger );
    bool join( AlleleGraph* ag, AlleleGraph* tar, unordered_set<AlleleGraph*>& joins, unordered_set<AlleleGraph*>& pathed );
    int32_t gap( AlleleGraph* base, int i, AlleleGraph* ext, int j, bool drxn );
    int multi( AlleleGraph* ag, int i, bool paired );
    bool output( vector<NodePath*> ext[3], bool drxn );
    void rescore();
    bool terminate( AlleleGraph* a, AlleleGraph* b, int i, int j );
    void trim( int i );
    vector< pair<AlleleGraph*, int> > align[2];
    vector<AlleleGraph*> path[2];
    vector<int32_t> gaps[2], lens;
    int32_t score;
};

class AllelePaths
{
    AllelePaths( NodePath* seed, unordered_set<NodePath*>& used );
    ~AllelePaths();
    void align( vector<NodePath*> path[2], int32_t len[2], int32_t ol[2], bool drxn );
    void align( bool drxn );
    void challenge( AllelePaths* ap );
    bool discard( AllelePaths* ap, unordered_map<NodePath*, int>& multis );
    bool dump();
    bool extend( bool drxn );
    unordered_map<NodePath*, int> getMultis();
    bool merge( vector<AllelePaths*>& alleles, bool drxn );
    void plot( AlleleGraph* a, int i, AlleleAlign& aligned, AlleleAlign& best, bool drxn );
    static void reduce( vector<AllelePaths*>& alleles );
    void setAligns( unordered_set<NodePath*> pathed[2], vector<AlleleGraph*> graphs[2], bool drxn );
    void setBranches( vector<NodePath*> branches[2][2], int len[2] );
    void setClones( vector<NodePath*>& paths, NodeRoll& nodes );
    void setDiffs();
    void setDiffs( NodePath* np, int32_t diff, unordered_set<NodePath*>& path, unordered_set<NodePath*>& pathed, int i, bool drxn );
//    bool setChallenged( AllelePaths* ap, vector<NodesPath*>& path, vector<AlleleNode*>& aligned, bool complete[2] );
//    bool setMerged( AlleleAlign& aligned, bool drxn );
    bool setMerged( AlleleTar tar[2], bool drxn );
    void setMultis( unordered_map<NodePath*, int>& multis );
    void setPathable( NodePath* np, unordered_set<NodePath*>& pathed, unordered_set<NodePath*>& merge, unordered_set<NodePath*>& base, int32_t limit, bool drxn );
    void setQuery( AlleleTar tar[2], int32_t len[2], int32_t diff[2], bool drxn );
    void setTarget( AlleleTar& tar, bool homo, bool drxn );
    void setTarget( NodePath* np, vector<NodePath*>& pathed, unordered_map<NodePath*, int32_t>& t, int32_t dist, int32_t limit, bool drxn );
    
    vector<AlleleNode*> path_;
    unordered_map<NodePath*, int32_t> diffs_[2];
    unordered_map<ReadId, int32_t> reads_[2];
    unordered_set<ReadId> unique_;
    int id_, score_;
    bool ambiguous_[2], contested_;
public:
    static void create( vector<NodePath*>& paths, NodeRoll& nodes );
};



#endif /* PATH_ALLELES_H */

