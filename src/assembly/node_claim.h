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

#ifndef NODE_CLAIM_H
#define NODE_CLAIM_H

#include "types.h"
#include "node_groups.h"
#include "node_structs.h"

class Node;
class ClaimNode;

struct ClaimEdge
{
    ClaimEdge( ClaimNode* cn, int ol, int32_t diff, bool isLeap ) : node( cn ), ol( ol ), diff( diff ), isLeap( isLeap ) {};
    ClaimNode* node;
    int ol, diff;
    bool isLeap;
};

struct ClaimJoin
{
    ClaimJoin( ClaimNode* fork, ClaimNode* branch, int32_t dist, bool drxn );
    ~ClaimJoin();
    int diff( bool drxn );
    void fill( Node* node, NodeDists& dists, bool drxn );
    ClaimNode* node[2];
    Nodes bridge;
    int dist;
};

struct ClaimRepair
{
    ClaimRepair( Node* l, Node* r );
    ClaimRepair( vector<Node*>& l, vector<Node*>& r );
    void dupe( vector< pair<ClaimNode*, ClaimNode*> >& cloned, NodeRoll& nodes );
    void fill( Node* node, int32_t off, bool drxn );
    unordered_map<Node*, int32_t> get( ClaimRepair& cr, bool drxn );
    void sever( Node* node, unordered_map<Node*, int32_t>& tar, bool drxn );
    void trim();
    vector<Node*> forks[2];
    unordered_map<Node*, int32_t> offs[2];
};

struct ClaimRedundant
{
    ClaimRedundant( vector<ClaimNode*> alts[2], Node* paired[2], int hits );
    ~ClaimRedundant();
    bool add( vector<ClaimNode*> alts[2], Node* paired[2], int hits );
    bool disregard( Node* paired[2] );
    int get( unordered_set<ClaimNode*>& lTar, unordered_set<ClaimNode*>& rTar, ClaimNode* l, ClaimNode* r, bool lEdge, bool rEdge );
    bool reach( ClaimNode* cn, ClaimNode* l, ClaimNode* r, vector<bool>& tars, unordered_set<ClaimNode*>& pathed, bool& reached, bool lEdge, bool rEdge );
    vector<ClaimNode*> nodes[2];
    vector< pair<Node*, Node*> > used;
    int score;
};

struct ClaimShared
{
    ClaimShared( ClaimNode* node ):node( node ){};
    ClaimNode* node;
    Nodes shared;
};

struct ClaimPairing
{
    ClaimPairing( ClaimNode* l, ClaimNode* r, int32_t dist );
    ~ClaimPairing();
    int hit( vector<int32_t>& hitDiffs );
    int get( ClaimNode* l, ClaimNode* r, bool lEdge, bool rEdge );
    int get( vector<ClaimNode*>& path, int mode );
    bool reach( ClaimNode* cn, vector<ClaimNode*>& path, int i, int mode, unordered_set<ClaimNode*>& pathed, int32_t diff, int& score, bool success, bool fail );
    bool reach( ClaimNode* cn, ClaimNode* l, ClaimNode* r, unordered_set<ClaimNode*>& pathed, int32_t diff, int& score, bool lEdge, bool rEdge );
    ClaimNode* node[2];
    vector<int32_t> diffs;
    vector<int> missed;
    int hits;
};

struct ClaimScore
{
    ClaimScore( Node* l, Node* r, ClaimPairing* cp, int32_t dist, int32_t est );
    void addRedundant( vector<ClaimNode*> alts[2], int hits );
    void cull( int32_t cutoff );
    bool redundant( vector<ClaimNode*>& claims, int32_t est, int32_t best, int32_t cutoff, int hits );
    void score( vector<ClaimNode*>& claims, int32_t est, int32_t cutoff, int hits );
    bool setAlts( int32_t est );
    bool setDistal( bool distal[2] );
    void setClaims( vector<ClaimNode*>& claims, bool drxn );
    Node* node[2];
    vector<ClaimPairing*> pairs;
    vector< vector<int32_t> > diffs;
};

struct ClaimBranch
{
    ClaimBranch( ClaimEdge& ce, bool shared, bool drxn );
    bool claim( unordered_set<ClaimBranch*> claims[2] );
    void fill( ClaimNode* cn, bool drxn );
    bool isShared( unordered_set<ClaimBranch*>& alts );
    ClaimEdge edge;
    vector<ClaimNode*> path;
    unordered_set<ClaimNode*> pathed;
    vector< pair<ClaimBranch*, int> > dumped[2];
    bool claimed, shared;
};

struct ClaimDupe
{
    static bool dupe( ClaimNode* fork, vector<ClaimNode*>& claims, NodeRoll& nodes );
private:
    ClaimDupe( ClaimNode* seed );
    ClaimDupe( ClaimBranch* seed, ClaimDupe* cd, bool drxn );
    ~ClaimDupe();
    bool advance( ClaimBranch* cb, vector<ClaimNode*>& claims, NodeRoll& nodes, bool& split, bool drxn );
    bool claim( unordered_set<ClaimBranch*> claimed[2], vector<ClaimNode*>& claims, NodeRoll& nodes );
    bool resolve( vector<ClaimNode*>& claims, NodeRoll& nodes );
    void setDettached( vector< pair<ClaimNode*, ClaimNode*> >& cloned, bool drxn );
    bool setDump( ClaimBranch* l, ClaimBranch* r );
    void setPath( vector< pair<ClaimNode*, ClaimNode*> >& cloned, vector<ClaimNode*>& claims, NodeRoll& nodes );
    int setScore( ClaimBranch* cd[2][2], unordered_set<ClaimNode*> alts[2], int l, int r );
    int setUniques( ClaimBranch* l, ClaimBranch* r );
    
    vector<ClaimNode*> path;
    unordered_set<ClaimNode*> edged[2];
    vector<ClaimBranch*> branches[2];
    ClaimNode* forks[2];
    ClaimDupe* edge[2];
    int gen;
};

struct ClaimTrim
{
    static bool trim( ClaimNode* fork, NodeRoll& nodes, bool drxn );
private:
    ClaimTrim( ClaimNode* l, ClaimNode* r, ClaimNode* bridge );
    void fill( ClaimNode* cn, unordered_set<ClaimNode*>& group, bool drxn );
    int score( int l, int r );
    void score( int l, int r, int& hits, int& unique );
    bool trim( int ol );
    ClaimNode* forks[2];
    vector<ClaimNode*> path;
    unordered_set<ClaimNode*> groups[2][2];
    int ols[2];
};

class ClaimNode
{
public:
    ClaimNode( ClaimNode* cn, vector<ClaimNode*>& claims, NodeRoll& nodes );
    ClaimNode( ClaimNode* cn, NodeRoll& cloned, Nodes& base );
    ~ClaimNode();
    void addEdge( ClaimNode* cn, int ol, int32_t diff, bool isLeap, bool drxn, bool reciprocate=true, bool nodeEdge=false );
    static bool claim( Node* fork, NodeRoll& nodes, bool drxn );
    void dismantle();
    int32_t* get( Node* node );
    ClaimEdge getEdge( ClaimNode* cn, bool drxn );
    Node* getFork( bool drxn );
//    vector<Node*> getForks( Node* fork, bool drxn );
    ClaimPairing* getPairing( ClaimNode* cn, bool drxn );
    bool isBridge();
    bool isForked( bool drxn );
    bool isInvalid( Node* node, int32_t off, ClaimNode* cn, bool drxn );
    bool isShared( Node* node, ClaimNode* cn, bool drxn );
    bool isValid( Node* node, bool drxn );
    bool reachViaEdgeOnly( ClaimNode* tar, ClaimNode* forks[2], unordered_set<ClaimNode*>& pathed, int32_t diff, vector<int32_t>& diffs, bool path[2], bool forked, bool drxn );
    bool reachViaForkOnly( ClaimNode* tar, ClaimNode* fork, unordered_set<ClaimNode*>& pathed, int32_t diff, vector<int32_t>& diffs, bool forked, bool drxn );
//    bool reach( ClaimNode* tar, int32_t diff, vector<int32_t>& diffs, bool drxn );
    bool removeEdge( ClaimNode* cn, bool drxn, bool reciprocate=true );
    void setBase( Nodes& base, bool distal );
    bool setState( bool drxn );
    
    ClaimNode* clone_;
    vector<Node*> path_;
    vector<ClaimEdge> edges_[2];
    vector<ClaimJoin*> joins_[2];
    vector<ClaimPairing*> pairs_[2];
    vector<ClaimRedundant*> redundant_[2];
    vector<ClaimShared> shared_;
    vector<ClaimNode*> alts_;
    Nodes distal_;
    bool edged_[2], split_;
    
private:
    ClaimNode( Node* node, vector<ClaimNode*>& claims, int32_t coord, bool orient, bool drxn );
    void addDistal( Node* node, bool drxn );
    bool addEdge( Edge& e, vector<ClaimNode*>& claims, bool orient, bool drxn );
    void addJoin( ClaimNode* cn, Node* hang, NodeDists& dists, int32_t dist, bool drxn );
//    void addJoin( ClaimNode* cn, int32_t diff, bool drxn, bool reciprocate=true );
    static bool branch( vector<ClaimNode*>& claims, vector<ClaimNode*>& forks, bool orient );
    bool branch( vector<ClaimNode*>& claims, unordered_set<ClaimNode*>& tested, bool branched, bool orient, bool drxn );
    static void complete( vector<ClaimNode*>& claims, bool drxn );
    bool create( vector<ClaimNode*>& claims, bool drxn );
    bool dupe( unordered_set<ClaimNode*> duped[2], vector<ClaimNode*>& claims, NodeRoll& nodes, int drxn );
    bool edge( vector<ClaimNode*>& claims, bool orient, bool drxn );
    void extend( vector<ClaimNode*>& claims, bool orient, bool drxn );
    void fill( Node* node, Nodes& block, int32_t dist, int32_t limit, bool orient, bool drxn );
    void fill( Node* node, Nodes& include, int32_t dist, bool orient, bool drxn );
    bool findFork( Node* q, bool drxn );
    int32_t getCoord( Node* fork, Edge& e, int32_t dist, bool orient, bool drxn );
    vector<Edge> getEdges( Node* fork, bool drxn, bool blunt=false );
    Nodes getExts();
//    vector<Node*> getForks( Node* fork );
    vector< pair<Node*, int32_t> > getHangs( bool orient, bool drxn );
    bool isBranched( bool drxn );
    bool isExtend( int32_t coord );
    bool reach( ClaimNode* tar, bool drxn );
//    bool removeJoin( ClaimNode* cn, bool drxn, bool reciprocate=true );
    void reset();
    static bool resplit( vector<ClaimNode*>& claims, bool orient );
    static void score( vector<ClaimNode*>& claims, bool orient );
    void setBlocked( Node* node, bool drxn );
//    void setBlocks( Nodes& block, bool drxn );
    void setBranched( vector<ClaimNode*>& claims, bool reversed, bool orient, bool drxn );
    bool setBranched( vector<ClaimNode*>& claims, unordered_set<ClaimNode*>& tested, bool branched, bool orient, bool drxn );
    void setDistal( bool drxn );
    static void setEnds( vector<ClaimNode*>& ends, Nodes& base, bool orient, bool drxn );
    bool setExts( Nodes& shared, int32_t limit, bool drxn );
    void setForked( vector<ClaimNode*>& claims, bool orient, bool drxn );
    void setForks( vector<ClaimNode*>& forks, bool drxn );
    static void setJoins( vector<ClaimNode*>& claims, Nodes& base, int32_t limit, bool orient, bool drxn );
    void setKeep( Node* node, Nodes& keep, bool drxn );
    void setPairs( ClaimNode* cn, unordered_set<ClaimNode*>& pathed, int32_t diff, bool lFork, bool rFork );
//    bool setRedundant( vector<ClaimNode*>& claims, ClaimPairing* cp, Node* paired[2], int32_t cutoff, int32_t est, int hits, bool orient );
    void setScores( vector<ClaimNode*>& claims, Nodes& base );
//    void setScores( vector<ClaimNode*>& claims, Nodes& base, bool orient );
    bool split( vector<ClaimNode*>& claims, bool orient, bool drxn );
    void testLoop( ClaimNode* cn );
    bool trim( unordered_set<ClaimNode*> trimmed[2], NodeRoll& nodes, int drxn );
    
    unordered_map<Node*, int32_t> offs_;
    int32_t ends_[2];
};

//struct Claim
//{
////    Claim( Node* fork, Edge& e, int32_t dist, bool branch, bool orient, bool drxn );
//    Claim( Node* fork, int32_t dist, bool orient, bool drxn );
//    Claim( Claim* fork, Edge& e, bool orient, bool drxn );
//    
//    void advance( Node* fork, bool orient, bool drxn );
//    static void block( vector<Claim*>& claims, bool drxn );
//    void block( Node* node, Nodesx& blocked, Nodesx& keep, bool drxn );
//    void cull( Nodesx& keep );
////    int32_t estimate( Claim* f, int32_t& best, int32_t est, bool drxn );
//    bool identical( Claim* rhs, bool drxn );
//    Node* place( int i, bool drxn );
//    bool redundant( Claim* rhs, bool drxn, bool reciprocate=true );
//    void retract();
//    vector<Claim*> split( bool init, bool orient, bool drxn );
//    static Nodesx target( vector<Claim*>& branches );
//    Node* terminus( bool drxn );
//    
//    vector<Node*> path;
//    vector<Edge> edges;
//    NodeOffsets offs;
//    int32_t coord;
//};
//
//struct ClaimPair
//{
//    ClaimPair( Claim* fwd, Edge e, int altCount );
//    int score();
//    Claim* path;
//    Edge edge;
//    int32_t work;
//    vector<int> alts, pairs;
//    int uniques;
//    bool drop;
//};
//
//struct ClaimBlocks
//{
//    ClaimBlocks( vector<Claim*> test[2] );
//    bool ambiguous( vector<Claim*> test[2], Node* bck, Node* fwd, bool drxn );
//    bool disregard( vector<Claim*> test[2], Node* bck, Node* fwd, bool drxn );
//    vector<Nodesx> blocks;
//    vector< pair<Nodesx, Nodesx> > used;
//};
//
//struct ClaimScores
//{
//    ClaimScores( Claim* path, vector<Claim*> test[2], bool drxn );
//    static void create( vector<Claim*> test[2], vector<ClaimScores*>& scores, bool orient, bool drxn );
//    bool contest();
//    static bool disregard( Node* bck, Node* fwd, Nodesx& used, bool drxn );
//    int get( ClaimScores* alt, Claim* branch );
////    int get( ClaimScores* alt, Claim* branch, int& ol );
//    bool paired( Claim* c, bool exists=false );
//    bool preferred( int hits[2][3], int ols[2], int altHits[2][3], int altOl, bool dual );
//    bool remove( Claim* branch );
//    static void score( vector<Claim*> test[2], vector<ClaimScores*>& scores, bool orient, bool drxn );
//    static bool score( vector<ClaimScores*>& scores, int32_t cutoff, bool drxn );
////    static bool score( vector<ClaimScores*>& scores, int32_t cutoff );
//    void score( Node* b, Node* f, int32_t est, int32_t& best, bool drxn );
//    void score( int hit, bool unique );
//    bool set( ClaimScores* alt, Claim* weak, Claim* strong, int hits[2][3], int ols[2], bool& paired );
//    void trim( bool drxn );
//    Claim* path;
//    vector<ClaimScores*> alts;
//    vector<ClaimPair> pairs;
//    bool claimed, dropped;
//};
//
//struct ClaimMap
//{
//    ClaimMap( ClaimScores* cs, Nodesx& bases, NodeRoll& cloned, bool drxn );
//    void claim( Claim* c, bool d );
//    bool confirm( bool drxn );
//    void connect( ClaimPair& cp, bool drxn );
////    void disconnect( bool drxn );
//    void disconnect( ClaimPair& cp, bool drxn );
//    Node* get( Node* node );
//    bool match( ClaimScores* cs, bool drxn );
//    void unmatch( Claim* path, bool drxn );
//    vector<ClaimScores*> scores[2];
//    Nodesx& bases, alts[2]; // alts[0] belong to different paths (must clone), alts[1] belong to the same paths, but shared by different pairs (must edge to)
//    NodeRoll& cloned, safe[2];
//    unordered_set<Claim*> claims[2];
//    unordered_map<Node*, Node*> used;
//    Node* fork;
//    bool claimed;
//};
//
//class ClaimFork
//{
//public:
//    ClaimFork( Node* fork, NodeRoll& nodes, bool drxn );
//    ClaimFork( Node* fork, Node* branch, bool drxn );
//    ClaimFork( Node* fork, vector<Claim*>& borrow, NodeRoll& nodes, bool orient, bool drxn );
//    ~ClaimFork();
//    
//    NodeRoll cloned_;
//    bool claimed_;
//private:
//    void add( Node* fork, Edge e, int32_t dist, bool drxn );
////    void add( Node* fork, Edge& e, int32_t dist, bool branch, bool drxn );
//    bool claim( NodeRoll& nodes );
//    bool complete( bool drxn );
//    bool confirm();
//    void contest( Claim* branch, vector<Claim*>& contested, bool drxn );
//    bool retry();
//    
//    Node* fork_;
//    vector<Claim*> paths_[2], added_;
////    vector<ClaimFork*> extends_[2];
//    vector<ClaimScores*> scores_;
//    bool orient_, drxn_;
//};

#endif /* NODE_CLAIM_H */

