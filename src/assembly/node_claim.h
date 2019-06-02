/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   node_claim.h
 * Author: glen
 *
 * Created on 23 January 2019, 9:10 PM
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
    void fill( Node* node, NodeOffsets& offs, bool drxn );
    ClaimNode* node[2];
    Nodes bridge;
    int dist;
};

struct ClaimRepair
{
    ClaimRepair( Node* l, Node* r );
    void dupe();
    void fill( Node* node, int32_t off, bool drxn );
    unordered_map<Node*, int32_t> get( ClaimRepair& cr, bool drxn );
    void trim();
    Node* fork[2];
    unordered_map<Node*, int32_t> offs[2];
};

struct ClaimRedundant
{
    ClaimRedundant( vector<ClaimNode*> alts[2], Node* paired[2], int hits );
    ~ClaimRedundant();
    bool add( vector<ClaimNode*> alts[2], Node* paired[2], int hits );
    bool disregard( Node* paired[2] );
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
    ClaimNode* node[2];
    vector<int32_t> diffs;
    vector<int> missed;
    int paths, score;
};

struct ClaimScore
{
    ClaimScore( Node* l, Node* r, ClaimPairing* cp, int32_t dist, int32_t est );
    void addRedundant( vector<ClaimNode*> alts[2], int hits );
    bool cull( int32_t cutoff );
    bool redundant( vector<ClaimNode*>& claims, int32_t est, int32_t best, int32_t cutoff, int hits );
    void score( vector<ClaimNode*>& claims, int32_t est, int32_t cutoff, int hits );
    bool setAlts( int32_t est );
    bool setDistal( bool distal[2] );
    void setClaims( vector<ClaimNode*>& claims, bool drxn );
    Node* node[2];
    vector<ClaimPairing*> pairs;
    vector< vector<int32_t> > diffs;
};

struct ClaimDupe
{
    static bool dupe( ClaimNode* fork, NodeRoll& nodes, bool drxn );
private:
    ClaimDupe( ClaimNode* seed, vector<ClaimDupe*>& paths, bool drxn );
    ClaimDupe( ClaimDupe* seed, ClaimNode* branch, vector<ClaimDupe*>& paths, bool drxn );
    static bool create( ClaimNode* fork, vector<ClaimDupe*> paths[2], bool drxn );
    static bool dump( ClaimNode* fork, ClaimDupe* l, ClaimDupe* r, vector<ClaimDupe*> paths[2], int& pairing );
    void dupeNodes( unordered_set<ClaimNode*>& dupes, NodeRoll& cloned, Nodes& base, bool drxn );
    bool dupePath( Node* fork, NodeRoll& cloned, Nodes& base, bool drxn );
    void dupeEdges( unordered_set<ClaimNode*>& trims, bool drxn );
    void extend( vector<ClaimDupe*>& paths, bool drxn );
    void fill( ClaimNode* cn, bool drxn );
    bool isDupe( vector<ClaimDupe*>& paths );
    bool isDupe( vector<ClaimDupe*> claim[2], unordered_set<ClaimNode*>& alts, bool drxn );
    static bool setDumps( ClaimNode* fork, vector<ClaimDupe*> paths[2] );
    static bool setDupes( Node* fork, vector<ClaimDupe*> paths[2], NodeRoll& nodes, bool drxn );
    static bool setIgnores( vector<ClaimDupe*> paths[2] );
    static void setSplits( vector<ClaimDupe*> paths[2] );
    vector<ClaimNode*> path;
    unordered_set<ClaimNode*> pathed;
    vector<ClaimDupe*> dumped, kept;
    vector<int> pairings;
    bool duped, ignored;
};

struct ClaimTrim
{
    static bool trim( ClaimNode* fork, bool drxn );
private:
    ClaimTrim( ClaimNode* fork, ClaimNode* branch, ClaimNode* bridge, bool drxn );
    void fill( ClaimNode* cn, unordered_set<ClaimNode*>& group, bool drxn );
    void score( int l, int r, int& hits, int& unique );
    bool trim( int ol );
    ClaimNode* forks[2];
    unordered_set<ClaimNode*> groups[2][2], shared[2];
    int ols[2];
};

class ClaimNode
{
public:
    ClaimNode( ClaimNode* cn, NodeRoll& cloned, Nodes& bases );
    ~ClaimNode();
    void addEdge( ClaimNode* cn, int ol, int32_t diff, bool isLeap, bool drxn, bool reciprocate=true );
    static bool claim( Node* fork, NodeRoll& nodes, bool drxn );
    int32_t* get( Node* node );
    Node* getFork( bool drxn );
    vector<Node*> getForks( Node* fork, bool drxn );
    ClaimPairing* getPairing( ClaimNode* cn, bool drxn );
    bool isBridge();
    bool isForked( bool drxn );
    bool isInvalid( Node* node, int32_t off, ClaimNode* cn, bool drxn );
    bool isShared( Node* node, ClaimNode* cn, bool drxn );
    bool isValid( Node* node, bool drxn );
    bool reach( ClaimNode* tar, ClaimNode* forks[2], int32_t diff, vector<int32_t>& diffs, bool path[2], bool forked, bool drxn );
    bool reach( ClaimNode* tar, ClaimNode* fork, int32_t diff, vector<int32_t>& diffs, bool forked, bool drxn );
    bool reach( ClaimNode* tar, int32_t diff, vector<int32_t>& diffs, bool drxn );
    bool removeEdge( ClaimNode* cn, bool drxn, bool reciprocate=true );
    void setBase( Nodes& base, bool distal );
    bool setState( bool drxn );
    ClaimNode* clone_;
    vector<ClaimEdge> edges_[2];
    vector<ClaimJoin*> joins_[2];
    vector<ClaimPairing*> pairs_[2];
    vector<ClaimRedundant*> redundant_[2];
    vector<ClaimShared> shared_;
    vector<ClaimNode*> alts_;
    Nodes distal_;
    bool split_;
    
private:
    ClaimNode( Node* node, vector<ClaimNode*>& claims, int32_t coord, bool orient, bool drxn );
    void addDistal( Node* node, bool drxn );
    bool addEdge( Edge& e, vector<ClaimNode*>& claims, bool orient, bool drxn );
    void addJoin( ClaimNode* cn, Node* hang, NodeOffsets& offs, int32_t dist, bool drxn );
//    void addJoin( ClaimNode* cn, int32_t diff, bool drxn, bool reciprocate=true );
    static bool branch( vector<ClaimNode*>& claims, vector<ClaimNode*>& forks, bool orient, bool test );
    static bool complete( vector<ClaimNode*>& claims, bool drxn );
    static bool create( Node* fork, vector<ClaimNode*>& claims, vector<ClaimNode*>& forks, bool drxn );
    bool dupe( NodeRoll& nodes, int drxn );
    void edge( vector<ClaimNode*>& claims, bool orient, bool drxn );
    void extend( vector<ClaimNode*>& claims, bool orient, bool drxn );
    void fill( Node* node, Nodes& block, int32_t dist, int32_t limit, bool orient, bool drxn );
    void fill( Node* node, Nodes& include, int32_t dist, bool orient, bool drxn );
    bool findFork( Node* q, bool drxn );
    int32_t getCoord( Node* fork, Edge& e, int32_t dist, bool orient, bool drxn );
    vector<Edge> getEdges( Node* fork, bool drxn, bool blunt=false );
    Nodes getExts();
    vector<Node*> getForks( Node* fork );
    vector< pair<Node*, int32_t> > getHangs( bool orient, bool drxn );
    bool isBranched( bool drxn );
    bool isExtend( int32_t coord );
    bool reach( ClaimNode* tar, bool drxn );
    bool removeJoin( ClaimNode* cn, bool drxn, bool reciprocate=true );
    void reset();
    static bool resplit( vector<ClaimNode*>& forks, vector<ClaimNode*>& claims, bool orient );
    static void score( vector<ClaimNode*>& claims, bool orient );
    void setBlocked( Node* node, bool drxn );
//    void setBlocks( Nodes& block, bool drxn );
    bool setBranched( vector<ClaimNode*>& claims, bool branched, bool orient, bool drxn );
    void setDistal( bool drxn );
    static void setEnds( vector<ClaimNode*>& ends, Nodes& base, bool orient, bool drxn );
    bool setExts( Nodes& shared, int32_t limit, bool drxn );
    void setForked( vector<ClaimNode*>& claims, bool orient, bool drxn );
    void setForks( vector<ClaimNode*>& forks, bool drxn );
    static void setJoins( vector<ClaimNode*>& claims, Nodes& base, int32_t limit, bool orient, bool drxn );
    void setKeep( Node* node, Nodes& keep, bool drxn );
    void setPairs( ClaimNode* cn, int32_t diff, bool lFork, bool rFork );
    bool setRedundant( vector<ClaimNode*>& claims, ClaimPairing* cp, Node* paired[2], int32_t cutoff, int32_t est, int hits, bool orient );
    void setScores( vector<ClaimNode*>& claims, Nodes& base );
    void setScores( vector<ClaimNode*>& claims, Nodes& base, bool orient );
    bool split( vector<ClaimNode*>& claims, bool orient, bool drxn );
    void testLoop( ClaimNode* cn );
    bool trim( int drxn );
    
    vector<Node*> path_;
    unordered_map<Node*, int32_t> offs_;
    int32_t ends_[2];
    bool edged_[2];
};

struct Claim
{
//    Claim( Node* fork, Edge& e, int32_t dist, bool branch, bool orient, bool drxn );
    Claim( Node* fork, int32_t dist, bool orient, bool drxn );
    Claim( Claim* fork, Edge& e, bool orient, bool drxn );
    
    void advance( Node* fork, bool orient, bool drxn );
    static void block( vector<Claim*>& claims, bool drxn );
    void block( Node* node, Nodes& blocked, Nodes& keep, bool drxn );
    void cull( Nodes& keep );
//    int32_t estimate( Claim* f, int32_t& best, int32_t est, bool drxn );
    bool identical( Claim* rhs, bool drxn );
    Node* place( int i, bool drxn );
    bool redundant( Claim* rhs, bool drxn, bool reciprocate=true );
    void retract();
    vector<Claim*> split( bool init, bool orient, bool drxn );
    static Nodes target( vector<Claim*>& branches );
    Node* terminus( bool drxn );
    
    vector<Node*> path;
    vector<Edge> edges;
    NodeOffsets offs;
    int32_t coord;
};

struct ClaimPair
{
    ClaimPair( Claim* fwd, Edge e, int altCount );
    int score();
    Claim* path;
    Edge edge;
    int32_t work;
    vector<int> alts, pairs;
    int uniques;
    bool drop;
};

struct ClaimBlocks
{
    ClaimBlocks( vector<Claim*> test[2] );
    bool ambiguous( vector<Claim*> test[2], Node* bck, Node* fwd, bool drxn );
    bool disregard( vector<Claim*> test[2], Node* bck, Node* fwd, bool drxn );
    vector<Nodes> blocks;
    vector< pair<Nodes, Nodes> > used;
};

struct ClaimScores
{
    ClaimScores( Claim* path, vector<Claim*> test[2], bool drxn );
    static void create( vector<Claim*> test[2], vector<ClaimScores*>& scores, bool orient, bool drxn );
    bool contest();
    static bool disregard( Node* bck, Node* fwd, Nodes& used, bool drxn );
    int get( ClaimScores* alt, Claim* branch );
//    int get( ClaimScores* alt, Claim* branch, int& ol );
    bool paired( Claim* c, bool exists=false );
    bool preferred( int hits[2][3], int ols[2], int altHits[2][3], int altOl, bool dual );
    bool remove( Claim* branch );
    static void score( vector<Claim*> test[2], vector<ClaimScores*>& scores, bool orient, bool drxn );
    static bool score( vector<ClaimScores*>& scores, int32_t cutoff, bool drxn );
//    static bool score( vector<ClaimScores*>& scores, int32_t cutoff );
    void score( Node* b, Node* f, int32_t est, int32_t& best, bool drxn );
    void score( int hit, bool unique );
    bool set( ClaimScores* alt, Claim* weak, Claim* strong, int hits[2][3], int ols[2], bool& paired );
    void trim( bool drxn );
    Claim* path;
    vector<ClaimScores*> alts;
    vector<ClaimPair> pairs;
    bool claimed, dropped;
};

struct ClaimMap
{
    ClaimMap( ClaimScores* cs, Nodes& bases, NodeRoll& cloned, bool drxn );
    void claim( Claim* c, bool d );
    bool confirm( bool drxn );
    void connect( ClaimPair& cp, bool drxn );
//    void disconnect( bool drxn );
    void disconnect( ClaimPair& cp, bool drxn );
    Node* get( Node* node );
    bool match( ClaimScores* cs, bool drxn );
    void unmatch( Claim* path, bool drxn );
    vector<ClaimScores*> scores[2];
    Nodes& bases, alts[2]; // alts[0] belong to different paths (must clone), alts[1] belong to the same paths, but shared by different pairs (must edge to)
    NodeRoll& cloned, safe[2];
    unordered_set<Claim*> claims[2];
    unordered_map<Node*, Node*> used;
    Node* fork;
    bool claimed;
};

class ClaimFork
{
public:
    ClaimFork( Node* fork, NodeRoll& nodes, bool drxn );
    ClaimFork( Node* fork, Node* branch, bool drxn );
    ClaimFork( Node* fork, vector<Claim*>& borrow, NodeRoll& nodes, bool orient, bool drxn );
    ~ClaimFork();
    
    NodeRoll cloned_;
    bool claimed_;
private:
    void add( Node* fork, Edge e, int32_t dist, bool drxn );
//    void add( Node* fork, Edge& e, int32_t dist, bool branch, bool drxn );
    bool claim( NodeRoll& nodes );
    bool complete( bool drxn );
    bool confirm();
    void contest( Claim* branch, vector<Claim*>& contested, bool drxn );
    bool retry();
    
    Node* fork_;
    vector<Claim*> paths_[2], added_;
//    vector<ClaimFork*> extends_[2];
    vector<ClaimScores*> scores_;
    bool orient_, drxn_;
};

#endif /* NODE_CLAIM_H */

