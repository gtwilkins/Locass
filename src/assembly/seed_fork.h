/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   seed_fork.h
 * Author: glen
 *
 * Created on 6 January 2019, 5:43 AM
 */

#ifndef SEED_FORK_H
#define SEED_FORK_H

#include "node.h"
#include "types.h"
#include "node_structs.h"
#include "node_path.h"

struct ExtendBranch
{
    ExtendBranch( Node* node, Nodes& fwd, bool drxn );
    ExtendBranch( Node* node, bool drxn );
    Node* node;
    Nodes branch;
    int score;
};

struct ExtendScores
{
    ExtendScores( Node* node, bool drxn );
    ExtendScores( Nodes& loop, bool drxn );
    void add( Node* node, bool drxn );
    vector<ExtendBranch> branch( Node* node, bool drxn );
    vector<ExtendBranch> branch( Nodes& loop, bool drxn );
    int get( Node* node );
    unordered_map<Node*, int> scores;
    Nodes added, pathed;
};

struct ExtendAlt
{
    ExtendAlt( Node* branch );
    static void add( vector<ExtendAlt>& alts, Node* node );
    static void confirm( vector<ExtendAlt>& alts, NodeRoll& nodes, Nodes& bck, bool drxn );
    bool get( Querier &bwt, NodeRoll &nodes, NodeRoll& ext, bool drxn );
    vector< pair<Node*, int> > get( Querier &bwt, NodeRoll &nodes, int32_t cutoff, bool drxn );
//    bool prune( Querier &bwt, NodeRoll &nodes, bool drxn );
    Node* branch;
    int lastExt, gen;
    bool stopped;
};

struct ExtendFork
{
    ExtendFork( Node* a, Node* b );
    bool advance( vector<ExtendAlt>& alts, bool drxn );
//    bool advance( vector<ExtendAlt>& alts, int i, ExtendScores& scores, bool advanced, bool drxn );
//    bool advance( Node*& node, Node*& alt, vector<ExtendAlt>& alts, NodeScores &scores, bool advanced, bool drxn );
//    bool cull( Querier &bwt, NodeRoll &nodes, Nodes& bck, bool priming, bool drxn );
//    static bool extend( vector<ExtendFork> &forks, vector<ExtendAlt>& alts, Querier &bwt, NodeRoll &nodes, bool priming, bool drxn );
//    bool getAlts( Querier& bwt, NodeRoll& nodes, NodeRoll &ext, Nodes& bck, bool drxn );
//    bool getAlt( Querier& bwt, NodeRoll& nodes, NodeRoll &ext, Nodesx& bck, bool drxn );
//    static Nodes getBack( vector<SeedFork> &forks, bool drxn );
    bool get( NodeRoll &ext, vector<ExtendAlt>& alts, Nodes& bck, bool priming, bool drxn, bool first=true );
    bool getLoop( NodeRoll &ext, vector<ExtendAlt>& alts, Nodes& bck, bool drxn );
//    bool getSide( NodeRoll &ext, Nodes& bck, bool drxn );
    bool isLoop( bool drxn );
//    bool restart( Querier& bwt, NodeRoll& nodes, vector<Node*> path[2], bool drxn );
//    void test( Querier &bwt, NodeRoll &nodes, bool drxn );
    Node* branch[2],* loop[2];
//    NodeRoll sides;
    int lastExt, impatience, gen;
};

class SeedExtend
{
    vector<ExtendFork> forks;
    vector<ExtendAlt> alts;
    int gen;
public:
    void addAlt( Node* node );
    void addFork( Node* node );
    bool empty();
    bool extend( Querier& bwt, NodeRoll& nodes, bool priming, bool drxn );
    void reset();
};


#endif /* SEED_FORK_H */

