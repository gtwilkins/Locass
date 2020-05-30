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

struct ExtendEdge
{
    bool get( NodeRoll& nodes, NodeRoll &ext, bool drxn );
    bool get( Node* node, vector<Node*>& adds, int readCount, bool drxn );
    Node* edge[2];
    ReadId id[2];
};

struct ExtendFork
{
    ExtendFork( Node* a, Node* b );
    bool advance( vector<ExtendAlt>& alts, bool drxn );
    bool get( NodeRoll &nodes, NodeRoll &ext, vector<ExtendAlt>& alts, Nodes& bck, bool priming, bool drxn, bool first=true );
    bool getLoop( NodeRoll &ext, vector<ExtendAlt>& alts, Nodes& bck, bool drxn );
    bool isLoop( bool drxn );
    Node* branch[2],* loop[2];
//    NodeRoll sides;
    int lastExt, impatience, gen;
};

class SeedExtend
{
    vector<ExtendFork> forks_;
    vector<ExtendEdge> edges_;
    vector<ExtendAlt> alts_;
    int gen_;
public:
    void addAlt( Node* node );
    void addEdge( Node* fork, Node* branch, bool drxn );
    void addFork( Node* node );
    void cull( Querier& bwt, NodeRoll& nodes, bool drxn );
    bool empty();
    bool extend( Querier& bwt, NodeRoll& nodes, bool priming, bool drxn );
    void reset();
};


#endif /* SEED_FORK_H */

