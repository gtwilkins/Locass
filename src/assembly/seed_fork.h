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

struct AltFork
{
    AltFork( Node* branch );
    static void add( vector<AltFork>& alts, Node* node );
    static void confirm( vector<AltFork>& alts, NodeRoll& nodes, Nodes& bck, bool drxn );
    bool get( Querier &bwt, NodeRoll &nodes, NodeRoll& ext, bool drxn );
    vector< pair<Node*, int> > get( Querier &bwt, NodeRoll &nodes, int32_t cutoff, bool drxn );
//    bool prune( Querier &bwt, NodeRoll &nodes, bool drxn );
    Node* branch;
    int lastExt, gen;
    bool stopped;
};

class SeedFork
{
public:
    SeedFork( Node* a, Node* b );
    bool advance( bool drxn );
    bool advance( Node*& node, Node*& alt, NodeScores &scores, bool advanced, bool drxn );
//    bool cull( Querier &bwt, NodeRoll &nodes, Nodes& bck, bool priming, bool drxn );
    static bool extend( vector<SeedFork> &forks, Querier &bwt, NodeRoll &nodes, bool priming, bool drxn );
//    bool getAlts( Querier& bwt, NodeRoll& nodes, NodeRoll &ext, Nodes& bck, bool drxn );
    bool getAlt( Querier& bwt, NodeRoll& nodes, NodeRoll &ext, Nodes& bck, bool drxn );
//    static Nodes getBack( vector<SeedFork> &forks, bool drxn );
    bool getExt( NodeRoll &ext, Nodes& bck, bool priming, bool drxn, bool first=true );
//    bool getSide( NodeRoll &ext, Nodes& bck, bool drxn );
    bool restart( Querier& bwt, NodeRoll& nodes, vector<Node*> path[2], bool drxn );
//    void test( Querier &bwt, NodeRoll &nodes, bool drxn );
    Node* branch[2];
//    NodeRoll sides;
    vector<AltFork> alts;
    int lastExt, impatience, gen;
};


#endif /* SEED_FORK_H */

