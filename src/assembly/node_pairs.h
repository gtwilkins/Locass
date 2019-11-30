/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   node_pairs.h
 * Author: glen
 *
 * Created on 23 January 2019, 1:27 AM
 */

#ifndef NODE_PAIRS_H
#define NODE_PAIRS_H

#include "types.h"
#include "node_groups.h"

class Node;

struct NodePair
{
    NodePair( ReadId id, int32_t est, bool mp );
    void add( ReadId id, int32_t est, bool mp );
    int32_t estimate();
    int getLib( ReadId id );
    int32_t margin();
    vector< pair<ReadId, int32_t> > pairs[2];
    vector<int> libs;
    int count;
};

struct NodePairs
{
    NodePairs( Node* self ) : self( self ), count( 0 ){};
    void add( Node* node, ReadId id, int32_t est, int32_t len, bool mp, bool drxn );
    void clean();
    bool empty();
    void erase( Node* node, bool drxn );
    int get();
    int get( bool drxn );
    int get( Node* node, bool drxn );
    int get( Node* node, NodeOffsets& selfOffs, NodeOffsets& tarOffs, bool drxn );
    int get( Nodes& nodes, bool drxn );
    int get( Nodes& nodes, NodeOffsets& selfOffs, NodeOffsets& tarOffs, bool drxn );
    NodePair* getPair( Node* node, bool drxn );
    void reset();
    void setRedundant( Nodes& q, Nodes& t, Node* fwd, bool drxn );
    void test();
    unordered_map<Node*, NodePair > pairs[2];
    Node* self;
    int count;
};

#endif /* NODE_PAIRS_H */

