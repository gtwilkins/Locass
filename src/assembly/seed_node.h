/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   seed_node.h
 * Author: glen
 *
 * Created on 20 December 2019, 3:26 PM
 */

#ifndef SEED_NODE_H
#define SEED_NODE_H

#include "node.h"

class SeedNode
{
    SeedNode( Node* node, Coords* coords, int i, int j );
    static bool setEdges( string& s, vector<Node*>& ends, SeedNode& seed, NodeRoll& nodes, int ol, bool drxn );
    static bool setEdges( vector<Node*>& ends, NodeRoll& nodes, int ol, bool drxn );
    static vector<SeedNode> setNodes( Querier& bwt, NodeRoll& nodes, string& s, int drxn );
    static void setPath( vector<SeedNode>& seeds, vector<Node*> ends[2], NodeRoll& nodes, int drxn );
    Node* node_;
    int32_t mapped_[2][2];
    bool created_, ended_;
public:
    static bool seed( Querier& bwt, NodeRoll& nodes, string s, int lOl, int rOl, int drxn );
};



#endif /* SEED_NODE_H */

