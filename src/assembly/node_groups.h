/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   node_groups.h
 * Author: glen
 *
 * Created on 23 January 2019, 1:42 AM
 */

#ifndef NODE_GROUPS_H
#define NODE_GROUPS_H

#include "types.h"

class Node;
struct Edge;

struct Nodes
{
    Nodes(){};
    Nodes( Node* node );
    Nodes( vector<Node*> base );
    Nodes( Node* node, bool drxn, bool inclNode, bool inclClone );
    Nodes( Node* node, int32_t limit, bool drxn, bool inclNode );
    Nodes( Node* node, int32_t limit, bool drxn, bool inclNode, bool inclClone );
//    Nodes( Node* base, NodeOffsets &offs, bool drxn, bool inclBase );
//    Nodes( Node* base, Node* avoid, NodeOffsets &offs, bool drxn, bool inclBase );
    Nodes( vector<Node*> &base, bool drxn, bool inclBase );
    void operator += ( Node* node );
    void operator += ( Nodes &rhs );
    void operator -= ( Node* node );
    void operator -= ( Nodes &rhs );
    bool add( Node* node );
    static Nodes between( Node* a, Node* b, bool inclNodes );
    static void cancel( Nodes &a, Nodes &b );
    void clear();
    static Nodes connected( Node* node );
    void dumpBad( bool bad );
    void dumpGraph( int graph, bool dump );
    bool empty();
    bool erase( Node* node );
    Nodes exclusive( Nodes& ignore );
    void fill( Node* node );
    void fill( Node* node, bool drxn, bool inclNode, bool inclClone );
//    bool fill( Node* node, NodeOffsets &offs, bool drxn );
    void fill( Node* node, int32_t limit, bool drxn, bool inclNode, bool inclClone );
    void fillBad( Node* node, bool bad, bool drxn );
    void fillBranch( Node* node );
    void fillIn( Node* node, Nodes &include, bool drxn, bool inclNode );
    void fillIn( Node* node, Nodes &include, int32_t limit, bool drxn, bool inclNode );
    void fillNext( Node* node, bool drxn );
    void fillNot( Node* node, Nodes &ignore, bool drxn, bool inclNode );
    void fillNot( Node* node, Nodes &ignore, int32_t limit, bool drxn, bool inclNode );
    bool find( Node* node );
    static void forkSets( Node* fork, Node* branch, Nodes fwds[2], Nodes bcks[2], int32_t limit, bool drxn );
    static Nodes inSet( Node* node, Nodes &include, bool drxn, bool inclNode );
    static Nodes inSet( Node* node, Nodes &include, int32_t limit, bool drxn, bool inclNode );
    static Nodes isBad( Node* node, bool bad, bool drxn );
    static Nodes isBranch( Node* node );
    static Nodes notSet( Node* node, Nodes &ignore, bool drxn, bool inclNode );
    static Nodes notSet( Node* node, Nodes &ignore, int32_t limit, bool drxn, bool inclNode );
    int size();
    unordered_set<Node*> nodes;
};

struct NodeRoll
{
    NodeRoll(){};
    NodeRoll( Node* node );
    NodeRoll( Nodes base );
    Node* operator[]( int i );
    void operator += ( NodeRoll &rhs );
    void operator += ( Nodes rhs );
    void operator -= ( Nodes& rhs );
    void operator += ( Node* node );
    void operator -= ( Node* node );
    bool add( Node* node );
    void add( NodeRoll &rhs );
    void clear();
    static NodeRoll clones( Node* node );
    bool find( Node* node );
    void erase( Node* node );
    void erase( Node* node, int &i );
    void drop( Node* node );
    bool empty();
    NodeRoll getGraph( int drxn );
    static NodeRoll next( Node* node, bool drxn );
    void remove( Node* node );
    void print( string fn, int32_t i, int32_t j );
    void print( ofstream& ofs, Node* node, Nodes used[2], int32_t limit );
    int size();
    void test( bool loop=false );
    vector<Node*> nodes;
};

struct NodeOffset
{
    NodeOffset( int dist ){ dists[0] = dists[1] = dists[2] = dist; }
    bool add( int32_t dist, bool minOnly=false );
    int32_t diff( int32_t est );
    int32_t diff( NodeOffset& rhs, int32_t est, bool drxn );
    int32_t& operator[]( int i ){ return dists[i]; }
    int32_t dists[3];
};

struct NodeOffsets
{
    NodeOffsets(){};
    NodeOffsets( Node* node, bool drxn, bool inclNode, bool inclClone=true );
    NodeOffsets( Node* node, int32_t limit, bool orient, bool drxn, bool inclNode );
    bool add( Node* node, int32_t dist );
    bool add( Node* fork, Edge& e, int32_t base, bool orient, bool drxn );
    void clear();
    bool doLoop( Node* clone, int32_t dist, bool drxn );
    bool empty();
    void erase( Node* node );
    static int32_t extend( Node* fork, Edge& e, int32_t dist, bool orient, bool drxn );
    bool find( Node* node );
    void fill( Node* node, int32_t dist, int32_t limit, bool orient, bool drxn, bool inclNode, bool inclClone, bool looped=false );
    void fillIn( Node* node, Nodes& include, int32_t dist, bool orient, bool drxn, bool inclNode, bool looped=false );
    void fillNot( Node* node, Nodes& ignore, int32_t dist, int32_t limit, bool orient, bool drxn, bool inclNode, bool looped=false );
//    void fillNot( Node* node, Nodes& ignore, int32_t dist, bool orient, bool drxn, bool inclNode, bool looped=false );
//    void fill( Node* node, Nodes& nodes, bool drxn, bool inclNode, bool inclClone );
//    void flip( int32_t dist, bool drxn );
//    static void forkSets( Node* fork, Edge& e, NodeOffsets offs[2], Nodes fwds[2], Nodes bcks[2], int32_t limit, bool drxn );
    NodeOffset* get( Node* node );
    void init( Node* node, int32_t base, int32_t limit, bool drxn, bool inclNode, bool rev );
//    static NodeOffsets offset( Node* node, int32_t base, int32_t limit, bool drxn, bool inclNode, bool rev );
//    void offset( int32_t dist );
//    void set( Node* node, int32_t dist, bool drxn, bool inclNode );
//    bool set( Node* node, int32_t dist, int32_t limit, bool drxn, bool inclNode, bool rev );
    unordered_map<Node*, NodeOffset> map;
};


#endif /* NODE_GROUPS_H */

