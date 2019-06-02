/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   node_path.h
 * Author: glen
 *
 * Created on 4 January 2019, 1:50 AM
 */

#ifndef NODE_PATH_H
#define NODE_PATH_H

#include "node.h"
#include "node_structs.h"

struct PathScores
{
    PathScores(){};
    PathScores( Node* node ){ add( node ); };
    void add( Node* node );
    void add( Node* node, float score, bool drxn );
    float get( Node* node, int drxn );
    float get( Nodes& nodes, int drxn );
    PathScores unadded();
    unordered_map<Node*, float> scores[2];
    Nodes added;
};

struct NodeBranch
{
    NodeBranch();
    NodeBranch( Node* fork, Node* branch );
    NodeBranch( Node* fork, Node* branch, Nodes& fwdIn, NodeScores& scores, bool drxn );
    NodeBranch( Node* fork, Node* branch, Nodes& fwdIn, PathScores& scores, bool drxn );
    static void add( vector<NodeBranch>& branches, Node* fork, Nodes& added, Nodes& btw, NodeScores& scores, int32_t dist, bool drxn );
    bool cull( Nodes &ignore );
    static vector<NodeBranch> get( Node* fork, Nodes fwd, NodeScores &scores, bool drxn );
    static vector<NodeBranch> get( Node* fork, Nodes& fwd, PathScores& scores, bool drxn );
    static void sort( vector<NodeBranch> &branches );
    Node* fork,* branch;
    Nodes fwd;
    float score;
};

struct PathChunk
{
    PathChunk( NodeList &pilot, Node* alt, int &i );
    void add( vector<PathChunk> &chunks, Node* alt, int &i );
    NodeBranch alt( NodeScores& scores, Nodes& added, Nodes& btw, int& cumul, bool drxn );
    static Nodes between( vector<PathChunk> &chunks );
//    void branch( NodeScores& scores, Nodes& added, Nodes& btw, bool drxn );
//    static vector<PathChunk> fill( NodeList &pilot, NodeScores &scores, Nodes &added );
    bool take( Node* fork, PathChunk& tar, Nodes &ignore, bool drxn  );
    bool setFork( Node* branch[2], bool drxn );
    NodeList path[2];
    vector<NodeBranch> branches[2];
    int hits[2];
    bool edges[2][2];
};

class NodePath
{
public:
    NodePath( vector<Node*> a, vector<Node*> b );
    NodePath( Node* seed, PathScores& scores );
    static void branch( vector<NodePath>& path, PathScores& scores, bool drxn );
    void branch( NodeBranch& best, PathScores& scores, PathScores& altScores, Nodes& btw, int& cumul, bool drxn );
    bool branch( NodeBranch& best, vector<NodePath>& path, PathScores& scores, vector<Node*>& ext, bool drxn );
    static vector< vector<NodePath> > create( NodeRoll& nodes );
    bool fill( vector<NodePath>& path, PathScores& scores, bool drxn );
    static bool path( NodeList& path, PathScores &scores, bool drxn );
    bool set( Node* branch[2], bool drxn );
    bool unpause( NodeList& path, PathScores& scores );
    vector<Node*> path_[2];
};


#endif /* NODE_PATH_H */

