/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   path_seed.h
 * Author: glen
 *
 * Created on 10 February 2018, 9:46 PM
 */

#ifndef PATH_SEED_H
#define PATH_SEED_H

#include "types.h"
#include "node.h"


//struct SeedBranch
//{
//    SeedBranch( Edge &e, NodeList &tNodes, NodeSet &shared, bool drxn );
//    bool operator>( SeedBranch &rhs );
//    Node* n;
//    int ol, score[2], reads;
//};

//struct SeedPath
//{
//    SeedPath( Node* node, NodeSet &usedSet, bool drxn );
//    
//    NodeList getTargets( bool drxn );
//    static NodeSet getUsed( vector<SeedPath> &paths, bool drxn );
//    static bool resolve( vector<SeedPath> &paths, NodeSet &altFwd, NodeSet &usedSet, int i, bool drxn );
//    void setBranches( Node* node, NodeSet &altFwd, NodeSet &usedSet, NodeList &tNodes, bool drxn );
//    bool setLoop( Node* node, NodeList &tNodes, bool drxn );
//    bool setNext( NodeSet &altFwd, NodeSet &usedSet, NodeSet &bckSet, NodeList &tNodes, bool drxn );
//    
//    
//    NodeList path, branches;
//    NodeSetList branchFwd;
//    NodeListList alts;
//    NodeSet starts, shared;
//    bool ended, doBranch;
//};

//class PathSeed
//{
//public:
//    PathSeed( NodeList &nodes );
//    void exportAlign( ofstream &fp );
//    void plot( bool drxn );
//    
//private:
//    void setEnds();
//    
//    vector<SeedPath> paths_;
//    NodeList &nodes_, ends_[2];
//    NodeSet usedSet_;
//};



#endif /* PATH_SEED_H */

