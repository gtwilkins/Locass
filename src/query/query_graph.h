/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   query_graph.h
 * Author: glen
 *
 * Created on 28 December 2018, 4:04 AM
 */

#ifndef QUERY_GRAPH_H
#define QUERY_GRAPH_H

#include "types.h"
#include "constants.h"
#include "query_binary.h"
#include "query_structs.h"
#include "query_state.h"

struct QueryRead
{
    QueryRead( ReadId id, int ol, int ext ): id( id ), ol( ol ), ext( ext ), redundant( false ){};
    void offset( int off );
    ReadId id;
    uint16_t ol, ext;
    bool redundant;
};

struct QueryNode
{
    QueryNode( string &ext, int maxOl );
    QueryNode( string &q, Overlap &ol, bool perfect );
    QueryNode( Overlap &ol, vector<QueryNode*> &hits, bool perfect );
    QueryNode( QueryNode* node, int j );
    ~QueryNode();
    void addEdge( QueryNode* qn );
//    static void cull( vector<QueryNode*> &nodes, int cutoff );
//    void cull( int cutoff );
    bool confirm();
//    void discard();
    void extend( Overlap &ol );
    string getSeq( bool drxn );
    static vector<QueryNode*> graph( QueryBinaries* qb, QState &qs, string seq, bool drxn );
    bool isBlunt();
    bool isIgnore( int len, vector<QueryNode*> &exts, bool drxn );
    bool match( Overlap &ol, vector<QueryNode*> &hits, int i, bool &added );
    bool merge();
    void removeEdge( QueryNode* node, int i );
    void resetBase();
    int setExt();
    int setReads(  int cutoff );
    bool setSeq( bool drxn );
//    void setNonAlt( vector<QueryNode*> &nodes );
    void trim();
    void trim( uint16_t ol );
    vector<QueryNode*> edges[2];
    vector<QueryRead> reads;
    string seq, ext;
    int len, base, maxOl, maxExt, readCount, merged;
    bool fixed, good, bad;
};


#endif /* QUERY_GRAPH_H */

