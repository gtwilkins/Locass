/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   locus_fill.h
 * Author: glen
 *
 * Created on 23 June 2019, 10:40 PM
 */

#ifndef LOCUS_FILL_H
#define LOCUS_FILL_H

#include "node.h"
#include "query.h"

struct FillNode;

struct FillPair
{
    FillPair( Node* node, ReadId id, int32_t est ): node( node ), id( id ), est( est ){};
    Node* node;
    ReadId id;
    int32_t est;
};

struct FillPaired
{
    FillPaired( Node* hit, Node* base, ReadId id, int32_t est ): node( hit ), pairs{ FillPair( base, id, est ) }{};
    Node* node;
    vector<FillPair> pairs;
};

struct FillRead
{
    FillRead( string seq, ReadId id, int32_t coord, bool unseeded ): seq( seq ), id( id ), coord( coord ), mapped( 0 ), unseeded( unseeded ){};
    string seq;
    ReadId id;
    int32_t coord, mapped;
    bool unseeded;
};

struct FillHit
{
    FillHit( Node* hit, int32_t coords[2], int32_t hang, int32_t diff );
    static void add( vector<FillHit>& hits, FillNode* fn, string& seq, int32_t coords[2], int32_t est, int32_t off, bool drxn );
    static bool confirm( vector<FillHit> hits[2], bool ext[2], int seqLen, bool& ignore  );
    void fill( Node* node, int32_t l, int32_t r );
    void set( ReadId id, bool ignore );
    vector<Node*> hits;
    vector<int32_t> mapped[2];
    int32_t len, diff, hang;
};

struct FillNode
{
    FillNode( Node* node, NodeRoll& nodes, vector<FillNode*>& fills, Nodes& used );
    void addRead( NodeRoll& nodes, Node* node, ReadId id, int32_t coord, int32_t dist, bool drxn );
    void addPair( Node* hit, Node* base, ReadId id, int32_t est, bool drxn );
    static vector<FillNode*> create( NodeRoll& nodes );
    void edge( NodeRoll& nodes, vector<FillNode*>& fills, Nodes& used, bool drxn );
    void getDist( unordered_map<FillNode*, int32_t>& dists, int32_t dist, int32_t limit, int32_t& best, bool drxn );
    Node* getFork( bool drxn );
    bool isBlunt( int readLimit, bool drxn );
    bool isStrong( int32_t dist, int32_t limit, bool drxn );
    int trim( bool drxn );
    vector< pair<Node*, int32_t> > path;
    vector< pair<FillNode*, int> > edges[2];
    vector< pair<ReadId, int32_t> > reads[2];
    vector<FillPaired> paired[2];
    string seq;
    int32_t len;
    bool bad;
};

class LocusFill
{
    LocusFill( Querier& bwt, NodeRoll& nodes );
    ~LocusFill();
    void clear();
    bool getDists( FillNode* fn, unordered_map<FillNode*, int32_t>& dists, int32_t& ext, int32_t& limit, bool drxn );
    void getSeeds( Querier& bwt, FillNode* fn, vector<FillRead>& seeds, bool drxn );
    string getSeq( Querier& bwt, ReadId id );
    bool join( Querier& bwt, NodeRoll& nodes );
    void map( Querier& bwt, NodeRoll& nodes, Nodes& repair, FillNode* fn, bool drxn );
    void refresh( Querier& bwt, NodeRoll& nodes );
    bool seed( Querier& bwt, NodeRoll& nodes );
    bool trim( NodeRoll& nodes );
    
    vector<FillNode*> fills_;
    unordered_set<ReadId> seeded_, mapped_, unmapped_;
    unordered_map<ReadId, string> reads_;
public:
    static bool fill( Querier& bwt, NodeRoll& nodes );
};

#endif /* LOCUS_FILL_H */

