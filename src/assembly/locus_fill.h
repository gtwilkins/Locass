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
    FillPaired( Node* hit, Node* base, Coords* coords, ReadId id, int32_t est ): node( hit ), pairs{ FillPair( base, id, est ) }{};
    Node* node;
//    int32_t ends[2];
    vector<FillPair> pairs;
};

struct FillRead
{
    FillRead( ReadId id, int32_t coord, bool paired ): id( id ), coord( coord ), paired( paired ){};
    ReadId id;
    int32_t coord;
    bool paired;
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

struct FillBranch
{
    FillBranch( string seq, FillRead& fr, int base );
    FillBranch( FillBranch* fb, int ext, bool drxn );
    ~FillBranch();
    bool add( string& seq, FillRead& fr, int ext, bool drxn );
    static bool confirm( vector<FillBranch*>& branches, bool drxn );
//    int confirm( Querier& bwt, int coverage, bool drxn );
    bool cull( Querier& bwt, string seq, int coverage, bool drxn );
    void get( string seq, vector<FillBranch*> path, vector<string>& seqs, vector< vector<FillBranch*> >& paths, bool refresh, bool drxn );
    void get( string seq, vector<string>& seqs, bool drxn );
    void merge( bool drxn );
    void setCounts();
    string seq_;
    vector<FillBranch*> branches_;
    vector< pair<FillRead, int> > exts_;
    int32_t base_, good_, unpaired_, total_;
    bool culled_, added_;
};

struct FillFork
{
    FillFork( string& seq, string& base, FillRead& fr, int32_t coords[2], bool drxn );
    ~FillFork();
    bool add( string& seq, FillRead& fr, int32_t coords[2], bool drxn );
    void add( string& seq, vector<FillBranch*>& path, FillRead& fr, int32_t coords[2], int len, bool drxn );
    int confirm( Querier& bwt, int unpairedCut, int totalCut, bool drxn );
    bool cull( Querier& bwt, int coverage, bool drxn );
    bool get( vector<string>& seqs, vector< vector<FillBranch*> >& paths, bool refresh, bool drxn );
    void get( vector<string>& seqs, bool drxn );
    string seq_;
    vector<FillBranch*> branches_;
    unordered_set<ReadId> used_;
    int32_t coord_, ol_;
};

struct FillNode
{
    FillNode( Node* node, NodeRoll& nodes, vector<FillNode*>& fills, Nodes& used );
    bool addBranch( string& seq, FillRead& fr, int32_t coords[2], bool drxn );
    void addRead( NodeRoll& nodes, Node* node, ReadId id, int32_t coord, int32_t dist, bool drxn );
//    void addPair( Node* hit, Node* base, ReadId id, int32_t est, bool drxn );
    static vector<FillNode*> create( NodeRoll& nodes );
    void edge( NodeRoll& nodes, vector<FillNode*>& fills, Nodes& used, bool drxn );
    void getDist( unordered_map<FillNode*, int32_t>& dists, int32_t dist, int32_t limit, bool drxn );
    void getDist( unordered_map<FillNode*, int32_t>& dists, int32_t dist, int32_t limit, int32_t& best, bool drxn );
    Node* getFork( bool drxn );
    bool isBlunt( int readLimit, bool drxn );
    bool isStrong( int32_t dist, int32_t limit, bool drxn );
    int trim( bool drxn );
    vector< pair<Node*, int32_t> > path_;
    vector< pair<FillNode*, int> > edges_[2];
    vector<FillRead> reads_[2];
    vector<FillPaired> paired[2];
    vector<FillFork*> forks_[2];
    string seq_;
    int32_t len;
    float cover_;
    bool bad;
};

class LocusFill
{
    LocusFill( Querier& bwt, NodeRoll& nodes );
    ~LocusFill();
    void clear();
    int extend( Querier& bwt, NodeRoll& nodes, FillNode* fn, bool drxn );
    bool getDists( FillNode* fn, unordered_map<FillNode*, int32_t>& dists, int32_t& ext, int32_t& limit, bool drxn );
//    void getSeeds( Querier& bwt, FillNode* fn, vector<FillRead>& seeds, bool drxn );
    string getSeq( Querier& bwt, ReadId id );
    bool join( Querier& bwt, NodeRoll& nodes );
    void map( Querier& bwt, NodeRoll& nodes, Nodes& repair, FillNode* fn, bool drxn );
    void refresh( Querier& bwt, NodeRoll& nodes );
    void seed( Querier& bwt, NodeRoll& nodes, FillNode* fn, unordered_set<ReadId>& added, bool drxn );
    bool seed( Querier& bwt, NodeRoll& nodes );
    bool trim( NodeRoll& nodes );
    
    vector<FillNode*> fills_;
    unordered_set<ReadId> seeded_, paired_, mapped_, unmapped_;
    unordered_map<ReadId, string> reads_;
public:
    static bool fill( Querier& bwt, NodeRoll& nodes );
};

#endif /* LOCUS_FILL_H */

