/*
 * Copyright (C) 2017 Glen T. Wilkins <glen.t.wilkins@gmail.com>
 * Written by Glen T. Wilkins
 * 
 * This file is part of the Locass software package <https://github.com/gtwilkins/Locass>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LOCUS_GAP_H
#define LOCUS_GAP_H

#include "node_path.h"

class Bud;

struct BudPair
{
    BudPair( string seq, ReadId id, int32_t dist, bool unique ): seq_( seq ), id_( id ), dist_( dist ), count_( 0 ), unique_( unique ){ best_[0] = best_[1] = 0; };
    void align( int32_t coords[2] );
    bool alignEnds( BudPair* bp, int i, int j, bool drxn );
    string seq_;
    ReadId id_;
    int32_t dist_, best_[2];
    int count_;
    bool unique_;
};

struct Kmer
{
    BudPair* bp_;
    int off_;
};

struct Kmers
{
    void add( BudPair* bp );
    vector<Kmer>* get( string& q );
    unordered_map<string, vector<Kmer> > kmers_;
};

struct BudExt
{
    BudExt( string& seq, ReadId id, int ol, bool drxn );
    string ext_;
    ReadId id_;
    int ol_, align_;
};

struct BudMapx
{
    BudMapx( BudPair* bp, int32_t q[2], int32_t t[2], bool drxn );
    BudPair* pair_;
    int32_t ext_, q_[2], t_[2];
};

struct BudAlign
{
    BudAlign( int32_t coord, int len, int base ): coord_( coord ), len_( len ), base_( base ), good_( 0 ), ol_( 0 ), forked_( false ){};
    ~BudAlign();
    bool addExt( BudExt* ext, bool drxn );
    void addMap( BudPair* bp, int32_t q[2], int32_t t[2], bool drxn );
    bool advance( int base, bool drxn );
    int align( string& a, string& b, int i, bool drxn );
    void fold( BudAlign* ba );
    bool isBranch( int limit );
    void print( int32_t base, bool drxn );
    void setExt( bool drxn );
    void setForks();
    bool setGood( int32_t good, bool drxn );
    bool used( ReadId id );
    string ext_;
    vector<BudExt*> exts_;
    vector<BudAlign*> edges_, folds_;
    vector<BudMapx*> maps_;
    int32_t coord_;
    int len_, base_, good_, ol_;
    bool forked_;
};

struct BudNode
{
    BudNode( BudPair* bp ): bp_( bp ), good_( 0 ), bad_( false ){};
    bool isChained( int chained, int good, int tar, bool drxn );
    void removeEdge( BudNode* bn, bool drxn );
    void setBad( bool drxn );
    void setEdges( bool drxn );
    void setGood( bool drxn );
    
    vector< pair<BudNode*, pair<int, int> > > edges_[2]; // len, excess
    BudPair* bp_;
    int good_;
    bool bad_;
};

struct BudIsle
{
    BudIsle( BudNode* bn, vector<BudIsle*>& isles, bool drxn );
    void addEdge( BudNode* bn, vector<BudIsle*>& isles, int ol, bool drxn );
    void extend( vector<BudIsle*>& isles, bool drxn );
    vector<BudNode*> path_;
    vector< pair<BudIsle*, int> > edges_[2];
    vector<BudPair*> mapped_;
};

struct BudGraph
{
    void clean( bool drxn );
    BudNode* get( BudPair* bp, bool create );
    unordered_map<BudPair*, BudNode*> nodes_;
};

struct BudRead
{
    BudRead( int32_t l, int32_t r, int32_t len, bool drxn );
    int32_t coord( int d );
    int32_t coords_[3];
    bool redundant_;
};

struct BudEdge
{
    BudEdge( Bud* l, Bud* r, int32_t diff, int ol );
    Bud* edge_[2];
    int32_t diff_, off_[2];
    int ol_;
};

class Bud
{
    Bud( NodePath* np, vector<Bud*>& buds, int32_t target, bool budding, bool drxn );
    Bud( Querier& bwt, vector<Bud*>& buds, Bud* b, BudAlign* ba, bool drxn );
    Bud( Querier& bwt, vector<Bud*>& buds, Bud* b, string s, bool drxn );
    Bud( Querier& bwt, vector<Bud*>& buds, string s, bool drxn );
    ~Bud();
    void addAlign( Querier& bwt, vector<Bud*>& buds, BudAlign* ba, bool drxn );
    void addBranch( vector<Bud*>& buds, string seq[2], ReadId id[2], int ol, int32_t dist, bool drxn );
    bool addExt( Querier& bwt, vector<Bud*>& buds, BudAlign* ba, bool drxn );
    void addPair( vector<Bud*>& buds, string& seq, ReadId id, int32_t coord, int32_t dist, bool drxn );
    static bool blockPair( Querier& bwt, vector<Bud*>& buds, string& seq, bool drxn );
    void build( vector<Bud*>& buds, bool drxn );
    void getDiffs( unordered_map<Bud*, int32_t>& diffs, int32_t diff, bool drxn );
    Node* getEdgeable( int32_t off, int32_t& diff, NodeRoll& nodes, bool drxn );
    Kmers getKmers( unordered_map<Bud*, int32_t>& diffs, bool drxn );
    pair<Node*, int> getNode( ReadId id, BudRead& br );
    bool isUnique( int32_t coords[2] );
    static void print( vector<Bud*>& buds, bool drxn );
    void print( unordered_set<ReadId>& used, unordered_set<ReadId>& mapped, int32_t limit, bool drxn );
    void rebuild( vector<Bud*>& buds, int32_t target, bool drxn );
    vector< pair<ReadId, int> > setAligned( BudAlign* ba, vector<Bud*>& buds, Kmers& kmers, bool drxn );
    void setEdge( PathEdge* pe, vector<Bud*>& buds, bool drxn );
    void setExt( Querier& bwt, vector<Bud*>& buds, bool drxn );
    static bool setForks( Querier& bwt, vector<Bud*>& buds, bool drxn );
    void setFresh( Querier& bwt, NodeRoll& nodes, bool drxn );
    static void setKept( Querier& bwt, vector<Bud*>& buds, bool drxn );
    static void setIslands( Querier& bwt, vector<Bud*>& buds, bool drxn );
    static void setMaps( Querier& bwt, vector<Bud*>& buds, bool drxn );
    static int setMatch( string& q, string& t, int qq[2], int tt[2] );
    void setNodes( NodeRoll& nodes, vector<Bud*>& buds, bool drxn );
    static void setReads( Querier& bwt, vector<Bud*>& buds, bool drxn );
    void setRedundant();
    void setTarget( NodePath* np, bool drxn );
    static void sprout( Querier& bwt, NodeRoll& nodes, vector<Bud*>& buds, bool drxn );
    void test( Querier& bwt, bool drxn );
    vector<Node*> path_;
    vector<NodePath*> pathed_;
    vector<BudEdge*> edges_[2];
    vector<BudAlign*> forks_;
    unordered_map<ReadId, BudRead> reads_;
    unordered_map<Node*, int32_t> offs_;
    unordered_map<ReadId, pair<BudPair*, int32_t> > pairs_;
    vector<bool> block_;
    string seq_;
    int32_t ends_[2], target_;
    bool bud_, mapped_;
    
public:
// Define homopolymers and non-homopolymers among pairs and exts
// Search for matches among non-homopolymers
    void amend( Querier& bwt, string s, bool drxn );
    void append( Querier& bwt, string s, bool drxn );
    static bool bud( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& ends, bool drxn );
    static void test( Querier& bwt, string s, int drxn );
};

#endif /* LOCUS_GAP_H */

