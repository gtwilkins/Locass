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

#ifndef LOCUS_TRAVERSE_H
#define LOCUS_TRAVERSE_H

#include "node.h"
#include "query.h"
#include "kmer_graph.h"

struct Vertex;

struct Link
{
    Link( Vertex* base, Vertex* edge, int ol, bool drxn );
    Vertex* node_[2];
    int coord_[2], ol_;
};

struct VRead
{
    VRead( ReadId id, int l, int r ): id_( id ){ coord_[0] = l; coord_[1] = r; };
    bool isRedundant( VRead& redundant );
    ReadId id_;
    int coord_[2];
};

struct VPairs
{
    VPairs( Vertex* a, Vertex* b, int aCoord, int bCoord, bool bDrxn );
    static int set( Vertex* l, Vertex* r );
    Vertex* node_[2];
    vector<VRead> pairs_;
    int coord_[2], dist_;
};

struct Vertex
{
    Vertex(){};
    Vertex( Ext* ext, bool drxn );
    Vertex( string& seq, ReadId id, bool drxn );
    void addPair( ReadId id, VRead& vr, bool drxn );
    void addRead( VRead& vr );
    void discardPaired();
    void get( vector<Vertex*>& nodes, int coord, int limit, bool drxn );
    VRead* getPair( ReadId id, bool drxn );
    void refill( int drxns );
//    int setPairs( Vertex* base, int coord, bool drxn, bool cont );
    string seq_;
    vector<VRead> reads_, redundant_;
    unordered_map<ReadId, VRead> pairs_[2];
    unordered_map<ReadId, vector<int> > alts_[2];
    vector<VPairs*> paired_[2];
    vector<Link*> edges_[2];
    int coord_[2];
};

struct VBranch
{
    VBranch( Link* l, Vertex* base, bool drxn );
    Vertex* node_;
    int end_, flay_, baseCoord_, baseHits_, altCoord_, altHits_;
};

struct TFork
{
    TFork( string seq, bool drxn );
    TFork( string& base, string seq, ReadId id, int ext, int ol, bool drxn );
    TFork( TFork* tf, string seq, ReadId id, int ext, int ol, bool drxn );
    ~TFork();
    bool add( string& seq, ReadId id, bool drxn );
    bool add( string& seq, ReadId id, int ol, bool drxn );
    bool update();
    string base_, ext_, seq_;
    vector<TFork*> exts_;
    vector< pair<ReadId, int> > reads_;
    int ol_;
    bool exted_;
};

class Traverse
{
    Traverse( Querier& bwt, Node* node, bool drxn );
    
    void print( TFork* fork, bool drxn );
    void setFlayed( Querier& bwt, int maxLen, bool flayDrxn, bool endDrxn );
    void setGraph( Vertex* v, vector<Vertex*>& base, int coord, bool pairDrxn, bool drxn );
    void setJumps( Vertex* v, vector<Vertex*>& base, int coord, bool drxn );
    bool setMatches( vector<Vertex*>& nodes, ReadId id, int minOl, bool drxn );
    void setPairs();
    void setReads( int len, int endDrxn );
    void setTraces( bool drxn );
    
    Querier* bwt_;
    vector<Node*> path_;
    unordered_map<Node*, int32_t> coords_;
    unordered_map<ReadId, string> reads_;
    vector<Vertex*> tips_;
    KmerGraph* graphs_[2];
    Vertex* base_;
    int32_t len_;
public:
    static bool leap( Querier& bwt, Node* node, bool drxn );
    static bool trace( Querier& bwt, Node* node, bool drxn );
};

#endif /* LOCUS_TRAVERSE_H */

