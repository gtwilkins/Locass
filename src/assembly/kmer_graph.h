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

#ifndef KMER_GRAPH_H
#define KMER_GRAPH_H

#include "types.h"

struct KmerNode;

struct KmerRead
{
    KmerRead( string seq, ReadId id ): node_( NULL ), seq_( seq ), id_( id ){};
    KmerNode* node_;
    string seq_;
    ReadId id_;
};

struct Kmer
{
    Kmer( KmerRead* kr, int coord ): coords_{ make_pair( kr, coord ) }{ for ( int d : { 0, 1 } ) for ( int i = 0; i < 4; i++ ) edges_[d][i] = NULL; };
    vector< pair<KmerRead*, int> > coords_;
    Kmer* edges_[2][4];
};

struct KmerCoord
{
    KmerCoord( KmerRead* kr, int base, int node, int len ): read_( kr ), base_( base ), node_( node ), len_( len ){};
    KmerRead* read_;
    int base_, node_, len_;
};

struct KmerNode
{
    KmerNode(){};
    KmerNode( string seq ): seq_( seq ){};
//    void add( KmerRead* read, int coord );
    void addEdge( KmerNode* node, int ol, bool drxn );
    void get( unordered_map<KmerNode*, int>& got, int dist, bool drxn );
    int getTerminalLen( bool drxn );
    bool isRedundant( KmerRead* read );
    string seq_;
    vector< pair<KmerNode*, int> > edges_[2];
    vector< pair<KmerRead*, int> > reads_, redundant_;
//    vector<KmerNode*> edges2_[2];
//    vector<KmerCoord> reads2_;
};

struct KmerExt
{
    KmerExt(): read_( NULL ), ext_( 0 ), ol_( 0 ){};
    KmerExt( KmerRead* read, int ol ): read_( read ), ext_( read->seq_.size() - ol ), ol_( ol ){};
    KmerExt( KmerExt& kc, int ext ): read_( kc.read_ ), ext_( kc.ext_-ext ), ol_( kc.ol_+ext ){};
    int match( KmerExt& ext, bool drxn );
    int match( KmerExt& ext, bool exact, bool drxn );
    void operator += ( int ol );
    KmerRead* read_;
    int ext_, ol_;
};

struct KmerOverlaps
{
    KmerOverlaps(): node_( NULL ){};
    void addExt( KmerRead* read, int ol, bool drxn );
    bool addOverlap( KmerRead* read, int base, int ol, bool drxn );
    void test( string& t, int ext, bool drxn );
    vector<KmerExt> exts_[2];
    vector< pair<KmerExt, int> > splints_[2];
    KmerNode* node_;
};

class KmerGraph
{
    void addMapped( KmerNode* node, KmerRead* read, int coord );
    void addSplints( KmerOverlaps* branch, KmerOverlaps* base, int ext, int ol, bool drxn );
    void advance( KmerOverlaps* ols, int ext, bool drxn );
    void create( KmerOverlaps* ols );
    void extend( KmerRead* read, bool drxn );
    bool getBranches( KmerOverlaps* ols, vector< pair<KmerOverlaps*, int> >& branches, bool branched, bool drxn );
    CharId getK( string& seq, int i );
    Kmer* getKmer( CharId k );
    string getSeq( Kmer* k );
    bool isBackOverlap( int base, int ol, int cut );
    void setBase( string seq, bool drxn );
    void setBranch( KmerRead* read, KmerNode* node, KmerOverlaps* branch, int ext, int ol, bool drxn );
    bool setExtend( KmerExt& ext, KmerOverlaps* ols, KmerOverlaps*& fork, bool drxn );
    KmerOverlaps* setOverlaps( string seq, KmerOverlaps* ols=NULL, int ext=0, bool map=true, bool drxn=true );
    void test( KmerOverlaps* ols );
    
    vector<KmerRead*> reads_;
    unordered_map<CharId, Kmer*> kmers_;
    unordered_set<KmerRead*> reused_;
    CharId mask_;
    int len_;
    
public:
    KmerGraph( int len );
    void add( ReadId id, string seq );
    void assemble( string seq, bool drxn );
};

#endif /* KMER_GRAPH_H */

