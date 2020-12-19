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

#ifndef REVIEWER_H
#define REVIEWER_H

#include "node.h"
#include "query.h"

struct ReviewMap
{
    ReviewMap( int i, int j ):multi( NULL ){ coord[0] = i; coord[1] = j; };
    int coord[2];
    vector<int>* multi;
};

struct ReviewPair
{
    ReviewPair( int i, int j ){ coord[0] = i; coord[1] = j; dupe = 0; };
    static vector<ReviewPair> add( ReviewMap* l, ReviewMap* r );
    static int add( vector<ReviewPair> pairs[4], ReviewMap* l, ReviewMap* r, int est );
    int dist();
    int getCut( int est );
    int getDiff( int est );
    ReadId id;
    int coord[2], dupe, para, homo; // dupe 0 = unduped, 1 = left, 2 = right, 3 = para, 4 = unpara, 5 = homo, 6 = unhomo, 7 = homopara, 8 = unhomopara
    bool dispute;
};

struct ReviewAlign
{
    static vector<ReviewAlign> align( string& a, string& b );
    static bool conflict( ReviewAlign& l, ReviewAlign& r );
    static bool queryable( string& q );
    int coords[2][2], len;
};

struct ReviewStretch
{
    int coords[2], add;
};

class Reviewer
{
    Reviewer( Querier& bwt, string s );
    void add( ReadId id, int32_t i, int32_t j );
    ReviewMap* get( ReadId id );
    void setAllelic( Reviewer* alt );
    static string getDecimal( float dec, int places );
    void setCover( vector<Reviewer*>& reviews );
    static void setFoci( vector<Reviewer*>& reviews, vector<string>& foci );
    void setGaps( string l, string r );
    static vector<int> setLibs( vector<Reviewer*>& reviews );
    static float setMedian( vector<Reviewer*>& reviews );
    void setMultis( vector<Reviewer*>& reviews );
    void setOut( Reviewer* allele, vector<int>& libs, string ofn );
    void setOut( vector<int>& libs, ofstream* ofs, float median, int spaced, int cut );
    void setPairs( vector<Reviewer*>& reviews );
    static void setReview( Querier& bwt, vector<Reviewer*>& reviews );
    void setReads( Querier& bwt );
    void setSeeds( vector<string> seeds[2] );
    void setSpan( vector<ReviewPair>& pairs, int base, int lib, int& cur );
    void setSpanned( vector<int>& counts, int l, int r );
    static void dumpReads( string ifn, vector<Reviewer*>& reviews );
    static void loadReads( string ifn, vector<Reviewer*>& reviews );
    
    string seq_;
    Reviewer* allele_;
    vector<Node*> path_;
    unordered_map<ReadId, ReviewMap> mapped_;
    vector<ReviewStretch> stretch_;
    vector< vector<ReviewPair> > pairs_;
    vector< vector<ReadId> > mispair_[2], unpair_[2];
    vector< vector< pair<ReadId, int> > > inverts_[2][2], shorts_[2][2], longs_[2][2];
    vector< pair<int, int> > foci_;
    vector< pair<int, int> > seeds_;
    vector< pair<int, int> > homology_[2];   // [0] = homologue, [1] = paralogue
    vector<int> coords_;
    vector<int> counts_;
    double* cover_,* multi_;
    int len_, base_, end_;
public:
    static void review( Querier& bwt, vector< vector<Node*> >& paths );
    static void review( Querier& bwt, vector<string>& seqs );
    static void review( Querier& bwt, vector<string>& seqs, vector<string>& foci, string ofn, bool diploid );
};

#endif /* REVIEW_H */

