/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   review.h
 * Author: glen
 *
 * Created on 24 January 2020, 5:19 AM
 */

#ifndef REVIEW_H
#define REVIEW_H

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

class Review
{
    Review( Querier& bwt, string s );
    void add( ReadId id, int32_t i, int32_t j );
    ReviewMap* get( ReadId id );
    static string getDecimal( float dec, int places );
    void setCover( vector<Review*>& reviews );
    void setGaps( string l, string r );
    static vector<int> setLibs( vector<Review*>& reviews );
    static float setMedian( vector<Review*>& reviews );
    void setMultis( vector<Review*>& reviews );
    void setOut( vector<int>& libs, ofstream* ofs, float median, int spaced, int cut );
    void setPairs( Querier& bwt, vector<Review*>& reviews );
    static void setReview( Querier& bwt, vector<Review*>& reviews );
    void setSeeds( vector<string> seeds[2] );
    void setSpan( vector<ReviewPair>& pairs, int base, int lib, int& cur );
    void setSpanned( vector<int>& counts, int l, int r );
    
    string seq_;
    Review* allele_;
    vector<Node*> path_;
    unordered_map<ReadId, ReviewMap> mapped_;
    vector< vector<ReviewPair> > pairs_;
    vector< vector<ReadId> > mispair_[2], unpair_[2];
    vector< vector< pair<ReadId, int> > > inverts_[2][2], shorts_[2][2], longs_[2][2];
    vector< pair<int, int> > seeds_;
    vector< pair<int, int> > homology_[2];   // [0] = homologue, [1] = paralogue
    vector<int> coords_;
    vector<int> counts_;
    float* cover_,* multi_;
    int len_;
public:
    static void review( Querier& bwt, vector< vector<Node*> >& paths );
    static void review( Querier& bwt, vector<string>& seqs );
};

#endif /* REVIEW_H */

