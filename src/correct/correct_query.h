/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   correct_query.h
 * Author: glen
 *
 * Created on 26 November 2018, 8:30 PM
 */

#ifndef CORRECT_QUERY_H
#define CORRECT_QUERY_H

#include "index_reader.h"
#include "types.h"
#include "query_binary.h"

struct CorrectExt
{
    CorrectExt( string &seq, int len, bool drxn );
    bool add( string &seq, bool drxn, bool doAdd=true );
//    static vector<CorrectExt> get( vector<Overlap> ols, int len, bool drxn );
    static vector<string> get( vector<Overlap> ols, int len, bool drxn );
    
    string ext;
    int maxLen;
    vector<int> lens;
};

struct CorrectAlign
{
    CorrectAlign( string &t, string &q );
    
    bool contendFull( CorrectAlign &full );
    bool contendPart( CorrectAlign &part, int &partLen );
    static vector<CorrectAlign> get( vector<string> &ts, vector<string> &qs, bool drxn );
    bool isBad();
    bool isFull( int maxMiss );
    bool isPart();
    bool operator >( CorrectAlign &rhs );
    void scoreGap( int d, int i, int gap );
    void write( string &seq, int len, bool drxn );
    
    static uint64_t nCorrectCount, errorCorrectCount, polyCorrectCount;
    string a[2];
    int lens[2], hit, miss, limit, poly, perf, score, maxScore;
};

struct CorrectBranch
{
    CorrectBranch( uint8_t i, int it, CharId rank, CharId count );
    void collect( vector<uint8_t> &ols, vector<CharId> &ranks, vector<CharId> &counts );
    bool beats( CorrectBranch &rhs, int len );
    int getNovel();
    bool steal( vector<uint8_t> &ols, vector<CharId> &ranks, vector<CharId> &counts, int len );
    string yield( bool drxn );
    static vector<string> yield( vector<CorrectBranch> &cbs, bool drxn );
    vector<CorrectBranch> branches;
    vector<uint8_t> q, endOverlaps;
    vector<CharId> endRanks, endCounts;
    CharId rank, count;
    int base, bad, good;
};

class CorrectQuery
{
public:
    CorrectQuery( IndexReader* ir, string &seq, int qLen, bool retract );
    CorrectQuery( IndexReader* ir, string &seq, bool drxn );
    CorrectQuery( IndexReader* ir, vector<uint8_t> &q, vector<CorrectBranch> &branches, int seqLen );
    
//    void correct( QueryBinaries* qb, string &seq, bool drxn );
    int correct( QueryBinaries* qb, string &seq, bool &trimmed );
    int correct( string &seq, vector<string> &seqs, bool &trimmed );
    int trim( QueryBinaries* qb, string &seq, bool &trimmed );
    int trim( QueryBinaries* qb, string &seq, vector<string> &seqs, bool &trimmed );
    
private:
    int correct( vector<CorrectAlign> &cas, string &seq, int maxMiss, bool partial, bool drxn );
    void query( uint8_t i, int it, CharId rank, CharId count );
    void query( CorrectBranch &cb, uint8_t i, int it );
    
    IndexReader* ir_;
    vector<CorrectBranch> branches_;
    vector<uint8_t> q_, endOverlaps_;
    vector<CharId> endRanks_, endCounts_;
    int base_, qLen_, seqLen_;
    bool collect_, trim_;
};


#endif /* CORRECT_QUERY_H */

