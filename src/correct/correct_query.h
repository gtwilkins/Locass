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
    void align( string seq, int len, bool drxn );
    void contend( CorrectAlign*& full, CorrectAlign*& part, int partLen, int imperf );
    void correct( CorrectAlign*& full, CorrectAlign*& part, bool &bad, int maxMiss );
    bool proceed( IndexReader* ir, int limit );
    void query( IndexReader* ir, int limit );
    
    vector<CorrectBranch> branches;
    vector<CorrectAlign> aligns;
    vector<uint8_t> q, endOverlaps;
    vector<CharId> endRanks, endCounts;
    CharId rank, count;
    int base, bad, good;
};

class CorrectQuery
{
public:
    CorrectQuery( IndexReader* ir, string &seq, int &len, bool &trimmed, bool dummy );
    CorrectQuery( IndexReader* ir, string &seq, vector<string> &seqs, vector<CorrectBranch> &branches, int &len, bool &trimmed );
    CorrectQuery( IndexReader* ir, string &seq, vector<string> &seqs, int &len, bool &trimmed );
    
private:
    int correct( IndexReader* ir, string &seq, bool drxn );
    int correct( IndexReader* ir, string &seq, vector<string> &seqs );
    int correct( string &seq, bool drxn );
    bool proceed( IndexReader* ir );
    ReadId query( IndexReader* ir, uint8_t i, int it, CharId rank, CharId count );
    void query( uint8_t i, int it, CharId rank, CharId count );
    bool setQuery( string &seq, bool drxn );
    
    vector<CorrectBranch> branches_;
    vector<uint8_t> q_, endOverlaps_;
    vector<CharId> endRanks_, endCounts_;
    int base_, qLen_, seqLen_, coords_[2];
    bool collect_, trim_, initial_;
};


#endif /* CORRECT_QUERY_H */

