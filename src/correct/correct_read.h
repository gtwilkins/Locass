/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   correct_read.h
 * Author: glen
 *
 * Created on 9 July 2018, 2:50 PM
 */

#ifndef CORRECT_READ_H
#define CORRECT_READ_H

#include "types.h"
#include "query.h"

struct Chunk
{
    Chunk( string &s, int coords[2], bool inCaps[2] );
    Chunk( string &s, int coords[2], bool inCaps[2], int confLen, bool bridged, bool d );
    bool match( ReadId id, string &q, int cutoff=0, bool doCutoff=false );
    bool matchPartial( ReadId id, string &q );
    static void reverse( vector<Chunk> &cs, int len );
    void setLimits( int* limits );
    
    string seq;
    unordered_map<ReadId, pair<int,int> > hits;
    int ends[2], conf[2], drxn;
    bool caps[2];
};

struct AlignStruct
{
    AlignStruct( CorrectionExt &ext, string &t, string &p, bool drxn );
    bool addLen( int len );
    void align( string &t, string &p );
    int charCount( char &c, int &i );
    bool charCount( char &c, int &i, int &copies );
//    bool confirm( bool initial, bool pyro, bool bridge );
    bool consolidate( AlignStruct &as, int polyLen );
//    void fold( vector<AlignStruct> (&alts)[2], unordered_map<ReadId, CorrectionRead> &reads, unordered_set<ReadId> &foldIds, bool drxn );
    int getAlignLimit( bool i );
    static AlignStruct* getExtension( AlignStruct*& best, vector<AlignStruct> &ass, int gapLen, int polyLen, bool initial );
    int getExtLen();
    int getNovelty( vector<AlignStruct>* ass );
    string getPhred( string &p, int seqLen, bool drxn );
    int getSharedIdCount( AlignStruct &as );
    string getSeq( int seqLen, bool drxn );
    int getValidGap( int gapLen );
    bool isBetter( AlignStruct &as, bool initial=false );
    bool isBetterErase( AlignStruct &as );
    void setCongruence( AlignStruct &as, int &congruence );
    bool setExtLens( int gapLen, int extLen );
    void setHits( vector<Chunk> &t, vector<bool> &tHits, int coord );
    bool setReads( vector<AlignStruct> &ass, unordered_map<ReadId, CorrectionRead> &reads );
    bool setScore( Chunk &q, vector<Chunk> (&t)[3], unordered_map<ReadId, CorrectionRead> &reads, int gap );
    string a[2], seq;
    vector<ReadId> ids;
    unordered_set<ReadId> hitIds;
    vector< pair<int,int> > errors;
    vector<bool> matched;
    int aligns[3], ends[3], lens[3], misses, polys, novelty, hitLen;
    bool bridged, valid, drxn;
};

struct Interval;
struct IntervalTargets
{
    IntervalTargets( vector<Interval*> &ils, Querier &bwt, unordered_map<ReadId, string> &reads ): ils( ils ), bwt( bwt ), reads( reads ), added( false ){};
    void addInterval( Interval* il, bool drxn );
    string& getRead( ReadId id );
    void query( vector<Chunk> &cs, unordered_set<ReadId> &ids, int cutoff=0, bool doCutoff=false );
    vector<Interval*> &ils;
    Querier &bwt;
    vector<Chunk> t[3];
    unordered_set<ReadId> ids;
    unordered_map<ReadId, string> &reads;
    bool added;
};

struct Interval
{
    Interval( int coord, bool good, bool pyro );
    Interval( string &s, string &p, int first, int last, Interval* edge, bool good, bool pyro );
    static bool clean( vector<Interval*> &ils );
    bool confirm( AlignStruct* best, bool initial, bool drxn );
    static void correct( Querier &bwt, vector<Interval*> (&ils)[4], bool pe, bool chim );
    bool correct( IntervalTargets &its, bool initial );
    bool correct( IntervalTargets &its, bool initial, bool drxn );
    bool correct( IntervalTargets &its, AlignStruct* best, bool drxn );
    static vector<Interval*> create( string &s, string &p, vector<bool> &marks, bool doSplit=false, bool isPyro=false );
    bool extend( IntervalTargets &its, AlignStruct* best, int minExt, bool drxn );
    void fold( AlignStruct* best, bool drxn );
    bool getAlignTargets( string &s, string &p, bool drxn );
    string getPhred( AlignStruct &as, bool drxn );
    string getSequence( AlignStruct &as, bool drxn );
    static Interval* getSplit( vector<Interval*> &ils );
    static void getTargets( IntervalTargets &its, vector<Interval*> &ils );
    void getTargets( unordered_set<ReadId> &ids, vector<Chunk> &t, int drxns, bool drxn );
    int getUniqueLen( int uniqueCount, bool drxn );
    bool isCap( bool isAlt, bool drxn );
    void merge( IntervalTargets &its, bool drxn );
    static void output( vector<Interval*> &ils, ofstream &ofsSeqs, string* lines, bool pyro );
//    static void output( vector<Interval*> (&ils)[3], ofstream &ofsIntervals, ofstream &ofsPairs, ofstream &ofsSingles );
    bool query( IntervalTargets &its, bool drxn );
    void query( CorrectionStruct &cs, bool drxn );
    void realign( bool drxn );
    void shift( int diff );
    static void shrink( vector<Interval*> &ils, int maxLen );
    bool shrink( int diff, bool exchange, bool excise, bool drxn );
    void split( Querier &bwt, vector<Interval*> &ils, vector<Interval*> &splitIls );
    bool split( Querier &bwt, bool d );
    
    string seq, phred;
    vector<AlignStruct> alts[2];
    unordered_set<ReadId> usedIds;
    unordered_map<ReadId, CorrectionRead> reads;
    Interval* edges[2];
    int len, coords[2];
    bool good, bad, pyro, extable[2], appendable[2];
    
    static ReadId capCount, cappedCount, bridgeCount, bridgedCount;
};

#endif /* CORRECT_READ_H */

