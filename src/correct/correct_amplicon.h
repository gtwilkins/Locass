/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   correct_amplicon.h
 * Author: glen
 *
 * Created on 29 May 2018, 1:58 PM
 */

#ifndef CORRECT_AMPLICON_H
#define CORRECT_AMPLICON_H

#include "types.h"
#include "correct_query.h"
#include "deadapter.h"

class AmpPile;

struct Amp
{
    Amp( uint32_t id, uint8_t len, uint8_t* buff );
    ~Amp();
    
    bool align( Amp* amp );
    uint8_t comp( int i );
    static uint64_t create( FILE* ifp, Amp* (&amps)[2] );
    int penalty( int i );
    
    AmpPile* pile;
    Amp* alt;
    uint32_t id;
    string ascii;
    vector<uint8_t> seq, phred;
    int phredLen, queryLen;
};

class AmpPile
{
public:
    AmpPile( vector<AmpPile*> &piles, vector<Amp*> &amps, IndexReader* ir );
    AmpPile( vector<AmpPile*> &piles, AmpPile* pile, int i, CharId rank, CharId count, IndexReader* ir );
    AmpPile( AmpPile* basis, AmpPile* pair );
    ~AmpPile();
    static vector<AmpPile*> compile( vector<Amp*> (&amps)[2], IndexReader* ir, QueryBinaries* qb );
    void output( ofstream &ofs, uint8_t* used, Deadapter* da, IndexReader* ir, QueryBinaries* qb );
    static uint64_t pileCount, ampCount, adpCount, trimCount;
private:
    void add( Amp* amp );
    void advance( uint8_t c, CharId rank, CharId count );
    void cannabalise( AmpPile* pile );
    static void clean( vector<AmpPile*> &piles );
    static void collapse( vector<AmpPile*> &piles );
    void confirm( vector<AmpPile*> &alts );
    void demate( IndexReader* ir, QueryBinaries* qb );
    void fill( int trimLen );
    Amp* getBest();
    int getDisagree( AmpPile* alt );
    void getPairs( vector<AmpPile*> &piles );
    string getPhred();
    int getPhred( int i );
    static AmpPile* getPile( Amp* amp, vector<AmpPile*> &piles );
    vector<string> getSeqs();
    bool inPile( Amp* amp );
    void mate();
    bool merge( AmpPile* alt );
    void query( vector<AmpPile*> &piles, IndexReader* ir  );
    void take( Amp* amp, vector<Amp*> &amps );
    
    vector<Amp*> amps_;
    vector<CorrectBranch> branches_;
    vector<uint8_t> seq_;
    string ascii_;
    AmpPile* pair_;
    CharId rank_, count_;
    int seqLen_, goodLen_;
};

//struct Amplicon
//{
//    Amplicon( uint8_t* buff, Amplicon* ampPair );
//    
//    void addScore( vector<uint8_t> &seq, vector<int> &scores );
//    static bool allSame( vector<Amplicon*> &amps, int i );
//    int getCharScore( int i );
//    int getPairScore();
////    int getRunScore( int i, int last );
//    void shift( uint8_t* bases, int base, int off );
//    
//    Amplicon* pair_;
//    string ascii_;
//    vector<uint8_t> seq_, phred_;
//    uint32_t id_;
//    int score_, strength_, miss_, perfect_, goodLen_;
//};
//
//struct Offset
//{
//    Offset(): off( 0 ){};
//    Offset( int off, Amplicon* amp ): off( off ){ amps.push_back( amp ); };
//    bool clean( vector<Amplicon*> &amp );
//    static bool set( vector<Offset> &offs, vector<Amplicon*> &good );
//    int off;
//    vector<Amplicon*> amps;
//};
//
//struct Shift
//{
//    Shift( int i, Amplicon* amp, int off ): pos( i ){ offs.push_back( Offset( off, amp ) ); }
//    bool add( int i, Amplicon* amp, int off );
//    int resolve( vector<uint8_t> &seq, vector<Amplicon*> &amplicons, int off, int &posLimit );
//    
//    int pos;
//    vector<Offset> offs;
//};
//
//class Pile
//{
//public:
//    
//    static vector<string> compile( vector<Amplicon*> (&amps)[2] );
//    
//    static int adapterCount_;
//    static uint32_t pileCount_, ampCount_;
//private:
//    Pile( Amplicon* amp );
//    Pile( vector<Amplicon*> &amps, int len );
//    ~Pile();
//    
//    void add( vector<Amplicon*> &amps, bool doSteal );
//    bool align( Amplicon* amp, int maxMiss, bool steal, bool force, bool add );
//    int alignExact( Amplicon* amp, int &score, int &miss, int &perfect );
//    bool alignInexact( int i, Amplicon* amp, int &score, int &miss, int &perfect, int maxMiss, vector< pair<int, int> > &offs );
//    void build();
//    bool close();
//    bool confirm();
//    bool couple( vector<Pile*> &piles );
//    bool couple( vector<Pile*> &curPiles, vector<Pile*> &altPiles );
//    bool disentagle( vector<Pile*> &altPiles, vector<Amplicon*> &curAmps, vector<Amplicon*> &altAmps );
//    void dismantle( vector<Amplicon*> &amps );
//    void extract( vector<Amplicon*> &curAmps, vector<Amplicon*> &altAmps );
//    bool fill( vector<Amplicon*> &amps );
//    static vector<Pile*> group( vector<Amplicon*> &amps );
//    static void group( vector<Pile*> &piles, vector<Amplicon*> &amps, int goodLen, int i );
//    bool inPile( Amplicon* amp );
//    bool merge( Pile* pile );
//    static void optimise( vector<Pile*> &piles );
//    bool optimise();
//    string output();
//    void output( string &seq, string &phred );
//    bool rebuild();
//    void scoreConsensus( double scores[][5], int len );
//    void setNextOffset( Amplicon* amp, int &pos, int &off );
//    void setPerfect();
////    bool shift( Amplicon* amp, int &iCurr, int &jCurr );
//    
//    vector<Amplicon*> amps_;
//    vector<uint8_t> seq_;
//    vector<Shift> shifts_;
//    Pile* pair_;
//    int score_, miss_, perfect_;
//    bool optimised_;
//};

#endif /* CORRECT_AMPLICON_H */

