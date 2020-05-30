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

class Pile;

struct Amp
{
    Amp( uint32_t id, uint8_t len, uint8_t* buff );
    ~Amp();
    
    bool align( Amp* amp );
    uint8_t comp( int i );
    static uint64_t create( FILE* ifp, Amp* (&amps)[2] );
    int penalty( int i );
    void setUsed( uint8_t* used );
    
    static uint64_t ampCount;
    Pile* pile;
    Amp* alt;
    uint32_t id;
    string ascii;
    vector<uint8_t> seq, phred;
    int phredLen, queryLen;
};

class Pile
{
public:
    Pile( vector<Pile*> &piles, vector<Amp*> &amps, IndexReader* ir );
    Pile( vector<Pile*> &piles, Pile* pile, int i, CharId rank, CharId count, IndexReader* ir );
    Pile( Pile* basis, Pile* pair );
    ~Pile();
    static vector<Pile*> compile( vector<Amp*> (&amps)[2], IndexReader* ir );
    void output( ofstream &ofs, uint8_t* used, Deadapter* da, IndexReader* ir );
    static uint64_t pileCount, adpCount, trimCount, created, deleted;
private:
    void add( Amp* amp );
    void advance( uint8_t c, CharId rank, CharId count );
    void cannabalise( Pile* pile );
    static void clean( vector<Pile*> &piles );
    static void collapse( vector<Pile*> &piles );
    void confirm( vector<Pile*> &alts );
    void demate( IndexReader* ir );
    void fill( int trimLen );
    Amp* getBest();
    int getDisagree( Pile* alt );
    void getPairs( vector<Pile*> &piles );
    string getPhred();
    int getPhred( int i );
    static Pile* getPile( Amp* amp, vector<Pile*> &piles );
    vector<string> getSeqs();
    bool inPile( Amp* amp );
    void mate();
    bool merge( Pile* alt );
    void query( vector<Pile*> &piles, IndexReader* ir  );
    void take( Amp* amp, vector<Amp*> &amps );
    
    vector<Amp*> amps_;
    vector<CorrectBranch> branches_;
    vector<uint8_t> q_;
    string seq_;
    Pile* pair_;
    CharId rank_, count_;
    int seqLen_, goodLen_;
};

#endif /* CORRECT_AMPLICON_H */

