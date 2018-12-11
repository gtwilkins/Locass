/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   correct.h
 * Author: glen
 *
 * Created on 19 January 2018, 8:43 AM
 */

#ifndef CORRECT_H
#define CORRECT_H

#include "types.h"
#include "query.h"
#include "query_binary.h"
#include "correct_kmers.h"
#include "correct_read.h"
#include "deadapter.h"
#include "query_state.h"
#include "correct_query.h"
#include "deamplify_structs.h"


class Correct
{
public:
    Correct( int argc, char** argv );
    ~Correct();
    
private:
    void correct();
    void deamplify();
    void deamplifyConvert();
    void deamplifyIdentify();
    void deamplifyScrutinise();
    
    
    bool readFastq( ifstream &ifs, string (&lines)[4] );
//    void censorFiles( string fn, ofstream &ofsOut, bool paired, bool chim, bool pyro );
//    void correct( Querier &bwt, string fn, bool paired, bool chim, bool pyro );
//    bool readSequence( ifstream &ifs, MarksFile &mf, string* lines, vector<Interval*> &ils, bool chim, bool pyro );
//    bool splitJoin( string &seq1, string &seq2, string &phred1, string &phred2 );
//    void trimEnd( string &seq, string &phred );
    
    QueryBinaries* qb_;
    IndexReader* ir_;
    Deadapter* da_;
    Filenames* fns_;
    DeamplifyFiles* df_;
    ofstream ofs_[2];
    uint8_t* used_;
    int64_t pairCount_, trimCount_, adpCount_, discardCount_;
    uint8_t readLen_;
    int32_t maxSeqMem_;
    string ifn_, ofns_[6], fnBase_, fnSeqs_, joiner_, ender_;
    ifstream ifsReads_;
    ofstream ofsIntervals_, ofsPairs_, ofsSingles_;
    uint64_t byteCount_, dontAddByte_;
    bool demate_, deadapter_, deamplify_;
};

#endif /* CORRECT_H */

