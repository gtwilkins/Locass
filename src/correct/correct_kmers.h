/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   correct_kmers.h
 * Author: glen
 *
 * Created on 1 August 2018, 7:23 PM
 */

#ifndef CORRECT_KMERS_H
#define CORRECT_KMERS_H

#include "types.h"

//struct MarksFile
//{
//    MarksFile( string fn, string mode );
//    ~MarksFile();
//    void advance( size_t posAdv, int pAdv );
//    void edit( bool mark );
//    void editSkip( int count );
//    vector< pair<int,int> > getIntervals( int len );
//    void increment();
//    bool read();
//    vector<bool> read( int count );
//    void write( bool mark );
//    
//    FILE* fp;
//    size_t pos;
//    uint8_t b[4096], x[8];
//    int i, p;
//    bool out;
//};
//
//struct KmerFileStruct
//{
//    KmerFileStruct( string &fn, int worth, bool seed );
//    void clean();
//    void close();
//    void open( size_t posStart=0 );
//    bool read();
//    
//    FILE* fp;
//    string fnSeq, fnMark;
//    vector<CharId> kmers;
//    CharId charCount, rChar[4];
//    ReadId seqCount, seqsCounted;
//    int kmerCount;
//    uint16_t seqLen, lastEnd, il[2];
//    uint8_t ilCount, worth, rKmer[4][256];
//    bool seed;
//};
//
//class KmerCount
//{
//public:
//    KmerCount( vector<string> &filenames, string fnBase );
//    void catalog();
//private:
//    bool read();
//    void write();
//    
//    string fnBase_;
//    vector<KmerFileStruct> kfss_;
//    MarksFile* mf_[2];
//    size_t posKfs_[2], posMarks_[2];
//    int iKfs_[2], pMarks_[2];
//};


#endif /* CORRECT_KMERS_H */

