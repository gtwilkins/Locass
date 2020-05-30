/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   damplify_structs.h
 * Author: glen
 *
 * Created on 7 November 2018, 3:37 PM
 */

#ifndef DAMPLIFY_STRUCTS_H
#define DAMPLIFY_STRUCTS_H

#include "types.h"

struct DeamplifyFiles
{
    DeamplifyFiles( string &ifnSeq, string &ifnStatus );
    ~DeamplifyFiles();
    
    void collate( FILE* ofpDump, string &fnKmers, uint64_t kmerSize, uint64_t maxMem, uint64_t &ampliconCount );
    void cycle( string &fnKmers, uint64_t kmerSize, bool &full );
    void edit( uint32_t id, uint8_t status );
    void edit( uint8_t status );
    void overwrite( string (&lines)[2][4] );
    bool read( uint8_t &status, uint64_t &kmer );
    int readSize();
    void reset( int j );
    void rewrite();
    void rewind();
    bool set( int j );
//    void test( uint64_t testKmer );
    
    FILE* ifpSeq,* ifpStatus,* ifpEdit;
    uint64_t kmersRead, kmersUnread, setPos[2], lineLen;
    uint32_t curPair, maxPair, setPair[2];
    uint8_t* bufStatus, bufSeq[1024], lens[2], mask[4], inTable[4][256], outTable[4][4], readLen;
    int bufSize, p, i;
    bool doRewrite;
};


#endif /* DAMPLIFY_STRUCTS_H */

