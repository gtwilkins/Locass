/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   correct_bitfiles.h
 * Author: glen
 *
 * Created on 19 November 2018, 6:54 PM
 */

#ifndef BITFILES_H
#define BITFILES_H

#include "types.h"

class BitFile
{
public:
    BitFile( string fn, int bits );
    ~BitFile();
    
    void close( bool ignore=false );
    void edit( uint8_t code );
    bool open( bool doRead, bool doWrite, bool overwrite=true );
    uint8_t read();
    void setPos( uint64_t pos, int i );
    void tellPos( uint64_t &pos, int &i  );
    void write( uint8_t code );
    
    string fn_;
    FILE* fp_;
    uint8_t flags_[8], readTable_[8][256], writeTable_[8][16], byte_;
    int i_, bits_, per_;
    bool edited_, writing_;
};



#endif /* CORRECT_BITFILES_H */

