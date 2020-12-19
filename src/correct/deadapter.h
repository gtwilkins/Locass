/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   deadaptor_identify.h
 * Author: glen
 *
 * Created on 2 November 2018, 3:37 PM
 */

#ifndef DEADAPTOR_H
#define DEADAPTOR_H

#include "types.h"

class Deadapter
{
public:
    Deadapter( string fn );
    bool isConnected( string &seq, string &phred );
    bool isOverlap( string &seq1, string &seq2, string &phred1, string &phred2, bool blank );
    
private:
    int align( string (&s)[2], string (&p)[2], int &hits, int &miss, int &diff, int trim );
    bool addEnds( vector< pair<string, uint32_t> > &ends, string &s );
    bool consolidateEnd( vector< pair<string, uint32_t> > &ends, int &len, string &adaptor );
    bool getPair( string (&s)[2] );
    int getDifference( string (&s)[2], int len[2] );
    bool isOverlap2( string &seq1, string &seq2, string &phred1, string &phred2, int i, bool blank );
    bool isOverlap3( string &seq, string &phred, int ol, int i );
    bool isOverlap4( string &seq, string &phred, int i );
    bool readLine( string &s, bool suppress );
    void rewind();
    void setAdapters( uint32_t kmerLen );
    bool setBlank( string &seq, string &phred, int start );
    bool setConnectors( vector<string> &connectors );
    
    string fn_, adapter_[2];
    vector<string> connectors_;
    ifstream ifs_;
    uint64_t lineCount_;
    int header_, footer_;
};

#endif /* DEADAPTOR_IDENTIFY_H */

