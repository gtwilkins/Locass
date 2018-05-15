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


class Correct
{
public:
    Correct( int argc, char** argv );
    
    void correct( Querier &bwt, ifstream &infile, string oprefix );
    void correct();
    void decontaminate();
    void trim( vector<string> &pes, vector<string> &mps);
    
private:
    bool correctRead( Querier &bwt, string &qSeq, string &finalSeq, int nCoords[2], int qCoords[2] );
    void trimPair( string &seq, string &rev );
    bool trimPrimers( string lines[4], string extra[2], bool isFastq, bool isMp );
    void writeFile( ofstream &fp, string lines[4], int lineCount );
};

#endif /* CORRECT_H */

