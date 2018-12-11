/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   local_alignment.h
 * Author: glen
 *
 * Created on 2 July 2018, 11:58 PM
 */

#ifndef LOCAL_ALIGNMENT_H
#define LOCAL_ALIGNMENT_H

#include <vector>
#include <string>

class LocalAlignment
{
public:
    LocalAlignment( std::string &a, std::string &b, bool glocal, bool freePolymer );
    
    static int isGapPoly( std::string (&a)[2], int d, int i, int gap );
    void print( int i, int j );
    void realign( std::string &a, std::string &b, bool conform, bool bluntStart=false, bool bluntEnd=false, bool trimEnd=false );
//    void setAlign( std::string &a, std::string &b );
//    int setCoords( int* coords );
private:
    void blunt( bool start, bool end );
    int score( int i, int j, char &c, int* coords );
    void setRuns( int* run, int i, int j );
    
    std::string& a_,& b_;
    std::vector< std::vector<int> > m_;
    bool freePolymer_;
};



#endif /* LOCAL_ALIGNMENT_H */

