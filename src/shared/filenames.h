/*
 * Copyright (C) 2017 Glen T. Wilkins <glen.t.wilkins@gmail.com>
 * Written by Glen T. Wilkins
 * 
 * This file is part of the Locass software package <https://github.com/gtwilkins/Locass>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FILENAMES_H
#define FILENAMES_H

#include "types.h"
#include <fstream>

struct Filenames
{
    Filenames( string inPrefix );
    
    FILE* getReadPointer( string &filename, bool doEdit );
    FILE* getWritePointer( string &filename );
    ifstream getReadStream( string &filename );
    ofstream getWriteStream( string &filename );
    
    FILE* getBinary( bool doRead, bool doEdit );
    void removeFile( string &filename );
    void setIndex( FILE* &inBin, FILE* &inBwt, FILE* &inIdx );
    
    string prefix;
    string bin;
    string bwt;
    string ids;
    string idx;
};

struct PreprocessFiles : public Filenames
{
    PreprocessFiles( string inPrefix );
    
    void clean();
    void setBinaryWrite( FILE* &outBin, FILE* &outBwt, FILE* &outEnd, FILE* (&outPos)[4], FILE* (&outSap)[4], FILE* (&outIds)[4][4] );
    void setCycler( FILE* &inBwt, FILE* &outBwt, FILE* &inEnd, FILE* &outEnd, FILE* (&outPos)[4], FILE* (&outSap)[4], FILE* (&outIds)[4][5], uint8_t cycle );
    void setCyclerIter( FILE* &inPos, FILE* &inSap, FILE* (&inIds)[5], uint8_t cycle, uint8_t i );
    void setCyclerFinal( FILE* &inBwt, FILE* &outBwt, FILE* &inEnd, FILE* &outEnd, uint8_t cycle );
    void setCyclerFinalIter( FILE* &inPos, FILE* &inSap, FILE* &inIds, uint8_t cycle, uint8_t i );
    void setCyclerUpdate( FILE* &outBwt,  FILE* &outEnd, FILE* (&outPos)[4], FILE* (&outSap)[4], FILE* (&outIds)[4][5], uint8_t cycle );
    void setIndexWrite( FILE* &inBwt, FILE* &outIdx );
    
    string tmpBwt[2];
    string tmpEnd[2];
    string tmpPos[2][4];
    string tmpSap[2][4];
    string tmpIds[2][4][5];
    string tmpSingles;
};


#endif /* FILENAMES_H */

