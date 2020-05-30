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

#ifndef QUERY_STRUCTS_H
#define QUERY_STRUCTS_H

#include "types.h"

struct ReadStruct
{
    string seq;
    int32_t tether[2], coords[2];
    ReadId id;
};

struct MappedSeqs
{
    MappedSeqs(){};
    void cull( vector<MappedSeqs>& alts );
    void setBest( string &seq );
    void sort();
    void updateTethers( string &seq );
    vector<int32_t> chunks;
    vector<ReadStruct> reads;
    unordered_set<ReadId> usedIds;
};

struct CorrectionRead
{
    CorrectionRead( ReadId id ) : id( id ){};
    string operator[]( int i ){ return i ? seq.substr( seq.length() - exts[1] ) : seq.substr( 0, exts[0] ); };
    bool congruent( string &ext, bool drxn );
    string seq;
    ReadId id;
    int exts[2], ol;
};

struct CorrectionExt
{
    CorrectionExt( string seq, CorrectionRead &read, bool drxn );
    CorrectionExt( string seq, vector<CorrectionRead> &reads, int i, bool drxn );
    bool addRead( CorrectionRead &read, bool drxn );
    static string getCongruent( vector<CorrectionExt> &exts, bool drxn );
    
    string seq;
    vector<ReadId> ids;
    int lens[2];
};

struct CorrectionStruct
{
    CorrectionStruct():overabundant( false ), error( false ), fork( false ){};
    void addReads( int i, bool drxn );
    void clear();
    void setExts( bool drxn );
    int validate( int len, int kmerLen, bool drxn );
    vector<CorrectionRead> reads;
    vector<CorrectionExt> exts[2];
    bool overabundant, error, fork;
};

struct QueryBranch
{
    QueryBranch( CharId rank, CharId count, int i, int ol ): rank( rank ), count( count ), i ( i ), ol( ol ){};
    CharId rank, count;
    int i, ol;
};

#endif /* QUERY_STRUCTS_H */

