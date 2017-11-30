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

#ifndef QUERY_EXTENSION_H
#define QUERY_EXTENSION_H

#include "types.h"
#include "parameters.h"

struct Overlap
{
    Overlap( SeqNum readId, uint16_t overlapLen, uint16_t extendLen ) : readId( readId ), overLen( overlapLen ), extLen( extendLen ), redundant( false ){}
    void offset( int offset, bool drxn );
    void truncate( bool drxn );
    SeqNum readId;
    string seq;
    uint16_t overLen, extLen;
    bool redundant;
};

struct Extension {
    Extension( Overlap &overlap, int readCount, bool drxn );
    Extension( Overlap &overlap, int readCount, bool drxn, vector<Extension> &fwdExts, bool doAdd );
    
    virtual ~Extension() {};
    
//    void addFwd( vector<Extension> &addFwdExts, bool drxn );
    void addOverlap( Overlap &overlap );
    bool canAdd( Overlap &overlap );
    bool checkLoop( Overlap &overlap );
    bool getExtendSum() const;
    bool getOverlapSum() const;
    bool isCongruent( Overlap &overlap );
    bool isRepeat( string &seq);
    bool isSuperior( Extension &ext, uint16_t cutoff );
    bool isSuperior( int minOver, int sumCutoff, int maxCount, bool valid=true );
//    bool addGuide( Overlap &overlap, bool drxn );
    bool isValid();
    void offsetFwdExts( int offset );
    bool operator<( const Extension &rhs ) const;
    void resetMaxOverLen();
    void setFullSeq( string &inSeq );
    void setScores( int minOver, int overCount );
    void trimToExtLen( uint16_t extCutoff );
    bool trimToOverlap( uint16_t &overlapCutoff );
    bool trimToSeeds( vector<Overlap> &seedOverlaps, unordered_set<SeqNum> &seeds );
    
    string seq;
    string fullSeq;
    uint16_t maxOverLen, maxExtLen, score;
    int readCount;
    vector<Overlap> overlaps; // sorted by descending extendLen_
    vector<Extension> fwdExts;
    bool drxn, valid, anyPe;
};


#endif /* QUERY_EXTENSION_H */

