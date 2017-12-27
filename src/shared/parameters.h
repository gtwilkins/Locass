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

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "types.h"

struct Lib {
    Lib() : size( 0 ), minDist( 0 ), maxDist( 0 ), orientation( 0 ), isPe( false ) {};
    bool doAddMarker( bool thisRev, int &pairRev, bool drxn );
    bool getPair( ReadId &readId, int32_t dist, int &pairDrxn );
    void setMinMax();
    
    int32_t size, minDist, maxDist;
    ReadId count, endCount;
    uint8_t orientation; // 0 = Unknown, 1= FF 2 = FR, 3 = RF
    bool isPe;
};

struct Parameters {
    Parameters();
    
    void checkReady();
    bool isReadPe( SeqNum &readId );
    int32_t getFurthestMpDist( int32_t coord, bool drxn );
    int32_t getFurthestPeDist( int32_t coord, bool drxn );
    int32_t getFurthestMpMean( int32_t coord, bool drxn );
    int32_t getFurthestPeMean( int32_t coord, bool drxn );
    Lib* getLib( const SeqNum &readId );
    int32_t getLibSize( SeqNum &readId );
    SeqNum getPairId( const SeqNum readId );
    void set();
    void setLimits( int32_t &limit );
    
    int readLen, queryMpCutoff;
    int32_t locusLimits[2];
    int32_t maxPeMax, maxPeMean, maxMpMax, maxMpMean, avgPeMean;
    vector<Lib> libs;
    float cover, readSpacing, branchMinHits, peCover, peRatio;
    SeqNum seqCount;
    uint8_t outMode;  // 0 = unitig, 1 = alleles
    bool isSet, isCalibrated;
    
};

#endif /* PARAMETERS_H */

