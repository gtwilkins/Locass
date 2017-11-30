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

#ifndef CALIBRATE_STRUCTS_H
#define CALIBRATE_STRUCTS_H

#include "types.h"
#include "parameters.h"

struct SeedLibraryCount
{
    SeedLibraryCount( Parameters &inParams );
    void add( ReadId id, int dist, bool leftFwd, bool rightFwd );
    
    Parameters &params;
    vector<ReadId> ffCounts, frCounts, rfCounts, libCounts;
    vector<CharId> ffSum, frSum, rfSum;
    vector<float> ratios;
    float ratioScore;
    int id;
};

struct SeedLibraryCounts
{
    SeedLibraryCounts( Parameters &params ): libRatios( params.libs.size() ){}
    vector<int> set( Parameters &params );
    vector<SeedLibraryCount> libs;
    vector<ReadId> ffCounts, frCounts, rfCounts;
    vector<CharId> ffSum, frSum, rfSum;
    vector< vector<float> > libRatios;
};

struct LocusLibraryCount
{
    LocusLibraryCount( Parameters &params ): lens( params.libs.size()), ffPairs( params.libs.size() ), frPairs( params.libs.size() ), rfPairs( params.libs.size() ) {}
    int32_t getMedian( int i, bool isFf, bool isFr, bool isRf );
    void setMedians( Parameters &params );
    vector<int32_t> lens;
    vector< vector<int32_t> > ffPairs, frPairs, rfPairs;
    float coverage;
    int len;
};

struct LocusLibraryCounts
{
    void set( Parameters &params );
    
    vector<LocusLibraryCount> libs;
};

#endif /* CALIBRATE_STRUCTS_H */

