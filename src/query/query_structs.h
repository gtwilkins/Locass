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
    SeqNum readId;
};

struct MappedSeqs
{
    MappedSeqs(){};
    void setBest( string &seq );
    void sort();
    vector<int32_t> chunks;
    vector<ReadStruct> reads;
    unordered_set<SeqNum> usedIds;
};

#endif /* QUERY_STRUCTS_H */

