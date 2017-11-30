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

#ifndef INDEX_WRITER_H
#define INDEX_WRITER_H

#define IDX_BUFFER (CharId)16384

#include "filenames.h"
#include "types.h"

class IndexWriter
{
public:
    IndexWriter( PreprocessFiles* fns );
    virtual ~IndexWriter();
    
private:
    void write();
    
    FILE* bwt,* idx;
    CharId id;
    
    uint8_t* buff;
    ReadId* marks;
    
    uint8_t decodeBaseChar[256], decodeBaseRun[256], maxBaseRun[5];
    uint8_t contFlag, contMask;
    uint8_t bwtBegin;
    
    ReadId bwtPerIndex;
    
    ReadId basePos[4], countsPerMark, markSizes[5];
    CharId indexSize, markSize;
    CharId currByte;
    CharId bwtSize;
    CharId charCounts[5], counts[5];
};

#endif /* INDEX_WRITER_H */

