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

#ifndef CORRECT_STRUCTS_H
#define CORRECT_STRUCTS_H

#include "types.h"
#include <fstream>

struct Fastq
{
    Fastq( string line, uint8_t qualCutoff );
    bool getSeq( string &seq, int nCoords[2], int qCoords[2], int i );
    bool setNext();
    
    ifstream fp;
    string seqs[2], quals[2];
    uint8_t cutoff;
    bool paired;
};


#endif /* CORRECT_STRUCTS_H */

