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

struct CorrectReader
{
    CorrectReader( string &fnReads, uint8_t* used );
    bool read( string (&lines)[2][4] );
    FILE* ifpReads;
    ifstream ifsReads;
    uint64_t id;
    uint8_t* used;
    uint8_t idTable[8];
};

#endif /* CORRECT_STRUCTS_H */

