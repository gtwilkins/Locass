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

#include "correct_structs.h"
#include <cassert>
#include "constants.h"

CorrectReader::CorrectReader( string &fnReads, uint8_t* used )
: ifsReads( fnReads ), used( used )
{
    ifpReads = NULL;
    id = 0;
    for ( int i = 0; i < 8; i++ ) idTable[i] = 1 << ( 7 - i );
}

bool CorrectReader::read( string (&lines)[2][4] )
{
    for ( int i = 0; i < 2; i++ ) for ( int j = 0; j < 4; j++ ) if ( !getline( ifsReads, lines[i][j] ) ) return false;
    if ( lines[0][0].empty() || lines[1][0].empty() ) return false;
    
    if ( !ifpReads ) return true;
    
    uint8_t lens[2], seqs[2][256];
    fread( lens, 1, 2, ifpReads );
    fread( seqs[0], 1, lens[0], ifpReads );
    fread( seqs[1], 1, lens[1], ifpReads );
    
    bool dump = used && ( used[id/8] & ( idTable[id%8] ) );
    id++;
}
