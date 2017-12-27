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

#include "shared_functions.h"
#include <cassert>

int getHomopolymerLen( string &s, bool drxn )
{
    int i = 0;
    while ( i+1 < s.length() &&
            ( drxn ? s[s.length()-i-1] == s[s.length()-i-2] : s[i] == s[i+1] ) ) i++;
    if ( i+2 < s.length() && ( drxn ? s[s.length()-i-2] == s[s.length()-i-3] : s[i+1] == s[i+2] ) )
    {
        int j = ++i;
        while ( i+1 < s.length() &&
                ( drxn ? s[s.length()-j-1] == s[s.length()-i-2] : s[j] == s[i+1] ) ) i++;
    }
    if ( i ) i += 2;
    return i;
}

int getHomopolymerScore( string &s )
{
    int score = 0;
    int i = 0, j = 0;
    while ( j++ < s.length() )
    {
        if ( s[i] != s[j] )
        {
            if ( j > i + 2 )
            {
                score += j - i - 2;
            }
            i = j;
        }
    }
    
    score += j - i - 1;
    return score;
}

bool mapSeqCongruent( string &left, string &right, int offset )
{
    int i = 0;
    while ( i + offset < left.length() && i < right.length() && left[i+offset] == right[i++] );
    return i + offset == left.length();
}

bool mapSeqEnd( string &q, string &t, int minLen, int32_t* coords, bool drxn )
{
    if ( minLen > q.length() ) return false;
    string qSub = drxn ? q.substr( q.length() - minLen ) : q.substr( 0, minLen );
    size_t it = t.find( qSub );
    bool didMap = false;
    while ( it != t.npos )
    {
        didMap = true;
        coords[0] = it;
        coords[1] = it + qSub.length();
        int i = drxn ? q.length() - qSub.length(): qSub.length();
        if ( drxn )
        {
            while ( coords[0]-1 >= 0 && i-1 >= 0 && q[i-1] == t[coords[0]-1] )
            {
                coords[0]--;
                i--;
            }
        }
        else
        {
            while ( i < q.length() && coords[1] < t.length() && q[i] == t[coords[1]] )
            {
                coords[1]++;
                i++;
            }
        }
        minLen = coords[1] - coords[0] + 1;
        if ( minLen > q.length() ) break;
        qSub = drxn ? q.substr( q.length() - minLen ) : q.substr( 0, minLen );
        it = t.find( qSub );
    }
    
    return didMap;
}

int mapSeqOverlap( string &left, string &right, int minLen )
{
    if ( minLen < right.length() )
    {
        string q = right.substr( 0, minLen );
        size_t it = left.find( q );
        while ( it != left.npos )
        {
            int i = it + minLen;
            int j = minLen;
            while ( i < left.length() && j < right.length() && left[i] == right[j] )
            {
                i++;
                j++;
            }
            if ( i == left.length() ) return j;
            it = left.find( q, it + 1 );
        }
    }
    return 0;
}
