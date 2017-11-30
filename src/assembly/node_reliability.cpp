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

#include "node.h"

bool Node::isMisassembled( bool drxn )
{
    if ( !assembled_[drxn] && !misassembled_[drxn] && !isContinue( drxn ) )
    {
        int i, j;
        assembled_[drxn] = true;
        if ( isMisassembled( i, j, drxn ) )
        {
            assembled_[drxn] = false;
            misassembled_[drxn] = true;
        }
    }
    return misassembled_[drxn];
}

bool Node::isMisassembled( Querier &bwt, NodeList &nodes, int32_t* validLimits, bool drxn )
{
    int i, j;
     if ( !assembled_[drxn] && !misassembled_[drxn] && !isContinue( 0 ) && !isContinue( 1 ) && isMisassembled( i, j, drxn ) )
    {
        misassembled_[drxn] = true;
        reliable_ = false;
//        int32_t coords[2];
//        coords[0] = drxn ? max( ends_[0], marks_[1][iBest].mark - ( 2 * params.readLen) ) : marks_[0][jBest].mark;
//        coords[1] = drxn ? marks_[1][jBest].mark : min( ends_[1], marks_[0][iBest].mark + ( 2 * params.readLen ) );
//        
//        string seq = seq_.substr( coords[0] - ends_[0], coords[1] - coords[0] );
//        bwt.mapBackForks( seq, drxn );
        int32_t splitCoord = drxn ? max( ends_[0], marks_[1][i].mark - int( 1.5 * params.readLen) )
                                  : min( ends_[1], marks_[0][i].mark + int( 1.5 * params.readLen ) );
        if ( splitCoord != ends_[!drxn] )
        {
            splitCoord = findNextRead( splitCoord, drxn );
            if ( splitCoord != ends_[drxn] )
            {
                Node* node = splitNode( splitCoord, nodes, drxn, drxn );
                node->reliable_ = false;
                node->unreliable_ = true;
                node->misassembled_[drxn] = true;
                misassembled_[drxn] = false;
                propagateValidation( validLimits, drxn );
            }
        }
    }
    
    return misassembled_[drxn];
}

bool Node::isMisassembled( int &iBest, int &jBest, bool drxn )
{
    int32_t peWindow = max( params.readLen, params.avgPeMean / 2 );
    int32_t peWindowCutoff = max( float( 3 ), float( peWindow * params.peCover ) / float( params.readLen * 6 ) );
    
    sortMarks( marks_[drxn], drxn );
    iBest = -1, jBest = -1;
    for ( int i ( 0 ); i + peWindowCutoff < marks_[drxn].size(); i++ )
    {
        if ( params.isReadPe( marks_[drxn][i].readId ) 
                && abs( marks_[drxn][i].mark - marks_[drxn][i+peWindowCutoff].mark ) < peWindow )
        {
            int j = i;
            int windowCount = 1;
            while ( j < marks_[drxn].size()
                    && abs( marks_[drxn][j+1].mark - marks_[drxn][i].mark ) < peWindow )
            {
                windowCount += params.isReadPe( marks_[drxn][j+1].readId );
                j++;
            }
            if ( windowCount >= peWindowCutoff )
            {
                iBest = i;
                jBest = j-1;
                return true;
            }
        }
    }
    
    return false;
}
