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

//void Node::complete( bool drxn )
//{
//    sortMarks( readMarks_, drxn );
//    vector<ReadMark> peMarks;
//    for ( ReadMark &mark : readMarks_ )
//    {
//        if ( params.isReadPe( mark.readId ) )
//        {
//            peMarks.push_back( mark );
//        }
//    }
//    
//    float peCover = coverage_ * params.peRatio;
//    int32_t len = params.maxPeMean - params.readLen;
//    int cutoff = max( 2, int(float(peCover * len) / (float)8) );
//    if ( ends_[1] - ends_[0] > 200 )
//    {
//        for ( int i( 0 ); i < peMarks.size() - cutoff; i++ )
//        {
//            int32_t dist = abs( peMarks[i].mark - peMarks[i + cutoff].mark );
//            if ( dist <= len )
//            {
//                
//            }
//        }
//    }
//}
//
