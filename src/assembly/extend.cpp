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

#include "extend.h"
#include "timer.h"

//Extend::Extend( string outStr )
//: outState_( outStr ), outFasta_( outStr ), outContig_( outStr+"out" ), outExtends_( outStr+"extends" ), locusCount_( 0 )
//{
//    outContig_.createFile();
//    outExtends_.createFile();
//}

vector<int32_t> Extend::close()
{
    sort( contigLens_.begin(), contigLens_.end() );
    return contigLens_;
}

void Extend::extend( vector<Locus*> loci )
{
    for ( Locus* locus : loci )
    {
        double startTime = clock();
        locus->extendLocus();
        locus->finalise();
//        ofstream align( "/media/glen/ssd/pdalign.fa" );
//        ofstream dump( "/media/glen/ssd/pddump" );
//        locus->exportLocus( align, dump );
        
        string header, origin, lft, rght;
        locus->getExtends( header, origin, lft, rght );
        int32_t len = origin.length() + lft.length() + rght.length();
        cout << "\tAssembled locus of length " << to_string( len ) << " in " << getDuration( startTime ) << endl;
        
        outFile_ << ">" << header << endl;
        outFile_ << origin << endl;
        outFile_ << ">" << header << "_left_extension" << endl;
        outFile_ << lft << endl;
        outFile_ << ">" << header << "_right_extension" << endl;
        outFile_ << rght << endl;
        
        contigLens_.push_back( len );
        delete locus;
    }
}

