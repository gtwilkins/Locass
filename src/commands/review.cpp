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

#include "review.h"
#include "parameters.h"
#include "filenames.h"
#include "error.h"
#include "query.h"
#include "query_binary.h"
#include "index_reader.h"
#include "reviewer.h"
#include <iostream>
#include <string.h>
#include <fstream>
#include <cassert>

extern Parameters params;

Review::Review( int argc, char** argv )
{
    Filenames* fns = NULL;
    string inFile[2], outFile, line;
    vector<string> seqs[2];
    bool diploid = false;
    for ( int i ( 2 ); i < argc; i++ )
    {
        if ( !strcmp( argv[i], "-h" ) ) printUsage( false );
        else if ( !strcmp( argv[i], "-i" ) )
        {
            if ( !inFile[0].empty() ) failure( "Multiple inputs provided." );
            inFile[0] = argv[++i];
        }
        else if ( !strcmp( argv[i], "-o" ) )
        {
            if ( !outFile.empty() ) failure( "Multiple outputs provided." );
            outFile = argv[++i];
        }
        else if ( !strcmp( argv[i], "-p" ) )
        {
            if ( fns ) failure( "More than one data prefix provided." );
            fns = new Filenames( argv[++i] );
        }
        else if ( !strcmp( argv[i], "-m" ) )
        {
            if ( !inFile[1].empty() ) failure( "Multiple mapping files provided." );
            inFile[1] = argv[++i];
        }
        else if ( !strcmp( argv[i], "--diploid" ) ) diploid = true;
        else printUsage( "Unrecognised argument: \"" + string( argv[i] ) + "\"" );
    }
    
    if ( argc <= 2 ) printUsage( "No arguements were supplied, see usage:" );
    if ( !fns ) failure( "No data prefix supplied." );
    if ( inFile[0].empty() ) failure( "No input filename given." );
    if ( outFile.empty() ) outFile = "locass-review";
    
    IndexReader* ir = new IndexReader( fns );
    QueryBinaries* qb = new QueryBinaries( fns );
    Querier bwt( ir, qb );
    params.checkReady();
    
    inFile[1] = "/home/glen/Thesis/Chapter 4/HtExons.fasta";
    for ( int i : { 0, 1 } ) if ( !inFile[i].empty() )
    {
        ifstream ifs( inFile[i] );
        while ( getline( ifs, line ) && !line.empty() && line[0] == '>' && getline( ifs, line ) ) seqs[i].push_back( line );
        ifs.close();
    }
    if ( diploid && seqs[0].size() != 2 ) failure( "Allelic flag expects 2 sequences as input." );
    Reviewer::review( bwt, seqs[0], seqs[1], outFile, diploid );
    assert( false );
}

void Review::printUsage( string msg )
{
    cerr << "Error: " << msg << endl << endl;
    printUsage( true );
}

void Review::printUsage( bool failed )
{
    cout << endl << "Locass version " << LOCASS_VERSION << endl;
    cout << endl << "Command:" << endl;
    cout << "    review" << endl;
    cout << endl << "Synopsis:" << endl;
    cout << "    Generates complete data on read coverage and pairing throughout input sequences." << endl;
    cout << endl << "Usage:" << endl;
    cout << "    locass assemble [args]" << endl;
    cout << endl << "Required arguments:" << endl;
    cout << "    -p    Prefix for BWT data files." << endl;
    cout << "    -i    Fasta file containing one or more sequences to review." << endl;
    cout << endl << "Optional arguments:" << endl;
    cout << "    -o    Output preview (default: ./locass-review" << endl;
//    cout << endl << "Notes:" << endl;
//    cout << "    - The data file prefix should be the same as was supplied to the preprocess command." << endl;
//    cout << "    - Accepted seed file formats are fasta, or a list of sequences, one per line." << endl;
    
    exit( failed ? EXIT_FAILURE : EXIT_SUCCESS );
}

