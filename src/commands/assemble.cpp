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

#include "assemble.h"
#include "parameters.h"
#include "timer.h"
#include "seed.h"
#include "query.h"
#include "query_binary.h"
#include "index_reader.h"
#include "extend.h"
#include <iostream>
#include <string.h>

extern Parameters params;

Assemble::Assemble( int argc, char** argv )
{
    Filenames* fns = NULL;
    string inFile, outFile;
    int32_t limit = 10000;
    int errorCount = 5;
    bool drxns[2] = { false, false };
    
    for ( int i ( 2 ); i < argc; )
    {
        if ( !strcmp( argv[i], "-h" ) )
        {
            printUsage();
            exit( EXIT_SUCCESS );
        }
        else if ( !strcmp( argv[i], "-i" ) )
        {
            if ( !inFile.empty() )
            {
                cerr << "Error: multiple inputs provided." << endl;
                exit( EXIT_FAILURE );
            }
            inFile = argv[i+1];
            i += 2;
        }
        else if ( !strcmp( argv[i], "-o" ) )
        {
            if ( !outFile.empty() )
            {
                cerr << "Error: multiple inputs provided." << endl;
                exit( EXIT_FAILURE );
            }
            outFile = argv[i+1];
            i += 2;
        }
        else if ( !strcmp( argv[i], "-p" ) )
        {
            if ( fns )
            {
                cerr << "Error: more than one output prefix provided." << endl;
                exit( EXIT_FAILURE );
            }
            fns = new Filenames( argv[i+1] );
            i += 2;
        }
        else if ( !strcmp( argv[i], "-l" ) )
        {
            for ( int j ( 0 ); j < strlen( argv[i+1] ); j++ )
            {
                if ( !isdigit( argv[i+1][j] ) )
                {
                    cerr << "Error: \"" << argv[i+1] << "\" is not recognised as a number." << endl;
                    exit( EXIT_FAILURE );
                }
            }
            limit = stoi( argv[i+1] );
            i += 2;
        }
        else if ( !strcmp( argv[i], "--haploid" ) )
        {
            params.haploid = true;
            i++;
        }
        else if ( !strcmp( argv[i], "--right" ) )
        {
            drxns[1] = true;
            params.drxns[0] = false;
            i++;
        }
        else if ( !strcmp( argv[i], "--left" ) )
        {
            drxns[0] = true;
            params.drxns[0] = true;
            i++;
        }
        else
        {
            cerr << "Unrecognised argument: \"" << argv[i] << "\"" << endl << endl;
            printUsage();
            exit( EXIT_FAILURE );
        }
    }
    
    if ( drxns[0] && drxns[1] ) params.drxns[0] = params.drxns[1] = true;
    
    if ( argc <= 2 )
    {
        cerr << "Error: no arguements supplied, see usage:" << endl << endl;
        printUsage();
        exit( EXIT_FAILURE );
    }
    
    if ( !fns )
    {
        cerr << "Error: no data prefix supplied." << endl;
        exit( EXIT_FAILURE );
    }
    
    if ( inFile.empty() )
    {
        cerr << "Error: no input filename given." << endl;
        exit( EXIT_FAILURE );
    }
    
    if ( outFile.empty() )
    {
        int i = 0;
        outFile = "locass-out.fa";
    }
    
    if ( limit > 200000 )
    {
        cerr << "Error: maximum extension limit is 200,000." << endl;
        exit( EXIT_FAILURE );
    }
    
    params.setLimits( limit );
    vector<string> headers, seqs;
    ifstream ifp = fns->getReadStream( inFile );
    readIn( ifp, headers, seqs );
    ifp.close();
    
    if ( seqs.size() != headers.size() )
    {
        cerr << "Unexpected error reading seed sequences." << endl;
        exit( EXIT_FAILURE );
    }
    
    IndexReader* ir = new IndexReader( fns );
    QueryBinaries* qb = new QueryBinaries( fns );
    Querier bwt( ir, qb );
    params.checkReady();
    
    cout << "Performing seeded locus assembly" << endl << endl;
    cout << "Error rate set to " << to_string( errorCount ) << "%" << endl;
    cout << "Extension limit set to " << to_string( limit ) << endl;
    cout << "Read in " << to_string( seqs.size() ) << ( seqs.size() > 1 ? " seed sequences." : " seed sequence." ) << endl << endl;
    
    ofstream ofp = fns->getWriteStream( outFile );
    
    Extend extend( ofp );
    double startTime = clock();
    int locusCount = 0, seedCount = 0;
    for ( int i = 0; i < seqs.size(); i++ )
    {
        cout << "Assembling locus " << to_string( i + 1 ) << " of " << to_string( seqs.size() ) << endl;
        cout << "    Searching for target loci... " << endl;
        Seed seed( headers[i], seqs[i], bwt, errorCount );
        seedCount++;
        if ( seed.warning() ) continue;
        seed.assemble();
        vector<Locus*> loci = seed.getLoci();
        cout << ( loci.empty() ? "    No target loci found." : "    tTarget loci found!" ) << endl;
        if ( !loci.empty() )
        {
            cout << "    Extending seeded loci... " << endl;
            extend.extend( loci );
            cout << "\tExtension complete!" << endl;
            locusCount+= loci.size();
        }
        cout << endl;
    }
    vector<int32_t> contigLens = extend.close();
    ofp.close();
    
    int32_t contigLenTotal = 0, cumulativeLen = 0, n50 = 0;
    for ( int32_t len : contigLens )
    {
        contigLenTotal += len;
    }
    int32_t halfTotalLen = contigLenTotal / 2;
    
    for ( int32_t len : contigLens )
    {
        cumulativeLen += len;
        n50 = len;
        if ( cumulativeLen >= halfTotalLen ) break;
    }
    
    
    cout << "Assembly completed!" << endl << endl;
    cout << "Assembled " << to_string( locusCount ) << ( locusCount > 1 ? " loci from " : " locus from " ) 
            << to_string( seedCount ) << ( seedCount > 1 ? " seed sequences" : "seed sequence." ) << endl;
    
    if ( locusCount > 1 )
    {
        cout << "Average extension length: " << to_string( contigLenTotal / ( locusCount * 2 ) ) << endl;
        cout << "Assembled length N50: " << to_string( n50 ) << endl;
    }
    else
    {
        cout << "Length assembled: " << to_string( contigLenTotal ) << endl;
    }
    
    cout << "Total time taken: " << getDuration( startTime ) << endl;
    cout << endl << "Output saved to \"" << outFile << "\"" << endl;
}

void Assemble::readIn( ifstream &fh, vector<string> &headers, vector<string> &seqs )
{
    string line;
    bool isFasta = false;
    int lineNum = 0;
    int seqNum = 0;
    string acceptedChars = "acgtACGT";
    while ( getline( fh, line ) && !line.empty() )
    {
        size_t it = line.find( '\r' );
        if ( it != line.npos ) line.erase( it, line.npos );
        if ( !lineNum )
        {
            isFasta = line[0] == '>';
        }
        ++lineNum;
        ++seqNum;
        string header;
        if ( isFasta )
        {
            if ( line[0] != '>' )
            {
                cerr << "Error: expected but did not find a fasta header on line " << to_string( lineNum ) << " of input file." << endl;
                exit( EXIT_FAILURE );
            }
            header = line.substr( 1 );
            getline( fh, line );
            it = line.find( '\r' );
            if ( it != line.npos ) line.erase( it, line.npos );
            ++lineNum;
        }
        
        if ( header.empty() )
        {
            header = "Input_sequence_" + to_string( seqNum );
        }
        if ( line.empty() ) break;
        
        it = line.find_first_not_of( acceptedChars );
        if ( it != line.npos )
        {
            cerr << "Error: invalid sequence character \"" << line[it] << "\" on line " << to_string( lineNum ) << " of input file." << endl;
            exit( EXIT_FAILURE );
        }
        
        headers.push_back( header );
        seqs.push_back( line );
    }
    
    if ( headers.empty() || seqs.empty() )
    {
        cerr << "Error: did not find input sequences in input file." << endl;
        exit( EXIT_FAILURE );
    }
    
    if ( seqs.size() > 200 )
    {
        cerr << "Error: too many seed sequences provided. Locass supports a maximum of 200 seed sequences in a batch." << endl;
        exit( EXIT_FAILURE );
    }
}

void Assemble::printUsage()
{
    cout << endl << "Locass version " << LOCASS_VERSION << endl;
    cout << endl << "Command:" << endl;
    cout << "\tassemble" << endl;
    cout << endl << "Synopsis:" << endl;
    cout << "\tAssembles the nucleotide sequence of one or more target genomic loci." << endl;
    cout << endl << "Usage:" << endl;
    cout << "\tlocass assemble [args]" << endl;
    cout << endl << "Required arguments:" << endl;
    cout << "\t-p\tPrefix for read data files." << endl;
    cout << "\t-i\tSeed sequence file containing one or more sequences." << endl;
    cout << endl << "Optional arguments:" << endl;
    cout << "\t-o\tOutput file name (default: ./locass-out.fa)." << endl;
    cout << "\t-l\tExtension length limit in nucleotides (default: 10,000, maximum: 200,000)." << endl;
    cout << endl << "Notes:" << endl;
    cout << "\t- The data file prefix should be the same as was supplied to the preprocess command." << endl;
    cout << "\t- Accepted seed file formats are fasta, or a list of sequences, one per line." << endl;
}

