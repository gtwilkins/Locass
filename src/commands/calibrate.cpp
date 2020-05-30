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

#include "calibrate.h"
#include <string.h>

Calibrate::Calibrate( int argc, char** argv )
{
    Filenames* fns = NULL;
    bool covered = false, forced = false;
    float cover;
    vector<int> libs, dists;
    vector<uint8_t> orients;
    for ( int i ( 2 ); i < argc; i++ )
    {
        if ( !strcmp( argv[i], "-h" ) )
        {
            printUsage();
            exit( EXIT_SUCCESS );
        }
        else if ( !strcmp( argv[i], "-p" ) )
        {
            if ( fns )
            {
                cerr << "Error: more than one prefix provided." << endl;
                exit( EXIT_FAILURE );
            }
            fns = new Filenames( argv[i+1] );
            i++;
        }
        else if ( isdigit( argv[i][1] ) )
        {
            if ( i+3 > argc )
            {
                cerr << "Error: not enough support information given for paired library update." << endl;
                exit( EXIT_FAILURE );
            }
            
            libs.push_back( atoi( &argv[i++][1] ) );
            orients.push_back( !strcmp( argv[i], "FR" ) ? 2 : ( !strcmp( argv[i], "RF" ) ? 3 : 0 ) );
            dists.push_back( atoi( argv[++i] ) );
            forced = true;
            
            if ( !orients.back() )
            {
                cerr << "Error: did not recognise user-defined library type. Expecting FR or RF, received: " << argv[i+1] << endl;
                exit( EXIT_FAILURE );
            }
        }
        else if ( !strcmp( argv[i], "-c" ) )
        {
            covered = forced = true;
            cover = atoi( argv[++i] );
        }
        else
        {
            cerr << "Unrecognised argument: \"" << argv[i] << "\"" << endl << endl;
            printUsage();
            exit( EXIT_FAILURE );
        }
    }
    
    if ( argc <= 2 )
    {
        cerr << "Error: no arguements supplied, see usage:" << endl << endl;
        printUsage();
        exit( EXIT_FAILURE );
    }
    
    if ( !fns )
    {
        cerr << "Error: no prefix supplied." << endl;
        exit( EXIT_FAILURE );
    }
    
    IndexReader* ir = new IndexReader( fns );
    QueryBinaries* qb = new QueryBinaries( fns );
    params.cover = 0;
    Querier bwt( ir, qb );
    
    if ( forced )
    {
        if ( covered ) forceCoverage( cover );
        for ( int i = 0; i < libs.size(); i++ ) forcePaired( libs[i], dists[i], orients[i] );
        write( fns );
        exit( EXIT_SUCCESS );
    }
    
    cout << "Calibrating parameters for dataset with:" << endl << endl;
    cout << "    " << to_string( params.seqCount / 2 ) << " reads" << endl;
    cout << "    " << to_string( params.libs.size() ) << " paired libraries" << endl << endl;
    CalibrateWriter cal ( bwt );
    
    cout << "Calibrating read coverage... " << endl;
    cal.coverage();
    cout << "Calibrating read coverage... complete!" << endl;
    
    cout << "Calibrating library specs... " << endl;
    cal.pairing();
    cout << "Calibrating library specs... complete!" << endl;
    
    cout << "Calibration complete!" << endl << endl;
    cal.write( fns );
    cout << endl << "Sequence dataset has been successfully calibrated and is ready to perform assemblies!" << endl;
}

void Calibrate::printUsage()
{
    cout << endl << "Locass version " << LOCASS_VERSION << endl;
    cout << endl << "Command:" << endl;
    cout << "\tcalibrate" << endl;
    cout << endl << "Synopsis:" << endl;
    cout << "\tEstimates sequencing coverage and pair read library specifications in preparation for assembly" << endl;
    cout << endl << "Usage:" << endl;
    cout << "\tlocass calibrate [args]" << endl;
    cout << endl << "Required arguments:" << endl;
    cout << "\t-p\tPrefix for sequence data files (as output from the transform command)." << endl;
}

void Calibrate::forceCoverage( float cover )
{
    params.cover = cover;
}

void Calibrate::forcePaired( int i, int dist, uint8_t orient )
{
    if ( params.libs.size() < i )
    {
        cerr << "Error: specified library number \"" << i << "\" is not a positive integer." << endl;
        exit( EXIT_FAILURE );
    }
    if ( params.libs.size() < i )
    {
        cerr << "Error: specified library number \"" << i << "\" is greater thatn the number of libraries in this dataset (" << params.libs.size() << " paired libraries)." << endl;
        exit( EXIT_FAILURE );
    }
    
    params.libs[i-1].size = dist;
    params.libs[i-1].orientation = orient;
    params.libs[i-1].setMinMax();
    params.libs[i-1].isPe = orient == 2;
}

void Calibrate::write( Filenames* fns )
{
    if ( !params.isCalibrated && params.cover && ( params.isCalibrated = true ) )
    {
        for ( Lib& lib : params.libs ) if ( !lib.size ) params.isCalibrated = false;
    }
    
    uint32_t coverage = params.cover * 100000;
    FILE* fp = fns->getBinary( true, true );
    fseek( fp, 11, SEEK_SET );
    fwrite( &params.isCalibrated, 1, 1, fp );
    fwrite( &coverage, 4, 1, fp );
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        uint16_t libMed = params.libs[i].size;
        uint16_t libMin = libMed ? params.libs[i].minDist : 0;
        uint16_t libMax = libMed ? params.libs[i].maxDist : 0;
        fseek( fp, 25 + ( i * 12 ), SEEK_SET );
        fwrite( &libMed, 2, 1, fp );
        fwrite( &libMin, 2, 1, fp );
        fwrite( &libMax, 2, 1, fp );
        fwrite( &params.libs[i].orientation, 1, 1, fp );
        fwrite( &params.libs[i].isPe, 2, 1, fp );
    }
    
    fclose( fp );
}
    