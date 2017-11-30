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
            i += 2;
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
    Querier bwt( ir, qb );
    cout << "Calibrating parameters for dataset with:" << endl << endl;
    cout << "\t" << to_string( params.seqCount / 2 ) << " reads" << endl;
    cout << "\t" << to_string( params.libs.size() ) << " paired libraries" << endl << endl;
    CalibrateWriter cal ( bwt );
    
    cout << "Calibrating read coverage... " << flush;
    cal.coverage();
    cout << "complete!" << endl;
    
    cout << "Calibrating library specs... " << flush;
    cal.pairing();
    cout << "complete!" << endl;
    
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

