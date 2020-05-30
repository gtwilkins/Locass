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

#include "export_file.h"
#include <algorithm>

ExportFile::ExportFile( string filePrefix )
{
    filePrefix_ = filePrefix;
    overWrite_ = true;
    fileCount_ = 0;
}

void ExportFile::createFile()
{
    setFileName();
    outFile_.open( filename_ );
}

void ExportFile::closeFile()
{
    outFile_.close();
}

ExportFasta::ExportFasta(string filePrefix)
        : ExportFile( filePrefix )
{}

void ExportFasta::exportContig( string header, string &seq )
{
    outFile_ << ">" << header << endl;
    outFile_ << seq << endl;
}

void ExportFasta::exportExtends( string header, string &origin, string &left, string &right )
{
    outFile_ << ">" << header << endl;
    outFile_ << origin << endl;
    outFile_ << ">" << "left_extension_0" << endl;
    outFile_ << left << endl;
    outFile_ << ">" << "right_extension_0" << endl;
    outFile_ << right << endl;
}

void ExportFasta::exportLoci( vector<Locus*> loci )
{
    createFile();
    
    for ( Locus* locus : loci )
    {
        NodeList nodes = locus->getAllNodes();
        sort( nodes.begin(), nodes.end(), [](Node* a, Node* b){ return a->ends_[0] < b->ends_[0]; } );
        int offset = nodes[0]->ends_[0];
        for ( Node* node : nodes )
        {
            outFile_ << node->getHeader( locus->header_ ) << endl;
            outFile_ << string( node->ends_[0] - offset, '-' ) << node->seq_ << endl;
        }
    }
    
    closeFile();
}

void ExportFasta::addNode( string header, string seq )
{
    outFile_ << ">" << header << "\n";
    outFile_ << seq << "\n";
}

void ExportFasta::setFileName()
{
    while ( true )
    {
        filename_ = filePrefix_ + ( fileCount_ > 0 ? "-" + to_string( fileCount_ ) : "" ) + ".fa";
        fileCount_++;
        ifstream checkFile( filename_ );
        if ( !checkFile.good() )
        {
            break;
        }
    }
}

ExportState::ExportState(string filePrefix)
        : ExportFile( filePrefix )
{}

void ExportState::setFileName()
{
    while ( true )
    {
        filename_ = filePrefix_ + ( fileCount_ > 0 ? "-" + to_string( fileCount_ ) : "" ) + ".dump";
        fileCount_++;
        ifstream checkFile( filename_ );
        if ( !checkFile.good() )
        {
            break;
        }
    }
}

void ExportState::exportLoci( vector<Locus*> loci )
{
    createFile();
    
    for ( int i( 0 ); i < loci.size(); i++ )
    {
        exportLocus( loci[i] );
    }
    closeFile();
}

void ExportState::exportLocus(Locus* locus)
{
    vector< vector<Node*> > nodes = locus->exportNodes();
    outFile_ << "LOCUS=" << locus->header_ << "\n";
    for ( int drxn : {2, 0, 1 } )
    {
        for ( Node* node : nodes[drxn] )
        {
            outFile_ << "\tNODE=" << node->id_ << "," << drxn <<"\n";
            outFile_ << "\t\tSEQ=" << node->seq_ << "\n";
            outFile_ << "\t\tGUIDED=" << 0 << "\n";
            outFile_ << "\t\tCOORDS=" << node->ends_[0] << "_" << node->ends_[1] << "\n";
            outFile_ << "\t\tSTOP=" << node->stop_[0] << "_" << node->stop_[1] << "\n";
            outFile_ << "\t\tL_EDGES=";
            for ( int i = 0; i < node->edges_[0].size(); i++ )
            {
                if ( i > 0 )outFile_ << ",";
                outFile_ << node->edges_[0][i].node->id_ << "_" << node->edges_[0][i].ol;
            }
            outFile_ << "\n";

            outFile_ << "\t\tR_EDGES=";
            for ( int i = 0; i < node->edges_[1].size(); i++ )
            {
                if ( i > 0 ) outFile_ << ",";
                outFile_ << node->edges_[1][i].node->id_ << "_" << node->edges_[1][i].ol;
            }
            outFile_ << "\n";

            outFile_ << "\t\tREADS=";
            int i = 0;
            for ( auto read : node->reads_ )
            {
                if ( i > 0 ) outFile_ << ",";
                outFile_ << read.first << "_" << read.second[0] << "_" << read.second[1] << "_" << (int)read.second.redundant;
                i++;
            }
            outFile_ << "\n";
        }
    }
}

void ExportState::setInFileName( string filename )
{
    filename_ = filename;
}

void ExportState::updateLocus( Locus* locus, bool del )
{
    ifstream inFile( filename_ );
    createFile();
    string line;
    bool doWrite = true;
    if ( inFile.is_open() )
    {
        while ( getline( inFile, line ) )
        {
            if ( line.substr( 0, 5 ) == "LOCUS" )
            {
                if ( line.substr( 6 ) == locus->header_ )
                {
                    doWrite = false;
                    if ( !del )
                    {
                        exportLocus( locus );
                    }
                }
                else
                {
                    doWrite = true;
                }
            }
            if ( doWrite )
            {
                outFile_ << line << "\n";
            }
        }
    }
    closeFile();
}

