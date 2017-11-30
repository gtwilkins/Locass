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

#include "import_file.h"

ImportFile::ImportFile( string filename )
{
    filename_ = filename;
}

vector<string> ImportFile::splitString(string s, char delim)
{
    vector<string> splitString;
    if ( !s.empty() )
    {
        int pos1 = 0;
        int pos2 = 0;
        do
        {
            pos2 = s.find( delim, pos1 );
            splitString.push_back( s.substr( pos1, pos2 - pos1 ) );
            pos1 = pos2 + 1;
        } while ( pos2 != s.npos );
    }
    return splitString;
}

ImportAlign::ImportAlign(string filename ) 
        : ImportFile( filename ) {}

vector<string> ImportAlign::getSeqs()
{
    ifstream inFile( filename_ );
    string line;
    vector<string> seqs;
    if ( inFile.is_open() )
    {
        while ( getline( inFile, line ) )
        {
            seqs.push_back( line );
        }
    }
    return seqs;
}

ImportState::ImportState( string filename )
        : ImportFile( filename ) {}

vector<Locus*> ImportState::getLoci()
{
    vector<Locus*> loci;
//    Locus* locus = NULL;
//    ifstream inFile( filename_ );
//    string line, tag;
//    NodeMap nodes;
//    Edges edges;
//    Node* node = NULL;
//    if ( inFile.is_open() )
//    {
//        vector<string> data;
//        vector<string> parts;
//        while ( getline( inFile, line ) )
//        {
//            if ( line.find_first_not_of(" \t") == string::npos )
//            {
//                cerr << "Error reading file." << endl;
//                exit(1);
//            }
//            line = line.substr( line.find_first_not_of(" \t") );
//            tag = line.substr( 0, line.find( '=' ) );
//            data = splitString( line.substr( tag.length() + 1 ), ',' );
//            if ( !locus && tag != "LOCUS")
//            {
//                cerr << "Error reading file." << endl;
//                exit(1);
//            }
//            if ( tag == "LOCUS")
//            {
//                node = NULL;
//                if ( locus )
//                {
//                    addEdgesToNodes( nodes, edges );
//                    loci.push_back( locus );
//                }
//                locus = new Locus();
//                locus->header_ = data[0];
//            }
//            else if ( !node && tag != "NODE" )
//            {
//                cerr << "Error reading file." << endl;
//                exit(1);
//            }
//            else if ( tag == "NODE" )
//            {
//                assert( data.size() == 2 );
//                node = new Node();
//                node->id_ = data[0];
//                nodes[node->id_] = node;
//                locus->addNode( node, stoi( data[1] ) );
//                
//            }
//            else if ( tag == "SEQ" )
//            {
//                node->seq_ = data[0];
//            }
//            else if ( tag == "GUIDED" )
//            {
////                node->guided_ = (bool)stoi( data[0] );
//            }
//            else if ( tag == "STOP" )
//            {
//                parts = this->splitString( data[0], '_' );
//                node->stop_[0] = stoi( parts[0] );
//                node->stop_[1] = stoi( parts[1] );
//            }
//            else if ( tag == "COORDS" )
//            {
//                parts = this->splitString( data[0], '_' );
//                node->ends_[0] = stoi( parts[0] );
//                node->ends_[1] = stoi( parts[1] );
//            }
//            else if ( tag == "L_EDGES" )
//            {
//                for ( string &datum : data )
//                {
//                    parts = this->splitString( datum, '_' );
//                    InEdge edge( parts[0], node->id_, stoi( parts[1] ) );
//                    if ( find( edges.begin(), edges.end(), edge ) == edges.end() )
//                    {
//                        edges.push_back( edge );
//                    }
//                }
//            }
//            else if ( tag == "R_EDGES" )
//            {
//                for ( string &datum : data )
//                {
//                    parts = this->splitString( datum, '_' );
//                    InEdge edge( node->id_, parts[0], stoi( parts[1] ) );
//                    if ( find( edges.begin(), edges.end(), edge ) == edges.end() )
//                    {
//                        edges.push_back( edge );
//                    }
//                }
//            }
//            else if ( tag == "READS" )
//            {
//                for ( string &datum : data )
//                {
//                    parts = this->splitString( datum, '_' );
//                    Coords coords( stoi( parts[1] ), stoi( parts[2] ), (bool)stoi( parts[3] ) );
//                    SeqNum readNum = stoul( parts[0] );
//                    node->reads_.insert( make_pair( readNum, coords ) );
//                }
//            }
//            else
//            {
//                assert( false );
//            }
//        }
//        addEdgesToNodes( nodes, edges );
//        loci.push_back( locus );
//    }
//    for ( Locus* locus : loci )
//    {
//        locus->rebootLocus();
//    }
    return loci;
}

void ImportState::addEdgesToNodes(NodeMap &nodes, Edges &edges)
{
    for ( InEdge &edge : edges )
    {
        nodes[get<0>(edge)]->edges_[1].push_back( Edge( nodes[get<1>(edge)], get<2>(edge), false ) );
        nodes[get<1>(edge)]->edges_[0].push_back( Edge( nodes[get<0>(edge)], get<2>(edge), false ) );
    }
    
    nodes.clear();
    edges.clear();
}
