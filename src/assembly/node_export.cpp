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
#include <algorithm>

void Node::exportNodeAlign( ofstream &fh, int32_t &coords )
{
    fh << ">Node" << id_ << endl;
    int32_t offset = max( 0, ends_[0] - coords );
    fh << string( offset, '-' ) << seq_ << endl;
}

void Node::exportNodeDump( ofstream &fh )
{
    fh << "NODE=" << id_ << endl;
    
    fh << "SUBGRAPH=" << drxn_ << endl;
    
    fh << "SEQ=" << seq_ << endl;
    
    fh << "COORDS=" << ends_[0] << "_" << ends_[1] << endl;
    
    fh << "CLONE=";
    bool first = true;
    for ( Node* clone : this->getCloneSet() )
    {
        if ( !first )
        {
            fh << "_";
        }
        fh << clone->id_;
        first = false;
    }
    fh << endl;
    
    fh << "EDGES=";
    first = true;
    for ( int drxn : { 0, 1 } )
    {
        for ( Edge &e : edges_[drxn] )
        {
            if ( !first )
            {
                fh << "_";
            }
            fh << drxn << "_" << e.node->id_ << "_" << e.ol << "_" << e.isLeap;
            first = false;
        }
    }
    fh << endl;
    
    fh << "READS=";
    first = true;
    for ( auto &read : reads_ )
    {
        if ( !first )
        {
            fh << "_";
        }
        fh << read.first << "_" << read.second[0] << "_" << read.second[1];
        first = false;
    }
    fh << endl;
}

NodeList Node::importNodes( ifstream &fh )
{
    string line, tag;
    vector<string> data;
    NodeList nodes;
    unordered_map<int, Node*> nodeMap;
    vector< unordered_set<int> > cloneSets;
    vector<int> edges[4];
    Node* node = NULL;
    while ( getline( fh, line ) )
    {
        tag = line.substr( 0, line.find( '=' ) );
        data = Node::importNodesGetSubStrings( line.substr( tag.length() + 1 ) );
        
        if ( tag == "NODE" )
        {
            node = new Node();
            node->id_ = data[0];
            nodeMap[stoi( data[0] )] = node;
            nodes.push_back( node );
        }
        else if ( tag == "SUBGRAPH" )
        {
            node->drxn_ = stoi( data[0] );
        }
        else if ( tag == "SEQ" )
        {
            node->seq_ = data[0];
        }
        else if ( tag == "COORDS" )
        {
            node->ends_[0] = stoi( data[0] );
            node->ends_[1] = stoi( data[1] );
        }
        else if ( tag == "CLONE" )
        {
            if ( !data.empty() )
            {
                unordered_set<int> ids = { stoi( node->id_ ) };
                for ( string id : data )
                {
                    ids.insert( stoi( id ) );
                }

                bool doAdd = true;
                for ( unordered_set<int> &cloneSet : cloneSets )
                {
                    doAdd = doAdd && cloneSet.find( stoi( node->id_ ) ) == cloneSet.end();
                }
                
                if ( doAdd )
                {
                    cloneSets.push_back( ids );
                }
            }
        }
        else if ( tag == "EDGES" )
        {
            assert( data.size() % 4 == 0 );
            for ( int i ( 0 ); i < data.size(); i += 4 )
            {
                if ( data[i] == "1" )
                {
                    edges[0].push_back( stoi( node->id_ ) );
                    edges[1].push_back( stoi( data[i+1] ) );
                    edges[2].push_back( stoi( data[i+2] ) );
                    edges[3].push_back( stoi( data[i+3] ) );
                }
            }
        }
        else if ( tag == "READS" )
        {
            assert( data.size() % 3 == 0 );
            for ( int i ( 0 ); i < data.size(); i += 3 )
            {
                node->addRead( stoull( data[i] ), stol( data[i+1] ), stol( data[i+2] ), false );
            }
        }
        else
        {
            assert( false );
        }
    }
    
    // Add edges
    for ( int i ( 0 ); i < edges[0].size(); i++ )
    {
        nodeMap[edges[0][i]]->addEdge( nodeMap[edges[1][i]], edges[2][i], 1, false, edges[3][i] );
    }
    
    // Add clones
    for ( unordered_set<int> &cloneSet : cloneSets )
    {
        int i = *cloneSet.begin();
        for ( int j : cloneSet )
        {
            if ( i != j )
            {
                nodeMap[i]->addClone( nodeMap[j] );
            }
        }
    }
    
//    ifstream afh( "/home/glen/PythonProjects/BioJunkyard/data/Export/align7.fa"  );
//    NodeList alignNodes = Node::importAlign( afh, 1 );
//    nodes.insert( nodes.end(), alignNodes.begin(), alignNodes.end() );
//    nodeMap[682]->dismantleNode();
//    nodeMap[674]->dismantleNode();
//    NodeSet delSet;
//    delSet.insert( nodeMap[682] );
//    delSet.insert( nodeMap[674] );
//    nodeMap[654]->addEdge( alignNodes[0], 94, 1, true );
//    alignNodes[0]->addEdge( alignNodes[1], 99, 1, true );
//    alignNodes[0]->addEdge( alignNodes[2], 76, 1, true );
//    alignNodes[1]->addEdge( alignNodes[3], 78, 1, true );
//    alignNodes[2]->addEdge( alignNodes[3], 37, 1, true );
//    alignNodes[3]->addEdge( alignNodes[4], 96, 1, true );
//    alignNodes[3]->addEdge( nodeMap[688], 80, 1, true );
//    alignNodes[4]->addEdge( nodeMap[676], 90, 1, true );
//    
//    for ( int i( 0 ); i < nodes.size(); )
//    {
//        if ( delSet.find( nodes[i] ) != delSet.end() )
//        {
//            delete nodes[i];
//            nodes.erase( nodes.begin() + i );
//            continue;
//        }
//        i++;
//    }
    
    for ( Node* n : nodes )
    {
        n->setCoverage();
        if ( n->drxn_ == 2 )
        {
            n->validLimits_[0] = n->validLimits_[1] = n->ends_[0];
            n->validLimits_[2] = n->validLimits_[3] = n->ends_[1];
        }
        
        n->stop_[0] = 0;
        n->stop_[1] = 0;
    }
    
//    node = nodeMap[12];
//    
//    node->trimSeq( -1476, 0, false );
////    vector< pair<SeqNum, Coords> > reads;
////    for ( auto &read : node->reads_ )
////    {
////        reads.push_back( read );
////    }
////    
////    bool drxn = false;
////    sort( reads.begin(), reads.end(), [&drxn]( pair<SeqNum, Coords> &a, pair<SeqNum, Coords> &b ) {
////        return ( drxn ? a.second[1] > b.second[1] 
////                      : a.second[0] < b.second[0] );
////    });
//    node->addEdge( nodeMap[28], 29, 0 );
    
//    assert( false );
    
    return nodes;
}

NodeList Node::importAlign( ifstream &fh, bool drxn )
{
    NodeList nodes;
    string line;
    Node* node;
    
    while ( getline( fh, line ) )
    {
        if ( line[0] = '>' )
        {
            line = line.substr( 1 );
            if ( line[0] == 'N' )
            {
                getline( fh, line );
                
                line = line.substr( 0, line.find_last_not_of( '-' ) + 1 );
                size_t it = line.find_first_not_of( '-' );
                line = line.substr( it );
                
                node = new Node();
                node->seq_ = line;
                node->ends_[0] = it;
                node->ends_[1] = it + node->seq_.length();
                node->drxn_ = drxn;
                nodes.push_back( node );
            }
            else
            {
                SeqNum readId = stoull( line );
                getline( fh, line );
                line = line.substr( 0, line.find_last_not_of( '-' ) + 1 );
                size_t it = line.find_first_not_of( '-' );
                assert( node );
                node->addRead( readId, it, line.length(), false );
            }
        }
    }
    
    for ( Node* node : nodes )
    {
        for ( auto &read : node->reads_ )
        {
            assert( node->ends_[0] <= read.second[0] && read.second[1] <= node->ends_[1] );
        }
    }
    
    return nodes;
}

void Node::importNodes( NodeList &nodes, unordered_map<int, Node*> &nodeMap, ifstream &fh, bool drxn )
{
    NodeList thisNodes;
    Node* node = NULL;
    string line;
    vector< vector<int> > edges;
    int i = 0;
    
    while ( getline( fh, line ) )
    {
        if ( line[0] = '>' )
        {
            line = line.substr( 1 );
            if ( line[0] == 'N' )
            {
                line = line.substr( 6 );
                auto it = line.find( ' ' );
                while ( it != line.npos )
                {
                    line = line.substr( it + 1 );
                    it = line.find( ' ' );
                }
                vector<int> thisEdges;
                thisEdges.push_back( stol( line ) );
                edges.push_back( thisEdges );
                getline( fh, line );
                line = line.substr( 0, line.find_last_not_of( '-' ) + 1 );
                node = new Node();
                node->seq_ = line;
                node->ends_[0] = 0;
                node->ends_[1] = node->seq_.length();
                node->drxn_ = drxn;
                thisNodes.push_back( node );
                i++;
            }
            else
            {
                SeqNum readId = stoull( line );
                getline( fh, line );
                line = line.substr( 0, line.find_last_not_of( '-' ) + 1 );
                assert( node );
                int i = 0;
                while ( line[i] == '-' )
                {
                    i++;
                }
                node->addRead( readId, i, line.length(), false );
            }
        }
    }
    
    for ( int i( 0 ); i < thisNodes.size(); i++ )
    {
        for ( int j : edges[i] )
        {
            Node* a = ( drxn ? nodeMap[j] : thisNodes[i] );
            Node* b = ( drxn ? thisNodes[i] : nodeMap[j] );
            int x = 0;
            int edge = 0;
            while ( edge == 0 && x < a->seq_.length() )
            {
                while ( a->seq_[x] != b->seq_[0] )
                {
                    x++;
                }
                int y = 1;
                while ( a->seq_[x + y] == b->seq_[y] )
                {
                    y++;
                }
                if ( x + y == a->seq_.length() )
                {
                    edge = y;
                }
                x++;
            }
            if ( drxn )
            {
                a->addEdge( b, edge, 1 );
            }
            else
            {
                b->addEdge( a, edge, 0 );
            }
        }
    }
    
    nodes.insert( nodes.end(), thisNodes.begin(), thisNodes.end() );
}

vector<string> Node::importNodesGetSubStrings( string line )
{
    vector<string> data;
    
    int i = 0, j = 0;
    
    while ( j < line.length() )
    {
        if ( line[j] == '_' )
        {
            data.push_back( line.substr( i, j - i ) );
            i = j + 1;
        }
        j++;
    }
    
    if ( j > i )
    {
        data.push_back( line.substr( i, j - i ) );
    }
    
    return data;
}
