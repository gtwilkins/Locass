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

#include "locus_port.h"
#include <fstream>
#include <string.h>

using namespace std;

LocusExport::LocusExport( NodeRoll& nodes, string fn )
{
    ofstream ofs( fn );
    int comma;
    for ( int i = 0; i < nodes.size(); i++ ) nodes[i]->id2_ = i;
    for ( int i = 0; i < nodes.size(); i++ )
    {
        ofs << "NODE=" << i << endl;
        
        ofs << "DRXN=" << nodes[i]->drxn_ << endl;
        
        ofs << "SEQ=" << nodes[i]->seq_ << endl;
        
        if ( nodes[i]->cloned_ )
        {
            int cloned = i;
            for ( Node* node : nodes[i]->cloned_->nodes ) cloned = min( cloned, node->id2_ );
            ofs << "CLONED=" << cloned << endl;
        }
        
        ofs << "STOP=" << nodes[i]->stop_[0] << "_" << nodes[i]->stop_[1] << endl;
        
        ofs << "ENDS=" << nodes[i]->ends_[0] << "_" << nodes[i]->ends_[1] << "_" << nodes[i]->ends_.origin[0] << "_" << nodes[i]->ends_.origin[1] << endl;
        
        ofs << "TAGS=";
        comma = 0;
        if ( nodes[i]->bad_ ) ofs << string( comma++ ? "," : "" ) + "BAD";
        if ( nodes[i]->verified_ ) ofs << string( comma++ ? "," : "" ) + "VERIFIED";
        if ( nodes[i]->mapped_ ) ofs << string( comma++ ? "," : "" ) + "MAPPED";
        if ( nodes[i]->culled_ ) ofs << string( comma++ ? "," : "" ) + "CULLED";
        ofs << endl;
        
        ofs << "READS=";
        comma = 0;
        for ( auto& read : nodes[i]->reads_ )
        {
            ofs << string( comma ? "," : "" ) << read.first << "_" << read.second[0] << "_" << read.second[1] << "_" <<  read.second.coords[2] << "_" << ( read.second.redundant ? "1" : "0" ) << "_" << ( read.second.ignore ? "1" : "0" );
            comma = 1;
        }
        ofs << endl;
    }
    
    ofs << "EDGES=";
    comma = 0;
    for ( int i = 0; i < nodes.size(); i++ ) for ( int d : { 0, 1 } ) for ( Edge& e : nodes[i]->edges_[d] )
    {
        ofs << string( comma ? "," : "" ) << i << "_" << e.node->id2_ << "_" << d << "_" << e.ol << "_" << e.leap;
        comma = 1;
    }
    ofs << endl;
    
//    ofs << "EDGES=";
//    comma = 0;
//    for ( int i = 0; i < nodes.size(); i++ ) for ( Edge& e : nodes[i]->edges_[1] )
//    {
//        ofs << string( comma ? "," : "" ) << i << "_" << e.node->id2_ << "_" << e.ol << "_" << e.isLeap;
//        comma = 1;
//    }
//    ofs << endl;
    
    ofs.close();
}

LocusImport::LocusImport( NodeRoll& nodes, string fn )
{
    assert( nodes.empty() );
    ifstream ifs( fn );
    string tag, line;
    Node* node = NULL;
    bool edged = false;
    
    unordered_map<int, Node*> ids;
    unordered_map<int, vector<Node*> > clones;
    
    while ( getline( ifs, line ) && !line.empty() )
    {
        size_t it = line.find( '=' );
        assert( it != string::npos );
        assert( !edged );
        string tag = line.substr( 0, it );
        if ( tag != "NODE" ) assert( node );
        
        if ( tag == "NODE" )
        {
            node = new Node();
            nodes += node;
            vector<int64_t> data = parse( parse( line, it ) );
            assert( data.size() == 1 );
            node->id2_ = data[0];
            assert( ids.insert( make_pair( node->id2_, node ) ).second );
        }
        else if ( tag == "DRXN" )
        {
            vector<int64_t> data = parse( parse( line, it ) );
            assert( data.size() == 1 );
            node->drxn_ = data[0];
        }
        else if ( tag == "SEQ" )
        {
            node->seq_ = parse( line, it );
            assert( node->seq_.find_first_not_of( "ACTGN" ) == string::npos && !node->seq_.empty() );
        }
        else if ( tag == "CLONED" )
        {
            auto ins = clones.insert( make_pair( stoi( parse( line, it ) ), vector<Node*>{ node } ) );
            if ( !ins.second ) ins.first->second.push_back( node );
        }
        else if ( tag == "STOP" )
        {
            vector<int64_t> data = parse( parse( line, it ) );
            assert( data.size() == 2 );
            node->stop_[0] = data[0];
            node->stop_[1] = data[1];
        }
        else if ( tag == "ENDS" )
        {
            vector<int64_t> data = parse( parse( line, it ) );
            assert( data.size() == 4 );
            node->ends_[0] = data[0];
            node->ends_[1] = data[1];
            node->ends_.origin[0] = data[2];
            node->ends_.origin[1] = data[3];
        }
        else if ( tag == "TAGS" )
        {
            for ( string mod = parse( line, it ); !mod.empty(); )
            {
                if ( mod == "BAD" ) node->bad_ = true;
                else if ( mod == "VERIFIED" ) node->verified_ = true;
                else if ( mod == "MAPPED" ) node->mapped_ = true;
                else if ( mod == "CULLED" ) node->culled_ = true;
                else assert( false );
                mod = parse( line, it );
            }
        }
        else if ( tag == "READS" )
        {
            for ( vector<int64_t> data = parse( parse( line, it ) ); !data.empty(); )
            {
                assert( data.size() == 6 );
                Coords coords( data[1], data[2] );
                coords.coords[2] = data[3];
                coords.redundant = data[4];
                coords.ignore = data[5];
                node->reads_.insert( make_pair( data[0], coords ) );
                data = parse( parse( line, it ) );
            }
        }
        else if ( tag == "EDGES" )
        {
            for ( vector<int64_t> data = parse( parse( line, it ) ); !data.empty(); )
            {
                assert( data.size() == 5 );
                auto a = ids.find( data[0] );
                auto b = ids.find( data[1] );
                assert( a != ids.end() && b != ids.end() );
                a->second->addEdge( Edge( b->second, data[3], data[4] ), data[2], false );
                data = parse( parse( line, it ) );
            }
//            for ( vector<int64_t> data = parse( parse( line, it ) ); !data.empty(); )
//            {
//                assert( data.size() == 4 );
//                auto l = ids.find( data[0] );
//                auto r = ids.find( data[1] );
//                assert( l != ids.end() && r != ids.end() );
//                l->second->addEdge( Edge( r->second, data[2], data[3] ), 1, true );
//                data = parse( parse( line, it ) );
//            }
            edged = true;
        }
        else assert( false );
    }
    
    for ( const pair<int, vector<Node*> >& cn : clones ) for ( int i = 1; i < cn.second.size(); i++ ) cn.second[0]->addCloned( cn.second[i] );
    for ( Node* n : nodes.nodes ) n->setCoverage();
}

string LocusImport::parse( string s, size_t& it )
{
    if ( it == string::npos ) return "";
    
    size_t it2 = s.find( ',', ++it );
    string r = s.substr( it, it2 == string::npos ? string::npos : it2 - it );
    it = it2;
    return r;
}

vector<int64_t> LocusImport::parse( string s )
{
    vector<int64_t> r;
    if ( s.empty() ) return r;
    size_t it[2]{ 0, 0 };
    while ( it[0] != string::npos && it[0] < s.size() )
    {
        it[1] = s.find( '_', it[0] );
        string tmp = s.substr( it[0], it[1] == string::npos ? string::npos : it[1] - it[0] );
        assert( !tmp.empty() );
        r.push_back( stoll( tmp ) );
        it[0] = it[1] == string::npos ? string::npos : it[1]+1;
    }
    assert( !r.empty() );
    return r;
}
