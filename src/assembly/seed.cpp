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

#include "seed.h"
#include "path_seed.h"
#include "node_path.h"
#include "leap.h"
#include "locus_path.h"
#include "path_alleles.h"
#include "path_cross.h"
#include "locus_port.h"
#include "locus_gap.h"
#include "locus_fill.h"
#include "seed_node.h"
#include "review.h"
#include "locus_traverse.h"
#include <algorithm>

Seed::Seed( string &header, string &seq, Querier &bwt, int errorCount, bool bestOnly )
: header_( header ), seq_( seq ), bwt_( bwt )
{
    validLimits_[0] = seq.length();
    validLimits_[1] = 0;
    tether_[0] = seq.length();
    tether_[1] = 0;
    ends_[0] = seq.length();
    ends_[1] = 0;
    assert( errorCount <= 15 );
    MappedSeqs ms = bwt_.mapSeed( seq, errorCount, bestOnly );
    readCount_ = ms.reads.size();
    
    for ( ReadStruct &read : ms.reads ) Node::addSeed( seed_, read );
    Node::trimSeed( bwt, seed_ );
}

Seed::Seed( string fn, Querier &bwt )
: bwt_( bwt )
{
//    SeedNode::seed( bwt_, seed_, "AATTCAAACTTACGGGACTGGAGTGAAACCCCAATACGTCAACATGGAGTTGAAAGCAACTGTGATCGTTGCTATCCTGGCTGCACTTGCTATGTCAGGTAAGCTCTCAAATTTCAATTTCGATCATCTAAAAGTATTTTACACAAGACGATATTTTGATTATCGCAATTCATTGATAGGAATCATAATTATTTCTATTTTGTTTGTTAAAGGAGGACTTCATATCATTATCAATCATCACTGCTGTATCAATGGGGAG", 0, 0, 2 );
//    SeedNode::seed( bwt_, seed_, "AATTCAAACTTACGGGACTGGAGTGAAACCCCCAATACGTTAACATGGAGTTGAAAGCAACTGTGATCGTTGCTATCCTGGCTGCACTTGCTATGTCAGGTAAGCTCTCAAATTTCAATTTCGATCATCTAAAAGTATTTTACACAAGACGATATTTTGATTATCGCAATTCATTGATAGGAATCATAATTATTTCTATTTTGTTTGTTAAAGGAGGACTTCATATCATTATCAATC", 0, 0, 2 );
//    SeedNode::seed( bwt_, seed_, "AATTCAAACTTACGGGACTGGAGTGAGACCCCCAATACGTCAACATGGAGTTGAAAGCAACTGTGATCGTTGCTATCCTGGCTGCACTTGCTATGTCAGGTAAGCTCTCAAATTTCAATTTCGATCATCTAAAAGTATTTTACACAAGACGATATTTTGATTATCGCAATTCATTGATAGGAATCATAATTATTTCTATTTTGTTTGTTAAAGGAGGACTTCATATCATTATCAATCATCACTGTTGTATCGATGGGG", 0, 0, 2 );
//    SeedNode::seed( bwt_, seed_, "AATTCAAACTTACGGGACTGGAGTGAGACCCCCAATACGTCAACATGGAGTTGAAAGCAACTGTGATCGTTGCTATCCTGGCTGCACTTGCTATGTCAGGTAAGCGCTCAAATTTCAATTTCGATCATCTAAAAGTATTTTACACAAGACGATATTTTGATTATCGCAATTCATTGATAGGAATCATAATTATTTCTATTTTGTTTGTTAAAGGAGGACTTCATATCATTATCAATCATCAC", 0, 0, 2 );
//    SeedNode::seed( bwt_, seed_, "AATTCAAACTTACGGGACTGGAGTGAAACCCCCAATACGTCAACATGGAGTTGAAAGCAACTGTGATAGTTGCTATCCTGGCTGCACTTGCTATGTCAGGTAAGCTCTCAAATTTCAATTTCGATCATCTAAAAGTATTTTACACAAGACGATATTTTGATTATCGCAATTCATTGATAGGAATCATAATTATTTCTATTTTGTTTGTTAAAGGAGGACTTCATATCATTATC", 0, 0, 2 );
//    SeedNode::seed( bwt_, seed_, "AATTCAAACTTACGGGACTGGAGTGAAACCCCCAATACGTCAACATGGAGTTGAAAGCAACTGTGATCGTTGCTATCCTGGCTGCACTTGCTATGTCAGGTAAGCTCTAAAATTTCAATTTCGATCATCTAAAAGTATTTTACACAAGATGATATTTTGATTATCGCAATTCATTGATAGAAATCATAATTATTTCTATTTTGTTTGTTAAAGGAGAACTTCATATCATTAT", 0, 0, 2 );
//    SeedNode::seed( bwt_, seed_, "AATTCAAACTTACGGGACCTGAGTGAAACCCCCAATACGACAACATGGAGTTGAAAGCAACTGTGATAGTTGCTATCCTGGCTGCACTTGCTGTGTCAGGTAAGCTCTCAAATTTCAATTTCGATCATCTAAAAGTATTTTACACAAGACGATATTTTGATTATCGCAATTCATTGATAGGAATCATAATTATTTCTATTTTGTTTGTTAAAGGAGGAC", 0, 0, 2 );
//    for ( Node* node : seed_.nodes ) node->setOrigin();
//    assemble( true );
//    assert( false );
    
    LocusImport( seed_, fn );
    seed_.summarise();
    Node::reverify( seed_ );
    assemble( false );
    assert( false );
    
    
//    for ( ReadStruct &read : ms.reads ) Node::addSeed( seed_, read );
//    Node::trimSeed( bwt, seed_ );
    
    
//    LocusImport( seed_, fn );
//    seed_.summarise();
    
//    string t = "AATGTTTTATTGGAAAATTAAGTTATCACATCTTAACGCAATATGACATCTATAGATTGAATTTGGGATGGAAAGTCGAGGTCGAGCAGCACAATTATTCTTCCTGACTTATTGCCTTTGTTCTGGTGCACAACCATCCTAGTCTCCTTCTGTGACGACGAATCGTTGATGGACTGCAATCAACACCGAGTTCTTTTTAGCTCATCCGGCCCGAAGGGCCCGATGAGCTTATGCCATGGCGTGGCGTCCGTCCTTCCGTCCGTCCTTCCGTCCACAGTTTTCATAGTGCTTCTTCTTCATTATTTTTCAACCGATTTTGATTCTGATTGTTTCATATGATAGGGCTAGGAGGGGGATTTAAAACTTCTACACAGAATCTATGCTAATTTATGCAAATTTATTCTGCCACTTGCGGAAGAGATATAACCAGAATAACCAGAACAATAATTGTTAACATTGTTCATTTGTTCAAAAGAACAATTAATCTAA";
//    int32_t coords[2];
//    Coords* c;
//    cout << ">Base" << endl << string( 150, '-' ) + t << endl;
//    unordered_set<ReadId> used;
//    Nodes found;
//    for ( int d : { 0, 1 } ) for ( QueryBranch& qb : bwt_.mapBranches( t, NULL, 80, d ) ) for ( ReadId id : bwt.getIds( qb.rank, qb.count, d ) )
//    {
//        if ( !used.insert( id ).second ) continue;
//        string s = bwt_.getSequence( id );
//        assert( mapSeqEnd( s, t, 80, coords, !d ) );
//        int base = 150 + coords[0] - ( d ? 0 : s.size() - ( coords[1] - coords[0] ) );
//        cout << ">" << id << endl << string( base, '-' ) + s << endl;
//        for ( Node* node : seed_.nodes ) if ( ( c = node->getRead( id ) ) && !c->redundant ) found += node;
//    }
//    assert( false );
//    Node::seedNode( bwt_, seed_, t, 2 );
    
//    Node::reverify( seed_ );
    assemble( false );
    assert( false );
//    SeedExtend seed[2];
//    restart( seed );
    assert( false );
}

void Seed::printExts( vector<string> exts[2], ofstream* ofs )
{
    Lib* lib;
    int d, base = params.maxPeMean;
    for ( int drxn : { 0, 1 } ) for ( string& s : exts[drxn] )
    {
        print( "Query", s, drxn ? 0 : base, ofs );
        
        for ( QueryBranch& qb : bwt_.mapBranches( s, NULL, 18, drxn ) ) for ( ReadId id : bwt_.getIds( qb.rank, qb.count, drxn ) )
        {
            int ol = qb.ol;
            string seq[2]{ bwt_.getSequence( id ), "" };
            int32_t coords[2];
            assert( ol = mapSeqEnd( seq[0], s, ol, coords, !drxn ) );
            
            bool used[2]{ false, false };
            for ( string& alt : exts[drxn] ) if ( alt.find( seq[0] ) != string::npos ) used[0] = true;
            
            if (!used[0] ) print( to_string( id ), seq[0], base + ( drxn ? coords[0] : coords[1] - seq[0].size() ), ofs );
            
//            if ( ( lib = params.getLib( id ) ) && lib->getPair( id, d ) && d == drxn )
//            {
//                seq[1] = bwt_.getSequence( id );
//                for ( string& alt : exts[drxn] ) if ( alt.find( seq[1] ) != string::npos ) used[1] = true;
//                
//                if ( !used[1] ) print( to_string( id ), seq[1], base + drxn ? coords[0] + lib->size - seq[1].size() : coords[1] - lib->size, ofs );
//            }
        }
    }
}

void Seed::printSeeds( vector<string>& seqs, int errors, ofstream* ofs )
{
    vector<ReadStruct> reads;
    unordered_set<ReadId> used;
    int32_t coords[2];
    Lib* lib;
    int matched = 0, qCount = 0, base = -params.readLen*3;
    vector< pair<int, pair<ReadId, string> > > prints;
    for ( string& s : seqs )
    {
        MappedSeqs ms = bwt_.mapSeed( s, errors, false );
        unordered_set<ReadId> ids;
        vector< pair<int, pair<ReadId, string> > > maps;
        for ( ReadStruct& rs : ms.reads ) if ( ids.find( rs.id ) == ids.end() )
        {
            bool mapped = ( s.find( rs.seq ) != string::npos ) || ( rs.seq.find( s ) != string::npos );
            if ( !mapped && mapSeqEnd( rs.seq, s, 20, coords, 1 ) && !coords[0] ) mapped = true;
            if ( !mapped && mapSeqEnd( rs.seq, s, 20, coords, 0 ) && coords[1] == s.size() ) mapped = true;
            if ( mapped ) used.insert( rs.id );
            
            ReadId id = rs.id;
            int d, est = rs.coords[0], ol = 0;
            string seq = rs.seq;
            if ( ( lib = params.getLib( rs.id ) ) && lib->getPair( id, d ) )
            {
                string alt = bwt_.getSequence( id );
                if ( ol = mapSeqOverlap( seq, alt, 15, d ) )
                {
                    matched++;
                    seq = d ? seq + alt.substr( ol ) : alt + seq.substr( ol );
                    if ( !d ) est -= alt.size() + ( ol ? -ol : 50 );
                    if ( mapped ) used.insert( id );
                    ids.insert( id );
                }
            }
            
            int i = 0;
            vector< pair<int, pair<ReadId, string> > >& pile = mapped ? maps : prints;
            while ( i < pile.size() && est > pile[i].first ) i++;
            pile.insert( pile.begin() + i, make_pair( est, make_pair( rs.id, seq ) ) );
        }
        
        if ( ofs ) ( *ofs ) << ">Query " << ++qCount << endl <<  string( -base, '-' ) << s << endl;
        else cout << ">Query " << ++qCount << endl <<  string( -base, '-' ) << s << endl;
        for ( int i = 0; i < maps.size(); i++ )
        {
            if ( ofs ) ( *ofs ) << ">" << maps[i].second.first << endl << string( maps[i].first - base, '-' ) + maps[i].second.second << endl;
            else cout << ">" << maps[i].second.first << endl << string( maps[i].first - base, '-' ) + maps[i].second.second << endl;
        }
    }
    
    for ( int i = 0; i < prints.size(); i++ ) if ( !used.insert( prints[i].second.first ).second ) prints.erase( prints.begin() + i-- );
    
    int len = params.readLen*.9;
    for ( int i = 0; i < prints.size(); i++ ) 
    {
        vector< pair<int, pair<ReadId, string> > > ols{ prints[i] };
        for ( int j = 0; j < ols.size(); j++ ) for ( int k = i+1; k < prints.size(); k++ )
        {
            if ( mapSeqOverlap( ols[j].second.second, prints[k].second.second, len )
                    || ols[j].second.second.find( prints[k].second.second ) != string::npos
                    || prints[k].second.second.find( ols[j].second.second ) != string::npos )
            {
                ols.push_back( prints[k] );
                prints.erase( prints.begin() + k-- );
            }
        }
        if ( ofs ) for ( pair<int, pair<ReadId, string> > p : ols ) (*ofs) << ">" << p.second.first << endl << string( p.first - base, '-' ) + p.second.second << endl;
        else for ( pair<int, pair<ReadId, string> > p : ols ) cout << ">" << p.second.first << endl << string( p.first - base, '-' ) + p.second.second << endl;
    }
}

void Seed::print( string header, string seq, int blanks, ofstream* ofs )
{
    if ( ofs ) ( *ofs ) << ">" << header << endl << string( max( 0, blanks ), '-' ) << seq << endl;
    else ( *ofs ) << ">" << header << endl << string( max( 0, blanks ), '-' ) << seq << endl;
}

void Seed::seed()
{
    vector<string> seqs;
    vector<MappedSeqs> maps;
    for ( string& s : seqs ) maps.push_back( bwt_.mapSeed( s, 10, false ) );
    for ( MappedSeqs& ms : maps ) ms.cull( maps );
    if ( !seed_.empty() ) for ( MappedSeqs& ms : maps ) cull( ms );
    assert( false );
}

void Seed::assembleRNA( bool seeded )
{
    SeedExtend seed[2];
    if ( seeded )
    {
        for ( Node* node : seed_.nodes ) for ( int d : { 0, 1 } ) if ( node->isContinue( d ) ) seed[d].addFork( node );
        for ( int d : { 0, 1 } ) while ( seed[d].extend( bwt_, seed_, true, d ) );
        Node::prune( bwt_, seed_ );
    }
//    vector<ReadId> ids{ 18635309, 24578001 };
//    for ( ReadId id : ids ) cout << ">" << id << endl << bwt_.getSequence( id ) << endl;
    
    vector< vector<NodePath*> > loci;
    vector<NodePath*> paths, seeds;
//    for ( int i = 0; i < seed_.size(); i++ ) if ( seed_[i]->seq_.find( "AGGAAAAGAGGGTGATACCTGACACAACCAGGTACGAAGTGGTTTTAATCTAGGGTAGAAAGATATAAATTCAGGAATTCAACTGATAGTTGTCAGTCAGCAGCTAATTCTCTCTTGGAAGGGGACGAAGAGAAACGTTAGGACCAGAGAGACTGTACAACCATGGAGTTGAAAGTGACACTGATCTTCGCCATTGTGGCTGCTATTACCATCTCAGTTCACGCGCGAAGAGATCATAACAAACGGGGAGGAGGAAGAGGGAATGGCAGAGAGAGGGGACAAGGACGTTTCGGAGGAGGACCAGGCTCTGATAGACAGCAGATTGGTGGTCCTAAACAAGGCGGTGGGCCAATGGGTGGAAGGAGGGCTGATGGTCCTGAATCTGGTGGCCCCCTGATGG" ) != string::npos ) seed_.erase( seed_[i], i );
//    for ( int i = 0; i < seed_.size(); i++ )
//    {
//        bool cull = false;
//        if ( seed_[i]->seq_ == "ACCGTAACAAGACCGAAGATGGTGGCCACGACCATCATCATAACAAGATCGAAGAAGGGGACCAGGATAGACCCATGTCAGAGATGAAGCCATTCCGGTTCAATCCTTTTGGTCGCAAGCCTTCCGGAGACCGTCCGTCTGGCAGATGCAACCATACCGAAGAGGGGCACCCCAGGCGCGATGGTCACCAACATCCCCATGG" ) cull = true;
//        if ( seed_[i]->seq_ == "CCTAAACAAGGCGGCGGACCAATGGGTGGACGGAGGTTTGATGGTCCTGAATCTGGCGGCCCAATAATGGATGGACGTAGACCAAATGGTGGGCCAAGGGGTGGACGGAGGTTTGATGGTCCTGAATCTGATCGCCCACGGATGGATGGACCTGGACATGGCGGACATAGACAAAATGGTGGACCAATGGGCGGTAGGAGGTACGATGGAACAGGATTTGCT" ) cull = true;
//        if ( seed_[i]->seq_ == "CCCTTATACAGTTCACGCACGAAGAGATCATAACGAACGGGGAGGAGGAAGAGGAAATGGCAGAGAGAGGGGACAAGGACGTTTCGGAGGAGGACCAGGCTCTGATAGACAGCAGATGGGTGGTCC" ) cull = true;
//        if ( seed_[i]->seq_ == "AGCTGATGGACAGCAGATGGGAGGTCCTAGACAAGGCGGTGGACCAATGGGTGGAAGGAGGTCTGATGGTCCTGAATCTGATGGACAGCAGATGGGAGGTCCTAGACAAGGCGGTGGACCAATGGGTGGAAGGAGGTCTGATGGTCCTGAATCTGGTGGCCA" ) cull = true;
//        if ( seed_[i]->seq_ == "GACTTTGGAATGAGAATGACAGTGAGGAACAAGCTACTACCACTAAAAGCGTAACTCCATCTACAGCTCCTGATATGCCTGAGGTGGTCGAGATTGATATTAACGAAATTGACAATAACATGATGG" ) cull = true;
////        if ( seed_[i]->seq_ == "GAGATGCCTGAGGTGGTCGAGATTGATATCAACGAAATCGACCAGGTTTAGAACTTGTCAAGACAATATCGCTGTTGTATATAATATTGACCACTGATGTCTATTCTCTCAACCATTTACCTGCGA" ) cull = true;
////        if ( seed_[i]->seq_ == "GAGATGCCTGAGGTGGTCGAGATTGATATCAACGAAATCGACCAGGTTTAGAACTTGTCAAGACAATATCGCTGTTGTATATAATATTGACCACTGATGTCTATTCTCTCAACCATTTACCTGCGG" ) cull = true;
////        if ( seed_[i]->seq_ == "GGAGATGCCTGAGGTGGTCGAGATTGATATCAACGAAATCGACCAGGTTTAGAACTTGTCAAGACAATATCGCTGTTGTATATAATATTGACCACTGATGTCTATTCTCTCAACCATTTACCTGCG" ) cull = true;
////        if ( seed_[i]->seq_ == "AACCATTTACCTGCAATAACCATTTGTGATAAAAACGTAGGCCTTCATGAGATTACAACTTTTACAGTTCTTCGTATATCATATTTTATTTATTTGTAAAATAAAGATATAACTCGAAACATAAAAAAAAA" ) cull = true;
//        if ( seed_[i]->seq_ == "" ) cull = true;
//        if ( cull )
//        {
//            seed_.erase( seed_[i], i );
//        }
//    }
    
    
//    Node::reverify( seed_ );
//    NodePath::create( bwt_, seed_, seeds, paths );
//    NodePath::print( paths );
//    loci = NodePath::resolve( bwt_, paths, seed_ );
//    assert( false );
//    
//    LocusFill::fill( bwt_, seed_ );
//    loci.clear();
    
    for ( int again = 1; again--; )
    {
        if ( !loci.empty() ) Node::prune( bwt_, seed_ );
        for ( Node* node : seed_.nodes ) node->setVerified();
        NodePath::create( bwt_, seed_, seeds, paths );
        loci = NodePath::resolve( bwt_, paths, seed_ );
        NodePath::setExtends( loci, seed );
        for ( int d : { 0, 1 } ) while ( seed[d].extend( bwt_, seed_, false, d ) ) again = 1;
    }
    NodePath::print( paths, loci );
    LocusExport( seed_, "/home/glen/LocassDump/dump48" );
    LocusFill::fill( bwt_, seed_ );
    loci.clear();
    assert( false );
}

void Seed::assemble( bool seeded )
{
    SeedExtend seed[2];
    if ( seeded )
    {
        for ( Node* node : seed_.nodes ) for ( int d : { 0, 1 } ) if ( node->isContinue( d ) ) seed[d].addFork( node );
        for ( int d : { 0, 1 } ) while ( seed[d].extend( bwt_, seed_, true, d ) );
    }
    else restart( seed );
    
    for ( bool good = true; good; )
    {
        for ( int d : { 0, 1 } ) while ( seed[d].extend( bwt_, seed_, false, d ) );
        Node::prune( bwt_, seed_ );
        if ( !restart( seed ) ) break;
    }
}

void Seed::cull( MappedSeqs& ms )
{
    Coords* coords;
    for ( int i = 0; i < ms.reads.size(); i++ ) 
    {
        bool bad = false;
        for ( Node* node : seed_.nodes ) if ( !node->bad_ && ( coords = node->getRead( ms.reads[i].id ) ) && ( bad = !coords->coords[2] ) ) break;
        if ( bad ) ms.reads.erase( ms.reads.begin() + i-- );
    }
}

void Seed::print( MappedSeqs& ms, ofstream* ofs )
{
    vector< pair<int, pair<ReadId, string> > > prints;
    
    Lib* lib;
    int matched = 0, len = 150;
    for ( ReadStruct& rs : ms.reads )
    {
        ReadId id = rs.id;
        int d, est = rs.coords[0];
        string seq = rs.seq;
        if ( ( lib = params.getLib( rs.id ) ) && lib->getPair( id, d ) )
        {
            string alt = bwt_.getSequence( id );
            if ( seq.size() > len )
            {
                if ( !d ) est += seq.size()-len;
                seq = d ? seq.substr( 0, len ) : seq.substr( seq.size()-len );
            }
            if ( alt.size() > len ) alt = d ? alt.substr( alt.size()-len ) : alt.substr( 0, len );
            int ol = mapSeqOverlap( seq, alt, 15, d );
            if ( ol )
            {
                matched++;
                seq = d ? seq + alt.substr( ol ) : alt + seq.substr( ol );
                if ( !d ) est -= alt.size() + ( ol ? -ol : 50 );
            }
        }
        int i = 0;
        while ( i < prints.size() && est > prints[i].first ) i++;
        prints.insert( prints.begin() + i, make_pair( est, make_pair( rs.id, seq ) ) );
    }
    
    cout << "Number of overlapping read pairs printed: " << matched << endl;
    
    if ( prints.empty() ) return;
    int base = prints[0].first;
    
    vector<string> seeds;
//    seeds.push_back( "AAAGGATATAAATTCAGCGTTCCAACTTGTAGTCGTCGGTAGACAGCTATCTGTGAAGAGCACGAAGAGGAAATAATTTAAGACCTGAGCGAATGCCATGGAATTGAAAGTAATATTGACCTTTGCCATTGTGGCTGCTGTTGCTATAACTATTACAGCTCATCCACGAAGAGGGGAAGACGAAGAG" );
//    seeds.push_back( "TTGTAGTCGTCGGTAGACAGCTATCTGTGAAGAGCACGAAGAGGAAATAATTTAAGACCTGAGCGAATGCCATGGAATTGAAAGTAATATTGACCTTTGCCATTGTGGCTGCTGTTGCTATAACTATTGCAGCTCATCCACGAAGAGGGGAAGACGAAGAGGG" );
//    seeds.push_back( "CAGCTATCTGTGAAGAGCACGAAGAGGAAATAATTTAAGACCTGACTGAATACCATGGAATTGAAAGTAATATTGACCTTTGCCATTGTGGCTGCTGTTGCTATAACTATTACAGCTCATCCACGAAGAGGGGAAGACGAAGAG" );
//    seeds.push_back( "CAACCATGGAGTTGAAAGTGACACTGATATTCGCCATTG" );
//    seeds.push_back( "CAACCATGGAGTTGAAAGTGACACTGATATTCGCCATTG" );
//    seeds.push_back( "CAACCATGGAGTTGAAAGTGACACTGATATTCGCCATTG" );
//    seeds.push_back( "CAACCATGGAGTTGAAAGTGACACTGATATTCGCCATTG" );
    int32_t coords[2];
    for ( string& s : seeds ) if ( ofs )
    {
        for ( int i = 0; i < prints.size(); i++ ) 
        {
            bool good = prints[i].second.second.find( s ) != string::npos;
            if ( s.find( prints[i].second.second ) != string::npos ) good = true;
            if ( mapSeqEnd( prints[i].second.second, s, 50, coords, 1 ) && !coords[0] ) good = true;
            if ( mapSeqEnd( prints[i].second.second, s, 50, coords, 0 ) && coords[1] == s.size() ) good = true;
            if ( !good ) continue;
            ( *ofs ) << ">" << prints[i].second.first << endl << string( prints[i].first - base, '-' ) + prints[i].second.second << endl;
            prints.erase( prints.begin() + i-- );
        }
    }
    
    len = ( len*4 ) / 5;
    for ( int i = 0; i < prints.size(); i++ ) 
    {
        vector< pair<int, pair<ReadId, string> > > ols{ prints[i] };
        for ( int j = 0; j < ols.size(); j++ ) for ( int k = i+1; k < prints.size(); k++ )
        {
            if ( mapSeqOverlap( ols[j].second.second, prints[k].second.second, len )
                    || ols[j].second.second.find( prints[k].second.second ) != string::npos
                    || prints[k].second.second.find( ols[j].second.second ) != string::npos )
            {
                ols.push_back( prints[k] );
                prints.erase( prints.begin() + k-- );
            }
        }
        if ( ofs ) for ( pair<int, pair<ReadId, string> > p : ols ) (*ofs) << ">" << p.second.first << endl << string( p.first - base, '-' ) + p.second.second << endl;
        else for ( pair<int, pair<ReadId, string> > p : ols ) cout << ">" << p.second.first << endl << string( p.first - base, '-' ) + p.second.second << endl;
        
    }
}

bool Seed::restart( SeedExtend seed[2] )
{
    static bool initial = true;
    
    for ( int d : { 0, 1 } ) seed[d].reset();
    
    vector<NodePath*> seeds, paths, stops[2], branches[2];
    
    vector<string> culls[3], vers, pures, adds, qq[2], solidify, ext[2], base[2], tar[2], blunts[2];
    for ( string& s : pures ) seed_.purify( s );
    for ( string& s : solidify ) seed_.solidify( s );
    
    vector<Node*> destroy;
    for ( string& s : adds ) for ( int i = 0; i < seed_.size(); i++ ) if ( s.find( seed_[i]->seq_ ) != string::npos ) destroy.push_back( seed_[i] );
    for ( Node* node : destroy ) seed_.erase( node );
    for ( string& s : adds ) SeedNode::seed( bwt_, seed_, s, 99, 99, 1 );
//    Node::prune( bwt_, seed_ );
    
    for ( int d : { 0, 1 } ) for ( string& s : qq[d] ) for ( Node* node : seed_.nodes ) if ( s.find( node->seq_ ) != string::npos  )
    {
        node->setAllPairs();
        cout << ">Node " << node->coverage_ << endl << node->seq_ << endl;
        for ( NodeMark& nm : node->pe_[d] ) if ( ( d ? node->ends_[1] - nm.coords[1] : nm.coords[0] - node->ends_[0] ) < 500 ) cout << ">" << nm.id << endl << bwt_.getSequence( nm.id ) << endl;
        int x = 0;
    }
    for ( int d = 0; d < 3; d++ ) for ( string& s : culls[d] ) for ( Node* node : seed_.nodes ) if ( s.find( node->seq_ ) != string::npos )
    {
        if ( d < 2 ) node->cleanEdges( seed_, d );
        seed_.erase( node );
        break;
    }
    Node::updateStates( seed_ );
//    assert( seed_.cull( "CGGAAAGCACCCTAAACACGTGATTTCCAAATGCTGAAAAAATGCACCCCTAACAAGTATTGGCGTGTGAAATCCTACCCTTAACAAGTATTGGAAACAAAACCGTTTTTGGCAAATATTTCGCCTTGAAATGAACCCCTAAACAAGTAA" ) );
//    Bud::test( bwt_, "AGCATTACTTTTTCATATAGATAAAGTTTGGTGACAGCCCGCTAGACCGAAGGTCCGCGAGACCGAGGGTCCGCTAGACCGAGGGTCCGCGAGACCGAGGGTCCGCTAGACCGAAGGTCCGCGAGACCGAGGGTCCGCGAGACCGAAGGCCCGCGAGGCCGAAGGTCCACGAGACCGACTTTTTTCCCTTCATCGGTTGCCGCGGTCGAGTGGTTTAAGGCGCCTCGTTGCAGCGGCGTAACAGGTACTTGTTATCAGAGGGGGCGGGGGGGGGGCAGATCCAAATATTTCCGAGTAATTCCGTCTAGAATTAATTGCTAGGATTTCAGTTTAGGTGGCAGAAGTGAGCACGGGAAGCGCGTTCGGGGGGGGGGGGTGCAGGGAGGGGAGTGTACCTACCCCCCCGCGAGATAAGAGGAAGGTCCGACGTTTCTTAATACCCCCATCTCACTAGATAGCGATTCCGCTGCGATCGTCAAAAATCGACAAAGCGTGGTATATCGTAGCTTGAACGTGGCTTGATCTTTTCAATCTTCTCCGTCGTACCTTGATCTTTTCTGAAGCGGCAAGGATCGCCACCATCTTCCTGGTCG", 1 );
//    Bud::test( bwt_, "AGCATTACTTTTTCATATAGATAAAGTTTGGTGACAGCCCGCTAGACCGAAGGTCCGCGAGACCGAGGGTCCGCTAGACCGAGGGTCCGCGAGACCGAGGGTCCGCTAGACCGAAGGTCCGCGAGACCGAGGGTCCGCGAGACCGAAGGCCCGCGAGGCCGAAGGTCCACGAGACCGACTTTTTTCCCTTCATCGGTTGCCGCGGTCGAGTGGTTTAAGGCGCCTCGTTGCAGCGGCGTAACAGGTACTTGTTATCAGAGGGGGCGGGGGGGGGGCAGATCCAAATATTTCCGAGTAATTCCGTCTAGAATTAATTGCTAGGATTTCAGTTTAGGTGGCAGAAGTGAGCACGGGAAGCGCGTTCGGGGGGGGGGGGTGCAGGGAGGGGAGTGTACCTACCCCCCCGCGAGATAAGAGGAAGGTCCGACGTTTCTTAATACCCCCATCTCACTAGATAGCGATTCCGCTGCGATCGTCAAAAATCGACAAAGCGTGGTATATCGTAGCTTGAACGTGGCTTGATCTTTTCAATCTTCTCCGTCGTACCTTGATCTTTTCTGAAGCGGCAAGGATCGCCACCATCTTCCTGGTCGCCACAAAATTTTGAACATGTTCAAAACTTACGTAGCGAGATCGCGGAGGTGCTAAAGCGCATCATAATCGAGTCAAGAGCGTATTACAAACGACATATTCGTAGCTGTAACGAAGCTCCAACGCACAAAGCGTAGTAGAAGCGCATTAAAAGCGGC", 2 );
    
//    assert( seed_.cull( "" ) );
//    assert( seed_.cull( "" ) );
//    assert( seed_.cull( "" ) );
//    SeedNode::seed( bwt_, seed_, "AGACCGAAGGTCCGCGAGACCGAGGGTCCGCGAGACCGAAGGCCCGCGAGGCCGAAGGTCCACGAGACCGACTTTTTTCCCTTCATCGGTTGCCGCGGTCGAGTGGTTTAAGGCGCCTCGTTGCAGCGGCGTAACAGGTACTTGTTATCAGAGGGGGCGGGGGGGGGGCAGATCCAAATATTTCCGAGTAATTCCGTCTAGAATTAATTGCTAGGATTTCAGTTTAGGTGGCAGAAGTGAGCACGGGAAGCGCGTTCGGGGGGGGGGGGTGCAGGGAGGGGAGTGTACCTACCCCCCCGCGAGATAAGAGGAAGGTCCGACGTTTCTTAATACCCCCATCTCACTAGATAGCGATTCCGCTGCGATCGTCAAAAATCGACAAAGCGTGGTATATCGTAGCTTGAACGTGGCTTGATCTTTTCAATCTTCTCCGTCGTACCTTGATCTTTTCTGAAGCGGCAAGGATCGCCACCATCTTCCTGGTCGCCACAAAATTTTGAACATGTTCAAAACTTACGTAGCGAGATCGCGGAGGTGCTAAAGCGCATCATAATCGAGTCAAGAGCGTATTACAAACGACATATTCGTAGCTGTAACGAAGCTCCAACGCACAAAGCGTAGTAGAAGCGCATTAAAAGCGGCACCGATCGCCCATGTTTTTCAGGCCCGTTCCCAAGACA", 149, 0, 1 );
    SeedNode::seed( bwt_, seed_, "AGACCGAAGGTCCGCGAGACCGAGGGTCCGCGAGACCGAAGGCCCGCGAGGCCGAAGGTCCACGAGACCGACTTTTTTCCCTTCATCGGTTGCCGCGGTCGAGTGGTTTAAGGCGCCTCGTTGCAGCGGCGTAACAGGTACTTGTTATCAGAGGGGGCGGGGGGGGGGCAGATCCAAATATTTCCGAGTAATTCCGTCTAGAATTAATTGCTAGGATTTCAGTTTAGGTGGCAGAAGTGAGCACGGGAAGCGCGTTCGGGGGGGGGGGGTGCAGGGAGGGGAGTGTACCTACCCCCCCGCGAGATAAGAGGAAGGTCCGACGTTTCTTAATACCCCCATCTCACTAGATAGCGATTCCGCTGCGATCGTCAAAAATCGACAAAGCGTGGTATATCGTAGCTTGAACGTGGCTTGATCTTTTCAATCTTCTCCGTCGTACCTTGATCTTTTCTGAAGCGGCAAGGATCGCCACCATCTTCCTGGTCGCCACAAAATTTTGAACATGTTCAAAACTTACGTAGCGAGATCGCGGAGGTGCTAAAGCGCATCATAATCGAGTCAAGAGCGTATTACAAACGACATATTCGTAGCTGTAACGAAGCTCCAACGCACAAAGCGTAGTAGAAGCGCATTAAAAGCGGCACCGATCGCCCATGTTTTTCAGGCCCGTTCCCAAGACAATTCTGCTACGCTGATTCCGTTTACGAGGCGACTGCCTTGCGCTCGACACGCTTCACACCGTTTCCACCACGACGCTTTCCCTTTGGGTGGCATGTGGCACA", 149, 0, 1 );
//    assert( seed_.purify( "" ) );
//    assert( seed_.purify( "" ) );
//    assert( seed_.purify( "" ) );
//    assert( seed_.purify( "" ) );
//    assert( seed_.purify( "" ) );
//    assert( seed_.purify( "" ) );
//    Node::remap( bwt_, seed_ );
//    Node::reverify( seed_ );
    seed_.verify( 1 );
    
    for ( int again = 1; again-- > 0; )
    {
        NodePath::create( bwt_, seed_, seeds, paths );
//        if ( again = PathCross::resolve( paths, seed_ ) )
//        {
//            Node::reverify( seed_ );
//            Node::prune( bwt_, seed_ );
//        }
    }
    NodePath::print( paths );
    Traverse::trace( bwt_, paths[0]->path_.back(), 1 );
    paths[0]->path_.back()->extendFork( bwt_, seed_, 1000, 1, 1 );
    for ( Node* node : paths[0]->path_ ) node->bad_ = false;
    for ( Node* node : paths[0]->path_ ) node->setVerified();
    for ( Node* f : Nodes( paths[0]->path_.back(), 1, true ).nodes ) f->bad_ = false;
    for ( Node* f : Nodes( paths[0]->path_.back(), 1, true ).nodes ) f->setVerified();
    NodePath::create( bwt_, seed_, seeds, paths );
    NodePath::print( paths );
//    vector<NodePath*> buds = { paths[33] };
//    Bud::bud( bwt_, seed_, buds, 0 );
    assert( false );
//    LocusExport( seed_, "/home/glen/Genomes/LvDump/dump10" );
//    seed_.erase( paths[183]->path_[0] );
//    assert( seed_.solidify( "" ) );
//    assert( seed_.solidify( "" ) );
//    assert( seed_.purify( "" ) );
//    vector<NodePath*> buds = { paths[156] };
//    Bud::bud( bwt_, seed_, buds, 0 );
//    assert( false );
    
//    for ( Node* node : paths[85]->path_ ) for ( NodeMark& nm : node->pe_[0] ) cout << ">" << nm.id << endl << bwt_.getSequence( nm.id ) << endl;
//    for ( NodeMark& nm : paths[30]->path_[8]->pe_[1] ) cout << ">" << nm.id << endl << bwt_.getSequence( nm.id ) << endl;

//    seed_.verify( 1 );
//    NodePath::create( bwt_, seed_, seeds, paths );
//    NodePath::print( paths );
//    for ( int i = 17; i < paths[162]->path_->size(); i++ ) for ( NodeMark& nm : paths[162]->path_[i]->pe_[0] ) cout << ">" << nm.id << endl << bwt_.getSequence( nm.id ) << endl;
//    LocusExport( seed_, "/home/glen/Genomes/LvDump/dump7" );

    if ( initial )
    {
        initial = false;
    }
    
    
    unordered_set<NodePath*> used;
    for ( int d : { 0, 1 } )
    {
        unordered_set<NodePath*> pathed;
        for ( NodePath* np : seeds ) np->setEnds( pathed, stops[d], d );
        used.insert( pathed.begin(), pathed.end() );
    }
    
    vector<Node*> ends[2];
    Nodes extable[2], leapable[2];
    for ( int d : { 0, 1 } ) for ( NodePath* np : branches[d] ) np->setBranch( used, extable[d], leapable[d], d );
    for ( int d : { 0, 1 } ) for ( NodePath* np : stops[d] ) ends[d].push_back( np->getEnd( bwt_, seed_, used, d ) );
    for ( int d : { 0, 1 } ) for ( Node* node : leapable[d].nodes ) Leap::leapBranch( bwt_, seed_, node, d );
    for ( int d : { 0, 1 } ) for ( Node* node : ends[d] ) if ( node ) seed[d].addFork( node );
    for ( int d : { 0, 1 } ) for ( Node* node : extable[d].nodes ) seed[d].addAlt( node );
    
    
    seed[0].reset();
    
    if ( !seed[0].extend( bwt_, seed_, false, 0 ) && !seed[1].extend( bwt_, seed_, false, 1 ) )
    {
        assert( false );
    }
    for ( NodePath* np : paths ) delete np;
    
    return !seed[0].empty() || !seed[1].empty();
}

void Seed::fill()
{
    unordered_set<ReadId> used;
    for ( Node* node : seed_.nodes ) for ( auto& read : node->reads_ ) used.insert( read.first );
    
    Nodes filled[3];
    for ( Node* node : seed_.nodes ) if ( !node->bad_ && ( node->size() > params.readLen * 1.5 ) && ( node->countReads( true ) > 15 ) )
    {
        for ( int d : { 0, 1 } ) filled[d].fill( node, d, true );
    }
    for ( Node* node : filled[0].nodes ) if ( filled[1].find( node ) ) filled[2] += node;
    for ( Node* node : filled[2].nodes ) for ( int d : { 0, 1 } )
    {
        vector< pair<ReadId, string> > reads;
        for ( NodeMark& nm : node->pe_[d] ) if ( used.find( nm.id ) == used.end() ) reads.push_back( make_pair( nm.id, bwt_.getSequence( nm.id ) ) );
        if ( reads.empty() ) continue;
        
        Nodes fwd( node, params.maxPeMean - params.readLen, d, true );
        for ( pair<ReadId, string> read : reads )
        {
            int ol = 0;
            int32_t coords[2];
            vector<string> exts;
            for ( Node* f : fwd.nodes ) if ( mapSeqEnd( read.second, f->seq_, max( 50, ol ), coords, !d ) && ( d ? coords[1] : f->size() - coords[0] ) >= params.readLen )
            {
                if ( coords[1] - coords[0] > ol ) exts.clear();
                ol = coords[1] - coords[0];
                if ( ol == params.readLen ) assert( false );
                if ( ol == params.readLen ) break;
                string ext = f->seq_.substr( d ? coords[1]: coords[1]-params.readLen+1, params.readLen-1-ol );
                for ( string e : exts ) if ( e == ext ) ext.clear();
                if ( !ext.empty() ) exts.push_back( ext );
            }
            if ( ol == params.readLen || exts.empty() || !bwt_.isExtendable( read.second, max( 50, 1+params.readLen-ol ), d ) ) continue;
            int ext = exts[0].size();
            for ( int i = 0; i < ext; i++ ) for ( int j = 1; j < exts.size(); j++ ) if ( exts[0][i] != exts[j][i] ) ext = i;
            ol += ext;
            string seq = d ? exts[0].substr( exts[0].size()-ext ) : exts[0].substr( 0, ext );
            SeedNode::seed( bwt_, seed_, seq, d?ol:0, d?ol:0, d );
            assert( false );
        }
    }
}

bool Seed::flip()
{
    for ( Node* f : seed_.getGraph( 2 ).nodes )  for ( auto& read : f->reads_ ) if ( !read.second.redundant )
    {
        ReadId id = params.getRevId( read.first );
        bool found = false;
        for ( Node* node : seed_.nodes ) if ( node->getRead( id ) )
        {
            Nodes alls[2], base[2], tar[2];
            for ( int d : { 0, 1 } ) base[d].fill( f, d, true );
            for ( int d : { 0, 1 } ) tar[d].fill( node, d, true );
            vector<Node*> forks[2];
            for ( int d : { 0, 1 } ) for ( Node* fork : base[d].nodes ) if ( tar[d].find( fork ) )
            {
                for ( Edge& e : fork->edges_[!d] ) if ( base[d].find( e.node ) && !tar[d].find( e.node ) ) forks[d].push_back( fork );
            }
            bool finds[2]{ base[1].find( node ), base[0].find( node ) };
            base[0] += base[1];
            vector<Node*> met = base[0].meet( node );
            if ( !met.empty() ) continue;
            alls[0].fillAll( f );
            alls[1].fillAll( node );
            vector<Node*> flips[2];
            for ( int d : { 0, 1 } ) for ( Node* node : alls[d].nodes ) if ( node->drxn_ == 2 ) flips[d].push_back( node );
            seed_.flip( alls[ alls[1].size() < alls[0].size()] );
            return true;
        }
        if ( found ) break;
    }
    return false;
}

//void Seed::assemble()
//{
//    for ( Node* node : nodes_ )
//    {
//        node->resetMarks();
//        node->setCoverage();
//    }
//    
//    NodeList extendNodes[2] = { nodes_, nodes_ }, dummy;
//    int32_t limits[2] = { (int32_t)seq_.length() - params.maxPeMean - params.readLen, params.maxPeMean + params.readLen };
//    
//    while ( !extendNodes[0].empty() || !extendNodes[1].empty() )
//    {
//        extendNodes[0].clear();
//        extendNodes[1].clear();
//        NodeSet seedSet, delSet;
//        for ( Node* node : nodes_ )
//        {
//            if ( node->isSeed( seq_.length() ) )
//            {
//                seedSet.insert( node );
//            }
//        }
//        
//        Node::seedGetExtend( extendNodes, seedSet, delSet, limits );
//        
//        for ( int drxn : { 0, 1 } )
//        {
//            for ( Node* node : extendNodes[drxn] )
//            {
//                if ( !node->clones_ )
//                {
//                    ExtVars ev( nodes_, dummy, validLimits_, bwt_, false, false );
//                    ev.ante = node->getDrxnNodes( !drxn );
//                    node->extendCount_ = 1 + ( ( params.maxPeMean * 3 ) / params.readLen );
//                    node->extendNode( ev, drxn );
//                    ends_[drxn] = drxn ? max( ends_[1], node->ends_[1] ) : min( ends_[0], node->ends_[0] );
//                }
//            }
//        }
//        
//        Node::seedValidate( seedSet, delSet, validLimits_, ends_ );
//        for ( Node* del : delSet )
//        {
//            nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
//            delete del;
//        }
//    }
//}

//void Seed::assembleGraph()
//{
//    vector<string> seqs, badSeqs, startSeqs, clearSeqs, clearEdges, truncate;
////    Sd
////    {
////        seqs.push_back( "GGAGTATATGATTATTTCATTTTTTTATAATATGCATTTTCAAATTGTTCGTTTCACAAAATAATTAATTTTATTTCTTAAGCCTACACCAATCCGTTTATTGAATGACAAAAGATGTGTATTATGTTTATGTTCAATCTGGTATTCAAATTTAATTCGGTAAGAATAAGCGTTTTGAACATCGAACCACACGTACA" );
////        seqs.push_back( "GGAGTATATGATTATTTCATTTTTTTATAATATGCATTTTCAAATTGTTCGTTTCACAAAATAATTAATTTTATTTCTTAAGCCTACACCAATCCGTTTATTGAATAACAAAAGATGTGTATTATGTTTATGTTCAATCTGGTATTCAAATTTAATTCGGTAAGAATAAGCGTTTTGAACATCGAACCACACGTACAGAG" );
////        seqs.push_back( "GAAACCGTAATGAAATCTTCACCACCTGAAGTGGTCGAGATCGCAGTCAATGAAGAAGACATCAATGTGATCGCCGAGGTGTAATTTTTTTTTTTTTTTTTTTTTTAATGTAAATATACTAATAGTCACCAATGTCATTCCTTGGACAGAATGCATAATTATCGTTCCAGAGACAACATCATGG" );
////        
////        badSeqs.push_back( "GGCGATATCTGAATTTTAAGGGGAAAAATGCATTTCTAG" );
////        badSeqs.push_back( "CGTCCCGGACCGTCCCGGACCATCCCGGACCCCAT" );
////        badSeqs.push_back( "GCAAGAGAGAGAGAGAGAGAGAGAGAGGGGGAGATGG" );
////        badSeqs.push_back( "GGTGGATAGAATAGAGATAGAGAGTGATAGATCTATACATCGAG" );
////        badSeqs.push_back( "TGATTTTGAATTTGAATGCCCGAATTC" );
////        badSeqs.push_back( "TGACACATATTTGTAAAGCT" );
////        badSeqs.push_back( "GATCCTCCTCTTTCATTTGAT" );
////        badSeqs.push_back( "GCTTGGAATAACAAAGTTCTGTTTGTTTAAAGTAGTGCTTGCAT" );
////        
////        startSeqs.push_back( "AAAAAAAATGCCCACGGCTAACCTAAAAAAGGTGACCTTACCCCACCCAACGAAGAAAAATGTCAAATTTTCAAATTTCGAGACCAAAAAAAAAAAAAAAGGGTCCGCTTAAGAAGGTTCGATGGACCTGGATTTGGTGCCCCACAGATGGATGGACGGAGACATAATAGCGGTCCGATGGGTGGTAGGAGATTCGACGGACCTGGATTTGGTGGCTCCAGACCAGATAGTGCTGGAGGAAGACCTTTCTTTGGCCAAGGAGGCAGGTATGCTGGTGGAGAAGAAGAAGAAACTGATGCTGCCCAACAAAGTGGTGATGGAGAAGAAGAAGAAACTGATGCTGCCCAACA" );
////    }
//    
////    Af
////    {
////        seqs.push_back( "ACAATATTCATTGATCGCTGTCCCGGTCGAAATCGAGTCTAAAGTTATGAAATTATAAATAATAAACTCTTTATGTTGTATAGCTCCCTCTCAAGGTCAAGATATCAATTTAGTTATGAAACAATGTTTCTTGACTAAAATCGACCGATTAGCTATT" );
////        seqs.push_back( "CACAACTACATATTTCTGAGACTTTCATGGTCGCAGCTACCATTGTACGGCCGTAGGCGCTATTTGAAATTTATCAACATTTACCGGGAGCACCGGTAAATTTTTATATTCAATTTGTCTCGAAAACAGAAGTGTTTTAGGGGTTAATTTTTTGTGTGCTTATAATACATCCACAGGGGCGTCCAAAT" );
////        seqs.push_back( "CGTAGGCGCTATTTGAAATTTATCAAAATTTACCGGGAGCACCGGTAAATTTTTATATTCAATTTTTCTCGAAAACAGAAGTTTTTTAGGGGTTATTTTTGTTGGTGCTTATAATACATCTACAGGGGCATCCAAATATATCAAAATTGTAAAGATTCACTATCCACCCAATTTTCAAAATTTCCCAAAACTTTGTATTGCCTTTAATTAATTCGATTATTTCATTACCCCTTACAGCTCACGCACGAAGAGATTTCAATGAACAGCGAGGAGAGGAGAATGGCAGAGAGA" );
////        seqs.push_back( "CGTAGGCGCTATTTGAAATTTATCAAAATTTACCGGGAGCACCGGTAAATTTTTATATTCAATTTTTCTCGAAAACAGAAGTTTTTTAGGGGTTATTTTTTTTGGTGCTTATAATACATCTACAGGGGCATCCAAATATATCAAAATTGTAAAGATTCACTATCCACCCAATTTTCAAAATTTCCCAAAACGTTGTATTGCCTTTAATTAATTCGATTATTTCATTACCCCTTACAGCTCACGCACGAAGAGATTTCAATGAACAGCGAGGAGAGGAGAATGGCAGAGAGA" );
////        
////        badSeqs.push_back( "TTGTCCCGGTCGAAATCGAGTCTAAAGTTATGAAATTATAAATAATAAACTCTTTATGTT" );
////        badSeqs.push_back( "GATGTTCTTCCCAGCAGCTTCTCTTGTCAGGATGCC" );
////        badSeqs.push_back( "AAGATACATACATGATCTAAGTTTGGCAGGAATCC" );
////        badSeqs.push_back( "TTTGAATTTGAATGCCCGAATTCACCTGCGAGCG" );
////        badSeqs.push_back( "TATTTAGAATATTTCGTTTGCATTAAATTAAGCTTTGGAAGAAAAAAAATCCCAAA" );
////        badSeqs.push_back( "ATGGGAGTGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGC" );
////        badSeqs.push_back( "AAGAGAGAGAGAGAGAGAGAGAGATTAACCATAAT" );
////        badSeqs.push_back( "GAAGATGAATAAGAGTGCCGATAAATGAAATTGATTGCTCGGATTTTGAATG" );
////        badSeqs.push_back( "ATTGTATTGAATAGTATACAAATATAC" );
////        badSeqs.push_back( "GTTTTAACCAGTACCTCGAATAGAAGGCATGGTATATATGTA" );
////        badSeqs.push_back( "CAGTACCTCGAATAGAAGGCATGGTATATTTC" );
////        badSeqs.push_back( "TTAACCAGTACCTCGAATAGAAGGCATGGTATATTTGC" );
////        badSeqs.push_back( "GTGACACCAATATGCAAATGGTACTTTTGGAGTCACCGAGGCCCATAGTTTTAACCAGTACCTCGAATAGAAGGCATGGTATATTTGTTTTATATCCT" );
////        badSeqs.push_back( "TTATAGTTCATTGGATTAGATCAAGGTTTAGGCCT" );
////        badSeqs.push_back( "TTCCGTTTAGAAATTATATATCACAACTACATATTTCTGAGACTTTCATGGTC" );
////        badSeqs.push_back( "GAATTGTAACTGGTATGAGTTGATACGAACGAGGGGCTGTTTT" );
////        
////        startSeqs.push_back( "AAACTCTTTATGTTGTATAGCTCCCTCTCAAGGTCAAGATATCAATTTAGTTATGAAACAATGTTTCTTGACTAAAATCGACCGATTAGCTATTCCCGTTTGGAAATTATATATCACAACTACATATTCTGAGATTTTCATGGTCGCAGCTACCATTGTACGGCCGTAGGCGCTATTTGAAATTTATCAAAATTTACCGGGAGCACCGGTAAATTTTTATATTCAATTTTTCTCGAAAACAGAAGTTTTTTAGGGGTTATTTT" );
////        
////        clearSeqs.push_back( "AATGAGAATGGTAGAGAGAGAGGACAAGGTCGCTTTGGAGGAAGGCCTGGTGGAATGCAGAATGGTGGACCAAGGCAAGATGGTGGACCAATGGGTGGA" );
////    }
//    
////    Sp
////    {
////        seqs.push_back( "TTTAATGAAAAAATGTTTCATGTTTCTTTATACCATCTGGTTTTCAAATTCAATTAGGTTAGAATTAGGCGTTTTGAATATCCAACCGCATGCATGCATTCTGGTAGGATTAAACATAACGGGACCAGAATGATATCGAAGAAGACAAATCTCCTTGCATAAACCAATGTAGAGCAAAAGATAATAATGCAGAATTTGATTAATTTGATTAATCAATTCGACTATTTCATTATCCCCTACAGCTCACGCACGAAGAGATTTCAATGAACGGCGAG" );
////        seqs.push_back( "AGGCGTTTTGAATATCGAACCGCATGTACAAAGTAGGGTTAAACATAACGGGAAAAAGCATCAAAGAAGGCGAGTTATCATTCTTATTTGTCACCTGCCATAAACCAATGTAGAGCTAAAGATAATAATGCAGAATTTGATAACTTAATTAATTCGATTTTTTTTCATTACCCCCTACAGCTCACGCACGGAGAGATTTCAATGAACGGCGAGGAAATGAGAATGGCAGAGAGAGAGGACAAGGTCGCTTTGGAGGAAGGCCTGGTGGAATGCAGATGGGTGGACC" );
////        seqs.push_back( "CCTAGCATAAACCAATGTAGAGCTAAAGATAATAATGCAGAATTTGATTAATTAGTTCGATGATTTCATTACCCCTTACAGCTCACGCACGGAGAGATTTCAATGAACGGCGAGGAAAGGAGAATGGTAGAGAGAGGGGACAAGGTCGCTTTGGAGGAAGGCCTGATGGAATGCAGATGGGTGGACCAAGGCAAGATGGTGGACCAATGGGTGGAAGGAGGTTCGATGGACCTGGATTCGGTGCCCCACAGATGGATGGACGGAGACAAAAT" );
////        
////    }
//    
////    Mn
////    {
////        seqs.push_back( "TATGTACTCTGGTAAGGGAAACAAACATCGAAGGAGCCGAGTTACCAATCAGGGGGGTGTTTCACAAAGATCCTAAGTTGAACTTATCTCTAAGTTGGACTTAACAATTACGGAAAGCCGTTAGCATCTATGAAAATATTTTGTCAGAGTTATTTTAAAAGGCATACGGTGTTGATTCAAATTATTCATATTTTCTATTATCAAGGAATCTTCATGTTCTTGATGTGGAATTTATAAAAAGTCTAAAATACTTGAATTTTCGCTTTTGAATATATTTTATTTTAAGGCTACAAATGGCTCTCCATAATGTTTAAGTCCAGCTTAACAGTTAAGTTTGACTTAGGATCTTTGTGAAACACCCCCCTGATGTGTCACCTACCATAAACCAATGTATAGCTCAAGATG" );
////    }
//    
////    Mf
//    {
//        seqs.push_back( "TTTCTTCTTAGGCCTACATCATTCCGTTTTTCGAATGAAAAAAAGAAATGGAAAATACGCAATGTCTTATGTTTCTGTACCATCTGGGGGTGTTTCACAAAGACTAAGATCGACTTTAGGTTGGACTTAACTATGACAGGCCATTCATGCCACTGGTTGCTCTCACTATTGGTCTCTCAATTATGATCTTTAGATAATATTTTTTCCTTCAATAATGAACTATCGAAAAACGTAAATTCTGTAGATTTGTGACCACGATTTATAAGTCCCTTGAGAGACCGATGCCAGTATTATGGTGAGAGCAAACGGTGGCATGAACAGCC" );
//        seqs.push_back( "GCTGTGCTGCCATTGCGGCGCTCTAGCATTTCAAGGAATTAATTCCTGCCAGGTACCCATTTACATCACCTGGGTTGAGTGTGGCAAATGTAGATCAATGCC" );
//        seqs.push_back( "GTGCTGCCATTGCGGCGCTCTAGCATTTCAAGGAATTAATTCCTGCCAGGTACCCATTTACATCACCTGGGTTGAGTGTGGCAAATGTAGATCAATGCCTTGCCAAAGGACATTAGTGCCGCGGTGGGACTCGAACCGCGGACCTTGTGATCGACAGTCAGGTGACG" );
//        seqs.push_back( "GTGCTGCCATTGCGGCGCTCTAGCATTTCAAGGAATTAATTCCTGCCAGGTACCCATTTACATCACCTGGGTTGAGTGTGGCAAATGTAGATCAATGCCTTGCCAAAGGATATTAGTGCTGCGGTGGGACTCGAACCGCGGACCTTGTGATCGACAGTCAGGTGACGCATCCAC" );
//        seqs.push_back( "GTGCTGCCATTGCGGCGCTCTAGCATTTCAAGGAATTAATTCCTGCCAGGTACCCATTTACATCACCTGGGTTGAGTGTGGCAAATGTAGATCAATGCCCTGCCAAAGGACATTAGTGCCGCGGTGGGACTCGAACCGCGGACCTTGTGATCGACAGTCAGGTGACGCATCCACT" );
//        
//        badSeqs.push_back( "ACCCCAACTCCTTTTTGTTCTTAATTGTATTAAATTTGTTG" );
//        badSeqs.push_back( "GGCTGACATTTTGATTGAATGGGGTAGGGA" );
//        badSeqs.push_back( "TAAACATCATCAGATTTTGAAATAAACTACGACTTAACTACGGAATTCAGACATTAGAAGTTTACTTATCTCGTATTCAACTCGCATTTTAGACTTTTAC" );
//        badSeqs.push_back( "TATTACTTGATAAGTGGTGAATATGAAGCCAACAAAGGGCTCACAGGAGTATATGATTTATTTCTTAAAATATGCATTATAATAATAATAATAACAACAT" );
//        badSeqs.push_back( "TAATAATAATAATAATAACAACATTTATAATGCGCACATATCCACCCCAAGAGATGCTCAAGGTGCCTAACAAGGTAAATACATTAAACATAAGAAAGTA" );
//        badSeqs.push_back( "TAATGCGCACATATCCACCCCAAGAGATGCTCAAGGTGCCTAACAAGGTAAATACATTAAACATAAGAAAGTACATGAGAAGTTACATAGAACGAAAATA" );
//        badSeqs.push_back( "CGGCGCTCTAGCATTTCAAGGAATTAATTCCTGCCAGGTACCCATTTACATCACCTGGGTTGAGTGTGGCAAATGTAGATCAATATCTTGCCAAAGGACA" );
//        badSeqs.push_back( "GGGGAACTTTATTCAAACAAGGATTAATAGTCCAGAACGGTACTACGAAACGATCTTTCTCTAAAAATTTCATGAAACAGGGATTAAGAGTCGAGAACGG" );
//        badSeqs.push_back( "TTCAAACAAGGATTAATAGTCCAGAACGGTACTACGAAACGATCTTTCTCTAAAAAATTCATGAAACAGGGATTAAGAGTCGAGAACGGTGCAACGAATG" );
//        badSeqs.push_back( "GGGAACTTTATTCAAACAAGGATTAATAGTCCAGAACGGTACTACGAAACGATCTTTCTCTAAAAATTTCATGAAACAGGGATTAAGAGTCGAGAAAGTA" );
//        badSeqs.push_back( "GAAACAGGGATTAAGAGTCGAGAACGGTACAACGAATCGGCTTCGTTCTCTAATAATTTAATGAAACGTCAAAGATAAAAAAAGACAATAAATGAACCGG" );
//        badSeqs.push_back( "AAACAGGGATTAAGAGTCGAGAACGGTACAACGAATCGGCTTCGTTCTCTAATAATTTAATGAAGCGTCAAAGATAAAAAAGACAATAAATGAACCGGTC" );
//        badSeqs.push_back( "AGAGTCGAGAACGGTACAACGAATCGGCTTCGTTCTCTAATAATTTAATGACACGTCACAGCTTCCGCGCTACGCGCGAGAGGACGGAAACTGAATTCAA" );
//        badSeqs.push_back( "GGCTTCGTTCTCTAATAATTTAATGAAACGTCACAGCTTCCGCGCTACGCGCGAGAGGCCGGAAACTGAATTCAAACAAGGATTAAGAGAGAGTCAAAAA" );
//        badSeqs.push_back( "CAACGAATCGGCTTCGTTCTCTAATAATTTAATGAAACGTCACAGCTTCCGCGCTA" );
//        badSeqs.push_back( "GTTAATTGAGAGCTGAGTTGATGAGATAACCACATGTCAGTGACCACGCTGGGAGGTCTAATTCCTCCTGGCCAAACTAGTGCAGACAGGCTTTAAAGGT" );
//        badSeqs.push_back( "CTAGTGCAGACAGGCTTTAAGGTAACCTTTCTCATTGAGTGATACAAAACTGATAAAAATTGGGACGTTATCTTCTTTATGAAAGTCTACTGTCATTTTT" );
//        badSeqs.push_back( "CTAGTGCAGACAGGCTTTAAGGTAACCTTTCTCATTGAGTGATACAAAACTGATAAAAATTGGGACGTTATCTTCTTTATGAAAGTCTACTGTCATTTTT" );
//        badSeqs.push_back( "TGTATACTTTTGAAGGGGTTGACGTTCATAAGGCCAATTGCCCTTTTTGTTAACCCAAATTTTGTTTGCTCCTTCCCTTTGTATTTTTCTTTTAAATCGA" );
//        
//        clearEdges.push_back( "ATTTCTTCTTAGGCCTACATCATTCCGTTTTTCGAATGAAAAAAAGAAATGGAAAATACGCAATGTCTTATGTTTCTGTACCATCTGGGGGTGTTTCACA" );
//        clearEdges.push_back( "AGCTGTGCTGCCATTGCGGCGCTCTAGCATTTCAAGGAATTAATTCCTGCCAGGTACCCATTTACATCACCTGGGTTGAGTGTGGCAAATGTAGATCAAT" );
//        clearEdges.push_back( "TTCTATATTCATAATTGTGCTGCCACGGAGAAGCAGGCTTTTCAAAAAGTATTGTCAAAATCCTCTTGTACCAAGTTTGGGCTGAAGCAATGAATGTGGG" );
//        clearEdges.push_back( "TATATACAAATGTATATATATATATATATATGTAACGCCCTGGAAAATAGGCCTCCCTCAAAATTTTGTATATATATCTTTTTATTATTGTCACATATTT" );
//        clearEdges.push_back( "GAAAAATTGAAGAAGGAAAATCGATATTCACACACACACACACACATATATATATATATATATATATATATATATATATATATATATATATATAGAGAGA" );
//        clearEdges.push_back( "TCTATAGTCGAGTTTATCGATTCATTACTATTGGCTTACAGACTTACGATATGCATTTTATATTTCTATAGCACAGTCCTTAATATTTAAATTTGTTAAA" );
//        
//        // Truncate
//        truncate.push_back( "TGGATTGATCATGATATACATTGATTGAGGAAGAC" );
//        truncate.push_back( "ACATATTGCTTTTGTTAAGATATAATATATATATATATATA" );
//        truncate.push_back( "GAAGCCAGTGAAAGATGCATGAAAATACTGAC" );
//        truncate.push_back( "ACTGAAATCATACAAATACTTGTTTAGCGTTTACTAATTTAGCCCCTAGTATTTTCCTAACGCAGCTGGATTGATCATGTTATATATTTATTGAAGAAGA" );
//        
//
//    }
//    
//    for ( string seq : startSeqs )
//    {
//        MapNode* mn = new MapNode();
//        mn->seq = seq;
//        bwt_.mapSequence( mn->seq, mn->ids, mn->coords );
//        mn->recoil();
//        Node* base = new Node( mn, 0, mn->ids.size()-1, 2 );
//        delete mn;
//        nodes_.push_back( base );
//    }
//    
//    bool finished = false, anyContinue = false;
//    int continueCount = 0;
//    bool drxn = true;
//    for ( Node* node : nodes_ ) if ( !node->edges_[!drxn].empty() ) node->drxn_ = drxn;
//    NodeList island;
//    ExtVars ev( nodes_, island, validLimits_, bwt_, false, false );
//    
//    while ( !finished )
//    {
//        NodeSet delSet;
//        finished = true;
//        for ( int i = 0; i < 16; i++ )
//        {
//            for ( int j = 0; j < clearEdges.size(); )
//            {
//                bool doErase = false;
//                for ( Node* node : nodes_ )
//                {
//                    if ( node->seq_.find( clearEdges[j] ) == node->seq_.npos ) continue;
//                    for ( Node* nxt : node->getNextNodes( drxn ) ) nxt->dismantleNode( delSet, drxn );
//                    doErase = true;
//                    clearEdges.erase( clearEdges.begin() + j );
//                    break;
//                }
//                if ( !doErase ) j++;
//            }
//            deleteNodes( delSet );
//            
//            for ( int j = 0; j < seqs.size(); )
//            {
//                bool doErase = false;
//                for ( Node* node : nodes_ )
//                {
//                    int ol = mapSeqOverlap( node->seq_, seqs[j], 98 );
//                    if ( ol )
//                    {
//                        MapNode* mn = new MapNode();
//                        mn->seq = seqs[j];bwt_.mapSequence( mn->seq, mn->ids, mn->coords );
//                        mn->recoil();
//                        Node* newNode = new Node( mn, 0, mn->ids.size()-1, drxn );
//                        nodes_.push_back( newNode );
//                        ol -= seqs[j].find( newNode->seq_ );
//                        node->addEdge( newNode, ol, drxn, true );
//                        delete mn;
//                        seqs.erase( seqs.begin() + j );
//                        doErase = true;
//                        break;
//                    }
//                }
//                if ( !doErase ) j++;
//            }
//            NodeList extNodes;
//            for ( Node* node : nodes_ )
//            {
////                if ( node->seq_.find( "TGTAGATTTGAGACCACGATTTGTACGCCCCTTGAAAGACAAATGCCAGTAAA" ) != node->seq_.npos ) node->stop_[drxn] = 1;
//                if ( node->isContinue( drxn ) ) extNodes.push_back( node );
//            }
//            for ( Node* node : extNodes )
//            {
//                node->extendSeed( ev, drxn );
//                node->seedLoop( drxn );
//                ends_[0] = min( ends_[0], node->ends_[0] );
//                ends_[1] = max( ends_[1], node->ends_[1] );
//                finished = false;
//            }
//            
//            for ( Node* node : nodes_ )
//            {
//                bool badSeq = false;
//                for ( string seq : badSeqs ) if ( node->seq_.find( seq ) != node->seq_.npos ) badSeq = true;
//                if ( badSeq )
//                {
//                    node->dismantleNode( delSet, drxn );
//                }
//                else if ( node->edges_[drxn].size() > 4 ) node->seedValidate( delSet, drxn );
//            }
//            for ( string &seq : truncate )
//            {
//                for ( Node* node : nodes_ )
//                {
//                    if ( node->seq_.find( seq ) == node->seq_.npos ) continue;
//                    for ( Node* nxt : node->getNextNodes( drxn ) ) nxt->dismantleNode( delSet, drxn );
//                }
//            }
//            deleteNodes( delSet );
//        }
//        
//        for ( Node* node : nodes_ )
//        {
//            ends_[0] = min( ends_[0], node->ends_[0] );
//            ends_[1] = max( ends_[1], node->ends_[1] );
//            node->stop_[!drxn] = 1;
//        }
//        
////        Node::seedLeap( bwt_, nodes_, 0, 2000 );
////        for ( bool d : { 0, 1 } )
////        {
////            for ( int i = 0; i < 3; i++ )
////            {
////                NodeList extNodes;
////                for ( Node* node : nodes_ )
////                {
////                    node->drxn_ = d;
////                    if ( node->isContinue( d ) ) extNodes.push_back( node );
////                }
////                for ( Node* node : extNodes )
////                {
////                    node->extendSeed( ev, d );
////                    node->seedLoop( d );
////                    ends_[0] = min( ends_[0], node->ends_[0] );
////                    ends_[1] = max( ends_[1], node->ends_[1] );
////                }
////            }
////        }
//        
//        for ( Node* node : nodes_ )
//        {
//            if ( node->isContinue( drxn ) )
//            {
//                anyContinue = true;
//                continueCount++;
//            }
//            for ( string seq : clearSeqs )
//            {
//                if ( node->seq_.find( seq ) != node->seq_.npos ) node->clearEdges( 0 );
//            }
//        }
//        
//        for ( Node* node : nodes_ )
//        {
//            bool canErase = node->ends_[0] < 100 ;
//            for ( Node* prv : node->getNextNodes( !drxn ) )
//            {
//                if ( prv->ends_[0] < 100 ) canErase = true;
//            }
//            for ( Node* prv : node->getNextNodes( !drxn ) )
//            {
//                if ( !canErase ) break;
//                if ( drxn ? node->ends_[1] + 1000 < prv->ends_[0] : prv->ends_[1] < node->ends_[0] - 1000 )
//                {
//                    prv->removeEdge( node, drxn );
//                    node->removeEdge( prv, !drxn );
//                }
//            }
//            if ( node->edges_[!drxn].empty() || node->drxn_ != 2 ) continue;
//            node->drxn_ = drxn;
//        }
//        Node::seedValidate( nodes_, delSet, drxn );
//        deleteNodes( delSet );
//        merge();
//        break;
//    }
//    
//    
//    drxn = true;
//    NodeList ends;
//    int id = 0;
//    for ( Node* node : nodes_ )
//    {
//        node->id2_ = -1;
//        if ( node->edges_[!drxn].empty() ) ends.push_back( node );
//    }
//    
//    ofstream fp( "/media/glen/ssd/test.fa" );
//    PathSeed path( nodes_ );
//    path.plot( drxn );
//    path.exportAlign( fp );
//    
//    cout << ( anyContinue ? "Not complete " : "Complete " ) + to_string( continueCount ) << endl;
//    
////    for ( Node* node : ends )
////    {
////        node->setIds( id, drxn );
////    }
//    sort( nodes_.begin(), nodes_.end(), []( Node* const &a, Node* const &b ){ return a->id2_ < b->id2_; } );
//    for ( Node* node : nodes_  )
//    {
//        unordered_set<int> hits[2];
//        bool loop = false;
//        for ( Node* fwd : node->getDrxnNodes( drxn ) ) if ( fwd == node ) loop = true;
//        for ( bool d : { 0, 1 } )
//        {
//            for ( ReadMark &mark : node->marks_[d] )
//            {
//                for ( Node* n : nodes_ )
//                {
//                    auto it = n->reads_.find( mark.id );
//                    if ( it == n->reads_.end() || n == node ) continue;
//                    hits[d].insert( n->id2_ );
//                }
//            }
//        }
//        node->setCoverage();
//        fp << ">" + to_string( node->id2_ ) + " Cover " + to_string( (int)node->coverage_ );
//        fp << " Reads " + to_string( node->reads_.size() );
//        if ( !hits[0].empty() || !hits[1].empty() )
//        {
//            fp << " |";
//            for ( int hit : hits[1] ) fp << " " << hit;
//            fp << " |";
//            for ( int hit : hits[0] ) fp << " " << hit;
//            fp << " |";
//        }
//        fp << ( node->edges_[drxn].empty() ? " End" : ( loop ? " Loop" : "" ) ) << endl;
//        fp << string( node->ends_[0] - ends_[0], '-' ) + node->seq_ << endl;
//    }
////    for ( Node* n1 : nodes_ )
////    {
////        for ( ReadMark &mark : n1->marks_[!drxn] )
////        {
////            bool found = false;
////            for ( Node* n2 : nodes_ ) if ( n2->reads_.find( mark.id ) != n2->reads_.end() ) found = true;
////            if ( found ) continue;
////            if ( drxn ) fp << ">" + to_string( mark.id ) << endl << string( mark.estimate - 100, '-' ) + bwt_.getSequence( mark.id ) << endl;
////            else fp << ">" + to_string( mark.id ) << endl << string( mark.estimate - ends_[0] + 100, '-' ) + bwt_.getSequence( mark.id ) << endl;
////        }
////    }
//    fp.close();
//    exit(0);
//}

//void Seed::assembleHaploid()
//{
//    {
//        MapNode* mn = new MapNode();
//        mn->seq = "CGAATTCCGCGGATCCTTCTATAGTGTCACCTAAATGTCGACGGCCAGGCGGCCGCCAGGCCTACCCACTAGTCAATTCGGGAGGATCGAAACGGCAGATCGCAAAAAACAGTACATACAGAAGGAGACATGAACATGAACATCAAAAAAATTGTAAAACAAGCCACAGTTCTGACTTTTACGACTGCACTTCTGGCAGGAGGAGCGACTCAAGCCTTCGCGAAAGAAAATAACCAAAAAGCATACAAAGAAACGTACGGCGTCTCTCATATTACACGCCATGATATGCTGCAGATCCCTAAACAGCAGCAAAACGAAAAATACCAAGTGCCTCAATTCGATCAATCAACGATTAAAAATATTGAGTCTGCAAAAGGACTTGATGTGTGGGACAGCTGGCCGCTGCAAAACGCTGACGGAACAGTAGCTGAATACAACGGCTATCACGTTGTGTTTGCTCTTGCGGGAAGCCCGAAAGACGCTGATGACACATCAATCTACATGTTTTATCAAAAGGTCGGCGACAACTCAATCGACAGCTGGAAAAACGCGGGCCGTGTCTTTAAAGACAGCGATAAGTTCGACGCCAACGATCCGATCCTGAAAGATCAGACGCAAGAATGGTCCGGTTCTGCAACCTTTACATCTGACGGAAAAATCCGTTTATTCTACACTGACTATTCCGGTAAACATTACGGCAAACAAAGCCTGACAACAGCGCAGGTAAATGTGTCAAAATCTGATGACACACTCAAAATCAACGGAGTGGAAGATCACAAAACGATTTTTGACGGAGACGGAAAAACATATCAGAACGTTCAGCAGTTTATCGATGAAGGCAATTATACATCCGGCGACAACCATACGCTGAGAGACCCTCACTACGTTGAAGACAAAGGCCATAAATACCTTGTATTCGAAGCCAACACGGGAACAGAAAACGGATACCAAGGCGAAGAATCTTTATTTAACAAAGCGTACTACGGCGGCGGCACGAACTTCTTCCGTAAAGAAAGCCAGAAGCTTCAGCAGAGCGCTAAAAAACGCGATGCTGAGTTAGCGAACGGCGCCCTCGGTATCATAGAGTTAAATAATGATTACACATTGAAAAAAGTAATGAAGCCGCTGATCACTTCAAACACGGTAACTGATGAAATCGAGCGCGCGAATGTTTTCAAAATGAACGGCAAATGGTACTTGTTCACTGATTCACGCGGTTCAAAAATGACGATCGATGGTATTAACTCAAACGATATTTACATGCTTGGTTATGTATCAAACTCTTTAACCGGCCCTTACAAGCCGCTGAACAAAACAGGGCTTGTGCTGCAAATGGGTCTTGATCCAAACGATGTGACATTCACTTACTCTCACTTCGCAGTGCCGCAAGCCAAAGGCAACAATGTGGTTATCACAAGCTACATGACAAACAGAGGCTTCTTCGAGGATAAAAAGGCAACATTTGCGCCAAGCTTCTTAATGAACATCAAAGGCAATAAAACATCCGTTGTCAAAAACAGCATCCTGGAGCAAGGACAGCTGACAGTCAACTAATAACAGCAAAAAGAAAATGCCGATACTTCATTGGCATTTTCTTTTATTTCTCAACAAGATGGTGAATTGACTAGTGGGTAGATCCACAGGACGGGTGTGGTCGCCATGATCGCGTAGTCGATAGTGGCTCCAAGTAGCGAAGCGAGCAGGACTGGGCGGCGGCCAAAGCGGTCGGACAGTGCTCCGAGAACGGGTGCGCATAGAAATTGCATCAACGCATATAGCGCTAGCAGCACGCCATAGTGACTGGCGATGCTGTCGGAATGGACGATATCCCGCAAGAGGCCCGGCAGTACCGGCATAACCAAGCCTATGCCTACAGCATCCAGGGTGACGGTGCCGAGGATGACGATGAGCGCATTGTTAGATTTCATACACGGTGCCTGACTGCGTTAGCAATTTAACTGTGATAAACTACCGCATTAAAGCTTATCGATGATAAGCTGTCAAACATGAGAATTGATCCGGAACCCTTAATATAACTTCGTATAATGTATGCTATACGAAGTTATTAGGTCCCTCGACTATAGGGTCACCGTCGACAGCGACACACTTGCATCGGATGCAGCCCGGTTAACGTGCCGGCACGGCCTGGGTAACCAGGTATTTTGTCCACATAACCGTGCGCAAAATGTTGTGGATAAGCAGGACACAGCAGCAATCCACAGCAGGCATACAACCGCACACCGAGGTTACTCCGTTCTACAGGTTACGACGACATGTCAATACTTGCCCTTGACAGGCATTGATGGAATCGTAGTCTCACGCTGATAGTCTGATCGACAATACAAGTGGGACCGTGGTCCCAGACCGATAATCAGACCGACAACACGAGTGGGATCGTGGTCCCAGACTAATAATCAGACCGACGATACGAGTGGGACCGTGGTCCCAGACTAATAATCAGACCGACGATACGAGTGGGACCGTGGTTCCAGACTAATAATCAGACCGACGATACGAGTGGGACCGTGGTCCCAGACTAATAATCAGACCGACGATACGAGTGGGACCATGGTCCCAGACTAATAATCAGACCGACGATACGAGTGGGACCGTGGTCCCAGTCTGATTATCAGACCGACGATACGAGTGGGACCGTGGTCCCAGACTAATAATCAGACCGACGATACGAGTGGGACCGTGGTCCCAGACTAATAATCAGACCGACGATACGAGTGGGACCGTGGTCCCAGTCTGATTATCAGACCGACGATACAAGTGGAACAGTGGGCCCAGAGAGAATATTCAGGCCAGTTATGCTTTCTGGCCTGTAACAAAGGACATTAAGTAAAGACAGATAAACGTAGACTAAAACGTGGTCGCATCAGGGTGCTGGCTTTTCAAGTTCCTTAAGAATGGCCTCAATTTTCTCTATACACTCAGTTGGAACACGAGACCTGTCCAGGTTAAGCACCATTTTATCGCCCTTATACAATACTGTCGCTCCAGGAGCAAACTGATGTCGTGAGCTTAAACTAGTTCTTGATGCAGATGACGTTTTAAGCACAGAAGTTAAAAGAGTGATAACTTCTTCAGCTTCAAATATCACCCCAGCTTTTTTCTGCTCATGAAGGTTAGATGCCTGCTGCTTAAGTAATTCCTCTTTATCTGTAAAGGCTTTTTGAAGTGCATCACCTGACCGGGCAGATAGTTCACCGGGGTGAGAAAAAAGAGCAACAACTGATTTAGGCAATTTGGCGGTGTTGATACAGCGGGTAATAATCTTACGTGAAATATTTTCCGCATCAGCCAGCGCAGAAATATTTCCAGCAAATTCATTCTGCAATCGGCTTGCATAACGCTGACCACGTTCATAAGCACTTGTTGGGCGATAATCGTTACCCAATCTGGATAATGCAGCCATCTGCTCATCATCCAGCTCGCCAACCAGAACACGATAATCACTTTCGGTAAGTGCAGCAGCTTTACGACGGCGACTCCCATCGGCAATTTCTATGACACCAGATACTCTTCGACCGAACGCCGGTGTCTGTTGACCAGTCAGTAGAAAAGAAGGGATGAGATCATCCAGTGCGTCCTCAGTAAGCAGCTCCTGGTCACGTTCATTACCTGACCATACCCGAGAGGTCTTCTCAACACTATCACCCCGGAGCACTTCAAGAGTAAACTTCACATCCCGACCACATACAGGCAAAGTAATGGCATTACCGCGAGCCATTACTCCTACGCGCGCAATTAACGAATCCACCATCGGGGCAGCTGGTGTCGATAACGAAGTATCTTCAACCGGTTGAGTATTGAGCGTATGTTTTGGAATAACAGGCGCACGCTTCATTATCTAATCTCCCAGCGTGGTTTAATCAGACGATCGAAAATTTCATTGCAGACAGGTTCCCAAATAGAAAGAGCATTTCTCCAGGCACCAGTTGAAGAGCGTTGATCAATGGCCTGTTCAAAAACAGTTCTCATCCGGATCTGACCTTTACCAACTTCATCCGTTTCACGTACAACATTTTTTAGAACCATGCTTCCCCAGGCATCCCGAATTTGCTCCTCCATCCACGGGGACTGAGAGCCATTACTATTGCTGTATTTGGTAAGCAAAATACGTACATCAGGCTCGAACCCTTTAAGATCAACGTTCTTGAGCAGATCACGAAGCATATCGAAAAACTGCAGTGCGGAGGTGTAGTCAAACAACTCAGCAGGCGTGGGAACAATCAGCACATCAGCAGCACATACGACATTAATCGTGCCGATACCCAGGTTAGGCGCGCTGTCAATAACTATGACATCATAGTCATGAGCAACAGTTTCAATGGCCAGTCGGAGCATCAGGTGTGGATCGGTGGGCAGTTTACCTTCATCAAATTTGCCCATTAACTCAGTTTCAATACGGTGCAGAGCCAGACAGGAAGGAATAATGTCAAGCCCCGGCCAGCAAGTGGGCTTTATTGCATAAGTGACATCGTCCTTTTCCCCAAGATAGAAAGGCAGGAGAGTGTCTTCTGCATGAATATGAAGATCTGGTACCCATCCGTGATACATTGAGGCTGTTCCCTGGGGGTCGTTACCTTCCACGAGCAAAACACGTAGCCCCTTCAGAGCCAGATCCTGAGCAAGATGAACAGAAACTGAGGTTTTGTAAACGCCACCTTTATGGGCAGCAACCCCGATCACCGGTGGAAATACGTCTTCAGCACGTCGCAATCGCGTACCAAACACATCACGCATATGATTAATTTGTTCAATTGTATAACCAACACGTTGCTCAACCCGTCCTCGAATTTCCATATCCGGGTGCGGTAGTCGCCCTGCTTTCTCGGCATCTCTGATAGCCTGAGAAGAAACCCCAACTAAATCCGCTGCTTCACCTATTCTCCAGCGCCGGGTTATTTTCCTCGCTTCCGGGCTGTCATCATTAAACTGTGCAATGGCGATAGCCTTCGTCATTTCATGACCAGCGTTTATGCACTGGTTAAGTGTTTCCATGAGTTTCATTCTGAACATCCTTTAATCATTGCTTTGCGTTTTTTTATTAAATCTTGCAATTTACTGCAAAGCAACAACAAAATCGCAAAGTCATCAAAAAACCGCAAAGTTGTTTAAAATAAGAGCAACACTACAAAAGGAGATAAGAAGAGCACATACCTCAGTCACTTATTATCACTAGCGCTCGCCGCAGCCGTGTAACCGAGCATAGCGAGCGAACTGGCGAGGAAGCAAAGAAGAACTGTTCTGTCAGATAGCTCTTACGCTCAGCGCAAGAAGAAATATCCACCGTGGGAAAAACTCCAGGTAGAGGTACACACGCGGATAGCCAATTCAGAGTAATAAACTGTGATAATCAACCCTCATCAATGATGACGAACTAACCCCCGATATCAGGTCACATGACGAAGGGAAAGAGAAGGAAATCAACTGTGACAAACTGCCCTCAAATTTGGCTTCCTTAAAAATTACAGTTCAAAAAGTATGAGAAAATCCATGCAGGCTGAAGGAAACAGCAAAACTGTGACAAATTACCCTCAGTAGGTCAGAACAAATGTGACGAACCACCCTCAAATCTGTGACAGATAACCCTCAGACTATCCTGTCGTCATGGAAGTGATATCGCGGAAGGAAAATACGATATGAGTCGTCTGGCGGCCTTTCTTTTTCTCAATGTATGAGAGGCGCATTGGAGTTCTGCTGTTGATCTCATTAACACAGACCTGCAGGAAGCGGCGGCGGAAGTCAGGCATACGCTGGTAACTTTGAGGCAGCTGGTAACGCTCTATGATCCAGTCGATTTTCAGAGAGACGATGCCTGAGCCATCCGGCTTACGATACTGACACAGGGATTCGTATAAACGCATGGCATACGGATTGGTGATTTCTTTTGTTTCACTAAGCCGAAACTGCGTAAACCGGTTCTGTAACCCGATAAAGAAGGGAATGAGATATGGGTTGATATGTACACTGTAAAGCCCTCTGGATGGACTGTGCGCACGTTTGATAAACCAAGGAAAAGATTCATAGCCTTTTTCATCGCCGGCATCCTCTTCAGGGCGATAAAAAACCACTTCCTTCCCCGCGAAACTCTTCAATGCCTGCCGTATATCCTTACTGGCTTCCGCAGAGGTCAATCCGAATATTTCAGCATATTTAGCAACATGGATCTCGCAGATACCGTCATGTTCCTGTAGGGTGCCATCAGATTTTCTGATCTGGTCAACGAACAGATACAGCATACGTTTTTGATCCCGGGAGAGACTATATGCCGCCTCAGTGAGGTCGTTTGACTGGACGATTCGCGGGCTATTTTTACGTTTCTTGTGATTGATAACCGCTGTTTCCGCCATGACAGATCCATGTGAAGTGTGACAAGTTTTTAGATTGTCACACTAAATAAAAAAGAGTCAATAAGCAGGGATAACTTTGTGAAAAAACAGCTTCTTCTGAGGGCAATTTGTCACAGGGTTAAGGGCAATTTGTCACAGACAGGACTGTCATTTGAGGGTGATTTGTCACACTGAAAGGGCAATTTGTCACAACACCTTCTCTAGAACCAGCATGGATAAAGGCCTACAAGGCGCTCTAAAAAAGAAGATCTAAAAACTATAAAAAAAATAATTATAAAAATATCCCCGTGGATAAGTGGATAACCCCAAGGGAAGTTTTTTCAGGCATCGTGTGTAAGCAGAATATATAAGTGCTGTTCCCTGGTGCTTCCTCGCTCACTCGAGGGCTTCGCCCTGTCGCTCAACTGCGGCGAGCACTACTGGCTGTAAAAGGACAGACCACATCATGGTTCTGTGTTCATTAGGTTGTTCTGTCCATTGCTGACATAATCCGCTCCACTTCAACGTAACACCGCACGAAGATTTCTATTGTTCCTGAAGGCATATTCAAATCGTTTTCGTTACCGCTTGCAGGCATCATGACAGAACACTACTTCCTATAAACGCTACACAGGCTCCTGAGATTAATAATGCGGATCTCTACGATAATGGGAGATTTTCCCGACTGTTTCGTTCGCTTCTCAGTGGATAACAGCCAGCTTCTCTGTTTAACAGACAAAAACAGCATATCCACTCAGTTCCACATTTCCATATAAAGGCCAAGGCATTTATTCTCAGGATAATTGTTTCAGCATCGCAACCGCATCAGACTCCGGCATCGCAAACTGCACCCGGTGCCGGGCAGCCACATCCAGCGCAAAAACCTTCGTGTAGACTTCCGTTGAACTGATGGACTTATGTCCCATCAGGCTTTGCAGAACTTTCAGCGGTATACCGGCATACAGCATGTGCATCGCATAGGAATGGCGGAACGTATGTGGTGTGACCGGAACAGAGAACGTCACACCGTCAGCAGCAGCGGCGGCAACCGCCTCCCCAATCCAGGTCCTGACCGTTCTGTCCGTCACTTCCCAGATCCGCGCTTTCTCTGTCCTTCCTGTGCGACGGTTACGCCGCTCCATGAGCTTATCGCGAATAAATACCTGTGACGGAAGATCACTTCGCAGAATAAATAAATCCTGGTGTCCCTGTTGATACCGGGAAGCCCTGGGCCAACTTTTGGCGAAAATGAGACGTTGATCGGCACGTAAGAGGTTCCAACTTTCACCATAATGAAATAAGATCACTACCGGGCGTATTTTTTGAGTTATCGAGATTTTCAGGAGCTAAGGAAGCTAAAATGGAGAAAAAAATCACTGGATATACCACCGTTGATATATCCCAATGGCATCGTAAAGAACATTTTGAGGCATTTCAGTCAGTTGCTCAATGTACCTATAACCAGACCGTTCAGCTGGATATTACGGCCTTTTTAAAGACCGTAAAGAAAAATAAGCACAAGTTTTATCCGGCCTTTATTCACATTCTTGCCCGCCTGATGAATGCTCATCCGGAGTTCCGTATGGCAATGAAAGACGGTGAGCTGGTGATATGGGATAGTGTTCACCCTTGTTACACCGTTTTCCATGAGCAAACTGAAACGTTTTCATCGCTCTGGAGTGAATACCACGACGATTTCCGGCAGTTTCTACACATATATTCGCAAGATGTGGCGTGTTACGGTGAAAACCTGGCCTATTTCCCTAAAGGGTTTATTGAGAATATGTTTTTCGTCTCAGCCAATCCCTGGGTGAGTTTCACCAGTTTTGATTTAAACGTGGCCAATATGGACAACTTCTTCGCCCCCGTTTTCACCATGGGCAAATATTATACGCAAGGCGACAAGGTGCTGATGCCGCTGGCGATTCAGGTTCATCATGCCGTTTGTGATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGCAGGGCGGGGCGTAATTTTTTTAAGGCAGTTATTGGTGCCCTTAAACGCCTGGTTGCTACGCCTGAATAAGTGATAATAAGCGGATGAATGGCAGAAATTCGATGATAAGCTGTCAAACATGAGAATTGGTCGACGGCGCGCCAAAGCTTGCATGCCTGCAGCCGCGTAACCTGGCAAAATCGGTTACGGTTGAGTAATAAATGGATGCCCTGCGTAAGCGGGGCACATTTCATTACCTCTTTCTCCGCACCCGACATAGATAATAACTTCGTATAGTATACATTATACGAAGTTATCTAGTAGACTTAATTAAGGATCGATCCGGCGCGCCAATAGTCATGCCCCGCGCCCACCGGAAGGAGCTGACTGGGTTGAAGGCTCTCAAGGGCATCGGTCGAGCTTGACATTGTAGGACTATATTGCTCTAATAAATTTGCGGCCGCTAATACGACTCACTATAGGGAGAGGATCCGCGGAATTC";
////        mn->seq = "GTTAGATGGTTAGAAATGAAATGAAATGTAAAAGTGCTCTGAAAGTATTGGAAGATGCATAAATTGCCCATATTTAAAGGCACATTCTTTTCTATAGAAAATCCCTAAAGAGAAACATTTGTATACCAAAATCTCAAATGAGATCAATGAATTATTCATATAAGATACCACTTTGGGGAGGGTGTTTCTTGGCAACATCTACAATGTATATGTAAAGATGTATTTTGCCTTTTTTCTCACATTCTCAGATTCCTTATTTTTTCATTTCATGATTGACAGCTGTCGTTGTGATGCAGACGGCAACCGTTCAGGATCTGAAGAGCGCAATCAAGAGGTACTTTGAGCTCAAGCAGGATCGAGAGGAAATCAAAAAGAAAATTAGCTGGTGAGAACAGAAACTGAACGATTCAAGTTTTCAAGTAGACTTGCATTAGATTTTAGTTTTATGATTTTTATTGGAAATGTTCTCAGCGCCTTAGATTGTTTTCCATGTTCTTACACTCTAATTTCTACTTTGATTATGCCAAATCAAACTTAAGATTGTAGTTCCTTAAACCTTTACAAGATGTTGATGAATATTGGTATGGAATAGATTTCCTACATGTAGTAATGTGTTTGGAAAAAAGGATAAAATTCTTCATAGTTTATAGTGTATTACAGTCAAGGGGGCTGATTAAGAATGGAGTTTAAGAATAAGGAATGAGATGAAAATTATAGACCTATTTACTGGGATTTAAGGCAGTCCTTTTTAATTTTCTAATTTGTTATGTAAATTTTTCCTTTAATGACTGTATTAGGAGATAAAATGGAAGATATTACCAGAGATGCCAGTCTTCAGGCTTTTGCCTGATTTCAGGCCCTGACAACTTCCCATTTCAAGCTTTTTCAGGCTTTATATTGTACAGCACATATGGGGCTTTTACAAGTTCTGGCTTCCTTTCAGGCCCAGTGAGTAAAACTGACTGGCATCCCTGTATTACAAATTGTGCTTCAAATTAGTTTGTAACTAGCATCAAAGATTGCAATGAATGACTGAGATGGGATATATATTGATGTAGCTGATGTATTTCAAAATACATTTAAATATTATTTCCTTCTGTCTTCCTGAATCTCTCAAGGCGATACGTCTGGAGGAGTTACTATCTTTGTTTTGAGGGTGAGAAGCTGACTGATGATTACAAGAAAATAAGAGAGTAAGTATTTTGATGCCTTTTTACTGAACAAAATGATCTATTTTTTTTATCCCCCCATATCTTGATATTATGATTTAGGACACATAATGCAAAGTTCCCATCGAAATCCTGGGTAATGAACGTTGTGAAGGATCAGAAGCTTTGTTGTGTCTTCTATGTTCAGTGTTCAAGATTGAATATACTGGATGTGATATTTCAAACAGTGAAAGATTATCACTATTTGGATGGTAGTTTGAGGGAAGGAGAAAAAGTATGTAATGTATTGTACATATACATGAGAAAGACATTTATTGAAATGAAACATTATTTTTCAATATTTCTCTTTTTTGTCTTCAGCTTTGGAATCAGGAATCGTGCCGAGGTGACCTTTGCCAAGAGACTGCATTCCAAAGGCGGAGGCTACTGATGATGATGATGAAATGACCAATCACGTCAAGATAAAGGGTTGCTATGGCAGCTGGATACGAGATGCCTATACCTGGGATTTGTCTTCTTGAGATGGAAAAAGAAGAGCAACCATGGCTCAAGGAGTTGTTTCTGTGAAGGAGGGGTGCTGGTTTGAGTGTCGACCTCCTGCCTGTATGAGCTGATGCTCTGGAACAGCGATGTGTGATAACATAGAGGGAGCTATTGCTTCGATGTCTTCATGAGCATGCAGGAAGTCGGAGGAATGTCTCAGATGAAGGACGAATGTTAATGGCATGAACCTAGAGGGCGCAATTGCTTCAGTGTTTTCATGAGCATGCAGGAAGCCAGAGGAATGTTCATGACATCAATGTAGTCTATGATTTGGTTCTTGCTTCAGATTGGTGTGACATAGTTATCCTGATGATTGGCAGCTGTGAAGGTTGAAGACACTTGGAAGGTTCCAATGCTTTCATTTCTTCACCTACATTGTAGATGTAGGACCAAGGAGTGTGGAGGCTGCAGGAAGTCTCTGGTGAAAGATACCCATCAGTATCGCCAATGTGAATTATGATCTTGGGGGTGTTTCACAAAGATCCTAAGTCAAACTTAACTCTAAGTTCGACTTAAAAACTATGGAGAGCCATCTGTAGCCTCAAAATAAAGTATATATATTTAAGTATTTTAGACTTTTTATAAATTCCATATCAAGAACATGAAATTTCTGTGATAATAGAAAATATGAGCAATTTGTATCAACACTTTATGCCTTTTAAAATAACTTTTACAAAATATTTTCATAGATGCTAACAGCTTTCCATAATTGTTAAGTCCAACTTAGACATAAGTTCAACTGGTGATCTTTGTGAAACACCCCCCAGCATGCCAGAAGTCGGGAGAAGTATGACTGGAAGAAGATTTCTATATCAATATGGTTTATGATCTGGGGGTTGTTTCACAAAACCAATGCGACGACAGTCATATATTTGTCTGCAACCTATGGAGATGTCAGGAAGCCTACAACCTGCTTCACAGAGCCATGACAAATCAACTGCTTCATGAAATCAATACTATTGGCTAAGAAACAGTCTCCATCACGATACATCACGAATGTGTGGGTCGTTGGCCTCATCAAGCGTTCATGCAACTTCTTTCTTCTTATTTTCTCTCTCTCTTCTTGTCTTCTTTACATTATGTTATCTTTATCTTTCCTCTGTCCTCTGTTGTTGTGAATGTATGTACTATGTTTTTATATAAAATAAGTTTGTATCTTGTTTATGCTTAGTATTTATTTGGGGATCATTCTCATCAAGCCAATTTTGGCTTTCATGATTCCCCACTTCCTTTTATCTTATGTTGATAATGAAGAAACATGTATACTTTTGATGTTTTTTAAAAGAGAGTGAAATGAATAAATGAAATGAAATGAAATCCAACATGCTCATCATACAGGCAGTATGTGACACTCAATAAATCAGAATTAACTTTAAATTGCGAGGATTTATTGTATTGTATTGAATTTATTATCAGGATACTTCTTGGGCTAAATCGTAAAGGCTTCTTCGCATAAGTTGCAGTAACCCACTGGCTTACCGCAGGAGCTTTTTAAGGGTACTTTCTTTTAAAAATACTTACTGCACTTCACATATGTATCTTGTAACACAAAAATACATTTTTTTGGTACGGGCGTGCAGCTCTAGGTTTTTCATCAACACTAGCCATTTATTTAATCACAAAAAAACTGCAGTTTACTCTGGTATATTTGAATGAATTTAAGTTTATGTTTTTATTGTGTACATGTATGTTTCATACAGCAGGTTACTAGATTCATTTATTTGAGACAAAATAAATTTCTTTTATTCGTATGTTGTGAAACAGAAGAATTTACATTCAAATATTGTCTAATTTGAAATCATTTTTTAACAACTGCATTCTAAACCAATCACCATAATCCTCATGACCAGACGAGTAATTTTATATGCATAATGACACGATTGGAACTATTAACATTTTTATCCCTGTGTGACCTTTGTTGACCCCTCACGAAAATGCCCGGGATTTTCAAGGGGGCAGCCCAAGTCAAAATTGTTCGGTATTCATCAGGGAACACAAAAACACTTTCAAACCTTTGTTCCCAGAGATTGGTCACGGGACCTACAGGCTGGTCTCCTACTAATAGACGGACCGTCTGGTTCTGCATGTTATATAACTTATATAACTTATCAACAATAGAATAGGAATTTTACAGTTTTGAGACGTATTTGAATACAAGCTACATTACAAAGAGTTTACAAATGACTAGAACCGTTCTGAAGCATCCTTATTATTGTATCCAGTCTAACGATCATGAAGTAGATATGTATTATGTATATAAAGAAGACTACTTGTAAAACATTACATAAATATAATCTAAATATTTCAATGTGAGCTCAAACAATAAATCATTAGTACCTCGGTCACACTACCCAACCGACTCAAGACCCTCTGCATGCACTCACCAAACGATTTTTGGCCAAAGGGGGTAATTTTAATTGTTTGTTTGGTCTGTAGTTGGTCAGTGGTATAGCTGACGTATAATAATATACAAGTGTAGATATTATTAAACATGATTTAAACAGTTCAGTGTTCTTATGATCTTTTTACACAAACTTAACTTTCCCTAAACTCTACTGTGCTACCAAATCAAGTTAATTCATCCACATCATTGGATAATTTTTTGAGAATTAAAGATGAAATACATGTATGTCCATACATTTTTTTACTAGTTTTATACTATCATATTTATAAAATATATTTCTCAATTACTTAAATGTAAACAACAACCCTAAGCAAACTCAGATCAGGCCAGGGTATTATAATATGTTATCAACCTTAATTAATAGTTGAACTCATGGAACAATCGTTGTTTCATTTCTCCAACCTGTGCATAAACAAAAGAAAAAGGAAATTATGTTTTTAACATGTCGGTTTCGATATTCACAGACTTGAAATAATTTCTCCATTTTGTAAGTGTTCATATACCTTCATGGTTTAGCCATGTTTTTTTTTTATTCAGAAACAAACACCCCCCCCCCCCTTTAAAAATAAATAATAAATTTAAAAATTAGAAAAAACTATGATTTTTGTCAAGAACAAAAAAAAAGGATATGATTATGGTGGGATAATGAACTTGAGTGAGAAATTATCAACAGAGAAAAATGGGGTTGAGAGAGAGAGAGAGAGAGAGAGGGGGGGGGGCAGAGGGGAGAGAGATAGGGAGAAAAGAGATGGATGGGACTAAAGAAAAGGAGAGATTTGATGAAAAGATGAAAGCTAAAGACCAAAAAGAAGAGAGTGTGATGGAGAAAAACAAAGATATAGCCTAATCATGTCAGTGATTGCTTACCTCTTCTGGTACACGGTGTCCTCTCTTTCGAGAAGAATCGGGCCTTCTTGAGCCATTCCTACGGATTGGTCCTCCACCCTGAAGGTACTGGCCGGGTAGCTATATAGAGAGGGAGAGGGGGAAGGGAAAACAGCCAGGCGGATAGATACAGAGAAAGAAGAACAGATTGAAAAAGCAGAGCGAGAGAGAGAGAGCGTGAAAAGAAAAGCAAAGACAATTGTGTGGTGTGGGAGAGAGTAAGAGAGATGCATACACATGAATGATTTGTGCAATGTGTATTTGATTTGAGAGAAAATAATAGAAAATATCTTTTTCTTTTTTAATACACAAAAAGTCATGCAAATACATGTACGGACACAGAATAAAGAGATAAAGTGAGAGATAAAATATTATCCATTTCATACAATTATTTGTGGATACTCGAGGATAAATTGTCTCTTACACCAGTTTATATATATATATATATATATATATATATATTTGTCTCTTATACACATCTGGTCAACCTCCAGATCGCGCTGGAATATGGCCAGGCTATGACGGTGATCGTCAGAAGAGCTGAAGGCGATGTCATGCGTAAGTAAACTGTTTTTTAATTATTTGGGAGTGCCTATATTACATCACCAGGGAGGCGTTACACAAAGATCCTAAGTCGAACTTAACTCTAATTGGACTTAAAAGTTATGGAGAGCTATCTGTATTCAAAAGCGAAAATTTGAATAT";
//        revComp( mn->seq );
//        bwt_.mapSequence( mn->seq, mn->ids, mn->coords );
//        mn->recoil();
//        Node* base = new Node( mn, 0, mn->ids.size()-1, 2 );
//        base->dontExtend_ = true;
//        delete mn;
//        nodes_.push_back( base );
//        mn = new MapNode();
////        mn->seq = "ATTACCTCTTTCTCCGCACCCGACATAGATAATAACTTCGTATAGTATACATTATACGAAGTTATCTAGTAGACTTAATTAAGGATCGATCCGGCGCGCCAATAGTCATGCCCCGCGCCCACCGGAAGGAGCTGACTGGGTTGAAGGCTCTCAAGGGCATCGGTCGAGCTTGACATTGTAGGACTATATTGCTCTAATAAATTTGCGGCCGCTAATACGACTCACTATAGGGAGAGGATCCGCGGAATTCTACAAGAGAGGAGAGATCCTACGGAAGATCAGTGAGCGCCGTGCCGCACCTCAGATCACTGCCATGCCCAGTATCCAACCAGATGACCTGCCTGAGAAAGGGGAACTAGGTAAGCATGGAAATACCCACCATGCTCTGTCTTTATATTTTCAATCAGCCAAAATATGAATACAAAATGTCAATATTACACATTCTGCAATATGAAACCTGGTGCACTTGCTTACTACAGCACCCTAGGTTTTCATACGGTACCGTCCTAAAACCAGTTGAAATATATTTGTTGTTACTACACATTTAGACTTGAGACATGTTTCTGCCTAGCCTATAATAGATGTCACAGCAAGGACATTTCTTGTGATATTTTTTCATATTCAAGGACTGAGGGAGTTTACCTATTTCATTTTGTATTACTGCAAAATCAAAAACTTTGAAGGCATGTAACTATTAGCGAATTACACAAATTCAGTGTTGTGTGAATTGTGATTCAATTATTTTTTCATGGCAAAGATCTTAGACATTTCTTGCCCACAAAAGTTTCTGATTTAACAGTATATAGTATTCAGATTTCTAGGTTTGAGAGAGTTTATCCACTGCATGCAGTATTTCCAAATTATTTCAATATAAATAATTGTTAGGCTCTAGGTTTGAGGGAGATTATTTCAATTTGTATTGTTTGGTACGGTATTTGTTTTTTTCCAAGTTGAGGGAGTTTATCCATTGCTTGCTGTCCTTTAAACATTGTTTTAATTTGTGTTGTTTTCTACCCTAGGTTTGAGAGAGTTTATCCACTGGATGCAGTACTTAGAGATAGCCATTGATGTAGATGTGATTGACAGGCTGTTCAATGCTTTATCAGCAGGATCCAACTGTGAGTAGTAATGCAAATATATATAGATATTAGCCTATTCTCTTCAAATGTAGTTTTGCACTATTGTTGAGAAGTTACTACGAGGTATATATGCTCCATTATATATTTCACAAAGCTGTCACCTTTGTTTTGAAGAAGTGAGAACCAATCTAATAGTGTCACATGCAAAATTTACTTTTATTTAACTTAGTCTCTTGCTACTTATTTTGTCATAACCAGACACTTTCACGGCAAGAAACCTACATTGCAAATAGCAAATTCTAGTGTTCAGATTCAGTCAATTTTGTGGCATGAAACTTTTGCAAATCTCATTATTTCTTCATGAAGTTGATGTTAAAGGTGTTCAGTACCATCATTTTATCAATTCTCTTGTTTCTGTAGTCCTTGATGCCGAGACTTTCATCTATTTCTATGAGGCATGGAAGGAAGGGTTAGAGGAGACAAAGTCCATACCACTCAGTATCTACAAGGAATACTTGGAGAGCACTGAGAGTGTTCTCAAGATGTCTTCCTTAATCAGAACAGATCATGGCATGGGAAGATTGGTCCTGACTGAGAAGAGGTACAGTCATTTGGAGATAATTTTAGTTATGGCTGTGACAAGGAAAATTGCAATGTTTGAATATCATGTCATATCATTCTACTCATATCATATATAACAAGTGCCTGATAAGTAAATGGTTTGCCCTTAGCTTCTGATGATTTTTAGTGTAATTTTCTGTATGGAATACCACATTAAAATTACATGTACTGCTTGCTTAAATCTGTCTCTGCCTAATTTCATAGCACCGCAGCCACTTTTCTGTGCACACACACCTTTCCCTTTATTTTACACAGTTTAATGGCACTGTTGAATTGACATTCCTCGCCAAATCTTTTAAGTGAAAGTGGAAAATGTCAATTTGTTTAGAAAGACAGGTTTTTGTCTGTGATTATATACAATCATATCAGAATATTCCTGCAATTTAAACTAAGATATGTGGATGTGAAGGATGTTTTTTTACCCGTCTTTCCCTTTGTCACTCCTTAGGTATAACTGCATGTCTTTGCATTCTCTGAAATGCAAATTTCCTGCAAGATGTATCCTGTATAAATTTCATATTTTCCTTATTTACTAGGCTCTTCTACATGGGCGAAGGCAGCAATATCTACAAGGAGATTGTGAGAGCCAGAGACATTGAGAAGTTAGAGAAGTTTGAATACTTTAACATCCTCCTCAGCTGCCAAGCCCTTAGGATATACAGCACCAGTAAGTAACACTAACCTCAGCATGTTATGTAACTTGACCGCAGAGCAATGGTCCCCGGAGAATGTGTCATAGATCCAGTCGGTGTTTTGTCACCTCATTCTGTACGTCGTTTACATGGTGGATAGACAATGTACTCTTTAGCTTATACAATGCAGGTTGTATCTGATTACCAGCTAATTTTTGAATTCCCACATGCGTTCAGTTTTAGATGTATGATAAACAGTTTGGCAAAATCTTTCATTTGTCATATAACATGCAAGAGGTAGCTAGTTTTATTACAAATATACAATTCATTAATTAAGTAGTAACAAAATGTACAAGCCATTTTTACTTTTCCTTTATGCTTTTATCATATTTTTGCACAGTATTGTTTGCCTGATTGTGTTTTGCAACATGGAAATGTGTATAATGTACACAAACATTGTATACAATGTTACGACATTACATCAAGCACTCATTTGATCACTCTTGTGTATATTTACTCTATATTTCTTTACCAGAGCTCAATACGAGTCCATACGTTGCCAATTTGAGGGCAGAGCGTTAGAGCTGGTTCACACTGATCACTGAGCTGTGACATCAAGCACTCATTTGACCACTCTTGTGTATATTAACTCCATATTTCCTTCCCAGAGCCCAATACCAGTCCATACGTTGCCAATTTGAAGGCAGAGCGTAACAGTTGGTTCACACTGATCACTGAGCTGTGGGCCGGGCGTGTCATTGCCGACGCTCAGAAGGATCCTCAGGTGGTGCAGCAAGCCGCAAGGAATGTGAAGCTGATAGACTCGGTCATACGGAGTGCAGAGAATGAGGATGCAACACACGCAAAGCATCTAGATAGTGCTGTATGGGTAGGTATTAGACAAGGGTCCTGTATTCAATTCTTATTAACAATTAAGGGAGGTTCGATCGCTCAGTGATAGAGCTGTGATCTCGTAAACAGGAGGTCCTGGGTTCAAAACCAAATCAATGCACCAGTGCCCTTTGGTAAGGCATTAATCCCCATTATTAAGTCCCTCAGAGACTCAAAGCCGTCATGCTTTCTTAGCAGTCGGGTAAAACAATCAAACAAACTAAATTTATTACCATTTTATTATTAAGATATTGAACATATTTACAAGGTATATTTATAAATCTCTCTATGCTCAGCAAAAGTTGCTTTGTGCCATTTTATGTCCCTCAACTTTTGTTACATTAAAACACCACAAATTGTTTGAGTACCACTTTCAAATGTGGCACATGAATGGCCTGAAAGAAACTGTCTATAATTTGTTTCAGTAACATCAAATTTTGTAATGATGAACTCTGATGTAACGTTTTTGTGTGTCTAACAAAGTGGAAGAGCATGACCTTTCCTTCAAGCAGTCTTCATTTCAAAATCTGGTAAAATTTCTTATTTCAGCATCTCTGTCACTTGTGACTGAGAGAGGAAGGCGTGGGGCACATTTAAAAAAATAATAATGGCTTCTTATATAGCGCACATGTCCGCCAGGTCGTGACGCTCCTGGCACTGCCCGATTATTACCCTGGCTTTAGCAAGGCGGCCATTATGGCGCTATAGCAATACACAATGAAGTATGGTAGAGACTATGATCATCCAATTACTGACTTTTTCCAACAAAATTTTCTTTTTCTTTCCTTTTGTATTTCAGCATCTGTGTCACTTCTCACGGTTGAGAGAGGAGGGCATGGGACGGGTTCCACCAGAGACAAGCTCAGCTTTGGTTCACAAGTTTAACCCATCCAGTAACGAGGCGCAGAGGACTACAGTGGAGGCTATGGTTTACACACCGGGCAGTAAGTCCACTCAGCTTGATGTCTGCATGGTTCTTTGGGTTTCTTTTATATCTTTCTTTTTTAATCTAAAAACTCATAGTAATACCATTTTCTAGAAAAGCTATGCATGATTAGCGATTTTCCTTTGTTTGTGACCACACTGTTTTCAAAATGCACATAATTCAACAGGTTAGATATTTCTTTTTCACCACACAATGCACGGTTTATCATCCTCTAAGGTCCTTGAAGTGAATAGTGAAATAGTTTTGTTTGTGTATTTTAGTTTATGAGGAATTTCAGCTAGATATTTCTTTGTGCTGGGTCAGGTAAGTATGCTGCCGTGTGGAATTTATTCTAACCGTGTAATAACTCAGAAAATTTTTACCAAAAATGTAATTACACTACAATAATAATTCATGAATTCTCCGGATCAGTGCTATTAGAATAGATAATACCTTTAATCAAGAGAGGCCTTGTATGTCAGCATCAATATCCCCTGATGCATGCTGCCCCCACCGACGGCACACACATCAGAACTTACCCCCTACTATACCCCTCTGTACATCCCTCAGCCTCTGTCTTTTCAGACCGGAGTATGAGCGAGGAGGAGGCCACGCCCAAACTCTGGTGCGCCATGGGGTCAGGCAAGGTCAAGGTCTTCGACGGGAGCAACTTCGTGCTGGAGGCGGAGTTTGCTGATGCCAAGGACAGAGTGGTAATTTTTTTGTTCTTCTCGCCTTTTTTTAAAGGTCTTTGTTGCTGAAGTTGAATACACCTATACACTGTGGAGGTTGATGTCATGATTAAATGTTCATTGTTAGTAGTTCAAGTTCAAGTTCAAGTTTCATTTCCACATTTTTAATTTTTTTTTAAAGTCAATTTATATGCATATATACAATAACAAAAGTTAAAATGAAACATCTACAAACATTCTTTTTGTTCAAAAGTCACATCAACTATTTTGATATCATTTATAACCTTCAAAGAGGGAAGAACGTAATGAGCAAGCTTAGTTGGAATTCAAATAATTTTACTTTCAAATGAAAATTGGGTAAAGTTCCTTTTATCATCAATCTAGCCTGTACTATTAAAATGAATTATTAAAACTCTTCATTTTATAAATTCTGGTATTTGTAATGGCCGAATAGCTGCATTTTATAATCAAACCCTTCATTTGTTTTACTAGTAATTATGTATCATTCATATATTTCTCTTTATGTGTAACAGTGCTGCCTGCTGAGTGTGAGAGGAGAGCAGGTATGGGCGGGGTCGTTCGATACAACCATCTACATCATCGACATACCCAGCTGCCTGTCAAACAAACAGCTTGTAGAACATAACGACATTGTCTCTGACATGACCATATCAGAAGATGGCAAGTAAGTATGATTTAACAACATACATCTTCAAAATTTATTAGTACCCTAACCCATTGCCTATGGGAAGTGATGGTCTCTGTATTAAGCCTATGGGCAATAGAGAATTCAGAAGCCAACACGTTATATAGATACTGTAGAATGTTTAATATGTTGTCATTCAAAGGATAGTAGATTCACTGTATCATGAATATTTTAGAAATGAAATTGCTCCTGAACTAATTTACAAATGTGAAACGACATATGTACATGTTGATTAATCTGAAGAAAAAAAATATAGAATCTGTGATTTACAGTGTATTACAAAAGAGAAAACATATTTGTTTGGTAGCACAGATAAATCTGACCAATGTTGTTCTCTAATAACCTTTAAAGCTTATTACAGAAGTGTGATTGATTTGGGTTTCGTAGATGCATAACGAAAATCACATTACCACCAAACCCTTGACCTTGTCCAATGCGATATTACATGTGATGTTGAGCCTTTAATTATTTTCCTTTTGCTGGTATTGTATACACATTGTTTCAATCAACATGTTGCATTTCTGTTTACAGGATTGCCTTCACATGCAGCTTGAATGGGCAGATCTTTGGTTGGGACACGCAGAGCTTGTCACAGAAACATCAGATACAGCTGAAGAATACCAAGACGCTCGTCTCAATGAGATGGTACAATGACAAGTTATGGTGCTGTAAGTTTCCTTGCAATCATAAACATCATGATTATCGTCATCATTTAGTACACCTGTCTCATATATTGTGATATTATTCATGAAGTCTTATGGCATGTTGATATCTTTTAATTTGCTCATTTTAGACATGTACTTCTTTTGCAGTAAATCATTATTATTTCAAGTTTCTTTATTATTATTGTCATTTGTAAAAGTGTGTAGAAACAGTCCTGTATTTGTCTATTATTTCTAGTATTGATATGATCATGATTATTATTGATGTATCTTTATCATTGTTATCGGTCTCAATGATTATCAGCAAACATTTGTTTACTTTTAAAAAAGACAAGAGAATGAGATTTTACTTTTTATATGATTATCATCTCATCATGCTGTGAACACATTTGCAATATGCAAAATGTTTTTCATTCCCTCGATCGTCAATGTCCCTGACTCTATGATAGCTAAACTCCACCTGGTATTGTATTTTGCAGGTACCAAGACTGATATCAAAGTCATAGGCCTGGATGGTAGTGAGCTGAGCAGCCTCCAGCATCTTGACAAAGATGGGGGTCCATCCCTCATTGAGTCGTTCCTTCTTCTAGGGAACAGGGTAAGTTATTTAAATATCAATATCATGAGTGCGACTATGGACCCAAGGCTGAGATGGAGGCAGGAGATGATGCAAATCGTGAGGATAAAGATTTTGGATGGTGGTGGTGATGAAAATGATGACGATGAGGAAGATGATGATGATGAGGAAGATGATGGTGGTGATGATAAAATTGAGGGGAATAATGCTGATGGTAGTGATGAAGATGATGATGATGAGGAGGAGGAGGATGTTTATTATGATAGGGAAGATGAAAATGATGATGGTGATGATGATGGTGGTGATGAGGAGGATGGGGATTTTTATATAATGAGGATGGGGAGGAGGAGGAGAGGGAGGATAATGAGAATGGTGCTTCATCTTATCATGATTGTGGGGATGACGGTGTTGAATAATGATGATCATTTTAAGTTTCATTTAGATTTTAATAACGTTATATAGAACGTTACTTACAGCTAATCTTTGTTCTTTGACAGATCTGGACAGGTTGTGGTAGGAGAGGAGAGGTAGCCTGCTGGAATGTGAACACATTCAAACAAGAGAAGCTTCTTACCATATCATGCCGGGGTATCAGCAAGCTTGTTGCCGTAGGAACCAGGGTCAGTATATTATACATTCTCCACTTAGAGCCCCTGCCTCCAAATGCTTTGTTATTATGAAAAGATTTTCTTCATAAATCATATGGAATGCAATAGCAGTCTTTAGTGGTTTGTTCTGTAACATATTCAATACATTTTTATCTCATAACATATATCAGGCTAGGTCTTAGATAAAATGTAAAGAAGCCATTTCAGAGGTCTACATTGATGTTGCCCTGTATCACTGATACTAATACCCAAAGCTGTCCCCCCCTCCCCCCCCCCCCTAGTGTAGCTATGATGTCATCCTCCTCTTTTCCCTTCCATCCTCTTTGGAGATCAGCTCATTATTCCAAAAAGAAAAGCCAATCCCGGCAGCCAGGACTCTCCGTAGTCTAGCTCAATGTGCGGTCAGCAATTTTATAATTCTTATTTTGATTACTTATTGTTTATATCGGTTTAGATATGGGCTGGTAGCAAACAAGGCAAGATTCACCTCTTCGATGCCAACACCTGTGAGTACGAGAAGGAGCTCGAGGCTCATGAAGATGCTATCAGGTCCATGTGCCAGGCAGAGTTACGCTACGTCATCACCGGGGCAGGCAGCAAGGATGGAAAGGTCGCCCTCTGGAGGGCAAACTTTGTGGCAAACTAAAGGAGTGTAGCTTGGCAAATACTGGTTAACGTCATCTACGACTTCTTTTGGACAACTGCACGTAAACGAATAGACTTCATTAACGGTAGTCAACGACTTAATCTCTACATGGTCAGTGCGCGCTAGTGTTTGCGAAACAGACCGAATTTCTCCTCACGCTTCAGTTTCGCGTATGTTGCATGCACATGATTACTGCGGTTCCGTTGTCCAGAGTAGTATTATACCTGATTGGCCCCAAGTCAGGCCACGGGAAAGACACCAAATTTCATAACTTCAACTTGAGACAGAACATGAAACCTAATCTTAAACCTGAAATCCAACTGTAATCCTAGTCCTAACCCTAACCATAAAAACTTTGATGAAATAAGTCCTGGAGGAGCCCATTAATGCCTCTCGTGTAGATGAACAGCTTGAGGATGGAAAATATCTGTCTTGAATATATTTAGAACAAGTTGTAACATCAAACTTTTCACCCAGATCATTTGAGAACAAGTTGCTGAAGCCTATGAGCTTAAGCTGCCCAATTTAACGAAGAAATTGAGCTTCTTGAAATTTAATAATGTTATGAACTGAAATATAGAGTTACAAGTAGGGATTGTTTACTTTTTTTGGTTGAAGTTTTCTACAGGGAAATGAAAGCTGCACAAAGTTGTCAATATTTATTGGTTATTCTTGTTTAAAATCCACCTTGGAAGGAACTACTCTTATCCTTCAATGTGCAATAAAGATCCTCATGATTTAAACGTAAAATCACAGTAAAATCACAATCAGTGCTGTATTATAAAGATATATATATACATCTGTCCATTTATTTAATCTTTGACATCTTTTGACATGTTTTTTTTAACCTGTGGCTAGTTTTATTCTTGTTCTTCTGATGTAAAGGCTTATTTAGAATTGTACATATACATGTGCCTGCACATCTTCTTTTTAGAAAGTCAGTTGGTATTTTCATTTCCTTTTTCCTCTTAAAACGTATATGAAGAAGAAAAAAAATCATTTAAATACCTTTTTCGTAAAAGATTAGAAAACAATCTTATCAAGAATATTTAATATACGGTACTTTCATCTAGAGTAGAGGCAATAGTAAATGCAGAACCAATACGTTTCTCACAATGTGTTTTTGACATTGAAATGAAGATGCAATCAATTTAGTGTCCTTTGTTAATGAGCCAATTGTCAACAGTGCCCAAATTGCAGACTGTTGAATGTTCTCAGAGGAACATAGATTCTTTTTATTTTTAATGCATTGTAACTTCAAATTTTATGCCGATATAATTTTTTCAAAAATGATCTGTGTTCTGGACCTGTTAAGAGTTCATAGTGAAATGTAATTTATATTAAGTAATGTATTATGTGCTATGTCTATGTGTAAATATGTCAGTGTATTATGTAAATACACTTGTCGTACTTAATGTTTTAATAGTCCAGCTATATATGTATATATATATCAACCACAGCTGCTTTTCATATATCGTATACTTCCATTTTTGATTTGCTTCCACAGATATGTCTTTCTTGAAAATATGCCAAATTACATGTACAATTATTTCCATAATAAAGATGCCTAGTAGTCAGGTTTTCGTTGGGAATTTGACACCATTCTGATGTTTACCCAATAACCGGTGTTGTCATTCACCAAGCGTTGACTTCACTTTTACGAGAGACATGTCTTGTGTAATGCGCTGTGCTCTAGCCTACCACCATATGGAATTCTCATCATGTACAGTGCAAATATTTGTAATAGCTGATGGCACTTGACAAACTAGACGTATGCATTAAATTTGGTTACATTTCAAAAATTACCAATTGGAATTGCCCTATTGGTACATAATCTGGTGATGCATTTGTTATTCAGGATAGTGCTTACATATAGTTTCTGGTGTCTAGTGTCTATGTATAAACTGTTCTTGTACAGTTATAGCATCAATATGTAGTGTTTACAACCTGTTTAATTGCAGGTGCTATCGTAGAATGGTGTATACATGACCTGAGTGTCGACAATTCTGTACTGTTTACTTTTTTTTGCTTCCATTCATTGAAACGTTGATTATAACCTTCCATGACCTGGTATGAATACACGTATGAATAAAGAGTGGCCTACTTGGGATACTGAATACACACTGTCTTCTATACAAACATTTCAAGGCTCCGACTCGCAGGCCCAGGTTTAAAACCAGCTGGCTGATCTTACTTTCACATTCTTACTGTCCGAACCTCGCGAGTCTTCTCATTCTATTGTATCCCAAGTTGGCCACTCAGTATTCAAATTTGCATTTCTTTTCTTATTTAAAATCTTTATATATCAATTTGAAATCCAACTACCACATCAGTCTGTGAATATTGATGTGCATCCCCTCTTTCTTCTCTTTTCTTTTTTCTTCTTTTTTCTGTCTTTTCCCCCTTTTACCAATATTTTGTTCCTCCATATTTGTATTATTATGTATCCTTTGAGGGTAGTCATACTATCAAGCTTTGCTTATTTTGTCAACCCTCCTGTTTGATGTTCTCGATTAGTTACACCCAGCCTTACTTTATTCTACATATATTCTTTTGTCATGTGTAAATTTACAATGTTTGTTTGTTATATATTTGTAAAAATGCAAAAGAAAAAGAAAATCAAACAGAAATCAATCAATCAATCAATCAATCAATCAATTCGTCCTTGTTTTCACGAAGTGACGGTAGTGAATGGAAATCGCATACACTTTGCGTGCAGATTTTCCTCACTACCGTCACTTCGTGAAAACAGGAACGAATTGCATATAATACATGTATTATTATTATGTCACTTTCATTCCATTGAGTTTTATTTTTAAAAGTTATTCTTGTGATACTTTTCAAGTTTCTATCAGTGTGTCAATTTCATACTTACTATGCGACCTAAACAATCATGATGTGCAATATTTGTTATGTTATCGTAGTCTTTTTCATGCTTTTCCCCATTTCCCAAAGTCTTCTGACTAAGAGGAAGAAATTGGGTTTCTTGTTGTTCAAATTTTGATTCCATCACATTTTTTTAAAATACACAGTGCAAGTTATTGGTCAAGGCCATTCATGTAACTTAAGTGAAAACATATTTTTGATGGGGGAGTCTAGTATATCAAGTACATGTGTATGCAACTTCTCAAAAATATAAACAAAATAAACCTGTGTGGTATTTTGTTCTACCCTAATTTTTATTGAGGTAATTGATGGAATTTTCTTTCATTCCTTAGACTTATATCAACACTGCATTACCTGTTCATATATTTCCATTTGGCATGCACTAAGAACATTCTCAGGATAATATTAATTATACATTTAAATGCAATCATGATGTGAACTTGTGGTTAGGAGAGAATTGAAAGTACATTGTACTTTAATATGTAGTTTCAAGCATTTGCGAGCCTTCACCACTCTAAAGGAATTTACAACAACAAAAAACTAGAGGTAGAACCTATGTTAGTATCGAAACTTTTACTTTTTTGTATAAAGTAACACAAATTAGTTAATTCTTTAGATTTAGTAGAGCAGGGCCCATTGCATAAAACTTACATTGATGGTAACTTTGCCATTAATTGTAATTTTCTAGGAAACCTTGATTTGATTGGCTGTTGAGTGATGTTGCCATGGTGACTAAGTTACATCAACGGTAAGTTTTATGCAATGGGGCCCTGGTGTCAAAACAAGGTTGTACATATCAAGTCTAGCAGGGTCTATACTTCAGTGTGTATTGTGTAGCACATATTGTCAGCAGGAATGCAGTTCAGTCATTTAGACTAGACTCTAATACCAAATACAATGCTAATCAGCTCTTAAAGGTAAAGACCAGTTATGGTCATGCGATTGCATTGTGAAAGAAACGTTGTGGAATCGATCAGAAAAGGTTGATCATCTAGCATCTAGCATCTGCAAGGTTCTGGAATCATCAAAATGATAAGAAATATGGCCCATTCTTCGTAAAGTTACAAGTATTATACACAGAAAGATGTATATGGGACTACCCTATATCTCTCACAGACCATGAGAGAAACAAATTTACCAAATATATGAATTAGTAGGCGGTTTCCCGGCAAATTTTTTATCATCATGTATTAGTCAAAACATATTCTTTCACATGATTGGTAATGAGAGCTTCCATTTCATTTTTTGAAATTGGGGCGTTTCCAATACTAGTCTTTACCTTTAACATATGGCAAAATATGTTCTCATAGTATATACTAAATATTCTCTAGGTCTCTCCAATGTTCCCCTTTGCGAGAAATATTTCATTGAAAGATTATCTCATATACCTCTAGCCTCTTGCACCATGGTTTGCACATGACGGACACAAAGTTTAACCCTCTAAGTAAAAGCCAAGAAGACTCTAATCATATTATTTAAAAAATTAAAAGATTGTTCAATAACAATGCACTTACTTGCAGTTCCAATTGGTCACATTGTCACTTTGCTAAATCTATGCATGTTTGTATATCAATAATTTATGCAATATATTTTCATATCTTCCAGCATAATGCAACCAAATATATTTTTTGCAAATTAAATTTACGTTATTATTAATTTGTTGATTATATGTTTATACTGTATATTGAACATTATTATTTTCATGTACAGGTATGTCTACAATAGATGAAACAAAAAAGCTATTACGAGACTATGTTGCAATAAAGTCTACATGTGAATATAGCCTAAAAATTATATTTGACTACAGGTATGTTGTCTTCTTTTCAAATATGTAAAATTGTCTTTAAAGCAAGCATTATTAATGATATACTGTAACACTTGATTATTATTTTTTTAATCACATCAACAATTCAGACAGAAAATTCTGAAACTGCATTTATTATCAAGCTGGACTTGGTACATAATGTATATTGTAAATCAAATCTTGTAGAATTTATCTACACTAAATCTTGAGACATATTAGTTGTAAAACAATGCTCAGAAATTCCTACTAGATTTTAACCCGTGGACAAATGAAGGCAAATAAATATGTATTTTTCTGATTGTATCATCAATACATTTGGAACATAACAATACCATATATTCCATTTCACAACATCTTTTGGAATGTGCTTCCTACTGTCAGTAGCTTGGTATGTGTGATACTCTTTCTCAATTAACCCGTCTTTAACAGCATGTGCTATCTATGTTTAGATTTCATGGTGTAAAACAACTTACATTGCATTATATATTCATGTTTTAGGTGTTGTATAATGATAACTAAATGCACAACAAAAGAAATTACTTTCACACTTGATATTTCCATTTATGTTTGTCTTTTATCTTTAAACAAAAATTTATTTCACATTATTTATAAACATGTTTTTTACACAAAATATACCAACTGTACATACATGTAGTACATTCTTGGCACATAATTCATATACAAACCTGTTTAGTGTCTAAGTAGTCTATAATTTACAAATAAAAACTGTTTTTCATATATATGTATACATGTACATGGGACACTCATACCCAAGAGTTAATAGGCAATCAATTACATCAAAAAAATTCTTTTACCAGCTGACCTCACCACCAAATTACATATATTGTGTATGCCCAATGACAATAACAATTTATGGAAACTTTAATCAATTATGAGCGTATTACGTGCATCTTAGAGAACTTTTACTCCTGTAACAGCTTGTCATAGCAATAAATAGAGGTATAACTACCGGTAACATTTCTGGACTATACTATATTATGATGCAAATTACAAATCAAAATGATCAGAGAATACAGAAGAAACGTCTCTAAGTGTCAACCAAAACTTCATAACAAATCACATTGGAGAGCGGGCTCGATCACACCTTTTGCACCAAGTATTTTTTACTTGGTCATATTAAAAAGGTATGTATTAAAGCTGTACATAACGAAGTAGATATGATGACAACATTATATACAATTACAATATTATTTCATCCACATTCATTCAGCTTCTTTGAAAGGTCCCAGGACATGTATACAAATAATGACTGTTTTGAGGGGCACAATGAAAACTACAGCCCAAGGACATCTCAAAACGATACTGTGCCCCGAGACAAAGTTGTATGCGTTAGTAGAACTATATCTGGGCAAACCTAAAAGCTAGATTCCCTTTAAGCAATAGCCTAAAAGCAATGTATTAATTACCATTTGCCTGTCAGGTGCCAACATGAAAGATGAAAGTTTACTTTGAATAATGTTATATATCCAAAAGTATTGTTTACAAACAATATACATGTATGCTTATGAAAAGAAATAGATATAATTTGTACTTGCACTGATGAGAAACTGGACGGGAAACACAAAAACATAAATACAATAAACACAATCCTAATATACACTGTACATGTACCGGAATATATCATCTGTTCGGAGCCAGAAATGCTCCAAACAGATTAAAAAAATTAGAAAATCTTCTTCTTTTTTTTTGTTTCCCTGAATCCCCCCTCTGCTGTTTGATACAAATATCGATTTCATTTAAAAAAGAAGATTAATCTTCGCTTTTTAAAATTCTATGCTTTTAGAGCATTCATGAATAATCTCTTGACACAACTTATTTCACAAAAGTGTAGATTCTGGACCTGAGAGCAATCAAATAAGCACTCTAAAATCATCTACAGTGTATACACCTGTATGTGTAATTAGTGGACACCAAACAGACATGTATATATATATATATATATAATATATATATATATATATATATACATGTATATATATGTAGGCCTATCATTTTAACATCAAACAGCAATTAATAAACAGTCACGTGTTTTGATAGACTATGTATGTACAAACATGAAAATGCAACCTAGCACAGCTTAAGAATTTAATATGGATTGTAACTAAACAGTTACTATATATATATCAAACTTTGCACATGGCATTTTGTTTTTCTAATGGGCTGATCAAAAGTTCCCACAAACTAATTTCAGGGGTTTCATATATCTACCAATTAATCAAAGTTCAAGTAGCTCAAACCTATTACTTCCAAGGAGCTGGTTTACTAGGTGTTGCTGAGCGCCTCCCCGAGGCAGACTGAGGCCTGCTTCCCATGGATGACCTTCGACCTCCAGTGGCCGACCCAGGCCGACTTCCACTGGCCGACCCAGGGCGTTCACCCATGCTTGAGTTTGACCTTGATTCAAATTCAAGACTACCCTGAACAGGGGACCCCCACAAAAGAGGTGTTTTTCTACGTTTATCGGCCTTCTTGTCCGTCCATGGTGCCTCAAAAGAGGGCTCCTCGAGAGGTCCCCCAAAGAGGGAGTGGTCGACGTAAGTCGGTGTGTGCTGCTGAATCCGATACTTATTCCTCGTCTTAGGCATCGACGGGGAGACAGTCAGCATTGGTTTGATTGCTGGAGATTTCCTTGTTGCCCCACTAGAGGGCGTTCTTTGTGGTGGATGATCGCTCCATTCAACCTGTGCTTTCCTATTCCTATCTCTATGTGATCCGAAGAGTGTCTCGTCCACTGAGGAGTCATTGGATATCACCCGATATCTGTAGTTTTTCTTCACGCTTCCATTTGAAGTGTCCCTGCCCATGGACTCGCCCTGTACGCTGAGGCTGGACTGGCGCCCACTGCCAGTATTACGAAAAGAATAATCCGACATTTTATCTCTTTTCAGATGATCTAATTATTGTAGAGTGCCATTGAAAACTGAAACCTCCTCCACATCCTAAGTAAGCAGCTACCTCATTAACCTTTCTCATGTACTTTTTAGAATGTAACAGTGAGCATGGTCGAAAATTTGAGTGGCTTGGGTAAACTTCTGCAGGTAGAAACTGTTACTAGACCTGGACAAGAATATCTCTGTCTAGCTTTTTTTGAAGTCCCAATACCTTGTGTCAACCGCTTGGTCAGGCTCCTGTGGACAAAAGAATAAAAGGGTACAACATTAGAGAATTGCAGCCTCATTTTCCTAAATTGGCGTCTACAAGGATTAAGGTCCCTTCTTATTGACTATCGCTGGCTTCCATTTGCTCAATGTTATAAGATATTGGTGCGAGGAGAAAGCTTGCGATTTTCACTTTGATCCGATACCAAAATCATGATATTTCAGTCATACCCTCTATATATTACAATGTTTTCAGAATCCTTAAAAGGTTAAGTATTTCACCTAAGTTTAAATCTTACGTAACTTTTCTTTAAATTTGAAAAATTCTTGTTTTCGTCATATCTAAGGTACAAATTGGAAATTTATTCTGAAATGTTTGTTAACTTTTCTATTTTATGTGTGTGCAGTGCCATTTGGGTGCGTCAGTTATGTGTGAGCGCGTCACCTCACTGTACCGCGCTCAATTCGCAAAGTTCTATGTGCTAGGCCAGGGAGCATGAAGAACTTTGCGCATTGAGCGCTGTACATGCAGTGAGCCTATGAGTTCCCATTGAATGATCGATCTGCTCGAGTGGCTTCAAATTTAAGTGAGAGAGGGCCTACTTAACTTCAGCGCAACTTAAATTCAATGGCAATGTCAATACGTAACATTTTGTCAAACCAAGCCAGATTTTGACATGTCAATCATGCTAACATTAATAACATTAACAGTGACAATACATTCTCGACGATAAGGCCTGATACTAGATCTAGAGTAAAGTACCCGGTTTTCGGCCACTTAAGACAAAATCCTATTTATTGGACTCTTTACTTCCTATTTATTATCATTTCTTGTTCAAATCATGAAGCAAATATTATTATTGTTTAAACAAAAAAAGAATCCCAAAAAAAAATCACGAATTGCCATGACAATATTCATTTCAATTTTCACTGTTTTGGAGCCGGTTATCGGCCGAGCGTTTCAGTTTACGGCCGCATGTACAAAGTCAATAGGAAGCGCTCGCGCGCGCCGTGCGTTTCCGCGCTGATGTATAGCACAAGCGGGCGGAGCAGGTATCCCACGCCGGAGCTCCCGGGTGCTCGCCTGACCGAATACAGGCACAGCGGTGATCGGAAGTTGGAATCTGAATTAGACGGACATATCGTGTGTGCGTCGCGGTGTGCGGCTGCATTGATTTAATACATGCGATGAGATACTTTTGAGCTTGCGACCAGCACAAAATCGGGCTCTATATCTGAAAAATACTTGAAATAAAGCTCATTTTATGGTACCCATATTTCTTACATACACGCGTAAGAAACTCGATTTTTCAGATTTGGGTAGCGGGACCGCCGATTACATCGCTAAAATGTGATATTTTAGGATATATTTTTCATTGCAAATGTGGTTCCTCAGTGACGCTGGCATGGTATGGTATGTGTGGGATATGTGTTATGCTATTATGCGATGTATTTTAAAGTGATTTTACCATGTATCCAACCAATTGGCCGAAAACCGGAGCTGACCGATTACCGGGTACTTTACTCTAGTATCAAATTCATTTTTCATTTCATTTCATTTTTATTTTTTCATTTCCAACAAAATACAATACATAATAACGAATACAATGTAAAAATGTATACACAACATAGTTACATTGGTAATAAACAAATCTTAATAATAATTTATACATTTTAGGAAGACAAAAAAATGAATTTTAAAATTTAATTGCATTTATATTTTGAGAAAATGGAGAGGTCCACTGAAAAGCAAAGCTTGTAAATTGTGGACCCCTCAGATAATAATGAGAATTCATTAGATGTCTAGCGAGTAAACTCTTCGTAAAACGGCATAGGAACATCAGAAACATCTCCGTATTAGGTGGACACAGCTACGACTGTCTCGTCATACGTTTCGAGAAACGAGTCCCTGGCGCCGGCTTTATCAGACAAGTTCAGATCGAGTAAGTGATAGTACGGAACGGTACGCAATTTGCGATACATGGGACCATACTATCAAAATTATTACAGTAAAACGTAAACAACGTAGTTAGACCATTACTTCATTCTACATCACTACAAAGAGATATGGAAACAAAATAAGAGGAAAATGTGAGCGATTCTTACCGGTAAATCAAAACAACTATATTGGACGTCGATATGGGATAAATTCGTCGAACATATTCCAGGAATATTGTTGTCACCAAAAAAAATGCGTCACCTACCAACGAACCGTGTCACGACGGATAGGATACGCAGGCGGATTTGGTACCGACGAGAGTACGGTGCGGCGCCAAATCCGTCCGGAACATTTGGTGCTCTTGATGAAATAAACCGTCTCCACCACTGATCTCTACGGAAATTTTTTTCCGTAGAGATCAGTGGTCTCCACTAATAATAAAAAAAGGCCAGCTTTAGTTGGTGTTATTGTTGTTGTCGTTGTTAATTGTTCGTTCAAAGAGCAGTAAGCCTAGGGCCTAGGGCCTTTATTTAATGTCATCATCATCATGTTATATTTGATATTTCTATCTTGGCCTGGATACCCTAAAATTATGATGAGGGCTCGCATGCCTCTGGGCAACAGTGAATTATTTTACTTTTGCAAACTCTGGCCCACTAGGTAACAATTTCTTTCCTCTTCTCCTTGCTTCTCATTTCCTTCTCTTTCTCTATTATTTTATATATAATCACTCTCTCTCTCTCTCTCTCCCTCTCTCTTATTACCTAGTGCTATATAGGCATGTAGGATATCCGATATGTAAAGCCATTCGAGGCGAAACTATATATCATTATACATCGTTTGAAATAAGAAACACATTCAAATAAAAACAAAATTAGTCTGTCTGTACCAGACCGTGAAGCACGGACATTTTGTCCCCATAAACCGGTATAATGATATAAATGTCCCCATAACCGGTGTATATATATATATATATATATATATATATATATATGTTCCCATAACCGGTGTATATAGATATACCCTATCAAGCCCAAGCTATTATGGTCGAAGGGCAGGGGCGGCCGAGCGACTTCCAAAGTGAGGGGGCACCAAGAAAAAAGGGCAATTTTCGGAAAAAGGGGCACCTACAAAGGAAAATCTCAAAGCAGAAATATACATTTTTGTTTACATTGGTATACATCCCTACATTGTGGTTACATTATCAGGGACGATTTCAATGTTTGTTTCACTTTATACACCATACAAGAAATTGCTTAAGCATCATACTTTGCAGCTGGGTGTGATGGTCTGACAATTTATTTGGGAGCGAGCAAGTTTAAGATTTGTGAGGTGGCAGTGCGAAAAATTACATCTTATTATGAGGAAAAAACCCAATCCCTGTGGGGTTATAAGCAAAGAAAAGACCTTTGATAAGTAACATAAGTATGTAATATCTTGTTTATTTTGAAAATCAGTGGTGAACATCAAATGAAGAATGCAGGTCACAGAAACCTTGAAAGCGGAAAAATACCTTTGACTTTGTAATGTTATGTTTCGTTAATACAGATTATGACTTCAATTCTATTTCATAAGATGGTGATAATTTAAACATAAGATACATTACAATTTATTCGATTTAATTCATTCATTGTACACCGTACTATACATATACCATATAAACAAATAATTTTAATATTGCTGTTGCTATTGTTGGAACTGTAATCTCCGTATTATAATAGCATTACTGTTTAATTTGCCATCTATACTTATATTTGTATGAATATAATCCAGCAAAGTTTTACTTAAAATGGATTAATATTAATATGCCAGTACTTAATTGATAGTGTTTCAGACTAGTTCAGCTAGTCTTTGGTATTTCAATTTATAACAACAGTATAGTGTACGAGTAAAATGTATTGTCATGTATGTTTTACAATGCTTTTATATAATATTATCGCTTTCAATTTAATTATTTTTGGTATGACATAATTTTGTCTTTCTGTAATGATTTTTTAACGCTATAATTATATAAACATATTTTAAAACGAAATAATCATCCTTGGAGAAGAAGGGTTATTATTTAACTGTTAGTTTGTATAAAATCTATATTCTGACCACTTGCACACTTTACTCCATTTTCTATCATATCTTGTACTATTTATTGTATTGGCTGACTGCTCCCTGGAGAATAGGATGAGAGGGGTTGTCTCTTTCTGTCTCTCACTCGTCCATGTGTCATGTTCACAATTTCTCCTCTCTTTTCATAAAGGATAGACTCGGTCTGCATCTCCACACTCTGTTTTTCCTATGTAATATATTAATTTCTTTTACAGGAAGTAACTCAAAATGTTGAACAGCAAATTTACCAACAAATTCTGATTTTACAAACTGTCATTAATTATAAATTCTAGCGTCATATCATAAGTAGAAACAATATTTTTTGTCTTGATAGTCTTGTTCCAATATGCGAACTTTATGAGGACAAAAAGTGACAAACACAAAAGTTACAAACCAAAAATCAATCTTACAAAATAGAGGGAAAATAATCAAAAAGACGAAACAAAAGATAAAGCAAGTTATCAGTCAATTTGAGGCATTATGTATTTTATTGCTTTATGATTTTCTTTCATACAATTTTTTTTGCAATTTTTTGTCGCAGTCATTACACACAATCTACATGATCTGTATTCATCAATTATGAGTTATACTACATACTCGTTCATTACAATAGCCCCCTTTGTTGATGCGGTAAAGATGGTTCAAAAGAATAATGTAATAATCGAAAGCGGAAAAAAGAAGAATGGAAAAAAATTGTTGAAAAATCTTTATCGAGATAAGGAGCTCGAGAATTGAAGGAATTTCTTCCTCTTCGGTATGCCCATTTACTTTACATGGGTGGAGAGTGGCAGAACGTATATTAATATAGTAGCAAAGGACATTAACTGTTGTGACAGGGATTCGAATCCCGTACCTTGTGCTCAATAGTTCAGGAACCTATTGACCAAACCGCGAGACCTTGCACTTGTCATGCTTAAGGTGAAAAGTAATAGTTTAAAGGGGCATGTTCAAAATATGTACGCGCTGGCGGACGGTGCTATAACTGCCTATTGACACCAACAGAAATCGGCTCTTACGCAGATGGCATGTTCCGATAACTACACTGGCGGGAACGTGCCCCTTTAAGGAGGTCAATGGTCACTAATCACAGGCATTGGTCGACGTCAATGGAAATCCAATCGACTATTTTCATTTAGAGTAACAGAAATATATTTTCTTGCCAGTTTTGATATCTGTTTTCATAAATGATGGCCCGTCTGAATTTAGATGTATGAAATCAAAGCGAATGTTTTAACCTTTGGTAGAACAGCGATGATTCATTTGACAAACTACATTAGGTAGTAAAGTAAACAGAATGACAATCAAAGTGATAAATGTAAAAACGAATCAGAGAGAAAAAATGAAGCTCAGCTATATTTGCTTGAAGACGGACTTCAGATTCAATACAAGACAATTCTTAGGGTATGAATCGTAAACAAAGAGAAAAAACAACGGAAAGTATCTATTTTTGTAGTGTAATAATATGTCATAAAGGTTTTCGTAAACTCCTTCGTCATGTTTTACAGACGTATACATGTACAACGCAAATTTTACTTTTACGGTCTAAAGAAAATTAAAACCATCATACTTGACCTCCATGTTAAACGTTCCACTACATATATGTGATGCCACCTTTAAACCTGCCTATCCTTTTTATATGTGTACTACCAACACCTGCACCATCTGCAGTTCTCCAAACTGGAGGAAGCCAATTGGCTTAAGCCAAGAAGAAGAAGATAAACAAAACAAGGAAAATTGGCAAATGAAGAATGCGATGTGTTTTAACTACAAATACCCACCACAACCAAATAACACAAGTCCTTTGTTATACTACTGTAATCAAAGCCTAAGAGCATTAAAACCTTTTAATATTGAAAGATAAATTGCTGGTATTTGTTCAAATTCATCTGTATTACGTCGGGGGAAATTGCTATAGCCATGACTTGTTGGACAAGCTTGAGCTTTTTCGTGTAACATTATTGATTCGAATATGGGTAAAAGCTCTCTTCTCGGCAATGATTTCTCTGTAATTTCGGATATTCATATTAAACTGTGGTAGAATCGCCATGACGACAGGAGGAAAACGGCATAACAGTGTGTGTTTATTAAAGGGATGCAGTCTACATCAGAACACATCAGGCGAAAGACTAATAACATGTTTAAAATACGTAAAGCTTCTACTATAATATTGCAAAATCCTCTATTCAATATAATGATTGCAAAGTTATGAAAATAGCTATAGCAAGAAAGCATGATAAATCCTGAAAACCAAAACATACGTATATTTAACAGCTAAATTATTCATACCCTTACAAAGTTTCGTAGAATTTTAAGTATTTAAAATGTTTATATGTCCTGTTCTATGTTTTTTTGAGAATGAGACATATTAAAGTTTACGTCTGACAATTTTCAACAAAGGCCCTAATTGAACCAGTCCACCTTTCATTCTTGAGTTGAGAAAAACACAAAAACATCAAATTTACAAGGGGTATCACTAAGTCTTAACTGTATAGCAACCATGAACCCCAAAGATAAGCAAACGAAGATACAACATAGTGAAGGTAGAATAAATTCGAGAATGCATAAAATCGAATATAATAATGAACGTTACATAGCGCTTATCACTCACGTAAAGAAAGTCCCTATGCGCTTAGAAAAGAAGGAAAAGGGCGGAATAGAAATGAACATTTATGCCTGGTACATGAAAATAATTAATTACAAACAATAATATGGTGTAATAATCATGCAAAAAGTTATTTCGTTTGCATGTTTTGTCGTTGTTGTATATATTTTGCAATTGCGAACAATAGCCAGTTCGGGTATCAACAAATACTTGTTGGTAGCTGCTCACATCAGCAGTCCCTGTAAGGGGAATCCCAAATAGGCCAATGAAATTGGTACATTGGATATGCCTGTCTCCAATAGATACATCAGATGTCAAATTAAATGCTTCACGCGATGACATCCTCAAAATGATAACACTTAGTAGACTTTAATGTCAATAATTTGCGCTCTTCACCATGTTCATTAGGCCTATATAATGATTGATGATTAATTTAGATTAATACAAATTCAAGTACCATCACAACCATGATAAGGTAGATTGACAAGTCCAAGGCATACGTTTTATATCGACATCATCAAAGTTTTTTTTAGTATAGCCGATGATGATGATGATGATGATGATAATAATACTACTAATAATAATTATAAGAATACTGATGATAATAATAATACAAATAGCTATAATAATAGTAAATTTAAAATGATAATGATAATAATAACAATAATACTCAGGTAGCCTCATCGGCATCAAGGCCCTTCATGAGGGCCCTGCATTCATTATTACTCCAGCAATTACCAGGTACCCATTTACACATGGTTGGAAAGTGAGAAATGCGGATTTAAACAAAATTAACTCACGGGGAATTTCAATACCCTGTTGTACATCTCATAGAAATCATGATCGAATTAAGTAACACGTCATTGCAATGTCCGCGCGTTGAAAAAGGTATATGACAGGGTAATACTGAAATTGCATTAATAAGTTGATTTTTCGAAAGTCTCTACGATATTATTGTACTAACATATTATTATATTAACACATCACTATTCCCATAACTCAACATTCTCTTTGATAACTGATATTTCAACTTTATTTATAATTTTCGTTGACCTTCACTCGATGCCTCGATGGTATAAATGTTGGACATGCGTCTCGCTGGTAATTTGTAGAGCACACCGGCTACCGTGGCAACAATCAAACACGTGATACACGCGAGCACTCCAATCAAAATGCTTGATTCCATACCGACGGGCGAACCATCTGCAGTTAAAGATCATAAATACCATGATTGAGAAAAATAAACAAAAATAATCCCCCATTCCAAAAACTAGGAAGAAGAAATGATTTATATAAAGGATTGTGCTAACTGTCATACATATTCTATTTTTTATAAGTTTACACGACTCTTTCTGAGAAACAGCCGCTGGCTCGTTGTAAGCGTCTTGAATGAAGAAAAGAGGACATTTCAGCTTAAAACTTCTAATCCATTAAGCTAATTTTATGTAATATTAAAACAATGTCATCATGAAGGACAATAATTAAAATACAAATACAACTTAAAAAGATTTGTACAAGAAAACAAGAGGTCAAATGTCAAAGGCAGGACAAGTAGGTGATACATTGAGGCTGCTCCCCATAACTGATGTAAATTTGGATATTTGACCATATAGAATACATACAAAGTGAAATACAAACAAACTACACACAAAAACAAGATCACAAATAGAAACAGTTAATTGGCATAGGCCTAGTCTGTAGGTAATAATAATAATAATAATAATGATGATTGGTATTTCTAAACCGCTCATTGCCCATCAGGCAAGGTGGTCAAGGCGCTCCAATATTTATACCCCGGCTCTGCTAGCTAGGCTACCGTTTTCAGCTCTCACATCATTCAAGGAATTCCTTCATACCGGTACCCATTTAAAACCTGGGCGGTGAGGGACAAATGTATATTAATATCTTGCCAACGTGGCGAGACTCGACTGGCGAACAACAATTGTGATCCACAGTCCAACGACGTTTTCACTAGACAAACACGACACCTTTAAAGTTAAGTGTTGTGTTCTTTAAAAAGGTAATTACCAGAGTTGTGTGCTTCTGCGGCATGCATAGTGTGGGCAGTTGAAAGAGGTCCATTAGAGATTAGATGAGGAGTGGAGCTTGCTCCACGGGTGTGACGGCTACCACGTCTGAATCTAGACTTACATCCCTGGGAACACCGAGACGATGGGTCATTATCATTACAAATGAGCACCTCACAGTGAATGAAAACCTGCATGGAAAAAAATGTATAACCGACCTCTATGACGATAACATATTCTTCAGTTGCGTTCAGTGGTATTCAGTAGTGGTTCTAGAGGGGGAAGTCTGCGGGTCTCCACCCTCCCCCCAAAATGGCACGAATGACAAGGTACCCCATAAGACAGGGTCGTATAACGTCATCTATCATTATAAACATCAAGACGGATGTTTCAAAATAAATCTAGTACTGGTCATGTGATTACGTTTTGAAAGAAATGGGCTAAAATTTACATAGAAAGGTTTAATATGTATGTAGTACCAGGCGCGTACGCAGGGGGGGGGGGGGGTGTTTAGGGGGTTAACCCCCCCCCCTTTGGTATTTTTATGTGTATAAACCCCCCTAAAAAAACAACCCCCCCCCCCTTTGTTTTGTTTTGTTTTTTGTTTTTTTGAATAATAAAAAAAAAAATAATAATAAAAAAACTATTCCCATGGCTGACCTAAAAAAGGCGACCTTAACCCCCCCCCCCCCCCCCCCCCACTTGAAAAAATCCTGCGTACGCGCCTGTTTAGTACAGAAGTGTAAGGTTTGAATTCATCAAAATTACAAGAAATTTAATCAGAATGTCTGAAAATGCTGACAAATGACTTTTTAGTGCACACATAAACGAAATTGGGTAGGCTAATATTAAACTACACATGTGGCATGGAATATTTATCAACAAATTTCCAGCTACCTCCGGTAATATTTCGAAAACGATGATCAATCCTTTTCAGGGGATATTTCTACTGTGATTTTTTTTTTCAGGCAGACGTGAACATGATACAGCTCCCTGATAACAATTATCTCAGCTTTACTCAGCCTTTTCCCTAAACCCTATTATTGTCATTCTTGCAGCTGTTCGCCTTTTACACCTTTTATTCATTATTTTGACCACTTTTGTACTTTTATCATTGTGAAAAAGAAATGGAAATACAAAAGAACGATGGCCTTATCATCCCACCTGTGTGTTTTCCTCGATGAAGGTGAAAGCGTCGATGACAAACTTGGAGAATGACGGATCAGAAGAGTACAATTTGCGCATAGTGGGATCCCCTGAGCATCTGAGATAAGAAAAGAACAACAAAACCGTCATGATAAATGTGAGTTGATCTAGTGGTTATGTCCGTGCATGGTCAGCTCTCTTGTCTGTCTTCGGCTACAAACCGACTCTGATCGAACCCAGGTTCGATACCCGGGGTTCAATACTGAGTCACCTCTTTCGGATGGTGACGTTAAAGGTCGGTCCTAGACGTACATAATCTTTATCTGATTGATACGCGATTGGTAAAAACTCATTACATACATACTCACAGATACTCACAGATACGCACGCACACACCCGCACACACCCACACATGCACGTCCGCACACACACGCCCACCGCGAATTGGAGAATTATGAGCAATTTCGTGCTATATGTCCCGTACGATTATGGCCATTTTCAATTAATCAAATCAATCAGTGATTACTTATATGACTTGAGATTTTCAACATTATTTATACAAAATGAAAAATCTCTTGTTTGACCACTGAATAATGGTGCTTTCTGGTATCACTACCAGTAGTCATAGTCAATAGATCCTGCAACTATGGATGAAACACTCGAATGAATGTTATTCTTTCGTATGAATCACGGAATCGCACTGAGTTTTTTTTTTTTTTTTTTACTATCCACCGAGAGGTGACGTTGCCTAGTTAAAATTTCAACAGACTGTGGTTCATATGGTTCGCGGTTCGAATCACAACCGCTGCACAACTGACCTTTGACGAAACATATTGATTTACATTTGTCACTCTCGACCCATGTGTAGTAAACGGGTATCTGGTAGGAAGTAACTCCTTGAATGCTTTCAGCGCCTTATGCAGCATTAATGCTTATGCCAGGGTAAAAATACTCTTTAGGTGCGTTGAGCATCCTCTAGGATGGATATGGACGCAATATAACTTGAGCTATTAAAAAAATATATATATATCATTACGGTAGCATTATGCAAACATTTCTTTATGCTTACTGATTTTGTGTATTGTTTTATCAAAAGCAAAAGTCATTTTCGTAATCTGGCAGACCTGGCAGCGGTCAGTAAATGTCCATAGAAGATGTGCGGGTGAACGTGTGTACATTTGTGTGTAAATTTGTGTGGTAAGGTGTGTGTAGTGTTTTTAGGTGTGTGTAGTTTGGAGGTTATGTTATCTTGCACATTTACAGAATGCATATCTACATGGGTTCCCATGGGGAAATTTGACAGGGAAGCACCCGAGAGTATGAATAGGTTACAAGTAGAATAGGTCTAAGTAAGGCGAATAATAGTGTTGATAATTAATAGTTTGGACTTTTTTAACCTAAGGTATACACATTAGTATTAGTACTAGCAGGACCCGTGGAAAATCCACGGGATTCTAGGGGGAGGGGAAGGGGGGGGGGGGGTAAATAGATAGAGGCCTATGTGCTATCGAAAATAACATTTGAGTATATTTCTATAGGGTCTAGACGATTTTTATAATATACGCTGTTTGACCATGAATGGCGCACCATAATACCCGTTTTAATGCATTGATTTGATGAAGAAACTGTATCATAGTAAAAGTATTCTATCCAGCATGGCTATGTTTTTAAACGACCATCTAAATTTTCTTGTGAAAACTTTACTTTGAAGAATGATTCTCAAAATTCAGTAAACATGGACAGGATAATACTAAAGTAACCCTGCTTGCATACCTCTATTCACGTAAACTATATACAAATCACGGCAATCATCATAGAATCAGAGTCGAAGCGACTGACCTGTGCGTGTGACGGTCTAGGTGGATCGAGTGAGGGTGGGGGGGGGGGGGGGGGGTAAATAGAAAGAGGCCTATGTGTATCGAAAATAAAAATTGTTTATTTTTATATAGGGTCTAGACGATTTTTATAATATACGCCTTTTGACCATGGCGCACCGTAATACCCGTTTTAATGCATTGACTTGATGCAAGACCTGTATCACAGTAAAAGTATTCTATCCAACTAATTGAGCTAGGCTATGTTTTTAAACGACCGTCTGATTTTCTTGTGAAAACTTTACTTTGAAAAATGATTATCGAAATTCAGCATACAGCTATACCTGGCCAGAACAATCACAATCTAAGCTCAGATTTCATACCTCTATAGTCACGAAACTAAATTGAAATCGGGGCAATAAGAGCAATTTAAAAGTATTGTCGAGTCTACTCAGCTCAGCTCACACACCTGCACAGTCACAGTTCCGACAGTGCTCGACAGCAGTACAGTGGAATGCACTGGCACTGTCCGTGTGTGGACTCAAGCCTCACTACAGCTACCTTGCAACTGCAACTGAACCGTCAGTTACGATGAATTTTTTTGTGAAAACTTTACTTTGAAGAATGATTGTCAAAATGCAGTAAACATGGACAGGATAATCACAAAGTAAGCCCTGCTTTCATACCTCTAATTCACGTAAACTATATATAAAAATCACGGCAATCATAGCATGAAAGCCGAATGAAAGTCGCAAACCGACTGTGCGTGTGAAAATCATCGCATACATTATACACGTACGTTCGTACACACAGCAGGTCCGTCTGGCCTACAAAGCATATCCTATAGAAAGACGCATAACACGGGTAAGCTCAACCCGTTGAGCGTACATATTGGCTACATCCCATCACCTGAGTGAGTCATCGGTTAATGTGTATTTCTCTGACACGCTACAAAGCATATCCTATAGTAAGAGGTCTAACTCGGGGTACGCCCCACCGTTTTAAGATGAGATTTGCTACACCCGTTGACGCGAGTGGATTCTCGTCTATTGTGTGTCTCTTCGTCTCGCTTTGATAATTAGACAAAGGGCGGGAGCACGCCTTGTCGACTCGCTCTCGCGGGTTCGTCGCCTCCTTGCAATCCACGACGCCCGCCGTGCTTTCGTCCGCCCGTCCTCCCCGGCACCCCGAAACCGCGTTTAACGGCATGGACGTCGATCGCGGGGTCGAGGTGGCACGCGGGTTACGGGAGACAGGCCTCGTGGAATTGTGCCACCGTCCTCGCGCCACGGCGGAGGCGAGCGGGGTCGGGGAAAGACGTTTTTGAGGAACCGTTTTGTACCCGTTTGGACCCCCCGCGGCGATTTTGCGTATTAATAATATAGACTAGGAGGACCCGTGGAAAATCCACGGGATTCTAGGGGAGGGGCTGATTTTAGGGTCTGAGAGGGGTTGTAGATAGATTTATGTGCGATGAAAATATAATTTGATAATTTTTATTTATTAATCTTGTAGTCTTTTATTTTTTTCAACTTTATTTATTATTATTAAATCGTTAAAAAATCACCCTAATGCCAACGTTTTATGTATTGATTTGAGCTGTTTTTATTTGTACCATTAATATATCTTGAAGTTTGGAAGAAAAAAGTGGGTTTAGTATATTGTAATAAATTGATGCATTTAGTGTATTCTTATAATTGTTTCATTTTTGTTATTTATGATCGTTTCATATTTATTCTTACTATTCTCTTCTCGTCCCCACTTTCACTATTTATTATGGTTGTTGTCATGCCTCTCTCATTCTCATTTTGCTAATTTTGCATATCTCTTTAGTTTTGCACAGCAATATCAAATTATGCAATCATTATTCTTCCTCAATCTCAGCATACTTTGAATGACAGTGACATTCATAATACATCATAAATGAAATTGATTTCAAAGATTGATTTTTAATACATGTATGAACAATTTCATATTTATTCATAACTTTTCCAAAACGAACTATGAATAACAGTTGATATATTTGACTTAACATATTCTATTAAAACTACATGTTTTAATAACAAAATTATCATCAAAAATATATTTTGTAAGTATGATACAATCAACTTGAAAGAATTACATTTATTGTATTACTTACTTATTTCATACTCTTACATATCATATAACATGATAGTTGTTATATTTAAAAATTACAATCTAATATACCTCAAAAGATTTATTGTTTGCTTTTAAAAGTAAATTTTGAAGTAAGAAAGTTGAAGTCAGGCAGATAATGAAAGATGAAGTTTACACTTACATATGTATTATTATTAAAGGTTAAAGGTTGTAGTGCAGAACACTTTAAATTGTATCATTTGCAATCAAAACATGAAAGGCATTTTACATATTGACAACATATGTATATTGTATGTATGTTTACGAACAACATTACAGTAGAACATGTAGTTTCATGGAATATAGTTTTGGTAATAAATAACATAACATTGTAAGTACAGAATCTTCTACATATGTTCTTTGGCAAGTGTTGCTGCATAATTAGTGATATGAGATCGTCTGATTTTTGGTAATAAATAACATAACATGTAAGTACAGAATCTTCTTCATTTGTTCTTTGGCAAGTGTTGTCTTTTTACTTAGAGATATGAGATCATCTGTGAGTTAAAGTCAGAAATATTTTATAATGATTTTTAAAATAATTGTAATAAGAAATCACATTTATTTTGTATACAAAAATTGAAAATGGGGTTGTGAAGAGTATCTAATTTTAAAATCAGATTATAATCAAACTTTCATATACTTGAATTTTATCACTTATTTCCATAAAGAAATAGCGATAAATATTTGAGCATTTTATTCTTTCAAAAAAACTTGTTTCAGGAAATTGAGGTAGATCAATTACTTAAAGGAATTAACGATCTGATTTATGAGAGGTGTCGGCGAGGTAGATTTAAATATTTAAAGGAAATAACGATTTGGAATATGAGAGGTGTCAGCGAATAATTGATTCATACATATAGTTAAATAAGAACATAACTTATCTACATCCACATATTTATAAAAACACTTATTTACATAGCTTTGGCAATATATGACGTATCCTCGTTTAAAAAAAGATTAAGTACTGGAATACTTTTTGATATTGCATGAATGAGTCAAGATAATCTTTGAACATAGAAATTTAATATCTATTTTAAAGGTTAAAGGCTGTAGTGCGGAACACTTTAAATTGTATCATTTGCAATCAAAACATGAAAAGGCATTTTACATATTGACAACATATGTACATTGATGTATGTTTACGAACAACATTACAGTAGAACATCTAGTTTCATGGAATACAGTTTTGGTAATCAATAACATAACATTGTAAGTACAGAAATTCTACATTTTTCTTTGGCAAGTGTTGCTGCATAATTAGTGATATGAGATCGTCTGTTTTTTGGTAATAAATAACATAACATGTAAGTACAGAATATTTTTCATTTGTTCTTCGGCAAGTGTCGTTTTTTTACTTAGAGATATGAGATCATTTGTGAGTTAAAGTCTGAAATATTTTTTAAATGATTGTTAAGATGTGATTTTTTATTTCATTTATTTTGTATACAAAAATTGAAAATGGGGTTGTGAAGAGAAATATAATTTTAAAATCAGATTATAATCAAACTTTCATATTCTTGAATTTTATCACTTATTTCCATAAAGAAATAGCGATAAATACTTGAGCATTTTATTCTTTCAAAAAAATTTGTTTCAGGAAATTGAGGTAGATCAATTTTTTAAAGGAATTAACGATCTGATATTCATGTATGAGAGGTGTCAGCGAGGTAGATTTAAATATTTAAAGGAAATAACGATTTGGAATATGAGAGGTGTCAGCGAATAATTGATTCATACATATAGTTAAATAAGAACATAACTTATCTACATCCACATATTTATAAAAACACTTATTTACATAGCTTTGGCAATATATGACATACCCTCATTAAAAAAGATTAAGTACTGGAATACTTTTTGATATTGCATGAATGAGTCAAGATATACTTTGAACATAGAAATTTAATATATATTTTCAGCGTTATGTAACATAGTCAACATACTTAAGAGTTTCAAGTATGGATTACGTCTATGTTATTGCATAATGTAGGCAAGATATCCAGTTCCCCATTAAATATTTGAAGGTAAAATATAAATGGAACTTTGATTATAAAAATACATGTAAAATACGTCATGTATGAAGTTGAAAAACCATATTTTTTCAGATCATTTTGTAGAGCTATGTTTTTTCTACTACTTTAAAAAAATGAATGCATAAAATTGTTTACCTTAATAAAAGATCAAACATTTTTGCATTTTTACAACAGGTCAAGTAGTTTTGCATACTGCAAGAGTTAAAAGAAATTATAGACATCAGTTTTTATCTCAAGCGTTATAAATGTATTTAGTTTTTAATTAGTTCTAAAAAAGAAAGAGAAATTGATATAGGTCTAAGTAATCATATTAAATAACTTTTCTTTTAATTTAAGTTAAATGAAAGTTTTAAAAAACATAGTTGAGAGTAACTTTATTCTTAAAAAACATTGGTGTACGGTACATGAGATCTTACATTATGGGCGATGATGAGGTGTCGTCGTTGTTGGAGGTTCCTTGTTCGTTTTCCTGGGGCGTTGATTGGGCTGTAGGTTGAGAAGCTATTATATTATCAGTACTGTTCATCACTTCATGTAGCTCATCTAATGCATTGCTGCATTTTAGCTGAGAGTCATTTAACGCAGATTTCATTTCTTTGAGCTGCTGCCGCAAGTTGAATCCCTCTTTGGAAAGTGCTGCATTCTTAATCAGTAATTTTTTGTTGGCTTCTTTCATACTTGATGTTTCATTTTTCATTTCTATAGCAGACTCTTCTATCTTTTCTGTCAGGTTTATTATTTGCTGCCATGCGTCCTTCGTGAGAGTCGGTGCCTTGACTTGCAGCCTTGTCATTGTACCTAAGATTTGTAGACTTGCTGCTTTGATGTTTTCCATGTTGCCTGTTTTAATGAAATGAAAGTTACATGTATGAGATACTAATTGAAGTAATTGGTATGATTTTCATGATTTTTTGGGTGGTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCATGCGTCTTCAAAGTAAATGGACGCAAAACACGTTTTTTAATTTAAAAAAAAAATAGAAATAATGACTTTAAATAATGTGATAGGTTATCGAGAAAATGTCCTTTTCATTTTTAAGTTGTGCTATTGTCCTTAATTCTTTCGTTTACAGCGGTTCTCAGAGTTAATGATTAGTTGTGCTGATATCGATTTTAAATTAAACTAAGTTAGTCATTCGCTTTTTTAGTTTGCTTTGGTGTTTTGCTCTACGGTTTGCGTTGCAGTCGCTGGCCGGGCGTTCTGGCTCGCTTGTGTTGCTCGGTATGGTAAATCGTTCAAAAATAAATTTAGGGTTCGAGTGGATTCATTTTAGATTTTTTAAGTTCAATTGTCTTAGCCTAAATGATATCAACTGACAAAATTCGTCACTCACCTACAAATGATGTTTATCGCGAAGTATTTATTCCAAAATTTCTTTTATAACAAAAGTTATCTCCAAATCCAAAGTACGGAAACTGTGAAGTGTATATCAAGAAATTAAATGGCGTGCGAAAGGCTTTATAAAACGCCGTGTGTGAAATACGGTACCATGCGTTGCATTTATTCGGCGGTCTTGTTGCCTTGCTTGTATTGCACGTGCTAGTGATTTTGTTTTAGATTTGTAAATTTAACTGTTTTAGTCTAAATTATTATATTCTAGAGAAGTTTTTATTCACCTAAAAATGATGATTATCGCGAAATATTTGTAAAGGAATTTTTTTCTCCCGCACAATATTAATTATCTCCGAAACCAAAGTACGAAAACCGTACAGTGGATATCAAGAAATTGATTGGCGTGCCAGATGCTTTATAAAACGCCATGTGTGAAATACCATGCGTTGCATTTATTCGGCGGTCTTGTTGCCTTGCTTGTATTGCACGTGACGATTATTTTGTTCAAATATATTCAAATTGCTAGTGAATTCGTAGATTTTTAAATTTTACTGTTTATCAAAGGTGCTAGTTTATTCGTTTTAGATTTTTAAATTCAACAGTTTTAGTCTAAAGTATAATAGTTGACAGAGTTTGTGCCTCACCTATAAATGATGTTTATGGTGAAATATTTATAAAGGAATTTTTCCTTCAGCACAATAATTTTCTTTAAAACCAAAGTACAAACAAAAAAGTTAAGTGGATATTATGAAATTCATTGGCATGCGAGATTTTTTTATATAAATTATTGTGTCTGCAATTATAAGTAAGTGAGCTTGGCAGAACAGGAAGCTTTAAAGCATAACAAAATAGTTTGTTTTGATCGAAAGCCATCATTTTAAAGACTTATGATAGAGATTGAATTTATTGTAGTTTTAGATTTTATAGTAATTGTTAGCTTTCATTTGATTCATGATTATATTTTTGACAGTTACAATTATTTTTAATTATTCATCAAAATATTTGCCATGTAATGGATTACGTTTTTGGTAGATGCGTGTATTTTTTTTTGCTATAGATTTTGGACAGTGTAGTTTTGAGGAATAAAAACAGTAAGGTTGAACCAAGATGTTTGACATAAATGATAGAGAACCTTTCTTTGATGGGGGGGGGGGGGGGGGGGGGGTCTAACATCTTTGTTTTTAGAAAAACTTTTTTTAACAACAATTTTTAAAATTATTTTGTTTAAATTTAAAACTGGTAGCTATGATTTGTTTTTTGTTTTGGTGGGAGGGGGGCTCACGATAAACATTACAAAGAAAGCGGACAAAATATACAGTTGCTCATTCTACCTTCAACAATTGCTACGGTAGACAAAGGAACATCCCTTCACTGCTATGTAATAATGTGTAGAACAGGTGAATATGCATGTGACACCTGTTCAACAGATTTTTACATAATAGCAGTGAAGGGATGTTCCTTTGTCTACCATAGCAATTGTTGAGGGTAGAATGAGGTACTGTATATATTTTTATTTTTTTTTGTCCGTTGTGTTCATAGTGGCTGGAGCCCTCTTGATTTTTTTTTTAGACAATTTAATTTTTTGTTAAATTGTTTGGTATGATACTTTTCAACACTTTACTTTGCATCAAATTTCCTTGAATATCTGGGTGATTGGGCCTTATTGTATTGTTATTTTTTTTAATATCTCCTTTTTTAAAAAAAAAAGGACATATTTAAATCAAATTACCTACCAAATTTATGAAAGAAAAATTGTTATTTTAACACAAATTTGTTTGGGCTTGTCCGCCACACATTTTTTTTGGACAATATGCATTATTTTTATGTTCAAGTAAGTTTAAATTAAGGACATCGGTTACTTTAGACATTTGTGGATTTCATAATGAAGGCAGGGAAACAATACATGGTTCTCATTTCATTCATTCTCAAAATGCCAAAAAGAGCCCCTGAAAGTCATAAAAGAGCCACAAAAATAATTGCATATAGAGTTCTAATTAAAGAGTGTCATTTTCATGCAAGCAAGCTTGAAAAATAGCCCTCAAAATGTTTTTATTTTTTTCCCTCCTGAATGAAGGGATTTTCATTTAATGCGTTTGCATTGTATGTATTATTTATATGCATTAACACTTTTGTAATGATGTGGTGTATCTTTTTTTAAAGGATAGCATAGTATTTTTAGTCAATATTACAATGACAAAGAGGAATTTTTTAATATAATATCTTTTTTAGTCAACCATTAAAATGACTAATGCATTCTTACTTACCCTTTAGAGGAAAAGAAATTTCAACCTAGGACTTCACGATACACGATGTTCTTTGTGTAGTATTTACTGTGTTTCTTTCCTTGGCCTGAAGTGTTGTCCATCATCACACGTATGCTGTTGAGTGATTTAGCCCTGGACAATGCTACGTATAGTTGTCCATGTGTGAAAACTGGGTACGGCAGATAGATCCCAACTCTGGAGAAGGTTTGTCCTTGAGCCTTATTGATTGTCATGCAATATGCTAAACGTATTGGAAATTGTACACGGTGTAAGTTAAATGGCATGGTGGTGTCAGATGGGGTCAACTTTATCCGTGGGATAAATACACGTTTGCCTGTTATCAACACTTGCGCATCAATGACATTGTTGTACAAATTGCAAATTCTCAATCTTGTACCATTACACAATCCATCATTGATATCAAGGTTACGAAGTAACATCACAATTGCTCCCTGTTTCAATTGCAGTTTATGGAGTGGTAATCCTGATGGTGTCATTGAGTTTACAAATTCTGTTGGATAGCTCAGAGCATCTGAAGGGTCATCTGTTGTGATTGTATCTGTACTGTAGTACACTGCTGTATTCTCGGGAAGTTTTTGAAGTATTTGCTGATTTATTATGAACGTGGCATCGTTCTTTGGTGTGAGGATAGCTTCTTTTTCGACATTATTGTTGTCATTGTCATCATACATGTCCTTCACAATGTCTTTTGTGATTACTGAGGATGGAAGTACAATGGAATCACCTGGAAGTGCGTCTGTTCCAGCATCATACGTTCCTTCTCCTACCTCAAGCAGCCATTTCGCAAACTGTTTCTCACTTTCATCTGCTCTCATGTTTGTGGTCAGCTTGAAGTTGTGAAACTCTTGGAAGGTGTTTGATTTCTTAATTGTGTTCTCCATCAGTACTGCCGGAGCAACTCTTGGAATTACTGGTAAGACTTGTCGAAAGTCTCCACCTAACAGAACAATCTTTCCTCCAAAGAGTTTGTCATTATCCATTACGTCTTGTAGCATCATGTTGATGGCATCCAGTGCATGTTTTGGCACCATTGAAGCTTCATCTATTATGAATAGATCAGTATCTCGCACCTTTTGCGCTTGTTGTGTATCTGGTGGAACATTGCATGTAGATCTTTCATCAATAGGTACTGGCAATTTGTAAAAGCTGTGAATTGTTTGACCCTTGTGTAGAAGTGTTGCAGCAATGCCAGTATATGCTGTAGCAGACACCTACAAAAGATAATAAATTAAATAAGGTTATTTATATTTTAGTAAAAGCAATGTTATTTACATGAAAATGTAATTTCAGATAACAGGAAAATTAGAAATATTGTTGTGTGTATTTGCTATTTTAAATCATATGTTTCAGTACAGAATGAATGACAGAGAGATAAATATGTAAATATTTTTTTGTAAGATATATGGAACGAAGTATCAAAATAAATTTTCTGAACAATGATTACATTTTAAATGAATCTAATATTTTACTATTTCACTAATAATAGCTTTTCGCAACAAATTTAAGTCAATAAATAATTGTAATAAGGATTGTTTATATGTAATTGTAAAAAAACAACTTACTTTCAGGTTTTTCTTGTGGAATGTAGCAATGAGTCGGTTATACATTGTAGTTTTTCCTGTACCTCCTGGACCATCAATGAAGAAGGCAGTAGTCGGTGGATTCATTCCTTGGGAGATGGTTTGTTGACATTGGAGGATTGTGTTGATGACATGTTCCTGTTCTTCATTCATTTGGACATCAGGGTTATCAAATTGGTTTGCCTCATGTAGAAGATCATTTTGTATTTGACCGTCTTGTAGGGTTTCCAAAGGTAATCCAAGGGATCTACATGAAAGTCCTAAATGAAAACAGCAAAGAAACAACTACGAATTAAATTAATGCTTCAAGTTTTAGAATTGATCACATGAATAAATGCAAAATAAGCAGGATTTCATACATTGAGTTAAAAAACATTATGTTGTATTTCCATAAACAATTGAGGTTTTTTTTTGTAGACGGTACATATTTGATAAGACAGTTGTGAAGTATATCTTTATCATTTTATTAACTAAATAATATTTCTGTATTTAAGGACTAAACAAATTTTAATTTATTAAAGGTTTTACAGGAATTTTTAAAATTTAATAGAAAACTAAATGTCATAGTTACTTACCATGCTGTTTTAGGATATCATTGATGTGTGCAAGGGCATGTTGGGTTGCTTGTGGTTGAGGAAATGTTCTTCTGAAATCCTCTGTGAGTGCACTTTTATGGTTTTCCCATAATGAGAGGACATCTGTTGGATGATTGTGTGCCACTATTGTTGCAAAGAGCAGACGCATCTGTCTTGGCATCTTGTATGATGATGCTTCTCTCAATGTGTCATCCCATACTTTGTCACTTTCAAGAAGGTTGAGTTTTAGGCATGCTTCTTTGAAAGTAGAGGCTTCAATTCCATCGACGGTTTTCATGTCTTTGTAGCTCTTGGCCCCAGGAACATGCAGTAGAAGTAGACGAAGAAAAAAATCTTTCTTCGTCTTTAGGACTCACACTGTACATTCTTTGCACAGATTTGACATTTTTTTTCCGTTTTTGCCAACGGCTTGTTTTATTATTGAAAACATAGTGGTTGGGTATTTCTGTGTAGAAGTATTGATTTGCATCTTCATCTTCTTCATTGAGTTGGAACCAAGCAGATAATTTGGATTGTTTATTGTTTGCTCTTTCAAGACCTTCCTCCTCCTGTCCTCGTGTGAAGTAGACTGGTTGTTTGTTGGCAAGATGAACTGGTAGACGATGGATTGTGTGTGAACGTCCTTGCATAGAGTACTCTGACAGTCGCCAAAATCCTTCTGGTGCACTAACATATCTGCAAAATATAATAAGCAAGTGCAAGAATGTTTTTTAAATGAATTTGGGATATGATAAGAGCAATCTATTCAATGAAACTTAGAATTGACATTAAAAAGATTTGCACAACTAATTATTAAAATACTGATTTTATTAAAAGAGTTCTCTGTGAAGTTAAATTTTTTTGCTATAATGTACATTTTGTAGATTCATTTCAAAATATAAGCAATTTGAAGTTGTGTATTTTAAAAGTTACCTTGCATCCGTGAAAGAGGAGATCTCATCGTGATTTACTGTAGTTTGTGAAGTTTGCAATGTTGCACAGTCGTATCCTTTATAGATGTATTTGTACAGGTACTTCACGCTTTTGATGGACATGCAGGCTTCCAAGTTTATGTGTGCTTGATATTTTTTGGATAGATATGGGCTGTAAGGAACAATCCATCTATAAACAGGAAGAAAGAAAAAAAGGAACAAAAATACTTATAAAATAAGTTGTAAAATATAAAAAGAAAGAGGAGTATCATTTATGTGAATTTATATACATGTATATGTAAAATCTGTTCATGATTATTTATTGTTTTACAAAATTAATTAATACATGAAACAATATTTCATTTTATCATTGCATTTATTGTTTATGTTTACAATGTAGTTTTAAGAACAATATTGGTAAAACATGAACACTTTCTCTGCATATTTAAGTTTGGTAGATCAACTTAAATGAAAGATTAGCAATATAGAGGTAATATTTGTTATAGAAAAGTTTATGCATTATGATAAAATTCTACATACCTGTTATCCAATTCCTTGTGTTTGACGGTGACGGTCCTGCTATTGTCTGGTCTGGCATATTGGGGATATCCATCAGTAGAGACTGCTGTTCTTGTTTGAAACTGTTTAGGATAATCCTTGCTGCACACCCCTTCTTTCATGCATACTGATTTAGGATTCAGTACCCCACATGGCCCATGCACCATAGTAGATTTGATGATTTCAAACAGCTGTGGGTCTTCTACTGGGTCAGGGATTGTTGCTTTGATGAGGCTGTCAATATCTTCAGGCTGACGTAATTTGCTGTCTTGTTCCAGAATAATCAACATGTGGCAATGTGGTAGACCTCTCTTTTGGAATTCGATGACATGAACCCATGCTACAGGTTTGCCAAAAATGTGATTTTTAGTTATGTCTTTTATGAGTTGGTCAAGCTTCAGTTTGAAAACTCTCGAGATTAAATCAGGCCTGTCTGAACTTTGTTGGTGTGGAAGAAGGTTGCTCTTGATTTCTTCCCATTTAGGATTGCACGTAAACGTCAGAAACAGGTCTGGCTTTCCAAAGTGAGCAACTATTGACATGGCATCCTGATAGTTCTGATTCATTGCACGTGGGGAACCTTGAAAGGATGAAGGCAATATCACAATTCTACCGAGTGTGCGTTCTGTTTCAATGGCTTGCTTATTGATGTGGTCCATTAAACCTTGGTAGTGTTCAACTCGCAATTGAGATTGATTATGTCGCAGATATGACAGGCGTTCAGATTCTGTTTTGACATAGGCGTCAACTATGTACTGTTGGAAGAGCTTGCCGGATTTGTGGATGGGAGAGAATTGTTGTCTTATGGCCAAACGAAATGAATAAAACTGGAGAAGAGTTACTGAATTCCTTTTTGCTGTTGCATGTGGTGCTGAATGTTGTATGCCAGGAATCCATCCACGTTCTCCATGTGGAAATAGGATTGGGTAAGTCATGGGGTCAACGTTGGGGGATATGGTGGAAATGTTCTGTAAAGGCATATCTCTGGGATAGACAATGAGGTCTCGATCGGGGGGGATGCCATCAGTTCCTACAAACACTGCAGCGACTTCATCGTGGGTTGGCAGGTTGTATCGCCTTTGATCTGTGCCTTTTTTGATTACCATGGTGACAGTGGGTGGCTGGCAATTTGATGTGATTGCTCTATTATTCTCTTCTGTCTCAACCTCATGCATGTGCTTGTATGCAGCAGCATAAGGGTTGATGGCAATTAGTTGATCATGAATGGAATCCATTATGCCCTCATGACAATGTTCATTTTGGGGCTGGGCCAGTCTCGTTTGCGTTGCTTCTGGTCCATCAAGGATATAGATTTGGCCATATCGTCGAGGTGTATTTGGTGGTGGGTGAAGGGTTCCAGTTCTATGATATAATTGCCCATGAATACGGAAACAGTATGGTCCATGGCCTGGTGGCGGAGCCATCATTGCGCCCATCGACGCAAACGCAAATGAACTGTTGTAGTTGCGAATATGCTTCAGGAAGTTTTTTGATTTGTGGTCTTGTCCCTCAAGAAGATGCTGCAGTAAAGTAGGGTATGGTGGTAGAGGCTGGAGGTCAACTTTCCCGTTGTGACAGCAGTTACTTGCTTCATTCTTGAATCGACTTGCATGGCAGTATTGGCATTCATGTGACATTGGTCCCAGAGAGAAACTGTTCAATGTTGAGTGATAGTTTTGAGCAGCATTCAAATTAGATCGATGTGCTACATTTCCTTGAAGTTCAGTGGAATCCTGAGCATCAGTTGGTAGTAGGTGAGCAGGAGAATCACTTTGTCTTTCTGTACAGTTTATGGAATCTTGAGCATCACTTGGTAGTAGGTGAGCAGGAGAATCATCTTGTCTTTCTGTACAGTTTGTGGAATCTTGAGCATCACTTGGTAGTAGGTGAGCAGGAAAATCATCTTGTCTTTCTGTACAGTTTGTGGAATCTTGAGCATCACTTGGTAGTAGGTGAGCAGGAGAATCACCTTGTCTTTCTGTACAGTTTGTGGAATCTTGAGCATCAGTTGCTGGTAGGTGAGCAGGAGAATCACCTTGTCGTTCTGTAAACTCTGTGAAATCCTGAGCATTAGTTGGTGGTAGTTGAGCCGTAGAGTTACTTTGTCTTTCTCTCGATTGTTGTTTTAAGAATTCTCTTCGTTGTTTTTTTTTCTCTGAGTTCCTTCTACCGTCAACAGTTTTGGCTTTTCTTTTCGGCATCTGAAAAAAAGAAAAAAGAAAAGAAAACATTTTTAAAATTTAAATATTCATCCCTTTAAAATGATTTTCAATGAAGCAAAGGGCAAATGTTATGTAAGGAGAAATTGAAGTATGTTGCATTATGAGTAGATAAGTAGTAAGTGTGTATGTAGCATGTATCTGTATAGGTCTACGACAGGGTGTAGAGAGAATTTCACGCCACCATTCTTTTTCTCCAACCGTCCCACCAACTCGTAATAGGAGACTCCTGCCCCCCCCCCCCCCTCACTCCCGCTCATGGCCTGGCAGACGATGTAGAATCTCTCGATTTATGATAAAAACGTATCATACACAAATCTTTTTTGCTCAATTGCTCAAGGGGGGGGGGGGGCATTTTTCTTCTCCATCATTTGGAAACATGTAAAAAAAAAAGTAGTGGAGCGCCACTCCTCTATGTAGAGTGAGGACGGGTACAACTTCGCCTTTATGAGTTCAGGCTTCAATTCAAACAAAAGAACACAATAATTTGTTTATATTTTAAAGGGCCACTCTAACCCAATTGTAAACATTAATAGGCTAGAATCTACACGTCTTCTTAGTAATATTTCAGCCTTTAAATTGACCAAAAAATGATGATAACCATTTATTTTATGCAGATCTCATGCTTGCTGGGCTACGGCTCTCAAAAAACTGTAATTCTGCAGACCTACAGTACAGTCGGCAGAAGACATGGTAGACTCTACACCGGGCCCTTTCACTTTAGTTTCAGTTTAAAAGATCTGCAATTAATTAATATCTAGATTTACTGGATAAGAATTAATTAATATAATTATGCTTACCTCTAATCTAATCGTAGACAGGATGCCAATGAGTTCATAGCATACACTGGGTCCGAAACTCCCGTGAGAAAAATAAGCCTAAAATCACATGTAGTCCACAGACAATTTTATGGCTTCATAATTTCATGCATCTAAGCTCCTTCTACCGTCAACAGTTTTGGCTTTTCTTTTCGGCATCTGAAAAAAATAAAAAAGAAAAGAAAACTTTTTAAAAATTTGAATATTCATCCCTTTAAAATGATTTTCAATGAAGCAAAGGGCAAATGTTATGCAAGGAGAAATTGAAGTATGTTGCATTATTAGTAGATAAGTAGTAAGTGTGTATGTAGCATGTATCTGTATAGGTCTACGACAGGGTGTAGAGAGAATTTCACGCCACCATTCTTTTTCTCCAACCGTCCCACCAACTCGTAATAGGAGACTCCTGCCCCCCCTTCACTCCCGCTCATGTCCTGGCAGACGATGTAGAATCTCTCGATTTACGATAAAAAAACGTATCATACGTCTACACAAATCTTTTTTGCTCAAGGGGGGGGGGGGGGCATTTTTCTTCTCCATCATTTGGAAACATGTAAAGAAAAGTAGTGGAGCGCCACTCCTCTATGTAGAGTGAGGACGGGTACAACTTCGCCTTTATGAGTTCAGGCTTCAATTCAAACAAAAGAACACAATAATTTGTTTATATTTTAAAGGGCCACTCTAACCCAATTGTAAACATTAATAATAATAGGCTAGAATCTAATCTTAGTAATATTTCAGCCTCTAAATTGACCAAAAAATGATAACCACTTATTTTATGCAGATCTCATGCTTGCTGGGCTACGGCCCTCAAAAAACTGTATTTCTGCAGACCTACAGTAGGTCTGGCCCTCAAAAAACTGTAATTCTGCAGACCTAAGTTACAGTAGGCCTAGGTCTGCAGTCAGCAGTAGACATGGTCATGGTAGACTCTAACTCCGGGAAGCTCCAACACCGGCGCCGCAGGCGTCGGACGCCTGCGGCGCCGGTGTTGGAGCTTCCCGGAGTTAGGTAGACTCTAACTTACACCGGGCCCTTTCACTTTGGTTTCAGTTTCAAAGAACTGCGATTAGTCTAATTGATATCTAGATTTACTGGATAAGAATTAATTAATATAATTATGCTTACCTCTAATCTAATCGTAGACAGGATGCCAATGAGTTCATAGTCATACACTGGGTCCGAAACTCCCGTGAGAAAAATAAGCCTAAAATCACATGTAGTCCACAGACAATTTTATGGCTTCATCATGCATGTGCGTCGCATACGTTATACACGTACACACATGTCCGCTTGGACTACAAAGCATATCCTATAGAAAGAAGCATAGCACGGGTAAGCTCAACCCGTTGAGGATACATATTGGCTACATCCCCTCACCCGGGTGAGTCGTCTGTTGATGTGTATTTCTCTTACACGCTACAAAGCGAACCCTATAGCAAGAGGTCTTACTCGGGATACGCCGCAGCCTTTTTCAATGAGACTTGCTACACCCGTTGACGCGAGTGGATTCTCGTCCATCGTGTGTCTCTTCGTCTCTCTTTGATAATTAGACAGAGGGCGGGTGCACGCCTTGCCGACTCGCTCTCGCGGGTTCGTCGCCTCCTCGCAATCCGCGACGCCCGCCGTGCTCTCGTCCGCCCGTCCTCCCCGGCCGCCCCGAAACCGCGTTTAACGGCGTGGACGTCGATCGCGGGGTCGACGTGGCACGCGGGTTACGGGAGACAGACCTCGTGGAATTGTGCCGCCCTCCCCGCGCCACGGCGGAGGCGGGCGGGGTCGGGGAAAGACGTTTTTGAGGAACCGTTTTGTACCCGTTTGGACCCCCCGCGGCGATTTTGCGTATTAGTAATATAGACTAGCAGGACCCGTGGAAAATCCACGGGATTCTAGGGGAGGGGGCCAAATAGATAGATAGATAGAGGCCTATGTGCTAACCAAAAAAAAATTAAGATTATCTTTATATAGAGAATTTTTATAATATATGCTGTTTGACCATGAATGGCGCATCGTAAAACAAACTATGTCAGATAGGATCTCGATAATTACTATTTATTTGCATTATTATTGTTGTGATTTCAATGTAGTTTCATGAAAAGACTACATATGCAGAGGCTTGGATTATGATTATGTTGGCATTGTTCACAAAATGTGTATAATTATTTTCAAAGTAAAGATTTTACGCGAAAAATCGGACGGTTTGTGCACAAATTGACTCGACTGTGTGTGGGACTGTGAGCGTAGAATCCAGTGAAAGATTATTTTTTTCGGGGACTCGATCTATATCAAATTTTGGGGAAATTTTTGGGGGCTACTTTTTTTATCTTTCGGGGGCTCATTTTTCTTGGGCTAATTTTTCATTAATTCATTTAAAGGTGCAGTTTGTTTAGTTTAGTAGCGTAGTGTTACGTAAAAAAAGGCGTACTTTCTGCGTCAAGTTTCTACGCCGAAAGTACGTGAATTTTGTTAAACATATTAAATTGTTTGCGTCTATCTTGCGTACTTTTTTGCGTCCCTCAGGGTATAGTAACGTAGGCCAGCGTACTCTCAGCTCAGCGTAGTTGGTACGCAACTCTACTATGGACCCCACGCCAATTAGTGATCGACCTGCCAGTGTCACAAATAAGGCATTAAAAATATGCTCATCTGACCGTAATTAAATGCCTAAAACCCCCGAAAGCTCAAATCGCTATGAGTTACCTGACCCCCCCCCCCATACTTTAAAAATTGATTTGATGCAAGAACTGTATCATACTAAAAGTATTATCATGTTTTTAAACGACCATGACGACCGTCTGATTTTCTTGTGAAAACTTTACTTTGAAAAATGATTGTCGAAATTCAGTATACAACTATACATGGCCAGAATGGAAAATCCACGGGATTCTAGGGGAGGGGGCCAAATAGATAGATAGAGGCGGCCTATGTGCTATCCAAAAAAAAATTTAGATTATCTTTATATAGAGAATTTTTATAATATATGCTGTTTGACCATGAATGGCGCATCGTAAAACAAACTATGTCAGATAGGATCTCGATAATTACGATTTATTTGCATTATTATTGTTGTGATTTCAATGTAGTTTCATGAAAAGACTACATATGCAGAGGCTTGGATTATGATTATGTTGGCATTGTTCACAAAATGTGTATAATTATTTTCAAAGTAAAGATTTTACGCGAAAAATCGGACGGTTTGTGCACAAATTGACTCGACTGTTTGTGGGACTGTGAGCGCAGAATCCAGTGAAAGATTTTTTTTTTTCGGGGACTCGATCTATATCAAATTTTGGGAAATTTTTTGGGGGCTACTTTTTTATCTTTCGGGGGCTCATTTTTTATTAATTCATTTAAAGATGCAGTTTGTTTAGTTTAGTAGCGTAGTGTTACGTAAAAAAGGCGTACTTTCTGCGTCAAGTTTCTACGCCGAAAGTACGCGAATTTTGTTAAACATATCAAATTGTTTGCGTCCATCTTGCGTACTTTTTTGCGTCCCTCAGGGTATAGTAACGTAGGCCAGCGTACTCTCAGCTCAGCGTAGTTGGTACGCAACTCTACTATGGACCCCACGCCAATTAGTGATCGACCTGCCAGTGTCACAAATAAGGCATTAGAAATATGCTCATCTGACGACTGACCGTAATTAAATGCCTAAAACCCCTGAGCTTTCGGGGGTTCAGCCCTTTGTACCCGTGGACACCACGCCTGCACTCAAATCGCTTTGAGTTACCTGAGCCCCCCCCCCCCCCATACTTCAAAAATTGATTTGATGCAAGAACTGTATCATACTAAAAGTATTATGTTTTTAAACGACCATGAGGACCGTCTGATTTTCTTGTGAAAACTTTACTTTGCTTGTGAAAACTTTACTTTGAAAAATAATTATCGAAATTCAGTATACAACTATACATGGCCAGAATAATCACAATCTAAGCCCAGGTTTCATACATCTATTCACGAAACTATATTGAAGTCGGGGCAATAATAGCAATTTCAAAGTGTTGTTTAGACTCTCCTCACACCACGGATAGGCCTGGCTACTCCGGTCACCGCTACTCCATGATGCTCACACACTCAGCTCCGACAGCACCGCCGACAGCAGCGCTCAACTGCCGTGGCAATACACTGTGGACTCAAGGCAGTTTGCAACTGAACTGTCACGTGAAGGCGTACTATATCACAGATCTATTTGCGCTATTCTTGTTGTGATTTCGATGTAGTTTTCATGAAAAGACCAATAAAACAGAAGCTTGGACTATGATTATGTTGGCATAGATCACAGAATTTGTATAATTATTTTCAGAGTAAAGATTTTACGTGAAAATAGGTCCATGGTTTGCGCCAGTTGTCGCAATCGCCTCGACTCGAATCGACTGTGGTGTGTGCATGTGTGTGACCTGAGAGCGTGTAATACACGTACACTATGGGCGCCGGGCCTACAAAGCATATATCCTATAGAAAGACGCATAACACGGGTCAGGTCAACCCGTTGAGCGTACATATTGGCTACATCCCATCACCCGGGTGAGTCGTCTGTTGATGTGTATTTCTCTTACACGCTACAAAGCGAATCCTATAGCAAGAGGTCTTACTCGGGATACGCCGCAGCCTTTTACAATGAGACTTGCTACACCCGTTGACGCGAGTGGATTCTCGTCCATCGTGTGTCTCTTCGTCTCTCTTTGATAATTAGACAGAGGGCGGGTGCACGCCTTGCGGACTCGCTCTCGCGGGTTCGTCGCCTCCTCGCAATCCACGACGCCCGCCGTGCTCTCGTCCGCCCGTCCTCCCCGGCACCCCGAAACCGCGTTTAACGGCGTGGACGTCGATCGCGGGGTCGACGTGGCACGCGGGTTACGGGAGACAGGCCTCGTGGAATTGTGCCGCCCTCCCCGCGCCACGGCGGAGGCGGGCGGGGTCGGGGAAAGACGTTTTTGAGGAACCGTTTTGGCCCGTGTTGGACCCCCCGCGGCGATTTTGCGTATTAGTAATATAGATTGTCTTTTCCTACTGGCCCTACAATTCTTATTACCCCATAATTGCGAAGTAGGCATGTACACCCGGGTGGAGTGGGACAAAATTGCGGGCAAAACACCTTGTCCAAGGATGTAAGGACCACGCAAGAATTATACCTGGGAGGCTGGGACCATGGAGTTACTCCATGCTGGGACCCTCGTATAAGAGCCAAGAGCCTTCCCCAATTAACCAAAGTGCTCCCACTACGGATAGGCGCATCATATGTTTTTCCGTGTATCTCATTTTCATTATATTTTCAATATGCTATTGTAATTTTTTAGCATTGTAATTAAATTAATCAATTTAATCAATTTCTGCATTAAGAGTGTTCTATTCAAATTGTTTTTCTGTTTTACAATTGAATGTCCTTAATCAGGAACTTGATGAACAACAGCTTTTGCTGATCTATTAAAAAGAAGAAAAAAATACTGTTTCAATAAATAAATGGATAAATGAAATGAGAGAGTGATTTCATTACCCATTGTCGAGGAAGTTGTAACTCTTTGCGTCATTCGGGTACTGAGATGGCGTGGCCCAGCAGTTCTCGAAAAAGACAGTCAAGCCAGCCACAGTTGTCAGATTAACCGCAAAGTAAAGATGTTCTCCAAGAATGATCTCGCCTTCGTCGTCACCTGCGAACTCGTCGAAGGTCGCATCTGTGTAGCGATCGAGTCGAAGAGAAAAGTCACCGTAGCCCGTTTCGTTGAACACTATCTCGCCGATCTTCGAGTCGAAGGAACTACCCAGGAGACGATGGCGATCAAGTTCACATTTGACAGGGATCATAAGGATGTTATCCCGCGTTATTTCCGTGCCTGGCAACGGCCGAGGTTTGTAGGACGTTACCACATTAGAATAGATAATGCTGGTCTCAGTTTCCTTAATAATAATGATGATGATAATAATAATGATAATGATGATGACGATGACAATGATGATAACAATAACAACGATAATAATAGAATAAGTAATGATAATACTGATAACAATAACAATAATAACAGCCTCTGGTGCTACAACTAATTATCACAATGATTTTGTTCATACAAAAACAATCCGCAGCATCAGTGCTATTCTTCCCGAGGGAGGTCCATTTACCATTAACTCAGTAATTGCCGGGTAGCCATTAATATGTTGACCGGGATTCGGATCCCCAACCCTTGTATTGAGAGCCAGCGAATAAACCACTACCACACAGGGCGTATCAATAAAAAGTTTACCCAACTTTGGTAGCTCCTAATTAAAAATTCATAACCAATATGACTTGGACATTGCTTTCATGATGAGGAGTAACTTATTTTCTGTTCAATGATACCCTTTTAATAATAATTTTTGCCACGGGTAAATGAGCAGCAGGGTCTTTAAAAAGTGTTGTCAACATCGCCTTTTTCTAAGATTGGGCTTTAGGAGCAAATGGGTCTTACCATAATTATGCATGTGCTTTTCTGGTCTTAAATCCTTTTTTAAGTTAAATTGTGCCCTGTGCTTCAGTGCTTCATTGTTCAAAATAATGCAATATAACTAGGTTTTCACAGAACCTTGTTTCAATTGACAGGATAACCCTTAATTATCTTAAACAGCCTACTACCTTCTTCAAGGCTTGTCCCGTAATATCTCAGTTCTCACCGCTCACGTGTCCACTGCGCCATTTCACACGCGTCAGGAGCATAAGGTGTTTAACTCGAATAGCGATTTTGATTTGGGTTACTAAACCCAAGTGCAAGAAAGAACATTTATGATGATTACATTTCAAAGAAAGTGATAACAGACTTTTATTTTTTAAGAAGGTAAAATCCTTGATAAACTTTATCCATTTTGTATCACTCAATGAGAAGGATATATTAAACAAGGGTTGCGGATGAAGGCAAGATTGCAATCATCATAACGAAGAGAATATCATAAATTTTGCTTGTAATAAGCATTTCTTTGAGTTTGAGTTGATGTTAACTATTTTTCTTCCAACTCCCTAATAGTCAATGTAAATAAATGAGCAATTGTTTTGGAATGACATAACACACACATAATGGCTTTTGACCCACATTCATTGTTCAGCCCAAACTTCGTATAGACGCGTACACGAGGATTTCGACAATACTTTTTGAAAAGCCTCCTGTTCAGTCATCGTTCGTAGCACAATTATGAAAATGGTATAATTAGACAGGGTATTAGTTTTTTGTACAATGAAATATAATTCAAGATCACACAGGTGAAAATTTTTTAATTAGAAGCTACCAAAGTTGTGTAAACTTTTATTTTGATACGCACTGTATACCTCAACACCCCGAATCACGCTACATTGAAACAGCGGTAAAATTATGCATTATGATGAAGATCGACCTCAAAAAAGAGTGTTACGAAAAGATTAAGTATGTCTCAGTTTTCTTTTTAACATGGTTAATATAGCGGTAAAATCATTTGGAAACACGACTCGAACAAAATACACAGCTCTATCTCCACTCTATAGGAACAATGGAATTGATATAATGATACATGTCAAACATGCTTTTGGTTTATATAGGCCTAGAAAATTATTCTAATTCATAGAATTTTATCAAAGGTGTGTTGTACGTCACCACGTTAGAATAGATGATGCTGGTCTCAGTTTCCTAAATAATCATAATAATAATAATAACAATAACAACAACGAAATAATAATAATAATAATGATAATGATGATGATGATGAAAATGATGATAACAATAACAAAAATAATAATAGTAGTAATGATGATAATGATAACAATAACAATAATAATAGCCTAGAGAATTATTCTGCTTCAAGGAATTTTATCAAACAGTAATAGTGCAACTGTTTAACGGCATATGACGTTTAGATCTTCTCACTGGTTCTAGGTGGGAAAGGCTCTATTAAAAGCACAGAAAAGAAGGAATCATGAATGCTTCGAAGCTCAAAGAATGATAGGCCAGTTTAACAGATTCACCTCAGTTCTGATTGGACAACGAAAATGACGGTACAAATTGTGCTTATGTTATTTACTTGCACTCTGGTTGAACACTTGTCGTATCTGGTGTATAGAACGATTTGAGTAGTATCATGATCGATCCCCACGCAATCTTCATCCTCGAAATGAATATCTTCTGCATTCGTTCCGTTTGGTAGCTGCCTTATGTCCAGGAAAACCGCCATGGATGTGCTGTCACACTCCAGGTTAACTTCTGTTAAAGATAATCATTTCCAGATGGTTTTCAAATATCGAGTGTACGCCGCATGCTCTTGCGACGATTAAACCTCGTCTATTTCCTATATAAACTAACATAATTTCTTAACTTGAACCCTGGACCCAACAATATTACTAATCCTAACCGTACACCTCACGTAAAACCCAATAGCAACCTTAACCATATTTGACGAAGTAAAGCAATTTCGTGAATACTATCTGCGCATGGGCAGATAAAGGTCATTTGTGATAAAATTGTAAAACCCACCGAGATGAATACGTGCTCATGTTATGTAATCATTAGGGACTTTTAGATTTGCATGGACTGATACGTGTGACGTCCATTCAAGGCAGTGCATTAAAAGGACGTCGATCGTAGCAATCCCCAATGTCTTCGTGTAAATCTAAAAGTCCCTAATATTCATGTGTACACCTACGGTACATTCTATGATACAAGATCAAGAGGGCAAACAAAAGCTTTTCATAAAGCACAAACATGACTTTGGCATTGATTATTCTTATACACACACAAAGTAGGACTTCTCCATGAGTAAAAATTCTGTGGATATTTGGGGTCATTTTACCCATGTGCTAACTACGAAATGGTCTATTGCCGCAGGAGCAAATGTCTGGTCACCAAAATAATATCAAAGAAATAAACTCTCAAGTGCACACAACCGTATGACGAACCTAAGCTACACAGGAAAGGTCTTCACCGATTCATGAATTTTGAAGTATAAAGCCGAAATTGATAGAATAAACAAGTTTGTGCAAATCCTACTCAACTCTTGCAAAATGGTGTTTTTTCAATGTTCGTTATCTGGCTCCAAACAGACAATGTAGTATAAAGTCCGATTAGAGAAAAAAAATTGTGTTTAAAAATTCACAGATTGCGAAATATCCTCAAAACAAATATATATAATATTCCTGGACAACTACCATTTGCAACTTCAGCATGGGGAGGGGGGTGGGCTAGGGTAGGCATCGACAGGTTTTCTAGGGGCTGGGAGGGGCGAGACTAATATCATTAATAATGCTTTTTTTAATTTTTATTCACAAAGCATTATTTTGCCTGAACTTGAGGTGGAGTGTGACCAATGCTCCGAGCCACTCACCAGTTGTGCCATACGATGTAGGTTGCGATGTAGCCGTTGCATAGTGATCATCATCGTGGGATATAGTTTAGATGATAACAAAGGGGGGAAAGAGAGAGAGAGAGAGAGAGAGAGAAAGAGAGAGATAGAGAAGGATGAGAGAGAAGGGGGAGACGGGAGTTAATGCAAGTTAAAATACATGTAAACAAACAATTACACCTTTGATCGTTTAACTTAAACGTAATAACTCGTATCACACAGAGGTATTCACATTGTGAATTTGAAGGATGTCACATTTATATTATACTTAAGTTAGCTTTTATGATATTTCTGTGAATGTTTTACTTAAAATGATGGTCTGAATATTCAATTAAAAATACGAAGTCCGACAAATTCAGAACTCAACATTTCCCAATACGCATTGCAGTTGAAATTCCACCTATTTAAAGGGATAGTCGAGTTCTGTTGCAGAGGGCGCTACTCAAAATATTTTTTTGGTGAGATAATGAAAAACTGACTGGCAAATTTGAAAGAACATACAATTTATGATAAATACGATTGGTCAGGTCGATCGTGCCTTTCATGCAGTGGGTACTTCGCAAAGCTGCGGAACAAGTAGTTTACTCCGCAACAAGGCTATCGATGTTCAGGAAAAAATTGCTTTTGAGAATAGCAGGTTGAATGAAATAAAAGTGTCTTTTGCAGGTGAAACTATCTGGTTTATTCATTTATATGCACCAAATGAAGTAAGCGACAGAATGGATTTTTTCCACGAATGTTTTCATTTATTTCATGAACATGTACTTCCAGGAGAAGGGATTGTTTCAGCCGGTGATTATAACACGGTCCTGAATTGTAAAAAGGACAAAAAAGGAGGAATAATCAAAGATGGTAATGATGCTTTGTATTTGCGTAACCAAATCTCTTTAAAAAAAAACTTAATTGATATATGGAGAGTGAGAAACCCTGAAGCAAGGGAGTGTACACATGAACAATCTTATCCATGCTTTGTACAGTCTCGAATAGACTTCGTACTGTTATCCGAATGGTTGAATAATTTTGCTTAAGTCATTATAGGTCGCTTTATGCAGACCATAATGCAGTGTATTTAAAATTCGAGTTTAAGGCCGCAAAAAAAGGGCAAATGCATTTGGAAACTGAATAACTCACTTTTGCAAGATGATGGTTATCAGAAAAATGTAATTTATCTTATTAATGACTTAGAAGCGAAATTACAATTAAAAACAGAGAATATATCAGTTTTGTGGGACATTTTCAAACAGAAAATAAGACACATGTCAATTAATGATTGTATTAATAAGGCAAGGGTGAAGAATTGTCTCATGTTAGAAGACAATTGTCAGAATTCTATGAATACGAAATACACGGTGCACAAATACGGTCAAAAATCAAATGGTATGAAGATGGTGGGCGAAATAATAAATAATTTTTAGGCTTAGAAAAAAGAAACTTTAAGTTAGAGAATATGTTTAAAGTAAATGTTAATGTAATAATATACAGAAATATGAAAGATATAATTGTAGCTGAAGTTAAATATTTTTCTCGTTTGTATTCGGCGACAGGAAAAGAGAATGAAGTGAATGATTATCTGAAAAACATTAAGGTTGATAATGTATTAAATGAAAAACAAAATAGATACTGTGATTTAGAAATAACACTTACAGAATGTAATGACGCCCTAAATGGGACTAAAAAAGAAAAAAAAACAGGTTCAGAGGGTTTGACTTTTGAGTTTAACCAATGTTTTTGGAATAATTTGCAACATATTCTGTTAAATTTATTCATAACATCATTGTTAATAGGAAATTTACCTATTTCACAAACTAGAGCTGTCATAACCCTTTTACATAAGAAAGGTGATAAAGAAAAATCTTTCAAATTGGAGGCCGATTAGTCTCCTAAATACCGATTATAAAAGCATCTGTGTTAGCAGGAATGTTGAAAAAAGTATTAAATGTTTGTATTGATTCAGATCAAAAGGCATATCTTAAAGGACGTTTTATTGGGAACAATATTCGTTTGATCGATGATCATATATGTTATTTTTAAGAAAAAGAGAATTCCGGAGCTGTTATCTTTTTGGACATAGAAAAGGTTTACGATCAAATAAGTTGGGATTTTCTTTATAAAATTTAGACAAATTTAACTTTGATTGCATATTTATTGGACATATAAAAACTCTTAAAAAAATGTTCAAAATTGTATTATTAATAATAATTGGCAATCATCTTTTTTTTCTCTGCAAAGAGGTCTTCGGCAAGGATGCCCAGTGTCTGCTCTATTGTTTTTAGTTGTCATCGAAATATTAGGAACCGAGATTAGATCAGATGAAGCTATTCAAGGTGTTCAATTACCGGGATTAGGAGATATGCCTGATACAATTTTCGTTGCAAACGAAATTCTTTTAAAGATGCCATGCGTAAAATAAAAACTTTTGGTTGTGCAGCAGGACCTACTTTAAATGTATCAAAATATAAGAGGTTTGAAAATGCCAGATTTTGTCTTGCAAAGTTTTGCTTTAAAAGTAAATTGGCTCTTTAAGTTAATGGATAATAAAGATCATAATTTTAGATGAAAGATTCTCCCATATTATTATTAAAACTATCTTGGGCCTAACTTTTTAGTGCTCAAGATGAGTTTCTTAAAAAACGAAGATTTAACAGTTCTAAAAGAAACACCTTTATTTTACCAACAATGTATTTCCGCATGGAAAAAATGTCAATGTGTAAATAAATCATGTAATGAGGTCAATAATATTAGAATGCAGCTGTCATGGGGAAACACATATATAAAAGAAAAAAAAATTATGTTTATGGTGGAAACACTGGCTAAAATCTGGAATATTATTTATTGGGGACATCATAGACGCCAAAGGGGATGTGGATGTAAATCGGGCAAAAGATTTACTTTGTATTAAAAATAATTATTTATGTGAGGTTAAATTACTTTGGAAGTGTTTACCTCGTGACTGGAAAAGAAAATTAAAATGTGAAACTGTTTTGTGTTCATATAAAAATTAAAAAAAGTTATCCTGGGAAAATATATGGAAAACGAATTTGACATACGTTAAGGAAAAAAGATTAAAATAGTTTAATTTCGAACTACTCTATAAGTTATTAGTATTTGAAGAATTTCTATATCGTTCAAAAAGATCTGATCATGATTATTGTTTGAAGTGCAAATGCACAGAAAAAAATATTATCATTTATTTTATGACTGTGAGTACTCAAAAATATTTTGGAACAAGGTCTTAATTTTCTCAGAAAAACAGAAGTATTGTTTGATAATGAAAATATTGAATACAGTGATATAATTGTAGGCAAGCAATGTGATAACATAAATTCAAAATATGTATTTTAAAAAAATCTATTAATTAATTTAGCTTCTTATGTTTTATATACGCAAAAAAGAAAAGGTAATTTTACAAACAAATTCGACATTTTTTTAAAGTAAATAATTACTATGAAGGGTAAATGTTGAAAAGCATCATCTCTGGAAAGAAATGTTGCCTATGTTAAACAAAAAAGTATAATTCTTTTTGTAATATAAGTATAATTCATGTGTTTTTCTTTGTAATCGAGATAATTTCGTGTCTCTTATGTACGTATTATTATTCAGTTAATAAATGATGGGAAAAAATACGAATTAAAACATACAATTTAAATGACAAAATCATCTTTGTTTGATAAAAATTGGTTCAGAGACGGCTGTAATATCCATAAATATAACGATTTCTAATACTTGTGTATCACCACAAATATTAGATAACTTGAAAAAACCGTTGAAAAAACGAAACAGCGCCGTTTACGTCATATTTCAGCCAATGCCGATGTTGGCTTTTCGAGTTGGTTATATTACGCACAACCAGCCGGCCCATGTATTGTTATTTTTAATTAAATAATTAATTAACGGTCAACGTGTTTTTTTTTAAAGTGGAGCTGAATTGGTATGATCGTATTGGCTAAATTTCACACATTTTGGGGTCTAAAATTACCTCCATAGCAATTTTAGGTTCTAAAAGGGCACTAGAAATTTCATAGATAAGAATAAGAGCAATAGCCATGCACTTTGAAGGCATTTGCGATCTGAAAACAGAACTCAACTATCCCTTTAACACAGAGTCATGAGTCATTAATAAATCGAGATAGCAAAGAGGTTAGAAAAAATGACAAGGCGATGCACATTGGCATGGGTGTCGATCTGTAATTGATACCTATTGTTCAGGGGGATGATCTATTCATTCATCCCCCTGAACAATAGGTATTACGTAAGGATCGACACCCATGCACATTGGTATTGTTGGTGTTGAATAAATAACTGAAAGACAGAGGCATTTTTAGAGTAAGATGAGGCCCTATTCATGGAAATGTCATAAAAGATACTAAGAAATAACTTTGAAATCTTTTAACCTTACAATATTCAAAACTCAGTTTTTTACGAGTTTAGGTTTTCAGCAATCATTAATAATTGGCAATTTCAAAACATTATTATCACTATCATTATCATTATTATCATAATAATCACTACAAAACATATCATCATCATCTTCATTATTTTGTTTAGCGTAGGCCTCCCTTTTGCGTCAGGTTTCTTCACATTTCATTGTTTCTCTTTATCAGTTATATAAATGATGAAGAAGACAAAGAAATTGATTCGATTTTTCTTTCTTCATAAAAAAACTATTTAAAAAAATAAATAAAAAAAGGATCTTTAATAACTTATTTCATAATTTAACTGAGGTCTGAAATATCTACGATATTACTTATATTTCATTTTTTTCTTTCTGTAATTTAGGCAGCACGGATAAACTATAAAGCCGAATAAAGGCTGTTCATGGCCGCGGTTTCCATTCTTTTTTTTGCGTACAGCGCTTGAAATATCCTATTTTAATCTTACAACAGGCCTACACGGATAATAGGTAGGCCTACGTGAATCCCCTCTAGACATTCTTTTAAAATTCTTTATTCATTGACGGGACTGTTTTTCTCATAATTTTCTATTTTTATGATTTTAGTTATGTGACTAAGAAGGTCGTTCAGCTACCTTCAAATAAACAGCTTATTGCGAAACTAGCTCCTTTGCCTGAGGTTCCGCGATTGTTGTAGTAGTAGACCGTCATGACGTTGCTGCTTGAATTCACACTGTCAGGTATGCTTAAATTACCAGTTGAATGTTGTATCACAGAATAACCAGTATTCGGGCCATCATACACGCGAAGATACCCGTTAGAATCAAGCCAGAAGTCGTTAAAATCGATTCGAACTGTTGTGTAGCTGGGACATTCTACGAACCACGTGCACTTGTAGTAGTCGCCATTTACATCCGGGTAATCAGGGGTTGTGAGCGTGATACGACTATCTTGGCCAGACAATGATACAGCGCTACCACATGAAGACATCGGTCCAACCTGCTCCCAAGTTGGTTCTCCTGAATATAACATACACCAAATATATGATAGATAATATATTAATAGTTGCAATAAATATAGCGTTTGGTACTTATGTTTCTAAGTGCTTTATAAGAATATTATTACACGGCATGCACCCCCTTCATCAGAATCTGCCCGGTTCTTTATAATTCCTTACACACTTTATTCAAGCCGTCGAAGATAGCGTTATCAACAAACTCACAAGCCCAGCCTGAGGAGGACCTATTTTTTTCGCCGGGCGTCCTTTTTTTTTTCATCACGCCATCATCGGCTGGAAGCAAGTAAGTAGGGGGGTAGCGTTGGTAACTAGCCTGGACTAGCCTACTTACGGCAAAACCACTTTGATAAGGCGCTAAAGCCATTTCCTTTACATGGGTGGAAAATGGCAAATGTAGATGATTGTCTTGCCAAAGGAATTTGCAAGGCCAGTCTTGTGATCCACACTCCGGAGACATACTCACTGGACCACGACATATCTTGACATGCATATTAATGTAGCTGCTTAGTAGTAAGATAGGAGAAGGCAAACAAAAACATAATTTTAAATAATGTGTCTGTTATAATTTAGATCATCCCATAATGCGGATCCGTTTTCATAATCGCAAAGTCCTGGGGGTGTTTCACAAAGATCCTAAGTAGAACTTAACACTAAGTTTGACTTAAAAATTATGGACAGCCATTGGAGCCTTAAAATGAAATAACCTCAAAAGCAGAAATTCAAGTATTTTAGACTTTTTATTAATTCCACATCCAGAACATGAAATCTCCTTGATGATAGAAAATATGAATGATTTCAATAAACAATGAATGCCTTTCAAAATAATGAATGCCTTTCAAAATAACTTTTACAATATATTTTCAAAGATGCTAATGGCTTTCCATAATTTTTAAGTCCAACTTAGAGATAAGTTCAACTTAGGATCTTTGTGAAACACCCCCAGGATATAAATTTGTCGTGCCTATTGATATCGCATTATCTTAGTTACCCAATCACATCCTCAGTCTATCAGGATTTGGATGAAGTGATTTTTAGAATACTTTGAAAGCAATGTCTAAAATAATTTTGGAATACTATATATTTATTTCATTTGTCTGTGTGTATTTGTGGAGGGGGGGGGGGGGAGTTCAGCCACTAGTGATCAAAATTACTGCATTTTAAATTACTGCAATTTTGGAATACTATATATTTATTTCATTTGTCTGTGTGTATTTGTGTAGGGGGGGGGAGTTCAGCCACTAGTGATCCAAATTACTGCATTTTAAAATTAAGCAAGTTTTGGTTTAACTAAAAATCTTAGAAACGAATTTATCAAGCGTCAGTGTGAAACCGAGTCATCCCTGAAATACAAACTTGGCGCTCGTACTACGTATTAATATTCATAACGCCAGTATGGATTATCAGAAATGATTCTTTTAGTCTAACTTCGAAAAACCTCGGTTGAGCCGAACCGCTGTTCTTCGAAAGGATCAACTCTTGGCTCGTATTCGGTTACTTCTATCACGGATATACTATACTGTCGCTAGTAGATCCATACTAATATTGTCCACTGAAATAAAAGAGCTATCGAGGGAGAGAAAGAAGGCAAAGAAGAAAAGAGAGGGGGAGAGGGGGAGACAGAGAGATAAAGAGAGAGGGAGAAAGAGAGGGGAGACAGAGAGAGAGCAGAGAGAAGGTATTAGCTGTAATATGGAAATTGTATTGTCTAGATATATAATTGATAATTTATGTACCTCCGTTTATTAAACATACCTATGCATTTGTATTACCATTATATATCTTGTGAGCCAAAATGACCTTTTTTAACGTGAAACTTGGAAATAACCCCGAACAGTACTGAACTGAATATATCAAACCTGACGAAACAATGGAGGCTTCAGCGAGTATTCCACGTTGCTTTGTTACGTCACTCCAGTCAGTATCATTGTAGTAAAGCGTCATGGTATTGCTGGATGACTGAATGGCTCCAGAACCTGGTGTGGAACTCGACGATATTTCGGCCAACTGATAGGTTGAAGAAACAGAAGAGGCATCGTAAACTCTGAGGCGGGAGCCAGAGAAAATGTCGAAATCGTTGAAAAGCAGCGAGATTCGGTATCCGTATGGGCAGGTGAAATACCAGTAGCAGATGTACCGGTCATTCCTAACCTCGACTGGATAAGTCGGTGTTCTAAAGACAATAGATGGTGATGAAAGGTCGAGGTGTACCGTTTGACCACAGGAATCAATACCCACACTATCCCACGACAGTTCGTCTGAAAAACAAGGTAAACACACAATATTATAGGTATATAGTTCAAAGATAACAAAGTTGAATATTGTTATCAATAGAGGAAAGAACGCTCGTTGAAGATAATTTGCGTTGGATTTTACATAAATTACGATCTTCGCAATTCCGATGAAAACCGACTCTGATTTTCCGTCAGATTTGCATTACTTTAAGCTTAACTCTCGACTTAAGTGCTTATATGAAATCAAAACATGTTAAGAGTCAAGTCCAAAACCACTCAGGTCATACCTGTTAATGCGACGTCCATGGCATAACCATGACCATTTCTACTACTAGATGTATCGTAGAAAACTCCAGAGAAACCTCCTTTTGTTGAAACCACCGTAGAATTACTGACAATACCTCGAAGATCGGAAATCTGGAAGAGCGATGACGATGGATAGTCGTACACCGTCACATACTCGTTTTCTGTCTGGAGGTCGTGTACACGCATCTCAACTTGCTTACCCGAGGGACTTCTGACGTACCACCAGCAGTCGATGGTCACATCAGATTCCGAATTTTTATGATACAAAGGAGTGGTGGTTGTGACCCGTGGTGTACTTGTACTCAAAGGAATCACTCCGCCACAGAAACCCGGCACCTCCGAGGGCATTGCATCTGAGGGCGACAATAAAGAACAAAATGCGGAGCTGTTCTTCGAATTGTCATGAGGGCGCACGGAAGAGTGATGGAAAGTAAAATTAATTGATGTTATGATGCTAATTGTTTATTTTTAGCTGCGTGTGTATGCTAAGGCTGGGGAGGTATTTCTTAGGATTTCTTACCACCTGTCTTCTCAAAGTATGAAGCATATACATGTATGCTACCCTCTTTCCTATGTGTTCATCCTTTACATTACATTAATTTGAAATTTTACCTATATACAGTTAAGAATTACTACTCATTCAGTGGTAGACCCACTTGTGTGTATTATGTTCGTTCATATGTCACATACAAGAGCATATAGAGATATAAAGTCATATAAATACTATTGTAAAAAATATTGTACAATTTAGCGCTTGGGCTGCCTTCTATATTTCTTAATAAAACAATTATCATTAGTATGATATTAGAACGAAAGCGTAATACATATTTTTGTTAATAAGAACTCAAACTAACGCCTTTTTTACGAACCATACGGGCCACAGTCCGATCACGTATCCACTACGTTATGACACATCTATATTGATTATAACTATTATGGGAAAAGCATATTTTTTAGTAGCCTAGAGCCCTAGCCTGCACTAATCAATCACCCTTCATCAAGTTGAATACACAAAAAAGGGGAAAAGGTGAATTGCAAAGGAAGACTTATATTATTTAAAAGCAGAGGGATTCCATCCCATCCACATGTTTTTGATGAAAATACCATTGAAACACCTGTAGGCACCGTAAAGGGTACAAGAGATGTCAGGAGGAGGAAACATATTTGCTTCAACTAACAACCACGCATATCTTTCTTGTATTATCAGTTTAACTTTTTCTACCAGCTCATCTAGAATGTAACATGTTATCAACGTTGGATTTCTCGTTCATACCTAGAATCCTCGCTTCGAACAGGATTCCTGGGCCGTTATAACCAGAGTTGTAATCGTAGTAACGCACTATCACTCCGCCCCGTGAAGACAAGACAGTAGCAGGTGATTCCAGGCGTCCAGACAAGTAAGATACCTGGTAACTTGGACTAATAGAAGTCGGATGGTCATAAATATATATCATTTCATAACTGCCTGTCTCTAGATCAAGTACTGTCAGTTGGACCTGTTTTCCTGAAGGACTTGTAATATACCAGTAGCAGTTCAACAAAACAGAAGATCTAGACGACTGATAGTAGTACGGGATGGAGAAGGTAGCATGTGGTGAAGACTGAGTCAGTGGAATCGTACCTCCGCATGAAGCCGCGCTCGAAGAATTGGGTGTGTCTAGGAAATGAAAATATTTATAACATACATGTGTGTAAAATTATCTCGATAGATAGACGGCAACTTGCTTACATTGCTAGGTCGGGTAAACATTATGGATTGCCCTACCTTGCAAACACGTAGGCCTAATATATTCATATGCACTCTGACTTTTGAAATAATTCGTAGCATAGGCGCCGGAAGTTGGGGGCGGGGGGGGGGCACTTGCGCCCCCCCCCCCCCACCCTTAAAAAAAATATTTAAAAAAGCAGTAAAGGAGTAAAGATGGTTTGTATTTTTGAAACAAAAAGAACATACAAATAAGAAAACCTGAATTGAATCACATACAAATTATATTCCATTCAAAGAATGGGATTTAAATTTCCTATACATTCCGCATGCCCCTATATTTATCAATCGAGGTTTTTATAATAGCGATACCCGCATATCGTAGAAATATTGTGAGAGGATGGGAAAGCGAGTGTGTTTTTTATTATAGAGGTGAGAAGACATGAGCCCCTGAACACTGACGTACTGATGGGGGGGGGGGGGTATTTCTAATTTGAAACCAGTGATTGATAAACGTTGCTGAAATGACATGATTTTCCGGTGTTTTTCGGGGTAAATCGTACCTTTTTATGCGTTTGTTGGAAAATTACAAATTCCCAAATTTACAGTCGGGCGAAAAAAAATGCTTATTTAGAAGCATTCATGTCAGGGCATGCAGAATGTCGGGGCTTGGGGGCTCTAACCCTCCTTCACGAATTTGAAAAAGTATTTAAAAGAGAAAAAGAAAGGAAGAGAGAGAAAAAAAAATGAGTGACCCGCCAATGTCACAAAGGAGGTATCAACATATGCTAATCTGGCCGTAATTAAATGCCTAAAACTCATGAGTTTACGGGGGTTCCGCCCCCTGGACCCCCCCCCCTCCCCAAAGGGTCCTAAGGCGGCCCCTGGATGCTCGTCTGACCGTAATTAAATGCCTAAAACCCATGAGCTTACGGGCGCTCTGCCCCTTGAACCCCACCTAAGGCCCCATGGCGGGCCCTTGGATGCACATCTGATCGTAATTAAATACCGAAAACCTCTGAGCTTCCGGGGGCTCTGCCCCATGGGCCCCTGGATGCTTATCTGACCGTAATTAAATACCTAAAACCCTTGAGCTTCAGAGGGCTCTGCCCCTGGACCCCGCGCCAATTAGTGAACTGCCAGTGTCACAAATAAGGTGTCACAAATAATCTAAATATGCTTATCTGACCGTAATATAATGCCTAAAAACGCTGAGCTTCCGGGGGCAAAGCCCCTGGACCCTACACCAATTTACGCTCGCTACACTCTACCTGATCGTTAATACTCCACGTTGAAAAACGTTCCGCGGGGCCTCCATGTATCATCTAAGTGTTGTGGTGATGCAGATTTTTCTTCTACCCCCCCCCCCAAAAAAAATGGAATGTTTGCCCCACCAATCCCAAACTGCTTCCAGCGCACCTGATTCGTAGTACTTGAATTAATTAATTCGTAACGAGGCATTAACAAAGATGAAGGCCTATGTATTCAACTTGATTCCGTTACAAGCATGGCTTACTTTGTCATAAATGTGAATTTCCCACCGTCAGCGCATGGGGAGGAGGGGAGGGAGGGGGCCGTGGGATTTGAGTGGGATTCGAGTTGGAAAATAAGTGCAGACAGAGACATTTGGGGGAATCTATTTTTACTTCTACCCCCCCCCCCAACCCCGTCGGGAAAAACCCGACATCAACATTATTTCCATACCATCCAAATCAATGCCATCAAATTCCTCGACAGAAATAAATCTCTGAAGTAGTCAATGATGCATTTCTTTGGATGGGGGAAGGGAAGAATGTACAAAGTAGCCGATGACAGTAAAGAGTTCTTGTTGAAAATTCGTAGGCGTTACAAATGGAGAGCCTCCTACTTTACCTCTTACCCTGATAGCAGACTGACGTCTTAAAGACGTCATAAATTGGTCGTAAAGTGGTAAACACGTCGTACGTCGTATGGACGACGTAATTATGACGTCATTATGACGTCTTCAAAAGACTATTTTTTCCAAGATAATTTTTGTCGCATTGACGACGTCATAACGACGTCGTAATGACGTCTTTAAAGGACTAATTTTTTCAAGATATTTTTTGTCTTCTCGACGACGTCATAATGACGTCTTTAAAAGACTCTTTTTTTCCAAGATAATTTTTGTCGTATTAACGACGTCATAAAGACGTCATTTTGACGTCTTTAAAAGACTATTTTTCCAAGCTAATTTTTGTCGTATTGACGACGTGATAACGACGTCATAATGACGTCTTTAAAGGACTAATTTTTTCAAGATATTTTTTGTCTTCTCGACGACGTCATTATGACGTCTTTAAAAGACTATTTTTTCCAAGATAATTTTTGTCGTATTAACGACGTCATAAAGACGTCATTTTGACGTCTTTAAAAGACTATTTTTTCCAAGCTAATTTTTGTCGTATCGACGACGTCATTATGACGTCTTTAAAAGACTAATTTTTCCAATGTAATGTTGTCGTATTGACGACGTGCTAACGACGTCATTATGACGACTTTAAAAGACTTATTTTTCCAATGTAATGTTGTCGTATTGAGGACGTGATAACGACGTCATTATGACGTCTTTAAAAGACTATTTTTTCCAAGCTAATTTTTGTTGTATCGACGACGTCATTATGACGTCTTTAAAAGACTAATTTTTCCAATGTAATGTTGTCGTTTGACGACGTGCTCACGACGTCATTATGACGACTTTAAAAGACTTATTTTTCCAATGTAATGTTGTCGTATTGACGACGTGCTAACGACGTCATTATGACGACTTTAAAAGACTTATTTTTCCAATGTAATGTTGTCGTATTGACGACGTGATAACGACGTCATTATGACGTCTTTAAAAGACTAATTTTTTCCAATGTAATTTTTGTCGTATTGACGACGTCATTACGACGTCATTATGACCTCTTTAAAATACTATTTTTTCCAAGATACTTGTTGTCGTATCGACGACGTCATTACGACGTCTTTAAAAGACTATTTTTTCCAATATATTTTTTGTCGTATTGACGACGTCATTACGACGTCATTATGACCTCTTTTAAAGACTATTTTTCCAAGATATTTTTTGTCGTAAGGACGACGTTATAACGACGTCATTATGACGTCTCTAATAACCTAATCATTTCCAAGATAATAATCGGGTAAAGACTACAACATCTAATATGAGCATAAGTCAAAGTAAATAGGGCAATTATCAGAATTGAAATTAGGGAACAACTCCCCCCCCCCCCCTCAAATGTCATTAGTCATTATAGATTGTGCAACATCTTAACACGCACCAAGGTCTGCAATATCGATCAAGACTAAACACGCTCGCATCTAGTTCAAACGGGAGTGTGTACGTGTGTGTACAGATGCATAGTACTCAAGGCATGCAAGGCGTAGGACTGACTACATGCACTATGCGATGTGTATGTAATGTATTTGAATTCAGAGTAAAAAAGCCGTTGCGTCAAGTATCATCCTCCTGAATACTCTCGTTTTCTCCTGCAAATTTACAGTATGGGTAAGTCCTAGAATTACCTAATTTATTTACAAAAAGTACATAAATCTATAGTTACTGTTTTTTGTTCTGAAAAGTTGATATGTAGACTCAATGTATGATCTTATGTCACTATAATGTCAGTGATAGTGAGCTCGGACTTTCAACACTTGGTGACAAGATTTGGAATGTATGATCCCATATATGTATTTTTATGTAGGGTCTTGATTGAAGGTATGCATGGCATGGGGCATCATTTAAAGCAAAACAGTTTACTCATTTGCCTACACATGGGCTGCATGGCTTTAAAGGGAAATCTAACTCGAATATACAACTGCTCTAAAAGACTAAAACTAAAAGAAAGCAGAAAAATAAGAGAAGCAAACTGATGAAAGCTTTGGCACTAGATCTAGGTAGACTCGTCTAGATTAGACGGCAATAACCATCTGATTTGGGTCATTTCTTTAATAAAACTGACATTTGTTGTCATTCCCCAATACAGATGGCAACCCAGAGATCTCGTACCGTGGGGGTGGCTATGTACATATTTTCGTCCACGACATTTTCTCCCCGATTTTCTTTGTTTAGATCTAATTTCCTGATTTCTAAATCAGATGTAAATTAACGCCATTACGAACGTGTAGATGGCGCACAGTTTTTCCTTCATAATCATTTTGCCTTCGGCGTAGACTAGAAAGGATGAGATTGATTTAGATGACTCCGGCAGCGAGAGTTTAACACGGCTTGCAAGATTTTCAACATGGATGTCGATTTTGACCCATAATTTCCTTAAAAGGGACCATTAAACGGATTATTTCAATCAAGTTAGTAAGTGGAACCGATGTATTATAGCCTAGGCCTAATTTACAATTGGATATGGACGAGGTAAAATGACATTTTGGGGGATTTTGGGGGCATCACCCTACGACCTGAAATGGCCGTGGCGTTATCAAGGGAATGTGTATAGACGGTTTCCGTCTCTAAAGAGAGAGCCGGCCGACCCTACCCCAATAAGAGGCTAACGAACGAAAGAGACGGAAAATTCGGTCTAGGCCTACCCCAAATCAAACAAATATTCATGAAATGTTCAAATAATGCCGGTCATGATTTAAAAGTCGTTATGTTTCTAACCTGTTTCTGGTCAATTTTCTTATTTTGGGGGGAGTATTTTTTCGGAAAAACATTGAACGTAGTATACACTTTCTCACAGGTTGCATGGCAAACCGCTGCACAATACAAGTGCATACAGATGGAGGGTTGCGGGGGCACGTAAGTCTGCGTGAAGTGATCTACCTCCTCATTGTACGGTAGATAACTTCACACAGGCGTGCGCCCGCCCGCGACCACCATCGGTATGCACTTGTAATGTGCAGCGGTTTGCAGTGAAAGTGAAACCTGTGAGAAAGTGTACTACCGGCCCGGGCCCAGATGTCTGGAAATTGAGCATTACGTGTCGATACGTAGAGGAATAAAGCAGTGACTTGTGTCACAAACATGTAACCTTAAACAAGTCGTTTCTATAGTCTTTAGACAGCAAATTTATCTTATAGGCTACCGCGTATACAGAAATCTTTATCTGGAAATAATCGGTTAGACCTGTAAATGGTAAACAAAATTCTTTCCCCACAGGCTGTTCAAAGAAGTTTGTAGTATATGTGTTGCACATGCATACATTATCAACCATTCATTTTGTTCTTTATTCATTTTCTCATTTATTTGTAGGCATCCTTTTTTAACATCTTTCGTCAAAGAAGTGGATTCAAAGCACATTGTGGCCTCAGTATATTGATGTGAAACCATGCAGAAAATGCTACTCTCTCTCAACCTTGTGGATGTTCCCTACTTCCCTAAGGATCAGCAGCCTCTTGAAACAGTTGAAGCCATTATCGGTAAATGTGTACATGTACTTGTATTATTTCTGTACATTTATTTGTCATATTGTGCATGATCACTGAGACACTGACCAATTAGACACTACATAATGGTTCATTTCATTTGAAAGGCTTATTTCAGTTATTTGTCATGAGTTTCAAATATATTTATAATTCAATTAATTTTAGATCACTAAAACCTGCAGGTATGTTGCAATAACCGACAATATTATGTGAGTTTCAAACAGTGTTTTTTTACCCTTCAGCATAGGTCAAGGAGCTTTCATTCTACCCTCATCAGCAAACTACAGGAATGCCAGTTTTCAGGCAATAGCCTGAATTCAGGCCCTGGCGAGTTATTTTCATACCATAGTTTAGGTTTTTTGTGTACAAGCTAGCACATTCTAGGTGATTTTTCAGGCCCTCTGATCAATTCTAACTGGCATCCCTGAAACTACAATTGTATTGATAATTTACCAAACAAAATCCTGAATCATGTTTTGCATTGCCCACAAATGTTTAAAATCCAAAGGTGTGTTGACTTTATGTTTGTGAGCAATGAGAAATTATTTCTGATAAATCTTACCAAGCAGTACTTTTCTGACATAATTCATTATTTAATTCATTGTTTACGTTCCAGTTTTAAAATATATAAAAGTTTGTCTTAAATTGTATTCCTTTCTTTTTGCTTTATTCATGCACTAGAATTCCTTGCGTCTATGTTTTCTGGATCTGAAAGCAAAGCATAATGCCAGAGATTCCACCCAAAAGCTAACCCATGAGTCCAAGCTGAAGAGGAATATTTCCAAATCGATGGCTCCACTAATCCGATTGGTGAGTCATGGCCCTCATCAGCAAACTACAATTGTATTGATAATTTACCAAACGAAATCCTGAATCACGTTTTGCATTGCCCACAAATGTTTAAAATCCAAAGGTGTGTTGACTTTATGTTTGTGAGCAATAAGAAATTATTTCTGATAAATCTTACCAAGCAGTACTTTTCTGACATAATTCATTATTTAATTCATTGTTTACGTTCCAGTTTTAAAATATATAAAAGTTTGTCTCAAATTGTATTCCTTTCTTTTTGCTTTATTCATGCACTAGAATTCCTTGCGTCTATGTTCCCTGGATCTGAAAGCAGTGCATAATGCCAGAGATTCCACCCAAAAGCTAACCCATGAGTCCAAGCTGAAGAGGAATATGTCCAAAATTACATGTAGATCTACACTCACCCGATAGTTAAACTGCTAGAAAAAAGGTTTTTTTTTCACTTCATTATAATTTATAGAGAGAGATTGGGGGAGGGGGAGAGGTGAGGTGTAGCCACAAACCTGAATGGCTAAGCTGTGGAACTGAAATTCCCATTCATGTCACTTACAATTTTAATGAATTACAGGGTGTACGAACTGACATAATAACAGTTTGAGTAAAAGATTCCAGTATTTTGCATATATTTGTAATTACTATTTATGAGGTTGACTTGTTTCATAATGCTTAAATGTAGTAGACTTGATATGCTCAGTCTCATGCTCACCCCCACCCCCCATCCGCTTCTTATCAGTTCTTAAAACCTGTTGCAAATTATATATAAGTTGGTAATAGAGAAAAATTGAGTGAGTTAAATACAATACAAATCGATTCCATTATTAAGCAACTAAGACTTGTTCAGATACGGAACGTGAACCGCCAACAATTTTCGAAAATGATGCAGTTTGGTTTTATGTGTAAAGCGAATGGCAACAAGTCTTTGTATGTTTTGGACTGCGCATATTCAGTGAGCTCTTTTTTAAAAAATAATATACAGAAATGTAAAGCATATTGAACTCTATGTCAACATTAAAAAAACTCTGAGATTTTGTTTTGCAGCATCATATCTGAACGAGCCTTAAAGGGACTTGGTAACATTGTTCATATGAGAGAAAATGATGCTAAGGTATTGACTTCAATGCTATACTGTATGTACACTTAAAGATAATCTGAACATTTGAGTTAAAAATAGACCACTATTGATGATTAAAATTGGGAAAACTGATCTAATCCATTCCTTAATTCTCCCAACTATATTTTTATAATTTGATTAAATTGAAGTGAGTGTAAATCATCACAAAGCAACGAATTAGCTCCAATACCTTTTCCAAATTTAGGCGTTGAAAGATTTTTAAAGGCCTGTTTTATATCTTATTTTCATGGAAAAACGGTTGGGGAAAATTGAAAATTTTGTTGGGAGAATTAACTGTTAAAGTGTGACACACAGTAGCGCCTTCAATGAAAACCATAAATACCTTGTGGTTTGAAAAAAAATTGGTAGTTAAAAGTATAAAATTCTTGGATTGTATCTGTGATTTTGAAAGCACAATATGGCTTACCTACATGTATATTGAATTTGCTGTTAAAATTTAAGTTAAATTTTAATTCTTTTGTGAAACTCCTGGGACTTGCCAAATGAATGTTGTGACATTAGCAATGGAATATATCTCTTTGTGCTAAACTTTTTCAGCATGTTTGCTATTTTATCTTAAACTGGGAAATATCTTCTGAATAATTATTCCACCGTCACCTTGGAATTGCCCAACTGAGCCAATCGTACTGGCGGGAATGTTATCTGATATAAATTGTATTTTATACTGGAACTTCAGCTTCTGTATTTGCAACAAAACTAATTTTTTGATATTTATAAAGATGAATGTTGCCACTAATGCCTTTTCTTATATTTAAGTCATGTGCCCGTAAGAAGCAGGTTTCTAGCACAGAGAATAGGTACAGTTGAGCTGATTTCTGCACAGTTTTCTAATATGCTGTATCGCCAAATGTGAACATCATGGTTTTGACGGGCAAGCAACCATCATCAGAAAACTCTGCAGGTATCAGTGAAGCTCTCCCTGTTCCCCTTGCTGCCAACTTGTTCCTTTTGAGTGTAAGAAATAACTTAGCTAAGCCATGTGTTTATCATTTGTGCCATTTTTCTCCTGTGTGTGATTTAGATTGTATAGTTTTATACTTTATACTGAATTCTTGCCATGAACTTAGATATTCCCCTTGTATGTAACAAGAATTGATCTTTCATGCTCAGAAATTCTTTCCAACAACCTTTTGTTTTATGTGATAGTCTTTATCTGCATTAAAGTCTTGGGGTTGTCAGGAGCGGTTTTGCAGCATGGAGACGAGGTAGAGCTGCACTGATGAATGCCTAGCTTTCCACTATGCTGTTTTGCTGAATATGAGGAACTTCCCCAGTGTGCTGTTTGCTTGTCTACTAAAATACTGCGGGTATCAGTATAGCTTTACCTAGTCTCTATGCTGCTAACCTGCTTCCCACAACCGTAAGACTTATAAGTTTAGCCAATATGTCATAGTTTCTTGTAGCACAATTTATAACGAAACAGATATCTTTGTGATGGCCTGTTAACTTTCTTAAAATAAATTATGAGACTGACAGAAATACACTAGGTATTTATACAGCATGTTTTGGAAGAAGTGTCAGGATAGTGATGAAAAGATTTTGTACCAATGGAAAACAAAGTTAGATAGAATAAGAAGGTATAGAAATACGTATTTTATATTTTTGTTATCCATTTTCATACCATATACGTGTACATTAAAAATATGAAAATTACTTTCTGTTTGGCATGTATATATATATTATGTGATGCAATATTTGTTGCTTTGATTTTACTGTGTGATTGTTCAGAATGTGTACATGAAGTGACTTTTACACTTTCTTCCTTGATCAGGTGGTCGTATAGATGTTATTTAGAGGACAACAAAGACGTAAAGATCATGTAAGGATAAAAGGTACATCGAAATGACCAGTAAGGTGTAAATATCACGTCAAAATTTGATGTCCATACAACGTATACAGGTCTTTCTTGGGTGAAAGTTGACGTCTTAAAGACGTCTTTATTTAACCAATAAGACCTAAAAAGAAAATATGACGTCCATACAACGTCTTTATATGGTACAAGTTGACGTCTTAAAGACGTCTTTCTTTGACTAATAAGACGTTAAAATCACGACAAAATATGACGTCCATACGACGTCTTCCTATGGTACAAGTTGACGTCTTACAGACGTCTTTCTTTGACTAAAAGACGTAAATATGACGTCCATACAACGTCTTTCTATGGCACAAGTTGACGTCTTAAAGACGTCTTTCTTTAACTAATAAGACGTAAAAATCACGTCAAAATATGACGTCCTCGTGACGACGTGCGCGGGTGCAAAATCAAGTCTTTTTGACGTCGTGTTATGTCGTCTTAAAGACGTCTTTCTTTAACTAAAAGACGTAAATATGACGTCCATATGACGTCTTCCTATAGTACAACTTGACGTCTTAAAGACGTCTTTCTTTGACTAATAAGACGTAAAAATCACGTCAAAATATGACGTCCTCGTGACGACGTGCGCGGGTGCAAAATCAAGTCTTTTTGACGTCGTGTTATGTCGTCTTAAAGACGTCTTTCTTTAACTAAAAGACGTAAATATGACGTCCATATGACGTCTTCCTATAGTACAAGTTGACGTCTTAAAGACGTCTTTCTTTGACTAATAAGACGTAAAAATCACGTCAAAATATGACGTCCTCGTGACGACGTGCGCGGGTGCAAAATCATGTCTTTTTGACGTCGTGTTATGTCGTCTTAAAGACGTAGTACAAGACGTAAATATGACGTGATATTTACCAGTCTCTGCTATCAGGGTAGACTGAACTCTGCCCGGAAGCCTTTTCCGTACCTCCCACAGTACCTGTCAGTGTGTTGTACTGTCCATCCACCACGTGAAGAAACCACAGGTGATATATCCACGTCAATGTCTCCGTTTATACTTGAAACTATGTAGTTGGTAGAGGATGGGTAGTCATAGATAATAACTTGATAAGTACTATAGTCTATGTCCAACCATGTCAGCTCAACCTCCTGTCCCGAGGGAGCATGTACACACCAGTAGCAGCGTATATATGTTTTAGAGTCATATGACGAGGTTTGGTAGTAGTGTGGGGAGTTGAGTGTAAAGTCTGGCATCGATTTGCTCAGCCTAATTGTCCCCCCGCATGATCCTGGGCCCACTACCGAAGTTCCATCTAGATCAAGAAGAGAACCGTTTCGATGTACTATGTAATGTAAACAATCATAATACATAGTACATGTATACGAATAGCGACGTTGATTTACCAGGGAATGCGTAAGGCGTCGATCCATACAGTTCTACTATCTTTCTGGCTTTAGCCAAGCTGCCCTAAGACGCTAAAGTATCCAAAGAATAGATATATACCTGATAACAATATACTGCTCTTGGGTGGAGAGAGGCACTGTAGATTGATGAAGAAAAACTAATGCCTTCGTTTGTCTTTTTAATAAGCATGGTTAAACTTTGTTTAAAAGATGATGATTTGTGAATCCTTATCATAAACATTCTTATTATCATTGTCATTGTGGCCTTTGGTGTTTTCATCATATTATCATACTGTACAATCGAAATTTTTACCTGAAAGGGAGACGTTGGCTCGGAACCCTCTTCCATTAGTATTAGATGAGAAGTCTTCTAGCATAAACGTGATGCCTCCTCGAGTAGACACGAACGGTGACAGGTCTGGATCAAATTCACCTCGGAATAGAGATACCGGGTAACTTTGAGATGATGGGTGATCGTACACCCTCACGTACTCATCCGATGTATCAATATCTAGCCATTCAAGCTGCACCTGTTTACCTGCCGGACTAGTCACATACCAATAGCAATCTGCTTCTCCTTTAACCGACTGACCGTAGTAGTACGGGGAGTTGAAGAAAAACTTATCTTGAAATTGATCCAGTTCTACTCTACCGCCACATGTCTTTTGCCCCACAACTGGGTATTCATCTGCAAGGATATTTTGCATAACAAAATCAAGGGCGCCGCTACAAGTGGGGGGGGGGGGGGGGGGNNNNNNNNNNCCCCACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTAAAATGAAATAAATAAATAAATAAACGCCCCCCATTTAAAAAAAAAAACAGCAGTAAAGGAATAAAAATGGTTTGTACGTTTTTAATCAAACAGAAAATAAATGATAATGAAACCTAAATTTAGTCACATGAAAATGATATTCCATTTTATAAAAAGTAGAGTAAGCATGCCCCGATATTTATCAATGGAGGATTTTATAATAACCATACGCGCATATCCTGGAAATATCGTGAGAGGATGGGAAAGCGTGTGTGTTTTTATTATAGGGGTGAGAAGACATGAGCCCCTGAACACTGGCGTACTGATTGGGGGAGGGACAGTTGCCCCCATCCAATTAAGATTCAAAACGTGGTTCGATTCGTTACCAGAGAAATACCAATGTTGCTAAAATAACATAATTCTTAGGGATTTTTCCTCACTTTTTAATAGATTATAATTTTTATGCGTTTTCGGAAAATTACATATTCCCACTTAACAGTCAGTCAAAAATATGTTTATTCAGAAGCCTTCGTGTCAGGGCCTTCAGAACGGCGAGCTCGGGGAGGCCCCCCTCACTTTTTCTGGAAAAAGTATGTTAAAGAGAAAAAGAAAAAAAGAGAGAAAGAAAGAATGAGTGGCCTGCCATTATCACAAATGAGGTAATAAAATATGGTCATCTGACCGAAGTTAAATGCCTAAAACCCCCGAGCTTTCTGGGGCTCCGCGCCCCTGGACCCATACCGGGGCCCTAAGACGGGCTCCTGAATGCTCATCTGACTGTAATTTAATGCCTAAAACCTCTGAGCTTCCGGGGGCTCCGTCCATTGACCCTCGGCTCCGCCCCTTGCCCCCGCTGGGGGAAAGCCCCTGGACCCTACGCCAATTTATGCTCCTCTCGCTCCGCTCTCTACTCTCTACCTGATCGTTGATACTCCACTTTGAAATACGTTCCGCGGGGCCTGCATATATTATGTAACATGTAACAGTGGTTGTGATGCAGATTTTCTTCTCCCCCCCCCCCCCCCACGTTTTTGCTCATTACTCTCTTTAAAAAAAAGAGTGAATTTTCAAAATCACGGACGTCACTCCTCTGGAGAGTGAAAATGTAGTGAAAATAGAGTGAACATTCACTCTTTTCTTCTTAGAGAGACCATAAAATCAATCGATTTTCACTCCTAAACAAGTGATCACGGACGTCACTCCCAAGAAAGTGAAAATAGAGTGAAAATTCACTCTTTTTTTTAAAGAGAGTACGCCACTGCCCCTGAAAGGTGAGCCCCCTAGGAGGTTTGCCCCCCCCCCCCCCCAATCTCAAACTGCTTTTAGCGCGCCGGGTTCAGACCAAATGCTCGATGTGGTCTTCCAGAAATACGCTTCCTGTCACGGTTGTGTCTCCACAGAACGCTCGTATGAACGCCTCTCTAACGAGGTGAGTATGAATAGGGCCTCAGTTCTACAGGGAGACCGTTTCGTTGAATAATGTTGTTAATTTCATTATATTTTCACTCTTAAACAATGTTCACCTGTAATGCTTATTGTAAATTGAAGGCCGACTCCTCCATAGGTTCCACTGAATTGGAAGTATAGTGTTCCTGAGTTCTGAGAGAAAACCCTATCTTCAGGCATGATAGATGGATCCGAGAAACGAGATACTATTCCACTTGTGTTAGAACTATCATATACATAGAAATGTCCTGCTTTAAGCTTATGGAATTCGACCTTTATACGTCTTCCGACTGGACTAGTAACAAACCAGTAGTTGTTCATATTTTCAACATAATATGTAGAAGTGAAATACATAGACTGAACGGTCAGCGATGTCGTATTTTCATTCAGTTCAATGACCCCACCGGATGATGGGAGGACGCTGGCATCGTTATGTCCATCTATATAGTATGAAAAAGATACGAAACACGAAAGAGACATAAGAGATAAGAAATGTTGCATGTCGCTTGCTGCAATAATTCGTTGGTAATGCTCTTTCTTGCAAAACGGATTCTATCTAATCGTTTATATTGATATGAAATCTACATTAGCTTGGATGTAAAATGTATCAGCGCCGATTTGGGACAAAGTCATACACTTCCTCTTTTAGCATATTTTTAGGGTTGGCGCAAAATAGGCAGGGGAATAACTTTGGGCTACGGCGATGGTCGGCAATTATGTAAATTCGAGAGTTTGAAGTTTCAAGTAATCTCTTATTTCTTGTCGATATCTCAAATTTCTAAAAACCTAAATAGGTTCATAAATCCGAAAGAAAATTGGCGTGCTATGTTGTATTGGGCTTATTCTGATAACCGACTTTCTAAAACTGCTCATTCGCTTTGCTTGCTCATATCATAATCCAGTTTAATCGGCGGGGACCAGCTTCCTGTTATTGGAACATCATGGAACATGAATTAGAAAAGTATTTACCTATGTAGCTGACTTTCACGCCGAAATAGGAATCCGACACATCATTATTGTAGTAAGCGTTGAAATAGACACGGGCACTGGAATCAGATATAACTGAGCGTGGGGCTTGTTGGAGTCTTCCATAAGCATAGAAAATTTGTCTTCCTCTTGAATCATATATACTCAGGTAATCGTTGTCATCAATTAATTGCCAATCATAGAACTCCACCACAATCATCGTGTTGGGAGGACTATTGATGTACCAGTTACAGGAGGCAGACACACCAGATGGATAGGAGTTGGTGACGTTGAAGTGAGATATCGGAGAAGTCTGGTCGAGATTGATAGTTCCTCCGCACTCTCCTGGAATAGCTCCTTCGATTGACAACTCACCTGCATGGAAAGATGAGTGAATGGGCATCCATAGTCAAGGAGAAAACCACAAAGACGTGTTTGATAAAATTCTAGGCATCTAATATGGGATAAAAATTCGTGCCTGATTTTATACGATACAAAACATTACTCGCTATACAGGTCTTTTACAAAGACTTGGGCCCCGTCTCATAAAGAGTTGCGATCGATTCAAATCATTCGCAAGTTGCAGTCTCCTATGTACATAGCTGTGATTGATTTCAATTACAAAAAAGTTGCGATTAATATCAACTCTTTATAAGACGGAGCCCTGGTTCCTGTCTTATAAACAGTTGCGATTAAACTTGCGATTGATTTAGATCGATCGCAACTCTTTATAAGACGGGGCCCTTGGGAAAACCTAATTCATAATATTATGAATTTAATCAGATGTCTAAGTCAGTCGGCAATCAAAATAATTCAGATGATCGTCGTTATTGCATTACAAAGAGATCAATGCGGTGCAGTGTGCTACACGCAACCTTGTTGGATCCACATTATATTGACAAATTAAAGGTAGCTGAAGCTGAATTGACAAAAGACTTATATCCGTTTCGTCGACAAACAGCCATGCATGAAGTCTTTCGTCGAAGGCTGATTGTCAACGAAACGGAGTTAAGTCTTTTGTTAATTCCGTTTGTTGCTACGGCTGAAAGCCAGTGAAAACGCATTCGTTTAATAAAATTATGGAAATGCATGCACTACTTTAAACAAACAGAAATTTGTCATTCTTTACTCATTTTTACTAATTGAGGTGTCATATGAAAAGACAGAAACTCAACTTTTTAGAAATGTGATTTTCTCTTTGAAATTCAGATACCGTTCGCCAAATGACTTTTTGGTTTCCTCTCTTGAAATTCATTTTAATCACTCATGCGTGTAATATGCATCATATAAGTTGCATCAAAATGAAAGAGGAGGATATAAGCTTTCCAAATATGTGTCATTTATTAATTATAAATGTCATCAGGTTATGATAAAACATCGATAATTCTCGGTTTCGTTTTTTCTGGGACGCACTGTATAGTCGAGTAAAAATACTTGACAGGAGATGCAAAAGTCGTTGGTTTGATCCTCGGCACGGCAGTCATATCCTTAAGCAAGGTATTTCGTCTACATTGCTTCTCATCACCCAGATGTTAACGGATAGCCCTACCCAGTAGGATGCAAAGTCATTGTGGTTTGCTTAGCATCGTGAGCGTCTTTATTATGGCTTAATACGGAATACAACCTTTAAGCGCTTTACAGATATTATTACCCCTGTCGTCGGATTCATTCCCGCACAATGTGCACCATCCCCACTCCCTGGGGAGTATTCCTACCATATGGTCACTGTGAGGGCATAGGCAATTACTAACACGAATGCGTTTGTCATACTACCGGGTCCCCATTTCTCTCTTGGGTATTGAGTGGCAAATGTAGATTAATGTCTTGCAAAAATACTTTGCCAAACCACTCTACTACGACACTTCCAAATACGACAAATACCGAATAAACCCCAACAAAATTATATTATTCTGTTATTATTATTATTTTTTATTAGGATTAGGATTTCAACCATCCCTTTATAAACACATTGTAAGCAATATGACCCACAGTTTAGTACCTTGTAGAGTTGCTTGCATGGTCACCGAAGTCTCTTGGTAGATGTCGTCTAATGTTGTGGGTAGTTTACACGTGAGTCCACGACCACTTGATACGAAACTGTTCCCATGACGTCCACCGAGTTGGATGCCTTGGTCATTATCAGTTTCGATCCCGAGTCCGTCATAGACAATGAAGGGTTCAGTCATGGGCTCCAAAGCCACTATGTCCAGTCGAATGGTCGTACCGGGTGGACCGCGGAAATACCAAAGGCATTCGTTTGTGGCATCGTCTTCAGTCCCTCCGAGACGGTAACGCTCGTCTGTGATTTGGGCAGTGGGTTGGTCCGCGTCCAGAAGGATGACTTTCCCGCATGGAGTCGGGTCTCCACCATCTAGTGATGATCAAGAAATAACAACAAAATATACAAACTATCAGGAATAAAAAATTGTTCAATTTGTCTCCATATTTCGCTCTCTTTACAATAAAATCAATTATTTTCTCACCTCAGGGGATTCCATGAAAGATAACATTACATTTTTTAATCATTTTGATATTTAAAAGTAATATATGATTTCATTGAGACGGGCTTAAAATAATATTACTTTCTACATTGAACACGATAGTTGACGGAAACGAAAAAACAAATATTAACAAGAGGCACAATTTCTATTTCTGCTTTTGTTACCTGTGAAAACCTTTCGACACTACCTAGCATGCTAGGGAAGGGACTTCATTGTAGATCGTACACATCAAGAAATATGAGGGCATCATGATTATGAAAATATAAAATTCATCGTATACGCAAAATGTCGGGAAATGATGAACGTTATGAAACAGGTATTGAAAATCATATCGTATGTAAATTTAAAACTTATATATATATATATATATATATAACAGTGCGATAATTCTCAACTCTGACATGACAACGGGTGTTCCAGAATATGTAGGCGTGAAGAAAATGGTCTAGTAATAAACTTGTAAAGTACAGCTACGTCTCACGTTCCTAGACCATTTTCTTCACGCCTATATATATATATATATATATATACATAAAAGGTATTGGAAGCAATAGGTCAATAAATAGACGTACATAGTTTGGTAGAATACCTAATGGACGGTGACAAAATCATTGGAAATGTGCTTCAATTTAAAATGGGAGATTCGTACATAACTGCTAAGTTAAAATCGGTGCATATCGGAAATACGTATACCGTAACAAAACAGACTAATTTTACATCTACAGTCGTTTTATGACTCAAGATATATTTGAGACTAATGGCTCATTATATTAACTGCATTTACAAGTTGGCCTTCCCTGCTGGAAAGAACACACTCATTCACAGTGTGCGCCACAATGTAAAATTGGGGCGTTTTCACATTGTGGTCTCGTCGTTCCACAGTGTGATTCACAGTATGGATAACCACAAGTTCACAATGTGACCACACTTTGAATCTGTAGTTATCCACACTGTGAAACGACGCGACCACAACACTTATTCAGACTCATTCACACTCATTCACAGTGTGAATCACAGTGTGCGCCACAATGTTAAATTGGGGCGTTTTCACATCGTGGTCGCGTCGTTCCACAGTGCGATTCACAGTGTGGATAACCACAAGTTCACAATGTGACCACACTGTGAATCTGTAGTTATCCACACTGTGAAACGACGCGACCACAACACTTATTCAGACTCATTCACACTCATTCACAGTGTGAATCACAGTGTGCGCCACAATGTTAAATTGGGGCGTTTTCACATTGTGGTCGCGACCAGACCAGACCACAACACTCATTCACAGTGTGAAACACAGTGTGAATCACACTGTGAAAACGCCCCAATTTAACATTGTGAACACTCACACAGTGAAATGCCCACCACACTGTGTAACGCCCACAATGTTAAATTGGGGCGTTTTCACATTGTGGTCGCGTCGTTTCATACTGTGAACCCAACTGTGAAACGAAGCGACCACAATGTGAAACGCCCCAATTTAACATTGTAAACACTCACACAGTGAAATGCCCAGCACACTTTGAAACAGCAAACACACTGTGAAACGCCCACAATGTTACATTGGGGTATTTTCACATTGTGGTCGCGTCGTTTCAGTGTAGATATCCACACTGTGAACCACACTGTAAAACGAAGCAACCACAATGTGAAAACGCCCCAATTTAATATTGTGAACACTCACACAGTGAAATGCCCACCACACTGTGAAACACCAAACACACTGTGAAACGCCCACTACGTTACATTGTGACACATTATATTGTATGAGATTACCACACTGTGAAACGCCCACCACACTGTGATACAAGTTTAGTTTTTGAATTCTGTGTGAAATATTTCAAACTATGAAAATGTAGTCATGTAGAGCTTAGTACGTATAAATACGAAGGAGAAAAGGTTTATAATCAATGTCTCACTTATAAACTCAACTTAAATAGAAATATACCCCCGAAACATCTGAAAATAATACAGTAATCTTAAACATATACATGTTATAGCCTTTATGAGGTAACTTAAAATGTATTTCAAAAGCTTATATTGATGTACTTTTCGCTAAAAGGACTTATGAACATATACCTCTTAGGTATACAAAGATTGCAGTAATACTATACTTACAGATAAACAAATTTTACTCTGTGAATGCCAACATGATGTGAACCAGAATGAAAACAAAGACAAATGAATCCAAACAGTTTAAATTAGATAGTAGAGATACTTCTGAAAACAGAACTTGCTATCTACAGTGAAACCTCATTATAATGAGCACATTGGTGTCCACTAAATTGACCTCTTTATATCAGATATCTTATTGTACAAGGGTTTAAAACAATGAATTACAGAAAAACGTGGACTTGTGAAATAGAGTTATAGAGGTTCCACTGTATAACATCGAGGAATAAATTACAGTGACGGCAAAAATTGCACATATCATCAATATTTCCCTGAACTCGTACTGACACTCGTAAAACATAATAAAAGCACATTAGCATGTTACGTCAGTAATAAAAATAGTTTTAGAGCATCTCTGGCAATATCTCTTAAACATTTCTGTAGGTATATTAAATTTATGAAGATACATTAGTATGGCTAGTCACATTTTCCAAGTATAGAATATCCACACATGACATGTAACTTTGGCAATTGCAATGTATCAGAAAATCTTTGGCACATACAGGTTCTTCATAAACCTATTAAAATGATAAAAAGTTTCAAAATATCTCTGTCAATACATCTTGAAAAAAAGTTGAAAACAAATCTATTTTAGGCAGACACGTTTTACTGTTATGGCTACATAAATTTCTACAACTAATACTAATAAATACCTGAAGATGAGAAAAAAAACAGATTTTGCTCATCTTTTTACAAATAAATAACCATACACTCTATATACGGTGGCACAGAAAAAAGCAGAATGGCATACAAAAATATATGTACTTCCTAAACTACATAAGTAACTTATGATTTTTTGTTCTAGCGATAATATACTGATGTTAGGATATTATAACTAGATGGAACTTCCGACTGTGCAGTCAAAAATAATACCTCCACTATCACACAATTTGTTAAAGCCCGAATTATCATCCCCCAGATCTTGAGATATCACTCCAGTCCACCGTATGGTGTATCACGTACGTCAGGTCTGGTGTCGCGAATATACAGTGCGAATTTAAAACTTTCTATAATATTCTATGAGCAGTTAAAAATCCTTCTTCCCCTTCTAACGGAGATGGAAATTTTGTACCTCTTGCTTAAAAGCTTTAGCGAAAGTACAGTTGACAAATCCAGCAGGTAGCAAGTTCCAAAGCATGATGGAATCGGGGAAGAATGACTTCTGGTATACGAGAGTCATTGTAAATTGTATTTGGAATTACTGGGTATGACCATGGTAGGATGTTAACACCGATAGAAGGTACATTTTGAGAAACTCTGAGTGCAGGTTGACAGTTCCGTACACTATCACCATCTTAAACCGTGCCCTGCGTTCTCGGAGAGATGGCCATTTGAGATGTTGCAGCATTGGTGCGATGCTGCTGGTATGATGGTGATCTCCAGTGATAAACTGGGCGTTGCGGAGCTGTACCATTTCCAGCCTGTGTATGTTGTCCTGGGTGTGAGGATCCCGCACTGTGCACGCACGCATACTCCAGGATTGGGCAAACAAGAGCCAAGTAGCACAATGCTTTCGTCCCTCTTGGGCAAGGCCTAATGTTCCGCAGCAAAAATGCACTGACAGAGCTGACTTTCTTTGTTGTTGTCTTTATGTGATCTGTCCAATTCAGGTTGCTCTTGAGGTCAGGGGCGTACGCAGGATTTTTTCAAGGGGAGGGGGGGGGGGGGTTTAAATCGGGCCCAAAATTTCGCAACGACTCAGCTTTATTTTTTTGCTCCCGAAAATTTGAAAATCGCGCCGGCCGAAAAGTCGCTCAGTGGCGGGGGGTGGGGGGAAAGGCACCCTTTTTAGGTTAGCCATGTGCATAATTTTTTTTATTTTATTTTTTTATAATTAAAACAAAAACAAAAACAAAAACAAAAGGGGGTGGGGGTTTTTTTTTTTTTTTAGGGGGGTTTACATACATAAAAATACCAAAGGGGGGGGGGGGGTTCAACCCCCTAACCCCCCCCCCCCCCATATATTATAGTTTCATATTTTGGTCTACATATAATTAGTCAAAGTCTGAAATCAACTAAATTTAAAAAGAAAACTTTCAAAACGACATACATTCAGCAATTTATTATGCTTTTCTAATCAAAAACAATTTCCCAATAGCATTGTTAGTATTATGAGCCAAGGGGCACTGTCAGTATATGTTATAGTAAAATGCCACACCAATGATGGCATCCACCTTGTCTTGGTCAAGCAACTGCTTCCCTTATTCTCCTCCTTTAATCTAAATGCTGAAACTATGCTGTCTGCCATCTGCTCATAGGTAAAGAGCATACCAATGAGTGCAAGAACAGCACGTGACCAGGATTGCTGGTGTAGCTCCAGTATTCCCAGCTGTACTTGGGTCAAGTATGTCCGAGAACCTCATGTATGGAACTGTGGAACAGAAAAGAAATGTCAAGAATTACTGAAAAGAAAAGCACGTTAACAGATTTTAATAAAATAATCAAACATATCCAGCCTAAATTTGACGAACTCCGTACGCCACACTTCGTACTTCTACGCAGATCGAGTCTAAAAGTAGTAGTACGACACAGTTCACGACCAAACATTATTCATACTACTCCTAGCTCAAGGAAACGCTAGGAAGACACATTTTCATACATTTTTAGACATTATTATTGTATATATGTCATTTCTAGTGACCTACATGTAGGCCTAGACTATATATTATGTTTATGCATTTTCGTCCGACTGTTTTTTTTTGCCGACAGAAAACATGAAATAAATTGTACCAACATTTATTTCTTACCTGATTTATTTCTAGACCTAGGGGCCGATATACACTACTTAACAGTTTGTGTGTTTAACATGCTCAATTATAATGTCCAAAAATGTTTGAAAAAGTGTCTTCATAGCGTTTCCCTGAGCGAAGTACTTATAACAAAGATGATACTCACCATTCAATGACAGCCCGGCAGCAGCCCGGCAACACCCTCGCCTGATAGCTATCCTAGCCGTTGGGAAACTTCACCAGGCCATGTCCGCCATTCCAATTATGTTCTCATCCGTCAACATGTAGTTGTACATATTCGTCTGTTATCACGGAAACTTGATTGGCCTCGATCCAAAACACTAGTATTCCAGCGGCCATACTGTGAATTCTATTTTTTCACAATGTGAAGTATTTTTCACAGAGTGGCCACACTGAGAAATTTACTGTTCACACCGTGAAATGTATGTCGTAAGGTATTCACACTGTGAAAATCTATTTATCACAGTTTGGTCACACTGTGAAATCATTTTCACACTGTGAGCAGGTACATTCACACTGTGAAATCGAATTTTACAATGTGAATACTGTGTGGTTTACCACATTGTGACAGTGTCCACACTGTGAAAAGGCTCCGTTCACAGTGTGAACTGTGTGTTCACACTGTGAATGAGTGTGTTTTTTTATACAAGTGTTACTGTTAGTGGCCACTATAAGGCCCCTTAGTGATTATAAAAAAAACATATGCTCTTCCTGGCTTCGGGAATGTGTCAGTTTGCAACAGAAATGGAGAAAAACGCATTTTGTTGAAAATCATGTGATATCACGAACATTTATACATCACTTGTATATGGAATTGAAAATATGTCATGGTTCAATTTTATTTGGGGGATGGCGCCCCCCCCCCCTCCCCCTAGCCCCCTACCAGTACGCCAGTGCTGACAGGTGTGCTCTTCGACCAAAGATCCGTCAGACTTCGTTGTTTCAATCGCTCATCGTACGACTATATACGAACTTGTGTTGCTTACATTTTTCCACGGCTGTGTTTACTACTCATGCTGTGTGCCTCTTTTAATGTCTAGTAAGTGCCTCTTTTCAGGTACGATACGTGCGCCTTTTCACTTGTAATAAGTGCCCCCTTCACCTCTGATAATTGCCTTATTAACATTTGAGAGGTGCCCCTTTTCCTGATGGTTAGTTCCCCCATTAAAACAAAGATCCCCTTCTCCTTTTCATAAATGCCTCTTTTCAGCTCTGGTGAATGACCCCTTTCAGGTAAAAAGAAGTTGCTCCTTTTCATTTCAGATAAGACTAGTGCTTCTTTCACGTCCGATAATTGCCCCTATTCAAATATGAGAGATGCTCGTTTTACTTTGACAAGTGCCCCTTTACCTTTCTAATAAATGCGCTTTTCATTAAAACAGAGTGCCCCTTCTCCTTTTTATAACTACTCCTTTTCTTCTCTGATAAGTGCCCCTATTCAAATGTGATCGGTGCCCCTTTTCATTCAATAAATGCACATTTTCCTTTATTATAAGTGCCCCTTTTTCATTCTTATAAGTGCCCCTTTTCATATCTGGTCAGTGTCCCTTCTCAGGTGCGATTAGTGTCCCTTTTCACCTCTGCTAATTTCCCCTTTTTCAATGATACGTGTCCCTTTTACTCTCTAATATTTTTCATTGCATATTTGTTGTGATATAAATATTCTTGTATTATGTTTTGCATGAGAAATAAAATCATTTGAATCAAACCTAGCCCCTTTTTACTTTTTATGACAAATGTCCATTTTTCCATTCGCATTAATGCCTTTTGTCATCCCTGGTGAGTGCCCTTTTTAAGGTGTAGTAGTGCCCCATTTTCAATTGTGTTAATTACCCTTAATAATTTGACACGTGCTACCTTTCTTTTCTAATAAGTGTCTAATCCCTATTGATAAGTGTCACCTTTTCCTCTCATCAGTAAACGTGCCCCTAATTTAAAGCACTTTCCCCCATAGCACGCTCGCTACGCTCGCTCACATTTATATATGGTATAATTATGTTCCTTATCAAAAATCTTGTTAGTCTTCCGTTATGAAAAAGTGTAATAACCCCATTAAGATTGGGGTTTTTCTGCTTTATAAGATGTAATTTTTTGGACTGCGACCTCACAAATCTTCAATTTTTACTCGCTCGCAACGCTTGCAAATGAATATTCAGACTATCACAGCCAACTGCAGAGTATGATGCTTAATGGGGCTTTCACATATACACTCCGGGGCCCGTACGGCGAGTAAAAACAGACATTTTGGGAATTTTCGGCCGGTGTCCGGCCGGTTGCCGTAAAGTTCGTACGGGCGCCGGGACAGTATTTTTAACCGTACGGGCGCCGTCGTGCTGAGTTGGCTGCTCCCAATTCTCCTGGCGTTTTTATCTGACCGTTAAGCCGTAGAAAGCCCGTAAAGGGCCAGTAGGGCACCTTACGGGTCCCATGCGGGTTCCGTACGGTGCCCGAACGGGCCCCTCACGATTTCGGGGGTCAATTTTCATGAATAATCTTTACGGGCTCCTTGCGGTGCTCAAGGGTTCCTTACGTGTACCCTACGGCGCCCTCACGGGCCCCTCAAGAGCCCTGGGTGATATGCACAACAACAGGCGTCTGCCGTCCGGGGGCCATTTACAATGAATATTGCCCGTGCGGGTCCCTCACGGCGGCATCTCGAAATGTACAACGCAGTTTTAAATACAGGAGGTGCGGATACCAGAGGACACATACAGAGTGAGTCAGAAAAAATGTCCCACTTTTTAAAAGTTGTATTGTTTGAAGTATGAAGAGTTAATCATGAGTTTCATGATATGTCTTGATAGAGTTATCCTCCCGGTACATTTAGGTACCATTTTCATTTTAATCTCTGCATGCATGACCGAATATCAGACAATCTTTCCAAGACTGCCAGATCTCACTTGCGCCAAGTTGCTAGGTTAGGAATAAAATAATGTTTTGAACAAAGACAGCCAGTTTTTTCAATTATCCAATTAGCCAAACAGACAAAAACACACGCACACACGCACACACGCACACACTCACGCACCCACCCACACACACACACACACACACACCCACACACACGCATACACAACACACACGGTGAGAATTGTTTCATAATTACCCCTTTTCTGATTGAGCATTGATTTTCATCAAGATCGCTGGCAAAACCCTAAAAATCTTTCCCTTTTGTACACATTCACTACCACCATCATGTCATCATCATCATCATTGTCATCATCTTAATCTTCATCATCATCATCATTATCATCATCATCATCATCATCATCATAAAGCACAACATCCGCACCAAGAGGACAGCGGCCACTCAACTTTTCGTACGACAAGAAGAATTCCCGACAAAAAGTAACAGTTTGGGTGGGCCTTGACGGGAATAGAATGCTCATTTGTGTGTGTGTGCGTGCGTGTATGCGTGTGTGTACGTGTGCGCGTGTCTGACTGTCTTTGTTTAATGCACTGTTTTATTCCTCACCCAGTAACTTGGCGCAAGTGAGATCTGGCAGTCTTGTAAAGATTGCCTTATATTCAGTGATGCATTCAAAGATTCAAATGAAAATGGTACCTAAATGTACTAGGAGGTTAACTCTATCAAGATATATCATGAGACCTATGATTAAATATTCATACTTCAAACAATACAACTTTTAAAAAGTGGGACATTTTTTCTGACTCACTCTGTAGTTCAATGCACCGCACGGCACCCGGCCAGGCCCCGGCAGGGCCCCGACAGGGCACCTTCCAGGAGCCATCGGTACACCTGCTCATCGGACGAAGCCCGTAAGAATACCGGCCGGCACTCGGCCGGCACCCTGCCGGCACCCTGCCGGCACCCTGCCGGCACCCTTACGGGCGCCTTGCATTTCACTCGCCCAATTTTCTCCGGGCACCTCCCGGGCGCCGCAGAATATGTGACCACACATGCCGTAGAAATGAAAATCGGACGGTTCCCGTACGGTCTCTAAATCCACGGCAAGCCGCCGGGTAGGGAAAATCGTACGGGGCCCGGTTTATATGTGACCTAGGCATAATAGCATTACCTTGTATAAAGAGAAACAACTATTGAAATCACCCCTGATAATGTAATCACTAATCTAGGGGTATATACCGATGTCTAAAAAATTGTATATTTTCGGTTTCCCTTTGAGGTGCCCTTTTTTCTCCAGAAAGTACCCTTTTGCGTTGGTCCCCTCCCCCCCCCCCCCACTTTGTAATTCGGCAATTATACAATAAGTTCAAAACTGATCATTTATCCAGTCGGTTGATTCACTGTACCTTATTGAAGTAAAATAAATTACAATTGTTAACATAATTCATTATGAAATGTAGGCCAAAATTTAAAAAAAAAAGAATTACGTAAATCTGAGAAACGTCTGAGCACAGAACAATAAAATTTCCAGATATAAAGACGCTGTCATGACAACTCACCTATTTCAGGCAGTAGAAGGAAAGCGACTGAAATGCAAAGCCAGGCACTTGACATCATAGCTGAAGAATTTCGCGCCATTTCGACGTCGCTGTCTGATCGCTTGTAGATAGATAAACTGTGGATCCTGTTTTACTCTCAGAACAATCTTAACTCTCCTCGAAATAATGATAGGAGTGTGCTAATTTCTCTGCTCTTATAATTAATGATTCTTTGAAGCTCTATGATAACGTAGCAACAGAGCCTAATCGATATGATTTCTATTGGCTATAAGGGTAAAATAATGTAAATTAATCATCGAATCTATCCAATCCAAAATGAGGTCAGTCTTTTTACGCAAATATATCTGCATAACTTATCAGCTTGTTCGGAGACTTGATTGGATACAAGTAAGTTAATATGGTTAATTTTTAACACAAGTCATAAGGGAGTGACGAGATGCAGAACATGTAATGTATGATTAAAGGTAGGAGGCAATTACGTTTTACCAGGCAAGCACTGATAAAAAGCAACGCTAACGGTTGTTATTTGCTAAAAATGGTCCCACAATTTCACACGCAATTATGATTGAGAAGGGTTTGCTAAAGTGTTGCAAAATGTTGAGTCACATTGAAGGTCATTGTGGGAATGAAAACCTATGGGCCATGGTTATTTGATTGACCGATTCCTTCTCAGGCCTCCTGGTCTGCTCTATAGTGTCCATGAAATTATTGCCATAAAGTTATTTGGTTACTTATTTTGTAGGCTCTTTTTTGACATTACTCCTATTTACAATGTCAAATGAATAATGCATGTAATAATAATAATGGTTTACTGACAGCTGTTAATGGGTCATTTCTCATGAATTGACATAAGCGTTGAGGATTGTAAGCCAACCAATTTAATATCCCAGGGATTCGGCTTAGAGGAATTCAGAAGAAAAAAAAATTTAAAAAAATAATATTTTGTTTCTTGACCAGAAAGGATTGGGGAAAGTAGTTTACTCTTAGGTGTACCGGTGTAGGTCTACCCGATGTTTAAGGGTTTTGGCCTAACCACCATGTTGCAGGTTTACTGTTGAACGCGAGCCATCTCTTTCTTCAAAATGTGGTTGAAAATTTAAAAAATAGTGAAAGGATCAGACACATGTACATTTTGTGGAGCATATCTTGGGTCTCATGGATTGATTATAATTATTAATACTGACACGGGTTCAAATATGACACCGTGGTGGTTAGGGGCATACGTGAGTCCTACTACAGCACACTCATTCAGAACATGCTCAACAGTTTGAAATTGACCACATGGACACAAGTCAGAGTCTTTGAGACCACAGAGCTTCATGCTTGCTGCAAAGCGGCCTGCACCTGATCGTATAGGCGGTTCAGCTGAACTCAGCCCTTTATTGGAAGGTTCAAAGTTTCCTATTTTACCAAAACAACGAAGGAAATAAAACATTCATGAGAAAAAAAGCACATGGCTATTTCATCAATGTTTGTTGCAGTTCAGGATTATATCTCATTACTTAATTATATTTATCTCATTATTCGTGATAAACTGATACAATGAATAATAGGTAACGCGGATAACTTAAATAAACATCATGTTCATGATTAATTATGATTGATGGCCAGTTGATTAGCAAAACTTTGTTTTTATCCCAGGCATGCAGCTGGCTGGAGTCGCATAGCTTTTCATTTGTAAAGTTGAAATAAACATTGATAAAAAGACATGATAATAATAAATGAGTAATCTTTTGATAACAATTAGAGTACACGGAGTGAACTGGGACGCGGGAGAAGTGGAACACTTCACACTTAAAATTGTAATCATCCTTTTATAATTAAAAAATACATATATTACAAGTGTTCTACTATCAATGTTTAGAAATCCATTTTTAGCTATAATACGATTTCGCACGTTGGTTTAATATCAACTTTTTCAAAACGATTTTTATTTATCGGCGTTCAGTGCCTAATGATTTACCAGTATACCAATTCATTGTTGTAACAATCTGGATATATAGGGCTGATATGACCAAATACTTTGAAATAATTTGTTTTGATTATGTAATTTTCTTTTCTAAAAAGTATCAACCTTTCTACACTACTATGTAATGCGCTTTATAAGAATATACTTATAATCATGATTAAGCATATGTCATGTTATCATTGATGAATTAGGGGAACGTGAAAGGAAAAAAGAGAACAATCCATCAATCCATCTTGGCAGGGCGCCCGAGGATAGGGGGCTCACGGGGTGCGACCCCTCCCCCTTCATTCACTCACACTTGATTATTGATGCCCACATAACAATACATACATGTATAGTGCTGGATACTTTGCTCACACTTTGATCATCTCACTTCTGTTAATATAGGCATTATAAAATTCGATCCACAGCTACCGTTGAGAAAATTGAAAAAATTTGTTTTAAGATTTGTATCCCTCTTTCTCCCTTCCTAGCTGTCACGCTATTTCGCTCAAACTTTCCATTCCATTCTCTTTCTCTCCAAATCCATCCCTCTGTCCTTATAATTTATCAAACTCCATCGTTATCTCTCTCTCTCTCTCTCTCAGTCTCTGTTTAGATTACACAAACCTATAATGTGAGGCGATTGTTGTTTGTCTGCATATGCTACGTTTACACTATGTTCCCGGCGGCCACGGCAACCACGTTTCAGGCCACCGTGGACGAAACGGGATTGACGTGACACAACGCGGGGGACCGTATGGGCAAACGTGACAGGCCGTGATTGCCGTGGATACCGTGTCAGATTTTTAAACTGTCAAAAAATTTGCCACGGCAGTCACGGTGCGGATGAGAAACCTAGTAGGCCGTAATAAACCGGGGTGAACGTAACAAAACGTGACAGTCCGTAATAGGCCATAACAAAACTTGGGGGTGGTCCGTAGTGGGACGGTTTAAGCGCCCGGCACGCTACGGCGCGTAGGGTTCTACTACGGTCTGACACAACGCGAGGTCTCAGAATACCATGATTGCTATTGCTTGGGAGACCTTCTTTTTATACCGCCGGCACAAACGTGCTCCAAAGCTGTCCTGCGTCCACGTTTCATCCGTACGAAGCCTTGACAGAACGTGAAAAACGTGTCAGACCTTGATCAGAACGTGGGAAACGTGAGGGACCCCATCAGAACGTGACGGGCCGGGAGGAAACGTTCTGGCATCCGTGGTAGAACGCGGGCCGAACCGTAGTGCGCCTTGGAGCCTCCGGGACTCAGCATACCCTTCCCCCCAACGGCCCGTCACGGATACCGTGATAGACCGTGAGAAAACTTGGCAAGACGTAACAAAACTTGAATCTAGAGATAATATTGGTCTAAATCCCACGGCGCACTACGTTTCGGGCAAACGGGGCTGCCGTGGCCGCCGGGAACATAGTGTAAACGTAGCAATATCGGCAGGATCGGAGGAGTGTGGGCGTTATACGTGCAGTGTGCGCGCGCTTGCATTGCTTGCTTGCATAAACGACTTTGGCGAAAGCGAGTATTTTTCACGTTGACATAATTCACCTCCCTCTCCCACTACATCTCTCCTTAATTACTAACCCCATCCCCTACATGCAGTCATACAGTATTCATGCACTAAAGGCTCTGTACTTTAGCCGATTTGACACAGCCTAACCAAATACTGACATTCTTAGGCTCTGCTCAAGGCTTTAGGAATTGTGGTACAGTTAGAGTCTGTGCATGGAGACTAAGTATTCTACGCTATTTATCCTCGTCTGTCTCTCTCTCTCTCTTTCTCTCTCTGTCTCTCTGTATCAAACTCTTTGTGGTAAATGTATGGCATGTTGTTTGCATTTTTGAGAGTCAAAGAATCGTTTTGAATTTCTTGAGATACTTTAATTAGCTGTCTCTCTGCCAGTGCCATAATTATGTCCTTAACCAGCCTTTATTCGAAATTCAAAATAAATTTTGAATTTTCTCATTGAATCATGTTTAAAAATAAAGAAGAAACAATCATCTTATGATTTGTCTATTATTAGTATGGTAATTCACACAATGAGAGGCGGAATTTCCGACACCGACATCCACGTAGGCAATTAGTAGGGAGTACAATTCCCCAAAATGTTATGGGGTTAGAAAACCCTACATGGAGGGGGAGCACTTTATACCCATCGATCAATGCATGTTTCTATCCCTAGATGTTTTTTGTTTTTTTATTGATGACGAGGAAAATGGGCGTCTCGCCCCACTTTGACCAACCTATGACAATAAAAAGCCTCCTGCTGAAGAAAAACAAAACTGACACCCCCCGCAAATGGCAGTACGCTGATGCTTGATATGATTTAAAGTCACGAGATCGTTGGCAAAAAGCCAAAGTTTCCACAGACGGCAGAGTTTCCATTTATATCCGTGTAAAGGTACATGCATTAGGAATGTCTTTAAAAATCCATTTATTCTACATATAACGTAGCTCTATACATGAATATAATTTCTTAGGTTCAAGGAAAATTACCCTCTAGAAATTCACACCCCGACTAATTACCAGTGGAACAAATAACCCACAAGGACAATTATCTCCAGAAAATTACCCCAAGGATATGTGTGCAGTGTATGTGCTTTCCAGCAGGACTGAGTTAATTTGTTATAATTTAATTGATTTAAAATCCTCATTATATACACAATTCCAATGCACTACATAATTATTGGTTCGTAAGAACCAATGACAATATGTGTAAATTACTATTCGCTTTGTTACTCTGAATTTTAGTTCGTAAGAACCAATGACAATATGTGTAAATTACTATTCGCTTGGTTACTAAGATCATTCTTCCCTTATTCGTCTGTTTTGTATTGTTTGTTACTTTGTTTACATTTCTGTGCACTCGTTTATAAAAACGCAATACAAAAATGCTGAACATTTTGCTTTTGATGTTTTTTATTAGTTGATCCGTGGACCTAGTAGGTCCATGGTTGGTTGTATCTGTCACACTGTTCTCTATTCGATAAAAGGACAGTACATGGCCTTCCACTCAAGCATTCACTTAATTGATTAGTAAAATGTGAATCATATATGAGAATATCCGCTGGTTACGAAAGGTGTTGTTTATCTATTAAACTTATTTCAATCTTTTTCATACAGGGTGGCCCTTTCAATGTCAAGTTTTGACTTTCCAATGGGACCTGTCTGTACAATTACAAAACCTAAAAACGTTACATATAGAACATGATGAGACACAAACATTTAAAACAAGGAATATAAAAATAATATAAAAATAAATGGATATACTATAAATATGTATAAAAATCAATTGATAAAAGTTGAATTTCAGAAGGTATTTATAATAGTTAGCTTTTTACATCGCCTCTTGAATGATAGGGAGGATAATGTACGTAGATCGGGTGATGGGGTGTTCCACAGGAGGTATCCTCAGTAGTGAAGGGATCTTTTGAAATAAAAAGTTCTTGGTTATTTTTCAATGAAGAGGGGGGAGTGTGTCGCAGGTCTCGAGTATTATAAAATGGAAACCTTTCTAGGAATAATAATTGGTTTCGGATACAAAGGAATGAGACCATTTCAGTTTGTATAAGTATCTGGACAAGTTTGGTATCACTTGACTCTTTGACGAACTGATTGCCAAGCAAATTCAGATGACTTCCAATTCTGATTTTTGTAACTACCATAATTATTTGAAATTCAGGTAGTATCTCGTTTGTGTGTGAAAATCAGTAGTAAACATTAAAATAAAATGCTGGTCACAGAAACCTTGAAAGTGGAAAAAATACCTTTGACTTTGTAATGTTATGTTTCGTTAATACAGATATATATGAAGATATATATGATTTAAATTATATTTCATAAGACGATGATAATTTAAAACGTTATTATTTAAGCGATTTAATTCATTCATTGAACACAGTACTATACATGTATACCATAAAAACAAATAAGTTTCATATTATTGTTGATATCGTTGGTACTGTACTCTGTATTATACTTGAATTACTGTTTAATTTGCCAGCTATACACGGTATATACAGTGCGTCCCAGAAAATAACGAAATCTAGAATTATCGATGTTTTATCATAACTTAGTCACATTTGAAAATAATAAATGACACATATTTGGAAAGCTTAGATCCTCCTCTTTCATTTTGAGTGATTAAAAGCAATTTGAAGAGAGGAAACCAAAAAGTAATTTGGCGGGCGGTATCTGAATTTCAAAAACAAAATCATTTTCAAGAAAGTTGAGTTTGTACCCTTTCATGTGACACCTCAATCAGGTAAATTAAGCAAGGAATGACAAAGTTCTGTTTGTTTAAAGTAGCGCTTGTATTTCCATAATTTTATCAAACAAATTCGTTTTCAATAGCTGCCGGAGCAACAGGCAGAATTGACAAAAGACTTATATCCGTTTCGTCGACAAACAGCCATGCATGAAGTCTTTCGTCTACGGCTAATTGTTAACGAAAGGGAGTTAAGTCTTTTGTTAATTTCGTTTGTTGCTCCGGCAGGAAGCCAGTGAAAACGAATTCGTTTTACCGGTAATACAATTATGGAAATGCAAGCACTACTTTAAACAAACAGAACTTTGTTATTCCAAGCTCATTTTCCTGATTGAGGTGTCATATGAAAGCGCAGAAATTCAATTTTCTAGAAATAGCATTTTCTCTTTGAAATTCAGATATCGCCCGCCAAATGACTTTTTTGTTTCCTCTCTTCAAATTGCTTTTAATCACTCATGCGTGAAATATGCATTATATTACTTGCATCAGATGAAGGAGGAGGATCTAAGCTTTCCAGATATGTGTCATTTATCAATTATAGTTGTGATTAAGTTATGATAAAACATCGATAATTCTCGGTTTCGTTTTTTCTGGAACGCACTGTATATATATATATATATATATATATATATATATATATATAATACTCCCTCACTCTCACTTTCTCACTCCCCCTCCCTTTCTTTCTCTCTCTCTCTCTCTCTCACTCTCCCTCTCTCTCTCCCTCTCTTTCTCTCTCAATGTCTCAATATCATCGCGTGCCTCTGGGCAACAGTGAAGTATTTCACTTTTGCGAACCCTGGCCATAGGTAACAATTTTTTTCCTCTTCTCCTCGCTTCTCATTTCCTTCTCTTTCTCTATTATTTGATATATTATCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTAATCTCAATAGATCAATATCATGGTCAATACAGCTGCAGAAGGAAAATAATAATCATTATTAGATGCTATAATATTGAACGCTTCACAAGTGTTTGTATGAGTTCAGTTATATACTTTTATTTAACGAGTCTTCAAAAATTCAGTATAAAATTAGATTTTAATAATTCGTCAGCATTTTCATGCACCTTTCACTGGCTTTCAGTCTCTTGTATCTGAATTGTATGACATTATTGTAGTCTCTGGTAAGATATGCTTATGTCCAAGGAAAGACATTGGTGACTATTATTTGTATATTTACACAGGCAATATTTTAAAAACAAATTCTACACCTCGGCGACATCATTGATCGCGATCTCGACCACGTCAGGTGGTGAAGATGTCGTTACGCTTTCCGTCGTTGGATGCTCCTTCTCCTCACTTTCATTCTCACCCCAACGTCCTCGGTTACCATGGGGCCGACGTTGGCCATCACGCCTGGGAGATCCTTCCTCGGTTCCGTTGCGTCTGCCGAAGGGACGGTCTCCGAAAGGCTTTCTACCGAAGGGGTTGAACCGAAAGGGCCTCATCTCGAACATTGGTCTGTCCTGATCTCCATCTCTGGTCTTGTTGTGATGACGATGGGGGTGGTCTCCCGTCTCATTATGACCCTGGTGACCTTCTGTATGGTTGCGTCTGCCGAAGGGACGGTCTCCGAAAGGCTTTCTACCATGATGACGACGTCCAGGACCATCAAACTGACCGGGCCCTCCTAGACCATCACCAATTTGTTGGGCAGCATCAGTTTCTTCTTCTCCATCACCACGCCTGCCTCCTTGGCCGAAGAAAGGTCTTCCTCCAGCACCATCTGGTCTGGAGCCACCAAATCCAGGTCCGTCGAATCTCCTACCACCCATCGGACCGCCATTTTGTCTGCGTCCATCCATTTGTGGGGCACCAGATTCAGGTCCATCGAACCTCATTCCACCCATTGGTCCACCATCTTGCCTAGGACTACCCGTCTGCATTCCACCAGGCCTTCCTCCAAAGCCACCTTGTCCTCTCTCTCTGCCATTCTCCTTTCCTCGTCGTTCATTGAAATCTCTTCGTGCGTGAGCTGTTGGGATAATGGAATAGTCGAATTTATTAATAAAATAATTAATCACATTCTGCATTATTATCTTTAGCTCTACATTGGTTTATGCTAGGAGATTTGTCTTCTTCGATGTCTGGTCCCGTTATGTTTAACTCTACCAGAATGCATGCATGCGGTTGGATATTCAAAACGCCTAATTTTAACCGAATTGAATTTGAAAACCAGATTGTAAAGTACCATGAAACATTTTTTTCATAAAAAAAACAGAATGCTGCAGGCCGCAGAATAAATGTTATTTTTTATTGTGTCACAAACAAGTTCAATATACCACTTATGAAAGAATACCGAGTAATAATTTGATTTCTTACCCGAGATAGCAAAAGCAGCCACAATGGCAACGATCAGTGTCACTTTCACCATGATAGTAATAGGTCTCTCCGATGCTATAGCTTTCTATAGGTTCGTAGCCTTCTATAAGAGAACGAGCTCCAAACTGATCAACTATCTGCTGGATCACTGAATTTATAGGCTTTCTACCATAGCTTGAAATCACTTTCGTACCAGACCGTGTCAGGAATTACCAGGTATACCCTGTTCTAGTTGATCATTGTTTTTATTTTGTTGTTGTTTTTTTTTCAATGTCAAACCTCGTATATATCTGTAACCTTTTTTGTACCAGGACGTATAACCAATTACATCATCTTTGGATATCCGGTAGAGCCAAGTATTGGGTGACACTGATATCGGGTAGTCACCCTTGGAATTTTCAAGGAAGCAACGGTCAACGATATTCTCTCTCTATCTCTCTCTCTCTCTTTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCATCCGCTTTCTAAACAAAATTACATTGATTCGAAAGCCATTAAGAAAGTTGTCAAAAAAATCTGGGCAATAACAACATCATTATAACAATCATAATCATCCCCATTATCATCATCATAATTAATGTATTGGTTATATAGATTAAGCTATGTCAATGATTTTGAACATGGATGTGATATTTCGAATCTTTTAGAAATAGAAACAGGGGATAAAGTGAATGTTTTTGTTATGAAAGTGAAAGATTACTATCGTTTGTTATTAAGTCAAGAGAACCAAATTTCATTTATTTTTGGAATTCCTTTTTAGATTGAGATAATGACTTTGTTTGGAAGGATGTATTTTCATTTAATTAAAGCAAATTCGTAATAATAAAGTTAAACAATTTAATTTCAAGATGATTCACAGGTTTGTTGCCTCAAAAGAAAACTTATATTAATGGCAAGTTGTGAATAATTTATGCAACTCTTGTGGACAAGTTGACTCAACTTTTCACTTTATGTTATATTGTAAGGATGTGACTTTGTTTTGGAAAATTATATTTAATTTGATAATTAACCAATATAAAAAAGATAAACCAAAAGCTATAAGTCGTAAATAAGGACATTGGAAACAAGAAATATTCTCTCCTGAACATTATTTTAAATTATGCGCAATATGCAAATTTATAAGTGTTAAGTTAAAAAGATTGTTAATGGTTCTGTTTATTACCCAAAGACTTTTTTAAAGTTTAAGTCGTTGCTAAGAGTGCAGCGTTTAGACAAATAAAAATGCAATAATCTTCTCGCTCGGGAGCTATGTCCCTCGCATAATATTCTTCAAAGTGTACAGTAAATATTCTAGAAAAGTGAAGTGTGAAAAAGATATATTGCTTGTTTTTATATTTTGTTAATACAACAAAACTTCAAAAACCTGCGGTGGGGGGGGGGGGATAGTCACTTCCGTCACCTTCACCCCTCTCGTTCACTATACTCCCTCGCCCTGGCGTAATGATGGGGGGATTGGGGGTAGTTGCCCCTTAATAAAGTTCAAACTTGATTTATTTCTAACTCGATACCAGTGATTTACAAATGTTTCTGAAATGGCATGGTTTTCCCTAATAAATGCCTAAAAACCCTGAGCTGAGCCCACGCCAATTTAAGCTCCACTTGCTTCGCTCGCTTTGCTCTACCTGATCGAAGCCCCCCCCCCCCCATACTTTGAAAAACGTTCCGCAGGGCGTGCATGTACCATGTAAGAGTGGGGATGCGGATTTTCTCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAACGTTTTTGATCATTATACGCCACCGTCCCTCGCATTTCGTCCCTCCCCCCCCCCCCCCCGGTCATGAAAAGAAATCGCCGCCCCTGGCAATCACCCTAGGAATTTGCGAGGTAGCTAACGGTCATCGGTTTCTCTCTCTCTTTCTCTGCCTGCCTGTCTCTCTCTCTCTCTCCTCTCCTTTTATCTGAAAATTGCATGGATTCGAAAGCCATTTAACAAAGTTATCAAAGAATCTAGGCATTAAAAACATTATTAAAACAATCATTACCATTACCATTATCATCATCATTATCATCATCACTACTAATCTCATTATCTTTTTCTCCTTCATCATGCAGGCATCATCATCATCAGAGGCGGCCCTCCTGGGGGCGGGGGTGGGCATAGTGCCCCCCCCCCACTTTTTCAAGTATCCCAAAAGTGCCCCTTTTTACATACGAAAAGTCCCCCTTAAACCTCGTAAGAAGTGCCCCTTCTACTTTGAAAAGTGCCCCTATTCCATTATAAACATGCCCCTCTCCCTTTCCACCTCAGATAAGTGCCCCTTTTTGGACACGTGCCCCTTTTCCTTCGACAGTTTTTCCTGTTTCCCTTTGCAATTATTGTCAATTTTTCCATTCTCATAAGTGCCTTTATTAACCTCAAGTAAGTGCCCCGTTTCAGGTACGAGAAATGCCCCTTTTCACTTTTGATTAAACGCTCGCTCGCATTTATAAATGATATAAATATGTTACTTATCAAAGGTCATTTCAGTCTTCTGTTATGAAAATGTGCTTATAACTCCATCAGAATTTGTTCTTCTTCTTTATAAGATGTAATTTTTCGCACTGCGACCTCACACCCAACTGCAAAGTACAATACTTAAGCAATATACTGTATAATGTATACAGAGAAACAAATATCGAAATCATCCCTGATAATGTAATCACAATGTAGGGGTATATACCAATGTCTAAAAAATGTGTATATTTTTGCTTTGAGATTTCCTATTGTAGGGGTGCCCTTTTTTTTCTCTAGAAAGTGCCCCTTTTTGTTGGTGCCCCCCCCCCCAACTTTGTAAATCACTCGGCCGCCCCTGATCATCATCATCATCATCATCAATGCAATCGATTCATATAGATGGTCATTATTCGAGGGATTTCTTTTCAAAATTTATTTCGTTTTTATTTGATTTCATAACTCGGAGCAAGAGGAATTTCGTAACACAACATTGCAGCTCATATTTTGTACGTAATAGGCCAACTGTTATGGAAATTTTGCGAACAAATGGAAAGTTACAAATGTTTTGATATGACGTTTTAAAACTAGATTTAAACAAAGGTAACTTCATTCTTTCTATTTCGGAAATTAAAAAAGTAATTCAGGACATCCGTTGAATCAGCCCCTATAGCTCGTATCATGTGATACCAGTTACAATTACAATAATTGCTTGTTTATATCCACATTCAAAAGTTAAGCTAATGGTATATAATGCAATATGAATGAATTGTAATTGCGTTCCTTGAAAATTATTTCGTCTCCCACAATTACAATTTGATATGTTACATACCTTTGCTCTTTGGCAGGTAGTCTGTTTGAGTGTTTGACTATATTGGACAAAACATTTCTGTGTCAGAACAGCATTTAGTTTTTGTTTTCCAGTATCAAATTTATTATTTTGAATATCTATTTATCAAAGGGTCGTGATGAATAAACCTCTGGGGTAGCGAGCCTTCGGAAAACTGAAAAATCAATAATACTTTTTACATGGACGGTATATGTAAATCTAGCAAGAGCCTAGCATTTGCCTTTCCCTTGGCTAATTGCCTAAACCAGCGTTTTTCAATCGATGTACTGCAAGACATCTCTAAGGTGCCACGAAAGAATTTAAATAATATATACAATTACAAAATAAATATATTGGAATTACGTCGATTTCTTGAAAATCTTAACAACCACAAATGTTTGTCTTCTTTCCATATGATGTAATATAGCTTTATTTCATTTTATCCGAAAATGATTAAATAATAGAGTTTTATCTGTTGTGCTTGTTAGTTTGATTCACTTTCAAATGTTCTTTTAATGAATTAAATAAAATAAAATTGTAATCATATCCTTAACTGTAATAATGCAGTATCTTGACCAGCCACATTTGATTTTTTTATAGATATATATCATTTCTGGCAGTTTATGCACTAGTAACAATTCTTATTACCAATTTTAATTTGAAGTATATTTACAGTACTCTTGTTTGGACACTCATGATTATTGTTATTCATATAAATATACTTTTTCGGCAGTACTATTCTATCGTTATTTTTTTAGATATTTATCTTGATTATTTCATATTTTGATGTACCATTATCATTTTTATGTTAATTCGATGTATTTACATTGTTTCATACGGCCCAATGGAAATCAGACTTGTAACTTTTGGGGCTTTTTTAACCGTGTATAAATAAAGATATCAAACTTAAGTATGTCATCTCATGAAAATAGCATTTTACTGGACAGAAAATTGTCTTTCTTCCTTGGTCACTTTGCTCCTTCGCAACACTTAGAAACATGATGATTTCTGGCAGCAAAGTGCCTCCATATTTGAGCATTGAGGTGCCCGAGTAGCCTATTTGATACCACGCTCGTTCTCTTCACTCGCTCGCAACATTTAAAAACTGGATAATTCTGGGCAATGAGATGCCCGAATTCCCTCATATTCGTGCTTCGAATTGGCCACACGTAATGATTTTTGTAATTATGTATTCGGCAGCGAGGTGCCTGCACGCTTAATTCTATATTTTATTTAGTGGGTTGACTATAAAGTCAATTATTATCATATTAAAAAGAGTTATTTTATTTTTTGCTCGCTCGCTCGCAAAATATAAAAACAAATGAGCTCATGTGTACAGTGATGGAGCAGGACGAGTGTCAGAGTGATGTTAAATGAATATTCTGTAAATTAATTTGGGATTAGGGTGATGCGAAACCTTGTTTTTTTTTTAAAAGGGTGCCTTGACTTAAAAAAAAAAGTTGAAAAACTTTGGCCTAAACCTTGATCTAATTCAATGAACTAAAATTTCTTACAAAAAGAAGCAAAACTTTGTGTTCTTGCTGGAAATGGCCTCCCCAGACAAGAAATTAGTCGAGTAAATGGGGAGGCATGGGAATCGACGGGACACCTCCTCCGCCAATCGAAAGCCAAAACATTATTTTTTAAATTAAAGAGTTGAAGAAGAAGAAGAGAAAAAAAAGGAAAAAAAAAAACAGGAGGTAACAGAATTTAGATAGAAAAGGAGTAGCAGAGCTGAAGAAGAACAAGACAGGAACTGCTTAACCCTTACCACAACCTCGGAGATACTTAAGTAGTAATATATAGCAGTAGTAGTAGTAGCAGTAGCAGCGGCAATATAATAGTAGTAGTAGTAGTATAAAAGCGTAGTAGCAGTAGTAGTAGTAGTATTAGCAGGCGCGCCGGAAGCAGTTTGGGATGGGGGGGGGGG";
//        mn->seq = "CTTTGTATGCTTTTTGGTTATTTTCTTTCGCGAAGGCTTGAGTCGCTCCTCCTGCCAGAAGTGCAGTCGTAAAAGTCAGAACTGTGGCTTGTTTTACAATTTTTTTGATGTTCATGTTCATGTCTCCTTCTGTATGTACTGTTTTTTGCGATCTGCCGTTTCGATCCTCCCGAATTGACTAGTGGGTAGGCCTGGCGGCCGCCTGGCCGTCGACATTTAGGTGACACTATAGAAGGATCCGCGGAATTCAAGAAACAGTTACAGTCTAACGTTATTATCGTGTAGCGATTTTGGGCGTATTTTCCTTTGCGCAATAATACACTAGACTAGATTGCGATCTGACTATTCGTCGTTGGAAATCTTCGCTTGTTCGAGCTGAGCAAGAAGGTGACGTAGAGACAGAAGAATTGTTAAGGTACGGCATGTGCACCCAAGCATGACTTGAACGCCGCATGTTGATTGGCCAAAGCACTAAACATTTTGAATGACGTAGGCTATTATATCAGTCTAGGGAGTAGAATGGCTAGGGTAAATATTAAGTTGTATTGTATATCGCATTGATAATTGTCTGCCAAACGGGCCTGCCCAAAAGAGAAGATTAAGTTACACGAAACACAAGAATATCTCATAAATTCATAATTTATATTTCGTCGGCAAATTATAAATGTTACGCATTGTGTTTCTTCGTTTTTCCCCAGGCTTATCTAATAGGTGGCGGGAGACATGTGCTTCAGCCTGTGATATTTTTTTTTATTTTCTCCTATTCTTTCTTTCCAAAGTTGGTATTTTCTTCCTTTTTTTCTTTCGATTTCTTTCGATTTGTTAAAGGCAGGCACAACATGCAATCGATATGACCAAAGAGATTTGTGGTTTCTCTTTGATATGACTATGTGGCTATGTTTTTTTTTTACGAAGGGGGTGACAGATGCTGTCGTACGTTGTACTTCTTTGCAAACCTTGCAAGGGCTTCTTGCATTGCTAGCATTGCATGGGCTACTTGCGCTGGTCTACGGAACCGCATGGTGCAGTATACGGAGTGGTGATGTTTCGCACGGCAGTTTTCTTTTCTTTCTCTCTCTCTCTTTCTCACACTCTTTTAATGTTGCCGTTCATACCTCAACTCCAATAAAGCCTCAATCACACCAGACTAGTTCTCACTACGGAAAATGTTGCCAAACTGGTTCAAACTGTGGGCGAACGGACTGAACGTATAATATAGCAACTGTTTAGGATCTTTTAAAGAACTTCATGGATGGGCTCGGACTGGCTTGGACGGGAAGGGCGGTCTAGAACTTCGATTGTCCCCGTCCAGCAATTAGACGGAGTAACAACGTACTCAGTGGCGGATCCAGGGGGGGGGGGTTGGGGGGGTTGCAACCCCCCCCCCTTGGGCGGCCAAAAAAAAAAAAAAAAAAAAAAAAAGGAAAAAAAAAGAAAAAAAGAAAGAAAAAAAAGAAAAAAAAAGAACGCACGCTGGAGAGAAACCCTTTCAAGTGCAAATATTGTGAGAAGATGTTCAACTCATCAGGTAACTTAACCAAACACATGATAACTCTGGATATAACTGTGTATTTAACTGTTACTTGATTTAAAAAGACGTTATCATGCACTGCTTTTTATGGTCTCCATTGTGTTAATAGATGCCAGGAAAACATCACTTACGGGGTTCAAATAAACACCCCAAAAAAAATCGGCTCGCCGCTTCGCGGCTCGTGGACGATGCTCCGCATCGTCATTATTTTGGTCAACCCCCCCCTTAGCAAAAAGCTAGATCCGCCCCTGGTACTTAAAACAAGATGGACGTACTACGAACGGAATTGCCGTTTTTCCATCCTTCTAGGAACGTGTTAAAACGTCGAGCATCAATTACCAGTGCGACCGGGTCAATTTGCTTGTTCCTTCCGTTAGCAAACATACCCAACTCGAATGTCGTGATCTTATTAGCATGAAGCGTTGAACAGGAATCAGGAATGACAAGAACGCGCTGAGTTTTTATAGCCTCGGACCATGCAGGGCGGGATTGCACACGCACTATAGGAGCGCGTTCGAACGCGCTAAGAACGAGATGAAGACGAGTTGAAATGTCTTTGAACAGAATGGTCGGGTTAAAAACTAGTTAGAACAAGCTGCGAACTAGTCGGAATGGGTTTGGACGTATTGGGAACGAGGTAGAACTAGTTAAGAACGGGTTTACTGACTGCGGACGGGCTTGGACTGCTTGGAAATGGCAAATTTACCCGATCTGAACTTGATATGAACTGACTGTGGACGGGATACAGACGAGATGGACGTATTGAGGACTGTTTATGAACGAAGTAGGATCGTATTATAGAAATTATTTACAAATTCTACCCCGACCATAATCCGTTCTACGGTTTTCTAAGACGTTTCCAGAACGTAGTGAGAACGAGTTTGATGTGATGAGGGCAACTTAATTGCAATTAGATATATTATCTCCAACGTAGCTCAACATGTTCAAAGTCGTTATTCTGTCACAGCGCCCCCTTCCAATATGTCGGCCAGATAATGAATAGGGACTTTTAGATTTGCGTGGACTGATACGCGTGACGTCCATTCAAGGCAGTGCATTAAAAGGACGTCGATCGTAGCGTGGACCTATTGTCTTCGTGCGAATCTGAAAGTCCCTAATATCCCGTGCTAGTCTGCATGTGAGCCTCTTCCAAAATAGTGTCCCTAACAAATCCGTATGTATTTGGCCTCTTAAAGTTACTTTTAGACCATCTAAACGTTGACAGTTATCAGTTCAGCATAACTGTCATGTTTAATCCATCATACCACGAAAAAGACACAAAATTCAGATCGATTACAACCTATAGCTTGCTAAATAGGGTTTATTATACAAAGAATAATAAACTGAAGGAAAACAAATACTTTCCGTCAGGATGAGCCCTGGATCAATTACCTCCTTCTGCCCATCGGCCATTTTCTTCCCCTTCAGAAATGATTATGCCTTTCCTTTACCAACAAAATCTTCCATCATGAAGGGGGGTGTTAGCAGGCGCGGATTCAGGATTTCTCGGGGGGGGGGGGGCAATTTTTTTGTTTGTTTTTTTTTGCTCCCGAAACTTGAAAATCGCGCCGGCCGAAAAGTCGCTCAGTGGGGGGGGGGGGGGGGGGGGGTATGGTACTGAGCACATTTTTTTTTTTTGGGGGGGGCCCTCAAAAATTAGAAAATTTGCACCCCCCCCCCCCACAAAAAAGCTAAAAATCACGTGTTTTTACACTTTATTTCCCAAGGGGGGGGGGGCCGAGCCCCCCGCCCCCTGGATCCGCCACTGGGTGTTAGACCCCCTCAAGATGATATCAAGGCTTCGTGTTCTCTGACGAAACGTTACCCCGGCCGAAATGCTCTGGATGTCAGATTTGTCTATACTTTTTGTAAATGTTGAATAATGTAGATACTATTTCAATTATACTGGCATAGACATATAGCATATACTTAGATTGTAGAAGTACTTAATGATTATATTGTGCTGCAAGTAGAATATGATTTCATTCTGTTTATTACTTGATTACCATTAGGAAGTTTAACAAAATATGTCTTCCATATACTAGGCTTTACATAGTTTAATATTTAGTATTGTACTCTTCCTGACCGGCAATTTAGTATAAAATGATTATCTTTTTCTTTCATAAGAACACAGATAATTTAAATCAGTGCTAGACATTATCTGAGAAGACCTCGTAATAATTAAACTAATATTTGTTATTGACTCCTCTTCAAGTTAAGATTCAGTGATCAGACAAGGCTTTATGTACAAGAACAGAGACTTTGCCTTTTATTGTTTGTAAATTGTGAATTATAAATTTAACCTTTTGCTAAGGGCATTGAAAGGTGGCGTAGGACGGGTCGATCCAACCATTGATGGGGCACTGGCTTCATGGGGGAGAGGGGTATGGGGGGGGGGGTCTTATTCGGCAAATGAAATACCAAAAGCCAAGGCAACCGATCAGTGGAGGTTGTGGGAGGGTCAGGTAGTTTCCTCGCTTGTCCTCTACCCGTACTCCCACTCCTTTTCAATCTAAATACGTTTACGTCTTTATTTGAGTTGTTTGCACTGGGTAGGGCGGGGGACGTGTCCCGTGCCTCCCCATTTCCACTTAACTGATTTCTAGTGTCGAGAAGCCATTTCCAGAAATACCACAAACTTTTGCTTCATCTCTTTTGGTTCCTTGGATTAGCTCAAGGTATAGGCAATTAGCCTAGGAAAAGGCAAGTGCTATAGGCTCTTGCTAGATTTAGAAATACCGTCCATTCACATAGTATTTTTGATTTTTTTTCGGTTTTCCAAAGGCTCATTTAACCCATCCGTTGCTTTATCACGAACCTTCAATTAATATTTTTTTTTCTTCAGAATAATAACTTCCGCACTGTAAACGAGACCAAAATTATCGTTTTTACAAAGAAACTTTGTAGTCAAACAGATTAACAAAAATTGTTGGAAAGTAATTACGTTAGATCTATTGACATACAGTTAGATCTAATTTTCATGATGTGAGTATAATTATGCAATCATTTTACATATTACTGGTTTGAATTGAAACGAACTATTGGCTGATTCAATAGACAACCTGAATTCATTTTCATTTCCGAAATAGAAAGAGTAATCTGACCTTTGCCTAGATCAAGATTTGCTACACCATATCAAAACATTTGTAACTTTCCAATTTGTTCTCTTTGTTCACAATTATGTTTACAAAAATAAAAATGCAACGATGTGTTATGAAAGTCCTCAAGCTTCGAGTTATGAAATCATTAAAAAAAAATACAATGGTTAAGATAATGATGATGAGGATGGTGATGATTATCATGTTTTTGATCATGATCATTATCTGATCAAGATAATGTTATTAGTGTATCATAGAGAGATAGACAGATAGACAGAGAAAGAGAGAGAGAGAGAGAGAAACCGATGAACGTTGCATCCTCCCAAATTTGAAGGATGACTACTACGCCGCGACATCGGTGTCACCCAATACCTGGCTCTACCGGATATCCGAAGAGGATGTAAATTGTTATACGTCTTGGTACAACACTCTCTCCAAAACAAAACATGATTTATAAGTCTAAGTAGGATAAAGGTCAAACACTTTGACAAAAACAGCAACATTGCAAAAACAGCAACATTGAAAAAACAGCAACAACACATTAAAAACAATGATCAACTATAGAACAGGGTATCCCTGGTAATTCCTGACACGGTCTGGTACGAAAGTGATTTCAATCTATGATAGAAAGCCTATAAATTCAGTGATCCAGCAGATAGTTATTCAGTTTGGAGCTCGTTCTCTTATAGAAGGCTACGAACCTATAGAAAGCTATAGCATCGGAGAGACCTATTACTAACATGGAGGTGAAAGTGACACTGATCGTTGCCATTGTGGCTGCTCTTGCTATCTCGGGTAAGAAATCAAATTATTACTTGGTATTACTTGATAAGTGGCAAATATTAAGCCAACAAAAGGCTCACAGGAGTATATTATTATTTCATTTATCATAATATGTGTTTCTTACCTGTTTGTTACACAATACACAAAATATTTCTTCTGTGGGCTTCAGCATTCCGTTTATTCGAATGAAAAACACGTTTTATGTTTCTTTACCATCTGGTTTTCAAATTCATTTCGGTTGGAATTAGGCGTTTTGAATATCCAACCGCATGTATTCTGATAGGGTTAAACTGAACAGGAACAGACATCAAAGAAGGCGAGTTACCATTCTTATTTGTCACCTGCAAAAAACCAATGTAGAGCCAAATATAAAAATGCAGAATGTGATTACTTAATTAACTCGATTTTTCATTACCCCTTACAGCTCACGCACAAAGCGATTTCAATGAACGACGAGGAAAGGAGAATGGCAGAGAGAGAGGACAAGATCGCTTTGGAGGAAGGCCTGATGGAATGCAGATGGGTGGACCTAGGCAAGATGGCGGTCCGATGGGTGGTAGGAGATTCGACGGACCTAGATTTGGTGCCCCGCAGATGGGTGGACCTAGGCAAAATGGTGGACCAATGGGTGGCAGAAGGTTCGATGGACCTGGATTTGGTGCCCCGCCGATGGGTGGACCAAGGCAAGATGGTGGACCAATGGGTGGAAGAAGGTTCGATGGACCTGGATTTGGTGCCCCGCAAATGGGTGGACCTAGGCAAAATGGCGGTCCGATGGGTGGTAGGAGATTCGACGGACCTCGATTTGGTGGCTCCAGACCAGATGGTGCTGGAGGGAGACCTTTCTTCGGCGAAGGAGGTAGGCGTGGTGATGGAGAAGAAGAAACTGATGCTGCCCGACAAATTGGGCCTGGTCGGTTTGATGGTCCTGGACATGGTCATTATGGTCATCATCAAGGTGCAGGAAGACCTTTCTTCGGCAATCCTCCTCCTTTTAACCCAGAACAGGAACCGCGCAACGACAGCAGCGAGGAGGATGGCCGTCATCACCGTCACCACGATCGCCACCACGCCCACCATGGCCACCATGGCCACCACGAACACCATCATCAACATCATAACCACACAGAAGGCCACCAAGATCATGACAGACCGATGTTTGAGATGAGGCCCTTCCGGTTCAACCCCTTGGGTAGAAAGCCTTTCGGAGACCATCCCTTCGGCAGACGCAATCACACAGAAGGTCACCAGGGTCATAATGAGACGGGAGATCACCCCCACCGTCATCACAGCAAAAACGTAGATGGAGATCAGGACACCGGCCACCACGGCCACCATGGCCACCACGAACACCATCATCATCAGCATGACCACAGAGAAGGCCACCAAGATCATGACAGACCGATGTTTGAGATGAGGCCCTTCCGGTTCAACCCCTTGGGTAGAAAGCCTTTCGGAGACCATCCCTTCGGCAGACGCAACCACACAGAAGGTCACCAGGGTCATAATGAGACGGGAGATCACCCCCACCGTCATCACAGCAAGACCGGAGATGGAGATCAGGACAGACCAATGTTTGAGACGAGGCCCTTCTGGGTCAACCCCTTCGGTAGAAAGCCTTTCGGAGACCGTCCCTTCGACAGACGCAACGGAACCGAAGAAGGATCTCCCAGGCGTGATGGCCACCCTCATCCCCATGGTAACCGCGGACGTTGGGGTGAGAATGAAAGTGAGGAGAAGGAGCATCCAACGACGGAAAGCGTAACGACATCTTCACCACTTAAAGTGATCGAGATCGCAATCAATGAAGTAGACACCAATGTGGTCGCCGAGGTGTAGAATTTGTATTAAAAAAAGAAGAAAACAAATATTGCTTTTGTTAAGATATAATATAAATATATATATATATATATACTTCACCAATAGCTTTCACTGGACAGAATGTTTTTATTAACGCGTCATAACCGTACCAAAGACAACACTGTGTTTCTACAAATTCATTTACAACAGACTAAAAGCCAGTGAAAGGTGCTTGAAAATACTGACGAATTTTAAATATACAATTTCATACTGAATACTAAAAAATTGTAAAATAAACATGTAGCTGAACTCATACGAAGGTAATGTGCAATAATAAAGCCCCGAATATTTCCCTCACGCAGATGAATTGATCATTGTTATTGAATCAATGAGGAAGACGAAATGAAGGGAAATATTGCTTTTGTTAAGATATAATATAAATATATATATACTTCACCAATAGCTTTCACTGGAGAGAATGTTTTATTAACGCGTCATAACCGTATCAAAGACAACACTGTGTTTCTACAAATTCATTTACAACAGACTGAAAGCCAGTGTAAGGTGCTTGAAAATACTGACGAATTGTTAAAATACAATTTCATACTGAATACTAAAAAATTGTAAAATAAAGATGTAGCTGAACTCATACGAAGGTAATGTGCAATAATAAAGCCCCGAATATTTTCCTGACGCAGATGAATTGATCATTGTTATTGAATCAATGAGGAAGACGAAATGAAGGGAAAAGGAATACAATGTATATAGAGATAGAGAGAGAGAGAGAGAGAGGGATAGAGAGGGATATATATATATATATATATATATATATAGAGAGAGAGAGAGAGAGAGAGAGAGAATGAGATAGAGAGAAATAAAAGAGAGAGAGAGAGAAAGAAAGAGAGAGGGAGAGGTGAGATAGATAGAAAGGGTAGAGAGAGATGAAGAGACAAAAGAAATAGAGATATGTAGAGATAGAGAGAGCGAGAGGGAGAGACAAAGAGACAGAGAGACAAAGAGAATATGAGAGTGAGATGGAGAGAGATATGAGTGTGAGAGAAAGAGGGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAGATAGATAAGATTATTTATATATACAATGCGTATAAAAAAAAGTGTGCCCAACTTTAGTAACCTCTACTTAAAAATTTATAACATATAAACTGATATCCTGTCTAATGATTTTATATTCATAATTGTACTGCCATGGATGATTGAGCAGCAGGCTTTTTATAAAGTTTTTTCAAAATCCTTTTTGAACCAAGTTTGGGCCGAAGCAATGGATGCGGGTCAAAAGTCATTATGTGTGGGTTATCTCATTCCATAAACAATTGTTCTCTTATGAATACTAACTATTAGGCTGTTGAAGATAAAAGGTGTGAATATCAACTCAAGCTATAAAGGAAATTTTATCACAAACAAAATGTATGATATTCTCTTTGTTATGATGAATAAAATCTTAAATTCATTCGTATCCCTTGTTTAAGGAAACCTTTCTCAGTGATAAAAAACCGATCGGGACGTTATCTTCTTTATTAAAGTCTCATCGTAATAATAATATGAAAAATTTATATAGCGCTTGTGACAAAAGTTTCAAAGCACTCGTGTGTTCGTTCCTGCATTTGGATGTAGTAACCTTTGAGTTTTCATTCGCTATTCAATATAAACACCATAATGTGCTCATTACGCGTGTGAAATTGTGTAGAGGACAAGTGAGCAGAGAGAAATTAGATATAATGAAACAAGAAGTTAGTTGGCTGTTCAAAATAACTAAAAGTTTTCCCCGGGGTCGTAACAAGGTTCTGCGAAACAGGGGCGGATCTAGCCGGCGGCGAGGGTGGGGGAGGGGGGGGGGGCAATTTAAGAAAATAGTTAGCGCCGAATTAGCCGGCGAAAAAGCAATAGGGGGGGGGTTAATGAGCAATAAATTGTTATTTCTTTTTTGCTCAACGGGGGGGGGGGGGGTGGGGTAGTGCATGCGCGGATCCAGGGGAGGCCCCCCGCCCCCCAAAAAAGTTTTTAATTTTTTGTTTTTTTAAATGGAGACAGATAAAAATTTAGTGCTCACTACCACCCCCCCCCCTCCCCCCCTACTGAGCAAAATTTATCGGCCGGCACGATTTTCGAATTTCACCACGCTAAATTAAAAATTGATTTCAAATTTGGGTCTCCCCTAAGGAATCCTGGACCCGCGCCAGTAGTAAGCACTATTTTTATTTACTGAACCCCCCCCCCCTTGAGCAAAAACAAAAAGAAGGTAAGATTGTGGATACTATCTTTGTTTAAATATAATTTTTTGCCAGGGGTGGGGGCTCGTCTTTTTTTTCTTCATTTTTTTATAGGTAAAATTGTGGTTACCGTCTTCGTTTTAAATTTTTGTTGTTGTTGCCAGGAGTGGGGTCCCGTGCGGTGCCCTCCATGCCCCCATGCTGGATCCGCCACTGCTGTGAACGCCTATAGACATTATTTTGAAAAATGAAGCACAGAAGCACGCGGCATAATCTAACTTTTAAAATGATTTGAGACCAAAGAAGCACATGCATATGGCCTGATAGTGTGGTACGGTATGTCACTGAAACAGCTCCCGTAGGACTTTTGTGCCAAAATTGATTTTTTGGCATTTTTGGATTGATTGGGGTCCATGCTATCAGAAATTGATGGGGTTCAAATTTTCCCGACCCTCAACCCCTCTCTGGGGGGGGGGGGGGGGTAGTTTCCCATCATGAGACCCCCTAAAATTGGGTGTGTCATTCATCTAATGACTGATTTCACAAAACTTGAGTCGTTTTCTAGAACTGTGATTATCCCTCTACAGCGAGAAAGGTAATTTGAAGACTTTTAGTTGGATAAACCCTCCAACCAACCCATCCCCATCCATACCCCACCCCGTCAGCCCAAATTGCGTACCTTTGTCATTGTTCGGCCAGCTGTTTTAAGGTGTAGTACTTTTATTACCATTATCATATATATATTTACTAGCCCATTTCCAAACACAAAGAGACAGAAACAAACATACAGACAACCTTCATGTATTCACATCAACCTCCCCACACACACACACACACACACACGCACGTTCACACACACACACACACCCTTCCTCCCACACGCACATACTCACACACACCCAAACATTCACAGGCAAGTGCGCTCGCGCGCACACAGCGCGCACACATTCCTACATAATTATGCTGCTCTGTGATCATTGTTTCCAGCTCTTCTGTGGTAAAAATTGTTTTGTTCCTTGGTCTTCGCAAGTTTGACAGAATAAAGACCCCTATTTTGCAGTATATTGCCTGATTTTCCAAGGTAAGCATAAACATATAGCCAATGGCAAATCTAGATTTTAAAATTTTAAAAATGTGCTCAAGCATTTTCGGTCATTAACAAAGGAGCAGTATCCCTCTCTTGGGTTTTGTGTGGATGGTGTGCTCTGTAAGAGTGTGTATTAAATAAAAAGAATTACCAATGGGGAGCAGCAAAGTTGAAAAGCAGACGATTTTGGTGATATAAAATTTATTATTACTGATTCGTACACTTTCGCACGCACCAATGACCTTATAAAACAATTAGAAGGCCATGATCTTTATAAGCCATGTGCGTTGTAATTCTGGAATATTATGAGGAAATTGGGGAAAACTCTCGAGTTCGTGACAAGAACCTTGTTCAAAACCTTTCAGAAGTCTACATTGGAAGTGACTATCACAAGATACTTGATCTAGAAAAGAGAGTTGAGCAGTCAGTTCTTCAGCGAATGAAGGACGCTGGGGGTTTTGTCTTCCTGATTTTGTCAAGAAGGGTGTGAACATCTGGTTTGCCATTGACAACATTGATCTCTTGGAAGACACATCCACAGGGCAAGGAACCTTTCACGGGACAGTGGTTGTAATCAACCAACAGGCTGTAGATGGAGAGCCAGTGAATCAACCACTTGTCATAACCGAGAAACTTTCTTCGCAAAACCTTCTAGCATTTGAGATGAATGTGCTTCCAGAGCCAGTCATCAGAACTAGCCCACTGAGATTTCAGGCCTACAAGAAGAGGAAACAAAATCTCATTTCTCAAGAATTCACTCATACCTGGGCTCTTGTGAACTACCTTACAGCAGATGACAACGGAGAAATCATCCAAACTGAACCTCAGCTGCACGATGAAGAAGCCCAGCGCAGTGAAGAGATGCCAAATGATGCAGAGTCTACAAATGATGAGGAAACTGAGCCTTCTAGCAACAAGGGGGAAGCTGGGTCTTCAAAGGACAGAGAGGAAACAGACTCTGTCGTGATCATCAAAAAGCAGGTGAAGAAATCAGATAAACTAGCAAAGAATGTTCTTCCAACATGGGCTGCAACTAGGTCCTTATTGTTGTCTGAATCCTCTCCCGCTAGCACGCCTACAAACACATCGGTGGTTGCCCCATTATTCAAGACATCGCCAACTGATTATGGGACCCTCTATACTGTTCTTCGGTTATCCCAGGGAATATCTGCGACTGTTGTTGGCCCTCACAGAAAGACATTGATAACGCTAGATCTCGATCTTTACTCCCGAGCGTTGAAGATTCAGCAGTCGGTGGGAAACGCCAACTGGATCTTAAGAGCGGGAGCTCTCCACATTGCTTTCGCTGCCTTACATGCCCTTGGCAAAACCATCGACGGTAGCGGTCTTGACACGTGTGCTATCGAATGTGGTGCGTACACCTCAGCATCTCTTCGCAGGATCTTTGGTGGTAAAGCTTATAAGCGTGGCCGTGAATTCCACATTACTGCAAGCCTTGCAATCATGATGTTGCGTTTTGATGCCATACTATCGGATCTTCCCAAGGGTCCAATTCGCATCCAGTGCAATTCTCTCAAAGAAAAACTTCATGGACGTGACCCAGAGATGGTGGAGATCTATGAAGAAATCCAATCCTGGTACTCAAGCAATGTCAAACCGCTTGAAGAAGCTGAAGACCTTGGCAAGTTTGCCCAGTTCCTGACTCAGTATCTTAACCAGGTGGAGAGTCTTCTACATCTCATTAGTTCATGTCGATCAGGGTACTGGGAGGGCTACTTGTCATCACTAGAAGACCTCACCAAGTACTTCTTTGCTCGTGATCTCCTGAACTATGCTCGCCTGATGCCAGTCCATCTCGCTCAAATGAATGCCCTAGAAGAGGATGATCCAGAAACATGGAATGCCCTCAAGTCTGGAGACTTTGTGGTGGCAAAGTCAGAGATACCCTTTTCTCTTCTCTTCACTGATCAGGCTCTAGAGCAAGAAATCAAGAAACTAAAGGGGAATGGTGGCATGGTTGGGCTCACAAGAAATGAAGCCGCTCTGGACCGACTAGTCACTACCACACCTCACCTCGCTGCTCTGGTAAATTACTACCTCAATGACTTCCCAAAAGCTACTGGAGCTTCGGTGAGGAAGGAGCATCATCAGCTCTCAGGAGACATTGCAGTGAGGTCAAAGAAGAATGCTCTGAAACTACGCCACCTAATCGAGCTGCAGAGCGGAGTCAATCCTTTCAAGGAGAAGACACAATTGAAGAGTCTGGTTTTTTTATTTTATTATTATTTTCTTTGGCAGAAAACACATGATATCAAAATACAGAGTTGAAATAAACAATTAATAGCAAATCACACACAAAAAAGGAACACTGGGTGCACCTTGCTGTGCTGCCAGACAATGGGCAAAGTTAACCAAATAACTATAGTCCGCAGACCGCCCTCACCTATATCCCAACCAGATGCAACCAGCGTTGCAGATTTATGTGACAAGGACAATATAAGAGATACATGTACCACATTCAAATAATTAAAATGCATCAAGACAAAAACTTCAAATTTCCCCAGGGTTGCATCACACGCCAGAGACTCATGACCATACCAGCAGAAAAAAGTGGCACAATATGAATTATAAGGGGAAATTTTGTATCCGGAGTGTGATGCAACCCTAGAAAATCATATGGACCTATACATACATTAAACATAATACATTATATATTAAAAAACATGCAAAATGTGACATTACAACAGTGATAATATTGTACAAATAAAGAAAATGCATATACGTGTAATTAAAGGAAATTGCATGAAAGTTGGAAGAATCAACTCAGTCGAGGTTAAAATCACTGAATAGTGAATATACATGTAGAACCATAAGCTGATAAAGGTCATCCAATAATCTATTGAACAAGTGAATAGTAAATTATTAAGTTCAGTTGAAATCTATGGAAAATAAGTATTAACAATTTATGAGTAAAAACTGCTGTGGCATGTTGGCAGGCCAGCTTATTAAAATATACATGTATTATTACAAAACATGCTGAGCTAAATATCTGAGATTTCATAATAGTATTGGCTCAACTAAGAGTGAACATTCAAATTTAATGGTGCCTTTATAATCAAACATTTTACTTTCCCCCTACAGAAATATTAATGTGTATATATAGATAAATACCATTAACTGGAATCAATACATAATTAAAATTCAGATAAATTGCATAAATATATAAACCAACATTGTCATATCAAACTCAAGATGCATCGGCAAGATTAAGATGGGATTTTAAGTGACTGCACCATGTGGAAAAAGACACAGTAAGGACCGAGTTCCTGAGTTCAGAAGGTAGGGTGTCCCACAGCCTGGGGAACTCATGTATTGTATTGCACTGAAATGCAGTCCATCTTTGCGTGCAAAGAGATGCCTGAACCTGGGTTCACTACAATGAGAGAGTTCACAACAGGGCATTGAATGTGTTTGAACAGACACTTAGTAATGAAGGCCATTGTACAGTAAGTCCTTCTCAACTTGAGAGAATTCCAATTAAGCACTTGTAGGCGTTCCTCATAGGGCCTTTGGCCTCGCCGTTGCTTCAAGGCAAGACGTGTTGCCCTGCGCTGTATTGACTCAAGTCTGTTGATTGATCCTACTTTGTGCGGTGTCCAGACAGGAAGACCATATTCCAGGATTGGGAGGACCAGACTCTTGTACAAGCTGAACACAGCCACTGATGAAAGACCCATGGCTAGGCTGGTTATCAGTCCAAGTAAGCGATTTGCTCGAGAGACGACATACCCAACATGGCTATCCCAAGACATGCGATTGTTTATGATGACACCCAAGTACTTGGTTTCAGTAACCACCTCCAGTTGTATATCATTCAGATAGTACTGCGGTTGGCCCGGTTGTCTTTTCCAGGATATTCTCATGACTTTGGTTTTGGAGGTATTCAGACTCATTTTGTTTCTGTCACACCAGTCCGCAAGAGCATTGACATCATCTTGAAGCGCTTGCTGATCCATGGGGTCATGGATAGGTCGTGAAATAACAGTGTCATCTGCAAACAATGCACATCTCGATAAAACAACCTCAGGGAGGTCATTAATGAATAGGTTAAACAGGCATGGTCCTAGGACACTCCCTTGTGGAATACCTGAGGTAGGTGTTGACCACTCAGAGTAATGTCCTTTGTGTTGGACGCGTTGTAGTCTGCCAGATAGAAATGCACGGATCCAGAACCAGATATCGCTTGAGATGTTGTAGTGAGAAAGCTTGTGCAGAAGGACATCATGAGGCATACGATCAAAAGCTTTCGCGAAGTCGAAGAATATCGCGTCAATAGAAGGGGGTCTACGCTTGTCCAGGTCCTGGAGCCAGGAGTGTGTCAAAGTCATTAGCGTGGTAACACATGACTTATGATGAACAAAACCATGTTGGAATGGAGTGAGGAGGCAATTTGATACCATGTGGTCAGAGATCGCTCTTGAAATCAGAGACTCCATCAGCTTGATGACAATAGATGTTAACGCGACCGGGCGATAATTTGAAACAATCAATCTATCCCCAGCTTTGAACACTGGGACAATATTTGCAGATTTCCATCCTTCAGGTAGTGATCGTTGTCTGAGTGACATATTAAAAAGACGTTGCAGGATAGGAGATATCACATCCACTGATTTCTTGAGAAGAACGGGAGTGATTAAATCAGGGCCAGGACTTTTCTTGTCCTGTAACTGCATGATGCCCCTGTGAATTTCCTGAGTTGAAATCGAAATGGAATTGAGTTGAGGTAGTGGATGGTCACATGCACTAGGACCATCGGGCAACACACATGGAGGAGGGTAGATTGACTGGAAGAATTCTCCAAAGCTGTTGGCTATGTCACTAGGCTTATTAATGACTGCCCCACCTACGAGTAGAACTGGAGGAGATGTACTCTTCTTCTTAGATCTGACGTACGCAAAGAAACGCTTTCTGTTGACTGGGAGAGCGAAGAGCTTATTTACATACATCCAGTAGCTTGCATTAATTGCATTCTTTGTATCATTGCGAACAGCTTTGAACTTTTCCCAATCTCGTGGATCATTTGAACACTTTGCATGCTTAAACGAATTGTGCTTTTTGTTGATCAGCTTCTTTATTTCCTTCGTCGTCTTTATTTCCTTCATCGGCACTTGTGTCAAATGAGGCGAAAAGTGACATTCTCCAATTTGCATAAAAGGGCCAGAAGCGCTTTGAGGAGTTTGTATCTGATCGTTTGCTCTCATCATCAACCCTCTCGGTGTGAAGATGAAGAAACTGAAGCTGAAATCCTTTTCGAACTGGATGGAGAAAAGGAAAGTGCGTGTTGGAGACAAGGTCATCAAGTTGCGGGAGGAACGTGAATTGCTGGGAAGATTCCTTATCATCCAAGGCAGCCGCCCAAGCTTAGTTCCTAAACTTGCGGAAACAATAGGTGAATATGAAATGTCAGTGGTTCCCCGTTCACTATGTGCTGTCGATAGTTCTTTGTACATCCCAACAGACAGAGCAAGTCTGATGCATGCAGTTGAAGATGCAAAGGCAGAACCTCCCGAAGCTGTAAAACAGCCTGATGTCGTGGAGGATAATCATTCAACTTCCCCACAGGTTGAAATGATGCAGGAAGAGGTTCCGGTTACCGTAGAACATGAGATAACGCAGGTTCCTTCTGGAGCAGCCCAACAACCTAACCTCTTGCATGATCTTCCTGTCAAGGTGCTGATAATTTACGCCATGGGTGTACTTCAAGGCATGAAAAAGACACCTGCCATGCAAAAGATGTCAGACCTGCAGAATGCATTCAACAGGCGCATCGAAGGGATGATGGCTAGCTATGATGAAGGTCGAGTTGTGTTTGATCGATATATGGAAGAGTCATTAAAGAACAAGACTCGACAAAAGAGAGCAACTACATCAGTAGAATACGAGATACACCCAGAAATGAAGCTCACAATGTCCATTAAGGAGCTTCTATCTTCATCGTCAACCAAGAAAAAGCTGACGTGCCTGTTAGGTCATGGGCTGCTTGATAACTTTTCACAGAACACGGACAGCCCTTTCAAGCTAGAGGTTGTGTATGATACCTTCATCAAGGGACATGATTTTGAAGAGGCGCACACACATGAAGAGGCTGATACCCTGATTCCTAATCAAGTCCTTGCCTCTATACTGCTAGTGGTGCCTTGCAAGAAATTACTGTCTGGTCCCCTGACACTGATGTCCTCCTTCTCTTAATTCATCTAGCATCTTGTGGAAACAATGCAGTGCCCACTTCTCTACAATTTTCTACAGGCAAGGGTACAAAGAAACGAGAAATAGATGTGCTGGAGCGGGTTCAAGTCATTGGACATCAAAAATGTCAAGGCCTCCTTGGACTTCATAATTTTTCTGGTGCTGACTGGGGAGGAAAGTTTGTTGGGATATGCAAGAAGACATGGATCAATGCCTACTTAAAGCTTGATGATGATGATCCTGCCATTATCTGCTTCAAAGAACTGGGTGAGGGTTCCATTCCAACTGAGCTTAGCAATGGAGAGCTTCCAACACAGGTGAAGGCACTAGAGCATTTTATATGTCGTGTATACTGTTCATCAGGTCCAACAACCCTACCATCGCTCAGATGGGAACTGTTCCGATCAAAGAACCTGGAAGGTGAGATGTTACCTCCAACTCGTGCCGCATTACTACCTCACATCCTCCGTGCCAATTACATCACGATGAGAGATAAATCATACAAGACTAACTACCCAGTGCTTCCTCCTATCGAAGAGAATGGATGGTATTCAGATAATCGAGGATGTCTTCCAGTCAAGTGTCTAGCGCTCCCTGCACCACGGGCAGTACTTGAACTCATCAAATGTGGCTGCAAATCAGGATGCAAGGGACGGTGTAGCTGCTCCAACAACGATTTGCCCTGCACTCCTCTCTGTAAATGCTATAGTGGAGACTGTGAAAACCGGACAAGGGAGGATGCTCCATATAGTGACAGCGATTGACTTTAACTGGAAGTCTGATTAACTTTCCTGTTGTCATGGTAACCACGAAGTAGGTGGTTACCCCTTTTCCGATTTTATTTTATGTTTTGTTTATAACTGGAAAATATTTCCAAATTTTTTTGATTTAGATCAAGATTGTGCCATATTTACAAGTATTTCGGCACAAAATAGCTCAATTTGTTTCTGTTGTTTTTTTTTTCTCCGGGGCGGGGGGACCATATCAACATTTTCGTATTATTTTGCATTTGAAAAAAAATATTACATATATAATTTGTCAACTTTCTGATTGACTGATGCAATTACAGTGCCTAGAAGTTGATATTTGAGAATTTTGTTATATTATCATCTGAAGGTGTTTGGCGCCATCTGGTGTTTAAGTAAGTTCTTCATTTGATAGAATACTAAAATAGGGTTGATTTTTATATTAACAGCCTAAAATAACAGTGTTAATGACATGTGATAGGTTTCACAATTATTCCACTGTATAATAACACACAACACACAATAATGCCATTTCCACTTATTATTTCAAAATTCTAGCCCAAAAGAAAGAAGAAGGGCACATTTTTTAGGAAATATTTCCCAGTTTTTCACAATTTTTCGCAAATTGCGCCCCCATAAGGTCAAAGTGAAAAAGTTTTATGACGCGATATCATTAAAAATGAATCCGTATATACCTCAGACAACTTCACAGTGAAATTGGCGAAAAAATATATGAATTAGGGTGGAAGACACCAAATTTCCCTATTGCGCTTTATGCAATTTTGGCGCCCAAAATGGACCCCCCAGGAAGTGAAGGAGGGGTCTCAAAATTTGAACCCCATCATTTTTGGCCAATTGGGACCCTTTTGAAGCCAAAACAGCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCAATTTTGGCAAAAAAATCCTACGGGATTACATACCATACCGCACTATGACCCATAATCGTTCGGTCAGCCTAAACCTGGTACAAGGCGACTTTGACCAAACTTTCAAAAGACCCTTATTGTGCCCAGTCACCCGTGGCAGAAATGTTATGAATAGGGTATCATTGGACAGAAGAGATGTTACTCGCCATCATGAAATAAATTTCAAAGTCATATTGTTTGTGATATTCTAATAAAGAGCTACCAAAGTTGGGTACATTTTTTTGTAAACTCACTGTATACATAGAGAGAGAGAGAGTGGGGGAGGAGAAAAGGGAAGAGGGGGGGGGGGGGAGAGAGAGAGAGAGATCCTGCTACTACTACTGTCGCTAGATTACTGCTATAAAAACTACTATTTCTTAAAAGGCAAGTACACTCCAAAAATACCTTACTTTCAATAGAAATCAGACAATATGTAATCAGACAATTATTTCCCTCACTTTCACTATTATATAAAATTTCTTTTGTCTATCTTCTTCTTCTCTCTCTCCCCCCCTCTCTCTCTCTCTCTCTCTCTCACACACACACACACAAACGTACACACTCACACACACACTCTCTCTCTCTCCCTCTCTCCATCTCACTCTCATCTCTGTCTTTGTCTCTTTGTCTCTCCCTCCCCCCCCCCTCTCTCTCTCTCTCTATATATAATCAGACAATATATAATACCTATATCCCTCACTCTCACTTTCTCACTCCCTCTCGTTTTTTTTTCTTTCTCTCTCTCTCTCTCTCTCTCTATCTATCTCTATCTCTACCTATCTCTCTCTCTCTCTTTCTTTTTCTCTCTCTCTCTCTCTATCAATCTCTATCTCTACCTATCTCTCTCTCTCTCTCTCTCTCTTTCTCTCTCTCTCTCTCTCTCTCTCTTATCTCAATCGATCAATATCGTGATCAATCCAGCTGCATAAGGAAAACAATTATTAGGGGCTATAATATTGAACGCTTCACAAGTATTAGTATGAGTTCAGTTATATATTTTTATTTGACAAAGTCTTGACAAATTCAGTATCAATAAGAGTTTAACAATTCGTCAACATTTCATGCACCTTTCACTGGCTAAGAGTCTCCTGTATCTGAATTTTAGACACAATAATGTTTTAGTCTCTGGTAAGATATGCATTCTGTCCAAGGAAAGACATTAGTGACTATTAGTATATTTACACAAGCAATCTTTTTTTATACAAATTGTACACCTCGGCGACCACATTGACGTCTTCTTCATTGACTGCGATCTCAACCACTTCAGGTGGTGAAGATGTCGTTACGCTTTCCGTCGTTGGATGCTCCTTCTCCTCACTTTCATTCTCACCCCAACGTCCTCGGTTACCATAGGGATGACGGTGGCCATCACGCCTGGGAGATCCTTCTTCGGTTCTGTTGCGTCTGCCGAAGGGACGGTCTCCGAAAGGCTTTCTACCGAAGGGGTTGACCCAGAAGGGCCTCATCTCAAACATTGGTCTGTCCTGATCTCCATCTCCGGTCTTGTTGTGATGACGATGGGGGTGATCTCCCGTCTCATTATGACCCTGGTGACCTTCTGTGTGGTTGCGTCTGCCGAAGTGGTTGAACCGGAAGGGCCTCGTATCATTTGGTTTGTCCTGATCTTGATCTCCCGTCTCATTATGTCCTTGGTGACCCTCTGTCTGGTTATGATGACGGTGGTGGTGGCCTTGGTGGCCATGGTGGCCATGATGGTGAGGGTGAGGACGGCCATCCTCCTCGTTGCGTTCATTGCGCTGACCGAACGGTTGTTCCTCTGGTCGGTCCTGAGGTGGACCTTGACGGTGACCATGATGACGACGTCCATGACCATCAAACTGACCGCGCCCTCCTAGACCATCACCAATTTGTTGGGAAGCATCAGTTTCTTCTTCTCCATCACCACGTCTTCCTCCTTGGCCGAAGAAAGGTCTTCCTCCAGCACCATCTGGTCTGGAGCCACCAAATCCAGGTCCGTCGAATCTCCTACCACCCATCGGACCGCCATTTTGTCTGCGTCCATCCATATGTGGGGCACCAAATCCAGGTCCATCGAACCTCCTTCCACCCATTGGTCCACCATCTTGCCTAGATCCACCCATCTGCATTCCACCAGGCCTTCCTCCAAAGCGACCTTGTCCTCTCTCTCTGCCATTCTCATTTCCTCGCCGTTCATTGAAATCTCTTCGTGCGTGAGCTGTAAGGGTTAATGAAATAATCGAATTAATTAAGTAATCACATTCTGCATTATTATATTTGGCTCTACATTGGTTTATAGCAGGTGACAAATAAGAATGGTAACTCGCCTTCTTATATTTTTTTCCCTTTATGTTTAATCCTAATCCGTAAGTGCGGTTCGATATTCAAAACGCCTAATTCTGACTGAATTAAATTTGATTACCAGATTGAACGTAAACATAATACACATCCTTTGTCATTCAGAAAACGGATTGGTGTAGGCTTAATAAATACAATTAATTATATTGTGTAACGAACAATTTAAAAATGCATATTACAAACTATGAAATAATCATATTTACTTCTGTGAGCCCTTCGTTGGCTTTATATTTAGCACTTATCAAGTAATACCGAGTAATAATTTGATTTCTTACCCGAGATAGCAAGAGCAGCCACAATGGCAACGATCAATGTTGCTTTCACCTCCATGTTTGTAAGGTCTCTCCGATGCTACAAGCTTTCTCTAGATTCGTTGCCTTCCAAGAGAGAACTAGCTCCAAACTTAACAACTACCTGCTGGGCCACTGAATTTATAGGTTTTCTACCCTAGATTGATATCTCACACTAGTACCAGACCGTGTCAGGAATTACCAGGCATACCCTGTTCTAGTTGATCTTTCCCCCTTTTTTTCTTCATTTTCACACTGTGTTAATATCTTTAACCTTTTATCAACCACTGATTTGGAGAGTGCGTTGTACGTACTAGTACCAAAACCAATCTTCGGATATCCGGTAGAGCCATGTATTGGGTGACACCGATGTCGCGGCGTAGTACAGTAGTCACCCTTTGAGTTTGCGAGGGAGCTATAACGGTCATCGGTTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCACTTTCTCTTTCTCTTTCTCTCCCCCTTTCCCCTATCCTTCTATCTGAAAAATTGCATTGATTCGAAAGCCATTTAACAAAGTTATCAAAAAATCTAGGCAATAACAACATCATCATTACAATCATTTTCATTACCATTATCATCATCAATATTCTAGTCATCATTTTCTTCGTTATCTTACAGTCGTCATCATCATCATTACCATTATCATCATTATCATCATCATCATAATCATAATCATCATTCGGAAATTCAAAAGTAATTCAGGACATCCTTTGAAACAGCCCCTCGTTCGTATCAACTCATATTAGTTACGATTCCAATAATTACATAATTATATCCACATTCTGAAATTAACCTAATGGTATATAATGCATTATGAATGAATTGCAATTGCGTCCCTTGAAAATTATTTCGTCTCCCACAATTACAAAATTTTTGTGTAAGAACCGCATTTCTGTTCTGTTTTCAGTACCGAATTTATTAGTCTAAATATGTATTTATCGAAGGGTCGTGATAAATAAACATCCGGGGTAGTGAGCCTTCGAAAAACTGAAAAATCAATAATTCTTTTGACATGGACGGCATATGTAAATCTAGCAAGAGCCTACCTTTACATTTTCCTAGGCTAATTGCCTAATCCAGCGTTTTTCAACCGATGTGCCCCAAGACATTTTTAGGGTACCACAAAAGAATTTGAAGAATATATAAGATTTTAAAAAAAATTCTTGGAATAACGTCGATTTCTTGAAAACCTTAACGACCACTGAAGTACAGTATGTAATTTCATGAGGAAATAGCAATTTACTGGATAGAAAATTGTCTTTCTTCCTTACTCGCTTCACTCCTTCGCAAGACTTAGAAACATGATGATTTCTATCAACGAAGTGCCTTAATATTTGAGCATTGAGGTGCCCGAATAGCCTTCTATTTGACAGTTAGGCGCCCACATAATTATGTGTTATATCTTTTATATGAAAACTTCAATTTTACCACGCTCGTTCCCTTCACCCGCTCGCAACATTTAGGGATAATTCTGGGCAATGAGATGCCTGAATCTCCTCATATTCGTGCTTCGAATTGACTGCACATATTGAGTTTTGTATCAATGTATTCGGCAGCGAGGTGCCTGGACATTTAATTCTATATTTTATCTAGTGTGTTGACTATAAAGTCAATTTTTATCGTATTAAAAGTGTCATGGAGCGGGACGAGTGTCAGAGTGATATCAAATGAATCTTTTGTAAATTAATTTGGAATTAGGGTGACGCGAAACCTTTTTGTTTCTTAAAAAGGTGCCTTGACTAAAACAAAAAGGGTAAAAAAACGCTGGCCTTCACCTTGATCTAATCCAATGAACTATAATTTCTTACAAGAAGAAACAAAAGCTTGTGTTCTTGCTGGAAATGGCTTCCCCAGACTAGCAATTAGTCAAGTAAATGGGGAGGCACGGGCATATCCTCCGCCCAATCGAAAGCCAAAACATTGTTTTATAACAAACAAAACAAAACACAAAAACCCAGAGGTAACCGTTTTTAGATAGAAAAGGAGTCGCAGAGCGGAAAGAGAACAAGAAAGCAGAATGGTTGAGAGCTTGTCACGTGACCAGAAACATTTCATCTATTAAAGTACTACTGTATATGAGTGTACTATTTGGGAGTTAGAAAAAAAATGTGGATTTTTTTTTTTTTTTAATATAATTTTTTCTTTTCTTTTTCAAATCTTAATTTAATGAAAACAAAATATTTCAAATCTTTTTAGATAAAAATTGATATTTTTGGGGTCTAAATTTTTGTTTTCTAAATGGTATTTAAAAAAAATTAAATAAATGTTTGATGCAGGATATAAAAAACAAATATGGACTGCTCTGTAATCGGATGATGGTGAGAACTACCACTCCCTCGGTGGTGTTCACACCGGCCCTTCTATCACCCCCTCGGCTTCGCCTCGGGGTGATAGAGGGGCCGGTGTGAGCACCACCTTGGGAGTGATAGTTCCCACCATCATCCTCATGAGCAGTCAATATTTGTATACTATCTTCTCCCTCTGTCTGTCTCTCTCTGTCTCTGTCTCTGTCTTTCTCTCTCTCTCTCACTCCCATCTCTTTATCTCTTTGTTTCTCTCTGTGTATGTCTCTCTCTCTCTCCCCCTTTCTCTATCTCAGTGATCTATCTCTATCCATCTCTATTTCTGTCTGTCCCTTGCTTTCTCTACCTTTTCTATCTATTTCACCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCTCTCTCTCTTTTCTCTCTTACTCCCATCTTTCTCTGTCTGTCTCCCTCTGCCTTTCTCTCTATCTATCTCTCTTTTACCTCTTCAACCAATCACAGCTGCATAGTGAAAATATAAGGGGCTTTAATATTGAACGCTTCACAAGTATTTGTACGATTTCAATTATATATATTTATTTAACAAGTCTTCAAAAATTCAGTATAAATAAGAGTTTGATAATTCGTCAGCATTTTCATGCACCTTTCACTGGCTTTCAGTCTCTAGCAATCCATGCTGCTGTCTCTGGTAAGATATGCATTCTGTCCAAGGAAAGACATTGATGACTATTATCAATATATTTACACAAGCAATATTTTAAAAACAAATTCTACACCTCGGCGACCACATTGATGTCTTCTTCATTGATTGCGACCACTTCAGGTGGTGAAGATGTCGTTACGCTTTCCGTCGTTGGATGCTCCTTCTCCTCACTTTCATTCTCACCCCATCGTCCTCGGTTACCATAGGGCCGACGATGGCCATCACGCCTGGGAGATCCTTCTTCGGTTCCGTTGCGTCTGTCGAAGGGACGGCCTCCGAAAGGCTTTCTACGGAAGGGGTTGAACCGGAAGGGCCTCATCTCGAACATTGGTCTGTCCTGATCTCCATCTCTGGTCTTGTTGTGATGACGATGGGGGTGGTCTCCCGTCTCATTATGACCCTGGTGACCTTCTGTATGGTTGCGTCTGCCGAAGGGACGGTCTCCGAAAGGCTTTCTACCGAAGTGGTTGTACCGGAAGGGCCTTGTATCATGTAGTTTGTCCTGATCTTGATCTCCCGTCTCATTATGACCTTGGTGACCTTCTGTGTGGTTATGATGACGGTGGTGGTGGCCATGGTGGCGATGGTGGTGAGGGTGAGGACGGCCATCCTCCTCGCTGCTTTCGTTGCGCTGACCAAACGGTTGTTCCTCTGCTTGGTCCTGAGGATGACCTTGACGGTGACCATGATGACGACGTCCAGGACCATCAAACTGACCGGGCCCTCCTAGACCATCACCAATTTGTTGGGCAGCATCAGTTTCTTCTTCTCCATCACCACGCCTGCCTCCTTGGCCGAAGAAAGGTCTTCCTCCAGCACCATCTGGTCTGGAGCCACCAAATCCAGGTCCGTCGAATCTCCTTCCACCCATCGGACCGCCATTTTGTCTCCGTCCATCCATCTCCGGGGCACCAAATCCAGGTCCATCGAACCTCCTTCCACCCATTGGTCCACCATCTTGTCTCCGTCCATCCATTTGTGGGGCACCAGAGTCAGGTCCATCGAACCTTCTCCCACCCATTGGTCCACCATCTTGCCTCGATCCACCCATCTGCATTCCACCCGGCCTTCCTCCAAAGCGACCTTGTCCTCTCTCTCTGCCATTCTTATTTCCTCGTAGTTCATTGTAATCTCTTTGTGCGTGAGCTGTAGGGGATAATGAAGAAAAAAATCGAATTAATTAAGTGATCACATTCTGCATTATTATCTTTAGCTCTACATTGGGTTATGGCAGGTGACAAATAAGAATGGTAACTCGCCTTCTTTGATGTTTTTTCCCTTTGTGTTTAACCCTATCGGTTCGATATTCAAAACGCCTAATTCTGACTGAATTGAACTTGAATACCAGATTGAACATAAAGATAATTCTTATCTTTTGCCACTCGCCCAACGGATTGGTGTAGGCTTAAGAAATAATATAAATTATATTTTGTAACGAACAATTTAGAAATGCATATTATAAAATATGAAATAATCATATATACTTCTGTGAGCTATTCGTTGGCTCTATATTTAGCACTTATCAAGTAATACCGAGTAATAATTTGATTTCTTACCTGAGATAGCAAGAGCAGCCACAATGGCAACGATCAGTGTCACTTTCACCATGTTAGTAATAGGTCTCTCCGATGCTATAGCTTTCTCTAGATTCGTTGTCTTCTAAGAGAGAACGAGCTCCAAACTGAATAACTATCTGCTGGATCACTGAATTTATAGGCTTTCTATCATAGCTTGAAATCACTTTCGTACCAGTCCGTGTCAGGAATTACCAGGGATACCCTGGTCCAGTTGATCTCCCCCCTTTTTATCTTTATTTTCATACTTTGTTAATATCTTTAACCTTTTATCAACCACTAATTTGGAGAGTGCGTTGTACGTACTAGTACCAAAACCAATCACTTCATCTTCGGATATCCGGTAGAGCCATGTATTGGGTGACACCGATGTCGCGGCGTAGTACAGTAGTCACCCTTTGAGTTTGCGAGGGAGCTAACGGTCATCGGTTTCTCTCTCTCTTTCTCAGTCTGCCTGCCTGCCTCCCTGTCCGTCTGAGTGTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCTTCTCCTTCTCCTTCTCTCCCCCGATCCTTCTATCTGAAAAATTGCATTGATTCGAAAGCCATTTAACAATGTTGTCAAAGAATCTAGGCAATAACAACATCATCATTACAATCGTTATCATTACCATTATCATCATCAATATTCTAGTCATCATTTTCTTCGTTATCTTACAGTCGTCATCATCATCATCATACCAATCATAATCATCACCATTATCAGCATCATCATCATCATCATTACCATTATCATCACTATCATCACCATCATCATTATCATCATCATCATAATCATTTGGAAATTAAAAAGTAATCCAGGACATCCGTTGAAACAACCCCTCGTTCGTATCAACTCATATTAGTTACAATTCCAATGATTACATAATTATATCCACATTCTGAAGTTAACCTAATGGTATATAATGCATTGTGAATGAATTGCAATTGCGTCCCTTGAAAATTATTTCGTCTCCCACAATTACAAAATTTTTGTGTAAGAACCGCATTTCTGTTCTGTTTTCAGTACCGAATTTATTAGTCTAAATATGTATTTATCGAAGGGTCGTGATAAATAAACATCCGGGGTAGTGAGCCTTCGAAAAACTGAAAAATCAATAATTCTTTTGACATGGACGGTATATGTAAATCTAGCAAGAGCCTACCTTTACATTTTCCTAGGCTAATTGCCTAATCCAGCGTTTTTCAACCGATGTGCCCCAAGACATTTTTAGGGTACCACAAAAGAATTTGAAGAATATATAAGATTTAAAAAAAAATTCTTGGAATAACGTCGATTTCTTGAAAACCTTAACGACCACTGAAGTACAGTATGTAATTTCATGAGGAAATAGCAATTTACTGGATAGAAAATTGTCTTTCTTCCTTACTCGCTTCACTCCTTCGCAAGACTTAGAAACATGATGATTTCTATCAACGAAGTGCCTTAATATTTGAGCATTGAGGTGCCCGAATAGCCTTCTATTTGACAGTTAGGCGCCCACATAATTATGTGTTATA";
//        bwt_.mapSequence( mn->seq, mn->ids, mn->coords );
//        mn->recoil();
//        Node* edge = new Node( mn, 0, mn->ids.size()-1, 1 );
//        edge->dontExtend_ = true;
//        base->addEdge( edge, 249, 1 );
//        delete mn;
//        nodes_.push_back( edge );
////        int len = edge->seq_.length();
////        ofstream fp( "/media/glen/ssd/test.fa" );
////        for ( ReadMark &mark : edge->marks_[0] )
////        {
////            int diff = mark.estimate - edge->ends_[1] + 150;
////            if ( !params.isReadPe( mark.id ) ) continue;
////            if ( diff < 0 ) continue;
////            fp << ">" << mark.id << endl << string( diff, '-' ) + bwt_.getSequence( mark.id ) << endl;
////            
////        }
////        fp.close();
////        exit( 0 );
//    }
//    
////    params.libs[1].size += 300;
//    
//    NodeList island;
//    ExtVars ev( nodes_, island, validLimits_, bwt_, false, false );
//    IslandVars iv( ev, 1 );
//    bool finished = false;
//    bool ended = false;
//    bool drxn = true;
//    int z = 0;
//    while ( !finished )
//    {
//        finished = true;
//        {
//            int minDiv = 0, maxDiv = 0;
//            bool didExtend = true;
//            while ( minDiv < 20 && didExtend )
//            {
//                minDiv = maxDiv;
//                NodeList extNodes;
//                vector<int> divCounts;
//                for ( Node* node : nodes_ )
//                {
//                    int thisDiv = 0;
//                    if ( !node->isContinue( drxn ) ) continue;
//                    extNodes.push_back( node );
//                    for ( Node* bck : node->getDrxnNodes( !drxn, false, false ) ) 
//                    {
//                        thisDiv += bck->edges_[drxn].size() - 1;
//                        thisDiv -= bck->edges_[!drxn].size() - 1;
//                    }
//                    divCounts.push_back( thisDiv );
//                    maxDiv = max( maxDiv, thisDiv );
//                    minDiv = min( minDiv, thisDiv );
//                }
//                didExtend = false;
//                for ( int i = 0; i < extNodes.size(); i++ )
//                {
//                    if ( divCounts[i] > minDiv + 5 ) continue;
//                    finished = false;
//                    didExtend = true;
//                    ev.ante = extNodes[i]->getDrxnNodes( !drxn );
//                    extNodes[i]->extendCount_ = 99;
//                    extNodes[i]->extendNode( ev, drxn );
//                    ends_[drxn] = drxn ? max( ends_[1], extNodes[i]->ends_[1] ) : min( ends_[0], extNodes[i]->ends_[0] );
//                }
//            }
//        }
//        {
//            NodeList endNodes;
//            for ( Node* node : nodes_ )
//            {
//                if ( node->drxn_ != 2 ) continue;
//                for ( Node* nxt : node->getNextNodes( drxn ) )
//                {
//                    if ( nxt->drxn_ != 2 )
//                    {
//                        endNodes.push_back( node );
//                        break;
//                    }
//                }
//            }
//            
//            NodeSet delSet;
//            for ( Node* node : endNodes )
//            {
//                ended = node->plotSeed( iv, delSet, drxn, finished );
//            }
//            for ( Node* del : delSet )
//            {
//                nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
//                delete del;
//            }
//        }
//        {
//            NodeSet delSet;
//            NodeList nodes[3];
//            for ( Node* node : nodes_ )
//            {
//                nodes[node->drxn_].push_back( node );
//            }
//            Node::mergeAll( nodes, delSet );
//            for ( Node* del : delSet )
//            {
//                nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
//                delete del;
//            }
//            if ( finished )
//            {
//                {
//                    ofstream fp( "/media/glen/ssd/test" + to_string( z++ ) + ".fa" );
//                    for ( int c = 0; c < nodes_.size(); c++ )
//                    {
//                        fp << ">" + to_string( c ) << endl << string( nodes_[c]->ends_[0], '-' )+ nodes_[c]->seq_ << endl;
//                    }
//                    fp.close();
//                }
//                nodes[drxn].clear();
//                for ( Node* node : nodes_ )
//                {
//                    if ( node->drxn_ == 2 ) continue;
//                    if ( !node->dontExtend_ ) node->remap( bwt_ );
//                    node->dontExtend_ = false;
//                    nodes[drxn].push_back( node );
//                }
//                delSet.clear();
//                Node::mergeAll( nodes, delSet );
//                for ( Node* del : delSet )
//                {
//                    nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
//                    delete del;
//                }
//                Node* seed = nodes_[0]->edges_[1][0].node;
//                while ( !seed->edges_[drxn].empty() )
//                {
//                    assert( seed->edges_[drxn].size() == 1 );
//                    seed = seed->edges_[drxn][0].node;
//                }
//                finished = ended || !seed->fixSeed( iv, drxn );
//                seed->dontExtend_ = true;
//                {
//                    ofstream fp( "/media/glen/ssd/test" + to_string( z++ ) + ".fa" );
//                    for ( int c = 0; c < nodes_.size(); c++ )
//                    {
//                        fp << ">" + to_string( c ) << endl << string( nodes_[c]->ends_[0], '-' )+ nodes_[c]->seq_ << endl;
//                    }
//                    fp.close();
//                }
//            }
//        }
//    }
//    
//    {
//        NodeSet delSet;
//        NodeList nodes[3];
//        for ( Node* node : nodes_ )
//        {
//            nodes[node->drxn_].push_back( node );
//        }
//        Node::mergeAll( nodes, delSet );
//        for ( Node* del : delSet )
//        {
//            nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
//            delete del;
//        }
//    }
//    
//    ofstream fp( "/media/glen/ssd/SpBACs/Clone42.fa" );
//    fp << ">Clone42" << endl << nodes_[1]->seq_ << endl;
//    fp.close();
//    exit(0);
//}

void Seed::checkDivergent( NodeList &path )
{
    NodeSet drxnSets[2] = {  path.back()->getDrxnNodes( 0, false, true ), path[0]->getDrxnNodes( 1, false, true ) };
    NodeSet pathSet = path[0]->getDrxnNodesInSet( drxnSets[0], 1, true );
    NodeList divNodes[2];
    
    for ( bool drxn : { 0, 1 } )
    {
        for ( Node* node : pathSet )
        {
            for ( Node* nxt : node->getNextNodes( drxn ) )
            {
                if ( pathSet.find( nxt ) == pathSet.end() )
                {
                    divNodes[drxn].push_back( nxt );
                }
            }
        }
        
        checkDivergentBack( divNodes[drxn], drxnSets[drxn], drxn );
    }
    
    NodeIntMap divScores;
    for ( bool drxn : { 0, 1 } )
    {
        for ( Node* div : divNodes[drxn] )
        {
            int score = div->getPairHitsTotal();
            for ( Node* fwd : div->getDrxnNodes( drxn ) )
            {
                score += fwd->getPairHitsTotal();
            }
            assert( score < 2 );
            divScores[div] = score;
        }
    }
}

void Seed::checkDivergentBack( NodeList &div, NodeSet &pathSet, bool drxn )
{
    NodeSet delSet;
    for ( int i ( 0 ); i < div.size(); )
    {
        NodeSet fwdSet = div[i]->getDrxnNodes( drxn, false, true );
        NodeList branches;
        for ( Node* fwd : fwdSet )
        {
            for ( Node* prv : fwd->getNextNodes( !drxn ) )
            {
                if ( find( branches.begin(), branches.end(), prv ) == branches.end()
                        && pathSet.find( fwd ) == pathSet.end() )
                {
                    branches.push_back( prv );
                }
            }
        }
        
        for ( Node* node : branches )
        {
            NodeSet tSet = node->getDrxnNodesInSet( pathSet, drxn );
            NodeSet qSet = node->getDrxnNodesNotInSet( pathSet, !drxn, true );
            int score[2] = { 0, 0 };
            for ( Node* t : tSet )
            {
                for ( auto &np : t->pairs_ )
                {
                    if ( qSet.find( np.first ) != qSet.end() )
                    {
                        score[0] += np.second;
                    }
                    else if ( tSet.find( np.first ) == tSet.end() && pathSet.find( np.first ) != pathSet.end() )
                    {
                        score[1] += np.second;
                    }
                }
            }
            
            if ( score[0] > score[1] && score[0] > 2 )
            {
                for ( Node* nxt : node->getNextNodes( drxn ) )
                {
                    for ( Node* prv : nxt->getNextNodes( !drxn ) )
                    {
                        if ( pathSet.find( nxt ) != pathSet.end() && pathSet.find( prv ) != pathSet.end() )
                        {
                            nxt->removeEdge( prv, !drxn );
                            prv->removeEdge( nxt, drxn );
                            if ( prv->edges_[drxn].empty() )
                            {
                                prv->dismantleNode( delSet, !drxn );
                            }
                        }
                    }
                }
            }
            else
            {
                for ( Node* nxt : node->getNextNodes( drxn ) )
                {
                    if ( pathSet.find( nxt ) != pathSet.end() )
                    {
                        node->removeEdge( nxt, drxn );
                        nxt->removeEdge( node, !drxn );
                    }
                }
                
                if ( node->edges_[drxn].empty() )
                {
                    node->dismantleNode( delSet, !drxn );
                }
            }
        }
        
        i++;
    }
    
    for ( Node* del : delSet )
    {
        delete del;
        div.erase( remove( div.begin(), div.end(), del ), div.end() );
        nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
    }
}

void Seed::deleteNodes( NodeSet &delSet )
{
    for ( Node* del : delSet )
    {
        nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
        delete del;
    }
    delSet.clear();
}

vector<Locus*> Seed::getLoci()
{
    resolveBackForks();
    vector<Locus*> loci;
    
    NodeIntMap scores;
    NodeList bestPath = getLociGetPath( loci.empty() );
    if ( bestPath.empty() ) return loci;
//    checkDivergent( bestPath );
    
    Node* forks[2];
    while ( getLociSetDivergent( bestPath, forks ) )
    {
        getLociResolveDivergent( scores, bestPath, forks );
    }
    
    while ( getLociSetConvergent( bestPath, forks ) )
    {
        if ( abs( forks[0]->ends_[0] + forks[0]->getBestOverlap( 1 ) ) < abs( forks[1]->ends_[1] - forks[1]->getBestOverlap( 0 ) ) )
        {
            bestPath.erase( bestPath.begin(), find( bestPath.begin(), bestPath.end(), forks[0] ) );
        }
        else
        {
            bestPath.erase( find( bestPath.begin(), bestPath.end(), forks[1] ) + 1, bestPath.end() );
        }
    }
    
    assert( !getLociSetDivergent( bestPath, forks ) );
    
    if ( !bestPath.empty() )
    {
        forks[0] = forks[0] ? forks[0] : bestPath[0];
        forks[1] = forks[1] ? forks[1] : bestPath.back();
        NodeList nodes[3];
        NodeList forkList( find( bestPath.begin(), bestPath.end(), forks[0] ), find( bestPath.begin(), bestPath.end(), forks[1] ) + 1 );
        
        Node* origin = NULL;
        int32_t bestTether = 0;
        for ( Node* node : forkList )
        {
            int32_t nodeTether = min( tether_[0] - node->ends_[0], node->ends_[1] - tether_[1] );
            if ( !origin || nodeTether > bestTether )
            {
                origin = node;
                bestTether = nodeTether;
            }
        }
        
        NodeSet locusSet = { origin };
        nodes[2].push_back( origin );
        for ( int drxn : { 0, 1 } )
        {
            for ( Node* node : origin->getDrxnNodes( drxn ) )
            {
                node->drxn_ = drxn;
                nodes[drxn].push_back( node );
                locusSet.insert( node );
            }
            origin->offsetForward( drxn, false, true );
        }
        
        Locus* locus = new Locus( bwt_, nodes );
        locus->header_ = header_;
        loci.push_back( locus );
        
        for ( Node* node : locusSet )
        {
            nodes_.erase( remove( nodes_.begin(), nodes_.end(), node ), nodes_.end() );
        }
    }
    
    for ( Node* node : nodes_ )
    {
        node->dismantleNode();
        delete node;
    }
    
    nodes_.clear();
    
    return loci;
}

NodeList Seed::getLociGetPath( bool doForce )
{
    NodeList startNodes = getLociGetPathGetStarts( doForce );
    NodeList bestPath;
    int bestPathScore = -1;
    
    for ( Node* node : startNodes )
    {
        NodeList path = { node };
        int score = getLociGetPathCurrScore( node, path, 1 );
        getLociGetPath( node, path, score, 1 );
        
        if ( score > bestPathScore )
        {
            bestPathScore = score;
            bestPath = path;
        }
    }
    
    if ( !bestPath.empty() && !bestPath[0]->edges_[0].empty() )
    {
        reverse( bestPath.begin(), bestPath.end() );
        getLociGetPath( bestPath.back(), bestPath, bestPathScore, 0 );
        reverse( bestPath.begin(), bestPath.end() );
    }
    
    return bestPath;
}

void Seed::getLociGetPath( Node* curr, NodeList &path, int &score, bool drxn )
{
    while ( !curr->edges_[drxn].empty() )
    {
        Node* bestNxt = NULL;
        int bestNxtScore = -1;
        for ( Node* nxt : curr->getNextNodes( drxn ) )
        {
            int nxtScore = 0;
            for ( Node* fwd : nxt->getDrxnNodes( drxn, false, true ) )
            {
                for ( Node* node : path )
                {
                    nxtScore += fwd->getPairHits( node );
                }
            }

            if ( nxtScore > bestNxtScore )
            {
                bestNxt = nxt;
                bestNxtScore = nxtScore;
            }
        }

        assert( bestNxt );
        curr = bestNxt;
        path.push_back( curr );
        score += getLociGetPathCurrScore( curr, path, drxn );
    }
    
}

int Seed::getLociGetPathCurrScore( Node* curr, NodeList &path, bool drxn )
{
    int score = 0;
    
    if ( drxn ? tether_[0] <= curr->ends_[0] || tether_[1] <= curr->ends_[1] 
              : curr->ends_[0] <= tether_[0] || curr->ends_[1] <= tether_[1] )
    {
        for ( Node* node : path )
        {
            if ( drxn ? node->ends_[0] <= tether_[0] || node->ends_[1] <= tether_[1]
                      : tether_[0] <= node->ends_[0] || tether_[1] <= curr->ends_[1] )
            {
                score += node->getPairHits( curr );
            }
        }
    }
    
    return score;
}

NodeList Seed::getLociGetPathGetStarts( bool doForce )
{
    NodeList spanNodes, startNodes, endNodes;
    NodeSet spanFwd;
    for ( Node* node : nodes_ )
    {
        if ( node->ends_[0] <= tether_[0] || node->ends_[1] <= tether_[1] )
        {
            bool isSpan = false;
            for ( auto &np : node->pairs_ )
            {
                isSpan = isSpan || tether_[0] <= np.first->ends_[0] || tether_[1] <= np.first->ends_[1];
            }
            if ( isSpan )
            {
                spanNodes.push_back( node );
                node->getDrxnNodes( spanFwd, 1 );
            }
        }
        
        if ( node->edges_[0].empty() )
        {
            endNodes.push_back( node );
        }
    }
    
    for ( Node* node : spanNodes )
    {
        if ( spanFwd.find( node ) == spanFwd.end() )
        {
            startNodes.push_back( node );
        }
    }
    
    return ( startNodes.empty() && doForce ? endNodes : startNodes );
}

void Seed::getLociResolveDivergent( NodeIntMap &scores, NodeList &path, Node** forks )
{
    NodeSet pathSet( path.begin(), path.end() );
    int divScores[2] = { 0, 0 };
    NodeSet divSets[2];
    
    for ( int drxn : { 0, 1 } )
    {
        divSets[drxn] = forks[drxn]->getDrxnNodesNotInSet( pathSet, drxn );
        for ( Node* fwd : divSets[drxn] )
        {
            divScores[drxn] += scores[fwd];
        }
    }
    
    NodeSet delSet;
    for ( int drxn : { 0, 1 } )
    {
        if ( divScores[drxn] < divScores[!drxn] || divScores[drxn] < 2 )
        {
            for ( Node* node : divSets[drxn] )
            {
                node->dismantleNode( delSet, drxn );
            }
        }
    }
    
    for ( Node* del : delSet )
    {
        delete del;
        nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
    }
}

bool Seed::getLociSetConvergent( NodeList &path, Node** forks )
{
    forks[0] = NULL;
    forks[1] = NULL;
    NodeSet pathSet = { path.begin(), path.end() };
    
    for ( Node* node : path )
    {
        if ( node->edges_[0].size() > 1 )
        {
            forks[0] = node;
        }
        if ( node->edges_[1].size() > 1 && !forks[1] )
        {
            forks[1] = node;
        }
    }
    
    if ( forks[0] && forks[1] )
    {
        NodeSet fwdSet = forks[1]->getDrxnNodesNotInSet( pathSet, 1 );
        NodeSet bckSet = forks[0]->getDrxnNodesInSet( fwdSet, 0 );
        return !bckSet.empty();
    }
    return false;
    
}

bool Seed::getLociSetDivergent( NodeList &path, Node** forks )
{
    forks[0] = NULL;
    forks[1] = NULL;
    bool problem = false;
    NodeSet pathSet = { path.begin(), path.end() };
    NodeSet drxnSets[2] = { path.back()->getDrxnNodes( 0 ), path[0]->getDrxnNodes( 1 ) };
    
    for ( Node* node : path )
    {
        for ( int drxn : { 0, 1 } )
        {
            for ( Node* nxt : node->getNextNodes( drxn ) )
            {
                if ( pathSet.find( nxt ) == pathSet.end()
                        && drxnSets[!drxn].find( nxt ) == drxnSets[!drxn].end() )
                {
                    problem = problem || ( !drxn && forks[1] );
                    forks[drxn] = !drxn || !forks[drxn] ? node : forks[drxn];
                }
            }
        }
    }
    
    return problem;
}

void Seed::merge()
{
    NodeSet delSet;
    NodeList nodes[3];
    for ( Node* node : nodes_ )
    {
        assert( node->drxn_ < 3 );
        nodes[node->drxn_].push_back( node );
    }
    Node::mergeAll( nodes, delSet );
    deleteNodes( delSet );
}

//void Seed::plot()
//{
//    NodeList ends[2];
//    for ( bool drxn : { 0, 1 } )
//    {
//        vector< pair<Node*, Node*> > forks;
//        for ( Node* node : nodes_ )
//        {
//            if ( !node->edges_[drxn].empty() ) continue;
//            Node* curr = node;
//            while ( curr->edges_[!drxn].size() == 1 && curr->edges_[drxn].size() < 2 ) curr = node->edges_[!drxn][0].node;
//            forks.push_back( make_pair( node, curr ) );
//        }
//        for ( int i = 0; i < forks.size(); i++ )
//        {
//            bool useFork = false;
//            for ( int j = i + 1; j < forks.size(); )
//            {
//                if ( forks[i].second == forks[j].second )
//                {
//                    useFork = true;
//                    forks.erase( forks.begin() + j );
//                }
//                else j++;
//            }
//            ends[drxn].push_back( ( useFork ? forks[i].second : forks[i].first ) );
//        }
//    }
//    
//    Node::plotSeed( ends[0], 1 );
//}

bool Seed::resolveBackFork( Node** forks, NodeSet &delSet )
{
    NodeSet qSets[2], tSets[2], forkSets[2];
    for ( bool drxn : { 0, 1 } )
    {
        forkSets[!drxn] = forks[!drxn]->getDrxnNodes( drxn, false, true );
        forks[!drxn]->getDrxnNodes( forkSets[!drxn], !drxn );
        qSets[drxn] = forks[drxn]->getDrxnNodesNotInSet( forkSets[!drxn], !drxn );
        tSets[drxn] = forks[drxn]->getDrxnNodes( drxn, false, true );
    }
    
    int midScore = 0;
    for ( Node* fwd : forks[0]->getDrxnNodesNotInSet( qSets[1], 1 ) )
    {
        bool isMid = tSets[1].find( fwd ) != tSets[1].end();
        for ( auto &np: fwd->pairs_ )
        {
            if ( tSets[0].find( np.first ) != tSets[0].end() || ( isMid && tSets[1].find( np.first ) != tSets[1].end() ) )
            {
                midScore += np.second;
            }
        }
    }
    
    if ( !qSets[0].empty() || !qSets[1].empty() )
    {
        int scores[2] = { 0, 0 };
        for ( bool drxn : { 0, 1 } )
        {
            for ( Node* q : qSets[drxn] )
            {
                for ( auto &np : q->pairs_ )
                {
                    if ( tSets[drxn].find( np.first ) != tSets[drxn].end() )
                    {
                        scores[drxn] += np.second;
                    }
                }
            }
        }
        
        bool delMid = midScore < scores[0] && midScore < scores[1];
        for ( bool drxn : { 0, 1 } )
        {
            NodeList midNodes;
            bool delQ = !delMid && scores[drxn] <= scores[!drxn];
            for ( Node* nxt : forks[drxn]->getNextNodes( !drxn ) )
            {
                bool isQ = qSets[drxn].find( nxt ) != qSets[drxn].end();
                if ( ( isQ && delQ ) || ( !isQ && delMid ) )
                {
                    forks[drxn]->removeEdge( nxt, !drxn );
                    nxt->removeEdge( forks[drxn], drxn );
                    if ( nxt->edges_[!drxn].empty() )
                    {
                        nxt->dismantleNode( delSet, drxn );
                    }
                }
                if ( !isQ ) midNodes.push_back( nxt );
            }
            if ( !delQ && !delMid && midNodes.size() > 1 )
            {
                Node* worst = NULL;
                int worstScore = 0;
                for ( Node* node : midNodes )
                {
                    int score = 0;
                    NodeSet midSet = node->getDrxnNodesNotInSet( tSets[!drxn], !drxn, true );
                    for ( Node* mid : midSet ) score += mid->getPairHitsTotal();
                    if ( !worst || score < worstScore )
                    {
                        worst = node;
                        worstScore = score;
                    }
                }
                worst->dismantleNode( delSet, !drxn );
            }
        }
        
        return true;
    }
    return false;
}

bool Seed::resolveBackForkBypass( Node* fork, NodeSet &delSet )
{
    NodeSet fwdSet = fork->getDrxnNodes( 1 );
    NodeSet tSet = fork->getDrxnNodes( 0, false, true );
    bool didErase = true;
    
    while ( didErase && fork->edges_[1].size() > 1 )
    {
        didErase = false;
        for ( Node* nxt : fork->getNextNodes( 1 ) )
        {
            NodeSet qSet = nxt->getDrxnNodesInSet( fwdSet, 0 );
            if ( !qSet.empty() )
            {
                int fwdScore = 0, qScore = 0, qReads = 0;

                for ( Node* q : qSet )
                {
                    qReads += q->reads_.size();
                    for ( auto &np : q->pairs_ )
                    {
                        if ( tSet.find( np.first ) != tSet.end() || fwdSet.find( np.first ) != fwdSet.end() )
                        {
                            qScore += np.second;
                        }
                    }
                }

                for ( Node* fwd : fwdSet )
                {
                    for ( auto &np : fwd->pairs_ )
                    {
                        if ( tSet.find( np.first ) != tSet.end() )
                        {
                            fwdScore += np.second;
                        }
                    }
                }
                
                if ( qScore == 0 && fwdScore > 0 && qReads > params.peCover )
                {
                    assert( false );
                    for ( Node* q : fork->getNextNodes( 1 ) )
                    {
                        if ( qSet.find( q ) != qSet.end() )
                        {
                            fork->removeEdge( q, 1 );
                            q->removeEdge( fork, 0 );
                            if ( q->edges_[0].empty() )
                            {
                                q->dismantleNode( delSet, 1 );
                            }
                        }
                    }
                }
                else
                {
                    fork->removeEdge( nxt, 1 );
                    nxt->removeEdge( fork, 0 );
                }
                
                didErase = true;
                break;
            }
        }
    }
    
    return fork->edges_[1].size() <= 1;
}

bool Seed::resolveBackForkDouble( Node* fork, NodeSet &delSet )
{
    bool didResolve = true;
    while ( didResolve && fork->edges_[1].size() > 1 )
    {
        didResolve = false;
        Node* nodes[4] = { fork, NULL, NULL, NULL };
        for ( int i ( 0 ); i < fork->edges_[1].size() - 1 && !didResolve; i++ )
        {
            nodes[2] = fork->edges_[1][i].node;
            NodeSet prvSet = nodes[2]->getNextNodes( 0 );
            for ( int j ( i + 1 ); j < fork->edges_[1].size() && !didResolve; j++ )
            {
                nodes[3] = fork->edges_[1][j].node;
                for ( Node* prv : nodes[3]->getNextNodes( 0 ) )
                {
                    if ( prv != nodes[0] && prvSet.find( prv ) != prvSet.end() )
                    {
                        nodes[1] = prv;
                        resolveBackForkDouble( nodes, delSet );
                        didResolve = true;
                    }
                }
            }
        }
    }
    
    return fork->edges_[1].size() <= 1;
}

void Seed::resolveBackForkDouble( Node** forks, NodeSet &delSet )
{
    NodeSetList bckSets = Node::getNodeSetsExclusive( forks[0], forks[1], 0 );
    NodeSetList fwdSets = Node::getNodeSetsExclusive( forks[2], forks[3], 1 );
    int pairScores[2][2];
    int pref[2] = { -1, -1 };
    for ( int i : { 0, 1 } )
    {
        for ( int j : { 0, 1 } )
        {
            pairScores[i][j] = 0;
            for ( Node* bck : bckSets[i] )
            {
                for ( Node* fwd : fwdSets[j] )
                {
                    auto it = bck->pairs_.find( fwd );
                    if ( it != bck->pairs_.end() )
                    {
                        pairScores[i][j] += it->second;
                    }
                }
            }
        }
    }
    
    bool didResolve = false;
    for ( int i ( 0 ); i < 2 && !didResolve; i++ )
    {
        for ( int j ( 0 ); j < 2 && !didResolve; j++ )
        {
            if ( pairScores[i][j] > pairScores[i][!j] )
            {
                if ( pairScores[!i][j] <= pairScores[!i][!j] )
                {
                    forks[i]->removeEdge( forks[2+!j], 1 );
                    forks[2+!j]->removeEdge( forks[i], 0 );
                    forks[!i]->removeEdge( forks[2+j], 1 );
                    forks[2+j]->removeEdge( forks[!i], 0 );
                    didResolve = true;
                }
            }
        }
    }
    
    if ( !didResolve )
    {
        int minPairs, minReads, iMinPairs = -1, iMinReads = -1;
        bool anyPairs = false;
        
        for ( int i ( 0 ); i < 4; i++ )
        {
            int pairCount = 0, readCount = 0;
            
            for ( Node* node : ( i < 2 ? bckSets[i] : fwdSets[i-2] ) )
            {
                pairCount += node->getPairHitsTotal();
                readCount += node->reads_.size();
            }
            
            anyPairs = anyPairs || pairCount > 0;
            
            if ( iMinPairs = -1 || pairCount < minPairs )
            {
                iMinPairs = i;
                minPairs = pairCount;
            }
            
            if ( iMinReads = -1 || readCount < minReads )
            {
                iMinReads = i;
                minReads = readCount;
            }
        }
        
        int iMinSet = ( anyPairs ? iMinPairs : iMinReads );
        if ( iMinSet == -1 ) assert( false );
        for ( Node* node : ( iMinSet < 2 ? bckSets[iMinSet] : fwdSets[iMinSet-2] ) )
        {
            node->dismantleNode();
            delSet.insert( node );
        }
    }
}

void Seed::resolveBackForks()
{
    NodeList forks;
    NodeSet delSet;
    for ( Node* node : nodes_ )
    {
        if ( node->edges_[1].size() > 1 
                && !resolveBackForkBypass( node, delSet )
                && !resolveBackForkDouble( node, delSet ) )
        {
            forks.push_back( node );
        }
    }
    
    for ( Node* del : delSet )
    {
        nodes_.erase( remove( nodes_.begin(), nodes_.end(), del ), nodes_.end() );
        forks.erase( remove( forks.begin(), forks.end(), del ), forks.end() );
    }
    
    for ( Node* fork : forks )
    {
        bool didResolve = true;
        while ( didResolve && fork->edges_[1].size() > 1 && delSet.find( fork ) == delSet.end() )
        {
            didResolve = false;
            
            for ( Node* branch : fork->getNextNodes( 1 ) )
            {
                Node* forkPair[2] = { fork, branch };

                while ( forkPair[1] && !didResolve )
                {
                    Node* nxt = NULL;
                    if ( forkPair[1]->edges_[0].size() > 1 )
                    {
                        didResolve = resolveBackFork( forkPair, delSet );
                    }
                    else if ( forkPair[1]->edges_[1].size() == 1 )
                    {
                        nxt = forkPair[1]->edges_[1][0].node;
                    }
                    forkPair[1] = nxt;
                }
            }
        }
    }
}

bool Seed::warning()
{
    if ( readCount_ > params.cover * 3 )
    {
        cout << "\tWarning: excessive reads mapped ( " << to_string( readCount_ ) << " )" << endl;
        cout << "\tEstimated multiplicity of " << to_string( readCount_ / params.cover ) << "." << endl;
        cout << "\tExtension will not be attempted. Please try a more refined query." << endl << endl;
        return true;
    }
    return false;
}

Seed::~Seed()
{
    for ( Node* node : nodes_ )
    {
        delete node;
    }
    nodes_.clear();
}

