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
        
//    vector<bool> redundant;
//    for ( int i = 0; i < ms.reads.size(); i++ )
//    {
//        bool isRedundant = false;
//        int iLen = ms.reads[i].coords[1] - ms.reads[i].coords[0];
//        for ( int j = 0; j < ms.reads.size(); j++ )
//        {
//            int jLen = ms.reads[j].coords[1] - ms.reads[j].coords[0];
//            if ( iLen >= params.readLen || isRedundant || jLen <= iLen ) continue;
//            isRedundant = ms.reads[j].seq.find( ms.reads[i].seq ) != ms.reads[j].seq.npos;
//        }
//        redundant.push_back( isRedundant );
//    }
//    
//    int k = 0;
//    for ( ReadStruct &read : ms.reads )
//    {
//        if ( redundant[k++] ) continue;
//        vector<int32_t> hits, coords;
//        NodeSet fwdSet;
//        tether_[0] = min( tether_[0], read.tether[0] );
//        tether_[1] = max( tether_[1], read.tether[1] );
//        
//        for ( int i( 0 ); i < nodes_.size(); i++ )
//        {
//            int32_t coord;
//            if ( nodes_[i]->seedCongruent( read, coord ) )
//            {
//                hits.push_back( i );
//                coords.push_back( coord );
//                nodes_[i]->getDrxnNodes( fwdSet, 0 );
//            }
//        }
//        
//        for ( int i ( 0 ); i < hits.size() && hits.size() > 1; )
//        {
//            if ( fwdSet.find( nodes_[hits[i]] ) != fwdSet.end() )
//            {
//                hits.erase( hits.begin() + i );
//                coords.erase( coords.begin() + i );
//            }
//            else
//            {
//                i++;
//            }
//        }
//        
//        if ( hits.size() > 1 || hits.size() == 1 && ( nodes_[ hits[0] ]->ends_[1] != coords[0] || !nodes_[ hits[0] ]->edges_[1].empty() ) )
//        {
//            Node* node = new Node( read );
//            nodes_.push_back( node );
//            for ( int i ( 0 ); i < hits.size(); i++ )
//            {
//                if ( coords[i] != nodes_[ hits[i] ]->ends_[1] )
//                {
//                    nodes_[ hits[i] ]->seedSplit( nodes_, coords[i] );
//                }
//                node->addEdge( nodes_[ hits[i] ], 0 );
//            }
//        }
//        else if ( hits.size() == 1 )
//        {
//            nodes_[ hits[0] ]->seedAdd( read );
//        }
//        else
//        {
//            Node* node = new Node( read );
//            nodes_.push_back( node );
//        }
//    }
//    
//    bool anyGood = false;
//    for ( Node* node : nodes_ ) if ( node->reads_.size() > 2 ) anyGood = true;
//    
//    for ( int i = 0; i < nodes_.size(); )
//    {
//        bool doErase = anyGood 
//                &&( nodes_[i]->seq_.length() < params.readLen || nodes_[i]->reads_.size() <= 2 )
//                && nodes_[i]->edges_[0].empty() 
//                && nodes_[i]->edges_[1].empty();
////        for ( bool drxn : { 0, 1 } )
////        {
////            if ( !nodes_[i]->edges_[drxn].empty() || nodes_[i]->edges_[!drxn].empty() || doErase ) continue;
////            int readCounts[2] = { (int)nodes_[i]->reads_.size(), 0 };
////            if ( readCounts[0] > 1 ) continue;
////            NodeSet fwdSet;
////            for ( Node* prv : nodes_[i]->getNextNodes( !drxn ) ) prv->getDrxnNodes( fwdSet, drxn );
////            for ( Node* fwd : fwdSet )
////            {
////                if ( fwd == nodes_[i] ) continue;
////                readCounts[1] += fwd->reads_.size();
////            }
////            doErase = readCounts[1] > readCounts[0] * 5;
////        }
//        if ( doErase )
//        {
//            nodes_[i]->dismantleNode();
//            delete nodes_[i];
//            nodes_.erase( nodes_.begin() + i );
//        }
//        else i++;
//    }
//    
//    for ( Node* node : nodes_ )
//    {
//        for ( int i = 0; i < redundant.size(); i++ )
//        {
//            if ( !redundant[i] ) continue;
//            size_t it = node->seq_.find( ms.reads[i].seq );
//            if ( it == node->seq_.npos ) continue;
//            node->addRead( ms.reads[i].readId, node->ends_[0] + it, node->ends_[0] + it + ms.reads[i].seq.length(), true );
//        }
//        node->resetMarks();
//        node->setCoverage();
////        cout << string( node->ends_[0] + params.readLen, '-' ) << node->seq_ << endl;
//        ends_[0] = min( ends_[0], node->ends_[0] );
//        ends_[1] = max( ends_[1], node->ends_[1] );
//        validLimits_[0] = min( validLimits_[0], node->validLimits_[1] );
//        validLimits_[1] = max( validLimits_[1], node->validLimits_[2] );
//    }
}

void Seed::assemble()
{
    vector<SeedFork> forks[2];
    for ( Node* node : seed_.nodes ) for ( int d : { 0, 1 } ) if ( node->isContinue( d ) ) forks[d].push_back( SeedFork( node, NULL ) );
    for ( int d = 0; d < 2; d++ ) while ( SeedFork::extend( forks[d], bwt_, seed_, true, d ) );
    
    for ( bool good = true; good; )
    {
        for ( int d : { 0, 1 } ) while ( SeedFork::extend( forks[d], bwt_, seed_, false, d ) );
        Node::prune( bwt_, seed_ );
        if ( !restart( forks ) ) break;
    }
}

bool Seed::extend( NodeRoll &exts, bool drxn )
{
    for ( Node* node : exts.nodes ) node->extendNode( bwt_, seed_, drxn );
    Nodes fwdSet( exts.nodes, drxn, true );
    exts.clear();
    for ( Node* fwd : fwdSet.nodes ) if ( fwd->isContinue( drxn ) ) exts.add( fwd );
    return !exts.empty() && exts.nodes.size() < 20;
}

bool Seed::restart( vector<SeedFork> forks[2] )
{
    forks[0].clear();
    forks[1].clear();
    vector< vector<NodePath> > paths = NodePath::create( seed_ );
    for ( vector<NodePath>& path : paths )
    {
        for ( int d = 0; d < 1; d++ )
        {
            SeedFork fork( NULL, NULL );
            if ( fork.restart( bwt_, seed_, ( d ? path.back() : path[0] ).path_, d ) ) forks[d].push_back( fork );
        }
    }
    
    seed_.test();
    assert( !forks[0].empty() || !forks[1].empty() );
    return !forks[0].empty() || !forks[1].empty();
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

