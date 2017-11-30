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
//
//#include "Locus.h"
//#include <algorithm>
//
//void Locus::fill()
//{
//    plot();
//    setIslandVars( bwt_, paths_[1][0].path.back(), 1 );
//    int32_t limits[3] = { 0, 4000, 0 };
//    vector< pair<SeqNum, int32_t> > reads;
//    vector<ReadMark> marks;
//    unordered_set<SeqNum> Ids;
//    for ( Node* node : getAllNodes() )
//    {
////        limits[0] = min( limits[0], node->ends_[0] );
////        limits[1] = max( limits[1], node->ends_[1] );
//        limits[2] = max( limits[2], node->ends_[1] );
//        for ( bool drxn : { 0, 1 } )
//        {
//            for ( ReadMark &mark : node->marks_[drxn] )
//            {
//                if ( limits[0] < mark.estimate )
//                {
//                    if ( params.isReadPe( mark.readId ) )
//                    {
//                        marks.push_back( mark );
//                        reads.push_back( make_pair( mark.readId, ( drxn ? mark.estimate : mark.estimate - params.readLen ) ) );
//                    }
//                    Ids.insert( mark.readId );
//                }
//            }
//        }
//    }
//    
//    for ( ReadMark &mark : marks )
//    {
//        Node::seedIslandsSingle( *ivs_[1], mark, Ids, 1 );
//    }
//    
//    leapCleanup( *ivs_[1] );
//    leapExtend( *ivs_[1] );
//    
////    leapReview( *ivs_[1] );
////    ofstream align( "/home/glen/PythonProjects/BioJunkyard/data/Export/align4.fa" );
////    ofstream dump( "/home/glen/PythonProjects/BioJunkyard/data/Export/dump4" );
////    exportLocus( align, dump );
////    align.close();
////    dump.close();
//    
////    return;
//    
//    
//    sort( reads.begin(), reads.end(), []( pair<SeqNum, int32_t> &a, pair<SeqNum, int32_t> &b ){
//        return a.second < b.second;
//    } );
//    
//    ofstream fh( "/home/glen/PythonProjects/BioJunkyard/data/Export/align4.fa" );
//    
//    NodeList nodes;
//    int i = 0;
//    for ( Node* node : nodes_[1] )
//    {
//        node->id_ = to_string( i );
//        nodes.push_back( node );
//        i++;
//    }
//    NodeSet islandSet;
//    for ( Node* node : nodes_[4] )
//    {
//        if ( islandSet.find( node ) == islandSet.end() )
//        {
//            NodeSet thisSet = node->getConnectedNodes( true );
//            int hits = 0;
//            NodeSet hitSet, thisHitSet;
//            for ( Node* n : thisSet )
//            {
//                for ( auto &np : n->pairs_ )
//                {
//                    if ( np.first->drxn_ <= 2 )
//                    {
//                        hits += np.second;
//                        hitSet.insert( np.first );
//                        thisHitSet.insert( n );
//                    }
//                }
//            }
//            if ( hitSet.size() > 1 || hits >= 4 )
//            {
//                for ( Node* n : thisSet )
//                {
//                    if ( islandSet.find( n ) == islandSet.end() && thisHitSet.find( n ) != thisHitSet.end() )
//                    {
//                        nodes.push_back( n );
//                        for ( Node* nxt : n->getNextNodes( 0 ) )
//                        {
//                            if ( islandSet.find( nxt ) == islandSet.end() )
//                            {
//                                nodes.push_back( nxt );
//                            }
//                        }
//                    }
//                    islandSet.insert( n );
//                }
//            }
//        }
//    }
//    
////    NodeSet delSet;
////    for ( Node* node : nodes )
////    {
////        if ( node->seq_ == "GAGAAAAAAAAAAAAGAAAAGGGGAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAGGGAGAGAGAGAGAGAGAGAGAGAGAGAAGGGGGAGGTGGGAAGAAGCCGAGAGCAAAAGGGGCAAAAATTAAATAAAATGAAATGATAAGAGGAAATCGTAAAGAGAAAGAGAGAGAGAGCGCGCGCACAAATCTGAAGAATAAGACTAGAGAGACAGAGAGGGGCATGTGTGTATGCGTGTGTGCGAAAACAAAAAGGAGAGAATGAGAGAGTGAAAGAGAGAGAAAGAGACAG" )
////        {
////            node->seq_ = "GAGAAGAAGAAAAGAGAAAAGGGGAATAAGAGAAAAGAAAGAAGAAACAGAAAAAGGGAGAGAGAGAGAGAGAGAGAGAGAGAAGGGGGAGGTGGGAAGAAGCCGAGAGCAAAAGGGGCAAAAATTAAATAAAATGAAATGATAAGAGGAAATCGTAAAGAGAAAGAGAGAGAGAGCGCGCGCACAAATCTGAAGAATAAGACTAGAGAGACAGAGAGGGGCATGTGTGTATGCGTGTGTGCGAAAACAAAAAGGAGAGAATGAGAGAGTGAAAGAGAGAGAAAGAGACAG";
////            node->drxn_ = 1;
////            for ( Node* n : nodes )
////            {
////                if ( n->id_ == "49" )
////                {
////                    n->addEdge( node, 95, 1, true );
////                }
////                if ( n->id_ == "13" )
////                {
////                    node->addEdge( n, 85, 1, false );
////                }
////                if ( n->id_ == "14" )
////                {
////                    node->addEdge( n, 85, 1, false );
////                }
////                if ( n->id_ == "47" )
////                {
////                    delSet.insert( n );
////                }
////                if ( n->id_ == "48" )
////                {
////                    delSet.insert( n );
////                }
////            }
////            node->addRead( SeqNum(383001974), node->ends_[0] + 15, node->ends_[0] + 15 + params.readLen, false );
////            node->addRead( SeqNum(168695518), node->ends_[0] + 16, node->ends_[0] + 16 + params.readLen, false );
////            nodes_[1].push_back( node );
////        }
////    }
////    
////    this->deleteNodes( delSet, 1 );
////    
////    ofstream align( "/home/glen/PythonProjects/BioJunkyard/data/Export/align5.fa" );
////    ofstream dump( "/home/glen/PythonProjects/BioJunkyard/data/Export/dump5" );
////    exportLocus( align, dump );
////    align.close();
////    dump.close();
//    
////    sort( nodes.begin(), nodes.end(), []( Node* a, Node* b ){
////        return a->ends_[0] < b->ends_[0];
////    } );
////    
////    for ( Node* node : nodes )
////    {
////        fh << ">Node" << node->id_ << endl;
////        fh << string( node->ends_[0] - nodes[0]->ends_[0], '-' ) << node->seq_ << endl;
////    }
////    for ( pair<SeqNum, int32_t> &read : reads )
////    {
////        string seq = bwt_.getSequence( read.first );
////        fh << ">" << read.first << endl;
////        fh << string( read.second - ends[0], '-' ) << seq << endl;
////    }
//    
//    fh.close();
//    assert( false );
//}
//
//void Locus::fill2()
//{
//    NodeSet nodes[2];
//    int32_t ends[2];
//    for ( Node* node : nodes_[1] )
//    {
//        if ( node->id_ == "606" )
//        {
//            ends[0] = node->ends_[0] - 500;
//            node->getDrxnNodes( nodes[0], 0, node->ends_[1] - params.maxMpMax );
//        }
//        if ( node->id_ == "614" )
//        {
//            ends[1] = node->ends_[1] + 500;
//            node->getDrxnNodes( nodes[1], 1, node->ends_[0] + params.maxMpMax );
//        }
//    }
//    
//    vector< pair<SeqNum, int32_t> > reads;
//    
//    for ( bool drxn : {0,1} )
//    {
//        for ( Node* node : nodes[drxn] )
//        {
//            if ( node->coverage_ < params.cover * 2 )
//            {
//                for ( ReadMark &mark : node->marks_[drxn] )
//                {
//                    int32_t est = ( drxn ? mark.estimate - params.readLen : mark.estimate );
//                    if ( ends[0] <= est && est <= ends[1] )
//                    reads.push_back( make_pair( mark.readId, est ) );
//                }
//            }
//        }
//    }
//    
//    sort( reads.begin(), reads.end(), []( pair<SeqNum, int32_t> &a, pair<SeqNum, int32_t> &b ){
//        return a.second < b.second;
//    } );
//    
//    ofstream fh( "/home/glen/PythonProjects/BioJunkyard/data/Export/align9.fa" );
//    for ( pair<SeqNum, int32_t> &read : reads )
//    {
//        string seq = bwt_.getSequence( read.first );
//        fh << ">" << read.first << endl;
//        fh << string( max( 0, read.second - reads[0].second ), '-' ) << seq << endl;
//    }
//    
//    assert( false );
//}
//
//void Locus::fill3()
//{
//    ofstream fh( "/home/glen/PythonProjects/BioJunkyard/data/Export/align10.fa" );
//    
//    vector<string> seqs;
////    seqs.push_back( "CATAATAAAAGAAAACGGACTGAAAATGATTGGACACACTCTTCTAAACAAGATTGTAATACTTCATTAGGTTACAATTTTATGAGAGAAATAGAGAGAGTGGGGGGGGGGGGAGAAGAAGTGAGGAGAATAAAAGAAGACTAATAGCATAAGAGGAAGGAAAAGAGCTCGGAAGACAAAGAAG" );
////    seqs.push_back( "GAGAAATAGAGAGAGTGGGGGGGGGGGGAGAAGAAGTGAGGAGAATAAAAGAAGACTAATAGCATAAGAGGAAGGAAAAGAGCTCGGAAGACAAAGAAGTGGGTGTACAAAGGAATAGGAGGAGAAGGAGAGAGAGAAAGAGAGAGGGAAATAGAGAGAGAGGGGGGGGGGGGGGGGGAGGTTCAATGATAAAAAGGGAGGAGATGTTTTACATACAAAGTACCTGATGCAGGCGTTAATATCCATACGAATAACGGCTTT" );
////    seqs.push_back( "GGGGGAGAAGAAGTGAGGAGAATAAAAGAAGACTAATAGCATAAGAGGAAGGAAAAGAGCTCGGAAGACAAAGAAGAGGGTGTACAAAGGAATAGGAGGAGAAGGAGAGAGAGAGAGAGAAATAGAGAGAGAGAGAGAAGGGGGGGGGGGGGGGGAGGTTCAATGATAAAAAGGGA" );
//    seqs.push_back( "CTTCTACATTTACTTCTTCTCATTGTACGACTGGCTGCTGCTACTACTACTACTACTACTACTACTATACTACTTCTACATTTACTTCTTCTATTTCTACT" );
////    seqs.push_back( "ATACTACTTCTACATTTACTTCTTCTATTTCTACTTCTTCTCATTGTACGACTGGCTGCTAC" );
//    for ( string &seq : seqs )
//    {
//        MappedSeqs ms = bwt_.mapSeed( seq, 1, false );
//        int x = 0;
//        for ( ReadStruct &read : ms.reads )
//        {
//            fh << ">" << read.readId << endl;
//            fh << string( read.coords[0] - ms.reads[0].coords[0], '-' ) << read.seq << endl;
//        }
//    }
//    
//    assert( false );
//    
//}
//
//void Locus::fill4()
//{
//    int32_t ends[2] = { -50000, 50000 };
//    vector<ReadPairStruct> reads, notReads;
//    
//    for ( Node* node : getAllNodes() )
//    {
//        bool isReliable = node->coverage_ <= params.cover * 1.3;
//        for ( bool drxn : {0,1} )
//        {
//            for ( ReadMark &mark : node->marks_[drxn] )
//            {
//                bool doAdd = true;
//                if ( ends_[0] <= mark.estimate && mark.estimate <= ends[1] )
//                {
//                    ReadPairStruct read;
//                    read.readId = mark.readId;
//                    read.seq = bwt_.getSequence( mark.readId );
//                    read.coord = drxn ? mark.estimate : mark.estimate - params.readLen;
//                    read.setMin();
//                    for ( ReadPairStruct &read2 : reads )
//                    {
//                        doAdd = doAdd && read.seq != read2.seq;
//                    }
//                    if ( doAdd )
//                    {
//                        ( isReliable && params.isReadPe( mark.readId ) ? reads : notReads ).push_back( read );
//                    }
//                }
//            }
//        }
//    }
//    
//    int iMax = reads.size();
//    for ( ReadPairStruct &read : notReads )
//    {
//        bool doAdd = true;
//        int overlaps[2] = { 0, 0 };
//        int diff = ( params.isReadPe( read.readId ) ? 1000 : 2000 );
//        for ( int i( 0 ); i < reads.size(); i++ )
//        {
//            doAdd = doAdd && reads[i].seq  != read.seq;
//            if ( abs( read.coord - reads[i].coord ) <= diff && i < iMax )
//            {
//                overlaps[0] = max( overlaps[0], ReadPairStruct::overlap( reads[i], read, 12 ) );
//                overlaps[1] = max( overlaps[1], ReadPairStruct::overlap( read, reads[i], 12 ) );
//            }
//        }
//        
//        if ( doAdd && overlaps[0] > 0 && overlaps[1] > 0 )
//        {
//            reads.push_back( read );
//        }
//        else if ( doAdd && params.isReadPe( read.readId ) && ( overlaps[0] > 18 || overlaps[1] > 18 ) )
//        {
//            reads.push_back( read );
//        }
//    }
//    
//    sort( reads.begin(), reads.end(), []( ReadPairStruct &a, ReadPairStruct &b ){
//        return a.coord < b.coord;
//    });
//    
//    unordered_map< int, pair<int, int> > edges[2];
//    for ( int i( 0 ); i < reads.size(); i++ )
//    {
//        for ( int j ( i + 1 ); j < reads.size(); j++ )
//        {
//            int diff = ( params.isReadPe( reads[i].readId ) && params.isReadPe( reads[j].readId ) ? 1000 : 2000 );
//            
//            if ( abs( reads[i].coord - reads[j].coord ) <= diff )
//            {
//                int overlaps[2] = { ReadPairStruct::overlap( reads[j], reads[i], 12 ), 
//                                    ReadPairStruct::overlap( reads[i], reads[j], 12 ) };
//
//                if ( overlaps[0] > 0 )
//                {
//                    auto it = edges[0].find( i );
//                    if ( it == edges[0].end() || overlaps[0] > it->second.second )
//                    {
//                        edges[0][i] = make_pair( j, overlaps[0] );
//                    }
//
//                    it = edges[1].find( j );
//                    if ( it == edges[1].end() || overlaps[0] > it->second.second )
//                    {
//                        edges[1][j] = make_pair( i, overlaps[0] );
//                    }
//                    
//                    assert( edges[0].find( i ) != edges[0].end() );
//                    assert( edges[1].find( j ) != edges[1].end() );
//                    assert( edges[0][i].second >= overlaps[0] );
//                    assert( edges[1][j].second >= overlaps[0] );
//                }
//
//                if ( overlaps[1] > 0 )
//                {
//                    auto it = edges[1].find( i );
//                    if ( it == edges[1].end() || overlaps[1] > it->second.second )
//                    {
//                        edges[1][i] = make_pair( j, overlaps[1] );
//                    }
//
//                    it = edges[0].find( j );
//                    if ( it == edges[0].end() || overlaps[1] > it->second.second )
//                    {
//                        edges[0][j] = make_pair( i, overlaps[1] );
//                    }
//                    
//                    assert( edges[0].find( j ) != edges[0].end() );
//                    assert( edges[1].find( i ) != edges[1].end() );
//                    assert( edges[0][j].second >= overlaps[1] );
//                    assert( edges[1][i].second >= overlaps[1] );
//                }
//            }
//        }
//    }
//    
//    for ( int i : {0,1} )
//    {
//        for ( pair< const int, pair<int, int> > &edge : edges[i] )
//        {
//            assert( edges[!i].find( edge.second.first ) != edges[!i].end() );
//        }
//    }
//    
//    vector< unordered_set<int> > clusters;
//    for ( int i : {0,1} )
//    {
//        for ( pair< const int, pair<int, int> > &edge : edges[i] )
//        {
//            int clusterHit = -1;
//            for ( int j ( 0 ); j < clusters.size(); )
//            {
//                if ( clusters[j].find( edge.first ) != clusters[j].end() || clusters[j].find( edge.second.first ) != clusters[j].end() )
//                {
//                    clusters[j].insert( edge.first );
//                    clusters[j].insert( edge.second.first );
//                    if ( clusterHit == -1 )
//                    {
//                        clusterHit = j;
//                    }
//                    else
//                    {
//                        clusters[clusterHit].insert( clusters[j].begin(), clusters[j].end() );
//                        clusters.erase( clusters.begin() + j );
//                        continue;
//                    }
//                }
//                j++;
//            }
//            
//            if ( clusterHit == -1 )
//            {
//                unordered_set<int> cluster = { edge.first, edge.second.first };
//                clusters.push_back( cluster );
//            }
//        }
//    }
//    
//    ofstream align( "/home/glen/PythonProjects/BioJunkyard/data/Export/align11.fa" );
//    
//    for ( unordered_set<int> &cluster : clusters )
//    {
//        int32_t sum = 0;
//        for ( int i : cluster )
//        {
//            sum += reads[i].coord;
//            auto it = edges[1].find( i );
//            if ( it != edges[1].end() )
//            {
//                auto it2 = edges[0].find( it->second.first );
//                if ( it2 != edges[0].end() && it2->second.first != i && it->second.second < it2->second.second )
//                {
//                    edges[1].erase( it );
//                }
//            }
//        }
//        int32_t coord = sum / cluster.size();
//        
//        for ( bool j : {1,0} )
//        {
//            unordered_set<int> notLoop;
//            for ( int i : cluster )
//            {
//                if ( notLoop.find( i ) == notLoop.end() )
//                {
//                    unordered_set<int> loop = { i };
//                    auto it = edges[j].find( i );
//                    auto weakest = it;
//                    while ( it != edges[j].end() )
//                    {
//                        i = it->second.first;
//                        it = edges[j].find( i );
//                        if ( it != edges[j].end() )
//                        {
//                            if ( it->second.second < weakest->second.second )
//                            {
//                                weakest = it;
//                            }
//
//                            if ( loop.find( i ) != loop.end() )
//                            {
//                                edges[j].erase( it );
//                                break;
//                            }
//                        }
//                        loop.insert( i );
//                    }
//                    notLoop.insert( loop.begin(), loop.end() );
//                }
//            }
//        }
//        
//        unordered_map<int, int> used;
//        for ( int i : cluster )
//        {
//            if ( edges[0].find( i ) == edges[0].end() )
//            {
//                int offset = 0, j = i;
//                do
//                {
//                    i = j;
//                    align << ">" << coord << "_" << reads[i].readId << endl;
//                    align << string( offset, '-' ) << reads[i].seq << endl;
//                    used.insert( make_pair( i, offset ) );
//                    auto it = edges[1].find( i );
//                    if ( it != edges[1].end() )
//                    {
//                        offset += reads[i].seq.length() - it->second.second;
//                        j = it->second.first;
//                    }
//                } while ( i != j && used.find( j ) == used.end() );
//            }
//        }
//        
//        bool didAdd = true;
//        while ( didAdd )
//        {
//            didAdd = false;
//            for ( int i : cluster )
//            {
//                auto it = edges[0].find( i );
//                if ( it != edges[0].end() 
//                        && used.find( it->second.first ) != used.end() 
//                        && used.find( i ) == used.end() )
//                {
//                    didAdd = true;
//                    auto it2 = used.find( it->second.first );
//                    int offset = it2->second + reads[it->second.first].seq.length() - it->second.second;
//                    int j = i;
//                    do
//                    {
//                        i = j;
//                        align << ">" << coord << "_" << reads[i].readId << endl;
//                        align << string( offset, '-' ) << reads[i].seq << endl;
//                        used.insert( make_pair( i, offset ) );
//                        auto it3 = edges[1].find( i );
//                        if ( it3 != edges[1].end() )
//                        {
//                            offset += reads[i].seq.length() - it3->second.second;
//                            j = it3->second.first;
//                        }
//                    } while ( i != j && used.find( j ) == used.end() );
//                }
//            }
//        }
//        
//        for ( int i : cluster )
//        {
//            if ( used.find( i ) == used.end() )
//            {
//                align << ">" << coord << "_" << reads[i].readId << endl;
//                align << reads[i].seq << endl;
//            }
//        }
//    }
//}
//
