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

#include "locus.h"
#include <algorithm>
extern struct Parameters params;

Locus::Locus( Querier &bwt )
: bwt_( bwt ), finalise_( false ), calibrate_( false )
{
    duration_ = leapTime_ = revTime_ =  0;
    multiplicity_ = 1;
    debugLog_ = NULL;
    for ( int drxn : { 0, 1 } )
    {
        ivs_[drxn] = NULL;
        reliable_[drxn] = params.locusLimits[!drxn];
        stopCodes_[drxn] = 0;
        forkCount_[drxn] = 0;
        loopCount_[drxn] = 0;
        ends_[drxn] = 0;
        limits_[drxn] = 0;
        leapFar_[drxn] = false;
        harsh_[drxn] = false;
        validLimits_[drxn] = 0;
        forkLimits_[drxn] = 0;
        desperation_[drxn] = 0;
        completed_[drxn] = !params.drxns[drxn];
        finished_[drxn] = !params.drxns[drxn];
        paths_[drxn].push_back( Path() );
    }
}

Locus::Locus( Querier &bwt, Node* origin )
: Locus( bwt )
{
    origin->drxn_ = 2;
    nodes_[2].push_back( origin );
    originEnds_[0].push_back( origin );
    originEnds_[1].push_back( origin );
}

Locus::Locus( Querier &bwt, NodeList* nodes )
: Locus( bwt )
{
    NodeSet seedSet, dummy, delSet;
    Node::mergeAll( nodes, delSet );
    
    for ( int i : { 0, 1, 2 } )
    {
        for ( Node* node : nodes[i] )
        {
            if ( delSet.find( node ) != delSet.end() )
            {
                node->dismantleNode();
                delete node;
                continue;
            }
            node->setCoverage();
            seedSet.insert( node );
            node->drxn_ = i;
            nodes_[i].push_back( node );
            node->stop_[0] = 0;
            node->stop_[1] = 0;
            ends_[0] = min( ends_[0], node->ends_[0] );
            ends_[1] = max( ends_[1], node->ends_[1] );
            validLimits_[0] = min( validLimits_[0], node->validLimits_[1] );
            validLimits_[1] = max( validLimits_[1], node->validLimits_[2] );
            limits_[0] = min( limits_[0], node->ends_[0] - params.maxPeMax );
            limits_[1] = max( limits_[1], node->ends_[1] + params.maxPeMax );
        }
    }
    
    
    setOriginEnds();
    
    Node::resetPairing( seedSet );
    for ( Node* node : nodes_[2] ) node->setValid();
    Node::seedValidate( seedSet, dummy, validLimits_, ends_, false );
}

//Locus::Locus( Querier &bwt, ifstream &fh )
//: Locus( bwt )
//{
//    int32_t limits[2];
//    NodeSet seedSet, dummy, delSet;
//    NodeSet nxtSets[2];
//    for ( Node* node : Node::importNodes( fh ) )
//    {
//        if ( delSet.find( node ) == delSet.end() && node->reads_.empty() )
//        {
//            node->dismantleNode();
//            delete node;
//            continue;
//        }
//        nodes_[node->drxn_].push_back( node );
//        seedSet.insert( node );
//        
//        node->stop_[0] = node->edges_[0].empty();
//        node->stop_[1] = node->edges_[1].empty();
////        node->clearReads();
//    }
////    deleteNodes( delSet, 1 );
//    deleteNodes( delSet );
//    delSet.clear();
//    nodes_[2][0]->offsetForward( 1, false, true );
//    nodes_[2][0]->offsetForward( 0, false, true );
//    
//    NodeList nodes = getAllNodes();
////    Node::graphCover( "/media/glen/ssd/HtCoverage.csv", nodes );
////    Node::graphPairs( "/media/glen/ssd/HtPairing.csv", nodes );
//    Node::remapGenes( bwt, nodes );
//    
////    setOriginEnds();
////    MapNode* mn1 = new MapNode();
////    mn1->seq = "ACAACAACAACAACAAACAACTGCACACGGATATGATTACAGCCTATACAGTTTCATTCAGTTTGATGATGCTACTGTTCTGACAGTCATTTTTATCAACACAAGTGAGGCAAATTATAAACAGTCACGTATTTTGGTGCAATATACCATGCTCCGTGATGAGATTTCATGACTTTGCCTAACAAAATTATATAACTTTATTACTCCCTTACAGTTCATGCACTAAGTGAACGGAGAGAAAGAGGAAATGGCAGGGGACTAGGTCGCTTCGGAGGAAGGCCAGGATCTGATAGACCCCAAATGATGGTCGGACCTAGGCAAGATGGACCACCAATGGGCGGAAGGAGGTTTGATGGCCTTGGACAAGGTGACCCACTGATCGGTGGACGTGGACAAAATGGCTGGCCGATGGGCGGTAGGAGGTTTGATGGCCCTGGACAAGATGACCCAGTGATCGGTGGACGTGGACCAAATGGCGGACCAGTGGGCGGTAGGAGGTTTGATGGACCTGGATTTGGTGGCTTCAGACCCGAAGGTATAGGGAGACCTTTCTTCGGTCACGGAGGAAGGCATGCTGACGGAGAAGGAGAAATGGAGGCTGCTCAACCAATCGGTGATGGTCAAGGATGGCCCGATCGTTTCGATGGTCCTGGACGATTTTCCGGACGTCCTCACCCAGGCCGTGGTGGCTATTATGGACATCACCATAGTCCTCCCCATGACCTGCAAACCGACGAACAACCGTTCGGTCAGTACAACGACAGCAGCAGCGAGGAGGATGGCCGACCTCACCCTCATCCTCACCCTCTCCATCACCCCCATCCTCCCCATCACCCTCCACCACTACATCACCCTCCCCATCACCCTCCCCATCACCCTCCACCACCCCATCACCCTCCCCATCACCCTCCACCTCCCCCTCACTCTTACCACCACCATCTCCACCACCATGACCATCATAACAAGACAGACGACCACCGTAATCATAATCACACTGGAGGTCACCTCCACCATCATCATAACCAGACAGAAGAGTGGGACCAGGACAGGCCAGATATGAGGCCATTCCGGTTCAATTCTTTCGATAGCCAGGAGGATGGCCGACCTTACCCTCACCCTTACCCTCCCCATCACCCACCCCATTACCC";
////    bwt.mapSequence( mn1->seq, mn1->ids, mn1->coords );
////    Node* node1 = new Node( mn1, 0, mn1->ids.size()-1, 1 );
////    MapNode* mn2 = new MapNode();
////    mn2->seq = "CCCCATCACCCTCCACCACCCCACCCCCCTCCCCATCACCCTCCACCTCCCCATCACCCTCCCCATCACCCTCACTCTCACCACCACCACCATGACCATCATAACAAGACAGACGACCACCGTAATCATAATCACACGGCAGGCCACATCCACCATCATCATAACCAGACAGAAGAGTGGGACCAGGACAGACCAGAGATGAGGCCATTCTGGTTCAATCCTTTCGGTCGCAAGCCTTTCGGACGACGTCCATTTGGCAGACGCTACCACACCGAAGAGGGATTTCCCAGGCGCGATGGACACCGTCATCCTCATGGCATGCATGGCAACCGAGGACGTTGGGATGAGAATGAAAATGAGGAGGAAGAACATCTCCCGACTGAAAGAATGACAACCTCTGCAATGCCTAAAGTGGTCGAAATCGATATCACCGAAATAGACAACAACGTCATCGACAAAGTGTAGACTTTGTCAAGAAATTTTTTGCTTTTGAGTAAACGAACTATTGTCAGGGGCGTCGGTCCGTTTTTCGGAGGGGGGGGGGGGGATAAAGAGAGCCTAATAATACCCCTCATGATTTGTTAACGCCTCTTAATTATGCTAATTACCAGCTTCTACCCGGAGCGAGTGGAGCGGAGCGAGCGGCAAAATTGAAGTCTGATCGTGTAATACATTTTATATTTCATTGTTTGGAAGGTTTCATTTCCCGAGTTTTGCATCATTTTCCCTTTTTTTATTGTACACCTCCGGGTTCATTGGTGGGGCAAAACGATATGTTTGCCTCCAAATTTTTTTTTGTGGGGGCAACTGCCCCCTACCACCCCACAATCGACGCCCCTGACTATTATCAATGTGCGGCTGGCAGATGCTAGGGCGTTTAACGATTTGCACCTTTCTGCTAAACGTAACGAGTACCATTTCGTTTAGCAAAAGGTCGTTGGCACTAATATCTTAAAAATGTATGATTATTTTTGATACAAAACTATATGTACTAAAATTGGGTAGTTACCATTTTCACTTTTTTATTAGTCGATAAAGAGAGAAATCCCTGTCATGTTTCGAAAATGGAATACAGAGCGATCAAAACCTTAGAAGCGAGTTACGAGAAAATTGATTTATGGCATGCTCGCGTGCCTAAAGTCATTCTGTCTTTTCCATTAACAGTATTCCAGAAATTCAATGCCATGTTCATAACATTCGAAAGTCTACCATTTTATTACAAACTCAAAGAGACTAGAAGGTGAAACCAACAAATTGAAATTATAATTTCCATTTGCACGGTATTGTATATGTGAAATAAAGATATATCGCAATTCATTACAAAATGATAATTATGTAGCTTTTATTTTATTATGTAGGCCTTGATATTTTCAATCGTGACTGAACTGAGTATACTCTTTTCTTTTTTTATTTTTTTTTTCATCCTTCTTAGTCATAAATACTAGGATCTATTAAATAATATATTATTAATAGATCCTAGTATTTATGACGAAGAAGAAGGATGAGAAAAATAAGAAAAAAAA";
////    bwt.mapSequence( mn2->seq, mn2->ids, mn2->coords );
////    Node* node2 = new Node( mn2, 0, mn2->ids.size()-1, 1 );
////    
////    node2->offset( 1200 );
////    cout << ">node1" << endl;
////    cout << node1->seq_ << endl;
////    cout << ">node2" << endl;
////    cout << node2->seq_ << endl;
////    node1->addEdge( node2, -80, 1, true, true );
//////    ExtVars ev( nodes_[1], nodes_[4], validLimits_, bwt_ );
//////    IslandVars iv( ev, 1 );
//////    iv.origin.push_back( nodes_[2][0] );
//////    node1->extendCount_ = 20;
//////    node1->extendIsland( iv, 0 );
////    node1->drxn_ = 1;
////    node2->drxn_ = 1;
////    nodes_[1].push_back( node1 );
////    nodes_[1].push_back( node2 );
//////    node2->extendCount_ = 20;
//////    node2->extendNode( ev, 1 );
////    for ( Node* prv : nxtSets[0] )
////    {
////        int ol = mapSeqOverlap( prv->seq_, node1->seq_, 5 );
////        prv->addEdge( node1, ol, 1, true, false );
////    }
////    for ( Node* nxt : nxtSets[1] )
////    {
////        int ol = mapSeqOverlap( node2->seq_, nxt->seq_, 5 );
////        node2->addEdge( nxt, ol, 1, true, false );
////    }
//    
////    for ( ReadMark &mark : node1->marks_[0] )
////    {
////        if ( !params.isReadPe( mark.id ) ) continue;
////        if ( node1->reads_.find( mark.id ) != node1->reads_.end() ) continue;
////        if ( node2->reads_.find( mark.id ) != node2->reads_.end() ) continue;
////        cout << ">" << mark.id << endl;
////        cout << string( mark.estimate - 101, '-' ) << bwt.getSequence( mark.id ) << endl;
////    }
////    for ( ReadMark &mark : node2->marks_[1] )
////    {
////        if ( !params.isReadPe( mark.id ) ) continue;
////        if ( node1->reads_.find( mark.id ) != node1->reads_.end() ) continue;
////        if ( node2->reads_.find( mark.id ) != node2->reads_.end() ) continue;
////        cout << ">" << mark.id << endl;
////        cout << string( mark.estimate, '-' ) << bwt.getSequence( mark.id ) << endl;
////    }
//    
////    vector<int> offsets;
////    for ( ReadMark &mark : node1->marks_[0] )
////    {
////        if ( !params.isReadPe( mark.id ) ) continue;
////        auto it = node2->reads_.find( mark.id );
////        if ( it != node2->reads_.end() )
////        {
////            int32_t coord = it->second[1];
////            offsets.push_back( mark.estimate - it->second[1] );
////        }
////    }
////    
////    assert( false );
////    
////    MapNode* mn = new MapNode();
////    mn->seq = "CCCTCCCCATCACCCTCCACCTCCCCATCACCCTCCCCATCACCCTCACTCTCACCACCACCACCATGACCATCATAACAAGACAGACGACCACCGTAATCATAAT";
//////    mn->seq = "CTTACCACCACCATCTCCACCACCATGACCATCATAACAAGACAGACGACCACCGTAATCATAATCACACTGGAGGTCACCTCCACCATCATCATAACCAGACAGAAGAGTGGGACCAGGACAGGCCAGATATGAGGCCATTCCGGTTCAATTCTTTCGATAGCCAGGAGGATGGCCGACCTTACCCTCACCCTTACCCTCCCCATCACCCACCCCATTACCC";
////    bwt.mapSequence( mn->seq, mn->ids, mn->coords );
////    mn->recoil();
////    Node* newNode = new Node( mn, 0, mn->ids.size()-1, 1 );
////    ExtVars ev( nodes_[1], nodes_[4], validLimits_, bwt_ );
////    IslandVars iv( ev, 1 );
////    iv.origin.push_back( nodes_[2][0] );
////    newNode->extendCount_ = 20;
////    newNode->extendIsland( iv, 0 );
////    nodes_[1].push_back( newNode );
////    for ( Node* n : iv.merged[0] ) n->drxn_ = 1;
////    newNode->extendNode( ev, 1 );
////    for ( bool drxn : { 0, 1 } )
////    {
////        for ( ReadMark &mark : newNode->marks_[drxn] )
////        {
////            if ( !params.isReadPe( mark.id ) ) continue;
////            cout << ">" << mark.id << endl;
////            cout << string( ( mark.estimate - newNode->ends_[0] )+600, '-' ) << bwt.getSequence( mark.id ) << endl;
////        }
////    }
////    assert( false );
////    for ( Node* nxt : newNode->getNextNodes( 1 ) )
////    {
////        nxt->extendCount_ = 20;
////        nxt->extendNode( ev, 1 );
////        if ( nxt->stop_[1] == 1 )
////        {
////            nxt->dismantleNode();
////            delSet.insert( nxt );
////            deleteNodes( delSet );
////            delSet.clear();
////        }
////    }
////    
////    assert( false );
//    
////    NodeList allNodes = getAllNodes();
//////    limits[0] -= 300;
//////    limits[1] += 300;
////    int j = 0;
////    for ( Node* n : allNodes )
////    {
////        for ( int i : { 0, 1 } )
////        {
////            for ( ReadMark &mark : n->marks_[i] )
////            {
////                if ( mark.estimate < limits[0] || limits[1] < mark.estimate ) continue;
////                if ( !params.isReadPe( mark.id ) ) continue;
////                bool doAdd = true;
////                for ( Node* n2 : allNodes )
////                {
////                    if ( n2->reads_.find( mark.id ) != n2->reads_.end() )
////                    {
////                        doAdd = false;
////                        break;
////                    }
////                }
////                if ( !doAdd ) continue;
////                string seq = bwt_.getSequence( mark.id );
////                cout << ">" << j << endl;
////                cout << string( mark.estimate - limits[0], '-' ) << seq << endl;
////                j++;
////            }
////        }
////    }
////    
////    assert( false );
//////    ExtVars ev( nodes_[0], nodes_[3], validLimits_, bwt_ );
//////    IslandVars iv( ev, 0 );
//////    iv.origin.push_back( nodes_[2][0] );
////    for ( Node* node : nodes_[1] )
////    {
////        if ( node->id_ != "589" ) continue;
////        int ol = mapSeqOverlap( node->seq_, newNode->seq_, 5 );
////        node->addEdge( newNode, ol, 1 );
////    }
////    for ( Node* node : nodes_[1] )
////    {
////        if ( node->id_ != "297" ) continue;
////        int ol = mapSeqOverlap( newNode->seq_, node->seq_, 5 );
////        newNode->addEdge( node, ol, 1 );
////    }
////    nodes_[1].push_back( newNode );
//    
////    newNode->drxn_ = 3;
////    newNode->extendCount_ = 20;
////    newNode->extendIsland( iv, 1 );
////    for ( Node* nxt : newNode->getNextNodes( 1 ) )
////    {
////        nxt->extendCount_ = 20;
////        nxt->extendIsland( iv, 1 );
////    }
////    newNode->drxn_ = 0;
////    newNode->extendCount_ = 20;
////    newNode->extendNode( ev, 0 );
////    for ( Node* nxt : newNode->getNextNodes( 0 ) )
////    {
////        if ( !nxt->isContinue( 0 ) ) continue;;
////        nxt->extendCount_ = 20;
////        nxt->extendNode( ev, 0 );
////    }
////    for ( Node* node : nodes_[0] )
////    {
////        if ( node->seq_.find( "TGGTCATACCCCCGCCCCCAGTTGAAAAGAAATTATTGCAT" ) == node->seq_.npos ) continue;
////        NodeSet nxtSet = node->getNextNodes( 1 );
//////        node->clearEdges( 1 );
////        int ol = mapSeqOverlap( node->seq_, newNode->seq_, 15 );
//////        if ( !ol ) ol = -1;
////        node->addEdge( newNode, -22, 0, true, true );
//////        for ( Node* nxt : nxtSet )
//////        {
//////            ol = mapSeqOverlap( newNode->seq_, nxt->seq_, 15 );
//////            newNode->addEdge( nxt, ol, 1 );
//////        }
//////        break;
////    }
//    
////    Node::seedValidate( seedSet, dummy, validLimits_, ends_, false );
////    int count = 0;
////    for ( Node* node : getAllNodes() )
////    {
////        node->mapMates( bwt_, count );
////    }
////    assert( false );
////    this->plot();
////    Node::mergeAll( nodes_, delSet );
////    deleteNodes( delSet );
////    delSet.clear();
////    ofstream dump( "/media/glen/ssd/dump15" );
////    ofstream align( "/media/glen/ssd/align2.fa" );
//    this->exportLocus( align, dump );
//    assert( false );
//}

Locus::~Locus()
{
    for ( int i( 0 ); i < 5; i++ )
    {
        for ( Node* node : nodes_[i] )
        {
            node->dismantleNode();
            delete node;
        }
        nodes_[i].clear();
    }
    for ( int i ( 0 ); i < 2; i++ )
    {
        if ( ivs_[i] )
        {
            delete ivs_[i];
            ivs_[i] = NULL;
        }
    }
}

//void Locus::locusTest()
//{
//    for ( int i = 0; i < 2; i++ )
//    {
//        for ( Node* node : forkNodes_[i] )
//        {
//            assert( node->drxn_ == 2 || find( nodes_[i].begin(), nodes_[i].end(), node ) != nodes_[i].end() );
//        }
//        
//        for ( Path &path : paths_[i] )
//        {
//            for ( Span &span : path.spans )
//            {
//                assert( span.node->drxn_ == 2 || find( nodes_[i].begin(), nodes_[i].end(), span.node ) != nodes_[i].end() );
//            }
//        }
//    }
//    for ( int i = 0; i < 5; i++ )
//    {
//        for ( Node* node : nodes_[i] )
//        {
//            int count = 0;
//            for ( Node* n2 : nodes_[i] ) count += node == n2;
//            assert( count == 1 );
//            
////            assert( !node->unreliable_ );
//            assert( node->drxn_ == i );
//            float cover = node->coverage_;
//            node->setCoverage();
//            assert( cover * 2 > node->coverage_ || node->edges_[1].empty() || node->edges_[0].empty() );
//            bool ended[2] = { false, false };
//            for ( auto &read : node->reads_ )
//            {
//                if ( read.second[0] == node->ends_[0] ) ended[0] = true;
//                if ( read.second[1] == node->ends_[1] ) ended[1] = true;
//                assert( node->ends_[0] <= read.second[0] && read.second[1] <= node->ends_[1] );
//            }
//            assert( ended[0] && ended[1] );
//            if ( i == 0 ) assert( !node->edges_[1].empty() );
//            if ( i == 1 ) assert( !node->edges_[0].empty() );
//            for ( bool drxn : { 0, 1 } )
//            {
//                int j = 0;
//                for ( Edge &e : node->edges_[drxn] )
//                {
////                    assert( drxn != node->drxn_ || node->drxn_ != 1 || e.node->drxn_ == 1 );
////                    assert( drxn != node->drxn_ || node->drxn_ != 0 || e.node->drxn_ == 0 );
//                    if ( i == 1 || i == 0 )
//                    {
//                        assert( e.node->drxn_ != i || find( nodes_[i].begin(), nodes_[i].end(), e.node ) != nodes_[i].end() );
//                    }
//                    bool found = false;
//                    int k = 0;
//                    for ( Edge &e2 : e.node->edges_[!drxn] )
//                    {
//                        if ( e2.node == node ) found = true;
//                        assert( k++ < 30 );
//                    }
//                    assert( j < 30 );
//                    assert( found );
//                    j++;
//                }
//            }
//            for ( bool drxn : { 0, 1 } )
//            {
//                for ( ReadMark &mark : node->marks_[drxn] )
//                {
//                    if ( drxn ) assert( mark.estimate < mark.mark );
//                    else ( assert( mark.mark < mark.estimate ) );
//                }
//                int seqLen = node->seq_.length();
//                for ( Edge &e : node->edges_[drxn] )
//                {
//                    assert( seqLen >= e.overlap );
//                }
//            }
//        }
//    }
//}

NodeListList Locus::exportNodes()
{
    vector< vector<Node*> > nodes( 3 );
    int nodeId( 0 );
    for ( Node* node : nodes_[2] )
    {
        node->id_ = to_string( nodeId );
        nodeId++;
        nodes[2].push_back( node );
    }
    for ( int drxn : {0, 1} )
    {
        for ( Node* node : nodes_[drxn] )
        {
            node->id_ = to_string( nodeId );
            nodeId++;
            nodes[drxn].push_back( node );
        }
    }
    return nodes;
}

void Locus::rebootLocus()
{
    setOriginEnds();
    vector<Node*> nodes = getAllNodes();
    bool first = true;
    for ( Node* node : nodes )
    {
        ends_[0] = first ? node->ends_[0] : min( ends_[0], node->ends_[0] );
        ends_[1] = first ? node->ends_[1] : max( ends_[1], node->ends_[0] );
        if ( node->drxn_ == 2 )
        {
            node->validate( validLimits_ );
        }
        node->resetMarks();
        first = false;
    }
    setExtend( 0 );
    setExtend( 1 );
}

void Locus::locusTest()
{
    for ( bool drxn : { 0, 1, 2 } )
    {
        for ( Node* node : nodes_[drxn] )
        {
            node->readTest();
//            assert( !node->farPairNodes_[0] || node->pairs_.find( node->farPairNodes_[0] ) != node->pairs_.end() );
            
//            assert( !node->reads_.empty() );
//            if ( node->clones_ )
//            {
//                NodeSet cloneSet = node->getCloneSet( true );
//                NodeListList cloneLists;
//                for ( Node* clone : cloneSet )
//                {
//                    cloneLists.push_back( *clone->clones_ );
//                }
//                for ( Node* clone : *node->clones_ )
//                {
//                    assert( clone->reads_.size() == node->reads_.size() );
//                }
//            }
        }
    }
    
    
    
//    for ( bool drxn : { 0, 1 } )
//    {
//        for ( Node* node : nodes_[drxn] )
//        {
////            node->farTest( drxn );
////            assert( !node->edges_[!drxn].empty() );
//            node->offsetTest( drxn );
//        }
//    }
    
//    for ( int i ( 0 ); i < 5; i++ )
//    {
//        for ( Node* node : nodes_[i] )
//        {
//            for ( bool drxn : {0,1} )
//            {
//                for ( Edge &f : node->edges_[drxn] )
//                {
//                    bool didFind = false;
//                    int j = 0;
//                    for ( Edge &r : f.node->edges_[!drxn] )
//                    {
//                        if ( r.node == node )
//                        {
//                            didFind = true;
//                            break;
//                        }
//                        assert( j < 20 );
//                        j++;
//                    }
//                    assert( didFind );
//                }
//            }
//        }
//    }
}

//void Locus::locusTest( bool drxn )
//{
//    for ( Node* node : nodes_[drxn] )
//    {
//        for ( Node* clone : node->getCloneSet() )
//        {
//            assert( clone->reads_.size() == node->reads_.size() );
//        }
//    }
//}

void Locus::addNode( Node* node, int drxn )
{
    node->drxn_ = drxn;
    nodes_[drxn].push_back( node );
}

void Locus::deleteNodes( NodeSet &delSet )
{
    for ( Node* del : delSet )
    {
        for ( int drxn : { 0, 1, 2, 3, 4 } )
        {
            nodes_[del->drxn_].erase( remove( nodes_[del->drxn_].begin(), nodes_[del->drxn_].end(), del ), nodes_[del->drxn_].end() );
            endNodes_[drxn].erase( remove( endNodes_[drxn].begin(), endNodes_[drxn].end(), del ), endNodes_[drxn].end() );
            sideNodes_[drxn].erase( remove( sideNodes_[drxn].begin(), sideNodes_[drxn].end(), del ), sideNodes_[drxn].end() );
            forkNodes_[drxn].erase( remove( forkNodes_[drxn].begin(), forkNodes_[drxn].end(), del ), forkNodes_[drxn].end() );
            toExtend_[drxn].erase( remove( toExtend_[drxn].begin(), toExtend_[drxn].end(), del ), toExtend_[drxn].end() );
        }
        delete del;
    }
}

void Locus::deleteNodes( NodeSet &delSet, bool drxn )
{
    NodeSet baseDelNodes( delSet.begin(), delSet.end() );
    for ( Node* node : baseDelNodes )
    {
        node->dismantleNode( delSet, drxn );
    }
    
    for ( Node* del : delSet )
    {
        assert( del->drxn_ != 2 );
        assert( del->drxn_ <= 2 );
        nodes_[del->drxn_].erase( remove( nodes_[del->drxn_].begin(), nodes_[del->drxn_].end(), del ), nodes_[del->drxn_].end() );
        endNodes_[drxn].erase( remove( endNodes_[drxn].begin(), endNodes_[drxn].end(), del ), endNodes_[drxn].end() );
        sideNodes_[drxn].erase( remove( sideNodes_[drxn].begin(), sideNodes_[drxn].end(), del ), sideNodes_[drxn].end() );
        forkNodes_[drxn].erase( remove( forkNodes_[drxn].begin(), forkNodes_[drxn].end(), del ), forkNodes_[drxn].end() );
        toExtend_[drxn].erase( remove( toExtend_[drxn].begin(), toExtend_[drxn].end(), del ), toExtend_[drxn].end() );
        delete del;
    }
}

void Locus::deleteNodes( NodeSet &delSet, NodeSet &goodSet, bool drxn )
{
    for ( auto it = delSet.begin(); it != delSet.end(); )
    {
        if ( goodSet.find( *it ) != goodSet.end() )
        {
            it = delSet.erase( it );
            continue;
        }
        it++;
    }
    if ( !delSet.empty() )
    {
        deleteNodes( delSet, drxn );
    }
}

NodeList Locus::getAllNodes()
{
    NodeList allNodes = nodes_[2];
    allNodes.insert( allNodes.end(), nodes_[0].begin(), nodes_[0].end() );
    allNodes.insert( allNodes.end(), nodes_[1].begin(), nodes_[1].end() );
    return allNodes;
}

NodeListList Locus::getEndLists( NodeSet &fwdSet, bool drxn )
{
    NodeListList endListList;
    for ( Node* fwd : fwdSet )
    {
        if ( fwd->isContinue( drxn ) )
        {
            NodeList endList = { fwd };
            NodeSet endSet;
            fwd->getDrxnNodesInSet( endSet, fwdSet, !drxn );
            endList.insert( endList.end(), endSet.begin(), endSet.end() );
            endListList.push_back( endList );
        }
    }
    
    return endListList;
}

uint Locus::getLen()
{
    return ends_[1] - ends_[0];
}

ScoreMap Locus::getScoreMap( Node* bgn, NodeList &tNodes, bool drxn )
{
    int32_t limits[2];
    Node::getLimits( limits, tNodes );
    ScoreMap scores;
    NodeSet qSet = { bgn };
    bgn->getDrxnNodes( qSet, drxn );
    for ( Node* q : qSet )
    {
        Score score = q->getPairScore( tNodes, limits, drxn );
        q->limitScore( score, true );
        scores.insert( make_pair( q, score ) );
    }
    return scores;
}

NodeSet Locus::getValidSet()
{
    NodeSet vSet;
    for ( NodeList &nodes : nodes_ )
    {
        for ( Node* node : nodes )
        {
            if ( node->isValidated() )
            {
                vSet.insert( node );
            }
        }
    }
    return vSet;
    
}

NodeList Locus::getValidNodes()
{
    NodeList vNodes;
    for ( NodeList &nodes : nodes_ )
    {
        for ( Node* node : nodes )
        {
            if ( node->isValidated() )
            {
                vNodes.push_back( node );
            }
        }
    }
    return vNodes;
}

int Locus::getWeakestEdge( Node* begin, Node* end, NodeSet &fwdSet, bool drxn )
{
    Node* thisNode = end;
    int weakest = params.readLen;
    while ( thisNode != begin && thisNode )
    {
        Node* nextNode = NULL;
        int thisBest = 0;
        for ( Edge &edge : thisNode->edges_[!drxn] )
        {
            if ( edge.overlap > thisBest && ( edge.node == begin || fwdSet.find( edge.node ) != fwdSet.end() ) )
            {
                thisBest = edge.overlap;
                nextNode = edge.node;
            }
        }
        thisNode = nextNode;
        weakest = min( weakest, thisBest );
    }
    return weakest;
}

void Locus::setForkLimits()
{
    forkLimits_[0] = forkLimits_[1] = 0;
    for ( Node* node : forkNodes_[0] ) forkLimits_[0] = min( forkLimits_[0], node->ends_[0] );
    for ( Node* node : forkNodes_[1] ) forkLimits_[1] = max( forkLimits_[1], node->ends_[1] );
}

void Locus::setHeader( string &header )
{
    header_ = header;
}

void Locus::setOriginEnds()
{
    for ( int drxn : { 0, 1 } )
    {
        originEnds_[drxn].clear();
        for ( Node* node : nodes_[2] )
        {
            bool originEnd = node->edges_[drxn].empty();
            for ( Edge edge : node->edges_[drxn] )
            {
                if ( find( nodes_[2].begin(), nodes_[2].end(), edge.node ) == nodes_[2].end() )
                {
                    originEnd = true;
                }
            }
            if( originEnd  )
            {
                originEnds_[drxn].push_back( node );
                node->extendCount_ = 8;
            }
            reliable_[!drxn] = drxn ? min( reliable_[0], node->ends_[1] ) : max( reliable_[1], node->ends_[0] );
        }
    }
}

void Locus::summarise( bool drxn )
{
    cout << endl << "NEXT EXTENSION" << endl;
    cout << "Number of nodes to extend:" << this->toExtend_[drxn].size() << endl;
    bool first = true;
    int leftest, rightest;
    for ( Node* node : toExtend_[drxn] )
    {
        if ( first )
        {
            leftest = node->ends_[drxn];
            rightest = node->ends_[drxn];
            first = false;
        }
        else
        {
            leftest = min( leftest, node->ends_[drxn] );
            rightest = max( rightest, node->ends_[drxn] );
        }
    }
    if ( !toExtend_[drxn].empty() )
    {
        cout << "To extend nodes range from " << leftest << " to " << rightest << endl << endl;
    }
}
