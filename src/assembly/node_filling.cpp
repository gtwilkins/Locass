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

void Node::checkClones()
{
    NodeSet eraseSet;
    for ( Node* clone : *clones_ )
    {
        if ( clone->seq_ != seq_ )
        {
            eraseSet.insert( clone );
        }
    }
    
    for ( Node* node : eraseSet )
    {
        node->clones_->erase( remove( node->clones_->begin(), node->clones_->end(), this ), node->clones_->end() );
        if ( node->clones_->empty() )
        {
            delete node->clones_;
            node->clones_ = NULL;
        }
        clones_->erase( remove( clones_->begin(), clones_->end(), node ), clones_->end() );
        if ( clones_->empty() )
        {
            delete clones_;
            clones_ = NULL;
        }
    }
}

void Node::clearReads()
{
    reads_.clear();
    marks_[0].clear();
    marks_[1].clear();
}

void Node::fillReads( Querier &bwt, NodeSet &delSet )
{
    int32_t offsets[2] = { 0, 0 };
    size_t it = seq_.find_first_not_of( 'N' );
    string seq = seq_;
    if ( it )
    {
        offsets[0] = it;
        seq = seq.substr( it );
    }
    it = seq.find_first_of( 'N' );
    if ( it != seq.npos )
    {
        offsets[1] = seq.length() - it;
        seq = seq.substr( 0, it );
    }
    clearReads();
    vector<ReadId> ids;
    vector<int32_t> coords[2];
    bwt.mapSequence( seq, ids, coords );
    int32_t ends[2] = { (int32_t)seq_.length(), 0 };
    int32_t cEnds[2] = { ends_[1], ends_[0] };
    for ( int i = 0; i < ids.size(); i++ )
    {
        coords[0][i] += offsets[0];
        coords[1][i] += offsets[0];
        ends[0] = min( ends[0], coords[0][i] );
        ends[1] = max( ends[1], coords[1][i] );
        coords[0][i] += ends_[0];
        coords[1][i] += ends_[0];
        cEnds[0] = min( cEnds[0], coords[0][i] );
        cEnds[1] = max( cEnds[1], coords[1][i] );
        addRead( ids[i], coords[0][i], coords[1][i], false );
    }
    int diffs[2] = { ends[0], (int32_t)seq_.length() - ends[1] };
    
    if ( reads_.empty() )
    {
        NodeSet blockSet = { this };
        for ( Edge &r : edges_[1] )
        {
            NodeSet reachSet = r.node->getDrxnNodesNotInSet( blockSet, 0 );
            vector<int> ols;
            NodeList eNodes;
            for ( Edge &l : edges_[0] )
            {
                if ( reachSet.find( l.node ) == reachSet.end() )
                {
                    ols.push_back( l.overlap - ( seq_.length() - r.overlap ) );
                    eNodes.push_back( l.node );
                }
            }
            for ( int i = 0; i < eNodes.size(); i++ )
            {
                r.node->addEdge( eNodes[i], ols[i], 0, false, ols[i] <= 0 );
            }
        }
        dismantleNode();
        delSet.insert( this );
        return;
    }
    
    if ( diffs[0] ) recoil( diffs[0], 0 );
    if ( diffs[1] ) recoil( diffs[1], 1 );
}

void Node::graphCover( string filename, NodeList &nodes )
{
    vector<int32_t> coords[2];
    coords[0].push_back( -6452 );
    coords[1].push_back( -4731 );
    coords[0].push_back( -52 );
    coords[1].push_back( 1269 );
    coords[0].push_back( 9202 );
    coords[1].push_back( 10326 );
    
    ofstream csv( filename );
    vector<int> allCovers[3];
    int maxLen = 0;
    for ( int i = 0; i < 3; i++ )
    {
        int len = coords[1][i]-coords[0][i];
        maxLen = max( maxLen, len );
        int covers[len]{0};
        for ( Node* n : nodes )
        {
            if ( n->ends_[1] < coords[0][i] || coords[1][i] < n->ends_[0] ) continue;
            for ( auto &read : n->reads_ )
            {
                int j = max( read.second[0], coords[0][i] );
                int k = min( read.second[1], coords[1][i] );
                j -= coords[0][i];
                k -= coords[0][i];
                while ( j < k ) covers[j++]++;
            }
        }
        for ( int j = 0; j < len; j++ )
        {
            allCovers[i].push_back( covers[j] );
        }
    }
    
    for ( int i = 0; i < maxLen; i++ )
    {
        for ( int j = 0; j < 3; j++ )
        {
            if ( i < allCovers[j].size() ) csv << to_string( allCovers[j][i] );
            if ( j < 2 ) csv << ",";
            else csv << endl;
        }
    }
    
    csv.close();
    assert( false );
}

void Node::graphPairs( string filename, NodeList &nodes )
{
    struct GraphCoords
    {
        ReadId id;
        int32_t coords[2];
    };
    
    sort( nodes.begin(), nodes.end(), []( Node* const &a, Node* const &b ) {
        return a->ends_[0] < b->ends_[0];
    });
    
    vector<GraphCoords> allGcs;
    vector<GraphCoords> gcs[params.libs.size()];
    int misCounts[3]{0};
    int32_t limits[2] = { 50000, -50000 };
    for ( Node* n : nodes )
    {
        size_t leader = n->seq_.find( "ATGGAGTTGAAAG" );
        if ( leader != n->seq_.npos )
        {
            limits[0] = min( limits[0], n->ends_[0] + (int)leader );
            limits[1] = max( limits[1], n->ends_[0] + (int)leader );
        }
    }
    limits[0] -= 1000;
    limits[1] += 1000;
    
    int orientCount = 0;
    int misorientCount = 0;
    for ( Node* n : nodes )
    {
        NodeOffsetMap offs = n->getDrxnNodesOffset( 1, 0, true );
        int nOff = 0;
        if ( n->edges_[1].size() == 1 )
        {
            nOff = n->getBiggestOffset( 1 );
        }
        for ( ReadMark &mark : n->marks_[0] )
        {
            int i = 0;
            while ( i < params.libs.size() && mark.id >= params.libs[i].endCount ) i++;
            assert( i < params.libs.size() );
            bool found = false, misoriented = false;
            GraphCoords gc;
            gc.coords[0] = mark.mark;
            gc.id = mark.id;
            for ( Node* t : nodes )
            {
                auto it = t->reads_.find( mark.id );
                if ( it == t->reads_.end() ) continue;
                int32_t pc[2] = { mark.mark, it->second[1] };
                auto it2 = offs.find( t );
                if ( t != n && it2 != offs.end() )
                {
                    int32_t diff = ( it2->second.first + it2->second.second ) / 2;
                    pc[0] += nOff;
                    pc[1] += diff - nOff;
                }
                if ( !found && it->second[0] < mark.mark )
                {
                    misoriented = true;
                    continue;
                }
                if ( found && abs( gc.coords[1] - gc.coords[0] - params.libs[i].size ) 
                            < abs( pc[1] - pc[0] - params.libs[i].size ) )
                {
                    continue;
                }
                found = true;
                misoriented = false;
                gc.coords[0] = pc[0];
                gc.coords[1] = pc[1];
            }
            
            if ( !found && !misoriented )
            {
                ReadId revId = params.getRevId( mark.id );
                for ( Node* t : nodes )
                {
                    auto it = t->reads_.find( revId);
                    if ( it == t->reads_.end() ) continue;
                    misoriented = true;
                    break;
                }
                
            }
            
            if ( ( limits[0] <= gc.coords[0] && gc.coords[0] <= limits[1] )
                    || ( limits[0] <= gc.coords[1] && gc.coords[1] <= limits[1] ) )
            {
                if ( found )
                {
                    gcs[i].push_back( gc );
                    allGcs.push_back( gc );
                    orientCount++;
                }
                else
                {
                    misCounts[i] += misoriented;
                    misorientCount += misoriented;
                }
            }
        }
    }
    
    vector<int32_t> markCoords;
    {
        int32_t coord = limits[0];
        while ( coord < limits[1] )
        {
            markCoords.push_back( coord );
            coord += 50;
        }
    }
    
    vector<int> graphMarks[params.libs.size()];
    vector<int> pairMarks;
    for ( int32_t coord : markCoords )
    {
        for ( int i = 0; i < params.libs.size(); i++ )
        {
            vector<int> dists;
            for ( GraphCoords &gc : gcs[i] )
            {
                if ( gc.coords[0] <= coord && coord <= gc.coords[1] ) dists.push_back( gc.coords[1] - gc.coords[0] );
            }
            sort( dists.begin(), dists.end() );
            int middle = dists.size() / 2;
            int median = dists[middle];
            if ( dists.size() % 2 && dists.size() > 1 )
            {
                median += dists[middle-1];
                median /= 2;
            }
            graphMarks[i].push_back( median );
        }
        int pairs = 0;
        for ( GraphCoords &gc : allGcs )
        {
            if ( gc.coords[0] <= coord && coord <= gc.coords[1] ) pairs++;
        }
        pairMarks.push_back( pairs );
    }
    
    ofstream csv( filename );
    csv << "Coordinates,550bp,3000bp,8000bp,PairCount\n"; 
    for ( int i = 0; i < markCoords.size(); i++ )
    {
        csv << to_string( markCoords[i] - limits[0] ) << "," 
                << to_string( graphMarks[0][i] ) << "," 
                << to_string( graphMarks[1][i] ) << "," 
                << to_string( graphMarks[2][i] ) << ","
                << to_string( pairMarks[i] ) << endl; 
    }
    
    csv.close();
    assert( false );
}

void Node::mapMates( Querier &bwt, int &count )
{
    for ( int i : { 0, 1 } )
    {
        NodeList tNodes = getTargetNodes( i, false );
        for ( ReadMark &mark : marks_[i] )
        {
            if ( params.isReadPe( mark.id ) ) continue;
            string seq = bwt.getSequence( mark.id );
//            if ( i )
//            {
//                size_t it = seq.find( "AGATGTGTATAAGAGACAG" );
//                if ( it == seq.npos ) continue;
//                seq = seq.substr( it + 19 );
//            }
//            else
//            {
//                size_t it = seq.find( "CTGTCTCTTATACACATCT" );
//                if ( it == seq.npos ) continue;
//                seq = seq.substr( 0, it );
//            }
//            if ( seq.size() < 45 ) continue;
            Node* best = NULL;
            int ol = 80;
            int32_t bestOffset;
            int32_t bestCoords[2];
            for ( Node* t : tNodes )
            {
                if ( t->reads_.find( mark.id ) != t->reads_.end() ) continue;
                
                int32_t coords[2];
                if ( mapSeqEnd( seq, t->seq_, ol, coords, i ) )
                {
                    coords[0] += t->ends_[0];
                    coords[1] += t->ends_[0];
                    if ( coords[0] < mark.coords[0] || mark.coords[1] < coords[1] ) continue;
                    int32_t offset = abs( coords[!i] - mark.estimate );
                    if ( best && ol == coords[1] - coords[0] )
                    {
                        best = NULL;
                        ol++;
                        continue;
                    }
                    best = t;
                    ol = coords[1] - coords[0];
                    bestCoords[0] = coords[0];
                    bestCoords[1] = coords[1];
                    bestOffset = offset;
                }
            }
            if ( best )
            {
                int ans = 1;
//                cout << i << " " << ol << " " << seq << endl;
//                cin >> ans;
                if ( ans )
                {
                    best->addRead( mark.id, bestCoords[0], bestCoords[1], true );
                    count++;
                }
            }
        }
    }
}

void Node::mergeAll( NodeList* nodes, NodeSet &delSet )
{
    for ( int i : { 2, 0, 1 } )
    {
        for ( Node* node : nodes[i] )
        {
            for ( int j : { 0, 1 } )
            {
                vector<int> ols;
                NodeList eNodes;
                for ( Edge &e : node->edges_[j] )
                {
                    if ( e.overlap <= 0 ) e.isLeap = true;
                    if ( e.isLeap ) continue;
                    Node* x[2];
                    x[0] = j ? node : e.node;
                    x[1] = j ? e.node : node;
                    string seqs[2] = { x[0]->seq_.substr( x[0]->seq_.length() - e.overlap ), x[1]->seq_.substr( 0, e.overlap ) };
                    if ( seqs[0] != seqs[1] )
                    {
                        int ol = mapSeqOverlap( x[0]->seq_, x[1]->seq_, 15 );
                        assert( false );
                        if ( !ol )
                        {
                            e.isLeap = true;
                            continue;
                        }
                        eNodes.push_back( e.node );
                        ols.push_back( ol );
                    }
                }
                for ( int k = 0; k < eNodes.size(); k++ )
                {
                    node->removeEdge( eNodes[k], j );
                    eNodes[k]->removeEdge( node, !j );
                    node->addEdge( eNodes[k], ols[k], j, false, false );
                }
            }
        }
    }
    
    for ( int i : { 2, 0, 1 } )
    {
        for ( Node* node : nodes[i] )
        {
            if ( delSet.find( node ) != delSet.end() ) continue;
            for ( int j : { 0, 1 } )
            {
                node->mergeDrxn( delSet, j );
            }
        }
    }
}

void Node::mergeDrxn( NodeSet &delSet, bool drxn )
{
    while ( edges_[drxn].size() == 1 )
    {
        if ( edges_[drxn][0].overlap <= 0 ) edges_[drxn][0].isLeap = true;
        if ( edges_[drxn][0].isLeap ) return;
        Node* node = edges_[drxn][0].node;
        if ( clones_ || node->clones_ || dontExtend_ || node->dontExtend_ ) return;
        int ol = edges_[drxn][0].overlap;
        if ( node->edges_[!drxn].size() != 1 ) return;
        int seqLen = node->seq_.length();
        
        for ( auto &read : node->reads_ )
        {
            if ( reads_.find( read.first ) != reads_.end() && !node->isRedundant( &read.second ) )
            {
                auto it = reads_.find( read.first );
                int32_t coords[2] = { it->second[0], it->second[1] };
                return;
            }
        }
        
        Node* l = drxn ? this : node;
        Node* r = drxn ? node : this;
        int32_t offset = drxn ? ends_[1] - ol - node->ends_[0]
                              : ends_[0] - node->ends_[1] + ol;
        assert( node->seq_.length() == node->ends_[1] - node->ends_[0] );
        node->clearPairs();
        string lSeq = l->seq_.substr( l->seq_.length() - ol );
        string rSeq = r->seq_.substr( 0, ol );
        assert( lSeq == rSeq );
        assert( l->seq_.substr( l->seq_.length() - ol ) == r->seq_.substr( 0, ol ) );
        string seq = l->seq_ + r->seq_.substr( ol, r->seq_.length() );
        seq_ = seq;
        if ( drxn ) ends_[1] += node->seq_.length() - ol;
        else ends_[0] -= ( node->seq_.length() - ol );
        clearEdges( drxn );
        for ( Edge &e : node->edges_[drxn] )
        {
            e.node->removeEdge( node, !drxn );
            addEdge( e.node, e.overlap, drxn, false, e.isLeap );
        }
        assert( ends_[1] - ends_[0] == seq_.length() ); 
        int32_t x[2] = { node->ends_[1], node->ends_[0] };
        for ( auto &read : node->reads_ )
        {
            x[0] = min( x[0], read.second[0] );
            x[1] = max( x[1], read.second[1] );
            addRead( read.first, read.second[0] + offset, read.second[1] + offset, read.second.redundant );
        }
        int32_t limits[2] = { ends_[1], ends_[0] };
        for ( auto &read : reads_ )
        {
            limits[0] = min( limits[0], read.second[0] );
            limits[1] = max( limits[1], read.second[1] );
        }
        
        assert( limits[0] == ends_[0] );
        if ( limits[1] < ends_[1] )
        {
            assert( edges_[1].empty() );
            ends_[1] = limits[1];
            seq_ = seq_.substr( 0, ends_[1] - ends_[0] );
        }
        node->dismantleNode();
        delSet.insert( node );
    }
}

void Node::recoil()
{
    int32_t ends[2] = { ends_[1], ends_[0] };
    assert( seq_.length() == ends_[1] - ends_[0] );
    for ( auto &read : reads_ )
    {
        ends[0] = min( ends[0], read.second[0] );
        ends[1] = max( ends[1], read.second[1] );
    }
    ends[0] = ends[0] - ends_[0];
    ends[1] = ends_[1] - ends[1];
    if ( ends[0] ) recoil( ends[0], 0 );
    if ( ends[1] ) recoil( ends[1], 1 );
    assert( seq_.length() == ends_[1] - ends_[0] );
}

void Node::recoil( int32_t diff, bool drxn )
{
    assert( diff > 0 );
    ends_[drxn] = drxn ? ends_[1] - diff : ends_[0] + diff;
    seq_ = drxn ? seq_.substr( 0, seq_.length() - diff ) : seq_.substr( diff );
    vector<int> ols;
    NodeList eNodes;
    for ( Edge &e : edges_[drxn] )
    {
        eNodes.push_back ( e.node ); 
        ols.push_back( e.overlap - diff );
    }
    clearEdges( drxn );
    for ( int i = 0; i < eNodes.size(); i++ )
    {
        this->addEdge( eNodes[i], ols[i], drxn, false, ols[i] <= 0 );
    }
}

void Node::remap( Querier &bwt )
{
    clearReads();
    vector<ReadId> ids;
    vector<int32_t> coords[2];
    bwt.mapSequence( seq_, ids, coords );
    for ( int i = 0; i < ids.size(); i++ )
    {
        addRead( ids[i], ends_[0] + coords[0][i], ends_[0] + coords[1][i], false );
    }
    setCoverage();
}

void Node::remapGenes( Querier &bwt, NodeList &nodes )
{
//    vector<ReadId> mapIds;
//    mapIds.push_back(342020973);
//    mapIds.push_back(427301160);
//    mapIds.push_back(517859060);
//    mapIds.push_back(610694375);
//    mapIds.push_back(484089616);
//    mapIds.push_back(649063290);
//    mapIds.push_back(325458350);
//    mapIds.push_back(500017413);
//    mapIds.push_back(432882647);
//    mapIds.push_back(673937720);
//    mapIds.push_back(445170277);
//    mapIds.push_back(414561562);
//    mapIds.push_back(625931473);
//    mapIds.push_back(390974329);
//    mapIds.push_back(647460820);
//    mapIds.push_back(390975725);
//    mapIds.push_back(346537812);
//    mapIds.push_back(424493340);
//    mapIds.push_back(491052479);
//    mapIds.push_back(547553316);
//    mapIds.push_back(483926524);
//    mapIds.push_back(455280565);
//    mapIds.push_back(519448264);
//    mapIds.push_back(510918341);
//    mapIds.push_back(724870756);
//    
//    for ( ReadId &id : mapIds )
//    {
//        string s = bwt.getSequence( id );
//        Node* bestNode = NULL;
//        int ol = 30;
//        int32_t bestCoords[2];
//        ReadId x = id % 4;
//        bool drxn = x == 1 || x == 2;
//        for ( Node* n : nodes )
//        {
//            int32_t coords[2];
//            if ( mapSeqEnd( s, n->seq_, ol, coords, drxn ) )
//            {
//                int thisOl = coords[1] - coords[0];
//                coords[0] += n->ends_[0];
//                coords[1] += n->ends_[0];
//                if ( thisOl == ol && bestNode->ends_[!drxn] == bestCoords[!drxn] )
//                {
//                    assert( false );
//                }
//                bestNode = n;
//                bestCoords[0] = coords[0];
//                bestCoords[1] = coords[1];
//                ol = thisOl;
//            }
//        }
//        bestNode->addRead( id, bestCoords[0], bestCoords[1], true );
//    }
//    
//    return;
//    assert( false );
    
    
    vector<string> seqs;
    seqs.push_back( "ATGGAGTTGAAAGGGACACTGATCGTTGTCATTGTGGCCTCTATTACCATTTCAG" );
    seqs.push_back( "ATGGAGTTGAAAGGGAAACTGATCGTTGTCATTGTGGCTGCTTTTACTATCTCAG" );
    seqs.push_back( "ATGGAGTTGAAAGGGACACTGATCGTTGTCATTGTGGCTGCTTTTACTATCTCAG" );
    seqs.push_back( "ATGGAGTTGAAAGTGAAACTGATCGTTGTTATTGTGGCTGCTATTGCAATCTCAG" );
    seqs.push_back( "ATGGAGTTGAAAGTGACACTGATCGTTGTTATTGTGGCTGCTATTGCAATCTCAG" );
    seqs.push_back( "TTCATGCACAAAGAGACGGGGGAGGAAGAGGAAATGGCAGAGAGAGGGGACAAGGCCGCTTTGGAGGAAGGCCACGACCTGATAGACCCCAGATGATGGGTGGACCTAGGCAAGGTGGTCCACCAATGGGCGGAAGGAGGTTTGATGTCCCTGGGCAAGGTGACCAACAGATGGATGGACGTGGACCGAATGGCGGGCCAATGGGTGGTAGGAGATTTGATAGACCAGGATTTGGTGGCTTCAGACCCGAAGGTGCCGGGAGACCTTTCTTCGGTCACGGAGGAAGGCATACTGATGGAGAAGGAGAAATGGAGGCTGCTCAACCAATTGGTGATGGTCAAGGATGGCTTGGTTTTTTCGATGGTACTGGAAGATTTTCCGGACGTCCTCACCCAGGCCGTGGTGATCATCATGGACACCACCATGGTCCTCCCCATGACCAGACCGACGAACACCCATTCGGTCAGCACAACGACAGCAACAGCGAGGAGGATGGCCGACCTCACCGTCACCACCACCACCACCACCACCATCATCACCATGACCGTCATAATGAGACAGACGACCACCGTCATCATAATCACACTGAAGGCCACCGCCACCATCATCATAATAAGACAGAAGAGGGTGACCAGGACAGACCAGAGATGAGGCCATTCCGGTTCAACCCTTTCGGTCGCAAGCCTTTCGGAGGACGTCCATTCGGCATGTTCGGCAGACGCAAACATACCGAAGAAGGATCTCACAGGCGCGATGGCCACCGTCATCCCCATGGCAACCGAGGACGTTGGGATGAGAATGAAGGTGAGGAGGAGGAGGAGGAACATCTTCCAACTGAAAACATGACAACATCGGCAGTGCCTGATGTGGTCGAGATCGACATCAACGAAATAGACAGCAACATTATCCCCGAGGTGTAG" );
    seqs.push_back( "TTCATGCACAAAGAGACGGGGGAGGAAGAGGAAATGGCAGAGAGAGGGGACAAGGCCGCTTTGGAGGAAGGCCACGACCTGATAGACCCCAGATGATGGGTGGACCTAGGCAAGGTGGTCCACCAATGGGCGGAAGGAGGTTTGATGTCCCTGGGCAAGGTGACCAACAGATGGATGGACGTGGACCGAATGGCAGGCCAATGGGTGGTAGGAGATTTGATAGACCAGGATTTGGTGGCTTCAGACCCGAAGGTGCCGGGAGACCTTTCTTCGGTCACGGAGGAAGGCATACTGATGGAGAAGGAGAAATGGAGGCTGCTCAACCAATTGGTGATGGTCAAGGATGGCTTGGTTTTTTCGATGGTACTGGAAGATTTTCCGGACGTCCTCACCCAGGCCGTGGTGATCATCATGGACACCACCATGGTCCTCCCCATGACCAGGCCGACGAACACCCATTCGGTCAGCACAACGACAGCAACAGCGAGGAGGATGGCCGACCTCACCGTCACCACCACCACCACCACCACCATCATCACCATGACCGTCATAATGAGACAGACGACCACCGTCATCATAATCACACTGAAGGCCACCGCCACCATCATCATAATAAGACAGAAGAGGGTGACCAGGACAGACCAGAGATGAGGCCATTCCGGTTCAACCCTTTCGGTCGCAAGCCTTTCGGAGGACGTCCATTCGGCATGTTCGGCAGACGCAAACATACCGAAGAAGGATCTCACAGGCGCGATGGCCACCGTCATCCCCATGGCAACCGAGGACGTTGGGATGAGAATGAAGGTGAGGAGGAGGAGGAGGAACATCTTCCAACTGAAAAAATGACAACATCGACAGTGCCTGATGTGGTCGAGATCGACATCAACGAAATAGACAGCAACATTATCCCCGAGGTGTAG" );
    seqs.push_back( "TTCATGCACAAAGAGAAAAGGCAAGAAGAGGGAATGGCAGAGAGAAGGAAGAGGGTCGCTTAAAAGGAAGGCAACGATCTGATAGACCCCAGATGATGGGTAGACCTAGGAAAGGTGGTCCACCAATGGGCGGAAGGGGGTTTGATGGCCCTGGACAAGGTGACCAACAGATGGGTGGACGTGTACCAAATGGCGGACCGATGGGCGGTAGGAGGTTTGATGGCCCTGGACAAGGTGACCAACAGATGGGTGGACGTGGACCAAATGGCGGACCGATGGGCGGTAGGAGGTTTGATGGCCCTGGACAAGGTGACCAACAGATGGGTTGACGTGGACCAAATGGATTTGGTGGCTTCAGACCCGAAGGTACAGGGAGACCTTTCTTCGGTCACGGAGGAAGGCACCGCCACCAAAATCACCACCATGACCGTCATAACAAGATAGACGATCACCGTCATCATAATCACACTGAAGGCCACCGTCATCATCATCATCATAACAAGACAGAAGAGGGTGACCAGGACAGACCAGAAATGTGGCCATTCCGGTTCAACAATACGGAAGAAGGATCTCCCAGGCGCGATGGACACCGTCATCCCAATGGCACTCGAGGACGTTGGGAGGAGAATGAAAGTGAGGAGGAAGGACATCTTTCGACTGAAAGCATGACAACATCTGCAGTGCCTGATGTGGTCGAGACCGACAGCAACGAAAAAGACAACATCATTATCCCTGAGGTGTAG" );
    seqs.push_back( "TTCATGCACAAAGAGAAAAGGCAAGAAGAGGGAATGGCAGAGAGAAGGAAGAGGGTCGCTTAAAAGGAAGGCAACGATCTGATAGACCCCAGATGATGGGTAGACCTAGGAAAGGTGGTCCACCAATGGGCGGAAGGGGGTTTGATGGCCCTGGACAAGGTGACCAACAGATGGGTGGACGTGTACCAAATGGCGGACCGATGGGCGGTAGGAGGTTAGATGGCCCTGGACAAGGTGACCAACAGATGGGTGGACGTGGACCAAATGGCGGACCGATGGGCGGTAGGAGGTTAGATGGCCCTGGACAAGGTGACCAACAGATGGGTGGACGTGGACCAAATGGATTTGGTGGCTTCAGACCCGAAGGTACAGGGAGACCTTTCTTCGGTCACGGAGGAAGGCACCGCCACCAAAATCACCACCATGACCGTCATAACAAGATAGACGATCACCGTCATCATAATCACACTGAAGGCCACCGTCATCATCATCATCATAACAAGACAGAAGAGGGTGACCAGGACAGACCAGAAATGTGGCCATTCCTGTTCAACCATACGGAAGAAGGATCTCCCAGGCGCGATGGACACCGTCATCCCAATGACACTCGAGGACGTTGGGATGAGAATGAAAGTGAGGAGGAAGGACATCTTTCGGCTGAAAGCATGACAACATCTGCAGTGCCTGATGTGGTCGAGACCGACAGCAACGAAAAAGACAACATCATTATCCCTGAAGTGTAG" );
    seqs.push_back( "TTCACGCCAAAGGAGAACGGAGAGGAAGAGGAAATGGCAGAGAGAGGGGAAAAGGTCGCGTCGGAGGAAGGCCAGGATCTGATAGACCCGAGATGATGGTTGGACCTATGCAAGGTGATCCACCAATGGGCGGAAGGAACTTTGATGGTTCTCCCCATGACCAGGCCGACAAACAACCATTCGGTCAACACAACGACAGCAGCAGCGAGGAGGATGGCCGACCTCACCGTCACCACCATGACCGTCATAATAAGACAGACGACCACCGTCATCATAATCACACTGAAGGCCACCGCCATCATCACCATAACCAGACAGAAGAGGGTGACCAGGACAGACCAGAGATGAGGCCATTCCGGTTCAACCCTTTCGGTCGCAAGCCCTTCGGAGGACGTCCATTCGGAAGACGCAACCATACCGAAGAAGGATCTCCTAGGCGCGGTGGACACCGTCAACGCAATGGCAACCGTGGACGTTGGGATGAGAATGAAAGCATGACAACATCTGCAGTGCCTGATGTGGTCGAGGTGTAG" );
    seqs.push_back( "TTCACGCCAAAGGAGAACGGAGAGGAAGAGGAAATGGCAGAGAGAGGGGAAAAGGTCGCGTCGGAGGAAGGCCAGGATCTGATAGACCCGAGATGATGGTTGGACCTAGGCAAGGTGATCCACCAATGGGCGGAAGGAACTTTGATGGTTCTCCCCATGATCAGGCCGACAAACAACCATTCGGTCAACACAACGACAGCAGCAGCGAGGAGGATGGCCGACCTCACCGTCACCACCATGACCGTCATAATAAGACAGACGACCACCGTCATCATAATCACACTGAAGGCCACCGCCATCATCACCATAACCAGACAGAAGAGGGTGACCAGGACAGACCAGAGATGAGGCCATTCCGGTTCAACCCTTTCGGTCGCAAGCCCTTCGGAGGACGTCCATTCGGAAGACGCAACCATACCGAAGAAGGATCTCCTAGGCGCGGTGGACACCGTCAACGCAATGGCAACCGTGGACGTTGGGATGAGAATGAAAGCATGACAACATCTGCAGTGCCTGATGTGGTCGAGGTGTAGTCACCCCC" );
    unordered_set<ReadId> ids;
    int i = 0;
    int geneCount = 0;
    int pseudoCount = 0;
    for ( string &seq : seqs )
    {
        MappedSeqs ms = bwt.mapSeed( seq, 10, false );
        cout << ">Seq" << to_string( i + 1 ) << endl;
        cout << string( 101, '-' ) << seq << endl;
        for ( ReadStruct &read : ms.reads )
        {
            bool docout = ids.find( read.readId ) == ids.end();
            if ( docout ) ids.insert( read.readId );
            for ( Node* n : nodes )
            {
                if ( n->reads_.find( read.readId ) == n->reads_.end() ) continue;
                docout = false;
            }
            docout = docout && read.seq.find( "TATAAGAGACAG" ) == read.seq.npos;
            docout = docout && read.seq.find( "CTGTCTCTTATA" ) == read.seq.npos;
            if ( docout )
            {
                cout << ">" << to_string( read.readId ) << endl;
                cout << string( 101 + read.coords[0], '-' ) << read.seq << endl;Node* best = NULL;
//                int ol = 50;
//                bool bestDrxn;
//                int32_t bestCoords[2];
//                Node* best = NULL;
//                for ( Node* n : nodes )
//                {
//                    for ( bool drxn : { 0, 1 } )
//                    {
//                        int32_t coords[2];
//                        if ( mapSeqEnd( read.seq, n->seq_, ol, coords, drxn ) )
//                        {
//                            coords[0] += n->ends_[0];
//                            coords[1] += n->ends_[0];
//                            if ( best && ol == coords[1] - coords[0] )
//                            {
//                                ol++;
//                                continue;
//                            }
//                            best = n;
//                            bestDrxn = drxn;
//                            ol = coords[1] - coords[0];
//                            bestCoords[0] = coords[0];
//                            bestCoords[1] = coords[1];
//                        }
//                    }
//                }
//                if ( best )
//                {
//                    ol = bestCoords[1] - bestCoords[0];
//                    string ext = bestDrxn ? read.seq.substr( 0, read.seq.length() - ol )
//                                          : read.seq.substr( ol );
//                    string base = "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG";
//                    if ( base.find( ext ) == base.npos && ext.find( "CTGTCTCTTATACACATCT" ) == ext.npos && ext.find( "AGATGTGTATAAGAGACAG" ) == ext.npos ) continue;
//                    if ( bestDrxn )
//                    {
//                        cout << read.seq.substr( 0, read.seq.length() - ol ) << endl;
//                        cout << read.seq << endl;
//                    }
//                    else
//                    {
//                        cout << string( ol, '-' ) << read.seq.substr( ol ) << endl;
//                        cout << read.seq << endl;
//                    }
//                    int ans = 1;
//                    cin >> ans;
//                    if ( ans )
//                    {
//                        best->addRead( read.readId, bestCoords[0], bestCoords[1], true );
//                    }
//                }
            }
        }
        i++;
    }
    
    int foundCount = 0;
    int notFoundCount = 0;
    int peBad = 0;
    int mpBad = 0;
    for ( ReadId id : ids )
    {
        bool found = false;
        for ( Node* node : nodes )
        {
            if ( node->reads_.find( id ) == node->reads_.end() ) continue;
            found = true;
            if ( node->ends_ [0] > 11000 ) pseudoCount++;
            else geneCount++;
            break;
        }
        if ( found ) foundCount++;
        else
        {
            string seq = bwt.getSequence( id );
            if ( seq.find( "TATAAGAGACAG" ) != seq.npos ) continue;
            if ( seq.find( "CTGTCTCTTATA" ) != seq.npos ) continue;
            if ( params.isReadPe( id ) ) peBad++;
            else mpBad++;
            notFoundCount++;
        }
    }
    assert( false );
}

