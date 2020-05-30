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
                    ols.push_back( l.ol - ( seq_.length() - r.ol ) );
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

//void Node::graphPairs( string filename, NodeList &nodes )
//{
//    struct GraphCoords
//    {
//        ReadId id;
//        int32_t coords[2];
//    };
//    
//    sort( nodes.begin(), nodes.end(), []( Node* const &a, Node* const &b ) {
//        return a->ends_[0] < b->ends_[0];
//    });
//    
//    vector<GraphCoords> allGcs;
//    vector<GraphCoords> gcs[params.libs.size()];
//    int misCounts[3]{0};
//    int32_t limits[2] = { 50000, -50000 };
//    for ( Node* n : nodes )
//    {
//        size_t leader = n->seq_.find( "ATGGAGTTGAAAG" );
//        if ( leader != n->seq_.npos )
//        {
//            limits[0] = min( limits[0], n->ends_[0] + (int)leader );
//            limits[1] = max( limits[1], n->ends_[0] + (int)leader );
//        }
//    }
//    limits[0] -= 1000;
//    limits[1] += 1000;
//    
//    int orientCount = 0;
//    int misorientCount = 0;
//    for ( Node* n : nodes )
//    {
//        NodeOffsetMap offs = n->getDrxnNodesOffset( 1, 0, true );
//        int nOff = 0;
//        if ( n->edges_[1].size() == 1 )
//        {
//            nOff = n->getBiggestOffset( 1 );
//        }
//        for ( ReadMark &mark : n->marks_[0] )
//        {
//            int i = 0;
//            while ( i < params.libs.size() && mark.id >= params.libs[i].endCount ) i++;
//            assert( i < params.libs.size() );
//            bool found = false, misoriented = false;
//            GraphCoords gc;
//            gc.coords[0] = mark.mark;
//            gc.id = mark.id;
//            for ( Node* t : nodes )
//            {
//                auto it = t->reads_.find( mark.id );
//                if ( it == t->reads_.end() ) continue;
//                int32_t pc[2] = { mark.mark, it->second[1] };
//                auto it2 = offs.find( t );
//                if ( t != n && it2 != offs.end() )
//                {
//                    int32_t diff = ( it2->second.first + it2->second.second ) / 2;
//                    pc[0] += nOff;
//                    pc[1] += diff - nOff;
//                }
//                if ( !found && it->second[0] < mark.mark )
//                {
//                    misoriented = true;
//                    continue;
//                }
//                if ( found && abs( gc.coords[1] - gc.coords[0] - params.libs[i].size ) 
//                            < abs( pc[1] - pc[0] - params.libs[i].size ) )
//                {
//                    continue;
//                }
//                found = true;
//                misoriented = false;
//                gc.coords[0] = pc[0];
//                gc.coords[1] = pc[1];
//            }
//            
//            if ( !found && !misoriented )
//            {
//                ReadId revId = params.getRevId( mark.id );
//                for ( Node* t : nodes )
//                {
//                    auto it = t->reads_.find( revId);
//                    if ( it == t->reads_.end() ) continue;
//                    misoriented = true;
//                    break;
//                }
//                
//            }
//            
//            if ( ( limits[0] <= gc.coords[0] && gc.coords[0] <= limits[1] )
//                    || ( limits[0] <= gc.coords[1] && gc.coords[1] <= limits[1] ) )
//            {
//                if ( found )
//                {
//                    gcs[i].push_back( gc );
//                    allGcs.push_back( gc );
//                    orientCount++;
//                }
//                else
//                {
//                    misCounts[i] += misoriented;
//                    misorientCount += misoriented;
//                }
//            }
//        }
//    }
//    
//    vector<int32_t> markCoords;
//    {
//        int32_t coord = limits[0];
//        while ( coord < limits[1] )
//        {
//            markCoords.push_back( coord );
//            coord += 50;
//        }
//    }
//    
//    vector<int> graphMarks[params.libs.size()];
//    vector<int> pairMarks;
//    for ( int32_t coord : markCoords )
//    {
//        for ( int i = 0; i < params.libs.size(); i++ )
//        {
//            vector<int> dists;
//            for ( GraphCoords &gc : gcs[i] )
//            {
//                if ( gc.coords[0] <= coord && coord <= gc.coords[1] ) dists.push_back( gc.coords[1] - gc.coords[0] );
//            }
//            sort( dists.begin(), dists.end() );
//            int middle = dists.size() / 2;
//            int median = dists[middle];
//            if ( dists.size() % 2 && dists.size() > 1 )
//            {
//                median += dists[middle-1];
//                median /= 2;
//            }
//            graphMarks[i].push_back( median );
//        }
//        int pairs = 0;
//        for ( GraphCoords &gc : allGcs )
//        {
//            if ( gc.coords[0] <= coord && coord <= gc.coords[1] ) pairs++;
//        }
//        pairMarks.push_back( pairs );
//    }
//    
//    ofstream csv( filename );
//    csv << "Coordinates,550bp,3000bp,8000bp,PairCount\n"; 
//    for ( int i = 0; i < markCoords.size(); i++ )
//    {
//        csv << to_string( markCoords[i] - limits[0] ) << "," 
//                << to_string( graphMarks[0][i] ) << "," 
//                << to_string( graphMarks[1][i] ) << "," 
//                << to_string( graphMarks[2][i] ) << ","
//                << to_string( pairMarks[i] ) << endl; 
//    }
//    
//    csv.close();
//    assert( false );
//}

//void Node::mapMates( Querier &bwt, int &count )
//{
//    for ( int i : { 0, 1 } )
//    {
//        NodeList tNodes = getTargetNodes( i, false );
//        for ( ReadMark &mark : marks_[i] )
//        {
//            if ( params.isReadPe( mark.id ) ) continue;
//            string seq = bwt.getSequence( mark.id );
////            if ( i )
////            {
////                size_t it = seq.find( "AGATGTGTATAAGAGACAG" );
////                if ( it == seq.npos ) continue;
////                seq = seq.substr( it + 19 );
////            }
////            else
////            {
////                size_t it = seq.find( "CTGTCTCTTATACACATCT" );
////                if ( it == seq.npos ) continue;
////                seq = seq.substr( 0, it );
////            }
////            if ( seq.size() < 45 ) continue;
//            Node* best = NULL;
//            int ol = 80;
//            int32_t bestOffset;
//            int32_t bestCoords[2];
//            for ( Node* t : tNodes )
//            {
//                if ( t->reads_.find( mark.id ) != t->reads_.end() ) continue;
//                
//                int32_t coords[2];
//                if ( mapSeqEnd( seq, t->seq_, ol, coords, i ) )
//                {
//                    coords[0] += t->ends_[0];
//                    coords[1] += t->ends_[0];
//                    if ( coords[0] < mark.coords[0] || mark.coords[1] < coords[1] ) continue;
//                    int32_t offset = abs( coords[!i] - mark.estimate );
//                    if ( best && ol == coords[1] - coords[0] )
//                    {
//                        best = NULL;
//                        ol++;
//                        continue;
//                    }
//                    best = t;
//                    ol = coords[1] - coords[0];
//                    bestCoords[0] = coords[0];
//                    bestCoords[1] = coords[1];
//                    bestOffset = offset;
//                }
//            }
//            if ( best )
//            {
//                int ans = 1;
////                cout << i << " " << ol << " " << seq << endl;
////                cin >> ans;
//                if ( ans )
//                {
//                    best->addRead( mark.id, bestCoords[0], bestCoords[1], true );
//                    count++;
//                }
//            }
//        }
//    }
//}

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
                    if ( e.ol <= 0 ) e.leap = true;
                    if ( e.leap ) continue;
                    Node* x[2];
                    x[0] = j ? node : e.node;
                    x[1] = j ? e.node : node;
                    string seqs[2] = { x[0]->seq_.substr( x[0]->seq_.length() - e.ol ), x[1]->seq_.substr( 0, e.ol ) };
                    if ( seqs[0] != seqs[1] )
                    {
                        int ol = mapSeqOverlap( x[0]->seq_, x[1]->seq_, 15 );
                        assert( false );
                        if ( !ol )
                        {
                            e.leap = true;
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

//Node* Node::merge( bool drxn )
//{
//    if ( edges_[drxn].size() != 1 ) return NULL;
//    Edge& e = edges_[drxn][0];
//    if ( e.ol <= 0 ) e.isLeap = true;
//    if ( e.isLeap || e.node->edges_[!drxn].size() != 1 ) return NULL;
//    if ( e.node->cloned_ ) return NULL;
//    if ( cloned_ ) for ( Node* c : cloned_->nodes ) if ( !c->edges_[drxn].empty() ) return NULL;
//    
//    Node* node = e.node,* l = drxn ? this : e.node, * r = drxn ? e.node : this;
//    if ( !bad_ && node->bad_ ) claimGood( drxn );
//    if ( bad_ && !node->bad_ ) node->claimGood( !drxn );
//    
//    int32_t off = drxn ? ends_[1] - e.ol - node->ends_[0] : ends_[0] - node->ends_[1] + e.ol;
//    node->offset( off );
//    if ( node->drxn_ >= 2 ) drxn_ = node->drxn_;
//    bad_ = bad_ && node->bad_;
//    mapped_ = mapped_ && node->mapped_;
//    ends_.merge( node->ends_ );
//    
//    clearPaired( false );
//    node->clearPaired( false );
//    verified_ = false;
//    string lSeq = l->seq_.substr( l->seq_.length() - e.ol );
//    string rSeq = r->seq_.substr( 0, e.ol );
//    assert( lSeq == rSeq );
//
//    string seq = l->seq_ + r->seq_.substr( e.ol, r->seq_.length() );
//    seq_ = seq;
//    stop_[drxn] = node->stop_[drxn];
//
//    clearEdges( drxn );
//    for ( Edge &edge : node->edges_[drxn] )
//    {
//        edge.node->removeEdge( node, !drxn );
//        addEdge( edge.node, edge.ol, drxn, false, edge.isLeap );
//    }
//    node->clearEdges( drxn );
//    assert( ends_[1] - ends_[0] == seq_.length() ); 
//
//    reads_.insert( node->reads_.begin(), node->reads_.end() );
//    setCoverage();
//    remark();
//    readTest();
//    
//    if ( !cloned_ ) return node;
//    
//    for ( Node* c : cloned_->nodes )
//    {
//        c->seq_ = seq_;
//        c->clearPaired( false );
//        c->verified_ = false;
//        c->ends_[drxn] = drxn ? c->ends_[0] + size() : c->ends_[1] - size();
//        c->reads_ = reads_;
//        off = c->ends_[0] - ends_[0];
//        for ( auto& read : c->reads_ ) read.second.offset( off );
//        c->setCoverage();
//        c->remark();
//        c->readTest();
//    }
//    
//    return node;
//}

void Node::mergeDrxn( NodeSet &delSet, bool drxn )
{
    while ( edges_[drxn].size() == 1 )
    {
        if ( edges_[drxn][0].ol <= 0 ) edges_[drxn][0].leap = true;
        if ( edges_[drxn][0].leap ) return;
        Node* node = edges_[drxn][0].node;
        if ( clones_ || node->clones_ || dontExtend_ || node->dontExtend_ ) return;
        int ol = edges_[drxn][0].ol;
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
            addEdge( e.node, e.ol, drxn, false, e.leap );
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
    assert( !reads_.empty() );
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
    vector<Edge> edges = edges_[drxn];
    clearEdges( drxn );
    for ( Edge& e : edges ) addEdge( e.node, e.ol-diff, drxn, false, e.leap || e.ol <= 0 );
}

void Node::recoordinate( NodeRoll& nodes )
{
    Nodes tested;
    for ( Node* node : nodes.getGraph( 2 ).nodes ) node->recoordinate( tested, 1 );
    
    int32_t limits[2]{0}, counts[2]{0}, verified[2]{0};
    for ( Node* node : tested.nodes )
    {
        if ( node->cloned_ ) continue;
        if ( node->verified_ ) verified[0] = min( verified[0], node->ends_[0] );
        if ( node->verified_ ) verified[1] = max( verified[1], node->ends_[1] );
        if ( node->isContinue( 0 ) ) limits[0] = min( limits[0], node->ends_[0] );
        if ( node->isContinue( 1 ) ) limits[1] = max( limits[1], node->ends_[1] );
        for ( int d : { 0, 1 } ) if ( node->edges_[d].empty() ) counts[d]++;
    }
    
    cout << "Left ends: " << counts[0] << ", right ends: " << counts[1] << endl;
    cout << "Left verified min: " << verified[0] << ", right verified max: " << verified[1] << endl;
    cout << "Left min: " << limits[0] << ", right max: " << limits[1] << endl;
}

void Node::recoordinate( Nodes& tested, bool drxn )
{
    if ( drxn_ < 2 ) for ( Edge& e : edges_[!drxn] ) if ( !e.node->bad_ && !tested.find( e.node ) ) return;
    if ( !tested.add( this ) ) return;
    assert( !bad_ );
    if ( drxn_ < 2 )
    {
        bool isset = false;
        int32_t best = 0;
        for ( Edge& e : edges_[!drxn] )
        {
            if ( e.node->bad_ ) continue;
            int32_t off = e.node->ends_[drxn] - ends_[!drxn] + ( drxn ? -e.ol : e.ol );
            if ( off )
            {
                int x = 0;
            }
            if ( !isset || ( drxn ? best < off : off < best ) ) best = off;
            isset = true;
        }
        if ( best )
        {
            int x = 0;
        }
        assert( isset );
        if ( best ) offset( best );
        
        for ( Edge& e : edges_[drxn] ) e.node->recoordinate( tested, drxn );
    }
    else for ( int d : { 0, 1 } ) for ( Edge& e : edges_[d] ) e.node->recoordinate( tested, d );
}

void Node::remap( Querier& bwt, NodeRoll& nodes )
{
    NodeRoll remapped;
    for ( Node* node : nodes.nodes )
    {
        if ( !node->verified_ || node->mapped_ ) continue;
        node->remap( bwt );
        remapped += node;
//        if ( node->cloned_ ) for ( Node* clone : node->cloned_->nodes ) remapped += clone;
    }
    if ( remapped.empty() ) return;
    NodeRoll tar, ignored;
    for ( Node* node : nodes.nodes ) if ( node->verified_ ) tar += node;
    for ( Node* node : remapped.nodes ) if ( node->remap( tar ) ) ignored += node;
    for ( Node* node : ignored.nodes ) node->clearPaired( true );
    for ( Node* node : remapped.nodes ) node->setVerified();
}

bool Node::remap( Querier &bwt )
{
    if ( mapped_ ) return false;
    vector<ReadId> ids;
    vector<int32_t> coords[2];
    bwt.mapSequence( seq_, ids, coords );
    int remapped = 0;
    for ( int i = 0; i < ids.size(); i++ ) if ( reads_.find( ids[i] ) == reads_.end() )
    {
        bool redundant = isRedundant( ends_[0] + coords[0][i], ends_[0] + coords[1][i] );
        add( ids[i], ends_[0] + coords[0][i], ends_[0] + coords[1][i], redundant );
        if ( cloned_ ) for ( Node* clone : cloned_->nodes ) clone->add( ids[i], clone->ends_[0] + coords[0][i], clone->ends_[0] + coords[1][i], redundant );
        remapped++;
    }
    if ( !remapped ) return false;
    setCoverage();
    mapped_ = true;
    if ( cloned_ ) for ( Node* clone : cloned_->nodes ) clone->mapped_ = true;
    cout << "Remapped " << reads_.size() << " " << remapped << endl;
    return true;
}

bool Node::remap( NodeRoll& tar )
{
    int ignored = 0;
    for ( auto& read : reads_ )
    {
        if ( !read.second.redundant || read.second.ignore ) continue;
        for ( Node* node : tar.nodes )
        {
            auto it = node->reads_.find( read.first );
            if ( it == node->reads_.end() || node == this ) continue;
            if ( cloned_ && cloned_->find( node ) ) continue;
            read.second.ignore = true;
            it->second.ignore = true;
            ignored++;
        }
        if ( !read.second.ignore || !cloned_ ) continue;
        for ( Node* clone : cloned_->nodes )
        {
            auto it = clone->reads_.find( read.first );
            assert( it != clone->reads_.end() );
            it->second.ignore = true;
        }
    }
    cout << "Ignored " << reads_.size() << " " << ignored << endl;
    return ignored;
}
