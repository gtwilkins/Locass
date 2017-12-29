/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "path_reassembly.h"
#include <algorithm>

Reassemble::Reassemble( Node* node, PathVars &pv, bool invalid )
: fork_( node ), invalid_( invalid )
{
    if ( invalid )
    {
        estLimits_[0] = node->ends_[!pv.drxn] - params.maxPeMean;
        estLimits_[1] = node->ends_[!pv.drxn] + params.maxPeMean;
        markLimits_[0] = node->ends_[!pv.drxn] + ( pv.drxn ? 0 : -params.readLen );
        markLimits_[1] = node->ends_[!pv.drxn] + ( pv.drxn ? params.readLen : 0 );
    }
    else
    {
        estLimits_[0] = pv.misassEst[0] + ( pv.drxn ? 0 : -params.readLen );
        estLimits_[1] = pv.misassEst[1] + ( pv.drxn ? params.readLen : 0 );
        markLimits_[0] = pv.misassMark[0] + ( pv.drxn ? -params.readLen : 0 );
        markLimits_[1] = pv.misassMark[1] + ( pv.drxn ? 0 : params.readLen );
    }
    
    int32_t tLimits[2] = { estLimits_[0] - params.readLen, estLimits_[1] + params.readLen };
    int32_t qLimits[2] = { tLimits[0] - params.maxPeMean, tLimits[1] + params.maxPeMean };
    int32_t pathLimits[2] = { node->ends_[0], node->ends_[1] };
    node->setOffsetMap( offsets_, pv.tSet, qLimits[0], 0 );
    node->setOffsetMap( offsets_, pv.tSet, qLimits[1], 1 );
    
    NodeSet usedSet;
    NodeListList paths( 1, { node } );
    for ( int i : { 0, 1 } )
    {
        int j = 0;
        while ( j < paths.size() )
        {
            pathLimits[0] = min( pathLimits[0], paths[j].back()->ends_[0] );
            pathLimits[1] = max( pathLimits[1], paths[j].back()->ends_[1] );
            auto it = offsets_.find( paths[j].back() );
            int32_t offset = it != offsets_.end() ? it->second : 0;
            bool doAdd = i ? paths[j].back()->ends_[1] - offset < tLimits[1]
                           : tLimits[0] < paths[j].back()->ends_[0] - offset;
            if ( usedSet.find( paths[j].back() ) != usedSet.end() ) doAdd = false;
            if ( paths[j].back()->drxn_ != pv.drxn ) doAdd = false;
            usedSet.insert( paths[j].back() );
            bool didAdd = false;
            for ( Node* nxt : paths[j].back()->getNextNodes( i ) )
            {
                if ( !doAdd || pv.tSet.find( nxt ) == pv.tSet.end() ) continue;
                if ( !didAdd ) paths[j].push_back( nxt );
                else
                {
                    NodeList path = { paths[j].end()[-2] };
                    path.push_back( nxt );
                    paths.push_back( path );
                }
                didAdd = true;
            }
            if ( !didAdd ) j++;
        }
        for ( int k = 0; k < paths.size(); k++ )
        {
            Node* anchor = paths[k][0];
            if ( !i ) reverse( paths[k].begin(), paths[k].end() );
            if ( !k ) continue;
            auto it = offsets_.find( anchor );
            assert( it != offsets_.end() );
            seqs_.push_back( new SeqPathReassemble( paths[k], anchor, ( it != offsets_.end() ? it->second : 0 ) ) );
        }
        paths.erase( paths.begin() + 1, paths.end() );
        usedSet.erase( node );
    }
    
    seqs_.push_back( new SeqPathReassemble( paths[0], node, 0 ) );
    estLimits_[0] = max( estLimits_[0], pathLimits[0] + params.readLen );
    estLimits_[1] = min( estLimits_[1], pathLimits[1] - params.readLen );
    assert( !seqs_.empty() );
    if ( estLimits_[1] < estLimits_[0] ) return;
    setWeakspot();
    map( pv );
}

Reassemble::~Reassemble()
{
    for ( int i = 0; i < seqs_.size(); i++ )
    {
        delete seqs_[i];
    }
}

bool Reassemble::doesSpanOrigin( PathVars &pv )
{
    int32_t limits[2] = { markLimits_[0] - params.readLen
                        , markLimits_[1] + params.readLen };
    limits[!pv.drxn] += pv.drxn ? -params.maxPeMean : params.maxPeMean;
    bool doesSpan = false;
    for ( SeqPathReassemble* s : seqs_ )
    {
        for ( Node* node : s->nodes )
        {
            if ( node->drxn_ == 2 )
            {
                if ( limits[0] < node->ends_[0] && node->ends_[1] < limits[1] ) doesSpan = true;
            }
        }
    }
    
    if ( doesSpan )
    {
        for ( Node* q : fork_->getDrxnNodes( pv.drxn , true, true ) )
        {
            for ( auto &np : q->pairs_ )
            {
                if ( np.first->drxn_ == 2 || np.first->drxn_ == !pv.drxn ) return true;
            }
        }
    }
    
    return false;
}

void Reassemble::map( PathVars &pv )
{
    NodeSet qSet = { fork_ };
    fork_->getDrxnNodes( qSet, 0, good_[0] - params.maxPeMean );
    fork_->getDrxnNodes( qSet, 1, good_[1] + params.maxPeMean );
    
    vector<string> seqs;
    vector<ReadMark> marks;
    vector< unordered_set<ReadId> > canIds( seqs_.size() );
    int32_t limits[2] = { estLimits_[0] - params.readLen * 2, estLimits_[1] + params.readLen * 2 };
    for ( Node* q : qSet )
    {
        if ( q->coverage_ > params.cover * 1.5 ) continue;
        auto it = offsets_.find( q );
        int32_t offset = it != offsets_.end() ? it->second : 0;
        for ( int i : { 0, 1 } )
        {
            // Set which sequences are reachable by this node
            vector<int> canSeq;
            NodeSet tSet = q->getDrxnNodesInSet( qSet, !i, true );
            for ( int j = 0; j < seqs_.size(); j++ )
            {
                bool thisCan = true;
                for ( Node* node : seqs_[j]->nodes )
                {
                    thisCan = thisCan && tSet.find( node ) != tSet.end();
                }
                thisCan = thisCan || find( seqs_[j]->nodes.begin(), seqs_[j]->nodes.end(), q ) != seqs_[j]->nodes.end();
                if ( thisCan ) canSeq.push_back( j ); 
            }
            
            for ( ReadMark mark : q->marks_[i] )
            {
                bool dontMap = false;
                if ( pv.usedIds.find( mark.id ) != pv.usedIds.end() ) continue;
                for ( Node* node : pv.nds[pv.drxn] ) if ( node->reads_.find( mark.id ) != node->reads_.end() ) dontMap = true;
                for ( Node* node : pv.nds[2] ) if ( node->reads_.find( mark.id ) != node->reads_.end() ) dontMap = true;
                if ( dontMap ) continue;
                mark.offset( -offset );
                if ( mark.estimate < limits[0] || limits[1] < mark.estimate ) continue;
                seqs.push_back( pv.bwt.getSequence( mark.id ) );
                marks.push_back( mark );
                for ( int j : canSeq ) canIds[j].insert( mark.id );
            }
        }
    }
    
    bool first = true;
    bool anyMap = true;
    while ( anyMap )
    {
        anyMap = false;
        for ( SeqPathReassemble* seq : seqs_ ) seq->complete = true;
        for ( int i = 0; i < seqs.size(); i++ )
        {
            int minOls[2] = { 16 + getHomopolymerLen( seqs[i], 0 ), 16 + getHomopolymerLen( seqs[i], 1 ) };
            for ( int j = 0; j < seqs_.size(); j++ )
            {
                if ( canIds[j].find( marks[i].id  ) == canIds[j].end() ) continue;
                if ( first ) seqs_[j]->map( seqs[i], marks[i], minOls );
                else seqs_[j]->remap( seqs[i], marks[i] );
            }
        }
        first = false;
        for ( SeqPathReassemble* seq : seqs_ ) anyMap = anyMap || !seq->complete;
    }
    
    for ( SeqPathReassemble* seq : seqs_ ) seq->setEdges();
}

bool Reassemble::reassemble( PathVars &pv, NodeSet &delSet )
{
    if ( !doesSpanOrigin( pv ) )
    {
        if ( estLimits_[1] <= estLimits_[0] )
        {
            fork_->dismantleNode( delSet, pv.drxn );
            return true;
        }
        if ( tryBridge( pv ) ) return true;
        if ( tryHalf( pv, delSet, !pv.drxn ) ) return true;
        if ( tryHalf( pv, delSet, pv.drxn ) ) return true;
        removeDubious( pv );
        if ( !invalid_ && tryGap( pv ) ) return true;
        if ( tryMap( pv ) ) return true;
        if ( trySlice( pv, delSet ) ) return true;
        assert( false );
    }
    
    return false;
}

void Reassemble::removeDubious( PathVars &pv )
{
    unordered_set<ReadId> badIds, usedIds;
    unordered_map<ReadId, int> ids[2];
    for ( SeqPathReassemble* s : seqs_ )
    {
        for ( int i : { 0, 1 } )
        {
            for ( ReadEndMap* read : s->reads[i] )
            {
                int ol = max( read->ol, read->edgeOl );
                auto r = ids[i].insert( make_pair( read->id, ol ) );
                if ( !r.second ) r.first->second = max( r.first->second, ol );
            }
        }
    }
    
    for ( SeqPathReassemble* s : seqs_ )
    {
        for ( int i : { 0, 1 } )
        {
            for ( ReadEndMap* read : s->reads[i] )
            {
                if ( ids[!i].find( read->id ) != ids[!i].end() ) continue;
                if ( read->edged || usedIds.find( read->id ) != usedIds.end() ) continue;
                usedIds.insert( read->id );
                int overlen = min( params.readLen * 0.8, params.readLen * 0.4 + getHomopolymerLen( read->seq, !i ) );
                if ( pv.bwt.isExtendable( read->seq, overlen, !i ) > 1 )
                {
                    badIds.insert( read->id );
                    while ( read->edge )
                    {
                        read = read->edge;
                        badIds.insert( read->id );
                    }
                }
            }
        }
    }
    
    for ( SeqPathReassemble* s : seqs_ )
    {
        for ( int i : { 0, 1 } )
        {
            s->paths[i].clear();
            for ( int j = 0; j < s->reads[i].size(); )
            {
                int ol = max( s->reads[i][j]->ol, s->reads[i][j]->edgeOl );
                int maxOl = ids[i][s->reads[i][j]->id ];
                if ( badIds.find( s->reads[i][j]->id ) != badIds.end()
                        || ol < maxOl )
                {
                    delete s->reads[i][j];
                    s->reads[i].erase( s->reads[i].begin() + j );
                } else j++;
            }
        }
    }
}

void Reassemble::removeRedundant( vector<ReadEndMap*> &reads )
{
    for ( int i = 0; i < reads.size(); i++ )
    {
        for ( int j = i+1; j < reads.size(); j++ )
        {
            if ( reads[i] == reads[j] ) assert( false );
            if ( reads[i]->id != reads[j]->id ) continue;
            if ( reads[i]->ol > reads[j]->ol ) reads[i]->doMap = false;
            else reads[j]->doMap = false;
        }
    }
}

void Reassemble::setWeakspot()
{
    NodeSet nodes;
    for ( SeqPathReassemble* seq : seqs_ ) nodes.insert( seq->nodes.begin(), seq->nodes.end() );
    
    int len = estLimits_[1] - estLimits_[0];
    int cover[len]{0};
    
    for ( Node* node : nodes )
    {
        int32_t offset = offsets_[node];
        if ( estLimits_[0] < node->ends_[1] + offset && node->ends_[0] + offset < estLimits_[1] )
        {
            for ( auto &read : node->reads_ )
            {
                int j = max( estLimits_[0], read.second[0] + offset ) - estLimits_[0];
                int k = min( estLimits_[1], read.second[1] + offset ) - estLimits_[0];
                while ( j < k )
                {
                    assert( j >= 0 && j < len );
                    cover[j++]++;
                }
            }
        }
    }
    
    int minCover = cover[0];
    int32_t coord;
    int32_t goodCover[3] = {0,len,0};
    int32_t highCover[3] = {0,len,0};
    int goodCutoff = params.cover * 0.6, highCutoff = params.cover * 1.2;
    for ( int i = 0; i < len; i++ )
    {
        if ( cover[i] > goodCutoff )
        {
            if ( i < goodCover[1] ) goodCover[1] = i;
            goodCover[2] = i;
        }
        
        if ( cover[i] > highCutoff )
        {
            if ( i < highCover[1] ) highCover[1] = i;
            highCover[2] = i;
        }
        
        if ( cover[i] < minCover )
        {
            if ( cover[i] < goodCutoff )
            {
                goodCover[0] = goodCover[2];
                goodCover[1] = len;
            }
            if ( cover[i] < highCutoff )
            {
                highCover[0] = highCover[2];
                highCover[1] = len;
            }
            minCover = cover[i];
            coord = i;
        }
    }
    
    est_ = coord + estLimits_[0];
    minCover_ = minCover;
    allHigh_ = !highCover[1];
    allGood_ = !goodCover[1];
    high_[0] = estLimits_[0] + highCover[0];
    high_[1] = estLimits_[0] + highCover[1];
    good_[0] = estLimits_[0] + goodCover[0];
    good_[1] = estLimits_[0] + goodCover[1];
}

bool Reassemble::tryBridge( PathVars &pv )
{
    MapNode* mn = new MapNode();
    
    float bestScore = 0;
    unordered_set<ReadId> usedIds; 
    for ( SeqPathReassemble* s : seqs_ )
    {
        s->setBridges( seqs_, mn, usedIds, bestScore, !pv.drxn, ( invalid_ ? 2 : pv.drxn ) );
        s->setBridges( seqs_, mn, usedIds, bestScore, pv.drxn, ( invalid_ ? 2 : pv.drxn ) );
    }
    
    bool didBridge = false;
    if ( mn->seq.length() >= params.readLen )
    {
        didBridge = Node::mapBridge( fork_, pv, mn );
        if ( didBridge ) pv.usedIds.insert( usedIds.begin(), usedIds.end() );
    }
    delete mn;
    return didBridge;
}

bool Reassemble::tryGap( PathVars &pv )
{
    int32_t limits[2] = { estLimits_[0], estLimits_[1] };
    limits[0] = max( estLimits_[0], good_[0] );
    limits[1] = min( estLimits_[1], good_[1] );
    if ( allGood_ || limits[1] < limits[0] ) return false;
    
    unordered_set<ReadId> ids[2];
    vector<int> anchors[2];
    int ante[2] = { 0 };
    int pro[2] = { 0 };
    int olTotal = 0;
    for ( int i : { !pv.drxn, pv.drxn } )
    {
        for ( SeqPathReassemble* s : seqs_ )
        {
            for ( ReadEndMap* read : s->reads[i] )
            {
                if ( ids[i].find( read->id ) != ids[i].end() ) continue;
                if ( read->coords[2] + s->ends[0] < limits[0] || limits[1] < read->coords[0] + s->ends[0] ) continue;
                anchors[i].push_back( read->coords[1] + s->ends[0] );
                ( read->drxn == i ? pro[i] : ante[i] )++;
                ids[i].insert( read->id );
                olTotal += max( read->ol, read->edgeOl );
            }
        }
    }
    
    if ( !ante[0] || !ante[1] 
            || olTotal < params.readLen 
            || ( pro[pv.drxn] + ante[!pv.drxn] < 3 ) ) return false;
    
    for ( int i : { 0, 1 } ) sort( anchors[i].begin(), anchors[i].end() );
    
    int mapAnchors[2] = { ( anchors[0][ anchors[0].size() / 2 ] 
                          + anchors[0][ ( anchors[0].size()-1 ) / 2 ] ) / 2
                        , ( anchors[1][ anchors[1].size() / 2 ] 
                          + anchors[1][ ( anchors[1].size()-1 ) / 2 ] ) / 2 };
    
    if ( mapAnchors[1] < mapAnchors[0] )
    {
        int tmpAnchors[2] = { mapAnchors[0], mapAnchors[1] };
        for ( int i = 0; i < anchors[0].size(); i++ )
        {
            if ( anchors[0][i] < tmpAnchors[1] ) mapAnchors[0] = anchors[0][i];
        }
        for ( int i = anchors[1].size(); --i >= 0; )
        {
            if ( tmpAnchors[0] < anchors[1][i] ) mapAnchors[1] = anchors[1][i];
        }
    }
    
    vector<ReadEndMap*> mapReads;
    for ( int i : { 0, 1 } )
    {
        for ( SeqPathReassemble* s : seqs_ )
        {
            for ( ReadEndMap* read : s->reads[i] )
            {
                if ( read->coords[2] + s->ends[0] < mapAnchors[0] 
                        || mapAnchors[1] < read->coords[0] + s->ends[0] ) continue;
                mapReads.push_back( read );
                read->doMap = true;
            }
        }
    }
    
    removeRedundant( mapReads );
    bool didMap = false;
    for ( SeqPathReassemble* s : seqs_ ) didMap = s->doMap( pv ) || didMap;
    
    return didMap;
}

bool Reassemble::tryHalf( PathVars &pv, NodeSet &delSet, bool drxn )
{
    MapNode* mn = new MapNode();
    
    int bestScore = 0;
    unordered_set<ReadId> usedIds; 
    for ( SeqPathReassemble* s : seqs_ )
    {
        s->setHalves( mn, usedIds, bestScore, drxn );
    }
    
    bool didBridge = false;
    if ( mn->seq.length() >= params.readLen )
    {
        pv.bwt.mapSequence( mn->seq, mn->ids, mn->coords );
        mn->recoil();
        Node* node = new Node( mn, 0, mn->ids.size()-1, pv.drxn );
        pv.usedIds.insert( usedIds.begin(), usedIds.end() );
        
        int32_t coords[2] = { mn->bridgeCoords[drxn][0], mn->bridgeCoords[drxn][0] };
        coords[drxn] = drxn ? coords[0] + mn->bridgeOverlaps[1][0] : coords[1] - mn->bridgeOverlaps[0][0];
        int32_t offset = coords[drxn] - node->ends_[drxn];
        node->offset( offset );
        
        if ( drxn == pv.drxn )
        {
            didBridge = node->sliceOrBridge( pv, mn->bridges[drxn][0], coords, delSet );
        }
        else
        {
            NodeList hitNodes;
            vector<int32_t> hitCoords[2];
            mn->bridges[drxn][0]->overlapExtend( pv.nds[pv.drxn], coords, hitNodes, hitCoords, pv.drxn, drxn );
            pv.nds[pv.drxn].push_back( node );
            pv.newSet.insert( node );
            int32_t dummy[2] = { params.locusLimits[0], params.locusLimits[1] };
            
            for ( int i = 0; i < hitNodes.size(); i++ )
            {
                hitNodes[i]->addEdge( node, hitCoords[1][i] - hitCoords[0][i], pv.drxn );
            }
            node->inheritEdges( drxn );
            node->propagateValidation( dummy, pv.drxn );
            didBridge = true;
        }
    }
    
    delete mn;
    return didBridge;
}

bool Reassemble::tryMap( PathVars &pv )
{
    if ( allHigh_ ) return false;
    
    vector<ReadEndMap*> mapReads;
    vector<int> anchors[2];
    unordered_set<ReadId> hitIds[2];
    for ( int i : { 0, 1 } )
        for ( SeqPathReassemble* s : seqs_ )
            for ( ReadEndMap* read : s->reads[i] )
                hitIds[i].insert( read->id );
        
    for ( int i : { 0, 1 } )
    {
        for ( SeqPathReassemble* s : seqs_ )
        {
            for ( ReadEndMap* read : s->reads[i] )
            {
                bool isBridge = hitIds[!i].find( read->id ) != hitIds[!i].end();
                if ( ( read->drxn == i && !isBridge )
                        || s->ends[0] + read->coords[2] < high_[0] 
                        || high_[1] < s->ends[0] + read->coords[0] ) continue;
                if ( !isBridge && read->ol < params.readLen * 0.4 ) continue;
                if ( abs( read->off ) > read->ol ) continue;
                anchors[i].push_back( s->ends[0] + read->coords[1] );
                read->doMap = true;
                mapReads.push_back( read );
            }
        }
    }
    
    int cutoff = params.maxPeMean - params.readLen;
    for ( int i = 0; i < anchors[0].size(); )
    {
        bool isGood = false;
        for ( int j = 0; j < anchors[1].size(); j++ )
        {
            if ( abs( anchors[0][i] - anchors[1][j] ) > cutoff ) continue;
            if ( anchors[0][i] <= anchors[1][j] ) isGood = true;
        }
        if ( !isGood ) anchors[0].erase( anchors[0].begin() + i );
        else i++;
    }
    
    for ( int j = 0; j < anchors[1].size(); )
    {
        bool isGood = false;
        for ( int i = 0; i < anchors[0].size(); i++ )
        {
            if ( abs( anchors[0][i] - anchors[1][j] ) > cutoff ) continue;
            if ( anchors[0][i] <= anchors[1][j] ) isGood = true;
        }
        if ( !isGood ) anchors[1].erase( anchors[1].begin() + j );
        else j++;
    }
    
    anchors[0].insert( anchors[0].end(), anchors[1].begin(), anchors[1].end() );
    if ( anchors[0].empty() ) return false;
    
    for ( int i : { 0, 1 } )
    {
        for ( SeqPathReassemble* s : seqs_ )
        {
            for ( ReadEndMap* read : s->reads[i] )
            {
                if ( read->doMap 
                        || s->ends[0] + read->coords[1] < high_[0] 
                        || high_[1] < s->ends[0] + read->coords[1] ) continue;
                int anchorSpans = 0;
                for ( int anchor : anchors[0] ) anchorSpans += ( s->ends[0] + read->coords[0] <= anchor
                                                            && anchor <= s->ends[0] + read->coords[2] );
                int score = read->ol + min( 15, anchorSpans * 3 ) + ( read->drxn != i ? 10 : 0 );
                if ( s->ends[0] + read->coords[1] < good_[0] 
                        || good_[1] < s->ends[0] + read->coords[1] ) score -= 10;
                if ( read->off > score || score < params.readLen * 0.45 ) continue;
                read->doMap = true;
                mapReads.push_back( read );
            }
        }
    }
    
    removeRedundant( mapReads );
    bool didMap = false;
    for ( SeqPathReassemble* s : seqs_ ) didMap = s->doMap( pv ) || didMap;
    
    return didMap;
}

bool Reassemble::trySlice( PathVars &pv, NodeSet &delSet )
{
//    bool misassembled = fork_->isMisassembledRev( markLimits_, pv.drxn );
//    if ( !misassembled )
//    {
//        if ( tryMap( pv ) ) return true;
//        assert( false );
//    }
//    if ( fork_->slice( pv, misassembled, pv.drxn ) ) return true;
    int32_t coord = markLimits_[!pv.drxn] + ( pv.drxn ? -params.readLen / 2 : params.readLen / 2 );
    coord = pv.drxn ? max( coord, fork_->ends_[0] ) : min( coord, fork_->ends_[1] );
    Node* node = fork_;
    if ( fork_->getNextReadCoord( coord, !pv.drxn, pv.drxn ) )
    {
        node = node->splitNode( coord, pv.nds[pv.drxn], pv.drxn, pv.drxn );
    }
    
    node->dismantleNode( delSet, pv.drxn );
//    else
//    {
//        node->setUnreliable();
//        assert( false );
//    }
    return true;
}
