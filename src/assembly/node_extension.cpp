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

void Node::addExtension( Extension &ext, ExtVars &ev, vector<MergeHit> &selfMerges, bool doesBranch, bool drxn )
{
    MergeHit merge;
    for ( Node* &node : ev.nodes )
    {
        if ( node->checkExtensionMerge( ext, merge ) ) break;
    }
    for ( Node* &node : ev.island )
    {
        if ( node->checkExtensionMerge( ext, merge ) ) break;
    }
    
    if ( merge.node )
    {
        if ( ev.ante.find( merge.node ) != ev.ante.end() && !ev.doLoop )
        {
            merge.node = NULL;
        }
        else if ( isMergeSelf( ev, merge ) )
        {
            if ( merge.overlap == ext.maxOverLen && ev.doLoop  )
            {
                selfMerges.push_back( merge );
            }
            else
            {
                merge.node = NULL;
            }
        }
        else if ( edges_[drxn].empty() 
                || merge.node->drxn_ > 2 
                || !isMergeFwd( ev, merge, drxn ) )
        {
            addExtensionMerge( merge, ext, ev, doesBranch, drxn );
        }
    }
    
    if ( !merge.node && !ext.overlaps.empty() )
    {
        if ( doesBranch )
        {
            Node* node = new Node( getSeqEnd( ext.maxOverLen, drxn ), ext, ends_[drxn], drxn );
            addEdge( node, ext.maxOverLen, drxn );
            ev.nodes.push_back( node );
            bool fwdBranch = false;
            if ( !ext.fwdExts.empty() )
            {
                ev.ante.insert( this );
                node->addExtensions( ext.fwdExts, ev, fwdBranch, drxn );
                ev.ante.erase( this );
            }
            node->setCoverage();
        }
        else 
        {
            appendNode( ext, drxn );
            bool fwdBranch = false;
            if ( !ext.fwdExts.empty() )
            {
                addExtensions( ext.fwdExts, ev, fwdBranch, drxn );
            }
        }
        setCoverage();
    }
}

void Node::addExtension( Extension &ext, ExtVars &ev, NodeSet &acceptable, NodeSet &ignore, bool doesBranch, bool drxn )
{
    MergeHit merge;
    for ( Node* node : ignore )
    {
        if ( node->checkExtensionMerge( ext, merge ) ) return;
    }
    for ( Node* &node : ev.nodes )
    {
        if ( node->checkExtensionMerge( ext, merge ) ) break;
    }
    for ( Node* &node : ev.island )
    {
        if ( node->checkExtensionMerge( ext, merge ) ) break;
    }
    
    
    if ( merge.node )
    {
        if ( ignore.find( merge.node ) != ignore.end() )
        {
            if ( merge.overlap == ext.maxOverLen )
            {
                return;
            }
            else
            {
                merge.node = NULL;
            }
        }
        else if ( acceptable.find( merge.node ) != acceptable.end() 
                && !isMergeFwd( ev, merge, drxn )
                && ev.ante.find( merge.node ) == ev.ante.end() )
        {
            addExtensionMerge( merge, ext, ev, doesBranch, drxn );
        }
        else
        {
            ev.rebranch.insert( this );
        }
    }
    
    if ( !merge.node )
    {
        if ( doesBranch )
        {
            Node* node = new Node( getSeqEnd( ext.maxOverLen, drxn ), ext, ends_[drxn], drxn );
            addEdge( node, ext.maxOverLen, drxn );
            ev.nodes.push_back( node );
            if ( !ext.fwdExts.empty() )
            {
                ev.ante.insert( this );
                ignore.insert( node );
                for ( Extension &fwdExt : ext.fwdExts )
                {
                    node->addExtension( fwdExt, ev, acceptable, ignore, true, drxn );
                }
                ev.ante.erase( this );
            }
            node->setCoverage();
        }
        else 
        {
            appendNode( ext, drxn );
            if ( !ext.fwdExts.empty() )
            {
                for ( Extension &fwdExt : ext.fwdExts )
                {
                    addExtension( fwdExt, ev, acceptable, ignore, true, drxn );
                }
            }
        }
        setCoverage();
    }
}

void Node::addExtensionMerge( MergeHit &merge, Extension &ext, ExtVars &ev, bool doesBranch, bool drxn )
{
    Node* thisNode = this;
    
    if ( !ext.overlaps.empty() )
    {
        if ( doesBranch )
        {
            thisNode = new Node( getSeqEnd( ext.maxOverLen, drxn ), ext, ends_[drxn], drxn );
            addEdge( thisNode, ext.maxOverLen, drxn, ev.doOffset );
            ev.nodes.push_back( thisNode );
            thisNode->setCoverage();
        }
        else
        {
            appendNode( ext, drxn );
        }
    }
    
    bool isBypass = ev.doOffset && isMergeBypass( ev, merge, drxn );
    bool isCis = isBypass || isMergeCis( ev, merge, drxn );
    Node* mergeNode = merge.node->splitNode( (*merge.coords)[!drxn], ( merge.node->drxn_ <= 2 ? ev.nodes : ev.island ), drxn, drxn );
    
    if ( mergeNode->drxn_ <= 2 && !isCis )
    {
        mergeNode = new Node( mergeNode, ev, drxn );
        ev.cloneSet.insert( mergeNode );
        ev.nodes.push_back( mergeNode );
    }

    thisNode->addEdge( mergeNode, merge.overlap, drxn, ev.doOffset );
    mergeNode->setCoverage();
    
    if ( mergeNode->drxn_ > 2 )
    {
        mergeNode->mergeIsland( ev, drxn, true );
        NodeSet mergeSet = { mergeNode };
        Node::reviewMerged( ev, mergeSet, drxn );
    }
    else if ( mergeNode->validated_ )
    {
        thisNode->validateMerge( mergeNode, drxn );
    }
    
    assert( ev.offset.find( mergeNode ) == ev.offset.end() || mergeNode->edges_[!drxn].size() > 1 || !ev.doOffset );
}

void Node::addExtensions( vector<Extension> &exts, ExtVars &ev, bool doesBranch, bool drxn )
{
    int32_t endCoord = ends_[drxn];
    
    doesBranch = doesBranch || exts.size() > 1 || !edges_[drxn].empty() || clones_;
    vector<MergeHit> selfMerges;
    for ( Extension &ext : exts )
    {
        addExtension( ext, ev, selfMerges, doesBranch, drxn );
    }
    
    if ( !selfMerges.empty() && !edges_[drxn].empty() && ev.doLoop )
    {
        addSelfLoop( selfMerges, ev, drxn );
    }
    
    stop_[drxn] = ( edges_[drxn].empty() && ends_[drxn] == endCoord ) ? 1 : stop_[drxn];
}

void Node::addSelfLoop( vector<MergeHit> &selfMerges, ExtVars &ev, bool drxn )
{
    if ( !edges_[drxn].empty() && ev.cloneSet.empty() )
    {
        sort( selfMerges.begin(), selfMerges.end(), [&drxn]( const MergeHit &a, const MergeHit &b ){
            return ( drxn ? (*a.coords)[0] < (*b.coords)[0] : (*a.coords)[1] > (*b.coords)[1] );
        });

        Node* thisNode = splitNode( (*selfMerges[0].coords)[!drxn], ev.nodes, drxn, drxn );
        Node* clone = new Node( thisNode, ev, drxn );
        ev.cloneSet.insert( clone );
        ev.nodes.push_back( clone );
        thisNode->addEdge( clone, selfMerges[0].overlap, drxn );
    }
}

void Node::appendNode( Extension &ext, bool drxn )
{
    int32_t anchor = ends_[drxn];

    // Add and map reads in order of ascending extend length
    for ( vector<Overlap>::reverse_iterator iter = ext.overlaps.rbegin(); iter != ext.overlaps.rend(); ++iter )
    {
        addRead( *iter, anchor, drxn );
    }
    
    assert( !clones_ );
    if ( !paired_ )
    {
        paired_ = new NodeSet;
        validated_ = false;
    }
    validated_ = false;
    paired_->clear();
    
    this->appendSeq( ext.seq, drxn );
}

void Node::appendSeq( string &seq, bool drxn )
{
    seq_ = drxn ? seq_ + seq : seq + seq_;
    ends_[drxn] += drxn ? seq.length() : -seq.length();
}

bool Node::checkExtensionMerge( Extension &ext, MergeHit &merge )
{
    for ( auto it = ext.overlaps.rbegin(); it != ext.overlaps.rend(); it++ )
    {
        if ( !(*it).redundant && findRead( (*it).readId, merge.coords, false ) )
        {
            merge.node = this;
            merge.read = (*it).readId;
            merge.overlap = (*it).overLen;
            return ext.trimToOverlap( merge.overlap );
        }
    }
    return false;
}

void Node::debriefExtension( ExtVars &ev, bool drxn )
{
    // End extension if beyond pre-defined limits
    if ( abs( ends_[drxn] ) > abs( params.locusLimits[drxn] ) * 1.1 )
    {
        stop_[drxn] = 4;
    }
    
//    setCoverage( ev, drxn, drxn );
    setCoverage();
    
    if ( !isContinue( 0 ) && !isContinue( 1 ) )
    {
        if ( drxn_ == 2 && !isContinue( !drxn ) )
        {
            reliable_ = true;
        }
    }
}

void Node::extendComplete( ExtVars &ev, bool drxn )
{
    while( isContinue( drxn ) && abs( ends_[drxn] ) < abs( params.locusLimits[drxn] * 1.1 ) && !clones_ )
    {
        vector<Extension> exts = ev.bwt.mapExtensions( seq_, drxn );
        if ( exts.size() == 1 )
        {
            MergeHit merge;
            bool doContinue = exts[0].fwdExts.empty();
            exts[0].fwdExts.clear();
            for ( Node* &node : ev.nodes )
            {
                if ( node->checkExtensionMerge( exts[0], merge ) ) break;
            }
            if ( !exts[0].overlaps.empty() )
            {
                appendNode( exts[0], drxn );
                if ( doContinue && !merge.node )
                {
                    continue;
                }
            }
        }
        break;
    }
}

bool Node::extendForward( ExtVars &ev, int base, int &best, int target, bool drxn )
{
    NodeSet fwdSet = { this };
    NodeSet currSet = { this };
    
    NodeIntMap scores, reads, edges;
    
    while ( !currSet.empty() && ( best < 1 || ( best + base < target ) ) )
    {
        NodeSet currFwd, nxtSet;
        for ( Node* curr : currSet )
        {
            curr->getDrxnNodes( currFwd, drxn );
        }
        
        for ( Node* curr : currSet )
        {
            if ( currFwd.find( curr ) == currFwd.end() )
            {
                // Compile previous scores, read counts and edge counts
                int prevScore = 0, prevReads = 0, prevEdges = 0;
                for ( Node* bck : curr->getNextNodes( !drxn ) )
                {
                    auto it = scores.find( bck );
                    if ( it != scores.end() ) prevScore = max( prevScore, it->second );
                    
                    it = reads.find( bck );
                    if ( it != reads.end() ) prevReads = max( prevReads, it->second );
                    
                    it = edges.find( bck );
                    if ( it != edges.end() ) prevEdges = max( prevEdges, it->second );
                }
                
                // Extend if possible and valid
                if ( curr->isContinue( drxn ) 
                        && prevEdges < 3
                        && ev.cloneSet.find( curr ) == ev.cloneSet.end()
                        && ( prevReads + curr->reads_.size() ) * max( 1, prevEdges ) - ( prevScore * 8 ) < (int)params.cover )
                {
                    curr->extendForward( ev, fwdSet, drxn );
                }
                
                // Updates scores, read counts and edge counts
                int score = prevScore + curr->getPairHitsTmp( drxn );
                scores[curr] = score;
                reads[curr] = prevReads + curr->reads_.size();
                edges[curr] = prevEdges + curr->edges_[drxn].size() - 1;
                
                best = max( best, score );
                
                curr->getNextNodes( nxtSet, drxn );
                curr->getDrxnNodes( fwdSet, drxn );
            }
            else
            {
                nxtSet.insert( curr );
            }
        }
        
        currSet = nxtSet;
    }
    
    return best > 0 && best + base > 0;
}

void Node::extendForward( ExtVars &ev, NodeSet &acceptable, bool drxn )
{
    ev.ante = getDrxnNodes( !drxn, true, true );
    NodeSet ignore = getDrxnNodesInSet( acceptable, !drxn, true );
    int32_t endCoord = ends_[!drxn];
    while ( isContinue( drxn ) && endCoord != ends_[drxn] && ends_[1] - ends_[0] < params.readLen * 3 )
    {
        endCoord = ends_[drxn];
        vector<Extension> exts  = ev.bwt.mapExtensions( seq_, drxn );
        bool doesBranch = exts.size() > 1 || clones_;
        for ( Extension &ext : exts )
        {
            addExtension( ext, ev, acceptable, ignore, doesBranch, drxn );
        }
    }
    setCoverage();
}

void Node::extendNode( ExtVars &ev, bool drxn )
{
    setExtVars( ev, drxn );
    while( isContinue( drxn ) && extendCount_ > 0 && abs( ends_[drxn] ) < abs( params.locusLimits[drxn] * 1.1 ) )
    {
        bool noMatches, doesBranch = false;
        vector<Extension> extensions = ev.bwt.mapExtensions( noMatches, seq_, drxn );
        if ( noMatches && reads_.size() > 1 )
        {
            reEnd( ev, drxn );
        }
        else
        {
            addExtensions( extensions, ev, doesBranch, drxn );
            extendCount_--;
        }
    }
    
    debriefExtension( ev, drxn );
}

bool Node::extendOrigin( ExtVars &ev, bool drxn )
{
    int endMisses = getEndMarks( drxn );   
//    PairingVars pv; 
//    pv.tNodes = getTargetNodes( drxn, true );
//    
//    cout << "Extending origin" << endl;
//    while( isContinue( drxn ) && endMisses < 30 && abs( ends_[drxn] ) < abs( params.locusLimits[drxn] ) )
//    {
//        int markCount = readMarks_.size();
//        
//        // Extend
//        bool doesBranch = false;
//        vector<Extension> extensions = ev.bwt.mapExtensions( seq_, drxn );
//        addExtensions( extensions, ev, doesBranch, drxn );
//        
//        // Check for pairing with new marks
//        vector<ReadMark> newMarks = { readMarks_.begin() + markCount, readMarks_.end() };
//        pv.marks = &newMarks;
//        setPairs( pv, drxn );
//        
//        // Update pairs and misses
//        if ( !pv.pairs.empty() )
//        {
//            addPairs( pv.pairs, pv.hitIds, drxn );
//            pv.pairs.clear();
//            pv.hitIds.clear();
//            endMisses = 0;
//        }
//        else
//        {
//            endMisses += newMarks.size();
//        }
//    }
//    debriefExtension( ev, drxn );
    return endMisses >= 30;
}

void Node::extendSeed( ExtVars &ev, bool drxn )
{
    setExtVars( ev, drxn );
    while( isContinue( drxn ) )
    {
        vector<Extension> extensions = ev.bwt.mapExtensions( seq_, drxn );
        if ( extensions.size() != 1 )
        {
            pause( drxn );
            break;
        }
        bool doesBranch = false;
        addExtensions( extensions, ev, doesBranch, drxn );
    }
}

int32_t Node::getBestMergeClone( MergeHit &merge, int32_t fromEnd, bool drxn )
{
    int32_t bestOffset = ( drxn ? fromEnd - (*merge.coords)[0] : (*merge.coords)[1] - fromEnd ) - merge.overlap ;
    for ( Node* clone : this->getCloneSet( true ) )
    {
        auto it = reads_.find( merge.read );
        int32_t cloneOffset = ( drxn ? fromEnd - it->second[0] : it->second[1] - fromEnd ) - merge.overlap;
        if ( cloneOffset < bestOffset )
        {
            merge.node = clone;
            merge.coords = &it->second;
            bestOffset = cloneOffset;
        }
    }
    return bestOffset;
}

bool Node::isMergeBypass( ExtVars &ev, MergeHit &merge, bool drxn )
{
    if ( (*merge.coords)[!drxn] == merge.node->ends_[!drxn] )
    {
        for ( Node* clone : merge.node->getCloneSet( true ) )
        {
            if ( ev.ante.find( clone ) == ev.ante.end() )
            {
                auto it = clone->reads_.find( merge.read );
                for ( Node* prv : clone->getNextNodes( !drxn ) )
                {
                    if ( ev.ante.find( prv ) != ev.ante.end() && it != clone->reads_.end() )
                    {
                        merge.node = clone;
                        merge.coords = &it->second;
                        ev.bypass.insert( clone );
                        merge.node->clearPairsForward( drxn );
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

bool Node::isMergeCis( ExtVars &ev, MergeHit &merge, bool drxn )
{
    for ( Node* clone : merge.node->getCloneSet() )
    {
        if ( ev.ante.find( clone ) == ev.ante.end() )
        {
            auto it = clone->reads_.find( merge.read );
            assert( it != clone->reads_.end() );
            merge.node = clone;
            merge.coords = &it->second;
        }
    }
    
    if ( ev.ante.find( merge.node ) == ev.ante.end() && !merge.node->edges_[!drxn].empty() )
    {
        if ( abs( ends_[drxn] - ( drxn ? merge.overlap : -merge.overlap ) - (*merge.coords)[!drxn] ) >= 200 )
        {
            ev.offset.insert( merge.node );
            merge.node->clearPairsForward( drxn );
        }
        
        return true;
    }
    
    return false;
}

//void Node::isMergeDualLocus( ExtVars &ev, MergeHit &merge, bool drxn )
//{
//    if ( drxn_ <= 2 && !( drxn ? drxn_ == 0 : drxn_ == 1 ) )
//    {
//        bool isDual = drxn ? merge.node->drxn_ == 0 : merge.node->drxn_ == 1;
//        if ( merge.node->drxn_ == 2 )
//        {
//            isDual = isDual || drxn_ != 2;
//            for ( Node* nxt : merge.node->getNextNodes( drxn ) )
//            {
//                isDual = isDual || nxt->drxn_ == 2;
//            }
//        }
//        
//        if ( isDual )
//        {
//            return !( merge.node, (*merge.coords)[drxn], drxn );
//        }
//    }
//    
//    return false;
//}

bool Node::isMergeFwd( ExtVars &ev, MergeHit &merge, bool drxn )
{
    NodeSet fwdSet = getDrxnNodes( drxn );
    return merge.node->anyCloneInSet( fwdSet );
}

bool Node::isMergeSelf( ExtVars &ev, MergeHit &merge )
{
    if ( this == merge.node )
    {
        return true;
    }
    for ( Node* clone : merge.node->getCloneSet() )
    {
        if ( clone == this )
        {
            auto read = reads_.find( merge.read );
            merge.node = clone;
            merge.coords = &read->second;
            return true;
        }
    }
    return false;
}

Node* Node::mergeNode( NodeList &nodes, int32_t coord, bool subGraph )
{
    if ( coord != ends_[subGraph] )
    {
        int32_t splitCoord = ends_[subGraph];
        for ( auto &read : reads_ )
        {
            if ( subGraph ? coord <= read.second[1] && read.second[0] < splitCoord
                          : read.second[0] <= coord && splitCoord < read.second[1] )
            {
                splitCoord = read.second[!subGraph];
            }
        }
        
        splitNode( splitCoord, nodes, subGraph, subGraph );
    }
    return this;
}

Node* Node::mergeNode( NodeList &nodes, Coords* coords, bool subGraph, bool drxn )
{
    if ( (*coords)[!drxn] != ends_[!drxn] )
    {
//        int32_t splitBgn = findNextRead( (*coords)[drxn], !drxn );
        int32_t splitCoord = ends_[!drxn];
        for ( auto &read : reads_ )
        {
            if ( drxn ? ( read.second[0] < (*coords)[0] && splitCoord < read.second[1] )
                      : ( (*coords)[1] < read.second[1] && read.second[0] < splitCoord ) )
            {
                splitCoord = read.second[drxn];
            }
        }
        splitNode( splitCoord, nodes, subGraph, !drxn );
    }
    return this;
}

void Node::reEnd( ExtVars &ev, bool drxn )
{
    int32_t limit = ends_[drxn];
    for ( int i : { 0, 1 } )
    {
        int32_t coords[2] = { ends_[1], ends_[0] };
        for ( auto &read : reads_ )
        {
            if ( drxn ? read.second[1] < limit && coords[drxn] < read.second[1] 
                      : limit < read.second[0] && read.second[0] < coords[drxn] )
            {
                coords[drxn] = read.second[drxn];
            }
            if ( read.second[drxn] == limit )
            {
                coords[!drxn] = drxn ? min( coords[!drxn], read.second[0] ) : max( coords[!drxn], read.second[1] );
            }
        }
        
        if ( coords[drxn] == ends_[!drxn] || coords[!drxn] == ends_[!drxn] || coords[!drxn] == ends_[drxn] ) break;
        
        string seq = drxn ? seq_.substr( 0, coords[1] - ends_[0] ) : seq_.substr( coords[0] - ends_[0] );
        vector<Extension> exts = ev.bwt.mapExtensions( seq, drxn );
        removeReExtensions( exts, drxn, true );
        if ( !exts.empty() )
        {
            splitNode( coords[!drxn], ( drxn_ <= 2 ? ev.nodes : ev.island ), drxn, drxn );
            addExtensions( exts, ev, true, drxn );
        }
        limit = coords[drxn];
    }
    
    stop_[drxn] = edges_[drxn].empty() ? 1 : 0;
}

void Node::rebranchNode( ExtVars &ev, bool drxn )
{
    ev.ante = getDrxnNodes( !drxn, true, true );
    stop_[drxn] = 0;
    vector<Extension> exts = ev.bwt.mapExtensions( seq_, drxn );
    removeReExtensions( exts, drxn );
    addExtensions( exts, ev, true, drxn );
}

void Node::removeReExtensions( vector<Extension> &exts, bool drxn, bool inclSelf )
{
    NodeSet fwdSet = getDrxnNodes( drxn, false, inclSelf );
    for ( auto it = exts.begin(); it != exts.end(); )
    {
        MergeHit merge;
        bool didMerge = false;
        for ( Node* fwd : fwdSet )
        {
            if ( fwd->checkExtensionMerge( *it, merge ) )
            {
                it = exts.erase( it );
                didMerge = true;
                break;
            }
        }
        if ( !didMerge )
        {
            it++;
        }
    }
    
}

void Node::setExtVars( ExtVars &ev, bool drxn )
{
    extendCount_ = max( extendCount_, 1 );
    getDrxnNodes( ev.ante, !drxn, drxn_ != 2 );
}

Node* Node::splitNode( int32_t splitBegin, NodeList &nodes, bool subGraph, bool drxn, bool isOriginal )
{
    // No need to split if cutoff is at start of node
    if ( splitBegin == ends_[!drxn] )
    {
        return this;
    }
    
    // Create split node and trim this node
    clearPairs();
    
    Node* node = new Node( getSeqEnd( abs( ends_[drxn] - splitBegin ), drxn ), splitBegin, stop_, drxn );
    node->drxn_ = drxn_ == 2 ? subGraph : drxn_;
    nodes.push_back( node );
    int overlap = splitReads( node, splitBegin, drxn );
    assert( !reads_.empty() && !node->reads_.empty() );
    
    seq_ = drxn ? seq_.substr( 0, ends_[1] - ends_[0] ) : seq_.substr( seq_.length() - ( ends_[1] - ends_[0] ) );
    
    // Transfer forward edges to forward split node
    for ( Edge &edge : edges_[drxn] )
    {
        edge.node->removeEdge( this, !drxn );
        node->addEdge( edge.node, edge.overlap, drxn, edge.isLeap );
    }
    
    // Add edge from this to split node
    node->edgeCount_[drxn] = edgeCount_[drxn];
    edges_[drxn].clear();
    edgeCount_[drxn] = 0;
    addEdge( node, overlap, drxn, false );
    
    // Split each clone of this node - condition only called on first node
    if ( clones_ && isOriginal )
    {
        vector<NodePair> splitClonePairs;
        
        // Split each clone
        for ( Node* clone : *clones_ )
        {
            assert( clone->reads_.size() == this->reads_.size() + node->reads_.size() );
            int32_t cloneSplitBegin = splitBegin + ( clone->ends_[!drxn] - ends_[!drxn] );
            Node* splitClone = clone->splitNode( cloneSplitBegin, nodes, subGraph, drxn, false );
            node->addClone( splitClone );
            splitClonePairs.push_back( make_pair( clone, splitClone ) );
            
            assert( splitClone->reads_.size() == node->reads_.size() );
            assert( clone->reads_.size() == reads_.size() );
        }
        
        // Re-validate valid clones after all nodes are set
        for ( NodePair &np : splitClonePairs )
        {
            if ( np.first->validated_ )
            {
                np.first->setValid();
                np.second->setValid();
            }
        }
    }
    
    // Re-validate if valid
    if ( isOriginal && validated_ )
    {
        this->setValid();
        node->setValid();
    }
    
    this->setCoverage();
    node->setCoverage();
    
    return node;
}

Node* Node::splitNodeDual( int32_t* coords, NodeList &nodes, int subGraph )
{
    int32_t splitCoords[2] = { ends_[bool(subGraph)],ends_[bool(subGraph)] };
    for ( auto &read : reads_ )
    {
        if ( subGraph )
        {
            if ( read.second[0] < splitCoords[0] && coords[0] <= read.second[0] )
            {
                splitCoords[0] = read.second[0];
            }
            if ( read.second[0] < splitCoords[1] && coords[1] < read.second[1] )
            {
                splitCoords[1] = read.second[0];
            }
        }
        else
        {
            if ( splitCoords[0] < read.second[1] && read.second[0] < coords[0] )
            {
                splitCoords[0] = read.second[1];
            }
            if ( splitCoords[1] < read.second[1] && read.second[1] <= coords[1] )
            {
                splitCoords[1] = read.second[1];
            }
        }
    }
    if ( splitCoords[bool(subGraph)] != ends_[bool(subGraph)] )
    {
        splitNode( splitCoords[bool(subGraph)], nodes, subGraph, bool(subGraph) );
    }
    return splitNode( splitCoords[!bool(subGraph)], nodes, subGraph, bool(subGraph) );;
}
