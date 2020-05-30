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

void Node::addAlts( vector<QueryNode*> &alts, NodeRoll &nodes, bool drxn, int graph )
{
    for ( QueryNode* alt : alts )
    {
        assert( !alt->reads.empty() );
        MergeHit mergeBck, mergeFwd;
        if ( checkMergeRev( alt, nodes, drxn ) ) break;
        for ( Node* node : nodes.nodes ) if ( node->checkMergeRev( alt, mergeBck, drxn ) ) break;
        if ( alt->merged ) continue;
        alt->merged = alt->reads.size();
        assert( !alt->reads.empty() );
        for ( Node* node : nodes.nodes ) if ( node->checkMerge( alt, mergeFwd, drxn ) ) break;
        assert( !alt->reads.empty() );
        assert( mergeBck.node || alt->merged );
        
        for ( QueryNode* edge : alt->edges[1] ) edge->base = alt->ext.size();
        vector<Node*> forks, branches;
        if ( mergeBck.node ) assert( !mergeBck.coords->coords[2] );
        if ( mergeBck.node ) forks = mergeBck.node->splitNode( (*mergeBck.coords)[drxn], nodes, !drxn );
        if ( mergeFwd.node ) branches = mergeFwd.node->splitNode( (*mergeFwd.coords)[!drxn], nodes, drxn );
        
        // Abort if merges back into a continuing fork
        bool fail = false;
        for ( Node* fork : forks ) if ( fork->isContinue( drxn ) ) fail = true;
        if ( fail ) continue;
        
        // Simple edge between pre-existing nodes
        if ( !forks.empty() && !alt->merged )
        {
            bool branched = false;
            for ( Node* fork : forks ) for ( Node* branch : branches ) if ( fork->isEdge( branch, drxn, true ) ) branched = true;
            if ( !branched ) for ( Node* fork : forks ) for ( Node* branch : branches ) fork->addEdge( branch, mergeBck.ol, drxn, false );
        }
        // Create novel node
        else
        {
            alt->base = 0;
            Node* node = new Node( alt->seq, alt, ends_[drxn], drxn, graph );
            nodes += node;
            node->bad_ = true;
            // New node joins back
            for ( Node* fork : forks ) fork->addEdge( node, mergeBck.ol, drxn, false );
            // New node joins forward
            for ( Node* branch : branches ) node->addEdge( branch, mergeFwd.ol, drxn, false );
            // Continue extending
            if ( branches.empty() ) node->addExtensions( alt->edges[1], nodes, drxn, graph );
        }
//        if ( fork && !alt->merged )
//        {
//            Node* branch = mergeFwd.node->splitNode( (*mergeBck.coords)[!drxn], nodes, drxn );
//            if ( fork->isEdge( branch, drxn, true ) ) continue;
//            fork->addEdge( branch, mergeBck.ol, drxn, false );
//        }
//        // Extend pre-existing node
//        else if ( fork && fork->edges_[drxn].empty() ) fork->appendNode( alt, drxn );
//        // Create novel node
//        else node = new Node( alt->seq, alt, ends_[drxn], drxn, graph );
//        
//        if ( node )
//        {
//            node->bad_ = fork ? fork->bad_ : true;
//            nodes.add( node );
//            // New node joins back
//            if ( fork ) fork->addEdge( node, mergeBck.ol, drxn, false );
//            // New node joins forward
//            if ( mergeFwd.node )
//            {
//                Node* branch = mergeFwd.node->splitNode( (*mergeFwd.coords)[!drxn], nodes, drxn );
//                node->addEdge( branch, mergeFwd.ol, drxn, false );
//            }
//            // Continue extending
//            else node->addExtensions( alt->edges[1], nodes, drxn, graph );
//        }
//        else if ( fork && !mergeFwd.node ) fork->addExtensions( alt->edges[1], nodes, drxn, graph );
        
    }
    for ( QueryNode* alt : alts ) alt->resetBase();
}

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
            if ( merge.ol == ext.maxOverLen && ev.doLoop  )
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
            if ( merge.ol == ext.maxOverLen )
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

    thisNode->addEdge( mergeNode, merge.ol, drxn, ev.doOffset );
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
}

void Node::addExtensions( vector<Extension> &exts, ExtVars &ev, bool doesBranch, bool drxn )
{
    int32_t endCoord = ends_[drxn];
    
    doesBranch = doesBranch || exts.size() > 1 || !edges_[drxn].empty() || clones_ || dontExtend_;
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

void Node::addExtensions( vector<QueryNode*>& exts, NodeRoll &nodes, bool drxn, int graph )
{
    assert( !cloned_ );
    for ( QueryNode* ext : exts )
    {
        Node* branch = this;
        for ( QueryNode* node : ext->edges[1] ) node->base = ext->ext.size();
        
        // Check for merging
        MergeHit merge;
        ext->merged = ext->reads.size();
        for ( Node* node : nodes.nodes ) if ( node->checkMerge( ext, merge, drxn ) ) break;
        assert( merge.node || ext->merged == ext->reads.size() );
        
        // Extend this node or branch
        if ( exts.size() == 1 && edges_[drxn].empty() && ext->edges[0].size() < 2 && ext->merged ) appendNode( ext, drxn );
        else if ( ext->base < ext->ext.size() && ext->merged )
        {
            branch = new Node( getSeqEnd( ext->maxOl + ext->base, drxn ), ext, ends_[drxn], drxn, graph );
            branch->bad_ = bad_ || branch->drxn_ != drxn || drxn_ == !branch->drxn_;
            addEdge( branch, ext->maxOl + ext->base, drxn, false );
            nodes += branch;
        }
        
        // Perform merge
        if ( merge.node )
        {
            if ( merge.coords->coords[2] && ext->isIgnore( (int)ext->reads[ ext->merged ].ext - abs( merge.coords->coords[2] ), exts, drxn ) ) continue;
            vector<Node*> split = merge.node->splitNode( (*merge.coords)[!drxn], nodes, drxn );
            for ( Node* node : split ) branch->addEdge( node, merge.ol, drxn, false );
        }
        // Continue forward branches
        else if ( !ext->edges[1].empty() ) branch->addExtensions( ext->edges[1], nodes, drxn, graph );
    }
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
        thisNode->addEdge( clone, selfMerges[0].ol, drxn );
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

void Node::appendNode( QueryNode* ext, bool drxn )
{
    assert( ext->merged );
    for ( int i = 0; i < ext->merged; i++ )
    {
        if ( drxn ) add( ext->reads[i].id, ends_[1] - ext->reads[i].ol - ext->base
                                         , ends_[1] + ext->reads[i].ext - ext->base, ext->reads[i].redundant );
        else        add( ext->reads[i].id, ends_[0] - ext->reads[i].ext + ext->base
                                         , ends_[0] + ext->reads[i].ol + ext->base, ext->reads[i].redundant );
    }
    
    if ( pairedNodes_ ) delete pairedNodes_;
    pairedNodes_ = NULL;
    verified_ = false;
    
    string extSeq = ext->getSeq( drxn );
    appendSeq( extSeq, drxn );
    readTest();
}

void Node::appendSeq( string &seq, bool drxn )
{
    seq_ = drxn ? seq_ + seq : seq + seq_;
    ends_[drxn] += drxn ? seq.size() : -seq.size();
}

bool Node::checkExtensionMerge( Extension &ext, MergeHit &merge )
{
    for ( auto it = ext.overlaps.rbegin(); it != ext.overlaps.rend(); it++ )
    {
        if ( !(*it).redundant && findRead( (*it).readId, merge.coords, false ) )
        {
            merge.node = this;
            merge.id = (*it).readId;
            merge.ol = (*it).overLen;
            return ext.trimToOverlap( merge.ol );
        }
    }
    return false;
}

bool Node::checkMerge( QueryNode* ext, MergeHit &merge, bool drxn )
{
    Coords* coords;
    for ( int i = 0; i < ext->merged; i++ )
    {
        if ( ext->reads[i].redundant || !findRead( ext->reads[i].id, coords, false ) ) continue;
        if ( coords->coords[2] )
        {
            if ( drxn ? coords->coords[2] < 0 : 0 < coords->coords[2] ) assert( false );
            bool bad = !ext->edges[1].empty();
            for ( int j = i+1; j < ext->reads.size(); j++ ) if ( !ext->reads[j].redundant ) bad = true;
            assert( !bad );
        }
        merge.coords = coords;
        merge.node = this;
        merge.id = ext->reads[i].id;
        merge.ol = ext->base;
        for ( int j = 0; j < i; j++ ) merge.ol = max( merge.ol, ext->reads[j].ext );
        merge.ol += ext->reads[i].ol;
        ext->merged = i;
        return !i;
    }
    
    return false;
}

bool Node::checkMergeRev( QueryNode* ext, MergeHit &merge, bool drxn )
{
    int i = 0;
    for ( ; i < ext->reads.size(); i++ )
    {
        if ( ext->reads[i].redundant ) continue;
        if ( !findRead( ext->reads[i].id, merge.coords ) ) break;
        ext->base = ext->reads[i].ext;
        merge.id = ext->reads[i].id;
    }
    if ( !i ) return false;
    if ( i == ext->reads.size() ) ext->merged = ext->reads.size();
    if ( ext->merged ) return true;
    merge.node = this;
    merge.ol = ext->reads[i].ol + ext->base;
    assert( merge.ol < params.readLen );
    ext->reads.erase( ext->reads.begin(), ext->reads.begin() + i );
    ext->seq = drxn ? ext->seq.substr( ext->seq.size() - ext->reads[0].ol ) : ext->seq.substr( 0, ext->reads[0].ol );
    for ( Edge &e : edges_[drxn] ) if ( e.node->checkMergeRev( ext, merge, drxn ) ) return true;
    return false;
}

bool Node::checkMergeRev( QueryNode* ext, NodeRoll& nodes, bool drxn )
{
    while ( !ext->reads.empty() && ( ext->reads[0].ol + ext->reads[0].ext ) < params.readLen )
    {
        string seq = drxn ? ext->seq + ext->ext.substr( 0, ext->reads[0].ext ) 
                          : ext->ext.substr( ext->ext.size() - ext->reads[0].ext ) + ext->seq;
        for ( Node* node : nodes.nodes )
        {
            int it = node->seq_.find( seq );
            if ( it == string::npos ) continue;
            int32_t coords[2]{ node->ends_[0] + it, node->ends_[0] + it + (int)seq.size() };
            bool redundant = node->isRedundant( coords[0], coords[1] );
            assert( redundant || node->reads_.find( ext->reads[0].id ) != node->reads_.end() );
            node->add( ext->reads[0].id, coords[0], coords[1], redundant );
            if ( redundant ) ext->reads[0].redundant = 1;
        }
        if ( !ext->reads[0].redundant ) break;
        while ( !ext->reads.empty() && ext->reads[0].redundant ) ext->reads.erase( ext->reads.begin() );
        if ( ext->reads.empty() ) return true;
        ext->maxOl = ext->reads[0].ol;
        ext->seq = drxn ? ext->seq.substr( ext->seq.size() - ext->maxOl ) : ext->seq.substr( 0, ext->maxOl );
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

bool Node::extendable( Querier& bwt, bool drxn )
{
    if ( !isContinue( drxn ) ) return false;
    return !( stop_[drxn] = !bwt.isExtendable( seq_, drxn ) );
    assert( false );
    return true;
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

void Node::extendEdge( Querier& bwt, NodeRoll &nodes, int minOl, bool drxn )
{
    Nodes block( this, params.readLen, drxn, true ), mapped;
    int32_t coords[2];
    vector< pair<Node*, pair<int, int> > > ols;
    for ( Node* node : nodes.nodes ) if ( !block.find( node ) && mapSeqEnd( seq_, node->seq_, minOl, coords, drxn ) )
    {
        ols.push_back( make_pair( node, make_pair( coords[1] - coords[0], node->ends_[0] + coords[!drxn] ) ) );
        mapped += node ;
    }
    for ( int i = 0; i < ols.size(); i++ )
    {
        bool bad = false;
        for ( Edge& e : ols[i].first->edges_[!drxn] ) if ( mapped.find( e.node ) ) bad = true;
        if ( bad ) ols.erase( ols.begin() + i-- );
    }
    if ( ols.empty() ) return;
    sort( ols.begin(), ols.end(), []( pair<Node*, pair<int, int> >& a, pair<Node*, pair<int, int> >& b ){ return a.second.first > b.second.first; } );
    
    int coord = ols[0].second.second, ol = ols[0].second.first;
    ols[0].first->getNextReadCoord( coord, !drxn, drxn );
    vector<Node*> edges = ols[0].first->splitNode( coord, nodes, drxn );
    ol -= abs( coord - ols[0].second.second );
    addEdge( edges[0], ol, drxn, false );
    edges[0]->setVerified();
    setVerified();
    Node::verify();
    
//    return;
//    MapNode mn;
//    mn.seq = drxn ? getSeqEnd( params.readLen-1, 1 ) + edges[0]->seq_.substr( ol, params.readLen-ol-1 )
//                  : edges[0]->seq_.substr( edges[0]->size()+1-params.readLen, params.readLen-1-ol ) + getSeqEnd( params.readLen-1, 0 );
//    
//    bwt.mapSequence( mn.seq, mn.ids, mn.coords );
//    Node* edge[2]{ drxn ? this : edges[0], drxn ? edges[0] : this };
//    int oled[2]{ mapSeqOverlap( edge[0]->seq_, mn.seq, 50, 1 ), mapSeqOverlap( edge[1]->seq_, mn.seq, 50, 0 ) };
//    assert( oled[0] == params.readLen-1 && oled[1] == params.readLen-1 );
//    for ( int i = 0; i < mn.ids.size(); i++ ) if ( mn.coords[1][i] <= oled[0] || mn.seq.size() - oled[1] <= mn.coords[0][i] )
//    {
//        for ( int d : { 0, 1 } ) mn.coords[d].erase( mn.coords[d].begin() + i );
//        mn.ids.erase( mn.ids.begin() + i-- );
//    }
//    assert( false );
}

bool Node::extendFork( Querier &bwt, NodeRoll &nodes, int32_t dist, int branchLimit, bool drxn )
{
    if ( dist <= 0 ) return true;
    for ( int i = 0; i < 20; i++ )
    {
        vector< pair<Node*, int32_t> > exts = getExtendable( dist, drxn );
        if ( exts.empty() ) return true;
        if ( exts.size() > branchLimit && exts[0].second >= dist  ) return true;
        if ( exts.size() > branchLimit && exts[branchLimit].second < dist / 2 ) return false;
        if ( exts.size() > branchLimit * 2 && exts[branchLimit*2].second < dist ) return false;
        for ( int j = 0; j < exts.size() && ( j < branchLimit || exts[j].second < dist ); j++ ) exts[j].first->extendNode( bwt, nodes, drxn );
    }
    
    return false;
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
        if ( noMatches && reads_.size() > 1 ) reEnd( ev, drxn );
        else
        {
            addExtensions( extensions, ev, doesBranch, drxn );
            extendCount_--;
        }
    }
    
    debriefExtension( ev, drxn );
}

void Node::extendNode( Querier &bwt, NodeRoll &nodes, bool drxn )
{
    while ( isContinue( drxn ) )
    {
        QueryJunction qj = bwt.mapJunction( seq_, drxn );
        stop( qj.failure_, drxn );
        if ( qj.failure_ && reEnd( bwt, nodes, drxn ) ) break;
        addAlts( qj.alts_, nodes, drxn, drxn_ < 2 ? drxn_ : drxn );
        addExtensions( qj.nodes_, nodes, drxn, drxn_ < 2 ? drxn_ : drxn );
    }
    
    setCoverage();
}

void Node::extendSeed( ExtVars &ev, bool drxn )
{
    extendCount_ = 100;
    while( isContinue( drxn ) )
    {
        bool noMatches, doesBranch = false;
        vector<Extension> extensions = ev.bwt.mapExtensions( noMatches, seq_, drxn );
        if ( noMatches && reads_.size() > 1 ) reEnd( ev, drxn );
        else
        {
            addExtensions( extensions, ev, doesBranch, drxn );
            extendCount_--;
        }
    }
    
    setCoverage();
}

int32_t Node::getBestMergeClone( MergeHit &merge, int32_t fromEnd, bool drxn )
{
    int32_t bestOffset = ( drxn ? fromEnd - (*merge.coords)[0] : (*merge.coords)[1] - fromEnd ) - merge.ol ;
    for ( Node* clone : this->getCloneSet( true ) )
    {
        auto it = reads_.find( merge.id );
        int32_t cloneOffset = ( drxn ? fromEnd - it->second[0] : it->second[1] - fromEnd ) - merge.ol;
        if ( cloneOffset < bestOffset )
        {
            merge.node = clone;
            merge.coords = &it->second;
            bestOffset = cloneOffset;
        }
    }
    return bestOffset;
}

vector< pair<Node*, int32_t> > Node::getExtendable( int32_t dist, bool drxn )
{
    vector< pair<Node*, int32_t> > nodes;
    for ( auto& nd : NodeDists( this, dist, drxn, drxn, true ).map ) if ( nd.first->isContinue( drxn ) ) nodes.push_back( nd );
    sort( nodes.begin(), nodes.end(), []( pair<Node*, int32_t>& a, pair<Node*, int32_t>& b ){ return a.second < b.second; } );
    return nodes;
}

bool Node::isMergeBypass( ExtVars &ev, MergeHit &merge, bool drxn )
{
    if ( (*merge.coords)[!drxn] == merge.node->ends_[!drxn] )
    {
        for ( Node* clone : merge.node->getCloneSet( true ) )
        {
            if ( ev.ante.find( clone ) == ev.ante.end() )
            {
                auto it = clone->reads_.find( merge.id );
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
            auto it = clone->reads_.find( merge.id );
            assert( it != clone->reads_.end() );
            merge.node = clone;
            merge.coords = &it->second;
        }
    }
    
    if ( ev.ante.find( merge.node ) == ev.ante.end() )
    {
        if ( abs( ends_[drxn] - ( drxn ? merge.ol : -merge.ol ) - (*merge.coords)[!drxn] ) >= 200 )
        {
            ev.offset.insert( merge.node );
            merge.node->clearPairsForward( drxn );
        }
        
        return true;
    }
    
    return false;
}

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
            auto read = reads_.find( merge.id );
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

//void Node::pauseBad( bool graph )
//{
//    for ( Edge &e : edges_[!graph] ) if ( bad_ && !e.node->bad_ ) setNotBad( graph );
//    if ( edges_[!graph].empty() ) stop( BACK_END, !graph );
//    for ( Edge &e : edges_[!graph] ) if ( e.node->bad_ && e.node->drxn_ == graph ) e.node->pauseBad( graph );
//}

bool Node::reEnd( Querier &bwt, NodeRoll &nodes, bool drxn )
{
    // Don't both re-ending if there were too many possible branches
    if ( stop_[drxn] == MANY_END || bad_ ) return true;
    
    // Find three most terminal reads, split at the biggest gap
    int32_t ends[3]{ ends_[!drxn], ends_[!drxn], ends_[!drxn] };
    for ( auto &read : reads_ )
    {
        if ( read.second.redundant ) continue;
        for ( int i = 0; i < 3; i++ )
        {
            if ( drxn ? read.second[1] < ends[i] : ends[i] < read.second[0] ) continue;
            for ( int j = 3; --j > i; ) ends[j] = ends[j-1];
            ends[i] = read.second[drxn];
            break;
        }
    }
    for ( int i = 1; i < 3 && ends[i] != ends_[!drxn]; i++ )
    {
        if ( ends[i] == ends[i-1] ) continue;
        string seq = drxn ? seq_.substr( 0, ends[i] - ends_[0] ) : seq_.substr( ends[i] - ends_[0] );
        QueryJunction qj = bwt.mapJunction( seq, drxn );
        if ( qj.failure_ ) continue;
        
        int32_t coords[2];
        getSplitCoords( coords, ends[i], !drxn );
        splitNode( coords[drxn], nodes, drxn );
        addAlts( qj.alts_, nodes, drxn, drxn );
        addExtensions( qj.nodes_, nodes, drxn, drxn );
        stop( 0, drxn );
    }
    
    return !edges_[drxn].empty();
}

void Node::reEnd( ExtVars &ev, bool drxn )
{
    int32_t limit = ends_[drxn];
    for ( int i : { 0, 1 } )
    {
        int32_t coords[2] = { ends_[1], ends_[0] };
        if ( !edges_[drxn].empty() ) coords[!drxn] += ( drxn ? -getBestOverlap( 1 ) : getBestOverlap( 0 ) );
        
        for ( auto &read : reads_ )
        {
            
            if ( read.second[drxn] != limit ) continue;
            coords[!drxn] = drxn ? min( coords[0], read.second[0] ) 
                                 : max( coords[1], read.second[1] );
        }
        if ( coords[!drxn] == ends_[!drxn] || coords[!drxn] == ends_[drxn] ) break;
        for ( auto &read : reads_ )
        {
            if ( drxn ? coords[0] <= read.second[0] : read.second[1] <= coords[1] ) continue;
            coords[drxn] = drxn ? max( read.second[1], coords[1] ) : min( read.second[0], coords[0] );
        }
        if ( coords[drxn] == ends_[!drxn] || coords[drxn] == ends_[drxn] ) break;
        
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
    
    assert( abs( ends_[drxn] - splitBegin ) >= this->getBestOverlap( drxn ) );
    
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
        node->addEdge( edge.node, edge.ol, drxn, edge.leap );
    }
    
    // Add edge from this to split node
    node->edgeCount_[drxn] = edgeCount_[drxn];
    edges_[drxn].clear();
    edgeCount_[drxn] = 0;
    addEdge( node, overlap, drxn, false );
    
    // Split each clone of this node - condition only called on first node
    if ( clones_ && isOriginal )
    {
        vector< pair<Node*,Node*> > splitClonePairs;
        
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
        for ( auto &np : splitClonePairs )
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

vector<Node*> Node::splitNode( int32_t cut, NodeRoll& nodes, bool drxn )
{
    if ( cut == ends_[!drxn] ) return vector<Node*>{ this };
    
    vector<Node*> split[2]{ clones(), vector<Node*>() };
    
    for ( Node* node : split[0] ) split[1].push_back( node->splitNode2( cut + ( node->ends_[!drxn] - ends_[!drxn] ), nodes, drxn ) );
    
    if ( split[0][0] != split[1][0] ) for ( int i = 1; i < split[1].size(); i++ )
    {
        split[1][0]->addCloned( split[1][i] );
    }
    
    for ( int i : { 0, 1 } ) for ( Node* node : split[i] )
    {
        if ( node->verified_ ) node->setVerified();
        else node->ends_.reset( node->drxn_ );
        node->setCoverage();
        node->readTest();
    }
    
    return split[1];
}

Node* Node::splitNode2( int32_t cut, NodeRoll& nodes, bool drxn )
{
    // No need to split if cutoff is at start of node
    if ( cut == ends_[!drxn] ) return this;
    
    assert( abs( ends_[drxn] - cut ) >= getBestOverlap( drxn ) );
    
    // Create split node and trim this node
    clearPaired( false );
    
    Node* node = new Node( getSeqEnd( abs( ends_[drxn] - cut ), drxn ), cut, stop_, drxn );
    node->bad_ = bad_;
    node->branch_ = branch_;
    node->drxn_ = drxn_;
    node->mapped_ = mapped_;
    node->stop( stop_[drxn], drxn );
    stop( 0, drxn );
    nodes += node;
    ends_[drxn] = cut;
    for ( auto it = reads_.begin(); it != reads_.end(); )
    {
        if ( drxn ? cut <= it->second[0] : it->second[1] <= cut )
        {
            node->reads_.insert( *it );
            it = reads_.erase( it );
            continue;
        }
        ends_[drxn] = drxn ? max( ends_[drxn], it->second[1] ) : min( ends_[drxn], it->second[0] );
        it++;
    }
    remark();
    node->remark();
    seq_ = drxn ? seq_.substr( 0, ends_[1] - ends_[0] ) : seq_.substr( seq_.size() - ( ends_[1] - ends_[0] ) );
    node->ends_.inherit( ends_.limits );
    ends_.recoil();
    if ( drxn_ >= 2 ) ends_.splitOrigin( node->ends_, drxn_, node->drxn_ );
    
    // Transfer forward edges to forward split node and add edge between the split
    for ( Edge &e : edges( drxn ) )
    {
        e.node->removeEdge( this, !drxn, true );
        node->addEdge( e, drxn, true );
    }
    addEdge( node, abs( ends_[drxn] - cut ), drxn, false );
    
    return node;
}

//Node* Node::splitNode( int32_t cut, NodeRoll& nodes, bool drxn )
//{
//    // No need to split if cutoff is at start of node
//    if ( cut == ends_[!drxn] ) return this;
//    
//    assert( abs( ends_[drxn] - cut ) >= getBestOverlap( drxn ) );
//    
//    // Create split node and trim this node
//    clearPaired( false );
//    
//    Node* node = new Node( getSeqEnd( abs( ends_[drxn] - cut ), drxn ), cut, stop_, drxn );
//    node->bad_ = bad_;
//    node->branch_ = branch_;
//    node->drxn_ = drxn_;
//    node->mapped_ = mapped_;
//    nodes += node;
//    ends_[drxn] = cut;
//    for ( auto it = reads_.begin(); it != reads_.end(); )
//    {
//        if ( drxn ? cut <= it->second[0] : it->second[1] <= cut )
//        {
//            node->reads_.insert( *it );
//            it = reads_.erase( it );
//            continue;
//        }
//        ends_[drxn] = drxn ? max( ends_[drxn], it->second[1] ) : min( ends_[drxn], it->second[0] );
//        it++;
//    }
//    remark();
//    node->remark();
//    seq_ = drxn ? seq_.substr( 0, ends_[1] - ends_[0] ) : seq_.substr( seq_.size() - ( ends_[1] - ends_[0] ) );
//    node->ends_.inherit( ends_.limits );
//    ends_.recoil();
//    if ( drxn_ >= 2 ) ends_.splitOrigin( node->ends_, drxn_, node->drxn_ );
//    
//    // Transfer forward edges to forward split node
//    for ( Edge &edge : edges_[drxn] )
//    {
//        edge.node->removeEdge( this, !drxn );
//        node->addEdge( edge.node, edge.ol, drxn, edge.isLeap );
//    }
//    
//    // Add edge from this to split node
//    int ol = abs( ends_[drxn] - cut );
//    edges_[drxn].clear();
//    addEdge( node, ol, drxn, false );
//    
//    if ( cloned_ ) for ( Node* clone : cloned_->nodes )
//    {
//        int32_t off = clone->ends_[!drxn] - ends_[!drxn];
//        clone->clearPaired( false );
//        assert( clone->drxn_ < 2 );
//        Node* edge = new Node( node, nodes, clone->drxn_, clone->bad_ );
//        edge->offset( off );
//        nodes += edge;
//        clone->seq_ = seq_;
//        clone->ends_[drxn] = ends_[drxn] + off;
//        for ( auto it = clone->reads_.begin(); it != clone->reads_.end(); )
//        {
//            if ( clone->ends_[0] <= it->second[0] && it->second[1] <= clone->ends_[1] ) it++;
//            else it = clone->reads_.erase( it );
//        }
//
//        for ( Edge& e : clone->edges_[drxn] )
//        {
//            e.node->removeEdge( clone, !drxn );
//            edge->addEdge( e.node, e.ol, drxn, false, e.isLeap );
//        }
//        clone->clearEdges( drxn );
//        clone->addEdge( edge, ol, drxn, false );
//        clone->remark();
//        edge->remark();
//        clone->ends_.init( clone->ends_[!clone->drxn_] );
//        edge->ends_.init( edge->ends_[!edge->drxn_] );
//        clone->readTest();
//        edge->readTest();
//    }
//    
//    // Re-validate if valid
//    for ( Node* n : { this, node } )
//    {
//        if ( n->verified_ || n->canSetVerified() ) n->setVerified();
//        else if ( n->drxn_ < 2 ) n->ends_.init( n->ends_[!n->drxn_] );
//        n->setCoverage();
//        n->readTest();
//    }
//    
//    return node;
//}

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
