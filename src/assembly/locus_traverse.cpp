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

#include "locus_traverse.h"
#include "shared_functions.h"
#include "query_flay.h"
#include <algorithm>

Link::Link( Vertex* base, Vertex* edge, int ol, bool drxn )
: ol_( ol )
{
    node_[!drxn] = base;
    node_[drxn] = edge;
    for ( int d : { 0, 1 } ) coord_[d] = node_[d] ? node_[d]->coord_[!d]: 0;
    for ( int d : { 0, 1 } ) if ( node_[d] ) for ( Link* l : node_[d]->edges_[!d] )
    {
        assert( l->node_[!d] != node_[!d] );
    }
    for ( int d : { 0, 1 } ) if ( node_[d] ) node_[d]->edges_[!d].push_back( this );
}

bool VRead::isRedundant( VRead& redundant )
{
    if ( redundant.coord_[1] - redundant.coord_[0] >= coord_[1]-coord_[0] ) return false;
    if ( redundant.coord_[0] < coord_[0] || coord_[1] < redundant.coord_[1] ) return false;
    assert( false );
    return true;
}

VPairs::VPairs( Vertex* a, Vertex* b, int aCoord, int bCoord, bool bDrxn )
{
    node_[!bDrxn] = a;
    node_[bDrxn] = b;
    coord_[!bDrxn] = aCoord;
    coord_[bDrxn] = bCoord;
    int d = node_[1]->pairs_[0].size() > node_[0]->pairs_[1].size();
    VRead* vr;
    for ( const pair<ReadId, VRead>& read : node_[!d]->pairs_[d] ) if ( d ? read.second.coord_[1] <= coord_[0] : coord_[1] <= read.second.coord_[0] )
    {
        if ( !( vr = node_[d]->getPair( read.second.id_, !d ) ) ) continue;
        if ( d ? vr->coord_[0] < coord_[1] : coord_[0] < vr->coord_[1]  ) continue;
        pairs_.push_back( VRead( d ? read.second.id_ : read.first, d ? read.second.coord_[0] : vr->coord_[0], d ? vr->coord_[1] : read.second.coord_[1] ) );
    }
    if ( !pairs_.empty() ) for ( int i : { 0, 1 } ) node_[i]->paired_[!i].push_back( this );
}

int VPairs::set( Vertex* l, Vertex* r )
{
    VPairs* vp = new VPairs( l, r, l->coord_[1], r->coord_[0], 1 );
    if ( !vp->pairs_.empty() ) return vp->pairs_.size();
    delete vp;
    return 0;
}

Vertex::Vertex( Ext* ext, bool drxn )
{
    coord_[0] = coord_[1] = 0;
    while ( ext )
    {
        seq_ = seq_.empty() ? ext->seq_ : ( drxn ? seq_ + ext->ext_ : ext->ext_ + seq_ );
        coord_[drxn] = drxn ? seq_.size() : -seq_.size();
        int base = coord_[drxn] + ( drxn ? -ext->ext_.size() : ext->ext_.size() );
        
        sort( ext->exts_.begin(), ext->exts_.end(), []( Ext* a, Ext* b ){ return a->count_ > b->count_; } );
        for ( int i = 1; i < ext->exts_.size(); i++ ) new Link( this, new Vertex( ext->exts_[i], drxn ), ext->exts_[i]->reads_[0].ol_, drxn );
        for ( ExtRead& er : ext->reads_ ) reads_.push_back( VRead( er.id_, base - ( drxn ? er.ol_ : er.ext_ ), base + ( drxn ? er.ext_ : er.ol_ ) ) );
        for ( ExtRead& er : ext->redundant_ ) redundant_.push_back( VRead( er.id_, base - ( drxn ? er.ol_ : er.ext_ ), base + ( drxn ? er.ext_ : er.ol_ ) ) );
        ext = ext->exts_.empty() ? NULL : ext->exts_[0];
    }
    refill( !drxn );
}

Vertex::Vertex( string& seq, ReadId id, bool drxn )
: seq_( seq )
{
    coord_[0] = drxn ? 0 : -seq.size();
    coord_[1] = drxn ? seq.size() : 0;
}

void Vertex::addPair( ReadId id, VRead& vr, bool drxn )
{
    auto ins = pairs_[drxn].insert( make_pair( id, vr ) );
    if ( ins.second || ins.first->second.coord_[0] == vr.coord_[0] ) return;
    assert( false );
    int diff = vr.coord_[0] - ins.first->second.coord_[0];
    auto it = alts_[drxn].find( id );
    if ( it != alts_[drxn].end() )
    {
        if ( find( it->second.begin(), it->second.end(), diff ) == it->second.end() ) it->second.push_back( diff );
    }
    else alts_[drxn].insert( make_pair( id, vector<int>{ diff } ) );
}

void Vertex::addRead( VRead& vr )
{
    for ( int i = 0; i < reads_.size(); i++ )
    {
        if ( vr.isRedundant( reads_[i] ) )
        {
            redundant_.push_back( reads_[i] );
            reads_[i] = vr;
            assert( false );
            return;
        }
        if ( reads_[i].isRedundant( vr ) )
        {
            redundant_.push_back( vr );
            return;
        }
        if ( reads_[i].coord_[0] <= vr.coord_[0] ) continue;
        reads_.insert( reads_.begin() + i, vr );
        return;
    }
    reads_.push_back( vr );
}

void Vertex::discardPaired()
{
    for ( auto it = pairs_[1].begin(); it != pairs_[1].end(); )
    {
        auto found = pairs_[0].find( it->second.id_ );
        if ( found != pairs_[0].end() )
        {
            pairs_[0].erase( found );
            it = pairs_[1].erase( it );
        }
        else it++;
    }
    for ( int d : { 0, 1 } ) for ( VPairs* vp : paired_[d] ) for ( VRead& vr : vp->pairs_ )
    {
        ReadId id[2]{ vr.id_, params.getPairId( vr.id_ ) };
        for ( int i : { 0, 1 } ) assert( vp->node_[i]->pairs_[!i].erase( id[!i] ) );
    }
}

void Vertex::get( vector<Vertex*>& nodes, int coord, int limit, bool drxn )
{
    if ( find ( nodes.begin(), nodes.end(), this ) != nodes.end() ) assert( false );
    if ( find ( nodes.begin(), nodes.end(), this ) != nodes.end() ) return;
    if ( limit < 0 ) return;
    nodes.push_back( this );
    for ( Link* l : edges_[drxn] ) if ( drxn ? coord <= l->coord_[0] : l->coord_[1] <= coord )
    {
        l->node_[drxn]->get( nodes, l->coord_[drxn], limit + l->ol_ - abs( coord - l->coord_[!drxn] ), drxn );
    }
}

VRead* Vertex::getPair( ReadId id, bool drxn )
{
    auto it = pairs_[drxn].find( id );
    return it != pairs_[drxn].end() ? &it->second : NULL;
}

void Vertex::refill( int drxns )
{
    for ( int d : { 0, 1 } ) if ( drxns == 2 || d == drxns ) pairs_[d].clear();
    for ( int d : { 0, 1 } ) if ( drxns == 2 || d == drxns ) alts_[d].clear();
    Lib* lib;
    int d;
    for ( vector<VRead>* reads : { &reads_, &redundant_ } ) for ( VRead vr : *reads ) if ( ( lib = params.getLib( vr.id_ ) ) && lib->isPe )
    {
        ReadId id = vr.id_;
        if ( lib->getPair( id, d ) && ( drxns == 2 || d == drxns ) ) addPair( id, vr, d );
    }
}

//int Vertex::setPairs( Vertex* base, int coord, bool drxn, bool cont )
//{
//    int hits = VPairs::set( drxn ? base : this, drxn ? this : base, drxn ? coord : coord_[1], drxn ? coord_[0] : coord );
//    if ( cont ) for ( Link* edge : edges_[drxn] ) hits += edge->node_[drxn]->setPairs( base, coord, drxn, true );
//    return hits;
//}

VBranch::VBranch( Link* l, Vertex* base, bool drxn )
: node_( l->node_[drxn] ), end_( l->coord_[!drxn] + ( drxn ? -l->ol_ : l->ol_ ) ), flay_( l->coord_[!drxn] ), baseHits_( 0 ), altHits_( 0 )
{
    baseCoord_ = altCoord_ = end_;
    vector<Vertex*> nodes;
    node_->get( nodes, node_->coord_[!drxn], params.maxPeMean, drxn );
    for ( Vertex* node : nodes )
    {
        VPairs* vp = new VPairs( base, node, l->coord_[!drxn], node->coord_[!drxn], drxn );
        if ( vp->pairs_.empty() ) delete vp;
        else for ( VRead& vr : vp->pairs_ )
        {
            baseHits_++;
            baseCoord_ = drxn ? min( baseCoord_, vr.coord_[0] ) : max( baseCoord_, vr.coord_[1] );
        }
    }
}

TFork::TFork( string seq, bool drxn )
: ol_( 0 ), exted_( false )
{
    int len = min( (int)seq.size(), params.readLen );
    base_ = seq_ = ( drxn ? seq.substr( seq.size() - len ) : seq.substr( 0, len ) );
}

TFork::TFork( string& base, string seq, ReadId id, int ext, int ol, bool drxn )
: base_( base ), reads_{ make_pair( id, seq.size()-max( 0, ol )-ext ) }, ol_( ol ), exted_( true )
{
    ext_ = drxn ? seq.substr( ol+ext ) : seq.substr( 0, seq.size()-ol-ext );
    seq_ = drxn ? base_ + ext_ : ext_ + base_;
}

TFork::TFork( TFork* tf, string seq, ReadId id, int ext, int ol, bool drxn )
: ol_( ol ), exted_( true )
{
    ext_ = drxn ? tf->ext_.substr( ext ) : tf->ext_.substr( 0, tf->ext_.size()-ext );
    tf->ext_ = drxn ? tf->ext_.substr( 0, ext ) : tf->ext_.substr( tf->ext_.size()-ext );
    base_ = tf->seq_ = drxn ? tf->base_ + tf->ext_ : tf->ext_ + tf->base_;
    seq_ = drxn ? base_ + ext_ : ext_ + base_;
    
    for ( int i = 0; i < tf->reads_.size(); i++ ) if ( tf->reads_[i].second > ext )
    {
        reads_.push_back( make_pair( tf->reads_[i].first, tf->reads_[i].second - ext ) );
        tf->reads_.erase( tf->reads_.begin() + i-- );
    }
    
    exts_ = tf->exts_;
    tf->exts_ = { this, new TFork( tf->seq_, seq, id, ext, ol, drxn ) };
}

TFork::~TFork()
{
    for ( TFork* tf : exts_ ) delete tf;
}

bool TFork::add( string& seq, ReadId id, bool drxn )
{
    bool added = false;
    int32_t coords[2];
    if ( !mapSeqEnd( seq, seq_, seq.size()/3, coords, !drxn ) )
    {
        for ( TFork* tf : exts_ ) if ( tf->add( seq, id, drxn ) ) added = true;
        return added;
    }
    int ext = drxn ? coords[1] - base_.size() : (int)ext_.size() - coords[0];
    int ol = drxn ? (int)base_.size() - coords[0] : coords[1] - (int)ext_.size();
    if ( coords[1] - coords[0] == seq.size() ? ext <= 0 : ext < 0 ) return false;
    
    ol_ = max( ol, ol_ );
    
    // Matches some or all of ext_
    if ( add( seq, id, ol, drxn ) ) return true;
    assert( ol > 0 );
    
    // Split at base of ext_
    if ( added = !ext_.empty() ) new TFork( this, seq, id, 0, ol, drxn );
    // Already split at base and matching one of the branches
    else for ( TFork* tf : exts_ ) if ( tf->add( seq, id, ol, drxn ) ) added = true;
    
    // New basal branch
    if ( !added ) exts_.push_back( new TFork( base_, seq, id, 0, ol, drxn ) );
    
    return true;
}

bool TFork::add( string& seq, ReadId id, int ol, bool drxn )
{
    int ext = max( 0, - ol );
    while ( ext+ol < seq.size() && ext < ext_.size() && ( drxn ? seq[ol+ext] == ext_[ext] : seq.end()[-ol-ext-1] == ext_.end()[-ext-1] ) ) ext++;
    assert( ext+ol > 0 );
    if ( !ext && ( ext < ext_.size() ) ) return false;
    ol_ = max( ol, ol_ );
    
    for ( int i = 0; i < reads_.size(); i++ ) if ( id == reads_[i].first ) return true;
    
    // Matching entire ext_ and one of the existing branches
    if ( ext == ext_.size() ) for ( TFork* tf : exts_ ) if ( tf->add( seq, id, ext+ol, drxn ) ) return true;
    
    if ( ext+ol == seq.size() ) reads_.push_back( make_pair( id, ext ) );
    else if ( ext == ext_.size() && exts_.empty() )
    {
        ext_ = drxn ? ext_ + seq.substr( ol+ext ) : seq.substr( 0, seq.size()-ol-ext ) + ext_;
        seq_ = drxn ? base_ + ext_ : ext_ + base_;
        reads_.push_back( make_pair( id, ext_.size() ) );
        exted_ = true;
    }
    else if ( ext == ext_.size() ) exts_.push_back( new TFork( base_, seq, id, ext, ol, drxn ) );
    else new TFork( this, seq, id, ext, ol, drxn );
    
    return true;
}

bool TFork::update()
{
    bool updated = exted_;
    for ( TFork* tf : exts_ ) if ( tf->update() ) updated = true;
    exted_ = false;
    return updated;
}

Traverse::Traverse( Querier& bwt, Node* node, bool drxn )
: bwt_( &bwt ), path_{ node }, coords_{ make_pair( node, drxn ? -node->size() : 0 ) }, base_( new Vertex() ), len_( node->size() )
{
    base_->coord_[0] = base_->coord_[1] = 0;
    for ( vector<Edge> edges = node->edges_[!drxn]; edges.size() == 1 && ( node = edges[0].node ) && node->edges_[drxn].size() == 1; )
    {
        len_ = len_ + node->size() - edges[0].ol;
        path_.insert( drxn ? path_.begin() : path_.end(), node );
        coords_.insert( make_pair( node, drxn ? -len_ : len_ - node->size() ) );
        edges = edges[0].node->edges_[!drxn];
    }
    base_->coord_[!drxn] = drxn ? -len_ : len_;
    base_->seq_ = Node::getSeq( path_ );
    assert( len_ == base_->seq_.size() );
    for ( int d : { 0, 1 } ) graphs_[d] = NULL;
}

bool Traverse::trace( Querier& bwt, Node* node, bool drxn )
{
    Traverse trav( bwt, node, drxn );
    trav.setReads( 500, drxn );
    trav.setFlayed( bwt, 500, !drxn, drxn );
    trav.setFlayed( bwt, 500, drxn, drxn );
    trav.setPairs();
    assert( false );
    return false;
}

//void Traverse::print( TFork* fork, bool drxn )
//{
//    cout << ">Ext: " << fork->ext_.size() << "  |  reads: " << fork->reads_.size() << ( fork->reads_.empty() ? "" : "  | " );
//    vector<int> bases;
//    for ( pair<ReadId, int> read : fork->reads_ )
//    {
//        auto it = reads_.find( params.getPairId( read.first ) );
//        assert( it != reads_.end() );
//        bases.push_back( it->second );
//    }
//    sort( bases.begin(), bases.end() );
//    for ( int i = 0; i < bases.size(); i++ ) cout << " " << bases[i];
//    cout << endl << fork->seq_ << endl;
//    for ( TFork* tf : fork->exts_ ) print( tf, drxn );
//}

void Traverse::setFlayed( Querier& bwt, int maxLen, bool flayDrxn, bool endDrxn )
{
    maxLen = !maxLen ? len_ : min( len_, maxLen );
    string seq = endDrxn ? base_->seq_.substr( base_->seq_.size()-maxLen ) : base_->seq_.substr( 0, maxLen );
    QueryFlay qf( seq, bwt.ir_, 120, flayDrxn );
    vector<Exts*> exts = qf.flay( seq, bwt.qb_, flayDrxn );
    if ( flayDrxn != endDrxn && maxLen < len_ ) while ( flayDrxn ? exts.back()->coord_ == maxLen : !exts.back()->coord_ )
    {
        delete exts.back();
        exts.pop_back();
    }
    int cutoff = 2, ol = 120;
    if ( !flayDrxn ) reverse( exts.begin(), exts.end() );
    for ( Exts* es : exts ) for ( Ext* e : es->exts_ ) if ( e->count( false, flayDrxn ) > cutoff && e->reads_[0].ol_ > ol )
    {
        if ( flayDrxn ? es->coord_ == maxLen : !es->coord_ ) e->sanitise( ol );
        Link* edge = new Link( NULL, new Vertex( e, flayDrxn ), e->reads_[0].ol_, flayDrxn );
        edge->coord_[!flayDrxn] = es->coord_ + ( endDrxn ? base_->coord_[1] - maxLen : base_->coord_[0] );
        base_->edges_[flayDrxn].push_back( edge );
    }
    for ( Exts* es : exts ) delete es;
}

void Traverse::setGraph( Vertex* v, vector<Vertex*>& base, int coord, bool pairDrxn, bool drxn )
{
    assert( graphs_[pairDrxn] );
    for ( auto& read : v->pairs_[drxn] )
    {
        graphs_[pairDrxn]->add( read.first, bwt_->getSequence( read.first ) );
    }
}

void Traverse::setJumps( Vertex* v, vector<Vertex*>& base, int coord, bool drxn )
{
    auto it = find( base.begin(), base.end(), v );
    if ( it == base.end() ) return;
    base.erase( it );
    
    
    vector<Vertex*> tar;
    v->get( tar, coord, params.maxPeMean*2+params.readLen, drxn );
    if ( coord == v->coord_[drxn] ) tar.erase( remove( tar.begin(), tar.end(), v ), tar.end() );
    
    coord = drxn ? coord - params.maxPeMean-params.readLen : coord + params.maxPeMean;
    for ( auto& read : v->pairs_[drxn] ) if ( drxn ? coord < read.second.coord_[1] : read.second.coord_[0] < coord )
    {
        setMatches( tar, read.first, 35, drxn );
    }
    
    for ( Link* l : v->edges_[drxn] ) setJumps( l->node_[drxn], base, l->node_[drxn]->coord_[!drxn], drxn );
}

bool Traverse::setMatches( vector<Vertex*>& nodes, ReadId id, int minOl, bool drxn )
{
    auto it = reads_.find( id );
    string& seq = it != reads_.end() ? it->second : reads_.insert( make_pair( id, bwt_->getSequence( id ) ) ).first->second;
    assert( !seq.empty() && seq.size() <= params.readLen );
    
    struct Coord{ int coords[2]; };
    Coord c;
    VRead vr( id, 0, 0 );
    Vertex* hit = NULL;
    vector<pair<Vertex*, Coord> > hits;
    for ( Vertex* v : nodes ) if ( mapSeqEnd( seq, v->seq_, minOl, c.coords, !drxn ) )
    {
        for ( int d : { 0, 1 } ) c.coords[d] += v->coord_[0];
        if ( drxn ? c.coords[1] < v->reads_[0].coord_[1] : v->reads_.back().coord_[0] < c.coords[0] ) continue;
        assert( !hit );
        int len = c.coords[1]-c.coords[0];
        bool add = true;
        if ( len == seq.size() )
        {
            hit = v;
            for ( int d : { 0, 1 } ) vr.coord_[d] = c.coords[d];
            hit->addRead( vr );
            continue;
        }
        for ( Link* l : v->edges_[drxn] ) if ( l->coord_[!drxn] == c.coords[drxn] )
        {
            int ext = 0, limit = min( seq.size() - len, l->node_[drxn]->seq_.size()-l->ol_ );
            while ( ext < limit && ( drxn ? seq[len+ext] == l->node_[1]->seq_[l->ol_+ext] : seq.end()[-len-ext-1] == l->node_[0]->seq_.end()[-l->ol_-ext-1] ) ) ext++;
            if ( ext < limit ) continue;
            if ( drxn ? c.coords[0] < l->coord_[0]-l->ol_ : l->coord_[1]+l->ol_ < c.coords[1] )
            {
                // NYI bridge
                assert( false );
            }
            else add = false;
        }
        if ( add ) hits.push_back( make_pair( v, c ) );
    }
    if ( hits.empty() ) return false;
    sort( hits.begin(), hits.end(), []( pair<Vertex*, Coord>& a, pair<Vertex*, Coord>& b){ return a.second.coords[1]-a.second.coords[0] > b.second.coords[1]-b.second.coords[0]; } );
    
    if ( !hit ) for ( auto& h : hits ) if ( h.second.coords[drxn] == h.first->coord_[drxn] )
    {
        hit = h.first;
        int ext = seq.size() - ( h.second.coords[1] - h.second.coords[0] );
        hit->seq_ = drxn ? hit->seq_ + seq.substr( seq.size()-ext ): seq.substr( 0, ext ) + hit->seq_;
        hit->coord_[drxn] += drxn ? ext : -ext;
        vr = VRead( id, drxn ? h.second.coords[0] : hit->coord_[0], drxn ? hit->coord_[1] : h.second.coords[1] );
        hit->reads_.insert( drxn ? hit->reads_.end() : hit->reads_.begin(), vr );
        break;
    }
    
    if ( !hit )
    {
        hit = new Vertex( seq, id, drxn );
        vr = VRead( id, hit->coord_[0], hit->coord_[1] );
        hit->reads_.push_back( vr );
        nodes.push_back( hit );
    }
    
    for ( auto& h : hits ) if ( h.first != hit )
    {
        Link* l = new Link( hit, h.first, h.second.coords[1]-h.second.coords[0], !drxn );
        l->coord_[!drxn] = h.second.coords[drxn];
    }
    
    if ( hits.size() > 1 )
    {
        int x = 0;
    }
    
    return true;
}

void Traverse::setPairs()
{
    int edges[2]{0};
    vector<VBranch*> branches[2], merged[2], unmerged[2];
    for ( int d : { 0, 1 } ) for ( Link* l : base_->edges_[d] ) branches[d].push_back( new VBranch( l, base_, d ) );
    for ( VBranch* l : branches[0] )
    {
        vector<Vertex*> nodes[2];
        l->node_->get( nodes[0], l->node_->coord_[1], params.maxPeMean, 0 );
        for ( VBranch* r : branches[1] ) if ( l->end_ <= r->flay_ || l->flay_ <= r->end_ )
        {
            r->node_->get( nodes[1], r->node_->coord_[0], params.maxPeMean, 1 );
            int hits = 0;
            for ( Vertex* ll : nodes[0] ) for ( Vertex* rr : nodes[1] ) hits += VPairs::set( ll, rr );
            if ( hits ) l->altCoord_ = max( l->altCoord_, r->flay_ );
            if ( hits ) r->altCoord_ = min( r->altCoord_, l->flay_ );
            if ( l->baseCoord_ <= r->flay_ ) l->altHits_ += hits;
            if ( l->flay_ <= r->baseCoord_ ) r->altHits_ += hits;
            nodes[1].clear();
        }
    }
    for ( int d : { 0, 1 } ) for ( VBranch* vb : branches[d] ) if ( vb->baseHits_ > 2 && !vb->altHits_ ) merged[d].push_back( vb );
    tips_ = { base_ };
    base_->discardPaired();
    // set blocked
    
    graphs_[1] = new KmerGraph( 32 );
    vector<Vertex*> base;
    setGraph( base_, base, 0, 1, 1 );
    graphs_[1]->assemble( base_->seq_, 1 );
    assert( false );
    int drxn = 1;
    for ( Vertex* v : tips_ )
    {
        vector<Vertex*> base;
        v->get( base, base_->coord_[drxn], params.maxPeMean*2+params.readLen, drxn );
        setJumps( v, base, v->coord_[drxn], drxn );
        assert( false );
    }
    for ( VBranch* vb : merged[1] );
    assert( false );
}

void Traverse::setReads( int len, int endDrxn )
{
    int coords[2]{ endDrxn == 1 ? base_->coord_[1]-len : base_->coord_[0], endDrxn == 0 ? base_->coord_[0]+len : base_->coord_[1] };
    for ( const pair<Node*, int32_t>& nc : coords_ ) if ( ( coords[0] <= nc.second+nc.first->size()-params.readLen ) && ( nc.second+params.readLen <= coords[1] ) )
    {
        Lib* lib;
        int off = nc.first->ends_[0] - nc.second, d;
        for ( auto& read : nc.first->reads_ ) if ( ( coords[0] <= read.second[0]-off ) && ( read.second[1]-off <= coords[1] ) )
        {
            ReadId id = read.first;
            if ( !( lib = params.getLib( id ) ) || !lib->isPe || !lib->getPair( id, d ) ) continue;
            VRead vr( read.first, read.second[0]-off, read.second[1]-off );
            base_->addPair( id, vr, d );
        }
    }
}

//void Traverse::setPairs( Querier& bwt, bool drxn )
//{
//    unordered_map<ReadId, int> reads;
//    for ( pair<Node*, int32_t> nc : coords_ ) if ( ( drxn ? -nc.second - nc.first->size() : nc.second ) < params.maxPeMean ) for ( auto& read : nc.first->reads_ )
//    {
//        int coord = nc.second + read.second[!drxn] - nc.first->ends_[0];
//        auto it = reads.insert( make_pair( read.first, coord ) );
//        if ( !it.second && abs( coord < it.first->second ) ) it.first->second = coord;
//    }
//    
//    Lib* lib;
//    int counted = 0;
//    for ( auto& read : reads ) if ( ( lib = params.getLib( read.first ) ) && lib->isPe )
//    {
//        ReadId id = read.first;
//        if ( ( lib->getPair( id ) != drxn ) || ( reads.find( id ) != reads.end() ) ) continue;
//        if ( ( drxn ? read.second + lib->size : lib->size - read.second ) < -300 ) continue;
//        reads_.insert( make_pair( read.first, read.second ) );
//        pairs_.insert( make_pair( id, bwt.getSequence( id ) ) );
//    }
//}
//
//void Traverse::setTraces( bool drxn )
//{
//    assert( forks_.empty() );
//    forks_.push_back( new TFork( Node::getSeq( path_ ), drxn ) );
//    vector<TFork*> forks = forks_;
//    for ( int again = 1; again-- > 0; )
//    {
//        for ( auto& read : pairs_ ) for ( TFork* tf : forks ) tf->add( read.second, read.first, drxn );
//        for ( TFork* tf : forks ) if ( tf->update() ) again = 1;
//    }
//    for ( TFork* tf : forks_ ) print( tf, drxn );
//    assert( false );
//}
