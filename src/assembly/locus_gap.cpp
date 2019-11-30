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

#include "locus_gap.h"
#include "shared_functions.h"
#include <algorithm>

void BudPair::align( int32_t coords[2] )
{
    if ( ( coords[1] - coords[0] ) > ( best_[1] - best_[0] ) ) for ( int d : { 0 , 1 } ) best_[d] = coords[d];
}

bool BudPair::alignEnds( BudPair* bp, int i, int j, bool drxn )
{
    while ( i-- > 0 && j-- > 0 ) if ( drxn ? seq_.end()[-i-1] != bp->seq_.end()[-j-1] : seq_[i] != bp->seq_[j] )
    {
        return false;
    }
    return true;
}

void Kmers::add( BudPair* bp )
{
    Kmer k;
    k.bp_ = bp;
    k.off_ = 0;
    for ( ; k.off_+31 < bp->seq_.size(); k.off_++ )
    {
        string q = bp->seq_.substr( k.off_, 32 );
        auto ins = kmers_.insert( make_pair( q, vector<Kmer>{ k } ) );
        if ( !ins.second ) ins.first->second.push_back( k );
    }
}

vector<Kmer>* Kmers::get( string& q )
{
    auto it = kmers_.find( q );
    return it != kmers_.end() ? &it->second : NULL;
}

BudExt::BudExt( string& seq, ReadId id, int ol, bool drxn )
: ext_( drxn ? seq.substr( ol ): seq.substr( 0, seq.size() - ol ) ), id_( id ), ol_( ol ), align_( 0 )
{}

BudMapx::BudMapx( BudPair* bp, int32_t q[2], int32_t t[2], bool drxn )
: pair_( bp ), ext_( drxn ? bp->seq_.size() - t[1] : t[0] )
{
    for ( int d : { 0, 1 } ) q_[d] = q[d];
    for ( int d : { 0, 1 } ) t_[d] = t[d];
}

//bool BudBranch::addExt( BudExt* ext, bool drxn )
//{
//    for ( int i = 0; i < ext->ext_.size(); i++ ) if ( drxn ? ext->ext_[i] != ext_->ext_[i] : ext->ext_.end()[-i-1] != ext_->ext_.end()[-i-1] ) return false;
//    exts_.push_back( ext );
//    return true;
//}
//
//void BudBranch::addMap( BudPair* bp, ReadId id, int32_t q[2], int32_t t[2], bool drxn )
//{
//    int i = 0;
//    for ( ; i < maps_.size() && ( q[!drxn] == maps_[i]->q_[!drxn] ? ( drxn ? q[1] <= maps_[i]->q_[1] : maps_[i]->q_[0] <= q[0] ) 
//                                                                  : ( drxn ? maps_[i]->q_[0] < q[0] : q[1] < maps_[i]->q_[1] ) ); i++ )
//    {
//        if ( q[0] == maps_[i]->q_[0] && q[1] == maps_[i]->q_[1] && id == maps_[i]->id_ ) assert( false );
//        if ( q[0] == maps_[i]->q_[0] && q[1] == maps_[i]->q_[1] && id == maps_[i]->id_ ) return;
//    }
//    maps_.insert( maps_.begin() + i, new BudMap( bp, id, q, t ) );
//}
//
//void BudBranch::setGood( bool drxn )
//{
//    int ext = 0, err = 0;
//    for ( BudMap* bm : maps_ )
//    {
//        int32_t coord[2]{ drxn ? bm->q_[0] - coord_ : coord_ - bm->q_[1], drxn ? bm->q_[1] - coord_ : coord_ - bm->q_[0] };
//    }
//}

BudAlign::~BudAlign()
{
    for ( BudAlign* ba : edges_ ) delete ba;
    for ( BudExt* be : exts_ ) delete be;
    for ( BudMapx* bm : maps_ ) delete bm;
}

bool BudAlign::addExt( BudExt* ext, bool drxn )
{
    if ( ( ext->align_ = align( ext->ext_, ext_, base_, drxn ) ) == base_ && base_ ) return false;
    
    if ( ext->align_ == ext->ext_.size() )
    {
        int i = 0;
        for ( ; i < exts_.size() && exts_[i]->align_ >= ext->align_; i++ );
        exts_.insert( exts_.begin() + i, ext );
        return true;
    }
    
    if ( ext->align_ == ext_.size() && !edges_.empty() )
    {
        for ( BudAlign* ba : edges_ ) if ( ba->addExt( ext, drxn ) ) return true;
        exts_.insert( exts_.begin(), ext );
        return true;
    }
    
    int i = 0;
    for ( ; i < exts_.size() && exts_[i]->align_ >= ext->align_; i++ ) if ( exts_[i]->align_ == ext->align_ )
    {
        if ( ( ext->align_ = align( exts_[i]->ext_, ext->ext_, ext->align_, drxn ) ) == exts_[i]->align_ ) continue;
        if ( exts_[i]->align_ == len_ )
        {
            len_ = exts_[i]->align_ = ext->align_;
            ext_ = drxn ? ext->ext_.substr( 0, len_ ) : ext->ext_.substr( ext->ext_.size() - len_ );
            if ( i ) swap( exts_[i], exts_[0] );
            exts_.insert( exts_.begin(), ext );
            return true;
        }
        vector<BudAlign*> edges{ new BudAlign( coord_, len_, exts_[i]->align_ ), new BudAlign( coord_, ext->align_, exts_[i]->align_ ) };
        len_ = exts_[i]->align_;
        exts_[i]->align_ = ext->align_;
        edges[0]->ext_ = ext_;
        edges[1]->ext_ = drxn ? ext->ext_.substr( 0, ext->align_ ) : ext->ext_.substr( ext->ext_.size() - ext->align_ );
        ext_ = drxn ? ext_.substr( 0, len_ ) : ext_.substr( ext_.size() - len_ );
        edges[0]->edges_ = edges_;
        edges[1]->exts_.push_back( exts_[i] );
        edges[1]->exts_.push_back( ext );
        exts_.erase( exts_.begin() + i );
        int j = 0;
        for ( ; j < exts_.size() && exts_[j]->align_ > len_; j++ );
        edges[0]->exts_.insert( edges[0]->exts_.begin(), exts_.begin(), exts_.begin() + j );
        exts_.erase( exts_.begin(), exts_.begin() + j );
        edges_ = edges;
        return true;
    }
    
    exts_.insert( exts_.begin()+i, ext );
    
    return true;
}

void BudAlign::addMap( BudPair* bp, int32_t q[2], int32_t t[2], bool drxn )
{
    bp->align( t );
    int i = 0;
    for ( ; i < maps_.size() && ( q[!drxn] == maps_[i]->q_[!drxn] ? ( drxn ? q[1] <= maps_[i]->q_[1] : maps_[i]->q_[0] <= q[0] ) 
                                                                  : ( drxn ? maps_[i]->q_[0] < q[0] : q[1] < maps_[i]->q_[1] ) ); i++ )
    {
        if ( q[0] == maps_[i]->q_[0] && q[1] == maps_[i]->q_[1] && bp->id_ == maps_[i]->pair_->id_ ) assert( false );
        if ( q[0] == maps_[i]->q_[0] && q[1] == maps_[i]->q_[1] && bp->id_ == maps_[i]->pair_->id_ ) return;
    }
    maps_.insert( maps_.begin() + i, new BudMapx( bp, q, t, drxn ) );
}

bool BudAlign::advance( int base, bool drxn )
{
    if ( !base_ )
    {
        if ( edges_.empty() || ext_.empty() ) return false;
        ext_.clear();
        base = len_;
        ol_ = good_ = len_ = 0;
        for ( BudAlign* ba : edges_ ) ol_ = max( ol_, ba->ol_ + base );
        for ( BudExt* be : exts_ ) delete be;
        exts_.clear();
    }
    else
    {
        for ( BudExt* be : exts_ )
        {
            be->ext_ = drxn ? be->ext_.substr( base ) : be->ext_.substr( 0, be->ext_.size() - base );
            be->align_ -= base;
            be->ol_ += base;
        }
        ext_ = drxn ? ext_.substr( base ) : ext_.substr( 0, ext_.size() - base );
        ol_ += base;
        len_ -= base;
        base_ -= base;
        good_ = 0;
    }
    
    for ( BudMapx* bm : maps_ ) delete bm;
    maps_.clear();
    coord_ += drxn ? base : -base;
    
    for ( BudAlign* ba : edges_ ) ba->advance( base, drxn );
    return true;
}

int BudAlign::align( string& a, string& b, int i, bool drxn )
{
    for ( ; i < min( a.size(), b.size() ) && ( drxn ? a[i] == b[i] : a.end()[-i-1] == b.end()[-i-1] ); i++ );
    return i;
}

void BudAlign::fold( BudAlign* ba )
{
    for ( BudAlign* edge : ba->edges_ ) fold( edge );
    ba->edges_.clear();
    for ( BudExt* be : ba->exts_ ) be->align_ = len_;
    ba->exts_.clear();
    delete ba;
}

bool BudAlign::isBranch( int limit )
{
    limit -= exts_.size();
    if ( limit <= 0 ) return true;
    for ( BudAlign* ba : edges_ ) if ( ba->isBranch( limit ) ) return true;
    return false;
}

void BudAlign::print( int32_t base, bool drxn )
{
    for ( BudAlign* ba : edges_ ) ba->print( base, drxn );
    for ( BudExt* be : exts_ )
    {
        cout << ">Ext " << be->id_ << ": " << coord_ + ( drxn ? (int)be->ext_.size() : -(int)be->ext_.size() ) << endl;
        cout << string( coord_ - base - ( drxn ? 0 : (int)be->ext_.size() ), '-' ) << be->ext_ << endl;
    }
}

void BudAlign::setExt( bool drxn )
{
    if ( good_ < len_ ) return;
    for ( BudAlign* ba : edges_ ) ba->setExt( drxn );
    if ( !edges_.empty() ) return;
    
    vector<BudMapx*> exts;
    for ( BudMapx* bm : maps_ ) if ( drxn ? bm->q_[1] >= coord_ + len_ : bm->q_[0] <= coord_ - len_ ) exts.push_back( bm );
    assert( !exts.empty() );
    sort( exts.begin(), exts.end(), [&]( BudMapx* a, BudMapx* b ){ return a->ext_ > b->ext_; } );
    vector<int> lens;
    int cut = exts[0]->ext_;
    for ( int i = 1; i < exts.size(); i++ )
    {
        int j = drxn ? exts[0]->t_[1] : exts[0]->t_[0]-1, k = drxn ? exts[i]->t_[1] : exts[i]->t_[0]-1, len = 0;
        for ( ; len < exts[i]->ext_ && exts[0]->pair_->seq_[j+(drxn?len:-len)] == exts[i]->pair_->seq_[k+(drxn?len:-len)]; len++ );
        if ( len != exts[i]->ext_ ) cut = min( cut, len );
        lens.push_back( len );
    }
    if ( !cut ) return;
    ext_ = drxn ? ext_ + exts[0]->pair_->seq_.substr( exts[0]->t_[1], cut ) : exts[0]->pair_->seq_.substr( exts[0]->t_[0]-cut, cut ) + ext_;
    for ( BudMapx* bm : exts )
    {
        bm->q_[drxn] = drxn ? bm->q_[1] + min( cut, bm->ext_ ) : bm->q_[0] - min( cut, bm->ext_ );
        bm->t_[drxn] = drxn ? bm->t_[1] + min( cut, bm->ext_ ) : bm->t_[0] - min( cut, bm->ext_ );
    }
    for ( BudExt* be : exts_ ) if ( be->align_ == len_ )
    {
        for ( int& i = be->align_; i < be->ext_.size() && ( drxn ? ext_[i] == be->ext_[i] : ext_.end()[-i-1] == be->ext_.end()[-i-1] ); i++ );
        if ( be->align_ > len_ ) lens.push_back( be->align_ - len_ );
    }
    len_ = ext_.size();
    sort( lens.rbegin(), lens.rend() );
    if ( lens.size() > 1 ) good_ += lens[1];
}

void BudAlign::setForks()
{
    forked_ = edges_.size() > 1;
    if ( edges_.size() > 1 && !ext_.empty() ) for ( BudAlign* ba : edges_ ) if ( !ba->isBranch( 5 ) ) forked_ = false;
    
    for ( BudAlign* ba : edges_ ) ba->setForks();
}

bool BudAlign::setGood( int32_t good, bool drxn )
{
    int32_t exted = coord_ + ( drxn ? good : -good );
    for ( int i = 0; i < maps_.size(); i++ )
    {
        int ext = abs( maps_[i]->q_[drxn] - coord_ );
        int ol = drxn ? exted - maps_[i]->q_[0] : maps_[i]->q_[1] - exted;
        bool bad = drxn ? maps_[i]->t_[1] < maps_[i]->pair_->best_[1] : maps_[i]->pair_->best_[0] < maps_[i]->t_[0];
        if ( bad )
        {
            delete maps_[i];
            maps_.erase( maps_.begin() + i-- );
        }
        else if ( ol > 0 && ext > good_ )
        {
            good_ = ext;
            exted = maps_[i]->q_[drxn];
        }
    }
    bool edged = false;
    for ( BudAlign* ba : edges_ ) if ( ba->setGood( good_, drxn ) ) edged = true;
    if ( edged ) good_ = len_;
    
    if ( edged ) for ( int i = 0; i < edges_.size(); i++ ) if ( !edges_[i]->good_ )
    {
        fold( edges_[i] );
        edges_.erase( edges_.begin() + i-- );
    }
    
    if ( edges_.size() > 1 ) for ( int i = 0; i < edges_.size(); i++ ) if ( edges_[i]->ol_ == 51 && edges_[i]->exts_[0]->id_ == 1094596642 )
    {
        delete edges_[i];
        edges_.erase( edges_.begin() + i-- );
    }
    if ( edges_.size() > 1 ) for ( int i = 0; i < edges_.size(); i++ ) if ( edges_[i]->ol_ == 42 && edges_[i]->exts_[0]->id_ == 1372544128 )
    {
        delete edges_[i];
        edges_.erase( edges_.begin() + i-- );
    }
    if ( edges_.size() > 1 ) for ( int i = 0; i < edges_.size(); i++ ) if ( edges_[i]->ol_ == 92 && edges_[i]->exts_[0]->id_ == 912682924 )
    {
        delete edges_[i];
        edges_.erase( edges_.begin() + i-- );
    }
    if ( edges_.size() > 1 ) for ( int i = 0; i < edges_.size(); i++ ) if ( edges_[i]->ol_ == 132 && edges_[i]->exts_[0]->id_ == 55861190 )
    {
        delete edges_[i];
        edges_.erase( edges_.begin() + i-- );
    }
    if ( edges_.size() > 1 ) for ( int i = 0; i < edges_.size(); i++ ) if ( edges_[i]->ol_ == 132 && edges_[i]->exts_[0]->id_ == 124544548 )
    {
        delete edges_[i];
        edges_.erase( edges_.begin() + i-- );
    }
    
    if ( edges_.size() == 1 )
    {
        vector<BudAlign*> edges = edges_[0]->edges_;
        edges_[0]->edges_.clear();
        exts_.insert( exts_.begin(), edges_[0]->exts_.begin(), edges_[0]->exts_.end() );
        edges_[0]->exts_.clear();
        maps_.insert( maps_.end(), edges_[0]->maps_.begin(), edges_[0]->maps_.end() );
        edges_[0]->maps_.clear();
        ext_ = edges_[0]->ext_;
        good_ = edges_[0]->good_;
        len_ = edges_[0]->len_;
        delete edges_[0];
        edges_ = edges;
    }
    
//    if ( edges_.empty() && good_ >= len_ ) extend( coord_ + ( drxn ? len_ : -len_ ), drxn );
    
    for ( BudExt* be : exts_ ) ol_ = max( ol_, be->ol_ );
    for ( BudMapx* bm : maps_ ) ol_ = max( ol_, drxn ? coord_ - bm->q_[0] : bm->q_[1] - coord_ );
    return good_;
}

bool BudAlign::used( ReadId id )
{
    for ( BudExt* be : exts_ ) if ( be->id_ == id ) return true;
    for ( BudAlign* ba : edges_ ) if ( ba->used( id ) ) return true;
    return false;
}

//BudFork::~BudFork()
//{
//    delete align_;
//}
//
//void BudFork::add( string& seq, ReadId id, int ol, bool drxn )
//{
//    if ( align_->used( id ) ) return;
//    assert( align_->addExt( new BudExt( seq, id, ol, drxn ), drxn ) );
//}

bool BudNode::isChained( int chained, int good, int tar, bool drxn )
{
    for ( pair<BudNode*, pair<int, int> >& edge : edges_[!drxn] )
    {
        if ( edge.second.second < good && edge.second.first + edge.second.second + chained >= tar ) return true;
    }
    for ( pair<BudNode*, pair<int, int> >& edge : edges_[!drxn] )
    {
        if ( edge.second.second < good && edge.first->isChained( chained + edge.second.second, edge.second.first, tar, drxn ) )
        {
            return true;
        }
    }
    return false;
}

void BudNode::removeEdge( BudNode* bn, bool drxn )
{
    for ( int i = 0; i < edges_[drxn].size(); i++ ) if ( edges_[drxn][i].first == bn )
    {
        edges_[drxn].erase( edges_[drxn].begin() + i );
        return;
    }
    assert( false );
}

void BudNode::setBad( bool drxn )
{
    bad_ = true;
    if ( !edges_[drxn].empty() ) good_ = edges_[drxn][0].second.first;
    bad_ = !good_ || !isChained( 0, good_, good_, drxn );
}

void BudNode::setEdges( bool drxn )
{
    if ( bad_ ) return;
    for ( int i = 0; i < edges_[drxn].size(); i++ ) for ( int j = i+1; j < edges_[drxn].size(); j++ )
    {
        if ( edges_[drxn][i].first->bad_ || edges_[drxn][j].first->bad_ ) continue;
        if ( edges_[drxn][i].second.first <= edges_[drxn][j].second.first ) continue;
        bool cull = false;
        for ( pair<BudNode*, pair<int, int> >& edge : edges_[drxn][j].first->edges_[!drxn] ) if ( !edge.first->bad_ )
        {
            if ( edge.second.first > edges_[drxn][j].second.first ) cull = true;
        }
        if ( !cull && edges_[drxn][j].first->good_ > edges_[drxn][j].second.first + edges_[drxn][j].second.second )
        {
            cull = true;
        }
        if ( !cull ) continue;
        edges_[drxn][j].first->removeEdge( this, !drxn );
        edges_[drxn].erase( edges_[drxn].begin() + j-- );
    }
    
}

void BudNode::setGood( bool drxn )
{
    sort( edges_[drxn].begin(), edges_[drxn].end(), [&]( pair<BudNode*, pair<int, int> >& a, pair<BudNode*, pair<int, int> >& b )
    { return a.second.first == b.second.first ? a.second.second < b.second.second : a.second.first > b.second.first; } );
    sort( edges_[!drxn].begin(), edges_[!drxn].end(), [&]( pair<BudNode*, pair<int, int> >& a, pair<BudNode*, pair<int, int> >& b )
    { return a.second.first + a.second.second > b.second.first + b.second.second; } );
    if ( edges_[drxn].size() > 1 )
    {
        int x = 0;
    }
    if ( edges_[!drxn].size() > 1 )
    {
        int x = 0;
    }
}

BudIsle::BudIsle( BudNode* bn, vector<BudIsle*>& isles, bool drxn )
: path_{ bn }
{
    isles.push_back( this );
}

void BudIsle::addEdge( BudNode* bn, vector<BudIsle*>& isles, int ol, bool drxn )
{
    BudIsle* edge = NULL;
    for ( BudIsle* bi : isles ) if ( ( drxn ? bi->path_[0] : bi->path_.back() ) == bn ) edge = bi;
    bool created = !edge;
    if ( created ) edge = new BudIsle( bn, isles, drxn );
    edges_[drxn].push_back( make_pair( edge, ol ) );
    edge->edges_[!drxn].push_back( make_pair( this, ol ) );
    if ( created ) edge->extend( isles, drxn );
}

void BudIsle::extend( vector<BudIsle*>& isles, bool drxn )
{
    for ( BudNode* bn = ( drxn ? path_.back() : path_[0] ); bn; )
    {
        vector< pair<BudNode*, pair<int, int> > > edges[2];
        for ( pair<BudNode*, pair<int, int> >& edge : bn->edges_[drxn] ) if ( !edge.first->bad_ ) edges[drxn].push_back( edge );
        bn = edges[drxn].size() == 1 ? edges[drxn][0].first : NULL;
        if ( bn ) for ( pair<BudNode*, pair<int, int> >& edge : bn->edges_[!drxn] ) if ( !edge.first->bad_ ) edges[!drxn].push_back( edge );
        if ( edges[!drxn].size() > 1 ) bn = NULL;
        
        if ( bn ) path_.insert( drxn ? path_.end() : path_.begin(), bn );
        else for ( pair<BudNode*, pair<int, int> >& edge : edges[drxn] ) addEdge( edge.first, isles, edge.second.first, drxn );
    }
}

BudNode* BudGraph::get( BudPair* bp, bool create )
{
    auto it = nodes_.find( bp );
    if ( it != nodes_.end() ) return it->second;
    if ( !create ) return NULL;
    BudNode* bn = new BudNode( bp );
    assert( nodes_.insert( make_pair( bp, bn ) ).second );
    return bn;
}

BudRead::BudRead( int32_t l, int32_t r, int32_t len, bool drxn )
{
    coords_[0] = l;
    coords_[1] = r;
    coords_[2] = drxn ? ( r - l ) - len : len - ( r - l );
}

BudEdge::BudEdge( Bud* l, Bud* r, int32_t diff, int ol )
: diff_( diff ), ol_( ol )
{
    edge_[0] = l;
    edge_[1] = r;
    off_[0] = off_[1] = 0;
}

Bud::Bud( NodePath* np, vector<Bud*>& buds, int32_t target, bool budding, bool drxn )
: pathed_{ np }, target_( target ), bud_( budding ), mapped_( false )
{
    ends_[0] = ends_[1] = np->ends_[drxn];
    buds.push_back( this );
    if ( budding ) setTarget( np, drxn );
}

Bud::Bud( Querier& bwt, vector<Bud*>& buds, Bud* b, BudAlign* ba, bool drxn )
: bud_( true ), mapped_( false )
{
    int ol = ba->ol_ + ba->base_, off = ba->base_ ? 0 : abs( b->ends_[drxn] - ba->coord_ );
    seq_ = drxn ? b->seq_.substr( b->seq_.size() - off - ol, ol ) + ba->ext_ : ba->ext_ + b->seq_.substr( off, ol );
    block_.resize( seq_.size(), false );
    ends_[!drxn] = ba->coord_ + ( drxn ? -ol : ol );
    ends_[drxn] = ba->coord_ + ( drxn ? ba->len_ : -ba->len_ );
    addAlign( bwt, buds, ba, drxn );
    BudEdge* be = new BudEdge( drxn ? b : this, drxn ? this : b, 0, ol );
    be->off_[!drxn] = off;
    b->edges_[drxn].push_back( be );
    edges_[!drxn].push_back( be );
    buds.push_back( this );
    assert( ba->edges_.empty() );
    for ( BudAlign* edge : ba->edges_ ) new Bud( bwt, buds, this, edge, drxn );
}

Bud::~Bud()
{
    for ( pair<ReadId, pair<BudPair*, int32_t> > bp : pairs_ ) if ( ( --bp.second.first->count_ ) < 1 ) delete bp.second.first;
    for ( int d : { 0, 1 } ) for ( BudEdge* be: edges_[d] )
    {
        be->edge_[d]->edges_[!d].erase( remove( be->edge_[d]->edges_[!d].begin(), be->edge_[d]->edges_[!d].end(), be ), be->edge_[d]->edges_[!d].end() );
        delete be;
    }
}

void Bud::addAlign( Querier& bwt, vector<Bud*>& buds, BudAlign* ba, bool drxn )
{
    target_ = ba->coord_ + ( drxn ? ba->good_ : -ba->good_ );
    for ( BudMapx* bm : ba->maps_ ) if ( ( drxn ? bm->t_[1] == bm->pair_->seq_.size() : !bm->t_[0] ) && reads_.find( bm->pair_->id_ ) == reads_.end() )
    {
        reads_.insert( make_pair( bm->pair_->id_, BudRead( bm->q_[0], bm->q_[1], bm->pair_->seq_.size(), drxn ) ) );
    }
    Lib* lib;
    for ( BudExt* be : ba->exts_ ) if ( reads_.find( be->id_ ) == reads_.end() )
    {
        ReadId id = be->id_;
        int32_t coords[2]{ ba->coord_ - ( drxn ? be->ol_ : be->align_ ), ba->coord_ + ( drxn ? be->align_ : be->ol_ ) };
        reads_.insert( make_pair( id, BudRead( coords[0], coords[1], be->ext_.size() + be->ol_, !drxn ) ) );
        if ( ( lib = params.getLib( id ) ) && lib->isPe && lib->getPair( id ) == drxn )
        {
            string seq = bwt.getSequence( id );
            addPair( buds, seq, id, coords[!drxn], drxn ? lib->size : -lib->size, drxn );
        }
    }
}

void Bud::addBranch( vector<Bud*>& buds, string seq[2], ReadId id[2], int ol, int32_t dist, bool drxn )
{
    int32_t coords[2];
    if ( !mapSeqEnd( seq[0], seq_, ol, coords, !drxn ) ) return;
    for ( int d : { 0, 1 } ) coords[d] += ends_[0];
//    if ( !seq[1].empty() ) pairs_.insert( make_pair( id[1], BudPair( seq[1], coords[!drxn] + dist, isUnique( coords ) ) ) );
    if ( !seq[1].empty() ) addPair( buds, seq[1], id[1], coords[!drxn], dist, drxn );
    if ( !bud_ || ( drxn ? coords[1] < target_ : target_ < coords[0] ) ) return;
    if ( id[0] == 1336325264 )
    {
        seq[0] = "TAGAGAGAGAGAGAGAGAGAGAGAAAGAAAAAAAAACGAGAGGGAGTGAGAAAGTGAGAGTGAGGGATATAGGTATTATATATTGTCTGATTATATATAGAGAGAGAGAGAGAGGGGGGGGGGAGGGAGAGACAAAGAGACAAAGACAGG";
    }
    if ( id[0] == 1664321330 )
    {
        seq[0] = "AGAGAGAGAGAGAGAGAAAGAAAAAAAAACGAGAGGGAGTGAGAAAGTGAGAGTGAGGGATATAGGTATTATATATTGTCTGATTATATATAGAGAGAGAGAGAGAGGGGGGGGGGAGGGAGAGACAAAGAGACAAAGACAGAGATGAGG";
    }
    
    BudAlign* fork = NULL;
    for ( BudAlign* ba : forks_ ) if ( ba->coord_ == coords[drxn] && ( fork = ba ) ) break;
    if ( !fork && ( fork = new BudAlign( coords[drxn], 0, 0 ) ) ) forks_.push_back( fork );
    if ( !fork->used( id[0] ) ) assert( fork->addExt( new BudExt( seq[0], id[0], coords[1] - coords[0], drxn ), drxn ) );
}

bool Bud::addExt( Querier& bwt, vector<Bud*>& buds, BudAlign* ba, bool drxn )
{
    if ( ba->coord_ != ends_[drxn] || ba->ext_.empty() || !path_.empty() || !edges_[drxn].empty() ) return false;
    ends_[drxn] += drxn ? ba->len_ : -ba->len_;
    seq_ = drxn ? seq_ + ba->ext_ : ba->ext_ + seq_;
    block_.insert( drxn ? block_.end() : block_.begin(), ba->len_, false );
    addAlign( bwt, buds, ba, drxn );
    if ( ba->advance( 0, drxn ) )
    {
        if ( ba->forked_ ) return false;
        forks_.push_back( ba );
    }
    else delete ba;
    bud_ = true;
    mapped_ = false;
    return true;
}

void Bud::addPair( vector<Bud*>& buds, string& seq, ReadId id, int32_t coord, int32_t dist, bool drxn )
{
    BudPair* bp = NULL;
    for ( Bud* b : buds )
    {
        auto it = b->pairs_.find( id );
        if ( it == b->pairs_.end() ) continue;
        if ( b == this )
        {
            if ( drxn ? it->second.second < coord : coord < it->second.second ) it->second.second = coord;
            return;
        }
        assert( !bp || bp == it->second.first );
        bp = it->second.first;
    }
    if ( !bp ) bp = new BudPair( seq, id, dist, false );
    bp->count_++;
    assert( pairs_.insert( make_pair( id, make_pair( bp, coord ) ) ).second );
}

bool Bud::blockPair( Querier& bwt, vector<Bud*>& buds, string& seq, bool drxn )
{
    int32_t coords[2];
    int ol = 0;
    for ( Bud* b : buds ) if ( !b->offs_.empty() && mapSeqEnd( seq, b->seq_, 32, coords, drxn ) ) ol = max( ol, coords[1] - coords[0] );
    return ol && ( ol == seq.size() || !bwt.isExtendable( seq, ol+1, drxn ) );
}

void Bud::build( vector<Bud*>& buds, bool drxn )
{
    assert( pathed_.size() == 1 && path_.empty() && offs_.empty() );
    ends_[!drxn] = ends_[drxn] + ( drxn ? -pathed_[0]->path_.back()->size() : pathed_[0]->path_[0]->size() );
    
    for ( NodePath* np = pathed_[0]; np && edges_[!drxn].empty(); )
    {
        for ( int i = 0; i < np->path_.size(); i++ )
        {
            Node* node = drxn ? np->path_.end()[-i-1] : np->path_[i];
            path_.insert( drxn ? path_.begin() : path_.end(), node );
            assert( offs_.insert( make_pair( node, ends_[!drxn] - node->ends_[!drxn] ) ).second );
            if ( abs( ends_[!drxn] - target_ ) > params.maxPeMean * 1.3 + 300 ) np = NULL;
            if ( !np || i+1 >= np->path_.size() ) break;
            Edge e = node->getEdge( drxn ? np->path_.end()[-i-2] : np->path_[i+1], !drxn );
            assert( e.node );
            ends_[!drxn] += ( e.node->size() - e.ol ) * ( drxn ? -1 : 1 );
        }
        
        PathEdge* edge = np && np->edges_[!drxn].size() == 1 ? np->edges_[!drxn][0] : NULL;
        if ( edge && edge->edge[!drxn]->edges_[drxn].size() == 1 && ( np = edge->edge[!drxn] ) )
        {
            assert( false );
            ends_[!drxn] += ( drxn ? -np->path_.back()->size() + edge->ol : np->path_[0]->size() - edge->ol );
            pathed_.insert( drxn ? pathed_.begin() : pathed_.end(), np );
        }
        else if ( np )
        {
            for ( PathEdge* pe : np->edges_[!drxn] ) setEdge( pe, buds, !drxn );
            np = NULL;
        }
    }
    
    seq_ = Node::getSeq( path_ );
    block_.resize( seq_.size(), false );
    for ( int i = params.maxPeMean * 1.3 + 300 - abs( target_ - ends_[drxn] ); i < block_.size(); i++ ) ( drxn ? block_.end()[-i-1] : block_[i] ) = true;
    assert( ends_[1] - ends_[0] == seq_.size() );
}

void Bud::getDiffs( unordered_map<Bud*, int32_t>& diffs, int32_t diff, bool drxn )
{
    auto ins = diffs.insert( make_pair( this, diff ) );
    if ( !ins.second )
    {
        if ( drxn ? ins.first->second <= diff : diff <= ins.first->second ) return;
        ins.first->second = diff;
    }
    for ( BudEdge* be : edges_[drxn] ) be->edge_[drxn]->getDiffs( diffs, diff + ( drxn ? -be->diff_ : be->diff_ ), drxn );
}

Kmers Bud::getKmers( unordered_map<Bud*, int32_t>& diffs, bool drxn )
{
    Kmers kmers;
    unordered_set<ReadId> used;
    for ( pair<Bud*, int32_t> bd : diffs ) for ( pair<ReadId, pair<BudPair*, int32_t> > bp : bd.first->pairs_ ) if ( used.insert( bp.first ).second )
    {
        kmers.add( bp.second.first );
    }
    return kmers;
}

bool Bud::isUnique( int32_t coords[2] )
{
    return false;
}

void Bud::print( vector<Bud*>& buds, bool drxn )
{
    unordered_set<ReadId> used;
    int32_t limit = buds[0]->ends_[0];
    int32_t target = drxn ? std::numeric_limits<int32_t>::max() : std::numeric_limits<int32_t>::min();
    for ( Bud* b : buds ) if ( b->bud_ ) target = drxn ? min( target, b->target_ ) : max( target, b->target_ );
    for ( Bud* b : buds )
    {
        limit = min( limit, b->ends_[0] );
        for ( pair<ReadId, pair<BudPair*, int32_t> > br : b->pairs_ ) limit = min( limit, br.second.first->dist_ + br.second.second );
        for ( pair<ReadId, BudRead> br : b->reads_ ) used.insert( br.first );
    }
    if ( !drxn ) limit -= params.readLen;
    int i = 1;
    for ( Bud* b : buds ) cout << ">" << i++ << "Node: " << b->ends_[0] << " - " << b->ends_[1] << endl << string( b->ends_[0] - limit, '-' ) << b->seq_ << endl;
    for ( Bud* b : buds ) for ( BudAlign* ba : b->forks_ ) ba->print( limit, drxn );
    for ( Bud* b : buds ) for ( pair<ReadId, pair<BudPair*, int32_t> > br : b->pairs_ ) if ( used.find( br.first ) == used.end() )
    {
        int32_t est = br.second.first->dist_ + br.second.second - ( drxn ? br.second.first->seq_.size() : 0 );
        if ( drxn ? est + 100 < target : target < est - 100 ) continue;
        cout << ">" << i++ << "Pair " << br.first << ": " << est << endl << string( est - limit, '-' ) << br.second.first->seq_ << endl;
    }
    int x = 0;
}

void Bud::rebuild( vector<Bud*>& buds, int32_t target, bool drxn )
{
    assert( false );
    if ( drxn ? target_ < target + 100 : target - 100 < target_ ) return;
    target_ = target;
    if ( !edges_[!drxn].empty() )
    {
        for ( BudEdge* be : edges_[!drxn] ) be->edge_[!drxn]->rebuild( buds, target_ + ( drxn ? be->diff_ : -be->diff_ ), drxn );
        return;
    }
    pathed_.erase( pathed_.begin() + drxn, pathed_.end() - ( !drxn ) );
    path_.clear();
    offs_.clear();
    seq_.clear();
    block_.clear();
}

vector< pair<ReadId, int> > Bud::setAligned( BudAlign* ba, vector<Bud*>& buds, Kmers& kmers, bool drxn )
{
    vector< pair<ReadId, int> > maps;
    int32_t off = abs( ends_[drxn] - ba->coord_ );
    for ( BudAlign* edge : ba->edges_ ) for ( pair<ReadId, int> m : setAligned( edge, buds, kmers, drxn ) ) maps.push_back( m );
    string seq = drxn ? seq_.substr( 0, seq_.size() - off ) + ba->ext_ : ba->ext_ + seq_.substr( off );
    assert( ba->ext_.size() == ba->len_ );
    
    vector<Kmer>* pk;
    for ( int i = 0; i < ba->len_ - ba->base_; i++ )
    {
        for ( int j = 0; j < maps.size(); j++ ) if ( ba->len_ - i < maps[j].second ) maps.erase( maps.begin() + j-- );
        
        string q = drxn ? seq.substr( seq.size()-i-32, 32 ) : seq.substr( i, 32 );
        if ( pk = kmers.get( q ) ) for ( Kmer& kmer : *pk )
        {
            bool bad = false;
            for ( pair<ReadId, int> m : maps ) if ( bad = ( m.first == kmer.bp_->id_ ) ) break;
            if ( bad ) continue;
            string t = kmer.bp_->seq_;
            
            int qq[2]{ drxn ? (int)seq.size()-i-32 : i, drxn ? (int)seq.size()-i : i+32 }, tt[2]{ kmer.off_, kmer.off_+32 };
            int len = setMatch( seq, t, qq, tt );
            int exted = ba->len_ - ( drxn ? seq.size() - qq[1] : qq[0] );
            bool breaks[2]{ qq[0] && tt[0], qq[1] < seq.size() && tt[1] < t.size() };
            
//            bad = bb->ext_->id_ == kmer.id_;
//            for ( int j = 0; bad && j < bb->exts_.size() && bb->exts_[i]->ext_.size() == bb->ext_->ext_.size(); j++ )
//            {
//                if ( bb->exts_[i]->id_ != kmer.id_ ) bad = false;
//            }
            
//            if ( !bad && len > exted+31 )
//            {
//                bb->good_ = max( bb->good_, exted );
//                if ( exted == bb->good_ && bb->bad_ != bb->good_ ) bb->bad_ = breaks[drxn] ? exted : -1;
//            }
            
            for ( int d : { 0, 1 } ) qq[d] = qq[d] + ends_[drxn] - ( drxn ? seq_ : ba->ext_ ).size();
            ba->addMap( kmer.bp_, qq, tt, drxn );
            maps.push_back( make_pair( kmer.bp_->id_, exted - len + 32 ) );
        }
    }
    
    return maps;
}

void Bud::setEdge( PathEdge* pe, vector<Bud*>& buds, bool drxn )
{
    BudEdge* edge = new BudEdge( this, this, drxn ? pe->edge[1]->ends_[0] + pe->ol - ends_[1] : ends_[0] + pe->ol - pe->edge[0]->ends_[1], pe->ol );
    
    bool built = false;
    int32_t target = target_ + ( drxn ? edge->diff_ : -edge->diff_ );
    for ( Bud* b : buds ) if ( ( drxn ? b->pathed_[0] == pe->edge[1] : b->pathed_.back() == pe->edge[0] ) && ( edge->edge_[drxn] = b ) ) built = true;
    if ( edge->edge_[0] == edge->edge_[1] ) edge->edge_[drxn] = new Bud( pe->edge[drxn], buds, target, false, !drxn );
    for ( int d : { 0, 1 } ) edge->edge_[d]->edges_[!d].push_back( edge );
    if ( built ) rebuild( buds, target, !drxn );
    else edge->edge_[drxn]->build( buds, !drxn );
    assert( edge->ol_ < edge->edge_[drxn]->block_.size() );
    for ( int i = 0; i < edge->ol_; i++ ) ( drxn ? edge->edge_[1]->block_[i] : edge->edge_[0]->block_.end()[-i-1] ) = true;
}

void Bud::setExt( Querier& bwt, vector<Bud*>& buds, bool drxn )
{
    assert( forks_.size() == 1 );
    
    vector<BudAlign*> forks = forks_;
    forks_.clear();
    bud_ = false;
    
    for ( int i = 0; i < forks.size(); i++ ) if ( addExt( bwt, buds, forks[i], drxn ) ) forks.erase( forks.begin() + i-- );
    
    for ( BudAlign* ba : forks )
    {
        if ( !ba->ext_.empty() ) new Bud( bwt, buds, this, ba, drxn );
        else for ( BudAlign* edge : ba->edges_ )
        {
            new Bud( bwt, buds, this, edge, drxn );
        }
        delete ba;
    }
}

bool Bud::setForks( Querier& bwt, vector<Bud*>& buds, bool drxn )
{
    vector<Bud*> budding;
    setIslands( bwt, buds, drxn );
    for ( Bud* b : buds ) if ( b->bud_ )
    {
        unordered_map<Bud*, int32_t> diffs;
        b->getDiffs( diffs, 0, drxn );
        Kmers kmers = b->getKmers( diffs, drxn );
        for ( int i = 0; i < b->forks_.size(); i++ ) if ( b->forks_[i]->coord_ != b->ends_[drxn] && !b->forks_[i]->isBranch( 5 ) )
        {
            delete b->forks_[i];
            b->forks_.erase( b->forks_.begin() + i-- );
        }
        for ( BudAlign* ba : b->forks_ ) b->setAligned( ba, buds, kmers, drxn );
        for ( BudAlign* ba : b->forks_ ) ba->setGood( 0, drxn );
        for ( BudAlign* ba : b->forks_ ) ba->setExt( drxn );
        for ( BudAlign* ba : b->forks_ ) ba->setForks();
        budding.push_back( b );
    }
    for ( Bud* b : budding ) if ( b->bud_ ) b->setExt( bwt, buds, drxn );
    for ( Bud* b : buds ) if ( !b->mapped_ && b->bud_ ) return true;
    assert( false );
    return false;
}

void Bud::setIslands( Querier& bwt, vector<Bud*>& buds, bool drxn )
{
    vector< pair<Bud*, unordered_map<Bud*, int32_t> > > budding;
    unordered_set<BudPair*> added;
    Kmers kmers;
    for ( Bud* b : buds ) if ( b->bud_ )
    {
        unordered_map<Bud*, int32_t> diffs;
        b->getDiffs( diffs, 0, drxn );
        unordered_set<ReadId> block;
        for ( pair<Bud*, int32_t> bd : diffs ) for ( const pair<ReadId, BudRead>& br : bd.first->reads_ ) block.insert( br.first );
        for ( pair<Bud*, int32_t> bd : diffs ) for ( pair<ReadId, pair<BudPair*, int32_t> > bp : bd.first->pairs_ ) if ( bp.second.first->seq_.size() > 32 )
        {
            if ( block.find( bp.first ) != block.end() )  continue;
            int32_t est = bp.second.second + ( bp.second.first->dist_ * 1.3 );
            if ( drxn ? est + bd.second + 200 < b->target_ : b->target_ < est - bd.second - 200 ) continue;
            if ( added.insert( bp.second.first ).second ) kmers.add( bp.second.first );
        }
        budding.push_back( make_pair( b, diffs ) );
    }
    
    vector<Kmer>* pk;
    BudGraph bg;
    for ( BudPair* bp : added )
    {
        string q = drxn ? bp->seq_.substr( bp->seq_.size() - 32 ) : bp->seq_.substr( 0, 32 );
        assert( pk = kmers.get( q ) );
        BudNode* bn[2]{ NULL, NULL };
        vector< pair<BudPair*, pair<int, int> > > edges;
        for ( Kmer& k : *pk ) if ( k.bp_ != bp )
        {
            int32_t qq[2]{ drxn ? (int)bp->seq_.size() - 32 : 0, drxn ? (int)bp->seq_.size() : 32 }, tt[2]{ k.off_, k.off_+32 };
            edges.push_back( make_pair( k.bp_, make_pair( setMatch( bp->seq_, k.bp_->seq_, qq, tt ), drxn ? (int)k.bp_->seq_.size()-tt[1] : tt[0] ) ) );
//            setMatch( bp->seq_, k.bp_->seq_, qq, tt );
//            if ( !bn[0] ) bn[0] = bg.get( bp, true );
//            bn[1] = bg.get( k.bp_, true );
//            int len = tt[1] - tt[0], excess = drxn ? (int)k.bp_->seq_.size()-tt[1] : tt[0];
//            bn[0]->edges_[drxn].push_back( make_pair( bn[1], make_pair( len, drxn ? (int)k.bp_->seq_.size()-tt[1] : tt[0] ) ) );
//            bn[1]->edges_[!drxn].push_back( make_pair( bn[0], make_pair( len, drxn ? qq[0] : (int)bp->seq_.size()-qq[1] ) ) );
        }
        sort( edges.begin(), edges.end(), [&]( pair<BudPair*, pair<int, int> >& a, pair<BudPair*, pair<int, int> >& b )
        { return a.second.first == b.second.first ? a.second.second < b.second.second : a.second.first > b.second.first; } );
        for ( int i = 0; i < edges.size(); i++ ) for ( int j = i+1; j < edges.size(); j++ )
        {
            if ( edges[i].first == edges[j].first || edges[i].first->alignEnds( edges[j].first, edges[i].second.second, edges[j].second.second, drxn ) ) edges.erase( edges.begin() + j-- );
        }
        for ( pair<BudPair*, pair<int, int> >& edge : edges )
        {
            if ( !bn[0] ) bn[0] = bg.get( bp, true );
            bn[1] = bg.get( edge.first, true );
            bn[0]->edges_[drxn].push_back( make_pair( bn[1], edge.second ) );
            bn[1]->edges_[!drxn].push_back( make_pair( bn[0], edge.second ) );
        }
    }
    
    vector<BudIsle*> seeds, isles;
    for ( pair<BudPair*, BudNode*> bpn : bg.nodes_ ) bpn.second->setBad( drxn );
    for ( pair<BudPair*, BudNode*> bpn : bg.nodes_ ) bpn.second->setEdges( drxn );
    for ( pair<BudPair*, BudNode*> bpn : bg.nodes_ ) if ( !bpn.second->bad_ )
    {
        bool bad = false;
        for ( pair<BudNode*, pair<int, int> >& edge : bpn.second->edges_[!drxn] ) if ( !edge.first->bad_ ) bad = true;
        if ( bad ) continue;
        BudIsle* bi = new BudIsle( bpn.second, isles, drxn );
        bi->extend( isles, drxn );
        seeds.push_back( bi );
    }
    for ( BudIsle* bi : isles ) delete bi;
}

void Bud::setMaps( Querier& bwt, vector<Bud*>& buds, bool drxn )
{
    Lib* lib;
    unordered_set<ReadId> used;
    for ( Bud* b : buds ) if ( !b->mapped_ ) for ( Node* node : b->path_ ) for ( auto& read : node->reads_ ) used.insert( read.first );
    for ( Bud* b : buds ) if ( !b->mapped_ ) for ( const pair<ReadId, BudRead>& br : b->reads_ ) used.insert( br.first );
    for ( int i = 0; i < buds.size(); i++ ) if ( !buds[i]->mapped_ )
    {
        for ( QueryBranch& qb : bwt.mapBranches( buds[i]->seq_, &buds[i]->block_, 32, drxn ) )
        {
            for ( ReadId id : bwt.getIds( qb.rank, qb.count, drxn ) ) if ( used.find( id ) == used.end() )
            {
                string seq[2]{ bwt.getSequence( id ), "" };
                int ol = qb.ol;
                int32_t coords[2];
                for ( Bud* b : buds ) if ( mapSeqEnd( seq[0], b->seq_, ol, coords, !drxn ) && ( ol = coords[1] - coords[0] ) == seq[0].size() ) break;
                if ( ol == seq[0].size() ) continue;
                
                ReadId pairs[2]{ id, id };
                if ( ( lib = params.getLib( pairs[0] ) ) && lib->isPe && lib->getPair( pairs[1] ) == drxn ) seq[1] = bwt.getSequence( pairs[1] );
                if ( !seq[1].empty() && blockPair( bwt, buds, seq[1], drxn ) ) seq[1].clear();
                for ( Bud* b : buds ) b->addBranch( buds, seq, pairs, ol, lib ? ( drxn ? lib->size : -lib->size ) : 0, drxn );
            }
            
            string q = buds[i]->seq_.substr( qb.i, qb.ol );
            for ( int j = 0; j+31 < qb.ol; j++ ) buds[i]->block_[ qb.i + ( drxn ? j : j+31 ) ] = true;
            for ( int j = i+1; j < buds.size(); j++ )
            {
                size_t it = buds[j]->seq_.find( q );
                if ( it != string::npos )
                {
                    int coords[2]{ (int)it, (int)it+qb.ol };
                    if ( drxn ) for ( int k = qb.i; coords[0] && --k > 0 && buds[j]->seq_[ coords[0]-1 ] == buds[i]->seq_[k]; ) coords[0]--;
                    if ( !drxn ) for ( int k = qb.i+qb.ol; coords[1] < buds[j]->seq_.size() && k < buds[i]->seq_.size() && buds[j]->seq_[ coords[1] ] == buds[i]->seq_[k]; k++ ) coords[1]++;
                    for ( int k = coords[0]; k+31 < coords[1]; k++ ) buds[j]->block_[ drxn ? k : k+31 ] = true;
                }
            }
        }
        int block = buds[i]->block_.size();
        for ( int j = 0; j < buds[i]->block_.size(); j++ ) if ( buds[i]->block_[j] )
        {
            while ( block < j ) buds[i]->block_[block++] = true;
            block = j+1;
        }
        buds[i]->mapped_ = true;
    }
}

int Bud::setMatch( string& q, string& t, int qq[2], int tt[2] )
{
    for ( ; qq[0] > 0 && tt[0] > 0 && q[ qq[0]-1 ] == t[ tt[0]-1 ]; qq[0]-- ) tt[0]--;
    for ( ; qq[1] < q.size() && tt[1] < t.size() && q[ qq[1] ] == t[ tt[1] ]; qq[1]++ ) tt[1]++;
    assert( qq[1] - qq[0] == tt[1] - tt[0] );
    return qq[1] - qq[0];
}

void Bud::setReads( Querier& bwt, vector<Bud*>& buds, bool drxn )
{
    unordered_set<ReadId> used;
    for ( Bud* b : buds ) for ( Node* node : b->path_ ) for ( auto& read : node->reads_ ) used.insert( read.first );
    Lib* lib;
    for ( Bud* b : buds ) for ( auto& no : b->offs_ ) for ( auto& read : no.first->reads_ )  if ( ( lib = params.getLib( read.first ) ) && lib->isPe )
    {
        ReadId id = read.first;
        int32_t est = read.second[!drxn] + no.second + ( drxn ? lib->size : -lib->size );
        if ( lib->getPair( id ) != drxn || b->pairs_.find( id ) != b->pairs_.end() ) continue;
        if ( used.find( id ) != used.end() ) continue;
        string seq = bwt.getSequence( id );
        b->addPair( buds, seq, id, read.second[!drxn], drxn ? lib->size : -lib->size, drxn );
    }
}

void Bud::setTarget( NodePath* np, bool drxn )
{
    target_ = np->ends_[1];
    for ( int32_t diff = 0; np; )
    {
        break;
    }
}

bool Bud::bud( Querier& bwt, NodeRoll& nodes, vector<NodePath*>& ends, bool drxn )
{
    vector<Bud*> budding, buds;
    for ( NodePath* np : ends ) budding.push_back( new Bud( np, buds, np->ends_[drxn], true, drxn ) );
    for ( Bud* b : budding ) b->build( buds, drxn );
    setMaps( bwt, buds, drxn );
    setReads( bwt, buds, drxn );
    print( buds, drxn );
    while ( setForks( bwt, buds, drxn ) )
    {
        setMaps( bwt, buds, drxn );
        print( buds, drxn );
    }
    assert( false );
}
