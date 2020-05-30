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

Vertex::Vertex( string seq, int coord, bool drxn )
: seq_( seq )
{
    coord_[!drxn] = coord;
    coord_[drxn] = drxn ? coord + seq.size() : coord - seq.size();
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

Traverse::Traverse( Node* node, bool drxn )
: path_{ node }, coords_{ make_pair( node, drxn ? -node->size() : 0 ) }, len_( node->size() )
{
    ends_[0] = ends_[1] = 0;
    for ( vector<Edge> edges = node->edges_[!drxn]; edges.size() == 1 && ( node = edges[0].node ) && node->edges_[drxn].size() == 1; )
    {
        len_ = len_ + node->size() - edges[0].ol;
        path_.insert( drxn ? path_.begin() : path_.end(), node );
        coords_.insert( make_pair( node, drxn ? -len_ : len_ - node->size() ) );
        edges = edges[0].node->edges_[!drxn];
    }
    ends_[!drxn] = drxn ? -len_ : len_;
}

bool Traverse::trace( Querier& bwt, Node* node, bool drxn )
{
    Traverse trav( node, drxn );
    trav.setFlayed( bwt, drxn );
    trav.setPairs( bwt, drxn );
    trav.setTraces( drxn );
    assert( false );
    return false;
}

void Traverse::print( TFork* fork, bool drxn )
{
    cout << ">Ext: " << fork->ext_.size() << "  |  reads: " << fork->reads_.size() << ( fork->reads_.empty() ? "" : "  | " );
    vector<int> bases;
    for ( pair<ReadId, int> read : fork->reads_ )
    {
        auto it = reads_.find( params.getPairId( read.first ) );
        assert( it != reads_.end() );
        bases.push_back( it->second );
    }
    sort( bases.begin(), bases.end() );
    for ( int i = 0; i < bases.size(); i++ ) cout << " " << bases[i];
    cout << endl << fork->seq_ << endl;
    for ( TFork* tf : fork->exts_ ) print( tf, drxn );
}

void Traverse::setFlayed( Querier& bwt, bool drxn )
{
    string seq = Node::getSeq( path_ );
    seq = seq.substr( seq.size()-500 );
    QueryFlay qf( seq, bwt.ir_, 120, !drxn );
    qf.flay( seq, bwt.qb_, drxn );
    assert( false );
}

void Traverse::setPairs( Querier& bwt, bool drxn )
{
    unordered_map<ReadId, int> reads;
    for ( pair<Node*, int32_t> nc : coords_ ) if ( ( drxn ? -nc.second - nc.first->size() : nc.second ) < params.maxPeMean ) for ( auto& read : nc.first->reads_ )
    {
        int coord = nc.second + read.second[!drxn] - nc.first->ends_[0];
        auto it = reads.insert( make_pair( read.first, coord ) );
        if ( !it.second && abs( coord < it.first->second ) ) it.first->second = coord;
    }
    
    Lib* lib;
    int counted = 0;
    for ( auto& read : reads ) if ( ( lib = params.getLib( read.first ) ) && lib->isPe )
    {
        ReadId id = read.first;
        if ( ( lib->getPair( id ) != drxn ) || ( reads.find( id ) != reads.end() ) ) continue;
        if ( ( drxn ? read.second + lib->size : lib->size - read.second ) < -300 ) continue;
        reads_.insert( make_pair( read.first, read.second ) );
        pairs_.insert( make_pair( id, bwt.getSequence( id ) ) );
    }
}

void Traverse::setTraces( bool drxn )
{
    assert( forks_.empty() );
    forks_.push_back( new TFork( Node::getSeq( path_ ), drxn ) );
    vector<TFork*> forks = forks_;
    for ( int again = 1; again-- > 0; )
    {
        for ( auto& read : pairs_ ) for ( TFork* tf : forks ) tf->add( read.second, read.first, drxn );
        for ( TFork* tf : forks ) if ( tf->update() ) again = 1;
    }
    for ( TFork* tf : forks_ ) print( tf, drxn );
    assert( false );
}
