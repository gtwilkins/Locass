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

#include "locus_fill.h"
#include "shared_functions.h"
#include "seed_node.h"
#include <algorithm>

FillHit::FillHit( Node* hit, int32_t coords[2], int32_t hang, int32_t diff )
: len( coords[1] - coords[0] ), hang( hang ), diff( diff )
{
    fill( hit, coords[0], coords[1] );
}

void FillHit::add( vector<FillHit>& hits, FillNode* fn, string& seq, int32_t coords[2], int32_t est, int32_t off, bool drxn )
{
    int32_t len = coords[1] - coords[0], diff = abs( coords[drxn] + off - est );
    if ( !hits.empty() && hits[0].len < len ) hits.clear();
    if ( !hits.empty() && hits[0].len > len ) return;
    
    Node* hit = NULL;
    for ( int i = 0; !hit && i < fn->path_.size(); i++ ) if ( drxn ? coords[1] <= fn->path_[i].second + fn->path_[i].first->size()
                                                                  : ( i+1 == fn->path_.size() || coords[0] < fn->path_[i+1].second ) )
    {
        hit = fn->path_[i].first;
        if ( drxn ) coords[0] = max( coords[0], fn->path_[i].second );
        if ( !drxn ) coords[1] = min( coords[1], fn->path_[i].second + fn->path_[i].first->size() );
        for ( int d : { 0, 1 } ) coords[d] += hit->ends_[0] - fn->path_[i].second;
        assert( coords[0] < coords[1] );
        assert( fn->path_[i].first->ends_[0] <= coords[0] && coords[1] <= fn->path_[i].first->ends_[1] );
    }
//    for ( int i = 0; !hit && i < fn->path.size(); i++ ) if ( fn->path[i].second <= coords[0] && coords[1] <= fn->path[i].second + fn->path[i].first->size() )
//    {
//        hit = fn->path[i].first;
//        for ( int d : { 0, 1 } ) coords[d] += hit->ends_[0] - fn->path[i].second;
//    }
    assert( hit );
    
    for ( FillHit& fh : hits ) for ( int i = 0; i < fh.hits.size(); i++ )
    {
        if ( fh.hits[i] == hit && fh.mapped[0][i] == coords[0] && fh.mapped[1][i] == coords[1] ) return;
    }
    hits.push_back( FillHit( hit, coords, drxn ? len - (int)seq.size() : seq.size() - len, diff ) );
}

bool FillHit::confirm( vector<FillHit> hits[2], bool ext[2], int seqLen, bool& ignore )
{
    int32_t lens[2]{0};
    for ( int d : { 0, 1 } ) for ( FillHit& fh : hits[d] ) lens[d] = max( lens[d], fh.len );
    for ( int d : { 0, 1 } ) if ( ext[d] && !lens[d] ) return false;
    for ( int d : { 0, 1 } ) if ( ext[d] && !ext[!d] && lens[d] ) hits[!d].clear();
    
    if ( lens[0] == seqLen || lens[1] == seqLen )
    {
        assert( lens[0] == lens[1] );
        for ( int d : { 0, 1 } ) for ( FillHit& fh : hits[d] ) if ( fh.hang ) assert( false );
        for ( int d : { 0, 1 } ) for ( FillHit& fh : hits[d] ) for ( Node* hit : fh.hits ) if ( !hits[0][0].hits[0]->isClone( hit ) ) ignore = true;
        return true;
    }
    
    if ( !hits[0].empty() && !hits[1].empty() )
    {
        if ( ext[0] && ext[1] ) return false;
        hits[ lens[1] < lens[0] ].clear();
    }
    
    for ( int d : { 0, 1 } ) for ( FillHit& fh : hits[d] ) for ( Node* hit : fh.hits ) if ( !hits[d][0].hits[0]->isClone( hit ) ) ignore = true;
    
    return !hits[0].empty() || !hits[1].empty();
}

void FillHit::fill( Node* node, int32_t l, int32_t r )
{
    if ( find( hits.begin(), hits.end(), node ) != hits.end() ) return;
    
    hits.push_back( node );
    mapped[0].push_back( l );
    mapped[1].push_back( r );
    for ( Node* clone : node->clones( false ) ) fill( clone, l + clone->ends_[0] - node->ends_[0], r + clone->ends_[0] - node->ends_[0] );
    
    for ( int d : { 0, 1 } ) for ( Edge& e : node->edges_[d] ) if ( d ? node->ends_[1] - e.ol <= l : r <= node->ends_[0] + e.ol )
    {
        int32_t off = d ? e.node->ends_[0] + e.ol - node->ends_[1] : e.node->ends_[1] - e.ol - node->ends_[0];
        fill( e.node, l + off, r + off );
    }
}

void FillHit::set( ReadId id, bool ignore )
{
    for ( int i = 1; i < hits.size(); i++ ) if ( !hits[0]->isClone( hits[i] ) ) ignore = true;
    for ( int i = 0; i < hits.size(); i++ )
    {
        hits[i]->add( id, mapped[0][i], mapped[1][i], hits[i]->isRedundant( mapped[0][i], mapped[1][i] ), ignore, hang );
    }
}

//FillBranch::FillBranch( string& seq, string& base, ReadId id, int32_t coords[2], bool drxn )
//: ext_( seq.size()-( coords[1]-coords[0] ) )
//{
//    int ol = coords[1] - coords[0];
//    coords[!drxn] = drxn ? max( 0, coords[1]-params.readLen+1 ) : min( (int)base.size(), coords[0]+params.readLen-1 );
//    ol_ = coords[1]-coords[0];
//    seq_ = drxn ? base.substr( coords[0], ol_ ) + seq.substr( ol ) : seq.substr( 0, ext_ ) + base.substr( coords[0], ol_ );
//    exts_.push_back( ext_ );
//}

FillBranch::FillBranch( string seq, FillRead& fr, int base )
: seq_( seq ), exts_{ make_pair( fr, base+seq.size() ) }, base_( base ), good_( 0 ), unpaired_( 0 ), total_( 0 ), culled_( false ), added_( true )
{ }

FillBranch::FillBranch( FillBranch* fb, int ext, bool drxn )
: seq_( drxn ? fb->seq_.substr( ext ) : fb->seq_.substr( 0, fb->seq_.size()-ext ) ), branches_( fb->branches_ ), base_( fb->base_+ext ), good_( 0 ), unpaired_( 0 ), total_( 0 ), culled_( false ), added_( true )
{
    for ( int i = 0; i < fb->exts_.size(); i++ ) if ( fb->exts_[i].second > base_ )
    {
        exts_.push_back( fb->exts_[i] );
        fb->exts_.erase( fb->exts_.begin() + i-- );
    }
    fb->seq_.erase( drxn ? fb->seq_.end()-( seq_.size() ) : fb->seq_.begin(), drxn ? fb->seq_.end() : fb->seq_.begin()+seq_.size() );
    fb->branches_.clear();
}

//FillBranch::FillBranch( FillBranch* fb, string& seq, int ext, bool drxn )
//: seq_( drxn ? fb->seq_ + seq.substr( seq.size()-ext ) : seq.substr( 0, ext ) + fb->seq_ ), ol_( fb->ol_ ), ext_( fb->ext_+ext )
//{
//    exts_.push_back( fb->ext_+ext );
//}

FillBranch::~FillBranch()
{
    for ( FillBranch* fb : branches_ ) delete fb;
}

bool FillBranch::add( string& seq, FillRead& fr, int ext, bool drxn )
{
    int i = 0;
    for ( ; i < seq_.size() && ext < seq.size() && ( drxn ? seq[ext] == seq_[i] : seq.end()[-ext-1] == seq_.end()[-i-1] ); i++ ) ext++;
    if ( !i ) return false;
    
    // Split
    if ( i < seq_.size() && ext < seq.size() )
    {
        branches_.push_back( new FillBranch( this, i, drxn ) );
        branches_.push_back( new FillBranch( drxn ? seq.substr( ext ) : seq.substr( 0, seq.size()-ext ), fr, base_ + i ) );
        assert( branches_.size() == 2 );
        return true;
    }
    
    // Branch
    if ( i == seq_.size() && ext < seq.size() && !branches_.empty() )
    {
        for ( FillBranch* fb : branches_ ) if ( fb->add( seq, fr, ext, drxn ) ) return true;
        branches_.push_back( new FillBranch( drxn ? seq.substr( ext ) : seq.substr( 0, seq.size()-ext ), fr, base_ + i ) );
        return true;
    }
    
    // Extend
    if ( i == seq_.size() && ext < seq.size() && branches_.empty() )
    {
        seq_ = drxn ? seq_ + seq.substr( ext ) : seq.substr( 0, seq.size()-ext ) + seq_;
        i = seq_.size();
    }
    
    exts_.push_back( make_pair( fr, base_ + i ) );
    
    return true;
}

bool FillBranch::confirm( vector<FillBranch*>& branches, bool drxn )
{
    for ( FillBranch* fb : branches ) if ( fb->total_ < 5 && fb->branches_.size() > 1 )
    {
        for ( FillBranch* ffb : fb->branches_ ) delete ffb;
        fb->branches_.clear();
    }
    for ( FillBranch* fb : branches ) confirm( fb->branches_, drxn );
    for ( FillBranch* fb : branches ) fb->merge( drxn );
    for ( int i = 0; i < branches.size(); i++ ) for ( int j = 0; j < branches.size(); j++ ) if ( i != j )
    {
        if ( branches[i]->unpaired_ > max( 2, branches[j]->unpaired_ / 2 ) ) continue;
        if ( branches[i]->total_ >= branches[j]->total_ / 2 ) continue;
        delete branches[i];
        branches.erase( branches.begin() + i-- );
        break;
    }
}

//int FillBranch::confirm( Querier& bwt, int coverage, bool drxn )
//{
//    int confirmed = 0;
//    for ( FillBranch* fb : branches_ ) confirmed = max( confirmed, fb->confirm( bwt, coverage, drxn ) );
//    return confirmed + exts_.size();
//}

bool FillBranch::cull( Querier& bwt, string seq, int coverage, bool drxn )
{
    if ( !good_ )
    {
        seq = drxn ? seq + seq_[0] : seq_.back() + seq;
        int diff = params.readLen / 10;
        int ol = params.readLen - diff;
        int exts = bwt.getExtendable( seq, ol, drxn );
        int cutoff = ( coverage * diff / params.readLen ) / 5;
        if ( exts < 2 ) return culled_ = true;
        if ( exts < cutoff ) return culled_ = true;
        seq = drxn ? seq.substr( 0, seq.size()-1 ) : seq.substr( 1 );
        good_ = 1;
    }
    for ( FillBranch* fb : branches_ ) fb->cull( bwt, drxn ? seq + seq_ : seq_ + seq, coverage, drxn );
    return culled_;
}

void FillBranch::get( string seq, vector<FillBranch*> path, vector<string>& seqs, vector< vector<FillBranch*> >& paths, bool refresh, bool drxn )
{
    path.insert( path.end(), this );
    bool advanced = false;
    seq = drxn ? seq + seq_ : seq_ + seq;
    for ( FillBranch* fb : branches_ ) if ( !fb->culled_ && ( refresh || fb->added_ ) && ( advanced = true ) ) fb->get( seq, path, seqs, paths, refresh, drxn );
    added_ = false;
    if ( advanced ) return;
    seqs.push_back( seq );
    paths.push_back( path );
}

void FillBranch::get( string seq, vector<string>& seqs, bool drxn )
{
    seq = drxn ? seq + seq_ : seq_ + seq;
    for ( FillBranch* fb : branches_ ) fb->get( seq, seqs, drxn );
    if ( branches_.empty() ) seqs.push_back( seq );
}

//bool FillBranch::cull( Querier& bwt, int coverage, int base, bool drxn )
//{
//    if ( !ext_ )
//    {
//        assert( branches_.size() < 5 );
//        for ( int i = 0; i < branches_.size(); i++ ) if ( branches_[i]->cull( bwt, coverage, ext_, drxn ) )
//        {
//            delete branches_[i];
//            branches_.erase( branches_.begin() + i-- );
//        }
//        if ( branches_.empty() ) return true;
//        merge();
//        return exts_.empty();
//    }
////    if ( setCounts() < max( 2, coverage / 1000 ) ) return true;
//    string seq = drxn ? seq_.substr( base+1, params.readLen ) : seq_.substr( ext_-base-1, params.readLen );
//    int diff = params.readLen / 10;
//    int ol = params.readLen - diff;
//    int exts = bwt.getExtendable( seq, ol, drxn );
//    if ( exts < 2 ) return true;
//    int cutoff = ( coverage * diff / params.readLen ) / 5;
//    if ( exts < cutoff ) return true;
//    
//    for ( int i = 0; i < branches_.size(); i++ ) if ( branches_[i]->cull( bwt, coverage, ext_, drxn ) )
//    {
//        delete branches_[i];
//        branches_.erase( branches_.begin() + i-- );
//    }
//    
//    merge();
//    
//    return exts_.empty();
//}

void FillBranch::merge( bool drxn )
{
    for ( int i = 0; i < branches_.size(); i++ ) if ( branches_[i]->culled_ )
    {
        delete branches_[i];
        branches_.erase( branches_.begin() + i-- );
    }
    if ( branches_.size() == 1 )
    {
        FillBranch* fb = branches_[0];
        if ( fb->good_ ) good_ = seq_.size() + fb->good_;
        seq_ = drxn ? seq_ + fb->seq_ : fb->seq_ + seq_;
        exts_.insert( exts_.end(), fb->exts_.begin(), fb->exts_.end() );
        branches_ = fb->branches_;
        fb->branches_.clear();
        delete fb;
        merge( drxn );
    }
    else for ( FillBranch* fb : branches_ ) fb->merge( drxn );
}

void FillBranch::setCounts()
{
    total_ = unpaired_ = 0;
    for ( FillBranch* fb : branches_ )
    {
        fb->setCounts();
        total_ = max( total_, fb->total_ );
        unpaired_ = max( unpaired_, fb->unpaired_ );
    }
    total_ += exts_.size();
    for ( pair<FillRead, int>& ext : exts_ ) if ( !ext.first.paired ) unpaired_++;
}

FillFork::FillFork( string& seq, string& base, FillRead& fr, int32_t coords[2], bool drxn )
: coord_( coords[drxn] ), used_{ fr.id }
{
    int ol = coords[1] - coords[0], ext = seq.size() - ( coords[1] - coords[0] );
    coords[!drxn] = drxn ? max( 0, coords[1]-params.readLen+1 ) : min( (int)base.size(), coords[0]+params.readLen-1 );
    ol_ = coords[1] - coords[0];
    seq_ = base.substr( coords[0], ol_ );
    branches_.push_back( new FillBranch( drxn ? seq.substr( ol ) : seq.substr( 0, ext ), fr, 0 ) );
}

bool FillFork::add( string& seq, FillRead& fr, int32_t coords[2], bool drxn )
{
    if ( coords[drxn] != coord_ ) return false;
    if ( !used_.insert( fr.id ).second ) return true;
    int ol = coords[1] - coords[0], ext = seq.size() - ( coords[1] - coords[0] );
    for ( FillBranch* fb : branches_ ) if ( fb->add( seq, fr, ol, drxn ) ) return true;
    branches_.push_back( new FillBranch( drxn ? seq.substr( ol ) : seq.substr( 0, ext ), fr, 0 ) );
    return true;
}

void FillFork::add( string& seq, vector<FillBranch*>& path, FillRead& fr, int32_t coords[2], int len, bool drxn )
{
    int ext = drxn ? coords[1] - seq_.size() : len - seq_.size() - coords[0];
    if ( ext <= 0 ) return;
    used_.insert( fr.id );
    assert( ext > 0 );
    for ( int i = 0; i < path.size(); i++ ) if ( i+1 == path.size() || path[i+1]->base_ >= ext )
    {
        int diff = ext - path[i]->base_;
        assert( diff > 0 && path[i]->add( seq, fr, coords[1]-coords[0]-diff, drxn ) );
        for ( int j = 0; j < i; j++ ) path[j]->added_ = true;
        return;
    }
    assert( false );
}

int FillFork::confirm( Querier& bwt, int unpairedCut, int totalCut, bool drxn )
{
    int unpaired = 0, total = 0;
    for ( FillBranch* fb : branches_ ) fb->merge( drxn );
    for ( FillBranch* fb : branches_ ) fb->setCounts();
    FillBranch::confirm( branches_, drxn );
//    for ( FillBranch* fb : branches_ ) confirmed = max( confirmed, fb->confirm( bwt, coverage, drxn) );
    for ( FillBranch* fb : branches_ ) unpaired = max( unpaired, fb->unpaired_ );
    for ( FillBranch* fb : branches_ ) total = max( total, fb->total_ );
    
    if ( total + unpaired > totalCut ) cout << "Confirmed total: " << total << ", unpaired: " << unpaired << " > " << totalCut << endl;
    
    if ( total > 5 )
    {
        int x = 0;
    }
    
    return total + unpaired > totalCut;
}

bool FillFork::cull( Querier& bwt, int coverage, bool drxn )
{
    for ( int i = 0; i < branches_.size(); i++ ) if ( branches_[i]->cull( bwt, seq_, coverage, drxn ) )
    {
        delete branches_[i];
        branches_.erase( branches_.begin() + i-- );
    }
    if ( !branches_.empty() ) return false;
    delete this;
    return true;
}

bool FillFork::get( vector<string>& seqs, vector< vector<FillBranch*> >& paths, bool refresh, bool drxn )
{
    seqs.clear();
    paths.clear();
    for ( FillBranch* fb : branches_ ) if ( !fb->culled_ && ( refresh || fb->added_ ) ) fb->get( seq_, vector<FillBranch*>(), seqs, paths, refresh, drxn );
    assert( seqs.size() == paths.size() );
    return !seqs.empty();
}

void FillFork::get( vector<string>& seqs, bool drxn )
{
    seqs.clear();
    for ( FillBranch* fb : branches_ ) fb->get( seq_, seqs, drxn );
}

FillFork::~FillFork()
{
    for ( FillBranch* fb : branches_ ) delete fb;
}

FillNode::FillNode( Node* node, NodeRoll& nodes, vector<FillNode*>& fills, Nodes& used )
: len( 0 ), bad( node->bad_ )
{
    for ( Node* init = node; node->edges_[0].size() == 1 && node->edges_[0][0].node->edges_[1].size() == 1; )
    {
        node = node->edges_[0][0].node;
        assert( node != init );
    }
    vector<Node*> path;
    for ( Edge* e = NULL; node; )
    {
        if ( e ) len -= e->ol;
        seq_ = e ? seq_ + node->seq_.substr( e->ol ) : node->seq_;
        path_.push_back( make_pair( node, len ) );
        path.push_back( node );
        for ( int d : { 0, 1 } ) for ( NodeMark& nm : node->pe_[d] ) addRead( nodes, node, nm.id, nm.coords[!d] - node->ends_[0] + len, nm.dist, d );
        len += node->size();
        used += node;
        e = node->edges_[1].size() == 1 ? &node->edges_[1][0] : NULL;
        node = e && !e->leap && e->node->edges_[0].size() == 1 ? e->node : NULL;
    }
    fills.push_back( this );
    assert( seq_.size() == len );
    assert( bad ==  path_[0].first->bad_ && bad == path_.back().first->bad_ );
    cover_ = Node::getCoverage( path );
}

bool FillNode::addBranch( string& seq, FillRead& fr, int32_t coords[2], bool drxn )
{
    int ol = coords[1]-coords[0];
    for ( pair<FillNode*, int>& e : edges_[!drxn] ) if ( drxn ? coords[1] <= e.second : len - coords[0] <= e.second ) return false;
    for ( pair<FillNode*, int>& e : edges_[drxn] ) if ( drxn ? coords[1] == len : !coords[0] )
    {
        if ( drxn ? seq[ol] == e.first->seq_[e.second] : seq.end()[-ol-1] == e.first->seq_.end()[-e.second-1] )
        {
            return false;
        }
    }
    for ( FillFork* ff : forks_[drxn] ) if ( ff->add( seq, fr, coords, drxn ) ) return true;
    forks_[drxn].push_back( new FillFork( seq, seq_, fr, coords, drxn ) );
    return true;
}

void FillNode::addRead( NodeRoll& nodes, Node* node, ReadId id, int32_t coord, int32_t dist, bool drxn )
{
    bool added = false;
    Coords* coords;
    for ( Node* tar : nodes.nodes ) if ( ( coords = tar->getRead( id ) ) && ( added = true ) && !coords->ignore )
    {
        bool pairs = false;
        int32_t est = coord + ( ( drxn ? 1 : -1 ) * ( dist + abs( tar->ends_[drxn] - coords->coords[drxn] ) ) );
        for ( FillPaired& fp : paired[drxn] ) if ( fp.node == tar && ( pairs = true ) ) fp.pairs.push_back( FillPair( node, id, est ) );
        if ( !pairs ) paired[drxn].push_back( FillPaired( tar, node, coords, id, est ) );
    }
    reads_[drxn].push_back( FillRead( id, coord + ( drxn ? dist : -dist ), added ) );
}

//void FillNode::addPair( Node* hit, Node* base, ReadId id, int32_t est, bool drxn )
//{
//    for ( FillPaired& fp : paired[drxn] ) if ( fp.node == hit )
//    {
//        fp.pairs.push_back( FillPair( base, id, est ) );
//        return;
//    }
//    paired[drxn].push_back( FillPaired( hit, base, id, est ) );
//}

vector<FillNode*> FillNode::create( NodeRoll& nodes )
{
    vector<FillNode*> fills;
    Nodes used;
    for ( Node* node : nodes.nodes ) if ( node->drxn_ == 2 && !used.find( node ) )
    {
        FillNode* fn = new FillNode( node, nodes, fills, used );
        for ( int d : { 0, 1 } ) fn->edge( nodes, fills, used, d );
    }
    
    return fills;
}

void FillNode::edge( NodeRoll& nodes, vector<FillNode*>& fills, Nodes& used, bool drxn )
{
    if ( bad || !getFork( drxn )->verified_ ) return;
    
    vector<FillNode*> added;
    for ( Edge& e : getFork( drxn )->edges_[drxn] )
    {
        FillNode* edged = NULL;
        bool ignore = false;
        for ( pair<FillNode*, int>& fe : edges_[drxn] ) if ( fe.first->getFork( !drxn ) == e.node ) ignore = true;
        if ( ignore ) continue;
        for ( FillNode* fn : fills ) if ( fn->getFork( !drxn ) == e.node && ( edged = fn ) ) break;
        assert( used.find( e.node ) == bool( edged ) );
        if ( !edged )
        {
            edged = new FillNode( e.node, nodes, fills, used );
            added.push_back( edged );
        }
        edges_[drxn].push_back( make_pair( edged, e.ol ) );
        edged->edges_[!drxn].push_back( make_pair( this, e.ol ) );
    }
    
    for ( int d : { 0, 1 } ) for ( FillNode* fn : added ) fn->edge( nodes, fills, used, d );
}

void FillNode::getDist( unordered_map<FillNode*, int32_t>& dists, int32_t dist, int32_t limit, bool drxn )
{
    auto it = dists.insert( make_pair( this, dist ) );
    if ( !it.second )
    {
        if ( drxn ? it.first->second <= dist : dist <= it.first->second ) return;
        it.first->second = dist;
    }
    
    if ( drxn ? dist < limit : limit < dist ) for ( pair<FillNode*, int>& edge : edges_[drxn] )
    {
        edge.first->getDist( dists, dist + ( ( drxn ? 1 : -1 ) * ( edge.first->len - edge.second ) ), limit, drxn );
    }
}

void FillNode::getDist( unordered_map<FillNode*, int32_t>& dists, int32_t dist, int32_t limit, int32_t& best, bool drxn )
{
    auto it = dists.insert( make_pair( this, dist ) );
    if ( !it.second )
    {
        if ( drxn ? it.first->second <= dist : dist <= it.first->second ) return;
        it.first->second = dist;
    }
    
    if ( drxn ? dist < limit : limit < dist ) for ( pair<FillNode*, int>& edge : edges_[drxn] )
    {
        edge.first->getDist( dists, dist + ( ( drxn ? 1 : -1 ) * ( edge.first->len - edge.second ) ), limit, best, drxn );
    }
    if ( !bad && edges_[drxn].empty() ) best = drxn ? min( best, dist ) : max( best, dist );
}

Node* FillNode::getFork( bool drxn )
{
    return ( drxn ? path_.back() : path_[0] ).first;
}

bool FillNode::isBlunt( int readLimit, bool drxn )
{
    return getFork( !drxn )->isBlunt( 0, readLimit, drxn );
}

bool FillNode::isStrong( int32_t dist, int32_t limit, bool drxn )
{
    for ( pair<Node*, int32_t>& nd : path_ ) if ( nd.first->verified_ )
    {
        if ( dist + ( drxn ? nd.first->size() + nd.second : len - nd.second ) > limit ) return true;
    }
    for ( pair<FillNode*, int>& e : edges_[drxn] ) if ( e.first->isStrong( dist + len - e.second, limit, drxn ) ) return true;
    return false;
}

int FillNode::trim( bool drxn )
{
    if ( bad || !getFork( drxn )->verified_ ) return false;
    
    vector<FillNode*> strong, blunt;
    for ( pair<FillNode*, int>& e : edges_[drxn] )
    {
        if ( e.first->isBlunt( 3, drxn ) ) blunt.push_back( e.first );
        else if ( e.first->isStrong( -e.second, 500, drxn ) ) strong.push_back( e.first );
    }
    
    if ( strong.empty() || blunt.empty() ) return false;
    
    for ( FillNode* fn : blunt )
    {
        getFork( drxn )->removeEdge( fn->getFork( !drxn ), drxn, true );
        fn->getFork( !drxn )->setState();
        for ( int i = 0; i < edges_[drxn].size(); i++ ) if ( edges_[drxn][i].first == fn ) edges_[drxn].erase( edges_[drxn].begin() + i-- );
        for ( int i = 0; i < fn->edges_[drxn].size(); i++ ) if ( fn->edges_[drxn][i].first == this ) fn->edges_[drxn].erase( fn->edges_[drxn].begin() + i-- );
    }
    return blunt.size();
}

LocusFill::LocusFill( Querier& bwt, NodeRoll& nodes )
: fills_( FillNode::create( nodes ) )
{ }

LocusFill::~LocusFill()
{
    clear();
}

void LocusFill::clear()
{
    for ( FillNode* fn : fills_ ) delete fn;
    fills_.clear();
}

int LocusFill::extend( Querier& bwt, NodeRoll& nodes, FillNode* fn, bool drxn )
{
    for ( int i = 0; i < fn->forks_[drxn].size(); i++ ) if ( fn->forks_[drxn][i]->cull( bwt, fn->cover_, drxn ) ) fn->forks_[drxn].erase( fn->forks_[drxn].begin() + i-- );
    
    if ( fn->forks_[drxn].empty() ) return 0;
    
    int32_t limit = drxn ? fn->len : 0;
    for ( FillFork* ff : fn->forks_[drxn] ) limit = drxn ? min( limit, ff->coord_ ) : max( limit, ff->coord_ );
    unordered_map<FillNode*, int32_t> dists;
    fn->getDist( dists, 0, limit + drxn ? -params.maxPeMean : params.maxPeMean, !drxn );
    
    int counted = 0;
    for ( FillFork* ff : fn->forks_[drxn] )
    {
        vector<string> seqs;
        vector< vector<FillBranch*> > paths;
        int32_t coords[2];
        while ( ff->get( seqs, paths, false, drxn ) )
        {
            for ( int i = 0; i < seqs.size(); i++ )
            {
                for ( auto& dist : dists ) for ( FillRead& fr : dist.first->reads_[drxn] ) if ( ff->used_.find( fr.id ) == ff->used_.end() )
                {
                    string q = getSeq( bwt, fr.id );
                    if ( mapSeqEnd( q, seqs[i], params.readLen / 3, coords, !drxn ) ) ff->add( q, paths[i], fr, coords, seqs[i].size(), drxn );
                }
            }
            assert( !ff->cull( bwt, fn->cover_, drxn ) );
        }
        if ( !ff->confirm( bwt, 3, 5 + ( fn->cover_ / 10 ), drxn ) ) continue;
        if ( fn->path_[0].first->drxn_ == 2 && !drxn )
        {
            int x = 0;
            continue;
        }
        ff->get( seqs, paths, true, drxn );
        for ( int i = 0; i < seqs.size(); i++ )
        {
            cout << "Seeding " << ( drxn ? " RIGHT: " : " LEFT: ") << ff->coord_ + fn->path_[0].first->ends_[0] << endl;
            SeedNode::seed( bwt, nodes, seqs[i], drxn ? ff->ol_ : 0, drxn ? 0 : ff->ol_, drxn );
            counted++;
        }
    }
    
    return counted;
}

bool LocusFill::fill( Querier& bwt, NodeRoll& nodes )
{
    LocusFill lf( bwt, nodes );
    if ( lf.seed( bwt, nodes ) ) return true;
    return false;
//    if ( lf.seed( bwt, nodes ) ) lf.refresh( bwt, nodes );
//    if ( lf.join( bwt, nodes ) ) lf.refresh( bwt, nodes );
//    if ( lf.trim( nodes ) ) lf.refresh( bwt, nodes );
//    Nodes repair;
//    for ( int d : { 0, 1 } ) for ( FillNode* fn : lf.fills_ ) lf.map( bwt, nodes, repair, fn, d );
//    for ( Node* node : repair.nodes ) node->setVerified();
//    return false;
}

bool LocusFill::getDists( FillNode* fn, unordered_map<FillNode*, int32_t>& dists, int32_t& ext, int32_t& limit, bool drxn )
{
    limit = drxn ? 0 : fn->len;
    for ( FillRead& fr : fn->reads_[drxn] ) if ( !fr.paired ) limit = drxn ? max( limit, fr.coord ) : min( limit, fr.coord );
    if ( limit == ( drxn ? 0 : fn->len ) ) return false;
    
    ext = drxn ? limit - fn->len + 300 : limit - 300;
    fn->getDist( dists, 0, ext, ext, drxn );
    limit = drxn ? fn->len + ext - 300 : ext + 300;
    return true;
}

//void LocusFill::getSeeds( Querier& bwt, FillNode* fn, vector<FillRead>& seeds, bool drxn )
//{
//    int32_t limit, ext;
//    unordered_map<FillNode*, int32_t> dists;
//    if ( fn->bad || !getDists( fn, dists, ext, limit, drxn ) ) return;
//    
//    vector<FillRead> reads;
//    
//    for ( const pair<FillNode*, int32_t>& dist : dists ) for ( pair<ReadId, int32_t>& read : dist.first->reads_[drxn] )
//    {
//        int32_t coord = read.second + dist.second;
//        bool unseeded = dist.first == fn && seeded_.find( read.first ) == seeded_.end();
//        if ( drxn ? coord < ext : ext < coord ) reads.push_back( FillRead( getSeq( bwt, read.first ), read.first, coord, unseeded ) );
//    }
//    
//    sort( reads.begin(), reads.end(), []( FillRead& a, FillRead& b ){ return a.coord < b.coord; } );
//    int x = 0;
//    for ( FillRead& fr : reads ) if ( fr.unseeded && drxn ? fr.coord < limit : limit < fr.coord )
//    {
//        int32_t coords[2];
//        for ( const pair<FillNode*, int32_t>& dist : dists )
//        {
//            if ( mapSeqEnd( fr.seq, dist.first->seq_, max( 50, fr.mapped ), coords, !drxn ) ) fr.mapped = coords[1] - coords[0];
//        }
//        if ( !fr.mapped || !bwt.isExtendable( fr.seq, drxn ) ) continue;
//        for ( FillRead& alt : reads ) if ( fr.id != alt.id && fr.seq != alt.seq && fr.coord - 300 < alt.coord )
//        {
//            if ( mapSeqOverlap( fr.seq, alt.seq, max( 50, (int)fr.seq.size() - fr.mapped ), drxn ) )
//            {
//                seeds.push_back( fr );
//                seeds.back().coord += fn->path_[0].first->ends_[0] + ( drxn ? -fr.seq.size() : 0 );
//                cout << seeds.back().coord << endl << fr.seq << endl;
//                break;
//            }
//            if ( fr.coord + 300 < alt.coord ) break;
//        }
//    }
//}

string LocusFill::getSeq( Querier& bwt, ReadId id )
{
    auto it = reads_.find( id );
    if ( it != reads_.end() ) return it->second;
    auto ins = reads_.insert( make_pair( id, bwt.getSequence( id ) ) );
    return ins.first->second;
}

bool LocusFill::join( Querier& bwt, NodeRoll& nodes )
{
    return false;
}

//void LocusFill::map( Querier& bwt, NodeRoll& nodes, Nodes& repair, FillNode* fn, bool drxn )
//{
//    int32_t limit, ext;
//    unordered_map<FillNode*, int32_t> dists;
//    if ( fn->bad || !getDists( fn, dists, ext, limit, drxn ) ) return;
//    
//    vector<ReadId> repaired;
//    for ( pair<ReadId, int32_t>& read : fn->reads_[drxn] ) if ( drxn ? read.second < limit : limit < read.second )
//    {
//        string seq = getSeq( bwt, read.first );
//        bool ext[2]{ bwt.isExtendable( seq, 0 ), bwt.isExtendable( seq, 1 ) };
//        if ( ext[0] && ext[1] ) continue;
//        
//        int32_t coords[2];
//        vector<FillHit> hits[2];
//        for ( int d : { 0, 1 } ) for ( const pair<FillNode*, int32_t>& dist : dists ) if ( mapSeqEnd( seq, dist.first->seq_, 50, coords, d ) )
//        {
//            FillHit::add( hits[d], dist.first, seq, coords, read.second, dist.second, drxn );
//        }
//        
//        bool ignore = false;
//        if ( FillHit::confirm( hits, ext, seq.size(), ignore ) )
//        {
//            for ( int d : { 0, 1 } ) for ( FillHit& fh : hits[d] ) fh.set( read.first, ignore );
//            repaired.push_back( params.getPairId( read.first ) );
//            mapped_.insert( read.first );
//        }
//        else if ( !ext[0] && !ext[1] )
//        {
//            ReadId id = params.getPairId( read.first );
//            for ( int i = 0; i < fn->path_.size(); i++ )
//            {
//                auto it = fn->path_[i].first->reads_.find( id );
//                if ( it == fn->path_[i].first->reads_.end() ) continue;
//                it->second.unpaired = true;
//                assert( fn->path_[i].first->rmvMark( read.first, drxn ) );
//            }
//        }
//        else ( unmapped_ ).insert( read.first );
//    }
//    for ( pair<Node*, int32_t>& nd : fn->path_ ) for ( ReadId id : repaired ) if ( nd.first->reads_.find( id ) != nd.first->reads_.end() )
//    {
//        repair += nd.first;
//        break;
//    }
//}

void LocusFill::refresh( Querier& bwt, NodeRoll& nodes )
{
    clear();
    Node::prunePaths( bwt, nodes );
    fills_ = FillNode::create( nodes );
}

bool LocusFill::seed( Querier& bwt, NodeRoll& nodes )
{
    unordered_set<ReadId> added;
    int counted = 0;
    for ( FillNode* fn : fills_ ) for ( int d : { 0, 1 } ) seed( bwt, nodes, fn, added, d );
    for ( FillNode* fn : fills_ ) for ( int d : { 0, 1 } ) counted += extend( bwt, nodes, fn, d );
    return counted;
}

void LocusFill::seed( Querier& bwt, NodeRoll& nodes, FillNode* fn, unordered_set<ReadId>& added, bool drxn )
{
    if ( fn->reads_[drxn].empty() ) return;
    int32_t limit = drxn ? 0 : fn->len;
    for ( FillRead& fr : fn->reads_[drxn] ) limit = drxn ? max( limit, fr.coord ) : min( limit, fr.coord );
    
    unordered_map<FillNode*, int32_t> dists;
    fn->getDist( dists, 0, limit + drxn ? params.readLen : -params.readLen, drxn );
    
    int minOl = params.readLen * .9;
    string seq;
    for ( FillRead& fr : fn->reads_[drxn] ) if ( ( seq = getSeq( bwt, fr.id ) ).size() == params.readLen )
    {
        int32_t coords[2];
        if ( seq.size() == params.readLen ) for ( const pair<FillNode*, int32_t>& dist : dists ) if ( mapSeqEnd( seq, dist.first->seq_, minOl, coords, !drxn ) )
        {
            if ( coords[1] - coords[0] < seq.size() ) dist.first->addBranch( seq, fr, coords, drxn );
        }
    }
}

bool LocusFill::trim( NodeRoll& nodes )
{
    int trimmed = 0;
    for ( int d : { 0, 1 } ) for ( FillNode* fn : fills_ ) trimmed += fn->trim( d );
    return trimmed;
}

