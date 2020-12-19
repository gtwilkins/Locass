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

#include "kmer_graph.h"
#include "constants.h"
#include "parameters.h"
#include "shared_functions.h"
#include <cassert>
#include <algorithm>
#include <iostream>

//void KmerNode::add( KmerRead* read, int coord )
//{
//    assert( !read->node_ );
//    read->node_ = this;
//    reads_.push_back( make_pair( read, coord ) );
//}

void KmerNode::addEdge( KmerNode* node, int ol, bool drxn )
{
    for ( pair<KmerNode*, int>& edge : edges_[drxn] ) if ( edge.first == node ) return;
    edges_[drxn].push_back( make_pair( node, ol ) );
    node->edges_[!drxn].push_back( make_pair( this, ol) );
}

void KmerNode::get( unordered_map<KmerNode*, int>& got, int dist, bool drxn )
{
    if ( got.insert( make_pair( this, dist ) ).second ) for ( pair<KmerNode*, int>& edge : edges_[drxn] )
    {
        edge.first->get( got, dist+seq_.size()-edge.second, drxn );
    }
}

int KmerNode::getTerminalLen( bool drxn )
{
    return ( reads_.empty() ? seq_: ( drxn ? reads_.back() : reads_[0] ).first->seq_ ).size();
}

bool KmerNode::isRedundant( KmerRead* read )
{
    for ( pair<KmerRead*, int>& kr : redundant_ ) if ( kr.first == read ) return true;
    return false;
}

int KmerExt::match( KmerExt& ext, bool exact, bool drxn )
{
    int i = 0;
    for ( ; i < min( ext.ext_, ext_ ); i++ )
    {
        if ( drxn ? ext.read_->seq_[ext.ol_+i] != read_->seq_[ol_+i] : ext.read_->seq_[ext.ext_-i-1] != read_->seq_[ext_-i-1] ) return exact ? 0 : i;
    }
    return i;
}

int KmerExt::match( KmerExt& ext, bool drxn )
{
    int i = 0;
    for ( ; i < min( ext.ext_, ext_ ); i++ )
    {
        if ( drxn ? ext.read_->seq_[ext.ol_+i] != read_->seq_[ol_+i] : ext.read_->seq_[ext.ext_-i-1] != read_->seq_[ext_-i-1] ) return i;
    }
    return i;
}

void KmerExt::operator += ( int ol )
{
    ext_ -= ol;
    ol_ += ol;
}

void KmerOverlaps::addExt( KmerRead* read, int ol, bool drxn )
{
    KmerExt ext( read, ol );
    int i = 0;
    while ( i < exts_[drxn].size() && ( exts_[drxn][i].ol_ == ext.ol_ ? exts_[drxn][i].ext_ <= ext.ext_ : ext.ol_ < exts_[drxn][i].ol_ ) )
    {
        assert( exts_[drxn][i].read_ != ext.read_ );
        i++;
    }
    exts_[drxn].insert( exts_[drxn].begin() + i, ext );
}

bool KmerOverlaps::addOverlap( KmerRead* read, int base, int ol, bool drxn )
{
    assert( node_ );
    int limit = min( (int)read->seq_.size(), drxn ? (int)node_->seq_.size() - base : base );
    while ( ol < limit && ( drxn ? read->seq_[ol] == node_->seq_[base+ol] : read->seq_.end()[-ol-1] == node_->seq_[base-ol-1] ) ) ol++;
    
    if ( ol == read->seq_.size() ) return false;
    else if ( ol == limit ) addExt( read, ol, drxn );
    else splints_[drxn].push_back( make_pair( KmerExt( read, ol ), drxn ? base + ol : node_->seq_.size() - base + ol ) );
    return true;
}

void KmerOverlaps::test( string& t, int ext, bool drxn )
{
    int coords[2];
    for ( int d : { 0, 1 } ) for ( pair<KmerExt, int>& ke : splints_[d] )
    {
        assert( mapSeqEnd( ke.first.read_->seq_, t, ke.first.ol_, coords, !d ) );
        assert( coords[1] - coords[0] == ke.first.ol_ );
        assert( d ? coords[1] == ke.second : coords[0] == t.size() - ke.second );
    }
    for ( KmerExt& ke : exts_[drxn] ) 
    {
        int ol = ke.ol_ + ext;
        if ( ke.ext_ >= ext ) assert( mapSeqOverlap( t, ke.read_->seq_, ol, drxn ) == ol );
        else assert( t.find( ke.read_->seq_ ) != string::npos );
    }
    for ( KmerExt& ke : exts_[!drxn] ) assert( mapSeqOverlap( t, ke.read_->seq_, ke.ol_, !drxn ) == ke.ol_ );
}

KmerGraph::KmerGraph( int len )
: len_( len ), mask_( 0 )
{
    assert( len_ <= 32 && len_ > 1 );
    for ( int i = 0; i < len_; i++ ) mask_ = ( mask_ << 2 ) + 3;
}

void KmerGraph::add( ReadId id, string seq )
{
    KmerRead* kr = new KmerRead( seq, id );
    reads_.push_back( kr );
    CharId k;
    Char c[2];
    Kmer* kmer[2]{ NULL, NULL };
    int coord = -len_;
    for ( int i = 0; i < seq.size(); i++ )
    {
        c[1] = charToInt[ seq[i] ];
        k = ( ( k << 2 ) + c[1] ) & mask_;
        if ( ++coord < 0 ) continue;
        if ( kmer[1] = getKmer( k ) ) kmer[1]->coords_.push_back( make_pair( kr, coord ) );
        else kmers_.insert( make_pair( k, ( kmer[1] = new Kmer( kr, coord ) ) ) );
        if ( kmer[0] )
        {
            assert( !kmer[0]->edges_[1][ c[1] ] || kmer[0]->edges_[1][ c[1] ] == kmer[1] );
            assert( !kmer[1]->edges_[0][ c[0] ] || kmer[1]->edges_[0][ c[0] ] == kmer[0] );
            kmer[0]->edges_[1][ c[1] ] = kmer[1];
            kmer[1]->edges_[0][ c[0] ] = kmer[0];
        }
        kmer[0] = kmer[1];
        c[0] = charToInt[ seq[coord] ];
    }
}

void KmerGraph::addMapped( KmerNode* node, KmerRead* read, int coord )
{
    bool redundant = false;
    int i = 0, len = read->seq_.size(), r = coord + read->seq_.size();
    for ( ; i < node->reads_.size() && ( node->reads_[i].second == coord ? node->reads_[i].first->seq_.size() > len : node->reads_[i].second < coord ); i++ )
    {
        if ( node->reads_[i].second+node->reads_[i].first->seq_.size() >= r && node->reads_[i].first->seq_.size() > len ) redundant = true;
    }
    while ( i < node->reads_.size() && node->reads_[i].first->seq_.size() < len && node->reads_[i].second+node->reads_[i].first->seq_.size() <= r )
    {
        assert( false );
    }
    if ( read->node_ ) reused_.insert( read );
    if ( read->node_ && !redundant )
    {
        assert( read->node_ != node );
        assert( read->node_->isRedundant( read ) );
    }
    if ( redundant ) node->redundant_.push_back( make_pair( read, coord ) );
    else node->reads_.insert( node->reads_.begin() + i, make_pair( read, coord ) );
    if ( !read->node_ || !redundant ) read->node_ = node;
    assert( node->seq_.find( read->seq_ ) == coord );
}

//void KmerGraph::addSplints( KmerExt ext, KmerOverlaps* branch, KmerOverlaps* base, bool drxn )
//{
//    for ( pair<KmerExt, int>& ke : base->splints_[!drxn] ) if ( ke.second - ke.first.ol_ + len_ <= ext.ol_ )
//    {
//        if ( ke.second < ext.ol_ ) branch->splints_[!drxn].push_back( make_pair( ke.first, ke.second+ext.ext_ ) );
//        else branch->addExt( ke.first.read_, ext.ol_ - ke.second + ke.first.ol_, !drxn );
//    }
//    
//    int unol = base->node_->seq_.size() - ext.ol_;
//    for ( pair<KmerExt, int>& ke : base->splints_[drxn] ) if ( unol <= ke.second - ke.first.ol_ )
//    {
//        branch->splints_[drxn].push_back( make_pair( ke.first, ke.second-unol ) );
//    }
//}

void KmerGraph::addSplints( KmerOverlaps* branch, KmerOverlaps* base, int ext, int ol, bool drxn )
{
    int len = base->node_->seq_.size();
    for ( pair<KmerRead*, int>& kr : base->node_->reads_ )
    {
        int dist = drxn ? len - kr.second - kr.first->seq_.size() : kr.second;
        if ( ol-dist >= len_ ) branch->addExt( kr.first, ol-dist, !drxn );
    }
    for ( pair<KmerRead*, int>& kr : base->node_->redundant_ )
    {
        int dist = drxn ? len - kr.second - kr.first->seq_.size() : kr.second;
        if ( ol-dist >= len_ ) branch->addExt( kr.first, ol-dist, !drxn );
    }
    
    for ( KmerExt& ke : base->exts_[!drxn] ) if ( ol + ke.ol_ - len >= len_ ) branch->addExt( ke.read_, ol + ke.ol_ - len, !drxn );

    for ( pair<KmerExt, int>& ke : base->splints_[!drxn] ) if ( ke.second - ke.first.ol_ + len_ <= ol )
    {
        if ( ke.second < ol ) branch->splints_[!drxn].push_back( make_pair( ke.first, ke.second+ext ) );
        else branch->addExt( ke.first.read_, ol - ke.second + ke.first.ol_, !drxn );
    }

    int unol = len - ol;
    for ( pair<KmerExt, int>& ke : base->splints_[drxn] ) if ( unol <= ke.second - ke.first.ol_ )
    {
        branch->splints_[drxn].push_back( make_pair( ke.first, ke.second-unol ) );
    }
}

void KmerGraph::advance( KmerOverlaps* ols, int ext, bool drxn )
{
    for ( int i = 0; i < ols->exts_[drxn].size(); i++ ) ols->exts_[drxn][i] += ext;
    for ( int i = 0; i < ols->exts_[drxn].size(); i++ ) if ( ols->exts_[drxn][i].ext_ <= 0 )
    {
        int coord = drxn ? ols->node_->seq_.size() - ols->exts_[drxn][i].ol_ : -ols->exts_[drxn][i].ext_;
        addMapped( ols->node_, ols->exts_[drxn][i].read_, coord );
        ols->exts_[drxn].erase( ols->exts_[drxn].begin() + i-- );
    }
}

void KmerGraph::assemble( string seq, bool drxn )
{
    setBase( seq, drxn );
    
//    assert( false );
//    vector<Kmer*> kmers;
//    for ( pair<CharId, Kmer*> k : kmers_ ) kmers.push_back( k.second );
//    kmers_.clear();
//    int assembled = 0, total = kmers.size();
//    for ( Kmer* k : kmers ) if ( k->coords_.size() > 1 ) assembled++;
//    
//    vector<KmerNode*> nodes;
//    for ( Kmer* k : kmers ) if ( k->coords_.size() > 1 )
//    {
//        unordered_map<Kmer*, KmerNode*> forks[2];
//        KmerNode* kn = create( k, getSeq( k ), nodes, forks );
//        sort( nodes.begin(), nodes.end(), []( KmerNode* a, KmerNode* b ){ return a->reads2_.size() > b->reads2_.size(); });
////        for ( KmerCoord& kc : nodes[200]->reads_ )
////        {
////            cout << ">" << kc.id_ << endl << string( 150 - kc.read_ + kc.node_, '-' ) + *getSeq( kc.id_ ) << endl;
////        }
//        assert( false );
//        int bads = 0;
//        for ( KmerNode* node : nodes ) if ( node->reads2_.size() == 1  && ( node->edges2_[0].empty() || node->edges2_[1].empty() ) )
//        {
//            bads++;
//        }
//        assert( false );
//    }
}

void KmerGraph::create( KmerOverlaps* ols )
{
    test( ols );
    KmerNode* node = ols->node_;
    assert( node );
    vector< pair<KmerOverlaps*, int> > branches[2];
    for ( int d : { 0, 1 } ) while ( !getBranches( ols, branches[d], !ols->node_->edges_[d].empty(), d ) )
    {
        KmerExt ext = ols->exts_[d][0];
        KmerOverlaps* fork = setOverlaps( ext.read_->seq_, NULL, ext.ext_, false, d );
        if ( !fork->exts_[!d].empty() )
        {
            vector<KmerExt> exts;
            for ( KmerExt& ke : ols->exts_[d] ) exts.push_back( KmerExt( ke, ext.ext_ ) );
            fork->exts_[d].insert( fork->exts_[d].begin(), exts.begin(), exts.end() );
            branches[d].push_back( make_pair( fork, ext.ol_ ) );
        }
        else
        {
            node->seq_ = d ? node->seq_ + ext.read_->seq_.substr( ext.ol_ ) : ext.read_->seq_.substr( 0, ext.ext_ ) + node->seq_;
            for ( int i = 0; i < ols->splints_[!d].size(); i++ ) ols->splints_[!d][i].second += ext.ext_;
            if ( !d ) for ( int i = 0; i < node->reads_.size(); i++ ) node->reads_[i].second += ext.ext_;
            advance( ols, ext.ext_, d );
            setOverlaps( node->seq_, ols, ext.ext_, true, d );
            delete fork->node_;
            delete fork;
        }
    }
    for ( int d : { 0, 1 } )
    {
        sort( ols->splints_[!d].begin(), ols->splints_[!d].end(), []( pair<KmerExt, int> a, pair<KmerExt, int> b ){
            return a.second == b.second ? a.first.ol_ > b.first.ol_ : a.second < b.second; }
        );
    }
    
    for ( int d : { 0, 1 } ) for ( pair<KmerOverlaps*, int> branch : branches[d] )
    {
        KmerExt ext = branch.first->exts_[d][0];
        addSplints( branch.first, ols, ext.read_->seq_.size()-branch.second, branch.second, d );
        branch.first->test( ext.read_->seq_, ext.ext_, d );
        setBranch( ext.read_, node, branch.first, ext.ext_, branch.second, d );
    }
    
    delete ols;
}

bool KmerGraph::getBranches( KmerOverlaps* ols, vector< pair<KmerOverlaps*, int> >& branches, bool branched, bool drxn )
{
    if ( !branches.empty() || ols->exts_[drxn].empty() ) return true;
    
    vector<KmerExt> exts = ols->exts_[drxn], redundant;
    int base = ols->node_->reads_.empty() ? ols->node_->seq_.size() : ( drxn ? ols->node_->reads_.back() : ols->node_->reads_[0] ).first->seq_.size();
    int maxol = exts[0].ol_;
    
    if ( exts[0].read_->node_ && !exts[0].read_->node_->isRedundant( exts[0].read_ ) ) branched = true;
    if ( !branched ) for ( pair<KmerExt, int>& ke : ols->splints_[!drxn] )
    {
        if ( branched = ( ke.second < base && maxol <= ke.second && ke.second - ke.first.ol_ + len_ <= maxol ) ) break;
    }
    
    while ( !exts.empty() )
    {
        vector<KmerExt> alts;
        KmerOverlaps* branch = new KmerOverlaps();
        KmerExt ext = exts[0];
        branch->exts_[drxn] = { ext };
        
        for ( pair<KmerOverlaps*, int> ko : branches ) for ( KmerExt& ke : ko.first->exts_[drxn] )
            if ( ke.ol_ <= ext.ol_ && ke.ext_ >= ko.first->exts_[drxn][0].ext_ )
        {
            int match = ext.match( ke, drxn );
            branch->splints_[drxn].push_back( make_pair( KmerExt( ke, match ), ext.ol_+match ) );
        }
        for ( KmerExt& ke : redundant ) if ( ke.ol_ <= ext.ol_ && ext.match( ke, drxn ) == ke.ext_ ) branch->exts_[drxn].push_back( ke );
        for ( int i = 1; i < exts.size(); i++ )
        {
            int match = ext.match( exts[i], drxn );
            if ( match == ext.ext_ ) branch->exts_[drxn].push_back( exts[i] );
            else if ( match == exts[i].ext_ )
            {
                branch->exts_[drxn].push_back( exts[i] );
                redundant.push_back( exts[i] );
            }
            else
            {
                branch->splints_[drxn].push_back( make_pair( KmerExt( exts[i], match ), ext.ol_+match ) );
                alts.push_back( exts[i] );
            }
        }
        
        if ( !branched && branches.empty() && alts.empty() )
        {
            delete branch;
            return false;
        }
        
        branches.push_back( make_pair( branch, ext.ol_ ) );
        exts = alts;
    }
    
//    for ( KmerOverlaps* ko : branches ) addSplints( ko->exts_[drxn][0], ko, ols, drxn );
//    for ( KmerOverlaps* ko : branches ) ko->test( ko->exts_[drxn][0].read_->seq_, ko->exts_[drxn][0].ext_, drxn );
    
    return true;
}

//bool KmerGraph::getBranches( KmerOverlaps* ols, vector< pair<KmerOverlaps*, int> >& branches, bool branched, bool drxn )
//{
//    
//    for ( int i = 1; i < ols->exts_[drxn].size(); i++ ) if ( !ols->exts_[drxn][0].match( ols->exts_[drxn][i], true, drxn ) )
//    {
//        int ext =  ols->exts_[drxn][i].ext_, ol = ols->exts_[drxn][i].ol_;
//        KmerOverlaps* alt = new KmerOverlaps();
//        alt->exts_[drxn].push_back( KmerExt( ols->exts_[drxn][i], ext ) );
//        for ( int j = i+1; j < ols->exts_[drxn].size(); j++ ) if ( ols->exts_[drxn][i].match( ols->exts_[drxn][j], true, drxn ) )
//        {
//            alt->exts_[drxn].push_back( KmerExt( ols->exts_[drxn][j], ext ) );
//            ols->exts_[drxn].erase( ols->exts_[drxn].begin() + j-- );
//        }
//        ols->exts_[drxn].erase( ols->exts_[drxn].begin() + i-- );
//        branches.push_back( make_pair( alt, ol ) );
//        branched = true;
//    }
//    
//    if ( !branched ) for ( pair<KmerExt, int>& ke : ols->splints_[!drxn] )
//        if ( branched = isBackOverlap( ke.second, ke.first.ol_, ols->exts_[drxn][0].ol_ ) ) break;
//    
//    // No branched means valid extension
//    if ( !branched ) return false;
//    
//    // Convert potential extension to branch
//    branches.push_back( make_pair( new KmerOverlaps(), ols->exts_[drxn][0].ol_ ) );
//    for ( int i = 0; i < ols->exts_[drxn].size(); i++ ) branches.back().first->exts_[drxn].push_back( KmerExt( ols->exts_[drxn][i], ols->exts_[drxn][0].ext_ ) );
//    ols->exts_[drxn].clear();
//    
//    // Set back overlaps from back splinters
//    int base = ols->node_->reads_.empty() ? ols->node_->seq_.size() : ( drxn ? ols->node_->reads_.back() : ols->node_->reads_[0] ).first->seq_.size();
//    for ( int d : { 0, 1 } ) for ( pair<KmerOverlaps*, int>& ko : branches )
//    {
//        int ext = ko.first->exts_[drxn][0].read_->seq_.size() - ko.second, unol = ols->node_->seq_.size() - ko.second;
//        if ( d != drxn ) for ( pair<KmerExt, int>& ke : ols->splints_[d] ) if ( ke.second < base && ke.second - ke.first.ol_ + len_ <= ko.second )
//        {
//            if ( isBackOverlap( ke.second, ke.first.ol_, ko.second ) )
//            {
//                ko.first->addExt( ke.first.read_, ko.second - ke.second + ke.first.ol_, d );
//                string l = drxn ? ke.first.read_->seq_ : ko.first->exts_[drxn][0].read_->seq_;
//                string r = drxn ? ko.first->exts_[drxn][0].read_->seq_ : ke.first.read_->seq_;
//                int ol = ko.second - ke.second + ke.first.ol_;
//                assert( l.substr( l.size() - ol ) == r.substr( 0, ol ) );
//            }
//            else
//            {
//                ko.first->splints_[d].push_back( make_pair( ke.first, ke.second+ext ) );
//                string q = ke.first.read_->seq_;
//                q = drxn ? q.substr( q.size() - ke.first.ol_ ) : q.substr( 0, ke.first.ol_ );
//                string t = ko.first->exts_[drxn][0].read_->seq_;
//                int coord = ke.second+ext;
//                assert( t.find( q ) == ( drxn ? t.size() - coord : coord - q.size() ) );
//            }
//            
//        }
//        
//        if ( d == drxn ) for ( pair<KmerExt, int>& ke : ols->splints_[d] ) if ( unol <= ke.second - ke.first.ol_ )
//        {
//            ko.first->splints_[d].push_back( make_pair( ke.first, ke.second-unol ) );
//            string q = ke.first.read_->seq_.substr( d ? 0 : ke.first.read_->seq_.size() - ke.first.ol_, ke.first.ol_ ) ;
//            string t = ko.first->exts_[drxn][0].read_->seq_;
//            assert( t.find( q ) == ( d ? ke.second-unol-ke.first.ol_ : t.size()-(ke.second-unol) ) );
//        }
//        
//        if ( d == drxn ) for ( pair<KmerOverlaps*, int>& alt : branches ) if ( alt.first != ko.first ) for ( KmerExt& ke : alt.first->exts_[d] )
//        {
//            if ( ke.ol_ <= ko.second )
//            {
//                int matched = ke.match( ko.first->exts_[d][0], false, d );
//                ko.first->splints_[d].push_back( make_pair( KmerExt( ke.read_, ke.ext_-matched ), ko.second+matched ) );
//                int ol = ke.ol_ + matched;
//                string q = d ? ke.read_->seq_.substr( 0, ol ) : ke.read_->seq_.substr( ke.read_->seq_.size()-ol );
//                string t = ko.first->exts_[drxn][0].read_->seq_;
//                assert( t.find( q ) == d ? ko.second - ke.ol_ : ko.first->exts_[drxn][0].ext_-matched );
//            }
//        }
//    }
//    
//    return true;
//}

CharId KmerGraph::getK( string& seq, int i )
{
    CharId k = 0;
    for ( int j = 0; j < len_; j++ ) k = ( k << 2 ) + charToInt[ seq[i+j] ];
    return k;
}

Kmer* KmerGraph::getKmer( CharId k )
{
    auto it = kmers_.find( k );
    return it != kmers_.end() ? it->second : NULL;
}

string KmerGraph::getSeq( Kmer* k )
{
    if ( !k->coords_.empty() ) return k->coords_[0].first->seq_.substr(  k->coords_[0].second, len_);
    return "";
}

bool KmerGraph::isBackOverlap( int dist, int ol, int cut )
{
    return cut <= dist && dist - ol + len_ <= cut;
}

void KmerGraph::setBase( string seq, bool drxn )
{
    KmerOverlaps* ols = setOverlaps( seq );
    KmerNode* node = ols->node_;
    vector< pair<KmerOverlaps*, int> > branches;
    getBranches( ols, branches, true, drxn );
    for ( pair<KmerOverlaps*, int> branch : branches )
    {
        KmerExt ext = branch.first->exts_[drxn][0];
        addSplints( branch.first, ols, ext.ext_, ext.ol_, drxn );
        setBranch( ext.read_, ols->node_, branch.first, ext.ext_, ext.ol_, drxn );
    }
    delete ols;
    unordered_map<KmerNode*, int> got;
    node->get( got, 0, drxn );
    vector< pair<KmerNode*, int> > nodes( got.begin(), got.end() );
    sort( nodes.begin(), nodes.end(), []( pair<KmerNode*, int> a, pair<KmerNode*, int> b ){ return a.second < b.second; } );
    for ( int i = 0; i < nodes.size(); i++ ) cout << ">" << i << endl << string( nodes[i].second, '-' ) << nodes[i].first->seq_ << endl;
    assert( false );
}

void KmerGraph::setBranch( KmerRead* read, KmerNode* node, KmerOverlaps* branch, int ext, int ol, bool drxn )
{
    if ( read->node_ && !read->node_->isRedundant( read ) )
    {
        size_t it = read->node_->seq_.find( read->seq_ );
        assert( drxn ? !it : it == read->node_->seq_.size() - read->seq_.size() );
        if ( branch->node_ ) delete branch->node_;
        node->addEdge( read->node_, ol, drxn );
        delete branch;
        return;
    }
    if ( !branch->node_ ) branch->node_ = new KmerNode( read->seq_ );
    advance( branch, ext, drxn );
    if ( ext ) setOverlaps( branch->node_->seq_, branch, ext, true, drxn );
    node->addEdge( branch->node_, ol, drxn );
    create( branch );
}

//bool KmerGraph::setExtend( KmerExt& ext, KmerOverlaps* ols, KmerOverlaps*& fork, bool drxn )
//{
//    fork = setOverlaps( ext.read_->seq_, NULL, ext.ext_, false, drxn );
//    if ( !fork->exts_[!drxn].empty() )
//    {
//        addSplints( ext, fork, ols, drxn );
//        for ( KmerExt& ke : ols->exts_[drxn] ) ke += ext.ext_;
//        fork->exts_[drxn].insert( fork->exts_[drxn].begin(), ols->exts_[drxn].begin(), ols->exts_[drxn].end() );
//        fork->test( ext.read_->seq_, 0, drxn );
//        return false;
//    }
//    delete fork->node_;
//    fork = NULL;
//    
////    node->seq_ = drxn ? node->seq_ + ext.read_->seq_.substr( ext.ol_ ) : ext.read_->seq_.substr( 0, ext.ext_ ) + node->seq_;
////    for ( int i = 0; i < ols->splints_[!drxn].size(); i++ ) ols->splints_[!drxn][i].second += ext.ext_;
////    if ( !drxn ) for ( int i = 0; i < node->reads_.size(); i++ ) node->reads_[i].second += ext.ext_;
////    advance( ols, ext.ext_, drxn );
//    
//    return true;
//}

KmerOverlaps* KmerGraph::setOverlaps( string seq, KmerOverlaps* ols, int ext,bool map, bool drxn )
{
    if ( !ext ) ext = seq.size() + 1 - len_;
    if ( !ols ) ols = new KmerOverlaps();
    if ( !ols->node_ ) ols->node_ = new KmerNode( seq );
    int coord = drxn ? seq.size() - len_ - ext : ext;
    Kmer* k = NULL;
    for ( int i = 0; i < ext; i++ )
    {
        coord += drxn ? 1 : -1;
        if ( k ) k = k->edges_[drxn][ charToInt[ seq[ drxn ? coord+len_-1 : coord ] ] ];
        if ( k ) assert( k == getKmer( getK( seq, coord ) ) );
        if ( !k ) k = getKmer( getK( seq, coord ) );
        if ( k ) for ( pair<KmerRead*, int>& kr : k->coords_ ) if ( !kr.second || kr.second+len_ == kr.first->seq_.size() )
        {
            if ( !ols->addOverlap( kr.first, kr.second ? coord+len_ : coord, len_, !kr.second ) && ( !kr.second == drxn ) )
            {
                if ( map ) addMapped( ols->node_, kr.first, kr.second ? coord + len_ - kr.first->seq_.size() : coord );
            }
        }
    }
    return ols;
}

void KmerGraph::test( KmerOverlaps* ols )
{
    unordered_set<KmerRead*> base[2], rebase[2];
    for ( int d : { 0, 1 } ) for ( KmerExt& ke : ols->exts_[d] ) base[d].insert( ke.read_ );
    
    string t = ols->node_->seq_;
    Kmer* k[2]{ getKmer( getK( t, 0 ) ), getKmer( getK( t, t.size()-len_ ) ) };
    if ( k[0] ) for( pair<KmerRead*, int>& kr : k[0]->coords_ )
    {
        int ol = len_, limit = kr.first->seq_.size()-kr.second;
        if ( limit > t.size() ) continue;
        while ( ol < limit && t[ol] == kr.first->seq_[ol+kr.second] ) ol++;
        if ( ol < limit || ol == kr.first->seq_.size() ) continue;
        assert( base[0].find( kr.first ) != base[0].end() );
        rebase[0].insert( kr.first );
    }
    if ( k[1] )for( pair<KmerRead*, int>& kr : k[1]->coords_ )
    {
        int ol = len_, limit = kr.second+len_;
        if ( limit > t.size() ) continue;
        while ( ol < limit && t.end()[-ol-1] == kr.first->seq_[kr.second+len_-ol-1] ) ol++;
        if ( ol < limit || ol == kr.first->seq_.size() ) continue;
        assert( base[1].find( kr.first ) != base[1].end() );
        rebase[1].insert( kr.first );
    }
    for ( int d : { 0, 1 } ) assert( base[d].size() == rebase[d].size() );
}

//void KmerGraph::setBranch( KmerNode* node, KmerOverlaps* branch, int ol, bool drxn )
//{
//    if ( branch->exts_[drxn][0].read_->node_ && !branch->exts_[drxn][0].read_->node_->isRedundant( branch->exts_[drxn][0].read_ ) )
//    {
//        size_t it = branch->exts_[drxn][0].read_->node_->seq_.find( branch->exts_[drxn][0].read_->seq_ );
//        assert( drxn ? !it : it == branch->exts_[drxn][0].read_->node_->seq_.size() - branch->exts_[drxn][0].read_->seq_.size() );
//        node->addEdge( branch->exts_[drxn][0].read_->node_, ol, drxn );
//        delete branch;
//        return;
//    }
//    assert( !branch->node_ );
//    if ( !branch->node_ ) branch->node_ = new KmerNode( branch->exts_[drxn][0].read_->seq_ );
//    while ( !branch->exts_[drxn].empty() && !branch->exts_[drxn][0].ext_ )
//    {
//        addMapped( branch->node_, branch->exts_[drxn][0].read_, drxn ? branch->node_->seq_.size() - branch->exts_[drxn][0].ol_ : 0 );
//        branch->exts_[drxn].erase( branch->exts_[drxn].begin() );
//    }
//    setOverlaps( branch->node_->seq_, branch, branch->node_->seq_.size()-ol, drxn );
//    node->addEdge( branch->node_, ol, drxn );
//    create( branch );
//}
