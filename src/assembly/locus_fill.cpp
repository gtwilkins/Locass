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
    for ( int i = 0; !hit && i < fn->path.size(); i++ ) if ( drxn ? coords[1] <= fn->path[i].second + fn->path[i].first->size()
                                                                  : ( i+1 == fn->path.size() || coords[0] < fn->path[i+1].second ) )
    {
        hit = fn->path[i].first;
        if ( drxn ) coords[0] = max( coords[0], fn->path[i].second );
        if ( !drxn ) coords[1] = min( coords[1], fn->path[i].second + fn->path[i].first->size() );
        for ( int d : { 0, 1 } ) coords[d] += hit->ends_[0] - fn->path[i].second;
        assert( coords[0] < coords[1] );
        assert( fn->path[i].first->ends_[0] <= coords[0] && coords[1] <= fn->path[i].first->ends_[1] );
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

FillNode::FillNode( Node* node, NodeRoll& nodes, vector<FillNode*>& fills, Nodes& used )
: len( 0 ), bad( node->bad_ )
{
    for ( Node* init = node; node->edges_[0].size() == 1 && node->edges_[0][0].node->edges_[1].size() == 1; )
    {
        node = node->edges_[0][0].node;
        assert( node != init );
    }
    for ( Edge* e = NULL; node; )
    {
        if ( e ) len -= e->ol;
        seq = e ? seq + node->seq_.substr( e->ol ) : node->seq_;
        path.push_back( make_pair( node, len ) );
        for ( int d : { 0, 1 } ) for ( NodeMark& nm : node->pe_[d] ) addRead( nodes, node, nm.id, nm.coords[!d] - node->ends_[0] + len, nm.dist, d );
        len += node->size();
        used += node;
        e = node->edges_[1].size() == 1 ? &node->edges_[1][0] : NULL;
        node = e && !e->leap && e->node->edges_[0].size() == 1 ? e->node : NULL;
    }
    fills.push_back( this );
    assert( seq.size() == len );
    assert( bad ==  path[0].first->bad_ && bad == path.back().first->bad_ );
}

void FillNode::addRead( NodeRoll& nodes, Node* node, ReadId id, int32_t coord, int32_t dist, bool drxn )
{
    bool added = false;
    for ( Node* tar : nodes.nodes )
    {
        auto it = tar->reads_.find( id );
        if ( it != tar->reads_.end() && ( added = true ) )
        {
            if ( it->second.ignore ) return;
            addPair( tar, node, id, coord + ( ( drxn ? 1 : -1 ) * ( dist + abs( tar->ends_[drxn] - it->second[drxn] ) ) ), drxn );
        }
    }
    if ( !added ) reads[drxn].push_back( make_pair( id, coord + ( drxn ? dist : -dist ) ) );
}

void FillNode::addPair( Node* hit, Node* base, ReadId id, int32_t est, bool drxn )
{
    for ( FillPaired& fp : paired[drxn] ) if ( fp.node == hit )
    {
        fp.pairs.push_back( FillPair( base, id, est ) );
        return;
    }
    paired[drxn].push_back( FillPaired( hit, base, id, est ) );
}

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
        for ( pair<FillNode*, int>& fe : edges[drxn] ) if ( fe.first->getFork( !drxn ) == e.node ) ignore = true;
        if ( ignore ) continue;
        for ( FillNode* fn : fills ) if ( fn->getFork( !drxn ) == e.node && ( edged = fn ) ) break;
        assert( used.find( e.node ) == bool( edged ) );
        if ( !edged )
        {
            edged = new FillNode( e.node, nodes, fills, used );
            added.push_back( edged );
        }
        edges[drxn].push_back( make_pair( edged, e.ol ) );
        edged->edges[!drxn].push_back( make_pair( this, e.ol ) );
    }
    
    for ( int d : { 0, 1 } ) for ( FillNode* fn : added ) fn->edge( nodes, fills, used, d );
}

void FillNode::getDist( unordered_map<FillNode*, int32_t>& dists, int32_t dist, int32_t limit, int32_t& best, bool drxn )
{
    auto it = dists.insert( make_pair( this, dist ) );
    if ( !it.second )
    {
        if ( drxn ? it.first->second <= dist : dist <= it.first->second ) return;
        it.first->second = dist;
    }
    
    if ( drxn ? dist < limit : limit < dist ) for ( pair<FillNode*, int>& edge : edges[drxn] )
    {
        edge.first->getDist( dists, dist + ( ( drxn ? 1 : -1 ) * ( edge.first->len - edge.second ) ), limit, best, drxn );
    }
    if ( !bad && edges[drxn].empty() ) best = drxn ? min( best, dist ) : max( best, dist );
}

Node* FillNode::getFork( bool drxn )
{
    return ( drxn ? path.back() : path[0] ).first;
}

bool FillNode::isBlunt( int readLimit, bool drxn )
{
    return getFork( !drxn )->isBlunt( 0, readLimit, drxn );
}

bool FillNode::isStrong( int32_t dist, int32_t limit, bool drxn )
{
    for ( pair<Node*, int32_t>& nd : path ) if ( nd.first->verified_ )
    {
        if ( dist + ( drxn ? nd.first->size() + nd.second : len - nd.second ) > limit ) return true;
    }
    for ( pair<FillNode*, int>& e : edges[drxn] ) if ( e.first->isStrong( dist + len - e.second, limit, drxn ) ) return true;
    return false;
}

int FillNode::trim( bool drxn )
{
    if ( bad || !getFork( drxn )->verified_ ) return false;
    
    vector<FillNode*> strong, blunt;
    for ( pair<FillNode*, int>& e : edges[drxn] )
    {
        if ( e.first->isBlunt( 3, drxn ) ) blunt.push_back( e.first );
        else if ( e.first->isStrong( -e.second, 500, drxn ) ) strong.push_back( e.first );
    }
    
    if ( strong.empty() || blunt.empty() ) return false;
    
    for ( FillNode* fn : blunt )
    {
        getFork( drxn )->removeEdge( fn->getFork( !drxn ), drxn, true );
        fn->getFork( !drxn )->setState();
        for ( int i = 0; i < edges[drxn].size(); i++ ) if ( edges[drxn][i].first == fn ) edges[drxn].erase( edges[drxn].begin() + i-- );
        for ( int i = 0; i < fn->edges[drxn].size(); i++ ) if ( fn->edges[drxn][i].first == this ) fn->edges[drxn].erase( fn->edges[drxn].begin() + i-- );
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

bool LocusFill::fill( Querier& bwt, NodeRoll& nodes )
{
    LocusFill lf( bwt, nodes );
    if ( lf.seed( bwt, nodes ) ) lf.refresh( bwt, nodes );
    if ( lf.join( bwt, nodes ) ) lf.refresh( bwt, nodes );
    if ( lf.trim( nodes ) ) lf.refresh( bwt, nodes );
    Nodes repair;
    for ( int d : { 0, 1 } ) for ( FillNode* fn : lf.fills_ ) lf.map( bwt, nodes, repair, fn, d );
    for ( Node* node : repair.nodes ) node->setVerified();
    return false;
}

bool LocusFill::getDists( FillNode* fn, unordered_map<FillNode*, int32_t>& dists, int32_t& ext, int32_t& limit, bool drxn )
{
    limit = drxn ? 0 : fn->len;
    for ( pair<ReadId, int32_t>& read : fn->reads[drxn] ) limit = drxn ? max( limit, read.second ) : min( limit, read.second );
    if ( limit == ( drxn ? 0 : fn->len ) ) return false;
    
    ext = drxn ? limit - fn->len + 300 : limit - 300;
    fn->getDist( dists, 0, ext, ext, drxn );
    limit = drxn ? fn->len + ext - 300 : ext + 300;
    return true;
}

void LocusFill::getSeeds( Querier& bwt, FillNode* fn, vector<FillRead>& seeds, bool drxn )
{
    int32_t limit, ext;
    unordered_map<FillNode*, int32_t> dists;
    if ( fn->bad || !getDists( fn, dists, ext, limit, drxn ) ) return;
    
    vector<FillRead> reads;
    
    for ( const pair<FillNode*, int32_t>& dist : dists ) for ( pair<ReadId, int32_t>& read : dist.first->reads[drxn] )
    {
        int32_t coord = read.second + dist.second;
        bool unseeded = dist.first == fn && seeded_.find( read.first ) == seeded_.end();
        if ( drxn ? coord < ext : ext < coord ) reads.push_back( FillRead( getSeq( bwt, read.first ), read.first, coord, unseeded ) );
    }
    
    sort( reads.begin(), reads.end(), []( FillRead& a, FillRead& b ){ return a.coord < b.coord; } );
    for ( FillRead& fr : reads ) if ( fr.unseeded && drxn ? fr.coord < limit : limit < fr.coord )
    {
        int32_t coords[2];
        for ( const pair<FillNode*, int32_t>& dist : dists )
        {
            if ( mapSeqEnd( fr.seq, dist.first->seq, max( 50, fr.mapped ), coords, !drxn ) ) fr.mapped = coords[1] - coords[0];
        }
        if ( !fr.mapped || !bwt.isExtendable( fr.seq, drxn ) ) continue;
        for ( FillRead& alt : reads ) if ( fr.id != alt.id && fr.seq != alt.seq && fr.coord - 300 < alt.coord )
        {
            if ( mapSeqOverlap( fr.seq, alt.seq, max( 50, (int)fr.seq.size() - fr.mapped ), drxn ) )
            {
                seeds.push_back( fr );
                seeds.back().coord += fn->path[0].first->ends_[0] + ( drxn ? -fr.seq.size() : 0 );
                cout << seeds.back().coord << endl << fr.seq << endl;
                break;
            }
            if ( fr.coord + 300 < alt.coord ) break;
        }
    }
}

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

void LocusFill::map( Querier& bwt, NodeRoll& nodes, Nodes& repair, FillNode* fn, bool drxn )
{
    int32_t limit, ext;
    unordered_map<FillNode*, int32_t> dists;
    if ( fn->bad || !getDists( fn, dists, ext, limit, drxn ) ) return;
    
    vector<ReadId> repaired;
    for ( pair<ReadId, int32_t>& read : fn->reads[drxn] ) if ( drxn ? read.second < limit : limit < read.second )
    {
        string seq = getSeq( bwt, read.first );
        bool ext[2]{ bwt.isExtendable( seq, 0 ), bwt.isExtendable( seq, 1 ) };
        if ( ext[0] && ext[1] ) continue;
        
        int32_t coords[2];
        vector<FillHit> hits[2];
        for ( int d : { 0, 1 } ) for ( const pair<FillNode*, int32_t>& dist : dists ) if ( mapSeqEnd( seq, dist.first->seq, 50, coords, d ) )
        {
            FillHit::add( hits[d], dist.first, seq, coords, read.second, dist.second, drxn );
        }
        
        bool ignore = false;
        if ( FillHit::confirm( hits, ext, seq.size(), ignore ) )
        {
            for ( int d : { 0, 1 } ) for ( FillHit& fh : hits[d] ) fh.set( read.first, ignore );
            repaired.push_back( params.getPairId( read.first ) );
            mapped_.insert( read.first );
        }
        else if ( !ext[0] && !ext[1] )
        {
            ReadId id = params.getPairId( read.first );
            for ( int i = 0; i < fn->path.size(); i++ )
            {
                auto it = fn->path[i].first->reads_.find( id );
                if ( it == fn->path[i].first->reads_.end() ) continue;
                it->second.unpaired = true;
                assert( fn->path[i].first->rmvMark( read.first, drxn ) );
            }
        }
        else ( unmapped_ ).insert( read.first );
    }
    for ( pair<Node*, int32_t>& nd : fn->path ) for ( ReadId id : repaired ) if ( nd.first->reads_.find( id ) != nd.first->reads_.end() )
    {
        repair += nd.first;
        break;
    }
}

void LocusFill::refresh( Querier& bwt, NodeRoll& nodes )
{
    clear();
    Node::prunePaths( bwt, nodes );
    fills_ = FillNode::create( nodes );
}

bool LocusFill::seed( Querier& bwt, NodeRoll& nodes )
{
    vector<FillRead> seeds[2];
    for ( FillNode* fn : fills_ ) for ( int d : { 0, 1 } ) getSeeds( bwt, fn, seeds[d], d );
    
    vector<Node*> added[2], seeded[2];
    for ( int d : { 0, 1 } ) 
    {
        for ( FillRead& fr : seeds[d] ) Node::seedNode( bwt, nodes, added[d], seeded[d], fr.seq, fr.id, fr.coord, fr.coord > 0 );
        Nodes ext( seeded[d] );
        for ( Node* node : added[d] ) ext += node;
        for ( Node* node : ext.nodes ) node->extendFork( bwt, nodes, params.maxPeMean, 8, !d );
    }
    return !seeded[0].empty() || !seeded[1].empty();
}

bool LocusFill::trim( NodeRoll& nodes )
{
    int trimmed = 0;
    for ( int d : { 0, 1 } ) for ( FillNode* fn : fills_ ) trimmed += fn->trim( d );
    return trimmed;
}

