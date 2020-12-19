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

#include "reviewer.h"
#include "shared_functions.h"
#include <algorithm>
#include <iomanip>

vector<ReviewPair> ReviewPair::add( ReviewMap* l, ReviewMap* r )
{
    vector<ReviewPair> pairs;
    if ( !l || !r ) return pairs;
    
    if ( l->coord[0] <= r->coord[0] )pairs.push_back( ReviewPair( l->coord[0], r->coord[1] ) );
    if ( l->multi ) for ( int i : *l->multi ) if ( i <= r->coord[0] ) pairs.push_back( ReviewPair( i, r->coord[1] ) );
    if ( r->multi ) for ( int j : *r->multi ) if ( l->coord[0] <= j ) pairs.push_back( ReviewPair( l->coord[0], j + r->coord[1] - r->coord[0] ) );
    if ( l->multi && r->multi ) for ( int i : *l->multi ) for ( int j : *r->multi ) if ( i <= j ) pairs.push_back( ReviewPair( i, j + r->coord[1] - r->coord[0] ) );
    
    sort( pairs.begin(), pairs.end(), [&]( ReviewPair& a, ReviewPair& b ){ return a.coord[1]-a.coord[0] < b.coord[1]-b.coord[0]; } );
    
    return pairs;
}
int ReviewPair::add( vector<ReviewPair> pairs[4], ReviewMap* l, ReviewMap* r, int est )
{
    if ( !l || !r ) return -1;
    
    pairs[ l->coord[0] <= r->coord[0] ? 0 : 3 ].push_back( ReviewPair( l->coord[0], r->coord[1] ) );
    if ( l->multi ) for ( int i : *l->multi ) pairs[ i <= r->coord[0] ? 0 : 3 ].push_back( ReviewPair( i, r->coord[1] ) );
    if ( r->multi ) for ( int j : *r->multi ) pairs[ l->coord[0] <= j ? 0 : 3 ].push_back( ReviewPair( l->coord[0], j + r->coord[1] - r->coord[0] ) );
    if ( l->multi && r->multi ) for ( int i : *l->multi ) for ( int j : *r->multi )
    {
        pairs[ i <= j ? 0 : 3 ].push_back( ReviewPair( i, j + r->coord[1] - r->coord[0] ) );
    }
    sort( pairs[0].begin(), pairs[0].end(), [&]( ReviewPair& a, ReviewPair& b ){ return abs( a.coord[1]-a.coord[0]-est ) < abs( b.coord[1]-b.coord[0]-est ); } );
    
    for ( int i = 0; i < pairs[0].size(); i++ ) if ( abs( pairs[0][i].getDiff( est ) ) > ( i ? pairs[0][i-1].getCut( est ) : est*.5 + 200 ) )
    {
        ReviewPair rp = pairs[0][i];
        int diff = pairs[0][i].getDiff( est );
        pairs[ pairs[0][i].dist() < est ? 1 : 2 ].push_back( pairs[0][i] );
        pairs[0].erase( pairs[0].begin() + i-- );
    }
    
    return l->multi && r->multi ? 3 : ( r->multi ? 2 : ( l->multi ? 1 : 0 ) );
}

int ReviewPair::dist()
{
    return coord[1] - coord[0];
}

int ReviewPair::getCut( int est )
{
    return min( abs( coord[1] - coord[0] - est ) * 2 + 100 + ( est * 0.2 ), est * 0.5 + 200 );
}

int ReviewPair::getDiff( int est )
{
    return coord[1] - coord[0] - est;
}

vector<ReviewAlign> ReviewAlign::align( string& a, string& b )
{
    int len = min( max( 33, (int)params.readLen / 3 ), 50 );
    int minAlign = len*2;
    vector<ReviewAlign> aligns;
    ReviewAlign align;
    if ( len < a.size() ) for ( int i = 0; i < a.size(); )
    {
        int j = min( i, (int)a.size()-len );
        i += len;
        string q = a.substr( j, len );
        if ( !queryable( q ) ) continue;
        size_t it = b.find( q );
        while ( it != string::npos )
        {
            int x[2][2]{ { j, j+len }, { (int)it, (int)it+len } };
            while ( x[0][0] && x[1][0] && a[ x[0][0]-1 ] == b[ x[1][0]-1 ] ){ x[0][0]--; x[1][0]--; };
            while ( x[0][1] < a.size() && x[1][1] < b.size() && a[ x[0][1] ] == b[ x[1][1] ] ){ x[0][1]++; x[1][1]++; };
            align.len = x[0][1]-x[0][0];
            bool good = align.len >= minAlign;
            for ( ReviewAlign& rab : aligns ) if ( rab.coords[0][0] == x[0][0] && rab.coords[1][0] == x[1][0] && align.len == rab.len ) good = false;
            for ( int k : { 0, 1 } ) for ( int d : { 0, 1 } ) align.coords[k][d] = x[k][d];
            if ( good ) aligns.push_back( align );
            it = b.find( q, it+1 );
        }
    }
    
    sort( aligns.begin(), aligns.end(), []( ReviewAlign& a, ReviewAlign& b ){ return a.coords[0][0] < b.coords[0][0]; } ); 
    for ( int i = 0; i < aligns.size(); i++ ) for ( int j = i+1; j < aligns.size(); j++ ) if ( conflict( aligns[i], aligns[j] ) )
    {
        int lens[2]{ aligns[i].len, aligns[j].len };
        for ( int k = j+1; k < aligns.size(); k++ )
        {
            bool conflicted[2]{ conflict( aligns[i], aligns[k] ), conflict( aligns[j], aligns[k] ) };
            if ( conflicted[0] && !conflicted[1] ) lens[1] += aligns[k].len;
            if ( !conflicted[0] && conflicted[1] ) lens[0] += aligns[k].len;
        }
        if ( lens[1] > lens[0] )
        {
            aligns.erase( aligns.begin() + i-- );
            break;
        }
        else aligns.erase( aligns.begin() + j-- );
    }
    return aligns;
}

bool ReviewAlign::conflict( ReviewAlign& l, ReviewAlign& r )
{
    return r.coords[0][0] < l.coords[0][1] || r.coords[1][0] < l.coords[1][1];
}

bool ReviewAlign::queryable( string& q )
{
    if ( q.find( 'N' ) != string::npos ) return false;
    for ( int i = 1; i < 4; i++ )
    {
        int hits = 0, len = q.size()-i;
        for ( int j = 0; j < len; j++ ) if ( q[j] == q[j+i] ) hits++;
        if ( len - hits < 3 ) return false;
    }
    return true;
}

Reviewer::Reviewer( Querier& bwt, string s )
: seq_( s ), allele_( NULL ), cover_( NULL ), multi_( NULL ), len_( s.size() ), base_( 0 ), end_( s.size() )
{}

void Reviewer::add( ReadId id, int32_t i, int32_t j )
{
    assert( i <= seq_.size() && j <= seq_.size() );
    ReviewMap rm( i, j );
    auto ins = mapped_.insert( make_pair( id, rm ) );
    if ( ins.second || i == ins.first->second.coord[0] ) return;
    bool added = !ins.first->second.multi;
    int place = added ? 0 : ins.first->second.multi->size();
    if ( added ) ins.first->second.multi = new vector<int>{ i };
    else for ( int k = 0; !( added = ( i == (*ins.first->second.multi)[k] ) ) && k < place; k++ ) if ( i < (*ins.first->second.multi)[k] ) place = k;
    
    if ( !added ) ins.first->second.multi->insert( ins.first->second.multi->begin() + place, i );
}

ReviewMap* Reviewer::get( ReadId id )
{
    auto it = mapped_.find( id );
    if ( it != mapped_.end() )
    {
        int coords[2]{ it->second.coord[0], it->second.coord[1] };
        assert( coords[0] >= 0 && coords[1] <= len_ );
    }
    return it != mapped_.end() ? &it->second : NULL;
}

void Reviewer::setAllelic( Reviewer* alt )
{
    allele_ = alt;
    alt->allele_ = this;
    vector<ReviewAlign> aligns = ReviewAlign::align( seq_, alt->seq_ );
    ReviewStretch stretch;
    for ( int i = 0; i+1 < aligns.size(); i++ )
    {
        int dist[2]{ aligns[i+1].coords[0][0]-aligns[i].coords[0][1], aligns[i+1].coords[1][0]- aligns[i].coords[1][1] };
        for ( int j : { 0, 1 } ) if ( dist[!j] > dist[j] )
        {
            stretch.add = dist[!j] - dist[j];
            stretch.coords[0] = aligns[i].coords[j][1];
            stretch.coords[1] = aligns[i+1].coords[j][0];
            ( j ? alt : this )->stretch_.push_back( stretch );
        }
        if ( dist[0] > dist[1] ) ;
        else if ( dist[1] > dist[0] );
    }
    assert( !aligns.empty() );
    int baseDiff = aligns[0].coords[0][0] - aligns[0].coords[1][0];
    int endDiff = ( len_ - aligns.back().coords[0][1] ) - ( alt->len_ - aligns.back().coords[1][1] );
    if ( baseDiff > 0 ) base_ = baseDiff;
    if ( baseDiff < 0 ) alt->base_ = -baseDiff;
    if ( endDiff > 0 ) end_ = len_ - endDiff;
    if ( endDiff < 0 ) alt->end_ = alt->len_ + endDiff;
}

string Reviewer::getDecimal( float dec, int places )
{
    string str = to_string( int( dec ) ) + ( places ? "." : "" );
    for ( int i = 0; i < places; i++ )
    {
        dec *= 10;
        str += to_string( int( dec ) % 10 );
    }
    return str;
}

void Reviewer::review( Querier& bwt, vector< vector<Node*> >& paths )
{
    vector<Reviewer*> reviews;
    for ( vector<Node*>& path : paths ) reviews.push_back( new Reviewer( bwt, Node::getSeq( path ) ) );
    for ( int i = 0; i < paths.size(); i++ )
    {
        int len = 0;
        ReviewMap* alt;
        reviews[i]->path_ = paths[i];
        for ( int j = 0; j < paths[i].size(); j++ )
        {
            reviews[i]->coords_.push_back( len );
            if ( j ) len -= paths[i][j]->getOverlap( paths[i][j-1], 0 );
//            for ( auto& read : paths[i][j]->reads_ ) if ( read.second[2] )
            for ( auto& read : paths[i][j]->reads_ )
            {
                bool good = true;
                for ( Reviewer* r : reviews ) if ( ( alt = r->get( read.first ) ) && alt->coord[1]-alt->coord[0] > read.second[1]-read.second[0] ) good = false;
                reviews[i]->add( read.first, len + read.second[0] - paths[i][j]->ends_[0], len + read.second[1] - paths[i][j]->ends_[0] );
            }
            len += paths[i][j]->size();
        }
    }
    reviews[0]->allele_ = reviews[1];
    reviews[1]->allele_ = reviews[0];
    for ( Reviewer* r : reviews ) r->setMultis( reviews );
    
//    vector< pair<Node*, ReadId> > hits[2];
//    cout << ">Sus" << endl << paths[0][17]->seq_ << endl;
//    for ( int d : { 0, 1 } ) for ( NodeMark& nm :  paths[0][17]->pe_[d] )
//    {
//        for ( vector<Node*>& path : paths ) for ( Node* node : path ) if ( node->getRead( nm.id ) ) hits[d].push_back( make_pair( node, nm.id ) );
//        string seq = bwt.getSequence( nm.id );
//        if ( bwt.isExtendable( seq, 0 ) && bwt.isExtendable( seq, 1 ) ) cout << ">" << nm.id << endl << seq << endl;
//    }
    setReview( bwt, reviews );
}

void Reviewer::review( Querier& bwt, vector<string>& seqs )
{
    vector<Reviewer*> reviews;
    for ( string& s : seqs ) reviews.push_back( new Reviewer( bwt, s ) );
    setReview( bwt, reviews );
}

void Reviewer::review( Querier& bwt, vector<string>& seqs, vector<string>& foci, string ofn, bool diploid )
{
    vector<Reviewer*> reviews;
    for ( string& s : seqs ) reviews.push_back( new Reviewer( bwt, s ) );
    if ( diploid && reviews.size() == 2 ) reviews[0]->setAllelic( reviews[1] );
    setFoci( reviews, foci );
    for ( Reviewer* r : reviews ) r->setReads( bwt );
//    string fn = "/home/glen/Thesis/HtReviewDump";
//    dumpReads( fn, reviews );
//    loadReads( fn, reviews );
    for ( Reviewer* r : reviews ) r->setCover( reviews );
    for ( Reviewer* r : reviews ) r->setPairs( reviews );
    float median = setMedian( reviews );
    vector<int> libs = setLibs( reviews );
    if ( diploid && reviews.size() == 2 ) reviews[0]->setOut( reviews[1], libs, ofn );
    assert( false );
}

void Reviewer::setCover( vector<Reviewer*>& reviews )
{
    vector<Reviewer*> alts;
    for ( Reviewer* r : reviews ) if ( r != this ) alts.push_back( r );
    if ( cover_ ) delete cover_;
    if ( multi_ ) delete multi_;
    cover_ = new double[ seq_.size() ]{0};
    multi_ = new double[ seq_.size() ]{0};
    int* counts = new int[ seq_.size() ]{0};
    
    ReviewMap* alt;
    for ( pair<ReadId, ReviewMap> m : mapped_ )
    {
        ReadId revId = params.getRevId( m.first );
        double multi = 1 + ( m.second.multi ? m.second.multi->size() : 0 );
        if ( alt = get( revId ) ) multi += 1 + ( alt->multi ? alt->multi->size() : 0 );
        for ( Reviewer* r : alts ) if ( alt = r->get( revId ) ) multi += 1 + ( alt->multi ? alt->multi->size() : 0 );
        for ( Reviewer* r : alts ) if ( alt = r->get( m.first ) ) multi += 1 + ( alt->multi ? alt->multi->size() : 0 );
        
        vector<int> offs{ 0 };
        if ( m.second.multi ) for ( int multi : *m.second.multi ) offs.push_back( multi - m.second.coord[0] );
        for ( int i = m.second.coord[0]; i < m.second.coord[1]; i++ ) for ( int off : offs )
        {
            cover_[i+off] += (double)1 / multi;
            multi_[i+off] += multi;
            counts[i+off]++;
        }
    }
    
    for ( int i = 0; i < seq_.size(); i++ ) multi_[i] = max( (double)1, multi_[i] / (double)max( 1, counts[i] ) );
    delete counts;
}

void Reviewer::setFoci( vector<Reviewer*>& reviews, vector<string>& foci )
{
    for ( int i = 0; i < foci.size(); i++ ) for ( int j = i+1; j < foci.size(); j++ ) if ( foci[i] == foci[j] ) foci.erase( foci.begin()+j-- );
    for ( string& s : foci )
    {
        string q[2]{ s, revCompNew( s ) };
        int found = 0;
        for ( int d : { 0, 1 } ) for ( Reviewer* r : reviews ) for ( size_t it = r->seq_.find( q[d] ); it != string::npos; )
        {
            r->foci_.push_back( make_pair( it, it+s.size() ) );
            it = r->seq_.find( q[d], it+1 );
            found++;
        }
        assert( found );
    }
    for ( Reviewer* r : reviews ) sort( r->foci_.begin(), r->foci_.end(), []( pair<int, int>& a, pair<int, int>& b ){ return a.first < b.first; } ); 
}

void Reviewer::setGaps( string l, string r )
{
    size_t it = seq_.find( l );
    size_t it2 = seq_.find( r, it+1 );
    assert( it != string::npos && it2 != string::npos );
    it += l.size();
    
    bool inserted = it2 > it;
    vector< pair<int,int> > coords;
    if ( inserted ) coords.push_back( make_pair( it, it2 ) );
    if ( inserted ) coords.push_back( make_pair( it, it ) );
    coords.push_back( make_pair( it2, it2 ) );
    
    cout << ( inserted ? "Insertion at: " : "Deletion at: " ) << ( it - seeds_[0].first ) << ( inserted ? " to " + to_string( it2 - seeds_[0].first ) : "" ) << endl;
    for ( pair<int, int> coord : coords )
    {
        cout << "   Coords:   " << coord.first - seeds_[0].first
                << ( coord.second == coord.first ? "" : "   to   " + to_string( coord.second-seeds_[0].first ) ) << endl;
        for ( int i = 0; i < pairs_.size(); i++ )
        {
            int counts[3]{0};
            vector<int> dists;
            for ( ReviewPair& rp : pairs_[i] ) if ( rp.coord[0] < coord.first && coord.second < rp.coord[1] )
            {
                counts[ rp.para ? 2 : ( rp.homo ? 1 : 0 ) ]++;
                if ( !rp.homo && !rp.para ) dists.push_back( rp.coord[1] - rp.coord[0] );
            }
            cout << "      Lib " << i+1 << ":   ";
            sort( dists.begin(), dists.end() );
            int median = dists.empty() ? 0 : ( dists[ ( dists.size()-1 ) / 2 ] + dists[ dists.size() / 2 ] ) / 2;
            if ( inserted && !counts[0] ) cout << "No pairs." << endl;
            else cout << counts[0] << " unambiguous,   " << counts[1] << " homologous ambiguous,   " << counts[2] << " paralogous ambiguous.   Median: " << median << endl;
            int x = 0;
            
        }
    }
    int x = 0;
}

vector<int> Reviewer::setLibs( vector<Reviewer*>& reviews )
{
    vector<int> libs;
    for ( int i = 0; i < params.libs.size(); i++ )
    {
        unordered_set<ReadId> used;
        vector<int> dists;
        for ( Reviewer* r : reviews ) for ( ReviewPair& rp : r->pairs_[i] ) if ( !rp.dispute && used.insert( rp.id ).second ) dists.push_back( rp.dist() );
        sort( dists.begin(), dists.end() );
        assert( !dists.empty() );
        libs.push_back( ( dists[ ( dists.size()-1 ) / 2 ] + dists[ dists.size() / 2 ] ) / 2 );
    }
    return libs;
}

float Reviewer::setMedian( vector<Reviewer*>& reviews )
{
    vector<float> covers;
    
    for ( Reviewer* r : reviews ) for ( pair<int, int> seed : r->seeds_ ) for ( int i = seed.first; i < seed.second; i++ ) covers.push_back( r->cover_[i] );
    for ( Reviewer* r : reviews ) for ( pair<int, int> seed : r->foci_ ) for ( int i = seed.first; i < seed.second; i++ ) covers.push_back( r->cover_[i] );
    
    sort( covers.begin(), covers.end() );
    
    if ( covers.empty() ) return 0;
    return ( covers[ ( covers.size()-1 ) / 2 ] + covers[ covers.size() / 2 ] ) / 2;
}

void Reviewer::setMultis( vector<Reviewer*>& reviews )
{
    int hCoord[2]{ len_, len_ }, pCoord[2]{ len_, len_ };
    for ( int i = 0; i < path_.size(); i++ )
    {
        int homo = 0;
        bool para = false;
        for ( Node* clone : path_[i]->clones( false ) ) for ( Reviewer* r : reviews ) if ( find( r->path_.begin(), r->path_.end(), clone ) != r->path_.end() )
        {
            if ( r == this || r != allele_ || homo++ ) para = true;
        }
        
        if ( para )
        {
            pCoord[0] = min( pCoord[0], coords_[i] );
            pCoord[1] = coords_[i] + path_[i]->size();
        }
        else if ( homo )
        {
            hCoord[0] = min( hCoord[0], coords_[i] );
            hCoord[1] = coords_[i] + path_[i]->size();
        }
        
        if ( pCoord[0] < len_ && ( path_[i] == path_.back() || !para ) )
        {
            homology_[1].push_back( make_pair( pCoord[0], pCoord[1] ) );
            pCoord[0] = pCoord[1] = len_;
        }
        
        if ( hCoord[0] < len_ && ( path_[i] == path_.back() || para || !homo ) )
        {
            homology_[0].push_back( make_pair( hCoord[0], hCoord[1] ) );
            hCoord[0] = hCoord[1] = len_;
        }
    }
}

void Reviewer::setOut( Reviewer* allele, vector<int>& libs, string ofn )
{
    assert( libs.size() == params.libs.size() );
    vector<Reviewer*> reviews{ this };
    if ( allele ) reviews.push_back( allele );
    int len = 0;
    for ( Reviewer* r : reviews )
    {
        int rLen = r->end_ - r->base_;
        for ( ReviewStretch rs : r->stretch_ ) rLen += rs.add;
        assert( !len || rLen == len );
        len = rLen;
    }
    
    vector<float> cover( len, 0 ), multi( len, 0 );
    vector<int> counts( len, 0 ), pairCounts( len, 0 );
    vector< vector<int> > pairMedians( params.libs.size(), vector<int>( len, 0 ) );
    vector<bool> foci( len, false );
    
    vector< vector< vector<ReviewPair> > > pairs( reviews.size(), vector< vector<ReviewPair> >( params.libs.size() ) );
    vector< vector<int> > cur( reviews.size(), vector<int>( params.libs.size(), 0 ) );
    vector<int> stretchLeft( reviews.size(), 0 ), stretchNext( reviews.size(), 0 ), stretchBlock( reviews.size(), 0 );
    vector<int> stretchDump( reviews.size(), 0 ), stretchDumpCount( reviews.size(), 0 ), stretch( reviews.size(), 0 ), base;
    for ( Reviewer* r : reviews ) base.push_back( r->base_ );
    
    for ( int i = 0; i < len; i++ )
    {
        for ( int j = 0; j < reviews.size(); j++ )
        {
            cover[i] += reviews[j]->cover_[ base[j] ];
            multi[i] += reviews[j]->multi_[ base[j] ];
            for ( pair<int, int> focus : reviews[j]->foci_ ) if ( focus.first <= base[j] && base[j] < focus.second ) foci[i] = true;
            for ( int k = 0; k < pairs[j].size(); k++ ) reviews[j]->setSpan( pairs[j][k], base[j], k, cur[j][k] );
            
            // Are we entering a segment that neeeds to be stretched?
            if ( !stretchLeft[j] && stretch[j] < reviews[j]->stretch_.size() && reviews[j]->stretch_[ stretch[j] ].coords[0] == base[j] )
            {
                int dist = reviews[j]->stretch_[ stretch[j] ].coords[1] - reviews[j]->stretch_[ stretch[j] ].coords[0];
                stretchDumpCount[j] = reviews[j]->stretch_[ stretch[j] ].add >= dist ? reviews[j]->stretch_[ stretch[j] ].add - max( 0, dist-1 ) : 0;
                stretchDump[j] = stretchDumpCount[j] ? reviews[j]->stretch_[stretch[j]].coords[0] + dist / 2 : 0;
                stretchLeft[j] = reviews[j]->stretch_[ stretch[j] ].add - stretchDumpCount[j];
                stretchBlock[j] = stretchNext[j] = stretchLeft[j] ? max( 1, dist / ( stretchLeft[j]+1 ) ) : 0;
                stretch[j]++;
            }
            
            if ( stretchNext[j] ) stretchNext[j]--;
            if ( stretchDumpCount[j] && stretchDump[j] == base[j] ) stretchDumpCount[j]--;
            else if ( stretchNext[j] || !stretchLeft[j] ) base[j]++;
            else if ( stretchLeft[j] && --stretchLeft[j] ) stretchNext[j] = stretchBlock[j];
        }
        
        for ( int j = 0; j < pairMedians.size(); j++ )
        {
            unordered_map<ReadId, int> used;
            for ( int k = 0; k < reviews.size(); k++ ) for ( ReviewPair& rp : pairs[k][j] )
            {
                auto ins = used.insert( make_pair( rp.id, rp.dist() ) );
                if ( !ins.second && abs( rp.dist()-libs[j] ) < abs( ins.first->second-libs[j] ) ) ins.first->second = rp.dist();
            }
            vector<int> dists;
            for ( pair<ReadId, int> p : used ) dists.push_back( p.second );
            sort( dists.begin(), dists.end() );
            pairMedians[j][i] = dists.empty() ? 0 : ( dists[ ( dists.size()-1 ) / 2 ] + dists[ dists.size() / 2 ] ) / 2;
            pairCounts[i] += used.size();
        }
    }
    
//    for ( Reviewer* r : reviews )
//    {
//        int i = r->base_, k = -1, stretchLeft = 0, stretchNext = 0, stretchBlock = 0, stretchDump = 0, stretchDumpCount = 0;
//        for ( int j = 0; j < len; j++ )
//        {
//            cover[j] += r->cover_[i];
//            multi[j] += r->multi_[i];
//            for ( pair<int, int> focus : r->foci_ ) if ( focus.first <= i && i < focus.second ) foci[j] = true;
//            
//            for ( int p = 0; p < pairs.size(); p++ ) setSpan( pairs[p], i, p, cur[p] );
//            
//            // Are we entering a segment that neeeds to be stretched?
//            if ( !stretchLeft && k+1 < r->stretch_.size() && r->stretch_[k+1].coords[0] == i )
//            {
//                int dist = r->stretch_[++k].coords[1] - r->stretch_[k].coords[0];
//                stretchDumpCount = r->stretch_[k].add >= dist ? r->stretch_[k].add - max( 0, dist-1 ) : 0;
//                stretchDump = stretchDumpCount ? r->stretch_[k].coords[0] + dist / 2 : 0;
//                stretchLeft = r->stretch_[k].add - stretchDumpCount;
//                stretchBlock = stretchNext = max( 1, dist / ( stretchLeft+1 ) );
//            }
//            
//            if ( stretchDumpCount && stretchDump == i ) stretchNext = 0;
//            if ( stretchNext || !stretchLeft ) i++;
//            if ( stretchNext ) stretchNext--;
//            else if ( stretchDumpCount && stretchDump == i ) stretchDumpCount--;
//            else if ( stretchLeft && --stretchLeft ) stretchNext = stretchBlock;
//        }
//        assert( k+1 == r->stretch_.size() );
//    }
    int limits[2]{ 0, len };
    if ( !foci[0] ) for ( int i = 1; i < len; i++ ) if ( foci[i] && ( limits[0] = i ) ) break;
    if ( !foci.back() ) for ( int i = len-1; i-- > 0; ) if ( foci[i] && ( limits[1] = i ) ) break;
    int coords[2]{ max( 0, limits[0]-1000 ), min( len, limits[1]+1000 ) };
    
    int f = 38;
    ofstream ofs( ofn );
    ofs << setprecision( 3 );
    ofs << "Coordinates,Coverage";
    if ( !foci.empty() ) ofs << ",Focus";
    for ( int i = 0; i < libs.size(); i++ ) ofs << ",Median paired dist " << i+1;
    ofs << ",Paired count\n";
    for ( int i = coords[0]; i < coords[1]; i++ )
    {
        ofs << i-coords[0]+1 << "," << cover[i];
        if ( !foci.empty() ) ofs << "," << ( foci[i] ? to_string( f ) : "" );
        if ( !foci.empty() ) ofs << "," << ( foci[i] ? to_string( 4500 ) : "" );
        for ( int j = 0; j < libs.size(); j++ ) ofs << "," << pairMedians[j][i] << "," << to_string( libs[j] );
        ofs << "," << pairCounts[i] << "," << ( foci[i] ? to_string( 200 ) : "" ) << "\n";
    }
    ofs.close();
    
    assert( false );
}

void Reviewer::setOut( vector<int>& libs, ofstream* ofs, float median, int spaced, int cut )
{
    for ( vector<ReviewPair>& p : pairs_ ) sort( p.begin(), p.end(), []( ReviewPair& a, ReviewPair& b ){ return a.coord[0] < b.coord[0]; } );
    vector< vector<ReviewPair> > pairs( params.libs.size() );
    vector<int> cur( pairs.size(), 0 );
    vector< pair<int, int> > cuts[2][cut];
    vector< pair<int, int> > paraAmbig[5], homoAmbig, unpair[20];
    
    auto blockUpdate = []( int32_t coord, vector< pair<int, int> >& blocks )
    {
        if ( blocks.empty() || blocks.back().second < coord ) blocks.push_back( make_pair( coord, coord ) );
        blocks.back().second = coord+1;
    };
    
    if ( ofs )
    {
        *ofs << setprecision( 3 );
        *ofs << "Coordinates,Coverage,Median Coverage,Exon,Multiplicity,Exon\n";
    }
    
    spaced = max( 1, spaced );
    int limits[2]{ seeds_.empty() ? 0 : max( 0, seeds_[0].first-1000 ), seeds_.empty() ? len_ : min( len_, seeds_.back().second+1000 ) };
    int limit = 0, rec = spaced / 2, coord = 1;
    int unambigMin[2]{0};
    int libMins[ pairs.size() ]{0};
    int totalMin = 0;
    
    string medstr = to_string( ( int( median + 10 ) / 10 ) * 10 );
    for ( int i = limits[0]; i < limits[1]; i++ )
    {
        int count[2]{0};
        int unambig[2]{0}, total = 0;
        for ( int j = 0; j < pairs.size(); j++ ) setSpan( pairs[j], i, j, cur[j] );
        for ( int j = 0; j < pairs.size(); j++ ) for ( ReviewPair& rp : pairs[j] ) if ( !rp.para ) unambig[ rp.homo ]++;
        for ( int j = 0; j < pairs.size(); j++ ) total += pairs[j].size();
        unambig[1] += unambig[0];
        
        if ( i == limits[0] )
        {
            for ( int j : { 0, 1 } ) unambigMin[j] = unambig[j];
            for ( int j = 0; j < pairs.size(); j++ ) libMins[j] = pairs[j].size();
            totalMin = total;
        }
        
        for ( int j : { 0, 1 } ) unambigMin[j] = min( unambigMin[j], unambig[j] );
        for ( int j = 0; j < pairs.size(); j++ ) libMins[j] = min( libMins[j], (int)pairs[j].size() );
        totalMin = min( totalMin, total );
        
        if ( !unambig[0] ) blockUpdate( i, homoAmbig );
        for ( int j = unambig[1]; j < 5; j++ ) blockUpdate( i, paraAmbig[j] );
        for ( int j = total; j < 20; j++ ) blockUpdate( i, unpair[j] );
        
        for ( int j = 0; j < pairs.size(); j++ ) count[!params.libs[j].isPe] += pairs[j].size();
        for ( int mp : { 0, 1 } ) for ( int j = count[mp]; j < cut; j++ )
        {
            if ( cuts[mp][j].empty() || cuts[mp][j].back().second < i ) cuts[mp][j].push_back( make_pair( i, i ) );
            cuts[mp][j].back().second = i+1;
        }
        
        if ( ( i % spaced ) == rec )
        {
            bool seeded = false;
            for ( pair<int, int> seed : seeds_ ) if ( seed.first <= i && i < seed.second ) seeded = true;
            if ( ofs ) *ofs << coord << "," << cover_[limit] << "," << median << "," << ( seeded ? medstr : "" ) << "," << multi_[limit] << "," << ( seeded ? "3" : "" ) << "\n";
            coord += spaced;
            limit = i;
        }
        else if ( cover_[i] > cover_[limit] ) limit = i;
    }
    
    for ( pair<int,int>& block : homoAmbig ) cout << "Homologous unspanned segment: " << block.first << " to " << block.second << endl;
    if ( !homoAmbig.empty() ) cout << endl;
    
    for ( pair<int,int>& block : paraAmbig[4] )
    {
        vector<int> counts( params.libs.size(), 0 );
        setSpanned( counts, block.first, block.second );
        cout << "Paralogous segment with less than 5 unambiguous pairs: " << block.first << " to " << block.second << endl;
        cout << "Spanned by lib of size:   ";
        for ( int i = 0; i < libs.size(); i++ ) cout << ( i ? ",   " : "" ) << libs[i] << " (x" << counts[i] << ")";
        cout << endl;
    }
    
    if ( ofs ) ofs->close();
}

void Reviewer::setPairs( vector<Reviewer*>& reviews )
{
    Lib* lib;
    ReadId id;
    int drxn;
    ReviewMap* paired[2],* alt[2];
    
    counts_.clear();
    counts_.resize( params.libs.size(), 0 );
    pairs_.resize( params.libs.size() );
    for ( vector<ReviewPair>& pairs : pairs_ ) pairs.clear();
    for ( int i : { 0, 1 } ) for ( int j : { 0, 1 } )
    {
        inverts_[i][j].clear();
        shorts_[i][j].clear();
        longs_[i][j].clear();
        if ( !j ) mispair_[i].clear();
        if ( !j ) unpair_[i].clear();
        inverts_[i][j].resize( params.libs.size() );
        shorts_[i][j].resize( params.libs.size() );
        longs_[i][j].resize( params.libs.size() );
        if ( !j ) mispair_[i].resize( params.libs.size() );
        if ( !j ) unpair_[i].resize( params.libs.size() );
    }
    
    for ( pair<ReadId, ReviewMap> m : mapped_ ) if ( ( id = m.first ) && ( lib = params.getLib( m.first ) ) && lib->getPair( id, drxn ) && drxn )
    {
        if ( !( paired[!drxn] = &m.second ) || !( paired[drxn] = get( id ) ) ) continue;
        
        int iLib = 0;
        while ( iLib < params.libs.size() && ( &params.libs[iLib] != lib ) ) iLib++;
        assert( iLib < params.libs.size() );
        counts_[iLib]++;
        
        bool para = paired[0]->multi || paired[1]->multi, homo = false, multi = false, invalid = false, dispute = false;
        vector<ReviewPair> pairs = ReviewPair::add( paired[0], paired[1] );
        for ( int i = 0; i < pairs.size(); i++ ) if ( pairs[i].dist() > lib->size*2 + 100 ) pairs.erase( pairs.begin() + i-- );
        if ( pairs.empty() )  continue;
        
        int best = pairs[0].dist();
        for ( ReviewPair& rp : pairs ) if ( abs( rp.dist() - lib->size ) < abs( best - lib->size ) ) best = rp.dist();
        bool bestLong = best > lib->size * 1.3 + 200;
        
        for ( Reviewer* r : reviews ) if ( r != this )
        {
            alt[0] = r->get( m.first );
            alt[1] = r->get( id );
            if ( !alt[0] && !alt[1] ) continue;
            for ( int d : { 0, 1 } ) if ( alt[d] && alt[d]->multi ) para = true;
            ( allele_ == r ? homo : para ) = true;
            
            if ( alt[0] && alt[1] ) for ( ReviewPair& rp :  ReviewPair::add( alt[0], alt[1] ) )
            {
                if ( bestLong && abs( rp.dist() - lib->size ) < abs( best - lib->size ) ) invalid = true;
                if ( abs( best - rp.dist() ) > 100 ) dispute = true;
            }
        }
        
        if ( invalid || ( bestLong && pairs[0].dist() < best ) ) continue;
        if ( pairs.size() > 1 && abs( best - lib->size ) > lib->size / 2 && abs( pairs[0].dist() - pairs.back().dist() ) > 200 )
        {
            assert( false );
            continue;
        }
        
        assert( pairs.size() < 9 );
        for ( ReviewPair& rp : pairs ) if ( abs( best - ( rp.dist() ) ) > 100 ) dispute = true;
        for ( ReviewPair& rp : pairs )
        {
            rp.homo = homo;
            rp.para = para;
            rp.dispute = dispute;
            rp.id = m.first;
            if ( rp.dist() < ( lib->size * 1.5 + 200 ) ) pairs_[iLib].push_back( rp );
        }
    }
    for ( vector<ReviewPair>& p : pairs_ ) sort( p.begin(), p.end(), []( ReviewPair& a, ReviewPair& b ){ return a.coord[0] < b.coord[0]; } );
}

void Reviewer::setReads( Querier& bwt )
{
    vector<ReadId> ids;
    vector<int32_t> coords[2];
    bwt.mapSequence( seq_, ids, coords );
    for ( int i = 0; i < ids.size(); i++ ) add( ids[i], coords[0][i], coords[1][i] );
}

void Reviewer::setReview( Querier& bwt, vector<Reviewer*>& reviews )
{
    string prefix = "/home/glen/Genomes/Hp/HpOut";
    vector<string> seeds[3];
    
    seeds[0].push_back( "CTCACGCACGAAGAGATTTCAATGAACGGCGAGGAGAGGAGAATGGT" );
    seeds[0].push_back( "CTCACGCACGACGAGAGTTCAATGAACGGCGAGGAAATGAGAATGGC" );
    seeds[0].push_back( "CTCACGCACAAAGAGATTTCAAAGAACGGCGAGGAGAGGAAAATGGT" );
    seeds[0].push_back( "CTCACGCACGAAGAGATTTCAACGAACGACGAGGAGAGGAAAATGGT" );
    seeds[1].push_back( "TCAATGAAGTAGACACCAATGCGATCGCCGAGGTGTAG" );
    seeds[1].push_back( "TCAATGGAGGAGACATCAATGTGGTCGCCGAGGTGTAG" );
    seeds[1].push_back( "TCAATGGAGAAGACATCAATGTGGTCGCCGAGGTGTAG" );
    seeds[1].push_back( "TCAATGAAGAAGAAATCAATGTGGTCGCCGAGGTGTAG" );
    seeds[2].push_back( "ATGAAGGTGAAAGTGACACTGATCGTTGCCATTGTGGCTGCTCTTGCTATCTCGG" );
    seeds[2].push_back( "ATGGAGGTGAAAGTGACACTGATCGTTGCCATTGTGACTGTTCTTGCTATCTCGG" );
    seeds[2].push_back( "ATGGAGGTGAAAGTGACACTGATCGTTGCCATTGTGGCTGCTCTTGCTATCTCGG" );
    seeds[2].push_back( "ATGGAGGTGAAAGTGACACTGATCGTTGCTGTTGTGGCTGCTCTTGCTATCTCGG" );
    seeds[2].push_back( "ATGGTGGTGAAAGTGACACTGATCGTTGCCATTGTGGCTGCTATTGCTATCTCGG" );
    seeds[2].push_back( "ATGGAGGTGAAAATGACACTGATTGTTGCCATTGTGGCTGTTCTAGCTATCTCGG" );
    
    for ( Reviewer* r : reviews ) r->setSeeds( seeds );
    for ( Reviewer* r : reviews ) r->setCover( reviews );
    for ( Reviewer* r : reviews ) r->setPairs( reviews );
    float median = setMedian( reviews );
    vector<int> libs = setLibs( reviews );
    
    reviews[1]->setGaps( "AATGCAGAAGTTTGAAGTCGGTGACATGTCTATT", "GTCGGTCGGCAAAACAAGTTTTTTTTCCCTATATACGAACTGAGA" );
    reviews[0]->setGaps( "AATGCAGAAGTTTGAAGTCGGTGACATGTATATT", "GTCGGTCGGCAAAACAATTTTTTTTTCCCTATATACGAACTGAGA" );
    reviews[1]->setGaps( "ATTCAGATAAAGAGTCGAAACGAATCAACGAATCGTCTACTTTCTAATAATTAAATGAAACGACAAAATTTC", "CGCGCTTCGCGCGGGAGGGGCAGAACTGAATTCAGATAAGGATTAAGAGTCG" );
    reviews[0]->setGaps( "ATTCAGATAAAGAGTCGAAACGAATCAACGAATCGTCTACTTTCTAATAATTAAATGAAACGACAAAATTTC", "CGCGCTTCCCAAGAGGGAGGGGCAGAACTGAATTCAGATAAGGATTAAGAGTCA" );
    reviews[1]->setGaps( "TACACCCCCAACTGAAAAAGGTCTAGCCTGGACAGATGCATCCTCGAAAATT", "AGGATAAAGGTTCAACACATTAAGAAATGCATACATTGAA" );
    reviews[0]->setGaps( "TACACCCCCCAACTGAAAAAGGTCTAGCGACACCCTGGACAGATGCATCCTCGAAAATT", "AGGATAAAGGTTAAACACATTAAGAAATGCATACATTGAA" );
    reviews[1]->setGaps( "GAGGGATAGCGTGTGTGTGTGTGTGTGTTTGTGTTTGTGT", "GTGTGAGAGAGAGGGAGGGAGAGAGAAAGAGAGGAG" );
    reviews[0]->setGaps( "GAGGGATAGCGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTTTGTGT", "GTGTGAGAGAGAGGGAGGGAGAGAGAAAGAGAGGAGAGAG" );
    reviews[1]->setGaps( "GATCTCCCGTCTCATTATGACCTTGATGACCTTCTGTGTGGTTGCGTCTGCCGAAGGGACGGTCTCCGAAAGGCCTTCTACCGAAGTGGTTAAACCTGAAGGGCCTTGTAGCATTATGTTTGTCCTGATCTTGATCTCCCGTCTCATTATGACCTTGGTG", "GCCAAAGCGACCTTGTCCTCTTTCTCTGCCATTCTCATTTCCTCGCCGTTCATTGAACTCTCGTCGTGCGTG" );
    reviews[0]->setGaps( "GATCTCCTGTCTCATTATGACCTTGGTGACCTTCTGTGTGGTTGCGTCTGCCGAAAGGACGGTCTCCGAAAGGCCTTCTACCGAAGTGGTTAAACCTGAAGGGCCTTGTAGCATTATGTTTGTCCTGATCTTGATCTCCCGTCTCATTATGACCTTGGTG", "GCCAAAGCGACCTTGTCCTCTTTCTCTGCCATTCTCATTTCCTCGCCGTTCATTGAACTCTCGTCGTGCGTG" );
    reviews[1]->setGaps( "ATCATATTTGGCTCTATACATTGGTTTATGACAGGTGACAAATAAGAATGGTATTACTCGCCTTCTTTGATATTCTTC", "CCTTTATGTTTAACCCTACTTCGTAATGCGGTGTAATGAACGCCTAATTCTGACTGAATTG" );
    reviews[0]->setGaps( "ATCATATTTGGCTCTATACATTGGTTTATGACAGGTGACAAATAAGAATGGTAACTCGCCTTCTTTGATATTCTTC", "CCTTTATGTTTAACCCTACTCCGTAATGCGGTTCAATGAACGCCTAATTCTGACTGAATTG" );
    reviews[1]->setGaps( "ATTCAAATTCTTTCGTGTTACCCTAGAGATGTTTTGCGGCACATCGGTTTGAAAACGCTTGTTTAGGAAATTAG", "CCTGCTACTACTACTTCTACTACAATAGGCCTACTCTTTATTGCCATAGTCTGATGGGAGTGAGAG" );
    reviews[0]->setGaps( "ATTCAAATTCTTTCGTGTTACCCTAGAGATGTTTTGCGGCACATCGGTTTGAAAACGCTTGTTTAGGAAATTAG", "CCTGCTACTACTACTTCTACTACAATAGGCCTACTCTTTATTGCCATAGTCTGATGGGAGTGAG" );
    reviews[1]->setGaps( "TCCAGTGTCACAAATTAGGTATCAAAATATGCTCATCTGACCCTAATTATATACCTATATTAAAAAAACCT", "AAGCTTGCGGGGGTTCCACCCCTGGATGCCCATCTGACCGTAATTAAATGCCTAAAACCCCTGAACTTCCGGGGCTCT" );
    reviews[0]->setGaps( "TCCAGTGTCACAAATTAGGTATCAAAATATGCTCATCTGACCCTAATTATATACCTATATAAAAAAAACCT", "GAGCTTGCGGGGGTTCCACCCCTGGATGCCCATCTGACCGTAATTAAATGCATAAAACCCATGAACTTCCGGGGCTCT" );
    

    
    int locus = 0;
    for ( Reviewer* r : reviews ) if ( !r->seeds_.empty() )
    {
        cout << "Locus: " << ++locus << ",   length: " << r->len_ << ",   seed limits: " << r->seeds_[0].first << " to " << r->seeds_.back().second << endl;
        for ( int i = 0; i*2 < r->seeds_.size(); i++ )
        {
            cout << "Gene " << i+1 << ":   " << r->seeds_[i*2].first << " - " << r->seeds_[i*2].second;
            cout <<" and " << r->seeds_[i*2+1].first << " - " << r->seeds_[i*2+1].second << endl;
        }
        cout << endl;
        ofstream ofs( prefix + "-" + to_string( locus ) + "_coverage.csv" );
        r->setOut( libs, &ofs, median, 50, 10 );
        cout << endl;
    }
    
    
    
//    for ( int i = 0; i < reviews.size(); i++ )
//    {
//        int limits[2]{ reviews[i]->seeds_[0].first-1000, reviews[i]->seeds_.back().second+1000 };
//        if ( i ) cout << endl;
//        cout << "Cluster " << i+1 << ", Length: " << limits[1]-limits[0]-2000 << " (" << limits[1]-limits[0] << ")" << endl;
//        
//        float intron = 0.0003, exon = 0.0009, len = limits[0];
//        int introns = 0, exons = 0;
//        
//        cout << "Exons:  ";
//        for ( int j = 0; j < reviews[i]->seeds_.size(); j++ )
//        {
//            cout << ( j ? "  |  " : "" ) << reviews[i]->seeds_[j].first-limits[0] << " - " << reviews[i]->seeds_[j].second-limits[0];
//        }
//        cout << endl;
//        
//        cout << "Measures:  ";
//        for ( int j = 0; j < reviews[i]->seeds_.size(); j++ )
//        {
//            introns += reviews[i]->seeds_[j].first - len;
//            float exLen = float( reviews[i]->seeds_[j].second - reviews[i]->seeds_[j].first ) * exon;
//            cout << ( j ? "  |  " : "" ) << getDecimal( introns * intron + exons * exon, 2 ) << " - " << getDecimal( exLen, 2 );
//            exons += reviews[i]->seeds_[j].second - reviews[i]->seeds_[j].first;
//            len = reviews[i]->seeds_[j].second;
//        }
//        introns += 1000;
//        cout << endl << "Total measure: " << getDecimal( introns * intron + exons * exon, 2 ) << endl << endl;
//    }
    
//    ofstream ofsSeqs( prefix + "_loci.fa" );
//    for ( int i = 0; i < reviews.size(); i++ ) ofsSeqs << ">Locus_" << i+1 << "\n" << reviews[i]->seq_ << "\n";
//    ofsSeqs.close();
//    for ( int i = 0; i < reviews.size(); i++ )
//    {
//        ofstream ofs( prefix + "-" + to_string( i+1 ) + "_coverage.csv" );
//        reviews[i]->setOut( reviews, &ofs, median, 50, 10 );
//    }
    assert( false );
//    Lib* lib;
//    ReadId id;
//    int d;
//    for ( int i = 0; i < reviews.size(); i++ )
//    {
//        for ( pair<ReadId, ReviewMap>& m : reviews[i].mapped_ ) if ( ( id = m.first ) && ( lib = params.getLib( m.first ) ) && lib->getPair( id, d ) )
//        {
//            bool found = false;
//            for ( int j = 0; j < reviews.size(); j++ )
//            {
//                auto it = reviews[j].mapped_.find( id );
//            }
//        }
//    }
}



void Reviewer::setSeeds( vector<string> seeds[3] )
{
    vector<int> coords[2][2];
    for ( int i = 0; i < 3; i++ ) for ( string& s : seeds[i] )
    {
        string q[2]{ s, revCompNew( s ) };
        for ( int d : { 0, 1 } ) for ( size_t it = seq_.find( q[d] ); it != string::npos; )
        {
            if ( i == 2 ) seeds_.push_back( make_pair( it, it + q[d].size() ) );
            else coords[d][ d ? !i : i ].push_back( i == d ? it : it + q[d].size() );
            it = seq_.find( q[d], it+1 );
        }
    }
    
    for ( int d : { 0, 1 } ) assert( coords[d][0].size() == coords[d][0].size() );
    for ( int d : { 0, 1 } ) for ( int i : { 0, 1 } ) sort( coords[d][i].begin(), coords[d][i].end() );
    for ( int d : { 0, 1 } ) for ( int i = 0; i < coords[d][0].size(); i++ ) assert( coords[d][0][i] < coords[d][1][i] );
    for ( int d : { 0, 1 } ) for ( int i = 0; i < coords[d][0].size(); i++ ) seeds_.push_back( make_pair( coords[d][0][i], coords[d][1][i] ) );
    
    sort( seeds_.begin(), seeds_.end(), []( pair<int, int>& a, pair<int, int>& b ){ return a.first < b.first; } ); 
}

void Reviewer::setSpan( vector<ReviewPair>& pairs, int base, int lib, int& cur )
{
    while ( !pairs.empty() && pairs.back().coord[1] <= base ) pairs.pop_back();
    for ( ; cur < pairs_[lib].size() && pairs_[lib][cur].coord[0] <= base; cur++ )
    {
        int i = pairs.size();
        for ( int j = 0; j < i; j++ ) if ( pairs_[lib][cur].coord[1] > pairs[j].coord[1] ) i = j;
        pairs.insert( pairs.begin()+i, pairs_[lib][cur] );
    }
}

void Reviewer::setSpanned( vector<int>& counts, int l, int r )
{
    assert( pairs_.size() == counts.size() );
    for ( int i = 0; i < pairs_.size(); i++ ) for ( ReviewPair& rp : pairs_[i] ) if ( rp.coord[0] < l && r < rp.coord[1] ) counts[i]++;
}

void Reviewer::dumpReads( string ofn, vector<Reviewer*>& reviews )
{
    ofstream ofs( ofn );
    for ( Reviewer* r : reviews )
    {
        int counted = 0;
        string reads;
        for ( pair<ReadId, ReviewMap> mapped : r->mapped_ )
        {
            reads += to_string( mapped.first ) + "," + to_string( mapped.second.coord[0] ) + "," + to_string( mapped.second.coord[1] ) + ";";
            if ( mapped.second.multi ) for ( int i : *mapped.second.multi )
            {
                reads += to_string( mapped.first ) + "," + to_string( i ) + "," + to_string( i + mapped.second.coord[1]-mapped.second.coord[0] ) + ";";
                counted++;
            }
            counted++;
        }
        ofs << reads << "\n";
        cout << "Dumped " << counted << " reads." << endl;
    }
    ofs.close();
}

void Reviewer::loadReads( string ifn, vector<Reviewer*>& reviews )
{
    ifstream ifs( ifn );
    string line;
    for ( Reviewer* r : reviews )
    {
        int counted = 0;
        assert( getline( ifs, line ) );
        size_t i = 0, j;
        while ( ( j = line.find( ',', i ) ) != string::npos )
        {
            ReadId id = stoul( line.substr( i, j-i ) );
            i = j+1;
            assert( ( j = line.find( ',', i ) ) != string::npos );
            int ii = stoi( line.substr( i, j-i ) );
            i = j+1;
            assert( ( j = line.find( ';', i ) ) != string::npos );
            int jj = stoi( line.substr( i, j-i ) );
            i = j+1;
            r->add( id, ii, jj );
            counted++;
        }
        cout << "Loaded " << counted << " reads." << endl;
    }
    ifs.close();
}

