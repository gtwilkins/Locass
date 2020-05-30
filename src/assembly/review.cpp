/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "review.h"
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

Review::Review( Querier& bwt, string s )
: seq_( s ), allele_( NULL ), cover_( NULL ), multi_( NULL ), len_( s.size() )
{
//    vector<ReadId> ids;
//    vector<int32_t> coords[2];
//    bwt.mapSequence( seq_, ids, coords );
//    for ( int i = 0; i < ids.size(); i++ ) add( ids[i], coords[0][i], coords[1][i] );
}

void Review::add( ReadId id, int32_t i, int32_t j )
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

ReviewMap* Review::get( ReadId id )
{
    auto it = mapped_.find( id );
    return it != mapped_.end() ? &it->second : NULL;
}

string Review::getDecimal( float dec, int places )
{
    string str = to_string( int( dec ) ) + ( places ? "." : "" );
    for ( int i = 0; i < places; i++ )
    {
        dec *= 10;
        str += to_string( int( dec ) % 10 );
    }
    return str;
}

void Review::review( Querier& bwt, vector< vector<Node*> >& paths )
{
    vector<Review*> reviews;
    for ( vector<Node*>& path : paths ) reviews.push_back( new Review( bwt, Node::getSeq( path ) ) );
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
                for ( Review* r : reviews ) if ( ( alt = r->get( read.first ) ) && alt->coord[1]-alt->coord[0] > read.second[1]-read.second[0] ) good = false;
                reviews[i]->add( read.first, len + read.second[0] - paths[i][j]->ends_[0], len + read.second[1] - paths[i][j]->ends_[0] );
            }
            len += paths[i][j]->size();
        }
    }
    reviews[0]->allele_ = reviews[1];
    reviews[1]->allele_ = reviews[0];
    for ( Review* r : reviews ) r->setMultis( reviews );
    
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

void Review::review( Querier& bwt, vector<string>& seqs )
{
    vector<Review*> reviews;
    for ( string& s : seqs ) reviews.push_back( new Review( bwt, s ) );
    setReview( bwt, reviews );
}

void Review::setCover( vector<Review*>& reviews )
{
    vector<Review*> alts;
    for ( Review* r : reviews ) if ( r != this ) alts.push_back( r );
    if ( cover_ ) delete cover_;
    if ( multi_ ) delete multi_;
    cover_ = new float[ seq_.size() ]{0};
    multi_ = new float[ seq_.size() ]{0};
    
    ReviewMap* alt;
    for ( pair<ReadId, ReviewMap> m : mapped_ )
    {
        ReadId revId = params.getRevId( m.first );
        int multi = 1 + ( m.second.multi ? m.second.multi->size() : 0 );
        if ( alt = get( revId ) ) multi += 1 + ( alt->multi ? alt->multi->size() : 0 );
        for ( Review* r : alts ) if ( alt = r->get( revId ) ) multi += 1 + ( alt->multi ? alt->multi->size() : 0 );
        for ( Review* r : alts ) if ( alt = r->get( m.first ) ) multi += 1 + ( alt->multi ? alt->multi->size() : 0 );
        
        vector<int> offs{ 0 };
        if ( m.second.multi ) for ( int multi : *m.second.multi ) offs.push_back( multi - m.second.coord[0] );
        for ( int i = m.second.coord[0]; i < m.second.coord[1]; i++ ) for ( int off : offs )
        {
            cover_[i+off]++;
            multi_[i+off] += multi;
        }
    }
    
    for ( int i = 0; i < seq_.size(); i++ ) multi_[i] = max( (float)1, multi_[i] / max( (float)1, cover_[i] ) );
    for ( int i = 0; i < seq_.size(); i++ ) cover_[i] /= max( (float)1, multi_[i] );
}

void Review::setGaps( string l, string r )
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

vector<int> Review::setLibs( vector<Review*>& reviews )
{
    vector<int> libs;
    for ( int i = 0; i < params.libs.size(); i++ )
    {
        unordered_set<ReadId> used;
        vector<int> dists;
        for ( Review* r : reviews ) for ( ReviewPair& rp : r->pairs_[i] ) if ( !rp.dispute && used.insert( rp.id ).second ) dists.push_back( rp.dist() );
        sort( dists.begin(), dists.end() );
        assert( !dists.empty() );
        libs.push_back( ( dists[ ( dists.size()-1 ) / 2 ] + dists[ dists.size() / 2 ] ) / 2 );
    }
    return libs;
}

float Review::setMedian( vector<Review*>& reviews )
{
    vector<float> covers;
    
    for ( Review* r : reviews ) for ( pair<int, int> seed : r->seeds_ ) for ( int i = seed.first; i < seed.second; i++ ) covers.push_back( r->cover_[i] );
    
    sort( covers.begin(), covers.end() );
    
    if ( covers.empty() ) return 0;
    return ( covers[ ( covers.size()-1 ) / 2 ] + covers[ covers.size() / 2 ] ) / 2;
}

void Review::setMultis( vector<Review*>& reviews )
{
    int hCoord[2]{ len_, len_ }, pCoord[2]{ len_, len_ };
    for ( int i = 0; i < path_.size(); i++ )
    {
        int homo = 0;
        bool para = false;
        for ( Node* clone : path_[i]->clones( false ) ) for ( Review* r : reviews ) if ( find( r->path_.begin(), r->path_.end(), clone ) != r->path_.end() )
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

void Review::setOut( vector<int>& libs, ofstream* ofs, float median, int spaced, int cut )
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

void Review::setPairs( Querier& bwt, vector<Review*>& reviews )
{
    Lib* lib;
    ReadId id;
    int drxn;
    ReviewMap* paired[2],* alt[2];
    
//    auto within = []( vector< pair<int, int> >& homology, ReviewMap& rm )
//    {
//        for ( pair<int, int>& h : homology ) if ( h.first <= rm.coord[0] && rm.coord[1] <= h.second ) return true;
//        return false;
//    };
    
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
        
        for ( Review* r : reviews ) if ( r != this && ( ( alt[0] = r->get( m.first ) ) || ( alt[1] = r->get( id ) ) ) )
        {
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
        
//        // Only pursue right-facing pairs
//        if ( ( paired[!drxn] = &m.second ) && ( paired[drxn] = get( id ) ) && !drxn ) continue;
//        
//        vector<ReviewPair> pairs[4]; // [0] = Correct, [1] = too short distance, [2] = too long distance, [3] inverted orientation
//        int dupe = ReviewPair::add( pairs, paired[0], paired[1], lib->size );
//        
//        bool recombine = false, alted = false;
//        for ( Review* r : reviews ) if ( r != this && ( ( alt[0] = r->get( m.first ) ) || ( alt[1] = r->get( id ) ) ) )
//        {
//            if ( alt[drxn] && !alt[!drxn] ) recombine = true;
//            vector<ReviewPair> alts[4];
//            int duped = ReviewPair::add( alts, alt[0], alt[1], lib->size );
//            if ( alts[0].empty() ) continue;
//            if ( !pairs[0].empty() && ( pairs[0][0].getCut( lib->size ) < abs( alts[0][0].getDiff( lib->size ) ) ) ) continue;
//            if ( !pairs[0].empty() && alts[0][0].getCut( lib->size ) < abs( pairs[0][0].getDiff( lib->size ) ) ) for ( int i = 0; i < pairs[0].size(); i++ )
//            {
//                assert( false );
//                pairs[ pairs[0][i].dist() < lib->size ? 1 : 2 ].push_back( pairs[0][i] );
//                pairs[0].erase( pairs[0].begin() + i-- );
//            }
//            alted = true;
//        }
//        
//        if ( recombine && !alted ) mispair_[!drxn][iLib].push_back( m.first );
//        if ( !recombine && !alted && !paired[drxn] ) unpair_[!drxn][iLib].push_back( m.first );
//        
//        for ( int d : { 0, 1 } ) for ( int i = 1; i < 4; i++ ) for ( int j = 0; j < pairs[i].size(); j++ )
//        {
//            bool used = false;
//            for ( int k = 0; k < j; k++ ) if ( pairs[i][j].coord[d] == pairs[i][k].coord[d] ) used = true;
//            if ( used ) continue;
//            ( i == 3 ? inverts_ : ( i == 2 ? longs_ : shorts_ ) )[d][ pairs[0].empty() ][iLib].push_back( make_pair( d == drxn ? id : m.first, pairs[i][j].coord[d] ) );
//        }
//        
//        pairs_[iLib].insert( pairs_[iLib].end(), pairs[0].begin(), pairs[0].end() );
    }
    
    return;
    
    vector<int> unpaired[2]{ vector<int>( params.libs.size(), 0 ), vector<int>( params.libs.size(), 0 ) };
    vector<int> repaired[2]{ vector<int>( params.libs.size(), 0 ), vector<int>( params.libs.size(), 0 ) };
    vector<int> inverts[2]{ vector<int>( params.libs.size(), 0 ), vector<int>( params.libs.size(), 0 ) };
    vector<int> errs[2]{ vector<int>( 3, 0 ),  vector<int>( 3, 0 ) };
    for ( int i = 2; i < params.libs.size(); i++ ) for ( int d : { 0, 1 } ) for ( ReadId base : unpair_[d][i] )
    {
        ReadId pairId = params.getPairId( base );
        int32_t coords[2];
        string seq = bwt.getSequence( pairId );
        ReviewMap* rm = get( base );
        int best = 12, off = 0;
        bool bad = true, inverted = false, good = false;
        for ( int j = 0; bad && j < 3 && j * 12 + best < seq.size(); j++ )
        {
            string q = d ? seq.substr( j ) : seq.substr( 0, seq.size() - j );
            if ( !mapSeqEnd( q, seq_, best, coords, !d ) ) continue;
            size_t qq[2]{ d ? j : seq.size() - j, ( d ? j : seq.size() - j ) + coords[1] - coords[0] };
            vector<int> offs;
            for ( size_t it = seq_.find( q ); it != string::npos; )
            {
                size_t tt[2]{ it, it+qq[1]-qq[0] };
                while ( qq[0] && tt[0] && seq[ qq[0] ] && seq_[ tt[0] ] && qq[0]-- && tt[0]-- ) offs.clear();
                while ( qq[1] < seq.size() && tt[1] < seq_.size() && seq[ qq[1] ] && seq_[ tt[1] ] && qq[1]++ && tt[1]++ ) offs.clear();
                offs.push_back( tt[0] );
                it = seq_.find( q, it+1 );
            }
            assert( !offs.empty() );
            for ( int32_t off : offs ) if ( off <= qq[0] ) inverted = true;
            bad = false;
//            repaired[d][i]++;
//            if ( d ? coords[0] < rm->coord[0] : rm->coord[1] <= coords[1] ) inverts[d][i]++;
//            errs[d][j]++;
        }
        
        if ( bad ) unpaired[d][i]++;
    }
    
    assert( false );
}

void Review::setReview( Querier& bwt, vector<Review*>& reviews )
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
    
    for ( Review* r : reviews ) r->setSeeds( seeds );
    for ( Review* r : reviews ) r->setCover( reviews );
    for ( Review* r : reviews ) r->setPairs( bwt, reviews );
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
    for ( Review* r : reviews ) if ( !r->seeds_.empty() )
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

void Review::setSeeds( vector<string> seeds[3] )
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

void Review::setSpan( vector<ReviewPair>& pairs, int base, int lib, int& cur )
{
    while ( !pairs.empty() && pairs.back().coord[1] <= base ) pairs.pop_back();
    for ( ; cur < pairs_[lib].size() && pairs_[lib][cur].coord[0] <= base; cur++ )
    {
        int i = pairs.size();
        for ( int j = 0; j < i; j++ ) if ( pairs_[lib][cur].coord[1] > pairs[j].coord[1] ) i = j;
        pairs.insert( pairs.begin()+i, pairs_[lib][cur] );
    }
}

void Review::setSpanned( vector<int>& counts, int l, int r )
{
    assert( pairs_.size() == counts.size() );
    for ( int i = 0; i < pairs_.size(); i++ ) for ( ReviewPair& rp : pairs_[i] ) if ( rp.coord[0] < l && r < rp.coord[1] ) counts[i]++;
}

