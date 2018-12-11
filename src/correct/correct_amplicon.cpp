/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "correct_amplicon.h"
#include <cassert>
#include <algorithm>
#include <string.h>
#include <valarray>
#include "constants.h"
#include "index_reader.h"
#include "local_alignment.h"
#include <iostream>

//int Pile::adapterCount_ = 0;
//uint32_t Pile::pileCount_ = 0;
//uint32_t Pile::ampCount_ = 0;
uint64_t AmpPile::pileCount;
uint64_t AmpPile::ampCount;
uint64_t AmpPile::adpCount;
uint64_t AmpPile::trimCount;

Amp::Amp( uint32_t id, uint8_t len, uint8_t* buff )
: id( id ), pile( NULL )
{
    queryLen = phredLen = 0;
    for ( int i = 0; i < len; i++ )
    {
        seq.push_back( buff[i] / 63 );
        phred.push_back( seq.back() > 3 ? 2 : buff[i] % 63 );
        ascii.push_back( intToChar[ seq.back() ] );
        if ( phred.back() > 4 ) phredLen = i + 1;
    }
}

Amp::~Amp()
{
    if ( alt && alt->alt == this ) alt->alt = NULL;
}

bool Amp::align( Amp* amp )
{
    string s[2] = { ascii.substr( 0, max( queryLen, phredLen ) )
                  , amp->ascii.substr( 0, max( amp->queryLen, amp->phredLen ) ) };
    string a[2];
    LocalAlignment la( s[0], s[1], false, true );
    la.realign( a[0], a[1], true, false, false, true );
    
    int j = 0, allMiss = 0, queryMiss = 0;
//    char c = 'X';
    for ( int i = 0; i < a[1].size(); i++ )
    {
        for ( int d = 0; d < 2; d++ )
        {
            if ( a[d][i] != '-' ) continue;
            int gap = 1;
            while ( i + gap < a[d].size() && a[d][i+gap] == '-' ) gap++;
            int poly = LocalAlignment::isGapPoly( a, d, i, gap );
            for ( ; --gap >= 0; i++ )
            {
                if ( poly < 1 ) allMiss += amp->penalty( j );
                if ( poly-- < 1 && j < amp->queryLen ) queryMiss += amp->penalty( j );
                if ( !d ) j++;
            }
            if ( i == a[d].size() ) break;
            d = -1;
        }
        if ( i >= a[1].size() ) break;
        if ( a[0][i] != a[1][i] )
        {
            allMiss += amp->penalty( j );
            if ( j < amp->queryLen ) queryMiss += amp->penalty( j );
        }
        j++;
    }
    
    return queryMiss < 40 && allMiss < 80;
}

uint8_t Amp::comp( int i )
{
    return i < ascii.size() ? charToIntComp[ ascii[i] ] : 4;
}

uint64_t Amp::create( FILE* ifp, Amp* (&amps)[2] )
{
    uint64_t kmer;
    uint32_t id;
    uint8_t lens[2], seqs[2][256];
    fread( &id, 4, 1, ifp );
    fread( lens, 1, 2, ifp );
    for ( int i = 0; i < 2; i++ )
    {
        fread( &seqs[i], 1, lens[i], ifp );
        for ( int j = 0; j < 16; j++ ) kmer = ( kmer << 2 ) + ( seqs[i][j] / 63 );
        amps[i] = new Amp( id, lens[i], seqs[i] );
    }
    amps[0]->alt = amps[1];
    amps[1]->alt = amps[0];
    
    return kmer;
}

int Amp::penalty( int i )
{
    if ( i >= phred.size() ) return 0;
    int pen = min( 20, (int) phred[i] );
    pen = ( pen * min( 40, i ) ) / 40;
    return pen;
}

AmpPile::AmpPile( vector<AmpPile*> &piles, vector<Amp*> &amps, IndexReader* ir )
: amps_( amps )
{
    seqLen_ = goodLen_ = 0;
    assert( ir );
    amps.clear();
    assert( !amps_.empty() && amps_[0]->seq.size() > 1 );
    for ( int i = 0; i < 12; i++ ) seq_.push_back( amps_[0]->comp( i ) );
    int it = ir->setBaseAll( seq_, rank_, count_ );
    assert( it == 12 );
    for ( Amp* amp : amps_ )
    {
        amp->queryLen = it;
        amp->pile = this;
    }
    goodLen_ = it;
    
//    seq_.push_back( amps_[0]->comp( 0 ) );
//    seq_.push_back( amps_[0]->comp( 1 ) );
//    ir->setBaseAll( seq_[0], seq_[1], rank_, count_ );
//    for ( Amp* amp : amps_ )
//    {
//        amp->queryLen = 2;
//        amp->pile = this;
//    }
//    goodLen_ = 2;
    pair_ = NULL;
    
    query( piles, ir );
    for ( Amp* amp : amps_ ) seqLen_ = max( seqLen_, (int)amp->ascii.size() );
}

AmpPile::AmpPile( vector<AmpPile*> &piles, AmpPile* pile, int i, CharId rank, CharId count, IndexReader* ir )
: rank_( rank ), count_( count ), seq_( pile->seq_.begin(), pile->seq_.end() - 1 )
{
    seqLen_ = goodLen_ = 0;
    seq_.push_back( pile->amps_[i]->comp( seq_.size() ) );
    add( pile->amps_[i] );
    for ( ; i < pile->amps_.size(); i++ )
    {
        if ( pile->amps_[i]->queryLen < seq_.size() - 1 || pile->amps_[i]->comp( seq_.size()-1 ) != seq_.back() ) continue;
        add( pile->amps_[i--] );
    }
    goodLen_ = count > 2 ? seq_.size() : min( pile->goodLen_, (int)seq_.size() - 1 );
    pair_ = NULL;
    
    query( piles, ir );
    for ( Amp* amp : amps_ ) seqLen_ = max( seqLen_, (int)amp->ascii.size() );
}

AmpPile::AmpPile( AmpPile* basis, AmpPile* pair )
: rank_( basis->rank_ ), count_( basis->count_ ), seq_( basis->seq_ ), goodLen_( basis->goodLen_ )
{
    seqLen_ = goodLen_ = 0;
    for ( Amp* amp : basis->amps_ ) if ( amp->alt->pile == pair ) add( amp );
}

AmpPile::~AmpPile()
{
    for ( Amp* amp : amps_ ) delete amp;
    amps_.clear();
    if ( pair_ && pair_->pair_ == this ) pair_->pair_ = NULL;
}

void AmpPile::add( Amp* amp )
{
    if ( amp->pile ) amp->pile->take( amp, amps_ );
    amp->pile = this;
    amp->queryLen = 0;
    for ( int &i = amp->queryLen; i < seq_.size(); i++ ) if ( amp->comp( i ) != seq_[i] ) break;
    seqLen_ = max( seqLen_, (int)amp->ascii.size() );
}

void AmpPile::advance( uint8_t c, CharId rank, CharId count )
{
    seq_.push_back( c );
    rank_ = rank;
    count_ = count;
    if ( count > 2 ) goodLen_ = seq_.size();
}

void AmpPile::cannabalise( AmpPile* pile )
{
    while ( !pile->amps_.empty() )
    {
        assert( pile->amps_[0]->pile == pile );
        add( pile->amps_[0] );
    }
}

void AmpPile::clean( vector<AmpPile*> &piles )
{
    for ( int i = 0; i < piles.size(); i++ )
    {
        if ( piles[i]->amps_.size() > 1 && piles[i]->pair_ && piles[i]->pair_->amps_.size() > 1 ) continue;
        delete piles[i];
        piles.erase( piles.begin() + i-- );
    }
}

void AmpPile::collapse( vector<AmpPile*> &piles )
{
    if ( piles.size() < 2 ) return;
    sort( piles.begin(), piles.end(), []( AmpPile* a, AmpPile* b ) { return a->goodLen_ > b->goodLen_; } );
    for ( int i = 0; i < piles.size(); i++ )
    {
        if ( piles[i]->amps_.empty() ) continue;
        for ( int j = i+1; j < piles.size(); j++ )
        {
            if ( piles[i]->getDisagree( piles[j] ) < piles[j]->goodLen_ ) continue;
            for ( Amp* amp : piles[j]->amps_ ) piles[i]->add( amp );
        }
    }
}

vector<AmpPile*> AmpPile::compile( vector<Amp*> (&amps)[2], IndexReader* ir, QueryBinaries* qb )
{
    vector<AmpPile*> piles[2];
    for ( int i : { 0, 1 } ) piles[i].push_back( new AmpPile( piles[i], amps[i], ir ) );
    
    for ( int i : { 0, 1 } ) AmpPile::collapse( piles[i] );
    
    for ( int i : { 0, 1 } ) for ( AmpPile* ap : piles[i] ) ap->mate();

    for ( int i : { 0, 1 } ) for ( AmpPile* ap : piles[i] ) ap->confirm( piles[!i]);
    
    for ( int i : { 0, 1 } ) AmpPile::clean( piles[i] );
    
    return piles[0];
}

void AmpPile::confirm( vector<AmpPile*> &alts )
{
    mate();
    if ( pair_ ) pair_->mate();
    if ( !pair_ || pair_->pair_ == this ) return;
    vector<AmpPile*> piles[2];
    piles[0].push_back( pair_->pair_ );
    piles[0].push_back( this );
    piles[1].push_back( pair_ );
    int sizes[2]{ 1, 1 };
    while ( sizes[0] < piles[0].size() && sizes[1] < piles[1].size() )
    {
        int cur[2]{ (int)piles[0].size(), (int)piles[1].size() };
        for ( int d = 0; d < 2; d++ )
        {
            for ( int i = sizes[d]; i < piles[d].size(); i++ ) piles[d][i]->getPairs( piles[!d] );
            sizes[d] = cur[d];
        }
    }
    bool good = true;
    for ( int i = 1; i < piles[0].size(); i++ ) if ( !piles[0][0]->merge( piles[0][i] ) ) good = false;
    for ( int i = 1; i < piles[1].size(); i++ ) if ( !piles[1][0]->merge( piles[1][i] ) ) good = false;
    
    if ( good )
    {
        for ( int i = 1; i < piles[0].size(); i++ ) piles[0][0]->cannabalise( piles[0][i] );
        for ( int i = 1; i < piles[1].size(); i++ ) piles[1][0]->cannabalise( piles[1][i] );
        return;
    }
    
    Amp* best = pair_->pair_->getBest();
    for ( int i = 0; i < amps_.size(); i++ )
    {
        if ( !best->align( amps_[i] ) ) continue;
        pair_->pair_->add( amps_[i--] );
        good = true;
    }
    
    if ( good )
    {
        mate();
        confirm( alts );
        return;
    }
    
    alts.push_back( new AmpPile( pair_, this ) );
    alts.back()->mate();
    mate();
}

void AmpPile::demate( IndexReader* ir, QueryBinaries* qb )
{
    if ( amps_.size() < 2 ) return;
    AmpPile* piles[2] = { goodLen_ > pair_->goodLen_ ? pair_ : this, goodLen_ > pair_->goodLen_ ? this : pair_ };
    vector<string> seqs[2] = { piles[0]->getSeqs(), piles[1]->getSeqs() };
    bool trimmed[2]{ seqs[0].empty(), seqs[1].empty() };
    for ( int i = 0; i < 2; i++ )
    {
        if ( trimmed[i] )
        {
            piles[i]->goodLen_ = 0;
            piles[i]->ascii_ = string( seqLen_, 'N' );
            for ( Amp* amp : piles[i]->amps_ )
            {
                if ( amp->queryLen <= piles[i]->goodLen_ ) continue;
                piles[i]->ascii_ = amp->ascii;
                piles[i]->goodLen_ = amp->queryLen;
            }
        }
        else
        {
            piles[i]->ascii_ = seqs[i][0];
            CorrectQuery cq( ir, piles[i]->seq_, piles[i]->branches_, piles[i]->seqLen_ );
            piles[i]->goodLen_ = cq.correct( piles[i]->ascii_, seqs[i], trimmed[i] );
        }
    }
    for ( int i = 0; i < 2; i++ )
    {
        if ( trimmed[i] || !piles[i]->goodLen_ || piles[i]->goodLen_ == piles[i]->ascii_.size() ) continue;
        if ( trimmed[!i] && ( piles[i]->goodLen_ > piles[i]->ascii_.size() * .9 ) && seqs[i].size() == 1 ) continue;
        CorrectQuery cq( ir, piles[i]->ascii_, piles[i]->goodLen_, true );
        piles[i]->goodLen_ = cq.trim( qb, piles[i]->ascii_, trimmed[i] );
    }
    for ( int i = 0; i < 2; i++ ) piles[i]->fill( trimmed[i] ? piles[i]->goodLen_ : piles[i]->seqLen_ );
    if ( trimmed[0] || trimmed[1] ) AmpPile::trimCount++;
}

void AmpPile::fill( int trimLen )
{
    if ( ascii_.size() > trimLen ) ascii_.erase( ascii_.begin() + trimLen, ascii_.end() );
    if ( ascii_.size() < seqLen_ ) ascii_ += string( seqLen_ - ascii_.size(), 'N' );
}

Amp* AmpPile::getBest()
{
    Amp* best = NULL;
    for ( Amp* amp : amps_ )
    {
        if ( best && best->queryLen > amp->queryLen ) continue;
        if ( !best || best->queryLen < amp->queryLen || best->phredLen < amp->phredLen ) best = amp;
    }
    return best;
}

int AmpPile::getDisagree( AmpPile* alt )
{
    int len = min( alt->seq_.size(), seq_.size() );
    for ( int i = 0; i < len; i++ ) if ( seq_[i] != alt->seq_[i] ) return i;
    return len;
}

void AmpPile::getPairs( vector<AmpPile*> &piles )
{
    for ( Amp* amp : amps_ )
    {
        if ( !amp->alt->pile ) continue;
        if ( find( piles.begin(), piles.end(), amp->alt->pile ) == piles.end() ) piles.push_back( amp->alt->pile );
    }
}

string AmpPile::getPhred()
{
    string phred;
    for ( int i = 0; i < ascii_.size(); i++ ) phred += ( i < goodLen_ ? char( 33 + 40 ) : char( 33 + 2 ) );
    return phred;
}

int AmpPile::getPhred( int i )
{
    uint8_t phred = 0;
    for ( Amp* amp : amps_ ) if ( amp->queryLen > i && amp->comp( i ) < 4 ) phred = max( phred, amp->phred[i] );
    return phred;
}

AmpPile* AmpPile::getPile( Amp* amp, vector<AmpPile*> &piles )
{
    for ( AmpPile* ap : piles ) if ( ap->inPile( amp ) ) return ap;
    assert( false );
    return NULL;
}

vector<string> AmpPile::getSeqs()
{
    vector<string> seqs;
    int maxCount = 0;
    vector<bool> added( amps_.size(), false );
    for ( int i = 0; i < amps_.size(); i++ )
    {
        if ( added[i] || amps_[i]->queryLen < seq_.size() ) continue;
        added[i] = true;
        int cloneCount = 1;
        for ( int j = i+1; j < amps_.size(); j++ )
        {
            if ( amps_[i]->ascii != amps_[j]->ascii ) continue;
            added[j] = true;
            cloneCount++;
        }
        if ( cloneCount < maxCount ) continue;
        if ( cloneCount > maxCount ) seqs.size();
        seqs.push_back( amps_[i]->ascii );
        maxCount = cloneCount;
    }
    
    return seqs;
}

bool AmpPile::inPile( Amp* amp )
{
    return find( amps_.begin(), amps_.end(), amp ) != amps_.end();
}

void AmpPile::mate()
{
    Amp* best = NULL;
    for ( Amp* amp : amps_ ) if ( !best || amp->queryLen + amp->alt->queryLen > best->queryLen + best->alt->queryLen ) best = amp;
    pair_ = best ? best->alt->pile : NULL;
}

bool AmpPile::merge( AmpPile* alt )
{
    Amp* bests[2] = { getBest(), alt->getBest() };
    assert( bests[0] && bests[1] );
    
    return bests[0]->align( bests[1] );
}

void AmpPile::output( ofstream &ofs, uint8_t* used, Deadapter* da, IndexReader* ir, QueryBinaries* qb )
{
    pileCount++;
    ampCount += amps_.size();
    bool adapter = false;
    for ( Amp* amp : amps_ )
    {
        if ( used ) used[amp->id/8] |= 1 << ( 7 - ( amp->id % 8 ) );
        if ( !amp->alt ) continue;
        string phred[2];
        for ( int i = 0; i < amp->phred.size(); i++ ) phred[0] += char( amp->phred[i] + 33 );
        for ( int i = 0; i < amp->alt->phred.size(); i++ ) phred[1] += char( amp->alt->phred[i] + 33 );
        if ( da && da->isOverlap( amp->ascii, amp->alt->ascii, phred[0], phred[1] ) ) adapter = true;
    }
    for ( Amp* amp : pair_->amps_ ) if ( used ) used[amp->id/8] |= 1 << ( 7 - ( amp->id % 8 ) );
    
    if ( adapter ) adpCount++;
    else
    {
        demate( ir, qb );
        assert( ofs );
        ofs << "@\n" << ascii_ << "\n+\n" << getPhred() << "\n";
        ofs << "@\n" << pair_->ascii_ << "\n+\n" << pair_->getPhred() << "\n";
    }
    delete pair_;
    delete this;
}

void AmpPile::query( vector<AmpPile*> &piles, IndexReader* ir  )
{
    CharCount ranks, counts;
    ir->countRange( seq_.back(), rank_, count_, ranks, counts );
    
    int len = seq_.size();
    for ( int i = 0; i < amps_.size(); i++ )
    {
        uint8_t c = amps_[i]->comp( len );
        if ( amps_[i]->queryLen < len || c > 3 || !counts[c] ) continue;
        amps_[i]->queryLen = len + 1;
        if ( seq_.size() == len ) advance( c, ranks[c], counts[c] );
        else if ( seq_.back() != c ) piles.push_back( new AmpPile( piles, this, i--, ranks[c], counts[c], ir ) );
    }
    
    if ( len < seq_.size() ) query( piles, ir );
    else for ( int i = 0; i < 4; i++ ) if ( counts[i] ) branches_.push_back( CorrectBranch( i, len, ranks[i], counts[i] ) );
}

void AmpPile::take( Amp* amp, vector<Amp*> &amps )
{
    assert( find( amps_.begin(), amps_.end(), amp ) != amps_.end() );
    amps_.erase( remove( amps_.begin(), amps_.end(), amp ), amps_.end() );
    amps.push_back( amp );
    amp->pile = NULL;
}

//Amplicon::Amplicon( uint8_t* buff, Amplicon* ampPair )
//: pair_( ampPair )
//{
//    score_ = strength_ = perfect_ = goodLen_ = 0;
//    memcpy( &id_, &buff[0], 4 );
//    int limit = pair_ ? 6 + buff[4] : 6 + buff[4] + buff[5];
//    for ( int i = pair_ ? 6 : 6 + buff[4]; i < limit; i++ )
//    {
//        seq_.push_back( buff[i] / 63 );
//        assert( seq_.back() < 5 );
//        phred_.push_back( seq_.back() > 3 ? 0 : buff[i] % 63 );
//        ascii_.push_back( !seq_.back() ? 'A' : ( seq_.back() == 1 ? 'C' : ( seq_.back() == 2 ? 'G' : ( seq_.back() == 3 ? 'T' : 'N' ) ) ) );
//        if ( phred_.back() > 4 ) goodLen_ = phred_.size();
//        if ( phred_.back() > 4 && phred_.size() >= 16 ) strength_ += min( 20, (int)phred_.back() );
//        assert( ascii_.size() > 16 || ascii_.back() != 'N' );
//    }
//    
//    if ( pair_ ) ampPair->pair_ = this;
//    assert( !ascii_.empty() );
//}
//
//bool Amplicon::allSame( vector<Amplicon*> &amps, int i )
//{
//    for ( int j = 1; j < amps.size(); j++ ) if ( i < amps[j]->seq_.size() && amps[0]->seq_[i] != amps[j]->seq_[i] ) return false;
//    return amps.size() > 1 && amps[1]->seq_.size() < i;
//}
//
//int Amplicon::getCharScore( int i )
//{
//    return i < phred_.size() ? min( 10, phred_[i] / 2 ) : 0;
//}
//
//int Amplicon::getPairScore()
//{
//    return score_ + pair_->score_ - ( miss_ * 2 );
//}
//
//void Amplicon::shift( uint8_t* bases, int pos, int off )
//{
//    assert( pos <= seq_.size() );
//    if ( off > 0 )
//    {
//        seq_.insert( seq_.begin() + pos, bases, bases + off );
//        phred_.insert( phred_.begin() + pos, off, 2 );
//    }
//    
//    if ( off < 0 )
//    {
//        seq_.erase( seq_.begin() + pos + off, seq_.begin() + pos );
//        phred_.erase( phred_.begin() + pos + off, phred_.begin() + pos );
//    }
//}
//
//bool Offset::clean( vector<Amplicon*> &good )
//{
//    for ( int i = 0; i < amps.size(); i++ ) if ( find( good.begin(), good.end(), amps[i] ) == good.end() ) amps.erase( amps.begin() + i-- );
//    return !amps.empty();
//}
//
//bool Offset::set( vector<Offset> &offs, vector<Amplicon*> &good )
//{
//    for ( int i = 0; i < offs.size(); i++ ) if ( !offs[i].clean( good ) ) offs.erase( offs.begin() + i-- );
//    if ( offs.empty() ) return false;
//    
//    unordered_set<Amplicon*> usedAmps;
//    for ( Offset &o : offs ) usedAmps.insert( o.amps.begin(), o.amps.end() );
//    offs.push_back( Offset() );
//    for ( Amplicon* amp : good ) if ( usedAmps.find( amp ) == usedAmps.end() ) offs.back().amps.push_back( amp );
//    
//    return true;
//}
//
//bool Shift::add( int i, Amplicon* amp, int off )
//{
//    if ( i != pos ) return false;
//    for ( Offset &o : offs ) o.amps.erase( remove( o.amps.begin() , o.amps.end(), amp ), o.amps.end() );
//    for ( Offset &o : offs )
//    {
//        if ( off != o.off ) continue;
//        o.amps.push_back( amp );
//        return true;
//    }
//    offs.push_back( Offset( off, amp ) );
//    return true;
//}
//
//int Shift::resolve( vector<uint8_t> &seq, vector<Amplicon*> &amps, int off, int &posLimit )
//{
//    // Clean up
//    pos += off;
//    off = 0;
//    if ( !Offset::set( offs, amps ) ) return 0;
//    
//    // Compile evidence for each offset
//    double numor = 0, denom = 0;
//    int totalNum = 0, lastNum = 0, limits[2]{0};
//    for ( Offset &o : offs )
//    {
//        limits[0] = min( limits[0], o.off );
//        limits[1] = max( limits[1], o.off );
//    }
//    int baseScores[4][ limits[1] - limits[0] ] = { 0 };
//    memset( baseScores, 0, sizeof( baseScores ) );
//    for ( Offset &o : offs )
//    {
//        vector< tuple< vector<uint8_t>, int, int> > scores;
//        int bestScore = 0, bestNum = 1;
//        for ( Amplicon* amp : o.amps )
//        {
//            int runScore = amp->phred_[pos];
//            vector<uint8_t> run( amp->seq_.begin() + pos + limits[0], amp->seq_.begin() + pos + o.off + 1 );
//            for ( int i = pos + limits[0]; i < pos + o.off + 1; i++ ) runScore = min( runScore, (int)amp->phred_[i] );
//            for ( int i = pos + limits[0]; i < pos + o.off; i++ ) if ( amp->seq_[i] < 4 ) baseScores[ amp->seq_[i] ][ i - pos - limits[0] ] += amp->phred_[i];
//            bool added = false;
//            for ( tuple< vector<uint8_t>, int, int> &score : scores )
//            {
//                if ( get<0>(score) != run ) continue;
//                get<1>(score) = max( get<1>(score), runScore );
//                bestNum = max( 1, ++get<2>(score) );
//                added = true;
//            }
//            if ( !added ) scores.push_back( make_tuple( run, runScore, 1 ) );
//        }
//        for ( tuple< vector<uint8_t>, int, int> &score : scores ) bestScore = max( bestScore, get<1>(score) * get<2>(score) );
//        numor += bestScore * o.off;
//        denom += bestScore;
//        totalNum += get<2>(scores[0]);
//        lastNum = get<2>(scores[0]);
//    }
//    double est = numor / denom;
//    
//    // Set the insertable bases
//    uint8_t bases[ limits[1] - limits[0] ]{0};
//    for ( int i = 0; i < limits[1] - limits[0]; i++ )  for ( int j = 1; j < 4; j++ ) if ( baseScores[i][j] > baseScores[i][ bases[i] ] ) bases[i] = j;
//    
//    for ( Offset &o : offs ) if ( abs( est - off ) > abs( est - o.off ) ) off = o.off;
//    if ( totalNum - lastNum <= lastNum || pos > posLimit ) off = 0;
//    for ( Offset &o : offs ) for ( Amplicon* amp : o.amps ) amp->shift( bases, pos + o.off, off - o.off );
//    
//    if ( off > 0 ) seq.insert( seq.begin() + pos, bases, bases + off );
//    if ( off < 0 ) seq.erase( seq.begin() + pos + off, seq.begin() + pos );
//    
//    posLimit += off;
//    return off;
//}
//
//Pile::Pile( Amplicon* amp )
//: optimised_( true )
//{
//    amps_.push_back( amp );
//    build();
//}
//
//Pile::Pile( vector<Amplicon*> &amps, int len )
//: amps_( amps.begin(), amps.end() ), optimised_( false )
//{
//    int bestScore = 0;
//    for ( int i = 0; i < amps_.size(); i++ )
//    {
//        int score = 0;
//        for ( int j = len; j < amps_[i]->seq_.size(); j++ ) score += amps_[i]->getCharScore( j );
//        if ( score <= bestScore ) continue;
//        bestScore = score;
//        if ( i ) iter_swap( amps_.begin(), amps_.begin() + i );
//    }
//    build();
//}
//
//Pile::~Pile()
//{
//    for ( Amplicon* amp : amps_ ) assert( false );
//    for ( Amplicon* amp : amps_ ) delete amp;
//}
//
//void Pile::add( vector<Amplicon*> &amps, bool doSteal )
//{
//    for ( int i = 0; i < amps.size(); i++ )
//    {
//        if ( inPile( amps[i] ) || align( amps[i], 40, doSteal, false, true ) ) amps.erase( amps.begin() + i-- );
//    }
//}
//
//bool Pile::align( Amplicon* amp, int maxMiss, bool steal, bool force, bool add )
//{
//    int score = 0, miss = 0, perfect = amp->seq_.size();
//    int len = alignExact( amp, score, miss, perfect );
//    if ( maxMiss >= 0 && miss > maxMiss ) return false;
//    
//    vector< pair<int, int> > offs;
//    if ( len < min( seq_.size(), amp->seq_.size() ) && !alignInexact( len, amp, score, miss, perfect, maxMiss, offs ) ) return false;
//    if ( !force && ( steal ? miss > amp->miss_ : score <= amp->score_ ) ) return false;
//    
//    amp->score_ = score;
//    amp->miss_ = miss;
//    amp->perfect_ = perfect;
//    if ( add ) amps_.push_back( amp );
//    if ( add && perfect_ < seq_.size() ) optimised_ = false;
//    assert( !amps_.empty() );
//    if ( amp != amps_[0] ) perfect_ = max( perfect_, perfect );
//    if ( amp != amps_[0] && score > score_ )
//    {
//        score_ = score;
//        miss_ = miss;
//    }
//    
//    for ( pair<int, int> &off : offs )
//    {
//        bool added = false;
//        for ( Shift &s : shifts_ ) if ( !added ) added = s.add( off.first, amp, off.second );
//        if ( !added ) continue;
//        int i = 0;
//        while ( i < shifts_.size() && shifts_[i].pos < off.first ) i++;
//        shifts_.insert( shifts_.begin() + i, Shift( off.first, amp, off.second ) );
//    }
//    
//    return true;
//}
//
//int Pile::alignExact( Amplicon* amp, int &score, int &miss, int &perfect )
//{
//    int limit = min( seq_.size(), amp->seq_.size() ), len = 16, tmpScore = 0, tmpMiss = 0, last[2]{ 16, 16 };
//    for ( int i = 16; i < limit; i++ )
//    {
//        if ( seq_[i] > 3 || amp->seq_[i] > 3 ) continue;
//        int phred = max( 4, min( 20, (int)amp->phred_[i] ) );
//        if ( seq_[i] == amp->seq_[i] )
//        {
//            if ( amp->seq_[i] != amp->seq_[i-1] )
//            {
//                last[0]++;
//                last[1]++;
//            }
//            tmpScore += phred;
//        }
//        else
//        {
//            if ( last[0] > 10 )
//            {
//                score = tmpScore;
//                miss = tmpMiss;
//                perfect = min( perfect, i );
//                len = i;
//            }
//            if ( last[1] <= 10 ) return len;
//            last[1] = last[0];
//            last[0] = 0;
//            tmpScore -= phred;
//            tmpMiss += phred;
//        }
//    }
//    score = tmpScore;
//    miss = tmpMiss;
//    return limit;
//}
//
//bool Pile::alignInexact( int i, Amplicon* amp, int &score, int &miss, int &perfect, int maxMiss, vector< pair<int, int> > &offs )
//{
//    int mid = 2, base = i-1;
//    while ( base > 0 && seq_[base] == seq_[base-1] ) base--;
//    for ( const vector<uint8_t> &s : { seq_, amp->seq_ } )
//    {
//        int run = 1;
//        for ( int i = base; i < s.size(); i++ )
//        {
//            run = s[i] == s[i-1] && s[i] < 4 ? run+1 : 0;
//            mid = min( 4, max( mid, run / 2 ) );
//        }
//    }
//    mid = max( 0, min( mid, (int)min( seq_.size(), amp->seq_.size() ) - i - 2 ) );
//    int maxGap = mid * 2 + 1, best = mid;
//    int curScores[maxGap]{0}, nxtScores[maxGap]{0}, curMisses[maxGap]{0}, nxtMisses[maxGap]{0}, curPerfect[maxGap], nxtPerfect[maxGap];
//    int polys[maxGap], poly = seq_[i-1] < 4;
//    vector< pair<int,int> > curOffs[maxGap], nxtOffs[maxGap];
//    for ( int j = i-1; --j >= 0 && seq_[j] == seq_[j+1]; ) poly++;
//    for ( int j = 0; j < maxGap; j++ )
//    {
//        int k = i - 1 + j - mid;
//        curScores[j] = nxtScores[j] = score - ( 80 * abs( j - mid ) );
//        curMisses[j] = nxtMisses[j] = miss + ( k < amp->goodLen_ ? 30 * abs( j - mid ) : 0 );
//        curPerfect[j] = nxtPerfect[j] = ( j == mid ? ( i == perfect ? seq_.size() : perfect ) : min( i, perfect ) );
//        if ( j != mid ) curOffs[j].push_back( make_pair( i, j - mid ) );
//        polys[j] = amp->seq_[k] < 4;
//        for ( ; --k >= 0 && amp->seq_[k] == amp->seq_[k+1]; ) polys[j]++;
//    }
//    
//    for ( ; i < seq_.size(); i++ )
//    {
//        for ( int j = 0; j < maxGap; j++ )
//        {
//            int k = i + j - mid, diffs[2]{0}, bestDiff = 0;
//            if ( k >= amp->seq_.size() ) continue;
//            
//            // Determine maximum free shift distances in either direction
//            if ( seq_[i-1] == amp->seq_[k-1] && ( seq_[i] != seq_[i-1] || amp->seq_[k] != amp->seq_[k-1] ) )
//            {
//                diffs[0] = min( j, polys[j] - 1 );
//                diffs[1] = min( maxGap - j, poly - 1 );
//            }
//            
//            // Score allowed shifts
//            for ( int m = j-diffs[0]; m < j + diffs[1]; m++ )
//            {
//                if ( m == j || curScores[m] - 40 <= nxtScores[j] ) continue;
//                nxtScores[j] = curScores[m] - 40;
//                nxtMisses[j] = curMisses[m];
//                nxtPerfect[j] = curPerfect[m];
//                bestDiff = m - j;
//            }
//            
//            // Score penalised gap
//            for ( int m = max( 0, j-1 ); m < min( j+2, maxGap ); m++ )
//            {
//                if ( m == j || curScores[m] - 80 <= nxtScores[j] ) continue;
//                nxtScores[j] = curScores[m] - 80;
//                nxtMisses[j] = curMisses[m] + ( k < amp->goodLen_ ? 30 : 0 );
//                nxtPerfect[j] = min( curPerfect[m], i );
//                bestDiff = m - j;
//            }
//            
//            // Add score and penalty for this base
//            if ( seq_[i] > 3 || seq_[i] != amp->seq_[k] ) nxtPerfect[j] = min( nxtPerfect[j], i );
//            if ( seq_[i] < 4 && amp->seq_[k] < 4 )
//            {
//                int match = ( seq_[i] == amp->seq_[k] ? 1 : - 1 ) * max( 4, min( 20, (int)amp->phred_[k] ) );
//                nxtScores[j] += match ;
//                if ( k < amp->goodLen_ && match < 0 ) nxtMisses[j] -= match;
//            }
//            
//            // Record a shift or gap
//            if ( bestDiff )
//            {
//                nxtOffs[j] = curOffs[j+bestDiff];
//                nxtOffs[j].push_back( make_pair( i, -bestDiff ) );
//            }
//            polys[j] = amp->seq_[k] < 4 && amp->seq_[k] == amp->seq_[k-1] ? polys[j]+1 : 0;
//        }
//        
//        bool good = maxMiss < 0;
//        for ( int j = 0; j < maxGap; j++ )
//        {
//            if ( nxtMisses[j] <= maxMiss ) good = true;
//            if ( nxtScores[j] > nxtScores[best] ) best = j;
//            if ( !nxtOffs[j].empty() ) curOffs[j] = nxtOffs[j];
//            nxtOffs[j].clear();
//        }
//        if ( !good ) return false;
//        memcpy( &curScores, &nxtScores, sizeof( int ) * maxGap );
//        memcpy( &curMisses, &nxtMisses, sizeof( int ) * maxGap );
//        memcpy( &curPerfect, &nxtPerfect, sizeof( int ) * maxGap );
//        poly = seq_[i] < 4 && seq_[i] == seq_[i-1] ? poly + 1 : 0;
//    }
//    
//    score = curScores[best];
//    miss = curMisses[best];
//    perfect = curPerfect[best];
//    offs = curOffs[best];
//    
//    return true;
//}
//
//void Pile::build()
//{
//    assert( !amps_.empty() );
//    seq_ = amps_[0]->seq_;
//    score_ = miss_ = perfect_ = 0;
//    pair_ = NULL;
//    shifts_.clear();
//    
//    for ( Amplicon* amp : amps_ )
//    {
//        amp->score_ = amp->miss_ = amp->perfect_ = 0;
//        align( amp, -1, false, true, false );
//        if ( amp == amps_[0] || amp->score_ < score_ ) continue;
//        score_ = amp->score_;
//        miss_ = amp->miss_;
//    }
//}
//
//bool Pile::close()
//{
//    if ( !amps_.empty() ) return false;
//    delete this;
//    return true;
//}
//
//vector<string> Pile::compile( vector<Amplicon*> (&amps)[2] )
//{
//    vector<Pile*> piles[2] = { Pile::group( amps[0] ), Pile::group( amps[1] ) };
//    
//    bool again = true;
//    while ( again )
//    {
//        again = false;
//        
//        for ( int i : { 0, 1 } ) Pile::optimise( piles[i] );
//        
//        // Set the best pair for each pile
//        for ( int i : { 0, 1 } ) for ( int j = 0; j < piles[i].size(); j++ ) piles[i][j]->couple( piles[!i] );
//
//        // Resolve any secondary pairings
//        for ( int i : { 0, 1 } ) for ( int j = 0; j < piles[i].size(); j++ ) if ( !piles[i][j]->couple( piles[i], piles[!i] ) ) again = true;
//        
//        // Resolve any four-way amplicon pairings
//        for ( int i : { 0, 1 } ) for ( int j = 0; !again && j < piles[i].size(); j++ ) again = piles[i][j]->disentagle( piles[!i], amps[i], amps[!i] );
//        
//        // Remove any defunct piles
//        for ( int i : { 0, 1 } ) for ( int j = 0; j < piles[i].size(); j++ ) if ( piles[i][j]->close() ) piles[i].erase( piles[i].begin() + j-- );
//    }
//    
//    // Remove any pile pairing that doesn't include two or more pairs in support
//    vector<string> lines;
//    for ( int i = 0; i < piles[0].size(); )
//    {
//        if ( piles[0][i]->confirm() ) lines.push_back( piles[0][i++]->output() );
//        else
//        {
//            piles[1].erase( remove( piles[1].begin(), piles[1].end(), piles[0][i]->pair_ ), piles[1].end() );
//            piles[0][i]->dismantle( amps[0] );
//            piles[0][i]->pair_->dismantle( amps[1] );
//            piles[0].erase( piles[0].begin() + i );
//        }
//    }
//    assert( piles[0].size() == piles[1].size() && amps[0].size() == amps[1].size() );
//    return lines;
//}
//
//bool Pile::confirm()
//{
//    assert( pair_ && pair_->pair_ == this && amps_.size() == pair_->amps_.size() );
//    if ( amps_.size() < 2 ) return false;
//    
//    Pile* pairs[2] = { this, pair_ };
//    vector<Amplicon*> candidates[2];
//    int good[2]{0}, perf[2]{0};
//    bool bad[2] = { false, false };
//    for ( int i = 0; i < 2; i++ )
//    {
//        for ( Amplicon* amp : pairs[i]->amps_ ) assert( pairs[!i]->inPile( amp->pair_ ) );
//        for ( Amplicon* amp : pairs[i]->amps_ ) good[i] = max( good[i], amp->goodLen_ );
//        good[i] = min( good[i], int( pairs[i]->seq_.size() *.8) );
//        for ( Amplicon* amp : pairs[i]->amps_ )
//        {
//            if ( amp->perfect_ >= good[i] ) candidates[i].push_back( amp );
//            else if ( amp->perfect_ < min( amp->goodLen_, good[i] ) ) bad[i] = true;
//        }
//        assert( !candidates[i].empty() );
//        sort( candidates[i].begin(), candidates[i].end(), []( Amplicon* a, Amplicon* b ){ return a->perfect_ > b->perfect_; } );
//        pairs[i]->perfect_ = candidates[i].size() > 1 ? candidates[i][1]->perfect_ : 0;
//        perf[i] = pairs[i]->perfect_;
//    }
//    
//    if ( ( bad[0] || bad[1] ) && candidates[0].size() < 2 && candidates[1].size() < 2 ) return false;
//    
//    for ( Amplicon* a : candidates[0] ) for ( Amplicon* b : candidates[1] ) if ( a->pair_ == b ) return true;
//    
//    // If at least one pile is perfect
//    for ( int i = 0; i < 2; i++ )
//    {
//        if ( pairs[!i]->perfect_ < pairs[!i]->seq_.size() ) continue;
//        for ( Amplicon* amp : candidates[i] )
//        {
//            if ( amp->pair_->miss_ <= ( 40 * amp->pair_->goodLen_ ) / pairs[!i]->perfect_ )
//            {
//                setPerfect();
//                return true;
//            }
//        }
//    }
//    
//    return false;
//}
//
//bool Pile::couple ( vector<Pile*> &piles )
//{
//    Amplicon* best = NULL;
//    pair_ = NULL;
//    for ( Amplicon* amp : amps_ ) if ( !best || amp->getPairScore() > best->getPairScore() ) best = amp;
//    assert( best && !pair_ );
//    for ( Pile* p : piles )
//    {
//        if ( p->inPile( best->pair_ ) )
//        {
//            assert( !pair_ );
//            pair_ = p;
//        }
//    }
//    assert( pair_ );
//    return pair_;
//}
//
//bool Pile::couple ( vector<Pile*> &curPiles, vector<Pile*> &altPiles )
//{
//    if ( amps_.empty() ) return true;
//    if ( !pair_ || !pair_->pair_ ) return false;
//    if ( pair_->pair_ == this ) return true;
//    
//    // Attempt to transfer all amplicons from this pile to the preferred pile
//    for ( int i = 0; i < amps_.size(); i++ )
//    {
//        if ( !pair_->inPile( amps_[i]->pair_ ) ) continue;
//        if ( pair_->pair_->align( amps_[i], 80, false, true, true ) ) amps_.erase( amps_.begin() + i-- );
//    }
//    
//    if ( amps_.empty() ) return false;
//    
//    if ( rebuild() ) return false;
//    
//    vector<Amplicon*> amps;
//    pair_->extract( amps, amps_ );
//    pair_->rebuild();
//    altPiles.push_back( new Pile( amps, 0 ) );
//    
//    return false;
//}
//
//bool Pile::disentagle( vector<Pile*> &altPiles, vector<Amplicon*> &curAmps, vector<Amplicon*> &altAmps )
//{
//    assert( pair_ );
//    for ( int i = 0; i < amps_.size(); i++ )
//    {
//        if ( pair_->inPile( amps_[i]->pair_ ) ) continue;
//        Pile* alt = NULL;
//        for ( Pile* pile : altPiles ) if ( pile->inPile( amps_[i]->pair_ ) ) alt = pile;
//        assert( alt );
//        if ( alt->merge( pair_ ) + merge( alt->pair_ ) ) return true;
//        bool altFirst = amps_[i]->score_ + pair_->score_ < amps_[i]->pair_->score_ + alt->score_;
//        int added = 2;
//        if ( ( altFirst ? alt->pair_ : pair_ )->align( altFirst ? amps_[i] : amps_[i]->pair_, 80, false, true, true ) ) added = altFirst;
//        else if ( ( altFirst ? pair_ : alt->pair_ )->align( altFirst ? amps_[i]->pair_ : amps_[i], 80, false, true, true ) ) added = !altFirst;
//        else
//        {
//            curAmps.push_back( amps_[i] );
//            altAmps.push_back( amps_[i]->pair_ );
//        }
//        if ( added != 1 )
//        {
//            alt->amps_.erase( remove( alt->amps_.begin(), alt->amps_.end(), amps_[i]->pair_ ), alt->amps_.end() );
//            alt->rebuild();
//        }
//        if ( added != 0 )
//        {
//            amps_.erase( amps_.begin() + i-- );
//            rebuild();
//        }
//        return true;
//    }
//    
//    return false;
//}
//
//void Pile::dismantle( vector<Amplicon*> &amps )
//{
//    amps.insert( amps.end(), amps_.begin(), amps_.end() );
//    amps_.clear();
//    delete this;
//}
//
//void Pile::extract( vector<Amplicon*> &curAmps, vector<Amplicon*> &altAmps )
//{
//    for ( int i = 0; i < amps_.size(); i++ )
//    {
//        if ( find( altAmps.begin(), altAmps.end(), amps_[i]->pair_ ) == altAmps.end() ) continue;
//        curAmps.push_back( amps_[i] );
//        amps_.erase( amps_.begin() + i-- );
//    }
//}
//
//bool Pile::fill( vector<Amplicon*> &amps )
//{
//    for ( int i = 0; i < amps_.size(); i++ )
//    {
//        auto it = find( amps.begin(), amps.end(), amps_[i] );
//        if ( it != amps.end() ) amps.erase( it );
//        else amps_.erase( amps.begin() + i-- );
//    }
//    
//    if ( !rebuild( ) ) return true;
//    assert( false );
//    amps.insert( amps.end(), amps_.begin(), amps_.end() );
//    sort( amps.begin(), amps.end(), []( Amplicon* a, Amplicon* b ){ return a->strength_ < b->strength_; } );
//    delete this;
//    return false;
//}
//
//vector<Pile*> Pile::group( vector<Amplicon*> &amps )
//{
//    vector<Pile*> piles;
//    int goodLen = 0;
//    for ( Amplicon* amp : amps ) goodLen = max( goodLen, amp->goodLen_ );
//    sort( amps.begin(), amps.end(), []( Amplicon* a, Amplicon* b ){ return a->seq_.size() > b->seq_.size(); } );
//    Pile::group( piles, amps, goodLen, 16 );
//    
//    sort( piles.begin(), piles.end(), []( Pile* a, Pile* b ) { return a->score_ > b->score_; } );
//    sort( amps.begin(), amps.end(), []( Amplicon* a, Amplicon* b ){ return a->strength_ < b->strength_; } );
//    
//    for ( int i = 0; i < piles.size() || !amps.empty(); i++ )
//    {
//        // Discard any piles that have been cannabalised
//        while ( i < piles.size() && !piles[i]->fill( amps ) ) piles.erase( piles.begin() + i );
//        
//        // Create a bew pile from the highest quality amplicon
//        if ( i == piles.size() )
//        {
//            if ( amps.empty() ) break;
//            piles.push_back( new Pile( amps.back() ) );
//            amps.pop_back();
//            for ( Pile* p : piles ) if ( p != piles.back() ) piles.back()->add( p->amps_, false );
//        }
//        
//        // Merge in any matching amplicons to this pile
//        piles[i]->add( amps, false );
//    }
//    
//    return piles;
//}
//
//void Pile::group( vector<Pile*> &piles, vector<Amplicon*> &amps, int goodLen, int i )
//{
//    assert( amps.size() > 1 );
//    while ( Amplicon::allSame( amps, i ) ) i++;
//    if ( i >= goodLen || i >= amps[1]->seq_.size() )
//    {
//        piles.push_back( new Pile( amps, i ) );
//        return;
//    }
//    vector<Amplicon*> tmpAmps[4];
//    for ( Amplicon* amp : amps ) if ( i < amp->seq_.size() && amp->seq_[i] < 4 ) tmpAmps[ amp->seq_[i] ].push_back( amp );
//    for ( int j = 0; j < 4; j++ ) if ( tmpAmps[j].size() > 1 ) Pile::group( piles, tmpAmps[j], goodLen, i+1 );
//}
//
//bool Pile::inPile( Amplicon* amp )
//{
//    return find( amps_.begin(), amps_.end(), amp ) != amps_.end();
//}
//
//bool Pile::merge( Pile* pile )
//{
//    if ( pile->amps_.size() == 1 && align( pile->amps_[0], 80, false, true, true ) )
//    {
//        pile->amps_.clear();
//        return true;
//    }
//    return false;
//}
//
//void Pile::optimise( vector<Pile*> &piles )
//{
//    for ( bool again = true; again; )
//    {
//        again = false;
//        for ( Pile* pile : piles )
//        {
//            if ( !pile->optimise() ) continue;
//            for ( int i = 0; i < piles.size(); i++ )
//            {
//                if ( piles[i] != pile ) pile->add( piles[i]->amps_, false );
//                if ( !piles[i]->close() ) continue;
//                assert( false );
//                piles.erase( piles.begin() + i-- );
//                again = true;
//            }
//        }
//        if ( piles.size() < 2 ) again = false;
//    }
//}
//
//bool Pile::optimise()
//{
//    if ( perfect_ == seq_.size() || amps_.size() < 3 || optimised_ ) return false;
//    int len = 0;
//    for ( Amplicon* amp : amps_ ) len = max( len, min( (int)seq_.size(), amp->goodLen_ ) );
//    
//    double scores[len][5];
//    scoreConsensus( scores, len );
//    
//    Amplicon* best = NULL;
//    double bestScore = -1;
//    for ( Amplicon* amp : amps_ )
//    {
//        double score = 0;
//        int j = 0, pos = 0, off = 0;
//        setNextOffset( amp, pos, off );
//        for ( int i = 0; i < len; i++ )
//        {
//            if ( off && i == min( pos, pos + off ) )
//            {
//                assert( false );
//                j += abs( off );
//                setNextOffset( amp, pos, off );
//            }
//            if ( j >= amp->seq_.size() ) break;
//            if ( amp->seq_[j] < 4 ) score += scores[i][ amp->seq_[j] ] * scores[i][ amp->seq_[j] ] / scores[i][4];
//            j++;
//        }
//        if ( score <= bestScore ) continue;
//        best = amp;
//        bestScore = score;
//    }
//    
//    assert( best );
//    if ( seq_ != best->seq_ )
//    {
//        assert( best && best->perfect_ < seq_.size() );
//        seq_ = best->seq_;
//        rebuild();
//        optimised_ = true;
//        return true;
//    }
//    assert( best && best->perfect_ == best->seq_.size() );
//    
//    return false;
//}
//
//string Pile::output()
//{
//    int ampCount = amps_.size();
//    Pile::pileCount_++;
//    Pile::ampCount_ += ampCount;
//    string line;
//    for ( int i = 0; i < 2; i++ )
//    {
//        string seq, phred;
//        ( i ? this : pair_ )->output( seq, phred );
//        line += "@" + to_string( Pile::pileCount_ ) + "." + to_string( i ) + " merged_" + to_string( ampCount ) + "_amplicons\n";
//        line += seq + "\n";
//        line += "+\n";
//        line += phred + "\n";
//    }
//    
//    delete pair_;
//    delete this;
//    return line;
//}
//
//void Pile::output( string &seq, string &phred )
//{
////    setPerfect();
//    int off = 0, len = seq_.size();
//    for ( Shift &shift : shifts_ ) off += shift.resolve( seq_, amps_, off, perfect_ );
//    
//    if ( !perfect_ ) perfect_ = seq_.size();
//    for ( int i = 0; i < len; i++ )
//    {
//        if ( i < seq_.size() && i < perfect_ ) seq.push_back( !seq_[i] ? 'A' : ( seq_[i] == 1 ? 'C' : ( seq_[i] == 2 ? 'G' : ( seq_[i] == 3 ? 'T' : 'N' ) ) ) );
//        else seq.push_back( 'N' );
//        uint8_t phredScore = 2;
//        for ( Amplicon* amp : amps_ ) if ( seq.back() != 'N' && i < amp->seq_.size() && seq_[i] == amp->seq_[i]  ) phredScore = max( phredScore, amp->phred_[i] );
//        phredScore += 33;
//        phred.push_back( char( phredScore ) );
//    }
//    
//    for ( Amplicon* amp : amps_ ) delete amp;
//    amps_.clear();
//}
//
//bool Pile::rebuild()
//{
//    for ( int i = 0; i < amps_.size(); i++ )
//    {
//        if ( amps_[i]->seq_ != seq_ ) continue;
//        if ( i )
//        {
//            iter_swap( amps_.begin(), amps_.begin() + i );
//            build();
//        }
//        return false;
//    }
//    
//    assert( !amps_.empty() );
//    
//    sort( amps_.begin(), amps_.end(), []( Amplicon* a, Amplicon* b ){ return a->strength_ > b->strength_; } );
//    build();
//    
//    return true;
//}
//
//void Pile::scoreConsensus( double scores[][5], int len )
//{
//    for ( int i = 0; i < len; i++ ) for ( int j = 0; j < 5; j++ ) scores[i][j] = 0;
//    
//    for ( Amplicon* amp : amps_ )
//    {
//        int j = 0, pos = 0, off = 0;
//        setNextOffset( amp, pos, off );
//        for ( int i = 0; i < len; i++ )
//        {
//            if ( off && i == min( pos, pos + off ) )
//            {
//                assert( false );
//                j += abs( off );
//                setNextOffset( amp, pos, off );
//            }
//            if ( j >= amp->seq_.size() ) break;
//            if ( amp->seq_[j] < 4 )
//            {
//                scores[i][ amp->seq_[j] ] += amp->phred_[j];
//                scores[i][ 4 ] += amp->phred_[j];
//            }
//            j++;
//        }
//    }
//}
//
//void Pile::setNextOffset( Amplicon* amp, int &pos, int &off )
//{
//    for ( Shift &s : shifts_ )
//    {
//        if ( s.pos <= pos ) continue;
//        for ( Offset &o : s.offs )
//        {
//            if ( find( o.amps.begin(), o.amps.end(), amp ) == o.amps.end() ) continue;
//            assert( false );
//            pos = s.pos;
//            off = o.off;
//            return;
//        }
//    }
//    pos = seq_.size();
//    off = 0;
//}
//
//void Pile::setPerfect()
//{
//    if ( perfect_ || amps_.empty() ) return;
//    
//    double scores[seq_.size()][5];
//    scoreConsensus( scores, seq_.size() );
//    int len = amps_[0]->goodLen_;
//    for ( int i = 0; i < seq_.size(); i++ )
//    {
//        if ( seq_[i] < 4 && ( i < len ? scores[i][ seq_[i] ] * 2 > scores[i][4] : scores[i][ seq_[i] ] == max( 1.0, scores[i][4] ) ) ) continue;
//        perfect_ = i;
//        return;
//    }
//}
