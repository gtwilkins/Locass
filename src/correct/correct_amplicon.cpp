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

uint64_t Amp::ampCount = 0;
uint64_t Pile::pileCount = 0;
uint64_t Pile::adpCount = 0;
uint64_t Pile::trimCount = 0;
uint64_t Pile::created = 0;
uint64_t Pile::deleted = 0;

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
    bool failed = !id--;
    if ( failed )
    {
        assert( false );
        fread( seqs[0], 1, 200, ifp );
        return 0;
    }
    for ( int i = 0; i < 2; i++ )
    {
        fread( &seqs[i], 1, lens[i], ifp );
        for ( int j = 0; j < 16; j++ ) kmer = ( kmer << 2 ) + ( seqs[i][j] / 63 );
        amps[i] = new Amp( id, lens[i], seqs[i] );
    }
    amps[0]->alt = amps[1];
    amps[1]->alt = amps[0];
//    if ( kmer == 433093095525958445 ) cout << id << endl;
    
    return kmer;
}

int Amp::penalty( int i )
{
    if ( i >= phred.size() ) return 0;
    int pen = min( 20, (int) phred[i] );
    pen = ( pen * min( 40, i ) ) / 40;
    return pen;
}

void Amp::setUsed( uint8_t* used )
{
    if ( !used ) return;
    uint8_t flag = 1 << ( 7 - ( id % 8 ) );
    if ( used[id/8] & flag ) return;
    used[id/8] |= flag;
    ampCount++;
}

Pile::Pile( vector<Pile*> &piles, vector<Amp*> &amps, IndexReader* ir )
: amps_( amps )
{
    created++;
    seqLen_ = goodLen_ = 0;
    assert( ir );
    amps.clear();
    assert( !amps_.empty() && amps_[0]->seq.size() > 1 );
    for ( int i = 0; i < 12; i++ ) q_.push_back( amps_[0]->comp( i ) );
    int it = ir->setBaseAll( q_, rank_, count_ );
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

Pile::Pile( vector<Pile*> &piles, Pile* pile, int i, CharId rank, CharId count, IndexReader* ir )
: rank_( rank ), count_( count ), q_( pile->q_.begin(), pile->q_.end() - 1 )
{
    created++;
    seqLen_ = goodLen_ = 0;
    q_.push_back( pile->amps_[i]->comp( q_.size() ) );
    add( pile->amps_[i] );
    for ( ; i < pile->amps_.size(); i++ )
    {
        if ( pile->amps_[i]->queryLen < q_.size() - 1 || pile->amps_[i]->comp( q_.size()-1 ) != q_.back() ) continue;
        add( pile->amps_[i--] );
    }
    goodLen_ = count > 2 ? q_.size() : min( pile->goodLen_, (int)q_.size() - 1 );
    pair_ = NULL;
    
    query( piles, ir );
    for ( Amp* amp : amps_ ) seqLen_ = max( seqLen_, (int)amp->ascii.size() );
}

Pile::Pile( Pile* basis, Pile* pair )
: rank_( basis->rank_ ), count_( basis->count_ ), q_( basis->q_ ), goodLen_( basis->goodLen_ )
{
    created++;
    seqLen_ = goodLen_ = 0;
    for ( int i = 0; i < basis->amps_.size(); i++ ) if ( basis->amps_[i]->alt->pile == pair ) add( basis->amps_[i--] );
}

Pile::~Pile()
{
    deleted++;
    for ( Amp* amp : amps_ ) delete amp;
    amps_.clear();
    if ( pair_ && pair_->pair_ == this ) pair_->pair_ = NULL;
}

void Pile::add( Amp* amp )
{
    if ( amp->pile ) amp->pile->take( amp, amps_ );
    amp->pile = this;
    amp->queryLen = 0;
    for ( int &i = amp->queryLen; i < q_.size(); i++ ) if ( amp->comp( i ) != q_[i] ) break;
    seqLen_ = max( seqLen_, (int)amp->ascii.size() );
}

void Pile::advance( uint8_t c, CharId rank, CharId count )
{
    q_.push_back( c );
    rank_ = rank;
    count_ = count;
    if ( count > 2 ) goodLen_ = q_.size();
}

void Pile::cannabalise( Pile* pile )
{
    while ( !pile->amps_.empty() )
    {
        assert( pile->amps_[0]->pile == pile );
        add( pile->amps_[0] );
    }
}

void Pile::clean( vector<Pile*> &piles )
{
    for ( int i = 0; i < piles.size(); i++ )
    {
        if ( piles[i]->amps_.size() > 1 && piles[i]->pair_ && piles[i]->pair_->amps_.size() > 1 ) continue;
        delete piles[i];
        piles.erase( piles.begin() + i-- );
    }
}

void Pile::collapse( vector<Pile*> &piles )
{
    if ( piles.size() < 2 ) return;
    sort( piles.begin(), piles.end(), []( Pile* a, Pile* b ) { return a->goodLen_ > b->goodLen_; } );
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

vector<Pile*> Pile::compile( vector<Amp*> (&amps)[2], IndexReader* ir )
{
    vector<Pile*> piles[2];
    for ( int i : { 0, 1 } ) piles[i].push_back( new Pile( piles[i], amps[i], ir ) );
    
    for ( int i : { 0, 1 } ) Pile::collapse( piles[i] );
    
    for ( int i : { 0, 1 } ) for ( Pile* ap : piles[i] ) ap->mate();

    for ( int i : { 0, 1 } ) for ( Pile* ap : piles[i] ) ap->confirm( piles[!i]);
    
    for ( int i : { 0, 1 } ) Pile::clean( piles[i] );
    
    return piles[0];
}

void Pile::confirm( vector<Pile*> &alts )
{
    mate();
    if ( pair_ ) pair_->mate();
    if ( !pair_ || pair_->pair_ == this ) return;
    vector<Pile*> piles[2];
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
    
    alts.push_back( new Pile( pair_, this ) );
    alts.back()->mate();
    mate();
}

void Pile::demate( IndexReader* ir )
{
    if ( amps_.size() < 2 ) return;
    Pile* piles[2] = { goodLen_ > pair_->goodLen_ ? pair_ : this, goodLen_ > pair_->goodLen_ ? this : pair_ };
    vector<string> seqs[2] = { piles[0]->getSeqs(), piles[1]->getSeqs() };
    bool trimmed[2]{ false, false };
    for ( int i = 0; i < 2; i++ ) CorrectQuery cq( ir, piles[i]->seq_, seqs[i], piles[i]->branches_, piles[i]->goodLen_, trimmed[i] );
    for ( int i = 0; i < 2; i++ )
    {
        if ( trimmed[i] || !piles[i]->goodLen_ || piles[i]->goodLen_ == piles[i]->seq_.size() ) continue;
        if ( trimmed[!i] && ( piles[i]->goodLen_ > piles[i]->seq_.size() * .9 ) && seqs[i].size() == 1 ) continue;
        CorrectQuery cq( ir, piles[i]->seq_, seqs[i], piles[i]->goodLen_, trimmed[i] );
    }
    
    for ( int i = 0; i < 2; i++ ) piles[i]->fill( trimmed[i] ? piles[i]->goodLen_ : piles[i]->seqLen_ );
    if ( trimmed[0] || trimmed[1] ) Pile::trimCount++;
}

void Pile::fill( int trimLen )
{
    if ( seq_.size() > trimLen ) seq_.erase( seq_.begin() + trimLen, seq_.end() );
    if ( seq_.size() < seqLen_ ) seq_ += string( seqLen_ - seq_.size(), 'N' );
}

Amp* Pile::getBest()
{
    Amp* best = NULL;
    for ( Amp* amp : amps_ )
    {
        if ( best && best->queryLen > amp->queryLen ) continue;
        if ( !best || best->queryLen < amp->queryLen || best->phredLen < amp->phredLen ) best = amp;
    }
    return best;
}

int Pile::getDisagree( Pile* alt )
{
    int len = min( alt->q_.size(), q_.size() );
    for ( int i = 0; i < len; i++ ) if ( q_[i] != alt->q_[i] ) return i;
    return len;
}

void Pile::getPairs( vector<Pile*> &piles )
{
    for ( Amp* amp : amps_ )
    {
        if ( !amp->alt->pile ) continue;
        if ( find( piles.begin(), piles.end(), amp->alt->pile ) == piles.end() ) piles.push_back( amp->alt->pile );
    }
}

string Pile::getPhred()
{
    string phred;
    for ( int i = 0; i < seq_.size(); i++ ) phred += ( i < goodLen_ ? char( 33 + 40 ) : char( 33 + 2 ) );
    return phred;
}

int Pile::getPhred( int i )
{
    uint8_t phred = 0;
    for ( Amp* amp : amps_ ) if ( amp->queryLen > i && amp->comp( i ) < 4 ) phred = max( phred, amp->phred[i] );
    return phred;
}

Pile* Pile::getPile( Amp* amp, vector<Pile*> &piles )
{
    for ( Pile* ap : piles ) if ( ap->inPile( amp ) ) return ap;
    assert( false );
    return NULL;
}

vector<string> Pile::getSeqs()
{
    vector<string> seqs;
    seq_.clear();
    for ( uint8_t q : q_ ) seq_ += intToCharComp[q];
    goodLen_ = seq_.size();
    
    for ( Amp* amp : amps_ )
    {
        bool added = amp->ascii.find( seq_ ) != 0;
        for ( string &seq : seqs ) if ( seq == amp->ascii ) added = true;
        if ( !added ) seqs.push_back( amp->ascii );
    }
    
    if ( seq_.size() < seqLen_ ) seq_ += string( seqLen_ - seq_.size(), 'N' );
    
    return seqs;
}

bool Pile::inPile( Amp* amp )
{
    return find( amps_.begin(), amps_.end(), amp ) != amps_.end();
}

void Pile::mate()
{
    Amp* best = NULL;
    for ( Amp* amp : amps_ ) if ( !best || amp->queryLen + amp->alt->queryLen > best->queryLen + best->alt->queryLen ) best = amp;
    pair_ = best ? best->alt->pile : NULL;
}

bool Pile::merge( Pile* alt )
{
    Amp* bests[2] = { getBest(), alt->getBest() };
    assert( bests[0] && bests[1] );
    
    return bests[0]->align( bests[1] );
}

void Pile::output( ofstream &ofs, uint8_t* used, Deadapter* da, IndexReader* ir )
{
    pileCount++;
    bool adapter = false;
    for ( Amp* amp : amps_ )
    {
        amp->setUsed( used );
        if ( !amp->alt ) continue;
        string phred[2];
        for ( int i = 0; i < amp->phred.size(); i++ ) phred[0] += char( amp->phred[i] + 33 );
        for ( int i = 0; i < amp->alt->phred.size(); i++ ) phred[1] += char( amp->alt->phred[i] + 33 );
        if ( da && da->isOverlap( amp->ascii, amp->alt->ascii, phred[0], phred[1], false ) ) adapter = true;
    }
    for ( Amp* amp : pair_->amps_ ) amp->setUsed( used );
    
    if ( adapter ) adpCount++;
    else
    {
        demate( ir );
        assert( ofs );
        ofs << "@\n" << seq_ << "\n+\n" << getPhred() << "\n";
        ofs << "@\n" << pair_->seq_ << "\n+\n" << pair_->getPhred() << "\n";
    }
    delete pair_;
    delete this;
}

void Pile::query( vector<Pile*> &piles, IndexReader* ir  )
{
    CharCount ranks, counts;
    ir->countRange( q_.back(), rank_, count_, ranks, counts );
    
    int len = q_.size();
    for ( int i = 0; i < amps_.size(); i++ )
    {
        uint8_t c = amps_[i]->comp( len );
        if ( amps_[i]->queryLen < len || c > 3 || !counts[c] ) continue;
        amps_[i]->queryLen = len + 1;
        if ( q_.size() == len ) advance( c, ranks[c], counts[c] );
        else if ( q_.back() != c ) piles.push_back( new Pile( piles, this, i--, ranks[c], counts[c], ir ) );
    }
    
    if ( len < q_.size() ) query( piles, ir );
    else for ( int i = 0; i < 4; i++ ) if ( counts[i] ) branches_.push_back( CorrectBranch( i, len, ranks[i], counts[i] ) );
}

void Pile::take( Amp* amp, vector<Amp*> &amps )
{
    assert( find( amps_.begin(), amps_.end(), amp ) != amps_.end() );
    amps_.erase( remove( amps_.begin(), amps_.end(), amp ), amps_.end() );
    amps.push_back( amp );
    amp->pile = NULL;
}

