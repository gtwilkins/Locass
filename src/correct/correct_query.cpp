/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "correct_query.h"
#include <cassert>
#include "constants.h"
#include "local_alignment.h"
#include <algorithm>

uint64_t CorrectAlign::nCorrectCount = 0;
uint64_t CorrectAlign::errorCorrectCount = 0;
uint64_t CorrectAlign::polyCorrectCount = 0;

CorrectAlign::CorrectAlign( string &t, string &q )
{
    LocalAlignment la( t, q, false, true );
    la.realign( a[0], a[1], false, true, true );
    lens[0] = t.size();
    lens[1] = q.size();
    hit = miss = limit = poly = perf = score = maxScore = 0;
    limit = a[0].size();
    for ( int d = 0; limit == a[0].size() && d < 2; d++ ) while ( limit > 0 && a[d][limit-1] == '-' ) limit--;
    for ( int i = 0; i < limit; i++ )
    {
        for ( int d = 0; d < 2 && i < limit; d++ )
        {
            if ( a[d][i] != '-' ) continue;
            int gap = 1;
            while ( ++i < limit && a[d][i] == '-' ) gap++;
            if ( i < limit ) scoreGap( !d, i - gap, gap );
            else miss += gap;
            d = -1;
        }
        if ( i >= limit ) break;
        if ( a[0][i] == a[1][i] ) hit++;
        else if ( a[0][i] != 'N' ) miss++;
        if ( !miss && !poly ) perf = i + 1;
        maxScore = max( maxScore, hit - ( miss * 4 ) - poly );
    }
    score = hit - ( miss * 4 ) - poly;
}

bool CorrectAlign::contendFull( CorrectAlign &full )
{
    if ( full > *this ) return false;
    int i = 0, j = 0, len = min( lens[0], min( lens[1], full.lens[1] ) );
    for ( int k = 0; k < len; k++ )
    {
        while ( i < full.a[1].size() && full.a[1][i] == '-' ) i++;
        while ( j < a[1].size() && a[1][j] == '-' ) j++;
        assert( i < full.a[1].size() && j < a[1].size() );
        if ( full.a[1][i] != a[1][j] ) return true;
    }
    
    return false;
}

bool CorrectAlign::contendPart( CorrectAlign &part, int &partLen )
{
    if ( !perf ) return *this > part;
    int len = min( partLen, perf );
    for ( int i = 0; i < min( partLen, len ); i++ ) if ( part.a[1][i] != a[1][i] ) partLen = i;
    return !partLen;
}

vector<CorrectAlign> CorrectAlign::get( vector<string> &ts, vector<string> &qs, bool drxn )
{
    vector<CorrectAlign> cas;
    if ( !drxn ) for ( string &seq : qs ) reverse( seq.begin(), seq.end() );
    if ( !drxn ) for ( string &seq : ts ) reverse( seq.begin(), seq.end() );
    for ( string &q : qs ) for ( string &t : ts ) cas.push_back( CorrectAlign( t, q ) );
    return cas;
}

bool CorrectAlign::isBad()
{
    return miss > 2 && score < 0;
}

bool CorrectAlign::isFull( int maxMiss )
{
    return lens[1] >= lens[0] && ( miss + poly <= maxMiss );
}

bool CorrectAlign::isPart()
{
    return perf;
}

bool CorrectAlign::operator >( CorrectAlign &rhs )
{
    if ( score > rhs.score && maxScore > rhs.maxScore ) return true;
    if ( score < rhs.score ) return false;
    int i = 0, j = 0;
    while ( i < a[1].size() && j < rhs.a[1].size() )
    {
        while ( i < a[1].size() && a[1][i] == '-' ) i++;
        while ( j < rhs.a[1].size() && rhs.a[1][j] == '-' ) j++;
        if ( i < a[1].size() && j < rhs.a[1].size() && a[1][i++] != rhs.a[1][j++] ) return false;
    }
    return true;
}

void CorrectAlign::scoreGap( int d, int i, int gap )
{
    int best = 0;
    for ( int len = 1; len <= gap && len*2 <= i; len++ )
    {
        bool same = len == 1;
        for ( int j = i-len+1; !same && j < i; j++ ) same = a[1][j] != a[1][j-1];
        for ( int j = i-len; same && j < i; j++ )
        {
            if ( a[0][j] != a[1][j] && a[0][j] != 'N' ) same = false;
            if ( len > 1 && a[1][j] != a[1][j-len] ) same = false;
            if ( a[1][j] != a[d][j+len] && a[d][j+len] != 'N' ) same = false;
        }
        if ( same ) best = max( best, len );
        for ( int k = len*2; same && k-len < gap; k += len )
        {
            for ( int j = i-len; same && j < i && a[!d][j+k] == '-'; j++ )
            {
                if ( a[1][j] != a[d][j+k] && a[d][j+k] != 'N' ) same = false;
                if ( same ) best = max( best, j+k-i+1 );
            }
        }
    }
    if ( best ) poly++;
    assert( best <= gap );
    miss += max( 0, gap - best );
}

void CorrectAlign::write( string &seq, int len, bool drxn )
{
    string base = ( drxn ? seq.substr( 0, seq.size() - lens[0] ) : seq.substr( lens[0] ) );
    string extend;
    for ( int i = 0; i < a[1].size() && extend.size() < lens[0]; i++ )
    {
        if ( extend.size() >= len && a[0][i] != '-' ) extend += a[0][i];
        if ( extend.size() < len && a[1][i] != '-' )
        {
            if ( a[0][i] == 'N' ) nCorrectCount++;
            else if ( a[0][i] == '-' ) polyCorrectCount++;
            else if ( a[0][i] != a[1][i]  ) errorCorrectCount++;
            extend += a[1][i];
        }
    }
    assert( extend.size() == lens[0] );
    seq = drxn ? base + extend : string( extend.rbegin(), extend.rend() ) + base;
}

CorrectBranch::CorrectBranch( uint8_t i, int it, CharId rank, CharId count )
: rank( rank ), count( count )
{
    q.push_back( i );
    base = it;
    bad = good = 0;
}

void CorrectBranch::align( string seq, int len, bool drxn )
{
    string ext[2]{ drxn ? seq.substr( seq.size() - len ) : seq.substr( 0, len ), "" };
    for ( int i = 0; i < q.size(); i++ ) ext[1] += drxn ? intToCharComp[ q[i] ] : intToChar[ q.end()[-i-1] ];
    aligns.push_back( CorrectAlign( ext[0], ext[1] ) );
}

void CorrectBranch::contend( CorrectAlign*& full, CorrectAlign*& part, int partLen, int imperf )
{
    for ( CorrectAlign &ca : aligns )
    {
        if ( full && full != &ca && ca.contendFull( *full ) ) full = NULL;
        if ( part && part != &ca && ca.contendPart( *part, partLen ) ) part = NULL;
        if ( !ca.perf ) imperf++;
    }
}

void CorrectBranch::correct( CorrectAlign*& full, CorrectAlign*& part, bool &bad, int maxMiss )
{
    for ( CorrectAlign &ca : aligns )
    {
        if ( ca.isBad() ) bad = true;
        if ( ( !full || ca > *full ) && ca.isFull( maxMiss ) ) full = &ca;
        if ( ( !part || ca.perf > part->perf ) && ca.isPart() ) part = &ca;
    }
}

bool CorrectBranch::proceed( IndexReader* ir, int limit )
{
    assert( branches.empty() && count );
    query( ir, limit );
    return true;
}

void CorrectBranch::query( IndexReader* ir, int limit )
{
    if ( q.size() >= limit ) return;
    
    CharCount ranks;
    CharCount counts;
    
    ir->countRange( q.back(), rank, count, ranks, counts );
    
    int i = 4;
    for ( int j = 0; j < 4; j++ ) if ( counts[j] && ( i == 4 || counts[j] > counts[i] ) ) i = j;
    for ( int j = 0; i < 4 && j < 4; j++ ) if ( j != i && ( counts[j] > min( (CharId)1, counts[i] / 8 ) ) ) j = 4;
    for ( int j = 0; j < 4; j++ ) if ( counts[j] && j != i ) branches.push_back( CorrectBranch( j, base+q.size(), ranks[j], counts[j] ) );
    if ( i > 3 ) return;
    
    q.push_back( i );
    rank = ranks[i];
    count = counts[i];
    query( ir, limit );
}

CorrectQuery::CorrectQuery( IndexReader* ir, string &seq, int &len, bool &trimmed, bool dummy )
: qLen_( 0 ), seqLen_( seq.size() ), trim_( false ), initial_( !len )
{
    // Set right query paramters
    coords_[0] = initial_ ? 0 : len - ( seqLen_ * .7 );
    coords_[1] = seqLen_;
    for ( int i = coords_[0]; i < seq.size() && i - coords_[0] < 32; i++ ) if ( seq[i] == 'N' ) coords_[0] = i + 1;
    if ( initial_ ? coords_[0] > 30 : coords_[0] <= 0 ) return;
    
    // Perform right query and correct or trim if necessary
    assert( setQuery( seq, 1 ) );
//    if ( !setQuery( seq, 0 ) ) return;
    CharId rank, count;
    int it = ir->setBaseAll( q_, rank, count );
    query( ir, q_[it-1], it, rank, count );
    len = correct( ir, seq, 1 ) + coords_[0];
    if ( trim_ ) trimmed = true;
    
    if ( !coords_[0] || !initial_ ) return;
    
    // Perform left query if undetermined bases present
    coords_[1] = min( coords_[0] + int( seqLen_ * .7 ), len );
    coords_[0] = 0;
    assert( setQuery( seq, 0 ) );
    it = ir->setBaseAll( q_, rank, count );
    query( ir, q_[it-1], it, rank, count );
    int left = coords_[1] - correct( ir, seq, 0 );
    for ( int i = 0; i < left; i++ ) if ( seq[i] != 'N' ) len = 0;
}

CorrectQuery::CorrectQuery( IndexReader* ir, string &seq, vector<string> &seqs, vector<CorrectBranch> &branches, int &len, bool &trimmed )
: branches_( branches ), qLen_( len ), seqLen_( seq.size() ), trim_( false ), initial_( false )
{
    coords_[0] = 0;
    coords_[1] = seqLen_;
    if ( branches_.empty() ) return;
    len = correct( ir, seq, seqs );
    if ( trim_ ) trimmed = true;
}

CorrectQuery::CorrectQuery( IndexReader* ir, string &seq, vector<string> &seqs, int &len, bool &trimmed )
: qLen_( len ), seqLen_( seq.size() ), trim_( false ), initial_( false )
{
    coords_[0] = len - ( seqLen_ * .7 );
    coords_[1] = seqLen_;
    if ( coords_[0] <= 0 ) return;
    
    assert( setQuery( seq, 1 ) );
    CharId rank, count;
    int it = ir->setBaseAll( q_, rank, count );
    query( ir, q_[it-1], it, rank, count );
    len = correct( ir, seq, seqs );
    if ( trim_ ) trimmed = true;
}

int CorrectQuery::correct( IndexReader* ir, string &seq, bool drxn )
{
    if ( proceed( ir ) )
    {
        for ( CorrectBranch &cb : branches_ ) cb.align( seq, q_.size() - cb.base, drxn );
        qLen_ += correct( seq, drxn );
    }
    
    branches_.clear();
    
    return qLen_;
}

int CorrectQuery::correct( IndexReader* ir, string &seq, vector<string> &seqs )
{
    if ( proceed( ir ) )
    {
        for ( string &s : seqs ) for ( CorrectBranch &cb : branches_ ) cb.align( s, seqLen_ - qLen_ - coords_[0], 1 );
        qLen_ += correct( seq, 1 );
    }
    if ( trim_ || qLen_ + coords_[0] >= seqLen_ ) seqs.clear();
    branches_.clear();
    
    string q = seq.substr( 0, coords_[0] + qLen_ );
    for ( int i = 0; i < seqs.size(); i++ ) if ( seqs[i].find( q ) != 0 ) seqs.erase( seqs.begin() + i-- );
    
    return coords_[0] + qLen_;
}

int CorrectQuery::correct( string &seq, bool drxn )
{
    CorrectAlign* full = NULL,* part = NULL;
    
    bool bad = false;
    for ( CorrectBranch &cb : branches_ ) cb.correct( full, part, bad, drxn && initial_ );
    
    int partLen = part ? part->perf : 0, imperf = 0;
    for ( CorrectBranch &cb : branches_ ) cb.contend( full, part, partLen, imperf );
    
    if ( full )
    {
        full->write( seq, full->lens[0], drxn );
        assert( seq.size() == seqLen_ );
        return full->lens[0];
    }
    else if ( drxn && part )
    {
        part->write( seq, partLen, drxn );
        if ( part->miss > 2 && ( part->score - part->perf < 0 ) ) trim_ = drxn;
        assert( seq.size() == seqLen_ );
        return partLen;
    }
    else if ( bad ) trim_ = drxn;
    
    return 0;
}

bool CorrectQuery::proceed( IndexReader* ir )
{
    if ( branches_.empty() || qLen_ >= min( coords_[1] - coords_[0], int( seqLen_ * .9 ) ) ) return false;
    int maxCount = 0;
    for ( CorrectBranch &cb : branches_ ) if ( cb.count > maxCount ) maxCount = cb.count;
    if ( maxCount < 10 ) return false;
    for ( int i = 0; i < branches_.size(); i++ )
    {
        if ( maxCount >= 30 && branches_[i].count < 2 ) branches_.erase( branches_.begin() + i-- );
        else branches_[i].query( ir, (int)q_.size() - qLen_ );
    }
    
    return true;
}

ReadId CorrectQuery::query( IndexReader* ir, uint8_t i, int it, CharId rank, CharId count )
{
    if ( count ) qLen_ = it;
    else return 0;
    
    CharCount ranks;
    CharCount counts;
    
    ir->countRange( i, rank, count, ranks, counts );
    ReadId ends = it < q_.size() && q_[it] < 4 ? query( ir, q_[it], it+1, ranks[ q_[it] ], counts[ q_[it] ] ) : 0;
    
    if ( !branches_.empty() || qLen_ == q_.size() || ends > 2 ) return ends;
    
    for ( int j = 0; j < 4; j++ ) if ( j != q_[it] && counts[j] ) branches_.push_back( CorrectBranch( j, it, ranks[j], counts[j] ) );
    
    return ends + counts.endCounts;
}

bool CorrectQuery::setQuery( string &seq, bool drxn )
{
    q_.clear();
    if ( coords_[1] - coords_[0] < 12 ) return false;
    if ( drxn ) for ( int i = coords_[0]; i < seq.size(); i++ ) q_.push_back( charToIntComp[ seq[i] ] );
    else for ( int i = coords_[1]; --i >= 0; ) q_.push_back( charToInt[ seq[i] ] );
    for ( int i = 0; i < 12; i++ ) if ( q_[i] > 3 ) return false;
    return true;
}
