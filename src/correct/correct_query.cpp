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

CorrectExt::CorrectExt( string &seq, int len, bool drxn )
{
    maxLen = seq.size();
    ext = len >= seq.size() ? seq : ( drxn ? seq.substr( 0, len ) : seq.substr( seq.size() - len ) );
    lens.push_back( min( seq.size(), ext.size() ) );
}

bool CorrectExt::add( string &seq, bool drxn, bool doAdd )
{
    for ( int i = 0; i < min( seq.size(), ext.size() ); i++ )
    {
        if ( seq[i] == 'N' ) continue;
        if ( drxn ? seq[i] != ext[i] : seq.end()[-i-1] != ext.end()[-i-1] ) return false;
    }
    if ( doAdd ) lens.push_back( min( seq.size(), ext.size() ) );
    return true;
}

vector<string> CorrectExt::get( vector<Overlap> ols, int len, bool drxn )
{
    vector<CorrectExt> exts;
    for ( Overlap &ol : ols )
    {
        bool added = false, strong = false;
        vector<int> weak;
        for ( int i = 0; i < exts.size(); i++ )
        {
            if ( !exts[i].add( ol.seq, drxn ) ) continue;
            if ( exts[i].lens.size() == 2 ) weak.push_back( i );
            if ( exts[i].lens.size() > 9 && exts[i].maxLen - ol.extLen > 29 ) strong = true;
            added = true;
        }
        if ( !added ) exts.push_back( CorrectExt( ol.seq, len, drxn ) );
        if ( !strong || weak.empty() ) continue;
        for ( auto it = weak.rbegin(); it != weak.rend(); it++ ) exts.erase( exts.begin() + *it );
    }
    bool strong = false;
    for ( CorrectExt &ext : exts ) if ( ext.lens.size() > 9 && ext.maxLen > 29 ) strong = true;
    
    vector<string> seqs;
    for ( CorrectExt &ext : exts ) if ( !strong || ext.lens.size() > 1 ) seqs.push_back( ext.ext );
    return seqs;
}

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
    assert( false );
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

void CorrectBranch::collect( vector<uint8_t> &ols, vector<CharId> &ranks, vector<CharId> &counts )
{
    int j = 0;
    for ( int i = 0; i < endOverlaps.size(); )
    {
        while ( j < ols.size() && endOverlaps[i] < ols[j] ) j++;
        int k = j < ols.size() ? i + 1 : endOverlaps.size();
        while ( k < endOverlaps.size() && ols[j] <= endOverlaps[k] ) k++;
        ols.insert( ols.begin(), endOverlaps.begin()+i, endOverlaps.begin()+k );
        ranks.insert( ranks.begin(), endRanks.begin()+i, endRanks.begin()+k );
        counts.insert( counts.begin(), endCounts.begin()+i, endCounts.begin()+k );
        i = k;
    }
    for ( CorrectBranch &cb : branches ) cb.collect( ols, ranks, counts );
}

bool CorrectBranch::beats( CorrectBranch &rhs, int len )
{
    return ( rhs.bad && ( good == len || max( good, bad ) - rhs.bad > 12 ) );
}

int CorrectBranch::getNovel()
{
    uint8_t last = -1;
    int novel = 0, len = min( (int)q.size(), max( good, bad ) );
    for ( int i = 0; i < len; i++ )
    {
        if ( q[i] != last ) novel++;
        last = q[i];
    }
    return novel;
}

bool CorrectBranch::steal( vector<uint8_t> &ols, vector<CharId> &ranks, vector<CharId> &counts, int len )
{
    for ( int i = 0; i < endOverlaps.size(); i++ )
    {
        if ( endOverlaps[i] > len ) continue;
        ols.insert( ols.begin(), endOverlaps.begin()+i, endOverlaps.end() );
        ranks.insert( ranks.begin(), endRanks.begin()+i, endRanks.end() );
        counts.insert( counts.begin(), endCounts.begin()+i, endCounts.end() );
        endOverlaps.erase( endOverlaps.begin()+i, endOverlaps.end() );
        endRanks.erase( endRanks.begin()+i, endRanks.end() );
        endCounts.erase( endCounts.begin()+i, endCounts.end() );
    }
    int cut = len - base;
    assert( cut <= q.size() );
    q.erase( q.begin(), q.begin() + len - base );
    base += cut;
    return q.empty();
}

string CorrectBranch::yield( bool drxn )
{
    string s;
    for ( int i = 0; i < q.size(); i++ ) s += drxn ? intToChar[ q.end()[-i-1] ] : intToCharComp[ q[i] ];
    return s;
}

vector<string> CorrectBranch::yield( vector<CorrectBranch> &cbs, bool drxn )
{
    vector<string> exts;
    for ( CorrectBranch &cb : cbs ) exts.push_back( cb.yield( drxn ) );
    return exts;
}

CorrectQuery::CorrectQuery( IndexReader* ir, string &seq, int qLen, bool retract )
: ir_( ir ), base_( 0 ), qLen_( 0 ), seqLen_( seq.size() ), collect_( true ), trim_( false )
{
    if ( retract ) while ( qLen > 1 && seq[qLen-1] == seq[qLen-2] ) qLen--;
    for ( int i = qLen; --i >= 0; ) q_.push_back( charToInt[ seq[i] ] );
    
    base_ = seqLen_ - q_.size();
    collect_ = true;
    if ( q_.size() < 16 ) return;
    CharId rank, count;
//    ir->setBaseAll( q_[0], q_[1], rank, count );
//    query( q_[1], 2, rank, count );
    int it = ir->setBaseAll( q_, rank, count );
    if ( count ) query( q_[it-1], it, rank, count );
}

CorrectQuery::CorrectQuery( IndexReader* ir, string &seq, bool drxn )
: ir_( ir ), base_( 0 ), qLen_( 0 ), seqLen_( seq.size() ), collect_( false ), trim_( false )
{
    for ( int i = 0; i < seq.size(); i++ )
    {
        if ( drxn ) q_.push_back( charToInt[ seq.end()[-i-1] ] );
        else q_.push_back( charToIntComp[ seq[i] ] );
        if ( q_.back() > 3 && i - base_ < 32 )
        {
            base_ = i + 1;
            q_.clear();
            collect_ = true;
        }
    }
    if ( !ir || q_.size() < 32 ) return;
    CharId rank, count;
//    ir->setBaseAll( q_[0], q_[1], rank, count );
//    query( q_[1], 2, rank, count );
    int it = ir->setBaseAll( q_, rank, count );
    if ( count ) query( q_[it-1], it, rank, count );
}

CorrectQuery::CorrectQuery( IndexReader* ir, vector<uint8_t> &q, vector<CorrectBranch> &branches, int seqLen )
: ir_( ir ), base_( 0 ), qLen_( q.size() ), seqLen_( seqLen ), q_( q ), branches_( branches ), collect_( false ), trim_( false )
{
    for ( CorrectBranch &cb : branches_ ) query( cb, cb.q.back(), cb.base+1 );
}

int CorrectQuery::correct( QueryBinaries* qb, string &seq, bool &trimmed )
{
    if ( qLen_ < seqLen_ * .9 && !branches_.empty() && base_ + qLen_ < seq.size() )
    {
        vector<string> base( 1, seq.substr( base_ + qLen_ ) );
        vector<string> exts = CorrectBranch::yield( branches_, 0 );
        vector<CorrectAlign> cas = CorrectAlign::get( base, exts, 1 );
        qLen_ += correct( cas, seq, 0, true, 1 );
        if ( trim_ ) trimmed = true;
    }
    
    if ( base_ )
    {
        for ( CorrectBranch &cb : branches_ ) cb.collect( endOverlaps_, endRanks_, endCounts_ );
        vector<string> base( 1, seq.substr( 0, base_ ) );
        vector<string> exts = CorrectExt::get( qb->getOverlaps( endOverlaps_, endRanks_, endCounts_, 80, 200, 0 ), base_, 0 );
        vector<CorrectAlign> cas = CorrectAlign::get( base, exts, 0 );
        int extLen = correct( cas, seq, 0, false, 0 );
        if ( !extLen && seq.find_first_not_of( "N" ) == base_ ) extLen = base_;
        base_ -= extLen;
        qLen_ += extLen;
        if ( !extLen ) qLen_ = 0;
    }
    
    return qLen_;
}

int CorrectQuery::correct( string &seq, vector<string> &seqs, bool &trimmed )
{
    assert( !seqs.empty() && !base_ );
    if ( branches_.empty() || qLen_ >= seqLen_ * .9 ) return qLen_;
    vector<string> base;
    for ( string &s : seqs ) base.push_back( s.substr( qLen_ + base_ ) );
    vector<string> exts = CorrectBranch::yield( branches_, 0 );
    vector<CorrectAlign> cas = CorrectAlign::get( base, exts, 1 );
    qLen_ += correct( cas, seq, 0, true, 1 );
    if ( trim_ ) trimmed = true;
    
    string q = seq.substr( 0, qLen_ );
    for ( int i = 0; i < seqs.size(); i++ ) if ( seqs[i].find( q ) != 0 ) seqs.erase( seqs.begin() + i-- );
    assert( !seqs.empty() );
    
    return qLen_;
}

int CorrectQuery::correct( vector<CorrectAlign> &cas, string &seq, int maxMiss, bool partial, bool drxn )
{
    CorrectAlign* full = NULL,* part = NULL,* bad = NULL;
    for ( CorrectAlign &ca : cas )
    {
        if ( ca.isBad() ) bad = &ca;
        if ( ( !full || ca > *full ) && ca.isFull( maxMiss ) ) full = &ca;
        if ( ( !part || ( ca.perf == part->perf ? ca > *part : ca.perf > part->perf ) ) && ca.isPart() ) part = &ca;
    }
    
    int partLen = part ? part->perf : 0, imperf = 0;
    for ( CorrectAlign &ca : cas )
    {
        if ( full && full != &ca && ca.contendFull( *full ) ) full = NULL;
        if ( part && part != &ca && ca.contendPart( *part, partLen ) ) part = NULL;
        if ( !ca.perf ) imperf++;
    }
    
    if ( full && ( full->miss || full->poly ) && cas.size() > 2 && full->hit < ( 5 + cas.size() / 4 ) ) full = NULL;
    if ( part && part->perf <= imperf ) part = NULL;
    
    if ( full )
    {
        full->write( seq, full->lens[0], drxn );
        assert( seq.size() == seqLen_ );
        return full->lens[0];
    }
    else if ( partial && part )
    {
        part->write( seq, partLen, drxn );
        if ( part->miss > 2 && ( part->score - part->perf < 0 ) ) trim_ = true;
        assert( seq.size() == seqLen_ );
        return partLen;
    }
    else if ( bad ) trim_ = true;
    
    return 0;
}

void CorrectQuery::query( uint8_t i, int it, CharId rank, CharId count )
{
    CharCount ranks;
    CharCount counts;
    
    ir_->countRange( i, rank, count, ranks, counts );
    if ( it < q_.size() && q_[it] < 4 && counts[ q_[it] ] )
    {
        qLen_ = it + 1;
        if ( qLen_ < q_.size() || collect_ ) query( q_[it], it+1, ranks[ q_[it] ], counts[ q_[it] ] );
    }
    
    for ( int j = 0; it >= qLen_ && j < 4; j++ )
    {
        if ( !counts[j] || ( it < q_.size() && j == q_[it] ) ) continue;
        if ( it < q_.size() && q_[it] < 4 && counts[j] <= counts[ q_[it] ] ) continue;
        branches_.push_back( CorrectBranch( j, it, ranks[j], counts[j] ) );
        query( branches_.back(), j, it+1 );
    }
    
    if ( !collect_ || it < 32 || it == seqLen_ || !counts.endCounts ) return;
    
    endOverlaps_.push_back( it );
    endRanks_.push_back( ranks.endCounts );
    endCounts_.push_back( counts.endCounts );
}

void CorrectQuery::query( CorrectBranch &cb, uint8_t i, int it )
{
    CharCount ranks;
    CharCount counts;
    
    ir_->countRange( i, cb.rank, cb.count, ranks, counts );
    bool free = cb.bad || it >= q_.size() || q_[it] > 3 || !counts[ q_[it] ];
    int jBest = free ? counts.getMaxBranch() : q_[it];
    int jCount = free ? counts.getBranchCount() : 1;
    if ( !cb.bad && it < q_.size() && jCount && q_[it] < 4 && jBest != q_[it] ) cb.bad = it;
    if ( !cb.bad && !cb.good && ( it >= q_.size() || !jCount ) ) cb.good = min( it, (int)q_.size() );
    if ( jBest < 4 && jCount == 1 )
    {
        cb.q.push_back( jBest );
        cb.rank = ranks[jBest];
        cb.count = counts[jBest];
        query( cb, jBest, it+1 );
    }
    else for ( int j = 0; j < 4; j++ ) if ( counts[j] ) cb.branches.push_back( CorrectBranch( j, it, ranks[j], counts[j] ) );
    
    if ( !collect_ || it < 32 || it == seqLen_ || !counts.endCounts ) return;
    
    cb.endOverlaps.push_back( it );
    cb.endRanks.push_back( ranks.endCounts );
    cb.endCounts.push_back( counts.endCounts );
}

int CorrectQuery::trim( QueryBinaries* qb, string &seq, bool &trimmed )
{
    for ( CorrectBranch &cb : branches_ ) cb.collect( endOverlaps_, endRanks_, endCounts_ );
    vector<string> base( 1, seq.substr( q_.size() ) );
    vector<string> exts = CorrectExt::get( qb->getOverlaps( endOverlaps_, endRanks_, endCounts_, 80, 200, 1 ), base_, 1 );
    vector<CorrectAlign> cas = CorrectAlign::get( base, exts, 1 );
    int extLen = correct( cas, seq, 1, true, 1 );
    if ( trim_ ) trimmed = true;
    return q_.size() + extLen;
}

int CorrectQuery::trim( QueryBinaries* qb, string &seq, vector<string> &seqs, bool &trimmed )
{
    assert( !seqs.empty() );
    for ( CorrectBranch &cb : branches_ ) cb.collect( endOverlaps_, endRanks_, endCounts_ );
    vector<string> base;
    for ( string &s : seqs ) base.push_back( s.substr( q_.size() ) );
    vector<string> exts = CorrectExt::get( qb->getOverlaps( endOverlaps_, endRanks_, endCounts_, 80, 200, 1 ), base_, 1 );
    vector<CorrectAlign> cas = CorrectAlign::get( base, exts, 1 );
    int extLen = correct( cas, seq, 1, true, 1 );
    if ( trim_ ) trimmed = true;
    return q_.size() + extLen;
}
