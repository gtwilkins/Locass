/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "correct_read.h"
#include "parameters.h"
#include "shared_functions.h"
#include "local_alignment.h"
#include <cassert>
#include <iostream>
#include <algorithm>

extern struct Parameters params;

Chunk::Chunk( string &s, int coords[2], bool inCaps[2] )
: seq( s )
{
    ends[0] = conf[0] = coords[0];
    ends[1] = conf[1] = coords[1];
    caps[0] = inCaps[0];
    caps[1] = inCaps[1];
    drxn = 2;
}

Chunk::Chunk( string &s, int coords[2], bool inCaps[2], int confLen, bool bridged, bool d )
: Chunk( s, coords, inCaps )
{
    ends[d] = d ? ends[0] + s.size() : ends[1] - s.size();
    conf[d] = bridged ? ends[d] : ( d ? conf[1] + confLen : conf[0] - confLen );
    drxn = d;
}

bool Chunk::match( ReadId id, string &q, int cutoff, bool doCutoff )
{
    int coords[2]{0};
    if ( !mapSeq( q, seq, coords, 32 ) ) return false;
    coords[0] += ends[0];
    coords[1] += ends[0];
    if ( doCutoff && drxn == 1 && coords[1] <= cutoff ) return true;
    if ( doCutoff && drxn == 0 && cutoff <= coords[0] ) return true;
    hits[id] = make_pair( coords[0], coords[1] );
    return true;
}

bool Chunk::matchPartial( ReadId id, string &q )
{
    if ( drxn == 2 || conf[drxn] == ends[drxn] ) return false;
    int coords[2]{0};
    if ( !mapSeqEnd( q, seq, 32, coords, drxn ) ) return false;
    coords[0] += ends[0];
    coords[1] += ends[0];
    if ( drxn ? coords[1] < conf[1] : conf[0] < coords[0] ) return false;
    hits[id] = make_pair( coords[0], coords[1] );
    return true;
}

void Chunk::reverse( vector<Chunk> &cs, int len )
{
    for ( Chunk &c : cs )
    {
        revComp( c.seq );
        int coords[2] = { len - c.ends[1], len - c.ends[0] };
        bool tmpCaps[2] = { c.caps[1], c.caps[0] };
        c.ends[0] = coords[0];
        c.ends[1] = coords[1];
        c.caps[0] = tmpCaps[0];
        c.caps[1] = tmpCaps[1];
    }
}

void Chunk::setLimits( int* limits )
{
    limits[0] = min( limits[0], ends[0] );
    limits[1] = max( limits[1], ends[1] );
}

AlignStruct::AlignStruct( CorrectionExt &ext, string &t, string &p, bool drxn )
: seq( ext.seq ), ids( ext.ids.begin(), ext.ids.end() ), drxn( drxn )
{
    lens[0] = seq.size();
    lens[1] = ext.lens[0];
    lens[2] = ext.lens[1];
    hitLen = 0;
    align( t, p );
}

bool AlignStruct::addLen( int len )
{
    len = min( len, (int)seq.size() );
    for ( int i = 0; i < 3; i++ )
    {
        if ( len <= lens[i] ) continue;
        for ( int j = 3; --j > i; ) lens[j] = lens[j-1];
        lens[i] = len;
        return !i;
    }
    return false;
}

void AlignStruct::align( string &t, string &p )
{
    ends[0] = ends[1] = -1;
    ends[2] = misses = polys = novelty = 0;
    a[0].clear();
    a[1].clear();
    errors.clear();
    
    bool unique = false;
    string q = t[0] + ( drxn ? seq : string( seq.rbegin(), seq.rend() ) );
    for ( int i = 0; i < min( q.size(), t.size() ) && !unique; i++ ) if ( t[i] != q[i] ) unique = true;
    if ( !unique )
    {
        a[0] = t + string( max( 0, int( q.size() - t.size() ) ), '-' );
        a[1] = q + string( max( 0, int( t.size() - q.size() ) ), '-' );
        ends[2] = min( q.size(), t.size() );
        ends[0] = ends[1] = ends[2]-1;
        return;
    }
    
    LocalAlignment lalign( t, q, false, true );
    lalign.realign( a[0], a[1], true );
    
    for ( int d = 0; d < 2; d++ )
    {
        int i = 0;
        while ( a[d][i] == '-' ) i++;
        for ( int j = 0; d && j+1 < i && a[0][j] == a[0][0]; j++ ) if ( a[0][j+1] == 'N' ) a[0][j+1] = a[0][0];
        for ( int j = i; i && j < a[d].size() && a[d][j-1] == '-' && a[d][j] == a[!d][0] && a[d][j] == a[!d][j-i]; j++ )
        {
            a[d][j-i] = a[d][j];
            a[d][j] = '-';
            if ( a[!d][j] != '-' ) continue;
            a[d].erase( j, 1 );
            a[!d].erase( j, 1 );
            j--; i--;
        }
    }
    
    assert( a[0].size() == a[1].size() );
    
    int limit = min( getAlignLimit( 0 ), getAlignLimit( 1 ) );
    assert( limit > 0 );
    
    for ( int& i = ends[2]; i < limit; i++ )
    {
        assert( a[0][i] != '-' || a[1][i] != '-' );
        for ( int d = 0; d < 2 && i < limit; d++ )
        {
            if ( a[d][i] != '-' ) continue;
            assert( i );
            bool poly = true;
            for ( ; i < limit && a[d][i] == '-'; i++ )
            {
                ends[!d]++;
                if ( i && a[!d][i] == a[!d][i-1] && ( a[d][i-1] == '-' || a[!d][i] == a[d][i-1] )  ) continue;
                errors.push_back( make_pair( ends[0], ends[1] ) );
                misses += ( p[ ends[0]-1 ] + p[ ends[0] ] - 70 ) / 2;
                poly = false;
            }
            if ( poly ) polys++;
            d = -1;
        }
        if ( i == limit ) break;
        assert( a[0][i] != '-' && a[1][i] != '-' );
        ends[0]++;
        ends[1]++;
        if ( a[0][i] != a[1][i] )
        {
            errors.push_back( make_pair( ends[0], ends[1] ) );
            misses += p[ ends[0] ]-35;
        }
    }
    
    for ( int i = 0; i < ends[1]; i++ ) if ( q[i] != q[i+1] ) novelty++;
}

int AlignStruct::charCount( char &c, int &i )
{
    assert( i <= seq.size() );
    c = drxn ? seq[i] : seq.end()[-i-1];
    int copies = 1;
    for ( ; ++i < seq.size() && ( drxn ? seq[i] : seq.end()[-i-1] ) == c; copies++ );
    return copies;
}

bool AlignStruct::charCount( char &c, int &i, int &copies )
{
    for ( ; i < a[1].size() && c == '-'; i++ ) c = a[1][i];
    bool diff = i && ( a[1][i-1] == c || a[0][i-1] == c ), gap = false;
    for ( int j = i-2; !diff && i >= 0 && a[1][j+1] == '-'; j++ ) if ( a[0][j] == c || a[1][j] == c ) diff = false;
    assert( ( i && ( a[1][i-1] == c || a[0][i-1] == c ) ) || !diff );
    
    for ( ; ( !gap || i < ends[2] ) && ( a[1][i] == '-' || a[1][i] == c ); i++ )
    {
        if ( gap ? ( a[0][i] == c || a[1][i] == c ) : ( i < ends[2] && a[0][i] != a[1][i] ) ) diff = true;
        if ( a[1][i] == c ) copies++;
        if ( i+1 == a[1].size() || a[1][i+1] == '-' ) gap = true;
        
    }
    if ( i < ends[2] && a[1][i] == c ) assert( false );
    if ( i < ends[2] && a[0][i] == c ) diff = true;
    return diff;
}

bool AlignStruct::consolidate( AlignStruct &as, int polyLen )
{
    char c[2];
    int i[2]{0}, polyCount = 0, diff = 0, polyNovel = 0;
    while ( i[0] < seq.size() && i[1] < as.seq.size() )
    {
        int copies[2]{ charCount( c[0], i[0] ), as.charCount( c[1], i[1] ) };
        for ( int j = 0; j < 2; j++ ) if ( i[j] - copies[j] == 0 ) copies[j] += polyLen - 1;
        if ( c[0] != c[1] ) return false;
        if ( polyCount ) polyNovel++;
        if ( copies[0] == copies[1] || i[0] == seq.size() || i[1] == as.seq.size() ) continue;
        int minCopy = min( copies[0], copies[1] );
        if ( polyCount++ || minCopy < 4 || i[1] < as.lens[1 + ( i[1] - copies[1] < 8 ) ] ) return false;
        diff = copies[1] - copies[0];
    }
    
    if ( polyNovel < 5 ) return false;
    addLen( min( lens[0], as.lens[0] + diff ) );
    for ( ReadId id : as.ids ) if ( find( ids.begin(), ids.end(), id ) == ids.end() ) ids.push_back( id );
    
    return true;
}

int AlignStruct::getAlignLimit( bool i )
{
    for ( int j = a[i].size(); j > 0; j-- ) if ( a[i][j-1] != '-' ) return j;
    return 0;
}

AlignStruct* AlignStruct::getExtension( AlignStruct*& best, vector<AlignStruct> &ass, int gapLen, int polyLen, bool initial )
{
    best = NULL;
    for ( AlignStruct &as : ass ) as.setExtLens( gapLen, as.seq.size() );
    if ( !initial )
    {
        for ( int i = 0; i < ass.size()-1; i++ )
        {
            for ( int j = i+1; j < ass.size(); j++ )
            {
                bool better[2] = { ass[i].isBetter( ass[j] ), ass[j].isBetter( ass[i] ) };
                if ( better[0] && !better[1] && ( ass[i].valid || !ass[j].bridged || !ass[j].valid ) ) ass.erase( ass.begin() + j-- );
                if ( !better[0] && better[1] && ( ass[j].valid || !ass[i].bridged || !ass[i].valid ) ) { ass.erase( ass.begin() + i-- ); break; }
            }
        }
    }
    for ( AlignStruct &as : ass )
    {
        if ( best )
        {
            if ( best->bridged && best->valid && !as.valid ) continue;
            else if ( as.bridged && as.valid && !best->valid );
            else if ( best->hitLen > as.hitLen ) continue;
            else if ( as.hitLen > best->hitLen );
            else if ( best->hitIds.size() > as.hitIds.size() ) continue;
            else if ( as.hitIds.size() > best->hitIds.size() );
            else if ( best->ids.size() >= as.ids.size() ) continue;
        }
        if ( as.lens[1] ) best = &as;
    }
    if ( !best ) return NULL;
    
    int iBest = -1;
    for ( int i = 0; i < ass.size(); i++ ) if ( &ass[i] == best ) iBest = i;
    assert( iBest >= 0 );
    for ( int i = 0; i < ass.size(); i++ )
    {
        if ( i == iBest || !best->consolidate( ass[i], polyLen ) ) continue;
        ass.erase( ass.begin() + i-- );
        if ( i < iBest ) best = &ass[--iBest];
    }
    if ( !best->lens[2] ) return NULL;
    
    
    int congruence = best->getExtLen();
    if ( !congruence ) return NULL;
    for ( AlignStruct &as : ass ) if ( &as != best && ( initial || !best->bridged || as.valid ) ) best->setCongruence( as, congruence );
    if ( best ) best->setExtLens( gapLen, congruence );
    if ( best && !best->bridged && best->aligns[0] < 1 && best->aligns[1] < 1 ) best = NULL;
    return best;
}

int AlignStruct::getExtLen()
{
    int anchor = errors.size() < 2 ? seq.size() : 0;
    for ( int i = errors.size(); i >= 0 && !anchor; i-- )
    {
        int cuts[2] = { i > 1 ? errors[i-1].second : -1, i > 2 ? errors[i-2].second : -1 };
        int last = i == errors.size() ? seq.size() : errors[i].second;
        int count = 1;
        for ( int j = last-1; --j > cuts[1] && !anchor; )
        {
            if ( drxn ? seq[j] != seq[j+1] : seq.end()[-j-1] != seq.end()[-j] ) count++;
            if ( count > ( j > cuts[0] ? 7 : 10 ) ) anchor = last;
        }
    }
    
    return min( anchor, lens[2] );
}

int AlignStruct::getNovelty( vector<AlignStruct>* ass )
{
    if ( !valid || ( !bridged && !ass ) ) return novelty;
    
    int ol = 0, score = novelty;
    if ( ass ) for ( AlignStruct &as : *ass ) ol = max( ol, mapSeqOverlap( ( drxn ? seq : as.seq ), ( drxn ? as.seq : seq ), 20 ) );
    if ( !bridged && !ol ) return novelty;
    
    char c;
    for ( int i = ( bridged ? aligns[1] : seq.size() - ol ); i < seq.size() && charCount( c, i ); score++ );
    
    return score;
}

string AlignStruct::getPhred( string &p, int seqLen, bool drxn )
{
    string phred;
    int phredLens[2] = { -1, -1 };
    for ( int i = 0; i < min( a[0].size(), a[1].size() ) && phred.size() < seqLen && phredLens[0] < (int)p.size(); i++ )
    {
        bool add[2] = { a[0][i] != '-', a[1][i] != '-' };
        if ( add[0] ) phredLens[0]++;
        if ( add[1] && ++phredLens[1] )
        {
            if ( phredLens[0] < 0 || a[0][i] != a[1][i] ) add[0] = false;
            phred += ( add[0] ? ( drxn ? p[ phredLens[0]-1 ] : p.end()[-phredLens[0] ] ) : '#' );
        }
    }
    
    assert( phred.size() <= seqLen );
    return drxn ? phred + string( seqLen - phred.size(), '#' ) : string( seqLen - phred.size(), '#' ) + phred;
}

int AlignStruct::getSharedIdCount( AlignStruct &as )
{
    int i = 0;
    while ( i < min( ids.size(), as.ids.size() ) && ids.end()[-i-1] == as.ids.end()[-i-1] ) i++;
    return i;
}

string AlignStruct::getSeq( int seqLen, bool drxn )
{
    return drxn ? seq.substr( 0, seqLen ) : seq.substr( seq.size() - seqLen );
}

int AlignStruct::getValidGap( int gapLen )
{
    int len = -gapLen -1, unique = 0, marks[2]{ (int)a[0].size(), -1 }, uniques[2]{0}, runs[2]{ -1, -1 }, gap = 2;
    char c = a[0][0];
    for ( int i = 0; i < a[0].size(); i++ )
    {
        if ( a[0][i] != '-' && ( !len++ || ( len > 0 && a[0][i] != a[0][i-1] ) ) ) unique++;
        if ( gap == 2 ) gap = a[0][i] == '-' ? 0 : ( a[1][i] == '-' ? 1 : 2 );
        if ( a[0][i] == ( gap == 0 ? '-' : c ) && a[1][i] == ( gap == 1 ? '-' : c ) )
        {
            for ( int d = 0; d < 2; d++ ) runs[d] += a[d][i] == c;
            if ( gap < 2 && i < ends[2] && ( a[gap][i+1] != '-' || a[!gap][i+1] != c ) && len > runs[gap] )
            {
                marks[1] = len - runs[gap];
                uniques[1] = unique;
                if ( marks[0] > len ) uniques[0] = unique;
                marks[0] = min( marks[0], len - abs( runs[0] - runs[1] ) );
            }
            continue;
        }
        gap = 2;
        c = a[0][i] == a[1][i] || a[1][i] == '-' ? a[0][i] : ( a[0][i] == '-' ? a[1][i] : '-' );
        for ( int d = 0; d < 2; d++ ) runs[d] = a[d][i] == c;
        if ( a[0][i] == a[1][i] || len < 1 || i >= ends[2] ) continue;
        marks[1] = len;
        uniques[1] = unique;
        if ( marks[0] < len ) continue;
        marks[0] = len - ( a[0][i] != '-' );
        uniques[0] = unique;
    }
    marks[0] = len - min( marks[0], len );
    uniques[0] = unique - uniques[0];
    assert( marks[1] > 0 );
    if ( uniques[0] < 3 && uniques[1] < 3 ) return len;
    if ( uniques[0] < 3 ) return -marks[0];
    if ( uniques[1] < 3 ) return marks[1];
    
    return 0;
}

bool AlignStruct::isBetter( AlignStruct &as, bool initial )
{
    if ( !initial )
    {
        for ( const ReadId &id : hitIds ) if ( as.hitIds.find( id ) == as.hitIds.end() ) return true;

        for ( int i = 0; i < hitLen && i < as.seq.size(); i++ )
        {
            if ( drxn ? seq[i] != as.seq[i] : seq.end()[-i-1] != as.seq.end()[-i-1] ) return true;
        }
        
        return false;
    }
    
    if ( errors.empty() && polys < 2 && novelty > 8 && lens[2] > 10 && !as.errors.empty() && as.errors.size() + as.polys > 1 )
    {
        if ( aligns[1] != lens[2] ) setExtLens( a[0].size(), lens[2] );
        if ( as.errors[0].first < lens[2] )
        {
            int shared = getSharedIdCount( as );
            if ( ( ids.size() - shared ) * 3 > as.ids.size() - shared ) return true;
        }
    }
    
    int errs[2]{0};
    for ( pair<int,int> &e : errors ) if ( e.first <= as.ends[0] || e.second <= as.ends[1] ) errs[0]++;
    for ( pair<int,int> &e : as.errors ) if ( e.first <= ends[0] || e.second <= ends[1]-4 ) errs[1]++;
    return ( errs[1] - errs[0] ) > ( max( 0, 2 - novelty / 10 ) + ( errs[0] > ( novelty / 8 ) ) );
}

bool AlignStruct::isBetterErase( AlignStruct &as )
{
    if ( !valid && as.bridged && as.valid ) return false;
    if ( !as.errors.empty() || !as.lens[2] || errors.empty() ) return true;
    return true;
}

void AlignStruct::setCongruence( AlignStruct &as, int &congruence )
{
    int i = 0, j = 0, limit = getAlignLimit( 0 );
    int aligns[2] = { -1, -1 }, error = seq.size(), poly = seq.size();
    int cutoff = as.lens[2] < lens[2] - 30 ? as.lens[2] : ( as.lens[1] < lens[2] - 20 ? as.lens[1] : as.lens[0] );
    int cutin = 10 + ( lens[0] - 10 ) / 2;
    
    bool anyDiff = hitLen < lens[2] && as.lens[2];
    while ( aligns[0] < error && aligns[1] < as.lens[0] )
    {
        char c[2] = { a[1][i], as.a[1][j] };
        int copies[2]{0};
        bool diffs[2] = { charCount( c[0], i, copies[0] ), as.charCount( c[1], j, copies[1] ) };
        if ( diffs[0] && !diffs[1] ) anyDiff = true;
        if ( copies[1] >= as.seq.size() ) return;
        
        if ( c[0] != c[1] ) error = min( aligns[0], error );
        else if ( copies[0] != copies[1] )
        {
            int err[2] = { aligns[0] + min( copies[0], copies[1] ), aligns[1] + copies[1] };
            if ( i >= limit ) error = err[0];
            else if ( diffs[0] && err[1] < cutoff ) anyDiff = true;
            poly = min( poly, err[0] );
        }

        aligns[0] += copies[0];
        aligns[1] += copies[1];
    }
    
    if ( ( anyDiff && poly < cutoff ) || poly > cutin ) congruence = min( congruence, poly );
    if ( poly < ends[1] && error > ends[1] ) poly = error;
    if ( error < cutoff || min( error, poly ) > cutin ) congruence = min( congruence, min( poly, error ) );
}

bool AlignStruct::setExtLens( int gapLen, int extLen )
{
    assert( a[0].size() == a[1].size() );
    bridged = false;
    valid = true;
    aligns[0] = aligns[1] = -1;
    
    bool bridge = gapLen < ends[0];
    int i = 0, limit = getAlignLimit( 1 ), last[2] = { -1, -1 };
    
    while ( i < ( bridged ? ends[2] : limit  ) && ( !bridged || valid ) )
    {
        assert( i < a[0].size() );
        int runs[2] = { a[0][i] != '-', a[1][i] != '-' };
        bool same = a[0][i] == a[1][i];
        for ( int j = 0; same && j < 2; j++ )
        {
            bool gap = i+1 == a[!j].size() || a[!j][i+1] == '-';
            for ( int &k  = runs[j]; i+k < a[j].size() && a[j][i] == a[j][i+k] && a[!j][i+k] == ( gap ? '-' : a[j][i] ); k++ ) if ( i+k+1 == a[j].size() || a[!j][i+k+1] == '-' ) gap = true;
        }
        i += max( runs[0], runs[1] );
        if ( ( aligns[0] == gapLen && !same && runs[0] ) || ( bridged && i <= ends[2] && runs[0] != runs[1] ) ) valid = !bridge;
        if ( bridged ) continue;
        if ( aligns[1] + runs[1] >= extLen )
        {
            int run = runs[0] ? min( extLen - aligns[1], gapLen - aligns[0] ) : extLen - aligns[1];
            assert( run >= 0 );
            aligns[0] += min( run, runs[0] );
            aligns[1] += min( run, runs[1] );
            if ( aligns[0] < gapLen && !same )
            {
                aligns[0] = last[0];
                aligns[1] = last[1];
            }
            valid = true;
            return valid;
        }
        if ( bridge && aligns[0] + runs[0] > gapLen )
        {
            bridged = true;
            int excess = aligns[0] + runs[0] - gapLen;
//            if ( aligns[1] + runs[1] - excess < 0 ) assert( false );
            if ( excess > runs[1] || aligns[1] + runs[1] - excess < 0 ) valid = false;
            runs[0] -= excess;
            runs[1] = max( 0, runs[1] - excess );
            assert( runs[1] >= 0 );
        }
        if ( same )
        {
            last[0] = aligns[0] + min( runs[0], runs[1] );
            last[1] = aligns[1] + min( runs[0], runs[1] );
        }
        aligns[0] += runs[0];
        aligns[1] += runs[1];
    }
    
    return valid;
}

void AlignStruct::setHits( vector<Chunk> &t, vector<bool> &tHits, int coord )
{
    if ( t.empty() ) return;
    for ( ReadId &id : ids )
    {
        if ( !params.isReadPe( id ) ) continue;
        ReadId pairId = params.getPairId( id );
        for ( Chunk &c : t )
        {
            auto it = c.hits.find( pairId );
            if ( it == c.hits.end() ) continue;
            for ( int i = max( 0, it->second.first - coord ); i < min( (int)tHits.size(), it->second.second - coord ); i++ ) tHits[i] = true;
            hitIds.insert( id );
        }
    }
}

bool AlignStruct::setReads( vector<AlignStruct> &ass, unordered_map<ReadId, CorrectionRead> &reads )
{
    lens[0] = lens[1] = lens[2] = 0;
    ReadId longest;
    for ( auto it = ids.begin(); it != ids.end(); it++ )
    {
        auto r = reads.find( *it );
        if ( r == reads.end() ) it = ids.erase( it )-1;
        else if ( addLen( r->second.exts[drxn] ) ) longest = *it;
    }
    
    if ( !lens[0] ) return false;
    if ( lens[0] == seq.size() ) return true;
    
    seq = drxn ? seq.substr( 0, lens[0] ) : seq.substr( seq.size() - lens[0] );
    for ( AlignStruct &as : ass ) if ( &as != this && find( as.ids.begin(), as.ids.end(), longest ) != as.ids.end() ) return false;
    
    return true;
}

bool AlignStruct::setScore( Chunk &q, vector<Chunk> (&t)[3], unordered_map<ReadId, CorrectionRead> &reads, int gap )
{
    assert( q.drxn < 2 );
    int hits = 0, hitsLimits[2], hitCoords[2], tLimits[4] = { q.ends[0], q.ends[1], 1000, -1000 };
    hitsLimits[q.drxn] = hitCoords[q.drxn] = q.ends[q.drxn];
    hitsLimits[!q.drxn] = hitCoords[!q.drxn] = gap;
    hitIds.clear();
    for ( int i = 0; i < 3; i++ ) for ( Chunk &c : t[i] ) c.setLimits( &tLimits[ i < 2 ? 0 : 2 ] );
    vector<bool> tHits[2] = { vector<bool>( max( 0, tLimits[1] - tLimits[0] ), false ), vector<bool>( max( 0, tLimits[3] - tLimits[2] ), false ) };
    for ( int i = 0; i < 3; i++ ) setHits( t[i], tHits[ i == 2 ], tLimits[ i == 2 ? 2 : 0 ] );
    
    if ( q.drxn == 1 && q.caps[1] )
    {
        for ( Chunk &c : t[2] ) if ( c.caps[0] ) hits = max( hits, mapSeqOverlap( q.seq, c.seq, 24 ) );
        hitCoords[1] = max( q.ends[1] - hits, gap );
        if ( hits > seq.size() ) hitCoords[0] = q.ends[1];
    }
    for ( int i = 0; i < 2; i++ ) for ( bool hit : tHits[i] ) if ( hit ) hits++;
    for ( const pair<ReadId, pair<int,int> > &hit : q.hits )
    {
        if ( hit.second.first <= hitsLimits[0] ) hitCoords[0] = max( hitCoords[0], hit.second.second );
        if ( hitsLimits[1] <= hit.second.second ) hitCoords[1] = min( hitCoords[1], hit.second.first );
        if ( find( ids.begin(), ids.end(), hit.first ) != ids.end() ) continue;
        if ( q.drxn ? hit.second.first < gap : gap < hit.second.second ) addLen( seq.size() );
    }
    
    hitLen = abs( hitCoords[!q.drxn] - gap );
    if ( hitCoords[1] < hitCoords[0] ) return true;
    if ( hits > 50 ) return true;
    return false;
}

void IntervalTargets::addInterval( Interval* il, bool drxn )
{
    auto it = find( ils.begin(), ils.end(), il );
    assert( it != ils.end() && !il->edges[drxn] && il->good );
    il->edges[drxn] = new Interval( il->coords[drxn], false, il->pyro );
    il->edges[drxn]->edges[!drxn] = il;
    il->appendable[drxn] = false;
    ils.insert( it+drxn, il->edges[drxn] );
    added = true;
}

string& IntervalTargets::getRead( ReadId id )
{
    auto it = reads.find( id );
    if ( it != reads.end() ) return it->second;
    auto r = reads.insert( make_pair( id, bwt.getSequence( id ) ) );
    return r.first->second;
}

void IntervalTargets::query( vector<Chunk> &cs, unordered_set<ReadId> &ids, int cutoff, bool doCutoff )
{
    if ( cs.empty() || ids.empty() ) return;
    for ( const ReadId &id : ids )
    {
        string s = getRead( id );
        bool added = false;
        for ( Chunk &c : cs ) if ( c.match( id, s, cutoff, doCutoff ) ) added = true;
        if ( !added ) for ( Chunk &c : cs ) c.matchPartial( id, s );
    }
}

ReadId Interval::capCount = 0;
ReadId Interval::cappedCount = 0;
ReadId Interval::bridgeCount = 0;
ReadId Interval::bridgedCount = 0;

Interval::Interval( int coord, bool good, bool pyro )
: good( good ), bad( false ), pyro( pyro )
{
    coords[0] = coords[1] = coord;
    len = 0;
    edges[0] = edges[1] = NULL;
    extable[0] = extable[1] = appendable[0] = appendable[1] = true;
}

Interval::Interval( string &s, string &p, int first, int last, Interval* edge, bool good, bool pyro )
: Interval( first, good, pyro )
{
    coords[1] = last;
    len = last - first;
    seq = s.substr( min( first, last ), abs( len ) );
    phred = p.substr( min( first, last ), abs( len ) );
    edges[0] = edge;
    if ( edge ) edge->edges[1] = this;
}

bool Interval::clean( vector<Interval*> &ils )
{
    bool incomplete = false;
    for ( int i = 0; i < ils.size(); i++ )
    {
        if ( ils[i]->bad )
        {
            delete ils[i];
            ils.erase( ils.begin() + i-- );
        }
        else if ( !ils[i]->good )
        {
            if ( ils[i]->edges[0] && ( ils[i]->edges[0]->extable[1] || !ils[i]->edges[0]->alts[1].empty() ) ) incomplete = true;
            if ( ils[i]->edges[1] && ( ils[i]->edges[1]->extable[1] || !ils[i]->edges[1]->alts[1].empty() ) ) incomplete = true;
        }
    }
    return incomplete;
}

bool Interval::confirm( AlignStruct* best, bool initial, bool drxn )
{
    if ( !best ) return false;
    
    if ( !pyro )
    {
        if ( initial ) return best->errors.size() < 2;
        return ( !edges[drxn] || best->errors.size() < 3 ) && best->misses < 40;
    }
    int score = best->getNovelty( ( edges[drxn] && !edges[drxn]->alts[!drxn].empty() ? &edges[drxn]->alts[!drxn] : NULL ) );
    if ( initial ) return best->errors.size() < 2 + ( drxn && !edges[1] ) + score / 20;
    score += best->hitLen / ( 1 + (bool)edges[drxn] );
//    assert( best->errors.size() < 3 + ( score ) / 10 );
    if ( best->errors.size() < 3 + ( score ) / 10 ) return true;
    int errorCount = -1;
    score = 1;
    for ( pair<int,int> &error : best->errors ) if ( error.second < best->aligns[1] ) errorCount++;
    for ( int i = 1; i < best->aligns[1]; i++ ) if ( drxn ? best->seq[i] != best->seq[i-1] : best->seq.end()[-i-1] != best->seq.end()[-i] ) score++;
    
    return score > errorCount * 7;
}

void Interval::correct( Querier &bwt, vector<Interval*> (&ils)[4], bool pe, bool chim )
{
    bool initial = true, incomplete[2] = { !ils[0].empty(), !ils[1].empty() }, updated[2]  = { !ils[0].empty(), !ils[1].empty() };
    unordered_map<ReadId, string> reads;
    
    while ( ( incomplete[0] || incomplete[1] ) && ( updated[0] || updated[1] ) )
    {
        updated[0] = updated[1] = initial;
        for ( int i = 0; i < 2; i++ )
        {
            if ( !incomplete[i] ) continue;
            IntervalTargets its( ils[i], bwt, reads );
            if ( pe && !initial ) Interval::getTargets( its, ils[!i] );
            for ( int j = 0; j < ils[i].size(); j++ ) if ( ils[i][j]->correct( its, initial ) ) updated[i] = true;
            if ( its.added ) updated[i] = true;
            incomplete[i] = Interval::clean( ils[i] );
        }
        initial = false;
    }
}

bool Interval::correct( IntervalTargets &its, bool initial )
{
    for ( int drxn = 0; drxn < 2; drxn++ )
    {
        if ( good || bad || !edges[drxn] || !edges[drxn]->query( its, !drxn ) ) continue;
        if ( edges[drxn]->correct( its, initial, !drxn ) )
        {
            if ( !good && !bad ) correct( its, initial );
            return true;
        }
        if ( initial ) return false;
    }
    
    return false;
}

bool Interval::correct( IntervalTargets &its, bool initial, bool drxn )
{
    // Attempt unchallenged extension
    int polyLen = drxn ? seq.size() - 1 - seq.find_last_not_of( seq.back() ) : seq.find_first_not_of( seq[0] );
    AlignStruct* best = AlignStruct::getExtension( best, alts[drxn], edges[drxn]->len, polyLen, true );
    
    if ( !edges[drxn]->confirm( best, true, drxn ) ) best = NULL;
    if ( extend( its, best, 20, drxn ) ) return true;
    if ( initial ) return false;
    
    // Attempt evidenced extension
    its.t[0].clear();
    its.t[1].clear();
    unordered_set<ReadId> qIds[2], tIds( its.ids.begin(), its.ids.end() );
    vector<Chunk> q;
    if ( edges[drxn] && edges[drxn]->edges[drxn] ) edges[drxn]->edges[drxn]->getTargets( tIds, its.t[drxn], drxn, drxn );
    getTargets( tIds, its.t[!drxn], !drxn, !drxn );
    bool caps[2] = { isCap( !drxn, 0 ), isCap( drxn, 1 ) };
    for ( AlignStruct &as : alts[drxn] )
    {
        if ( !as.bridged && as.aligns[1] < as.lens[0] ) as.setExtLens( edges[drxn]->len, as.lens[0] );
        string s = getSequence( as, drxn );
        q.push_back( Chunk( s, coords, caps, as.lens[1], as.bridged && as.valid, drxn ) );
        for ( ReadId id : as.ids )
        {
            if ( params.setPairId( id, 0 ) ) qIds[0].insert( id );
            else if ( params.setPairId( id, 1 ) ) qIds[1].insert( id );
        }
    }
    its.query( q, tIds, coords[drxn], true );
    its.query( its.t[0], qIds[0] );
    its.query( its.t[1], qIds[1] );
    its.query( its.t[2], qIds[1] );
    
    bool anyMatched = false;
    for ( int i = 0; i < q.size(); i++ ) if ( alts[drxn][i].setScore( q[i], its.t, reads, coords[drxn] ) ) anyMatched = true;
    best = AlignStruct::getExtension( best, alts[drxn], edges[drxn]->len, polyLen, !anyMatched );
    if ( !edges[drxn]->confirm( best, !anyMatched, drxn ) ) best = NULL;
    
    return extend( its, best, 1, drxn );
}

bool Interval::correct( IntervalTargets &its, AlignStruct* best, bool drxn )
{
    int shrinkLen = best->getValidGap( len );
    
    if ( !shrinkLen ) return false;
    
    if ( edges[drxn]->len - abs( shrinkLen ) < 32 )
    {
        merge( its, drxn );
        return true;
    }
    
    Interval* il = edges[drxn];
    if ( shrinkLen > 0 ) drxn = !drxn;
    bool shrinkable = il->edges[drxn] ? il->edges[drxn]->extable[!drxn] : il->appendable[drxn];
    for ( AlignStruct &as : il->alts[drxn] ) if ( as.lens[2] ) shrinkable = false;
    
    if ( shrinkable )
    {
        if ( !il->edges[drxn] ) its.addInterval( il, drxn );
        return il->shrink( abs( shrinkLen ), true, false, drxn );
    }
    
    return false;
}

vector<Interval*> Interval::create( string &s, string &p, vector<bool> &marks, bool doSplit, bool isPyro )
{
    vector<Interval*> intervals;
    vector< pair<int,int> > ils;
    if ( s.empty() ) return intervals;
    if ( marks.empty() )
    {
        intervals.push_back( new Interval( s, p, 0, s.size(), NULL, false, isPyro ) );
        return intervals;
    }
    
    int firstGood = 0, lastEnd = 0;
    for ( int i = 0; i <= marks.size(); i++ )
    {
        if ( i < marks.size() && marks[i] ) continue;
        if ( firstGood < i )ils.push_back( make_pair( firstGood, i+31 ) );
        firstGood = i + 1;
    }
    
    for ( int i = 0; i+1 < ils.size(); i++ )
    {
        bool polymer = true;
        char c = s[ ils[i+1].first ];
        for ( int j = ils[i+1].first+1; j < ils[i].second; j++ ) if ( s[j] != c ) polymer = false;
        for ( int j = ils[i].second; j < ils[i+1].first; j++ ) if ( s[j] != c ) polymer = false;
        if ( !polymer ) continue;
        while ( ils[i].second > 0 && s[ ils[i].second-1 ] == c ) ils[i].second--;
        while ( ils[i+1].first < s.size() && s[ ils[i+1].first ] == c ) ils[i+1].first++;
        
        for ( int j = i; --j >= 0; ) ils[j].second = min( ils[j].second, ils[i].second );
        for ( int j = i+2; j < ils.size(); j++ ) ils[j].first = max( ils[j].first, ils[i+1].first );
    }
    for ( int i = 0; i < ils.size(); i++ ) if ( ils[i].second - ils[i].first < 32 ) ils.erase( ils.begin() + i-- );
    
    for ( int i = 0; doSplit && i+1 < ils.size(); i++ )
    {
        int ol = ils[i].second - ils[i+1].first + 1;
        if ( ol > 0 ) ils[i].second -= ol;
        if ( ol > 0 ) ils[i+1].first += ol;
    }
    for ( int i = 0; i < ils.size(); i++ ) if ( ils[i].second - ils[i].first < 32 ) ils.erase( ils.begin() + i-- );
    for ( int i = 0; i < ils.size(); i++ )
    {
        if ( ils[i].first ) intervals.push_back( new Interval( s, p, lastEnd, ils[i].first, ( intervals.empty() ? NULL : intervals.back() ), false, isPyro ) );
        intervals.push_back( new Interval( s, p, ils[i].first, ils[i].second, ( intervals.empty() ? NULL : intervals.back() ), true, isPyro ) );
        lastEnd = ils[i].second;
    }
    if ( intervals.empty()|| lastEnd < s.size() ) intervals.push_back( new Interval( s, p, lastEnd, s.size(), ( intervals.empty() ? NULL : intervals.back() ), false, isPyro ) );
    assert( !intervals.empty() && !intervals[0]->coords[0] && intervals.back()->coords[1] == s.size() );
    
    return intervals;
}

bool Interval::extend( IntervalTargets &its, AlignStruct* best, int minExt, bool drxn )
{
    if ( !best || ( best->aligns[1] < minExt && best->aligns[0] < edges[drxn]->len ) ) return false;
    if ( !best->valid ) return edges[drxn]->correct( its, best, drxn );
    
    bool bridged = best->bridged;
    int diff = best->aligns[1] - best->aligns[0];
    seq = getSequence( *best, drxn );
    phred = getPhred( *best, drxn );
    coords[drxn] = drxn ? coords[0] + seq.size() : coords[1] - seq.size();
    len = seq.size();
    assert( seq.size() == phred.size() && len == seq.size() && coords[1] - coords[0] == len );
    Interval* il = edges[drxn]->edges[drxn]; 
    if ( !best->bridged ) edges[drxn]->shrink( best->aligns[0], false, false, !drxn );
    if ( il && !best->bridged ) il->realign( !drxn );
    fold( best, drxn );
    if ( bridged )
    {
        string ext = drxn ? seq.substr( 0, len - il->len ) : seq.substr( il->len );
        for ( auto read : il->reads )
        {
            if ( !read.second.congruent( ext, !drxn ) ) continue;
            read.second.exts[!drxn] = max( read.second.exts[!drxn] - len + il->len, 0 );
            if ( read.second.exts[0] || read.second.exts[1] )
            {
                auto it = reads.insert( read );
                if ( !it.second ) it.first->second = read.second;
            }
            else if ( params.isReadPe( read.second.id ) ) usedIds.insert( read.first );
        }
        for ( AlignStruct &alt : il->alts[drxn] ) if ( alt.setReads( il->alts[drxn], reads ) ) alts[drxn].push_back( alt );
        edges[drxn]->bad = il->bad = true;
        edges[drxn] = il->edges[drxn];
        if ( edges[drxn] ) edges[drxn]->edges[!drxn] = this;
        Interval::bridgedCount++;
    }
    else if ( il ) Interval::bridgeCount++;
    else if ( edges[drxn] ) Interval::capCount++;
    else Interval::cappedCount++;
    
    if ( diff && edges[drxn] ) ( drxn ? edges[drxn] : this )->shift( diff );
    extable[drxn] = true;
    
    return true;
}

void Interval::fold( AlignStruct* best, bool drxn )
{
    for ( auto it = reads.begin(); it != reads.end(); )
    {
        bool doErase = it->second.exts[drxn] && find( best->ids.begin(), best->ids.end(), it->first ) == best->ids.end();
        it->second.exts[drxn] = max( 0 , it->second.exts[drxn] - ( best->bridged ? best->ends[1] : best->aligns[1] ) );
        if ( !doErase && !it->second.exts[0] && !it->second.exts[1] && params.isReadPe( it->second.id ) ) usedIds.insert( it->first );
        if ( doErase || !it->second.exts[!drxn] ) it = reads.erase( it );
        else it++;
    }
    for ( int i = 0; i < alts[!drxn].size(); i++ )
    {
        if ( !alts[!drxn][i].setReads( alts[!drxn], reads ) ) alts[!drxn].erase( alts[!drxn].begin() + i-- );
    }
    alts[drxn].clear();
}

bool Interval::getAlignTargets( string &s, string &p, bool drxn )
{
    s = drxn ? seq.back() : seq[0];
    p = drxn ? phred.back() : phred[0];
    Interval* il = edges[drxn];
    int count = 0;
    while ( il && count++ < 2 )
    {
        s += drxn ? il->seq : string( il->seq.rbegin(), il->seq.rend() );
        p += drxn ? il->phred : string( il->phred.rbegin(), il->phred.rend() );
        il = il->edges[drxn];
    }
    
    return count == 2;
}

string Interval::getPhred( AlignStruct &as, bool drxn )
{
    string parts[2];
    parts[0] = edges[drxn] ? as.getPhred( edges[drxn]->phred, as.aligns[1], drxn ) : string( as.aligns[1], '#');
    parts[1] = as.bridged ? edges[drxn]->edges[drxn]->phred : "";
    return drxn ? phred + parts[0] + parts[1] : parts[1] + parts[0] + phred;
}

string Interval::getSequence( AlignStruct &as, bool drxn )
{
    string parts[2];
    if ( as.bridged && !as.valid )
    {
        as.bridged = false;
        as.aligns[1] = as.lens[0];
    }
    parts[0] = as.getSeq( as.aligns[1], drxn );
    if ( as.bridged ) parts[1] = edges[drxn]->edges[drxn]->seq;
    return drxn ? seq + parts[0] + parts[1] : parts[1] + parts[0] + seq;
}

Interval* Interval::getSplit( vector<Interval*> &ils )
{
    if ( ils.empty() ) return NULL;
    if ( ils[0]->seq.empty() ) { delete ils[0]; ils.erase( ils.begin() ); }
    if ( ils.back()->seq.empty() ) { delete ils.back(); ils.erase( ils.end()-1 ); }
    if ( ils.size() == 1 && ils[0]->good ) return NULL;
    return ils[0];
}

void Interval::getTargets( IntervalTargets &its, vector<Interval*> &ils )
{
    if ( ils.empty() ) return;
    ils[0]->getTargets( its.ids, its.t[2], 2, 1 );
    Chunk::reverse( its.t[2], ils.back()->coords[1] );
    vector<ReadId> revIds;
    for ( ReadId id : its.ids ) revIds.push_back( id & 0x1 ? id-1 : id+1 );
    its.ids.clear();
    its.ids.insert( revIds.begin(), revIds.end() );
}

void Interval::getTargets( unordered_set<ReadId> &ids, vector<Chunk> &t, int drxns, bool drxn )
{
    if ( !good && !edges[drxn] ) return;
    if ( !good ) edges[drxn]->getTargets( ids, t, 2, drxn );
    bool added = false;
    for ( int d = drxns == 1; d < ( drxns ? 2 : 1 ); d++ )
    {
        if ( !edges[d] ) continue;
        for ( AlignStruct &as : alts[d] )
        {
            if ( !as.setExtLens( edges[d]->len, as.seq.size() ) ) continue;
            added = true;
            string s = getSequence( as, d );
            bool caps[2] = { isCap( !d, 0 ), isCap( d, 1 ) };
            if ( as.bridged ) caps[d] = edges[d]->edges[d]->isCap( false, d );
            t.push_back( Chunk( s, coords, caps, as.lens[1], as.bridged, d ) );
            for ( ReadId id : as.ids ) if ( params.setPairId( id, !drxn ) ) ids.insert( id );
        }
    }
    
    bool caps[2] = { !edges[0] || !edges[0]->edges[0], !edges[1] || !edges[1]->edges[1] };
    if ( !added ) t.push_back( Chunk( seq, coords, caps ) );
    for ( ReadId id : usedIds ) if ( params.setPairId( id, !drxn ) ) ids.insert( id );
    
    if ( edges[drxn] ) edges[drxn]->getTargets( ids, t, 2, drxn );
}

int Interval::getUniqueLen( int uniqueCount, bool drxn )
{
    assert( uniqueCount );
    for ( int i = 1; i < seq.size(); i++ )
    {
        if ( drxn ? seq.end()[-i-1] == seq.end()[-i] : seq[i] == seq[i-1] ) continue;
        if ( --uniqueCount < 1 ) return i;
    }
    return len;
}

bool Interval::isCap( bool isAlt, bool drxn )
{
    return !edges[drxn] || ( !edges[drxn]->edges[drxn] && ( isAlt || alts[drxn].empty() ) );
}

void Interval::merge( IntervalTargets &its, bool drxn )
{
    its.added = true;
    assert( edges[drxn] && !bad );
    extable[drxn] = true;
    Interval* il = edges[drxn];
    il->bad = true;
    usedIds.insert( il->usedIds.begin(), il->usedIds.end() );
    coords[drxn] = il->coords[drxn];
    len += il->len;
    if ( il->len > 0 )
    {
        seq = drxn ? seq + il->seq : il->seq + seq;
        phred = drxn ? phred + il->phred : il->phred + phred;
    }
    else if ( il->len < 0 )
    {
        assert( false );
        seq = drxn ? seq.substr( 0, len ) : seq.substr( seq.size()-len );
        phred = drxn ? phred.substr( 0, len ) : phred.substr( phred.size()-len );
    }
    if ( il->edges[drxn] ) il->edges[drxn]->edges[!drxn] = this;
    edges[drxn] = il->edges[drxn];
    if ( edges[drxn] && edges[drxn]->good == good ) merge( its, drxn );
    else for ( int d = 0; d < 2; d++ ) if ( edges[d] && edges[d]->good ) edges[d]->realign( !d );
}

void Interval::output( vector<Interval*> &ils, ofstream &ofsSeqs, string* lines, bool pyro )
{
    assert ( !ils.empty() || ( lines[1].empty() && lines[3].empty() ) );
    if ( !ils.empty() && pyro )  ils[0]->shift( -ils[0]->coords[0] );
    else if ( !ils.empty() ) Interval::shrink( ils, 95 );
    lines[1].clear();
    lines[2] = "+";
    lines[3].clear();
    
    int expLen = ils.empty() ? 0 : ils.back()->coords[1];
    assert( lines[0].find( '|' ) == lines[0].npos );
    for ( Interval* il : ils )
    {
        assert( pyro || il->len >= 0 );
        if ( il->len < 0 ) lines[1].erase( lines[1].end()-il->len, lines[1].end() );
        if ( il->len < 0 ) lines[3].erase( lines[3].end()-il->len, lines[3].end() );
        if ( il->len > 0 ) lines[1] += il->seq;
        if ( il->len > 0 ) lines[3] += il->phred;
        if ( pyro && il->good ) lines[0] += "|" + to_string( il->coords[0] ) + "-" + to_string( il->coords[1] );
        delete il;
    }
    ils.clear();
    assert( lines[1].size() == lines[3].size() && lines[1].size() == expLen );
    
    for ( int i = 0; i < 4; i++ ) ofsSeqs << lines[i] + "\n";
}

bool Interval::query( IntervalTargets &its, bool drxn )
{
    if ( !alts[drxn].empty() ) return true;
    if ( !extable[drxn] ) return false;
    extable[drxn] = false;
    CorrectionStruct cs = its.bwt.mapCorrection( seq, len, drxn );
    if ( len < 50 || ( cs.overabundant && len < params.readLen ) ) extable[!drxn] = false;
    
    if ( cs.error )
    {
        if ( len > 64 )
        {
            if ( !shrink( 16, true, false, drxn ) ) return false;
            cs = its.bwt.mapCorrection( seq, len, drxn );
            its.added = true;
        }
        else if ( pyro )
        {
            edges[drxn]->merge( its, !drxn );
            if ( edges[drxn]->edges[!drxn] ) return edges[drxn]->edges[!drxn]->query( its, drxn );
            return false;
        }
    }
    
    if ( cs.overabundant || cs.error ) return false;
    
    int validated = cs.validate( len, 32, drxn );
    int diff = len - validated;
    query( cs, drxn );
    
    if ( diff )
    {
        bool requery = true;
        for ( AlignStruct &as : alts[drxn] ) if ( as.errors.empty() && as.ids.size() > 4 && as.novelty > 20 && as.lens[0] - as.lens[2] < 10 ) requery = false;
        if ( requery ) alts[drxn].clear();
        if ( requery || !cs.exts[!drxn].empty() ) alts[!drxn].clear();
        if ( diff * 2 < validated )
        {
            if ( !edges[!drxn] ) its.addInterval( this, !drxn );
            if ( shrink( diff, true, false, !drxn ) )
            {
                extable[!drxn] = true;
                its.added = true;
                if ( requery ) return query( its, drxn );
            }
        }
        assert( good );
        
    }
    
    return !alts[drxn].empty();
}

void Interval::query( CorrectionStruct &cs, bool drxn )
{
    assert( alts[drxn].empty() );
    bool prequeried = !alts[!drxn].empty();
    unordered_set<ReadId> goodIds[2], badIds[2];
    for ( bool d : { 0, 1 } )
    {
        if ( !alts[d].empty() ) continue;
        string s, p;
        getAlignTargets( s, p, d );
        for ( CorrectionExt &ext : cs.exts[d] )
        {
            alts[d].emplace_back( AlignStruct( ext, s, p, d ) );
            for ( int i = 0; i < alts[d].size()-1; i++ )
            {
                if ( alts[d].back().isBetter( alts[d][i], true ) )
                {
                    for ( ReadId &id : alts[d][i].ids ) badIds[d].insert( id );
                    alts[d].erase( alts[d].begin() + i-- );
                }
                else if ( alts[d][i].isBetter( alts[d].back(), true ) )
                {
                    for ( ReadId &id : alts[d].back().ids ) badIds[d].insert( id );
                    alts[d].pop_back();
                    break;
                }
            }
        }
    }
    for ( bool d : { 0, 1 } ) for ( AlignStruct &as : alts[d] ) goodIds[d].insert( as.ids.begin(), as.ids.end() );
    
    for ( CorrectionRead &read : cs.reads )
    {
        if ( goodIds[0].find( read.id ) == goodIds[0].end() && badIds[0].find( read.id ) != badIds[0].end() ) continue;
        if ( goodIds[1].find( read.id ) == goodIds[1].end() && badIds[1].find( read.id ) != badIds[1].end() ) continue;
        reads.insert( make_pair( read.id, read ) );
    }
    
    for ( bool d : { 0, 1 } )
    {
        for ( int i = 0; i < alts[d].size(); i++ )
        {
            if ( !alts[d][i].setReads( alts[d], reads ) ) alts[d].erase( alts[d].begin() + i-- );
        }
    }
    
    if ( len < 50 || prequeried ) return;
    if ( !edges[!drxn] ) alts[!drxn].clear();
    bool bridged = false;
    for ( AlignStruct &as : alts[!drxn] )
    {
        if ( !as.setExtLens( edges[!drxn]->len, as.lens[2] ) ) continue;
        if ( as.aligns[0] == edges[!drxn]->len && as.errors.empty() && as.lens[2] >= edges[drxn]->len + 15 ) bridged = true;
    }
    if ( !bridged ) alts[!drxn].clear();
}

void Interval::realign( bool drxn )
{
    if ( alts[drxn].empty() ) return;
    string s, p;
    getAlignTargets( s, p, drxn );
    for ( AlignStruct &as : alts[drxn] ) as.align( s, p );
}

void Interval::shift( int diff )
{
    coords[0] += diff;
    coords[1] += diff;
    if ( edges[1] ) edges[1]->shift( diff );
}

void Interval::shrink( vector<Interval*> &ils, int maxLen )
{
    if ( ils.empty() ) return;
    int diff = ils.back()->coords[1] - ils[0]->coords[0] - maxLen;
    while ( diff > 0 && ( !ils.back()->good || ( ils[0]->good && ils.back()->coords[1] > maxLen ) ) )
    {
        ils.back()->shrink( ( ils.back()->good ? min( diff, ils.back()->coords[1] - maxLen ) : diff ), false, true, 1 );
        if ( ils.back()->bad ) Interval::clean( ils );
        diff = ils.back()->coords[1] - ils[0]->coords[0] - maxLen;
    }
    while ( diff > 0 && ( !ils[0]->good || ( ils.back()->good && ils[0]->coords[0] < 0 ) ) )
    {
        ils[0]->shrink( ( ils[0]->good ? min( diff, -ils[0]->coords[0] ) : diff ), false, true, 0 );
        if ( ils[0] ) Interval::clean( ils );
        diff = ils.back()->coords[1] - ils[0]->coords[0] - maxLen;
    }
    ils[0]->shift( -ils[0]->len );
    while ( ils.back()->coords[1] > maxLen )
    {
        ils.back()->shrink( ils.back()->coords[1] - maxLen, false, true, 1 );
        if ( ils.back()->bad ) Interval::clean( ils );
    }
}

bool Interval::shrink( int diff, bool exchange, bool excise, bool drxn )
{
    if ( exchange && !edges[drxn] ) assert( false );
    if ( exchange && !edges[drxn]->extable[!drxn] ) return false;
    
    if ( ( excise || !edges[0] || !edges[1] ) && diff >= len )
    {
        for ( int i = 0; i < 2; i ++ ) if ( edges[i] && !edges[!i] ) edges[i]->edges[!i] = NULL;
        diff = len;
        bad = true;
    }
    
    assert( diff <= len );
    if ( exchange )
    {
        Interval* il = edges[drxn];
        il->extable[!drxn] = false;

        assert( il && diff > 0 ); 
        il->seq = drxn ? seq.substr( seq.size() - diff ) + il->seq : il->seq + seq.substr( 0, diff );
        il->phred = drxn ? phred.substr( phred.size() - diff ) + il->phred : il->phred + phred.substr( 0, diff );
        il->coords[!drxn] += drxn ? -diff : diff;
        il->len += diff;
    }
    
    seq = drxn ? seq.substr( 0, seq.size()-diff ) : seq.substr( diff );
    phred = drxn ? phred.substr( 0, phred.size()-diff ) : phred.substr( diff );
    coords[drxn] += drxn ? -diff : diff;
    len -= diff;
    
    if ( exchange && good && edges[!drxn] && edges[!drxn]->edges[!drxn] )
    {
        edges[!drxn]->edges[!drxn]->realign( drxn );
    }
    
    return true;
}

void Interval::split( Querier &bwt, vector<Interval*> &ils, vector<Interval*> &splitIls )
{
    if ( !split( bwt, 1 ) ) return;
    
    splitIls.insert( splitIls.end(), ils.begin() + 1, ils.end() );
    ils.erase( ils.begin() + 1, ils.end() );
    if ( splitIls.size() > 1 && !splitIls[0]->good ){ delete splitIls[0]; splitIls.erase( splitIls.begin() ); }
    
    splitIls[0]->split( bwt, 0 );
}

bool Interval::split( Querier &bwt, bool d )
{
    CorrectionStruct cs = bwt.mapCorrection( seq, len, d );
    if ( cs.reads.empty() ) return false;
    
    string s[2] = { CorrectionExt::getCongruent( cs.exts[0], 0 ), CorrectionExt::getCongruent( cs.exts[1], 1 ) };
    if ( s[d].empty() ) return false;
    if ( edges[!d] ) s[!d].clear();
    
    seq = s[0] + seq + s[1];
    phred = string( s[0].size(), '#' ) + phred + string( s[1].size(), '#' );
    len += s[0].size() + s[1].size();
    coords[0] -= s[0].size();
    coords[1] += s[1].size();
    
    if ( edges[d] && edges[d]->edges[!d] == this ) edges[d]->edges[!d] = NULL;
    edges[d] = NULL;
    
    return true;
}
