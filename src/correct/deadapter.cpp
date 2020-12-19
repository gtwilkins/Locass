/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include "deadapter.h"
#include "shared_functions.h"
#include "local_alignment.h"
#include "error.h"
#include "timer.h"

Deadapter::Deadapter( string fn )
: fn_( fn )
{
    header_ = footer_ = 0;
    lineCount_ = 0;
    rewind();
    
    string line;
    bool fileTypeSet = false;
    while ( !fileTypeSet )
    {
        for ( int i = 0; i < footer_+1; i++ ) readLine( line, false );
        if ( !line.empty() && line[0] == '@' )
        {
            header_ = 1;
            footer_ = 2;
            readLine( line, false );
        }
        else if ( !line.empty() && line[0] == '>' )
        {
            header_ = 1;
            readLine( line, false );
        }
        if ( !isSequence( line ) )
        {
            cerr << "Unexpected characters on line: " << lineCount_ << ", of file " << fn_ << endl;
            exit( 1 );
        }
        fileTypeSet = !line.empty();
    }
    
    ifs_.close();
    
    cout << "Attempting to find adaptors in read file \"" << fn << "\"..." << endl << endl;
    setAdapters( 13 );
    
    if ( !adapter_[0].empty() && !adapter_[1].empty() )
    {
        cout << "Found terminal adaptors:" << endl;
        cout << "    " << adapter_[0] << " - Pair 1" << endl;
        cout << "    " << adapter_[1] << " - Pair 2" << endl;
    }
    else
    {
        cout << "Did not find terminal adaptors." << endl;
        adapter_[0].clear();
        adapter_[1].clear();
    }
    
    if ( !connectors_.empty() )
    {
        cout << "Found mate-connecting adaptor:" << endl;
        for ( string &adp : connectors_ ) cout << "    " << adp << endl;
    }
    else cout << "Did not find mate-connecting adaptors." << endl;
    cout << endl;
}

int Deadapter::align( string (&s)[2], string (&p)[2], int &hits, int &miss, int &diff, int trim )
{
    if ( s[0].empty() || s[1].empty() ) return 0;
    string a[2];
    LocalAlignment la( s[0], s[1], false, true );
    la.realign( a[0], a[1], false, true );
    
    for ( int d = 0; d < 2; d++ )
    {
        if ( trim != 2 && trim != d ) continue;
        int len = a[d].end()[-1] == '-';
        if ( !len ) continue;
        while ( len < a[d].size() && a[d].end()[-1-len] == '-' ) len++;
        a[0].erase( a[0].end() - len, a[0].end() );
        a[1].erase( a[1].end() - len, a[1].end() );
        break;
    }
    
    int j[2]{0};
    bool polys = false;
    char c = 'X';
    for ( int i = 0; i < a[0].size() && j[0] < s[0].size() && j[1] < s[1].size(); i++ )
    {
        for ( int d = 0; d < 2 && i < a[0].size(); d++ )
        {
            if ( a[d][i] != '-' ) continue;
            bool poly = c != 'X';
            int polyLen = 0;
            for ( ; a[d][i] == '-' && i < a[0].size(); i++ )
            {
                if ( c == 'N' ) c = a[!d][i];
                if ( a[!d][i] != 'N' && a[!d][i] != c ) poly = false;
                if ( poly ) polyLen++;
                else miss += 20;
                j[!d]++;
            }
            miss += ( max( 0, polyLen - 3 ) ) * 5 + ( polys ? 10 : 0 );
            if ( poly ) polys = true;
            d = -1;
        }
        if ( i == a[0].size() ) break;
        
        if ( a[0][i] != 'N' && a[1][i] != 'N' )
        {
            if ( a[0][i] == a[1][i] ) hits++;
            else miss += min( 20, min( p[0][ j[0] ], p[1][ j[1] ] ) - 33 );
            c = a[0][i] == a[1][i] ? a[0][i] : 'X';
        }
        else if ( a[0][i] != 'N' ) c = a[0][i];
        else if ( a[1][i] != 'N' ) c = a[1][i];
        if ( a[0][i] != a[1][i] ) diff += max( -20, min( 20, p[0][ j[0] ] - p[1][ j[1] ] ) );
        j[0]++;
        j[1]++;
    }
    
    assert( j[0] <= s[0].size() );
    
    return j[0];
}

bool Deadapter::addEnds( vector< pair<string, uint32_t> > &ends, string &s )
{
    for ( int i = 0; i <= ends.size(); i++ )
    {
        if ( i == ends.size() ) ends.push_back( make_pair( s, 1 ) );
        else if ( ends[i].first == s ) ends[i].second++;
        else if ( ends[i].first.size() >= s.size() ) continue;
        else ends.insert( ends.begin()+i, make_pair( s, 1 ) );
        break;
    }
}

bool Deadapter::consolidateEnd( vector< pair<string, uint32_t> > &ends, int &len, string &adaptor )
{
    vector< pair<string, uint32_t> > kmers;
    for ( pair<string, uint32_t> &a : ends )
    {
        if ( a.first.size() < len ) break;
        string q = a.first.substr( 0, len );
        for ( pair<string, uint32_t> &b : kmers )
        {
            if ( q.empty() || b.first.size() < len ) break;
            if ( q != b.first ) continue;
            b.second += a.second;
            q.clear();
        }
        if ( !q.empty() ) kmers.push_back( make_pair( q, a.second ) );
    }
    sort( kmers.begin(), kmers.end(), []( pair<string, uint32_t> &a, pair<string, uint32_t> &b ){ return a.second > b.second; } );
    
    if ( kmers.empty() ) return false;
    if ( kmers.size() > 1 && kmers[0].second < 50 ) return false;
    
    int kmersCount = max( 2, (int)kmers.size() );
    if ( kmers[0].second / 4 > kmers[1].second ) kmers.erase( kmers.begin()+1, kmers.end() );
    for ( int i = 1; i < kmers.size()-1; i++ ) if ( ( kmers[i].second / 2 ) > kmers[i+1].second ) kmers.erase( kmers.begin()+i+1, kmers.end() );
    if ( kmersCount == kmers.size() ) return false;
    
    for ( int i = 0; i < ends.size(); i++ )
    {
        bool good = false;
        for ( int j = 0; j < kmers.size(); j++ ) if ( !ends[i].first.find( kmers[j].first ) ) good = true;
        if ( !good ) ends.erase( ends.begin() + i-- );
    }
    
    if ( kmers.size() > 1 ) return true;
    len++;
    adaptor = kmers[0].first;
    
    if ( kmers[0].second > 50 ) consolidateEnd( ends, len, adaptor );
    
    return true;
}

bool Deadapter::getPair( string (&s)[2] )
{
    string line;
    for ( int i = 0; i < 2; i++ )
    {
        for ( int j = 0; j < header_; j++ ) if ( !readLine( line, true ) ) return false;
        if ( !readLine( s[i], true ) ) return false;
        if ( !isSequence( s[i] ) )
        {
            cerr << "Unexpected characters on line: " << lineCount_ << ", of file " << fn_ << endl;
            exit( 1 );
        }
        for ( int j = 0; j < footer_; j++ ) if ( !readLine( line, true ) ) return false;
    }
    if ( s[0].empty() && s[1].empty() ) return false;
    return true;
}

int Deadapter::getDifference( string (&s)[2], int len[2] )
{
    string a[2];
    LocalAlignment la( s[0], s[1], false, false );
    la.realign( a[0], a[1], false );
    int diff = 0;
    for ( int i = 0; i < a[0].size() && len[0] < s[0].size() && len[1] < s[1].size(); i++ )
    {
        if ( a[0][i] != a[1][i] ) diff++;
        if ( a[0][i] != '-' ) len[0]++;
        if ( a[1][i] != '-' ) len[1]++;
    }
    return diff;
}

bool Deadapter::isConnected( string &seq, string &phred )
{
    for ( string& q : connectors_ )
    {
        size_t it = seq.find( q );
        if ( it != string::npos ) return setBlank( seq, phred, it );
        if ( it = mapSeqOverlap( seq, q, 12 ) ) return setBlank( seq, phred, seq.size() - it );
        if ( it = mapSeqOverlap( q, seq, 13 ) ) return setBlank( seq, phred, 0 );
        if ( q.size() <= 13 ) continue;
        if ( ( it = seq.find( q.substr( 0, 13 ) ) ) != string::npos )
        {
            string qq[2]{ q.substr( 13 ), seq.substr( it+13 ) };
            int exact = 0, len[2]{0}, limit = 2 + min( qq[0].size(), qq[1].size() ) / 10;
            while ( exact < min( qq[0].size(), qq[1].size() ) && qq[0][exact] == qq[1][exact] ) exact++;
            if ( exact+13 >= q.size() / 2 || ( getDifference( qq, len ) < limit ) ) return setBlank( seq, phred, it );
        }
        if ( ( it = seq.find( q.substr( q.size()-13 ) ) ) != string::npos )
        {
            string qq[2]{ string( q.rbegin()+13, q.rend() ), string( seq.rbegin()+seq.size()-it, seq.rend() ) };
            int exact = 0, len[2]{0}, limit = 2 + min( qq[0].size(), qq[1].size() ) / 10;
            while ( exact < min( qq[0].size(), qq[1].size() ) && qq[0][exact] == qq[1][exact] ) exact++;
            if ( getDifference( qq, len ) < limit ) return setBlank( seq, phred, it-len[1] );
            if ( exact+13 >= q.size() / 2 ) return setBlank( seq, phred, it+13-( q.size()/2 ) );
        }
    }
    
    return false;
}

bool Deadapter::isOverlap( string &seq1, string &seq2, string &phred1, string &phred2, bool blank )
{
    if ( seq1.size() < 24 || seq2.size() < 24 || adapter_[0].empty() || adapter_[1].empty()  ) return false;
    
    // This whole pair is junk
    if ( adapter_[0].size() > 21 && mapSeqOverlap( adapter_[0], seq1, 21 ) ||
            adapter_[1].size() > 21 && mapSeqOverlap( adapter_[1], seq2, 21 ) )
    {
        setBlank( seq1, phred1, 0 );
        setBlank( seq2, phred2, 0 );
        return true;
    }
    
    // Simple reverse complement overlap
    if ( isOverlap2( seq1, seq2, phred1, phred2, 0, blank ) ) return true;
    if ( isOverlap2( seq2, seq1, phred2, phred1, 1, blank ) ) return true;
    
    // Adaptor clearly overlaps with end
    bool ols[2]{ isOverlap4( seq1, phred1, 0 ), isOverlap4( seq2, phred2, 1 ) };
    if ( !ols[0] && !ols[1] ) return false;
    
    if ( !ols[0] ) setBlank( seq1, phred1, 0 );
    if ( !ols[1] ) setBlank( seq2, phred2, 0 );
    
    return true;
}

bool Deadapter::isOverlap2( string &seq1, string &seq2, string &phred1, string &phred2, int i, bool blank )
{
    string q = seq1.substr( 0, 24 ), rev = revCompNew( seq2 );
    size_t it = rev.find( q );
    if ( it == string::npos ) return false;
    
    int altOl = seq2.size() - it, ol = 24, hits = 0, miss = 0, diff = 0;
    while ( ol < seq1.size() && it + ol < rev.size() && seq1[ol] == rev[it+ol] ) ol++;
    hits = ol;
    string s[2] = { seq1.substr( ol ), rev.substr( ol + it ) };
    string p[2] = { phred1.size() > ol ? phred1.substr( ol ) : "", phred1.size() > it + ol ? string( phred2.rbegin() + it + ol, phred2.rend() ) : "" };
    ol += align( s, p, hits, miss, diff, 1 );
    bool confirmed = hits == seq1.size() && !miss;
    if ( ol * 2 - hits + ( miss / 10 ) > seq1.size() ) return false;
    
    assert( ol <= seq1.size() );
    if ( !confirmed && miss <= 40 && hits >= ol * .95  ) confirmed = isOverlap3( seq1, phred1, ol, i );
    if ( !confirmed ) diff = -1;
    if ( !confirmed ) confirmed = isOverlap3( seq2, phred2, altOl, !i );
    if ( !confirmed ) return false;
//        cout << lines[i][1] << endl << rev.substr( it ) << endl << endl;
    setBlank( seq1, phred1, blank ? 0 : ol );
    setBlank( seq2, phred2, altOl );
    return true;
}

bool Deadapter::isOverlap3( string &seq, string &phred, int ol, int i )
{
    string s[2] = { seq.substr( ol ), adapter_[i] };
    string p[2] = { phred.substr( ol ), string( s[1].size(), char( 53 ) ) };
    
    int hits = 0, miss = 0, diff = 0, len = min( s[0].size(), s[1].size() );
    align( s, p, hits, miss, diff, 2 );
    
    return hits + ( hits > 3 ) + ( hits / 10 ) + ( max( 0, hits - 5 ) / 10 ) - ( miss / 10 ) >= len;
}

bool Deadapter::isOverlap4( string &seq, string &phred, int i )
{
    size_t it = seq.find( adapter_[i] );
    if ( it == string::npos && adapter_[i].size() > 21 )
    {
        int ol = mapSeqOverlap( seq, adapter_[i], 21 );
        if ( ol ) it = seq.size() - ol;
    }
    if ( it == string::npos ) return false;
    
    for ( ; it < seq.size(); it++ )
    {
        seq[it] = 'N';
        if ( it < phred.size() ) phred[it] = char( 35 );
    }
    return true;
}

bool Deadapter::readLine( string &s, bool suppress )
{
    lineCount_++;
    if ( getline( ifs_, s ) ) return true;
    if ( !suppress )
    {
        cerr << "Failed to read line: " << lineCount_ << ", from file: " << fn_ << endl;
        exit( 1 );
    }
    return false;
}

void Deadapter::rewind()
{
    ifs_.close();
    ifs_.open( fn_ );
    if ( !ifs_.good() )
    {
        cerr << "Failed to open file." << endl;
        exit(1);
    }
}

void Deadapter::setAdapters( uint32_t kmerLen )
{
    rewind();
    if ( kmerLen > 13 )
    {
        cerr << "Kmer length set too large in adaptor identification." << endl;
        exit(1);
    }
    
    uint32_t kmerCount = pow( 4, kmerLen );
    uint32_t* kmers = new uint32_t[kmerCount]{0};
    
    lineCount_ = 0;
    adapter_[0] = adapter_[1] = "";
    int pairCount = 0, overlapCount = 0;
    string seqs[2];
    uint32_t kmer, kmerMask = 0;
    for ( int i = 0; i < kmerLen; i++ ) kmerMask = ( kmerMask << 2 ) + 3;
    
    // Set end adaptors
    vector< pair<string, uint32_t> > ends[2];
    bool finished[2] = { false, false };
    int endLens[2] = { 12, 12 };
    while ( pairCount < 500000 && getPair( seqs ) )
    {
        pairCount++;
        string revs[2] = { revCompNew( seqs[0] ), revCompNew( seqs[1] ) };
        int ols[2] = { mapSeqOverlap( revs[1], seqs[0], ( seqs[0].size() + 2 ) / 2 )
                     , mapSeqOverlap( revs[0], seqs[1], ( seqs[1].size() + 2 ) / 2 ) };
        for ( int i = 0; i < 2 && ols[0] && ols[1]; i++ )
        {
            if ( finished[i] || ( seqs[i].size() - ols[i] ) < endLens[i] ) continue;
            string adp = seqs[i].substr( ols[i] );
            if ( !adapter_[i].empty() && adp.find( adp ) ) continue;
            addEnds( ends[i], adp );
            if ( ends[i].size() >= 300 ) finished[i] = !consolidateEnd( ends[i], endLens[i], adapter_[i] );
        }
        if ( ols[0] || ols[1] ) overlapCount++;
        if ( ols[0] || ols[1] ) continue;
        
        for ( int i = 0; i < 2; i++ )
        {
            int n = -1, cLen = 0;
            char c = 'N';
            for ( int j = 0; j < seqs[i].size(); j++ )
            {
                kmer <<= 2;
                if ( seqs[i][j] == 'A' );
                else if ( seqs[i][j] == 'C' ) kmer += 1;
                else if ( seqs[i][j] == 'G' ) kmer += 2;
                else if ( seqs[i][j] == 'T' ) kmer += 3;
                else n = j;
                cLen = seqs[i][j] == c ? cLen + 1 : 1;
                c = seqs[i][j];
                if ( cLen > 3 ) n = j;
                if ( j - n < kmerLen ) continue;
                kmers[kmer&kmerMask]++;
                assert( kmers[kmer&kmerMask] );
            }
        }
    }
    for ( int i = 0; i < 2; i++ ) consolidateEnd( ends[i], endLens[i], adapter_[i] );
    
    // Find the most common k-mer that is likely to be the start of the mate-pair adaptor
    vector<uint32_t> highKmers;
    int listLimit = 4, minCutoff = pairCount / 10, maxKmerCount = 0;
    for ( uint32_t i = 0; i < kmerCount; i++ )
    {
        maxKmerCount = max( maxKmerCount, (int)kmers[i] );
        if ( kmers[i] < minCutoff ) continue;
        if ( highKmers.size() == listLimit && kmers[i] <= kmers[ highKmers.back() ] ) continue;
        uint32_t q[2] = { i >> 2, ( i << 2 ) & kmerMask }, qBest[2]{0}, qTotal[2]{0}, qInc[2] = { (uint32_t)1 << ( ( kmerLen - 1 ) * 2 ), 1 };
        for ( int j = 0; j < 4; j++ )
        {
            for ( int d = 0; d < 2; d++ )
            {
                if ( kmers[ q[d] ] > qBest[d] ) qBest[d] = kmers[ q[d] ];
                qTotal[d] += kmers[ q[d] ];
                q[d] += qInc[d];
            }
        }
        
        // This heuristic should remove k-mers that are part of repetitive elements or not the start of the connector
        if ( qBest[0] > qTotal[0] / 2 )
        {
            cout << "A   " << qBest[0] << "   " << qTotal[0] << endl;
            continue;
        }
        
        // This heuristic should remove k-mers that precede the start of the connector
        if ( qBest[1] * 3 / 4 > kmers[i] )
        {
            cout << "B   " << qBest[1] << "   " << kmers[i] << endl;
            continue;
        }
        
        if ( highKmers.size() == listLimit ) highKmers.erase( highKmers.end()-1 );
        int j = 0;
        while ( j < highKmers.size() && kmers[i] < kmers[ highKmers[j] ] ) j++;
        highKmers.insert( highKmers.begin()+j, i );
    }
    
    vector<string> connectors;
    for ( int i = 0; i < highKmers.size(); i++ )
    {
        if ( i && kmers[ highKmers[0] ] > kmers[ highKmers[i] ]*2 ) continue;
        string kmerStr( kmerLen, 'N' );
        uint32_t c = highKmers[i];
        for ( int j = 0; j < kmerLen; j++ )
        {
            kmerStr[kmerLen-j-1] = !(c&0x3) ? 'A' : ( (c&0x3) == 1 ? 'C' : ( (c&0x3) == 2 ? 'G' : 'T' ) );
            c >>= 2;
        }
        connectors.push_back( kmerStr );
    }
    
    delete kmers;
    
    if ( !connectors.empty() ) setConnectors( connectors );
}

bool Deadapter::setBlank( string &seq, string &phred, int start )
{
    for ( int i = start; i < seq.size(); i++ )
    {
        seq[i] = 'N';
        if ( i < phred.size() ) phred[i] = char( 35 );
    }
    return true;
}

bool Deadapter::setConnectors( vector<string> &connectors )
{
    rewind();
    int pairCount = 0;
    string seqs[2];
    vector< vector< pair<string, uint32_t> > > ends( connectors.size() );
    vector<bool> finished( connectors.size(), false );
    vector<int> lens;
    bool complete = false;
    for ( int i = 0; i < connectors.size(); i++ ) lens.push_back( connectors[i].size()+1 );
    while ( !complete && pairCount < 1000000 && getPair( seqs ) )
    {
        pairCount++;
        for ( int i = 0; i < 2; i++ )
        {
            for ( int j = 0; j < connectors.size(); j++ )
            {
                if ( finished[j] ) continue;
                auto it = seqs[i].find( connectors[j] );
                if ( it == seqs[i].npos ) continue;
                string adp = seqs[i].substr( it );
                if ( adp.find( 'N' ) != adp.npos ) continue;
                addEnds( ends[j], adp );
                if ( ends[j].size() < 300 ) continue;
                finished[j] = !consolidateEnd( ends[j], lens[j], connectors[j] );
                complete = finished[j];
                for ( int k = 0; k < finished.size(); k++ ) if ( !finished[k] ) complete = false;
            }
        }
    }
    
    for ( int i = 0; i < connectors.size(); i++ ) if ( !finished[i] ) consolidateEnd( ends[i], lens[i], connectors[i] );
    
    if ( !connectors.empty() )
    {
        string rev = revCompNew( connectors[0] );
        connectors_.push_back( connectors[0] );
        bool doAdd = false;
        if ( rev != connectors[0] ) for ( int i = 1; i < connectors.size(); i++ ) if ( doAdd = ( rev == connectors[i] ) ) break;
        if ( doAdd ) connectors_.push_back( rev );
    }
    
    return !connectors_.empty();
}
