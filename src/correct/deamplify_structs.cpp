/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "deamplify_structs.h"
#include <string.h>
#include <cassert>
#include <iostream>
#include "constants.h"

DeamplifyFiles::DeamplifyFiles( string &ifnSeq, string &ifnStatus )
{
    bufSize = 65536;
    bufStatus = new uint8_t[bufSize];
    ifpSeq = fopen( ifnSeq.c_str(), "rb" );
    ifpStatus = fopen( ifnStatus.c_str(), "rb+" );
    ifpEdit = fopen( ifnStatus.c_str(), "rb+" );
    curPair = 0;
    fread( &maxPair, 4, 1, ifpSeq );
    fread( &readLen, 1, 1, ifpSeq );
    fread( bufStatus, 1, readSize(), ifpStatus );
    lineLen = 4 + 2 + ( readLen * 2 );
    p = -1;
    i = 3;
    doRewrite = false;
    
    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < 4; j++ ) outTable[i][j] = j << ( ( 3-i ) * 2 );
        mask[i] = ~outTable[i][3];
        for ( int j = 0; j < 256; j++ ) inTable[i][j] = ( j >> ( ( 3-i ) * 2 ) ) & 0x3;
    }
    
    set( 0 );
}

DeamplifyFiles::~DeamplifyFiles()
{
    delete bufStatus;
    fclose( ifpSeq );
    fclose( ifpStatus );
}

void DeamplifyFiles::collate( FILE* ofpDump, string &fnKmers, uint64_t kmerSize, uint64_t maxMem, uint64_t &ampliconCount  )
{
    FILE* ifpKmers = fopen( fnKmers.c_str(), "rb" );
    while ( kmersRead < kmersUnread )
    {
        uint8_t* bufSeqs = new uint8_t[maxMem]{0};
        unordered_map<uint64_t, uint32_t>* kmers = new unordered_map<uint64_t, uint32_t> ( kmerSize );
        uint64_t seqCount = 0, added = 0, kmerLimit = kmers->bucket_count() / 4, kmer;
        uint32_t pairs, test;
        uint8_t status;
        while ( kmersRead < kmersUnread )
        {
            assert( fread( &kmer, 8, 1, ifpKmers ) == 1 );
            assert( fread( &pairs, 4, 1, ifpKmers ) == 1 );
            if ( ( seqCount + pairs ) * lineLen > maxMem || kmers->size() >= kmerLimit ) break;
            kmers->insert( make_pair( kmer, seqCount ) );
            seqCount += pairs;
            kmersRead++;
        }
        if ( kmersRead < kmersUnread ) fseek( ifpKmers, -12, SEEK_CUR );

        reset( 0 );
        while ( read( status, kmer ) )
        {
            if ( curPair < setPair[1] ? status : status != 1 ) continue;
            auto it = kmers->find( kmer );
            if ( it == kmers->end() ) continue;
            if ( curPair >= setPair[1] ) edit( 2 );
            uint64_t p = it->second * lineLen;
            memcpy( &test, &bufSeqs[p], 4 );
            assert( !test );
            memcpy( &bufSeqs[p], &curPair, 4 );
            memcpy( &bufSeqs[p+4], lens, 2 );
            memcpy( &bufSeqs[p+6], bufSeq, lens[0] + lens[1] );
            it->second++;
            added++;
        }
        
        uint8_t blank[1024];
        memcpy( &blank, bufSeqs, 1024 );
        assert( blank[0] || blank[1] || blank[2] || blank[3] || blank[4] );
        fwrite( bufSeqs, 1, seqCount*lineLen, ofpDump );
        delete kmers;
        delete bufSeqs;
        ampliconCount += seqCount;
    }
    fclose( ifpKmers );
}

void DeamplifyFiles::cycle( string &fnKmers, uint64_t kmerSize, bool &full )
{
    unordered_map<uint64_t, uint32_t>* kmers = new unordered_map<uint64_t, uint32_t> ( kmerSize );
    uint64_t kmerLimit = kmers->bucket_count() / 4, kmer;
    uint8_t status;
    full = false;
    set( 0 );
    while ( read( status, kmer ) )
    {
        if ( full && status == 1 ) edit( 2 );
        if ( status ) continue;
        auto it = kmers->find( kmer );
        if ( it != kmers->end() ) it->second++;
        else if ( full ) continue;
        else if ( kmers->size() < kmerLimit ) kmers->insert( make_pair( kmer, 1 ) );
        
        if ( full ) edit( 1 );
        else if ( kmers->size() >= kmerLimit ) full = set( 1 );
    }
    if ( !full ) set( 1 );
        
    // Write out all multiple k-mers
    FILE* ofpKmers = fopen( fnKmers.c_str(), "wb" );
    uint64_t ampCount = 0;
    kmersUnread = kmersRead = 0;
    for ( const pair<uint64_t, uint32_t> &kmerCount : *kmers )
    {
        if ( kmerCount.second < 2 ) continue;
        fwrite( &kmerCount.first, 8, 1, ofpKmers );
        fwrite( &kmerCount.second, 4, 1, ofpKmers );
        kmersUnread++;
        ampCount += kmerCount.second;
    }
    cout << to_string( ampCount ) << " of " << to_string( ( full ? setPair[1] : maxPair ) - setPair[0] ) << " were candidate amplicons in iteration." << endl;
    fclose( ofpKmers );
    delete kmers;
}

void DeamplifyFiles::edit( uint32_t id, uint8_t status )
{
    uint64_t pos = id / 4;
    int place = id % 4;
    uint8_t byte;
    
    fseek( ifpEdit, pos, SEEK_SET );
    assert( fread( &byte, 1, 1, ifpEdit ) == 1 );
    byte = ( byte & mask[place] ) | outTable[place][status];
    fseek( ifpEdit, -1, SEEK_CUR );
    assert( fwrite( &byte, 1, 1, ifpEdit ) == 1 );
}

void DeamplifyFiles::edit( uint8_t status )
{
    bufStatus[p] = ( bufStatus[p] & mask[i] ) | outTable[i][status];
    doRewrite = true;
}

void DeamplifyFiles::overwrite( string (&lines)[2][4] )
{
    uint8_t status;
    uint64_t kmer;
    if ( !read( status, kmer ) ) return;
    for ( int i = 0; i < 2; i++ )
    {
        if ( lines[i][1].find( 'N' ) == string::npos || lens[i] != lines[i][1].size() ) continue;
        for ( int j = 0; j < lines[i][1].size(); j++ ) lines[i][1][j] = intToChar[ bufSeq[ i ? j+lens[0] : j ] ];
    }
}

bool DeamplifyFiles::read( uint8_t &status, uint64_t &kmer )
{
    i = i == 3 ? 0 : i+1;
    if ( i == 0 ) p++;
    if ( p == bufSize || curPair == maxPair )
    {
        if ( doRewrite ) rewrite();
        doRewrite = false;
        if ( curPair == maxPair ) return false;
        fread( bufStatus, 1, readSize(), ifpStatus );
        p = 0;
    }
    
    status = inTable[i][ bufStatus[p] ];
    curPair++;
    if ( status == 3 ) return true;
    
    fread( &lens, 1, 2, ifpSeq );
    fread( &kmer, 8, 1, ifpSeq );
    fread( bufSeq, 1, lens[0] + lens[1], ifpSeq );
    assert( lens[0] > 15 && lens[1] > 15 && lens[0] <= readLen && lens[1] <= readLen );
    for ( int i = 0; i < 16; i++ )
    {
        uint8_t c[2] = { bufSeq[i], bufSeq[ i+lens[0] ] };
        assert( c[0] < 252 && c[1] < 252 );
    }
    return true;
}

int DeamplifyFiles::readSize()
{
    return min( (uint32_t)bufSize, 1 + ( maxPair-curPair-1 ) / 4 );
}

void DeamplifyFiles::reset( int j )
{
    curPair = setPair[j];
    i = ( setPair[j] - 1 ) % 4;
    setPair[j] /= 4;
    p = i == 3 ? -1 : 0;
    fseek( ifpStatus, setPair[j], SEEK_SET );
    fseek( ifpSeq, setPos[j], SEEK_SET );
    clearerr( ifpStatus );
    clearerr( ifpSeq );
    fread( bufStatus, 1, readSize(), ifpStatus );
}

void DeamplifyFiles::rewrite()
{
    uint8_t bufTest[bufSize];
    memcpy( bufTest, bufStatus, 65536 );
    long int dist = p + bool( i );
    fseek( ifpStatus, -dist, SEEK_CUR );
    fwrite( bufStatus, 1, dist, ifpStatus );
}

void DeamplifyFiles::rewind()
{
    fseek( ifpSeq, 0, SEEK_SET );
    fseek( ifpStatus, 0, SEEK_SET );
    fread( &maxPair, 4, 1, ifpSeq );
    fread( &readLen, 1, 1, ifpSeq );
    fread( bufStatus, 1, readSize(), ifpStatus );
    curPair = 0;
}

bool DeamplifyFiles::set( int j )
{
    setPair[j] = curPair;
    setPos[j] = ftell( ifpSeq );
    return true;
}
