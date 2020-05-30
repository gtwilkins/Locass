/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "deamplify.h"
#include "correct_amplicon.h"
#include "error.h"
#include "timer.h"
#include <cassert>
#include <iostream>
#include <string.h>

//Deamplify::Deamplify( string &ifn, string &ofn, string &fnTmp1, string &fnTmp2 )
//{
//    fns_[0] = ifn;
//    fns_[1] = ofn;
//    fns_[2] = fnTmp1;
//    fns_[3] = fnTmp2;
//}
//
//void Deamplify::convert( ifstream* pIfs[2] )
//{
//    double startTime = clock();
//    cout << "Reading in sequences from file... " << endl;
//    
//    FILE* ofpSeq = fopen( fns_[1].c_str(), "wb" ),* ofpStatus = fopen( fns_[2].c_str(), "wb" );
//    string seqs[2], phreds[2];
//    int p = 0, b = 0, bufSize = 65536;
//    uint64_t kmers[2];
//    uint32_t pairCount = 0;
//    uint8_t readLen = 0, bufStatus[bufSize]{0}, bufSeq[2][1024], status, lens[2];
//    fwrite( &pairCount, 4, 1, ofpSeq );
//    fwrite( &readLen, 1, 1, ofpSeq );
//    
//    while ( readFastq( *pIfs[0][i], seqs[0], phreds[0] ) )
//    {
//        if ( !readFastq( *pIfs[1][i], seqs[1], phreds[1] ) ) failure( "Unequal number of read pairs." );
//        pairCount++;
//        status = seqs[0].size() < 16 || seqs[1].size() < 16 ? 3 : 0;
//        for ( int j = 0; j < 2; j++ )
//        {
//            lens[j] = seqs[j].size();
//            if ( lens[j] > readLen ) readLen = lens[j];
//            kmers[j] = 0;
//            for ( int k = 0; !status && k < seqs[j].size(); k++ )
//            {
//                uint8_t c = 4;
//                if ( seqs[j][k] == 'A' || seqs[j][k] == 'a' ) c = 0;
//                else if ( seqs[j][k] == 'C' || seqs[j][k] == 'c' ) c = 1;
//                else if ( seqs[j][k] == 'G' || seqs[j][k] == 'g' ) c = 2;
//                else if ( seqs[j][k] == 'T' || seqs[j][k] == 't' ) c = 3;
//                else if ( k < 16 ) status = 3;
//                bufSeq[j][k] = c > 3 ? 255 : c * 63 + phreds[j][k] - 33;
//                if ( k < 16 ) kmers[j] = ( kmers[j] << 2 ) + c;
//            }
//        }
//        if ( p == bufSize )
//        {
//            fwrite( bufStatus, 1, p, ofpStatus );
//            memset( bufStatus, 0, p );
//            p = 0;
//        }
//        bufStatus[p] += status << ( ( 3 - b ) * 2 );
//        b = b == 3 ? 0 : b+1;
//        if ( b == 0 ) p++;
//
//        if ( status ) lens[0] = lens[1] = 0;
//        fwrite( &lens, 1, 2, ofpSeq );
//        if ( status ) continue;
//        int d = kmers[1] < kmers[0];
//        kmers[0] = ( kmers[d] << 32 ) + kmers[!d];
//        fwrite( &kmers[0], 8, 1, ofpSeq );
//        fwrite( bufSeq[d], 1, seqs[d].size(), ofpSeq );
//        fwrite( bufSeq[!d], 1, seqs[!d].size(), ofpSeq );
//    }
//    if ( readFastq( *pIfs[1][i], seqs[1], phreds[1] ) ) failure( "Error: unequal number of read pairs." );
//    if ( readLen > 256 ) failure( "Error: maximum read length is too high. Max supported is 256." );
//    
//    if ( p || b ) fwrite( bufStatus, 1, p + bool( b ), ofpStatus );
//    
//    fclose( ofpSeq );
//    fclose( ofpStatus );
//    ofpSeq = fopen( ofns_[1].c_str(), "rb+" );
//    fwrite( &pairCount, 4, 1, ofpSeq );
//    fwrite( &readLen, 1, 1, ofpSeq );
//    fclose( ofpSeq );
//    
//    cout << "Read " << to_string( pairCount ) << " read pairs in: " << getDuration( startTime ) << endl;
//}
//
//void Deamplify::parse( DamplifyReadStruct &drs )
//{
//    double startTime = clock();
//    cout << endl << "Checking read pairs for amplicons... " << endl;
//    
//    uint64_t pairCount[2]{0}, pos[2]{0}, kmer;
//    uint32_t pairs;
//    uint8_t status;
//    drs.set( pos[0], pairCount[0] );
//    bool full = drs.maxPair;
//    lineLen_ = 4 + 2 + ( drs.readLen * 2 );
//    
//    int iteration = 0;
//    while ( full )
//    {
//        iteration++;
//        kmers_ = new unordered_map<uint64_t, uint32_t>( 500 * 1000 * 1000 );
//        uint64_t kmerLimit = kmers_->bucket_count() / 4;
//        full = false;
//        
//        // Parse as many k-mers as possible
////        while ( drs.read( status, kmer ) )
////        {
//////            if ( status == 1 ) drs.edit( 2 );
//////            if ( status ) continue;
////            if ( status > 1 ) continue;
////            auto it = kmers_->find( kmer );
////            if ( it != kmers_->end() ) it->second++;
////            else if ( full ) continue;
////            else if ( kmers_->size() >= kmerLimit ) full = drs.set( pos[1], pairCount[1] );
////            else kmers_->insert( make_pair( kmer, 1 ) );
////            drs.edit( 1 );
////        }
////        
////        // Write out all multiple k-mers
////        FILE* ofpKmers = fopen( ofns_[3].c_str(), "wb" );
////        uint32_t kmersUnread = 0, kmersRead = 0, ampCount = 0;
////        for ( const pair<uint64_t, uint32_t> &kmerCount : *kmers_ )
////        {
////            if ( kmerCount.second < 2 ) continue;
////            fwrite( &kmerCount.first, 8, 1, ofpKmers );
////            fwrite( &kmerCount.second, 4, 1, ofpKmers );
////            kmersUnread++;
////            ampCount += kmerCount.second;
////        }
////        cout << to_string( ampCount ) << " of " << to_string( ( full ? pairCount[1] : drs.maxPair ) - pairCount[0] ) << " were candidate amplicons in iteration " << iteration << endl;
////        kmers_->clear();
////        fclose( ofpKmers );
//        
//        // Collate and de-amplify all read for multiple k-mers
//        uint64_t kmersUnread = 1451007, kmersRead = 0;
//        FILE* ifpKmers = fopen( ofns_[3].c_str(), "rb" );
//        bufSeqs_ = new uint8_t[maxSeqMem_]{0};
//        uint64_t res = 0;
//        drs.reset( pos[0], pairCount[0] );
//        while ( kmersRead++ < kmersUnread )
//        {
//            fread( &kmer, 8, 1, ifpKmers );
//            fread( &pairs, 4, 1, ifpKmers );
//            if ( ( res + pairs ) * lineLen_ > maxSeqMem_ ) parseKmerSet( drs, res, pos[0], pairCount[0] );
//            kmers_->insert( make_pair( kmer, res ) );
//            res += pairs;
//        }
//        parseKmerSet( drs, res, pos[0], pairCount[0] );
//        
//        delete kmers_;
//        delete bufSeqs_;
//        
//        
//        drs.reset( pos[1], pairCount[1] );
//    }
//}
//
//void  Deamplify::parseKmerSet( DamplifyReadStruct &drs, uint64_t &res, uint64_t setPos, uint32_t setPair )
//{
//    uint64_t kmer, lastKmer = -1;
//    uint32_t test, last;
//    uint8_t status;
//    drs.reset( setPos, setPair );
//    
//    while ( drs.read( status, kmer ) )
//    {
//        assert( status != 0 );
//        if ( status > 2 ) continue;
//        last = drs.curPair;
//        auto it = kmers_->find( kmer );
//        if ( it == kmers_->end() ) continue;
////        drs.edit( 2 );
//        uint64_t p = it->second * lineLen_;
//        memcpy( &test, &bufSeqs_[p], 4 );
//        assert( !test );
//        memcpy( &bufSeqs_[p], &drs.curPair, 4 );
//        memcpy( &bufSeqs_[p+4], drs.lens, 2 );
//        memcpy( &bufSeqs_[p+6], drs.bufSeq, drs.lens[0] + drs.lens[1] );
//        it->second++;
//    }
//    for ( uint64_t i = 0; i < res; i++ )
//    {
//        uint64_t p = i * lineLen_;
//        memcpy( &test, &bufSeqs_[p], 4 );
//        assert( test && test <= drs.maxPair );
//    }
//    
//    vector<Amplicon*> amps[2];
//    for ( uint64_t i = 0; i < res; i++ )
//    {
//        uint64_t p = i * lineLen_;
//        for ( uint64_t j = p+6; j < p+22; j++ ) kmer = ( kmer << 2 ) + ( bufSeqs_[j] / 63 );
//        for ( uint64_t j = p+6+bufSeqs_[p+4]; j < p+22+bufSeqs_[p+4]; j++ ) kmer = ( kmer << 2 ) + ( bufSeqs_[j] / 63 );
//        if ( kmer != lastKmer && !amps[0].empty() )
//        {
//            assert( amps[0].size() > 1 && amps[1].size() > 1 );
//            for ( string line : Pile::compile( amps ) );
////            for ( Amplicon* amp : amps[0] ) drs.edit( amp->id_, 3 );
//            for ( int i = 0; i < 2; i++ )
//            {
//                for ( Amplicon* amp : amps[i] ) delete amp;
//                amps[i].clear();
//            }
//        }
//        amps[0].push_back( new Amplicon( &bufSeqs_[p], NULL ) );
//        amps[1].push_back( new Amplicon( &bufSeqs_[p], amps[0].back() ) );
//        lastKmer = kmer;
//    }
//    
//    memset( bufSeqs_, 0, maxSeqMem_ );
//    res = 0;
//}
//
//bool Deamplify::readFastq( ifstream &ifs, string &seq, string &phred )
//{
//    if ( !ifs.good() || !getline( ifs, seq ) || !getline( ifs, seq ) || !getline( ifs, phred ) || !getline( ifs, phred ) ) return false;
//    return true;
//}
